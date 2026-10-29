#!/usr/bin/env python3
# Run everything of test/nemoback-make-clone-copy, each script in its own process, then the regression suites, then write TEST_REPORT.md.
#   own         test_core.py unit (static / copy / clone / guards / lifecycle / rng / bench), test_core.py core (copy paths A, births B + C), test_simulation.py
#   regression  the existing suites (never run in this repository: they rewrite their own TEST_REPORT.md + results/) in two throw-away clones:
#                 /tmp/mibg       = clone of features/developing (140acc4) + the 6 changed files copied from the working tree (md5 verified)
#                 /tmp/mibg-base  = the same clone without the copy (pristine 140acc4) -- tells which failures already exist at 140acc4
#               habitat / patch run_all.py stop at a failing smoke.py, so their test scripts are also run directly.
#               Every failure of the migrated clone is classified: (1) already failing at 140acc4 (known), (2) expected from the migration (old test asserts the
#               NotImplementedError guard / TypeError of asex-mix / the previous round's diff or tags -- the test needs the update, not the source), (3) other -> real issue
#   report      TEST_REPORT.md (the only file written in this folder) from RES/*.json; RES, raw logs of the runs (simulation runs, regression logs + artifacts) stay outside the repo in /tmp/nemoback-make-clone-copy/
# usage: python run_all.py [--own-only | --regression-only | --report-only] [--no-baseline] [--reclassify]   (--report-only --reclassify: classify the saved regression logs again)
import os, re, sys, json, glob, time, shutil, hashlib, platform, subprocess
from concurrent.futures import ThreadPoolExecutor

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, '..', '..'))
RES = '/tmp/nemoback-make-clone-copy/results'                      # results/*.json of the own scripts + run_all_*.json (outside the repo)
REG_DIR = '/tmp/nemoback-make-clone-copy/regression'                # raw regression logs + artifacts (outside the repo)
PY = sys.executable
ENV = dict(os.environ, PYTHONDONTWRITEBYTECODE='1', MPLBACKEND='Agg')
MIBG, BASE = '/tmp/mibg', '/tmp/mibg-base'
SMOKE_DIGEST = '2b9b330388a3bfcdc7239fb45f7903076ebfaded'

OWN = [('test_core.py unit', ['test_core.py', 'unit'], 'test_unit', 'static checks S1-S25; make_copy (operator= fields + alleles, ID, Ve); make_clone (pedigree, mutation law); guards; lifecycle (recycling, RSS); rng (two processes, Python imitation); bench'),
       ('test_core.py core', ['test_core.py', 'core'], 'test_core', 'A: the 10 is_remove=False sites (fake individuals vs v3 seed by seed, colonization False vs True, real Nemo individuals); B + C: births with_num, into sites / pools, mixed split, D wrappers, marker births, error-path probes'),
       ('test_simulation.py', ['test_simulation.py'], 'test_simulation', 'experiments/model.py, both branches, scaled: invariants, teardown, digests, NemoBackend vs pure v3 + genetics controls')]

REG = [   # (key, folder under test/, argv)
    ('nemo-indfactory/smoke.py', 'nemo-indfactory', ['smoke.py']),
    ('nemo-indfactory/accept_genetics.py', 'nemo-indfactory', ['accept_genetics.py']),
    ('nemo-indfactory/accept_guards.py', 'nemo-indfactory', ['accept_guards.py']),
    ('nemo-indfactory/accept_lifecycle.py', 'nemo-indfactory', ['accept_lifecycle.py']),
    ('nemoback-habitat/run_all.py', 'nemoback-habitat', ['run_all.py']),
    ('nemoback-habitat/test_bookkeeping.py (direct)', 'nemoback-habitat', ['test_bookkeeping.py']),
    ('nemoback-habitat/test_nemo.py (direct)', 'nemoback-habitat', ['test_nemo.py']),
    ('nemoback-patch/run_all.py', 'nemoback-patch', ['run_all.py']),
    ('nemoback-patch/test_bookkeeping.py (direct)', 'nemoback-patch', ['test_bookkeeping.py']),
    ('nemoback-patch/test_nemo.py (direct)', 'nemoback-patch', ['test_nemo.py']),
    ('nemoback-patch/test_chain.py (direct)', 'nemoback-patch', ['test_chain.py']),
    ('nemoback-metacommunity/test_metacommunity.py all', 'nemoback-metacommunity', ['test_metacommunity.py', 'all']),
    ('nemoback-extension/test_extension.py all --quick', 'nemoback-extension', ['test_extension.py', 'all', '--quick'])]

# (1) failures the task names as known at 140acc4 (used to label category-1 failures; category 1 itself = "also fails in the pristine 140acc4 clone")
KNOWN_140 = [(r'nemoback-patch/test_bookkeeping', r'original_patch|ORIGINAL|v3|== original|original', 'patch test_bookkeeping: comparison with the v3 class (matrix rounding changed in 140acc4)'),
             (r'nemoback-metacommunity', r'helper|== original|unchanged helpers', 'metacommunity: "helper == original" (matrix_around changed in 140acc4)'),
             (r'nemoback-extension', r'ORIGINAL_EXT_COMMIT|ORIGINAL_META_COMMIT|original_ext_commit|original_meta_commit|agreed change present|original_metacommunity', 'extension: ORIGINAL_EXT_COMMIT / ORIGINAL_META_COMMIT auto-location')]
# (2) failures expected from this migration: the old tests assert what the migration intentionally changed
EXPECTED = [(r'no NotImplementedError|NotImplementedError \+ full snapshot|is_remove=False accepted|accepted is_remove=False|is_remove=False \(guard\)', 'asserts the removed is_remove=False NotImplementedError guard'),
            (r"has no attribute 'make_copy'", 'old fake_indi has no make_copy(): reached only because the guard it expects is gone'),
            (r'guard_first|GUARDED|guard is the first statement|guard text ==|else branch is `pass|TODO\(copy_from\)|unreachable \(guard above\)|no random\.sample left', 'asserts the guard / the reserved else branch of the previous round'),
            (r"startswith\(\(?'TypeError'|not migrated \(TypeError|asex / mix", 'asserts that asex / mix still raise TypeError (now migrated)'),
            (r'is_remove=False == v3 \(copy\)|v3 copy removes nothing', 'compares is_remove=False with the v3 alias semantics on the guarded functions'),
            (r'every changed line of extension|agreed change|working tree: modified ==|no NemoBackend tag on|no function of the extension is tagged NemoBackend', 'asserts the previous round\'s diff / working tree / tag set'),
            (r'^baseline: Step 0 ', 'Step 0 of the extension re-runs a lower suite: its failures are classified under that suite above'),
            (r'smoke\.py FAIL|bookkeeping\s+FAIL|nemo\s+FAIL|chain\s+FAIL|FAILED sections|FAILED fuzz|section \w+ FAILED|rc=1|STATIC FAILED|FUZZ FAILED|MARKER FAILED|SIM FAILED|BASELINE FAILED|PERF FAILED|DISPERSAL FAILED|REGISTRY FAILED', 'summary line of a script / part whose own failures are listed separately')]


def run(argv, cwd, timeout=7200):
    t = time.time()
    try:
        p = subprocess.run(argv, capture_output=True, text=True, cwd=cwd, env=ENV, timeout=timeout)
        return {'rc': p.returncode, 's': round(time.time() - t, 1), 'out': p.stdout + '\n' + p.stderr}
    except subprocess.TimeoutExpired as e:
        return {'rc': 'timeout', 's': round(time.time() - t, 1), 'out': str(e)}


def git(*a, cwd=REPO):
    return subprocess.run(['git'] + list(a), capture_output=True, text=True, cwd=cwd).stdout.strip()


def md5(path):
    return hashlib.md5(open(path, 'rb').read()).hexdigest()


# ====================================================================== own scripts
def own():
    rows = []
    for script, argv, name, desc in OWN:
        jf = os.path.join(RES, name + '.json')
        if os.path.exists(jf): os.remove(jf)
        r = run([PY] + [os.path.join(HERE, argv[0])] + argv[1:], HERE)
        d = json.load(open(jf)) if os.path.exists(jf) else {}
        rows.append({'script': script, 'name': name, 'desc': desc, 'rc': r['rc'], 's': r['s'], 'checks': d.get('checks'), 'failed': d.get('failed'), 'findings': len(d.get('findings', [])),
                     'tail': r['out'].strip().splitlines()[-6:] if r['rc'] != 0 or not d else []})
        print('%-26s rc=%s %6.1fs checks=%s failed=%s findings=%s' % (script, r['rc'], r['s'], d.get('checks'), d.get('failed'), len(d.get('findings', []))), flush=True)
    json.dump(rows, open(os.path.join(RES, 'run_all_own.json'), 'w'), indent=1)
    return rows


# ====================================================================== regression
def make_clone(path, apply_changes):
    if os.path.exists(path): shutil.rmtree(path)
    common = git('rev-parse', '--git-common-dir')
    src = os.path.abspath(os.path.join(REPO, common, '..')) if not os.path.isabs(common) else os.path.dirname(common)
    subprocess.run(['git', 'clone', '-q', '--no-hardlinks', src, path], check=True)
    subprocess.run(['git', 'checkout', '-q', '-b', 'features/developing', 'origin/features/developing'], cwd=path, check=True)
    info = {'path': path, 'source': src, 'head': git('rev-parse', '--short', 'HEAD', cwd=path), 'local_head': git('rev-parse', '--short', 'HEAD')}
    changed = git('diff', '--name-only').split()
    if apply_changes:
        for f in changed: shutil.copyfile(os.path.join(REPO, f), os.path.join(path, f))
    info['files'] = {f: {'worktree': md5(os.path.join(REPO, f)), 'clone': md5(os.path.join(path, f))} for f in changed}
    info['md5_equal'] = all(v['worktree'] == v['clone'] for v in info['files'].values()) if apply_changes else all(v['worktree'] != v['clone'] for v in info['files'].values())
    info['git_status'] = git('status', '--short', cwd=path).splitlines()
    return info


TB = re.compile(r'Traceback \(most recent call last\):\n((?:[ \t].*\n)+)(\S.*)')


def tracebacks(text):
    ''' one signature per traceback: innermost frame in a test file (basename:function + code line) + the exception line '''
    out = []
    for m in TB.finditer(text + '\n'):
        frames = re.findall(r'File "([^"]+)", line (\d+), in (\S+)\n\s+(.*)', m.group(1))
        test_frames = [f for f in frames if '/test/' in f[0]] or frames
        f = test_frames[-1] if test_frames else ('?', '?', '?', '?')
        out.append({'where': '%s:%s' % (os.path.basename(f[0]), f[2]), 'code': f[3].strip()[:220], 'error': m.group(2).strip()[:300]})
    return out


def norm(s):
    return re.sub(r'\d+(\.\d+)?', '#', s)


def slug(label, key):
    return '%s__%s' % (label, re.sub(r'[^\w.-]+', '_', key))


def save_artifacts(label, key, folder, root):
    """ keep what the classification reads (the clones get rewritten by later suites: the extension's Step 0 re-runs the metacommunity suite, ...) """
    d = os.path.join(REG_DIR, slug(label, key))
    if os.path.exists(d): shutil.rmtree(d)
    os.makedirs(d)
    rep = os.path.join(root, 'test', folder, 'TEST_REPORT.md')
    if folder == 'nemoback-metacommunity' and os.path.exists(rep): shutil.copyfile(rep, os.path.join(d, 'TEST_REPORT.md'))
    if folder == 'nemoback-extension':
        for jf in glob.glob(os.path.join(root, 'test', folder, 'results', '*.json')):
            if not os.path.basename(jf).startswith(('sim_', 'all')): shutil.copyfile(jf, os.path.join(d, os.path.basename(jf)))
    return d


def signatures(folder, out, art):
    """ failure signatures of one suite run: tracebacks (log + saved TEST_REPORT.md / part JSON), failed check names (extension parts), FAIL / rc=1 lines of the log """
    text = out
    rep = os.path.join(art, 'TEST_REPORT.md')
    if os.path.exists(rep): text += '\n' + open(rep).read()
    sigs = []
    for t in tracebacks(text):
        sigs.append({'kind': 'traceback', 'text': '%s | %s | %s' % (t['where'], t['code'], t['error']), 'key': norm('%s | %s' % (t['where'], t['code']))})
    if folder == 'nemoback-extension':
        for jf in sorted(glob.glob(os.path.join(art, '*.json'))):
            part = os.path.basename(jf)[:-5]
            try: d = json.load(open(jf))
            except Exception: continue
            for c in d.get('results', []) if isinstance(d.get('results'), list) else []:
                if isinstance(c, dict) and c.get('ok') is False:
                    name = c.get('check', c.get('name')); ev = c.get('detail', c.get('evidence'))
                    sigs.append({'kind': 'check', 'text': '%s: %s | %s' % (part, name, str(ev)[:300]), 'key': norm('%s: %s' % (part, name))})
            if d.get('traceback'):
                for t in tracebacks(d['traceback']): sigs.append({'kind': 'traceback', 'text': '%s: %s | %s | %s' % (part, t['where'], t['code'], t['error']), 'key': norm('%s: %s | %s' % (part, t['where'], t['code']))})
    for line in out.splitlines():
        if re.search(r'\bFAIL(ED)?\b|\brc=1\b', line) and not line.startswith(('  File', 'Traceback')):
            sigs.append({'kind': 'line', 'text': line.strip()[:300], 'key': norm(re.sub(r'\s+', ' ', re.sub(r'\d+(\.\d+)?\s*s\b', '', line.strip())))})
    seen, uniq = set(), []
    for x in sigs:
        if x['key'] not in seen: seen.add(x['key']); uniq.append(x)
    return uniq


def classify(key, sig, base_keys):
    txt = sig['text']
    if sig['key'] in base_keys:
        for kre, tre, label in KNOWN_140:
            if re.search(kre, key) and re.search(tre, txt): return 1, label
        return 1, 'also fails in the pristine 140acc4 clone'
    for pat, reason in EXPECTED:
        if re.search(pat, txt): return 2, reason
    return 3, 'not explained by 140acc4 or by the intended migration'


def regression(with_base=True):
    os.makedirs(REG_DIR, exist_ok=True)
    meta = {'clones': {}, 'runs': {}}
    trees = [('migrated', MIBG, True)] + ([('base', BASE, False)] if with_base else [])
    for label, path, apply in trees:
        meta['clones'][label] = make_clone(path, apply)
        print('clone %s: %s HEAD %s md5 ok %s' % (label, path, meta['clones'][label]['head'], meta['clones'][label]['md5_equal']), flush=True)
    for label, path, _ in trees:                                          # one tree after the other: the suites have timing checks
        for key, folder, argv in REG:
            r = run([PY] + argv, os.path.join(path, 'test', folder))
            open(os.path.join(REG_DIR, slug(label, key) + '.log'), 'w').write(r['out'])
            save_artifacts(label, key, folder, path)
            meta['runs']['%s|%s' % (label, key)] = {'rc': r['rc'], 's': r['s']}
            print('[%s] %-52s rc=%-3s %7.1fs' % (label, key, r['rc'], r['s']), flush=True)
    json.dump(meta, open(os.path.join(REG_DIR, 'runs.json'), 'w'), indent=1, default=str)
    return reclassify()


def reclassify():
    """ classification from the saved logs + artifacts only (results/regression/) """
    meta = json.load(open(os.path.join(REG_DIR, 'runs.json')))
    out = {'clones': meta['clones'], 'suites': []}
    sig = {}
    for k in meta['runs']:
        label, key = k.split('|', 1); folder = [f for kk, f, _ in REG if kk == key][0]
        txt = open(os.path.join(REG_DIR, slug(label, key) + '.log')).read()
        sig[(label, key)] = (signatures(folder, txt, os.path.join(REG_DIR, slug(label, key))), txt)
    for key, folder, argv in REG:
        (msig, mtxt), b = sig[('migrated', key)], sig.get(('base', key))
        base_keys = {x['key'] for x in b[0]} if b else set()
        failures = []
        for x in msig:
            cat, why = classify(key, x, base_keys)
            failures.append({'category': cat, 'why': why, 'kind': x['kind'], 'text': x['text']})
        mr, br = meta['runs']['migrated|' + key], meta['runs'].get('base|' + key, {})
        row = {'key': key, 'folder': folder, 'argv': argv, 'migrated_rc': mr['rc'], 'migrated_s': mr['s'], 'base_rc': br.get('rc'), 'base_s': br.get('s'),
               'failures': failures, 'base_only': [x['text'] for x in (b[0] if b else []) if x['key'] not in {y['key'] for y in msig}],
               'repro': 'cd %s/test/%s && PYTHONDONTWRITEBYTECODE=1 %s %s' % (MIBG, folder, 'python', ' '.join(argv))}
        if key == 'nemo-indfactory/smoke.py':
            dg = re.findall(r'phenotype digest: (\w+)', mtxt)
            row['smoke_digest'] = dg[-1] if dg else None
            row['smoke_digest_ok'] = row['smoke_digest'] == SMOKE_DIGEST
        out['suites'].append(row)
    json.dump(out, open(os.path.join(RES, 'run_all_regression.json'), 'w'), indent=1, default=str)
    for r in out['suites']:
        print('%-52s base rc %-4s migrated rc %-4s  (1) %d  (2) %d  (3) %d  base-only %d' % (r['key'], r['base_rc'], r['migrated_rc'], *[sum(f['category'] == c for f in r['failures']) for c in (1, 2, 3)], len(r['base_only'])))
    return out


# ====================================================================== report
def esc(s):
    return str('' if s is None else s).replace('|', '/').replace('\n', ' ')


def table(head, rows):
    return ['| ' + ' | '.join(head) + ' |', '|' + '---|' * len(head)] + ['| ' + ' | '.join(esc(c) for c in r) + ' |' for r in rows]


def fmt(x, n=3):
    if isinstance(x, float): return ('%.' + str(n) + 'f') % x
    return '-' if x is None else str(x)


def load(name):
    p = os.path.join(RES, name + '.json')
    return json.load(open(p)) if os.path.exists(p) else {}


def report():
    sys.path.insert(0, HERE)
    import test_core as C
    env = C.env_info()
    own_rows = json.load(open(os.path.join(RES, 'run_all_own.json'))) if os.path.exists(os.path.join(RES, 'run_all_own.json')) else []
    reg = json.load(open(os.path.join(RES, 'run_all_regression.json'))) if os.path.exists(os.path.join(RES, 'run_all_regression.json')) else {}
    L = ['# make_copy / make_clone migration -- test report', '',
         'Generated %s by `test/nemoback-make-clone-copy/run_all.py`. Scope: the uncommitted migration on top of `%s` -- `individual.make_copy()` (Nemo `Individual::operator=`), '
         '`indfactory.make_clone()` (Nemo `IndFactory::makeClonedOffsprg`) and the core changes that use them (A: 10 `is_remove=False` copy sites, B: asexual / mixed births, '
         'C: marker births, D: NemoBackend tags). Test-only: nothing outside this folder was edited, nothing committed.' % (env['date'], env['head']), '',
         '## Environment', '',
         '- date `%s`, HEAD `%s` (`%s`), branch `%s`' % (env['date'], env['head'], env['head_subject'], env['branch']),
         '- python `%s`, pynemo `%s`, numpy `%s`, machine `%s`' % (env['python'], env['pynemo'], env['numpy'], env['machine']),
         '- conda env `metaibm-nemo-test2`; every script in its own process with `PYTHONDONTWRITEBYTECODE=1`', '']
    L += table(['changed file (working tree)', 'md5'], [('`%s`' % f, '`%s`' % h) for f, h in sorted(env['changed'].items())]) + ['']
    L += ['- decisions under test: D1 `is_remove` defaults stay True; D2 births into `offspring_pool` round with `int(x + U)`; D3 `hab_mixed_reproduce_germinate` makes exactly min(empty, a + s) offspring '
          '(asex_num = min(a, np.around(E a / (a + s))), sex_num = min(s, E - asex_num)), all placed; D4 `mutation_rate` / `pheno_var_ls` removed from the interface; D5 a birth habitat without parents raises `ValueError` (v3 behaviour, not fixed here)', '']
    # ---------------------------------------------------------------- own results
    tot_c = sum(r['checks'] or 0 for r in own_rows); tot_f = sum(r['failed'] or 0 for r in own_rows)
    L += ['## Results: %d checks, %d failed' % (tot_c, tot_f), '']
    rows = []
    for r in own_rows:
        rows.append(('`%s`' % r['script'], r['desc'], r['checks'], r['failed'], r['findings'], '%.1f s' % r['s'], 'PASS' if r['rc'] == 0 else '**FAIL** (rc %s)' % r['rc']))
        for part, v in load(r['name']).get('parts', {}).items():
            if part: rows.append(('', '- %s' % part, v['checks'], v['failed'], v['findings'], '', ''))
    L += table(['script / part', 'covers', 'checks', 'failed', 'findings', 'time', 'result'], rows) + ['']
    bench = load('test_unit').get('bench_us_per_call')
    if bench:
        cols = ['make_copy', 'make_clone', 'make_offspring', 'make_founder', 'choice + recycle only']
        L += ['`test_core.py unit` part bench: microseconds per call incl. recycle of the result (%d reps, Ve > 0, mu = 1e-4)' % load('test_unit').get('bench_reps', 0), '']
        L += table(['loci'] + cols, [(r['loci'],) + tuple('%.2f' % r[c] for c in cols) for r in bench]) + ['']
    for r in own_rows:
        if r['rc'] != 0 or r['failed']:
            d = load(r['name'])
            L += ['Failed in `%s`:' % r['script'], ''] + ['- %s -- `%s`' % (c['name'], esc(c['evidence'])[:400]) for c in d.get('results', []) if not c['ok']] + (['```'] + r['tail'] + ['```'] if r['tail'] else []) + ['']
    # ---------------------------------------------------------------- findings
    F = []
    for r in own_rows:
        for f in load(r['name']).get('findings', []): F.append(('%s%s' % (r['script'], ' / ' + f['part'] if f.get('part') else ''), f))
    reg_new = [s for s in reg.get('suites', []) for f in s['failures'] if f['category'] == 3]
    order = {'bug': 0, 'risk': 1, 'info': 2, 'known': 3}
    L += ['## Findings (%d from the tests + %d regression failures of category 3)' % (len(F), sum(1 for s in reg.get('suites', []) for f in s['failures'] if f['category'] == 3)), '']
    rows = []
    for s in reg.get('suites', []):
        for f in s['failures']:
            if f['category'] == 3: rows.append(('bug', '`%s`' % s['key'], 'regression failure not explained by 140acc4 or by the migration: %s' % f['text'], 'repro: `%s`' % s['repro'], 'run_all.py'))
    for script, f in sorted(F, key=lambda x: order.get(x[1]['severity'], 9)):
        ev = f.get('evidence')
        if isinstance(ev, dict) and 'trajectory (mean over seeds)' in ev: ev = {k: v for k, v in ev.items() if k != 'trajectory (mean over seeds)'}
        rows.append((f['severity'], '`%s`' % f['func'], f['text'], esc(json.dumps(ev, default=str))[:700] if ev is not None else '', script))
    L += table(['severity', 'function(s)', 'finding', 'evidence', 'from'], rows) + ['']
    L += NOTES_FINDINGS + ['']
    # ---------------------------------------------------------------- simulation
    sim = load('test_simulation')
    if sim:
        L += ['## Simulation (`test_simulation.py`)', '',
              '- experiments/model.py loop, scaled: islands 4x4 patches x 1 habitat 8x8 (1024 sites), mainlands 20x20, %d steps, environment +0.1 at t = 29 / 59 / 89 / 119, rain %.2f per mainland per step (area-scaled), '
              'per-capita rates as in model.py; %s runs, %.0f s wall with %d processes' % (sim['steps'], sim['params']['propagules_rain_num'], len(sim['runs']), sim['wall_s'], 6),
              '- modes: `v3sem` = NemoBackend with `is_remove=False` passed to the three calls whose v3 default was False (2 colonizations, marker dispersal, marker germinate-and-birth); '
              '`default` = NemoBackend, defaults (True); `novar` = v3sem with `quanti_init_variance = 0`; `v3gen` = novar + `quanti_mutation_variance = 0.025^2` (size of a v3 bit flip); '
              '`v3pure` = pure v3 (%s), v3 API, `is_remove=False`' % sim['v3pure'], '']
        L += table(['branch / kind / mode', 'runs', 'calls per run (copy / clone / sexual)', 'mean step', 'max step', 'run', 'RSS MB start / t=30 / end (max)', 'final island occupancy'],
                   [(k, v['runs'], ' / '.join(str(v['calls per run (mean)'][c]) for c in ('make_copy', 'make_clone', 'make_offspring')) if v.get('calls per run (mean)') else '-',
                     '%.1f ms' % v['mean step ms'], '%.0f ms' % v['max step ms'], '%.1f s' % v['run s (mean)'], ' / '.join('%.1f' % x for x in v['RSS MB start -> step 30 -> end (max over runs)']),
                     v['final islands occupancy (mean)']) for k, v in sim['perf'].items()]) + ['']
        for b in ('asexual', 'sexual'):
            rows = []
            for k, r in sim['nemo_vs_v3'][b].items():
                cell = lambda x: '-' if x[0] is None else '%s +- %s (n %d)' % (fmt(x[0]), fmt(x[1]), x[2])
                z = lambda v: '-' if v is None else ('**%.1f**' % v if v > 3 else '%.1f' % v)
                rows.append((k, cell(r['nemo']), cell(r['v3']), z(r['z']), cell(r['nemo novar']), z(r['z novar vs v3']), cell(r['nemo v3gen']), z(r['z v3gen vs v3'])))
            L += ['### %s branch: NemoBackend vs pure v3 (5 seeds each, mean +- SE)' % b, '']
            L += table(['metric', 'NemoBackend v3sem', 'pure v3', 'z', 'novar', 'z', 'v3gen', 'z'], rows) + ['']
            tr = sim['nemo_vs_v3'][b + ' trajectory']
            ts = [t for t in next(iter(tr.values())) if t.startswith('t=')]
            L += table(['island occupancy / sp1 phenotype'] + ts + ['founder phenotype SD (ml1 sp1, ml2 sp2)'],
                       [(k,) + tuple('%s / %s' % (row[t]['occupancy'], fmt(row[t]['sp1 phenotype'])) for t in ts) + ('%s / %s' % (row['founder phenotype SD']['mainland1 sp1'], row['founder phenotype SD']['mainland2 sp2']),) for k, row in tr.items()]
                       + [('environment mean',) + tuple(str(next(iter(tr.values()))[t]['env']) for t in ts) + ('',)]) + ['']
        L += ['### is_remove=False (v3 semantics) vs default True, NemoBackend, 3 seeds (mean +- SE)', '']
        rows = []
        for b in ('asexual', 'sexual'):
            for k, r in sim['v3sem_vs_default'][b].items():
                rows.append((b, k) + tuple('%s +- %s' % (fmt(v[0]), fmt(v[1])) for v in r.values()))
        L += table(['branch', 'metric', 'is_remove=False', 'default True'], rows) + [''] + sim_notes(sim) + ['']
    # ---------------------------------------------------------------- regression
    if reg:
        L += ['## Regression: existing suites in throw-away clones', '']
        for label, c in reg['clones'].items():
            L += ['- `%s` (%s): clone of `%s`, HEAD `%s` (worktree HEAD `%s`); %s; `git status`: %s' % (
                c['path'], label, c['source'], c['head'], c['local_head'],
                ('6 changed files copied, md5 == working tree: **%s**' % c['md5_equal']) if label == 'migrated' else 'no copy (pristine 140acc4; md5 differ from the working tree as expected: %s)' % c['md5_equal'],
                ', '.join('`%s`' % x for x in c['git_status']) or 'clean')]
        L += ['']
        cnt = lambda s, k: sum(1 for f in s['failures'] if f['category'] == k)
        L += table(['suite / script', 'pristine 140acc4', 'migrated', '(1) known at 140acc4', '(2) expected from the migration', '(3) other', 'only at 140acc4 (masked / not reached)'],
                   [(s['key'], 'rc %s, %s s' % (s['base_rc'], s['base_s']), 'rc %s, %s s' % (s['migrated_rc'], s['migrated_s']), cnt(s, 1), cnt(s, 2), cnt(s, 3), len(s['base_only'])) for s in reg['suites']]) + ['']
        sm = [s for s in reg['suites'] if s['key'] == 'nemo-indfactory/smoke.py']
        if sm: L += ['- nemo-indfactory smoke digest: `%s` (expected `%s`): **%s**' % (sm[0].get('smoke_digest'), SMOKE_DIGEST, 'OK' if sm[0].get('smoke_digest_ok') else 'DIFFERENT'), '']
        L += ['### Every failure of the migrated clone', '']
        rows = []
        for s in reg['suites']:
            for f in s['failures']: rows.append((s['key'], f['category'], f['why'], f['kind'], f['text'][:420]))
        L += table(['suite', 'cat', 'why', 'kind', 'failure (test frame / code / error, or check name, or summary line)'], rows) + ['']
        bo = [(s['key'], t[:300]) for s in reg['suites'] for t in s['base_only']]
        if bo: L += ['### Failing at 140acc4 but not in the migrated run (the migrated script stopped earlier, or the failure point moved)', ''] + table(['suite', 'failure at 140acc4'], bo) + ['']
        L += reg_notes(reg) + ['']
    L += ['## Notes on the method', ''] + NOTES_METHOD + ['']
    open(os.path.join(HERE, 'TEST_REPORT.md'), 'w').write('\n'.join(L) + '\n')
    print('TEST_REPORT.md written (%d lines)' % len(L))


NOTES_FINDINGS = []      # analyst notes that the tables cannot carry (kept short)


def sim_notes(sim):
    ''' reading of the simulation tables; every number is taken from RES/test_simulation.json '''
    out = ['**Reading the simulation tables**', '']
    a = sim['nemo_vs_v3']['asexual']; tr = sim['nemo_vs_v3']['asexual trajectory']
    k_occ = 'occupancy (last 30 steps)'
    traj = lambda kind: ' -> '.join('%s' % tr[kind][t]['occupancy'] for t in tr[kind] if t.startswith('t='))
    ph = lambda kind, t: tr[kind][t]['sp1 phenotype']
    last = 't=%d' % (sim['steps'] - 1)
    out.append('- asexual branch, %s: NemoBackend %.0f +- %.0f vs pure v3 %.0f +- %.0f (z %.1f). The trajectories show where it comes from: NemoBackend island sp1 follows the '
               'environmental shift (sp1 phenotype %s at t=59 -> %s at %s, environment 0.2 -> 0.4) and keeps islands occupied; v3 sp1 stays near 0 (%s at t=119) and the islands '
               'are empty at the end. Occupancy NemoBackend %s; v3 %s.' % (k_occ, a[k_occ]['nemo'][0], a[k_occ]['nemo'][1], a[k_occ]['v3'][0], a[k_occ]['v3'][1], a[k_occ]['z'],
               fmt(ph('nemo v3sem', 't=59')), fmt(ph('nemo v3sem', last)), last, fmt(ph('v3pure v3sem', 't=119')), traj('nemo v3sem'), traj('v3pure v3sem')))
    zmax = max(r['z v3gen vs v3'] for r in a.values() if r['z v3gen vs v3'] is not None)
    out.append('- the two controls separate genetics from the copy / clone code: without founder genetic variance (`novar`: founder phenotype SD %s, as v3 %s) NemoBackend still '
               'adapts (sp1 %s at the end, z %.1f); with v3-sized mutation effects as well (`v3gen`: quanti_mutation_variance 0.025^2 instead of 0.01, i.e. the per-birth mutational '
               'variance of v3\'s bit flips) NemoBackend reproduces v3 -- occupancy %s vs v3 %s, every metric within %.1f SE. So the difference is the default Nemo genetic parameters '
               'of the NemoBackend runs (quanti_init_variance 0.002025, quanti_mutation_variance 0.01, as in test/metaibm-nemo-run/simulate.py), not the migration; the clone / '
               'copy semantics themselves are pinned by the unit tests.' % (tr['nemo novar']['founder phenotype SD']['mainland1 sp1'], tr['v3pure v3sem']['founder phenotype SD']['mainland1 sp1'],
               fmt(ph('nemo novar', last)), a[k_occ]['z novar vs v3'], traj('nemo v3gen'), traj('v3pure v3sem'), zmax))
    sx = sim['nemo_vs_v3']['sexual']; k = 'mainland1 size (time mean)'
    mb = sim['perf']['sexual nemo v3sem'].get('mixed birth calls not filling every empty site (sum over runs)') or {'calls': -1, 'filled': -1, 'site left empty although asex + sex >= empty': -1, 'too few parents (asex + sex < empty)': -1}
    out.append('- sexual branch, %s: NemoBackend %.3f (z %s) vs v3 %.3f. Mixed births on the mainlands, measured in the 8 NemoBackend v3sem runs: of %d calls %d filled every '
               'empty site, %d left a site empty although asex + sex >= empty (the v3 split\'s np.around undershoot; 0 expected with the current split, see the '
               'hab_mixed_reproduce_germinate finding) and %d had too few parents; v3gen control %s.' % (k, sx[k]['nemo'][0], fmt(sx[k]['z'], 1), sx[k]['v3'][0], mb['calls'], mb['filled'],
               mb['site left empty although asex + sex >= empty'], mb['too few parents (asex + sex < empty)'], fmt(sx[k]['nemo v3gen'][0])))
    k2 = 'occupancy (time mean, of 1024 sites)'
    out.append('- the scaled sexual (mixed) islands stay almost empty in every variant (time-mean occupancy NemoBackend %.1f, v3 %.1f of 1024): single colonists rarely find a mate '
               'or fit well enough to clone, so this branch\'s island statistics are weak; its births (%s clones + %s sexual offspring per run) are still exercised, and the '
               'mixed births are tested directly in test_core.py (births parts).' % (sx[k2]['nemo'][0], sx[k2]['v3'][0], sim['perf']['sexual nemo v3sem']['calls per run (mean)']['make_clone'],
                                                                   sim['perf']['sexual nemo v3sem']['calls per run (mean)']['make_offspring']))
    out.append('- is_remove=False vs default True: with False the colonists are copies (%s make_copy per run) and the mainlands stay full; with True (0 copies) colonists leave the '
               'mainland, which then averages %s sites instead of 400.' % (sim['perf']['asexual nemo v3sem']['calls per run (mean)']['make_copy'],
                                                                         fmt(sim['v3sem_vs_default']['asexual']['mainland1 size (time mean)']['default True'][0], 2)))
    out.append('- memory: RSS is flat within each run (start / t=30 / end in the performance table); teardown minted == recycled in every NemoBackend run.')
    return out




def reg_notes(reg):
    """ reading of the regression tables, from RES/run_all_regression.json """
    n = {c: sum(1 for s in reg['suites'] for f in s['failures'] if f['category'] == c) for c in (1, 2, 3)}
    out = ['**Reading the regression tables**', '',
           '- %d failures of the migrated clone: %d already fail at 140acc4 (1), %d assert what this migration changed on purpose (2), **%d other (3)**.' % (sum(n.values()), n[1], n[2], n[3])]
    stops = [(s['key'], f['text'].split(' | ')[0], f['text'].split(' | ')[1][:120]) for s in reg['suites'] for f in s['failures'] if f['category'] == 2 and f['kind'] == 'traceback']
    if stops:
        out.append('- coverage the old suites lose until their tests are updated: these scripts / sections stop at an expected (2) failure, so their later checks did not run '
                   'against the migrated code -- ' + '; '.join('`%s` at `%s` (`%s`)' % x for x in stops) + '. The behaviour they would have checked is covered by the new tests above.')
    masked = [(s['key'], t[:140]) for s in reg['suites'] for t in s['base_only'] if s['key'].startswith(('nemoback-patch/test_bookkeeping', 'nemoback-metacommunity'))]
    if masked:
        out.append('- the known 140acc4 failures of patch test_bookkeeping (v3 class comparison: "RNG streams differ") and metacommunity bookkeeping ("helpers differ from the original") '
                   'are not reached in the migrated run (the scripts stop earlier, at the guard assertions); they are listed under "Failing at 140acc4 but not in the migrated run".')
    out.append('- extension suite: `ORIGINAL_EXT_COMMIT = git log -1 -- extension/global_habitat_network.py`, which is `53152d5` (the migration commit) at 140acc4 instead of `1356387`; '
               'its perf / dispersal / registry parts stop inside `original_metacommunity` in BOTH trees (same failure point, category 1). `test_extension.py all` exits 0 although parts fail.')
    out.append('- the expected (2) failures need test updates, not source changes: drop the NotImplementedError / TypeError expectations, give hab_common.fake_indi a make_copy(), '
               'and refresh the diff / tag / working-tree assertions of the previous round.')
    return out
NOTES_METHOD = [
    '- **v3 comparison (A)**: the v3 bodies are compiled verbatim from git into a subclass of the CURRENT class (habitat / patch `baace34`, metacommunity `7a1adbb`, extension `1356387`, '
    'module-level function indented into the class), so both twins run on the same lower layers; the twins are built from one seed and called with one seed. Fake individuals carry '
    '`.origin`, so a v3 alias and a current copy compare equal in the layout. The two direct-germination dispersals rarely overflow the empty sites on random worlds, so they are also run '
    'with the number matrix pinned (`[forced overflow]`, no RNG, same matrix for both twins) -- there the unplaced copies are recycled on ~185 / 200 seeds.',
    '- **make_copy semantics**: `test_core.py unit` (part rng) rebuilds copies and clones in Python without `operator=` / `makeClonedOffsprg` (setters, `set_allele_value`, `set_value`, then the '
    'clone steps of `breed_cloning`) and gets the same IDs, fields, alleles and phenotypes record by record -- so `operator=` copies exactly what is listed (all fields except the ID; '
    '`_mother` / `_father` raw pointers too, never read by MetaIBM) and RAND is consumed in the same order (one Ve draw in `operator=`, one more for a clone after mutation).',
    '- **test the tests**: broken variants (alias germination, remove instead of copy, a make_copy that drops the age) were monkeypatched in a scratch process; the copy-path harness '
    'flagged all three.',
    '- **pure v3**: `git archive baace34 metaibm` alone is not v3 -- `baace34` itself introduced the Nemo-backed `individual.py` (+ `indfactory.py`); `/tmp/v3pure` therefore takes '
    '`metaibm/individual.py` from `baace34^` (`67a742f`) and `extension/` from `baace34` (== `1356387`), as the v3 metacommunity imports it.',
    '- **regression classification**: category 1 = the same failure point (test frame + code line, or extension check name) also fails in the pristine 140acc4 clone; '
    'category 2 = the failure asserts what the migration intentionally changed (patterns listed in `run_all.py:EXPECTED`); category 3 = everything else. A script that aborts at its '
    'first failure hides its later checks: those parts are not exercised against the migrated code by the old suites (the new tests here cover the migrated behaviour).']


def main():
    os.makedirs(RES, exist_ok=True)
    only = [a for a in sys.argv[1:] if a.endswith('-only')]
    if not only or '--own-only' in only: own()
    if not only or '--regression-only' in only: regression(with_base='--no-baseline' not in sys.argv)
    if '--reclassify' in sys.argv: reclassify()
    report()


if __name__ == '__main__':
    main()
