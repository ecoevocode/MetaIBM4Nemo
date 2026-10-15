# Step 0 baseline (test/nemo-indfactory/*.py + test/nemoback-habitat/run_all.py, unchanged; the habitat folder is backed up and restored because its run_all rewrites
# TEST_REPORT.md + results/), then smoke.py (stop if it fails), then every test in its own process; method coverage; call-chain matrix; write TEST_REPORT.md.
#   python run_all.py [--perf] [--no-baseline]
import os, re, ast, sys, json, glob, time, shutil, hashlib, inspect, platform, tempfile, subprocess, datetime

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, '..', '..'))
TESTS = [('test_bookkeeping.py', 'six clear-ups, disturbance, is_remove=False guards, equivalence with the original (bodies on the current habitat + full original stack), forced / natural extras, sampling uniformity, edge shapes (fake individuals, no Nemo)'),
         ('test_nemo.py', 'real Nemo individuals, one fresh process per section: initialize / reproduction / two_species / lifecycle / getters / memory / reproducibility / fuzz / edges_known'),
         ('test_chain.py', 'call chain patch -> habitat -> indfactory -> individual (method wrapping), ordering rule, offspring homepatch, static scan of callers')]
NOISE = ('setting random seed from input value', 'outputs: ', 'the microsite in the habitat is occupied.', 'the microsite in the habitat is empty.')
ENV = dict(os.environ, PYTHONDONTWRITEBYTECODE='1')
BASE = [None]


def run(script, cwd=HERE, keep=25):
    t = time.time()
    p = subprocess.run([sys.executable, os.path.join(cwd, script)], capture_output=True, text=True, cwd=cwd, env=ENV)
    lines = [l for l in (p.stdout + '\n' + p.stderr).splitlines() if l.strip() and not any(n in l for n in NOISE) and not l.startswith('[FINDING')]
    return {'rc': p.returncode, 's': time.time() - t, 'tail': lines[-keep:]}


def load(name):
    try: return json.load(open(os.path.join(HERE, 'results', name + '.json')))
    except Exception: return {}


def git(*a):
    return subprocess.run(['git'] + list(a), capture_output=True, text=True, cwd=REPO).stdout.rstrip()


def md5(rel):
    return hashlib.md5(open(os.path.join(REPO, rel), 'rb').read()).hexdigest()


# ------------------------------------------------------------------ Step 0 baseline
def baseline():
    out = {'indfactory': {}, 'habitat': None}
    d = os.path.join(REPO, 'test', 'nemo-indfactory')
    for s in sorted(os.path.basename(x) for x in glob.glob(os.path.join(d, '*.py')) if not x.endswith('common.py')):
        r = run(s, cwd=d, keep=3); out['indfactory'][s] = (r['rc'], round(r['s'], 1))
    hd = os.path.join(REPO, 'test', 'nemoback-habitat')
    tmp = tempfile.mkdtemp(prefix='hab_suite_backup_')
    shutil.copytree(hd, os.path.join(tmp, 'suite'), copy_function=shutil.copy2)
    try:
        r = run('run_all.py', cwd=hd, keep=12)
        rep = open(os.path.join(hd, 'TEST_REPORT.md')).read()
        os.makedirs(os.path.join(HERE, 'results'), exist_ok=True)
        open(os.path.join(HERE, 'results', 'baseline_habitat_TEST_REPORT.md'), 'w').write(rep)
        fails = re.findall(r'\| (\S+\.py) \| .*? \| \*\*FAIL\*\* \|', rep)
        tb = re.search(r'test_bookkeeping.py</summary>\n\n```\n(.*?)```', rep, re.S)
        out['habitat'] = {'rc': r['rc'], 's': round(r['s'], 1), 'lines': [l for l in r['tail'] if re.match(r'(smoke|test_\S*)\.py', l)], 'failed': fails,
                          'bookkeeping_tail': (tb.group(1).strip().splitlines()[-4:] if tb else [])}
    finally:
        shutil.rmtree(hd); shutil.copytree(os.path.join(tmp, 'suite'), hd, copy_function=shutil.copy2); shutil.rmtree(tmp)
    out['habitat_restored_clean'] = git('status', '--short', '--', 'test/nemoback-habitat') == ''
    hab_dirty = git('status', '--short', '--', 'metaibm/habitat.py') != ''
    out['habitat_py_committed'] = not hab_dirty
    out['habitat_py_last_commit'] = git('log', '-1', '--format=%h %s', '--', 'metaibm/habitat.py')
    return out


# ------------------------------------------------------------------ method coverage
def _strip_doc(fn):
    fn = ast.parse(ast.unparse(fn)).body[0]
    if fn.body and isinstance(fn.body[0], ast.Expr) and isinstance(getattr(fn.body[0], 'value', None), ast.Constant) and isinstance(fn.body[0].value.value, str): fn.body = fn.body[1:] or [ast.Pass()]
    return ast.dump(fn)


def coverage(orig_commit, chain):
    src = open(os.path.join(REPO, 'metaibm', 'patch.py')).read(); osrc = git('show', orig_commit + ':metaibm/patch.py')
    cls = [n for n in ast.parse(src).body if isinstance(n, ast.ClassDef) and n.name == 'patch'][0]
    ocls = [n for n in ast.parse(osrc).body if isinstance(n, ast.ClassDef) and n.name == 'patch'][0]
    lines, olines = src.splitlines(), osrc.splitlines()
    orig = {n.name: n for n in ocls.body if isinstance(n, ast.FunctionDef)}
    texts = {os.path.basename(p): open(p).read() for p in glob.glob(os.path.join(HERE, '*.py')) if not p.endswith(('run_all.py', 'common.py'))}
    internal = {}
    for root, row in chain.get('matrix', {}).items():
        for m in row.get('patch (internal)', []): internal.setdefault(m, []).append(root)
    rows, missing = [], []
    cur = [n for n in cls.body if isinstance(n, ast.FunctionDef)]
    for n in cur:
        m = n.name
        if m not in orig: status = 'new'
        elif ast.dump(orig[m]) == ast.dump(n) and '\n'.join(lines[n.lineno - 1:n.end_lineno]) == '\n'.join(olines[orig[m].lineno - 1:orig[m].end_lineno]): status = 'unchanged'
        elif _strip_doc(orig[m]) == _strip_doc(n): status = 'docstring / comment only'
        else: status = 'changed'
        scope = 'asex/mix (probe only)' if re.search(r'asex|mix', m) else 'in scope'
        hits = sorted(f for f, t in texts.items() if re.search(r'\b%s\b' % re.escape(m), t))
        via = internal.get(m, [])
        via_s = ('%d patch methods, e.g. %s' % (len(via), ', '.join(sorted(via)[:2]))) if len(via) > 2 else ', '.join(sorted(via))
        if m == '__init__': via_s = 'constructor'
        rows.append((m, status, scope, ', '.join(h.replace('.py', '') for h in hits) or '-', via_s))
        if scope == 'in scope' and not hits and not via and m != '__init__': missing.append(m)
    for m in orig:
        if m not in {n.name for n in cur}: rows.append((m, 'removed', 'in scope', '-', 'replaced by the six patch_clear_up_* (callers: metacommunity.py, next round)'))
    return rows, missing


# ------------------------------------------------------------------ perf (record only)
def perf():
    sys.path.insert(0, HERE)
    import random, warnings
    import numpy as np
    warnings.filterwarnings('ignore')
    from common import make_patch, recycle_everything, nemo_species, seed_all
    seed_all(91)
    F = {'sp1': nemo_species('sp1', 0.2, seed=91), 'sp2': nemo_species('sp2', 0.5, seed=92)}
    T = [[0.2, 0.2], [0.5, 0.5]]
    out = {}
    for L in (10, 20, 40):
        np.random.seed(L); p = make_patch(4, L=L, K=L * L // 2, gidx0=0, means=[(0.2, 0.2), (0.5, 0.5)] * 2)
        p.patch_initialize(F, T)
        calls = [('patch_initialize', lambda: p.patch_initialize(F, T)), ('patch_dead_selection', lambda: p.patch_dead_selection(0.2, 0.5, 'niche_gaussian')),
                 ('patch_sex_reproduce_mutate_into_offspring_pool', lambda: p.patch_sex_reproduce_mutate_into_offspring_pool(1.0, None, None)),
                 ('patch_dispersal_within_from_offspring_pool_to_immigrant_pool', lambda: p.patch_dispersal_within_from_offspring_pool_to_immigrant_pool(0.3)),
                 ('patch_local_germinate_from_offspring_and_immigrant_pool', p.patch_local_germinate_from_offspring_and_immigrant_pool),
                 ('patch_dormancy_process_from_offspring_pool_and_immigrant_pool', p.patch_dormancy_process_from_offspring_pool_and_immigrant_pool),
                 ('patch_dead_selection ', lambda: p.patch_dead_selection(0.2, 0.5, 'niche_gaussian')),
                 ('patch_sex_reproduce_mutate_into_offspring_pool ', lambda: p.patch_sex_reproduce_mutate_into_offspring_pool(1.0, None, None)),
                 ('patch_dispersal_within_from_offspring_pool_and_dormancy_pool_to_immigrant_pool', lambda: p.patch_dispersal_within_from_offspring_pool_and_dormancy_pool_to_immigrant_pool(0.3)),
                 ('patch_dipersal_within_from_offspring_and_dormancy_pool', lambda: p.patch_dipersal_within_from_offspring_and_dormancy_pool(0.3)),
                 ('patch_local_germinate_from_offspring_immigrant_and_dormancy_pool', p.patch_local_germinate_from_offspring_immigrant_and_dormancy_pool),
                 ('patch_dormancy_process_from_offspring_pool_to_dormancy_pool', p.patch_dormancy_process_from_offspring_pool_to_dormancy_pool),
                 ('patch_sexual_birth_germinate', lambda: p.patch_sexual_birth_germinate(1.0, None, None)),
                 ('patch_sex_reproduce_calculation_into_offspring_marker_pool', lambda: p.patch_sex_reproduce_calculation_into_offspring_marker_pool(1.0)),
                 ('patch_dispersal_within_from_offspring_marker_pool_to_immigrant_marker_pool', lambda: p.patch_dispersal_within_from_offspring_marker_pool_to_immigrant_marker_pool(0.3, True)),
                 ('get_patch_microsites_individals_phenotype_values', lambda: p.get_patch_microsites_individals_phenotype_values('pheno1')),
                 ('get_patch_microsites_individals_sp_id_values', p.get_patch_microsites_individals_sp_id_values),
                 ('patch_disturbance_process', p.patch_disturbance_process), ('patch_clear_up_dormancy_pool', p.patch_clear_up_dormancy_pool)]
        tm = {}
        for rep in range(3):
            for name, fn in calls:
                t = time.perf_counter(); fn(); tm[name.strip()] = tm.get(name.strip(), 0) + (time.perf_counter() - t) / 3
            p.patch_initialize(F, T)
        out['4 habitats x %dx%d' % (L, L)] = {k: round(v, 4) for k, v in tm.items()}
        recycle_everything(p)
    return out


def main():
    for f in glob.glob(os.path.join(HERE, 'results', '*.json')): os.remove(f)
    sys.path.insert(0, HERE)
    try:
        import importlib.metadata as md; pv = md.version('pynemo')
    except Exception: pv = '?'
    orig_commit = subprocess.run([sys.executable, '-c', 'import common; print(common.ORIGINAL_COMMIT)'], capture_output=True, text=True, cwd=HERE, env=ENV).stdout.split()[-1]
    files = ['metaibm/patch.py', 'metaibm/habitat.py', 'metaibm/indfactory.py', 'metaibm/individual.py']
    L = ['# patch.py test report (NemoBackend, sexual)', '',
         'Generated %s by `test/nemoback-patch/run_all.py`. Scope: `metaibm/patch.py` and its call chain patch -> habitat -> indfactory -> individual (metacommunity layer only scanned statically; asex / mix only probed).' % datetime.datetime.now().strftime('%Y-%m-%d %H:%M'), '',
         '## Environment', '', '- python `%s`, pynemo `%s`, numpy `%s`, `%s`' % (sys.version.split()[0], pv, __import__('numpy').__version__, platform.platform()),
         '- HEAD `%s` (`%s`), branch `%s`' % (git('rev-parse', '--short', 'HEAD'), git('log', '-1', '--format=%s'), git('rev-parse', '--abbrev-ref', 'HEAD')), '']
    L += ['| file | md5 | state |', '|---|---|---|'] + ['| `%s` | `%s` | %s |' % (f, md5(f), 'modified (working tree)' if git('status', '--short', '--', f) else 'committed') for f in files]
    L += ['', '- original reference: `git show %s:metaibm/patch.py` (blob `%s`; verified to still contain `hab_clear_up_offspring_and_immigrant_pool`). `%s` did not touch patch.py itself: '
          'the last commit that changed it is `%s`, so `%s:patch.py` == `%s:patch.py` == `HEAD:patch.py`. Full original stack = patch.py + habitat.py from `%s`.' %
          (orig_commit, git('rev-parse', '--short', orig_commit + ':metaibm/patch.py'), orig_commit, git('log', '-1', '--format=%h', '--', 'metaibm/patch.py'), orig_commit,
           git('log', '-1', '--format=%h', '--', 'metaibm/patch.py'), orig_commit), '']
    # ---- Step 0
    if '--no-baseline' not in sys.argv:
        t = time.time(); b = baseline(); BASE[0] = b; print('baseline %.1fs' % (time.time() - t))
        L += ['## Step 0 -- baseline (existing suites, run unchanged)', '', '| suite | script | result | time |', '|---|---|---|---|']
        L += ['| test/nemo-indfactory | %s | %s | %.1f s |' % (s, 'PASS' if rc == 0 else '**FAIL**', sec) for s, (rc, sec) in b['indfactory'].items()]
        h = b['habitat']
        for l in h['lines']:
            m = re.match(r'(\S+\.py)\s+(PASS|FAIL)\s+([\d.]+)s', l)
            if m: L.append('| test/nemoback-habitat (run_all.py) | %s | %s | %s s |' % (m.group(1), 'PASS' if m.group(2) == 'PASS' else '**FAIL**', m.group(3)))
        hab_orig = re.search(r"^ORIGINAL_COMMIT = '([^']+)'", open(os.path.join(REPO, 'test', 'nemoback-habitat', 'common.py')).read(), re.M).group(1)
        L += ['', '- habitat suite folder backed up before and restored after its run (its run_all.py rewrites TEST_REPORT.md + results/): restored clean = `%s`; its fresh report is kept as `results/baseline_habitat_TEST_REPORT.md`.' % b['habitat_restored_clean'],
              '- habitat suite original reference: `ORIGINAL_COMMIT = \'%s\'`; `metaibm/habitat.py` is %s (last commit `%s`).' % (hab_orig, 'committed' if b['habitat_py_committed'] else 'modified', b['habitat_py_last_commit'])]
        if hab_orig == 'HEAD' and b['habitat_py_committed']:
            L += ['- `HEAD` now holds the **migrated** habitat.py, so that suite\'s "== original" checks compare the new code with itself (dormancy "eliminate diff {}"; coverage table all "unchanged") and `test_bookkeeping.py` stops at '
                  '`assert o.habitat_disturbance_process() is None`. The suite is stale, not the code; fix: set its `ORIGINAL_COMMIT` to `%s` (= c3711cd^).' % orig_commit]
        if h['bookkeeping_tail'] and h['failed']: L += ['', '```'] + h['bookkeeping_tail'] + ['```']
        L.append('')
    # ---- smoke
    smoke = run('smoke.py'); sj = load('smoke')
    print('smoke.py %s %.1fs' % ('PASS' if smoke['rc'] == 0 else 'FAIL', smoke['s']))
    L += ['## Smoke: %s' % ('PASS' if smoke['rc'] == 0 else '**FAIL**'), '', '| step | return | invariants |', '|---|---|---|'] + ['| %s | `%s` | %s |' % (s['step'], s['return'], s['invariants']) for s in sj.get('steps', [])] + ['']
    if smoke['rc'] != 0:
        L += ['Failed step `%s`:' % sj.get('failed_step'), '```', sj.get('traceback', '\n'.join(smoke['tail'])), '```', 'Other tests not run.']
        open(os.path.join(HERE, 'TEST_REPORT.md'), 'w').write('\n'.join(L) + '\n'); sys.exit(1)
    res = {}
    for f, _ in TESTS:
        res[f] = run(f, keep=60); print('%-22s %s %.1fs' % (f, 'PASS' if res[f]['rc'] == 0 else 'FAIL', res[f]['s']))
    chain = load('chain')
    rows, missing = coverage(orig_commit, chain)
    pf = None
    if '--perf' in sys.argv:
        t = time.time(); pf = perf(); print('perf (record only) %.1fs' % (time.time() - t))
    n_fail = sum(r['rc'] != 0 for r in res.values()) + bool(missing)
    L += ['## Tests: %d / %d PASS' % (sum(r['rc'] == 0 for r in res.values()), len(res)), '', '| script | covers | result | time |', '|---|---|---|---|']
    L += ['| %s | %s | %s | %.1f s |' % (f, d, 'PASS' if res[f]['rc'] == 0 else '**FAIL**', res[f]['s']) for f, d in TESTS]
    L += ['', '## Key numbers', '']
    bk, nm = load('bookkeeping'), load('nemo')
    kn = []
    if bk:
        kn += ['- equivalence with `%s` (fake individuals): %s' % (orig_commit, bk.get('equivalence_cases')),
               '- OD->sites extras recycled == dropped by the original: %s over the sweep; forced extras %s' % (bk.get('OD->sites_extras_recycled_total'), bk.get('forced_extras')),
               '- natural extras (P = 0.25): %s' % bk.get('natural_extras'), '- uniformity chi2 (loose bound df + 6 sqrt(2 df)): %s' % bk.get('uniformity_chi2')]
    if nm:
        mem = (nm.get('memory') or {})
        kn += ['- memory: RSS %s MB, address reuse %s, teardown (minted, recycled) %s' % (mem.get('rss_mb'), mem.get('address_reuse'), mem.get('teardown')),
               '- fuzz: %s operations' % (nm.get('fuzz') or {}).get('operations'), '- reproducibility: %s' % (nm.get('reproducibility') or {}).get('seed7')]
    if chain: kn += ['- chain: %s recycle() calls checked for the ordering rule, %s make_offspring calls checked for homepatch' % (chain.get('recycle_calls'), chain.get('make_offspring_calls'))]
    L += kn + ['']
    for f, _ in TESTS:
        L.append('<details><summary>%s output</summary>' % f); L += ['', '```'] + res[f]['tail'] + ['```', '</details>']
    if pf:
        names = list(next(iter(pf.values())).keys())
        L += ['', '## Performance (record only, `--perf`): seconds per call, mean of 3', '', '| patch function | ' + ' | '.join(pf) + ' |', '|---|' + '---|' * len(pf)]
        L += ['| `%s` | %s |' % (n, ' | '.join('%.4f' % pf[k][n] for k in pf)) for n in names]
    fl = []
    for n in ('bookkeeping', 'nemo', 'chain'):
        for x in load(n).get('findings', []): fl.append(x)
    if BASE[0] and BASE[0]['habitat_py_committed'] and 'test_bookkeeping.py' in BASE[0]['habitat']['failed']:
        fl.insert(0, {'severity': 'info', 'text': 'Step 0: test/nemoback-habitat is stale since habitat.py was committed (%s): its ORIGINAL_COMMIT = "HEAD" is now the migrated habitat.py, so its "== original" checks compare the code with itself and test_bookkeeping.py fails at "assert o.habitat_disturbance_process() is None". Not a code defect; that suite was not edited.' % BASE[0]['habitat_py_last_commit'],
                      'repro': 'cd test/nemoback-habitat && python run_all.py', 'fix': 'set ORIGINAL_COMMIT = "baace34" (c3711cd^) in test/nemoback-habitat/common.py and run_all.py (git show HEAD:... in coverage())'})
    seen, uniq = set(), []
    for x in fl:
        if x['text'] not in seen: seen.add(x['text']); uniq.append(x)
    order = {'bug': 0, 'risk': 1, 'info': 2, 'known': 3}
    uniq.sort(key=lambda x: order.get(x['severity'], 9))
    cnt = {s: sum(x['severity'] == s for x in uniq) for s in order}
    L += ['', '## Findings (%d: %s)' % (len(uniq), ', '.join('%s %d' % kv for kv in cnt.items())), '', 'Severity: **bug** = wrong behaviour in the code under test; **risk** = correct today, fragile / surprising; **info** = observation; **known** = expected / next migration round / pre-existing.', '',
          '| severity | finding | minimal reproducer | suggested fix (not applied) |', '|---|---|---|---|']
    esc = lambda s: (s or '').replace('|', '/').replace('\n', ' ')
    L += ['| %s | %s | %s | %s |' % (x['severity'], esc(x['text']), ('`%s`' % esc(x['repro'])) if x.get('repro') else '', esc(x.get('fix'))) for x in uniq]
    L += ['', '## Method coverage%s' % ('' if not missing else ' -- **not exercised: %s**' % missing), '',
          '`vs original` compares each method with `%s:metaibm/patch.py` (AST + source text). `scripts` = test scripts that name the method; `reached via` = patch methods that call it internally (from the chain test).' % orig_commit, '',
          '| method | vs original | scope | scripts | reached via |', '|---|---|---|---|---|']
    L += ['| `%s` | %s | %s | %s | %s |' % r for r in rows]
    mx = chain.get('matrix', {})
    if mx:
        L += ['', '## Call-chain matrix (observed)', '', 'What each public patch method reached while it was the outermost call (union over all calls in `test_chain.py`).', '',
              '| patch method | habitat | indfactory | individual | patch (internal) |', '|---|---|---|---|---|']
        L += ['| `%s` | %s | %s | %s | %s |' % (m, ', '.join(v['habitat']) or '-', ', '.join(v['indfactory']) or '-', ', '.join(v['individual']) or '-', ', '.join(v['patch (internal)']) or '-')
              for m, v in mx.items() if v['habitat'] or v['indfactory'] or v['individual'] or v['patch (internal)']]
        called = set(chain.get('call_counts', {}))
        quiet = [m for m, v in mx.items() if m in called and not (v['habitat'] or v['indfactory'] or v['individual'] or v['patch (internal)'])]
        L += ['', 'Pure patch-level (no habitat / indfactory / individual call): %s' % ', '.join('`%s`' % m for m in quiet),
              '', 'Not called in the chain process (asex / mix, probed in separate processes by `test_nemo.py` edges_known): %s' % ', '.join('`%s`' % m for m in mx if m not in called)]
        sc = chain.get('static_scan', [])
        L += ['', '### Static scan: patch calls in metacommunity.py / simulator.py / extension/*.py', '', '| location | caller | patch method | issue |', '|---|---|---|---|'] + ['| %s | %s | `%s` | %s |' % tuple(r) for r in sc]
        if not any('extension' in r[0] for r in sc): L += ['', '`extension/global_habitat_network.py` reads `patch_object.set` / habitat pools directly and calls no patch method: nothing to flag.']
    L += ['', '## Repository state', '', '```', git('status', '--short'), git('diff', '--stat'), '```']
    open(os.path.join(HERE, 'TEST_REPORT.md'), 'w').write('\n'.join(L) + '\n')
    print('coverage: %s' % ('all in-scope methods exercised' if not missing else 'MISSING %s' % missing))
    print('report:', os.path.join(HERE, 'TEST_REPORT.md'))
    sys.exit(1 if n_fail else 0)


if __name__ == '__main__':
    main()
