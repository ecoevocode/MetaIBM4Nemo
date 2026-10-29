# Run smoke.py (stop if it fails), then every test, each in its own process; check method coverage; write TEST_REPORT.md.
# python run_all.py --perf  also records timings (perf(), record only, no assertions).
import os, re, ast, sys, json, glob, time, hashlib, platform, subprocess, datetime

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, '..', '..'))
TESTS = [('test_bookkeeping.py', 'pool / microsite bookkeeping vs original, dormancy branches, clear-up, sampling uniformity (fake individuals, no Nemo)'),
         ('test_nemo.py', 'real Nemo individuals: reproduction | lifecycle | memory | edges_known (one section each)')]
ORIGINAL_COMMIT = 'baace34'                                     # same as common.ORIGINAL_COMMIT: pre-NemoBackend habitat.py (c3711cd^)
NOISE = ('setting random seed from input value', 'outputs: ', 'the microsite in the habitat is occupied.', 'the microsite in the habitat is empty.')


def run(script):
    t = time.time()
    p = subprocess.run([sys.executable, os.path.join(HERE, script)], capture_output=True, text=True, cwd=HERE, env=dict(os.environ, PYTHONDONTWRITEBYTECODE='1'))
    lines = [l for l in (p.stdout + '\n' + p.stderr).splitlines() if l.strip() and not any(n in l for n in NOISE) and not l.startswith('[FINDING')]
    return {'rc': p.returncode, 's': time.time() - t, 'tail': lines[-25:]}


def load(name):
    try: return json.load(open(os.path.join(HERE, 'results', name + '.json')))
    except Exception: return {}


def git(*a):
    return subprocess.run(['git'] + list(a), capture_output=True, text=True, cwd=REPO).stdout.rstrip()


def coverage():
    ''' every habitat method -> number of test scripts that call it; in-scope methods (not asex / mix) must be called somewhere '''
    src = open(os.path.join(REPO, 'metaibm', 'habitat.py')).read()
    cls = [n for n in ast.parse(src).body if isinstance(n, ast.ClassDef) and n.name == 'habitat'][0]
    orig = {n.name: ast.dump(n) for n in ast.walk(ast.parse(git('show', ORIGINAL_COMMIT + ':metaibm/habitat.py'))) if isinstance(n, ast.FunctionDef)}
    texts = [open(p).read() for p in glob.glob(os.path.join(HERE, '*.py')) if not p.endswith('run_all.py')]
    via = {'__init__': 'constructor', '__str__': 'str(hab)', '_empty_microsite_grid': '__init__', 'get_microsite_env_val_ls': 'dead selection',
           'survival_rate': 'dead selection', 'hab_sexual_pairwise_parents_ls': '..._num()', '_recycle_dropped': 'dormancy processes', '_germinate_from_pools': 'the 3 germination wrappers'}
    rows, missing = [], []
    for n in cls.body:
        if not isinstance(n, ast.FunctionDef): continue
        m = n.name
        status = 'new' if m not in orig else 'changed' if orig[m] != ast.dump(n) else 'unchanged'
        scope = 'asex/mix (probe only)' if re.match(r'hab_(asex|mix|mixed)_|hab_asexual_reprodece', m) else 'in scope'
        hits = sum(bool(re.search(r'\.%s\b' % re.escape(m), t)) for t in texts)
        rows.append((m, status, scope, hits, via.get(m, '')))
        if scope == 'in scope' and not hits and m not in via: missing.append(m)
    return rows, missing


def perf():
    ''' record only: per-phase time of one v3 step vs habitat size; germination vs candidates (list.remove -> O(n^2)) '''
    sys.path.insert(0, HERE)
    import random
    import numpy as np
    from common import habitat, fake_indi, make_species, make_hab, recycle_everything, seed_all
    out = {'per_step_s': {}, 'germinate_s': {}}
    seed_all(91)
    sp = make_species('sp1', 0.2, seed=91)
    for L in (25, 50, 100):
        h = make_hab(L=L)
        h.hab_initialize({'sp1': sp}, [[0.2, 0.2]])
        tm = {}
        for _ in range(3):
            for name, fn in (('dead_selection', lambda: h.hab_dead_selection(0.1, 0.5, 'niche_gaussian')),
                             ('reproduce', lambda: h.hab_sex_reproduce_mutate_into_offspring_pool(1.0)),
                             ('germinate', h.hab_local_germinate_from_offspring_and_immigrant_pool),
                             ('clear_up', lambda: h.hab_clear_up_offspring_pool() + h.hab_clear_up_immigrant_pool())):
                t = time.perf_counter(); fn(); tm[name] = tm.get(name, 0) + (time.perf_counter() - t) / 3
        out['per_step_s']['L=%d' % L] = {k: round(v, 4) for k, v in tm.items()}
        recycle_everything(h)
    for n in (2500, 5000, 10000, 20000):
        side = int(n ** 0.5) + 1
        np.random.seed(0)
        h = habitat('h', 0, (0, 0), 2, ['env1', 'env2'], [0.2, 0.2], [0.025, 0.025], side, side)
        h.offspring_pool = [fake_indi('o%d' % i) for i in range(2 * n)]
        random.seed(0)
        t = time.perf_counter(); h.hab_local_germinate_from_offspring_and_immigrant_pool()
        out['germinate_s']['%d candidates -> %d sites' % (2 * n, side * side)] = round(time.perf_counter() - t, 3)
    return out


def main():
    for f in glob.glob(os.path.join(HERE, 'results', '*.json')): os.remove(f)
    md5 = hashlib.md5(open(os.path.join(REPO, 'metaibm', 'habitat.py'), 'rb').read()).hexdigest()
    try:
        import importlib.metadata as md; pv = md.version('pynemo')
    except Exception: pv = '?'
    L = ['# habitat.py test report (NemoBackend, sexual)', '', 'Generated %s by `test/nemoback-habitat/run_all.py`. Scope: `metaibm/habitat.py` (patch / metacommunity layers not tested; asex / mix only probed).' % datetime.datetime.now().strftime('%Y-%m-%d %H:%M'), '',
         '- python `%s`, pynemo `%s`, `%s`' % (sys.version.split()[0], pv, platform.platform()),
         '- habitat.py md5 `%s`; original reference = `git show %s:metaibm/habitat.py` (HEAD `%s`)' % (md5, ORIGINAL_COMMIT, git('rev-parse', '--short', 'HEAD')), '']
    smoke = run('smoke.py'); sj = load('smoke')
    print('smoke.py %s %.1fs' % ('PASS' if smoke['rc'] == 0 else 'FAIL', smoke['s']))
    L += ['## Smoke: %s' % ('PASS' if smoke['rc'] == 0 else '**FAIL**'), '', '| step | return |', '|---|---|'] + ['| %s | `%s` |' % (s['step'], s['return']) for s in sj.get('steps', [])] + ['']
    if smoke['rc'] != 0:
        L += ['Failed step `%s`:' % sj.get('failed_step'), '```', sj.get('traceback', '\n'.join(smoke['tail'])), '```', 'Other tests not run.']
        open(os.path.join(HERE, 'TEST_REPORT.md'), 'w').write('\n'.join(L) + '\n'); sys.exit(1)
    res = {}
    for f, _ in TESTS:
        res[f] = run(f); print('%-22s %s %.1fs' % (f, 'PASS' if res[f]['rc'] == 0 else 'FAIL', res[f]['s']))
    rows, missing = coverage()
    pf = None
    if '--perf' in sys.argv:
        t = time.time(); pf = perf(); print('perf (record only) %.1fs' % (time.time() - t))
    n_fail = sum(r['rc'] != 0 for r in res.values()) + bool(missing)
    L += ['## Tests: %d / %d PASS' % (sum(r['rc'] == 0 for r in res.values()), len(res)), '', '| script | covers | result | time |', '|---|---|---|---|']
    L += ['| %s | %s | %s | %.1f s |' % (f, d, 'PASS' if res[f]['rc'] == 0 else '**FAIL**', res[f]['s']) for f, d in TESTS]
    L += ['', '## Key numbers', '']
    for f, _ in TESTS:
        L.append('<details><summary>%s</summary>' % f); L += ['', '```'] + res[f]['tail'] + ['```', '</details>']
    if pf: L += ['', '## Performance (record only, `--perf`)', '', '```', json.dumps(pf, indent=1), '```']
    fl = [(x['severity'], n, x['text']) for n in ('bookkeeping', 'nemo') for x in load(n).get('findings', [])]
    order = {'high': 0, 'medium': 1, 'low': 2, 'info': 3}
    L += ['', '## Findings / known behaviour (%d)' % len(fl), '', '| severity | finding |', '|---|---|'] + ['| %s | %s |' % (s, t.replace('|', '/')) for s, n, t in sorted(fl, key=lambda x: order.get(x[0], 9))]
    L += ['', '## Method coverage%s' % ('' if not missing else ' -- **not exercised: %s**' % missing), '', '| method | vs original | scope | scripts | reached via |', '|---|---|---|---|---|']
    L += ['| `%s` | %s | %s | %d | %s |' % r for r in rows]
    L += ['', '## Repository state', '', '```', git('status', '--short'), git('diff', '--stat', 'metaibm/', 'extension/'), '```']
    open(os.path.join(HERE, 'TEST_REPORT.md'), 'w').write('\n'.join(L) + '\n')
    print('coverage: %s' % ('all in-scope methods exercised' if not missing else 'MISSING %s' % missing))
    print('report:', os.path.join(HERE, 'TEST_REPORT.md'))
    sys.exit(1 if n_fail else 0)


if __name__ == '__main__':
    main()
