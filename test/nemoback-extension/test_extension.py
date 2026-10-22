# test/nemoback-extension: tests of the migrated extension/global_habitat_network.py (NemoBackend, sexual). Every part runs in its own process (Nemo's RNG and
# Individual::currentID are process-wide); helpers in common.py, full simulations in simulate_ext.py. Nothing in metaibm/, extension/, experiments/ is edited.
# Originals = the git bodies (extension at its last commit, add_patch at HEAD) compiled into a subclass of the CURRENT metacommunity.
#   python test_extension.py all [--quick|--full]          Step 0 baseline, smoke, static, perf, then registry / dispersal / marker / fuzz / sim in parallel, then TEST_REPORT.md
#   python test_extension.py <part> [--quick|--full]       one part -> results/<part>.json
#     baseline   Step 0: the other suites, folders backed up and restored byte-identical
#     smoke      one step of each flow (marker True / marker False / pool) with invariants after every call
#     static     A: git diff whitelist, signatures, pool_name tuple, guard / else, TYPE_CHECKING import, old-name residue, callers, NemoBackend tags
#     registry   B: add order, global_index record, distance matrix (brute force, v3, one recompute per add_patch), edges, rate-matrix properties x 5 kernels, kernel underflow
#     dispersal  C: offspring network dispersal -- guard, patch_num < 2, equivalence vs v3 + removal, conservation, v3 copy statistics, capping, distance law, Nemo level, chains
#     marker     D: marker network dispersal -- both is_remove modes vs v3, germinate-and-birth (habitat, home), known ValueError risk, default False -> True on GRFE births, chains
#     perf       F: one call at 16 / 64 / 256 / 1024 habitats vs v3, shares kernel / matrix_around / sampling, a zero-skip variant (test-only) for the fix suggestion
#     fuzz       G: random mixes of every meta operation + both network dispersals + guarded calls, invariants after each
#     sim        E: simulate_ext.py runs (3 groups + patch-level comparator, 2 scenarios), reproducibility, whole-run equivalence vs v3, ecology, high disturbance
#     report     TEST_REPORT.md from results/*.json
import os, re, ast, sys, json, time, glob, shutil, random, hashlib, inspect, platform, tempfile, subprocess, traceback, warnings, datetime
from collections import Counter, defaultdict

MODE = 'full' if '--full' in sys.argv else ('quick' if '--quick' in sys.argv else os.environ.get('EXT_MODE', 'quick'))
sys.argv = [a for a in sys.argv if a not in ('--full', '--quick')]
os.environ['EXT_MODE'] = MODE

import numpy as np                                                                  # noqa: E402
import common as C                                                                  # noqa: E402
from common import metacommunity, ext, build_world, random_spec, brute_distance, md5, original_metacommunity, finding, save, rng_state, same_rng_state, seed_all, KERNELS, HERE, REPO   # noqa: E402

sys.dont_write_bytecode = True
warnings.filterwarnings('ignore', category=RuntimeWarning); warnings.filterwarnings('ignore', category=PendingDeprecationWarning)
T0 = time.time()
RESULTS = []
OFF, MARK = C.OFF, C.MARK
DISP = OFF
REG = 'global_habitat_network_registry'
SIM_PY = os.path.join(HERE, 'simulate_ext.py')
PY_ENV = dict(os.environ, PYTHONDONTWRITEBYTECODE='1', MPLBACKEND='Agg', EXT_MODE=MODE)


def Q(quick, full): return full if MODE == 'full' else quick


def check(name, ok, detail=''):
    RESULTS.append({'check': name, 'ok': bool(ok), 'detail': str(detail)[:1500]})
    print('%s  %s  %s' % ('PASS' if ok else 'FAIL', name, str(detail)[:300]), flush=True)
    return ok


def _spawn(args, timeout=None):
    ''' run this file with args in a fresh process; the worker prints 'JSON <payload>' '''
    pr = subprocess.run([sys.executable, os.path.abspath(__file__)] + list(args), capture_output=True, text=True, cwd=HERE, env=PY_ENV, timeout=timeout)
    line = [l for l in pr.stdout.splitlines() if l.startswith('JSON ')]
    if not line: return {'error': 'rc=%d %s' % (pr.returncode, (pr.stdout + pr.stderr)[-2500:])}
    return json.loads(line[-1][5:])


def _pmap(jobs, max_workers=4):
    ''' run worker arg-lists in parallel processes; returns results in order '''
    from concurrent.futures import ThreadPoolExecutor
    with ThreadPoolExecutor(max_workers=max_workers) as ex:
        return list(ex.map(_spawn, jobs))


def _strip_comment(line):
    ''' code part of one source line (tokenize: a '#' inside a string literal is kept) '''
    import io, tokenize
    try:
        toks = list(tokenize.generate_tokens(io.StringIO(line + '\n').readline))
    except (tokenize.TokenError, IndentationError, SyntaxError):
        return line.split('#')[0].rstrip()
    cut = [t.start[1] for t in toks if t.type == tokenize.COMMENT]
    return line[:cut[0]].rstrip() if cut else line


def _code_lines(src, name):
    ''' stripped non-empty lines of a module-level def, docstring removed '''
    fn = C._functions(src)[name]
    lines = src.splitlines()[fn.lineno - 1:fn.end_lineno]
    d = fn.body[0] if fn.body and isinstance(fn.body[0], ast.Expr) and isinstance(getattr(fn.body[0], 'value', None), ast.Constant) and isinstance(fn.body[0].value.value, str) else None
    skip = set(range(d.lineno - fn.lineno, d.end_lineno - fn.lineno + 1)) if d else set()
    return [l.strip() for k, l in enumerate(lines) if k not in skip and l.strip()]


def _counter(ret):
    mt = re.search(r'there are (\d+) individuals disperse into immigrant_pool', ret or '')
    return int(mt.group(1)) if mt else None


def _mcounter(ret):
    mt = re.search(r'there are (\d+) offspring markers disperse into immigrant_marker_pool', ret or '')
    return int(mt.group(1)) if mt else None


def _meta_guard_message():
    fn = C._functions(open(os.path.join(REPO, C.META_FILE)).read(), 'metacommunity')['dispersal_aomng_patches_from_offspring_pool_to_immigrant_pool']
    for n in ast.walk(fn):
        if isinstance(n, ast.Raise) and isinstance(n.exc, ast.Call) and getattr(n.exc.func, 'id', '') == 'NotImplementedError': return n.exc.args[0].value


def nemo_counts(F): return (C.minted(F), C.recycled(F), C.indfactory.get_current_individual_id())


def part_cli(name, fn):
    out = {'mode': MODE, 'original_ext_commit': C.ORIGINAL_EXT_COMMIT, 'original_meta_commit': C.ORIGINAL_META_COMMIT, 'gidx_dict_name': C.GIDX_DICT}
    try:
        out.update(fn() or {}); out['ok'] = all(r['ok'] for r in RESULTS) and bool(RESULTS)
    except BaseException:
        out.update(ok=False, traceback=traceback.format_exc()); print(out['traceback'])
    out['results'] = RESULTS; out['seconds'] = round(time.time() - T0, 2)
    save(name, out)
    print('%s %s: %d/%d checks, %d findings (%.2f s)' % (name.upper(), 'PASSED' if out['ok'] else 'FAILED', sum(r['ok'] for r in RESULTS), len(RESULTS), len(C.FINDINGS), out['seconds']))
    sys.stdout.flush(); os._exit(0 if out['ok'] else 1)


# ############################################################################
# PART: smoke
# ############################################################################
def smoke_main():
    seed_all(11); F = C.species(11); fs = list(F.values())
    ml = C.make_mainland(L=8); isl = C.make_islands(3, 2, L=4, K=5, gidx0=100)
    ml.meta_initialize(F, C.TABLE); isl.meta_initialize(F, C.TABLE)
    rain = C.rain_for(ml, 0.01); rows = []; t = time.time()
    flows = (('marker, network is_remove=True (default)', C.marker_ops(ml, isl, rain, 0.5, 'exponential', {'rho': 2}, net_call=lambda: getattr(isl, MARK)(0.5, method='exponential', rho=2))),
             ('marker, network is_remove=False', C.marker_ops(ml, isl, rain, 0.5, 'exponential', {'rho': 2}, is_remove=False)),
             ('pool, network offspring (default)', C.pool_ops(ml, isl, rain, 0.5, 'exponential', {'rho': 2})))
    ok = True
    for flow, ops in flows:
        for name, op in ops:
            ret = op()
            try: C.check_world([ml, isl], fs, 'smoke %s %s' % (flow, name)); inv = True
            except AssertionError as e: inv = False; ok = False; print(e)
            rows.append([flow, name, (ret if isinstance(ret, str) else repr(ret)).strip()[:130], inv])
    C.teardown([ml, isl], fs)
    check('smoke: 3 flows x one step (mainland + 3x2 islands), invariants after every call, teardown minted == recycled', ok, '%d calls' % len(rows))
    return {'rows': rows, 'smoke_s': round(time.time() - t, 2)}


# ############################################################################
# PART: static  (A)
# ############################################################################
EXT, META, DOCS = C.EXT_FILE, C.META_FILE, 'docs-developer/extension-global-habitat-network.md'
STEP0_FILES = ['test/nemoback-metacommunity/test_metacommunity.py', 'test/nemoback-patch/test_chain.py']
_DN = '(?:%s)' % '|'.join(map(re.escape, C.GIDX_DICT_NAMES))
_FN = '(?:%s|%s)' % (OFF, MARK)
AGREED = {
    EXT: [('TYPE_CHECKING import', r'^from typing import TYPE_CHECKING$'),
          ('TYPE_CHECKING import', r'^if TYPE_CHECKING: from metaibm import patch, metacommunity$'),
          ('registry rename', r'^def global_habitat_id_idx_registry\(self, patch_object\):$'),
          ('registry rename', r"^def global_habitat_network_registry\(self, patch_object: ?'patch'\):$"),
          ('record stores habitat.global_index', r'^self\.global_habitat_id_2_index_dir\[\(patch_id, h_id\)\] = len\(self\.global_habitat_id_ls\) - 1$'),
          ('record stores habitat.global_index', r'^self\.%s\[\(patch_id, h_id\)\] = h_object\.global_index$' % _DN),
          ('distance matrix recomputed inside the registry', r'^self\.update_global_habitat_distance_matrix\(\)$'),
          ('is_remove in the signature', r"^def %s\(self, total_disp_among_rate, method='uniform', \*\*kwargs\):$" % _FN),
          ('is_remove in the signature', r"^def %s\(self(: ?'metacommunity')?, total_disp_among_rate, method='uniform', is_remove=True, \*\*kwargs\):$" % _FN),
          ('guard', r'^if not is_remove: raise NotImplementedError\('),
          ('sampling', r'^migrants_to_j_indi_object_ls \+= random\.sample\(h_i_object\.offspring_pool, migrants_i_j_num\)$'),
          ('sampling', r'^if is_remove: migrants_to_j_indi_object_ls \+= h_i_object\.sample_offspring_without_replacement\(migrants_i_j_num\)$'),
          ('else branch', r'^else: pass\s+# unreachable \(guard above\); TODO\(copy_from\)'),
          ('sampling', r'^migrants_to_j_offspring_marker_ls \+= random\.sample\(h_i_object\.offspring_marker_pool, migrants_i_j_num\)$'),
          ('sampling', r"^if is_remove: migrants_to_j_offspring_marker_ls \+= h_i_object\.sample_offspring_without_replacement\(migrants_i_j_num, pool_name=\('offspring_marker_pool',\)\)$"),
          ('sampling', r'^else: migrants_to_j_offspring_marker_ls \+= random\.sample\(h_i_object\.offspring_marker_pool, migrants_i_j_num\)$'),
          ('install name', r'^cls\.global_habitat_id_idx_registry = global_habitat_id_idx_registry$'),
          ('install name', r'^cls\.global_habitat_network_registry = global_habitat_network_registry$')],
    META: [('add_patch: registry call renamed', r'^self\.global_habitat_id_idx_registry\(patch_object\)'),
           ('add_patch: registry call renamed', r'^self\.global_habitat_network_registry\(patch_object\)'),
           ('add_patch: distance-matrix call moved into the registry', r'^self\.update_global_habitat_distance_matrix\(\)'),
           ('add_patch: distance-matrix call moved into the registry', r'^#self\.incremental_update_global_habitat_distance_matrix\(\)$')],
    DOCS: [('docs heading renamed', r'^### `global_habitat_(id_idx_registry|network_registry)`$')]}
LAYOUT = [r'^if self\.patch_num < 2:$', r"^log_info = '\[Dispersal among patches in global habitat network\] in %s: patch_num < 2", r'^return log_info$',
          r"^if self\.patch_num < 2: return '\[Dispersal among patches in global habitat network\]"]


def _diff_lines(path):
    out = []
    for l in C.git('diff', '-U0', 'HEAD', '--', path).splitlines():
        if l.startswith(('+++', '---')) or not l.startswith(('+', '-')): continue
        out.append((l[0], l[1:]))
    return out


def classify_diff():
    ''' every +/- line of the three changed files: agreed label | comment only | blank | layout | UNEXPECTED '''
    rows = []
    for path in (EXT, META, DOCS):
        lines = _diff_lines(path)
        code = {s: Counter(_strip_comment(t.strip()) for k, t in lines if k == s) for s in '+-'}
        for sign, text in lines:
            t = text.strip()
            lab = next((lab for lab, rx in AGREED[path] if re.match(rx, t)), None)
            if lab is None:
                if not t: lab = 'blank line'
                elif path != DOCS and _strip_comment(t) and code['+' if sign == '-' else '-'][_strip_comment(t)]: lab = 'comment only'
                elif path == EXT and any(re.match(rx, t) for rx in LAYOUT): lab = 'layout: patch_num < 2 early return folded into one line'
                else: lab = 'UNEXPECTED'
            rows.append([path, sign, t[:160], lab])
    return rows


def _import_order(first):
    code = ('import sys; sys.path.insert(0, %r); %s; import metaibm, extension.global_habitat_network as e; '
            'print("OK", hasattr(metaibm.metacommunity, "global_habitat_network_registry"), "patch" in vars(e), "metacommunity" in vars(e))') % (REPO, first)
    pr = subprocess.run([sys.executable, '-c', code], capture_output=True, text=True, cwd=HERE, env=PY_ENV)
    return pr.returncode, (pr.stdout.strip().splitlines() or [''])[-1], pr.stderr[-400:]


def _git_grep(pattern, extended=False):
    args = ['grep', '-n', '-I', '--untracked'] + (['-E'] if extended else []) + [pattern]
    return [l for l in C.git(*args).splitlines() if l]


def _sig_ok(fn):
    sig = inspect.signature(fn); ps = list(sig.parameters.values()); names = [p.name for p in ps]
    return ('is_remove' in names and sig.parameters['is_remove'].default is True and names.index('is_remove') == names.index('method') + 1
            and ps[-1].kind == inspect.Parameter.VAR_KEYWORD and names.index('is_remove') == len(names) - 2), str(sig)


def static_main():
    out = {}
    src = C.current_ext_source()
    # A1 git diff whitelist + repository state
    rows = classify_diff(); out['diff'] = rows
    bad = [r for r in rows if r[3] == 'UNEXPECTED']
    check('A1 every changed line of extension / add_patch / docs is one of the agreed changes (or blank / comment / layout)', not bad, bad[:5])
    labels = Counter(r[3] for r in rows); out['diff_labels'] = dict(labels)
    for want in ('registry rename', 'record stores habitat.global_index', 'distance matrix recomputed inside the registry', 'is_remove in the signature', 'guard', 'else branch', 'sampling', 'install name', 'TYPE_CHECKING import'):
        check('A1 agreed change present: %s' % want, labels[want] > 0, labels[want])
    tolerated = [r for r in rows if r[3] in ('comment only', 'layout: patch_num < 2 early return folded into one line')]
    if tolerated:
        finding('info', '%s, %s' % (OFF, MARK) if any('patch_num' in r[2] for r in tolerated) else EXT,
                'diff lines outside the agreed list, non-functional: %s (behaviour pinned by C2 / MB1 patch_num < 2 and the equivalence tests)' % ', '.join('%s x%d' % (k, v) for k, v in labels.items() if k in ('comment only', 'layout: patch_num < 2 early return folded into one line')),
                evidence=['%s %s' % (r[1], r[2][:110]) for r in tolerated])
    st = [l for l in C.git('status', '--porcelain').splitlines() if l]
    mod = {l[3:] for l in st if not l.startswith('??')}; unt = {l[3:] for l in st if l.startswith('??')}
    check('A1 working tree: modified == extension, metacommunity.py, docs (+ the two Step 0 test files); untracked == test/nemoback-extension/ only',
          mod <= {EXT, META, DOCS} | set(STEP0_FILES) and {EXT, META, DOCS} <= mod and unt <= {'test/nemoback-extension/'}, st)
    # A2 signatures
    for name in (OFF, MARK):
        ok, s = _sig_ok(getattr(metacommunity, name))
        check('A2 %s: is_remove=True right after method, before **kwargs' % name[-40:], ok, s)
    # A3 sampling arguments, guard, else branch
    fnm = C._functions(src)[MARK]
    calls = [n for n in ast.walk(fnm) if isinstance(n, ast.Call) and getattr(n.func, 'attr', '') == 'sample_offspring_without_replacement']
    pn = [kw.value for c in calls for kw in c.keywords if kw.arg == 'pool_name']
    check("A3 marker: pool_name=('offspring_marker_pool',) is a TUPLE (a bare string would iterate its characters)",
          len(calls) == 1 and len(pn) == 1 and isinstance(pn[0], ast.Tuple) and [getattr(e, 'value', None) for e in pn[0].elts] == ['offspring_marker_pool'], '; '.join(ast.unparse(c) for c in calls))
    check('A3 marker: no guard (is_remove=False stays valid for markers)', not [n for n in ast.walk(fnm) if isinstance(n, ast.Raise)])
    fno = C._functions(src)[OFF]
    callso = [n for n in ast.walk(fno) if isinstance(n, ast.Call) and getattr(n.func, 'attr', '') == 'sample_offspring_without_replacement']
    check('A3 offspring: sample_offspring_without_replacement(migrants_i_j_num) with the default pool (no pool_name)', len(callso) == 1 and not callso[0].keywords and len(callso[0].args) == 1, '; '.join(ast.unparse(c) for c in callso))
    body = fno.body[1:] if fno.body and isinstance(fno.body[0], ast.Expr) and isinstance(getattr(fno.body[0], 'value', None), ast.Constant) else fno.body
    first = body[0] if body else None
    guard_first = isinstance(first, ast.If) and ast.unparse(first.test) == 'not is_remove' and isinstance(first.body[0], ast.Raise) and getattr(first.body[0].exc.func, 'id', '') == 'NotImplementedError'
    check('A3 offspring: guard is the first statement (before patch_num < 2 and any RNG / state change)', guard_first, ast.unparse(first).splitlines()[0] if first else '')
    msg = first.body[0].exc.args[0].value if guard_first and first.body[0].exc.args and isinstance(first.body[0].exc.args[0], ast.Constant) else None
    check('A3 offspring: guard text == metacommunity.dispersal_aomng_patches_from_offspring_pool_to_immigrant_pool', msg is not None and msg == _meta_guard_message(), repr(msg))
    code = _code_lines(src, OFF)
    check('A3 offspring: else branch is `pass  # unreachable (guard above); TODO(copy_from) ...`', any(re.match(r'^else: pass\s+# unreachable \(guard above\); TODO\(copy_from\)', l) for l in code))
    check('A3 offspring: no random.sample left in code (# comments stripped first)', not any('random.sample' in _strip_comment(l) for l in code))
    # A4 TYPE_CHECKING
    tree = ast.parse(src)
    tc = [n for n in tree.body if isinstance(n, ast.If) and ast.unparse(n.test) == 'TYPE_CHECKING']
    check('A4 `from metaibm import patch, metacommunity` sits under `if TYPE_CHECKING:` (module level)', len(tc) == 1 and any(isinstance(b, ast.ImportFrom) and b.module == 'metaibm' for b in tc[0].body))
    check('A4 runtime: the names are not bound in the extension module (TYPE_CHECKING is False)', 'patch' not in vars(ext) and 'metacommunity' not in vars(ext))
    for first_imp in ('import metaibm', 'import extension.global_habitat_network'):
        rc, line, err = _import_order(first_imp)
        check('A4 fresh process, %s first: imports work, nothing bound at runtime' % first_imp, rc == 0 and line == 'OK True False False', 'rc=%d %r %s' % (rc, line, err.strip()[-200:]))
    anns = {n: {k: p.annotation for k, p in inspect.signature(getattr(ext, n)).parameters.items() if p.annotation is not inspect._empty} for n in (REG, OFF, MARK)}
    check('A4 annotations are strings (never evaluated)', all(isinstance(v, str) for d in anns.values() for v in d.values()), anns)
    import typing
    try: typing.get_type_hints(getattr(ext, REG)); gth = 'resolved'
    except NameError as e: gth = 'NameError: %s' % e
    out['get_type_hints'] = gth
    # A5 old name residue + A6 callers
    old = _git_grep(C.ORIGINAL_REGISTRY)
    def kind(l):
        p = l.split(':', 1)[0]
        if p.startswith('test/nemoback-extension/'): return 'this suite (v3 reference name)'
        if p.startswith(('docs', 'docs-developer')): return 'docs'
        if p.startswith('test/'): return 'other test'
        return 'code'
    out['old_name'] = [[l[:200], kind(l)] for l in old]
    live = [l for l, k in out['old_name'] if k in ('code', 'docs')]
    check('A5 old name global_habitat_id_idx_registry: no residue outside this suite (repo-wide incl. experiments, examples, .ipynb, docs, test)', not live, live[:5])
    callers = []
    for l in _git_grep('dispersal_among_patches_in_global_habitat_network_from_(offspring_pool_to_immigrant_pool|offspring_marker_pool_to_immigrant_marker_pool)', extended=True):
        p, ln, txt = l.split(':', 2)
        if p.startswith('test/nemoback-extension/'): continue
        which = 'offspring' if 'from_offspring_pool_to_immigrant_pool' in txt else 'marker'
        t = txt.strip(); t_code = t[1:] if t.startswith('"') else t                  # .ipynb: JSON string
        if p == EXT: k = 'definition / install'
        elif p.endswith('.md'): k = 'docs / report text'
        elif t_code.lstrip().startswith('#'): k = 'commented out'
        elif p.startswith('test/lecacy'): k = 'legacy v3 test'
        elif p.startswith('test/'): k = 'test suite'
        elif "'method':" in t or '"method":' in t: k = 'simulator config (method name)'
        else: k = 'live call'
        callers.append({'file': p, 'line': int(ln), 'function': which, 'kind': k, 'passes_is_remove': 'is_remove' in t, 'text': t[:170]})
    out['callers'] = callers
    false_calls = [c for c in callers if c['function'] == 'offspring' and c['kind'] in ('live call', 'simulator config (method name)') and 'is_remove=False' in c['text']]
    check('A6 no live caller passes is_remove=False to the offspring version (it raises)', not false_calls, false_calls)
    legacy = [c for c in callers if c['kind'] == 'legacy v3 test']
    if legacy:
        finding('info', '%s, %s' % (OFF, MARK), 'legacy v3 test folders still call the network dispersal functions (not part of any suite run here; the offspring calls now MOVE instead of copy)',
                evidence=['%s:%d' % (c['file'], c['line']) for c in legacy])
    exp_calls = [c for c in callers if c['file'].startswith(('experiments/', 'examples/', 'playgrounds/')) and c['kind'] in ('live call', 'simulator config (method name)')]
    finding('info', MARK, 'live callers in experiments / playgrounds use the marker version with the DEFAULT is_remove -> they now MOVE markers (v3 copied); quantified in D4 / E',
            evidence=['%s:%d (%s)' % (c['file'], c['line'], c['kind']) for c in exp_calls])
    # A7 NemoBackend tags; todo(accepted)
    tags = {n: 'NemoBackend' in (getattr(ext, n).__doc__ or '') for n in (REG, OFF, MARK)}
    check('A7 no NemoBackend tag on the registry / the two network dispersals (they touch no Nemo object)', not any(tags.values()), tags)
    all_tagged = [n for n, f in C._functions(src).items() if 'NemoBackend' in (ast.get_docstring(f) or '')]
    check('A7 no function of the extension is tagged NemoBackend', not all_tagged, all_tagged)
    for n in (REG, OFF, MARK):
        doc = [l.strip() for l in (getattr(ext, n).__doc__ or '').strip().splitlines()]
        if not (len(doc) == 3 and doc[1].startswith('INPUT') and doc[2].startswith(('RETURN', 'OUTPUT'))):
            finding('todo', n, 'docstring kept as is (no three-line INPUT / RETURN format) -- accepted by the author', evidence=doc)
    spac = [m.group(0) for m in re.finditer(r"\b(?:self|patch_object):'\w+'", src)]
    if spac: finding('todo', '%s, %s, %s' % (REG, OFF, MARK), 'annotation spacing kept as is (no space after the colon) -- accepted by the author', evidence=sorted(set(spac)))
    out['tags'] = tags
    out['gidx_dict'] = {'name': C.GIDX_DICT, 'known_names': C.GIDX_DICT_NAMES}
    return out


# ############################################################################
# PART: registry  (B)
# ############################################################################
FIXED_SPEC = [   # unsorted names, shuffled non-contiguous global_index incl. both ends of 0..65535, unequal habitat counts
    ('pB', 0, (0.0, 0.0), [('hz', 0, (0.0, 0.0), 907), ('ha', 1, (0.7, -0.3), 12)]),
    ('pA', 1, (3.0, 1.0), [('h1', 0, (3.0, 1.0), 65535), ('h0', 1, (3.7, 0.7), 0), ('h2', 2, (4.4, 0.4), 400)]),
    ('pC', 2, (6.0, 2.0), [('x', 0, (6.0, 2.0), 33), ('b', 1, (6.7, 1.7), 7), ('a', 2, (7.4, 1.4), 1000), ('m', 3, (8.1, 1.1), 5)]),
]
RATE_GRID = (0.0, 0.01, 0.1, 0.5, 0.9, 1.0)


def rate_props(M, hl, m):
    ''' diagonal == 1 - m, same-patch off-diagonal == 0, rows sum to 1; NaN rows counted separately '''
    A = np.asarray(M); H = len(hl)
    nan_rows = [i for i in range(H) if np.isnan(A[i]).any()]
    ok_rows = [i for i in range(H) if i not in nan_rows]
    diag = all(A[i, i] == 1 - m for i in range(H))
    same = all(A[i, j] == 0 for i in range(H) for j in range(H) if i != j and hl[i][0] == hl[j][0] and i in ok_rows)
    rows = all(abs(A[i].sum() - 1) < 1e-12 for i in ok_rows)
    return {'nan_rows': len(nan_rows), 'diag': diag, 'same_patch_zero': same, 'row_sum_1': rows}


def compare_world(spec, seed, rates=(0.3,)):
    ''' build the same world with the current and the original registry + add_patch; returns (mismatches, rate-property failures) '''
    Orig = original_metacommunity(ext_names=[C.ORIGINAL_REGISTRY], meta_names=['add_patch'])
    cur, cur_calls = build_world(spec, metacommunity, seed); cur_rng = rng_state()
    org, org_calls = build_world(spec, Orig, seed); org_rng = rng_state()
    bad, props_bad = [], []
    n_habs = sum(len(h) for *_, h in spec)
    if cur.global_habitat_id_ls != [(pn, hn) for pn, _, _, habs in spec for hn, *_ in habs]: bad.append('ls != add_patch x add_habitat order')
    if cur.global_habitat_id_ls != org.global_habitat_id_ls: bad.append('ls != original')
    gi = {(pn, hn): g for pn, _, _, habs in spec for hn, _, _, g in habs}
    rec = C.gidx_dict(cur)
    if rec != gi: bad.append('record != global_index passed to add_habitat')
    if any(cur.set[a].set[b].global_index != v for (a, b), v in rec.items()): bad.append('record != habitat.global_index')
    if org.global_habitat_id_2_index_dir != {k: i for i, k in enumerate(org.global_habitat_id_ls)}: bad.append('original record is not the row number (reference changed?)')
    if list(rec) != cur.global_habitat_id_ls: bad.append('record key order != ls')
    if cur_calls != [1] * len(spec): bad.append('distance-matrix calls per add_patch %s' % cur_calls)
    if cur_calls != org_calls: bad.append('calls %s != original %s' % (cur_calls, org_calls))
    Dc, Do = cur.global_habitat_distance_matrix, org.global_habitat_distance_matrix
    if not isinstance(Dc, np.matrix): bad.append('distance matrix type %s' % type(Dc).__name__)
    if n_habs and Dc.shape != (n_habs, n_habs): bad.append('shape %s' % (Dc.shape,))
    if Dc.shape != Do.shape or md5(Dc) != md5(Do): bad.append('distance matrix != original')
    if n_habs and not np.array_equal(np.asarray(Dc), brute_distance(cur)): bad.append('distance matrix != brute force')
    if not same_rng_state(cur_rng, org_rng): bad.append('RNG state after build != original')
    if len([p for p in spec if p[3]]) >= 2:
        for kern, kw in KERNELS:
            for m in rates:
                with np.errstate(all='ignore'):
                    a = cur.global_habitat_dispersal_among_rate_matrix(m, method=kern, **kw); b = org.global_habitat_dispersal_among_rate_matrix(m, method=kern, **kw)
                if md5(a) != md5(b): bad.append('rate matrix (%s, m=%s) != original' % (kern, m))
                pr = rate_props(a, cur.global_habitat_id_ls, m)
                if pr['nan_rows'] or not (pr['diag'] and pr['same_patch_zero'] and pr['row_sum_1']): props_bad.append((kern, m, pr))
    return bad, props_bad


def kernel_threshold(m, kern, kw, hi):
    ''' smallest distance where dispersal_kernel_strength underflows to exactly 0.0 (log-scale bisection, then linear); None if not below hi '''
    f = lambda d: m.dispersal_kernel_strength(np.float64(d), method=kern, **kw)          # the pipeline passes np.float64 (a Python float would raise OverflowError)
    with np.errstate(all='ignore'):
        if f(hi) != 0.0: return None
        lo = 1e-6
        if f(lo) == 0.0: return lo
        while hi / lo > 2:
            mid = float(np.exp((np.log(lo) + np.log(hi)) / 2))                           # lo * hi would overflow
            if f(mid) == 0.0: hi = mid
            else: lo = mid
        for _ in range(200):
            mid = (lo + hi) / 2
            if f(mid) == 0.0: hi = mid
            else: lo = mid
            if hi - lo <= 1e-12 * hi: break
    return float(hi)


def underflow_part():
    ''' B6: kernel underflow -> NaN row -> the habitat silently neither emits nor (through the column) receives '''
    out = {}
    m0 = metacommunity('probe')
    grid = [('exponential', {'rho': r}) for r in (0.2, 1, 2, 5, 10)] + [('gaussian', {'sigma': s}) for s in (0.5, 1.5, 5)] + \
           [('cauchy', {'gamma': 1.0}), ('power_law', {'alpha': 2.0, 'r0': 1.0}), ('power_law', {'alpha': 8.0, 'r0': 1.0}), ('uniform', {})]
    th = []
    for kern, kw in grid:
        d = kernel_threshold(m0, kern, kw, 1e300)
        th.append({'kernel': kern, 'params': kw, 'zero_from_distance': d})
    out['thresholds'] = th
    exp2 = next(r['zero_from_distance'] for r in th if r['kernel'] == 'exponential' and r['params']['rho'] == 2)
    check('B6 exponential(rho=2) underflows to exactly 0.0 from d = %.6g (expected 745.13 / rho = 372.57: exp(-rho d) underflows first)' % (exp2 or float('nan')), exp2 is not None and abs(exp2 - 745.1332 / 2) < 0.01, exp2)
    # experiments: largest habitat distance vs the threshold of their kernel
    import pandas as pd
    lay = pd.read_csv(os.path.join(REPO, 'experiments', 'patch_habitat_layouts.csv')); lay = lay[lay['patch_num'] == 16]
    xy = lay[['habitat_x_location', 'habitat_y_location']].to_numpy(dtype=float)
    dmax_grfe = float(np.sqrt(((xy.max(0) - xy.min(0)) ** 2).sum()))
    th02 = next(r['zero_from_distance'] for r in th if r['kernel'] == 'exponential' and r['params']['rho'] == 0.2)
    out['experiments'] = {'grfe: max habitat distance (patch_num 16)': round(dmax_grfe, 2), 'grfe: exponential rho=0.2 zero from': th02,
                          'net: max habitat distance (16 patches x 4, coordinates 0..7)': round(float(np.sqrt(2) * 7), 2), 'net: exponential rho=2 zero from': exp2}
    check('B6 not reachable in the experiments (max habitat distance << underflow distance)', dmax_grfe < th02 and np.sqrt(2) * 7 < exp2, out['experiments'])
    # pipeline: 3 patches, A far away (d = 400 > 372.6 at rho=2), B and C close to each other; control: A at d = 300 (kernel ~1e-261, no underflow)
    def pipeline(xa):
        spec = [('pA', 0, (0.0, 0.0), [('h0', 0, (0.0, 0.0), 11)]), ('pB', 1, (xa, 0.0), [('h0', 0, (xa, 0.0), 12), ('h1', 1, (xa + 0.5, 0.0), 13)]),
                ('pC', 2, (xa + 1.0, 0.0), [('h0', 0, (xa + 1.0, 0.0), 14)])]
        m, _ = build_world(spec, metacommunity, 5)
        hobj = [m.set[p].set[h] for p, h in m.global_habitat_id_ls]
        for (p, hn), h in zip(m.global_habitat_id_ls, hobj):
            h.offspring_marker_pool = [(p, hn, 'sexual')] * 40
            h.offspring_pool = [object() for _ in range(40)]
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter('always')
            M = np.asarray(m.global_habitat_dispersal_among_rate_matrix(0.5, method='exponential', rho=2))
            OM = np.asarray(m.get_global_habitat_network_offspring_marker_pool_num_matrix() * m.global_habitat_dispersal_among_rate_matrix(0.5, method='exponential', rho=2))
        warn = sorted({'%s: %s' % (x.category.__name__, str(x.message)[:80]) for x in w if x.category is not PendingDeprecationWarning})
        sent = Counter(); recv = Counter(); logs = set()
        base = [(list(h.offspring_pool), list(h.offspring_marker_pool)) for h in hobj]
        for r in range(Q(100, 400)):
            seed_all(r)
            for fn, kw in ((OFF, {}), (MARK, {'is_remove': True}), (MARK, {'is_remove': False})):
                for h, (op, omp) in zip(hobj, base): h.offspring_pool = list(op); h.offspring_marker_pool = list(omp); h.immigrant_pool = []; h.immigrant_marker_pool = []
                origin = {id(x): i for i, h in enumerate(hobj) for x in h.offspring_pool}
                ret = getattr(m, fn)(0.5, method='exponential', rho=2, **kw); logs.add(re.sub(r'\d+', 'N', ret.strip()))
                for j, h in enumerate(hobj):
                    src = [origin.get(id(x)) for x in h.immigrant_pool] if fn == OFF else [m.global_habitat_id_ls.index((mk[0], mk[1])) for mk in h.immigrant_marker_pool]
                    for i in src: sent[i] += 1; recv[j] += 1
        return {'distance_A_to_B': xa, 'nan_rows_rate_matrix': [list(m.global_habitat_id_ls[i]) for i in range(M.shape[0]) if np.isnan(M[i]).any()],
                'nan_cells_O_times_M': int(np.isnan(OM).sum()), 'cells': int(OM.size), 'warnings': warn, 'A_sent': sent[0], 'A_received': recv[0],
                'B_C_exchanged': sum(v for i, v in sent.items() if i != 0), 'logs': sorted(logs)}
    far, ctrl = pipeline(400.0), pipeline(300.0)
    out['pipeline'] = {'far (d=400)': far, 'control (d=300)': ctrl}
    check('B6 control (A at d=300, no underflow): no NaN, A emits and B <-> C exchange', not ctrl['nan_rows_rate_matrix'] and ctrl['A_sent'] > 0 and ctrl['B_C_exchanged'] > 0, ctrl)
    check('B6 far habitat (A at d=400): its rate-matrix row is NaN', far['nan_rows_rate_matrix'] == [['pA', 'h0']], far['nan_rows_rate_matrix'])
    H4 = int(round(np.sqrt(far['cells'])))
    check('B6 far habitat: O * M is a MATRIX PRODUCT, 0 * NaN = NaN spreads the NaN row into every cell outside A\'s own column (%d / %d NaN; A\'s diagonal is reset to 1 - m) '
          '-> after matrix_around nobody emits: A 0, B <-> C 0 (3 dispersal modes x %d seeds), the log just says 0' % (far['nan_cells_O_times_M'], far['cells'], Q(100, 400)),
          far['nan_cells_O_times_M'] == H4 * (H4 - 1) and far['A_sent'] == 0 and far['A_received'] == 0 and far['B_C_exchanged'] == 0, far)
    finding('risk', 'normalized_calculate_global_habitat_network_dispersal_kernel_strength_matrix -> both network dispersals',
            'kernel underflow switches the WHOLE network dispersal off: when every other-patch kernel value of ONE habitat underflows to 0.0, its row of D / D.sum() is 0/0 = NaN; '
            'the dispersal number matrix is the matrix PRODUCT O * M (O diagonal), so 0 * NaN = NaN reaches every cell outside that habitat\'s own column; matrix_around turns NaN into 0 -> no habitat emits anything, '
            'silently (a numpy RuntimeWarning only, the log reports 0). Not reachable in the experiments (GRFE max habitat distance %.1f vs underflow from %.0f at rho=0.2; '
            'net 9.9 vs 372.6 at rho=2); reachable with large coordinates, a large rho or a small sigma (exponential: d > 745.13 / rho; gaussian: d > 38.6 sigma; '
            'cauchy / power_law: no practical limit), and whenever a habitat has no habitat in another patch (e.g. patches without habitats, B6 last check). Pre-existing: the matrix code is unchanged from v3 (B2: rate matrices identical).' % (dmax_grfe, th02),
            evidence={'thresholds': [(r['kernel'], r['params'], None if r['zero_from_distance'] is None else float('%.6g' % r['zero_from_distance'])) for r in th],
                      'far': {k: far[k] for k in ('nan_rows_rate_matrix', 'nan_cells_O_times_M', 'cells', 'A_sent', 'B_C_exchanged', 'warnings')},
                      'control': {k: ctrl[k] for k in ('A_sent', 'B_C_exchanged')}},
            fix='scale rows with np.multiply(pool_sizes[:, None], M) (element-wise, not O * M) and treat rows with D.sum() == 0 as "no dispersal" (log it); '
                'optionally compute the exponential kernel relative to the row minimum distance so it cannot underflow')
    # other patches without habitats: patch_num >= 2 but nobody to disperse to
    spec2 = [('pA', 0, (0.0, 0.0), [('h0', 0, (0.0, 0.0), 21), ('h1', 1, (0.5, 0.0), 22)]), ('pB', 1, (3.0, 0.0), [])]
    m2, _ = build_world(spec2, metacommunity, 6)
    with np.errstate(all='ignore'): M2 = np.asarray(m2.global_habitat_dispersal_among_rate_matrix(0.5))
    for h in (m2.set['pA'].set['h0'], m2.set['pA'].set['h1']): h.offspring_marker_pool = [('pA', h.name, 'sexual')] * 10
    ret = getattr(m2, MARK)(0.5)
    out['empty_other_patch'] = {'nan_rows': int(np.isnan(M2).any(axis=1).sum()), 'log': ret.strip()}
    check('B6 patch_num == 2 but the other patch has no habitat: all rows NaN, 0 markers disperse, no error (no warning in the log)', out['empty_other_patch']['nan_rows'] == 2 and _mcounter(ret) == 0, out['empty_other_patch'])
    return out


def registry_main():
    out = {}
    # G install: rename only
    check('B0 new name installed', hasattr(metacommunity, REG))
    check('B0 old name gone', not hasattr(metacommunity, C.ORIGINAL_REGISTRY))
    def installed(src):
        fn = C._functions(src)['install_global_habitat_network_methods']
        return [t.attr for n in ast.walk(fn) if isinstance(n, ast.Assign) for t in n.targets if isinstance(t, ast.Attribute) and isinstance(t.value, ast.Name) and t.value.id == 'cls']
    cur_inst, org_inst = installed(C.current_ext_source()), installed(C.original_ext_source())
    renamed = [x if x != C.ORIGINAL_REGISTRY else REG for x in org_inst]
    check('B0 install list == original with the one rename', cur_inst == renamed, 'diff %s' % sorted(set(cur_inst) ^ set(renamed)))
    calls = [l.strip() for l in open(os.path.join(REPO, C.META_FILE)).read().split('def add_patch')[1].split('\n    def ')[0].splitlines() if 'global_habitat' in l]
    check('B0 add_patch calls only the registry (matrix recomputed inside it)', len(calls) == 1 and calls[0].startswith('self.global_habitat_network_registry(patch_object)'), calls)
    check('B0 record attribute name read from _init_global_habitat_network_fields: %s' % C.GIDX_DICT, C.GIDX_DICT in C.GIDX_DICT_NAMES, C.GIDX_DICT)
    # B1 fixed world
    bad, pb = compare_world(FIXED_SPEC, 3, RATE_GRID)
    check('B1 fixed world (3 patches, 2/3/4 habitats, unsorted names, global_index incl. 0 and 65535) == original and brute force', not bad, '; '.join(bad))
    check('B1 fixed world: rate matrix x 5 kernels x m in %s: diag == 1 - m, same patch == 0, rows sum to 1, no NaN' % (RATE_GRID,), not pb, pb[:3])
    # B2 random worlds
    rng = random.Random(20261021); fails, pfails = [], []
    n = Q(200, 1000)
    for s in range(n):
        bad, pb = compare_world(random_spec(rng), s, RATE_GRID if s % 4 == 0 else (0.3,))
        if bad: fails.append((s, bad))
        if pb: pfails.append((s, pb[:2]))
    check('B2 %d random worlds (2-5 patches, 1-4 habitats, unsorted names, random global_index): order, record, distance (v3 + brute force), one recompute per add_patch' % n, not fails, fails[:3])
    check('B2 %d random worlds: rate matrix == v3 (md5) and diag / same patch / row sum properties (every 4th world: 6 rates x 5 kernels)' % n, not pfails, pfails[:3])
    # B3 edges
    pn, _, l, h = FIXED_SPEC[1]
    bad, _ = compare_world([(pn, 0, l, h)], 1)
    check('B3 single patch == original', not bad, '; '.join(bad))
    for where in ('first', 'middle', 'last'):
        sp = [(pn_, i, l_, h_) for i, (pn_, _, l_, h_) in enumerate(FIXED_SPEC)]
        k = {'first': 0, 'middle': 1, 'last': 3}[where]
        sp.insert(k, ('p_empty', 0, (9.0, 9.0), []))
        sp = [(a, i, b, c) for i, (a, _, b, c) in enumerate(sp)]
        try: bad, pb = compare_world(sp, 2, (0.3,))
        except Exception as e: bad, pb = ['%s: %s' % (type(e).__name__, e)], []
        check('B3 patch without habitats (%s) == original, contributes no row' % where, not bad and not pb, '; '.join(bad) + str(pb[:1]))
    dup = [('pB', 0, (0.0, 0.0), [('h0', 0, (0.0, 0.0), 7), ('h1', 1, (0.5, 0.0), 7)]), ('pA', 1, (3.0, 0.0), [('h0', 0, (3.0, 0.0), 7)])]
    try:
        bad, _ = compare_world(dup, 4); m, _ = build_world(dup, metacommunity, 4)
        rec = C.gidx_dict(m); ok = not bad and sorted(rec.values()) == [7, 7, 7]
    except Exception as e: ok, rec, bad = False, None, ['%s: %s' % (type(e).__name__, e)]
    check('B3 duplicate global_index: recorded as passed, no validation, no error (agreed: record only)', ok, '%s %s' % (rec, bad))
    finding('info', REG, 'duplicate global_index values are recorded without validation (agreed: the registry records, it does not assign or check); several habitats then share one Nemo homepatch', evidence=[[list(k), v] for k, v in (rec or {}).items()])
    m, _ = build_world(FIXED_SPEC, metacommunity, 4)
    m.set['pA'].add_habitat('late', 9, (5.0, 5.0), 2, ['env1', 'env2'], [0.2, 0.2], [0.025, 0.025], 3, 3, 5, global_index=4242)
    late = ('pA', 'late') in m.global_habitat_id_ls or ('pA', 'late') in C.gidx_dict(m)
    check('B3 habitat added after add_patch is not registered and the matrices keep their size (same as v3)', not late and m.global_habitat_distance_matrix.shape == (9, 9))
    finding('known', REG, 'a habitat added to a patch after add_patch is not registered (no row / column in the network, distance matrix not updated); unchanged from v3',
            evidence='pA.add_habitat("late") after add_patch -> not in global_habitat_id_ls, distance matrix stays 9 x 9', fix='validation module: forbid add_habitat after add_patch, or re-register the patch')
    # B4 the record is never used as a matrix index
    uses = []
    for d in ('metaibm', 'extension', 'experiments', 'examples', 'playgrounds'):
        for root, _, files in os.walk(os.path.join(REPO, d)):
            for f in files:
                if f.endswith(('.py', '.ipynb')):
                    for i, ln in enumerate(open(os.path.join(root, f), errors='ignore'), 1):
                        if any(nm in ln for nm in C.GIDX_DICT_NAMES): uses.append('%s:%d' % (os.path.relpath(os.path.join(root, f), REPO), i))
    check('B4 the record (%s) is only touched in the extension (init + registry), never read as a matrix index' % '/'.join(C.GIDX_DICT_NAMES), len(uses) == 2 and all(u.startswith(C.EXT_FILE) for u in uses), uses)
    # B5 RNG: building draws exactly what v3 draws (in compare_world); B6 underflow
    out['underflow'] = underflow_part()
    return out


# ############################################################################
# PART: dispersal  (C)  offspring_pool -> immigrant_pool
# ############################################################################
V3_SAMPLE = 'migrants_to_j_indi_object_ls += random.sample(h_i_object.offspring_pool, migrants_i_j_num)'
# reference = v3 body + "the picked ones leave their source pool". random.sample draws the same indices from any list of the same length, and
# habitat.sample_offspring_without_replacement samples (pool_name, obj) pairs of the same pool in the same order -> the migrated function must be bit-identical to it
REF_PATCH = {DISP: [(V3_SAMPLE, 'picked = random.sample(h_i_object.offspring_pool, min(migrants_i_j_num, len(h_i_object.offspring_pool))); '
                                '[h_i_object.offspring_pool.remove(x) for x in picked]; migrants_to_j_indi_object_ls += picked')]}
RATES = (0.1, 0.37, 0.5, 0.9, 1.0)
N_EQUIV = Q(60, 200)
CHAIN = [(1, Q(200, 500)), (2, Q(200, 500)), (3, Q(200, 500))]


def dispersal_bookkeeping():
    ''' C1-C2, C4, C8 in this process '''
    Ref = original_metacommunity(ext_names=[DISP], patches=REF_PATCH)
    seed_all(1); F = C.species(1); fs = list(F.values())
    def world(s, n_p, n_h):
        rng = random.Random(s); seed_all(s)
        m = C.make_islands(n_p, n_h, L=rng.choice((3, 4, 5)), W=rng.choice((3, 4)), K=rng.choice((0, 5, 20)), gidx0=100 + 50 * n_p, name='isl%d' % n_p)
        m.meta_initialize(F, C.TABLE); C.fill_pools(m); return m
    m3, m1 = world(5, 3, 2), world(6, 1, 2)
    msg = _meta_guard_message()
    for m in (m3, m1):
        for style, call in (('keyword', lambda m: getattr(m, DISP)(0.3, is_remove=False)), ('positional', lambda m: getattr(m, DISP)(0.3, 'uniform', False)),
                            ('kernel kwargs', lambda m: getattr(m, DISP)(0.3, method='exponential', is_remove=False, rho=2))):
            snap = C.full_snapshot([m3, m1], fs); err = None
            try: call(m)
            except NotImplementedError as e: err = str(e)
            except Exception as e: err = '%s: %s' % (type(e).__name__, e)
            check('C1 is_remove=False (%s, patch_num=%d): NotImplementedError, text as metacommunity, state / random / np.random / minted / recycled / currentID unchanged' % (style, m.patch_num),
                  err == msg and C.same_snapshot(snap, C.full_snapshot([m3, m1], fs)), repr(err)[:120])
    snap = C.full_snapshot([m3, m1], fs)
    a = getattr(m1, DISP)(0.3); b = getattr(Ref, DISP)(m1, 0.3)
    check('C2 patch_num < 2 (default is_remove): returns the log at once (== v3 text), nothing changed', a == b and C.same_snapshot(snap, C.full_snapshot([m3, m1], fs)), a.strip()[:100])
    hobj = [m3.set[p].set[h] for p, h in m3.global_habitat_id_ls]; hl = m3.global_habitat_id_ls
    origin = {id(x): i for i, h in enumerate(hobj) for x in h.offspring_pool}
    home0 = {id(x): x.nemo_indi_obj.getHome() for h in hobj for x in h.offspring_pool}
    imm0 = {id(x) for h in hobj for x in h.immigrant_pool}
    off0, imm_n0 = sum(len(h.offspring_pool) for h in hobj), sum(len(h.immigrant_pool) for h in hobj)
    nemo0 = nemo_counts(F)
    ret = getattr(m3, DISP)(0.9, method='exponential', rho=2)
    new = [(origin.get(id(x)), j, x) for j, h in enumerate(hobj) for x in h.immigrant_pool if id(x) not in imm0]
    moved = {id(x) for _, _, x in new}
    off1, imm_n1 = sum(len(h.offspring_pool) for h in hobj), sum(len(h.immigrant_pool) for h in hobj)
    check('C4 every new immigrant came from an offspring_pool of ANOTHER patch', len(new) > 0 and all(i is not None and hl[i][0] != hl[j][0] for i, j, _ in new), 'moved %d' % len(new))
    check('C4 migrants left their source pool (one owner per individual)', not any(id(x) in moved for h in hobj for x in h.offspring_pool))
    check('C4 offspring -n, immigrant +n, n == logged counter', off0 - off1 == imm_n1 - imm_n0 == len(new) == _counter(ret), 'off %d->%d imm %d->%d log %s' % (off0, off1, imm_n0, imm_n1, _counter(ret)))
    check('C8 Nemo level: the dispersal call itself mints nothing, recycles nothing, currentID unchanged', nemo0 == nemo_counts(F), str(nemo0))
    gidx = [h.global_index for h in hobj]
    check('C8 migrants keep their birth habitat: getHome() unchanged == global_index of the SOURCE habitat (!= destination)',
          all(x.nemo_indi_obj.getHome() == home0[id(x)] == gidx[i] != gidx[j] for i, j, x in new))
    # germination from the immigrant pool after the network dispersal
    rec_ids = []
    nemo1 = nemo_counts(F)
    gret = m3.meta_local_germinate_from_offspring_and_immigrant_pool()
    placed = {id(x): j for j, h in enumerate(hobj) for x in C.meta_common.micro(h)}
    germ = [(i, j, x) for i, j, x in new if id(x) in placed]
    still = [(i, j, x) for i, j, x in new if id(x) in {id(y) for h in hobj for y in h.immigrant_pool}]
    check('C8 germination after the network dispersal: %d immigrants germinate, each in its destination habitat, the rest stay in that habitat\'s immigrant_pool' % len(germ),
          len(germ) > 0 and all(placed[id(x)] == j for i, j, x in germ) and len(germ) + len(still) == len(new) and all(x.nemo_indi_obj.getHome() == gidx[i] for i, j, x in germ), gret.strip()[:120])
    check('C8 germination moves, it does not mint / recycle', nemo_counts(F)[:2] == nemo1[:2])
    try: C.check_world([m3, m1], fs, 'C8'); inv = True
    except AssertionError as e: inv = False; print(e)
    check('C8 world invariants after dispersal + germination (one owner, alive, minted - recycled == held, microsites)', inv)
    # clean-up recycles every individual exactly once
    for f in fs:
        real = f.recycle
        def counted(indi, real=real): rec_ids.append(id(indi)); return real(indi)
        f.recycle = counted
    held = C.held([m3, m1]); nheld = len(held)
    C.teardown([m3, m1], fs)
    check('C8 clean-up: every held individual recycled exactly once (%d), none twice, minted == recycled' % nheld,
          len(rec_ids) == nheld and len(set(rec_ids)) == nheld and set(rec_ids) == {id(x) for x in held} and all(f.minted_num == f.recycled_num for f in fs))
    for f in fs: del f.recycle


def equiv_worker(variant, n_seeds):
    ''' one process, variant 'new' (current code) or 'ref' (v3 + removal): N seeds of random worlds (2-4 patches x 1-3 habitats), 5 kernels x 5 rates;
        every 3rd seed forces the rate matrix to 0.9 off-diagonal (requests exceed the pool -> capping). Stops at the first broken invariant. '''
    seed_all(1); F = C.species(1); fs = list(F.values())
    cls = original_metacommunity(ext_names=[DISP], patches=REF_PATCH) if variant == 'ref' else metacommunity
    recs = []
    for s in range(n_seeds):
        rng = random.Random(s); seed_all(s)
        isl = C.make_islands(rng.choice((2, 3, 4)), rng.choice((1, 2, 3)), L=rng.choice((3, 4, 5)), W=rng.choice((3, 4)), K=rng.choice((0, 5, 20)), gidx0=100)
        isl.meta_initialize(F, C.TABLE); C.fill_pools(isl)
        isl.__class__ = cls
        hl = isl.global_habitat_id_ls; H = len(hl); hobj = [isl.set[p].set[h] for p, h in hl]
        forced = s % 3 == 2
        if forced:
            R = np.matrix(np.full((H, H), 0.9)); np.fill_diagonal(R, 0.1)
            isl.global_habitat_dispersal_among_rate_matrix = lambda *a, **k: R.copy()
        captured, real_ma = [], isl.matrix_around
        def ma(M):
            out = real_ma(M); captured.append(np.asarray(out).copy()); return out
        isl.matrix_around = ma
        origin = {id(x): i for i, h in enumerate(hobj) for x in h.offspring_pool}
        imm0 = {id(x) for h in hobj for x in h.immigrant_pool}
        nemo0 = nemo_counts(F)
        seed_all(10000 + s)
        k, kw = KERNELS[s % len(KERNELS)]
        ret = getattr(isl, DISP)(RATES[(s // 5) % 5], method=k, **kw)
        new = [(origin.get(id(x)), j, x) for j, h in enumerate(hobj) for x in h.immigrant_pool if id(x) not in imm0]
        moved = {id(x) for _, _, x in new}
        Rm = captured[0] if captured else np.zeros((H, H))
        req = int(sum(Rm[i, j] for i in range(H) for j in range(H) if i != j and hl[i][0] != hl[j][0]))
        rec = {'seed': s, 'forced': forced, 'H': H, 'digest': C.digest([isl]), 'rng': C.rng_digest(), 'ret': ret, 'counter': _counter(ret),
               'moved': [[i, j, x.individual_id] for i, j, x in new], 'requested': req, 'capped': req - len(new),
               'same_patch_or_unknown': sum(1 for i, j, _ in new if i is None or hl[i][0] == hl[j][0]),
               'in_two_pools': sum(1 for h in hobj for x in h.offspring_pool if id(x) in moved),
               'nemo_delta': [a - b for a, b in zip(nemo_counts(F), nemo0)]}
        try: C.check_world([isl], fs, 'equiv %s seed=%d' % (variant, s)); rec['invariants'] = 'ok'
        except AssertionError as e: rec['invariants'] = str(e)[:300]
        recs.append(rec)
        if rec['invariants'] != 'ok' or rec['in_two_pools']: break          # double ownership: a teardown would recycle twice
        C.teardown([isl], fs)
    return recs


def chain_worker(seed, steps):
    ''' mainland + 3x3 islands, pool schedule with the network offspring dispersal on every step (m_among 0.2 / 0.6 / 1.0, 5 kernels), invariants after every call '''
    seed_all(seed); F = C.species(seed); fs = list(F.values())
    ml = C.make_mainland(L=8); isl = C.make_islands(3, 3, L=4, K=10, gidx0=100)
    ml.meta_initialize(F, C.TABLE); isl.meta_initialize(F, C.TABLE)
    rain = C.rain_for(ml, 0.01)
    moved, n_ops, t0, name, t = 0, 0, time.time(), '', -1
    try:
        for t in range(steps):
            k, kw = KERNELS[t % len(KERNELS)]; m_among = (0.2, 0.6, 1.0)[t % 3]
            for name, op in C.pool_ops(ml, isl, rain, m_among, k, kw):
                ret = op(); n_ops += 1
                if name == DISP: moved += _counter(ret) or 0
                C.check_world([ml, isl], fs, 'chain seed=%d t=%d %s' % (seed, t, name))
        C.teardown([ml, isl], fs)
        return {'seed': seed, 'steps': steps, 'ok': True, 'ops': n_ops, 'network_migrants': moved, 'seconds': round(time.time() - t0, 1)}
    except BaseException as e:
        return {'seed': seed, 'steps': steps, 'ok': False, 'at': 't=%d %s' % (t, name), 'error': ('%s: %s' % (type(e).__name__, e))[:400], 'network_migrants': moved}


def _restore(hobj, snap):
    for h, (op, ip) in zip(hobj, snap): h.offspring_pool = list(op); h.immigrant_pool = list(ip)


def cstats_worker():
    ''' C5: expected immigrants per destination, current (move) vs v3 (copy), same world, pools restored after every call; C5b: which source individual is picked (uniform) '''
    seed_all(3); F = C.species(3); fs = list(F.values())
    m = C.make_islands(3, 2, L=6, K=0, gidx0=100); m.meta_initialize(F, C.TABLE)
    m.meta_dead_selection(0.2, 0.5, 'niche_gaussian'); m.meta_sex_reproduce_mutate_into_offspring_pool(1.0, None, None)
    hobj = [m.set[p].set[h] for p, h in m.global_habitat_id_ls]; H = len(hobj)
    snap = [(list(h.offspring_pool), list(h.immigrant_pool)) for h in hobj]
    pools = [len(h.offspring_pool) for h in hobj]
    V3 = original_metacommunity(ext_names=[DISP])
    origin = {id(x): (i, k) for i, h in enumerate(hobj) for k, x in enumerate(h.offspring_pool)}
    R = Q(1500, 5000); out = {'pools': pools, 'R': R, 'rates': {}}
    for rate in (0.1, 0.5, 0.9, 1.0):
        res = {}
        for variant, fn, s0 in (('current', lambda: getattr(m, DISP)(rate, method='exponential', rho=1), 100000), ('v3 copy', lambda: getattr(V3, DISP)(m, rate, method='exponential', rho=1), 200000),
                                ('current b', lambda: getattr(m, DISP)(rate, method='exponential', rho=1), 300000)):       # independent seeds per sample
            per_j = [[] for _ in range(H)]; pick = Counter()
            for r in range(R):
                seed_all(s0 + r)
                fn()
                for j, h in enumerate(hobj):
                    imm = h.immigrant_pool[len(snap[j][1]):]
                    per_j[j].append(len(imm))
                    for x in imm:
                        i, k = origin[id(x)]
                        if i == 0: pick[k] += 1
                _restore(hobj, snap)
            res[variant] = {'mean': [float(np.mean(v)) for v in per_j], 'per_j': per_j, 'pick': [pick[k] for k in range(pools[0])]}
        pv = [C.welch_p(res['current']['per_j'][j], res['v3 copy']['per_j'][j]) for j in range(H)]
        pc = [C.welch_p(res['current']['per_j'][j], res['current b']['per_j'][j]) for j in range(H)]
        st, df, pu = C.chi2_gof(res['current']['pick'], [sum(res['current']['pick']) / pools[0]] * pools[0])
        out['rates'][str(rate)] = {'mean_current': [round(x, 3) for x in res['current']['mean']], 'mean_v3': [round(x, 3) for x in res['v3 copy']['mean']],
                                   'p_current_vs_v3': [round(x, 4) for x in pv], 'p_current_vs_current': [round(x, 4) for x in pc], 'pick_uniform_p': round(pu, 4),
                                   'total_current': round(sum(res['current']['mean']), 3), 'total_v3': round(sum(res['v3 copy']['mean']), 3)}
    C.teardown([m], fs)
    return out


def ctrunc_worker():
    ''' C6: one source habitat (pool of n plain objects) and 8 destination habitats (4 patches x 2) at the same distance (uniform kernel): requested (matrix_around row)
        vs realized per destination, in loop order j = 1..8 '''
    spec = [('src', 0, (0.0, 0.0), [('h0', 0, (0.0, 0.0), 1)])] + [('d%d' % k, k + 1, (1.0, 0.0), [('h0', 0, (1.0, 0.0), 10 + 2 * k), ('h1', 1, (1.0, 0.0), 11 + 2 * k)]) for k in range(4)]
    m, _ = build_world(spec, metacommunity, 1)
    hobj = [m.set[p].set[h] for p, h in m.global_habitat_id_ls]; H = len(hobj)
    cap, real_ma = [], m.matrix_around
    def ma(M):
        out = real_ma(M); cap.append(np.asarray(out)[0].copy()); return out
    m.matrix_around = ma
    R = Q(1500, 6000); res = []
    for n in (1, 3, 10, 50):
        items = [object() for _ in range(n)]
        for rate in (0.01, 0.1, 0.5, 0.9, 1.0):
            exc = trunc = 0; short = [0] * (H - 1); req_tot = real_tot = 0
            for r in range(R):
                seed_all(r)
                hobj[0].offspring_pool = list(items)
                for h in hobj: h.immigrant_pool = []
                cap.clear(); getattr(m, DISP)(rate, method='uniform')
                req = cap[0][1:].astype(int); got = [len(hobj[j].immigrant_pool) for j in range(1, H)]
                req_tot += int(req.sum()); real_tot += sum(got)
                if req.sum() > n: exc += 1; trunc += int(req.sum()) - n
                for j in range(H - 1): short[j] += int(req[j]) - got[j]
            res.append({'pool': n, 'm': rate, 'P_excess': round(exc / R, 4), 'mean_truncated': round(trunc / R, 4), 'requested_mean': round(req_tot / R, 3), 'realized_mean': round(real_tot / R, 3),
                        'truncated_share_of_request': round(trunc / req_tot, 4) if req_tot else 0.0, 'shortfall_by_destination_order': [round(x / R, 4) for x in short]})
    return {'R': R, 'rows': res}


def cdist_worker():
    ''' C7: one source habitat (200 plain objects) and one destination habitat per patch at distances 1..6; per kernel the realized counts vs the rate matrix (chi2),
        exponential: slope of log(count) vs distance ~ -rho '''
    spec = [('src', 0, (0.0, 0.0), [('h0', 0, (0.0, 0.0), 1)])] + [('d%d' % k, k + 1, (float(k + 1), 0.0), [('h0', 0, (float(k + 1), 0.0), 10 + k)]) for k in range(6)]
    m, _ = build_world(spec, metacommunity, 1)
    hobj = [m.set[p].set[h] for p, h in m.global_habitat_id_ls]
    items = [object() for _ in range(200)]; R = Q(300, 1000); rate = 0.5
    kernels = KERNELS + [('exponential', {'rho': 0.5}), ('exponential', {'rho': 1.0})]
    rows = []
    for k, kw in kernels:
        cnt = np.zeros(6)
        for r in range(R):
            seed_all(r)
            hobj[0].offspring_pool = list(items)
            for h in hobj: h.immigrant_pool = []
            getattr(m, DISP)(rate, method=k, **kw)
            cnt += [len(h.immigrant_pool) for h in hobj[1:]]
        M = np.asarray(m.global_habitat_dispersal_among_rate_matrix(rate, method=k, **kw))[0, 1:]
        exp_ = R * 200 * M
        st, df, p = C.chi2_gof(cnt.tolist(), exp_.tolist())
        row = {'kernel': k, 'params': kw, 'observed': cnt.astype(int).tolist(), 'expected': [round(x, 1) for x in exp_], 'chi2': round(st, 2), 'df': df, 'p': round(p, 4),
               'total_obs': int(cnt.sum()), 'total_exp': round(float(exp_.sum()), 1)}
        if k == 'exponential':
            ok = cnt > 0
            d = np.arange(1, 7)[ok]; y = np.log(cnt[ok])
            row['slope'] = round(float(np.polyfit(d, y, 1)[0]), 4) if ok.sum() >= 2 else None
        rows.append(row)
    return {'R': R, 'pool': 200, 'm': rate, 'rows': rows}


def dispersal_main():
    out = {}
    dispersal_bookkeeping()
    jobs = [['worker', 'equiv', 'new', str(N_EQUIV)], ['worker', 'equiv', 'ref', str(N_EQUIV)]] + [['worker', 'chain', str(s), str(n)] for s, n in CHAIN] + \
           [['worker', 'cstats'], ['worker', 'ctrunc'], ['worker', 'cdist']]
    res = _pmap(jobs, max_workers=4)
    new, ref, chains, cst, ctr, cdi = res[0], res[1], res[2:5], res[5], res[6], res[7]
    # C3 equivalence
    if isinstance(ref, dict) or isinstance(new, dict):
        check('C3 worker run', False, str(new if isinstance(new, dict) else ref)[:500])
    else:
        check('C3 reference (v3 + removal): %d seeds, invariants hold, never two pools' % len(ref), len(ref) == N_EQUIV and all(r['invariants'] == 'ok' and not r['in_two_pools'] for r in ref))
        bad_new = [r for r in new if r['invariants'] != 'ok' or r['in_two_pools']]
        check('C3 current: %d seeds, invariants hold, never two pools' % len(new), not bad_new and len(new) == N_EQUIV, bad_new[:1])
        n = min(len(new), len(ref))
        same = {k: sum(a[k] == b[k] for a, b in zip(new[:n], ref[:n])) for k in ('digest', 'rng', 'ret', 'moved')}
        check('C3 current == v3 + removal on every one of %d seeds (state digest, random + np.random, log, who moved where in which order)' % N_EQUIV,
              len(new) == len(ref) == N_EQUIV and all(v == N_EQUIV for v in same.values()), same)
        check('C4 current: no same-patch migrant, counter == moved, no Nemo birth / recycle / currentID change (every seed)',
              all(r['same_patch_or_unknown'] == 0 and r['counter'] == len(r['moved']) and r['nemo_delta'] == [0, 0, 0] for r in new))
        nat = [r for r in ref if not r['forced']]; frc = [r for r in ref if r['forced']]
        out['capping_equiv'] = {'natural_seeds': len(nat), 'natural_seeds_capped': sum(r['capped'] > 0 for r in nat), 'natural_capped_total': sum(r['capped'] for r in nat),
                                'natural_requested_total': sum(r['requested'] for r in nat), 'forced_seeds': len(frc), 'forced_capped_total': sum(r['capped'] for r in frc)}
        check('C6 forced over-request in the equivalence worlds: capping exercised', out['capping_equiv']['forced_capped_total'] > 0, out['capping_equiv'])
        out['equiv_summary'] = {'seeds': N_EQUIV, 'moved_total': sum(len(r['moved']) for r in new)}
    # C5 v3 copy statistics
    if 'error' in cst: check('C5 worker run', False, cst['error'][:500])
    else:
        out['stats_vs_v3'] = cst
        alpha = 0.001 / len(cst['pools'])
        for rate in ('0.1', '0.5'):
            r = cst['rates'][rate]
            check('C5 m=%s: expected immigrants per destination, current (move) == v3 (copy) (Welch, Bonferroni alpha %.4g; %d calls each)' % (rate, alpha, cst['R']), min(r['p_current_vs_v3']) > alpha, r)
        for rate in ('0.9', '1.0'):
            r = cst['rates'][rate]
            check('C5 m=%s: current total <= v3 total (removal caps what v3 over-copied)' % rate, r['total_current'] <= r['total_v3'] + 1e-9, r)
        check('C5 sanity: the test detects nothing between two independent current samples (m=0.5)', min(cst['rates']['0.5']['p_current_vs_current']) > alpha)
        check('C5 which source individual migrates is uniform (chi2, habitat 0, m=0.5)', cst['rates']['0.5']['pick_uniform_p'] > 0.001, cst['rates']['0.5']['pick_uniform_p'])
    # C6 truncation
    if 'error' in ctr: check('C6 worker run', False, ctr['error'][:500])
    else:
        out['truncation'] = ctr
        r1 = [r for r in ctr['rows'] if r['m'] == 1.0]
        check('C6 truncation quantified (m x pool): requests beyond the pool are cut; realized <= pool always', all(r['realized_mean'] <= r['pool'] + 1e-9 for r in ctr['rows']), [(r['pool'], r['m'], r['P_excess']) for r in r1])
        late = [r for r in ctr['rows'] if r['mean_truncated'] > 0]
        bias = all(r['shortfall_by_destination_order'][-1] >= r['shortfall_by_destination_order'][0] for r in late)
        check('C6 shortfall falls on the LATER destinations of the j loop (order bias)', bias, [(r['pool'], r['m'], r['shortfall_by_destination_order']) for r in late[:3]])
        tab = {(r['pool'], r['m']): r for r in ctr['rows']}
        finding('known', DISP, 'matrix_around rounds each (i, j) separately: with high m the rounded row can exceed the offspring pool; sample_offspring_without_replacement caps with min() '
                'and the excess is cut from the destinations visited LAST in the j loop (v3 copied, never short). Agreed: not fixed, quantified.',
                evidence={'P(excess) m=1.0 by pool 1/3/10/50': [tab[(n, 1.0)]['P_excess'] for n in (1, 3, 10, 50)], 'P(excess) m=0.9': [tab[(n, 0.9)]['P_excess'] for n in (1, 3, 10, 50)],
                          'P(excess) m=0.5': [tab[(n, 0.5)]['P_excess'] for n in (1, 3, 10, 50)], 'mean truncated m=1.0': [tab[(n, 1.0)]['mean_truncated'] for n in (1, 3, 10, 50)],
                          'share of request truncated m=1.0': [tab[(n, 1.0)]['truncated_share_of_request'] for n in (1, 3, 10, 50)],
                          'shortfall by destination order (pool 3, m=1.0)': tab[(3, 1.0)]['shortfall_by_destination_order']},
                fix='(not requested) sample the whole row once (e.g. multinomial of min(pool, sum) over destinations, or shuffle the j order)')
    # C7 distance law
    if 'error' in cdi: check('C7 worker run', False, cdi['error'][:500])
    else:
        out['distance'] = cdi
        for r in cdi['rows']:
            check('C7 %s %s: migrants per distance ~ rate matrix (chi2 p=%.3g, %d migrants)' % (r['kernel'], r['params'], r['p'], r['total_obs']), r['p'] > 0.001, {k: r[k] for k in ('observed', 'expected')})
        for r in cdi['rows']:
            if r['kernel'] == 'exponential':
                rho = r['params']['rho']
                check('C7 exponential rho=%s: slope of log(migrants) vs distance = %s ~ -rho' % (rho, r.get('slope')), r.get('slope') is not None and abs(r['slope'] + rho) < 0.05 * rho + 0.05, r.get('slope'))
    for c in chains:
        check('C9 chain seed=%s x %s steps (pool schedule, network offspring dispersal every step): invariants after every call, teardown minted == recycled' % (c.get('seed'), c.get('steps')),
              c.get('ok'), json.dumps({k: v for k, v in c.items() if k not in ('seed', 'steps')})[:300])
    out['chain'] = chains
    return out


# ############################################################################
# PART: marker  (D)  offspring_marker_pool -> immigrant_marker_pool (both is_remove modes valid, no guard)
# ############################################################################
V3_MARK_SAMPLE = 'migrants_to_j_offspring_marker_ls += random.sample(h_i_object.offspring_marker_pool, migrants_i_j_num)'
MARK_REF_PATCH = {MARK: [(V3_MARK_SAMPLE, 'picked = random.sample(h_i_object.offspring_marker_pool, min(migrants_i_j_num, len(h_i_object.offspring_marker_pool))); '
                                          '[h_i_object.offspring_marker_pool.remove(x) for x in picked]; migrants_to_j_offspring_marker_ls += picked')]}
MARK_CHAIN = [(1, Q(200, 300), 'true'), (2, Q(200, 300), 'true'), (3, Q(200, 300), 'false')]


def _mark_cls(variant):
    ''' new -> current class; ref_true -> v3 + removal; v3 -> the v3 body unchanged (copy) '''
    if variant == 'ref_true': return original_metacommunity(ext_names=[MARK], patches=MARK_REF_PATCH)
    if variant == 'v3': return original_metacommunity(ext_names=[MARK])
    return metacommunity


def _mark_call(m, variant, mode, rate, k, kw):
    if variant == 'new': return getattr(m, MARK)(rate, method=k, is_remove=(mode == 'true'), **kw)
    return getattr(m, MARK)(rate, method=k, **kw)                                       # originals have no is_remove


def _marker_world(s, F, n_p=None, n_h=None, gidx0=100):
    rng = random.Random(s); seed_all(s)
    m = C.make_islands(n_p or rng.choice((2, 3, 4)), n_h or rng.choice((1, 2, 3)), L=rng.choice((3, 4, 5)), W=rng.choice((3, 4)), K=0, gidx0=gidx0)
    m.meta_initialize(F, C.TABLE)
    m.meta_dead_selection(0.3, 0.5, 'niche_gaussian')
    m.meta_sex_reproduce_calculation_with_offspring_marker_pool(RATES[(s // 5) % 5])
    m.meta_dispersal_within_patch_from_offspring_marker_to_immigrant_marker_pool(0.5)
    return m


def _marker_pools(m):
    hl = m.global_habitat_id_ls
    return [list(m.set[p].set[h].offspring_marker_pool) for p, h in hl], [list(m.set[p].set[h].immigrant_marker_pool) for p, h in hl]


def marker_bookkeeping():
    ''' D1 in-process: patch_num < 2 == v3 in both modes; True moves (multiset conserved), False copies; D2 germinate-and-birth after the network dispersal; D3 known ValueError risk '''
    V3 = _mark_cls('v3')
    seed_all(1); F = C.species(1); fs = list(F.values())
    m1 = _marker_world(7, F, 1, 2)
    for mode in ('true', 'false'):
        snap = C.full_snapshot([m1], fs)
        a = getattr(m1, MARK)(0.5, is_remove=(mode == 'true')); b = getattr(V3, MARK)(m1, 0.5)
        check('D1 patch_num < 2, is_remove=%s: same log as v3, nothing changed' % (mode == 'true'), a == b and C.same_snapshot(snap, C.full_snapshot([m1], fs)), a.strip()[:90])
    for mode, s in (('true', 8), ('false', 9)):
        m = _marker_world(s, F, 3, 2, gidx0=300)
        off0, imm0 = _marker_pools(m)
        nemo0 = nemo_counts(F)
        ret = getattr(m, MARK)(0.9, method='exponential', is_remove=(mode == 'true'), rho=2)
        off1, imm1 = _marker_pools(m)
        n = _mcounter(ret); new = sum(len(b) - len(a) for a, b in zip(imm0, imm1))
        appended = all(b[:len(a)] == a for a, b in zip(imm0, imm1))
        if mode == 'true':
            cons = Counter(x for p in off0 + imm0 for x in p) == Counter(x for p in off1 + imm1 for x in p)
            check('D1 is_remove=True: markers move (offspring -n, immigrant +n, n == log, multiset of markers conserved, old immigrants kept in order)',
                  n and new == n == sum(map(len, off0)) - sum(map(len, off1)) and cons and appended, 'n=%s new=%d' % (n, new))
        else:
            check('D1 is_remove=False: markers copy (offspring_marker pools untouched, immigrant +n, n == log)', n and off0 == off1 and new == n and appended, 'n=%s new=%d' % (n, new))
        check('D1 is_remove=%s: no Nemo birth / recycle / currentID change' % (mode == 'true'), nemo0 == nemo_counts(F))
        if m1 is not None: C.release([m1]); m1 = None              # m1 shares the factories: recycle it first (teardown without its balance assert, m is alive), else held(m) != minted - recycled
        C.check_world([m], fs, 'MB %s' % mode); C.teardown([m], fs)
    # D2 germinate-and-birth after the network marker dispersal: newborns in the right habitat, home == birth habitat's global_index
    stats = Counter(); bad = []
    for s in range(Q(12, 40)):
        m = _marker_world(100 + s, F, 3, 2, gidx0=500)
        getattr(m, MARK)(0.9, method='exponential', rho=1)
        hl = m.global_habitat_id_ls; hobj = [m.set[p].set[h] for p, h in hl]; gidx = {(p, h): m.set[p].set[h].global_index for p, h in hl}
        before = [Counter(h.offspring_marker_pool + h.immigrant_marker_pool) for h in hobj]
        micro0 = [{id(x) for x in C.meta_common.micro(h)} for h in hobj]
        births = {}
        for (p, hn), h in zip(hl, hobj):
            real = h.hab_sex_reproduce_mutate_with_num
            def wrapped(mr, pv, num, real=real, key=(p, hn)):
                out = real(mr, pv, num)                     # parents are drawn from this (birth) habitat's microsites
                for x in out: births[id(x)] = (key,)
                return out
            h.hab_sex_reproduce_mutate_with_num = wrapped
        try:
            m.meta_local_germinate_and_birth_from_offspring_marker_and_immigrant_marker_pool(None, None)
        except ValueError as e:
            stats['ValueError'] += 1
        after = [Counter(h.offspring_marker_pool + h.immigrant_marker_pool) for h in hobj]
        for j, h in enumerate(hobj):
            newborn = [x for x in C.meta_common.micro(h) if id(x) not in micro0[j]]
            consumed = before[j] - after[j]
            born_from = Counter((births[id(x)][0]) for x in newborn if id(x) in births)
            if born_from != Counter({(mk[0], mk[1]): c for mk, c in consumed.items()}): bad.append(('habitat', s, hl[j], dict(born_from), dict(consumed)))
            for x in newborn:
                key = births.get(id(x), (None,))[0]
                if key is None: bad.append(('not born by the marker loop', s, hl[j])); continue
                if x.nemo_indi_obj.getHome() != gidx[key]: bad.append(('home', s, hl[j], key, x.nemo_indi_obj.getHome(), gidx[key]))
                stats['newborn'] += 1; stats['newborn_from_other_patch'] += key[0] != hl[j][0]
        for h in hobj: del h.hab_sex_reproduce_mutate_with_num
        C.check_world([m], fs, 'D2'); C.teardown([m], fs)
    check('D2 germinate-and-birth after the network marker dispersal: each newborn sits in the habitat whose pool held its marker; consumed markers == newborns per habitat; '
          'getHome() == global_index of the PARENTS\' (birth) habitat (%d newborns, %d from another patch)' % (stats['newborn'], stats['newborn_from_other_patch']),
          not bad and stats['newborn'] > 0 and stats['newborn_from_other_patch'] > 0, bad[:3])
    # D3 known risk: the birth habitat of a marker has no female-male pair any more
    raised, births_before, n = 0, [], Q(15, 40)
    for s in range(n):
        m = _marker_world(200 + s, F, 3, 2, gidx0=700)
        getattr(m, MARK)(0.9, method='uniform')
        moved_from_p0 = sum(1 for p, h in m.global_habitat_id_ls if p != 'p0' for mk in m.set[p].set[h].immigrant_marker_pool if mk[0] == 'p0')
        m.set['p0'].patch_disturbance_process()
        a0 = m.get_meta_individual_num()
        try: m.meta_local_germinate_and_birth_from_offspring_marker_and_immigrant_marker_pool(None, None)
        except ValueError as e: raised += 1; births_before.append(m.get_meta_individual_num() - a0); err = str(e)
        C.check_world([m], fs, 'D3'); C.teardown([m], fs)
    check('D3 known risk reproduced: network marker dispersal -> disturbance of the birth patch -> germinate-and-birth raises ValueError part-way (%d / %d seeds), invariants hold afterwards' % (raised, n), raised > 0)
    finding('risk', 'meta_local_germinate_and_birth_from_offspring_marker_and_immigrant_marker_pool (+ network marker dispersal)',
            'a marker whose birth habitat has no female-male pair any more raises ValueError (random.sample on an empty pair list) part-way through the loop -- earlier markers already born. '
            'Reachable through the network dispersal when a disturbance / clear-up of the birth patch runs between the network marker dispersal and germination; not in the experiments\' order (disturbance after germination). Pre-existing (v3).',
            evidence='network marker dispersal (m=0.9) -> p0.patch_disturbance_process() -> germinate-and-birth: ValueError in %d / %d seeds, births before the raise %s' % (raised, n, births_before[:8]),
            fix='skip (or drop) a marker whose birth habitat has no parent pair; or clear the markers of a disturbed habitat in every immigrant_marker_pool')


def marker_equiv_worker(variant, mode, n_seeds):
    ''' one process: N marker worlds (2-4 patches x 1-3 habitats), 5 kernels x 5 rates, every 3rd seed the rate matrix forced to 0.9 off-diagonal '''
    seed_all(1); F = C.species(1); fs = list(F.values())
    cls = _mark_cls(variant); recs = []
    for s in range(n_seeds):
        m = _marker_world(s, F); m.__class__ = cls
        H = len(m.global_habitat_id_ls)
        if s % 3 == 2:
            R = np.matrix(np.full((H, H), 0.9)); np.fill_diagonal(R, 0.1)
            m.global_habitat_dispersal_among_rate_matrix = lambda *a, **k: R.copy()
        off0, imm0 = _marker_pools(m)
        seed_all(20000 + s)
        k, kw = KERNELS[s % len(KERNELS)]
        ret = _mark_call(m, variant, mode, RATES[(s // 5) % 5], k, kw)
        off1, imm1 = _marker_pools(m)
        recs.append({'seed': s, 'digest': C.digest([m]), 'rng': C.rng_digest(), 'ret': ret, 'counter': _mcounter(ret),
                     'new': [[list(x) for x in b[len(a):]] for a, b in zip(imm0, imm1)], 'off_removed': sum(map(len, off0)) - sum(map(len, off1)),
                     'conserved': Counter(x for p in off0 + imm0 for x in p) == Counter(x for p in off1 + imm1 for x in p)})
        C.check_world([m], fs, 'mequiv %s %s seed=%d' % (variant, mode, s)); C.teardown([m], fs)
    return recs


def marker_chain_worker(seed, steps, mode):
    ''' mainland + 3x3 islands, marker schedule with the network marker dispersal on every step (m_among 0.2 / 0.6 / 1.0, 5 kernels), invariants after every call '''
    seed_all(seed); F = C.species(seed); fs = list(F.values())
    ml = C.make_mainland(L=8); isl = C.make_islands(3, 3, L=4, K=0, gidx0=100)
    ml.meta_initialize(F, C.TABLE); isl.meta_initialize(F, C.TABLE)
    rain = C.rain_for(ml, 0.01)
    moved, n_ops, t0, name, t = 0, 0, time.time(), '', -1
    try:
        for t in range(steps):
            k, kw = KERNELS[t % len(KERNELS)]; m_among = (0.2, 0.6, 1.0)[t % 3]
            for name, op in C.marker_ops(ml, isl, rain, m_among, k, kw, is_remove=(mode == 'true')):
                ret = op(); n_ops += 1
                if name == MARK: moved += _mcounter(ret) or 0
                C.check_world([ml, isl], fs, 'mchain seed=%d t=%d %s' % (seed, t, name))
        C.teardown([ml, isl], fs)
        return {'seed': seed, 'steps': steps, 'mode': mode, 'ok': True, 'ops': n_ops, 'network_markers': moved, 'seconds': round(time.time() - t0, 1)}
    except BaseException as e:
        return {'seed': seed, 'steps': steps, 'mode': mode, 'ok': False, 'at': 't=%d %s' % (t, name), 'error': ('%s: %s' % (type(e).__name__, e))[:400], 'network_markers': moved}


def _onestep(seed, mode):
    f = os.path.join(tempfile.mkdtemp(prefix='ext_onestep_'), 'o.json')
    pr = subprocess.run([sys.executable, SIM_PY, 'onestep', '--scenario', 'grfe', '--scale', MODE, '--seed', str(seed), '--mode', mode, '--out', f], capture_output=True, text=True, cwd=HERE, env=PY_ENV)
    try: return json.load(open(f))
    except Exception: return {'error': (pr.stdout + pr.stderr)[-1500:]}


def marker_main():
    out = {}
    marker_bookkeeping()
    jobs = [('new', 'true'), ('ref_true', 'true'), ('new', 'false'), ('v3', 'false')]
    args = [['worker', 'mequiv', a, b, str(N_EQUIV)] for a, b in jobs] + [['worker', 'mchain', str(s), str(n), mode] for s, n, mode in MARK_CHAIN]
    res = _pmap(args, max_workers=4)
    eq = dict(zip(jobs, res[:4])); chains = res[4:]
    for (a, b), label in (((('new', 'true'), ('ref_true', 'true')), 'True == v3 + removal'), ((('new', 'false'), ('v3', 'false')), 'False == v3 (copy)')):
        x, y = eq[a], eq[b]
        if isinstance(x, dict) or isinstance(y, dict):
            check('D1 is_remove=%s worker run' % label, False, str(x if isinstance(x, dict) else y)[:500]); continue
        same = {k: sum(p[k] == q[k] for p, q in zip(x, y)) for k in ('digest', 'rng', 'ret', 'new')}
        check('D1 is_remove=%s on every one of %d seeds (state, random + np.random, log, which markers reach which pool in which order)' % (label, N_EQUIV),
              len(x) == len(y) == N_EQUIV and all(v == N_EQUIV for v in same.values()), same)
    if not isinstance(eq[('new', 'true')], dict):
        r = eq[('new', 'true')]
        check('D1 is_remove=True: total markers conserved and counter == removed from the offspring pools (every seed)', all(q['conserved'] and q['counter'] == q['off_removed'] for q in r))
    if not isinstance(eq[('v3', 'false')], dict):
        check('D1 reference sanity: v3 copy removes nothing', all(q['off_removed'] == 0 for q in eq[('v3', 'false')]))
    for c in chains:
        check('D5 chain seed=%s x %s steps, is_remove=%s: invariants after every call, teardown minted == recycled' % (c.get('seed'), c.get('steps'), c.get('mode')), c.get('ok'),
              json.dumps({k: v for k, v in c.items() if k not in ('seed', 'steps', 'mode')})[:300])
    out['chain'] = chains; out['equiv_summary'] = {'/'.join(j): (len(v) if isinstance(v, list) else v) for j, v in eq.items()}
    # D4 default False -> True on GRFE births (step 0 from the same state, one process per seed and mode)
    seeds = list(range(1, Q(5, 9)))
    from concurrent.futures import ThreadPoolExecutor
    with ThreadPoolExecutor(max_workers=4) as ex:
        one = {(s, md): f for s in seeds for md in ('true', 'false') for f in [ex.submit(_onestep, s, md)]}
        one = {k: f.result() for k, f in one.items()}
    errs = [v['error'] for v in one.values() if 'error' in v]
    if errs: check('D4 onestep runs', False, errs[0])
    else:
        same_state = all(one[(s, 'true')]['before'] == one[(s, 'false')]['before'] for s in seeds)
        check('D4 GRFE (%s scale, 1024 habitats, filled with nearest-species founders): identical state before the network marker dispersal in both modes (%d seeds)' % (MODE, len(seeds)), same_state)
        rows = [{'seed': s, 'births_true': one[(s, 'true')]['births'], 'births_false': one[(s, 'false')]['births'], 'network_markers': one[(s, 'true')]['network_markers'],
                 'offspring_markers_after_true': one[(s, 'true')]['offspring_markers_after'], 'offspring_markers_after_false': one[(s, 'false')]['offspring_markers_after'],
                 'immigrant_markers_after_true': one[(s, 'true')]['immigrant_markers_after'], 'immigrant_markers_after_false': one[(s, 'false')]['immigrant_markers_after'],
                 'empty_sites': one[(s, 'true')]['empty_sites']} for s in seeds]
        bt, bf = sum(r['births_true'] for r in rows), sum(r['births_false'] for r in rows)
        out['default_effect'] = {'rows': rows, 'births_true': bt, 'births_false': bf, 'ratio_true_over_false': round(bt / bf, 4) if bf else None}
        check('D4 births with is_remove=True <= births with False (fewer markers available) on every seed', all(r['births_true'] <= r['births_false'] for r in rows), out['default_effect']['ratio_true_over_false'])
        finding('info', MARK, 'default is_remove False -> True: with the default the network marker dispersal now MOVES markers (v3 copied: the source pool kept every marker and a copy travelled). '
                'GRFE step 0 from a metacommunity filled with founders (%s scale, %d seeds): births %d (True) vs %d (False), ratio %s; network markers per call %s. Long-run effect: see E (mT vs mF).' % (MODE, len(seeds), bt, bf, out['default_effect']['ratio_true_over_false'], [r['network_markers'] for r in rows]),
                evidence=rows[:4])
    return out


# ############################################################################
# PART: perf  (F)
# ############################################################################
def _variant_cls(name, inserts):
    ''' test-only subclass with the CURRENT body of `name`, plus one line inserted after each anchor line (same indentation); never written to the repo '''
    src = C.current_ext_source(); fn = C._functions(src)[name]
    lines = C._cut(src, fn, '    ').split('\n')
    for anchor, extra in inserts:
        k = [i for i, l in enumerate(lines) if l.strip() == anchor]
        assert len(k) == 1, anchor
        lines.insert(k[0] + 1, lines[k[0]][:len(lines[k[0]]) - len(lines[k[0]].lstrip())] + extra)
    ns = {'metacommunity': metacommunity, 'random': random, 'np': np}
    exec(compile('class variant_metacommunity(metacommunity):\n' + '\n'.join(lines) + '\n', '<test-only variant %s>' % name, 'exec'), ns)
    return ns['variant_metacommunity']


SKIP_ZERO = {OFF: [('migrants_i_j_num = int(global_hab_net_offspring_dispersal_matrix[i, j])', 'if migrants_i_j_num == 0: continue')],
             MARK: [('migrants_i_j_num = int(global_hab_net_offspring_marker_dispersal_matrix[i, j])', 'if migrants_i_j_num == 0: continue')]}


def perf_world(P, hp, pool, seed=1):
    side = int(np.ceil(np.sqrt(P))); hs = int(np.ceil(np.sqrt(hp)))
    spec = []; g = 100
    for k in range(P):
        px, py = (k // side) * hs, (k % side) * hs
        habs = []
        for h in range(hp):
            habs.append(('h%d' % h, h, (px + h // hs, py + h % hs), g)); g += 1
        spec.append(('p%d' % k, k, (px, py), habs))
    t = time.perf_counter(); m, _ = build_world(spec, metacommunity, seed); build_s = time.perf_counter() - t
    for p, h in m.global_habitat_id_ls:
        hb = m.set[p].set[h]
        hb.offspring_marker_pool = [(p, h, 'sexual')] * pool
        hb.offspring_pool = [object() for _ in range(pool)]
    return m, build_s, spec


def perf_main():
    out = {'rows': [], 'build': []}
    configs = [(4, 4), (16, 4), (64, 4), (16, 64)]
    V3 = original_metacommunity(ext_names=[OFF, MARK])
    REF = original_metacommunity(ext_names=[OFF], patches=REF_PATCH)
    SZ = {OFF: _variant_cls(OFF, SKIP_ZERO[OFF]), MARK: _variant_cls(MARK, SKIP_ZERO[MARK])}
    for P, hp in configs:
        H = P * hp
        m, build_s, spec = perf_world(P, hp, pool=20)
        Orig = original_metacommunity(ext_names=[C.ORIGINAL_REGISTRY], meta_names=['add_patch'])
        t = time.perf_counter(); build_world(spec, Orig, 1); build_v3 = time.perf_counter() - t
        out['build'].append({'H': H, 'patches': P, 'current_s': round(build_s, 3), 'v3_s': round(build_v3, 3)})
        hobj = [m.set[p].set[h] for p, h in m.global_habitat_id_ls]
        base = [(list(h.offspring_pool), list(h.offspring_marker_pool)) for h in hobj]
        def restore():
            for h, (op, mp) in zip(hobj, base): h.offspring_pool = list(op); h.offspring_marker_pool = list(mp); h.immigrant_pool = []; h.immigrant_marker_pool = []
        reps = {16: 20, 64: 8, 256: 3, 1024: 1}[H] * Q(1, 2)
        variants = [('marker current is_remove=True', lambda: getattr(m, MARK)(0.1, method='exponential', rho=2)),
                    ('marker current is_remove=False', lambda: getattr(m, MARK)(0.1, method='exponential', is_remove=False, rho=2)),
                    ('marker v3 (copy)', lambda: getattr(V3, MARK)(m, 0.1, method='exponential', rho=2)),
                    ('marker zero-skip variant (test-only)', lambda: getattr(SZ[MARK], MARK)(m, 0.1, method='exponential', rho=2)),
                    ('offspring current', lambda: getattr(m, OFF)(0.1, method='exponential', rho=2)),
                    ('offspring v3 + removal', lambda: getattr(REF, OFF)(m, 0.1, method='exponential', rho=2)),
                    ('offspring zero-skip variant (test-only)', lambda: getattr(SZ[OFF], OFF)(m, 0.1, method='exponential', rho=2))]
        res = {}
        for name, fn in variants:
            ts = []
            for r in range(reps):
                restore(); seed_all(r); t = time.perf_counter(); fn(); ts.append(time.perf_counter() - t)
            res[name] = round(float(np.median(ts)), 5)
        # zero-skip variant is behaviour-neutral (random.sample(x, 0) draws nothing)
        same = {}
        for fn_name, cur_fn, sz_fn in ((MARK, lambda: getattr(m, MARK)(0.1, method='exponential', rho=2), lambda: getattr(SZ[MARK], MARK)(m, 0.1, method='exponential', rho=2)),
                                       (OFF, lambda: getattr(m, OFF)(0.1, method='exponential', rho=2), lambda: getattr(SZ[OFF], OFF)(m, 0.1, method='exponential', rho=2))):
            sig = []
            for fn in (cur_fn, sz_fn):
                restore(); seed_all(77); ret = fn()
                sig.append((ret, C.rng_digest(), [len(h.offspring_pool) for h in hobj], [[id(x) for x in h.immigrant_pool] for h in hobj], [list(h.immigrant_marker_pool) for h in hobj]))
            same[fn_name[-30:]] = sig[0] == sig[1]
        # breakdown of the current marker call
        parts = defaultdict(float); calls = Counter()
        wrap = {'kernel strength matrix (dispersal_kernel_strength x H^2)': 'calculate_global_habitat_network_dispersal_kernel_strength_matrix',
                'normalize + rate matrix (incl. kernel)': 'global_habitat_dispersal_among_rate_matrix', 'matrix_around (H^2 np.random draws)': 'matrix_around',
                'pool-size matrix': 'get_global_habitat_network_offspring_marker_pool_num_matrix'}
        for label, attr in wrap.items():
            real = getattr(m, attr)
            def timed(*a, real=real, label=label, **k):
                t = time.perf_counter(); r = real(*a, **k); parts[label] += time.perf_counter() - t; return r
            setattr(m, attr, timed)
        nsample = [0, 0]
        for h in hobj:
            real = h.sample_offspring_without_replacement
            def counted(num, pool_name=('offspring_pool',), real=real):
                nsample[0] += 1; nsample[1] += num == 0; return real(num, pool_name)
            h.sample_offspring_without_replacement = counted
        restore(); seed_all(5); t = time.perf_counter(); getattr(m, MARK)(0.1, method='exponential', rho=2); tot = time.perf_counter() - t
        for attr in wrap.values(): delattr(m, attr)
        for h in hobj: del h.sample_offspring_without_replacement
        kern = parts['kernel strength matrix (dispersal_kernel_strength x H^2)']; rate_all = parts['normalize + rate matrix (incl. kernel)']
        brk = {'kernel strength matrix': kern, 'normalize + rate (rest)': rate_all - kern, 'matrix_around': parts['matrix_around (H^2 np.random draws)'],
               'pool-size matrix': parts['pool-size matrix']}
        brk['loop + sampling (rest)'] = tot - sum(brk.values())
        out['rows'].append({'H': H, 'patches': P, 'habitats_per_patch': hp, 'reps': reps, 'median_s': res, 'zero_skip_identical': same,
                            'breakdown_s': {k: round(v, 4) for k, v in brk.items()}, 'breakdown_share': {k: round(v / tot, 3) for k, v in brk.items()}, 'total_s': round(tot, 4),
                            'sample_calls': nsample[0], 'sample_calls_with_zero': nsample[1]})
        print('perf H=%d %s' % (H, res), flush=True)
    # pool-size sweep at H=256: the per-(i, j) rebuild of the source pool list makes the current call O(H^2 x pool)
    m, _, _ = perf_world(64, 4, pool=1)
    hobj = [m.set[p].set[h] for p, h in m.global_habitat_id_ls]
    out['pool_sweep'] = []
    for pool in (20, 100, 400):
        base = [[(p, h, 'sexual')] * pool for p, h in m.global_habitat_id_ls]
        res = {}
        for name, fn in (('marker current is_remove=True', lambda: getattr(m, MARK)(0.1, method='exponential', rho=2)), ('marker v3 (copy)', lambda: getattr(V3, MARK)(m, 0.1, method='exponential', rho=2)),
                         ('marker zero-skip variant (test-only)', lambda: getattr(SZ[MARK], MARK)(m, 0.1, method='exponential', rho=2))):
            ts = []
            for r in range(Q(2, 4)):
                for h, b0 in zip(hobj, base): h.offspring_marker_pool = list(b0); h.immigrant_marker_pool = []
                seed_all(r); t = time.perf_counter(); fn(); ts.append(time.perf_counter() - t)
            res[name] = round(float(np.median(ts)), 4)
        out['pool_sweep'].append({'H': 256, 'pool': pool, 'median_s': res, 'ratio_current_over_v3': round(res['marker current is_remove=True'] / res['marker v3 (copy)'], 2)})
        print('perf pool sweep', out['pool_sweep'][-1], flush=True)
    for r in out['rows']:
        check('F H=%d: zero-skip variant gives the identical result (state, RNG, log) -- the fix suggestion is behaviour-neutral' % r['H'], all(r['zero_skip_identical'].values()), r['zero_skip_identical'])
    big = out['rows'][-1]
    ratio = big['median_s']['marker current is_remove=True'] / big['median_s']['marker v3 (copy)']
    out['ratio_current_over_v3_H1024'] = round(ratio, 2)
    check('F one call recorded at H = 16 / 64 / 256 / 1024 for current, v3, reference and the test-only zero-skip variant', len(out['rows']) == 4)
    sw = out['pool_sweep']
    finding('risk' if sw[-1]['ratio_current_over_v3'] > 1.5 else 'info', '%s, %s' % (OFF, MARK),
            'performance: the migrated loop calls habitat.sample_offspring_without_replacement for EVERY (i, j) pair of different patches, also when migrants_i_j_num == 0, and each call rebuilds the '
            'whole (pool_name, obj) list of the source pool: O(H^2 x pool) per call (v3: random.sample(pool, 0) is O(1)). H=1024, pools of 20: %.2f s vs %.2f s for v3 (x%.2f); '
            'H=256, pools of 20 / 100 / 400: x%s vs v3. Skipping zero counts is bit-identical (random.sample(x, 0) draws nothing) and is faster than v3 (H=1024: %.2f s).'
            % (big['median_s']['marker current is_remove=True'], big['median_s']['marker v3 (copy)'], ratio, ' / x'.join(str(x['ratio_current_over_v3']) for x in sw), big['median_s']['marker zero-skip variant (test-only)']),
            evidence={'H=%d' % r['H']: {k: r['median_s'][k] for k in ('marker current is_remove=True', 'marker v3 (copy)', 'marker zero-skip variant (test-only)', 'offspring current', 'offspring zero-skip variant (test-only)')} for r in out['rows']} | {'pool sweep H=256': sw},
            fix='in both network dispersals: `if migrants_i_j_num == 0: continue` before sampling (identical results), or sample once per source row')
    shares = big['breakdown_share']
    finding('info', 'global habitat network (H=1024)', 'cost of one marker call at H=1024: %s; the kernel matrix is a Python double loop of H^2 dispersal_kernel_strength calls and matrix_around draws H^2 single uniforms; both are recomputed on every call although the distance matrix is fixed.'
            % ', '.join('%s %.0f%%' % (k, 100 * v) for k, v in shares.items()),
            evidence=big['breakdown_s'], fix='cache the normalized kernel per (method, kwargs); vectorize the kernel (numpy on the distance matrix) and matrix_around (one uniform array)')
    b = out['build'][-1]
    finding('info', 'global_habitat_network_registry -> update_global_habitat_distance_matrix', 'every add_patch recomputes the full H x H distance matrix in a Python double loop (same as v3): building GRFE (16 patches, 1024 habitats) costs %.1f s; incremental_update_global_habitat_distance_matrix exists but is unused.' % b['current_s'],
            evidence=out['build'], fix='use the incremental update in the registry (or one numpy computation after the last add_patch)')
    return out


# ############################################################################
# PART: fuzz  (G)
# ############################################################################
ALLOWED = {'colonize: rain > mainland population': ('ValueError', 'Sample larger than population'),
           'germinate-and-birth: marker birth habitat without a parent pair (known risk D3)': ('ValueError', 'Sample larger than population')}


def fuzz_worker(seed, n_runs, n_ops):
    rng = random.Random(seed); out = {'ops': 0, 'kinds': Counter(), 'allowed': Counter(), 'guards': 0, 'runs': 0, 'errors': []}
    seed_all(seed); F = C.species(seed); fs = list(F.values())
    for run in range(n_runs):
        seed_all(seed * 1000 + run)
        ml = C.make_mainland(L=rng.choice((3, 6))); isl = C.make_islands(rng.choice((2, 3, 4)), rng.choice((1, 2, 3)), L=rng.choice((3, 4)), K=rng.choice((0, 5)), gidx0=100)
        ml.meta_initialize(F, C.TABLE); isl.meta_initialize(F, C.TABLE)
        def kern(): return KERNELS[rng.randrange(len(KERNELS))]
        def net_off():
            k, kw = kern(); return getattr(isl, OFF)(rng.choice(RATES), method=k, **kw)
        def net_mark(flag):
            def f():
                k, kw = kern(); return getattr(isl, MARK)(rng.choice(RATES), method=k, is_remove=flag, **kw)
            return f
        def guard():
            snap = C.full_snapshot([ml, isl], fs)
            k, kw = kern()
            try: getattr(isl, OFF)(rng.choice(RATES), method=k, is_remove=False, **kw); raise AssertionError('no NotImplementedError')
            except NotImplementedError: pass
            assert C.same_snapshot(snap, C.full_snapshot([ml, isl], fs)), 'guarded call changed something'
            out['guards'] += 1
        def getters():
            k, kw = kern(); r = rng.choice(RATES)
            with np.errstate(all='ignore'):
                for g in ('get_global_habitat_network_offspring_pool_num_matrix', 'get_global_habitat_network_offspring_marker_pool_num_matrix', 'get_global_habitat_network_dormancy_pool_num_matrix', 'get_global_habitat_network_empty_sites_num_matrix'):
                    getattr(isl, g)()
                isl.get_global_habitat_network_emigrants_matrix(r, method=k, **kw); isl.get_global_habitat_network_immigrants_matrix(r, method=k, **kw)
                isl.get_global_habitat_network_disp_among_num_matrix(r, method=k, **kw)
        def among(fn):
            def f():
                k, kw = kern(); return getattr(isl, fn)(rng.choice(RATES), method=k, **kw)
            return f
        ops = {'mainland dead + birth': lambda: (ml.meta_dead_selection(0.2, 0.5, 'niche_gaussian'), ml.meta_mainland_sexual_birth_mutate_germinate(1.0, None, None)),
               'dead selection': lambda: isl.meta_dead_selection(rng.choice((0.05, 0.2, 0.5)), 0.5, 'niche_gaussian'),
               'reproduce into offspring_pool': lambda: isl.meta_sex_reproduce_mutate_into_offspring_pool(rng.choice((0.5, 1.0)), None, None),
               'offspring markers': lambda: isl.meta_sex_reproduce_calculation_with_offspring_marker_pool(rng.choice((0.5, 1.0))),
               'colonize': lambda: isl.meta_colonize_from_propagules_rains(ml, rng.choice((0.5, 2.0, 6.4))),
               'pairwise colonize': lambda: isl.pairwise_sexual_colonization_from_prpagules_rains(ml, rng.choice((0.5, 2.0))),
               'within offspring->immigrant': lambda: isl.meta_dispersal_within_patch_from_offspring_to_immigrant_pool(rng.random()),
               'within offspring+dormancy->immigrant': lambda: isl.meta_dispersal_within_patch_from_offspring_and_dormancy_to_immigrant_pool(rng.random()),
               'within offspring+dormancy (sites)': lambda: isl.meta_dispersal_within_patch_from_offspring_and_dormancy_pool(rng.random()),
               'within markers': lambda: isl.meta_dispersal_within_patch_from_offspring_marker_to_immigrant_marker_pool(rng.random()),
               'patch-level among offspring->immigrant': among('dispersal_aomng_patches_from_offspring_pool_to_immigrant_pool'),
               'patch-level among markers': among('dispersal_aomng_patches_from_offspring_marker_pool_to_immigrant_marker_pool'),
               'patch-level among offspring+dormancy->immigrant': among('dispersal_among_patches_from_offsrping_pool_and_dormancy_pool_to_immigrant_pool'),
               'patch-level among offspring+dormancy (sites)': among('dispersal_among_patches_from_offspring_pool_and_dormancy_pool'),
               'NETWORK offspring (is_remove=True)': net_off, 'NETWORK markers is_remove=True': net_mark(True), 'NETWORK markers is_remove=False': net_mark(False),
               'NETWORK offspring is_remove=False (guard)': guard, 'network getters': getters,
               'germinate offspring+immigrant': lambda: isl.meta_local_germinate_from_offspring_and_immigrant_pool(),
               'germinate offspring+dormancy': lambda: isl.meta_local_germinate_from_offspring_and_dormancy_pool(),
               'germinate offspring+immigrant+dormancy': lambda: isl.meta_local_germinate_from_offspring_immigrant_and_dormancy_pool(),
               'germinate-and-birth from markers': lambda: isl.meta_local_germinate_and_birth_from_offspring_marker_and_immigrant_marker_pool(None, None),
               'dormancy offspring->dormancy': lambda: isl.meta_dormancy_process_from_offspring_pool_to_dormancy_pool(),
               'dormancy offspring+immigrant': lambda: isl.meta_dormancy_process_from_offspring_pool_and_immigrant_pool(),
               'disturbance patches': lambda: isl.meta_disturbance_process_in_patches(0.3),
               'disturbance habitats': lambda: isl.meta_disturbance_process_in_habitat(0.3),
               'clear offspring': lambda: isl.meta_clear_up_offspring_pool(), 'clear immigrant': lambda: isl.meta_clear_up_immigrant_pool(),
               'clear dormancy': lambda: isl.meta_clear_up_dormancy_pool(), 'clear offspring+immigrant': lambda: isl.meta_clear_up_offspring_and_immigrant_pool(),
               'clear markers': lambda: isl.meta_clear_up_offspring_marker_and_immigrant_marker_pool()}
        weights = {k: (4 if k.startswith('NETWORK') else 1) for k in ops}
        names = list(ops)
        for k in range(n_ops):
            name = rng.choices(names, weights=[weights[n] for n in names])[0]
            try:
                ops[name]()
            except ValueError as e:
                if name in ('colonize', 'pairwise colonize') and 'Sample larger than population' in str(e): out['allowed']['colonize: rain > mainland population'] += 1
                elif name == 'germinate-and-birth from markers' and 'Sample larger than population' in str(e): out['allowed']['germinate-and-birth: marker birth habitat without a parent pair (known risk D3)'] += 1
                else: out['errors'].append('run %d op %d %s: %s: %s' % (run, k, name, type(e).__name__, e)); break
            except Exception as e:
                out['errors'].append('run %d op %d %s: %s: %s' % (run, k, name, type(e).__name__, str(e)[:200])); break
            out['ops'] += 1; out['kinds'][name] += 1
            try: C.check_world([ml, isl], fs, 'fuzz run %d op %d %s' % (run, k, name))
            except AssertionError as e: out['errors'].append(str(e)[:300]); break
        if out['errors']: break
        C.teardown([ml, isl], fs); out['runs'] += 1
    return out


def fuzz_main():
    n_proc, runs, ops = 3, Q(10, 30), 110
    res = _pmap([['worker', 'fuzz', str(s), str(runs), str(ops)] for s in range(1, n_proc + 1)], max_workers=n_proc)
    errs = [r['error'] for r in res if 'error' in r and isinstance(r['error'], str)]
    if errs: check('G worker run', False, errs[0][:500]); return {'fuzz': res}
    tot = sum(r['ops'] for r in res); errors = [e for r in res for e in r['errors']]
    kinds = Counter(); allowed = Counter()
    for r in res: kinds.update(r['kinds']); allowed.update(r['allowed'])
    check('G %d random operations (%d kinds, %d processes x %d runs x %d ops): invariants after each, teardown minted == recycled after every run, no unexpected exception' % (tot, len(kinds), n_proc, runs, ops),
          not errors and tot >= Q(3000, 9000), errors[:3])
    check('G guarded is_remove=False calls (offspring network): NotImplementedError + full snapshot unchanged', sum(r['guards'] for r in res) > 0, sum(r['guards'] for r in res))
    return {'ops': tot, 'kinds': dict(kinds), 'allowed': dict(allowed), 'guards': sum(r['guards'] for r in res), 'errors': errors}


# ############################################################################
# PART: sim  (E)  -> simulate_ext.py
# ############################################################################
def sim_jobs():
    J = []
    def add(name, *a): J.append((name, list(a)))
    if MODE == 'quick':
        for g in ('mT', 'mF', 'pool', 'mPatch'):
            for s in (1, 2, 3): add('net_%s_s%d' % (g, s), '--scenario', 'net', '--group', g, '--steps', '300', '--seed', str(s), '--scale', 'quick')
            add('grfe_%s_s1' % g, '--scenario', 'grfe', '--group', g, '--steps', '20', '--seed', '1', '--scale', 'quick', '--check-every', '10')
        rep = (('net', '100', '50'), ('grfe', '5', '5')); eqv = (('net', '300', '100'), ('grfe', '10', '5')); dist = ('net', 'quick', '300')
    else:
        for g in ('mT', 'mF', 'pool', 'mPatch'):
            for sc in ('net', 'grfe'): add('%s_%s_s1' % (sc, g), '--scenario', sc, '--group', g, '--steps', '2000', '--seed', '1', '--scale', 'full')
        rep = (('net', '500', '100'), ('grfe', '30', '10')); eqv = (('net', '1000', '100'), ('grfe', '50', '10')); dist = ('net', 'full', '300')
    for sc, st, ce in rep:
        for k in ('a', 'b'): add('repro_%s_%s' % (sc, k), '--scenario', sc, '--group', 'mT', '--steps', st, '--seed', '7', '--scale', 'quick', '--digest-every', '1', '--check-every', ce)
    for sc, st, ce in eqv:
        for v in ('cur', 'v3'): add('equiv_%s_%s' % (sc, v), '--scenario', sc, '--group', 'mF', '--steps', st, '--seed', '8', '--scale', 'quick', '--variant', v, '--digest-every', '1', '--check-every', ce)
    for g in ('mT', 'pool'):
        add('disturb_%s' % g, '--scenario', dist[0], '--group', g, '--steps', dist[2], '--seed', '9', '--scale', dist[1], '--dist', '0.05', '--watch-dist', '--check-every', '50')
    main = [j for j in J if not j[0].startswith(('repro_', 'equiv_', 'disturb_'))]
    return [j for j in J if j not in main] + main                       # short runs first: the long ones would otherwise hold every slot for hours


def sim_launch(max_procs):
    os.makedirs(os.path.join(HERE, 'results'), exist_ok=True)
    jobs = sim_jobs(); procs = {}; pending = list(jobs); done = {}
    for f in glob.glob(os.path.join(HERE, 'results', 'sim_*.json')): os.remove(f)
    logd = tempfile.mkdtemp(prefix='ext_sim_logs_')
    while pending or procs:
        while pending and len(procs) < max_procs:
            name, a = pending.pop(0)
            outf = os.path.join(HERE, 'results', 'sim_%s.json' % name)
            lf = open(os.path.join(logd, name + '.log'), 'w')
            procs[name] = (subprocess.Popen([sys.executable, SIM_PY, 'run'] + a + ['--out', outf], cwd=HERE, env=PY_ENV, stdout=lf, stderr=subprocess.STDOUT), lf, time.time())
            print('sim start %s' % name, flush=True)
        for name in list(procs):
            p, lf, t = procs[name]
            if p.poll() is not None:
                lf.close(); done[name] = {'rc': p.returncode, 's': round(time.time() - t, 1), 'log_tail': open(lf.name).read()[-1500:]}
                print('sim done %s rc=%d %.0fs' % (name, p.returncode, time.time() - t), flush=True)
                del procs[name]
        time.sleep(2)
    return done


def _load_sim(name):
    try: return json.load(open(os.path.join(HERE, 'results', 'sim_%s.json' % name)))
    except Exception: return None


def _half_mean(xs): xs = xs[len(xs) // 2:]; return float(np.mean(xs)) if xs else float('nan')


def _decay(r):
    ''' per-pair immigrant flow by distance bin + slope of log(flow) vs distance (bins with flow) '''
    fb = {int(k): v for k, v in r.get('flow_bins', {}).items()}; pairs = {int(k): v for k, v in r.get('pairs_by_distance', {}).items()}
    pts = sorted((d, fb[d] / pairs[d]) for d in fb if pairs.get(d) and fb[d] > 0)
    slope = float(np.polyfit([d for d, _ in pts], np.log([f for _, f in pts]), 1)[0]) if len(pts) >= 2 else None
    return pts, slope


def sim_analyze(launch):
    out = {'launch': launch, 'runs': {}}
    names = [n for n, _ in sim_jobs()]
    runs = {n: _load_sim(n) for n in names}
    for n, r in runs.items():
        if r is None: check('E run %s produced its result file' % n, False, launch.get(n, {}).get('log_tail', '')[-400:]); continue
        out['runs'][n] = {k: r.get(k) for k in ('status', 'error', 'at', 'H', 'patches', 'sites', 'build_s', 'total_s', 'steps_done', 'teardown_ok', 'teardown_balance', 'phase_seconds')}
        out['runs'][n]['rows'] = r.get('rows', [])
        s = r.get('series', {})
        out['runs'][n]['summary'] = {'alive_end': (s.get('alive') or [None])[-1], 'occupancy_end': r['rows'][-1]['occupancy'] if r.get('rows') else None,
                                     'richness_end': (s.get('richness') or [None])[-1], 'births_mean_2nd_half': _half_mean(s.get('births', [])),
                                     'migrants_mean_2nd_half': _half_mean(s.get('migrants', [])), 'alive_mean_2nd_half': _half_mean(s.get('alive', [])),
                                     'richness_mean_2nd_half': _half_mean(s.get('richness', [])), 'step_s_mean': float(np.mean(s['step_s'])) if s.get('step_s') else None,
                                     'step_s_max': float(np.max(s['step_s'])) if s.get('step_s') else None,
                                     'rss_first_last': [r['rows'][0]['rss_mb'], r['rows'][-1]['rss_mb']] if r.get('rows') else None}
        pts, slope = _decay(r); out['runs'][n]['decay'] = {'points': [[d, round(f, 5)] for d, f in pts], 'slope': slope}
        ok = r.get('status') == 'done' and r.get('teardown_ok') and r.get('gidx_distinct') and all(x.get('invariants') == 'ok' and x['minted'] - x['recycled'] == x['held'] for x in r.get('rows', []))
        check('E1 %s: %s steps, distinct global_index per habitat, invariants every %s steps (one owner, alive, minted - recycled == held, microsites), teardown minted == recycled' % (n, r.get('steps_done'), r['args'].get('check_every')),
              ok, r.get('error') or r.get('at') or '')
    # E2 reproducibility
    for sc in ('net', 'grfe'):
        a, b = runs.get('repro_%s_a' % sc), runs.get('repro_%s_b' % sc)
        if a and b:
            check('E2 reproducible: %s mT seed 7, two fresh processes, digest of the whole state + RNG after each of %d steps identical' % (sc, len(a['digests'])),
                  a['digests'] == b['digests'] and len(a['digests']) == a['args']['steps'], 'first mismatch %s' % next((x[0] for x, y in zip(a['digests'], b['digests']) if x != y), None))
    # E3 whole-run equivalence with v3 (is_remove=False)
    for sc in ('net', 'grfe'):
        a, b = runs.get('equiv_%s_cur' % sc), runs.get('equiv_%s_v3' % sc)
        if a and b:
            mm = next((x[0] for x, y in zip(a['digests'], b['digests']) if x != y), None)
            check('E3 whole run, %s marker flow, is_remove=False: current == v3 extension bodies (registry + add_patch + marker dispersal) after every one of %d steps' % (sc, len(a['digests'])),
                  a['digests'] == b['digests'] and len(a['digests']) == a['args']['steps'], 'first mismatch t=%s' % mm)
    # E4 ecology: network vs patch-level comparator; E5 decay; E6 mT vs mF
    eco = {}
    for sc in ('net', 'grfe'):
        for g in ('mT', 'mF', 'pool', 'mPatch'):
            rs = [out['runs'][n] for n in out['runs'] if n.startswith('%s_%s_s' % (sc, g)) and out['runs'][n]['status'] == 'done']
            if not rs: continue
            sm = lambda k: [x['summary'][k] for x in rs]
            eco['%s/%s' % (sc, g)] = {'runs': len(rs), 'alive_mean_2nd_half': round(float(np.mean(sm('alive_mean_2nd_half'))), 1), 'occupancy_end': round(float(np.mean(sm('occupancy_end'))), 4),
                                      'richness_mean_2nd_half': round(float(np.mean(sm('richness_mean_2nd_half'))), 2), 'births_mean_2nd_half': round(float(np.mean(sm('births_mean_2nd_half'))), 2),
                                      'among_patch_migrants_mean_2nd_half': round(float(np.mean(sm('migrants_mean_2nd_half'))), 3),
                                      'decay_slope': [x['decay']['slope'] for x in rs], 'step_s_mean': round(float(np.mean(sm('step_s_mean'))), 4)}
    out['ecology'] = eco
    for sc in ('net', 'grfe'):
        a, b = eco.get('%s/mT' % sc), eco.get('%s/mPatch' % sc)
        if a and b:
            rel = {k: (round(abs(a[k] - b[k]) / max(abs(b[k]), 1e-9), 3)) for k in ('alive_mean_2nd_half', 'richness_mean_2nd_half', 'among_patch_migrants_mean_2nd_half')}
            out.setdefault('network_vs_patch', {})[sc] = rel
            steps = min(runs[n]['steps_done'] for n in runs if runs[n] and n.startswith(('%s_mT_s' % sc, '%s_mPatch_s' % sc)))
            if steps < 200:
                out['network_vs_patch'][sc + ' (not assessed)'] = 'runs of %d steps from an empty metacommunity (%.0f alive): too short to compare' % (steps, a['alive_mean_2nd_half'])
            else:
                if a['decay_slope'] and b['decay_slope'] and None not in a['decay_slope'] + b['decay_slope']:
                    finding('info', '%s vs dispersal_aomng_patches_from_offspring_marker_pool_to_immigrant_marker_pool' % MARK,
                            '%s: same totals, different spatial pattern -- with the same rho the network kernel acts on habitat-to-habitat distances, the patch-level one on patch-centre distances: '
                            'log-slope of immigrants per habitat pair vs distance %s (network) vs %s (patch-level); alive / richness / among-patch migrants differ by %s' %
                            (sc, [round(x, 2) for x in a['decay_slope']], [round(x, 2) for x in b['decay_slope']], rel))
                close = rel['alive_mean_2nd_half'] < 0.15 and rel['richness_mean_2nd_half'] < 0.15
                check('E4 %s: network marker dispersal vs the commented patch-level alternative -- alive / richness within 15%% (2nd half means)' % sc, close, rel)
                if not close or rel['among_patch_migrants_mean_2nd_half'] > 0.25:
                    finding('info', '%s vs dispersal_aomng_patches_from_offspring_marker_pool_to_immigrant_marker_pool' % MARK,
                            '%s: network vs patch-level marker dispersal differ by %s (relative; 2nd-half means, network = habitat-to-habitat distances with the same rho, patch-level = patch-centre distances)' % (sc, rel),
                            evidence={'network': a, 'patch': b})
        for g in ('mT', 'mF', 'pool'):
            e = eco.get('%s/%s' % (sc, g))
            if e and any(s is not None for s in e['decay_slope']):
                rho = 2 if sc == 'net' else 0.2
                sl = [s for s in e['decay_slope'] if s is not None]
                check('E5 %s/%s: immigrants per habitat pair decay with distance (slope of log flow %s, kernel rho %s)' % (sc, g, [round(s, 3) for s in sl], rho), all(s < 0 for s in sl), sl)
        a, b = eco.get('%s/mT' % sc), eco.get('%s/mF' % sc)
        if a and b:
            out.setdefault('true_vs_false', {})[sc] = {'births_ratio': round(a['births_mean_2nd_half'] / b['births_mean_2nd_half'], 4) if b['births_mean_2nd_half'] else None,
                                                       'alive_ratio': round(a['alive_mean_2nd_half'] / b['alive_mean_2nd_half'], 4) if b['alive_mean_2nd_half'] else None,
                                                       'migrants_ratio': round(a['among_patch_migrants_mean_2nd_half'] / b['among_patch_migrants_mean_2nd_half'], 4) if b['among_patch_migrants_mean_2nd_half'] else None}
    # E7 disturbance
    for g in ('mT', 'pool'):
        r = runs.get('disturb_%s' % g)
        if not r: continue
        ev = r.get('disturbances', [])
        rec = []
        for e in ev:
            target = 0.5 * e['pre_alive']
            k = next((i + 1 for i, v in enumerate(e['recolonize']) if v >= target), None) if e['pre_alive'] else 0
            rec.append(k)
        out.setdefault('disturbance', {})[g] = {'events': len(ev), 'recycled_ok': sum(e['recycled_ok'] for e in ev), 'empty_after': sum(e['empty_after'] for e in ev),
                                                'others_untouched': sum(e['others_untouched'] for e in ev), 'dormancy_kept': sum(e['dormancy_kept'] for e in ev),
                                                'recycled_total': sum(e['recycled'] for e in ev), 'steps_to_50pct_recolonized': rec,
                                                'not_recolonized_by_end': sum(1 for k in rec if k is None)}
        check('E7 high disturbance (%s, patch_dist_rate 0.05): %d events; each recycles exactly what the patch held (microsites + pools), leaves it empty, keeps dormancy, other patches untouched'
              % (g, len(ev)), len(ev) > 0 and all(e['recycled_ok'] and e['empty_after'] and e['others_untouched'] and e['dormancy_kept'] for e in ev), out['disturbance'][g])
    return out


def sim_main():
    launch = sim_launch(Q(6, 8))
    return sim_analyze(launch)


# ############################################################################
# PART: baseline  (Step 0) -- the other suites, run with their folders backed up and restored byte-identical
# ############################################################################
RA_NOISE = ('setting random seed from input value', 'outputs: ', 'the microsite in the habitat is occupied.', 'the microsite in the habitat is empty.')
SUITE_ENV = dict(os.environ, PYTHONDONTWRITEBYTECODE='1', MPLBACKEND='Agg')
SUITE_ENV.pop('EXT_MODE', None)


def ra_run(script, cwd, args=(), keep=40):
    t = time.time()
    p = subprocess.run([sys.executable, os.path.join(cwd, script)] + list(args), capture_output=True, text=True, cwd=cwd, env=SUITE_ENV)
    lines = [l for l in (p.stdout + '\n' + p.stderr).splitlines() if l.strip() and not any(n in l for n in RA_NOISE) and not l.startswith('[FINDING')]
    return {'rc': p.returncode, 's': round(time.time() - t, 1), 'tail': lines[-keep:]}


def _manifest(d):
    ''' {relative path: md5} of every file of the folder (results/ and __pycache__ are untracked: git status alone cannot prove the restore) '''
    return {os.path.relpath(os.path.join(r, f), d): hashlib.md5(open(os.path.join(r, f), 'rb').read()).hexdigest() for r, _, fs in os.walk(d) for f in fs}


def _report_rows(rep, tail):
    ''' per-script results from a suite's fresh TEST_REPORT.md (table rows '| x.py | ... | PASS | 1.2 s |'), else from its stdout '''
    rows = re.findall(r'^\| (\S+\.py) \| .*? \| (PASS|\*\*FAIL\*\*) \| ([\d.]+) s \|$', rep, re.M)
    sm = re.search(r'^## Smoke: (PASS|\*\*FAIL\*\*)', rep, re.M)
    if sm: rows = [('smoke.py', sm.group(1), '')] + rows
    tests = re.search(r'^## Tests: .*?\n\n(.*?)\n\n', rep, re.M | re.S)              # metacommunity style: '| bookkeeping | covers | PASS | 2.1 s |'
    if tests and not re.findall(r'^\| (\S+\.py) \| ', tests.group(1), re.M):
        rows = [(m.group(1), m.group(2), m.group(3)) for m in re.finditer(r'^\| (\w+) \| .*? \| (PASS|\*\*FAIL\*\*) \| ([\d.]+) s \|$', tests.group(1), re.M)]
        if sm: rows = [('smoke', sm.group(1), '')] + rows
    if not rows:
        rows = [(m.group(1), m.group(2), m.group(3)) for l in tail for m in [re.match(r'(\S+\.py)\s+(PASS|FAIL)\s+([\d.]+)s', l)] if m]
    return [(a, 'FAIL' if 'FAIL' in b else 'PASS', c) for a, b, c in rows]


def _run_suite_restored(rel, script, args, tag):
    ''' run <rel>/<script> with a backup of the whole folder (taken AFTER the Step 0 edits), restore it afterwards; keep its fresh report as results/baseline_<tag>_TEST_REPORT.md '''
    d = os.path.join(REPO, rel)
    before, git_before = _manifest(d), C.git('status', '--short', '--', rel)
    tmp = tempfile.mkdtemp(prefix='suite_backup_')
    shutil.copytree(d, os.path.join(tmp, 'suite'), copy_function=shutil.copy2)
    try:
        r = ra_run(script, d, args)
        rep = open(os.path.join(d, 'TEST_REPORT.md')).read() if os.path.exists(os.path.join(d, 'TEST_REPORT.md')) else ''
    finally:
        shutil.rmtree(d); shutil.copytree(os.path.join(tmp, 'suite'), d, copy_function=shutil.copy2); shutil.rmtree(tmp)
    os.makedirs(os.path.join(HERE, 'results'), exist_ok=True)
    open(os.path.join(HERE, 'results', 'baseline_%s_TEST_REPORT.md' % tag), 'w').write(rep)
    rows = _report_rows(rep, r['tail'])
    fails = {}
    for sc, res, _ in rows:
        if res == 'FAIL':
            blk = re.search(r'%s output</summary>\n\n```\n(.*?)```' % re.escape(sc), rep, re.S)
            fails[sc] = blk.group(1).strip().splitlines()[-8:] if blk else r['tail'][-8:]
    if r['rc'] and not fails: fails['(run)'] = r['tail'][-8:]
    inner = re.findall(r'restored byte-identical \(md5 of every file \+ git status\) = (.*?); their fresh', rep)
    return {'suite': rel, 'cmd': ' '.join([script] + list(args)), 'rc': r['rc'], 's': r['s'], 'rows': rows, 'fails': fails, 'inner_restore_note': inner[0] if inner else None,
            'restored_clean': C.git('status', '--short', '--', rel) == git_before and _manifest(d) == before, 'files': len(before)}


STEP0_EDITS = [
    ('test/nemoback-extension/test_extension.py', 'S3 (no random.sample left): strip # comments (tokenize) before scanning; the agreed else comment mentions random.sample'),
    ('test/nemoback-extension/test_extension.py', 'marker_bookkeeping: world m1 shares the factories with the later worlds -> recycle m1 (teardown clear-ups; C.release: without teardown\'s minted == recycled assert, which cannot hold while m is alive) before check_world([m])'),
    ('test/nemoback-metacommunity/test_metacommunity.py', "R1: skip methods with fn.__module__ == 'extension.global_habitat_network' when collecting is_remove defaults"),
    ('test/nemoback-patch/test_chain.py', 'lines 279-280: removed == set() and changed == set() (metacommunity.py is migrated); the stale "next migration round" finding deleted')]


def baseline_main():
    out = {'edits': STEP0_EDITS, 'indfactory': []}
    d = os.path.join(REPO, 'test', 'nemo-indfactory')
    before = _manifest(d)
    for s in sorted(os.path.basename(x) for x in glob.glob(os.path.join(d, '*.py')) if not x.endswith('common.py')):
        r = ra_run(s, d, keep=4); out['indfactory'].append({'script': s, 'rc': r['rc'], 's': r['s'], 'tail': r['tail']})
        print('baseline nemo-indfactory/%s rc=%d %.1fs' % (s, r['rc'], r['s']), flush=True)
    out['indfactory_clean'] = _manifest(d) == before
    out['suites'] = []
    for rel, script, args, tag in (('test/nemoback-habitat', 'run_all.py', [], 'habitat'), ('test/nemoback-patch', 'run_all.py', ['--no-baseline'], 'patch'),
                                   ('test/nemoback-metacommunity', 'test_metacommunity.py', ['all'], 'metacommunity')):
        r = _run_suite_restored(rel, script, args, tag); out['suites'].append(r)
        print('baseline %s %s rc=%d %.1fs restored=%s %s' % (rel, r['cmd'], r['rc'], r['s'], r['restored_clean'], r['rows']), flush=True)
    for x in out['indfactory']: check('Step 0 test/nemo-indfactory/%s' % x['script'], x['rc'] == 0, x['tail'][-1:] if x['rc'] else '')
    check('Step 0 test/nemo-indfactory untouched', out['indfactory_clean'])
    for x in out['suites']:
        check('Step 0 %s %s (%s)' % (x['suite'], x['cmd'], ', '.join('%s %s' % (a, b) for a, b, _ in x['rows'])), x['rc'] == 0, x['fails'])
        check('Step 0 %s restored byte-identical (md5 of every file + git status)' % x['suite'], x['restored_clean'])
    mc = next(x for x in out['suites'] if x['suite'].endswith('metacommunity'))
    if mc.get('inner_restore_note') and 'False' in mc['inner_restore_note']:
        finding('info', 'test/nemoback-metacommunity (Step 0)', "its own nested Step 0 reports a restore flag False (%s): its check requires `git status --short -- <suite>` to be EMPTY, and test/nemoback-patch/test_chain.py now carries "
                "the Step 0 edit (modified). The md5 manifests match (this suite's own restore check of the same folder is True). Not a defect." % mc['inner_restore_note'],
                fix='compare git status before / after instead of requiring it to be empty (test-side only)')
    return out


# ############################################################################
# PART: report
# ############################################################################
PARTS_INFO = [('smoke', 'one step of each flow (marker True / False, pool) with invariants after every call'),
              ('static', 'A: diff whitelist, signatures, pool_name tuple, guard / else, TYPE_CHECKING, old name, callers, tags'),
              ('registry', 'B: order, global_index record, distance matrix (brute force + v3 + one recompute), edges, rate-matrix properties x 5 kernels, underflow'),
              ('dispersal', 'C: offspring network dispersal -- guard, patch_num < 2, %d-seed equivalence, conservation, v3 statistics, capping, distance law, Nemo level, chains' % N_EQUIV),
              ('marker', 'D: marker network dispersal -- both modes vs v3 (%d seeds each), germinate-and-birth homes, ValueError risk, default effect on GRFE births, chains' % N_EQUIV),
              ('sim', 'E: full simulations (2 scenarios x 3 groups + patch-level comparator), reproducibility, whole-run equivalence vs v3, ecology, disturbance'),
              ('perf', 'F: one call at H = 16 / 64 / 256 / 1024 vs v3, cost breakdown, zero-skip variant'),
              ('fuzz', 'G: random mixes of every meta operation + both network dispersals + guarded calls')]
FILES = ['extension/global_habitat_network.py', 'metaibm/metacommunity.py', 'metaibm/patch.py', 'metaibm/habitat.py', 'metaibm/indfactory.py', 'metaibm/individual.py', DOCS]


def esc(s): return str(s if s is not None else '').replace('|', '/').replace('\n', ' ')
def table(head, rows): return ['| ' + ' | '.join(head) + ' |', '|' + '---|' * len(head)] + ['| ' + ' | '.join(esc(c) for c in r) + ' |' for r in rows]
def load(name):
    try: return json.load(open(os.path.join(HERE, 'results', name + '.json')))
    except Exception: return {}
def res_word(ok): return 'PASS' if ok else '**FAIL**'


def checks_table(d):
    return table(['check', 'result', 'detail'], [(r['check'], res_word(r['ok']), r['detail'][:260]) for r in d.get('results', [])])


def coverage_rows():
    src, osrc = C.current_ext_source(), C.original_ext_source()
    cur, org = C._functions(src), C._functions(osrc)
    def strip_doc(fn):
        fn = ast.parse(ast.unparse(fn)).body[0]
        if fn.body and isinstance(fn.body[0], ast.Expr) and isinstance(getattr(fn.body[0], 'value', None), ast.Constant): fn.body = fn.body[1:] or [ast.Pass()]
        for a in fn.args.args + [fn.args.vararg, fn.args.kwarg]:
            if a is not None: a.annotation = None
        return ast.dump(fn)
    texts = {}
    me = open(os.path.abspath(__file__)).read()
    for part, body in re.findall(r'# PART: (\w+).*?\n(.*?)(?=\n# #+\n# PART: |\Z)', me, re.S): texts[part] = body
    texts['simulate_ext (E)'] = open(SIM_PY).read(); texts['common (schedules)'] = open(os.path.join(HERE, 'common.py')).read()
    rows = []
    for name, fn in cur.items():
        oname = C.ORIGINAL_REGISTRY if name == REG else name
        o = org.get(oname)
        if o is None: st = 'new'
        elif ast.get_source_segment(src, fn) == ast.get_source_segment(osrc, o): st = 'unchanged'
        elif name != oname and strip_doc(fn).replace(name, oname) == strip_doc(o): st = 'renamed only'
        elif strip_doc(fn) == strip_doc(o): st = 'docstring / comment / annotation only'
        else: st = 'changed' + (' (renamed)' if name != oname else '')
        alias = {OFF: ['OFF', 'DISP'], MARK: ['MARK'], REG: ['REG']}.get(name, [])          # the tests name these through constants
        rx = r'\b(%s)\b' % '|'.join([re.escape(name)] + alias)
        hits = [p for p, t in texts.items() if p != 'report' and re.search(rx, t)]
        tag = 'NemoBackend' if 'NemoBackend' in (ast.get_docstring(fn) or '') else ''
        rows.append((name, st, tag, ', '.join(hits) or '(via the dispersal calls)'))
    for oname in org:
        if oname not in cur and oname != C.ORIGINAL_REGISTRY: rows.append((oname, 'removed', '', '-'))
    return rows


def report_main():
    import importlib.metadata as md
    allr = load('all'); L = []
    now = datetime.datetime.now().strftime('%Y-%m-%d %H:%M')
    L += ['# extension/global_habitat_network.py test report (NemoBackend, sexual)', '',
          'Generated %s by `test/nemoback-extension/test_extension.py %s` (mode `%s`). Scope: `extension/global_habitat_network.py` (the registry and the two network dispersals migrated to NemoBackend) '
          'and `metacommunity.add_patch`, on the NemoBackend call chain metacommunity -> patch -> habitat -> indfactory -> individual; sexual only.' % (now, 'all --' + MODE if allr else 'report', MODE), '']
    # Environment
    L += ['## Environment', '']
    try: pv = md.version('pynemo')
    except Exception: pv = '?'
    L += ['- python `%s`, pynemo `%s`, numpy `%s`, `%s`' % (platform.python_version(), pv, np.__version__, platform.platform()),
          '- HEAD `%s` (`%s`), branch `%s`' % (C.git('rev-parse', '--short', 'HEAD').strip(), C.git('log', '-1', '--format=%s').strip(), C.git('rev-parse', '--abbrev-ref', 'HEAD').strip()), '']
    rows = []
    for f in FILES:
        h = hashlib.md5(open(os.path.join(REPO, f), 'rb').read()).hexdigest()
        rows.append(('`%s`' % f, '`%s`' % h, 'modified (working tree)' if C.git('status', '--porcelain', '--', f).strip() else 'committed'))
    L += table(['file', 'md5', 'state'], rows) + ['']
    L += ['- originals: extension = `git show %s:%s` (last commit that touched it), add_patch = `git show %s:%s` (HEAD). Their bodies are compiled into a subclass of the CURRENT metacommunity '
          '(they run on the NemoBackend patch / habitat). The record attribute is read from `_init_global_habitat_network_fields`: `%s` (both `global_habitat_id_2_index_dir` and `global_habitat_id_2_gidx_dir` accepted).'
          % (C.ORIGINAL_EXT_COMMIT, C.EXT_FILE, C.ORIGINAL_META_COMMIT, C.META_FILE, C.GIDX_DICT), '- mode `%s`: %s' % (MODE, Q('reduced simulation landscapes (net habitats 10x10, GRFE 3x3, 300 / 20 steps), fewer seeds and repetitions', 'experiment-size simulations (2000 steps, net habitats 40x40, GRFE 10x10), more seeds and repetitions')), '']
    # Step 0
    b = load('baseline')
    L += ['## Step 0 -- test-script fixes and baseline', '', 'Four test-script edits (no source code touched):', '']
    L += table(['file', 'edit'], [('`%s`' % f, e) for f, e in STEP0_EDITS]) + ['']
    if b:
        rows = [('test/nemo-indfactory', x['script'], res_word(x['rc'] == 0), '%.1f s' % x['s']) for x in b.get('indfactory', [])]
        for s in b.get('suites', []):
            for sc, res, t in s['rows']: rows.append(('%s (%s)' % (s['suite'], s['cmd']), sc, res_word(res == 'PASS'), (t + ' s') if t else ''))
            rows.append((s['suite'], '**%s total**' % s['cmd'], res_word(s['rc'] == 0), '%.1f s' % s['s']))
        L += table(['suite', 'script', 'result', 'time'], rows) + ['']
        L += ['- every folder was backed up before and restored after its run (their runners rewrite TEST_REPORT.md + results/): restored byte-identical (md5 of every file + git status vs before) = %s; '
              'test/nemo-indfactory untouched = %s. Fresh reports kept as `results/baseline_{habitat,patch,metacommunity}_TEST_REPORT.md`.'
              % (', '.join('%s `%s`' % (s['suite'].split('/')[-1], s['restored_clean']) for s in b.get('suites', [])), b.get('indfactory_clean')), '']
        for s in b.get('suites', []):
            for sc, tail in s['fails'].items(): L += ['`%s/%s` failed; last lines:' % (s['suite'], sc), '```'] + tail + ['```', '']
    # Smoke
    sm = load('smoke')
    L += ['## Smoke: %s (%s s)' % (res_word(sm.get('ok')), sm.get('seconds')), '']
    if sm.get('rows'): L += table(['flow', 'step', 'return', 'invariants'], sm['rows']) + ['']
    # Tests
    L += ['## Tests', '']
    rows = []
    for p, cov in PARTS_INFO:
        d = load(p); t = allr.get('parts', {}).get(p, {}).get('s', d.get('seconds'))
        rows.append((p, cov, res_word(d.get('ok')) if d else 'not run', '%d / %d' % (sum(r['ok'] for r in d.get('results', [])), len(d.get('results', []))) if d else '', '%s s' % t if t is not None else ''))
    L += table(['part', 'covers', 'result', 'checks', 'time'], rows) + ['']
    if allr: L += ['- `all --%s` wall clock %s s (parts run in parallel processes; perf ran alone)' % (MODE, allr.get('seconds')), '']
    for md_ in ('quick', 'full'):
        x = load('all_%s' % md_)
        if x and md_ != MODE:
            L += ['- last `all --%s` run (%s, %s s): ' % (md_, x.get('finished'), x.get('seconds')) + ', '.join('%s %s (%d/%d)' % (p, 'PASS' if v['rc'] == 0 else 'FAIL', *x['checks'].get(p, [0, 0])) for p, v in x['parts'].items()), '']
    # A
    d = load('static')
    L += ['## A static: %s' % res_word(d.get('ok')), ''] + checks_table(d) + ['']
    if d.get('diff'):
        lab = Counter(r[3] for r in d['diff'])
        L += ['Diff lines by label: ' + ', '.join('%s x%d' % kv for kv in lab.items()), '', '<details><summary>every +/- line</summary>', ''] + table(['file', '+/-', 'line', 'label'], d['diff']) + ['', '</details>', '']
    if d.get('callers'):
        L += ['Callers of the two network dispersals (repo-wide incl. untracked, excluding this suite):', ''] + table(['file:line', 'function', 'kind', 'passes is_remove', 'text'],
              [('%s:%d' % (c['file'], c['line']), c['function'], c['kind'], c['passes_is_remove'], c['text'][:110]) for c in d['callers']]) + ['']
    if d.get('old_name') is not None: L += ['Old name `%s`: %s' % (C.ORIGINAL_REGISTRY, '; '.join('%s (%s)' % (a[:90], k) for a, k in d['old_name']) or 'no hit'), '']
    if d.get('get_type_hints'): L += ['- `typing.get_type_hints(global_habitat_network_registry)` at runtime: %s (string annotations under TYPE_CHECKING are for static tools only)' % d['get_type_hints'], '']
    # B
    d = load('registry')
    L += ['## B registry and matrices: %s' % res_word(d.get('ok')), ''] + checks_table(d) + ['']
    uf = d.get('underflow', {})
    if uf.get('thresholds'):
        L += ['Kernel underflow -- smallest distance where `dispersal_kernel_strength` returns exactly 0.0 (bisection):', ''] + \
             table(['kernel', 'params', 'zero from distance'], [(r['kernel'], r['params'], 'never (< 1e300)' if r['zero_from_distance'] is None else '%.4g' % r['zero_from_distance']) for r in uf['thresholds']]) + \
             ['', '- experiments: %s' % uf.get('experiments'), '- pipeline (A at d=400 from B/C, exponential rho=2): %s' % uf.get('pipeline'), '- other patch without habitats: %s' % uf.get('empty_other_patch'), '']
    # C
    d = load('dispersal')
    L += ['## C offspring network dispersal: %s' % res_word(d.get('ok')), ''] + checks_table(d) + ['']
    if d.get('stats_vs_v3'):
        s = d['stats_vs_v3']
        L += ['Expected immigrants per destination habitat, current (move) vs v3 (copy), %d calls each, offspring pools %s, exponential rho=1:' % (s['R'], s['pools']), '']
        L += table(['m', 'mean current', 'mean v3', 'p (Welch) per destination', 'total current', 'total v3', 'pick uniform p'],
                   [(m, r['mean_current'], r['mean_v3'], r['p_current_vs_v3'], r['total_current'], r['total_v3'], r['pick_uniform_p']) for m, r in s['rates'].items()]) + ['']
    if d.get('truncation'):
        t = d['truncation']
        L += ['Capping (C6): one source habitat with `pool` offspring, 8 destinations at equal distance (uniform kernel), %d calls per cell. `P(excess)` = P(rounded request > pool); '
              '`shortfall by destination order` = mean requested - realized for j = 1..8 in loop order:' % t['R'], '']
        L += table(['pool', 'm', 'P(excess)', 'mean truncated', 'requested', 'realized', 'share truncated', 'shortfall by destination order'],
                   [(r['pool'], r['m'], r['P_excess'], r['mean_truncated'], r['requested_mean'], r['realized_mean'], r['truncated_share_of_request'], r['shortfall_by_destination_order']) for r in t['rows']]) + ['']
    if d.get('distance'):
        t = d['distance']
        L += ['Distance law (C7): source pool %d, m=%s, destinations at d = 1..6, %d calls per kernel:' % (t['pool'], t['m'], t['R']), '']
        L += table(['kernel', 'params', 'observed', 'expected', 'chi2 / df', 'p', 'log-slope'],
                   [(r['kernel'], r['params'], r['observed'], r['expected'], '%s / %s' % (r['chi2'], r['df']), r['p'], r.get('slope', '')) for r in t['rows']]) + ['']
    if d.get('capping_equiv'): L += ['- capping in the equivalence worlds: %s' % d['capping_equiv'], '']
    # D
    d = load('marker')
    L += ['## D marker network dispersal: %s' % res_word(d.get('ok')), ''] + checks_table(d) + ['']
    if d.get('default_effect'):
        de = d['default_effect']
        L += ['Default False -> True, GRFE step 0 from the same state (%s scale): births True %d vs False %d (ratio %s)' % (MODE, de['births_true'], de['births_false'], de['ratio_true_over_false']), '']
        L += table(list(de['rows'][0].keys()), [list(r.values()) for r in de['rows']]) + ['']
    # E
    d = load('sim')
    L += ['## E full simulations: %s' % res_word(d.get('ok')), '',
          'Scenarios copied from `experiments/model-sloss-global-habitat-network.py` (`net`: 16 patches x 4 habitats, rho 2, m_among 0.1) and `experiments/model-sloss-GRFE.py` (`grfe`: 16 patches x 64 habitats from the layout CSV, '
          'rho 0.2, m_among 0.01); sexual functions; 16 Nemo indfactories; distinct global_index per habitat (mainland 0.., metacommunity 100..). The metacommunity starts empty and is colonized from the mainland, as in the experiments. '
          '`net`\'s sexual branch calls the patch-level marker dispersal in the experiment (network line commented out): here the network line is used. Groups: mT = marker + network is_remove default (True), '
          'mF = marker + is_remove=False, pool = individual pools + network offspring dispersal, mPatch = marker + the patch-level dispersal (the comparator). Scale %s.' % MODE, ''] + checks_table(d) + ['']
    if d.get('runs'):
        rows = []
        for n, r in d['runs'].items():
            s = r.get('summary', {})
            ph = r.get('phase_seconds') or {}
            net_s = next((v for k, v in ph.items() if k in ('network_marker', 'network_offspring', 'patch_marker')), None)
            rows.append((n, r['status'], r.get('H'), r.get('sites'), r.get('steps_done'), r.get('total_s'), '%.3f' % s['step_s_mean'] if s.get('step_s_mean') else '',
                         net_s, s.get('rss_first_last'), s.get('alive_end'), s.get('occupancy_end'), s.get('richness_end'), r.get('teardown_ok')))
        L += table(['run', 'status', 'H', 'sites', 'steps', 'total s', 's/step', 'among-dispersal s (total)', 'RSS MB first/last', 'alive end', 'occupancy end', 'species end', 'teardown ok'], rows) + ['']
    if d.get('ecology'):
        L += ['Ecology (2nd-half means over steps, averaged over seeds):', ''] + table(['scenario/group', 'runs', 'alive', 'occupancy end', 'species', 'births / step', 'among-patch migrants / step', 'decay slope', 's/step'],
              [(k, v['runs'], v['alive_mean_2nd_half'], v['occupancy_end'], v['richness_mean_2nd_half'], v['births_mean_2nd_half'], v['among_patch_migrants_mean_2nd_half'],
                [None if s is None else round(s, 3) for s in v['decay_slope']], v['step_s_mean']) for k, v in d['ecology'].items()]) + ['']
        L += ['- network vs patch-level (relative differences): %s' % d.get('network_vs_patch'), '- is_remove True vs False (ratios mT / mF): %s' % d.get('true_vs_false'), '']
    if d.get('runs'):
        L += ['Distance decay (immigrants per habitat pair by distance bin):', '']
        for n, r in d['runs'].items():
            if r.get('decay', {}).get('points') and not n.startswith(('repro', 'equiv')): L += ['- `%s`: slope %s; %s' % (n, None if r['decay']['slope'] is None else round(r['decay']['slope'], 3), r['decay']['points'][:12])]
        L += ['']
    if d.get('disturbance'): L += ['Disturbance (patch_dist_rate 0.05): %s' % d['disturbance'], '']
    if d.get('runs'):
        L += ['<details><summary>phase seconds per run</summary>', '']
        for n, r in d['runs'].items(): L += ['- `%s`: %s' % (n, dict(list((r.get('phase_seconds') or {}).items())[:8]))]
        L += ['', '</details>', '']
    # F
    d = load('perf')
    L += ['## F performance: %s' % res_word(d.get('ok')), ''] + checks_table(d) + ['']
    if d.get('rows'):
        keys = list(d['rows'][0]['median_s'].keys())
        L += ['Median seconds per call (pools of 20 per habitat, m=0.1, exponential rho=2; perf ran alone):', ''] + table(['H (patches x habitats)'] + keys, [('%d (%dx%d)' % (r['H'], r['patches'], r['habitats_per_patch']),) + tuple(r['median_s'][k] for k in keys) for r in d['rows']]) + ['']
        L += ['Breakdown of one current marker call (share of the total):', ''] + table(['H', 'total s'] + list(d['rows'][0]['breakdown_share'].keys()) + ['sample calls', 'with n=0'],
              [(r['H'], r['total_s']) + tuple('%.0f%%' % (100 * v) for v in r['breakdown_share'].values()) + (r['sample_calls'], r['sample_calls_with_zero']) for r in d['rows']]) + ['']
        if d.get('pool_sweep'): L += ['Pool-size sweep (H=256, marker call, median s):', ''] + table(['pool per habitat'] + list(d['pool_sweep'][0]['median_s']) + ['current / v3'], [(x['pool'],) + tuple(x['median_s'].values()) + (x['ratio_current_over_v3'],) for x in d['pool_sweep']]) + ['']
        L += ['Build (all add_patch calls, distance matrix recomputed each time): %s' % d.get('build'), '']
    # G
    d = load('fuzz')
    L += ['## G fuzz: %s' % res_word(d.get('ok')), ''] + checks_table(d) + ['']
    if d.get('kinds'): L += ['- operations by kind: %s' % d['kinds'], '- allowed pre-existing exceptions: %s' % d.get('allowed'), '']
    # Findings
    F = []; seen = set()
    for p in ['baseline', 'smoke', 'static', 'registry', 'dispersal', 'marker', 'sim', 'perf', 'fuzz']:
        for f in load(p).get('findings', []):
            k = (f['severity'], f['func'], f['text'])
            if k not in seen: seen.add(k); F.append(f)
    order = {'bug': 0, 'risk': 1, 'info': 2, 'known': 3, 'todo': 4}
    F.sort(key=lambda f: order.get(f['severity'], 9))
    cnt = Counter(f['severity'] for f in F)
    L += ['## Findings (%d: bug %d, risk %d, info %d, known %d, todo(accepted) %d)' % (len(F), cnt['bug'], cnt['risk'], cnt['info'], cnt['known'], cnt['todo']), '',
          'Severity: **bug** = wrong behaviour of the code under test; **risk** = correct today, fragile / surprising / costly; **info** = observation or intended change vs v3; '
          '**known** = agreed / pre-existing; **todo(accepted)** = style items the author decided to keep. Every suggested fix is a suggestion only -- **not applied**.', '']
    L += table(['severity', 'function', 'finding', 'evidence', 'suggested fix'],
               [('todo(accepted)' if f['severity'] == 'todo' else f['severity'], f['func'], f['text'], json.dumps(f['evidence'], default=str)[:700] if f['evidence'] is not None else '',
                 ('%s (not applied)' % f['fix']) if f.get('fix') else '') for f in F]) + ['']
    # Coverage
    L += ['## Method coverage', '', '`vs v3` compares each function with `%s:%s` (AST without docstrings / annotations, and source text). `tests` = parts of this suite that name the function.' % (C.ORIGINAL_EXT_COMMIT, C.EXT_FILE), '']
    L += table(['function', 'vs v3', 'tag', 'tests'], coverage_rows()) + ['', '- `metacommunity.add_patch`: changed (one registry call; the distance matrix is recomputed inside the registry) -- static A1 / B0, every world build in B-G.', '']
    # Repo state
    L += ['## Repository state', '', '```', C.git('status', '--short').rstrip(), C.git('diff', '--stat').rstrip(), '```', '']
    open(os.path.join(HERE, 'TEST_REPORT.md'), 'w').write('\n'.join(L))
    print('TEST_REPORT.md written (%d lines)' % len(L))
    return {}


# ############################################################################
# all
# ############################################################################
def run_part(p):
    t = time.time()
    pr = subprocess.run([sys.executable, os.path.abspath(__file__), p], capture_output=True, text=True, cwd=HERE, env=PY_ENV)
    tail = (pr.stdout + pr.stderr).strip().splitlines()[-3:]
    print('[all] %s rc=%d %.1fs %s' % (p, pr.returncode, time.time() - t, tail[-1:] if tail else ''), flush=True)
    return {'rc': pr.returncode, 's': round(time.time() - t, 1), 'tail': tail}


def all_main():
    t0 = time.time(); parts = {}
    if '--no-baseline' not in sys.argv: parts['baseline'] = run_part('baseline')
    for p in ('smoke', 'static', 'perf'): parts[p] = run_part(p)
    from concurrent.futures import ThreadPoolExecutor
    with ThreadPoolExecutor(max_workers=5) as ex:
        fut = {p: ex.submit(run_part, p) for p in ('sim', 'registry', 'dispersal', 'marker', 'fuzz')}
        for p, f in fut.items(): parts[p] = f.result()
    summ = {'mode': MODE, 'parts': parts, 'seconds': round(time.time() - t0, 1), 'finished': datetime.datetime.now().strftime('%Y-%m-%d %H:%M'),
            'checks': {p: [sum(r['ok'] for r in load(p).get('results', [])), len(load(p).get('results', []))] for p in parts}}
    json.dump(summ, open(os.path.join(HERE, 'results', 'all.json'), 'w'), indent=1)
    json.dump(summ, open(os.path.join(HERE, 'results', 'all_%s.json' % MODE), 'w'), indent=1)          # kept per mode: the full report also shows the last quick run
    run_part('report')
    print('ALL %s in %.0f s: %s' % (MODE, time.time() - t0, {p: v['rc'] for p, v in parts.items()}))


# ############################################################################
# command line
# ############################################################################
def worker_cli():
    kind = sys.argv[1]; a = sys.argv[2:]
    if kind == 'equiv': r = equiv_worker(a[0], int(a[1]))
    elif kind == 'chain': r = chain_worker(int(a[0]), int(a[1]))
    elif kind == 'cstats': r = cstats_worker()
    elif kind == 'ctrunc': r = ctrunc_worker()
    elif kind == 'cdist': r = cdist_worker()
    elif kind == 'mequiv': r = marker_equiv_worker(a[0], a[1], int(a[2]))
    elif kind == 'mchain': r = marker_chain_worker(int(a[0]), int(a[1]), a[2])
    elif kind == 'fuzz': r = fuzz_worker(int(a[0]), int(a[1]), int(a[2]))
    print('JSON ' + json.dumps(r, default=str))
    sys.stdout.flush(); os._exit(0)                     # skip interpreter teardown: a broken run may hold individuals owned twice


MAINS = {'baseline': baseline_main, 'smoke': smoke_main, 'static': static_main, 'registry': registry_main, 'dispersal': dispersal_main, 'marker': marker_main,
         'perf': perf_main, 'fuzz': fuzz_main, 'sim': sim_main}

if __name__ == '__main__':
    cmds = list(MAINS) + ['report', 'all', 'worker']
    if len(sys.argv) < 2 or sys.argv[1] not in cmds: sys.exit('usage: python test_extension.py {%s} [--quick|--full]' % '|'.join(cmds))
    part = sys.argv.pop(1)
    if part == 'worker': worker_cli()
    elif part == 'all': all_main()
    elif part == 'report': report_main()
    else: part_cli(part, MAINS[part])
