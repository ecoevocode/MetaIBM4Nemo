import os, sys, tempfile
REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))   # test/nemo-indfactory/ -> repo root
sys.path.insert(0, REPO)                      # metaibm + extension packages; pynemo comes from the conda env
from pynemo.params import SimParameters
from metaibm.indfactory import indfactory

NEMO_CACHE = tempfile.mkdtemp(prefix='nemo_cache_')

def nemo_params(sp_log_name, init_value, mutation_rate=0.0001, env_variance=0.002025, extra=()):
    params = SimParameters()
    rows = [('replicates', 1), ('generations', 1), ('patch_number', 1), ('patch_capacity', 1), ('aging', 1),
            ('run_mode', 'overwrite'), ('root_dir', NEMO_CACHE), ('logfile', NEMO_CACHE + '/nemo.log'), ('filename', sp_log_name),
            ('quanti_traits', 2), ('quanti_loci', 40), ('quanti_allele_model', 'continuous'),
            ('quanti_init_model', 2), ('quanti_init_value', '{{%g, %g}}' % (init_value, init_value)), ('quanti_init_variance', '{{0.002025, 0.002025}}'),
            ('quanti_mutation_rate', mutation_rate), ('quanti_mutation_variance', 0.01)]
    if env_variance is not None: rows.append(('quanti_environmental_variance', '{{%g, %g}}' % (env_variance, env_variance)))
    for k, v in list(rows) + list(extra):
        params.add(k, v)
    return params
