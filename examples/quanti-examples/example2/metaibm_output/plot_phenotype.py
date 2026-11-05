''' plot_phenotype -- example2's islands over time: one panel per island habitat (the islands start empty, so a panel follows a place).
    example2 has one species (sp1): every panel shows all the individuals of its habitat.
    abundance.png:         every reproduce_mode overlaid; per mode, the occupied microsites of the habitat (line) against its capacity (dashed line).
    phenotype-1.png, phenotype-2.png: every reproduce_mode overlaid; per mode, the mean phenotype of the habitat's individuals (line)
                           and their 5-95% range (band), against the habitat's mean environment (dashed line): phenotype-1 <-> env1, phenotype-2 <-> env2.
    survival_rate.png:     each individual's expected survival under example2's selection, same layout as the phenotypes.
    Lives in metaibm_output, next to its data, one set per reproduce_mode (example2.main writes one per run):
    Reads  {sp,pheno1,pheno2}_dis_over_time-<mode>.csv.gz (example2.save_snapshot)
    Usage: python plot_phenotype.py [mode ...]   (no mode: every mode whose sp file is present)
'''
import inspect
import os
import re
import sys
import warnings

import numpy as np
import pandas as pd
import matplotlib
IN_KERNEL = 'ipykernel' in sys.modules              # run as a Jupyter / VS Code Interactive Window cell: no __file__, sys.argv holds the kernel's own flags
if not IN_KERNEL:
    matplotlib.use('Agg')                           # a script only writes pngs; inside a kernel, leave its backend alone
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

OUT_DIR = os.path.dirname(os.path.abspath(__file__)) if '__file__' in globals() else os.getcwd()   # metaibm_output: the csv.gz files are read from here, the pngs are written here
if not os.path.exists(os.path.join(OUT_DIR, 'plot_phenotype.py')):     # a kernel's cwd is the file's folder only by default (VS Code jupyter.notebookFileRoot)
    raise SystemExit('cannot locate metaibm_output from %s -- run this from metaibm_output (or %%cd there first)' % OUT_DIR)
sys.path.insert(0, os.path.dirname(OUT_DIR))                   # example2/: example2.py (file names, landscape) and bootstrap_metaibm.py
import example2
from example2 import build_metacommunity

MODES = ('asexual', 'sexual')                       # example2.main's plotted reproduce_mode values
TRAIT_2_ENV = {'phenotype-1': 'env1', 'phenotype-2': 'env2'}
BAND_Q = (5, 95)                                    # percentiles of the band

# colors: reference palette, light mode (dataviz skill references/palette.md)
SURFACE, INK, INK_2, MUTED, GRID, AXIS = '#fcfcfb', '#0b0b0b', '#52514e', '#898781', '#e1e0d9', '#c3c2b7'
MODE_COLOR = {'asexual': '#2a78d6', 'sexual': '#eb6834', 'mix': '#1baf7a'}   # categorical slots 1-3, fixed per mode (never by plotting order)


def dis_files(reproduce_mode):
    ''' the file names example2.main(reproduce_mode) writes; RETURN: sp_file, {trait_name: pheno_file} '''
    sp_file = 'sp_dis_over_time-%s.csv.gz' % reproduce_mode
    pheno_files = {'phenotype-1': 'pheno1_dis_over_time-%s.csv.gz' % reproduce_mode,
                   'phenotype-2': 'pheno2_dis_over_time-%s.csv.gz' % reproduce_mode}
    return sp_file, pheno_files


def read_dis(file_name):
    ''' one row per recorded step, one column per microsite (patch, habitat, microsite); RETURN: time steps, columns, values '''
    df = pd.read_csv(os.path.join(OUT_DIR, os.path.basename(file_name)), header=[0, 1, 2], index_col=0)
    time_steps = np.array([int(i.replace('time_step', '')) for i in df.index])
    return time_steps, df.columns, df.to_numpy(dtype=float)


def island_habitats():
    ''' {(patch_id, h_id): {env_name: mean}} from example2's islands (mean_env_ls as configured, not the drawn microsite values), in landscape order '''
    meta_obj = build_metacommunity()
    return {(patch_id, h_id): dict(zip(h_object.env_types_name, h_object.mean_env_ls))
            for patch_id, patch_obj in meta_obj.set.items() for h_id, h_object in patch_obj.set.items()}


def habitat_masks(sp_columns, habitats):
    ''' {(patch_id, h_id): bool mask over the microsite columns of that habitat} '''
    habitat_of_col = np.array(['%s|%s' % (p, h) for p, h, _ in sp_columns])
    return {key: habitat_of_col == '%s|%s' % key for key in habitats}


def example2_selection_params():
    ''' base_dead_rate, fitness_wid as example2.main sets them (locals of main, so read off its source); RETURN: (d, w) '''
    src = inspect.getsource(example2.main)
    found = re.search(r'base_dead_rate,\s*fitness_wid\s*=\s*([0-9.eE+-]+),\s*([0-9.eE+-]+)', src)
    if found is None:
        raise SystemExit("example2.main no longer sets 'base_dead_rate, fitness_wid = <d>, <w>' -- update example2_selection_params()")
    return float(found.group(1)), float(found.group(2))


def final_text(finals, fmt='%.3f'):
    ''' finals = [(reproduce_mode, last value of its series)] -> 'final  asexual 0.402 · sexual empty' (nan = the habitat is empty at the last step) '''
    return 'final  ' + ' · '.join('%s %s' % (m, 'empty' if np.isnan(v) else fmt % v) for m, v in finals)


def style_axes(ax):
    ''' recessive grid and axes '''
    ax.set_facecolor(SURFACE)
    ax.grid(axis='y', color=GRID, linewidth=0.8)
    ax.set_axisbelow(True)
    for side in ('top', 'right', 'left'):
        ax.spines[side].set_visible(False)
    ax.spines['bottom'].set_color(AXIS)
    ax.tick_params(colors=MUTED, labelcolor=INK_2, labelsize=9, length=0)


def style_panel(ax, key, finals, right_title='', fmt='%.3f'):
    ''' bold habitat title + right_title on the top line, the per-mode final values on a line of their own below them '''
    ax.set_title(key[0], loc='left', fontsize=11, color=INK, fontweight='bold', pad=16)
    ax.set_title(right_title, loc='right', fontsize=9, color=INK_2, pad=16)
    ax.text(0, 1.015, final_text(finals, fmt), transform=ax.transAxes, ha='left', va='bottom', fontsize=8.5, color=INK_2)
    style_axes(ax)


def finish_figure(fig, axes, runs, ylabel, title, subtitle, ref_label, out_png):
    ''' shared axes labels, title block, legend (one line per mode + the dashed reference), save '''
    for ax in axes[-1]:
        ax.set_xlabel('time step', color=INK_2, fontsize=10)
    for ax in axes[:, 0]:
        ax.set_ylabel(ylabel, color=INK_2, fontsize=10)
    all_steps = np.concatenate([run[0] for run in runs.values()])
    axes[0, 0].set_xlim(all_steps.min(), all_steps.max())

    fig.suptitle(title, x=0.06, y=0.975, ha='left', fontsize=14, color=INK, fontweight='bold')
    fig.text(0.06, 0.928, subtitle, ha='left', fontsize=10, color=INK_2)
    handles = [Line2D([], [], color=MODE_COLOR[m], linewidth=1.8, label=m) for m in runs] + \
              [Line2D([], [], color=INK_2, linewidth=1.2, linestyle=(0, (5, 3)), label=ref_label)]
    fig.legend(handles=handles, loc='upper right', bbox_to_anchor=(0.985, 0.985), ncol=len(handles), frameon=False, fontsize=9, labelcolor=INK_2)
    fig.tight_layout(rect=(0.03, 0, 1, 0.935), h_pad=2.5, w_pad=1.5)
    fig.savefig(out_png, dpi=150, facecolor=SURFACE)
    plt.close(fig)
    return out_png


def habitat_grid(habitats):
    ''' one panel per island habitat, two per row; RETURN: fig, axes (2-D) '''
    nrows = int(np.ceil(len(habitats) / 2))
    fig, axes = plt.subplots(nrows, 2, figsize=(10, 3.5 * nrows), sharex=True, sharey=True, facecolor=SURFACE, squeeze=False)
    for ax in axes.ravel()[len(habitats):]:
        ax.set_visible(False)
    return fig, axes


def draw_band(ax, time_steps, values, reproduce_mode):
    ''' one mode in one panel: mean line + 5-95% band of the individuals (gaps while the habitat is empty); RETURN: the mean series '''
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', RuntimeWarning)     # an empty habitat is an all-nan row: nan mean, nan band -- a gap, not an error
        mean = np.nanmean(values, axis=1)
        lo, hi = np.nanpercentile(values, BAND_Q, axis=1)
    ax.fill_between(time_steps, lo, hi, color=MODE_COLOR[reproduce_mode], alpha=0.15, linewidth=0)
    ax.plot(time_steps, mean, color=MODE_COLOR[reproduce_mode], linewidth=1.8)
    return mean


def plot_abundance(runs, habitats, out_png):
    ''' runs = {reproduce_mode: (time_steps, sp_values, {trait_name: pheno_values}, masks)} -- one line per mode in every panel: the habitat's occupied microsites '''
    fig, axes = habitat_grid(habitats)
    capacity = 0
    for ax, (key, env_means) in zip(axes.ravel(), habitats.items()):
        finals = []
        for reproduce_mode, (time_steps, sp_values, _, masks) in runs.items():
            own = sp_values[:, masks[key]]                                     # this habitat's microsites only; an empty microsite is nan
            counts = (~np.isnan(own)).sum(axis=1)
            ax.plot(time_steps, counts, color=MODE_COLOR[reproduce_mode], linewidth=1.8)
            finals.append((reproduce_mode, counts[-1]))
        capacity = max(capacity, own.shape[1])
        ax.axhline(own.shape[1], color=INK_2, linewidth=1.2, linestyle=(0, (5, 3)))
        style_panel(ax, key, finals, ' '.join('%s %.2f' % kv for kv in env_means.items()), fmt='%d')
    axes[0, 0].set_ylim(0, 1.05 * capacity)

    return finish_figure(fig, axes, runs, 'individuals', 'abundance on the islands',
                         'one panel per island habitat, one species (sp1); line = occupied microsites',
                         'capacity (%d microsites)' % capacity, out_png)


def plot_trait(trait_name, runs, habitats, out_png):
    ''' one line + band per mode in every panel: the habitat's individuals '''
    env_name = TRAIT_2_ENV[trait_name]
    fig, axes = habitat_grid(habitats)
    for ax, (key, env_means) in zip(axes.ravel(), habitats.items()):
        env = env_means[env_name]
        finals = []
        for reproduce_mode, (time_steps, sp_values, pheno_values, masks) in runs.items():
            mean = draw_band(ax, time_steps, pheno_values[trait_name][:, masks[key]], reproduce_mode)
            finals.append((reproduce_mode, mean[-1]))
        ax.axhline(env, color=INK_2, linewidth=1.2, linestyle=(0, (5, 3)))
        style_panel(ax, key, finals, '%s %.2f' % (env_name, env))

    return finish_figure(fig, axes, runs, trait_name, '%s vs %s on the islands' % (trait_name, env_name),
                         'one panel per island habitat, one species (sp1); line = mean phenotype, band = %d-%d%% of individuals' % BAND_Q,
                         '%s mean of the habitat' % env_name, out_png)


def plot_survival(runs, habitats, out_png):
    ''' expected survival rate of every recorded individual, from example2's niche_gaussian selection (habitat.survival_rate):
        s = (1 - d) * prod_k exp(-((z_k - e_k) / w)^2), with e_k = the habitat's mean env (the microsite draws are not recorded) '''
    d, w = example2_selection_params()
    fig, axes = habitat_grid(habitats)
    for ax, (key, env_means) in zip(axes.ravel(), habitats.items()):
        finals = []
        for reproduce_mode, (time_steps, sp_values, pheno_values, masks) in runs.items():
            mismatch = sum(((pheno_values[trait_name][:, masks[key]] - env_means[env_name]) / w) ** 2 for trait_name, env_name in TRAIT_2_ENV.items())
            mean = draw_band(ax, time_steps, (1 - d) * np.exp(-mismatch), reproduce_mode)
            finals.append((reproduce_mode, mean[-1]))
        ax.axhline(1 - d, color=INK_2, linewidth=1.2, linestyle=(0, (5, 3)))
        style_panel(ax, key, finals)
    axes[0, 0].set_ylim(0, 1)

    return finish_figure(fig, axes, runs, 'survival rate', 'expected survival rate on the islands',
                         'one panel per island habitat, one species (sp1); line = mean over individuals, band = %d-%d%%; selection d = %g, w = %g' % (BAND_Q + (d, w)),
                         '1 - d (phenotype on the optimum)', out_png)


def load_mode(reproduce_mode, habitats):
    ''' one example2 run; RETURN: time_steps, sp_values, {trait_name: pheno_values}, {(patch_id, h_id): column mask} '''
    sp_file, pheno_files = dis_files(reproduce_mode)
    time_steps, sp_columns, sp_values = read_dis(sp_file)
    pheno = {}
    for trait_name, pheno_file in pheno_files.items():
        pheno_steps, pheno_columns, pheno[trait_name] = read_dis(pheno_file)
        assert np.array_equal(pheno_steps, time_steps) and pheno_columns.equals(sp_columns), '%s and %s are not aligned' % (pheno_file, sp_file)
    return time_steps, sp_values, pheno, habitat_masks(sp_columns, habitats)


def main(modes):
    modes = modes or [m for m in MODES if os.path.exists(os.path.join(OUT_DIR, dis_files(m)[0]))]
    if not modes:
        raise SystemExit('no sp_dis_over_time-<mode>.csv.gz in %s -- run example2.py first' % OUT_DIR)
    habitats = island_habitats()
    runs = {m: load_mode(m, habitats) for m in modes}
    print('wrote', plot_abundance(runs, habitats, os.path.join(OUT_DIR, 'abundance.png')))
    for trait_name in TRAIT_2_ENV:
        print('wrote', plot_trait(trait_name, runs, habitats, os.path.join(OUT_DIR, '%s.png' % trait_name)))
    print('wrote', plot_survival(runs, habitats, os.path.join(OUT_DIR, 'survival_rate.png')))


if __name__ == '__main__':
    main([] if IN_KERNEL else sys.argv[1:])
