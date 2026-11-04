''' plot_phenotype -- example1's phenotype trajectories: one png per trait, one panel per species, every reproduce_mode overlaid.
    Each panel: per mode, the species' mean phenotype (line) and the 5-95% range of its individuals (band) over time,
    against the mean environment of the habitat the species lives in (dashed line): phenotype-1 <-> env1, phenotype-2 <-> env2.
    Lives in metaibm_output, next to its data, one set per reproduce_mode (example1.main writes one per run):
    Reads  {sp,pheno1,pheno2}_dis_over_time-<mode>.csv.gz (example1.save_snapshot)
    Writes phenotype-1.png, phenotype-2.png, and survival_rate.png (each individual's expected survival under example1's selection)
    Usage: python plot_phenotype.py [mode ...]   (no mode: every mode whose sp file is present)
'''
import inspect
import os
import re
import sys

import numpy as np
import pandas as pd
import matplotlib
IN_KERNEL = 'ipykernel' in sys.modules              # run as a Jupyter / VS Code Interactive Window cell: no __file__, sys.argv holds the kernel's own flags
if not IN_KERNEL:
    matplotlib.use('Agg')                           # a script only writes pngs; inside a kernel, leave its backend alone
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

OUT_DIR = os.path.dirname(os.path.abspath(__file__)) if '__file__' in globals() else os.getcwd()   # metaibm_output: the csv.gz files are read from here, the pngs are written here
if not os.path.exists(os.path.join(OUT_DIR, 'plot_phenotype.py')):     # a kernel's cwd is the file's folder only by default (VS Code jupyter.notebookFileRoot)
    raise SystemExit('cannot locate metaibm_output from %s -- run this from metaibm_output (or %%cd there first)' % OUT_DIR)
sys.path.insert(0, os.path.dirname(OUT_DIR))                   # example1/: example1.py (file names, landscape) and bootstrap_metaibm.py
import example1
from example1 import build_metacommunity

MODES = ('asexual', 'sexual')                       # example1.main's plotted reproduce_mode values
TRAIT_2_ENV = {'phenotype-1': 'env1', 'phenotype-2': 'env2'}
BAND_Q = (5, 95)                                    # percentiles of the band

# colors: reference palette, light mode (dataviz skill references/palette.md)
SURFACE, INK, INK_2, MUTED, GRID, AXIS = '#fcfcfb', '#0b0b0b', '#52514e', '#898781', '#e1e0d9', '#c3c2b7'
MODE_COLOR = {'asexual': '#2a78d6', 'sexual': '#eb6834', 'mix': '#1baf7a'}   # categorical slots 1-3, fixed per mode (never by plotting order)


def dis_files(reproduce_mode):
    ''' the file names example1.main(reproduce_mode) writes; RETURN: sp_file, {trait_name: pheno_file} '''
    sp_file = 'sp_dis_over_time-%s.csv.gz' % reproduce_mode
    pheno_files = {'phenotype-1': 'pheno1_dis_over_time-%s.csv.gz' % reproduce_mode,
                   'phenotype-2': 'pheno2_dis_over_time-%s.csv.gz' % reproduce_mode}
    return sp_file, pheno_files


def read_dis(file_name):
    ''' one row per recorded step, one column per microsite (patch, habitat, microsite); RETURN: time steps, columns, values '''
    df = pd.read_csv(os.path.join(OUT_DIR, os.path.basename(file_name)), header=[0, 1, 2], index_col=0)
    time_steps = np.array([int(i.replace('time_step', '')) for i in df.index])
    return time_steps, df.columns, df.to_numpy(dtype=float)


def habitat_env_means():
    ''' {(patch_id, h_id): {env_name: mean}} from example1's landscape (mean_env_ls as configured, not the drawn microsite values) '''
    meta_obj = build_metacommunity()
    return {(patch_id, h_id): dict(zip(h_object.env_types_name, h_object.mean_env_ls))
            for patch_id, patch_obj in meta_obj.set.items() for h_id, h_object in patch_obj.set.items()}


def species_homes(sp_columns, sp_values):
    ''' {species_id: (patch_id, h_id)} -- the habitat each species occupies most, over all recorded steps '''
    habitat_of_col = np.array(['%s|%s' % (p, h) for p, h, _ in sp_columns])
    homes = {}
    for sp_id in sorted(int(s) for s in np.unique(sp_values[~np.isnan(sp_values)])):
        occupied = np.where(sp_values == sp_id, habitat_of_col[None, :], '')
        labels, counts = np.unique(occupied[occupied != ''], return_counts=True)
        homes[sp_id] = tuple(labels[np.argmax(counts)].split('|'))
    return homes


def example1_selection_params():
    ''' base_dead_rate, fitness_wid as example1.main sets them (locals of main, so read off its source); RETURN: (d, w) '''
    src = inspect.getsource(example1.main)
    found = re.search(r'base_dead_rate,\s*fitness_wid\s*=\s*([0-9.eE+-]+),\s*([0-9.eE+-]+)', src)
    if found is None:
        raise SystemExit("example1.main no longer sets 'base_dead_rate, fitness_wid = <d>, <w>' -- update example1_selection_params()")
    return float(found.group(1)), float(found.group(2))


def all_homes(runs):
    ''' {species_id: (patch_id, h_id)} over every mode '''
    homes = {}
    for _, _, _, mode_homes in runs.values():
        homes.update({sp_id: home for sp_id, home in mode_homes.items() if sp_id not in homes})
    return dict(sorted(homes.items()))


def final_text(finals):
    ''' finals = [(reproduce_mode, last value of its mean series)] -> 'final  asexual 0.402 · sexual extinct' (nan = no individual left) '''
    return 'final  ' + ' · '.join('%s %s' % (m, 'extinct' if np.isnan(v) else '%.3f' % v) for m, v in finals)


def style_panel(ax, sp_id, home, finals, right_title=''):
    ''' bold species title + right_title on the top line, the per-mode final values on a line of their own below them '''
    ax.set_title('sp%d  ·  %s' % (sp_id, home[0]), loc='left', fontsize=11, color=INK, fontweight='bold', pad=16)
    ax.set_title(right_title, loc='right', fontsize=9, color=INK_2, pad=16)
    ax.text(0, 1.015, final_text(finals), transform=ax.transAxes, ha='left', va='bottom', fontsize=8.5, color=INK_2)
    ax.grid(axis='y', color=GRID, linewidth=0.8)
    ax.set_axisbelow(True)
    for side in ('top', 'right', 'left'):
        ax.spines[side].set_visible(False)
    ax.spines['bottom'].set_color(AXIS)
    ax.tick_params(colors=MUTED, labelcolor=INK_2, labelsize=9, length=0)


def finish_figure(fig, axes, runs, ylabel, title, subtitle, ref_label, out_png):
    ''' shared axes labels, title block, legend (one line per mode + the dashed reference), save '''
    for ax in axes[1]:
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


def draw_band(ax, time_steps, values, reproduce_mode):
    ''' one mode in one panel: mean line + 5-95% band of the individuals; RETURN: the mean series '''
    mean = np.nanmean(values, axis=1)
    lo, hi = np.nanpercentile(values, BAND_Q, axis=1)
    ax.fill_between(time_steps, lo, hi, color=MODE_COLOR[reproduce_mode], alpha=0.15, linewidth=0)
    ax.plot(time_steps, mean, color=MODE_COLOR[reproduce_mode], linewidth=1.8)
    return mean


def plot_trait(trait_name, runs, env_means, out_png):
    ''' runs = {reproduce_mode: (time_steps, sp_values, {trait_name: pheno_values}, homes)} -- one line + band per mode in every panel '''
    env_name = TRAIT_2_ENV[trait_name]
    fig, axes = plt.subplots(2, 2, figsize=(10, 7), sharex=True, sharey=True, facecolor=SURFACE)
    for ax, (sp_id, home) in zip(axes.ravel(), all_homes(runs).items()):
        env = env_means[home][env_name]
        ax.set_facecolor(SURFACE)
        finals = []
        for reproduce_mode, (time_steps, sp_values, pheno_values, _) in runs.items():
            own = np.where(sp_values == sp_id, pheno_values[trait_name], np.nan)   # this species' individuals only
            finals.append((reproduce_mode, draw_band(ax, time_steps, own, reproduce_mode)[-1]))
        ax.axhline(env, color=INK_2, linewidth=1.2, linestyle=(0, (5, 3)))
        style_panel(ax, sp_id, home, finals, '%s %.2f' % (env_name, env))

    return finish_figure(fig, axes, runs, trait_name, '%s vs %s' % (trait_name, env_name),
                         'one panel per species; line = mean phenotype, band = %d-%d%% of individuals' % BAND_Q,
                         '%s mean of its habitat' % env_name, out_png)


def plot_survival(runs, env_means, out_png):
    ''' expected survival rate of every recorded individual, from example1's niche_gaussian selection (habitat.survival_rate):
        s = (1 - d) * prod_k exp(-((z_k - e_k) / w)^2), with e_k = the habitat's mean env (the microsite draws are not recorded) '''
    d, w = example1_selection_params()
    fig, axes = plt.subplots(2, 2, figsize=(10, 7), sharex=True, sharey=True, facecolor=SURFACE)
    for ax, (sp_id, home) in zip(axes.ravel(), all_homes(runs).items()):
        ax.set_facecolor(SURFACE)
        finals = []
        for reproduce_mode, (time_steps, sp_values, pheno_values, _) in runs.items():
            mismatch = sum(((pheno_values[trait_name] - env_means[home][env_name]) / w) ** 2 for trait_name, env_name in TRAIT_2_ENV.items())
            survival = np.where(sp_values == sp_id, (1 - d) * np.exp(-mismatch), np.nan)   # this species' individuals only
            finals.append((reproduce_mode, draw_band(ax, time_steps, survival, reproduce_mode)[-1]))
        ax.axhline(1 - d, color=INK_2, linewidth=1.2, linestyle=(0, (5, 3)))
        style_panel(ax, sp_id, home, finals)
    axes[0, 0].set_ylim(0, 1)

    return finish_figure(fig, axes, runs, 'survival rate', 'expected survival rate',
                         'one panel per species; line = mean over individuals, band = %d-%d%%; selection d = %g, w = %g' % (BAND_Q + (d, w)),
                         '1 - d (phenotype on the optimum)', out_png)


def load_mode(reproduce_mode):
    ''' one example1 run; RETURN: time_steps, sp_values, {trait_name: pheno_values}, homes '''
    sp_file, pheno_files = dis_files(reproduce_mode)
    time_steps, sp_columns, sp_values = read_dis(sp_file)
    pheno = {}
    for trait_name, pheno_file in pheno_files.items():
        pheno_steps, pheno_columns, pheno[trait_name] = read_dis(pheno_file)
        assert np.array_equal(pheno_steps, time_steps) and pheno_columns.equals(sp_columns), '%s and %s are not aligned' % (pheno_file, sp_file)
    return time_steps, sp_values, pheno, species_homes(sp_columns, sp_values)


def main(modes):
    modes = modes or [m for m in MODES if os.path.exists(os.path.join(OUT_DIR, dis_files(m)[0]))]
    if not modes:
        raise SystemExit('no sp_dis_over_time-<mode>.csv.gz in %s -- run example1.py first' % OUT_DIR)
    runs, env_means = {m: load_mode(m) for m in modes}, habitat_env_means()
    for trait_name in TRAIT_2_ENV:
        print('wrote', plot_trait(trait_name, runs, env_means, os.path.join(OUT_DIR, '%s.png' % trait_name)))
    print('wrote', plot_survival(runs, env_means, os.path.join(OUT_DIR, 'survival_rate.png')))


if __name__ == '__main__':
    main([] if IN_KERNEL else sys.argv[1:])
