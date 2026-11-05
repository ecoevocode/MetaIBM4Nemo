# example1 - a bilingual tutorial 

Four species, four patches: MetaIBM's life cycle on Nemo's genetics (metaibm-nemo); each species lives alone in its own patch.
Each species = one `spN.json` (its values, then a `description` of the species). What every row means in general: `../USER-MANUAL-quanti-parameters.md`;
This file explains the four genetically varied species in a simulation and what their values do, with real individuals printed out.

1. Genetics basics
2. Views of the genetic model: the genetic map, genotype -> phenotype, parents -> offspring
3. Implementation: what MetaIBM and Nemo each do
4. The rows of the four species: shared and not shared
5. The four species: each one's rows, real individuals and mutations

Every number marked "measured" comes from a real run of example1, not from theory.

## 1. Genetics basics

Every individual is diploid and carries 2 quantitative traits: `phenotype-1` and `phenotype-2`.

| Term | Meaning in example1 |
|---|---|
| locus | one gene position; the genome has 40 (`quanti_loci`) |
| allele copy | each locus is carried twice: copy 0 comes from the father, copy 1 from the mother (a clone gets both from its mother) |
| allele value $a$ | the number one allele copy adds to its trait |
| allele model | which values an allele can take and what a mutation does to it: di-allelic (two values; a mutation flips one into the other) or continuous (any real number; a mutation adds to it, or replaces it) |
| locus value $g$ | the two copies of one locus together: $g = a_1 + a_2$; with dominance, plus $k\,\lvert a_2 - a_1 \rvert$ |
| genotypic value $G$ | one trait's sum of locus values over all its loci, plus the epistatic terms; inherited |
| phenotype $P$ | $P = G + \varepsilon$, $\varepsilon \sim \mathcal{N}(0, V_e)$, drawn once when the individual is made; what MetaIBM's selection reads; $\varepsilon$ is not inherited |
| $V_g$, $V_e$, $h^2$ | the variance of $G$ and of $\varepsilon$ among individuals; heritability $h^2 = V_g / (V_g + V_e)$ |
| pleiotropy | one locus acting on several traits: off in sp1-sp3 (each locus belongs to one trait), full in sp4 (every locus acts on both traits, with one allele value per trait) |
| dominance | at a heterozygous locus $g$ differs from $a_1 + a_2$ (sp2) |
| epistasis | the locus values of different loci multiply into interaction terms (sp3) |
| mutation | at each birth, each allele copy mutates with probability $\mu$ (`quanti_mutation_rate`) |
| recombination | which of a parent's two copies each locus passes on; rate 0.5 = an independent coin flip per locus (no linkage) |

## 2. Views of the genetic model

### 2.1 The genetic map: where each locus sits

Nemo stores the genome as two sequences (copy 0, copy 1), with one entry per (locus, trait). Each position can be looked up with `get_locus_seq_pos(locus, trait)`:

```
sp1, sp2, sp3 -- no pleiotropy: 40 loci, 20 per trait (quanti_loci_per_trait {{20, 20}})

  sequence position    0   1   2   3  ...  19 | 20  21  22  23  ...  39
  trait                1   1   1   1  ...   1 |  2   2   2   2  ...   2
  locus (in trait)     1   2   3   4  ...  20 |  1   2   3   4  ...  20
  copy 0 (father)      a   a   a   a  ...   a |  a   a   a   a  ...   a
  copy 1 (mother)      a   a   a   a  ...   a |  a   a   a   a  ...   a
                       '-----------'             '-----------'
                       sp3: loci 1-4 of each trait interact (epistasis)

sp4 -- full pleiotropy: 40 loci, each holds one allele value per trait

  sequence position    0   1 |  2   3 |  4   5 | ... | 78  79
  locus                1   1 |  2   2 |  3   3 | ... | 40  40
  trait                1   2 |  1   2 |  1   2 | ... |  1   2
```

There are no map distances. At `quanti_recombination_rate 0.5` every locus is inherited on its own, so the order above is only storage order.
Note: Nemo keeps one genetic map per Python process, shared by every species started in that process, so all four species have `quanti_loci` 40.

### 2.2 Genotype -> phenotype

```
                  copy 0 (father)              copy 1 (mother)
  locus l, trait d     a_l,0,d                      a_l,1,d
                         \                         /
  locus value  (1) g_l,d = a_l,0,d + a_l,1,d + k * |a_l,1,d - a_l,0,d|
                   k = quanti_dominance_effects (sp2: 0.5, others: 0)
                                    |
  trait d      (2) G_d = sum_l g_l,d + sum_((i,j) in P) e_ij * g_i,d * g_j,d
                   l: the loci of trait d (20; sp4: all 40); P, e_ij: sp3 only
                   sp4, full pleiotropy: each copy of a locus holds one value per trait,
                   so one locus gives trait 1 its g_l,1 and trait 2 its own g_l,2
                                    |
  phenotype    (3) P_d = G_d + envNoise,  envNoise ~ N(0, Ve = 0.025^2)
                   envNoise drawn once, at birth
                                    |
  MetaIBM          selection (niche_gaussian) reads (P_1, P_2)
                   against the habitat's (env1, env2)
```

$$
g_{l,d} = a_{l,0,d} + a_{l,1,d} + k\,\lvert a_{l,1,d} - a_{l,0,d} \rvert \tag{1}
$$

$$
G_d = \sum_{l \in L_d} g_{l,d} + \sum_{(i,j) \in \mathcal{P}} e_{ij}\, g_{i,d}\, g_{j,d} \tag{2}
$$

$$
P_d = G_d + \varepsilon, \qquad \varepsilon \sim \mathcal{N}(0, V_e) \tag{3}
$$

| Symbol | Meaning | In example1 |
|---|---|---|
| $a_{l,0,d}$, $a_{l,1,d}$ | allele values of copy 0 (father) and copy 1 (mother) at locus $l$ for trait $d$. Without pleiotropy each locus belongs to one trait and each copy holds one value. Under full pleiotropy each copy holds one value per trait (sequence position $2(l-1)$ for trait 1, $2(l-1)+1$ for trait 2); each trait uses its own value, computed separately by (1) and (2) | see section 5 |
| $k$ | dominance coefficient, `quanti_dominance_effects` | sp2: 0.5; others: 0, so $g_{l,d} = a_{l,0,d} + a_{l,1,d}$ |
| $L_d$ | the loci of trait $d$ | 20 per trait; sp4: all 40, the same loci in both $L_1$ and $L_2$, but each trait gets its own $g_{l,1}$ or $g_{l,2}$ |
| $\mathcal{P}$ | the interacting locus pairs, each counted once: every pair within `quanti_epistatic_network` | sp3: the 6 pairs of loci 1-4; others: none |
| $e_{ij}$ | epistatic coefficient of the pair $(i,j)$, `quanti_epistatic_effects` | sp3: $\pm 20$ |
| $\varepsilon$ | environmental deviation, not inherited | drawn once when the individual is made |
| $V_e$ | `quanti_environmental_variance` | $0.025^2$ |



### 2.3 Parents -> offspring

| Birth | MetaIBM call | The newborn's alleles |
|---|---|---|
| founder | `indfactory.make_founder(gender, home_habitat)` | drawn from the species' init rows (see section 5); sp1: then redrawn by MetaIBM from `quanti_init_freq` |
| clone (asexual) | `indfactory.make_clone(mother, home_habitat)` | a copy of the mother's two copies, then mutation |
| sexual | `indfactory.make_offspring(mother, father, gender, home_habitat)` | copy 1 from the mother, copy 0 from the father; which of the parent's two copies each locus passes on is a coin flip; then mutation |

Every birth ends with $P$ recomputed from the new alleles by (1)-(3), with a new $\varepsilon$.

What one mutation does:

| Allele model (species) | A mutation ... |
|---|---|
| `diallelic` (sp1) | flips the allele: $+a \leftrightarrow -a$ |
| `diallelic_asymmetrical` (sp2) | flips the allele: $0.0075 \leftrightarrow 0.0225$ |
| `continuous` (sp3) | **adds** a draw of $\mathcal{N}(0, \sigma^2)$ to the allele, $\sigma^2$ = `quanti_mutation_variance` |
| `continuous_in_place` + full pleiotropy (sp4) | **replaces** the allele's values on both traits with a pair drawn from $\mathcal{N}(0, \sigma^2)$, whose correlation is `quanti_mutation_correlation` |

## 3. Implementation

One species = one `metaibm.indfactory` = one Nemo simulation (used only to make individuals, never run) + one Nemo `IndFactory`.
Nemo keeps the genes and does every genetic operation in C++. MetaIBM keeps the landscape, decides who dies and who breeds, and reads back only the phenotypes.

| Step | MetaIBM (Python) | Nemo (C++, through pynemo) |
|---|---|---|
| boot | `indfactory(species_id, pheno_name_ls, nemo_params, trait_keys=('quant',), seed, metaibm_nemo_params)` | builds the species' genetic model from the `quanti_*` rows: map, allele model, mutation model, $V_e$ |
| founder | `make_founder(gender, home_habitat)` | makes the alleles by `quanti_init_model`, then computes $P = G + \varepsilon$ |
| founder, sp1 only | done inside `make_founder()`: every copy is set to bit 0 with probability `quanti_init_freq`, else to bit 1; then $P$ is recomputed | -- |
| clone | `make_clone(mother, home_habitat)` | copies the mother's alleles, mutates, computes $P$ |
| sexual | `make_offspring(mother, father, gender, home_habitat)` | inherits with recombination, mutates, computes $P$ |
| selection | `individual.get_indi_phenotype_ls()` | returns the stored $P$ |
| death | `individual.recycle()` | puts the individual back into Nemo's recycling pool, for reuse at a later birth |

## 4. The rows of the four species

### 4.1 Rows shared by the four species

| Row | Value | What it does here |
|---|---|---|
| `replicates`, `generations`, `patch_number`, `patch_capacity`, `aging` | 1 | Nemo needs them to set up its simulation; example1 never runs that simulation, so they have no effect on the individuals |
| `run_mode` | 'silent' | no Nemo console output |
| `root_dir`, `logfile` | 'nemo_cache', 'nemo_cache/nemo.log' | where Nemo puts its log files |
| `quanti_traits` | 2 | two traits, named by `pheno_name_ls=('phenotype-1', 'phenotype-2')` |
| `quanti_loci` | 40 | loci per genome |
| `quanti_environmental_variance` | {{$0.025^2$, $0.025^2$}} | $V_e$ per trait: $\varepsilon \sim \mathcal{N}(0, 0.025^2)$ |
| `quanti_recombination_rate` | 0.5 | free recombination |

Every row outside the two tables of this section is `"_UNSET"` in all four species. `"_UNSET"` is example1's placeholder: the row stays in the json file but is never passed to Nemo, so Nemo uses its own default.

### 4.2 Rows not shared by the four species

"--" means the species leaves the row unset (`"_UNSET"`). What each row does in each species: see the matching subsection of section 5.

| Row | sp1 | sp2 | sp3 | sp4 | What it does |
|---|---|---|---|---|---|
| `filename` | 'sp1' | 'sp2' | 'sp3' | 'sp4' | names each species' Nemo log; must differ between species |
| `quanti_loci_per_trait` | {{20, 20}} | {{20, 20}} | {{20, 20}} | -- | 20 loci of its own per trait; cannot be set under full pleiotropy |
| `quanti_pleiotropy` | -- | -- | -- | 'full' | every locus acts on both traits |
| `quanti_allele_model` | 'diallelic' | 'diallelic_asymmetrical' | 'continuous' | 'continuous_in_place' | which values an allele can take and what a mutation does to it (see section 2.3) |
| `quanti_allele_value` | 0.008 | {{0.0075}{0.0225}} | -- | -- | allele values of the di-allelic models: $\pm a$ for sp1, {{allele 0}{allele 1}} for sp2 |
| `quanti_init_model` | -- | 0 | 2 | 1 | how Nemo makes the founder alleles |
| `quanti_init_value` | -- | -- | {{0.5, 0.5}} | {{0.7, 0.7}} | founder $E[G]$ per trait |
| `quanti_init_variance` | -- | -- | {{$0.025^2$, $0.025^2$}} | -- | founder additive $V_g$ per trait |
| `metaibm_nemo_params` `quanti_init_freq` | 0.65625 | -- | -- | -- | founder alleles redrawn by MetaIBM: bit 0 ($+a$) with probability $p$ |
| `quanti_mutation_rate` | 1e-4 | 1e-3 | 1e-4 | 1e-4 | $\mu$, the mutation probability per allele copy per birth |
| `quanti_mutation_variance` | -- | -- | $0.025^2$ | $0.01^2$ | $\sigma^2$ of the mutation effects; only the continuous models (sp3, sp4) use it, a di-allelic mutation is just a flip |
| `quanti_mutation_correlation` | -- | -- | -- | 0.5 | correlation between the two traits' effects of one mutation |
| `quanti_dominance_effects` | -- | 0.5 | -- | -- | dominance coefficient $k$ in (1) |
| `quanti_epistatic_network` | -- | -- | {{1, 2, 3, 4}} | -- | the interacting loci in (2): loci 1-4 interact pairwise |
| `quanti_epistatic_effects` | -- | -- | {{20, -20, 20, -20, 20, -20}} | -- | $e_{ij}$ in (2), one per locus pair |

The four species share one design:
- At most one of pleiotropy / dominance / epistasis per species, so each one can be watched on its own.
- Founder $V_g$ is kept small next to the founders' distance to the optimum, so at first, selection on this variation is equally fast in the sexual and the asexual mode.
- After that, asexual stalls near its best founder clone and waits for mutations, while sexual recombination keeps making new genotypes. That is what makes sexual faster.
- `quanti_init_model` sets how Nemo makes the founder alleles (0 when unset). Under the four `quanti_allele_model`s it means the following, with $n$ the loci per trait (20 without pleiotropy, 40 under sp4's full pleiotropy). Measured values in brackets: 900 founders, seed 1, no pleiotropy, $n = 20$, no dominance, no epistasis; $a = 0.008$, {{0.0075}{0.0225}}, $G_0 = 0.5$, $V_0 = 0.025^2$, $V_m = 0.01^2$.

  | `quanti_allele_model` | bit 0 / bit 1 value | 0 | 1 | 2 |
  |---|---|---|---|---|
  | 'diallelic' | bit 0 = $+a$, bit 1 = $-a$ (`quanti_allele_value` = $a$) | every copy is bit 0 ($+a$): $G = 2na$ (0.32), $V_g = 0$ | every copy draws bit 0 or bit 1 on its own ($E[G]$ = -0.009, $V_g$ = 0.0030) | loci alternate, homozygous: odd-numbered loci bit 0 on both copies ($+a$), even-numbered loci bit 1 on both ($-a$): $G = 0$, $V_g = 0$ |
  | 'diallelic_asymmetrical' | bit 0 = $a_0$, bit 1 = $a_1$ (`quanti_allele_value` = {{$a_0$}{$a_1$}}) | every copy is bit 0 ($a_0$): $G = 2n\,a_0$ (0.3), $V_g = 0$ | every copy draws $a_0$ or $a_1$ on its own ($E[G]$ = 0.61, $V_g$ = 0.0027) | odd-numbered loci $a_0$ on both copies, even-numbered loci $a_1$ on both: $G = n\,(a_0 + a_1)$ (0.6), $V_g = 0$ |
  | 'continuous' | no bits: the allele is any real value (its bit always reads 0) | every copy $= G_0 / 2n$: $G = G_0$, $V_g = 0$ | every copy $= G_0 / 2n$, then at each locus one random copy gets a mutation effect $\mathcal{N}(0, V_m)$ **added**: $V_g = n\,V_m$ (0.0021 / 0.0020) | every copy $= G_0 / 2n + \mathcal{N}(0,\, V_0 / 2n)$: $V_g = V_0$ (0.00071 / 0.00059); `quanti_init_variance` must be set |
  | 'continuous_in_place' | as 'continuous' | as 'continuous' | as 'continuous' (added too, not replaced) | as 'continuous' |

## 5. The four species

Founders, theoretical values (the same for both traits):

| Species | Habitat optimum | mean $P$ | mean $G$ | $V_g$ | $\mathrm{sd}(P - G)$ | $h^2 = V_g / (V_g + V_e)$ |
|---|---|---|---|---|---|---|
| sp1 | 0.3 | $E[G] = 0.1$ | $40\,a\,(2p - 1) = 0.1$ | $40\,(2a)^2\,p\,(1 - p) = 0.0023$ | $\sqrt{V_e} = 0.025$ | 0.79 |
| sp2 | 0.5 | $E[G] = 0.3$ | $40\,a_0 = 0.3$ | 0 | $\sqrt{V_e} = 0.025$ | 0 |
| sp3 | 0.7 | $E[G] = 0.5$ | $G_0 = 0.5$ | $\frac{V_0}{20}\left[16 + \sum_{i=1}^{4}(1 + \bar g\, c_i)^2\right] + \sum e_{ij}^2 \left(\frac{V_0}{20}\right)^2 = 0.00072$ | $\sqrt{V_e} = 0.025$ | 0.54 |
| sp4 | 0.9 | $E[G] = 0.7$ | $G_0 = 0.7$ | $40\,V_m = 0.004$ | $\sqrt{V_e} = 0.025$ | 0.86 |

sp3: $V_0 / 20$ is the variance of one locus's $g$, $\bar g = 0.025$, and $c_i = \sum_j e_{ij}$ = 20, 20, -60, 20, so the bracket is $16 + 3 \times 1.5^2 + 0.5^2 = 23$; $\sum e_{ij}^2 = 6 \times 20^2$.

### 5.1 sp1 -- diallelic, founders from quanti_init_freq (the baseline species: no pleiotropy, dominance or epistasis)

| Row | Value | Meaning |
|---|---|---|
| `quanti_allele_model` | 'diallelic' | two alleles per locus, symmetric around 0: bit 0 = $+a$, bit 1 = $-a$ |
| `quanti_allele_value` | 0.008 | $a$ |
| `quanti_init_model` | unset | Nemo's own founder draw, overwritten by quanti_init_freq below |
| `metaibm_nemo_params` `quanti_init_freq` | 0.65625 | re-implemented by MetaIBM: every founder allele copy is redrawn on its own, bit 0 ($+a$) with probability $p$, else bit 1 ($-a$) |

Expected, where $40 = 20 \times 2$ is the number of allele copies per trait (20 loci, 2 copies each):

$$
E[G] = 40\,a\,(2p - 1) = 0.1, \qquad V_g = 40\,(2a)^2\,p\,(1 - p) = 0.0023, \qquad G_{\max} = 40\,a = 0.32
$$

$G_{\max}$ lies past the optimum 0.3. Measured: $E[G] = 0.100$, $V_g$ = 0.0024 / 0.0023 (seed 1, 900 founders).

The founder at patch1 microsite (0, 0):

```
>>> print(ind)
individual
  metaibm-info : species_id=sp1, age=0, pheno_names=('phenotype-1', 'phenotype-2')
  nemo-info    : id=1, gender=male, phenotype=(0.0800, 0.0662), founder
>>> home habitat 0
trait 1 'phenotype-1': 20 loci at sequence positions 0-19
  bit 0 = allele value +0.0080, bit 1 = allele value -0.0080
  locus         1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16  17  18  19  20
  copy0 bit     1   1   1   0   1   0   1   1   0   0   0   1   1   0   0   0   0   1   0   0
  copy1 bit     0   0   0   0   0   0   1   1   0   1   0   0   0   1   0   1   1   0   0   0
  g (1e-3)      0   0   0  16   0  16 -16 -16  16   0  16   0   0   0  16   0   0   0  16  16
  G = 0.08000   P = 0.08001   P - G = +0.00001
trait 2 'phenotype-2': 20 loci at sequence positions 20-39
  bit 0 = allele value +0.0080, bit 1 = allele value -0.0080
  locus         1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16  17  18  19  20
  copy0 bit     0   0   0   0   1   0   1   0   0   0   0   0   1   0   0   1   1   0   1   1
  copy1 bit     1   0   1   0   1   0   1   0   0   0   0   1   1   0   1   0   0   0   1   1
  g (1e-3)      0  16   0  16 -16  16 -16  16  16  16  16   0 -16  16   0   0   0  16 -16 -16
  G = 0.06400   P = 0.06618   P - G = +0.00218
```

Reading it:
- $g = +2a = +0.016$ when both copies are bit 0, $g = 0$ when the locus is heterozygous, and $g = -2a = -0.016$ when both copies are bit 1.
- Trait 1: 7 loci at +16 and 2 at -16, so $G_1 = (7 - 2) \times 0.016 = 0.080$.
- Trait 2: 9 loci at +16 and 5 at -16, so $G_2 = (9 - 5) \times 0.016 = 0.064$.



### 5.2 sp2 -- diallelic_asymmetrical, founders from quanti_init_model 0, dominance

| Row | Value | Meaning |
|---|---|---|
| `quanti_allele_model` | 'diallelic_asymmetrical' | two arbitrary allele values per locus |
| `quanti_allele_value` | {{0.0075}{0.0225}} | two rows {{allele 0}{allele 1}}: bit 0 = 0.0075 (ancestral), bit 1 = 0.0225. One row {{b}} would mean $0$ / $b$, and a plain number the symmetric $\pm b$ |
| `quanti_init_model` | 0 | every founder allele copy is bit 0 |
| `quanti_dominance_effects` | 0.5 | $k$ in (1): $g = a_1 + a_2 + k\,\lvert a_2 - a_1 \rvert$ |
| `quanti_mutation_rate` | 1e-3 | 10x the others: with no standing variation, the default 1e-4 often lets the sexual run go extinct before mutations can rescue it |

Locus values by (1), by the number of mutant copies:
- 0 copies: $g = 2 \times 0.0075 = 0.015$.
- 1 copy (heterozygous): $g = 0.0075 + 0.0225 + 0.5 \times (0.0225 - 0.0075) = 0.0375$.
- 2 copies: $g = 2 \times 0.0225 = 0.045$.

So every founder has $G = 20 \times 0.015 = 0.3$ and $V_g = 0$: the founders are genetically identical. All variation comes from mutation, and every mutation raises $G$ ($0 \to 1$ copy: $+0.0225$; $1 \to 2$ copies: $+0.0075$).

The founder at patch2 microsite (0, 0):

```
>>> print(ind)
individual
  metaibm-info : species_id=sp2, age=0, pheno_names=('phenotype-1', 'phenotype-2')
  nemo-info    : id=901, gender=female, phenotype=(0.3294, 0.3413), founder
>>> home habitat 1
trait 1 'phenotype-1': 20 loci at sequence positions 0-19
  bit 0 = allele value +0.0075
  locus          1    2    3    4    5    6    7    8    9   10   11   12   13   14   15   16   17   18   19   20
  copy0 bit      0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0
  copy1 bit      0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0
  g (1e-3)    15.0 15.0 15.0 15.0 15.0 15.0 15.0 15.0 15.0 15.0 15.0 15.0 15.0 15.0 15.0 15.0 15.0 15.0 15.0 15.0
  G = 0.30000   P = 0.32939   P - G = +0.02939
trait 2 'phenotype-2': 20 loci at sequence positions 20-39
  bit 0 = allele value +0.0075
  locus          1    2    3    4    5    6    7    8    9   10   11   12   13   14   15   16   17   18   19   20
  copy0 bit      0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0
  copy1 bit      0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0
  g (1e-3)    15.0 15.0 15.0 15.0 15.0 15.0 15.0 15.0 15.0 15.0 15.0 15.0 15.0 15.0 15.0 15.0 15.0 15.0 15.0 15.0
  G = 0.30000   P = 0.34134   P - G = +0.04134
```

No dominance shows here: every locus is homozygous for bit 0, so $G = 0.3$, and founders differ only in $\varepsilon$.

A real mutation, from clones of this founder:

```
trait1 locus 17 copy0: 0.00750 -> 0.02250 (bit 0 -> 1)
```

This makes locus 17 heterozygous: $g$ goes from 0.015 to 0.0375.

Dominance shows only after mutations have built up. The individual below is from patch2 after example1 has run 200 steps in sexual mode (example1's default settings: base_dead_rate 0.1, fitness_wid 0.5, sexual_birth_rate 1). Patch2 then holds 900 individuals with 16.5 mutant allele copies each on average (max 30). This one sits at microsite (0, 0) and shows all three genotype classes:

```
>>> print(ind)
individual
  metaibm-info : species_id=sp2, age=0, pheno_names=('phenotype-1', 'phenotype-2')
  nemo-info    : id=347924, gender=male, phenotype=(0.4418, 0.5533), mother 336859, father 344773
>>> home habitat 1
trait 1 'phenotype-1': 20 loci at sequence positions 0-19
  bit 0 = allele value +0.0075, bit 1 = allele value +0.0225
  locus          1    2    3    4    5    6    7    8    9   10   11   12   13   14   15   16   17   18   19   20
  copy0 bit      0    1    0    0    0    0    0    0    0    1    0    0    0    0    0    0    0    0    1    0
  copy1 bit      0    1    0    0    1    1    0    1    0    0    0    0    0    0    0    1    0    0    0    0
  g (1e-3)    15.0 45.0 15.0 15.0 37.5 37.5 15.0 37.5 15.0 37.5 15.0 15.0 15.0 15.0 15.0 37.5 15.0 15.0 37.5 15.0
  G = 0.46500   P = 0.44185   P - G = -0.02315
trait 2 'phenotype-2': 20 loci at sequence positions 20-39
  bit 0 = allele value +0.0075, bit 1 = allele value +0.0225
  locus          1    2    3    4    5    6    7    8    9   10   11   12   13   14   15   16   17   18   19   20
  copy0 bit      0    1    0    0    0    1    0    0    0    0    0    1    0    0    1    1    0    0    1    0
  copy1 bit      0    0    0    1    1    0    0    0    0    1    1    1    0    0    0    1    0    0    1    0
  g (1e-3)    15.0 37.5 15.0 37.5 37.5 37.5 15.0 15.0 15.0 37.5 37.5 45.0 15.0 15.0 37.5 45.0 15.0 15.0 45.0 15.0
  G = 0.54750   P = 0.55335   P - G = +0.00585
```

Reading it:
- In units of $10^{-3}$: $g = 15.0$ at loci with no mutant copy, $g = 37.5$ at heterozygous loci, and $g = 45.0$ at loci with both copies bit 1 (loci 12, 16 and 19 of trait 2).
- Trait 2 has 7 heterozygous loci, 10 ancestral and 3 with both copies mutant: $G_2 = 7 \times 0.0375 + 10 \times 0.015 + 3 \times 0.045 = 0.5475$.
- Without dominance ($k = 0$), each heterozygous locus would give only $0.0075 + 0.0225 = 0.030$, and $G_2$ would be $7 \times 0.030 + 10 \times 0.015 + 3 \times 0.045 = 0.495$.

### 5.3 sp3 -- continuous, founders from quanti_init_model 2, epistasis

| Row | Value | Meaning |
|---|---|---|
| `quanti_allele_model` | 'continuous' | any real value; a mutation **adds** a draw of $\mathcal{N}(0, \sigma^2)$ to the allele |
| `quanti_init_model` | 2 | every founder allele copy $= G_0 / 40 + \mathcal{N}(0,\, V_0 / 40)$, with $G_0$, $V_0$ the next two rows |
| `quanti_init_value` | {{0.5, 0.5}} | $G_0$: founder $E[G]$ per trait |
| `quanti_init_variance` | {{$0.025^2$, $0.025^2$}} | $V_0$: founder additive $V_g$ per trait |
| `quanti_epistatic_network` | {{1, 2, 3, 4}} | the set $\mathcal{P}$ of (2): loci 1-4 interact pairwise. They are counted within a trait (1-20), and the same pairs act on both traits |
| `quanti_epistatic_effects` | {{20, -20, 20, -20, 20, -20}} | $e_{ij}$ of (2), one per pair, in the order (1,2) (1,3) (1,4) (2,3) (2,4) (3,4) |
| `quanti_mutation_variance` | $0.025^2$ | $\sigma^2$: a mutation adds a draw of $\mathcal{N}(0, 0.025^2)$ to the allele |

Expected and measured:
- An allele copy is about $0.5 / 40 = 0.0125$, so $g \approx 0.025$.
- One pair therefore adds $e_{ij}\, g_i\, g_j \approx \pm 20 \times 0.025^2 = \pm 0.0125$.
- Three positive and three negative pairs cancel at the mean, so $E[G]$ stays 0.5. The founders only gain non-additive variance: measured $V_g$ = 0.00074 / 0.00081, against $0.025^2$ additive.

The founder at patch3 microsite (0, 0):

```
>>> print(ind)
individual
  metaibm-info : species_id=sp3, age=0, pheno_names=('phenotype-1', 'phenotype-2')
  nemo-info    : id=1801, gender=male, phenotype=(0.5726, 0.5373), founder
>>> home habitat 2
trait 1 'phenotype-1': 20 loci at sequence positions 0-19
  locus     1     2     3     4     5     6     7     8     9    10    11    12    13    14    15    16    17    18    19    20
  copy0 11.04 17.67 11.45 12.42 15.71 13.39 11.61 26.60  7.20 20.25 13.28 14.57 12.67 17.85 14.69 13.62 13.33 11.06 16.35  9.86
  copy1  5.03 13.16 12.43 15.16 15.79  9.75 12.26 13.41 14.15 14.98 18.29 15.03 14.59  3.47 16.41 10.84 13.21 20.09 11.61 14.35
  g     16.07 30.83 23.87 27.58 31.50 23.14 23.88 40.01 21.35 35.23 31.57 29.60 27.26 21.32 31.10 24.47 26.54 31.15 27.95 24.21
  (allele values and g in units of 1e-3)
  epistasis e_ij*g_i*g_j:
    (1,2) +20*0.0161*0.0308 = +0.00991,  (1,3) -20*0.0161*0.0239 = -0.00767,  (1,4) +20*0.0161*0.0276 = +0.00886
    (2,3) -20*0.0308*0.0239 = -0.01472,  (2,4) +20*0.0308*0.0276 = +0.01700,  (3,4) -20*0.0239*0.0276 = -0.01317
  sum g = 0.54863, sum epistasis = +0.00021
  G = 0.54883   P = 0.57259   P - G = +0.02376
trait 2 'phenotype-2': 20 loci at sequence positions 20-39
  locus     1     2     3     4     5     6     7     8     9    10    11    12    13    14    15    16    17    18    19    20
  copy0 14.25 13.27 11.48 16.11 16.23  9.06 12.20 13.58 13.18 11.13 11.12 16.11 11.01 14.16 14.48  9.15 16.48 16.65 15.09  9.51
  copy1  7.51  9.36 18.87 10.84 10.87 14.90 10.90  6.09 15.73  8.06 17.11  5.55 12.96  3.26 12.59 16.22  7.01 19.12 13.65 11.87
  g     21.76 22.63 30.36 26.95 27.10 23.96 23.10 19.67 28.92 19.19 28.22 21.66 23.96 17.42 27.06 25.37 23.49 35.77 28.74 21.38
  (allele values and g in units of 1e-3)
  epistasis e_ij*g_i*g_j:
    (1,2) +20*0.0218*0.0226 = +0.00985,  (1,3) -20*0.0218*0.0304 = -0.01321,  (1,4) +20*0.0218*0.0270 = +0.01173
    (2,3) -20*0.0226*0.0304 = -0.01374,  (2,4) +20*0.0226*0.0270 = +0.01220,  (3,4) -20*0.0304*0.0270 = -0.01636
  sum g = 0.49672, sum epistasis = -0.00953
  G = 0.48719   P = 0.53728   P - G = +0.05009
```

Reading it: the two sums of (2) are `sum g` and `sum epistasis` in the printout.
- Trait 1: $G_1 = 0.54863 + 0.00021 = 0.54883$. The six pair terms almost cancel.
- Trait 2: $G_2 = 0.49672 - 0.00953 = 0.48719$.

A real mutation, from clones of this founder:

```
trait2 locus 9 copy0: 0.01318 -> 0.01066 (bit 0 -> 0)
```

The allele value moved by -0.00253, one draw of $\mathcal{N}(0, 0.025^2)$ added to it. A continuous model has no allele index, so the bit is always 0.

### 5.4 sp4 -- continuous_in_place, founders from quanti_init_model 1, full pleiotropy

| Row | Value | Meaning |
|---|---|---|
| `quanti_allele_model` | 'continuous_in_place' | a mutation **replaces** the allele with its effect, drawn from $\mathcal{N}(0, \sigma^2)$ |
| `quanti_loci_per_trait` | unset | not with full pleiotropy: all 40 loci act on both traits |
| `quanti_pleiotropy` | 'full' | each locus holds one allele value per trait |
| `quanti_mutation_correlation` | 0.5 | correlation between the two traits' effects of one mutation |
| `quanti_mutation_variance` | $0.01^2$ | $\sigma^2$ per trait, which also sets founder $V_g = 40\,\sigma^2 \approx 0.0632^2$ |
| `quanti_init_model` | 1 | every founder allele value $= G_0 / 80$ ($G_0$ the next row), then each locus gets one mutation effect **added** (added, even in this in-place model) |
| `quanti_init_value` | {{0.7, 0.7}} | $G_0$: founder $E[G]$ per trait |

- Equivalent matrix form: drop the `quanti_mutation_correlation` and `quanti_mutation_variance` rows and set `quanti_mutation_matrix` {{$\sigma_1^2$, $\rho\,\sigma_1\sigma_2$}{$\rho\,\sigma_1\sigma_2$, $\sigma_2^2$}} = {{0.0001, 0.00005}{0.00005, 0.0001}} (here $\sigma_1 = \sigma_2 = 0.01$, $\rho = 0.5$). The diagonal holds each trait's variance; the off-diagonal holds the covariance $\rho\,\sigma_1\sigma_2 = 0.5 \times 0.01 \times 0.01$, not the correlation. Measured: both forms give exactly the same results.

How a founder is made:
- Each locus draws one correlated effect pair (trait 1, trait 2).
- Each trait's part lands on a randomly chosen copy, independently of the other trait.
- $G_1$ and $G_2$ are therefore correlated (measured $r = 0.52$). But on one copy, the two traits' values need not come from the same draw.
- Mutations during the run are different: one mutation replaces both traits' values on the same copy (see below).

The founder at patch4 microsite (0, 0):

```
>>> print(ind)
individual
  metaibm-info : species_id=sp4, age=0, pheno_names=('phenotype-1', 'phenotype-2')
  nemo-info    : id=2701, gender=female, phenotype=(0.6291, 0.5764), founder
>>> home habitat 3
trait 1 'phenotype-1': 40 loci at sequence positions 0, 2, ..., 78 (step 2)
  locus     1     2     3     4     5     6     7     8     9    10    11    12    13    14    15    16    17    18    19    20    21    22    23    24    25    26     27    28    29    30    31     32    33    34    35    36    37    38    39    40
  copy0  3.73  8.75  0.96  8.75  8.45  8.75 -1.05  8.75  8.75 18.82  8.75  9.81  2.72  8.75  6.46  8.75 19.09  8.75 24.84  8.75  9.87  8.75  4.80  8.75  8.75 15.88   8.75 12.56 -4.04  8.75  3.08   8.75 -7.59  8.75 16.29 19.39  8.75  3.79 29.08  8.75
  copy1  8.75 10.76  8.75 19.99  8.75 -7.03  8.75  0.03 -1.59  8.75 -1.92  8.75  8.75  1.89  8.75 18.38  8.75 -2.58  8.75  2.63  8.75 19.79  8.75 -9.39  8.19  8.75  -9.76  8.75  8.75 -0.37  8.75   7.64  8.75 28.68  8.75  8.75 15.91  8.75  8.75  5.18
  g     12.48 19.51  9.71 28.74 17.20  1.72  7.70  8.78  7.16 27.57  6.83 18.56 11.47 10.64 15.21 27.13 27.84  6.17 33.59 11.38 18.62 28.54 13.55 -0.64 16.94 24.63  -1.01 21.31  4.71  8.38 11.83  16.39  1.16 37.43 25.04 28.14 24.66 12.54 37.83 13.93
  (allele values and g in units of 1e-3)
  G = 0.65337   P = 0.62907   P - G = -0.02430
trait 2 'phenotype-2': 40 loci at sequence positions 1, 3, ..., 79 (step 2)
  locus     1     2     3     4     5     6     7     8     9    10    11    12    13    14    15    16    17    18    19    20    21    22    23    24    25    26     27    28    29    30    31     32    33    34    35    36    37    38    39    40
  copy0 10.69 -3.07  7.91 14.98 16.58  8.75 -1.86  1.12  1.26 20.54 -8.57 10.12  1.24  8.75 -1.56  7.39 13.18 -5.21 16.49  6.59  3.13  8.75  5.25  1.05 -3.05  3.69 -15.65 18.38 15.58  8.75  0.21 -11.98  8.75  8.75 11.57  2.80  8.75  8.75 23.64  8.75
  copy1  8.75  8.75  8.75  8.75  8.75 14.46  8.75  8.75  8.75  8.75  8.75  8.75  8.75 13.96  8.75  8.75  8.75  8.75  8.75  8.75  8.75 14.44  8.75  8.75  8.75  8.75   8.75  8.75  8.75  6.48  8.75   8.75  6.65 15.67  8.75  8.75  7.52 -5.18  8.75 18.59
  g     19.44  5.68 16.66 23.73 25.33 23.21  6.89  9.87 10.01 29.29  0.18 18.87  9.99 22.71  7.19 16.14 21.93  3.54 25.24 15.34 11.88 23.19 14.00  9.80  5.70 12.44  -6.90 27.13 24.33 15.23  8.96  -3.23 15.40 24.42 20.32 11.55 16.27  3.57 32.39 27.34
  (allele values and g in units of 1e-3)
  G = 0.60505   P = 0.57639   P - G = -0.02866
```

Reading it:
- 8.75 ($0.7 / 80 = 8.75 \times 10^{-3}$) is an allele copy that received no init effect.
- At every locus, each trait has exactly one copy that differs from 8.75.
- At locus 2, the trait-1 part sits on copy 1 (10.76) and the trait-2 part on copy 0 (-3.07): different copies.
- Adding up the 40 loci by (2): $G_1 = 0.65337$, $G_2 = 0.60505$.

A real mutation, from clones of this founder:

```
trait1 locus 28 copy0: 0.01256 -> 0.00427 (bit 0 -> 0)
trait2 locus 28 copy0: 0.01838 -> -0.00091 (bit 0 -> 0)
```

One mutation hits one locus on one copy, and replaces both traits' values at once.
The new values are centred on 0, while founder allele copies sit near 0.00875.
So in sp4 a mutation lowers that copy's value by about 0.00875 on average. Selection toward the optimum 0.9 has to use the rarer mutations that go up.

---

# example1

四个物种，四个 patch：在 Nemo 的遗传模型上运行 MetaIBM 的生命周期（metaibm-nemo）；每个物种单独生活在自己的 patch 里。
每个物种 = 一个 `spN.json`（先是参数值，最后是该物种的 `description`）。每一行参数的通用含义见 `../USER-MANUAL-quanti-parameters.md`（中文在后半部分）；
本文件讲解这四个遗传模型以及它们的参数值起什么作用，并打印真实个体作为例子。英文版在本文件前半部分。

1. 遗传学基础
2. 遗传模型的视图：基因图、基因型 -> 表型、亲代 -> 子代
3. 实现方式：MetaIBM 和 Nemo 各负责什么
4. 四个物种的参数行：共用的、不共用的
5. 四个物种：各自的参数行、真实的个体和突变

凡是标"实测"的数字都来自 example1 的真实运行，不是理论推算。

## 1. 遗传学基础

每个个体都是二倍体，带有 2 个数量性状：`phenotype-1` 和 `phenotype-2`。

| 术语 | 在 example1 中的含义 |
|---|---|
| 位点（locus） | 一个基因位置；基因组共 40 个（`quanti_loci`） |
| 等位基因拷贝（allele copy） | 每个位点有两份：拷贝 0 来自父亲，拷贝 1 来自母亲（克隆的两份都来自母亲） |
| 等位基因值 $a$ | 一份等位基因拷贝给它所属性状贡献的数值 |
| 等位基因模型 | 等位基因能取哪些值、突变对它做什么：双等位（两个值；突变把一个翻成另一个）或连续（任意实数；突变加到原值上，或替换原值） |
| 位点值 $g$ | 一个位点两份拷贝合起来的值：$g = a_1 + a_2$；有显性时再加 $k\,\lvert a_2 - a_1 \rvert$ |
| 基因型值 $G$ | 一个性状在其所有位点上的位点值之和，再加上上位性项；可遗传 |
| 表型 $P$ | $P = G + \varepsilon$，$\varepsilon \sim \mathcal{N}(0, V_e)$，在个体生成时抽一次；MetaIBM 的选择读的就是它；$\varepsilon$ 不遗传 |
| $V_g$、$V_e$、$h^2$ | 个体间 $G$ 的方差和 $\varepsilon$ 的方差；遗传力 $h^2 = V_g / (V_g + V_e)$ |
| 多效性（pleiotropy） | 一个位点同时作用于多个性状：sp1-sp3 关闭（每个位点只属于一个性状），sp4 为 full（每个位点作用于两个性状，每个性状各有一个等位基因值） |
| 显性（dominance） | 在杂合位点上 $g$ 不等于 $a_1 + a_2$（sp2） |
| 上位性（epistasis） | 不同位点的位点值相乘产生交互项（sp3） |
| 突变 | 每次出生时，每份等位基因拷贝以概率 $\mu$ 突变（`quanti_mutation_rate`） |
| 重组 | 每个位点传下去的是亲本两份拷贝中的哪一份；重组率 0.5 = 每个位点各自独立地掷一次硬币（无连锁） |

## 2. 遗传模型的视图

### 2.1 基因图：每个位点在哪里

Nemo 把基因组存成两条序列（拷贝 0、拷贝 1），每个（位点，性状）占一个位置。每个位置可以用 `get_locus_seq_pos(locus, trait)` 查到：

```
sp1, sp2, sp3 -- no pleiotropy: 40 loci, 20 per trait (quanti_loci_per_trait {{20, 20}})

  sequence position    0   1   2   3  ...  19 | 20  21  22  23  ...  39
  trait                1   1   1   1  ...   1 |  2   2   2   2  ...   2
  locus (in trait)     1   2   3   4  ...  20 |  1   2   3   4  ...  20
  copy 0 (father)      a   a   a   a  ...   a |  a   a   a   a  ...   a
  copy 1 (mother)      a   a   a   a  ...   a |  a   a   a   a  ...   a
                       '-----------'             '-----------'
                       sp3: loci 1-4 of each trait interact (epistasis)

sp4 -- full pleiotropy: 40 loci, each holds one allele value per trait

  sequence position    0   1 |  2   3 |  4   5 | ... | 78  79
  locus                1   1 |  2   2 |  3   3 | ... | 40  40
  trait                1   2 |  1   2 |  1   2 | ... |  1   2
```

这里没有图距。在 `quanti_recombination_rate 0.5` 下，每个位点各自独立遗传，所以上图的顺序只是存储顺序。
注意：Nemo 的基因图在一个 Python 进程里只有一份，同一进程里启动的所有物种共用它，所以四个物种的 `quanti_loci` 都是 40。

### 2.2 基因型 -> 表型

```
                  copy 0 (father)              copy 1 (mother)
  locus l, trait d     a_l,0,d                      a_l,1,d
                         \                         /
  locus value  (1) g_l,d = a_l,0,d + a_l,1,d + k * |a_l,1,d - a_l,0,d|
                   k = quanti_dominance_effects (sp2: 0.5, others: 0)
                                    |
  trait d      (2) G_d = sum_l g_l,d + sum_((i,j) in P) e_ij * g_i,d * g_j,d
                   l: the loci of trait d (20; sp4: all 40); P, e_ij: sp3 only
                   sp4, full pleiotropy: each copy of a locus holds one value per trait,
                   so one locus gives trait 1 its g_l,1 and trait 2 its own g_l,2
                                    |
  phenotype    (3) P_d = G_d + envNoise,  envNoise ~ N(0, Ve = 0.025^2)
                   envNoise drawn once, at birth
                                    |
  MetaIBM          selection (niche_gaussian) reads (P_1, P_2)
                   against the habitat's (env1, env2)
```

$$
g_{l,d} = a_{l,0,d} + a_{l,1,d} + k\,\lvert a_{l,1,d} - a_{l,0,d} \rvert \tag{1}
$$

$$
G_d = \sum_{l \in L_d} g_{l,d} + \sum_{(i,j) \in \mathcal{P}} e_{ij}\, g_{i,d}\, g_{j,d} \tag{2}
$$

$$
P_d = G_d + \varepsilon, \qquad \varepsilon \sim \mathcal{N}(0, V_e) \tag{3}
$$

| 符号 | 含义 | 在 example1 中 |
|---|---|---|
| $a_{l,0,d}$、$a_{l,1,d}$ | 位点 $l$ 上拷贝 0（父亲）和拷贝 1（母亲）在性状 $d$ 上的等位基因值。无多效性：每个位点只属于一个性状，每份拷贝只有一个值。full 多效性：每份拷贝在两个性状上各存一个值（序列位置 $2(l-1)$ 是性状 1，$2(l-1)+1$ 是性状 2），两个性状各用自己的那个值，按 (1)(2) 分别计算 | 见第 5 节 |
| $k$ | 显性系数，`quanti_dominance_effects` | sp2：0.5；其余：0，此时 $g_{l,d} = a_{l,0,d} + a_{l,1,d}$ |
| $L_d$ | 性状 $d$ 的位点集合 | 每个性状 20 个；sp4：全部 40 个，同一批位点同时在 $L_1$ 和 $L_2$ 里，但给两个性状的是各自的 $g_{l,1}$、$g_{l,2}$ |
| $\mathcal{P}$ | 相互作用的位点对，每对只算一次：`quanti_epistatic_network` 内的所有位点对 | sp3：位点 1-4 的 6 对；其余：没有 |
| $e_{ij}$ | 位点对 $(i,j)$ 的上位性系数，`quanti_epistatic_effects` | sp3：$\pm 20$ |
| $\varepsilon$ | 环境偏差，不遗传 | 个体生成时抽一次 |
| $V_e$ | `quanti_environmental_variance` | $0.025^2$ |



### 2.3 亲代 -> 子代

| 出生方式 | MetaIBM 调用 | 新个体的等位基因 |
|---|---|---|
| 建群个体 | `indfactory.make_founder(gender, home_habitat)` | 按该物种的初始化参数行抽取（见第 5 节）；sp1 之后还会由 MetaIBM 按 `quanti_init_freq` 重抽 |
| 克隆（无性） | `indfactory.make_clone(mother, home_habitat)` | 复制母亲的两份拷贝，然后突变 |
| 有性 | `indfactory.make_offspring(mother, father, gender, home_habitat)` | 拷贝 1 来自母亲，拷贝 0 来自父亲；每个位点传下亲本两份拷贝中的哪一份，掷硬币决定；然后突变 |

每次出生最后都会用新的等位基因按 (1)-(3) 重新计算 $P$，并重新抽一个 $\varepsilon$。

一次突变做了什么：

| 等位基因模型（物种） | 一次突变…… |
|---|---|
| `diallelic`（sp1） | 翻转等位基因：$+a \leftrightarrow -a$ |
| `diallelic_asymmetrical`（sp2） | 翻转等位基因：$0.0075 \leftrightarrow 0.0225$ |
| `continuous`（sp3） | 在原值上**加**一个 $\mathcal{N}(0, \sigma^2)$ 的抽样，$\sigma^2$ = `quanti_mutation_variance` |
| `continuous_in_place` + full 多效性（sp4） | 用一对从 $\mathcal{N}(0, \sigma^2)$ 抽出的值**替换**该等位基因在两个性状上的值，两者的相关系数为 `quanti_mutation_correlation` |

## 3. 实现方式

一个物种 = 一个 `metaibm.indfactory` = 一个 Nemo 模拟（只用来生成个体，从不运行）+ 一个 Nemo `IndFactory`。
Nemo 保存基因，并用 C++ 完成所有遗传操作。MetaIBM 负责景观、决定谁死谁繁殖，只把表型读回来。

| 步骤 | MetaIBM（Python） | Nemo（C++，经 pynemo） |
|---|---|---|
| 启动 | `indfactory(species_id, pheno_name_ls, nemo_params, trait_keys=('quant',), seed, metaibm_nemo_params)` | 由 `quanti_*` 行建立该物种的遗传模型：基因图、等位基因模型、突变模型、$V_e$ |
| 建群个体 | `make_founder(gender, home_habitat)` | 按 `quanti_init_model` 生成等位基因，然后计算 $P = G + \varepsilon$ |
| 建群个体，仅 sp1 | `make_founder()` 内部完成：每份拷贝以概率 `quanti_init_freq` 设为 bit 0，否则设为 bit 1；然后重新计算 $P$ | -- |
| 克隆 | `make_clone(mother, home_habitat)` | 复制母亲的等位基因，突变，计算 $P$ |
| 有性 | `make_offspring(mother, father, gender, home_habitat)` | 带重组的遗传，突变，计算 $P$ |
| 选择 | `individual.get_indi_phenotype_ls()` | 返回已存的 $P$ |
| 死亡 | `individual.recycle()` | 把个体放回 Nemo 的回收池，供之后的出生重用 |

## 4. 四个物种的参数行

### 4.1 四个物种共用的参数行

| 参数行 | 值 | 在这里的作用 |
|---|---|---|
| `replicates`、`generations`、`patch_number`、`patch_capacity`、`aging` | 1 | Nemo 建立模拟时需要这些行；example1 从不运行这个模拟，所以它们对个体没有影响 |
| `run_mode` | 'silent' | Nemo 不在控制台输出 |
| `root_dir`、`logfile` | 'nemo_cache'、'nemo_cache/nemo.log' | Nemo 日志文件的存放位置 |
| `quanti_traits` | 2 | 两个性状，名字由 `pheno_name_ls=('phenotype-1', 'phenotype-2')` 给出 |
| `quanti_loci` | 40 | 每个基因组的位点数 |
| `quanti_environmental_variance` | {{$0.025^2$, $0.025^2$}} | 每个性状的 $V_e$：$\varepsilon \sim \mathcal{N}(0, 0.025^2)$ |
| `quanti_recombination_rate` | 0.5 | 自由重组 |

本节两张表之外的参数行，四个物种都是 `"_UNSET"`。`"_UNSET"` 是 example1 的占位值：这一行留在 json 文件里，但不会传给 Nemo，Nemo 使用自己的默认值。

### 4.2 四个物种不共用的参数行

"--" 表示该物种未设置（`"_UNSET"`）。每一行在各物种里的具体作用见第 5 节对应的小节。

| 参数行 | sp1 | sp2 | sp3 | sp4 | 作用 |
|---|---|---|---|---|---|
| `filename` | 'sp1' | 'sp2' | 'sp3' | 'sp4' | 给每个物种的 Nemo 日志命名；各物种必须不同 |
| `quanti_loci_per_trait` | {{20, 20}} | {{20, 20}} | {{20, 20}} | -- | 每个性状各有 20 个自己的位点；full 多效性下不能设 |
| `quanti_pleiotropy` | -- | -- | -- | 'full' | 每个位点同时作用于两个性状 |
| `quanti_allele_model` | 'diallelic' | 'diallelic_asymmetrical' | 'continuous' | 'continuous_in_place' | 等位基因能取哪些值、突变对它做什么（见第 2.3 节） |
| `quanti_allele_value` | 0.008 | {{0.0075}{0.0225}} | -- | -- | 双等位模型的等位基因值：sp1 为 $\pm a$，sp2 为 {{allele 0}{allele 1}} |
| `quanti_init_model` | -- | 0 | 2 | 1 | Nemo 生成建群等位基因的方式 |
| `quanti_init_value` | -- | -- | {{0.5, 0.5}} | {{0.7, 0.7}} | 建群个体每个性状的 $E[G]$ |
| `quanti_init_variance` | -- | -- | {{$0.025^2$, $0.025^2$}} | -- | 建群个体每个性状的加性 $V_g$ |
| `metaibm_nemo_params` `quanti_init_freq` | 0.65625 | -- | -- | -- | 由 MetaIBM 重抽建群等位基因：以概率 $p$ 取 bit 0（$+a$） |
| `quanti_mutation_rate` | 1e-4 | 1e-3 | 1e-4 | 1e-4 | $\mu$，每份等位基因拷贝每次出生的突变概率 |
| `quanti_mutation_variance` | -- | -- | $0.025^2$ | $0.01^2$ | 突变效应的 $\sigma^2$；只有连续模型（sp3、sp4）用它，双等位的突变只是翻转 |
| `quanti_mutation_correlation` | -- | -- | -- | 0.5 | 同一次突变在两个性状上的效应之间的相关系数 |
| `quanti_dominance_effects` | -- | 0.5 | -- | -- | (1) 式中的显性系数 $k$ |
| `quanti_epistatic_network` | -- | -- | {{1, 2, 3, 4}} | -- | (2) 式中相互作用的位点：位点 1-4 两两交互 |
| `quanti_epistatic_effects` | -- | -- | {{20, -20, 20, -20, 20, -20}} | -- | (2) 式中的 $e_{ij}$，每个位点对一个 |

四个物种共用同一套设计：
- 每个物种最多只开启多效性 / 显性 / 上位性中的一个，以便单独观察每一种效应。
- 建群个体的 $V_g$ 相对于它们到最优值的距离较小，所以起初利用这部分变异的选择，在有性和无性两种模式下一样快。
- 之后无性繁殖停滞在最好的建群克隆附近、等待突变，而有性重组不断产生新的基因型。这就是有性繁殖更快的原因。
- `quanti_init_model` 决定 Nemo 怎么生成建群等位基因（未设置时为 0）。它在四种 `quanti_allele_model` 下的含义如下，$n$ 是每个性状的位点数（无多效性为 20，sp4 的 full 多效性为 40）。括号里是实测值：900 个建群个体，seed 1，无多效性，$n = 20$，无显性、无上位性；$a = 0.008$，{{0.0075}{0.0225}}，$G_0 = 0.5$，$V_0 = 0.025^2$，$V_m = 0.01^2$。

  | `quanti_allele_model` | bit 0 / bit 1 的值 | 0 | 1 | 2 |
  |---|---|---|---|---|
  | 'diallelic' | bit 0 = $+a$，bit 1 = $-a$（`quanti_allele_value` = $a$） | 每份拷贝都是 bit 0（$+a$）：$G = 2na$（0.32），$V_g = 0$ | 每份拷贝各自随机取 bit 0 或 bit 1（$E[G]$ = -0.009，$V_g$ = 0.0030） | 位点交替纯合：奇数号位点两份都是 bit 0（$+a$），偶数号位点两份都是 bit 1（$-a$）：$G = 0$，$V_g = 0$ |
  | 'diallelic_asymmetrical' | bit 0 = $a_0$，bit 1 = $a_1$（`quanti_allele_value` = {{$a_0$}{$a_1$}}） | 每份拷贝都是 bit 0（$a_0$）：$G = 2n\,a_0$（0.3），$V_g = 0$ | 每份拷贝各自随机取 $a_0$ 或 $a_1$（$E[G]$ = 0.61，$V_g$ = 0.0027） | 奇数号位点两份都是 $a_0$，偶数号位点两份都是 $a_1$：$G = n\,(a_0 + a_1)$（0.6），$V_g = 0$ |
  | 'continuous' | 没有 bit：等位基因是任意实数（读出来的 bit 总是 0） | 每份拷贝 $= G_0 / 2n$：$G = G_0$，$V_g = 0$ | 每份拷贝 $= G_0 / 2n$，再在每个位点随机挑一份拷贝**加**一个突变效应 $\mathcal{N}(0, V_m)$：$V_g = n\,V_m$（0.0021 / 0.0020） | 每份拷贝 $= G_0 / 2n + \mathcal{N}(0,\, V_0 / 2n)$：$V_g = V_0$（0.00071 / 0.00059）；必须设 `quanti_init_variance` |
  | 'continuous_in_place' | 同 'continuous' | 同 'continuous' | 同 'continuous'（也是加，不是替换） | 同 'continuous' |

## 5. 四个物种

建群个体的理论值（两个性状相同）：

| 物种 | 栖息地最优值 | mean $P$ | mean $G$ | $V_g$ | $\mathrm{sd}(P - G)$ | $h^2 = V_g / (V_g + V_e)$ |
|---|---|---|---|---|---|---|
| sp1 | 0.3 | $E[G] = 0.1$ | $40\,a\,(2p - 1) = 0.1$ | $40\,(2a)^2\,p\,(1 - p) = 0.0023$ | $\sqrt{V_e} = 0.025$ | 0.79 |
| sp2 | 0.5 | $E[G] = 0.3$ | $40\,a_0 = 0.3$ | 0 | $\sqrt{V_e} = 0.025$ | 0 |
| sp3 | 0.7 | $E[G] = 0.5$ | $G_0 = 0.5$ | $\frac{V_0}{20}\left[16 + \sum_{i=1}^{4}(1 + \bar g\, c_i)^2\right] + \sum e_{ij}^2 \left(\frac{V_0}{20}\right)^2 = 0.00072$ | $\sqrt{V_e} = 0.025$ | 0.54 |
| sp4 | 0.9 | $E[G] = 0.7$ | $G_0 = 0.7$ | $40\,V_m = 0.004$ | $\sqrt{V_e} = 0.025$ | 0.86 |

sp3：$V_0 / 20$ 是一个位点 $g$ 的方差，$\bar g = 0.025$，$c_i = \sum_j e_{ij}$ = 20, 20, -60, 20，所以方括号里是 $16 + 3 \times 1.5^2 + 0.5^2 = 23$；$\sum e_{ij}^2 = 6 \times 20^2$。

### 5.1 sp1 -- diallelic，建群个体来自 quanti_init_freq（基准物种：多效性、显性、上位性都不开）

| 参数行 | 值 | 含义 |
|---|---|---|
| `quanti_allele_model` | 'diallelic' | 每个位点两个等位基因，关于 0 对称：bit 0 = $+a$，bit 1 = $-a$ |
| `quanti_allele_value` | 0.008 | $a$ |
| `quanti_init_model` | 未设置 | Nemo 自己抽的建群等位基因，会被下面的 quanti_init_freq 覆盖 |
| `metaibm_nemo_params` `quanti_init_freq` | 0.65625 | 由 MetaIBM 重新实现：每份建群等位基因拷贝各自重抽，以概率 $p$ 取 bit 0（$+a$），否则取 bit 1（$-a$） |

期望值，其中 $40 = 20 \times 2$ 是每个性状的等位基因拷贝数（20 个位点，每个 2 份拷贝）：

$$
E[G] = 40\,a\,(2p - 1) = 0.1, \qquad V_g = 40\,(2a)^2\,p\,(1 - p) = 0.0023, \qquad G_{\max} = 40\,a = 0.32
$$

$G_{\max}$ 超过了最优值 0.3。实测：$E[G] = 0.100$，$V_g$ = 0.0024 / 0.0023（seed 1，900 个建群个体）。

patch1 微位点 (0, 0) 上的建群个体：

```
>>> print(ind)
individual
  metaibm-info : species_id=sp1, age=0, pheno_names=('phenotype-1', 'phenotype-2')
  nemo-info    : id=1, gender=male, phenotype=(0.0800, 0.0662), founder
>>> home habitat 0
trait 1 'phenotype-1': 20 loci at sequence positions 0-19
  bit 0 = allele value +0.0080, bit 1 = allele value -0.0080
  locus         1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16  17  18  19  20
  copy0 bit     1   1   1   0   1   0   1   1   0   0   0   1   1   0   0   0   0   1   0   0
  copy1 bit     0   0   0   0   0   0   1   1   0   1   0   0   0   1   0   1   1   0   0   0
  g (1e-3)      0   0   0  16   0  16 -16 -16  16   0  16   0   0   0  16   0   0   0  16  16
  G = 0.08000   P = 0.08001   P - G = +0.00001
trait 2 'phenotype-2': 20 loci at sequence positions 20-39
  bit 0 = allele value +0.0080, bit 1 = allele value -0.0080
  locus         1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16  17  18  19  20
  copy0 bit     0   0   0   0   1   0   1   0   0   0   0   0   1   0   0   1   1   0   1   1
  copy1 bit     1   0   1   0   1   0   1   0   0   0   0   1   1   0   1   0   0   0   1   1
  g (1e-3)      0  16   0  16 -16  16 -16  16  16  16  16   0 -16  16   0   0   0  16 -16 -16
  G = 0.06400   P = 0.06618   P - G = +0.00218
```

怎么读：
- 两份拷贝都是 bit 0 时 $g = +2a = +0.016$；杂合时 $g = 0$；两份都是 bit 1 时 $g = -2a = -0.016$。
- 性状 1：7 个位点为 +16、2 个为 -16，所以 $G_1 = (7 - 2) \times 0.016 = 0.080$。
- 性状 2：9 个位点为 +16、5 个为 -16，所以 $G_2 = (9 - 5) \times 0.016 = 0.064$。



### 5.2 sp2 -- diallelic_asymmetrical，建群个体来自 quanti_init_model 0，显性

| 参数行 | 值 | 含义 |
|---|---|---|
| `quanti_allele_model` | 'diallelic_asymmetrical' | 每个位点两个任意的等位基因值 |
| `quanti_allele_value` | {{0.0075}{0.0225}} | 两行 {{allele 0}{allele 1}}：bit 0 = 0.0075（祖先型），bit 1 = 0.0225。只写一行 {{b}} 表示 $0$ / $b$，写一个普通数字则表示对称的 $\pm b$ |
| `quanti_init_model` | 0 | 每份建群等位基因拷贝都是 bit 0 |
| `quanti_dominance_effects` | 0.5 | (1) 式 $g = a_1 + a_2 + k\,\lvert a_2 - a_1 \rvert$ 中的 $k$ |
| `quanti_mutation_rate` | 1e-3 | 是其他物种的 10 倍：没有现存变异时，默认的 1e-4 常常让有性运行在突变救援之前就灭绝 |

按 (1) 式，位点值随突变拷贝数为：
- 0 份：$g = 2 \times 0.0075 = 0.015$。
- 1 份（杂合）：$g = 0.0075 + 0.0225 + 0.5 \times (0.0225 - 0.0075) = 0.0375$。
- 2 份：$g = 2 \times 0.0225 = 0.045$。

所以每个建群个体都是 $G = 20 \times 0.015 = 0.3$，$V_g = 0$：建群个体在遗传上完全相同。所有变异都来自突变，而且每次突变都让 $G$ 变大（$0 \to 1$ 份：$+0.0225$；$1 \to 2$ 份：$+0.0075$）。

patch2 微位点 (0, 0) 上的建群个体：

```
>>> print(ind)
individual
  metaibm-info : species_id=sp2, age=0, pheno_names=('phenotype-1', 'phenotype-2')
  nemo-info    : id=901, gender=female, phenotype=(0.3294, 0.3413), founder
>>> home habitat 1
trait 1 'phenotype-1': 20 loci at sequence positions 0-19
  bit 0 = allele value +0.0075
  locus          1    2    3    4    5    6    7    8    9   10   11   12   13   14   15   16   17   18   19   20
  copy0 bit      0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0
  copy1 bit      0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0
  g (1e-3)    15.0 15.0 15.0 15.0 15.0 15.0 15.0 15.0 15.0 15.0 15.0 15.0 15.0 15.0 15.0 15.0 15.0 15.0 15.0 15.0
  G = 0.30000   P = 0.32939   P - G = +0.02939
trait 2 'phenotype-2': 20 loci at sequence positions 20-39
  bit 0 = allele value +0.0075
  locus          1    2    3    4    5    6    7    8    9   10   11   12   13   14   15   16   17   18   19   20
  copy0 bit      0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0
  copy1 bit      0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0
  g (1e-3)    15.0 15.0 15.0 15.0 15.0 15.0 15.0 15.0 15.0 15.0 15.0 15.0 15.0 15.0 15.0 15.0 15.0 15.0 15.0 15.0
  G = 0.30000   P = 0.34134   P - G = +0.04134
```

这里看不出显性：每个位点都是 bit 0 纯合，所以 $G = 0.3$，建群个体之间只有 $\varepsilon$ 不同。

一次真实的突变，来自这个建群个体的克隆：

```
trait1 locus 17 copy0: 0.00750 -> 0.02250 (bit 0 -> 1)
```

这使位点 17 变成杂合：$g$ 从 0.015 变为 0.0375。

显性要等突变积累之后才看得出来。下面这个个体来自 example1 以有性模式运行 200 步之后的 patch2（example1 的默认设置：base_dead_rate 0.1，fitness_wid 0.5，sexual_birth_rate 1）。此时 patch2 有 900 个个体，平均每个个体带 16.5 份突变等位基因拷贝（最多 30 份）。它在微位点 (0, 0)，三种基因型它都有：

```
>>> print(ind)
individual
  metaibm-info : species_id=sp2, age=0, pheno_names=('phenotype-1', 'phenotype-2')
  nemo-info    : id=347924, gender=male, phenotype=(0.4418, 0.5533), mother 336859, father 344773
>>> home habitat 1
trait 1 'phenotype-1': 20 loci at sequence positions 0-19
  bit 0 = allele value +0.0075, bit 1 = allele value +0.0225
  locus          1    2    3    4    5    6    7    8    9   10   11   12   13   14   15   16   17   18   19   20
  copy0 bit      0    1    0    0    0    0    0    0    0    1    0    0    0    0    0    0    0    0    1    0
  copy1 bit      0    1    0    0    1    1    0    1    0    0    0    0    0    0    0    1    0    0    0    0
  g (1e-3)    15.0 45.0 15.0 15.0 37.5 37.5 15.0 37.5 15.0 37.5 15.0 15.0 15.0 15.0 15.0 37.5 15.0 15.0 37.5 15.0
  G = 0.46500   P = 0.44185   P - G = -0.02315
trait 2 'phenotype-2': 20 loci at sequence positions 20-39
  bit 0 = allele value +0.0075, bit 1 = allele value +0.0225
  locus          1    2    3    4    5    6    7    8    9   10   11   12   13   14   15   16   17   18   19   20
  copy0 bit      0    1    0    0    0    1    0    0    0    0    0    1    0    0    1    1    0    0    1    0
  copy1 bit      0    0    0    1    1    0    0    0    0    1    1    1    0    0    0    1    0    0    1    0
  g (1e-3)    15.0 37.5 15.0 37.5 37.5 37.5 15.0 15.0 15.0 37.5 37.5 45.0 15.0 15.0 37.5 45.0 15.0 15.0 45.0 15.0
  G = 0.54750   P = 0.55335   P - G = +0.00585
```

怎么读：
- 以 $10^{-3}$ 为单位：没有突变拷贝的位点 $g = 15.0$，杂合位点 $g = 37.5$，两份拷贝都是 bit 1 的位点 $g = 45.0$（性状 2 的位点 12、16、19）。
- 性状 2 有 7 个杂合位点、10 个祖先型位点、3 个两份拷贝都突变的位点：$G_2 = 7 \times 0.0375 + 10 \times 0.015 + 3 \times 0.045 = 0.5475$。
- 如果没有显性（$k = 0$），每个杂合位点只贡献 $0.0075 + 0.0225 = 0.030$，$G_2$ 会是 $7 \times 0.030 + 10 \times 0.015 + 3 \times 0.045 = 0.495$。

### 5.3 sp3 -- continuous，建群个体来自 quanti_init_model 2，上位性

| 参数行 | 值 | 含义 |
|---|---|---|
| `quanti_allele_model` | 'continuous' | 任意实数值；突变在原值上**加**一个 $\mathcal{N}(0, \sigma^2)$ 的抽样 |
| `quanti_init_model` | 2 | 每份建群等位基因拷贝 $= G_0 / 40 + \mathcal{N}(0,\, V_0 / 40)$，$G_0$、$V_0$ 见下面两行 |
| `quanti_init_value` | {{0.5, 0.5}} | $G_0$：每个性状建群个体的 $E[G]$ |
| `quanti_init_variance` | {{$0.025^2$, $0.025^2$}} | $V_0$：每个性状建群个体的加性 $V_g$ |
| `quanti_epistatic_network` | {{1, 2, 3, 4}} | (2) 式中的 $\mathcal{P}$：位点 1-4 两两交互。位点编号在性状内部计数（1-20），同样的位点对作用于两个性状 |
| `quanti_epistatic_effects` | {{20, -20, 20, -20, 20, -20}} | (2) 式中的 $e_{ij}$，每对一个，顺序为 (1,2) (1,3) (1,4) (2,3) (2,4) (3,4) |
| `quanti_mutation_variance` | $0.025^2$ | $\sigma^2$：一次突变在原值上加一个 $\mathcal{N}(0, 0.025^2)$ 的抽样 |

期望值和实测：
- 一份等位基因拷贝约为 $0.5 / 40 = 0.0125$，所以 $g \approx 0.025$。
- 因此一对位点贡献 $e_{ij}\, g_i\, g_j \approx \pm 20 \times 0.025^2 = \pm 0.0125$。
- 三正三负在均值处相互抵消，所以 $E[G]$ 仍是 0.5。建群个体只是多了非加性方差：实测 $V_g$ = 0.00074 / 0.00081，而加性部分是 $0.025^2$。

patch3 微位点 (0, 0) 上的建群个体：

```
>>> print(ind)
individual
  metaibm-info : species_id=sp3, age=0, pheno_names=('phenotype-1', 'phenotype-2')
  nemo-info    : id=1801, gender=male, phenotype=(0.5726, 0.5373), founder
>>> home habitat 2
trait 1 'phenotype-1': 20 loci at sequence positions 0-19
  locus     1     2     3     4     5     6     7     8     9    10    11    12    13    14    15    16    17    18    19    20
  copy0 11.04 17.67 11.45 12.42 15.71 13.39 11.61 26.60  7.20 20.25 13.28 14.57 12.67 17.85 14.69 13.62 13.33 11.06 16.35  9.86
  copy1  5.03 13.16 12.43 15.16 15.79  9.75 12.26 13.41 14.15 14.98 18.29 15.03 14.59  3.47 16.41 10.84 13.21 20.09 11.61 14.35
  g     16.07 30.83 23.87 27.58 31.50 23.14 23.88 40.01 21.35 35.23 31.57 29.60 27.26 21.32 31.10 24.47 26.54 31.15 27.95 24.21
  (allele values and g in units of 1e-3)
  epistasis e_ij*g_i*g_j:
    (1,2) +20*0.0161*0.0308 = +0.00991,  (1,3) -20*0.0161*0.0239 = -0.00767,  (1,4) +20*0.0161*0.0276 = +0.00886
    (2,3) -20*0.0308*0.0239 = -0.01472,  (2,4) +20*0.0308*0.0276 = +0.01700,  (3,4) -20*0.0239*0.0276 = -0.01317
  sum g = 0.54863, sum epistasis = +0.00021
  G = 0.54883   P = 0.57259   P - G = +0.02376
trait 2 'phenotype-2': 20 loci at sequence positions 20-39
  locus     1     2     3     4     5     6     7     8     9    10    11    12    13    14    15    16    17    18    19    20
  copy0 14.25 13.27 11.48 16.11 16.23  9.06 12.20 13.58 13.18 11.13 11.12 16.11 11.01 14.16 14.48  9.15 16.48 16.65 15.09  9.51
  copy1  7.51  9.36 18.87 10.84 10.87 14.90 10.90  6.09 15.73  8.06 17.11  5.55 12.96  3.26 12.59 16.22  7.01 19.12 13.65 11.87
  g     21.76 22.63 30.36 26.95 27.10 23.96 23.10 19.67 28.92 19.19 28.22 21.66 23.96 17.42 27.06 25.37 23.49 35.77 28.74 21.38
  (allele values and g in units of 1e-3)
  epistasis e_ij*g_i*g_j:
    (1,2) +20*0.0218*0.0226 = +0.00985,  (1,3) -20*0.0218*0.0304 = -0.01321,  (1,4) +20*0.0218*0.0270 = +0.01173
    (2,3) -20*0.0226*0.0304 = -0.01374,  (2,4) +20*0.0226*0.0270 = +0.01220,  (3,4) -20*0.0304*0.0270 = -0.01636
  sum g = 0.49672, sum epistasis = -0.00953
  G = 0.48719   P = 0.53728   P - G = +0.05009
```

怎么读：(2) 式的两个求和就是打印输出里的 `sum g` 和 `sum epistasis`。
- 性状 1：$G_1 = 0.54863 + 0.00021 = 0.54883$。六个位点对的项几乎完全抵消。
- 性状 2：$G_2 = 0.49672 - 0.00953 = 0.48719$。

一次真实的突变，来自这个建群个体的克隆：

```
trait2 locus 9 copy0: 0.01318 -> 0.01066 (bit 0 -> 0)
```

等位基因值变了 -0.00253，即在原值上加了一个 $\mathcal{N}(0, 0.025^2)$ 的抽样。连续模型没有等位基因编号，所以 bit 总是 0。

### 5.4 sp4 -- continuous_in_place，建群个体来自 quanti_init_model 1，full 多效性

| 参数行 | 值 | 含义 |
|---|---|---|
| `quanti_allele_model` | 'continuous_in_place' | 突变把等位基因**替换**成它的效应值，效应从 $\mathcal{N}(0, \sigma^2)$ 抽取 |
| `quanti_loci_per_trait` | 未设置 | full 多效性下不能设：全部 40 个位点都作用于两个性状 |
| `quanti_pleiotropy` | 'full' | 每个位点为每个性状各存一个等位基因值 |
| `quanti_mutation_correlation` | 0.5 | 同一次突变在两个性状上的效应之间的相关系数 |
| `quanti_mutation_variance` | $0.01^2$ | 每个性状的 $\sigma^2$，它同时决定建群 $V_g = 40\,\sigma^2 \approx 0.0632^2$ |
| `quanti_init_model` | 1 | 每个建群等位基因值 $= G_0 / 80$（$G_0$ 见下一行），然后每个位点**加上**一个突变效应（即使在这个 in-place 模型里也是加） |
| `quanti_init_value` | {{0.7, 0.7}} | $G_0$：每个性状建群个体的 $E[G]$ |

- 等价的矩阵写法：去掉 `quanti_mutation_correlation` 和 `quanti_mutation_variance` 两行，改设 `quanti_mutation_matrix` {{$\sigma_1^2$, $\rho\,\sigma_1\sigma_2$}{$\rho\,\sigma_1\sigma_2$, $\sigma_2^2$}} = {{0.0001, 0.00005}{0.00005, 0.0001}}（本例 $\sigma_1 = \sigma_2 = 0.01$，$\rho = 0.5$）。对角线是两个性状各自的方差，非对角线是协方差 $\rho\,\sigma_1\sigma_2 = 0.5 \times 0.01 \times 0.01$，不是相关系数。实测：两种写法的结果完全相同。

建群个体是怎么生成的：
- 每个位点抽一对相关的效应值（性状 1，性状 2）。
- 每个性状的那一半各自随机落在一份拷贝上，与另一个性状无关。
- 所以 $G_1$ 和 $G_2$ 是相关的（实测 $r = 0.52$）。但在同一份拷贝上，两个性状的值不一定来自同一次抽样。
- 运行中的突变不一样：一次突变会在同一份拷贝上同时替换两个性状的值（见下）。

patch4 微位点 (0, 0) 上的建群个体：

```
>>> print(ind)
individual
  metaibm-info : species_id=sp4, age=0, pheno_names=('phenotype-1', 'phenotype-2')
  nemo-info    : id=2701, gender=female, phenotype=(0.6291, 0.5764), founder
>>> home habitat 3
trait 1 'phenotype-1': 40 loci at sequence positions 0, 2, ..., 78 (step 2)
  locus     1     2     3     4     5     6     7     8     9    10    11    12    13    14    15    16    17    18    19    20    21    22    23    24    25    26     27    28    29    30    31     32    33    34    35    36    37    38    39    40
  copy0  3.73  8.75  0.96  8.75  8.45  8.75 -1.05  8.75  8.75 18.82  8.75  9.81  2.72  8.75  6.46  8.75 19.09  8.75 24.84  8.75  9.87  8.75  4.80  8.75  8.75 15.88   8.75 12.56 -4.04  8.75  3.08   8.75 -7.59  8.75 16.29 19.39  8.75  3.79 29.08  8.75
  copy1  8.75 10.76  8.75 19.99  8.75 -7.03  8.75  0.03 -1.59  8.75 -1.92  8.75  8.75  1.89  8.75 18.38  8.75 -2.58  8.75  2.63  8.75 19.79  8.75 -9.39  8.19  8.75  -9.76  8.75  8.75 -0.37  8.75   7.64  8.75 28.68  8.75  8.75 15.91  8.75  8.75  5.18
  g     12.48 19.51  9.71 28.74 17.20  1.72  7.70  8.78  7.16 27.57  6.83 18.56 11.47 10.64 15.21 27.13 27.84  6.17 33.59 11.38 18.62 28.54 13.55 -0.64 16.94 24.63  -1.01 21.31  4.71  8.38 11.83  16.39  1.16 37.43 25.04 28.14 24.66 12.54 37.83 13.93
  (allele values and g in units of 1e-3)
  G = 0.65337   P = 0.62907   P - G = -0.02430
trait 2 'phenotype-2': 40 loci at sequence positions 1, 3, ..., 79 (step 2)
  locus     1     2     3     4     5     6     7     8     9    10    11    12    13    14    15    16    17    18    19    20    21    22    23    24    25    26     27    28    29    30    31     32    33    34    35    36    37    38    39    40
  copy0 10.69 -3.07  7.91 14.98 16.58  8.75 -1.86  1.12  1.26 20.54 -8.57 10.12  1.24  8.75 -1.56  7.39 13.18 -5.21 16.49  6.59  3.13  8.75  5.25  1.05 -3.05  3.69 -15.65 18.38 15.58  8.75  0.21 -11.98  8.75  8.75 11.57  2.80  8.75  8.75 23.64  8.75
  copy1  8.75  8.75  8.75  8.75  8.75 14.46  8.75  8.75  8.75  8.75  8.75  8.75  8.75 13.96  8.75  8.75  8.75  8.75  8.75  8.75  8.75 14.44  8.75  8.75  8.75  8.75   8.75  8.75  8.75  6.48  8.75   8.75  6.65 15.67  8.75  8.75  7.52 -5.18  8.75 18.59
  g     19.44  5.68 16.66 23.73 25.33 23.21  6.89  9.87 10.01 29.29  0.18 18.87  9.99 22.71  7.19 16.14 21.93  3.54 25.24 15.34 11.88 23.19 14.00  9.80  5.70 12.44  -6.90 27.13 24.33 15.23  8.96  -3.23 15.40 24.42 20.32 11.55 16.27  3.57 32.39 27.34
  (allele values and g in units of 1e-3)
  G = 0.60505   P = 0.57639   P - G = -0.02866
```

怎么读：
- 8.75（$0.7 / 80 = 8.75 \times 10^{-3}$）表示这份等位基因拷贝没有收到初始化效应。
- 每个位点上，每个性状恰好有一份拷贝不等于 8.75。
- 位点 2：性状 1 的那一半在拷贝 1 上（10.76），性状 2 的那一半在拷贝 0 上（-3.07），是不同的拷贝。
- 按 (2) 式把 40 个位点加起来：$G_1 = 0.65337$，$G_2 = 0.60505$。

一次真实的突变，来自这个建群个体的克隆：

```
trait1 locus 28 copy0: 0.01256 -> 0.00427 (bit 0 -> 0)
trait2 locus 28 copy0: 0.01838 -> -0.00091 (bit 0 -> 0)
```

一次突变落在一个位点的一份拷贝上，同时替换两个性状的值。
新值以 0 为中心，而建群等位基因拷贝在 0.00875 附近。
所以在 sp4 里，一次突变平均会让那份拷贝的值下降约 0.00875。朝最优值 0.9 的选择只能利用较少见的、往上走的突变。
