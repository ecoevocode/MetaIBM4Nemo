# example2

One species, sp1, lives on a mainland and colonizes four empty islands whose environments differ from the mainland's. Every individual carries two quantitative traits built from many genes.
`example2.py` runs the same model twice, once with asexual and once with sexual reproduction, and shows that sex adapts the islands much faster when the founders carry hidden genetic variation.
Two programs work together: MetaIBM (Python) runs the landscape and the life cycle; Nemo (C++, through pynemo) holds every individual's genes and makes every birth. The Chinese version is in the second half of this file.

After this tutorial you know:
- how genes make a phenotype in Nemo's `continuous` allele model, and how clones and sexual offspring inherit them;
- when Nemo's `quanti_init_model` does not meet your needs, how to customize the founders' initial genes: `TTrait.set_allele_value(pos, strand, float)` or `TTrait.set_allele_bit(pos, strand, bool)`;
- what each step of the MetaIBM-Nemo life cycle does, and how to run it (colonization from the mainland, selection, reproduction and mutation, dispersal among islands, random germination into empty microsites).

1. Genetics basics
2. MetaIBM and Nemo: who does what?
3. The landscape: four islands and a mainland
4. One time step: the life cycle (life cycle events, LCE)
5. Custom heterozygous founders: variation that recombination in sexual reproduction quickly releases in the offspring
6. Results: asexual vs sexual

Every number marked "measured" comes from a real run of example2 (seed 1), not from theory.

## 1. Genetics basics

Every individual is diploid and carries 2 quantitative traits: `phenotype-1` and `phenotype-2`.

| Term | Meaning in example2 |
|---|---|
| locus | one gene position; each trait has 20 loci of its own (40 in all) |
| allele copy | every locus is carried twice: copy 0 comes from the father, copy 1 from the mother (a clone gets both from its mother) |
| allele value $a$ | the number one allele copy adds to its trait; in the `continuous` model any real number |
| locus value $g$ | the two copies of one locus together: $g = a_0 + a_1$ |
| genotypic value $G$ | the sum of a trait's 20 locus values; inherited |
| phenotype $P$ | $P = G + \varepsilon$, with $\varepsilon \sim \mathcal{N}(0, V_e)$ drawn once at birth; selection reads $P$; $\varepsilon$ is not inherited |
| $V_g$, $V_e$, $h^2$ | the variance of $G$ and of $\varepsilon$ among individuals; heritability $h^2 = V_g / (V_g + V_e)$ |
| mutation | at each birth, each allele copy mutates with probability $\mu$; a mutation **adds** a draw of $\mathcal{N}(0, \sigma^2)$ to the allele value |
| recombination | at each locus, a parent passes on copy 0 or copy 1 by its own coin flip (no linkage between loci) |

For trait $d$ (1 or 2) and its loci $l = 1, \dots, 20$:

$$
g_{l,d} = a_{l,0,d} + a_{l,1,d}, \qquad G_d = \sum_{l=1}^{20} g_{l,d}, \qquad P_d = G_d + \varepsilon_d, \quad \varepsilon_d \sim \mathcal{N}(0, V_e) \tag{1}
$$

Three kinds of birth:

| Birth | Call | The newborn's alleles |
|---|---|---|
| founder | `indfactory.make_founder(gender, home_habitat)`, then `configure_heterozygous_founder()` | set by example2's founder recipe (section 5) |
| clone (asexual) | `indfactory.make_clone(mother, home_habitat)` | the mother's two copies at every locus, then mutation |
| sexual offspring | `indfactory.make_offspring(mother, father, gender, home_habitat)` | at every locus one of the father's two copies (-> copy 0) and one of the mother's (-> copy 1), each by a coin flip; then mutation |

Every birth ends with $P$ computed from the new alleles by (1), with a new $\varepsilon$.

## 2. MetaIBM and Nemo: who does what

- One species = one `metaibm.indfactory`. It boots one Nemo simulation from the species' parameters (used only to make individuals, never run) and holds one Nemo `IndFactory`.
- Nemo keeps the genes and does every genetic step in C++: founders, inheritance, recombination, mutation, $P = G + \varepsilon$.
- MetaIBM keeps the landscape and decides who dies, who breeds and who moves. It reads only $P$.
- Nemo's own rule: a dead or unused individual is not deleted but put into the `IndFactory`'s recycling pool (`RecyclingPOOL`); every later birth first takes one individual from the pool, resets it and reuses it, and makes a new one only when the pool is empty. MetaIBM4Nemo already handles this (nothing for the user to do).

## 3. The landscape: four islands and a mainland

| | Patch | Habitat | Size | (env1, env2) mean | env sd | At t = 0 |
|---|---|---|---|---|---|---|
| island | patch0 | h0 | 30 x 30 | (0.1, 0.1) | 0.025 | empty |
| island | patch1 | h1 | 30 x 30 | (0.3, 0.3) | 0.025 | empty |
| island | patch2 | h2 | 30 x 30 | (0.7, 0.7) | 0.025 | empty |
| island | patch3 | h3 | 30 x 30 | (0.9, 0.9) | 0.025 | empty |
| mainland | mainland | h0-h3 | 4 x (30 x 30) | (0.5, 0.5) | 0.025 | 3600 founders of sp1 |

```
┌─────────────────────────────────────────────────────────────────────────────┐
│ MAINLAND: build_mainland(), full at t = 0                                   │
│ ┌────────────────┐ ┌────────────────┐ ┌────────────────┐ ┌────────────────┐ │
│ │ h0             │ │ h1             │ │ h2             │ │ h3             │ │
│ │ hab_index 0    │ │ hab_index 1    │ │ hab_index 2    │ │ hab_index 3    │ │
│ │ global_index   │ │ global_index   │ │ global_index   │ │ global_index   │ │
│ │ 65535          │ │ 65534          │ │ 65533          │ │ 65532          │ │
│ │ 900 founders   │ │ 900 founders   │ │ 900 founders   │ │ 900 founders   │ │
│ └────────────────┘ └────────────────┘ └────────────────┘ └────────────────┘ │
│                                                                             │
└─────────────────────────────────────────────────────────────────────────────┘
                                       │
                                       │  colonization: 10 individuals / step, straight into empty island microsites
                                       │ 
                                       ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│ METACOMMUNITY: build_metacommunity(), 4 islands, empty at t = 0             │
│                                                                             │
│ ┌────────────────┐ ┌────────────────┐ ┌────────────────┐ ┌────────────────┐ │
│ │ patch0         │ │ patch1         │ │ patch2         │ │ patch3         │ │
│ │ patch_index 0  │ │ patch_index 1  │ │ patch_index 2  │ │ patch_index 3  │ │
│ │ h0             │ │ h1             │ │ h2             │ │ h3             │ │
│ │ hab_index 0    │ │ hab_index 0    │ │ hab_index 0    │ │ hab_index 0    │ │
│ │ global_index 0 │ │ global_index 1 │ │ global_index 2 │ │ global_index 3 │ │
│ │ env (0.1, 0.1) │ │ env (0.3, 0.3) │ │ env (0.7, 0.7) │ │ env (0.9, 0.9) │ │
│ │ far: 0.4 away  │ │ near: 0.2 away │ │ near: 0.2 away │ │ far: 0.4 away  │ │
│ └────────────────┘ └────────────────┘ └────────────────┘ └────────────────┘ │
│    ◀──── dispersal: 10% of each island's offspring, to the other 3 ────▶    │
│                                                                             │
└─────────────────────────────────────────────────────────────────────────────┘
```

- `global_index` is a habitat attribute (`add_habitat(..., global_index=...)`), unique over the whole landscape (islands 0-3, mainland 65535-65532). It has two roles:
  - Nemo: the homepatch of every individual born in that habitat (`unsigned short`, 0-65535). It marks the individual's birth habitat; read it with `ind.nemo_indi_obj.getHome()`.
  - MetaIBM: the habitat's global ID as a node of the habitat network, recorded in `global_habitat_id_2_index_dir` ((patch, habitat) -> `global_index`). To find a habitat's row/column in the network matrices for dispersal, look up its (patch, habitat) in that dict, then take `global_habitat_id_ls.index((patch, habitat))`.

## 4. One time step: the life cycle

| # | Where | `example2.py` -> MetaIBM call | What happens |
|---|---|---|---|
| 1 | mainland | `meta_dead_selection(0.1, 0.5, 'niche_gaussian')` | each individual survives with probability $s$ of (2), against its own microsite; the dead are recycled |
| 2 | mainland | `mainland_birth()` -> `meta_mainland_asexual_birth_mutate_germinate(0.5)` / `meta_mainland_sexual_birth_mutate_germinate(1)` | newborns go straight into the habitat's empty microsites. Asexual: $\min(\text{empty}, \lfloor 0.5\,N \rfloor)$ clones of distinct random mothers. Sexual: $\min(\text{empty}, \text{pairs})$ offspring, at most one per pair; the pairs are the habitat's females and males, shuffled and zipped |
| 3 | islands | `meta_dead_selection(0.1, 0.5, 'niche_gaussian')` | as step 1 |
| 4 | islands | `meta_reproduce()` -> `meta_asex_reproduce_mutate_into_offspring_pool(0.5)` / `meta_sex_reproduce_mutate_into_offspring_pool(1)` | newborns go into each habitat's offspring pool, not yet into microsites. Asexual: $0.5\,N$ clones; sexual: one offspring per pair; a fractional part is born by chance |
| 5 | islands | `disperse_among_islands()` -> `dispersal_aomng_patches_from_offspring_pool_to_immigrant_pool(0.1, method='uniform')` | 10% of every island's offspring pool moves into the immigrant pools of the other three islands, a third to each |
| 6 | mainland -> islands | `colonize()` -> asexual: `meta_colonize_from_propagules_rains(mainland, 10)`; sexual: `pairwise_sexual_colonization_from_prpagules_rains(mainland, 10)` | 10 individuals leave the mainland and land straight in empty island microsites. Asexual: 10 random individuals, each to a random empty microsite on any island. Sexual: 5 (female, male) pairs, each pair from one mainland habitat into two empty microsites of one island habitat, since a lone sexual colonizer could never breed. The mainland refills those sites at its next birth |
| 7 | islands | `meta_local_germinate_from_offspring_and_immigrant_pool(is_remove=True)` | offspring and immigrants are shuffled together and fill the remaining empty microsites at random: a lottery, fitness plays no part |
| 8 | islands | `meta_clear_up_offspring_and_immigrant_pool()` | what did not germinate is recycled |
| 9 | islands | `save_snapshot()`, every 10 steps | species id and both phenotypes of every island microsite -> `metaibm_output/*.csv.gz` |

Selection (steps 1 and 3) is the only step that reads the phenotype. An individual with phenotype $(P_1, P_2)$ on a microsite with environment $(e_1, e_2)$ survives with probability

$$
s = (1 - d)\, \exp-\left[\left(\frac{P_1 - e_1}{w}\right)^2 + \left(\frac{P_2 - e_2}{w}\right)^2\right], \qquad d = 0.1,\ w = 0.5 \tag{2}
$$

On the optimum $s = 1 - d = 0.9$. A founder at (0.5, 0.5):

| Habitat | Distance per trait | $s$ |
|---|---|---|
| mainland | 0 | 0.90 |
| patch1, patch2 | 0.2 | $0.9\,e^{-0.4^2 - 0.4^2} = 0.65$ |
| patch0, patch3 | 0.4 | $0.9\,e^{-0.8^2 - 0.8^2} = 0.25$ |

## 5. Custom heterozygous founders

### 5.1 The recipe

`fill_founders()` makes each mainland founder with `make_founder()` and hands it to `configure_heterozygous_founder()`, which rewrites all its alleles (2 traits x 20 loci x 2 copies) with `TTrait.set_allele_value()`; `set_value()` then recomputes the phenotype $P$.

The recipe: every locus is heterozygous, its two copies $c + a$ / $c - a$ ($c = 0.0125$, $a = 0.02$):

$$
a_{l,0} = c + s_l\,a + \eta_{l,0}, \qquad a_{l,1} = c - s_l\,a + \eta_{l,1}
$$

where $s_l = \pm 1$; $\eta \sim \mathcal{N}(0,\ 0.00125/40)$

The effect:

- **Founders differ little.** The $+a$ and $-a$ of a locus cancel: $G = 0.5 + \sum \eta$, $V_g(F_0) = 0.00125$. A clone copies both copies, so asexual offspring also have $V_g = 0.00125$; only mutations add more.
- **Sex releases the differences.** A parent passes one copy per locus at random, $c + a$ or $c - a$ with equal chance. Each of the 40 inherited copies adds $a^2$: $V_g(F_1) = 40\,a^2 + 0.00125 = 0.01725$, about 14 times more after one generation.

| | $V_g$ theory | $V_g$ measured (trait 1 / 2) | $G$ range measured (trait 1 / 2) |
|---|---|---|---|
| founders $F_0$ | 0.00125 | 0.00127 / 0.00125 | 0.376-0.617 / 0.363-0.633 |
| clones of $F_0$ | 0.00125 | 0.00128 / 0.00125 | |
| sexual $F_1$ | 0.01725 | 0.0176 / 0.0171 | 0.080-0.925 / 0.081-0.935 |

Measured: the 3600 founders of the script's asexual run and one clone of each; $F_1$: 3600 more founders with random sexes, paired at random into 1780 (female, male) pairs, one offspring each.

So the asexual islands start from $G$ 0.36-0.63 and must wait for mutations, while one generation of sex already spans 0.08-0.93, past the far optima 0.1 and 0.9.

### 5.2 A real founder

The founder at mainland h0, microsite (0, 0), from the script's asexual run, read with the individual's own methods:

```python
>>> ind = mainland_obj.set['mainland'].set['h0'].set['microsite_individuals'][0][0]
>>> print(ind)
individual
  metaibm-info : species_id=sp1, age=0, pheno_names=('phenotype-1', 'phenotype-2')
  nemo-info    : id=1, gender=female, phenotype=(0.4340, 0.4799), founder
>>> ind.get_trait_alleles('quant', 'get_allele_value', trait_dim=0)   # trait 1: array (20 loci, 2 copies)
>>> ind.get_trait_value('quant', trait_dim=0)                         # trait 1: P = 0.43402
```

The two arrays, with $g$, $G$ and $P$ (allele values and $g$ in units of $10^{-3}$):

```
trait 1 'phenotype-1'
  locus        1      2      3      4      5      6      7      8      9     10     11     12     13     14     15     16     17     18     19     20
  copy0   -11.85  34.37  -9.27  40.41  32.30  29.46 -14.67 -16.59  -5.20 -20.68  35.91  25.46 -12.58  32.97  34.67  -7.39 -14.35 -12.43  38.87   0.92
  copy1    23.80  -4.52  31.41  -4.11  -2.40  -8.94  30.12  32.45  25.53  30.34  -8.12  -4.87  26.58   2.86 -15.23  33.22  37.07  26.41  -8.79  26.34
  g        11.95  29.86  22.14  36.30  29.90  20.52  15.45  15.86  20.33   9.66  27.79  20.59  14.00  35.82  19.43  25.83  22.72  13.98  30.07  27.26
  G = 0.44948   P = 0.43402   P - G = -0.01546
trait 2 'phenotype-2'
  locus        1      2      3      4      5      6      7      8      9     10     11     12     13     14     15     16     17     18     19     20
  copy0   -10.91 -10.09  33.98  34.59   4.61  31.74  -7.94  -8.99  23.01 -13.50  10.10 -13.43  -9.73 -23.16  31.24  34.79   0.48  -1.96  -3.56  -9.30
  copy1    27.93  28.67  -7.90 -13.82  25.25  -9.27  27.56  38.01  -7.13  44.45  30.96  46.00  39.85  21.43 -14.65 -11.93  24.33  39.56  28.70  26.09
  g        17.01  18.58  26.08  20.77  29.86  22.47  19.62  29.02  15.87  30.95  41.06  32.57  30.11  -1.73  16.59  22.86  24.81  37.60  25.14  16.79
  G = 0.47605   P = 0.47987   P - G = +0.00382
```

Reading it (units of $10^{-3}$):
- At every locus one copy sits near $c + a = 32.5$ and the other near $c - a = -7.5$. Trait 1, locus 1: copy 0 = -11.85 (low, $\eta = -4.35$), copy 1 = 23.80 (high, $\eta = -8.70$). At locus 2 the high copy is copy 0: $s_l$ is drawn per locus.
- $g$ stays near $2c = 25$: locus 1 has $g = 11.95 = 25 - 4.35 - 8.70$.
- Adding up the 20 loci by (1): $G_1 = 0.44948$, $G_2 = 0.47605$. $P - G$ is $\varepsilon$.

### 5.3 A real sexual offspring

A pair of founders and their offspring, trait 1 (units of $10^{-3}$). `f0` / `f1`: the child's copy 0 is the father's copy 0 / copy 1; `m0` / `m1`: the child's copy 1 is the mother's copy 0 / copy 1.

```
mother: id=8916, gender=female, phenotype=(0.5123, 0.4640), founder
father: id=8361, gender=male,   phenotype=(0.4163, 0.4637), founder
child:  id=10801, gender=male,  phenotype=(0.5543, 0.5860), mother 8916, father 8361

  locus                1      2      3      4      5      6      7      8      9     10     11     12     13     14     15     16     17     18     19     20
  mother copy0     31.02  31.13  26.29  21.31 -11.61  -4.40  32.14  33.95  32.45  32.67 -13.82  -2.86  38.09  -6.71 -11.17  31.39  32.12 -10.21  -1.83  35.31
  mother copy1     -3.99  -8.20   7.27  -1.36  29.18  29.49  -5.98 -13.79  -7.66 -19.78  39.28  27.76   2.05  28.48  26.83  -6.42   0.15  44.31  29.52 -15.39
  father copy0     33.69  31.42  -5.77  31.05  29.46 -11.31  -6.86  -1.34 -11.17  36.92 -14.64 -10.86  20.20 -14.45  39.64 -12.48 -10.91  31.76  -6.95 -17.93
  father copy1      4.97 -11.11  40.22  -5.33  -6.10  38.24  23.00  31.26  34.52  -9.26  30.32  28.40  -3.67  33.56 -12.50  19.80  26.72  -5.25  27.08  23.97
  child  copy0      4.97 -11.11  40.22  31.05  29.46 -11.31  23.00  31.26  34.52  -9.26  30.32 -10.86  -3.67 -14.45  39.64 -12.48 -10.91  31.76  27.08  23.97
  child  copy1     -3.99  31.13   7.27  -1.36 -11.61  -4.40  -5.98  33.95  -7.66 -19.78  39.28  -2.86  38.09  -6.71  26.83  31.39  32.12  44.31  29.52  35.31
  child copy0 <-      f1     f1     f1     f0     f0     f0     f1     f1     f1     f1     f1     f0     f1     f0     f0     f0     f0     f0     f1     f1
  child copy1 <-      m1     m0     m1     m1     m0     m0     m1     m0     m1     m1     m1     m0     m0     m0     m1     m0     m0     m1     m1     m0
  child  g          0.98  20.02  47.49  29.69  17.85 -15.71  17.02  65.22  26.86 -29.04  69.60 -13.72  34.42 -21.16  66.47  18.91  21.21  76.06  56.60  59.29

  mother G1 = 0.49702  G2 = 0.46289
  father G1 = 0.43831  G2 = 0.44923
  child  G1 = 0.54805  G2 = 0.56047
```

Reading it:
- Every child copy equals one of the parent's two copies: inheritance only picks, by a coin flip per locus (no mutation here).
- The parents' high/low pairs are broken up. Locus 18 got two high copies ($g = 31.76 + 44.31 = 76.06$), locus 10 two low ones ($g = -9.26 - 19.78 = -29.04$), while a founder's $g$ stays near 25.
- The child's $G_1 = 0.548$ and $G_2 = 0.560$ lie above both parents: one generation of sex makes genotypes outside the parents' range.

### 5.4 A real mutation

From clones of the founder of section 5.2 (clone 74 was the first that differed):

```
trait1 locus 12 copy0: 0.02546 -> 0.01273 (change -0.01273)
G1: 0.44948 -> 0.43675
```

The allele value moved by -0.01273, one draw of $\mathcal{N}(0, 0.05^2)$ added to it. Most draws are larger: the sd is 0.05, more than the founders' $\mathrm{sd}(G) = 0.035$, but a mutation hits only about one birth in 125 (80 copies x $\mu = 10^{-4}$).

## 6. Results: asexual vs sexual

Measured, seed 1. "Within 0.05": the first recorded step at which the island's mean phenotype is within 0.05 of the optimum on both traits. "Full": the first recorded step with all 900 microsites occupied.

| Island (optimum) | Mean $P$ at t = 990, asexual | sexual | Within 0.05, asexual | sexual | Full, asexual | sexual |
|---|---|---|---|---|---|---|
| patch0 (0.1) | (0.139, 0.167) | (0.104, 0.105) | not by t = 990 | t = 80 | t = 460 | t = 70 |
| patch1 (0.3) | (0.284, 0.305) | (0.312, 0.309) | t = 310 | t = 40 | t = 90 | t = 40 |
| patch2 (0.7) | (0.671, 0.683) | (0.699, 0.711) | t = 640 | t = 40 | t = 130 | t = 40 |
| patch3 (0.9) | (0.774, 0.751) | (0.887, 0.882) | not by t = 990 | t = 80 | t = 910 | t = 70 |

- Sexual: every island is full by t = 70 and within 0.05 of its optimum by t = 80, the far islands too.
- Asexual: the near islands patch1 and patch2 fill by t = 90 / 130, but reach the 0.05 band only at t = 310 / 640. The far islands patch0 and patch3 hold 29 individuals each at t = 100, fill only at t = 460 / 910, and are still 0.04-0.15 off at t = 990.
- The asexual curves move in jumps: the 5-95% band narrows abruptly when one genotype takes over the island (in `phenotype-1.png`: patch0 near t = 430, patch3 near t = 880).
- A clone copies its mother's alleles, so the asexual islands' phenotypes outside the founders' $G$ range (0.36-0.63), such as patch0's 0.139, can only come from mutations.

The figures in `metaibm_output/`: `abundance.png` (occupied microsites), `phenotype-1.png` and `phenotype-2.png` (mean and 5-95% band of the phenotype), `survival_rate.png` (the survival probability (2) of every individual, against the habitat's mean environment).

---

# example2

一个物种 sp1 生活在大陆上，向四个空岛殖民，岛上的环境和大陆不同。每个个体带有两个由许多基因决定的数量性状。
`example2.py` 把同一个模型跑两遍，一遍无性生殖、一遍有性生殖，展示当建群个体携带隐藏的遗传变异时，有性生殖让岛屿适应得快得多。
两个程序协同工作：MetaIBM（Python）负责景观和生命周期；Nemo（C++，经 pynemo）保存每个个体的基因并完成每一次出生。英文版在本文件前半部分。

读完本教程，你会知道：
- 在 Nemo 的 `continuous` 等位基因模型里，基因怎样决定表型，克隆和有性后代怎样继承基因；
- 当 Nemo 的 `quanti_init_model` 不满足需求时，用户如何定制建群个体（founder）的初始基因：`TTrait.set_allele_value(pos, strand, float)` 或 `TTrait.set_allele_bit(pos, strand, bool)`；
- MetaIBM-Nemo生命周期的每一步做了什么，以及怎样运行（从大陆殖民、选择、繁殖与突变、岛屿间扩散、空微生境里的随机萌发）。

1. 遗传学基础
2. MetaIBM 和 Nemo：谁负责什么？
3. 景观：四个岛和一个大陆
4. 一个时间步：生命周期（life cycle events, LCE）
5. 定制杂合的建群个体：有性生殖中的基因重组过程能使其子代个体快速释放的变异
6. 结果：无性 vs 有性

凡是标"实测"的数字都来自 example2 的真实运行（seed 1），不是理论推算。

## 1. 遗传学基础

每个个体都是二倍体，带有 2 个数量性状：`phenotype-1` 和 `phenotype-2`。

| 术语 | 在 example2 中的含义 |
|---|---|
| 位点（locus） | 一个基因位置；每个性状有 20 个自己的位点（共 40 个） |
| 等位基因拷贝（allele copy） | 每个位点有两份：拷贝 0 来自父亲，拷贝 1 来自母亲（克隆的两份都来自母亲） |
| 等位基因值 $a$ | 一份等位基因拷贝给它所属性状贡献的数值；`continuous` 模型里可以是任意实数 |
| 位点值 $g$ | 一个位点两份拷贝合起来的值：$g = a_0 + a_1$ |
| 基因型值 $G$ | 一个性状 20 个位点值之和；可遗传 |
| 表型 $P$ | $P = G + \varepsilon$，$\varepsilon \sim \mathcal{N}(0, V_e)$ 在出生时抽一次；选择读的是 $P$；$\varepsilon$ 不遗传 |
| $V_g$、$V_e$、$h^2$ | 个体间 $G$ 的方差和 $\varepsilon$ 的方差；遗传力 $h^2 = V_g / (V_g + V_e)$ |
| 突变 | 每次出生时，每份等位基因拷贝以概率 $\mu$ 突变；一次突变在等位基因值上**加**一个 $\mathcal{N}(0, \sigma^2)$ 的抽样 |
| 重组 | 在每个位点，亲本传下拷贝 0 还是拷贝 1 由各自独立的一次掷硬币决定（位点之间无连锁） |

对性状 $d$（1 或 2）及其位点 $l = 1, \dots, 20$：

$$
g_{l,d} = a_{l,0,d} + a_{l,1,d}, \qquad G_d = \sum_{l=1}^{20} g_{l,d}, \qquad P_d = G_d + \varepsilon_d, \quad \varepsilon_d \sim \mathcal{N}(0, V_e) \tag{1}
$$

三种出生方式：

| 出生方式 | 调用 | 新个体的等位基因 |
|---|---|---|
| 建群个体 | `indfactory.make_founder(gender, home_habitat)`，然后 `configure_heterozygous_founder()` | 由 example2 的建群配方设定（第 5 节） |
| 克隆（无性） | `indfactory.make_clone(mother, home_habitat)` | 每个位点都是母亲的两份拷贝，然后突变 |
| 有性后代 | `indfactory.make_offspring(mother, father, gender, home_habitat)` | 每个位点从父亲的两份拷贝中取一份（-> 拷贝 0），从母亲的两份中取一份（-> 拷贝 1），各掷一次硬币；然后突变 |

每次出生最后都会用新的等位基因按 (1) 计算 $P$，并重新抽一个 $\varepsilon$。

## 2. MetaIBM 和 Nemo：谁负责什么

- 一个物种 = 一个 `metaibm.indfactory`。它用该物种的参数启动一个 Nemo 模拟（只用来生成个体，从不运行），并持有一个 Nemo `IndFactory`。
- Nemo 保存基因，并用 C++ 完成所有遗传步骤：建群个体、遗传、重组、突变、$P = G + \varepsilon$。
- MetaIBM 负责景观，决定谁死、谁繁殖、谁移动。它只读取 $P$。
- Nemo 本身的规则：死亡或没用上的个体不会被删除，而是放进 `IndFactory` 的回收池（`RecyclingPOOL`）；之后每次出生先从池里取出一个个体、清空重置后重用，池空了才新建。MetaIBM4Nemo 已适配（无需用户额外操作）。

## 3. 景观：四个岛和一个大陆

| | patch | 栖息地 | 大小 | (env1, env2) 均值 | 环境 sd | t = 0 时 |
|---|---|---|---|---|---|---|
| 岛 | patch0 | h0 | 30 x 30 | (0.1, 0.1) | 0.025 | 空 |
| 岛 | patch1 | h1 | 30 x 30 | (0.3, 0.3) | 0.025 | 空 |
| 岛 | patch2 | h2 | 30 x 30 | (0.7, 0.7) | 0.025 | 空 |
| 岛 | patch3 | h3 | 30 x 30 | (0.9, 0.9) | 0.025 | 空 |
| 大陆 | mainland | h0-h3 | 4 x (30 x 30) | (0.5, 0.5) | 0.025 | 3600 个 sp1 建群个体 |

```
┌─────────────────────────────────────────────────────────────────────────────┐
│ MAINLAND: build_mainland(), full at t = 0                                   │
│ ┌────────────────┐ ┌────────────────┐ ┌────────────────┐ ┌────────────────┐ │
│ │ h0             │ │ h1             │ │ h2             │ │ h3             │ │
│ │ hab_index 0    │ │ hab_index 1    │ │ hab_index 2    │ │ hab_index 3    │ │
│ │ global_index   │ │ global_index   │ │ global_index   │ │ global_index   │ │
│ │ 65535          │ │ 65534          │ │ 65533          │ │ 65532          │ │
│ │ 900 founders   │ │ 900 founders   │ │ 900 founders   │ │ 900 founders   │ │
│ └────────────────┘ └────────────────┘ └────────────────┘ └────────────────┘ │
│                                                                             │
└─────────────────────────────────────────────────────────────────────────────┘
                                       │
                                       │  colonization: 10 individuals / step, straight into empty island microsites
                                       │ 
                                       ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│ METACOMMUNITY: build_metacommunity(), 4 islands, empty at t = 0             │
│                                                                             │
│ ┌────────────────┐ ┌────────────────┐ ┌────────────────┐ ┌────────────────┐ │
│ │ patch0         │ │ patch1         │ │ patch2         │ │ patch3         │ │
│ │ patch_index 0  │ │ patch_index 1  │ │ patch_index 2  │ │ patch_index 3  │ │
│ │ h0             │ │ h1             │ │ h2             │ │ h3             │ │
│ │ hab_index 0    │ │ hab_index 0    │ │ hab_index 0    │ │ hab_index 0    │ │
│ │ global_index 0 │ │ global_index 1 │ │ global_index 2 │ │ global_index 3 │ │
│ │ env (0.1, 0.1) │ │ env (0.3, 0.3) │ │ env (0.7, 0.7) │ │ env (0.9, 0.9) │ │
│ │ far: 0.4 away  │ │ near: 0.2 away │ │ near: 0.2 away │ │ far: 0.4 away  │ │
│ └────────────────┘ └────────────────┘ └────────────────┘ └────────────────┘ │
│    ◀──── dispersal: 10% of each island's offspring, to the other 3 ────▶    │
│                                                                             │
└─────────────────────────────────────────────────────────────────────────────┘
```

- `global_index` 是 habitat 的属性（`add_habitat(..., global_index=...)`），在整个景观里唯一（岛屿 0-3，大陆 65535-65532），有两个作用：
  - Nemo：在这个 habitat 出生的个体，homepatch 就是它（`unsigned short`，取值 0-65535），标记个体的出生 habitat，用 `ind.nemo_indi_obj.getHome()` 读出。
  - MetaIBM：habitat 网络里节点的全局唯一编号，登记在 `global_habitat_id_2_index_dir`（(patch, habitat) -> `global_index`）。扩散时要找某个 habitat 在网络矩阵里的行列号：先在这张表里反查出它的 (patch, habitat)，再用 `global_habitat_id_ls.index((patch, habitat))` 得到行列号。

## 4. 一个时间步：生命周期

| # | 位置 | `example2.py` -> MetaIBM 调用 | 发生什么 |
|---|---|---|---|
| 1 | 大陆 | `meta_dead_selection(0.1, 0.5, 'niche_gaussian')` | 每个个体以 (2) 式的概率 $s$ 存活，对照的是它自己的微生境；死亡个体被回收 |
| 2 | 大陆 | `mainland_birth()` -> `meta_mainland_asexual_birth_mutate_germinate(0.5)` / `meta_mainland_sexual_birth_mutate_germinate(1)` | 新生个体直接进入栖息地的空微生境。无性：$\min(\text{空位}, \lfloor 0.5\,N \rfloor)$ 个克隆，母亲随机抽取、互不重复。有性：$\min(\text{空位}, \text{配对数})$ 个后代，每对最多一个；配对是把栖息地里的雌性和雄性各自打乱后一一配上 |
| 3 | 岛屿 | `meta_dead_selection(0.1, 0.5, 'niche_gaussian')` | 同第 1 步 |
| 4 | 岛屿 | `meta_reproduce()` -> `meta_asex_reproduce_mutate_into_offspring_pool(0.5)` / `meta_sex_reproduce_mutate_into_offspring_pool(1)` | 新生个体进入各栖息地的 offspring pool，还没进入微生境。无性：$0.5\,N$ 个克隆；有性：每对一个后代；小数部分按概率出生 |
| 5 | 岛屿 | `disperse_among_islands()` -> `dispersal_aomng_patches_from_offspring_pool_to_immigrant_pool(0.1, method='uniform')` | 每个岛 offspring pool 的 10% 移到另外三个岛的 immigrant pool，每个岛各三分之一 |
| 6 | 大陆 -> 岛屿 | `colonize()` -> 无性：`meta_colonize_from_propagules_rains(mainland, 10)`；有性：`pairwise_sexual_colonization_from_prpagules_rains(mainland, 10)` | 10 个个体离开大陆，直接落进岛上的空微生境。无性：随机 10 个个体，各自落到任意岛上一个随机空位。有性：5 对（雌，雄），每对来自大陆同一个栖息地，落进同一个岛屿栖息地的两个空位，因为单独一个有性殖民者永远无法繁殖。大陆在下一次出生时补上这些空位 |
| 7 | 岛屿 | `meta_local_germinate_from_offspring_and_immigrant_pool(is_remove=True)` | 后代和移入者混在一起打乱，随机填满剩下的空微生境：一种抽签，与适合度无关 |
| 8 | 岛屿 | `meta_clear_up_offspring_and_immigrant_pool()` | 没有萌发的个体被回收 |
| 9 | 岛屿 | `save_snapshot()`，每 10 步一次 | 每个岛屿微生境的物种编号和两个表型 -> `metaibm_output/*.csv.gz` |

选择（第 1、3 步）是唯一读取表型的步骤。表型为 $(P_1, P_2)$ 的个体，在环境为 $(e_1, e_2)$ 的微生境上的存活概率是

$$
s = (1 - d)\, \exp-\left[\left(\frac{P_1 - e_1}{w}\right)^2 + \left(\frac{P_2 - e_2}{w}\right)^2\right], \qquad d = 0.1,\ w = 0.5 \tag{2}
$$

在最优值上 $s = 1 - d = 0.9$。一个位于 (0.5, 0.5) 的建群个体：

| 栖息地 | 每个性状的距离 | $s$ |
|---|---|---|
| 大陆 | 0 | 0.90 |
| patch1、patch2 | 0.2 | $0.9\,e^{-0.4^2 - 0.4^2} = 0.65$ |
| patch0、patch3 | 0.4 | $0.9\,e^{-0.8^2 - 0.8^2} = 0.25$ |

## 5. 定制杂合的建群个体

### 5.1 配方

`fill_founders()` 先用 `make_founder()` 生成大陆建群个体，再交给 `configure_heterozygous_founder()`：用 `TTrait.set_allele_value()` 重写它的全部等位基因（2 个性状 x 20 个位点 x 2 份拷贝），最后用 `set_value()` 重新计算表型 $P$。

配方：每个位点都杂合，两份拷贝分别为 $c + a$ / $c - a$（$c = 0.0125$，$a = 0.02$)：

$$
a_{l,0} = c + s_l\,a + \eta_{l,0}, \qquad a_{l,1} = c - s_l\,a + \eta_{l,1}
$$

其中，$s_l = \pm 1$；$\eta \sim \mathcal{N}(0,\ 0.00125/40)$

效果：

- **建群个体之间差别很小。**同一位点的 $+a$ 和 $-a$ 互相抵消，$G = 0.5 + \sum \eta$，$V_g(F_0) = 0.00125$。克隆原样复制两份拷贝，所以无性后代的 $V_g$ 也只有 0.00125，之后只能靠突变增加。
- **有性生殖把差别释放出来。**亲本在每个位点随机传下一份拷贝，是 $c + a$ 还是 $c - a$ 各占一半。传下的 40 份拷贝各贡献 $a^2$：$V_g(F_1) = 40\,a^2 + 0.00125 = 0.01725$，一代就变成约 14 倍。

| | $V_g$ 理论 | $V_g$ 实测（性状 1 / 2） | $G$ 范围实测（性状 1 / 2） |
|---|---|---|---|
| 建群个体 $F_0$ | 0.00125 | 0.00127 / 0.00125 | 0.376-0.617 / 0.363-0.633 |
| $F_0$ 的克隆 | 0.00125 | 0.00128 / 0.00125 | |
| 有性 $F_1$ | 0.01725 | 0.0176 / 0.0171 | 0.080-0.925 / 0.081-0.935 |

实测：脚本无性运行的 3600 个建群个体，以及它们各自的一个克隆；$F_1$ 是另外生成的 3600 个随机性别建群个体，随机配成 1780 对，每对一个后代。

所以无性的岛只能从 $G$ 0.36-0.63 出发、等待突变；一代有性生殖就覆盖了 0.08-0.93，越过了远岛的最优值 0.1 和 0.9。

### 5.2 一个真实的建群个体

脚本无性运行中，大陆 h0 微生境 (0, 0) 上的建群个体，用个体自己的方法读取：

```python
>>> ind = mainland_obj.set['mainland'].set['h0'].set['microsite_individuals'][0][0]
>>> print(ind)
individual
  metaibm-info : species_id=sp1, age=0, pheno_names=('phenotype-1', 'phenotype-2')
  nemo-info    : id=1, gender=female, phenotype=(0.4340, 0.4799), founder
>>> ind.get_trait_alleles('quant', 'get_allele_value', trait_dim=0)   # trait 1: array (20 loci, 2 copies)
>>> ind.get_trait_value('quant', trait_dim=0)                         # trait 1: P = 0.43402
```

两个数组，加上 $g$、$G$ 和 $P$（等位基因值和 $g$ 的单位为 $10^{-3}$）：

```
trait 1 'phenotype-1'
  locus        1      2      3      4      5      6      7      8      9     10     11     12     13     14     15     16     17     18     19     20
  copy0   -11.85  34.37  -9.27  40.41  32.30  29.46 -14.67 -16.59  -5.20 -20.68  35.91  25.46 -12.58  32.97  34.67  -7.39 -14.35 -12.43  38.87   0.92
  copy1    23.80  -4.52  31.41  -4.11  -2.40  -8.94  30.12  32.45  25.53  30.34  -8.12  -4.87  26.58   2.86 -15.23  33.22  37.07  26.41  -8.79  26.34
  g        11.95  29.86  22.14  36.30  29.90  20.52  15.45  15.86  20.33   9.66  27.79  20.59  14.00  35.82  19.43  25.83  22.72  13.98  30.07  27.26
  G = 0.44948   P = 0.43402   P - G = -0.01546
trait 2 'phenotype-2'
  locus        1      2      3      4      5      6      7      8      9     10     11     12     13     14     15     16     17     18     19     20
  copy0   -10.91 -10.09  33.98  34.59   4.61  31.74  -7.94  -8.99  23.01 -13.50  10.10 -13.43  -9.73 -23.16  31.24  34.79   0.48  -1.96  -3.56  -9.30
  copy1    27.93  28.67  -7.90 -13.82  25.25  -9.27  27.56  38.01  -7.13  44.45  30.96  46.00  39.85  21.43 -14.65 -11.93  24.33  39.56  28.70  26.09
  g        17.01  18.58  26.08  20.77  29.86  22.47  19.62  29.02  15.87  30.95  41.06  32.57  30.11  -1.73  16.59  22.86  24.81  37.60  25.14  16.79
  G = 0.47605   P = 0.47987   P - G = +0.00382
```

解读（单位 $10^{-3}$）：
- 每个位点都是一份拷贝在 $c + a = 32.5$ 附近，另一份在 $c - a = -7.5$ 附近。性状 1 位点 1：拷贝 0 = -11.85（低，$\eta = -4.35$），拷贝 1 = 23.80（高，$\eta = -8.70$）。位点 2 的高拷贝则是拷贝 0：$s_l$ 每个位点单独抽。
- $g$ 保持在 $2c = 25$ 附近：位点 1 的 $g = 11.95 = 25 - 4.35 - 8.70$。
- 按 (1) 把 20 个位点加起来：$G_1 = 0.44948$，$G_2 = 0.47605$。$P - G$ 就是 $\varepsilon$。

### 5.3 一个真实的有性后代

一对建群个体和它们的后代，性状 1（单位 $10^{-3}$）。`f0` / `f1`：后代的拷贝 0 是父亲的拷贝 0 / 拷贝 1；`m0` / `m1`：后代的拷贝 1 是母亲的拷贝 0 / 拷贝 1。

```
mother: id=8916, gender=female, phenotype=(0.5123, 0.4640), founder
father: id=8361, gender=male,   phenotype=(0.4163, 0.4637), founder
child:  id=10801, gender=male,  phenotype=(0.5543, 0.5860), mother 8916, father 8361

  locus                1      2      3      4      5      6      7      8      9     10     11     12     13     14     15     16     17     18     19     20
  mother copy0     31.02  31.13  26.29  21.31 -11.61  -4.40  32.14  33.95  32.45  32.67 -13.82  -2.86  38.09  -6.71 -11.17  31.39  32.12 -10.21  -1.83  35.31
  mother copy1     -3.99  -8.20   7.27  -1.36  29.18  29.49  -5.98 -13.79  -7.66 -19.78  39.28  27.76   2.05  28.48  26.83  -6.42   0.15  44.31  29.52 -15.39
  father copy0     33.69  31.42  -5.77  31.05  29.46 -11.31  -6.86  -1.34 -11.17  36.92 -14.64 -10.86  20.20 -14.45  39.64 -12.48 -10.91  31.76  -6.95 -17.93
  father copy1      4.97 -11.11  40.22  -5.33  -6.10  38.24  23.00  31.26  34.52  -9.26  30.32  28.40  -3.67  33.56 -12.50  19.80  26.72  -5.25  27.08  23.97
  child  copy0      4.97 -11.11  40.22  31.05  29.46 -11.31  23.00  31.26  34.52  -9.26  30.32 -10.86  -3.67 -14.45  39.64 -12.48 -10.91  31.76  27.08  23.97
  child  copy1     -3.99  31.13   7.27  -1.36 -11.61  -4.40  -5.98  33.95  -7.66 -19.78  39.28  -2.86  38.09  -6.71  26.83  31.39  32.12  44.31  29.52  35.31
  child copy0 <-      f1     f1     f1     f0     f0     f0     f1     f1     f1     f1     f1     f0     f1     f0     f0     f0     f0     f0     f1     f1
  child copy1 <-      m1     m0     m1     m1     m0     m0     m1     m0     m1     m1     m1     m0     m0     m0     m1     m0     m0     m1     m1     m0
  child  g          0.98  20.02  47.49  29.69  17.85 -15.71  17.02  65.22  26.86 -29.04  69.60 -13.72  34.42 -21.16  66.47  18.91  21.21  76.06  56.60  59.29

  mother G1 = 0.49702  G2 = 0.46289
  father G1 = 0.43831  G2 = 0.44923
  child  G1 = 0.54805  G2 = 0.56047
```

解读：
- 后代的每份拷贝都等于亲本两份拷贝中的一份：遗传只是挑选，每个位点掷一次硬币（这里没有突变）。
- 亲本的高低配对被拆开了。位点 18 拿到了两份高拷贝（$g = 31.76 + 44.31 = 76.06$），位点 10 拿到两份低拷贝（$g = -9.26 - 19.78 = -29.04$），而建群个体的 $g$ 一直在 25 附近。
- 后代的 $G_1 = 0.548$、$G_2 = 0.560$ 都高于父母双方：一代有性生殖就能产生超出亲本范围的基因型。

### 5.4 一次真实的突变

来自第 5.2 节那个建群个体的克隆（第 74 个克隆是第一个出现差异的）：

```
trait1 locus 12 copy0: 0.02546 -> 0.01273 (change -0.01273)
G1: 0.44948 -> 0.43675
```

等位基因值变了 -0.01273，即在原值上加了一个 $\mathcal{N}(0, 0.05^2)$ 的抽样。多数抽样比这个大：sd 是 0.05，比建群个体的 $\mathrm{sd}(G) = 0.035$ 还大；但大约每 125 次出生才有一次突变（80 份拷贝 x $\mu = 10^{-4}$）。

## 6. 结果：无性 vs 有性

实测，seed 1。"进入 0.05"：第一个两个性状上岛屿平均表型都离最优值不到 0.05 的记录步。"占满"：第一个 900 个微生境全部被占的记录步。

| 岛（最优值） | t = 990 平均 $P$，无性 | 有性 | 进入 0.05，无性 | 有性 | 占满，无性 | 有性 |
|---|---|---|---|---|---|---|
| patch0 (0.1) | (0.139, 0.167) | (0.104, 0.105) | t = 990 前未达到 | t = 80 | t = 460 | t = 70 |
| patch1 (0.3) | (0.284, 0.305) | (0.312, 0.309) | t = 310 | t = 40 | t = 90 | t = 40 |
| patch2 (0.7) | (0.671, 0.683) | (0.699, 0.711) | t = 640 | t = 40 | t = 130 | t = 40 |
| patch3 (0.9) | (0.774, 0.751) | (0.887, 0.882) | t = 990 前未达到 | t = 80 | t = 910 | t = 70 |

- 有性：每个岛在 t = 70 前占满，t = 80 前进入最优值 0.05 以内，远处的岛也一样。
- 无性：近处的 patch1、patch2 在 t = 90 / 130 占满，但到 t = 310 / 640 才进入 0.05 以内。远处的 patch0、patch3 在 t = 100 时各只有 29 个个体，到 t = 460 / 910 才占满，t = 990 时仍差 0.04-0.15。
- 无性曲线是跳跃式移动的：某个基因型占领全岛时，5-95% 色带会突然变窄（`phenotype-1.png` 中 patch0 在 t = 430 附近，patch3 在 t = 880 附近）。
- 克隆复制母亲的等位基因，所以无性岛屿上超出建群个体 $G$ 范围（0.36-0.63）的表型，例如 patch0 的 0.139，只能来自突变。

`metaibm_output/` 里的图：`abundance.png`（被占的微生境数）、`phenotype-1.png` 和 `phenotype-2.png`（表型均值和 5-95% 色带）、`survival_rate.png`（每个个体按 (2) 式计算的存活概率，环境取栖息地均值）。
