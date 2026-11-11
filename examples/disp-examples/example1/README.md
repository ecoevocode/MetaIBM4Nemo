# MetaIBM: dispersal driven by the disp trait

Once a species carries the fdisp / mdisp traits, each individual's disp value v can take part in dispersal. MetaIBM computes the rates and the allocation as before (rate matrix, kernel, `matrix_around`); only "who leaves, how many leave, and whether they die on the way" changes. The new parameters default to the old behaviour. The Chinese version is in the second half of this file.

- Section 1: how it is computed
- Section 2: usage
- Section 3: TODO

## 1. How it is computed

### 1.1 The four schemes

| Scheme | In plain words | Probability that individual k goes from i to j |
|---|---|---|
| 1: old MetaIBM | m sets how many leave; every individual is equally likely to be picked | m_ij (number m_ij × N_i) |
| 2: number unchanged | m sets the number of slots; individuals with higher v are more likely to be picked | ≈ m_ij × N_i × v_k / Σv (number m_ij × N_i unchanged) |
| 3: v decides leaving, m allocates | v decides whether an individual leaves; leavers are allocated by m | among: v_k × m_a × Normalize(D)_ij; within: v_k × m_w × Normalize(D^w)_hh′ |
| 4: 3 + death on the way | as scheme 3; each leaver dies on the way with probability c (recycled) | among: v_k × m_a × Normalize(D)_ij × (1 − c_a); within: v_k × m_w × Normalize(D^w)_hh′ × (1 − c_w) |

`sex_specific` can be used in schemes 2–4. Section 2.2 gives the parameters of each scheme.

### 1.2 Notation

- i, j: source and destination: patches, or habitats in the network version / within; h, h′: habitats; H: the number of habitats in a patch.
- k: an individual; v_k: the disp value of individual k, read by `get_trait_value`, in 0–1, the mean of its two alleles.
- N_i: the number of individuals in source i's pool (offspring_pool; some functions also include dormancy_pool).
- m: the dispersal rate parameter. among: `total_disp_among_rate` (m_a); within: `disp_within_rate` (m_w).
- d_ij: the distance between i and j.
- D: the dispersal kernel matrix; D_ij is the relative strength of i→j.
- Normalize(D): D after normalization; the off-diagonal elements of each row sum to 1.
- M: the rate matrix; element m_ij is the dispersal rate i→j.
- E = [e_ij]: the matrix of migrant numbers; e_ij is the number of migrants i→j (rounded to an integer).
- P(k: i → j): the probability that individual k in source i's pool is moved to j in one dispersal call. Summed over all j ≠ i, it is the probability that k leaves i; 1 minus this sum is the probability that k stays in i.

### 1.3 Rates and numbers: distance matrix → D → Normalize(D) → M → E

All four schemes share this step. among, the network version and within use the same formulas; only the "units", the kernel and m differ. $\mathbf{x}_i$ is the coordinate of unit i.

$$
\mathrm{Distance}_{ij} = \lVert \mathbf{x}_i - \mathbf{x}_j \rVert_2
$$

$$
D = \mathrm{kernel}(\mathrm{Distance}), \qquad D_{ii} = 0
$$

$$
\mathrm{Normalize(D)}_{ij} = \frac{D_{ij}}{\sum_{j'} D_{ij'}}, \quad \text{i.e. } \texttt{D / D.sum(axis=1)}
$$

$$
M = m \cdot \mathrm{Normalize(D)} , M_{ii}=1-m
$$

The three cases:

| | among | network version | within |
|---|---|---|---|
| units i, j | patches | the habitats of all patches | the H habitats of this patch |
| $\mathbf{x}_i$ | patch coordinates | habitat coordinates | habitat coordinates |
| kernel | chosen by `method` | chosen by `method` | `'uniform'` |
| elements of D set to 0 | $D_{ii}$ | $D_{ii}$, and $D_{ij}$ where i and j belong to the same patch | $D_{ii}$ |
| m | m_a (`total_disp_among_rate`) | m_a (`total_disp_among_rate`) | m_w (`disp_within_rate`) |

The network version's rate matrix is 0 between habitats of the same patch and is normalized over the habitats of the other patches, so each habitat's total rate to other patches is m_a, the same meaning as in the patch-level among.

**Migrant numbers**: the number of individuals in each source pool (P is the number of units)

$$
\mathbf{N} = \begin{pmatrix} N_1 \\ N_2 \\ \vdots \\ N_P \end{pmatrix}, \qquad \mathrm{diag}(\mathbf{N}) = \begin{pmatrix} N_1 & 0 & \cdots & 0 \\ 0 & N_2 & \cdots & 0 \\ \vdots & \vdots & \ddots & \vdots \\ 0 & 0 & \cdots & N_P \end{pmatrix}
$$

$$
E = \mathrm{matrix\_around}\big(\mathrm{diag}(\mathbf{N})\, M\big)
$$

`matrix_around` rounds each row randomly to integers and keeps each row sum. The functions that move migrants straight into microsites also take the element-wise minimum with the number of empty sites. The dispersal loop skips i = j (the network version also skips habitats of the same patch), so the diagonal of M and the network version's same-patch elements are never used.

The functions:

| Quantity | among | network version | within |
|---|---|---|---|
| Distance | `pairwise_patch_distance_matrix` | `global_habitat_distance_matrix` | — |
| D | `calculate_dispersal_kernel_strength_matrix` | `calculate_global_habitat_network_dispersal_kernel_strength_matrix` | — |
| Normalize(D) | `normalized_calculate_dispersal_kernel_strength_matrix` | `normalized_calculate_global_habitat_network_dispersal_kernel_strength_matrix` | — |
| M | `get_disp_among_rate_matrix(m)` | `global_habitat_dispersal_among_rate_matrix(m)` | `get_patch_dispersal_within_rate_matrix(m)` |

The $\mathrm{kernel}$ is chosen by `method`:

| `method` | D_ij | Parameters |
|---|---|---|
| `'uniform'` | 1 | none |
| `'gaussian'` | exp(−d_ij² / (2σ²)) | sigma |
| `'exponential'` | ρ × exp(−ρ × d_ij) | rho |
| `'cauchy'` | 1 / (1 + (d_ij / γ)²) | gamma |
| `'power_law'` | (1 + d_ij / r0)^(−α) | alpha, r0 |

After normalization the prefactor of D (e.g. ρ in exponential) cancels; only the shape matters.

### 1.4 Scheme 1: old MetaIBM

1. **Count**: E = [e_ij] from the rate matrix and the kernel (section 1.3).
2. **Pick**: from source i's pool, `random.sample` picks e_ij individuals, all equally likely.

The loop runs over destinations j (outer) and sources i (inner), one draw per pair (i, j). fdisp / mdisp are only inherited and mutated; they do not take part in dispersal. Schemes 2–4 keep this loop and replace only the "pick" step.

### 1.5 Scheme 2: number unchanged

$$
E = \mathrm{matrix\_around}\big(\mathrm{diag}(\mathbf{N})\, M\big)
$$

For each $e_{ij}$, draw $e_{ij}$ individuals from i's pool; individual k is drawn with probability $v_k / \sum_{l \in i} v_l$: one at a time, renormalized over the v of the remaining individuals after each draw. Within one $e_{ij}$ nobody is drawn twice; across destinations, `is_remove` decides whether the drawn individuals are put back (section 1.9). Individuals with v = 0 are never drawn; if fewer than $e_{ij}$ individuals have v > 0, only those are drawn.

$$
P(k: i \to j) \approx \underset{\substack{\big\downarrow \\ \text{rate i→j}}}{m_{ij}} \cdot \underset{\substack{\big\downarrow \\ \text{chance that one slot picks k}}}{\frac{v_k}{\sum_{l \in i} v_l}} \cdot \underset{\substack{\big\downarrow \\ \text{individuals in source i}}}{N_i}
$$

### 1.6 Scheme 3: v decides leaving, m allocates

Take $m_a + m_w = 1$. Each individual k in i's pool leaves with probability $v_k \cdot m$ and is then allocated to some j by $\mathrm{Normalize(D)}_{i\cdot}$; otherwise it stays in i.

$$
P(k: i \to j) = \underset{\substack{\big\downarrow \\ \text{k's own chance to leave}}}{v_k} \cdot \underset{\substack{\big\downarrow \\ \text{share of leavers going among / within: } m_a \text{ / } m_w}}{m} \cdot \underset{\substack{\big\downarrow \\ \text{share of leavers allocated to j}}}{\mathrm{Normalize(D)}_{ij}}
$$

$P(k \text{ stays in } i) = 1 - v_k m_w - v_k m_a = 1 - v_k$.

In practice: one test per pair (i, j). Each individual in i's pool draws a uniform random number $u_k$ and goes to j if $u_k < v_k \cdot m_{ij}$; individuals with $v_k = 0$ always stay. With `is_remove=True` the destinations j are tested in order, and individuals that have already left take no part in later j, so the realized P(k: i → j) is slightly below the formula above (section 1.10).

### 1.7 Scheme 4: scheme 3 + death on the way

Take $m_a + m_w = 1$. Each individual k in i's pool leaves with probability $v_k \cdot m$, is allocated to some j by $\mathrm{Normalize(D)}_{i\cdot}$, and then dies on the way with probability c (removed from its pool and recycled); otherwise it stays in i.

$$
P(k: i \to j) = \underset{\substack{\big\downarrow \\ \text{k's own chance to leave}}}{v_k} \cdot \underset{\substack{\big\downarrow \\ \text{share of leavers going among / within: } m_a \text{ / } m_w}}{m} \cdot \underset{\substack{\big\downarrow \\ \text{share of leavers allocated to j}}}{\mathrm{Normalize(D)}_{ij}} \cdot \underset{\substack{\big\downarrow \\ \text{chance to survive the journey}}}{(1 - c)}
$$

$P(k \text{ dies on the way}) = v_k m_a c_a + v_k m_w c_w$, $P(k \text{ stays in } i) = 1 - v_k$.

- c: among takes $c_a$, within takes $c_w$; they can be equal. By sex: `disp_cost=(c_fem, c_mal)`.
- As in Nemo's `evoldisp`: first test whether k goes to j; an individual going to j draws one more random number and dies on the way with probability c.
- `disp_cost=None` (default): no death on the way; `disp_cost=0` gives the same result.
- With `is_remove=False`, a copy that dies is never made, and the original is not affected.

### 1.8 `sex_specific`

- True: females use fdisp and males use mdisp, as in Nemo's default `exec_evol2sex`.
- False: every individual uses fdisp; Nemo has no such mode.
- Under asexual reproduction every individual is female, so this parameter has no effect.

### 1.9 The two cases of `is_remove`

| | `is_remove=True` (without replacement) | `is_remove=False` (with replacement; copies made by `make_copy()` disperse) |
|---|---|---|
| The original | leaves its source pool | stays in its source pool |
| Drawn from | the remaining pool, for each pair (i, j) | the full pool, for each pair (i, j) |
| Scheme 2 | an individual disperses at most once | no repeats within one j; can be copied to several j |
| Schemes 3/4 | an individual goes to at most one j | copied to each j independently with $v_k \cdot m_{ij}$; can be copied to several j |
| Death in scheme 4 | the dead are removed from the pool and recycled | the dead copy is never made; the original is not affected |

### 1.10 Effects of handling the destinations in order (`is_remove=True`)

The loop runs over destinations j (outer) and sources i (inner), so one source is drawn for each j in turn. With `is_remove=True`, a later destination faces a pool from which the earlier destinations have already taken individuals. One cause, two effects:

- **Destination-order bias**: scheme 1 draws uniformly and is not affected. Schemes 2–4 weight by v, so the individuals taken first are biased towards high v: earlier destinations receive migrants with higher v and later ones with lower v; the larger m, the stronger the bias.
  - Scheme 2: the fixed slots compete; who is drawn earlier changes every individual's chance later.
  - Schemes 3/4: $P(k: i \to j) = v_k\, m_{ij} \prod_{j' \text{ before } j} (1 - v_k\, m_{ij'})$, slightly below $v_k m_{ij}$.
- **Leave-rate gap** (schemes 3/4): $P(k \text{ leaves } i) = 1 - \prod_j (1 - v_k\, m_{ij}) < v_k\, m$; the larger m, the larger the gap (largest at m = 1). The product does not depend on the order, so neither does the gap; it equals the sum of what the later destinations miss.

Among first, then within, is the same kind of sequential handling: within faces the pool left after among, and the within leave probability falls short by $m_a m_w v_k^2$ (the two-scale error).

With `is_remove=False`, every pair (i, j) draws from the full pool, and none of these effects occur. Section 3 lists the fixes.

## 2. Usage

### 2.1 Supported functions

The 7 dispersal functions all take the same new parameters and are all called on `meta.`:

| Scale | Function | Source pool → target |
|---|---|---|
| among | `dispersal_among_patches_from_offspring_pool_and_dormancy_pool` | offspring + dormancy → empty microsites of the destination |
| among | `dispersal_aomng_patches_from_offspring_pool_to_immigrant_pool` | offspring → immigrant pool |
| among | `dispersal_among_patches_from_offsrping_pool_and_dormancy_pool_to_immigrant_pool` | offspring + dormancy → immigrant pool |
| among (network version) | `dispersal_among_patches_in_global_habitat_network_from_offspring_pool_to_immigrant_pool` | offspring → immigrant pool |
| within | `meta_dispersal_within_patch_from_offspring_to_immigrant_pool` | offspring → immigrant pool |
| within | `meta_dispersal_within_patch_from_offspring_and_dormancy_to_immigrant_pool` | offspring + dormancy → immigrant pool |
| within | `meta_dispersal_within_patch_from_offspring_and_dormancy_pool` | offspring + dormancy → empty microsites of the destination |

- Signature of among (the network version included): `(total_disp_among_rate, method='uniform', is_remove=True, use_disp_trait=False, sex_specific=False, is_rescale=True, disp_cost=None, **kwargs)`; the kernel parameters go in `**kwargs`.
- Signature of within: `(disp_within_rate, is_remove=True, use_disp_trait=False, sex_specific=False, is_rescale=True, disp_cost=None)`; the `meta_` functions pass the parameters unchanged to each patch's corresponding function.
- The network version is installed on `metacommunity` by `extension/global_habitat_network.py`; importing metaibm is enough.
- Functions that move migrants straight into microsites: migrants left without an empty site are recycled.

### 2.2 Parameters

| Parameter | Default | Effect |
|---|---|---|
| `use_disp_trait` | `False` | False: scheme 1, the three parameters below have no effect; True: sample by v (schemes 2–4) |
| `is_rescale` | `True` | True: scheme 2, the migrant numbers E are unchanged; False: schemes 3/4, each individual is tested by M |
| `disp_cost` | `None` | Only with `use_disp_trait=True, is_rescale=False`. None: no death on the way (scheme 3); c or `(c_fem, c_mal)`, in 0–1: the probability of dying on the way (scheme 4) |
| `sex_specific` | `False` | True: males use mdisp, females fdisp; False: all use fdisp (section 1.8) |
| `is_remove` | `True` | Existing parameter. True: migrants leave their source pool; False: copies disperse, the originals stay (section 1.9) |

| Scheme | `use_disp_trait` | `is_rescale` | `disp_cost` | Rates |
|---|---|---|---|---|
| 1 | `False` | — | — | as before |
| 2 | `True` | `True` | — (no effect) | as before |
| 3 | `True` | `False` | `None` | m_a + m_w = 1; one scale only: m = 1 |
| 4 | `True` | `False` | c or `(c_fem, c_mal)` | as in scheme 3 |

### 2.3 Example

```python
# among with patches as nodes; the global habitat network version (habitats as nodes, skips habitats of the same patch) takes the same parameters:
#     meta.dispersal_among_patches_in_global_habitat_network_from_offspring_pool_to_immigrant_pool(total_disp_among_rate, method='uniform', is_remove=True, use_disp_trait=False, sex_specific=False, is_rescale=True, disp_cost=None, **kwargs)
# kernel kwargs by method: gaussian sigma; exponential rho; cauchy gamma; power_law alpha, r0; uniform none

# scheme 1 (old MetaIBM): the new kwargs keep their defaults
meta.dispersal_aomng_patches_from_offspring_pool_to_immigrant_pool(total_disp_among_rate=m_among, method='gaussian', is_remove=True, use_disp_trait=False, sex_specific=False, is_rescale=True, disp_cost=None, sigma=1.0)
meta.meta_dispersal_within_patch_from_offspring_to_immigrant_pool(disp_within_rate=m_within, is_remove=True, use_disp_trait=False, sex_specific=False, is_rescale=True, disp_cost=None)

# scheme 2: count unchanged, v picks who
meta.dispersal_aomng_patches_from_offspring_pool_to_immigrant_pool(total_disp_among_rate=m_among, method='gaussian', is_remove=True, use_disp_trait=True, sex_specific=True, is_rescale=True, disp_cost=None, sigma=1.0)
meta.meta_dispersal_within_patch_from_offspring_to_immigrant_pool(disp_within_rate=m_within, is_remove=True, use_disp_trait=True, sex_specific=True, is_rescale=True, disp_cost=None)

# scheme 3: v decides leaving, m allocates (m_a + m_w = 1; one scale: m = 1)
meta.dispersal_aomng_patches_from_offspring_pool_to_immigrant_pool(total_disp_among_rate=m_a, method='gaussian', is_remove=True, use_disp_trait=True, sex_specific=True, is_rescale=False, disp_cost=None, sigma=1.0)
meta.meta_dispersal_within_patch_from_offspring_to_immigrant_pool(disp_within_rate=m_w, is_remove=True, use_disp_trait=True, sex_specific=True, is_rescale=False, disp_cost=None)

# scheme 4: scheme 3 + death on the way
meta.dispersal_aomng_patches_from_offspring_pool_to_immigrant_pool(total_disp_among_rate=m_a, method='gaussian', is_remove=True, use_disp_trait=True, sex_specific=True, is_rescale=False, disp_cost=c_a, sigma=1.0)
meta.meta_dispersal_within_patch_from_offspring_to_immigrant_pool(disp_within_rate=m_w, is_remove=True, use_disp_trait=True, sex_specific=True, is_rescale=False, disp_cost=c_w)
```

### 2.4 Notes

- among and within are two calls, each with its own parameters; in scheme 2, `use_disp_trait` can be switched on or off for among and within separately.
- With `use_disp_trait=True`, every species in the pool must have `'fdisp'` and `'mdisp'` in its `trait_keys`; otherwise reading v raises KeyError.
- With `is_remove=True`, schemes 2–4 have the effects of section 1.10.

## 3. TODO

- **Destination-order bias** (`is_remove=True`, schemes 2–4): later destinations receive migrants with lower v; drawing the destinations of one source in random order, instead of in sequence, removes it.
- **Leave-rate gap** (schemes 3/4, `is_remove=True`): $P(k \text{ leaves } i) = 1 - \prod_j (1 - v_k\, m_{ij}) < v_k\, m$; the larger m, the larger the gap (largest at m = 1). The gap does not depend on the order, so random order cannot remove it. The fix is one draw per source: each individual leaves with $v_k\, m$ and is then allocated to j by $m_{ij} / m$. This removes both effects.
- **Two-scale error** (among first, then within): the within leave probability falls short by $m_a m_w v_k^2$.
  - Correction set by the user: within takes $m_w' = \min\big(1,\ m_w \cdot V_{\text{before}} / V_{\text{after}}\big)$, $V = \sum_k v_k$, computed before and after among.
  - Exact approach: a single two-scale step: leaving is tested once, then the leavers are allocated among all within and among destinations.
- **No individual moved on error** (to be decided): a species without fdisp / mdisp raises KeyError when v is read, and earlier pairs (i, j) may already have moved individuals. To guarantee that no individual moves, check before the j loop that every species in every source pool has fdisp / mdisp.
- **Non-uniform kernel for within**: within uses only uniform now (kernel ≡ 1); the formulas of section 1.3 are written for any kernel, so only the within kernel needs to change.
- **Distance-dependent cost** $c_{ij}$.
- **Weighting the marker pool and the mainland propagule rain by v**: a marker is only `(patch_id, h_id, mode)` and carries no traits.

---

# MetaIBM：用 disp 性状控制 dispersal

物种设置了 fdisp / mdisp 后，可以让个体的 disp 值 v 参与 dispersal。MetaIBM 计算迁移率和分配的算法（rate 矩阵、kernel、`matrix_around`）不变，只改“谁走、走几个、路上是否死亡”。新参数的默认值等于原来的行为。英文版在本文件前半部分。

- 第 1 节：计算原理
- 第 2 节：用法
- 第 3 节：TODO

## 1. 计算原理

### 1.1 四个方案

| 方案 | 大白话 | 个体 k 从 i 去 j 的概率 |
|---|---|---|
| 一：原 MetaIBM | 走几个由 m 决定，等概率挑 | m_ij（数量 m_ij × N_i） |
| 二：数量不变 | 名额由 m 决定，v 高的优先 | ≈ m_ij × N_i × v_k / Σv（数量 m_ij × N_i 不变） |
| 三：v 决定离开，m 决定分配 | v 决定走不走，走出来的按 m 分配 | among：v_k × m_a × Normalize(D)_ij；within：v_k × m_w × Normalize(D^w)_hh′ |
| 四：三 + 路上死亡 | 同方案三，走出来的个体以 c 的概率死在路上（recycle） | among：v_k × m_a × Normalize(D)_ij × (1 − c_a)；within：v_k × m_w × Normalize(D^w)_hh′ × (1 − c_w) |

`sex_specific` 在方案二到四里都可以用。各方案对应的参数见 2.2 节。

### 1.2 符号

- i、j：来源和目的地，指 patch，或网络版 / within 里的 habitat；h、h′：habitat；H：一个 patch 的 habitat 数。
- k：个体；v_k：个体 k 的 disp 值，由 `get_trait_value` 读取，取值 0–1，等于两个等位基因的平均。
- N_i：来源 i 的 pool（offspring_pool，有的函数还包括 dormancy_pool）里的个体数。
- m：迁移率参数。among 是 `total_disp_among_rate`（m_a），within 是 `disp_within_rate`（m_w）。
- d_ij：i 与 j 之间的距离。
- D：dispersal kernel 矩阵，D_ij 是 i→j 的相对强度。
- Normalize(D)：归一化（normalize）后的 D，每一行的非对角元素加起来等于 1。
- M：rate 矩阵，元素 m_ij 是 i→j 的迁移率。
- E = [e_ij]：迁移数矩阵，e_ij 是 i→j 的迁移数（取整后的整数）。
- P(k: i → j)：一次 dispersal 中，来源 i 的 pool 里的个体 k 被迁到 j 的概率。对所有 j ≠ i 求和，就是 k 离开 i 的概率；1 减去这个和，就是 k 留在 i 的概率。

### 1.3 迁移率与迁移数：distance 矩阵 → D → Normalize(D) → M → E

四个方案共用这一步。among、网络版、within 用同一组公式，只是“单元”、kernel 和 m 不同。记 $\mathbf{x}_i$ 为单元 i 的坐标。

$$
\mathrm{Distance}_{ij} = \lVert \mathbf{x}_i - \mathbf{x}_j \rVert_2
$$

$$
D = \mathrm{kernel}(\mathrm{Distance}), \qquad D_{ii} = 0
$$

$$
\mathrm{Normalize(D)}_{ij} = \frac{D_{ij}}{\sum_{j'} D_{ij'}}, \quad \text{i.e. } \texttt{D / D.sum(axis=1)}
$$

$$
M = m \cdot \mathrm{Normalize(D)} , M_{ii}=1-m
$$

三种情况代入：

| | among | 网络版 | within |
|---|---|---|---|
| 单元 i、j | patch | 所有 patch 的 habitat | 本 patch 内的 H 个 habitat |
| $\mathbf{x}_i$ | patch 坐标 | habitat 坐标 | habitat 坐标 |
| kernel | 由 `method` 选择 | 由 `method` 选择 | `'uniform'` |
| D 中置 0 的元素 | $D_{ii}$ | $D_{ii}$，以及 i、j 属于同一 patch 的 $D_{ij}$ | $D_{ii}$ |
| m | m_a（`total_disp_among_rate`） | m_a（`total_disp_among_rate`） | m_w（`disp_within_rate`） |

网络版的 rate 矩阵在同一 patch 的 habitat 之间为 0，按其他 patch 的 habitat 归一化，所以每个 habitat 去其他 patch 的总迁移率就是 m_a，和 patch 版 among 的含义一样。

**迁移数**：各来源 pool 的个体数（P 为单元的个数）

$$
\mathbf{N} = \begin{pmatrix} N_1 \\ N_2 \\ \vdots \\ N_P \end{pmatrix}, \qquad \mathrm{diag}(\mathbf{N}) = \begin{pmatrix} N_1 & 0 & \cdots & 0 \\ 0 & N_2 & \cdots & 0 \\ \vdots & \vdots & \ddots & \vdots \\ 0 & 0 & \cdots & N_P \end{pmatrix}
$$

$$
E = \mathrm{matrix\_around}\big(\mathrm{diag}(\mathbf{N})\, M\big)
$$

`matrix_around` 按行随机取整，保持每一行的和不变。直接进 microsite 的函数还要再和空位上限逐元素取较小值。dispersal 循环跳过 i = j（网络版还跳过同一 patch 的 habitat），所以 M 的对角线和网络版同一 patch 的元素都不会被用到。

对应的函数：

| 量 | among | 网络版 | within |
|---|---|---|---|
| Distance | `pairwise_patch_distance_matrix` | `global_habitat_distance_matrix` | — |
| D | `calculate_dispersal_kernel_strength_matrix` | `calculate_global_habitat_network_dispersal_kernel_strength_matrix` | — |
| Normalize(D) | `normalized_calculate_dispersal_kernel_strength_matrix` | `normalized_calculate_global_habitat_network_dispersal_kernel_strength_matrix` | — |
| M | `get_disp_among_rate_matrix(m)` | `global_habitat_dispersal_among_rate_matrix(m)` | `get_patch_dispersal_within_rate_matrix(m)` |

$\mathrm{kernel}$ 由 `method` 选择：

| `method` | D_ij | 参数 |
|---|---|---|
| `'uniform'` | 1 | 无 |
| `'gaussian'` | exp(−d_ij² / (2σ²)) | sigma |
| `'exponential'` | ρ × exp(−ρ × d_ij) | rho |
| `'cauchy'` | 1 / (1 + (d_ij / γ)²) | gamma |
| `'power_law'` | (1 + d_ij / r0)^(−α) | alpha、r0 |

归一化之后，D 的前置系数（例如 exponential 的 ρ）会被约掉，只有形状起作用。

### 1.4 方案一：原 MetaIBM

1. **算数量**：由 rate 矩阵和 kernel 算出 E = [e_ij]（1.3 节）。
2. **挑个体**：从来源 i 的 pool 里，用 `random.sample` 等概率挑出 e_ij 个。

循环是外层目的地 j、内层来源 i，每一对 (i, j) 抽一次。fdisp / mdisp 只被遗传、突变，不参与 dispersal。方案二到四保留这个循环，只换掉“挑个体”这一步。

### 1.5 方案二：数量不变

$$
E = \mathrm{matrix\_around}\big(\mathrm{diag}(\mathbf{N})\, M\big)
$$

对每个 $e_{ij}$，从 i 的 pool 里抽 $e_{ij}$ 个个体，个体 k 被抽中的概率为 $v_k / \sum_{l \in i} v_l$：一个一个抽，每次按剩下个体的 v 重新归一化。同一个 $e_{ij}$ 里不重复；不同目的地之间放不放回由 `is_remove` 控制（1.9 节）。v = 0 的个体不会被抽到；v > 0 的个体不够 $e_{ij}$ 个时，只抽这些。

$$
P(k: i \to j) \approx \underset{\substack{\big\downarrow \\ \text{i→j 的迁移率}}}{m_{ij}} \cdot \underset{\substack{\big\downarrow \\ \text{每个名额抽中 k 的概率}}}{\frac{v_k}{\sum_{l \in i} v_l}} \cdot \underset{\substack{\big\downarrow \\ \text{来源 i 的个体数}}}{N_i}
$$

### 1.6 方案三：v 决定离开，m 决定分配

取 $m_a + m_w = 1$。对 i 的 pool 里每个个体 k：以概率 $v_k \cdot m$ 离开，离开后按 $\mathrm{Normalize(D)}_{i\cdot}$ 分配到某个 j；否则留在 i。

$$
P(k: i \to j) = \underset{\substack{\big\downarrow \\ \text{k 自己离开的概率}}}{v_k} \cdot \underset{\substack{\big\downarrow \\ \text{离开者里去 among / within 的比例：} m_a \text{ / } m_w}}{m} \cdot \underset{\substack{\big\downarrow \\ \text{离开者分到 j 的比例}}}{\mathrm{Normalize(D)}_{ij}}
$$

$P(k \text{ 留在 } i) = 1 - v_k m_w - v_k m_a = 1 - v_k$。

实际计算：每对 (i, j) 判断一次，i 的 pool 里每个个体取一个均匀随机数 $u_k$，$u_k < v_k \cdot m_{ij}$ 就去 j；$v_k = 0$ 的个体一定留下。`is_remove=True` 时按 j 的顺序依次判断，已经走掉的个体不再参与后面的 j，所以实际的 P(k: i → j) 略小于上式（1.10 节）。

### 1.7 方案四：方案三 + 路上死亡

取 $m_a + m_w = 1$。对 i 的 pool 里每个个体 k：以概率 $v_k \cdot m$ 离开，离开后按 $\mathrm{Normalize(D)}_{i\cdot}$ 分配到某个 j，再以概率 c 死在路上（移出 pool 并 recycle）；否则留在 i。

$$
P(k: i \to j) = \underset{\substack{\big\downarrow \\ \text{k 自己离开的概率}}}{v_k} \cdot \underset{\substack{\big\downarrow \\ \text{离开者里去 among / within 的比例：} m_a \text{ / } m_w}}{m} \cdot \underset{\substack{\big\downarrow \\ \text{离开者分到 j 的比例}}}{\mathrm{Normalize(D)}_{ij}} \cdot \underset{\substack{\big\downarrow \\ \text{路上存活的概率}}}{(1 - c)}
$$

$P(k \text{ 死在路上}) = v_k m_a c_a + v_k m_w c_w$，$P(k \text{ 留在 } i) = 1 - v_k$。

- c：among 传 $c_a$，within 传 $c_w$，可以相同；按性别写成 `disp_cost=(c_fem, c_mal)`。
- 和 Nemo 的 `evoldisp` 一样：先判断去不去 j，去 j 的个体再取一个随机数，以 c 判定死在路上。
- `disp_cost=None`（默认）时没有路上死亡；`disp_cost=0` 结果相同。
- `is_remove=False` 时，死亡的复制品不生成，原个体不受影响。

### 1.8 `sex_specific`

- True：雌性用 fdisp，雄性用 mdisp，等于 Nemo 默认的 `exec_evol2sex`。
- False：所有个体都用 fdisp，Nemo 没有这种模式。
- 无性繁殖时所有个体都是雌性，这个参数不起作用。

### 1.9 `is_remove` 的两种情况

| | `is_remove=True`（不放回） | `is_remove=False`（放回，迁移的是 `make_copy()`） |
|---|---|---|
| 原个体 | 离开来源 pool | 留在来源 pool |
| 怎么抽 | 每对 (i, j) 从剩下的 pool 里抽 | 每对 (i, j) 从完整的 pool 里抽 |
| 方案二 | 一个个体最多迁移一次 | 同一个 j 不重复，可以被复制到多个 j |
| 方案三、四 | 一个个体最多去一个 j | 每个 j 独立以 $v_k \cdot m_{ij}$ 复制，可以被复制到多个 j |
| 方案四的死亡 | 死亡个体移出 pool 并 recycle | 死亡的复制品不生成，原个体不受影响 |

### 1.10 依次处理目的地的影响（`is_remove=True`）

循环外层是目的地 j、内层是来源 i，同一个来源按 j 的顺序依次被抽。`is_remove=True` 时，后轮到的目的地面对的是已经被前面拿走一部分的 pool。同一个原因有两个表现：

- **目的地顺序偏差**：方案一等概率抽样，没有影响。方案二到四按 v 加权，被拿走的偏向 v 高的个体，所以先轮到的目的地收到的移民 v 偏高，后轮到的偏低，m 越大越明显。
  - 方案二：固定名额之间的竞争，前面抽走谁会改变后面每个个体的概率。
  - 方案三、四：$P(k: i \to j) = v_k\, m_{ij} \prod_{j' \text{ 在 } j \text{ 之前}} (1 - v_k\, m_{ij'})$，略小于 $v_k m_{ij}$。
- **离开率缺口**（方案三、四）：$P(k \text{ 离开 } i) = 1 - \prod_j (1 - v_k\, m_{ij}) < v_k\, m$，m 越大越明显（m = 1 时最明显）。乘积和顺序无关，所以缺口和顺序无关；它等于后轮到的目的地少收的部分之和。

先 among、后 within 也是这样依次处理：within 面对的是 among 之后剩下的 pool，within 的离开概率少了 $m_a m_w v_k^2$（两尺度误差）。

`is_remove=False` 时，每对 (i, j) 都从完整的 pool 里抽，没有这些偏差。处理办法见第 3 节。

## 2. 用法

### 2.1 支持的函数

7 个 dispersal 函数都加了同样的新参数，都用 `meta.` 调用：

| 尺度 | 函数 | 来源 pool → 去处 |
|---|---|---|
| among | `dispersal_among_patches_from_offspring_pool_and_dormancy_pool` | offspring + dormancy → 目的地的空 microsite |
| among | `dispersal_aomng_patches_from_offspring_pool_to_immigrant_pool` | offspring → immigrant pool |
| among | `dispersal_among_patches_from_offsrping_pool_and_dormancy_pool_to_immigrant_pool` | offspring + dormancy → immigrant pool |
| among（网络版） | `dispersal_among_patches_in_global_habitat_network_from_offspring_pool_to_immigrant_pool` | offspring → immigrant pool |
| within | `meta_dispersal_within_patch_from_offspring_to_immigrant_pool` | offspring → immigrant pool |
| within | `meta_dispersal_within_patch_from_offspring_and_dormancy_to_immigrant_pool` | offspring + dormancy → immigrant pool |
| within | `meta_dispersal_within_patch_from_offspring_and_dormancy_pool` | offspring + dormancy → 目的地的空 microsite |

- among（含网络版）的签名：`(total_disp_among_rate, method='uniform', is_remove=True, use_disp_trait=False, sex_specific=False, is_rescale=True, disp_cost=None, **kwargs)`，kernel 参数放在 `**kwargs`。
- within 的签名：`(disp_within_rate, is_remove=True, use_disp_trait=False, sex_specific=False, is_rescale=True, disp_cost=None)`，`meta_` 函数把参数原样传给每个 patch 的对应函数。
- 网络版函数由 `extension/global_habitat_network.py` 装到 `metacommunity` 上，导入 metaibm 即可使用。
- 直接进 microsite 的函数：空位不够时，多出的移民 recycle。

### 2.2 参数

| 参数 | 默认值 | 作用 |
|---|---|---|
| `use_disp_trait` | `False` | False：方案一，下面三个参数都不起作用；True：用 v 抽样（方案二到四） |
| `is_rescale` | `True` | True：方案二，迁移数 E 不变；False：方案三、四，按 M 对每个个体判断 |
| `disp_cost` | `None` | 只在 `use_disp_trait=True, is_rescale=False` 时起作用。None：没有路上死亡（方案三）；c 或 `(c_fem, c_mal)`，取值 0–1：路上死亡的概率（方案四） |
| `sex_specific` | `False` | True：雄性用 mdisp，雌性用 fdisp；False：都用 fdisp（1.8 节） |
| `is_remove` | `True` | 原有参数。True：移民离开来源 pool；False：迁出复制品，原个体留下（1.9 节） |

| 方案 | `use_disp_trait` | `is_rescale` | `disp_cost` | 迁移率 |
|---|---|---|---|---|
| 一 | `False` | — | — | 同原来 |
| 二 | `True` | `True` | —（不起作用） | 同原来 |
| 三 | `True` | `False` | `None` | m_a + m_w = 1；只用一个尺度时 m = 1 |
| 四 | `True` | `False` | c 或 `(c_fem, c_mal)` | 同方案三 |

### 2.3 例子

```python
# among with patches as nodes; the global habitat network version (habitats as nodes, skips habitats of the same patch) takes the same parameters:
#     meta.dispersal_among_patches_in_global_habitat_network_from_offspring_pool_to_immigrant_pool(total_disp_among_rate, method='uniform', is_remove=True, use_disp_trait=False, sex_specific=False, is_rescale=True, disp_cost=None, **kwargs)
# kernel kwargs by method: gaussian sigma; exponential rho; cauchy gamma; power_law alpha, r0; uniform none

# scheme 1 (old MetaIBM): the new kwargs keep their defaults
meta.dispersal_aomng_patches_from_offspring_pool_to_immigrant_pool(total_disp_among_rate=m_among, method='gaussian', is_remove=True, use_disp_trait=False, sex_specific=False, is_rescale=True, disp_cost=None, sigma=1.0)
meta.meta_dispersal_within_patch_from_offspring_to_immigrant_pool(disp_within_rate=m_within, is_remove=True, use_disp_trait=False, sex_specific=False, is_rescale=True, disp_cost=None)

# scheme 2: count unchanged, v picks who
meta.dispersal_aomng_patches_from_offspring_pool_to_immigrant_pool(total_disp_among_rate=m_among, method='gaussian', is_remove=True, use_disp_trait=True, sex_specific=True, is_rescale=True, disp_cost=None, sigma=1.0)
meta.meta_dispersal_within_patch_from_offspring_to_immigrant_pool(disp_within_rate=m_within, is_remove=True, use_disp_trait=True, sex_specific=True, is_rescale=True, disp_cost=None)

# scheme 3: v decides leaving, m allocates (m_a + m_w = 1; one scale: m = 1)
meta.dispersal_aomng_patches_from_offspring_pool_to_immigrant_pool(total_disp_among_rate=m_a, method='gaussian', is_remove=True, use_disp_trait=True, sex_specific=True, is_rescale=False, disp_cost=None, sigma=1.0)
meta.meta_dispersal_within_patch_from_offspring_to_immigrant_pool(disp_within_rate=m_w, is_remove=True, use_disp_trait=True, sex_specific=True, is_rescale=False, disp_cost=None)

# scheme 4: scheme 3 + death on the way
meta.dispersal_aomng_patches_from_offspring_pool_to_immigrant_pool(total_disp_among_rate=m_a, method='gaussian', is_remove=True, use_disp_trait=True, sex_specific=True, is_rescale=False, disp_cost=c_a, sigma=1.0)
meta.meta_dispersal_within_patch_from_offspring_to_immigrant_pool(disp_within_rate=m_w, is_remove=True, use_disp_trait=True, sex_specific=True, is_rescale=False, disp_cost=c_w)
```

### 2.4 注意

- among 和 within 是两次调用，参数各自设置；方案二里 among 和 within 可以分别打开或关闭 `use_disp_trait`。
- `use_disp_trait=True` 时，pool 里每个物种的 `trait_keys` 都必须包含 `'fdisp'` 和 `'mdisp'`，否则读 v 时会抛出 KeyError。
- 方案二到四在 `is_remove=True` 时有 1.10 节的偏差。

## 3. TODO

- **目的地顺序偏差**（`is_remove=True`，方案二到四）：后轮到的目的地收到的移民 v 偏低；把同一来源按 j 的顺序依次抽改为乱序执行可以消除。
- **离开率缺口**（方案三、四，`is_remove=True`）：$P(k \text{ 离开 } i) = 1 - \prod_j (1 - v_k\, m_{ij}) < v_k\, m$，m 越大越明显（m = 1 时最明显）。缺口和顺序无关，乱序执行消不掉；要改成每个来源只抽一次签：每个个体以 $v_k\, m$ 离开，再按 $m_{ij} / m$ 分配到 j。这样两项都能消除。
- **两尺度误差**（先 among、后 within）：within 的离开概率少了 $m_a m_w v_k^2$。
  - 用户赋值矫正：within 传 $m_w' = \min\big(1,\ m_w \cdot V_{\text{before}} / V_{\text{after}}\big)$，$V = \sum_k v_k$，分别在 among 前、后计算。
  - 精确做法：两尺度单步，离开只判断一次，再在 within 和 among 的所有目的地之间分配。
- **报错时不移动个体**（待确认）：缺 fdisp / mdisp 时读 v 抛 KeyError，排在前面的 (i, j) 可能已经移动了个体。要保证没有个体被移动，需要在 j 循环之前先确认所有来源 pool 的物种都有 fdisp / mdisp。
- **within 的非 uniform kernel**：within 现在只用 uniform（kernel ≡ 1），1.3 节的式子已按通用 kernel 写好，以后只需换掉 within 的 kernel。
- **随距离变化的代价** $c_{ij}$。
- **marker pool 和 mainland 繁殖体雨按 v 加权**：marker 只是 `(patch_id, h_id, mode)`，没有性状。
