# MetaIBM-Nemo quanti parameters -- bilingual manual

The Chinese version is in the second half of this file.

How the parameters get in: every Nemo row below goes through `pynemo.params.SimParameters` -- one `nemo_params.add(key, value)` per row -- and `nemo_params` goes to
`metaibm.indfactory.indfactory(species_id, pheno_name_ls, nemo_params, ...)`: one species = one indfactory = one Nemo boot sim. A row never added is unset.
Never add a key with `None`: `SimParameters.add(k, None)` writes a bare key, which Nemo reads as the BOOL flag `k 1` -- leave the row out instead.
The rows of section (3) never reach Nemo: they go to `indfactory(..., metaibm_nemo_params={key: value})`.
For example, example1 keeps one `spN.json` per species, every row below in this order, and writes `"_UNSET"` for a row that stays in the file but is never added (`load_species()` in example1.py).

## Tags

- (1) Nemo's Effect ON MetaIBM: ON / OFF / FORBIDDEN
- (2) Requirement of MetaIBM-NEMO: REQ / OPT
- (3) Parameters across species: FIX / VAR (static attributes?)
- Note: `*` = conditional: the tag shown is the safe default, and the line names when it relaxes
- Note: `!` = must (VAR! = must differ across species)
- FORBIDDEN (ON / OFF): ON = setting it does something (registers a writer, or kills the boot); OFF = no effect on its own

| Tag shape (quanti table) | Rows | Meaning |
|---|---|---|
| ON REQ FIX | 1 | Nemo needs it, it shapes the individuals, every species must agree |
| ON REQ FIX\* | 1 | ... and VAR is allowed under the condition that line names |
| ON REQ VAR | 3 | Nemo needs it, it shapes the individuals, each species may set its own |
| ON OPT FIX | 4 | optional; once set it shapes the individuals and every species must agree |
| ON OPT FIX\* | 1 | ... and VAR is allowed under the condition that line names |
| ON OPT VAR | 14 | optional per-species knob -- most of the genetics lives here |
| OFF | 10 | purely inert -- parsed, validated, then never read; the line says why |
| FORBIDDEN BY METAIBM | 7 | Nemo's FileServices: no writer here can ever see MetaIBM's individuals -- leave the whole family unset |

## (1) non_genetic

Builds the boot sim, which nothing here steps -- measured bit-identical either way.

| Parameter | (1) Effect | (2) Req. | (3) Species | Note |
|---|---|---|---|---|
| `replicates` | OFF | REQ | FIX | boot sim only; without it the consistency check fails, no record is built, and no prototype comes out; e.g. 1 |
| `generations` | OFF | REQ | FIX | boot sim only; same failure mode as replicates; e.g. 1 |
| `patch_number` | OFF | REQ | FIX | boot sim only; without it 'param "patch_number" is missing!'; e.g. 1 |
| `patch_capacity` | OFF | REQ | FIX | boot sim only; without it 'population parameters are not properly set!' + TMatrix overflow; e.g. 1 |
| `aging` | OFF | REQ | FIX | boot sim only; To avoid NEMO's warning >>>> WARNING::Metapop::setPopulation: generation 0 is empty!! |
| `run_mode` | OFF | OPT | FIX | file placement + console only; 'silent' turns SILENT_RUN on (no 'setting random seed' / 'outputs:' lines; warnings muted too, \*\*\*ERROR\*\*\* still shows), the .log files are still written. Measured vs 'overwrite' with seed 42: founders identical, logs differ only in start time and this line; e.g. 'overwrite' to see Nemo's warnings again |
| `root_dir` | OFF | OPT | FIX | file placement only, but the directory IS created at boot; e.g. 'nemo_cache' |
| `logfile` | OFF | OPT | FIX | file placement only; No effect in metaibm-nemo; e.g. 'nemo_cache/nemo.log' |
| `filename` | OFF | REQ | VAR! | names this species' .log -- the only one that is per species; sharing it does not error, the species just overwrite each other's log; e.g. 'sp1' |

## (2) quanti

The genetics, which is what MetaIBM reads back off the prototype (trait_keys=('quant',)).

### architecture

| Parameter | (1) Effect | (2) Req. | (3) Species | Note |
|---|---|---|---|---|
| `quanti_traits` | ON | REQ | FIX\* | quantitative traits per individual. |
| `quanti_loci` | ON | REQ | FIX | loci summed over all traits (full pleiotropy: loci that each affect every trait); the genetic map is process-global, so every species keeps the same count -- even at recombination 0.5 |
| `quanti_loci_per_trait` | ON | OPT | VAR | loci of each trait, one row of quanti_traits columns (not with full pleiotropy) |

### allele model

| Parameter | (1) Effect | (2) Req. | (3) Species | Note |
|---|---|---|---|---|
| `quanti_allele_model` | ON | REQ | VAR | per species: diallelic / diallelic_asymmetrical / continuous / continuous_in_place |
| `quanti_allele_value` | ON | REQ | VAR | di-allelic models only: the allele values |
| `quanti_diallele_datatype` | ON | OPT | VAR | di-allelic only: how Nemo stores the sequence, 'bit' (default) \| 'byte'; same individuals either way, except next to quanti_dominance_effects, where 'bit' gives a different phenotype -- dominance on -> 'byte'; e.g. 'byte' |

### founder initialisation

| Parameter | (1) Effect | (2) Req. | (3) Species | Note |
|---|---|---|---|---|
| `quanti_init_model` | ON | OPT | VAR | unset = 0. di-allelic: 0 every allele is allele 0, 1 random alleles (avoid: drawn with RAND::RandBool, biased in this build), 2 alternating homozygous loci; continuous: 0 every allele init_value / 2L, 1 + one mutation per locus, 2 + Gaussian jitter of quanti_init_variance |
| `quanti_init_value` | ON | OPT | VAR | continuous only: genotypic-value centre per trait (Nemo prints an ERROR for the di-allelic models) |
| `quanti_init_variance` | ON | OPT | VAR | continuous only: initial genotypic VARIANCE per trait (init model 2) |
| `quanti_init_freq` | OFF | - | - | but MetaIBM rebuilt it for diallelic(-asymmetric) models without pleiotropy in indfactory.py: a parameter of Nemo's LCE quanti_init, and MetaIBM4Nemo never reach Nemo's LCE. |
| `quanti_init_trait_values` | OFF | - | - | the other half of LCE quanti_init: per-patch trait centres, never read in metaibm-nemo; the per-species centre is quanti_init_value, which IS read at founding |

### genotype to phenotype

| Parameter | (1) Effect | (2) Req. | (3) Species | Note |
|---|---|---|---|---|
| `quanti_environmental_variance` | ON | OPT | VAR | Ve, non-heritable part: genotypic value + N(0, sqrt(Ve)) |
| `quanti_heritability` | OFF | - | - | a target heritability; only Nemo's LCEs (quanti_modifier, breed_quanti) use it, to reset Ve, and no LCE runs in metaibm-nemo -- set Ve itself with quanti_environmental_variance |
| `quanti_heritability_setTime` | OFF | - | - | how often (in generations) LCE breed_quanti re-applies quanti_heritability; never read in metaibm-nemo |
| `quanti_heritability_isBroad` | OFF | - | - | switches the heritability target from narrow (Va) to broad sense (Vg); read only by the quanti stat handler, which never runs in metaibm-nemo |
| `quanti_modifier_at_generation` | OFF | - | - | LCE quanti_modifier: the generation(s), one row, at which it resets Ve to reach quanti_heritability; no LCE runs in metaibm-nemo |
| `quanti_modifier_at_age` | OFF | - | - | LCE quanti_modifier: the age class Va is measured on, 'adults' \| 'offspring' \| 'all' |
| `quanti_modifier_is_permanent` | OFF | - | - | LCE quanti_modifier: turns quanti_modifier_at_generation into a period (1 = every generation) |

### mutation

| Parameter | (1) Effect | (2) Req. | (3) Species | Note |
|---|---|---|---|---|
| `quanti_mutation_rate` | ON | REQ | VAR | per ALLELE COPY per generation |
| `quanti_mutation_variance` | ON | OPT | VAR | VARIANCE of one mutation's effect (continuous adds N(0, sqrt(v)), continuous_in_place replaces with it); the di-allelic models never read it (a mutation just flips the allele) |

### pleiotropy

| Parameter | (1) Effect | (2) Req. | (3) Species | Note |
|---|---|---|---|---|
| `quanti_pleiotropy` | ON | OPT | VAR | unset/'no' (default) one locus one trait; 'full' every locus affects every trait |
| `quanti_mutation_correlation` | ON | OPT | VAR | needs pleiotropy |
| `quanti_mutation_matrix` | ON | OPT | VAR | variance-COVARIANCE matrix of mutation effects |
| `quanti_pleio_matrix` | ON | OPT | VAR | loci-by-traits 0/1 connectivity ('var' pleiotropy; not with epistasis) |

### dominance

| Parameter | (1) Effect | (2) Req. | (3) Species | Note |
|---|---|---|---|---|
| `quanti_dominance_effects` | ON | OPT | VAR | k: locus value a1 + a2 + k \|a2 - a1\| |
| `quanti_dominance_mean` | OFF | - | - | dead in Nemo: declared, read nowhere; draw random k with a config macro instead, quanti_dominance_effects {{rnorm(n, mean, sd)}} |
| `quanti_dominance_sd` | OFF | - | - | dead in Nemo: the other half of quanti_dominance_mean |

### epistasis

| Parameter | (1) Effect | (2) Req. | (3) Species | Note |
|---|---|---|---|---|
| `quanti_epistatic_effects` | ON | OPT | VAR | one coefficient per locus pair: G += e_ij g_i g_j |
| `quanti_epistatic_network` | ON | OPT | VAR | the interacting loci, one network per row |

### genetic map: process-global, keyed by trait type

| Parameter | (1) Effect | (2) Req. | (3) Species | Note |
|---|---|---|---|---|
| `quanti_recombination_rate` | ON | OPT | FIX\* | free recombination (= 0.5) |
| `quanti_genetic_map_resolution` | ON | OPT | FIX | map grain in cM (default 1), read only once a map is active; process-global: Nemo keeps the smallest over all species and rescales the whole map, so keep it identical; e.g. 0.1 |
| `quanti_chromosome_num_locus` | ON | OPT | FIX | loci per chromosome, one row; cannot switch the map on alone -- read next to a matrix quanti_recombination_rate or quanti_random_genetic_map; e.g. {{10, 10, 10, 10}} |
| `quanti_genetic_map` | ON | OPT | FIX | locus positions in cM, one row per chromosome, totalling quanti_loci; switches the map on, not together with quanti_recombination_rate; every species recombines on the FIRST registrant's map; e.g. 40 loci on 2 chromosomes: '{' + '{0, 5, ... , 95}' \* 2 + '}' |
| `quanti_random_genetic_map` | ON | OPT | FIX | chromosome lengths in cM, loci scattered uniformly; switches the map on, shared like quanti_genetic_map (the first registrant's map wins); e.g. {{100, 100}} |

### output: Nemo's FileServices, forbiddened by MetaIBM

| Parameter | (1) Effect | (2) Req. | (3) Species | Note |
|---|---|---|---|---|
| `quanti_output` | FORBIDDEN (ON) | - | - | Nemo's trait file: it only ever sees the boot sim's empty Metapop, never MetaIBM's individuals, and a bad option or a missing quanti_logtime kills the boot. Options: genotypes \| snp_genotypes \| snp_id \| plink |
| `quanti_logtime` | FORBIDDEN (OFF) | - | - | how often quanti_output writes (generations); read only next to it, where it is mandatory |
| `quanti_dir` | FORBIDDEN (OFF) | - | - | directory of the quanti output files, relative to root_dir |
| `quanti_freq_output` | FORBIDDEN (ON) | - | - | per-locus allele frequencies, di-allelic only (inert on continuous models); same empty Metapop, and fatal without quanti_freq_logtime |
| `quanti_freq_logtime` | FORBIDDEN (OFF) | - | - | how often quanti_freq_output writes; read only next to it, where it is mandatory |
| `quanti_ohta_output` | FORBIDDEN (ON) | - | - | Ohta's LD statistics per locus pair; fatal on anything but a di-allelic, single-trait, multi-patch setup |
| `quanti_ohta_logtime` | FORBIDDEN (OFF) | - | - | how often quanti_ohta_output writes; read only next to it, where it is mandatory |

## (3) metaibm_nemo_params

Nemo parameters that MetaIBM re-implements itself (indfactory(metaibm_nemo_params=...)).

| Parameter | (1) Effect | (2) Req. | (3) Species | Note |
|---|---|---|---|---|
| `quanti_init_freq` | ON | OPT | VAR | quanti_init_freq = P(allele 0) per locus, di-allelic models without pleiotropy only |

---

# MetaIBM-Nemo quanti 参数双语手册

英文版在本文件前半部分。

参数怎么传进去：下面每一行 Nemo 参数都经 `pynemo.params.SimParameters` 传入——每行调用一次 `nemo_params.add(key, value)`——再把 `nemo_params` 交给
`metaibm.indfactory.indfactory(species_id, pheno_name_ls, nemo_params, ...)`：一个物种 = 一个 indfactory = 一个 Nemo 启动模拟。从没 add 过的行就是未设置。
不要用 `None` 去 add：`SimParameters.add(k, None)` 会写出一个没有值的裸键，Nemo 会把它读成 BOOL 开关 `k 1`——不想设就直接不 add。
第 (3) 节的参数根本不会传给 Nemo：它们交给 `indfactory(..., metaibm_nemo_params={key: value})`。
举例：example1 每个物种一个 `spN.json`，按下表顺序列出每一行；不 add 的行写 `"_UNSET"` 占位，留在文件里但不会被 add（见 example1.py 的 `load_species()`）。

## 标签

- (1) Nemo 对 MetaIBM 的作用：ON / OFF / FORBIDDEN
- (2) MetaIBM-NEMO 是否需要：REQ / OPT
- (3) 各物种之间：FIX / VAR（静态属性？）
- 注：`*` = 有条件：所示标签是安全的默认值，该行写明何时可以放宽
- 注：`!` = 必须（VAR! = 各物种必须不同）
- FORBIDDEN (ON / OFF)：ON = 设置它会产生作用（注册一个 writer，或让启动失败）；OFF = 单独设置没有作用

| 标签组合（quanti 表） | 行数 | 含义 |
|---|---|---|
| ON REQ FIX | 1 | Nemo 需要它，它会影响生成的个体，所有物种必须一致 |
| ON REQ FIX\* | 1 | ……在该行写明的条件下允许 VAR |
| ON REQ VAR | 3 | Nemo 需要它，它会影响生成的个体，每个物种可以各自设置 |
| ON OPT FIX | 4 | 可选；一旦设置就会影响生成的个体，所有物种必须一致 |
| ON OPT FIX\* | 1 | ……在该行写明的条件下允许 VAR |
| ON OPT VAR | 14 | 可选的物种级旋钮——大部分遗传设定都在这里 |
| OFF | 10 | 完全无效——会被解析、校验，之后再也不读；该行写明原因 |
| FORBIDDEN BY METAIBM | 7 | Nemo 的 FileServices：这里的 writer 永远看不到 MetaIBM 的个体——整组都不要设置 |

## (1) non_genetic

用来构建启动模拟（boot sim），这里没有任何东西推进它——实测设与不设，结果逐位相同。

| 参数 | (1) 作用 | (2) 需要 | (3) 物种间 | 说明 |
|---|---|---|---|---|
| `replicates` | OFF | REQ | FIX | 只用于启动模拟；没有它，一致性检查失败，不会建立 record，也不会产出 prototype；例：1 |
| `generations` | OFF | REQ | FIX | 只用于启动模拟；缺失时的失败方式同 replicates；例：1 |
| `patch_number` | OFF | REQ | FIX | 只用于启动模拟；没有它会报 'param "patch_number" is missing!'；例：1 |
| `patch_capacity` | OFF | REQ | FIX | 只用于启动模拟；没有它会报 'population parameters are not properly set!'，并出现 TMatrix 溢出；例：1 |
| `aging` | OFF | REQ | FIX | 只用于启动模拟；用来避免 Nemo 的警告 >>>> WARNING::Metapop::setPopulation: generation 0 is empty!! |
| `run_mode` | OFF | OPT | FIX | 只影响文件位置和控制台输出；'silent' 会打开 SILENT_RUN（不打印 'setting random seed' / 'outputs:' 行；警告也被静音，\*\*\*ERROR\*\*\* 仍会显示），.log 文件照常写出。用种子 42 与 'overwrite' 对比实测：founder 完全相同，日志只在开始时间和这一行上不同；例：想重新看到 Nemo 的警告就用 'overwrite' |
| `root_dir` | OFF | OPT | FIX | 只影响文件位置，但该目录确实会在启动时创建；例：'nemo_cache' |
| `logfile` | OFF | OPT | FIX | 只影响文件位置；在 metaibm-nemo 中没有作用；例：'nemo_cache/nemo.log' |
| `filename` | OFF | REQ | VAR! | 为该物种的 .log 命名——唯一一个按物种设置的参数；共用同一个名字不会报错，只是各物种会互相覆盖日志；例：'sp1' |

## (2) quanti

遗传部分，也就是 MetaIBM 从 prototype 上读回的内容（trait_keys=('quant',)）。

### 遗传结构（architecture）

| 参数 | (1) 作用 | (2) 需要 | (3) 物种间 | 说明 |
|---|---|---|---|---|
| `quanti_traits` | ON | REQ | FIX\* | 每个个体的数量性状数。 |
| `quanti_loci` | ON | REQ | FIX | 所有性状的位点总数（完全多效性时：每个位点都影响所有性状）；遗传图谱是进程全局的，所以每个物种保持相同的位点数——即使重组率为 0.5 也是如此 |
| `quanti_loci_per_trait` | ON | OPT | VAR | 每个性状的位点数，一行、quanti_traits 列（不能与完全多效性同用） |

### 等位基因模型（allele model）

| 参数 | (1) 作用 | (2) 需要 | (3) 物种间 | 说明 |
|---|---|---|---|---|
| `quanti_allele_model` | ON | REQ | VAR | 按物种设置：diallelic / diallelic_asymmetrical / continuous / continuous_in_place |
| `quanti_allele_value` | ON | REQ | VAR | 仅限双等位模型：等位基因的值 |
| `quanti_diallele_datatype` | ON | OPT | VAR | 仅限双等位：Nemo 如何存储序列，'bit'（默认）\| 'byte'；两种方式得到相同的个体，但与 quanti_dominance_effects 同用时例外，此时 'bit' 会给出不同的表型——开启显性时用 'byte'；例：'byte' |

### 建群个体初始化（founder initialisation）

| 参数 | (1) 作用 | (2) 需要 | (3) 物种间 | 说明 |
|---|---|---|---|---|
| `quanti_init_model` | ON | OPT | VAR | 不设 = 0。双等位：0 所有等位基因都是 allele 0，1 随机等位基因（避免使用：由 RAND::RandBool 抽取，在此构建中有偏），2 交替的纯合位点；连续模型：0 每个等位基因都是 init_value / 2L，1 在此基础上每个位点加一次突变，2 在此基础上加方差为 quanti_init_variance 的高斯扰动 |
| `quanti_init_value` | ON | OPT | VAR | 仅限连续模型：每个性状的基因型值中心（双等位模型下 Nemo 会打印 ERROR） |
| `quanti_init_variance` | ON | OPT | VAR | 仅限连续模型：每个性状的初始基因型方差（init model 2） |
| `quanti_init_freq` | OFF | - | - | 但 MetaIBM 在 indfactory.py 中为无多效性的 diallelic(-asymmetric) 模型重新实现了它：它是 Nemo 的 LCE quanti_init 的参数，而 MetaIBM4Nemo 永远不会走到 Nemo 的 LCE。 |
| `quanti_init_trait_values` | OFF | - | - | LCE quanti_init 的另一半：按 patch 设置的性状中心，在 metaibm-nemo 中从不读取；按物种的中心是 quanti_init_value，它在建群时确实会被读取 |

### 基因型到表型（genotype to phenotype）

| 参数 | (1) 作用 | (2) 需要 | (3) 物种间 | 说明 |
|---|---|---|---|---|
| `quanti_environmental_variance` | ON | OPT | VAR | Ve，不可遗传的部分：基因型值 + N(0, sqrt(Ve)) |
| `quanti_heritability` | OFF | - | - | 目标遗传力；只有 Nemo 的 LCE（quanti_modifier、breed_quanti）用它来重设 Ve，而 metaibm-nemo 中没有 LCE 运行——直接用 quanti_environmental_variance 设定 Ve |
| `quanti_heritability_setTime` | OFF | - | - | LCE breed_quanti 每隔多少代重新应用一次 quanti_heritability；在 metaibm-nemo 中从不读取 |
| `quanti_heritability_isBroad` | OFF | - | - | 把遗传力目标从狭义（Va）切换为广义（Vg）；只有 quanti 统计处理器会读它，而它在 metaibm-nemo 中从不运行 |
| `quanti_modifier_at_generation` | OFF | - | - | LCE quanti_modifier：在哪一代（或哪几代，一行）重设 Ve 以达到 quanti_heritability；metaibm-nemo 中没有 LCE 运行 |
| `quanti_modifier_at_age` | OFF | - | - | LCE quanti_modifier：在哪个年龄类上测量 Va，'adults' \| 'offspring' \| 'all' |
| `quanti_modifier_is_permanent` | OFF | - | - | LCE quanti_modifier：把 quanti_modifier_at_generation 变成周期（1 = 每代） |

### 突变（mutation）

| 参数 | (1) 作用 | (2) 需要 | (3) 物种间 | 说明 |
|---|---|---|---|---|
| `quanti_mutation_rate` | ON | REQ | VAR | 每个等位基因拷贝、每代 |
| `quanti_mutation_variance` | ON | OPT | VAR | 单次突变效应的方差（continuous 加上 N(0, sqrt(v))，continuous_in_place 用它替换原值）；双等位模型不读它（突变只翻转等位基因） |

### 多效性（pleiotropy）

| 参数 | (1) 作用 | (2) 需要 | (3) 物种间 | 说明 |
|---|---|---|---|---|
| `quanti_pleiotropy` | ON | OPT | VAR | 不设/'no'（默认）一个位点只影响一个性状；'full' 每个位点影响所有性状 |
| `quanti_mutation_correlation` | ON | OPT | VAR | 需要多效性 |
| `quanti_mutation_matrix` | ON | OPT | VAR | 突变效应的方差-协方差矩阵 |
| `quanti_pleio_matrix` | ON | OPT | VAR | 位点 × 性状的 0/1 连接矩阵（'var' 多效性；不能与上位性同用） |

### 显性（dominance）

| 参数 | (1) 作用 | (2) 需要 | (3) 物种间 | 说明 |
|---|---|---|---|---|
| `quanti_dominance_effects` | ON | OPT | VAR | k：位点值 a1 + a2 + k \|a2 - a1\| |
| `quanti_dominance_mean` | OFF | - | - | Nemo 中的死参数：声明了，但没有任何地方读取；如需随机 k，改用配置宏 quanti_dominance_effects {{rnorm(n, mean, sd)}} |
| `quanti_dominance_sd` | OFF | - | - | Nemo 中的死参数：quanti_dominance_mean 的另一半 |

### 上位性（epistasis）

| 参数 | (1) 作用 | (2) 需要 | (3) 物种间 | 说明 |
|---|---|---|---|---|
| `quanti_epistatic_effects` | ON | OPT | VAR | 每对位点一个系数：G += e_ij g_i g_j |
| `quanti_epistatic_network` | ON | OPT | VAR | 相互作用的位点，每行一个网络 |

### 遗传图谱：进程全局，按性状类型区分（genetic map: process-global, keyed by trait type）

| 参数 | (1) 作用 | (2) 需要 | (3) 物种间 | 说明 |
|---|---|---|---|---|
| `quanti_recombination_rate` | ON | OPT | FIX\* | 自由重组（= 0.5） |
| `quanti_genetic_map_resolution` | ON | OPT | FIX | 图谱粒度，单位 cM（默认 1），只有图谱启用后才会读取；进程全局：Nemo 取所有物种中最小的值并重新缩放整张图谱，所以要保持一致；例：0.1 |
| `quanti_chromosome_num_locus` | ON | OPT | FIX | 每条染色体的位点数，一行；不能单独启用图谱——要和矩阵形式的 quanti_recombination_rate 或 quanti_random_genetic_map 一起读取；例：{{10, 10, 10, 10}} |
| `quanti_genetic_map` | ON | OPT | FIX | 位点位置，单位 cM，每条染色体一行，总数等于 quanti_loci；会启用图谱，不能与 quanti_recombination_rate 同用；所有物种都在第一个注册者的图谱上重组；例：40 个位点分在 2 条染色体上：'{' + '{0, 5, ... , 95}' \* 2 + '}' |
| `quanti_random_genetic_map` | ON | OPT | FIX | 染色体长度，单位 cM，位点均匀随机散布；会启用图谱，共享方式同 quanti_genetic_map（第一个注册者的图谱生效）；例：{{100, 100}} |

### 输出：Nemo 的 FileServices，被 MetaIBM 禁用（output: Nemo's FileServices, forbiddened by MetaIBM）

| 参数 | (1) 作用 | (2) 需要 | (3) 物种间 | 说明 |
|---|---|---|---|---|
| `quanti_output` | FORBIDDEN (ON) | - | - | Nemo 的性状文件：它只能看到启动模拟中空的 Metapop，永远看不到 MetaIBM 的个体；选项写错或缺少 quanti_logtime 都会让启动失败。选项：genotypes \| snp_genotypes \| snp_id \| plink |
| `quanti_logtime` | FORBIDDEN (OFF) | - | - | quanti_output 每隔多少代写一次；只在设了 quanti_output 时读取，此时必须设置 |
| `quanti_dir` | FORBIDDEN (OFF) | - | - | quanti 输出文件的目录，相对于 root_dir |
| `quanti_freq_output` | FORBIDDEN (ON) | - | - | 每个位点的等位基因频率，仅限双等位（连续模型下无效）；同样只看到空的 Metapop，且缺少 quanti_freq_logtime 时致命 |
| `quanti_freq_logtime` | FORBIDDEN (OFF) | - | - | quanti_freq_output 每隔多少代写一次；只在设了它时读取，此时必须设置 |
| `quanti_ohta_output` | FORBIDDEN (ON) | - | - | 每对位点的 Ohta 连锁不平衡统计量；除了双等位、单性状、多 patch 的设置以外都会致命 |
| `quanti_ohta_logtime` | FORBIDDEN (OFF) | - | - | quanti_ohta_output 每隔多少代写一次；只在设了它时读取，此时必须设置 |

## (3) metaibm_nemo_params

MetaIBM 自己重新实现的 Nemo 参数（indfactory(metaibm_nemo_params=...)）。

| 参数 | (1) 作用 | (2) 需要 | (3) 物种间 | 说明 |
|---|---|---|---|---|
| `quanti_init_freq` | ON | OPT | VAR | quanti_init_freq = 每个位点的 P(allele 0)，仅限无多效性的双等位模型 |
