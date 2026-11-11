# MetaIBM-Nemo disp parameters -- bilingual manual

How the parameters get in: every Nemo row below goes through `pynemo.params.SimParameters` -- one `nemo_params.add(key, value)` per row -- and `nemo_params` goes to
`metaibm.indfactory.indfactory(species_id, pheno_name_ls, nemo_params, trait_keys=('fdisp', 'mdisp'), ...)`: one species = one indfactory = one Nemo boot sim. A row never added is unset.
Never add a key with `None`: `SimParameters.add(k, None)` writes a bare key, which Nemo reads as `k 1` (e.g. `add('disp_init_rate', None)` gives every founder allele 1.0) -- leave the row out instead.
Nemo's dispersal life-cycle rows (`disperse`, `disperse_evoldisp`, `dispersal_*`, ...) are not listed: MetaIBM never runs them, so leave them unset; section (3) shows which of them MetaIBM's own dispersal provides.
For example, example1/use-disp-traits.ipynb writes the rows inline: `non_genetic` + `filename` + `mutation` + one `init` dict per species (`disp_params()`).

## Tags

- (1) Nemo's Effect ON MetaIBM: ON / OFF
- (2) Requirement of MetaIBM-NEMO: REQ / OPT
- (3) Parameters across species: FIX / VAR (static attributes?)
- Note: `*` = conditional: the tag shown is the safe default, and the line names when it relaxes
- Note: `!` = must (VAR! = must differ across species)

| Tag shape (disp table) | Rows | Meaning |
|---|---|---|
| ON REQ VAR | 2 | Nemo needs it to build the disp traits, it shapes the individuals, each species may set its own |
| ON OPT VAR | 5 | optional per-species knob: how the founders are drawn |

## (1) non_genetic

Builds the boot sim, which nothing here steps -- the individuals are the same whatever these rows hold, except where the line says.

| Parameter | (1) Effect | (2) Req. | (3) Species | Note |
|---|---|---|---|---|
| `replicates` | OFF | REQ | FIX | boot sim only; without it Nemo waits 10 s, builds no disp trait, and indfactory raises ValueError; e.g. 1 |
| `generations` | OFF | REQ | FIX | boot sim only; same failure mode as replicates; e.g. 1 |
| `patch_number` | OFF | REQ | FIX | boot sim only; without it 'param "patch_number" is missing!' and indfactory raises ValueError; e.g. 1 |
| `patch_capacity` | OFF\* | REQ | FIX | keep 1: with 2 or more, next to `aging` and random founders (disp_init_distribution, or no init row), the same seed no longer gives the same individuals (their distribution is unchanged); without it Python is killed ('population parameters are not properly set!'); e.g. 1 |
| `aging` | OFF | REQ | FIX | boot sim only; To avoid NEMO's warning >>>> WARNING::Metapop::setPopulation: generation 0 is empty!! |
| `run_mode` | OFF | OPT | FIX | console and files only; 'silent': no 'outputs:' line, the .log files are still written; 'create_init' also writes `<filename>.ini` into the working directory; an unknown value: indfactory raises ValueError; e.g. 'overwrite' to see Nemo's warnings |
| `root_dir` | OFF | OPT | FIX | where the .log files go; the directory IS created at boot; unset: the working directory; e.g. 'nemo_cache' |
| `logfile` | OFF | OPT | FIX | Nemo's run log; unset: nemo.log in the working directory; its directory must already exist, otherwise '\*\*\*ERROR\*\*\* could not create simulation logfile' and no log (the boot goes on); e.g. 'nemo_cache/nemo.log' |
| `filename` | OFF | REQ | VAR! | names this species' .log -- the only one that is per species; without it the same failure as replicates; sharing it does not error, the species just overwrite each other's log; e.g. 'sp1' |

## (2) disp

The genetics, which is what MetaIBM reads back off the prototype (trait_keys=('fdisp', 'mdisp')).
Every individual carries both traits, whatever its own sex: fdisp (a female's dispersal rate) and mdisp (a male's). One locus each, two alleles in [0, 1] (maternal, paternal); the trait value = their mean.
The rows always build both traits, so name both in trait_keys (one alone raises ValueError). Each species may set its own rows: two species in one run keep their own founders and mutation.

### mutation

| Parameter | (1) Effect | (2) Req. | (3) Species | Note |
|---|---|---|---|---|
| `disp_mutation_rate` | ON | REQ | VAR | μ, per ALLELE COPY per birth (sexual offspring or clone; founders and make_copy() never mutate): a newborn's trait changes with probability 1 - e^(-2μ), μ 0.02 -> about 4 %; one rate for both traits; without it (or without disp_mutation_mean) Nemo waits 10 s per trait, builds no disp trait, and indfactory raises ValueError; outside [0, 1] Python is killed ('argument is out of bounds'); e.g. 0.02 |
| `disp_mutation_mean` | ON | REQ | VAR | m: one mutation moves one allele up or down by a step ~ Exponential(mean m), so the mean step is m; a step past 0 or 1 stops at 0 or 1; m 0: no allele ever changes; steps go up a little more often than down (0.50 - 0.54, a known Nemo bias); same failures as disp_mutation_rate; e.g. 0.05 |

### founder initialisation

| Parameter | (1) Effect | (2) Req. | (3) Species | Note |
|---|---|---|---|---|
| `disp_init_rate` | ON | OPT | VAR | every founder allele of both traits = this value; outranks disp_init_rate_fem / _mal (no warning), loses to disp_init_distribution; outside [0, 1] Python is killed; e.g. 0.5 |
| `disp_init_rate_fem` | ON | OPT | VAR | every founder fdisp allele = this value; only together with disp_init_rate_mal: one alone kills Python at boot; e.g. 0.2 |
| `disp_init_rate_mal` | ON | OPT | VAR | every founder mdisp allele = this value; only together with disp_init_rate_fem; e.g. 0.8 |
| `disp_init_distribution` | ON | OPT | VAR | random founders, 'normal' \| 'uniform'; outranks every disp_init_rate\* row (Nemo warns 'discarding "disp_init_rate[...]" parameters'); no init row at all = 'uniform'. 'uniform': each allele U(0, 1). 'normal': each allele mean + N(0, sd), drawn again until it falls in [0, 1], so near 0 or 1 the founders' mean moves inward ({0.05, 0.3} gives 0.25); a mean outside [0, 1] (e.g. 1.5, sd 0.01) hangs make_founder() forever; another name ('gamma') kills Python at the first make_founder(); e.g. 'normal' |
| `disp_init_dist_params` | ON | OPT | VAR | {{mean, sd}}, one row per trait; required next to disp_init_distribution (missing: Python is killed at boot), even for 'uniform', which never reads it; 2 rows: **row 0 -> mdisp, row 1 -> fdisp** (Nemo's order: male = 0, female = 1); any other number of rows: row 0 for both; one column under 'normal' kills Python at the first make_founder(); without disp_init_distribution it is ignored; e.g. '{{0.7, 0.1}{0.3, 0.1}}' |

## (3) dispersal: Nemo's dispersal rows in MetaIBM

MetaIBM's dispersal functions read fdisp / mdisp when `use_disp_trait=True`; every dispersal function that has `use_disp_trait` (among patches, within patches, global habitat network) takes the same keywords, and one call applies them to every species in its pools.
Example numbers below: founders with fdisp 0.2 and mdisp 0.8, `method='uniform'`, two patches.

| Nemo row | MetaIBM | (1) Effect | (2) Req. | (3) Species | Note |
|---|---|---|---|---|---|
| `disperse_evoldisp` (Nemo's dispersal that reads fdisp / mdisp) | `use_disp_trait=True, is_rescale=False, sex_specific=True` | ON | OPT | FIX | an offspring goes from patch i to patch j with probability v × m_ij, v = fdisp for a female, mdisp for a male, m = total_disp_among_rate; with m = 1 it leaves with probability v, as in Nemo: about 20 % of the females and 80 % of the males leave; sex_specific=False: males read fdisp too; is_rescale=True (default): as many leave as without disp traits, v only decides who; every species in the call needs both traits, otherwise KeyError('fdisp') |
| `dispersal_cost` | `disp_cost=c` | ON | OPT | FIX | each leaver dies on the way with probability c (c 0.5: half of the leavers die) and goes back to Nemo; None (default) = 0; ignored when is_rescale=True or use_disp_trait=False |
| `dispersal_cost_fem` + `dispersal_cost_mal` | `disp_cost=(c_fem, c_mal)` | ON | OPT | FIX | c_fem for a female leaver, c_mal for a male one: (0.5, 0) kills half of the female leavers and no male |
| `dispersal_rate_fem` + `dispersal_rate_mal` (a fixed rate per sex) | `disp_init_rate_fem` + `disp_init_rate_mal` + `disp_mutation_rate 0`, then the first row with m = 1 | ON | OPT | VAR | no new keyword needed: with μ 0 the alleles never change, so a female always leaves with probability disp_init_rate_fem and a male with disp_init_rate_mal |
| `dispersal_fixed_trait` + `dispersal_fixed_rate` | `metaibm_nemo_params = {dispersal_fixed_trait, dispersal_fixed_rate}` + `use_disp_trait = True` | ON | OPT | VAR | same names, per species: `indfactory(..., metaibm_nemo_params={'dispersal_fixed_trait': 'female', 'dispersal_fixed_rate': 0.6})`; with use_disp_trait=True the named sex leaves with the fixed rate instead of its own trait, the other sex keeps reading (and evolving) its own: at m = 1 about 60 % of the females and 80 % of the males leave; the fixed sex pays disp_cost too; unset (default): no change; a trait other than 'female' / 'male', or a rate outside [0, 1]: ValueError; set both together (a trait alone: KeyError at the first dispersal; rate None: the fixed sex never leaves) |

---

# MetaIBM-Nemo disp 参数双语手册

参数怎么传进去：下面每一行 Nemo 参数都经 `pynemo.params.SimParameters` 传入——每行调用一次 `nemo_params.add(key, value)`——再把 `nemo_params` 交给
`metaibm.indfactory.indfactory(species_id, pheno_name_ls, nemo_params, trait_keys=('fdisp', 'mdisp'), ...)`：一个物种 = 一个 indfactory = 一个 Nemo 启动模拟。从没 add 过的行就是未设置。
不要用 `None` 去 add：`SimParameters.add(k, None)` 会写出一个没有值的裸键，Nemo 会把它读成 `k 1`（例如 `add('disp_init_rate', None)` 让每个 founder 等位基因都是 1.0）——不想设就直接不 add。
这里不列 Nemo 的扩散生命周期参数（`disperse`、`disperse_evoldisp`、`dispersal_*` 等）：MetaIBM 从不运行它们，都不要设置；第 (3) 节说明其中哪些由 MetaIBM 自己的扩散提供。
举例：example1/use-disp-traits.ipynb 把参数行直接写在代码里：`non_genetic` + `filename` + `mutation` + 每个物种一个 `init` 字典（`disp_params()`）。

## 标签

- (1) Nemo 对 MetaIBM 的作用：ON / OFF
- (2) MetaIBM-NEMO 是否需要：REQ / OPT
- (3) 各物种之间：FIX / VAR（静态属性？）
- 注：`*` = 有条件：所示标签是安全的默认值，该行写明何时可以放宽
- 注：`!` = 必须（VAR! = 各物种必须不同）

| 标签组合（disp 表） | 行数 | 含义 |
|---|---|---|
| ON REQ VAR | 2 | Nemo 需要它来构建 disp 性状，它会影响生成的个体，每个物种可以各自设置 |
| ON OPT VAR | 5 | 可选的物种级旋钮：founder 怎么抽取 |

## (1) non_genetic

用来构建启动模拟（boot sim），这里没有任何东西推进它——这些行取什么值，生成的个体都一样，该行另有说明的除外。

| 参数 | (1) 作用 | (2) 需要 | (3) 物种间 | 说明 |
|---|---|---|---|---|
| `replicates` | OFF | REQ | FIX | 只用于启动模拟；没有它 Nemo 会等待 10 秒，不构建 disp 性状，indfactory 抛出 ValueError；例：1 |
| `generations` | OFF | REQ | FIX | 只用于启动模拟；缺失时的失败方式同 replicates；例：1 |
| `patch_number` | OFF | REQ | FIX | 只用于启动模拟；没有它会报 'param "patch_number" is missing!'，indfactory 抛出 ValueError；例：1 |
| `patch_capacity` | OFF\* | REQ | FIX | 保持 1：取 2 及以上时，若同时设了 `aging` 且 founder 是随机的（disp_init_distribution，或不设任何初始化行），同一个种子不再得到同样的个体（分布不变）；没有它 Python 会被杀掉（'population parameters are not properly set!'）；例：1 |
| `aging` | OFF | REQ | FIX | 只用于启动模拟；用来避免 Nemo 的警告 >>>> WARNING::Metapop::setPopulation: generation 0 is empty!! |
| `run_mode` | OFF | OPT | FIX | 只影响控制台输出和文件；'silent'：不打印 'outputs:' 行，.log 文件照常写出；'create_init' 还会在工作目录写出 `<filename>.ini`；未知的值：indfactory 抛出 ValueError；例：想看到 Nemo 的警告就用 'overwrite' |
| `root_dir` | OFF | OPT | FIX | .log 文件放在哪里；该目录确实会在启动时创建；不设：工作目录；例：'nemo_cache' |
| `logfile` | OFF | OPT | FIX | Nemo 的运行日志；不设：工作目录下的 nemo.log；它所在的目录必须事先存在，否则报 '\*\*\*ERROR\*\*\* could not create simulation logfile'，不写日志（启动继续）；例：'nemo_cache/nemo.log' |
| `filename` | OFF | REQ | VAR! | 为该物种的 .log 命名——唯一一个按物种设置的参数；没有它的失败方式同 replicates；共用同一个名字不会报错，只是各物种会互相覆盖日志；例：'sp1' |

## (2) disp

遗传部分，也就是 MetaIBM 从 prototype 上读回的内容（trait_keys=('fdisp', 'mdisp')）。
每个个体都带两个性状，与它自己的性别无关：fdisp（雌性的扩散率）和 mdisp（雄性的扩散率）。每个性状一个位点、两个 [0, 1] 内的等位基因（母源、父源）；性状值 = 两者的平均。
这些参数总是同时构建两个性状，所以 trait_keys 要两个都写（只写一个会抛出 ValueError）。每个物种可以各自设置：同一次运行里的两个物种各自保持自己的 founder 和突变。

### 突变（mutation）

| 参数 | (1) 作用 | (2) 需要 | (3) 物种间 | 说明 |
|---|---|---|---|---|
| `disp_mutation_rate` | ON | REQ | VAR | μ，每个等位基因拷贝、每次出生（有性后代或克隆；founder 和 make_copy() 从不突变）：新生个体的一个性状改变的概率为 1 - e^(-2μ)，μ 0.02 -> 约 4 %；两个性状共用这一个突变率；没有它（或没有 disp_mutation_mean）Nemo 每个性状等待 10 秒，不构建 disp 性状，indfactory 抛出 ValueError；超出 [0, 1] 时 Python 被杀掉（'argument is out of bounds'）；例：0.02 |
| `disp_mutation_mean` | ON | REQ | VAR | m：一次突变把一个等位基因向上或向下移动一步，步长 ~ Exponential(均值 m)，所以平均步长就是 m；越过 0 或 1 的一步停在 0 或 1；m 0：等位基因永远不变；向上的步略多于向下（0.50 - 0.54，Nemo 已知的偏差）；缺失和越界的失败方式同 disp_mutation_rate；例：0.05 |

### 建群个体初始化（founder initialisation）

| 参数 | (1) 作用 | (2) 需要 | (3) 物种间 | 说明 |
|---|---|---|---|---|
| `disp_init_rate` | ON | OPT | VAR | 两个性状的每个 founder 等位基因都等于这个值；优先于 disp_init_rate_fem / _mal（无警告），但让位于 disp_init_distribution；超出 [0, 1] 时 Python 被杀掉；例：0.5 |
| `disp_init_rate_fem` | ON | OPT | VAR | 每个 founder 的 fdisp 等位基因都等于这个值；只能与 disp_init_rate_mal 一起设：只设一个会在启动时杀掉 Python；例：0.2 |
| `disp_init_rate_mal` | ON | OPT | VAR | 每个 founder 的 mdisp 等位基因都等于这个值；只能与 disp_init_rate_fem 一起设；例：0.8 |
| `disp_init_distribution` | ON | OPT | VAR | 随机 founder，'normal' \| 'uniform'；优先于所有 disp_init_rate\* 行（Nemo 警告 'discarding "disp_init_rate[...]" parameters'）；一个初始化行都不设 = 'uniform'。'uniform'：每个等位基因 U(0, 1)。'normal'：每个等位基因 = mean + N(0, sd)，落在 [0, 1] 之外就重抽，所以靠近 0 或 1 时 founder 的均值会向内偏（{0.05, 0.3} 得到 0.25）；均值在 [0, 1] 之外（如 1.5，sd 0.01）会让 make_founder() 永远卡住；其他名字（'gamma'）会在第一次 make_founder() 时杀掉 Python；例：'normal' |
| `disp_init_dist_params` | ON | OPT | VAR | {{mean, sd}}，每个性状一行；设了 disp_init_distribution 就必须设（缺失：启动时 Python 被杀掉），'uniform' 也一样，尽管它从不读取；两行：**第 0 行 -> mdisp，第 1 行 -> fdisp**（Nemo 的顺序：雄 = 0，雌 = 1）；其他行数：两个性状都用第 0 行；'normal' 下只有一列会在第一次 make_founder() 时杀掉 Python；没有 disp_init_distribution 时它被忽略；例：'{{0.7, 0.1}{0.3, 0.1}}' |

## (3) 扩散：Nemo 的扩散参数在 MetaIBM 中的对应（dispersal: Nemo's dispersal rows in MetaIBM）

`use_disp_trait=True` 时，MetaIBM 的扩散函数读取 fdisp / mdisp；所有带 `use_disp_trait` 的扩散函数（patch 之间、patch 内部、全局 habitat 网络）都接受同样的关键字，一次调用对其 pool 里的所有物种生效。
下面的示例数字：founder 的 fdisp 0.2、mdisp 0.8，`method='uniform'`，两个 patch。

| Nemo 参数 | MetaIBM | (1) 作用 | (2) 需要 | (3) 物种间 | 说明 |
|---|---|---|---|---|---|
| `disperse_evoldisp`（Nemo 读取 fdisp / mdisp 的扩散） | `use_disp_trait=True, is_rescale=False, sex_specific=True` | ON | OPT | FIX | 后代从 patch i 去 patch j 的概率为 v × m_ij，雌性 v = fdisp，雄性 v = mdisp，m = total_disp_among_rate；m = 1 时后代以概率 v 离开，与 Nemo 相同：约 20 % 的雌性和 80 % 的雄性离开；sex_specific=False：雄性也读 fdisp；is_rescale=True（默认）：离开的数量与不用 disp 性状时相同，v 只决定谁走；调用里的每个物种都必须带两个性状，否则 KeyError('fdisp') |
| `dispersal_cost` | `disp_cost=c` | ON | OPT | FIX | 每个离开者以概率 c 死在路上（c 0.5：一半离开者死亡），并回到 Nemo；None（默认）= 0；is_rescale=True 或 use_disp_trait=False 时被忽略 |
| `dispersal_cost_fem` + `dispersal_cost_mal` | `disp_cost=(c_fem, c_mal)` | ON | OPT | FIX | 雌性离开者用 c_fem，雄性用 c_mal：(0.5, 0) 让一半雌性离开者死亡，雄性不死 |
| `dispersal_rate_fem` + `dispersal_rate_mal`（每个性别一个固定扩散率） | `disp_init_rate_fem` + `disp_init_rate_mal` + `disp_mutation_rate 0`，再用第一行，m = 1 | ON | OPT | VAR | 不需要新关键字：μ 0 时等位基因永远不变，所以雌性总是以 disp_init_rate_fem 的概率离开，雄性以 disp_init_rate_mal |
| `dispersal_fixed_trait` + `dispersal_fixed_rate` | `metaibm_nemo_params = {dispersal_fixed_trait, dispersal_fixed_rate}` + `use_disp_trait = True` | ON | OPT | VAR | 同名，按物种设置：`indfactory(..., metaibm_nemo_params={'dispersal_fixed_trait': 'female', 'dispersal_fixed_rate': 0.6})`；use_disp_trait=True 时，指定的性别以固定扩散率离开，不看自己的性状，另一个性别照旧读自己的性状并演化：m = 1 时约 60 % 的雌性和 80 % 的雄性离开；被固定的性别同样承担 disp_cost；不设（默认）：没有任何变化；trait 不是 'female' / 'male'，或 rate 超出 [0, 1]：ValueError；两个要一起设（只设 trait：第一次扩散时 KeyError；rate 为 None：被固定的性别永远不离开） |
