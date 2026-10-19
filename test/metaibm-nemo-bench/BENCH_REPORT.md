# Benchmark: MetaIBM-Nemo vs pure MetaIBM (sexual)

Generated 2026-10-19 02:41 by `test/metaibm-nemo-run/bench.py` on `macOS-15.7.3-arm64-arm-64bit`, python `3.11.17`.

- nemo = `/Users/jianhaolin/Library/CloudStorage/OneDrive-UniversityofHelsinki/MetaIBM4Nemo_workspace/MetaIBM4Nemo-features-developing` (NemoBackend); pure = `/Users/jianhaolin/Library/CloudStorage/OneDrive-UniversityofHelsinki/MetaIBM_workspaces/Macbook_Air/MetaIBM-bug-fixing` (pure-Python MetaIBM, read only)
- same landscape / parameters / seed / call sequence on both sides; `is_remove=True` passed wherever the function takes it (same semantics: things move); one case at a time on an idle machine (load at the end of each case is listed)
- genome size is identical on both sides: Nemo quanti_loci = loci per trait x traits (asserted per case via get_num_loci; an earlier run that passed loci per trait as the Nemo total was discarded)
- speed-up = pure / nemo (> 1: MetaIBM-Nemo faster). Population sizes evolve stochastically and differ between the two runs: compare also "s per 1000 island individuals per step".

## Interfaces (experiment setups: islands start empty, marker schedule, 100 steps, CSV every 10 steps)

| case | pure s/step | nemo s/step | speed-up | islands alive mean pure / nemo | per 1000 individuals pure / nemo | setup pure / nemo s | RSS end pure / nemo MB |
|---|---|---|---|---|---|---|---|
| model | 1.775 | 1.122 | **1.58x** | 33757 / 31144 | 0.0526 / 0.0360 (1.46x) | 0.6 / 0.6 | 742 / 1232 |
| sloss-p16 | 1.349 | 1.474 | **0.92x** | 14819 / 28843 | 0.0910 / 0.0511 (1.78x) | 0.1 / 0.1 | 778 / 1284 |
| sloss-p64-hetero | 1.033 | 0.497 | **2.08x** | 16658 / 13144 | 0.0620 / 0.0378 (1.64x) | 1.1 / 1.0 | 793 / 1079 |
| sloss-p1 | 4.312 | 1.691 | **2.55x** | 20143 / 9498 | 0.2141 / 0.1780 (1.20x) | 0.1 / 0.1 | 762 / 1058 |
| grfe-p16 | 0.733 | 0.545 | **1.34x** | 4892 / 3431 | 0.1498 / 0.1588 (0.94x) | 3.8 / 3.8 | 773 / 880 |

## Boundaries (islands fully occupied from step 0, 102400 sites in 16 patches, 30 steps): habitat size x schedule

| habitat | schedule | pure s/step | nemo s/step | speed-up | init 102400 + 6400 founders pure / nemo s | RSS end pure / nemo MB |
|---|---|---|---|---|---|---|
| 10x10 | marker | 2.314 | 1.183 | **1.96x** | 2.1 / 0.5 | 685 / 1031 |
| 10x10 | pool | 5.064 | 1.149 | **4.41x** | 2.1 / 0.5 | 744 / 1238 |
| 20x20 | marker | 2.609 | 1.437 | **1.82x** | 2.2 / 0.5 | 650 / 994 |
| 20x20 | pool | 5.056 | 0.949 | **5.33x** | 2.2 / 0.5 | 686 / 1190 |
| 40x40 | marker | 4.506 | 3.567 | **1.26x** | 2.2 / 0.5 | 647 / 984 |
| 40x40 | pool | 5.279 | 1.087 | **4.86x** | 2.2 / 0.5 | 679 / 1191 |
| 80x80 | marker | 13.055 | 12.890 | **1.01x** | 2.2 / 0.6 | 591 / 940 |
| 80x80 | pool | 6.919 | 1.703 | **4.06x** | 2.2 / 0.6 | 627 / 1152 |

## Parameter sweeps (islands 16 x 16 habitats 20x20 fully occupied, 20 steps; one parameter changed from {'dead': 0.1, 'wid': 0.5, 'birth': 1.0, 'rain': 6.4, 'm_within': 0.1, 'm_among': 0.1, 'dist': 1e-05, 'loci': 20, 'traits': 2, 'mut': 0.0001})

| variant | change | schedule | pure s/step | nemo s/step | speed-up | init pure / nemo s | RSS end pure / nemo MB |
|---|---|---|---|---|---|---|---|
| base | baseline | marker | 2.622 | 1.430 | **1.83x** | 2.2 / 0.5 | 648 / 995 |
| base | baseline | pool | 5.700 | 0.932 | **6.11x** | 2.2 / 0.5 | 687 / 1188 |
| dead0.05 | dead=0.05 | marker | 1.548 | 1.006 | **1.54x** | 2.2 / 0.5 | 642 / 994 |
| dead0.05 | dead=0.05 | pool | 5.731 | 0.829 | **6.92x** | 2.2 / 0.5 | 687 / 1230 |
| dead0.3 | dead=0.3 | marker | 3.939 | 1.322 | **2.98x** | 2.2 / 0.5 | 642 / 988 |
| dead0.3 | dead=0.3 | pool | 3.445 | 0.698 | **4.93x** | 2.2 / 0.5 | 641 / 998 |
| wid0.2 | wid=0.2 | marker | 3.355 | 1.924 | **1.74x** | 2.2 / 0.5 | 654 / 993 |
| wid0.2 | wid=0.2 | pool | 5.659 | 0.953 | **5.94x** | 2.2 / 0.5 | 676 / 1136 |
| wid2.0 | wid=2.0 | marker | 2.185 | 1.170 | **1.87x** | 2.2 / 0.5 | 647 / 996 |
| wid2.0 | wid=2.0 | pool | 5.484 | 0.860 | **6.38x** | 2.2 / 0.5 | 686 / 1212 |
| birth0.5 | birth=0.5 | marker | 2.545 | 1.403 | **1.81x** | 2.2 / 0.5 | 644 / 992 |
| birth0.5 | birth=0.5 | pool | 2.966 | 0.769 | **3.86x** | 2.2 / 0.5 | 668 / 1083 |
| rain64 | rain=64.0 | marker | 2.634 | 1.419 | **1.86x** | 2.2 / 0.5 | 648 / 993 |
| rain64 | rain=64.0 | pool | 5.595 | 0.911 | **6.14x** | 2.2 / 0.5 | 682 / 1186 |
| disp0.5 | m_within=0.5, m_among=0.5 | marker | 4.295 | 2.043 | **2.10x** | 2.2 / 0.5 | 647 / 1031 |
| disp0.5 | m_within=0.5, m_among=0.5 | pool | 4.884 | 0.933 | **5.23x** | 2.2 / 0.5 | 674 / 1167 |
| dist0.01 | dist=0.01 | marker | 2.464 | 1.325 | **1.86x** | 2.2 / 0.5 | 642 / 990 |
| dist0.01 | dist=0.01 | pool | 5.562 | 0.867 | **6.42x** | 2.2 / 0.5 | 677 / 1187 |
| loci5 | loci=5 | marker | 1.757 | 1.403 | **1.25x** | 1.1 / 0.5 | 604 / 948 |
| loci5 | loci=5 | pool | 2.525 | 0.907 | **2.78x** | 1.1 / 0.5 | 641 / 1125 |
| loci80 | loci=80 | marker | 5.617 | 1.436 | **3.91x** | 10.2 / 0.7 | 847 / 1247 |
| loci80 | loci=80 | pool | 15.816 | 0.960 | **16.48x** | 10.2 / 0.7 | 884 / 1536 |
| loci320 | loci=320 | marker | 18.181 | 1.500 | **12.12x** | 103.0 / 1.5 | 1685 / 2003 |
| loci320 | loci=320 | pool | 54.565 | 1.097 | **49.75x** | 102.2 / 1.4 | 1719 / 2552 |
| traits1 | traits=1 | marker | 1.683 | 1.187 | **1.42x** | 1.2 / 0.5 | 581 / 989 |
| traits1 | traits=1 | pool | 3.061 | 0.834 | **3.67x** | 1.2 / 0.5 | 609 / 1148 |
| traits4 | traits=4 | marker | 4.109 | 1.706 | **2.41x** | 4.1 / 0.6 | 833 / 1074 |
| traits4 | traits=4 | pool | 8.969 | 1.083 | **8.29x** | 4.1 / 0.6 | 829 / 1277 |
| mut0.01 | mut=0.01 | marker | 2.633 | 1.524 | **1.73x** | 2.2 / 0.5 | 648 / 993 |
| mut0.01 | mut=0.01 | pool | 5.109 | 0.939 | **5.44x** | 2.2 / 0.5 | 681 / 1183 |

## Per phase (seconds per step)

<details><summary>model</summary>

| phase | pure | nemo | speed-up |
|---|---|---|---|
| germinate + birth from markers | 1.0637 | 0.4813 | 2.21x |
| mainland dead selection | 0.2151 | 0.2509 | 0.86x |
| islands dead selection | 0.1516 | 0.1628 | 0.93x |
| mainland sexual birth | 0.1312 | 0.0171 | 7.68x |
| marker dispersal among | 0.0690 | 0.0684 | 1.01x |
| colonization | 0.0649 | 0.0616 | 1.05x |
| data: phenotype csv.gz | 0.0397 | 0.0413 | 0.96x |
| data: sp_id csv.gz | 0.0335 | 0.0317 | 1.06x |
| marker calculation | 0.0055 | 0.0055 | 0.99x |
| marker clear-up | 0.0002 | 0.0003 | 0.61x |
| disturbance | 0.0001 | 0.0001 | 1.00x |
| environment offset | 0.0000 | 0.0000 | 1.00x |

- load at end: pure (3.5, 3.7, 3.2), nemo (3.3, 3.6, 3.2)
</details>

<details><summary>sloss-p16</summary>

| phase | pure | nemo | speed-up |
|---|---|---|---|
| germinate + birth from markers | 0.9239 | 0.9888 | 0.93x |
| colonization | 0.1247 | 0.1084 | 1.15x |
| islands dead selection | 0.1094 | 0.2120 | 0.52x |
| mainland sexual birth | 0.0741 | 0.0049 | 15.24x |
| data: phenotype csv.gz | 0.0542 | 0.0809 | 0.67x |
| data: sp_id csv.gz | 0.0365 | 0.0376 | 0.97x |
| mainland dead selection | 0.0183 | 0.0287 | 0.64x |
| marker calculation | 0.0027 | 0.0047 | 0.59x |
| marker dispersal among | 0.0026 | 0.0044 | 0.60x |
| marker dispersal within | 0.0017 | 0.0025 | 0.70x |
| marker clear-up | 0.0000 | 0.0002 | 0.24x |
| disturbance | 0.0000 | 0.0000 | 0.67x |

- load at end: pure (3.1, 3.4, 3.2), nemo (3.3, 3.5, 3.3)
</details>

<details><summary>sloss-p64-hetero</summary>

| phase | pure | nemo | speed-up |
|---|---|---|---|
| germinate + birth from markers | 0.6024 | 0.1318 | 4.57x |
| colonization | 0.1255 | 0.1200 | 1.05x |
| islands dead selection | 0.0913 | 0.1002 | 0.91x |
| mainland sexual birth | 0.0746 | 0.0046 | 16.08x |
| data: phenotype csv.gz | 0.0566 | 0.0543 | 1.04x |
| data: sp_id csv.gz | 0.0322 | 0.0275 | 1.17x |
| marker dispersal among | 0.0230 | 0.0215 | 1.07x |
| mainland dead selection | 0.0184 | 0.0286 | 0.64x |
| marker dispersal within | 0.0050 | 0.0047 | 1.05x |
| marker calculation | 0.0029 | 0.0024 | 1.20x |
| marker clear-up | 0.0001 | 0.0001 | 0.67x |
| disturbance | 0.0001 | 0.0001 | 1.00x |

- load at end: pure (3.3, 3.3, 3.2), nemo (3.3, 3.4, 3.2)
</details>

<details><summary>sloss-p1</summary>

| phase | pure | nemo | speed-up |
|---|---|---|---|
| germinate + birth from markers | 3.6423 | 1.2192 | 2.99x |
| islands dead selection | 0.2465 | 0.1278 | 1.93x |
| data: sp_id csv.gz | 0.1356 | 0.1244 | 1.09x |
| colonization | 0.1232 | 0.1313 | 0.94x |
| mainland sexual birth | 0.0768 | 0.0058 | 13.16x |
| data: phenotype csv.gz | 0.0614 | 0.0487 | 1.26x |
| mainland dead selection | 0.0186 | 0.0284 | 0.65x |
| marker dispersal within | 0.0034 | 0.0024 | 1.40x |
| marker calculation | 0.0029 | 0.0015 | 1.97x |
| marker clear-up | 0.0000 | 0.0000 | 1.33x |
| disturbance | 0.0000 | 0.0000 | 0.50x |
| marker dispersal among | 0.0000 | 0.0000 | - |

- load at end: pure (3.3, 3.4, 3.3), nemo (4.3, 3.9, 3.5)
</details>

<details><summary>grfe-p16</summary>

| phase | pure | nemo | speed-up |
|---|---|---|---|
| marker dispersal within | 0.2113 | 0.2100 | 1.01x |
| germinate + birth from markers | 0.1584 | 0.0330 | 4.80x |
| colonization | 0.1332 | 0.1332 | 1.00x |
| mainland sexual birth | 0.0748 | 0.0047 | 15.85x |
| islands dead selection | 0.0607 | 0.0625 | 0.97x |
| data: phenotype csv.gz | 0.0379 | 0.0359 | 1.06x |
| data: sp_id csv.gz | 0.0317 | 0.0302 | 1.05x |
| mainland dead selection | 0.0183 | 0.0287 | 0.64x |
| marker dispersal among | 0.0027 | 0.0026 | 1.04x |
| marker calculation | 0.0025 | 0.0027 | 0.94x |
| marker clear-up | 0.0001 | 0.0001 | 0.64x |
| disturbance | 0.0000 | 0.0000 | 0.67x |

- load at end: pure (3.9, 3.9, 3.6), nemo (3.2, 3.7, 3.5)
</details>

<details><summary>L10-marker</summary>

| phase | pure | nemo | speed-up |
|---|---|---|---|
| germinate + birth from markers | 1.5752 | 0.3414 | 4.61x |
| islands dead selection | 0.3268 | 0.4725 | 0.69x |
| marker dispersal within | 0.2547 | 0.2515 | 1.01x |
| mainland sexual birth | 0.0721 | 0.0097 | 7.41x |
| colonization | 0.0262 | 0.0361 | 0.72x |
| marker calculation | 0.0240 | 0.0238 | 1.01x |
| mainland dead selection | 0.0183 | 0.0300 | 0.61x |
| marker dispersal among | 0.0147 | 0.0142 | 1.03x |
| marker clear-up | 0.0012 | 0.0012 | 0.97x |
| disturbance | 0.0000 | 0.0000 | 1.00x |

- load at end: pure (3.8, 3.8, 3.6), nemo (3.6, 3.8, 3.6)
</details>

<details><summary>L10-pool</summary>

| phase | pure | nemo | speed-up |
|---|---|---|---|
| reproduce into offspring pool | 4.3133 | 0.2159 | 19.98x |
| islands dead selection | 0.3215 | 0.5470 | 0.59x |
| dispersal within -> immigrant pool | 0.2541 | 0.2525 | 1.01x |
| mainland sexual birth | 0.0728 | 0.0056 | 13.07x |
| colonization | 0.0255 | 0.0267 | 0.96x |
| dormancy | 0.0220 | 0.0255 | 0.86x |
| germinate O+I+D | 0.0197 | 0.0293 | 0.67x |
| mainland dead selection | 0.0185 | 0.0299 | 0.62x |
| dispersal among -> immigrant pool | 0.0152 | 0.0146 | 1.04x |
| clear-up O+I | 0.0001 | 0.0002 | 0.67x |
| disturbance | 0.0000 | 0.0000 | 1.00x |

- load at end: pure (3.1, 3.6, 3.5), nemo (3.1, 3.5, 3.5)
</details>

<details><summary>L20-marker</summary>

| phase | pure | nemo | speed-up |
|---|---|---|---|
| germinate + birth from markers | 2.0182 | 0.8114 | 2.49x |
| islands dead selection | 0.3986 | 0.5012 | 0.80x |
| mainland sexual birth | 0.0803 | 0.0074 | 10.80x |
| marker dispersal within | 0.0273 | 0.0264 | 1.03x |
| marker calculation | 0.0250 | 0.0180 | 1.38x |
| colonization | 0.0249 | 0.0262 | 0.95x |
| mainland dead selection | 0.0183 | 0.0295 | 0.62x |
| marker dispersal among | 0.0143 | 0.0139 | 1.03x |
| marker clear-up | 0.0009 | 0.0009 | 1.00x |
| disturbance | 0.0000 | 0.0000 | 1.00x |

- load at end: pure (3.1, 3.5, 3.5), nemo (3.0, 3.4, 3.5)
</details>

<details><summary>L20-pool</summary>

| phase | pure | nemo | speed-up |
|---|---|---|---|
| reproduce into offspring pool | 4.5141 | 0.2761 | 16.35x |
| islands dead selection | 0.3356 | 0.4988 | 0.67x |
| mainland sexual birth | 0.0718 | 0.0054 | 13.39x |
| dispersal within -> immigrant pool | 0.0299 | 0.0293 | 1.02x |
| colonization | 0.0243 | 0.0340 | 0.71x |
| dormancy | 0.0223 | 0.0257 | 0.87x |
| germinate O+I+D | 0.0218 | 0.0319 | 0.68x |
| mainland dead selection | 0.0184 | 0.0299 | 0.61x |
| dispersal among -> immigrant pool | 0.0165 | 0.0159 | 1.04x |
| disturbance | 0.0000 | 0.0000 | 1.00x |
| clear-up O+I | 0.0000 | 0.0000 | 1.00x |

- load at end: pure (3.2, 3.5, 3.5), nemo (4.4, 3.8, 3.6)
</details>

<details><summary>L40-marker</summary>

| phase | pure | nemo | speed-up |
|---|---|---|---|
| germinate + birth from markers | 3.8900 | 2.8348 | 1.37x |
| islands dead selection | 0.4527 | 0.6295 | 0.72x |
| mainland sexual birth | 0.0804 | 0.0051 | 15.86x |
| colonization | 0.0246 | 0.0274 | 0.90x |
| marker calculation | 0.0185 | 0.0182 | 1.01x |
| mainland dead selection | 0.0183 | 0.0296 | 0.62x |
| marker dispersal among | 0.0143 | 0.0140 | 1.02x |
| marker dispersal within | 0.0057 | 0.0058 | 0.99x |
| marker clear-up | 0.0008 | 0.0008 | 0.92x |
| disturbance | 0.0000 | 0.0000 | 1.00x |

- load at end: pure (3.4, 3.6, 3.6), nemo (3.4, 3.5, 3.5)
</details>

<details><summary>L40-pool</summary>

| phase | pure | nemo | speed-up |
|---|---|---|---|
| reproduce into offspring pool | 4.6086 | 0.2328 | 19.80x |
| islands dead selection | 0.4562 | 0.6499 | 0.70x |
| mainland sexual birth | 0.0714 | 0.0053 | 13.56x |
| germinate O+I+D | 0.0366 | 0.0471 | 0.78x |
| colonization | 0.0239 | 0.0549 | 0.43x |
| dispersal among -> immigrant pool | 0.0236 | 0.0226 | 1.04x |
| dormancy | 0.0231 | 0.0259 | 0.89x |
| mainland dead selection | 0.0184 | 0.0301 | 0.61x |
| dispersal within -> immigrant pool | 0.0167 | 0.0162 | 1.03x |
| disturbance | 0.0000 | 0.0000 | 1.00x |
| clear-up O+I | 0.0000 | 0.0000 | - |

- load at end: pure (4.6, 3.8, 3.6), nemo (3.2, 3.5, 3.5)
</details>

<details><summary>L80-marker</summary>

| phase | pure | nemo | speed-up |
|---|---|---|---|
| germinate + birth from markers | 11.9728 | 11.6619 | 1.03x |
| islands dead selection | 0.8791 | 1.1040 | 0.80x |
| mainland sexual birth | 0.0777 | 0.0048 | 16.31x |
| marker dispersal among | 0.0614 | 0.0355 | 1.73x |
| colonization | 0.0233 | 0.0315 | 0.74x |
| marker calculation | 0.0205 | 0.0193 | 1.06x |
| mainland dead selection | 0.0181 | 0.0303 | 0.60x |
| marker clear-up | 0.0006 | 0.0006 | 1.06x |
| marker dispersal within | 0.0004 | 0.0004 | 0.92x |
| disturbance | 0.0000 | 0.0000 | 1.00x |

- load at end: pure (1.6, 2.0, 2.8), nemo (1.6, 1.6, 2.2)
</details>

<details><summary>L80-pool</summary>

| phase | pure | nemo | speed-up |
|---|---|---|---|
| reproduce into offspring pool | 5.6868 | 0.2756 | 20.63x |
| islands dead selection | 0.8683 | 1.0748 | 0.81x |
| germinate O+I+D | 0.1143 | 0.1589 | 0.72x |
| dispersal among -> immigrant pool | 0.0919 | 0.0960 | 0.96x |
| mainland sexual birth | 0.0810 | 0.0050 | 16.09x |
| colonization | 0.0340 | 0.0338 | 1.01x |
| dormancy | 0.0224 | 0.0261 | 0.86x |
| mainland dead selection | 0.0190 | 0.0298 | 0.64x |
| dispersal within -> immigrant pool | 0.0004 | 0.0004 | 1.00x |
| disturbance | 0.0000 | 0.0000 | 1.00x |
| clear-up O+I | 0.0000 | 0.0000 | - |

- load at end: pure (1.5, 1.6, 2.1), nemo (1.7, 1.7, 2.1)
</details>

<details><summary>sweep-marker-base</summary>

| phase | pure | nemo | speed-up |
|---|---|---|---|
| germinate + birth from markers | 2.0476 | 0.8049 | 2.54x |
| islands dead selection | 0.3695 | 0.4989 | 0.74x |
| mainland sexual birth | 0.0907 | 0.0083 | 10.87x |
| marker dispersal within | 0.0275 | 0.0269 | 1.02x |
| marker calculation | 0.0272 | 0.0183 | 1.48x |
| colonization | 0.0245 | 0.0264 | 0.93x |
| mainland dead selection | 0.0182 | 0.0297 | 0.61x |
| marker dispersal among | 0.0142 | 0.0138 | 1.03x |
| marker clear-up | 0.0009 | 0.0009 | 0.94x |
| disturbance | 0.0000 | 0.0001 | 0.00x |

- load at end: pure (1.4, 1.6, 2.1), nemo (1.5, 1.6, 2.0)
</details>

<details><summary>sweep-pool-base</summary>

| phase | pure | nemo | speed-up |
|---|---|---|---|
| reproduce into offspring pool | 5.1496 | 0.2626 | 19.61x |
| islands dead selection | 0.3352 | 0.4984 | 0.67x |
| mainland sexual birth | 0.0805 | 0.0052 | 15.34x |
| dispersal within -> immigrant pool | 0.0303 | 0.0298 | 1.02x |
| colonization | 0.0249 | 0.0309 | 0.81x |
| dormancy | 0.0220 | 0.0253 | 0.87x |
| germinate O+I+D | 0.0215 | 0.0316 | 0.68x |
| mainland dead selection | 0.0182 | 0.0299 | 0.61x |
| dispersal among -> immigrant pool | 0.0163 | 0.0158 | 1.03x |
| disturbance | 0.0001 | 0.0001 | 1.00x |
| clear-up O+I | 0.0001 | 0.0001 | 1.00x |

- load at end: pure (1.6, 1.6, 2.0), nemo (1.5, 1.6, 2.0)
</details>

<details><summary>sweep-marker-dead0.05</summary>

| phase | pure | nemo | speed-up |
|---|---|---|---|
| germinate + birth from markers | 1.1159 | 0.4804 | 2.32x |
| islands dead selection | 0.2752 | 0.4138 | 0.66x |
| mainland sexual birth | 0.0437 | 0.0036 | 11.99x |
| marker calculation | 0.0345 | 0.0190 | 1.82x |
| marker dispersal within | 0.0284 | 0.0276 | 1.03x |
| colonization | 0.0173 | 0.0190 | 0.91x |
| mainland dead selection | 0.0158 | 0.0244 | 0.65x |
| marker dispersal among | 0.0154 | 0.0149 | 1.03x |
| marker clear-up | 0.0010 | 0.0011 | 0.95x |
| disturbance | 0.0000 | 0.0001 | 0.00x |

- load at end: pure (1.4, 1.5, 1.9), nemo (1.3, 1.5, 1.9)
</details>

<details><summary>sweep-pool-dead0.05</summary>

| phase | pure | nemo | speed-up |
|---|---|---|---|
| reproduce into offspring pool | 5.2912 | 0.2515 | 21.03x |
| islands dead selection | 0.2722 | 0.4101 | 0.66x |
| mainland sexual birth | 0.0437 | 0.0038 | 11.50x |
| dispersal within -> immigrant pool | 0.0314 | 0.0304 | 1.03x |
| dormancy | 0.0248 | 0.0307 | 0.81x |
| dispersal among -> immigrant pool | 0.0176 | 0.0171 | 1.03x |
| colonization | 0.0168 | 0.0288 | 0.58x |
| germinate O+I+D | 0.0161 | 0.0290 | 0.55x |
| mainland dead selection | 0.0157 | 0.0246 | 0.64x |
| disturbance | 0.0001 | 0.0001 | 1.00x |
| clear-up O+I | 0.0001 | 0.0001 | 1.00x |

- load at end: pure (1.3, 1.4, 1.8), nemo (1.4, 1.4, 1.8)
</details>

<details><summary>sweep-marker-dead0.3</summary>

| phase | pure | nemo | speed-up |
|---|---|---|---|
| germinate + birth from markers | 3.2647 | 0.7759 | 4.21x |
| islands dead selection | 0.3054 | 0.3745 | 0.82x |
| mainland sexual birth | 0.2271 | 0.0091 | 24.95x |
| colonization | 0.0785 | 0.0922 | 0.85x |
| mainland dead selection | 0.0219 | 0.0324 | 0.68x |
| marker dispersal within | 0.0209 | 0.0198 | 1.06x |
| marker calculation | 0.0112 | 0.0093 | 1.20x |
| marker dispersal among | 0.0080 | 0.0071 | 1.13x |
| marker clear-up | 0.0001 | 0.0001 | 1.00x |
| disturbance | 0.0000 | 0.0001 | 0.00x |

- load at end: pure (1.4, 1.4, 1.8), nemo (1.6, 1.5, 1.8)
</details>

<details><summary>sweep-pool-dead0.3</summary>

| phase | pure | nemo | speed-up |
|---|---|---|---|
| reproduce into offspring pool | 2.7414 | 0.1044 | 26.26x |
| islands dead selection | 0.3121 | 0.3943 | 0.79x |
| mainland sexual birth | 0.2124 | 0.0090 | 23.48x |
| colonization | 0.0944 | 0.0860 | 1.10x |
| germinate O+I+D | 0.0308 | 0.0413 | 0.75x |
| dispersal within -> immigrant pool | 0.0220 | 0.0209 | 1.05x |
| mainland dead selection | 0.0217 | 0.0321 | 0.68x |
| dispersal among -> immigrant pool | 0.0089 | 0.0080 | 1.13x |
| dormancy | 0.0002 | 0.0002 | 0.75x |
| disturbance | 0.0001 | 0.0001 | 1.00x |
| clear-up O+I | 0.0000 | 0.0001 | 0.00x |

- load at end: pure (2.0, 1.6, 1.8), nemo (1.8, 1.6, 1.8)
</details>

<details><summary>sweep-marker-wid0.2</summary>

| phase | pure | nemo | speed-up |
|---|---|---|---|
| germinate + birth from markers | 2.7941 | 1.2521 | 2.23x |
| islands dead selection | 0.3386 | 0.5395 | 0.63x |
| mainland sexual birth | 0.1136 | 0.0064 | 17.75x |
| colonization | 0.0304 | 0.0371 | 0.82x |
| marker dispersal within | 0.0265 | 0.0260 | 1.02x |
| mainland dead selection | 0.0192 | 0.0319 | 0.60x |
| marker calculation | 0.0176 | 0.0166 | 1.06x |
| marker dispersal among | 0.0135 | 0.0123 | 1.10x |
| marker clear-up | 0.0006 | 0.0005 | 1.18x |
| disturbance | 0.0000 | 0.0001 | 0.00x |

- load at end: pure (1.3, 1.5, 1.7), nemo (1.2, 1.5, 1.7)
</details>

<details><summary>sweep-pool-wid0.2</summary>

| phase | pure | nemo | speed-up |
|---|---|---|---|
| reproduce into offspring pool | 5.0656 | 0.2375 | 21.33x |
| islands dead selection | 0.3455 | 0.5329 | 0.65x |
| mainland sexual birth | 0.1083 | 0.0070 | 15.59x |
| colonization | 0.0298 | 0.0405 | 0.74x |
| dispersal within -> immigrant pool | 0.0294 | 0.0279 | 1.06x |
| germinate O+I+D | 0.0259 | 0.0431 | 0.60x |
| mainland dead selection | 0.0195 | 0.0316 | 0.61x |
| dormancy | 0.0187 | 0.0163 | 1.15x |
| dispersal among -> immigrant pool | 0.0155 | 0.0140 | 1.11x |
| disturbance | 0.0001 | 0.0001 | 1.00x |
| clear-up O+I | 0.0001 | 0.0001 | 1.00x |

- load at end: pure (1.5, 1.5, 1.7), nemo (1.8, 1.6, 1.7)
</details>

<details><summary>sweep-marker-wid2.0</summary>

| phase | pure | nemo | speed-up |
|---|---|---|---|
| germinate + birth from markers | 1.7062 | 0.6081 | 2.81x |
| islands dead selection | 0.3050 | 0.4486 | 0.68x |
| mainland sexual birth | 0.0733 | 0.0044 | 16.67x |
| marker dispersal within | 0.0273 | 0.0270 | 1.01x |
| colonization | 0.0219 | 0.0219 | 1.00x |
| marker calculation | 0.0186 | 0.0177 | 1.05x |
| mainland dead selection | 0.0171 | 0.0258 | 0.66x |
| marker dispersal among | 0.0141 | 0.0139 | 1.01x |
| marker clear-up | 0.0009 | 0.0009 | 1.00x |
| disturbance | 0.0000 | 0.0001 | 0.00x |

- load at end: pure (1.6, 1.6, 1.7), nemo (1.5, 1.5, 1.7)
</details>

<details><summary>sweep-pool-wid2.0</summary>

| phase | pure | nemo | speed-up |
|---|---|---|---|
| reproduce into offspring pool | 4.9800 | 0.2044 | 24.36x |
| islands dead selection | 0.3035 | 0.5001 | 0.61x |
| mainland sexual birth | 0.0728 | 0.0046 | 16.00x |
| dispersal within -> immigrant pool | 0.0299 | 0.0301 | 0.99x |
| dormancy | 0.0222 | 0.0273 | 0.81x |
| colonization | 0.0219 | 0.0224 | 0.98x |
| germinate O+I+D | 0.0197 | 0.0264 | 0.75x |
| mainland dead selection | 0.0172 | 0.0262 | 0.66x |
| dispersal among -> immigrant pool | 0.0162 | 0.0161 | 1.00x |
| disturbance | 0.0001 | 0.0001 | 1.00x |
| clear-up O+I | 0.0001 | 0.0001 | 1.00x |

- load at end: pure (1.5, 1.5, 1.7), nemo (1.4, 1.5, 1.6)
</details>

<details><summary>sweep-marker-birth0.5</summary>

| phase | pure | nemo | speed-up |
|---|---|---|---|
| germinate + birth from markers | 2.0503 | 0.7962 | 2.58x |
| islands dead selection | 0.3167 | 0.4999 | 0.63x |
| mainland sexual birth | 0.0896 | 0.0049 | 18.29x |
| colonization | 0.0245 | 0.0262 | 0.94x |
| marker dispersal within | 0.0204 | 0.0202 | 1.01x |
| mainland dead selection | 0.0180 | 0.0294 | 0.61x |
| marker calculation | 0.0166 | 0.0163 | 1.02x |
| marker dispersal among | 0.0077 | 0.0075 | 1.03x |
| marker clear-up | 0.0003 | 0.0003 | 1.17x |
| disturbance | 0.0000 | 0.0001 | 0.00x |

- load at end: pure (1.3, 1.4, 1.6), nemo (1.5, 1.5, 1.6)
</details>

<details><summary>sweep-pool-birth0.5</summary>

| phase | pure | nemo | speed-up |
|---|---|---|---|
| reproduce into offspring pool | 2.4418 | 0.1470 | 16.61x |
| islands dead selection | 0.3265 | 0.4971 | 0.66x |
| mainland sexual birth | 0.0771 | 0.0051 | 15.12x |
| colonization | 0.0485 | 0.0260 | 1.87x |
| dispersal within -> immigrant pool | 0.0210 | 0.0211 | 1.00x |
| mainland dead selection | 0.0179 | 0.0294 | 0.61x |
| germinate O+I+D | 0.0162 | 0.0261 | 0.62x |
| dormancy | 0.0083 | 0.0075 | 1.10x |
| dispersal among -> immigrant pool | 0.0081 | 0.0079 | 1.03x |
| disturbance | 0.0001 | 0.0001 | 1.00x |
| clear-up O+I | 0.0001 | 0.0001 | 1.00x |

- load at end: pure (1.7, 1.5, 1.6), nemo (1.6, 1.5, 1.6)
</details>

<details><summary>sweep-marker-rain64</summary>

| phase | pure | nemo | speed-up |
|---|---|---|---|
| germinate + birth from markers | 2.0529 | 0.7987 | 2.57x |
| islands dead selection | 0.3773 | 0.4944 | 0.76x |
| mainland sexual birth | 0.0987 | 0.0084 | 11.74x |
| marker dispersal within | 0.0277 | 0.0267 | 1.04x |
| colonization | 0.0254 | 0.0267 | 0.95x |
| marker calculation | 0.0183 | 0.0181 | 1.01x |
| mainland dead selection | 0.0181 | 0.0293 | 0.62x |
| marker dispersal among | 0.0141 | 0.0137 | 1.03x |
| marker clear-up | 0.0009 | 0.0009 | 0.94x |
| disturbance | 0.0000 | 0.0001 | 0.00x |

- load at end: pure (1.4, 1.5, 1.6), nemo (1.2, 1.4, 1.6)
</details>

<details><summary>sweep-pool-rain64</summary>

| phase | pure | nemo | speed-up |
|---|---|---|---|
| reproduce into offspring pool | 5.0584 | 0.2598 | 19.47x |
| islands dead selection | 0.3205 | 0.4828 | 0.66x |
| mainland sexual birth | 0.0832 | 0.0054 | 15.41x |
| dispersal within -> immigrant pool | 0.0301 | 0.0295 | 1.02x |
| colonization | 0.0244 | 0.0308 | 0.79x |
| germinate O+I+D | 0.0215 | 0.0311 | 0.69x |
| dormancy | 0.0215 | 0.0249 | 0.86x |
| mainland dead selection | 0.0181 | 0.0287 | 0.63x |
| dispersal among -> immigrant pool | 0.0163 | 0.0157 | 1.04x |
| disturbance | 0.0001 | 0.0001 | 1.00x |
| clear-up O+I | 0.0001 | 0.0001 | 1.00x |

- load at end: pure (1.3, 1.4, 1.5), nemo (1.4, 1.4, 1.5)
</details>

<details><summary>sweep-marker-disp0.5</summary>

| phase | pure | nemo | speed-up |
|---|---|---|---|
| germinate + birth from markers | 3.7408 | 1.3598 | 2.75x |
| islands dead selection | 0.3603 | 0.5432 | 0.66x |
| mainland sexual birth | 0.0833 | 0.0103 | 8.09x |
| colonization | 0.0379 | 0.0454 | 0.83x |
| marker dispersal within | 0.0283 | 0.0280 | 1.01x |
| mainland dead selection | 0.0184 | 0.0293 | 0.63x |
| marker calculation | 0.0162 | 0.0161 | 1.01x |
| marker dispersal among | 0.0090 | 0.0087 | 1.03x |
| marker clear-up | 0.0004 | 0.0004 | 1.12x |
| disturbance | 0.0000 | 0.0001 | 0.00x |

- load at end: pure (1.5, 1.4, 1.5), nemo (1.8, 1.6, 1.6)
</details>

<details><summary>sweep-pool-disp0.5</summary>

| phase | pure | nemo | speed-up |
|---|---|---|---|
| reproduce into offspring pool | 4.1818 | 0.1590 | 26.29x |
| islands dead selection | 0.3741 | 0.5921 | 0.63x |
| colonization | 0.1364 | 0.0421 | 3.24x |
| mainland sexual birth | 0.0793 | 0.0052 | 15.40x |
| dispersal within -> immigrant pool | 0.0364 | 0.0348 | 1.04x |
| germinate O+I+D | 0.0296 | 0.0455 | 0.65x |
| mainland dead selection | 0.0186 | 0.0292 | 0.63x |
| dormancy | 0.0152 | 0.0123 | 1.24x |
| dispersal among -> immigrant pool | 0.0111 | 0.0104 | 1.06x |
| disturbance | 0.0001 | 0.0001 | 1.00x |
| clear-up O+I | 0.0001 | 0.0001 | 1.00x |

- load at end: pure (2.0, 1.7, 1.6), nemo (2.1, 1.8, 1.7)
</details>

<details><summary>sweep-marker-dist0.01</summary>

| phase | pure | nemo | speed-up |
|---|---|---|---|
| germinate + birth from markers | 1.9235 | 0.7306 | 2.63x |
| islands dead selection | 0.3373 | 0.4503 | 0.75x |
| mainland sexual birth | 0.0864 | 0.0083 | 10.41x |
| colonization | 0.0302 | 0.0479 | 0.63x |
| marker dispersal within | 0.0268 | 0.0258 | 1.04x |
| marker calculation | 0.0259 | 0.0162 | 1.60x |
| mainland dead selection | 0.0180 | 0.0292 | 0.62x |
| marker dispersal among | 0.0136 | 0.0126 | 1.08x |
| marker clear-up | 0.0008 | 0.0008 | 1.00x |
| disturbance | 0.0007 | 0.0009 | 0.82x |

- load at end: pure (2.0, 1.8, 1.7), nemo (1.9, 1.8, 1.7)
</details>

<details><summary>sweep-pool-dist0.01</summary>

| phase | pure | nemo | speed-up |
|---|---|---|---|
| reproduce into offspring pool | 5.0058 | 0.1789 | 27.99x |
| islands dead selection | 0.3249 | 0.5161 | 0.63x |
| mainland sexual birth | 0.0815 | 0.0052 | 15.83x |
| colonization | 0.0425 | 0.0369 | 1.15x |
| dispersal within -> immigrant pool | 0.0294 | 0.0286 | 1.03x |
| germinate O+I+D | 0.0220 | 0.0316 | 0.70x |
| dormancy | 0.0203 | 0.0221 | 0.92x |
| mainland dead selection | 0.0181 | 0.0289 | 0.63x |
| dispersal among -> immigrant pool | 0.0155 | 0.0149 | 1.04x |
| disturbance | 0.0009 | 0.0011 | 0.83x |
| clear-up O+I | 0.0001 | 0.0001 | 1.00x |

- load at end: pure (1.6, 1.8, 1.7), nemo (1.8, 1.8, 1.7)
</details>

<details><summary>sweep-marker-loci5</summary>

| phase | pure | nemo | speed-up |
|---|---|---|---|
| germinate + birth from markers | 1.2402 | 0.7970 | 1.56x |
| islands dead selection | 0.3800 | 0.4852 | 0.78x |
| mainland sexual birth | 0.0320 | 0.0047 | 6.81x |
| marker dispersal within | 0.0277 | 0.0267 | 1.04x |
| colonization | 0.0245 | 0.0262 | 0.94x |
| marker calculation | 0.0186 | 0.0180 | 1.03x |
| mainland dead selection | 0.0183 | 0.0286 | 0.64x |
| marker dispersal among | 0.0142 | 0.0138 | 1.03x |
| marker clear-up | 0.0008 | 0.0009 | 0.94x |
| disturbance | 0.0000 | 0.0001 | 0.00x |

- load at end: pure (2.1, 1.9, 1.7), nemo (1.8, 1.9, 1.7)
</details>

<details><summary>sweep-pool-loci5</summary>

| phase | pure | nemo | speed-up |
|---|---|---|---|
| reproduce into offspring pool | 2.0149 | 0.2495 | 8.07x |
| islands dead selection | 0.3176 | 0.4886 | 0.65x |
| mainland sexual birth | 0.0439 | 0.0050 | 8.86x |
| colonization | 0.0399 | 0.0307 | 1.30x |
| dispersal within -> immigrant pool | 0.0299 | 0.0296 | 1.01x |
| germinate O+I+D | 0.0216 | 0.0314 | 0.69x |
| dormancy | 0.0216 | 0.0251 | 0.86x |
| mainland dead selection | 0.0181 | 0.0294 | 0.61x |
| dispersal among -> immigrant pool | 0.0162 | 0.0158 | 1.03x |
| disturbance | 0.0001 | 0.0001 | 1.00x |
| clear-up O+I | 0.0001 | 0.0001 | 1.00x |

- load at end: pure (1.5, 1.8, 1.7), nemo (1.3, 1.7, 1.7)
</details>

<details><summary>sweep-marker-loci80</summary>

| phase | pure | nemo | speed-up |
|---|---|---|---|
| germinate + birth from markers | 4.9190 | 0.8265 | 5.95x |
| islands dead selection | 0.3432 | 0.4862 | 0.71x |
| mainland sexual birth | 0.2497 | 0.0060 | 41.61x |
| marker dispersal within | 0.0273 | 0.0268 | 1.02x |
| colonization | 0.0253 | 0.0266 | 0.95x |
| marker calculation | 0.0187 | 0.0181 | 1.03x |
| mainland dead selection | 0.0179 | 0.0285 | 0.63x |
| marker dispersal among | 0.0141 | 0.0137 | 1.03x |
| marker clear-up | 0.0009 | 0.0009 | 0.94x |
| disturbance | 0.0000 | 0.0001 | 0.00x |

- load at end: pure (1.6, 1.6, 1.6), nemo (1.6, 1.6, 1.6)
</details>

<details><summary>sweep-pool-loci80</summary>

| phase | pure | nemo | speed-up |
|---|---|---|---|
| reproduce into offspring pool | 15.1142 | 0.3099 | 48.78x |
| islands dead selection | 0.3236 | 0.4792 | 0.68x |
| mainland sexual birth | 0.2429 | 0.0062 | 39.19x |
| dispersal within -> immigrant pool | 0.0303 | 0.0295 | 1.03x |
| colonization | 0.0245 | 0.0311 | 0.79x |
| dormancy | 0.0222 | 0.0245 | 0.91x |
| germinate O+I+D | 0.0217 | 0.0312 | 0.70x |
| mainland dead selection | 0.0182 | 0.0288 | 0.63x |
| dispersal among -> immigrant pool | 0.0164 | 0.0158 | 1.04x |
| disturbance | 0.0001 | 0.0001 | 1.00x |
| clear-up O+I | 0.0001 | 0.0001 | 1.00x |

- load at end: pure (1.6, 1.5, 1.6), nemo (1.5, 1.5, 1.6)
</details>

<details><summary>sweep-marker-loci320</summary>

| phase | pure | nemo | speed-up |
|---|---|---|---|
| germinate + birth from markers | 16.7863 | 0.8716 | 19.26x |
| mainland sexual birth | 0.8908 | 0.0085 | 105.42x |
| islands dead selection | 0.3943 | 0.4889 | 0.81x |
| marker dispersal within | 0.0275 | 0.0267 | 1.03x |
| colonization | 0.0250 | 0.0377 | 0.66x |
| marker calculation | 0.0192 | 0.0181 | 1.06x |
| mainland dead selection | 0.0189 | 0.0286 | 0.66x |
| marker dispersal among | 0.0149 | 0.0138 | 1.08x |
| marker clear-up | 0.0009 | 0.0009 | 1.00x |
| disturbance | 0.0000 | 0.0001 | 0.00x |

- load at end: pure (2.8, 1.9, 1.7), nemo (3.0, 2.0, 1.7)
</details>

<details><summary>sweep-pool-loci320</summary>

| phase | pure | nemo | speed-up |
|---|---|---|---|
| reproduce into offspring pool | 53.2208 | 0.4363 | 121.98x |
| mainland sexual birth | 0.8769 | 0.0087 | 100.22x |
| islands dead selection | 0.3281 | 0.4838 | 0.68x |
| dispersal within -> immigrant pool | 0.0303 | 0.0292 | 1.04x |
| colonization | 0.0245 | 0.0312 | 0.79x |
| dormancy | 0.0229 | 0.0244 | 0.94x |
| germinate O+I+D | 0.0219 | 0.0316 | 0.69x |
| mainland dead selection | 0.0186 | 0.0290 | 0.64x |
| dispersal among -> immigrant pool | 0.0172 | 0.0157 | 1.10x |
| disturbance | 0.0001 | 0.0001 | 1.00x |
| clear-up O+I | 0.0001 | 0.0001 | 1.00x |

- load at end: pure (3.4, 3.5, 3.1), nemo (3.3, 3.5, 3.2)
</details>

<details><summary>sweep-marker-traits1</summary>

| phase | pure | nemo | speed-up |
|---|---|---|---|
| germinate + birth from markers | 1.2743 | 0.6587 | 1.93x |
| islands dead selection | 0.2718 | 0.4136 | 0.66x |
| mainland sexual birth | 0.0373 | 0.0046 | 8.20x |
| marker dispersal within | 0.0275 | 0.0271 | 1.01x |
| colonization | 0.0234 | 0.0239 | 0.98x |
| marker calculation | 0.0170 | 0.0172 | 0.99x |
| mainland dead selection | 0.0153 | 0.0245 | 0.63x |
| marker dispersal among | 0.0145 | 0.0143 | 1.01x |
| marker clear-up | 0.0009 | 0.0009 | 0.89x |
| disturbance | 0.0000 | 0.0001 | 0.00x |

- load at end: pure (3.4, 3.5, 3.2), nemo (3.6, 3.5, 3.2)
</details>

<details><summary>sweep-pool-traits1</summary>

| phase | pure | nemo | speed-up |
|---|---|---|---|
| reproduce into offspring pool | 2.5978 | 0.2106 | 12.34x |
| islands dead selection | 0.2870 | 0.4329 | 0.66x |
| mainland sexual birth | 0.0505 | 0.0048 | 10.53x |
| dispersal within -> immigrant pool | 0.0301 | 0.0298 | 1.01x |
| colonization | 0.0230 | 0.0565 | 0.41x |
| germinate O+I+D | 0.0207 | 0.0291 | 0.71x |
| dormancy | 0.0181 | 0.0273 | 0.66x |
| dispersal among -> immigrant pool | 0.0167 | 0.0166 | 1.01x |
| mainland dead selection | 0.0158 | 0.0244 | 0.65x |
| disturbance | 0.0001 | 0.0001 | 1.00x |
| clear-up O+I | 0.0001 | 0.0001 | 1.00x |

- load at end: pure (3.0, 3.4, 3.2), nemo (3.0, 3.4, 3.1)
</details>

<details><summary>sweep-marker-traits4</summary>

| phase | pure | nemo | speed-up |
|---|---|---|---|
| germinate + birth from markers | 3.4453 | 0.9466 | 3.64x |
| islands dead selection | 0.3974 | 0.6253 | 0.64x |
| mainland sexual birth | 0.1556 | 0.0057 | 27.05x |
| marker dispersal within | 0.0274 | 0.0267 | 1.03x |
| colonization | 0.0263 | 0.0295 | 0.89x |
| mainland dead selection | 0.0228 | 0.0379 | 0.60x |
| marker calculation | 0.0185 | 0.0180 | 1.03x |
| marker dispersal among | 0.0141 | 0.0134 | 1.06x |
| marker clear-up | 0.0008 | 0.0009 | 0.94x |
| disturbance | 0.0000 | 0.0001 | 0.00x |

- load at end: pure (2.9, 3.3, 3.1), nemo (3.5, 3.4, 3.2)
</details>

<details><summary>sweep-pool-traits4</summary>

| phase | pure | nemo | speed-up |
|---|---|---|---|
| reproduce into offspring pool | 8.2663 | 0.2660 | 31.08x |
| islands dead selection | 0.3971 | 0.6278 | 0.63x |
| mainland sexual birth | 0.1428 | 0.0060 | 23.60x |
| colonization | 0.0438 | 0.0338 | 1.30x |
| dispersal within -> immigrant pool | 0.0299 | 0.0290 | 1.03x |
| dormancy | 0.0264 | 0.0230 | 1.15x |
| germinate O+I+D | 0.0230 | 0.0414 | 0.56x |
| mainland dead selection | 0.0224 | 0.0378 | 0.59x |
| dispersal among -> immigrant pool | 0.0163 | 0.0152 | 1.07x |
| disturbance | 0.0001 | 0.0001 | 1.00x |
| clear-up O+I | 0.0001 | 0.0001 | 1.00x |

- load at end: pure (3.6, 3.5, 3.3), nemo (3.4, 3.5, 3.3)
</details>

<details><summary>sweep-marker-mut0.01</summary>

| phase | pure | nemo | speed-up |
|---|---|---|---|
| germinate + birth from markers | 2.0807 | 0.8722 | 2.39x |
| islands dead selection | 0.3508 | 0.5252 | 0.67x |
| mainland sexual birth | 0.0968 | 0.0090 | 10.70x |
| marker dispersal within | 0.0271 | 0.0264 | 1.03x |
| colonization | 0.0255 | 0.0286 | 0.89x |
| marker calculation | 0.0185 | 0.0175 | 1.06x |
| mainland dead selection | 0.0181 | 0.0293 | 0.61x |
| marker dispersal among | 0.0141 | 0.0136 | 1.04x |
| marker clear-up | 0.0008 | 0.0009 | 0.94x |
| disturbance | 0.0000 | 0.0001 | 0.00x |

- load at end: pure (3.4, 3.5, 3.3), nemo (3.4, 3.4, 3.3)
</details>

<details><summary>sweep-pool-mut0.01</summary>

| phase | pure | nemo | speed-up |
|---|---|---|---|
| reproduce into offspring pool | 4.5732 | 0.1919 | 23.83x |
| islands dead selection | 0.3251 | 0.5725 | 0.57x |
| mainland sexual birth | 0.0762 | 0.0120 | 6.33x |
| dispersal within -> immigrant pool | 0.0297 | 0.0289 | 1.03x |
| colonization | 0.0262 | 0.0278 | 0.94x |
| germinate O+I+D | 0.0227 | 0.0343 | 0.66x |
| dormancy | 0.0211 | 0.0234 | 0.90x |
| mainland dead selection | 0.0181 | 0.0303 | 0.59x |
| dispersal among -> immigrant pool | 0.0162 | 0.0155 | 1.05x |
| disturbance | 0.0001 | 0.0001 | 1.00x |
| clear-up O+I | 0.0001 | 0.0001 | 1.00x |

- load at end: pure (3.0, 3.3, 3.2), nemo (3.1, 3.3, 3.2)
</details>

