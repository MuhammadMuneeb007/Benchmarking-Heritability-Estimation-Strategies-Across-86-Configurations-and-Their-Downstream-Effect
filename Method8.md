
1. **clump_p1**: A threshold p-value for the primary SNP (single nucleotide polymorphism) in linkage disequilibrium (LD) clumping. A primary SNP is retained if its p-value is below this threshold.
2. **clump_r2**: The LD r-squared threshold for clumping. SNPs that have an LD r-squared value greater than this threshold with the primary SNP are clumped together.
3. **clump_kb**: The distance in kilobases within which SNPs are clumped around the primary SNP.
4. **p_window_size**: Pruning
5. **p_slide_size**: Pruning
6. **p_LD_threshold**: Pruning
7. **numberofpca**: The number of principal components used in the analysis.
8. **h2**: The heritability estimate of the trait being studied (BMI in this case).
9. **numberofvariants**: The number of genetic variants considered in the analysis.
10. **alphamodelvalue**: Ignore it
11. **relatedmatrix**: Ignore it
12. **data**: Ignore it
13. **h2model**: Specifies the heritability model used in the analysis, such as "LDAK_Human/LDAK_GCTA/LDAK_ALPHA".
14. **clumpprune**: Indicates whether clumping/pruning was performed (yes or no).



Phenotype:  body_mass_index_bmi
|    |   clump_p1 |   clump_r2 |   clump_kb |   p_window_size |   p_slide_size |   p_LD_threshold |   numberofpca |       h2 |   numberofvariants |   alphamodelvalue |   relatedmatrix |   data | h2model             | clumpprune   |
|---:|-----------:|-----------:|-----------:|----------------:|---------------:|-----------------:|--------------:|---------:|-------------------:|------------------:|----------------:|-------:|:--------------------|:-------------|
|  0 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.161787 |            10445.4 |                 0 |               0 |      0 | LDAK_human_yes      | yes          |
|  1 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.158402 |            10445.4 |                 0 |               0 |      0 | LDAK_GCTA_yes       | yes          |
|  2 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.158402 |            10445.4 |                 0 |               0 |      0 | LDAK_BLD-LDAK_yes   | yes          |
|  3 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.158402 |            10445.4 |                 0 |               0 |      0 | LDAK_alpha_1.0_yes  | yes          |
|  4 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.158402 |            10445.4 |                 0 |               0 |      0 | LDAK_alpha_0.5_yes  | yes          |
|  5 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.158402 |            10445.4 |                 0 |               0 |      0 | LDAK_alpha_0.0_yes  | yes          |
|  6 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.158402 |            10445.4 |                 0 |               0 |      0 | LDAK_alpha_-0.5_yes | yes          |
|  7 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.158402 |            10445.4 |                 0 |               0 |      0 | LDAK_alpha_yes      | yes          |



Phenotype:  hypertension
|    |   clump_p1 |   clump_r2 |   clump_kb |   p_window_size |   p_slide_size |   p_LD_threshold |   numberofpca |       h2 |   numberofvariants |   alphamodelvalue |   relatedmatrix |   data | h2model             | clumpprune   |
|---:|-----------:|-----------:|-----------:|----------------:|---------------:|-----------------:|--------------:|---------:|-------------------:|------------------:|----------------:|-------:|:--------------------|:-------------|
|  0 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.11167  |             9153.4 |                 0 |               0 |      0 | LDAK_human_yes      | yes          |
|  1 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.117757 |             9153.4 |                 0 |               0 |      0 | LDAK_GCTA_yes       | yes          |
|  2 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.117757 |             9153.4 |                 0 |               0 |      0 | LDAK_BLD-LDAK_yes   | yes          |
|  3 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.117757 |             9153.4 |                 0 |               0 |      0 | LDAK_alpha_1.0_yes  | yes          |
|  4 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.117757 |             9153.4 |                 0 |               0 |      0 | LDAK_alpha_0.5_yes  | yes          |
|  5 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.117757 |             9153.4 |                 0 |               0 |      0 | LDAK_alpha_0.0_yes  | yes          |
|  6 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.117757 |             9153.4 |                 0 |               0 |      0 | LDAK_alpha_-0.5_yes | yes          |
|  7 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.117757 |             9153.4 |                 0 |               0 |      0 | LDAK_alpha_yes      | yes          |



Phenotype:  depression
|    |   clump_p1 |   clump_r2 |   clump_kb |   p_window_size |   p_slide_size |   p_LD_threshold |   numberofpca |       h2 |   numberofvariants |   alphamodelvalue |   relatedmatrix |   data | h2model             | clumpprune   |
|---:|-----------:|-----------:|-----------:|----------------:|---------------:|-----------------:|--------------:|---------:|-------------------:|------------------:|----------------:|-------:|:--------------------|:-------------|
|  0 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.149579 |             173904 |                 0 |               0 |      0 | LDAK_human_yes      | yes          |
|  1 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.156389 |             173904 |                 0 |               0 |      0 | LDAK_GCTA_yes       | yes          |
|  2 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.156389 |             173904 |                 0 |               0 |      0 | LDAK_BLD-LDAK_yes   | yes          |
|  3 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.156389 |             173904 |                 0 |               0 |      0 | LDAK_alpha_1.0_yes  | yes          |
|  4 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.156389 |             173904 |                 0 |               0 |      0 | LDAK_alpha_0.5_yes  | yes          |
|  5 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.156389 |             173904 |                 0 |               0 |      0 | LDAK_alpha_0.0_yes  | yes          |
|  6 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.156389 |             173904 |                 0 |               0 |      0 | LDAK_alpha_-0.5_yes | yes          |
|  7 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.156389 |             173904 |                 0 |               0 |      0 | LDAK_alpha_yes      | yes          |



Phenotype:  asthma
|    |   clump_p1 |   clump_r2 |   clump_kb |   p_window_size |   p_slide_size |   p_LD_threshold |   numberofpca |        h2 |   numberofvariants |   alphamodelvalue |   relatedmatrix |   data | h2model             | clumpprune   |
|---:|-----------:|-----------:|-----------:|----------------:|---------------:|-----------------:|--------------:|----------:|-------------------:|------------------:|----------------:|-------:|:--------------------|:-------------|
|  0 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.0343506 |            10321.2 |                 0 |               0 |      0 | LDAK_human_yes      | yes          |
|  1 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.036129  |            10321.2 |                 0 |               0 |      0 | LDAK_GCTA_yes       | yes          |
|  2 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.036129  |            10321.2 |                 0 |               0 |      0 | LDAK_BLD-LDAK_yes   | yes          |
|  3 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.036129  |            10321.2 |                 0 |               0 |      0 | LDAK_alpha_1.0_yes  | yes          |
|  4 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.036129  |            10321.2 |                 0 |               0 |      0 | LDAK_alpha_0.5_yes  | yes          |
|  5 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.036129  |            10321.2 |                 0 |               0 |      0 | LDAK_alpha_0.0_yes  | yes          |
|  6 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.036129  |            10321.2 |                 0 |               0 |      0 | LDAK_alpha_-0.5_yes | yes          |
|  7 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.036129  |            10321.2 |                 0 |               0 |      0 | LDAK_alpha_yes      | yes          |



Phenotype:  osteoarthritis
|    |   clump_p1 |   clump_r2 |   clump_kb |   p_window_size |   p_slide_size |   p_LD_threshold |   numberofpca |        h2 |   numberofvariants |   alphamodelvalue |   relatedmatrix |   data | h2model             | clumpprune   |
|---:|-----------:|-----------:|-----------:|----------------:|---------------:|-----------------:|--------------:|----------:|-------------------:|------------------:|----------------:|-------:|:--------------------|:-------------|
|  0 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.0012338 |              329.4 |                 0 |               0 |      0 | LDAK_human_yes      | yes          |
|  1 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.0012358 |              329.4 |                 0 |               0 |      0 | LDAK_GCTA_yes       | yes          |
|  2 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.0012358 |              329.4 |                 0 |               0 |      0 | LDAK_BLD-LDAK_yes   | yes          |
|  3 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.0012358 |              329.4 |                 0 |               0 |      0 | LDAK_alpha_1.0_yes  | yes          |
|  4 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.0012358 |              329.4 |                 0 |               0 |      0 | LDAK_alpha_0.5_yes  | yes          |
|  5 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.0012358 |              329.4 |                 0 |               0 |      0 | LDAK_alpha_0.0_yes  | yes          |
|  6 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.0012358 |              329.4 |                 0 |               0 |      0 | LDAK_alpha_-0.5_yes | yes          |
|  7 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.0012358 |              329.4 |                 0 |               0 |      0 | LDAK_alpha_yes      | yes          |



Phenotype:  high_cholesterol
|    |   clump_p1 |   clump_r2 |   clump_kb |   p_window_size |   p_slide_size |   p_LD_threshold |   numberofpca |      h2 |   numberofvariants |   alphamodelvalue |   relatedmatrix |   data | h2model             | clumpprune   |
|---:|-----------:|-----------:|-----------:|----------------:|---------------:|-----------------:|--------------:|--------:|-------------------:|------------------:|----------------:|-------:|:--------------------|:-------------|
|  0 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.11929 |             181530 |                 0 |               0 |      0 | LDAK_human_yes      | yes          |
|  1 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.12515 |             181530 |                 0 |               0 |      0 | LDAK_GCTA_yes       | yes          |
|  2 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.12515 |             181530 |                 0 |               0 |      0 | LDAK_BLD-LDAK_yes   | yes          |
|  3 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.12515 |             181530 |                 0 |               0 |      0 | LDAK_alpha_1.0_yes  | yes          |
|  4 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.12515 |             181530 |                 0 |               0 |      0 | LDAK_alpha_0.5_yes  | yes          |
|  5 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.12515 |             181530 |                 0 |               0 |      0 | LDAK_alpha_0.0_yes  | yes          |
|  6 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.12515 |             181530 |                 0 |               0 |      0 | LDAK_alpha_-0.5_yes | yes          |
|  7 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.12515 |             181530 |                 0 |               0 |      0 | LDAK_alpha_yes      | yes          |



Phenotype:  irritable_bowel_syndrome
|    |   clump_p1 |   clump_r2 |   clump_kb |   p_window_size |   p_slide_size |   p_LD_threshold |   numberofpca |        h2 |   numberofvariants |   alphamodelvalue |   relatedmatrix |   data | h2model             | clumpprune   |
|---:|-----------:|-----------:|-----------:|----------------:|---------------:|-----------------:|--------------:|----------:|-------------------:|------------------:|----------------:|-------:|:--------------------|:-------------|
|  0 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.0820036 |             180928 |                 0 |               0 |      0 | LDAK_human_yes      | yes          |
|  1 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.0868732 |             180928 |                 0 |               0 |      0 | LDAK_GCTA_yes       | yes          |
|  2 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.0868732 |             180928 |                 0 |               0 |      0 | LDAK_BLD-LDAK_yes   | yes          |
|  3 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.0868732 |             180928 |                 0 |               0 |      0 | LDAK_alpha_1.0_yes  | yes          |
|  4 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.0868732 |             180928 |                 0 |               0 |      0 | LDAK_alpha_0.5_yes  | yes          |
|  5 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.0868732 |             180928 |                 0 |               0 |      0 | LDAK_alpha_0.0_yes  | yes          |
|  6 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.0868732 |             180928 |                 0 |               0 |      0 | LDAK_alpha_-0.5_yes | yes          |
|  7 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.0868732 |             180928 |                 0 |               0 |      0 | LDAK_alpha_yes      | yes          |



Phenotype:  hypothyroidism_myxoedema
|    |   clump_p1 |   clump_r2 |   clump_kb |   p_window_size |   p_slide_size |   p_LD_threshold |   numberofpca |        h2 |   numberofvariants |   alphamodelvalue |   relatedmatrix |   data | h2model             | clumpprune   |
|---:|-----------:|-----------:|-----------:|----------------:|---------------:|-----------------:|--------------:|----------:|-------------------:|------------------:|----------------:|-------:|:--------------------|:-------------|
|  0 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.0600962 |              10590 |                 0 |               0 |      0 | LDAK_human_yes      | yes          |
|  1 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.0622988 |              10590 |                 0 |               0 |      0 | LDAK_GCTA_yes       | yes          |
|  2 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.0622988 |              10590 |                 0 |               0 |      0 | LDAK_BLD-LDAK_yes   | yes          |
|  3 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.0622988 |              10590 |                 0 |               0 |      0 | LDAK_alpha_1.0_yes  | yes          |
|  4 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.0622988 |              10590 |                 0 |               0 |      0 | LDAK_alpha_0.5_yes  | yes          |
|  5 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.0622988 |              10590 |                 0 |               0 |      0 | LDAK_alpha_0.0_yes  | yes          |
|  6 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.0622988 |              10590 |                 0 |               0 |      0 | LDAK_alpha_-0.5_yes | yes          |
|  7 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.0622988 |              10590 |                 0 |               0 |      0 | LDAK_alpha_yes      | yes          |



Phenotype:  gastro_oesophageal_reflux_gord_gastric_reflux
|    |   clump_p1 |   clump_r2 |   clump_kb |   p_window_size |   p_slide_size |   p_LD_threshold |   numberofpca |        h2 |   numberofvariants |   alphamodelvalue |   relatedmatrix |   data | h2model             | clumpprune   |
|---:|-----------:|-----------:|-----------:|----------------:|---------------:|-----------------:|--------------:|----------:|-------------------:|------------------:|----------------:|-------:|:--------------------|:-------------|
|  0 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.017674  |            10576.4 |                 0 |               0 |      0 | LDAK_human_yes      | yes          |
|  1 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.0189228 |            10576.4 |                 0 |               0 |      0 | LDAK_GCTA_yes       | yes          |
|  2 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.0189228 |            10576.4 |                 0 |               0 |      0 | LDAK_BLD-LDAK_yes   | yes          |
|  3 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.0189228 |            10576.4 |                 0 |               0 |      0 | LDAK_alpha_1.0_yes  | yes          |
|  4 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.0189228 |            10576.4 |                 0 |               0 |      0 | LDAK_alpha_0.5_yes  | yes          |
|  5 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.0189228 |            10576.4 |                 0 |               0 |      0 | LDAK_alpha_0.0_yes  | yes          |
|  6 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.0189228 |            10576.4 |                 0 |               0 |      0 | LDAK_alpha_-0.5_yes | yes          |
|  7 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.0189228 |            10576.4 |                 0 |               0 |      0 | LDAK_alpha_yes      | yes          |



Phenotype:  migraine
|    |   clump_p1 |   clump_r2 |   clump_kb |   p_window_size |   p_slide_size |   p_LD_threshold |   numberofpca |        h2 |   numberofvariants |   alphamodelvalue |   relatedmatrix |   data | h2model             | clumpprune   |
|---:|-----------:|-----------:|-----------:|----------------:|---------------:|-----------------:|--------------:|----------:|-------------------:|------------------:|----------------:|-------:|:--------------------|:-------------|
|  0 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.0620786 |             174699 |                 0 |               0 |      0 | LDAK_human_yes      | yes          |
|  1 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.0668028 |             174699 |                 0 |               0 |      0 | LDAK_GCTA_yes       | yes          |
|  2 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.0668028 |             174699 |                 0 |               0 |      0 | LDAK_BLD-LDAK_yes   | yes          |
|  3 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.0668028 |             174699 |                 0 |               0 |      0 | LDAK_alpha_1.0_yes  | yes          |
|  4 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.0668028 |             174699 |                 0 |               0 |      0 | LDAK_alpha_0.5_yes  | yes          |
|  5 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.0668028 |             174699 |                 0 |               0 |      0 | LDAK_alpha_0.0_yes  | yes          |
|  6 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.0668028 |             174699 |                 0 |               0 |      0 | LDAK_alpha_-0.5_yes | yes          |
|  7 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.0668028 |             174699 |                 0 |               0 |      0 | LDAK_alpha_yes      | yes          |

