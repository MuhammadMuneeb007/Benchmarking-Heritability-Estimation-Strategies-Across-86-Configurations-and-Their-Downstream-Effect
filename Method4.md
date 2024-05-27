
Table Explanation

- **clump_p1**: Significance threshold for clumping. Indicates the p-value threshold for selecting index SNPs during clumping.
- **clump_r2**: R-squared threshold for clumping. Represents the linkage disequilibrium (LD) threshold for pruning variants.
- **clump_kb**: Distance in kilobases (kb) for clumping. Defines the window size around each index SNP to include variants in LD.
- **p_window_size**: Size of the sliding window in base pairs for calculating heritability.
- **p_slide_size**: Size of the step in base pairs for sliding the window.
- **p_LD_threshold**: LD threshold for pruning.  
- **numberofpca**: Number of principal components (PCA) included in the analysis.
- **h2**: Estimated heritability (h²). Indicates the proportion of variance in the trait that can be attributed to genetic factors.
- **numberofvariants**: Number of genetic variants (SNPs) considered in the analysis.
- **h2model**: Model used for heritability estimation (e.g., HE regression with genotype or genotype + covariate).
- **clumpprune**: Indicates whether clumping and pruning were applied (yes) or not (no).



Phenotype:  body_mass_index_bmi
|    |   clump_p1 |   clump_r2 |   clump_kb |   p_window_size |   p_slide_size |   p_LD_threshold |   numberofpca |        h2 |   numberofvariants |   relatedmatrix |   data | h2model                                  | clumpprune   |
|---:|-----------:|-----------:|-----------:|----------------:|---------------:|-----------------:|--------------:|----------:|-------------------:|----------------:|-------:|:-----------------------------------------|:-------------|
|  0 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.155912  |            27515   |               0 |      0 | HE regression_genotype                   | yes          |
|  1 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.148626  |            27182.2 |               0 |      0 | HE regression_genotype_covariate         | yes          |
|  2 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.156778  |            27513.2 |               0 |      0 | HE regression_genotype_covariate_pca     | yes          |
|  3 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.129021  |            82362   |               0 |      0 | HE regression_genotype                   | no           |
|  4 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.114299  |            81650   |               0 |      0 | HE regression_genotype_covariate         | no           |
|  5 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.129763  |            82362   |               0 |      0 | HE regression_genotype_covariate_pca     | no           |
|  6 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.0771655 |           119641   |               0 |      0 | REML AI algorithm_genotype               | yes          |
|  7 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.0677847 |            86715.2 |               0 |      0 | REML AI algorithm_genotype_covariate     | yes          |
|  8 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.0588795 |           101356   |               0 |      0 | REML AI algorithm_genotype_covariate_pca | yes          |
|  9 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.117045  |           171793   |               0 |      0 | REML AI algorithm_genotype               | no           |
| 10 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.075998  |           211400   |               0 |      0 | REML AI algorithm_genotype_covariate     | no           |
| 11 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.0962617 |           322300   |               0 |      0 | REML AI algorithm_genotype_covariate_pca | no           |



Phenotype:  hypertension
|    |   clump_p1 |   clump_r2 |   clump_kb |   p_window_size |   p_slide_size |   p_LD_threshold |   numberofpca |        h2 |   numberofvariants |   relatedmatrix |   data | h2model                                  | clumpprune   |
|---:|-----------:|-----------:|-----------:|----------------:|---------------:|-----------------:|--------------:|----------:|-------------------:|----------------:|-------:|:-----------------------------------------|:-------------|
|  0 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.127176  |            25824.2 |               0 |      0 | HE regression_genotype                   | yes          |
|  1 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.125457  |            25836   |               0 |      0 | HE regression_genotype_covariate         | yes          |
|  2 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.12838   |            20654.2 |               0 |      0 | HE regression_genotype_covariate_pca     | yes          |
|  3 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.0609309 |            78742.2 |               0 |      0 | HE regression_genotype                   | no           |
|  4 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.0606503 |            78742.2 |               0 |      0 | HE regression_genotype_covariate         | no           |
|  5 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.0612456 |            78742.2 |               0 |      0 | HE regression_genotype_covariate_pca     | no           |
|  6 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.0629598 |           120997   |               0 |      0 | REML AI algorithm_genotype               | yes          |
|  7 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.0413074 |            70264.4 |               0 |      0 | REML AI algorithm_genotype_covariate     | yes          |
|  8 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.035073  |            67618   |               0 |      0 | REML AI algorithm_genotype_covariate_pca | yes          |
|  9 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.105678  |           135269   |               0 |      0 | REML AI algorithm_genotype               | no           |
| 10 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.0689101 |           140297   |               0 |      0 | REML AI algorithm_genotype_covariate     | no           |
| 11 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.12375   |           176021   |               0 |      0 | REML AI algorithm_genotype_covariate_pca | no           |



Phenotype:  depression
|    |   clump_p1 |   clump_r2 |   clump_kb |   p_window_size |   p_slide_size |   p_LD_threshold |   numberofpca |        h2 |   numberofvariants |   relatedmatrix |   data | h2model                                  | clumpprune   |
|---:|-----------:|-----------:|-----------:|----------------:|---------------:|-----------------:|--------------:|----------:|-------------------:|----------------:|-------:|:-----------------------------------------|:-------------|
|  0 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.155898  |             168684 |               0 |      0 | HE regression_genotype                   | yes          |
|  1 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.121622  |             134883 |               0 |      0 | HE regression_genotype_covariate         | yes          |
|  2 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.124088  |             134932 |               0 |      0 | HE regression_genotype_covariate_pca     | yes          |
|  3 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.081918  |             531018 |               0 |      0 | HE regression_genotype                   | no           |
|  4 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.0808639 |             531018 |               0 |      0 | HE regression_genotype_covariate         | no           |
|  5 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.0817495 |             531018 |               0 |      0 | HE regression_genotype_covariate_pca     | no           |
|  6 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.0816398 |             530484 |               0 |      0 | REML AI algorithm_genotype               | yes          |
|  7 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.0548996 |             428789 |               0 |      0 | REML AI algorithm_genotype_covariate     | yes          |
|  8 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.0297016 |             552259 |               0 |      0 | REML AI algorithm_genotype_covariate_pca | yes          |
|  9 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.0370766 |             553038 |               0 |      0 | REML AI algorithm_genotype               | no           |
| 10 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.066976  |             553330 |               0 |      0 | REML AI algorithm_genotype_covariate     | no           |
| 11 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.066976  |             553330 |               0 |      0 | REML AI algorithm_genotype_covariate_pca | no           |



Phenotype:  asthma
|    |   clump_p1 |   clump_r2 |   clump_kb |   p_window_size |   p_slide_size |   p_LD_threshold |   numberofpca |        h2 |   numberofvariants |   relatedmatrix |   data | h2model                                  | clumpprune   |
|---:|-----------:|-----------:|-----------:|----------------:|---------------:|-----------------:|--------------:|----------:|-------------------:|----------------:|-------:|:-----------------------------------------|:-------------|
|  0 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.0215642 |            16519.8 |               0 |      0 | HE regression_genotype                   | yes          |
|  1 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.0356364 |            27198.2 |               0 |      0 | HE regression_genotype_covariate         | yes          |
|  2 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.036281  |            27198.2 |               0 |      0 | HE regression_genotype_covariate_pca     | yes          |
|  3 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.0260743 |            81159.8 |               0 |      0 | HE regression_genotype                   | no           |
|  4 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.0260693 |            81159.8 |               0 |      0 | HE regression_genotype_covariate         | no           |
|  5 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.0263832 |            81159.8 |               0 |      0 | HE regression_genotype_covariate_pca     | no           |
|  6 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.0289468 |            82657.6 |               0 |      0 | REML AI algorithm_genotype               | yes          |
|  7 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.0987152 |           117329   |               0 |      0 | REML AI algorithm_genotype_covariate     | yes          |
|  8 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.0794329 |           102460   |               0 |      0 | REML AI algorithm_genotype_covariate_pca | yes          |
|  9 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.0860216 |           175425   |               0 |      0 | REML AI algorithm_genotype               | no           |
| 10 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.0995396 |           140797   |               0 |      0 | REML AI algorithm_genotype_covariate     | no           |
| 11 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.0817641 |           530558   |               0 |      0 | REML AI algorithm_genotype_covariate_pca | no           |



Phenotype:  osteoarthritis
|    |   clump_p1 |   clump_r2 |   clump_kb |   p_window_size |   p_slide_size |   p_LD_threshold |   numberofpca |          h2 |   numberofvariants |   relatedmatrix |   data | h2model                                  | clumpprune   |
|---:|-----------:|-----------:|-----------:|----------------:|---------------:|-----------------:|--------------:|------------:|-------------------:|----------------:|-------:|:-----------------------------------------|:-------------|
|  0 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.00080969  |              949.4 |               0 |      0 | HE regression_genotype                   | yes          |
|  1 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.000792697 |              944.8 |               0 |      0 | HE regression_genotype_covariate         | yes          |
|  2 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.000643442 |             5926.2 |               0 |      0 | HE regression_genotype_covariate_pca     | yes          |
|  3 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.0941081   |            22251.6 |               0 |      0 | HE regression_genotype                   | no           |
|  4 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.000822791 |             1245.8 |               0 |      0 | HE regression_genotype_covariate         | no           |
|  5 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.000828757 |             1245.8 |               0 |      0 | HE regression_genotype_covariate_pca     | no           |
|  6 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.038467    |            32955.6 |               0 |      0 | REML AI algorithm_genotype               | yes          |
|  7 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.189347    |           149123   |               0 |      0 | REML AI algorithm_genotype_covariate     | yes          |
|  8 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.050928    |            66209   |               0 |      0 | REML AI algorithm_genotype_covariate_pca | yes          |
|  9 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.0414835   |            66586.8 |               0 |      0 | REML AI algorithm_genotype               | no           |
| 10 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.0221295   |            51308.2 |               0 |      0 | REML AI algorithm_genotype_covariate     | no           |
| 11 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.110537    |           136236   |               0 |      0 | REML AI algorithm_genotype_covariate_pca | no           |



Phenotype:  high_cholesterol
|    |   clump_p1 |   clump_r2 |   clump_kb |   p_window_size |   p_slide_size |   p_LD_threshold |   numberofpca |        h2 |   numberofvariants |   relatedmatrix |   data | h2model                                  | clumpprune   |
|---:|-----------:|-----------:|-----------:|----------------:|---------------:|-----------------:|--------------:|----------:|-------------------:|----------------:|-------:|:-----------------------------------------|:-------------|
|  0 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.100186  |             140721 |               0 |      0 | HE regression_genotype                   | yes          |
|  1 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.0973384 |             140705 |               0 |      0 | HE regression_genotype_covariate         | yes          |
|  2 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.124326  |             175896 |               0 |      0 | HE regression_genotype_covariate_pca     | yes          |
|  3 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.0674643 |             553320 |               0 |      0 | HE regression_genotype                   | no           |
|  4 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.0664132 |             553320 |               0 |      0 | HE regression_genotype_covariate         | no           |
|  5 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.0672279 |             553320 |               0 |      0 | HE regression_genotype_covariate_pca     | no           |
|  6 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.066976  |             553330 |               0 |      0 | REML AI algorithm_genotype               | yes          |
|  7 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.066976  |             553330 |               0 |      0 | REML AI algorithm_genotype_covariate     | yes          |
|  8 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.066976  |             553330 |               0 |      0 | REML AI algorithm_genotype_covariate_pca | yes          |
|  9 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.066976  |             553330 |               0 |      0 | REML AI algorithm_genotype               | no           |
| 10 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.066976  |             553330 |               0 |      0 | REML AI algorithm_genotype_covariate     | no           |
| 11 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.066976  |             553330 |               0 |      0 | REML AI algorithm_genotype_covariate_pca | no           |



Phenotype:  irritable_bowel_syndrome
|    |   clump_p1 |   clump_r2 |   clump_kb |   p_window_size |   p_slide_size |   p_LD_threshold |   numberofpca |        h2 |   numberofvariants |   relatedmatrix |   data | h2model                                  | clumpprune   |
|---:|-----------:|-----------:|-----------:|----------------:|---------------:|-----------------:|--------------:|----------:|-------------------:|----------------:|-------:|:-----------------------------------------|:-------------|
|  0 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.0866639 |             175387 |               0 |      0 | HE regression_genotype                   | yes          |
|  1 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.0846073 |             175387 |               0 |      0 | HE regression_genotype_covariate         | yes          |
|  2 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.0861963 |             175387 |               0 |      0 | HE regression_genotype_covariate_pca     | yes          |
|  3 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.0295101 |             552363 |               0 |      0 | HE regression_genotype                   | no           |
|  4 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.0290978 |             552363 |               0 |      0 | HE regression_genotype_covariate         | no           |
|  5 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.0295235 |             552363 |               0 |      0 | HE regression_genotype_covariate_pca     | no           |
|  6 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.023759  |             442073 |               0 |      0 | REML AI algorithm_genotype               | yes          |
|  7 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.038626  |             442494 |               0 |      0 | REML AI algorithm_genotype_covariate     | yes          |
|  8 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.066976  |             553330 |               0 |      0 | REML AI algorithm_genotype_covariate_pca | yes          |
|  9 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.066976  |             553330 |               0 |      0 | REML AI algorithm_genotype               | no           |
| 10 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.066976  |             553330 |               0 |      0 | REML AI algorithm_genotype_covariate     | no           |
| 11 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.066976  |             553330 |               0 |      0 | REML AI algorithm_genotype_covariate_pca | no           |



Phenotype:  hypothyroidism_myxoedema
|    |   clump_p1 |   clump_r2 |   clump_kb |   p_window_size |   p_slide_size |   p_LD_threshold |   numberofpca |        h2 |   numberofvariants |   relatedmatrix |   data | h2model                                  | clumpprune   |
|---:|-----------:|-----------:|-----------:|----------------:|---------------:|-----------------:|--------------:|----------:|-------------------:|----------------:|-------:|:-----------------------------------------|:-------------|
|  0 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.0636228 |            27895.6 |               0 |      0 | HE regression_genotype                   | yes          |
|  1 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.0743282 |            27908.8 |               0 |      0 | HE regression_genotype_covariate         | yes          |
|  2 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.0757672 |            27908.8 |               0 |      0 | HE regression_genotype_covariate_pca     | yes          |
|  3 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.0481431 |            66736.2 |               0 |      0 | HE regression_genotype                   | no           |
|  4 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.0598492 |            83411.6 |               0 |      0 | HE regression_genotype_covariate         | no           |
|  5 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.0605275 |            83411.6 |               0 |      0 | HE regression_genotype_covariate_pca     | no           |
|  6 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.0672329 |           140242   |               0 |      0 | REML AI algorithm_genotype               | yes          |
|  7 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.0608334 |            78638   |               0 |      0 | REML AI algorithm_genotype_covariate     | yes          |
|  8 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.0845161 |           175446   |               0 |      0 | REML AI algorithm_genotype_covariate_pca | yes          |
|  9 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.0643402 |           119092   |               0 |      0 | REML AI algorithm_genotype               | no           |
| 10 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.155812  |           168499   |               0 |      0 | REML AI algorithm_genotype_covariate     | no           |
| 11 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0         |                0   |               0 |      0 | REML AI algorithm_genotype_covariate_pca | no           |



Phenotype:  gastro_oesophageal_reflux_gord_gastric_reflux
|    |   clump_p1 |   clump_r2 |   clump_kb |   p_window_size |   p_slide_size |   p_LD_threshold |   numberofpca |         h2 |   numberofvariants |   relatedmatrix |   data | h2model                                  | clumpprune   |
|---:|-----------:|-----------:|-----------:|----------------:|---------------:|-----------------:|--------------:|-----------:|-------------------:|----------------:|-------:|:-----------------------------------------|:-------------|
|  0 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.0111209  |            16724.6 |               0 |      0 | HE regression_genotype                   | yes          |
|  1 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.0184215  |            27865.6 |               0 |      0 | HE regression_genotype_covariate         | yes          |
|  2 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.0186843  |            27866   |               0 |      0 | HE regression_genotype_covariate_pca     | yes          |
|  3 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.00755422 |            66681.4 |               0 |      0 | HE regression_genotype                   | no           |
|  4 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.00750755 |            66699.8 |               0 |      0 | HE regression_genotype_covariate         | no           |
|  5 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.00951608 |            83369.2 |               0 |      0 | HE regression_genotype_covariate_pca     | no           |
|  6 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.0672329  |           140242   |               0 |      0 | REML AI algorithm_genotype               | yes          |
|  7 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.0608334  |            78638   |               0 |      0 | REML AI algorithm_genotype_covariate     | yes          |
|  8 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.0845161  |           175446   |               0 |      0 | REML AI algorithm_genotype_covariate_pca | yes          |
|  9 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.0831802  |           138073   |               0 |      0 | REML AI algorithm_genotype               | no           |
| 10 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.155812   |           168499   |               0 |      0 | REML AI algorithm_genotype_covariate     | no           |
| 11 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0          |                0   |               0 |      0 | REML AI algorithm_genotype_covariate_pca | no           |



Phenotype:  migraine
|    |   clump_p1 |   clump_r2 |   clump_kb |   p_window_size |   p_slide_size |   p_LD_threshold |   numberofpca |         h2 |   numberofvariants |   relatedmatrix |   data | h2model                                  | clumpprune   |
|---:|-----------:|-----------:|-----------:|----------------:|---------------:|-----------------:|--------------:|-----------:|-------------------:|----------------:|-------:|:-----------------------------------------|:-------------|
|  0 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.0654348  |             169812 |               0 |      0 | HE regression_genotype                   | yes          |
|  1 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.0638261  |             169812 |               0 |      0 | HE regression_genotype_covariate         | yes          |
|  2 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.0650382  |             169812 |               0 |      0 | HE regression_genotype_covariate_pca     | yes          |
|  3 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.001976   |             535177 |               0 |      0 | HE regression_genotype                   | no           |
|  4 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.00194932 |             535177 |               0 |      0 | HE regression_genotype_covariate         | no           |
|  5 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.00197174 |             535175 |               0 |      0 | HE regression_genotype_covariate_pca     | no           |
|  6 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.065332   |             424424 |               0 |      0 | REML AI algorithm_genotype               | yes          |
|  7 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.0816398  |             530484 |               0 |      0 | REML AI algorithm_genotype_covariate     | yes          |
|  8 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.017839   |             331410 |               0 |      0 | REML AI algorithm_genotype_covariate_pca | yes          |
|  9 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.0177535  |             331742 |               0 |      0 | REML AI algorithm_genotype               | no           |
| 10 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.066976   |             553330 |               0 |      0 | REML AI algorithm_genotype_covariate     | no           |
| 11 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.066976   |             553330 |               0 |      0 | REML AI algorithm_genotype_covariate_pca | no           |
