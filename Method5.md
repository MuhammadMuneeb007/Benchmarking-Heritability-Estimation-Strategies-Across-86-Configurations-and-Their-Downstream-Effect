
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
- **h2model**: Model used for heritability estimation (e.g., REML AI algorithm with genotype or genotype + covariate).
- **clumpprune**: Indicates whether clumping and pruning were applied (yes) or not (no).




Phenotype:  body_mass_index_bmi
|    |   clump_p1 |   clump_r2 |   clump_kb |   p_window_size |   p_slide_size |   p_LD_threshold |   numberofpca |        h2 |   numberofvariants |   relatedmatrix |   data | h2model                                  | clumpprune   |
|---:|-----------:|-----------:|-----------:|----------------:|---------------:|-----------------:|--------------:|----------:|-------------------:|----------------:|-------:|:-----------------------------------------|:-------------|
|  0 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | -0.125561 |              27515 |               0 |      0 | REML AI algorithm_genotype               | yes          |
|  1 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 |  0.705895 |              27515 |               0 |      0 | REML AI algorithm_genotype_covariate     | yes          |
|  2 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 |  2.73531  |              27515 |               0 |      0 | REML AI algorithm_genotype_covariate_pca | yes          |



Phenotype:  hypertension
|    |   clump_p1 |   clump_r2 |   clump_kb |   p_window_size |   p_slide_size |   p_LD_threshold |   numberofpca |       h2 |   numberofvariants |   relatedmatrix |   data | h2model                                  | clumpprune   |
|---:|-----------:|-----------:|-----------:|----------------:|---------------:|-----------------:|--------------:|---------:|-------------------:|----------------:|-------:|:-----------------------------------------|:-------------|
|  0 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.135893 |            25824.2 |               0 |      0 | REML AI algorithm_genotype               | yes          |
|  1 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.134684 |            25824.2 |               0 |      0 | REML AI algorithm_genotype_covariate     | yes          |
|  2 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.137224 |            25824.2 |               0 |      0 | REML AI algorithm_genotype_covariate_pca | yes          |



Phenotype:  depression
|    |   clump_p1 |   clump_r2 |   clump_kb |   p_window_size |   p_slide_size |   p_LD_threshold |   numberofpca |         h2 |   numberofvariants |   relatedmatrix |   data | h2model                                  | clumpprune   |
|---:|-----------:|-----------:|-----------:|----------------:|---------------:|-----------------:|--------------:|-----------:|-------------------:|----------------:|-------:|:-----------------------------------------|:-------------|
|  0 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 |  0.0727557 |             168684 |               0 |      0 | REML AI algorithm_genotype               | yes          |
|  1 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | -0.200843  |             168684 |               0 |      0 | REML AI algorithm_genotype_covariate     | yes          |
|  2 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 |  0.0943977 |             168684 |               0 |      0 | REML AI algorithm_genotype_covariate_pca | yes          |



Phenotype:  asthma
|    |   clump_p1 |   clump_r2 |   clump_kb |   p_window_size |   p_slide_size |   p_LD_threshold |   numberofpca |        h2 |   numberofvariants |   relatedmatrix |   data | h2model                                  | clumpprune   |
|---:|-----------:|-----------:|-----------:|----------------:|---------------:|-----------------:|--------------:|----------:|-------------------:|----------------:|-------:|:-----------------------------------------|:-------------|
|  0 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.0410093 |            27198.2 |               0 |      0 | REML AI algorithm_genotype               | yes          |
|  1 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.0407034 |            27198.2 |               0 |      0 | REML AI algorithm_genotype_covariate     | yes          |
|  2 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.0414508 |            27198.2 |               0 |      0 | REML AI algorithm_genotype_covariate_pca | yes          |



Phenotype:  osteoarthritis
|    |   clump_p1 |   clump_r2 |   clump_kb |   p_window_size |   p_slide_size |   p_LD_threshold |   numberofpca |          h2 |   numberofvariants |   relatedmatrix |   data | h2model                                  | clumpprune   |
|---:|-----------:|-----------:|-----------:|----------------:|---------------:|-----------------:|--------------:|------------:|-------------------:|----------------:|-------:|:-----------------------------------------|:-------------|
|  0 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.000786219 |              947.4 |               0 |      0 | REML AI algorithm_genotype               | yes          |
|  1 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.000786059 |              947.4 |               0 |      0 | REML AI algorithm_genotype_covariate     | yes          |
|  2 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.00079394  |              947.4 |               0 |      0 | REML AI algorithm_genotype_covariate_pca | yes          |



Phenotype:  high_cholesterol
|    |   clump_p1 |   clump_r2 |   clump_kb |   p_window_size |   p_slide_size |   p_LD_threshold |   numberofpca |         h2 |   numberofvariants |   relatedmatrix |   data | h2model                                  | clumpprune   |
|---:|-----------:|-----------:|-----------:|----------------:|---------------:|-----------------:|--------------:|-----------:|-------------------:|----------------:|-------:|:-----------------------------------------|:-------------|
|  0 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | -0.259972  |              70347 |               0 |      0 | REML AI algorithm_genotype               | yes          |
|  1 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | -0.125289  |             175896 |               0 |      0 | REML AI algorithm_genotype_covariate     | yes          |
|  2 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | -0.0770651 |             175896 |               0 |      0 | REML AI algorithm_genotype_covariate_pca | yes          |



Phenotype:  irritable_bowel_syndrome
|    |   clump_p1 |   clump_r2 |   clump_kb |   p_window_size |   p_slide_size |   p_LD_threshold |   numberofpca |        h2 |   numberofvariants |   relatedmatrix |   data | h2model                                  | clumpprune   |
|---:|-----------:|-----------:|-----------:|----------------:|---------------:|-----------------:|--------------:|----------:|-------------------:|----------------:|-------:|:-----------------------------------------|:-------------|
|  0 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 |  0.123336 |             112484 |               0 |      0 | REML AI algorithm_genotype               | yes          |
|  1 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | -0.721122 |             140263 |               0 |      0 | REML AI algorithm_genotype_covariate     | yes          |
|  2 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | -0.296956 |             175387 |               0 |      0 | REML AI algorithm_genotype_covariate_pca | yes          |



Phenotype:  hypothyroidism_myxoedema
|    |   clump_p1 |   clump_r2 |   clump_kb |   p_window_size |   p_slide_size |   p_LD_threshold |   numberofpca |        h2 |   numberofvariants |   relatedmatrix |   data | h2model                                  | clumpprune   |
|---:|-----------:|-----------:|-----------:|----------------:|---------------:|-----------------:|--------------:|----------:|-------------------:|----------------:|-------:|:-----------------------------------------|:-------------|
|  0 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.0883189 |            27908.8 |               0 |      0 | REML AI algorithm_genotype               | yes          |
|  1 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.0877105 |            27908.8 |               0 |      0 | REML AI algorithm_genotype_covariate     | yes          |
|  2 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.0894546 |            27908.8 |               0 |      0 | REML AI algorithm_genotype_covariate_pca | yes          |



Phenotype:  gastro_oesophageal_reflux_gord_gastric_reflux
|    |   clump_p1 |   clump_r2 |   clump_kb |   p_window_size |   p_slide_size |   p_LD_threshold |   numberofpca |        h2 |   numberofvariants |   relatedmatrix |   data | h2model                                  | clumpprune   |
|---:|-----------:|-----------:|-----------:|----------------:|---------------:|-----------------:|--------------:|----------:|-------------------:|----------------:|-------:|:-----------------------------------------|:-------------|
|  0 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.0200909 |              27866 |               0 |      0 | REML AI algorithm_genotype               | yes          |
|  1 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.0199449 |              27866 |               0 |      0 | REML AI algorithm_genotype_covariate     | yes          |
|  2 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.0203122 |              27866 |               0 |      0 | REML AI algorithm_genotype_covariate_pca | yes          |



Phenotype:  migraine
|    |   clump_p1 |   clump_r2 |   clump_kb |   p_window_size |   p_slide_size |   p_LD_threshold |   numberofpca |        h2 |   numberofvariants |   relatedmatrix |   data | h2model                                  | clumpprune   |
|---:|-----------:|-----------:|-----------:|----------------:|---------------:|-----------------:|--------------:|----------:|-------------------:|----------------:|-------:|:-----------------------------------------|:-------------|
|  0 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.0940761 |             169812 |               0 |      0 | REML AI algorithm_genotype               | yes          |
|  1 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.0979021 |             169812 |               0 |      0 | REML AI algorithm_genotype_covariate     | yes          |
|  2 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.10044   |             169812 |               0 |      0 | REML AI algorithm_genotype_covariate_pca | yes          |
