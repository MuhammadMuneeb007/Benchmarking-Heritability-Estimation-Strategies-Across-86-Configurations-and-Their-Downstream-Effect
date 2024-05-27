
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
|    |   clump_p1 |   clump_r2 |   clump_kb |   p_window_size |   p_slide_size |   p_LD_threshold |   numberofpca |       h2 |   numberofvariants |   relatedmatrix |   data | h2model                              | clumpprune   |
|---:|-----------:|-----------:|-----------:|----------------:|---------------:|-----------------:|--------------:|---------:|-------------------:|----------------:|-------:|:-------------------------------------|:-------------|
|  0 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.155912 |              27515 |               0 |      0 | HE regression_genotype               | yes          |
|  1 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.154453 |              27515 |               0 |      0 | HE regression_genotype_covariate     | yes          |
|  2 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.157302 |              27515 |               0 |      0 | HE regression_genotype_covariate_pca | yes          |
|  3 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.129021 |              82362 |               0 |      0 | HE regression_genotype               | no           |
|  4 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.128277 |              82362 |               0 |      0 | HE regression_genotype_covariate     | no           |
|  5 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.129763 |              82362 |               0 |      0 | HE regression_genotype_covariate_pca | no           |



Phenotype:  hypertension
|    |   clump_p1 |   clump_r2 |   clump_kb |   p_window_size |   p_slide_size |   p_LD_threshold |   numberofpca |        h2 |   numberofvariants |   relatedmatrix |   data | h2model                              | clumpprune   |
|---:|-----------:|-----------:|-----------:|----------------:|---------------:|-----------------:|--------------:|----------:|-------------------:|----------------:|-------:|:-------------------------------------|:-------------|
|  0 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.127176  |            25824.2 |               0 |      0 | HE regression_genotype               | yes          |
|  1 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.126026  |            25824.2 |               0 |      0 | HE regression_genotype_covariate     | yes          |
|  2 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.12838   |            25824.2 |               0 |      0 | HE regression_genotype_covariate_pca | yes          |
|  3 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.0609309 |            78742.2 |               0 |      0 | HE regression_genotype               | no           |
|  4 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.0606503 |            78742.2 |               0 |      0 | HE regression_genotype_covariate     | no           |
|  5 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.0612456 |            78742.2 |               0 |      0 | HE regression_genotype_covariate_pca | no           |



Phenotype:  depression
|    |   clump_p1 |   clump_r2 |   clump_kb |   p_window_size |   p_slide_size |   p_LD_threshold |   numberofpca |        h2 |   numberofvariants |   relatedmatrix |   data | h2model                              | clumpprune   |
|---:|-----------:|-----------:|-----------:|----------------:|---------------:|-----------------:|--------------:|----------:|-------------------:|----------------:|-------:|:-------------------------------------|:-------------|
|  0 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.155898  |             168684 |               0 |      0 | HE regression_genotype               | yes          |
|  1 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.152228  |             168684 |               0 |      0 | HE regression_genotype_covariate     | yes          |
|  2 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.155135  |             168684 |               0 |      0 | HE regression_genotype_covariate_pca | yes          |
|  3 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.081918  |             531018 |               0 |      0 | HE regression_genotype               | no           |
|  4 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.0808639 |             531018 |               0 |      0 | HE regression_genotype_covariate     | no           |
|  5 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.0817495 |             531018 |               0 |      0 | HE regression_genotype_covariate_pca | no           |



Phenotype:  asthma
|    |   clump_p1 |   clump_r2 |   clump_kb |   p_window_size |   p_slide_size |   p_LD_threshold |   numberofpca |        h2 |   numberofvariants |   relatedmatrix |   data | h2model                              | clumpprune   |
|---:|-----------:|-----------:|-----------:|----------------:|---------------:|-----------------:|--------------:|----------:|-------------------:|----------------:|-------:|:-------------------------------------|:-------------|
|  0 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.0359066 |            27198.2 |               0 |      0 | HE regression_genotype               | yes          |
|  1 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.0356364 |            27198.2 |               0 |      0 | HE regression_genotype_covariate     | yes          |
|  2 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.036281  |            27198.2 |               0 |      0 | HE regression_genotype_covariate_pca | yes          |
|  3 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.0260743 |            81159.8 |               0 |      0 | HE regression_genotype               | no           |
|  4 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.0260693 |            81159.8 |               0 |      0 | HE regression_genotype_covariate     | no           |
|  5 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.0263832 |            81159.8 |               0 |      0 | HE regression_genotype_covariate_pca | no           |



Phenotype:  osteoarthritis
|    |   clump_p1 |   clump_r2 |   clump_kb |   p_window_size |   p_slide_size |   p_LD_threshold |   numberofpca |          h2 |   numberofvariants |   relatedmatrix |   data | h2model                              | clumpprune   |
|---:|-----------:|-----------:|-----------:|----------------:|---------------:|-----------------:|--------------:|------------:|-------------------:|----------------:|-------:|:-------------------------------------|:-------------|
|  0 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.000797282 |              947.4 |               0 |      0 | HE regression_genotype               | yes          |
|  1 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.000797158 |              947.4 |               0 |      0 | HE regression_genotype_covariate     | yes          |
|  2 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.000805351 |              947.4 |               0 |      0 | HE regression_genotype_covariate_pca | yes          |
|  3 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.000821994 |             1245.8 |               0 |      0 | HE regression_genotype               | no           |
|  4 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.000822791 |             1245.8 |               0 |      0 | HE regression_genotype_covariate     | no           |
|  5 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.000828757 |             1245.8 |               0 |      0 | HE regression_genotype_covariate_pca | no           |



Phenotype:  high_cholesterol
|    |   clump_p1 |   clump_r2 |   clump_kb |   p_window_size |   p_slide_size |   p_LD_threshold |   numberofpca |        h2 |   numberofvariants |   relatedmatrix |   data | h2model                              | clumpprune   |
|---:|-----------:|-----------:|-----------:|----------------:|---------------:|-----------------:|--------------:|----------:|-------------------:|----------------:|-------:|:-------------------------------------|:-------------|
|  0 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.125216  |             175896 |               0 |      0 | HE regression_genotype               | yes          |
|  1 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.122009  |             175896 |               0 |      0 | HE regression_genotype_covariate     | yes          |
|  2 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.124326  |             175896 |               0 |      0 | HE regression_genotype_covariate_pca | yes          |
|  3 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.0674643 |             553320 |               0 |      0 | HE regression_genotype               | no           |
|  4 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.0664132 |             553320 |               0 |      0 | HE regression_genotype_covariate     | no           |
|  5 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.0672279 |             553320 |               0 |      0 | HE regression_genotype_covariate_pca | no           |



Phenotype:  irritable_bowel_syndrome
|    |   clump_p1 |   clump_r2 |   clump_kb |   p_window_size |   p_slide_size |   p_LD_threshold |   numberofpca |        h2 |   numberofvariants |   relatedmatrix |   data | h2model                              | clumpprune   |
|---:|-----------:|-----------:|-----------:|----------------:|---------------:|-----------------:|--------------:|----------:|-------------------:|----------------:|-------:|:-------------------------------------|:-------------|
|  0 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.0866639 |             175387 |               0 |      0 | HE regression_genotype               | yes          |
|  1 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.0846073 |             175387 |               0 |      0 | HE regression_genotype_covariate     | yes          |
|  2 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.0861963 |             175387 |               0 |      0 | HE regression_genotype_covariate_pca | yes          |
|  3 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.0295101 |             552363 |               0 |      0 | HE regression_genotype               | no           |
|  4 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.0290978 |             552363 |               0 |      0 | HE regression_genotype_covariate     | no           |
|  5 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.0295235 |             552363 |               0 |      0 | HE regression_genotype_covariate_pca | no           |



Phenotype:  hypothyroidism_myxoedema
|    |   clump_p1 |   clump_r2 |   clump_kb |   p_window_size |   p_slide_size |   p_LD_threshold |   numberofpca |        h2 |   numberofvariants |   relatedmatrix |   data | h2model                              | clumpprune   |
|---:|-----------:|-----------:|-----------:|----------------:|---------------:|-----------------:|--------------:|----------:|-------------------:|----------------:|-------:|:-------------------------------------|:-------------|
|  0 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.0749919 |            27908.8 |               0 |      0 | HE regression_genotype               | yes          |
|  1 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.0743282 |            27908.8 |               0 |      0 | HE regression_genotype_covariate     | yes          |
|  2 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.0757672 |            27908.8 |               0 |      0 | HE regression_genotype_covariate_pca | yes          |
|  3 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.0601052 |            83411.6 |               0 |      0 | HE regression_genotype               | no           |
|  4 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.0598492 |            83411.6 |               0 |      0 | HE regression_genotype_covariate     | no           |
|  5 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.0605275 |            83411.6 |               0 |      0 | HE regression_genotype_covariate_pca | no           |



Phenotype:  gastro_oesophageal_reflux_gord_gastric_reflux
|    |   clump_p1 |   clump_r2 |   clump_kb |   p_window_size |   p_slide_size |   p_LD_threshold |   numberofpca |         h2 |   numberofvariants |   relatedmatrix |   data | h2model                              | clumpprune   |
|---:|-----------:|-----------:|-----------:|----------------:|---------------:|-----------------:|--------------:|-----------:|-------------------:|----------------:|-------:|:-------------------------------------|:-------------|
|  0 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.0184918  |            27866   |               0 |      0 | HE regression_genotype               | yes          |
|  1 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.0183522  |            27866   |               0 |      0 | HE regression_genotype_covariate     | yes          |
|  2 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.0186843  |            27866   |               0 |      0 | HE regression_genotype_covariate_pca | yes          |
|  3 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.00943993 |            83369.2 |               0 |      0 | HE regression_genotype               | no           |
|  4 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.00939441 |            83369.2 |               0 |      0 | HE regression_genotype_covariate     | no           |
|  5 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.00951608 |            83369.2 |               0 |      0 | HE regression_genotype_covariate_pca | no           |



Phenotype:  migraine
|    |   clump_p1 |   clump_r2 |   clump_kb |   p_window_size |   p_slide_size |   p_LD_threshold |   numberofpca |         h2 |   numberofvariants |   relatedmatrix |   data | h2model                              | clumpprune   |
|---:|-----------:|-----------:|-----------:|----------------:|---------------:|-----------------:|--------------:|-----------:|-------------------:|----------------:|-------:|:-------------------------------------|:-------------|
|  0 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.0654348  |             169812 |               0 |      0 | HE regression_genotype               | yes          |
|  1 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.0638261  |             169812 |               0 |      0 | HE regression_genotype_covariate     | yes          |
|  2 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.0650382  |             169812 |               0 |      0 | HE regression_genotype_covariate_pca | yes          |
|  3 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.001976   |             535177 |               0 |      0 | HE regression_genotype               | no           |
|  4 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.00194932 |             535177 |               0 |      0 | HE regression_genotype_covariate     | no           |
|  5 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.00197048 |             535177 |               0 |      0 | HE regression_genotype_covariate_pca | no           |

