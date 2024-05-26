
These tables present various parameters and metrics related to clumping and pruning, as well as heritability estimates and model details for a specific scenario. 

| Parameter          | Description                                             |
|--------------------|---------------------------------------------------------|
| clump_p1           | Proportion of SNPs retained during clumping.            |
| clump_r2           | LD threshold for clumping.                              |
| clump_kb           | Physical distance threshold (in kilobases) for clumping.|
| p_window_size      | Size of the sliding window used during pruning.         |
| p_slide_size       | Amount the sliding window moves during pruning.         |
| p_LD_threshold     | LD threshold used during pruning.                       |
| numberofpca        | Number of principal components used.                    |
| h2                 | Estimated heritability.                                 |
| numberofvariants   | Number of variants used in the analysis.                |
| h2model            | Model used for estimating heritability (e.g., "GCTA_genotype"). |
| clumpprune         | Whether clumping and pruning were applied (yes/no).     |


Phenotype:  body_mass_index_bmi
|    |   clump_p1 |   clump_r2 |   clump_kb |   p_window_size |   p_slide_size |   p_LD_threshold |   numberofpca |        h2 |   numberofvariants | h2model                     | clumpprune   |
|---:|-----------:|-----------:|-----------:|----------------:|---------------:|-----------------:|--------------:|----------:|-------------------:|:----------------------------|:-------------|
|  0 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.393054  |            28606.8 | GCTA_genotype               | yes          |
|  1 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.659147  |           619652   | GCTA_genotype               | no           |
|  2 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.0038016 |            28606.8 | GCTA_genotype_covariate     | yes          |
|  3 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.330665  |           619652   | GCTA_genotype_covariate     | no           |
|  4 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.0463022 |            28606.8 | GCTA_genotype_covariate_pca | yes          |
|  5 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.36103   |           619652   | GCTA_genotype_covariate_pca | no           |



Phenotype:  hypertension
|    |   clump_p1 |   clump_r2 |   clump_kb |   p_window_size |   p_slide_size |   p_LD_threshold |   numberofpca |        h2 |   numberofvariants | h2model                     | clumpprune   |
|---:|-----------:|-----------:|-----------:|----------------:|---------------:|-----------------:|--------------:|----------:|-------------------:|:----------------------------|:-------------|
|  0 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.0713942 |            26593.4 | GCTA_genotype               | yes          |
|  1 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.805235  |           619652   | GCTA_genotype               | no           |
|  2 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 2e-06     |            26593.4 | GCTA_genotype_covariate     | yes          |
|  3 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.759073  |           619652   | GCTA_genotype_covariate     | no           |
|  4 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 2e-06     |            26593.4 | GCTA_genotype_covariate_pca | yes          |
|  5 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.761765  |           619652   | GCTA_genotype_covariate_pca | no           |



Phenotype:  depression
|    |   clump_p1 |   clump_r2 |   clump_kb |   p_window_size |   p_slide_size |   p_LD_threshold |   numberofpca |       h2 |   numberofvariants | h2model                     | clumpprune   |
|---:|-----------:|-----------:|-----------:|----------------:|---------------:|-----------------:|--------------:|---------:|-------------------:|:----------------------------|:-------------|
|  0 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.033539 |             174078 | GCTA_genotype               | yes          |
|  1 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 1e-06    |             619652 | GCTA_genotype               | no           |
|  2 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.291195 |             174078 | GCTA_genotype_covariate     | yes          |
|  3 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.285443 |             619652 | GCTA_genotype_covariate     | no           |
|  4 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.248543 |             174078 | GCTA_genotype_covariate_pca | yes          |
|  5 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.283982 |             619652 | GCTA_genotype_covariate_pca | no           |



Phenotype:  asthma
|    |   clump_p1 |   clump_r2 |   clump_kb |   p_window_size |   p_slide_size |   p_LD_threshold |   numberofpca |       h2 |   numberofvariants | h2model                     | clumpprune   |
|---:|-----------:|-----------:|-----------:|----------------:|---------------:|-----------------:|--------------:|---------:|-------------------:|:----------------------------|:-------------|
|  0 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 1e-06    |              28076 | GCTA_genotype               | yes          |
|  1 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.201014 |             619652 | GCTA_genotype               | no           |
|  2 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 1e-06    |              28076 | GCTA_genotype_covariate     | yes          |
|  3 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.54633  |             619652 | GCTA_genotype_covariate     | no           |
|  4 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 1e-06    |              28076 | GCTA_genotype_covariate_pca | yes          |
|  5 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.628824 |             619652 | GCTA_genotype_covariate_pca | no           |



Phenotype:  osteoarthritis
|    |   clump_p1 |   clump_r2 |   clump_kb |   p_window_size |   p_slide_size |   p_LD_threshold |   numberofpca |       h2 |   numberofvariants | h2model                     | clumpprune   |
|---:|-----------:|-----------:|-----------:|----------------:|---------------:|-----------------:|--------------:|---------:|-------------------:|:----------------------------|:-------------|
|  0 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.130503 |              975.6 | GCTA_genotype               | yes          |
|  1 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 1e-06    |           619652   | GCTA_genotype               | no           |
|  2 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.177943 |              975.6 | GCTA_genotype_covariate     | yes          |
|  3 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.544966 |           619652   | GCTA_genotype_covariate     | no           |
|  4 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.185921 |              975.6 | GCTA_genotype_covariate_pca | yes          |
|  5 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.509423 |           619652   | GCTA_genotype_covariate_pca | no           |



Phenotype:  high_cholesterol
|    |   clump_p1 |   clump_r2 |   clump_kb |   p_window_size |   p_slide_size |   p_LD_threshold |   numberofpca |        h2 |   numberofvariants | h2model                     | clumpprune   |
|---:|-----------:|-----------:|-----------:|----------------:|---------------:|-----------------:|--------------:|----------:|-------------------:|:----------------------------|:-------------|
|  0 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.401021  |             145373 | GCTA_genotype               | yes          |
|  1 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.308373  |             619652 | GCTA_genotype               | no           |
|  2 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.618465  |             181748 | GCTA_genotype_covariate     | yes          |
|  3 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.0585888 |             619652 | GCTA_genotype_covariate     | no           |
|  4 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.57904   |             181748 | GCTA_genotype_covariate_pca | yes          |
|  5 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.0524474 |             619652 | GCTA_genotype_covariate_pca | no           |



Phenotype:  irritable_bowel_syndrome
|    |   clump_p1 |   clump_r2 |   clump_kb |   p_window_size |   p_slide_size |   p_LD_threshold |   numberofpca |       h2 |   numberofvariants | h2model                     | clumpprune   |
|---:|-----------:|-----------:|-----------:|----------------:|---------------:|-----------------:|--------------:|---------:|-------------------:|:----------------------------|:-------------|
|  0 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.295308 |             181402 | GCTA_genotype               | yes          |
|  1 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.393782 |             619652 | GCTA_genotype               | no           |
|  2 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 1e-06    |             150508 | GCTA_genotype_covariate     | yes          |
|  3 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.170137 |             619652 | GCTA_genotype_covariate     | no           |
|  4 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.052254 |             181402 | GCTA_genotype_covariate_pca | yes          |
|  5 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.389177 |             619652 | GCTA_genotype_covariate_pca | no           |



Phenotype:  hypothyroidism_myxoedema
|    |   clump_p1 |   clump_r2 |   clump_kb |   p_window_size |   p_slide_size |   p_LD_threshold |   numberofpca |       h2 |   numberofvariants | h2model                     | clumpprune   |
|---:|-----------:|-----------:|-----------:|----------------:|---------------:|-----------------:|--------------:|---------:|-------------------:|:----------------------------|:-------------|
|  0 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.16943  |            28789.2 | GCTA_genotype               | yes          |
|  1 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.414696 |           619652   | GCTA_genotype               | no           |
|  2 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.426782 |            28789.2 | GCTA_genotype_covariate     | yes          |
|  3 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.597121 |           619652   | GCTA_genotype_covariate     | no           |
|  4 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.413135 |            28789.2 | GCTA_genotype_covariate_pca | yes          |
|  5 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.645272 |           619652   | GCTA_genotype_covariate_pca | no           |



Phenotype:  gastro_oesophageal_reflux_gord_gastric_reflux
|    |   clump_p1 |   clump_r2 |   clump_kb |   p_window_size |   p_slide_size |   p_LD_threshold |   numberofpca |        h2 |   numberofvariants | h2model                     | clumpprune   |
|---:|-----------:|-----------:|-----------:|----------------:|---------------:|-----------------:|--------------:|----------:|-------------------:|:----------------------------|:-------------|
|  0 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.0776772 |            28777.8 | GCTA_genotype               | yes          |
|  1 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.35827   |           619652   | GCTA_genotype               | no           |
|  2 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.0038396 |            28777.8 | GCTA_genotype_covariate     | yes          |
|  3 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.349698  |           619652   | GCTA_genotype_covariate     | no           |
|  4 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 1e-06     |            28777.8 | GCTA_genotype_covariate_pca | yes          |
|  5 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.293561  |           619652   | GCTA_genotype_covariate_pca | no           |



Phenotype:  migraine
|    |   clump_p1 |   clump_r2 |   clump_kb |   p_window_size |   p_slide_size |   p_LD_threshold |   numberofpca |        h2 |   numberofvariants | h2model                     | clumpprune   |
|---:|-----------:|-----------:|-----------:|----------------:|---------------:|-----------------:|--------------:|----------:|-------------------:|:----------------------------|:-------------|
|  0 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.55832   |             149213 | GCTA_genotype               | yes          |
|  1 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.0882328 |             619652 | GCTA_genotype               | no           |
|  2 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.481936  |             174870 | GCTA_genotype_covariate     | yes          |
|  3 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.236098  |             619652 | GCTA_genotype_covariate     | no           |
|  4 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.507982  |             141274 | GCTA_genotype_covariate_pca | yes          |
|  5 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.186043  |             619652 | GCTA_genotype_covariate_pca | no           |

