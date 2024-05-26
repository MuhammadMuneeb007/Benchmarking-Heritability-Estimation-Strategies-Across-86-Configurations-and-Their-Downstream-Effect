These tables summarize various parameters and metrics related to clumping and pruning, as well as heritability estimates and model details for different scenarios.  

| Column            | Description                                                                                     |
|-------------------|-------------------------------------------------------------------------------------------------|
| clump_p1          | The proportion of SNPs to retain during clumping.                                               |
| clump_r2          | The LD threshold for clumping.                                                                  |
| clump_kb          | The physical distance threshold (in kilobases) for clumping.                                    |
| p_window_size     | The size of the sliding window used during pruning.                                              |
| p_slide_size      | The amount the sliding window moves during pruning.                                             |
| p_LD_threshold    | The LD threshold used during pruning.                                                           |
| numberofpca       | The number of principal components used.                                                        |
| h2                | The estimated heritability.                                                                     |
| numberofvariants  | The number of variants used in the analysis.                                                    |
| h2model           | The model used for estimating heritability (e.g., LDpred-2_full, LDpred-2_hapmap).               |
| clumpprune        | Whether clumping and pruning were applied (yes/no).                                             |



Phenotype:  body_mass_index_bmi
|    |   clump_p1 |   clump_r2 |   clump_kb |   p_window_size |   p_slide_size |   p_LD_threshold |   numberofpca |        h2 |   numberofvariants | h2model         | clumpprune   |
|---:|-----------:|-----------:|-----------:|----------------:|---------------:|-----------------:|--------------:|----------:|-------------------:|:----------------|:-------------|
|  0 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.50282   |            28607.8 | LDpred-2_full   | yes          |
|  1 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.0338087 |            86190   | LDpred-2_full   | no           |
|  2 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.206038  |            13432.6 | LDpred-2_hapmap | yes          |
|  3 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.101754  |            40257   | LDpred-2_hapmap | no           |



Phenotype:  hypertension
|    |   clump_p1 |   clump_r2 |   clump_kb |   p_window_size |   p_slide_size |   p_LD_threshold |   numberofpca |          h2 |   numberofvariants | h2model         | clumpprune   |
|---:|-----------:|-----------:|-----------:|----------------:|---------------:|-----------------:|--------------:|------------:|-------------------:|:----------------|:-------------|
|  0 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 |  0.54183    |            26594.4 | LDpred-2_full   | yes          |
|  1 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | -0.00434896 |            81475   | LDpred-2_full   | no           |
|  2 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 |  0.145374   |            12614.2 | LDpred-2_hapmap | yes          |
|  3 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 |  0.0576817  |            38277   | LDpred-2_hapmap | no           |



Phenotype:  depression
|    |   clump_p1 |   clump_r2 |   clump_kb |   p_window_size |   p_slide_size |   p_LD_threshold |   numberofpca |         h2 |   numberofvariants | h2model       | clumpprune   |
|---:|-----------:|-----------:|-----------:|----------------:|---------------:|-----------------:|--------------:|-----------:|-------------------:|:--------------|:-------------|
|  0 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.124767   |             150146 | LDpred-2_full | yes          |
|  1 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.00503778 |             551925 | LDpred-2_full | no           |



Phenotype:  asthma
|    |   clump_p1 |   clump_r2 |   clump_kb |   p_window_size |   p_slide_size |   p_LD_threshold |   numberofpca |        h2 |   numberofvariants | h2model         | clumpprune   |
|---:|-----------:|-----------:|-----------:|----------------:|---------------:|-----------------:|--------------:|----------:|-------------------:|:----------------|:-------------|
|  0 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.0801001 |            28052.6 | LDpred-2_full   | yes          |
|  1 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.0801001 |            28052.6 | LDpred-2_full   | no           |
|  2 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.0488669 |             6828   | LDpred-2_hapmap | yes          |
|  3 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.0488669 |             6828   | LDpred-2_hapmap | no           |



Phenotype:  osteoarthritis
|    |   clump_p1 |   clump_r2 |   clump_kb |   p_window_size |   p_slide_size |   p_LD_threshold |   numberofpca |           h2 |   numberofvariants | h2model         | clumpprune   |
|---:|-----------:|-----------:|-----------:|----------------:|---------------:|-----------------:|--------------:|-------------:|-------------------:|:----------------|:-------------|
|  0 |        1   |       0.1  |        200 |             200 |             50 |             0.25 |           6   | -0.0117809   |              976.6 | LDpred-2_full   | yes          |
|  1 |        1   |       0.1  |        200 |             200 |             50 |             0.25 |           6   | -3.90756e-05 |             1294   | LDpred-2_full   | no           |
|  2 |        0.4 |       0.04 |         80 |              80 |             20 |             0.1  |           2.4 |  4.26297e-05 |              218.8 | LDpred-2_hapmap | yes          |
|  3 |        0.4 |       0.04 |         80 |              80 |             20 |             0.1  |           2.4 |  3.76102e-06 |              254.4 | LDpred-2_hapmap | no           |



Phenotype:  high_cholesterol
|    |   clump_p1 |   clump_r2 |   clump_kb |   p_window_size |   p_slide_size |   p_LD_threshold |   numberofpca |        h2 |   numberofvariants | h2model       | clumpprune   |
|---:|-----------:|-----------:|-----------:|----------------:|---------------:|-----------------:|--------------:|----------:|-------------------:|:--------------|:-------------|
|  0 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.0268799 |             181529 | LDpred-2_full | yes          |
|  1 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.0398526 |             575564 | LDpred-2_full | no           |



Phenotype:  irritable_bowel_syndrome
|    |   clump_p1 |   clump_r2 |   clump_kb |   p_window_size |   p_slide_size |   p_LD_threshold |   numberofpca |           h2 |   numberofvariants | h2model       | clumpprune   |
|---:|-----------:|-----------:|-----------:|----------------:|---------------:|-----------------:|--------------:|-------------:|-------------------:|:--------------|:-------------|
|  0 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 |  0.0195793   |             181200 | LDpred-2_full | yes          |
|  1 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | -0.000515425 |             576983 | LDpred-2_full | no           |



Phenotype:  hypothyroidism_myxoedema
|    |   clump_p1 |   clump_r2 |   clump_kb |   p_window_size |   p_slide_size |   p_LD_threshold |   numberofpca |        h2 |   numberofvariants | h2model         | clumpprune   |
|---:|-----------:|-----------:|-----------:|----------------:|---------------:|-----------------:|--------------:|----------:|-------------------:|:----------------|:-------------|
|  0 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.146424  |            28790.2 | LDpred-2_full   | yes          |
|  1 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.0517344 |            86545   | LDpred-2_full   | no           |
|  2 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.0751322 |            13464.6 | LDpred-2_hapmap | yes          |
|  3 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.0276094 |            40322   | LDpred-2_hapmap | no           |



Phenotype:  gastro_oesophageal_reflux_gord_gastric_reflux
|    |   clump_p1 |   clump_r2 |   clump_kb |   p_window_size |   p_slide_size |   p_LD_threshold |   numberofpca |          h2 |   numberofvariants | h2model         | clumpprune   |
|---:|-----------:|-----------:|-----------:|----------------:|---------------:|-----------------:|--------------:|------------:|-------------------:|:----------------|:-------------|
|  0 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.0625676   |            28778.8 | LDpred-2_full   | yes          |
|  1 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.000330644 |            86545   | LDpred-2_full   | no           |
|  2 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.0340128   |            13473.8 | LDpred-2_hapmap | yes          |
|  3 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | 0.00785986  |            40322   | LDpred-2_hapmap | no           |



Phenotype:  migraine
|    |   clump_p1 |   clump_r2 |   clump_kb |   p_window_size |   p_slide_size |   p_LD_threshold |   numberofpca |           h2 |   numberofvariants | h2model       | clumpprune   |
|---:|-----------:|-----------:|-----------:|----------------:|---------------:|-----------------:|--------------:|-------------:|-------------------:|:--------------|:-------------|
|  0 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | -0.00706787  |             148978 | LDpred-2_full | yes          |
|  1 |          1 |        0.1 |        200 |             200 |             50 |             0.25 |             6 | -0.000786406 |             555483 | LDpred-2_full | no           |


