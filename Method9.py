#!/usr/bin/env python
# coding: utf-8

# # LDpred-2 
# 
# In this notebook, we will use LDpred-2 to calculate the Polygenic Risk Score (PRS).
# 
# **Note:** LDpred-2 needs to be installed using R.
# 
# We will use the same flow of the code, and when we have to invoke LDpred-2 functions, we will call the script using Python. When LDpred-2 finishes executions, results will be stored in the specific directories for each fold and will be retrieved in Python for further calculations.
# 
# ## Install LDpred-2
# 
# Install LDpred-2 using the following commands:
# 
# 1. Activate the conda Environment.
# 2. Type R
# 3. Run the following commands:
# 
# ```R
# install.packages("remotes")
# library(remotes)
# remotes::install_github("https://github.com/privefl/bigsnpr.git")
# ```
# 
# The content in this notebook has been taken from the [PRS Tutorial](https://choishingwan.github.io/PRS-Tutorial/ldpred/), [LDpred-2 Documentation](https://privefl.github.io/bigsnpr/articles/LDpred2.html), and [Polygenic Scores (PGS)](https://privefl.github.io/bigsnpr-extdoc/polygenic-scores-pgs.html).
# 
# Thanks to Florian Privé for creating one of the best documentations explaining each aspect of the calculation. We will be using LDpred-2 for cross-validation design and hyper-parameter optimization, and this tutorial complements the existing one.
# 
# 

# ## LDpred-2  Hyperparameters
# 
# ### Hyperparameters for LDpred-2  performed using PLINK
# 
# #### Pruning Parameters
# 
# Informs Plink that we wish to perform pruning with a window size of 200 variants, sliding across the genome with a step size of 50 variants at a time, and filter out any SNPs with LD \( r^2 \) higher than 0.25.
# 
# 
# ```python
# 1. p_window_size = [200]
# 2. p_slide_size = [50]
# 3. p_LD_threshold = [0.25]
# ```
# 
# #### Clumping Parameters
# 
# The P-value threshold for an SNP to be included. 1 means to include all SNPs for clumping. SNPs having \( r^2 \) higher than 0.1 with the index SNPs will be removed. SNPs within 200k of the index SNP are considered for clumping.
# 
# ```python
# 1. clump_p1 = [1]
# 2. clump_r2 = [0.1]
# 3. clump_kb = [200]
# ```
# #### PCA
# Pca also affects the results evident from the initial analysis; however, including more PCA overfits the model. 
# 
# ### Hyperparameters for LDpred-2 
#  
# 
# 
# #### Heritability
# 
# To calculate h2 using LDpred-2, we followed the [LDpred-2 tutorial](https://choishingwan.github.io/PRS-Tutorial/ldpred/).
# The code for this part is in R, as LDpred-2 is in R.
# 
# 1. Using SNPs from the HapMap as preferred by the authors. The name in the code is `LDpred-2_hapmap`.
# 2. Using all the SNPs. The name in the code is `LDpred-2_full`.
# 
# This approach is computationally expensive, as the correlation between all the SNPs in a specific chromosome is being calculated.
# 
# ## LDpred-2 Models
# 
# Each model has its own set of hyperparameters, so first, we will generate those hyperparameters for those models and call them when invoking LDpred-2.
# 
# Use the following R code to know more about the hyperparameters:
# - `help(snp_ldpred2_inf)`
# - `help(snp_ldpred2_grid)`
# - `help(snp_ldpred2_auto)`
# - `help(snp_lassosum2)`
# 
# ### 1. LDpred2-inf: Infinitesimal Model
# 
# No Hyperparameters
# 
# ### 2. LDpred2(-grid): Grid of Models
# 
# Each model has its own set of parameters:
# - `burn_in = 100`
# - `num_iter = 10`
# - `p = c(1e-4, 1e-5)`
# - `h2 = c(0.1, 0.2)`
# - `sparse = c(FALSE, TRUE)`
# 
# ### 3. LDpred2-auto: Automatic Model
# 
# - `burn_in = 100`
# - `num_iter = 100`
# - `sparse = c(FALSE, TRUE)`
# - `alpha_bounds = c(1, 0.5)`
# - `use_MLE = c(FALSE, TRUE)`
# - `p = c(1e-4, 1e-5)`
# - `shrink_corr = 0.7`
# - `allow_jump_sign = c(FALSE, TRUE)`
# 
# ### 4. lassosum2: Grid of Models
# 
# - `lambda = c(0.001, 0.01, 0.1, 1)`
# - `delta = 30`
# 
# We will execute the code for each mkodel seperatly.
# 

# ### GWAS file processing for Plink.
# When the effect size relates to disease risk and is thus given as an odds ratio (OR) rather than BETA (for continuous traits), the PRS is computed as a product of ORs. To simplify this calculation, take the natural logarithm of the OR so that the PRS can be computed using summation instead.

# In[ ]:

import os
import sys
import pandas as pd
import numpy as np

filedirec = sys.argv[1]

GWAS = filedirec + os.sep + filedirec+".gz"
df = pd.read_csv(GWAS,compression= "gzip",sep="\s+")

import pandas as pd
if "BETA" in df.columns.to_list():
    df = df[['CHR', 'BP', 'SNP', 'A1', 'A2', 'N', 'SE', 'P', 'BETA', 'INFO', 'MAF']]
    df['Z'] = df['BETA'] / df['SE']
    #df = df[['SNP','N','Z','A1', 'A2']]


else:
    df["BETA"] = np.log(df["OR"])
    df = df[['CHR', 'BP', 'SNP', 'A1', 'A2', 'N', 'SE', 'P', 'BETA', 'INFO', 'MAF']]
    df['Z'] = df['BETA'] / df['SE']
    #df = df[['SNP','N','Z','A1', 'A2']]

df.to_csv(filedirec + os.sep +filedirec+".txt",sep="\t",index=False)
import subprocess
df_transformed = pd.DataFrame({
    'Predictor': df['CHR'].astype(str) + ":" + df['BP'].astype(str),
    'A1': df['A1'],
    'A2': df['A2'],
    'n': df['N'],
    'Z': df['BETA']/df['SE'],
    'SNP':df['SNP']
}) 

columns_to_check = ['A1', 'A2']
 
def specified_single_characters(row, cols):
    return all(len(str(row[col])) == 1 for col in cols)
 
df_transformed = df_transformed[df_transformed.apply(lambda row: specified_single_characters(row, columns_to_check), axis=1)]

df_transformed.to_csv(filedirec + os.sep +filedirec+".ldak",sep="\t",index=False)


commands = [
    "wget https://genetics.ghpc.au.dk/doug/ldak.thin.hapmap.gbr.tagging.gz",
    "wget https://genetics.ghpc.au.dk/doug/bld.ldak.hapmap.gbr.tagging.gz",
    "wget https://genetics.ghpc.au.dk/doug/bld.ldak.lite.alpha.hapmap.gbr.tagging.gz",
    "gunzip ldak.thin.hapmap.gbr.tagging.gz",
    "gunzip bld.ldak.hapmap.gbr.tagging.gz",
    "gunzip bld.ldak.lite.alpha.hapmap.gbr.tagging.gz",
    "wget https://www.dropbox.com/s/o7xphugm4mln9xa/pow.txt"
]

# Execute each command
for command in commands:
    #subprocess.run(command, shell=True)
    pass



 

from operator import index
import pandas as pd
import numpy as np
import os
import subprocess
import sys
import pandas as pd
import statsmodels.api as sm
import pandas as pd
from sklearn.metrics import roc_auc_score, confusion_matrix
from statsmodels.stats.contingency_tables import mcnemar

def create_directory(directory):
    """Function to create a directory if it doesn't exist."""
    if not os.path.exists(directory):  # Checking if the directory doesn't exist
        os.makedirs(directory)  # Creating the directory if it doesn't exist
    return directory  # Returning the created or existing directory

filedirec = sys.argv[1]  # Setting the base directory path

foldnumber = sys.argv[2]
#foldnumber = "1"  # Setting 'foldnumber' to "0"

folddirec = filedirec + os.sep + "Fold_" + foldnumber  # Creating a directory path for the specific fold
trainfilename = "train_data"  # Setting the name of the training data file
newtrainfilename = "train_data.QC"  # Setting the name of the new training data file

testfilename = "test_data"  # Setting the name of the test data file
newtestfilename = "test_data.QC"  # Setting the name of the new test data file

# Number of PCA to be included as a covariate.
numberofpca = ["6"]  # Setting the number of PCA components to be included

# Clumping parameters.
clump_p1 = [1]  # List containing clump parameter 'p1'
clump_r2 = [0.1]  # List containing clump parameter 'r2'
clump_kb = [200]  # List containing clump parameter 'kb'

# Pruning parameters.
p_window_size = [200]  # List containing pruning parameter 'window_size'
p_slide_size = [50]  # List containing pruning parameter 'slide_size'
p_LD_threshold = [0.25]  # List containing pruning parameter 'LD_threshold'



# Initializing an empty DataFrame with specified column names
prs_result = pd.DataFrame(columns=["clump_p1", "clump_r2", "clump_kb", "p_window_size", "p_slide_size", "p_LD_threshold",
                                "numberofpca","h2model","h2","numberofvariants","clumpprune","relatedmatrix","precomputedfile"])

 
def transform_ldak_data(traindirec, newtrainfilename,p,precomputedfile,p1_val, p2_val, p3_val, c1_val, c2_val, c3_val,Name):
 
    def readpve(file):
        file_path = file+os.sep+ "snpher1.hers"  
        with open(file_path, 'r') as file:
            for line in file:
                if "Her_All" in line:
                    pve_estimate = line.split(" ")[1].strip()
                    print("Heritability", pve_estimate)
                    return pve_estimate

    def readvariants(file):
        file_path = file + os.sep+ "snpher1.overlap"  
        with open(file_path, 'r') as file:
            for line in file:
                if "Summary_Statistic_Predictors" in line:
                    pve_estimate = line.split(" ")[1].strip()
                    print("Variants", pve_estimate)
                    return pve_estimate

    #if precomputedfile == "bld.ldak.hapmap.gbr.tagging":
    command = [
        "./ldak",
        "--sum-hers", traindirec+os.sep+"snpher1",
        "--summary", filedirec + os.sep +filedirec+".ldak",
        "--tagfile", "LDAKFILES"+os.sep+precomputedfile,
        "--check-sums","NO"
    ]


    print(" ".join(command))
    subprocess.run(command)

     
    pve_estimate = readpve(traindirec)
    variants = readvariants(traindirec)



    global prs_result 
    prs_result = prs_result._append({
  
        "h2":pve_estimate,
        "h2model":"LDAK"+"_"+precomputedfile,
        "clumpprune":"no",
        "numberofvariants":variants,
        "numberofpca":"-"
         
    }, ignore_index=True)

    prs_result.to_csv(traindirec+os.sep+Name+os.sep+"Results.csv",index=False)
    
 
    return

 

precomputedfiles = ["bld.ldak.hapmap.gbr.tagging","ldak.thin.hapmap.gbr.tagging",
"bld.ldak.lite.alpha.hapmap.gbr.tagging",
"bld.ldak.genotyped.gbr.tagging","ldak.thin.genotyped.gbr.tagging",
"bld.ldak.lite.alpha.genotyped.gbr.tagging"]
 
clumpprunes = ["yes" ]
# Nested loops to iterate over different parameter values
create_directory(folddirec+os.sep+"LDAK-Precomputed-Heritability")
for p1_val in p_window_size:
 for p2_val in p_slide_size: 
  for p3_val in p_LD_threshold:
   for c1_val in clump_p1:
    for c2_val in clump_r2:
     for c3_val in clump_kb:
      for p in numberofpca:
          for precomputedfile  in precomputedfiles:
            transform_ldak_data(folddirec, newtrainfilename, p,precomputedfile,str(p1_val), str(p2_val), str(p3_val), str(c1_val), str(c2_val), str(c3_val), "LDAK-Precomputed-Heritability")
            #exit(0)
exit(0)




