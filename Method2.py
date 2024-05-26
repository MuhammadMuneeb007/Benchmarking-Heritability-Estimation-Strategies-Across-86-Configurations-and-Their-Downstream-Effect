
# 
# #### Calculate h2 Heritability
# 
# There are multiple ways to calculate heritability, and we considered two methods required for lambda = m * (1 / h2SNP - 1).
# 
# ##### Using LDpred-2
# 
# To calculate h2 using LDpred-2, we followed the [LDpred-2 tutorial](https://choishingwan.github.io/PRS-Tutorial/ldpred/).
# The code for this part is in R, as LDpred-2 is in R.
# 
# 1. Using SNPs from the HapMap as preferred by the authors. The name in the code is `LDpred-2_hapmap`.
# 2. Using all the SNPs. The name in the code is `LDpred-2_full`.
# 
# This approach is computationally expensive, as the correlation between all the SNPs in a specific chromosome is being calculated.
# 
# ##### Using GCTA
# 
# To calculate h2 using GCTA, we followed the [GCTA tutorial](https://yanglab.westlake.edu.cn/software/gcta/#Tutorial).
# 
# 1. Using only genotype data and phenotype. The name in the code is `GCTA_genotype`.
# 
# 1. `gcta --bfile FILE  --make-grm --out FILE`
# 2. `gcta --grm FILE --pheno FILE.PHENP --reml --out FILE`
# 
# 
# 2. Using genotype data, covariates, and phenotype. The name in the code is `GCTA_genotype_covariate`.
# 
# 1. `gcta --bfile FILE  --make-grm --out FILE`
# 2. `gcta --grm FILE --pheno FILE.PHENP  --qcovar  FILE.cov --reml --out FILE`
# 
# 
# 3. Using genotype data, covariates, PCA, and phenotype. The name in the code is `GCTA_genotype_covariate_pca`.
# 
# 1. `gcta --bfile FILE  --make-grm --out FILE`
# 2. `gcta --grm FILE --pheno FILE.PHENP --reml --qcovar  FILE.cov_pca --out FILE`
# 
# OR
# 
# 2. `gcta --grm FILE --pheno FILE.PHENP --reml-no-constrain --qcovar  FILE.cov_pca --out FILE`
# 
# **Handling REML Non-convergence:**
# 
# In some cases, the REML calculation may not converge ([source](https://gcta.freeforums.net/thread/366/error-log-likelihood-converged)). To address this, the `--reml-no-constrain` flag can be used. However, it's important to note that in such cases, the heritability value might exceed 1, often attributed to sample relatedness ([source](https://www.researchgate.net/post/What_does_it_mean_when_heritability_is_larger_than_1#:~:text=1)%20Heritability%20can%20be%20greater,can%20also%20cause%20this%20result.)).
# 
# An alternative approach involves using a flag like `--grm-cutoff 0.025`, similar to Plink's `--rel-cutoff 0.125`. This allows for the exclusion of specific samples, facilitating the recalculation of heritability.
# 
# ```bash
# gcta --grm FILE --pheno FILE.PHENP --reml --grm-cutoff 0.025 --qcovar FILE.cov_pca --out FILE
# ```
# 
# 
# 
# Heritability calculated using GCTA and using genotype data, covariates, PCA, and phenotype, and LDpred-2 method using all the methods were almost the same.
# 
#  
# 
# All GCTA functions generate the following file:
# 
# Summary result of REML analysis:
# 
# | Source | Variance | SE      |
# |--------|----------|---------|
# | V(G)   | 0.869695 | 0.790524|
# | V(e)   | 0.000001 | 0.787625|
# | Vp     | 0.869696 | 0.063776|
# | V(G)/Vp| 0.999999 | 0.905632|
# 
# h2 = V(G)/Vp
# 
 
  

# ### GWAS File Processing for GCTA
# 
# GCTA requires GWAS file in a specific format as specified in this document https://yanglab.westlake.edu.cn/software/gcta/#COJO. The GWAS we have contains all the required columns for processing.
 


import os
import pandas as pd
import numpy as np
import sys

filedirec = sys.argv[1]
GWAS = filedirec + os.sep + filedirec+".gz"

df = pd.read_csv(GWAS,compression= "gzip",sep="\s+")

import pandas as pd
if "BETA" in df.columns.to_list():
    df = df[['CHR', 'BP', 'SNP', 'A1', 'A2', 'N', 'SE', 'P', 'BETA', 'INFO', 'MAF']]
    print("yes")
else:
    df["BETA"] = np.log(df["OR"])
    df = df[['CHR', 'BP', 'SNP', 'A1', 'A2', 'N', 'SE', 'P', 'BETA', 'INFO', 'MAF']]

df_transformed = pd.DataFrame({
    'SNP': df['SNP'],
    'A1': df['A1'],
    'A2': df['A2'],
    'freq': df['MAF'],
    'b': df['BETA'],
    'se': df['SE'],
    'p': df['P'],
    'N': df['N']
})
 
output_file = filedirec+os.sep+"SBLUPGWAS.ma"
df_transformed.to_csv(output_file,sep="\t",index=False)
  
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

 
foldnumber = sys.argv[2]
#foldnumber = "4"  # Setting 'foldnumber' to "0"

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

prs_result = pd.DataFrame(columns=["clump_p1", "clump_r2", "clump_kb", "p_window_size", "p_slide_size", "p_LD_threshold",
                                "numberofpca","h2model","h2","numberofvariants","clumpprune"])


 


def transform_gcta_data(traindirec, newtrainfilename,models,p, clumpprune,p1_val, p2_val, p3_val, c1_val, c2_val, c3_val,Name):

    ### First perform clumping on the file and save the clumpled file.
    
    #clupmedfile = traindirec+os.sep+newtrainfilename+".clump"
    #prunedfile = traindirec+os.sep+newtrainfilename+".clump.prune"

    command = [
    "./plink",
    "--bfile", traindirec+os.sep+newtrainfilename,
    "--indep-pairwise", p1_val, p2_val, p3_val,
    "--out", traindirec+os.sep+trainfilename
    ]
    subprocess.run(command)
    # First perform pruning and then clumping and the pruning.
 
    command = [
    "./plink",
    "--bfile", traindirec+os.sep+newtrainfilename,
    "--clump-p1", c1_val,
    "--extract", traindirec+os.sep+trainfilename+".prune.in",
    "--clump-r2", c2_val,
    "--clump-kb", c3_val,
    "--clump", filedirec+os.sep+filedirec+".txt",
    "--clump-snp-field", "SNP",
    "--clump-field", "P",
    "--out", traindirec+os.sep+trainfilename
    ]    
    subprocess.run(command)

    # Extract the valid SNPs from th clumped file.
    # For windows download gwak for linux awk commmand is sufficient.
    ### For windows require GWAK.
    ### https://sourceforge.net/projects/gnuwin32/
    ##3 Get it and place it in the same direc.
    #os.system("gawk "+"\""+"NR!=1{print $3}"+"\"  "+ traindirec+os.sep+trainfilename+".clumped >  "+traindirec+os.sep+trainfilename+".valid.snp")
    #print("gawk "+"\""+"NR!=1{print $3}"+"\"  "+ traindirec+os.sep+trainfilename+".clumped >  "+traindirec+os.sep+trainfilename+".valid.snp")

    #Linux:
    os.system("awk "+"\""+"NR!=1{print $3}"+"\"  "+ traindirec+os.sep+trainfilename+".clumped >  "+traindirec+os.sep+trainfilename+".valid.snp")
    #print("awk "+"\""+"NR!=1{print $3}"+"\"  "+ traindirec+os.sep+trainfilename+".clumped >  "+traindirec+os.sep+trainfilename+".valid.snp")


    command = [
    "./plink",
    "--make-bed",
    "--bfile", traindirec+os.sep+newtrainfilename,
    "--indep-pairwise", p1_val, p2_val, p3_val,
    "--extract", traindirec+os.sep+trainfilename+".valid.snp",
    "--out", traindirec+os.sep+newtrainfilename+".clumped.pruned"
    ]
    subprocess.run(command)
    
    # Also extract the PCA at this point.
    command = [
        "./plink",
        "--bfile", traindirec+os.sep+newtrainfilename+".clumped.pruned",
        "--extract", traindirec+os.sep+trainfilename+".valid.snp",
        "--pca", p,
        "--out", traindirec+os.sep+trainfilename
    ]
    subprocess.run(command)
    
    
    # At this stage, we will merge the PCA and COV file. 
    tempphenotype_train = pd.read_table(traindirec+os.sep+newtrainfilename+".clumped.pruned"+".fam", sep="\s+",header=None)
    phenotype = pd.DataFrame()
    phenotype = tempphenotype_train[[0,1,5]]
    phenotype.to_csv(traindirec+os.sep+trainfilename+".PHENO",sep="\t",header=['FID', 'IID', 'PHENO'],index=False)
 
    pcs_train = pd.read_table(traindirec+os.sep+trainfilename+".eigenvec", sep="\s+",header=None, names=["FID", "IID"] + [f"PC{str(i)}" for i in range(1, int(p)+1)])
    covariate_train = pd.read_table(traindirec+os.sep+trainfilename+".cov",sep="\s+")
    covariate_train.fillna(0, inplace=True)
    covariate_train.to_csv(traindirec+os.sep+trainfilename+".cov",sep="\t",index=False)


    print(covariate_train.head())
    print(len(covariate_train))
    covariate_train = covariate_train[covariate_train["FID"].isin(pcs_train["FID"].values) & covariate_train["IID"].isin(pcs_train["IID"].values)]
    print(len(covariate_train))
 
    covariate_train['FID'] = covariate_train['FID'].astype(str)
    pcs_train['FID'] = pcs_train['FID'].astype(str)
    covariate_train['IID'] = covariate_train['IID'].astype(str)
    pcs_train['IID'] = pcs_train['IID'].astype(str)
    covandpcs_train = pd.merge(covariate_train, pcs_train, on=["FID","IID"])
    covandpcs_train.to_csv(traindirec+os.sep+trainfilename+".COV_PCA",sep="\t",index=False)

    if clumpprune=="yes":
        BFILE = traindirec+os.sep+newtrainfilename+".clumped.pruned"
        variants = len(pd.read_csv(BFILE+".bim"))
    else:
        BFILE = traindirec+os.sep+newtrainfilename
        variants = len(pd.read_csv(BFILE+".bim"))

    
    if models=="GCTA_genotype":
        if os.path.exists(traindirec+os.sep+"train_data.hsq"):
            os.remove(traindirec+os.sep+"train_data.hsq")
        command1 = [
            './gcta',
            '--bfile', BFILE,
            '--make-grm',
            '--out', traindirec+os.sep+trainfilename
        ]
        subprocess.run(command1)
 
        command2 = [
            './gcta',
            '--grm', traindirec+os.sep+trainfilename,
            '--pheno', traindirec+os.sep+trainfilename+".PHENO",
            #'--reml-no-constrain',
            '--reml',
            '--out', traindirec+os.sep+trainfilename
        ]
        subprocess.run(command2)
        try:
            tempdata = pd.read_csv(traindirec+os.sep+"train_data.hsq",sep="\t")
            tempdata = tempdata[tempdata["Source"]=="V(G)/Vp"]["Variance"].values[0]
        except:
            print("Heritibility not found using '--reml' ")
            command3 = [
            './gcta',
            '--grm', traindirec+os.sep+trainfilename,
            '--pheno', traindirec+os.sep+trainfilename+".PHENO",
            #'--qcovar', traindirec+os.sep+trainfilename+".cov",
            #'--reml-no-constrain',
            '--reml',
            '--grm-cutoff', str(0.01),
            '--out', traindirec+os.sep+trainfilename+"a"
            ]
            subprocess.run(command3)
            tempdata = pd.read_csv(traindirec+os.sep+"train_dataa.hsq",sep="\t")
            tempdata = tempdata[tempdata["Source"]=="V(G)/Vp"]["Variance"].values[0]
  
 
        print(tempdata)
      
    if models=="GCTA_genotype_covariate":
        if os.path.exists(traindirec+os.sep+"train_data.hsq"):
            os.remove(traindirec+os.sep+"train_data.hsq")
        
        command1 = [
            './gcta',
            '--bfile', BFILE,
            '--make-grm',
            '--out', traindirec+os.sep+trainfilename
        ]
        subprocess.run(command1)
        print(" ".join(command1))
        
        command3 = [
            './gcta',
            '--grm', traindirec+os.sep+trainfilename,
            '--pheno', traindirec+os.sep+trainfilename+".PHENO",
            '--qcovar', traindirec+os.sep+trainfilename+".cov",
            #'--reml-no-constrain',
            '--reml',
            '--out', traindirec+os.sep+trainfilename
        ]
        subprocess.run(command3)
        print(" ".join(command3))
      
        try:
            tempdata = pd.read_csv(traindirec+os.sep+"train_data.hsq",sep="\t")
            tempdata = tempdata[tempdata["Source"]=="V(G)/Vp"]["Variance"].values[0]
        except:
            print("Heritibility not found using '--reml' ")
            command3 = [
            './gcta',
            '--grm', traindirec+os.sep+trainfilename,
            '--pheno', traindirec+os.sep+trainfilename+".PHENO",
            '--qcovar', traindirec+os.sep+trainfilename+".cov",
            #'--reml-no-constrain',
            
            '--grm-cutoff', str(0.01),
            '--reml',
            '--out', traindirec+os.sep+trainfilename+"a"
            ]
            print(" ".join(command3))
            subprocess.run(command3)
            #exit(0)
            tempdata = pd.read_csv(traindirec+os.sep+"train_dataa.hsq",sep="\t")
            tempdata = tempdata[tempdata["Source"]=="V(G)/Vp"]["Variance"].values[0]

    
        
    if models=="GCTA_genotype_covariate_pca":
        if os.path.exists(traindirec+os.sep+"train_data.hsq"):
            os.remove(traindirec+os.sep+"train_data.hsq")
        
        command1 = [
            './gcta',
            '--bfile', BFILE,
            '--make-grm',
            '--out', traindirec+os.sep+trainfilename
        ]
        subprocess.run(command1)

        command3 = [
            './gcta',
            '--grm', traindirec+os.sep+trainfilename,
            '--pheno', traindirec+os.sep+trainfilename+".PHENO",
            '--qcovar', traindirec+os.sep+trainfilename+".COV_PCA",
            #'--reml-no-constrain',
            '--reml',
            '--out', traindirec+os.sep+trainfilename
        ]
        subprocess.run(command3)
        try:
            tempdata = pd.read_csv(traindirec+os.sep+"train_data.hsq",sep="\t")
            tempdata = tempdata[tempdata["Source"]=="V(G)/Vp"]["Variance"].values[0]
        except:
            print("Heritibility not found using '--reml' ")
            command3 = [
            './gcta',
            '--grm', traindirec+os.sep+trainfilename,
            '--pheno', traindirec+os.sep+trainfilename+".PHENO",
             '--qcovar', traindirec+os.sep+trainfilename+".COV_PCA",
            #'--reml-no-constrain',
            '--reml',
            '--grm-cutoff', str(0.01),
             '--out', traindirec+os.sep+trainfilename+"a"
            ]
            subprocess.run(command3)
            tempdata = pd.read_csv(traindirec+os.sep+"train_dataa.hsq",sep="\t")
            tempdata = tempdata[tempdata["Source"]=="V(G)/Vp"]["Variance"].values[0]
 

    global prs_result 
    prs_result = prs_result._append({
        "clump_p1": c1_val,
        "clump_r2": c2_val,
        "clump_kb": c3_val,
        "p_window_size": p1_val,
        "p_slide_size": p2_val,
        "p_LD_threshold": p3_val,
        
        "h2model":models,
        "h2":tempdata,
        "clumpprune":clumpprune,
        "numberofvariants":variants,
        "numberofpca":p,
          
    }, ignore_index=True)

 
    prs_result.to_csv(traindirec+os.sep+Name+os.sep+"Results.csv",index=False)
    
 
 
h2models = ["GCTA_genotype","GCTA_genotype_covariate","GCTA_genotype_covariate_pca"]
clumpprunes = ["yes","no"]
create_directory(folddirec+os.sep+"GCTA-Heritability")
for p1_val in p_window_size:
 for p2_val in p_slide_size: 
  for p3_val in p_LD_threshold:
   for c1_val in clump_p1:
    for c2_val in clump_r2:
     for c3_val in clump_kb:
      for p in numberofpca:
       for model in  h2models:
        for clumpprune in clumpprunes:
         transform_gcta_data(folddirec, newtrainfilename,model, p, clumpprune,str(p1_val), str(p2_val), str(p3_val), str(c1_val), str(c2_val), str(c3_val), "GCTA-Heritability")

exit(0)
