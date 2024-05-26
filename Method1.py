# Importing necessary libraries
import os
import sys
import pandas as pd
import numpy as np
import statsmodels.api as sm
from sklearn.metrics import roc_auc_score, confusion_matrix
from statsmodels.stats.contingency_tables import mcnemar
import subprocess

# Retrieving command line arguments
filedirec = sys.argv[1]  # Base directory path
foldnumber = sys.argv[2]  # Fold number

# Creating directory if it doesn't exist
def create_directory(directory):
    if not os.path.exists(directory):
        os.makedirs(directory)
    return directory

#"""
# Reading GWAS data
GWAS = filedirec + os.sep + filedirec + ".gz"
df = pd.read_csv(GWAS, compression="gzip", sep="\s+")

# Checking for 'BETA' column in the GWAS data
# If not present, converting 'OR' to 'BETA' and updating DataFrame
if "BETA" in df.columns.to_list():
    df = df[['CHR', 'BP', 'SNP', 'A1', 'A2', 'N', 'SE', 'P', 'BETA', 'INFO', 'MAF']]
else:
    df["BETA"] = np.log(df["OR"])
    df = df[['CHR', 'BP', 'SNP', 'A1', 'A2', 'N', 'SE', 'P', 'BETA', 'INFO', 'MAF']]

# Saving processed GWAS data to file
df.to_csv(filedirec + os.sep + filedirec + ".txt", sep="\t", index=False)
#"""

# Setting up parameters and variables for further analysis
folddirec = filedirec + os.sep + "Fold_" + foldnumber  # Directory path for the specific fold
trainfilename = "train_data"  # Training data file name
newtrainfilename = "train_data.QC"  # New training data file name
testfilename = "test_data"  # Test data file name
newtestfilename = "test_data.QC"  # New test data file name
numberofpca = ["6"]  # Number of PCA components to be included
clump_p1 = [1]  # Clumping parameter 'p1'
clump_r2 = [0.1]  # Clumping parameter 'r2'
clump_kb = [200]  # Clumping parameter 'kb'
p_window_size = [200]  # Pruning parameter 'window_size'
p_slide_size = [50]  # Pruning parameter 'slide_size'
p_LD_threshold = [0.25]  # Pruning parameter 'LD_threshold'

# Initializing DataFrame to store results
prs_result = pd.DataFrame(columns=["clump_p1", "clump_r2", "clump_kb", "p_window_size", "p_slide_size", 
                                   "p_LD_threshold", "numberofpca", "h2model", "h2", "numberofvariants", "clumpprune"])


def transform_ldpred2_data(traindirec, newtrainfilename,models,p, clumpprune,p1_val, p2_val, p3_val, c1_val, c2_val, c3_val,Name):

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


 
    if models=="LDpred-2_full":
        if clumpprune=="yes":
            os.system("Rscript LDpred-2-Heritability.R "+os.path.join(filedirec)+"  "+traindirec+" "+trainfilename+" "+newtrainfilename+".clumped.pruned"+ " "+"3"+" "+c3_val+" "+c1_val+" "+c2_val+" "+p)
            print("Rscript LDpred-2-Heritability.R "+os.path.join(filedirec)+"  "+traindirec+" "+trainfilename+" "+newtrainfilename+".clumped.pruned"+ " "+"3"+" "+c3_val+" "+c1_val+" "+c2_val+" "+p)
            tempdata1 = pd.read_csv(traindirec+os.sep+"ldpred_h2_full.txt",sep="\s+",header=None)[1].values[0]
            tempdata2 = pd.read_csv(traindirec+os.sep+"ldpred_h2_full_variants.txt",sep="\s+",header=None)[1].values[0]
        
        else:
            os.system("Rscript LDpred-2-Heritability.R "+os.path.join(filedirec)+"  "+traindirec+" "+trainfilename+" "+newtrainfilename+ " "+"3"+" "+c3_val+" "+c1_val+" "+c2_val+" "+p)
            tempdata1 = pd.read_csv(traindirec+os.sep+"ldpred_h2_full.txt",sep="\s+",header=None)[1].values[0]
            tempdata2 = pd.read_csv(traindirec+os.sep+"ldpred_h2_full_variants.txt",sep="\s+",header=None)[1].values[0]
        
        
       

    if models=="LDpred-2_hapmap":
        if clumpprune=="yes":
            os.system("Rscript LDpred-2-Heritability.R "+os.path.join(filedirec)+"  "+traindirec+" "+trainfilename+" "+newtrainfilename+".clumped.pruned"+ " "+"2"+" "+c3_val+" "+c1_val+" "+c2_val+" "+p)
            tempdata1 = pd.read_csv(traindirec+os.sep+"ldpred_h2_hapmap.txt",sep="\s+",header=None)[1].values[0]
            tempdata2 = pd.read_csv(traindirec+os.sep+"ldpred_h2_hapmap_variants.txt",sep="\s+",header=None)[1].values[0]
        else:
            os.system("Rscript LDpred-2-Heritability.R "+os.path.join(filedirec)+"  "+traindirec+" "+trainfilename+" "+newtrainfilename+ " "+"2"+" "+c3_val+" "+c1_val+" "+c2_val+" "+p)
            tempdata1 = pd.read_csv(traindirec+os.sep+"ldpred_h2_hapmap.txt",sep="\s+",header=None)[1].values[0]
            tempdata2 = pd.read_csv(traindirec+os.sep+"ldpred_h2_hapmap_variants.txt",sep="\s+",header=None)[1].values[0]



    global prs_result 
    prs_result = prs_result._append({
        "clump_p1": c1_val,
        "clump_r2": c2_val,
        "clump_kb": c3_val,
        "p_window_size": p1_val,
        "p_slide_size": p2_val,
        "p_LD_threshold": p3_val,
        "numberofpca":p,
        "h2":tempdata1,
        "h2model":models,
        "clumpprune":clumpprune,
        "numberofvariants":tempdata2,
        
    }, ignore_index=True)

    prs_result.to_csv(traindirec+os.sep+Name+os.sep+"Results.csv",index=False)
    
        
    print("Heritability Value",tempdata1)
    print("Heritability Variants",tempdata2)
    #exit(0)
    return

 


h2models = ["LDpred-2_full","LDpred-2_hapmap"]
clumpprunes = ["yes","no"]

create_directory(folddirec+os.sep+"LDpred-2-Heritability")
for p1_val in p_window_size:
 for p2_val in p_slide_size: 
  for p3_val in p_LD_threshold:
   for c1_val in clump_p1:
    for c2_val in clump_r2:
     for c3_val in clump_kb:
      for p in numberofpca:
       for model in  h2models:
        for clumpprune in clumpprunes:
         transform_ldpred2_data(folddirec, newtrainfilename,model, p, clumpprune,str(p1_val), str(p2_val), str(p3_val), str(c1_val), str(c2_val), str(c3_val), "LDpred-2-Heritability")

exit(0)

