
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

df.to_csv(filedirec + os.sep +filedirec+".ldsc.txt",sep="\t",index=False)

 

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
                                "numberofpca","h2model","h2","numberofvariants","clumpprune","relatedmatrix","data"])

 
def transform_gemma_data(traindirec, newtrainfilename,p, clumpprune,p1_val, p2_val, p3_val, c1_val, c2_val, c3_val,Name):

    ### First perform clumping on the file and save the clumpled file.
    
    #clupmedfile = traindirec+os.sep+newtrainfilename+".clump"
    #prunedfile = traindirec+os.sep+newtrainfilename+".clump.prune"
 

    command = [
    "./plink",
    "--bfile", traindirec+os.sep+newtrainfilename,
    "--indep-pairwise", p1_val, p2_val, p3_val,
    "--out", traindirec+os.sep+trainfilename
    ]
    subprocess.call(command)
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
    subprocess.call(command)

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
    subprocess.call(command)
    
    # Also extract the PCA at this point.
    command = [
        "./plink",
        "--bfile", traindirec+os.sep+newtrainfilename+".clumped.pruned",
        "--extract", traindirec+os.sep+trainfilename+".valid.snp",
        "--pca", p,
        "--out", traindirec+os.sep+trainfilename
    ]
    subprocess.call(command)
    
    # At this stage, we will merge the PCA and COV file. 
    tempphenotype_train = pd.read_table(traindirec+os.sep+newtrainfilename+".clumped.pruned"+".fam", sep="\s+",header=None)
    phenotype = pd.DataFrame()
    phenotype = tempphenotype_train[[0,1,5]]
    phenotype.to_csv(traindirec+os.sep+trainfilename+".PHENO",sep="\t",header=['FID', 'IID', 'PHENO'],index=False)
 
    #pcs_train = pd.read_table(traindirec+os.sep+trainfilename+".eigenvec", sep="\s+",header=None, names=["FID", "IID"] + [f"PC{str(i)}" for i in range(1, int(p)+1)])
    pcs_train = pd.read_table(traindirec+os.sep+trainfilename+".eigenvec", sep="\s+", header=None, names=["FID", "IID"] + ["PC{}".format(str(i)) for i in range(1, int(p) + 1)])

    covariate_train = pd.read_table(traindirec+os.sep+trainfilename+".cov",sep="\s+")
    covariate_train.iloc[:, 2:].to_csv(traindirec+os.sep+trainfilename+".covgemma", header=False, index=False,sep="\t")
 
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
    covandpcs_train.iloc[:, 2:].to_csv(traindirec+os.sep+trainfilename+".COV_PCAgemma", header=False, index=False,sep="\t")
 
  

    if clumpprune=="yes":
        BFILE = traindirec+os.sep+newtrainfilename+".clumped.pruned"
    else:
        BFILE = traindirec+os.sep+newtrainfilename
 

    def readpve(file): 
        file_path = file
        with open(file_path, 'r') as file:
            for line in file:
                if "Total Observed scale h2" in line:
                    pve_estimate = line.split(":")[-1].split(" ")[1]
                    print("PVE Estimate:", pve_estimate)
                    return pve_estimate

    def readvariants(file): 
        file_path = file
        with open(file_path, 'r') as file:
            for line in file:
                if "After merging with regression SNP LD" in line:
                    pve_estimate = line.split(",")[-1].strip().split(" ")[0]
                    print("Variants:", pve_estimate)
                    return pve_estimate
    try:
        os.mkdir("LDSCFILES/"+filedirec)
    except:
        pass    

    
    try:
        os.mkdir("LDSCFILES/"+traindirec)
    except:
        pass    
    


    chromosomes = range(1, 23)   
 
    for chr_num in chromosomes:
        command = [
            "./plink",
            "--bfile", BFILE,
            "--chr", str(chr_num),
            "--make-bed",
            "--out", BFILE+"_"+str(chr_num)
        ]
        command_str = ' '.join(command)
        subprocess.call(command_str, shell=True)

        command = [
            "python", "ldsc.py",
            "--bfile", BFILE+"_"+str(chr_num),
            "--l2",
            "--yes-really",
            "--ld-wind-cm", "1",
            "--out", "LDSCFILES"+os.sep+traindirec+os.sep+str(chr_num)
        ]
        command_str = ' '.join(command)
        subprocess.call(command_str, shell=True)
    
    command = [
        "python ldsc.py",
        "--h2", filedirec + os.sep +filedirec+".ldsc.txt",
        "--ref-ld-chr", "LDSCFILES"+os.sep+traindirec+"/",
        "--w-ld-chr", "LDSCFILES"+os.sep+traindirec+"/",
        "--out", "LDSC"+os.sep+filedirec
    ]
    
    subprocess.call(" ".join(command), shell=True)
    h2 = readpve("LDSC"+os.sep+filedirec+".log")
    variants = readvariants("LDSC"+os.sep+filedirec+".log")

    global prs_result

    prs_result = prs_result.append({
            "clump_p1": c1_val,
            "clump_r2": c2_val,
            "clump_kb": c3_val,
            "p_window_size": p1_val,
            "p_slide_size": p2_val,
            "p_LD_threshold": p3_val,
            "h2":h2,
            "h2model":"LDSC"+"_"+"Genotype_Reference",
            "clumpprune":clumpprune,
            "numberofvariants":variants,
            "numberofpca":"-"
            }, ignore_index=True)

    prs_result.to_csv(traindirec+os.sep+Name+os.sep+"Results.csv")
 
 
    return

 
clumpprunes = ["yes","no"]
clumpprunes = ["yes"]

# Nested loops to iterate over different parameter values
create_directory(folddirec+os.sep+"LDSC-Genotype-Heritability")
for p1_val in p_window_size:
 for p2_val in p_slide_size: 
  for p3_val in p_LD_threshold:
   for c1_val in clump_p1:
    for c2_val in clump_r2:
     for c3_val in clump_kb:
      for p in numberofpca:
        for clumpprune in clumpprunes:
            transform_gemma_data(folddirec, newtrainfilename, p, clumpprune,str(p1_val), str(p2_val), str(p3_val), str(c1_val), str(c2_val), str(c3_val), "LDSC-Genotype-Heritability")
            #exit(0)
exit(0)




