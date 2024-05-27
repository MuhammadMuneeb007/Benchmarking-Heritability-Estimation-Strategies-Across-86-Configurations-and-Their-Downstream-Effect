
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
    print("yes")
else:
    df["BETA"] = np.log(df["OR"])
    df = df[['CHR', 'BP', 'SNP', 'A1', 'A2', 'N', 'SE', 'P', 'BETA', 'INFO', 'MAF']]

df.to_csv(filedirec + os.sep +filedirec+".txt",sep="\t",index=False)
 
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

 
def transform_gemma_data(traindirec, newtrainfilename,models,p, clumpprune,relatedmatrix,data,p1_val, p2_val, p3_val, c1_val, c2_val, c3_val,Name):

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
 
 
    relatedmatrixname = ""

    if relatedmatrix=="1":
        relatedmatrixname = "centered"
        outputfilename = "output/"+traindirec+".cXX.txt"
    else: 
        relatedmatrixname = "standardized"
        outputfilename = "output/"+traindirec+".sXX.txt"
    
    h2modelname = ""


    

    if models=="1":
        h2modelname = "HE regression"
    else: 
        h2modelname = "REML AI algorithm"

    temppve = ""

    if clumpprune=="yes":
        BFILE = traindirec+os.sep+newtrainfilename+".clumped.pruned"
        variants = len(pd.read_csv(BFILE+".bim"))
    else:
        BFILE = traindirec+os.sep+newtrainfilename
        variants = len(pd.read_csv(BFILE+".bim"))


    def readpve(tempfile):
        file_path = "output"+os.sep+tempfile+".log.txt" 
        with open(file_path, 'r') as file:
            for line in file:
                if "pve estimates" in line:
                    pve_estimate = line.split("=")[-1].strip()
                    print("PVE Estimate:", pve_estimate)
                    return pve_estimate

    def readvariants(tempfile):
        file_path = "output"+os.sep+tempfile+".log.txt"  
        with open(file_path, 'r') as file:
            for line in file:
                if "## number of analyzed SNPs" in line:
                    variants = line.split("=")[-1].strip()
                    print("variants:", variants)
                    return variants
            #exit(0)
    try:
        os.mkdir("output/"+filedirec)
    except:
        pass

    try:
        os.remove("output/"+traindirec+".sXX.txt")
    except:
        pass

    try:
        os.remove("output/"+traindirec+".cXX.txt")
    except:
        pass
    try:
        os.remove("output"+os.sep+traindirec+".log.txt" )
    except:
        pass
    

    if data == "genotype":
         
        print("Processing genotype data")
        subprocess.run(["./gemma",
                "--bfile", BFILE,
                "-gk", relatedmatrix,
                "-o", traindirec])
 
        tempvariants = readvariants(traindirec)        
        print(BFILE)
        # Command 2: ./gemma -p train_data.fam -k output/output.cXX.txt -n 6 -vc 1 -o output
        subprocess.run(["./gemma",
                        "-p", BFILE+".fam",
                        "-k", outputfilename,
                        "-n", "6",
                        "-vc", model,
                        "-o", traindirec])
        temppve =  readpve(traindirec)  


    elif data == "genotype_covariate":
        subprocess.run(["./gemma",
                "--bfile", BFILE,
                "-gk", relatedmatrix,
                "-o", traindirec])      
                 
 
        tempvariants = readvariants(traindirec)
        subprocess.run([
                './gemma',
                '-p', BFILE+".fam",
                '-k', outputfilename,
                '-n', '6',
                '-vc', model,
                '-c', traindirec+os.sep+trainfilename+".covgemma",
                '-o', traindirec
            ])
 
        temppve =  readpve(traindirec)    

    elif data == "genotype_covariate_pca":
        subprocess.run(["./gemma",
                "--bfile", BFILE,
                "-gk", relatedmatrix,
                "-o", traindirec]) 
                   
 
        tempvariants = readvariants(traindirec)         
        subprocess.run([
                './gemma',
                '-p', BFILE+".fam",
                '-k', outputfilename,
                '-n', '6',
                '-vc', model,
                '-c', traindirec+os.sep+trainfilename+".COV_PCAgemma",
                '-o', traindirec
            ])
        temppve =  readpve(traindirec)  


    

    global prs_result 
    prs_result = prs_result._append({
        "clump_p1": c1_val,
        "clump_r2": c2_val,
        "clump_kb": c3_val,
        "p_window_size": p1_val,
        "p_slide_size": p2_val,
        "p_LD_threshold": p3_val,

        "h2":temppve,
        "h2model":h2modelname+"_"+relatedmatrixname+"_"+data,
        "clumpprune":clumpprune,
        "numberofvariants":tempvariants,
        "numberofpca":p,
        #"relatedmatrix":relatedmatrixname,
        #"data":data,
         
    }, ignore_index=True)

    prs_result.to_csv(traindirec+os.sep+Name+os.sep+"Results.csv",index=False)
    
 
    return


# `-gk 1` calculates the centered relatedness matrix, while `-gk 2` calculates the standardized relatedness matrix.
relatedmatrixs = ["1", "2"]
# `-vc 1` (default) uses HE regression, and `-vc 2` uses REML AI algorithm.
h2models = ["1", "2"]
datas = ["genotype", "genotype_covariate", "genotype_covariate_pca"]

clumpprunes = ["yes","no"]
# Nested loops to iterate over different parameter values

create_directory(folddirec+os.sep+"Gemma1-Heritability")
for p1_val in p_window_size:
 for p2_val in p_slide_size: 
  for p3_val in p_LD_threshold:
   for c1_val in clump_p1:
    for c2_val in clump_r2:
     for c3_val in clump_kb:
      for p in numberofpca:
       for model in  h2models:
        for clumpprune in clumpprunes:
         for relatedmatrix in relatedmatrixs:
          for data  in datas:
            transform_gemma_data(folddirec, newtrainfilename,model, p, clumpprune,relatedmatrix,data,str(p1_val), str(p2_val), str(p3_val), str(c1_val), str(c2_val), str(c3_val), "Gemma1-Heritability")
             
exit(0)




