
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

variants = len(df)

import subprocess
def create_directory(directory):
    """Function to create a directory if it doesn't exist."""
    if not os.path.exists(directory):  # Checking if the directory doesn't exist
        os.makedirs(directory)  # Creating the directory if it doesn't exist
    return directory  # Returning the created or existing directory



prs_result = pd.DataFrame(columns=["refset","weightset","h2"])

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
    os.mkdir("LDSC")
except:
    pass    

try:
    os.remove("LDSC"+os.sep+filedirec)
except:
    pass
 


direc = "LDSC-Heritability"
create_directory(filedirec+os.sep+"Fold_0"+os.sep+direc)
create_directory(filedirec+os.sep+"Fold_1"+os.sep+direc)
create_directory(filedirec+os.sep+"Fold_2"+os.sep+direc)
create_directory(filedirec+os.sep+"Fold_3"+os.sep+direc)
create_directory(filedirec+os.sep+"Fold_4"+os.sep+direc)


refs = ["eur_ref_ld_chr","eur_w_ld_chr"]
weights = ["eur_ref_ld_chr","eur_w_ld_chr"]

for ref in refs:
    for weight in weights:
        command = [
            "python ldsc.py",
            "--h2", filedirec + os.sep +filedirec+".ldsc.txt",
            "--ref-ld-chr", "LDSCFILES"+os.sep+ref+"/",
            "--w-ld-chr", "LDSCFILES"+os.sep+weight+"/",
            "--out", "LDSC"+os.sep+filedirec
        ]
        
        subprocess.call(" ".join(command), shell=True)

        h2 = readpve("LDSC"+os.sep+filedirec+".log")
        variants = readvariants("LDSC"+os.sep+filedirec+".log")


        prs_result = prs_result.append({
            "h2":h2,
            "h2model":"LDSC"+"_"+ref+"_"+weight,
            "clumpprune":"no",
            "numberofvariants":variants,
            "numberofpca":"-"
            }, ignore_index=True)
        

        prs_result.to_csv(filedirec+os.sep+"Fold_0"+os.sep+direc+os.sep+"Results.csv")
        prs_result.to_csv(filedirec+os.sep+"Fold_1"+os.sep+direc+os.sep+"Results.csv")
        prs_result.to_csv(filedirec+os.sep+"Fold_2"+os.sep+direc+os.sep+"Results.csv")
        prs_result.to_csv(filedirec+os.sep+"Fold_3"+os.sep+direc+os.sep+"Results.csv")
        prs_result.to_csv(filedirec+os.sep+"Fold_4"+os.sep+direc+os.sep+"Results.csv")
exit(0)

