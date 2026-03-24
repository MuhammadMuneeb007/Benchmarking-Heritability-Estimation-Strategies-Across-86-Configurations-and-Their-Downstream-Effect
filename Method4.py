import os
import sys
import pandas as pd
import numpy as np
import subprocess
import statsmodels.api as sm
from sklearn.metrics import roc_auc_score, confusion_matrix
from statsmodels.stats.contingency_tables import mcnemar

filedirec = sys.argv[1]
foldnumber = sys.argv[2]

GWAS = filedirec + os.sep + filedirec + ".gz"
df = pd.read_csv(GWAS, compression="gzip", sep=r"\s+")

if "BETA" in df.columns.to_list():
    df = df[['CHR', 'BP', 'SNP', 'A1', 'A2', 'N', 'SE', 'P', 'BETA', 'INFO', 'MAF']]
    df['Z'] = df['BETA'] / df['SE']
    df = df[['SNP', 'N', 'Z', 'A1', 'A2']]
else:
    df["BETA"] = np.log(df["OR"])
    df = df[['CHR', 'BP', 'SNP', 'A1', 'A2', 'N', 'SE', 'P', 'BETA', 'INFO', 'MAF']]
    df['Z'] = df['BETA'] / df['SE']
    df = df[['SNP', 'N', 'Z', 'A1', 'A2']]

df.to_csv(filedirec + os.sep + filedirec + ".gemma.txt", sep="\t", index=False)


def create_directory(directory):
    if not os.path.exists(directory):
        os.makedirs(directory)
    return directory


folddirec = filedirec + os.sep + "Fold_" + foldnumber
trainfilename = "train_data"
newtrainfilename = "train_data.QC"
testfilename = "test_data"
newtestfilename = "test_data.QC"

numberofpca = ["6"]
clump_p1 = [1]
clump_r2 = [0.1]
clump_kb = [200]
p_window_size = [200]
p_slide_size = [50]
p_LD_threshold = [0.25]

prs_result = pd.DataFrame(columns=[
    "clump_p1", "clump_r2", "clump_kb",
    "p_window_size", "p_slide_size", "p_LD_threshold",
    "numberofpca", "h2model", "h2", "h2_se",
    "numberofvariants", "clumpprune"
])


def readvariants(tempfile):
    file_path = "output" + os.sep + tempfile + ".log.txt"
    with open(file_path, 'r') as file:
        for line in file:
            if "## number of analyzed SNPs" in line:
                variants = line.split("=")[-1].strip()
                print("variants:", variants)
                return variants
    return np.nan


def readpve(tempfile):
    file_path = "output" + os.sep + tempfile + ".log.txt"
    with open(file_path, 'r') as file:
        for line in file:
            if "pve estimates" in line:
                pve_estimate = line.split("=")[-1].strip()
                print("PVE Estimate:", pve_estimate)
                return pve_estimate
    return np.nan


def readpve_se(tempfile):
    file_path = "output" + os.sep + tempfile + ".log.txt"
    with open(file_path, 'r') as file:
        for line in file:
            if "se(pve)" in line:
                pve_se = line.split("=")[-1].strip()
                print("PVE SE:", pve_se)
                return pve_se
    return np.nan


def transform_gemma_data(traindirec, newtrainfilename, models, p, clumpprune, relatedmatrix, data,
                          p1_val, p2_val, p3_val, c1_val, c2_val, c3_val, Name):

    results_file = traindirec + os.sep + Name + os.sep + "Results.csv"
    if os.path.exists(results_file):
        os.remove(results_file)

    command = [
        "./plink",
        "--bfile", traindirec + os.sep + newtrainfilename,
        "--indep-pairwise", p1_val, p2_val, p3_val,
        "--out", traindirec + os.sep + trainfilename
    ]
    subprocess.call(command)

    command = [
        "./plink",
        "--bfile", traindirec + os.sep + newtrainfilename,
        "--clump-p1", c1_val,
        "--extract", traindirec + os.sep + trainfilename + ".prune.in",
        "--clump-r2", c2_val,
        "--clump-kb", c3_val,
        "--clump", filedirec + os.sep + filedirec + ".txt",
        "--clump-snp-field", "SNP",
        "--clump-field", "P",
        "--out", traindirec + os.sep + trainfilename
    ]
    subprocess.call(command)

    os.system("awk " + "\"" + "NR!=1{print $3}" + "\"  " +
              traindirec + os.sep + trainfilename + ".clumped >  " +
              traindirec + os.sep + trainfilename + ".valid.snp")

    command = [
        "./plink",
        "--make-bed",
        "--bfile", traindirec + os.sep + newtrainfilename,
        "--indep-pairwise", p1_val, p2_val, p3_val,
        "--extract", traindirec + os.sep + trainfilename + ".valid.snp",
        "--out", traindirec + os.sep + newtrainfilename + ".clumped.pruned"
    ]
    subprocess.call(command)

    command = [
        "./plink",
        "--bfile", traindirec + os.sep + newtrainfilename + ".clumped.pruned",
        "--extract", traindirec + os.sep + trainfilename + ".valid.snp",
        "--pca", p,
        "--out", traindirec + os.sep + trainfilename
    ]
    subprocess.call(command)

    tempphenotype_train = pd.read_table(
        traindirec + os.sep + newtrainfilename + ".clumped.pruned" + ".fam",
        sep=r"\s+", header=None
    )
    phenotype = tempphenotype_train[[0, 1, 5]]
    phenotype.to_csv(traindirec + os.sep + trainfilename + ".PHENO",
                     sep="\t", header=['FID', 'IID', 'PHENO'], index=False)

    pcs_train = pd.read_table(
        traindirec + os.sep + trainfilename + ".eigenvec",
        sep=r"\s+", header=None,
        names=["FID", "IID"] + ["PC{}".format(str(i)) for i in range(1, int(p) + 1)]
    )
    covariate_train = pd.read_table(
        traindirec + os.sep + trainfilename + ".cov", sep=r"\s+"
    )
    covariate_train.iloc[:, 2:].to_csv(
        traindirec + os.sep + trainfilename + ".covgemma",
        header=False, index=False, sep="\t"
    )
    covariate_train.fillna(0, inplace=True)
    print(covariate_train.head())
    print(len(covariate_train))
    covariate_train = covariate_train[
        covariate_train["FID"].isin(pcs_train["FID"].values) &
        covariate_train["IID"].isin(pcs_train["IID"].values)
    ]
    print(len(covariate_train))

    covariate_train['FID'] = covariate_train['FID'].astype(str)
    pcs_train['FID'] = pcs_train['FID'].astype(str)
    covariate_train['IID'] = covariate_train['IID'].astype(str)
    pcs_train['IID'] = pcs_train['IID'].astype(str)
    covandpcs_train = pd.merge(covariate_train, pcs_train, on=["FID", "IID"])
    covandpcs_train.to_csv(traindirec + os.sep + trainfilename + ".COV_PCA",
                           sep="\t", index=False)
    covandpcs_train.iloc[:, 2:].to_csv(
        traindirec + os.sep + trainfilename + ".COV_PCAgemma",
        header=False, index=False, sep="\t"
    )

    if relatedmatrix == "1":
        relatedmatrixname = "centered"
    else:
        relatedmatrixname = "standardized"

    if models == "1":
        h2modelname = "HE regression"
    else:
        h2modelname = "REML AI algorithm"

    if clumpprune == "yes":
        BFILE = traindirec + os.sep + newtrainfilename + ".clumped.pruned"
    else:
        BFILE = traindirec + os.sep + newtrainfilename

    temppve = np.nan
    temppve_se = np.nan
    tempvariants = np.nan

    try:
        os.mkdir("output/" + filedirec)
    except:
        pass

    try:
        os.remove("output/" + traindirec + ".sXX.txt")
    except:
        pass

    try:
        os.remove("output/" + traindirec + ".cXX.txt")
    except:
        pass

    try:
        os.remove("output" + os.sep + traindirec + ".log.txt")
    except:
        pass

    if data == "genotype":
        print("Processing genotype data")
        subprocess.call(["./gemma",
                         "-beta", filedirec + os.sep + filedirec + ".gemma.txt",
                         "-bfile", BFILE,
                         "-vc", models,
                         "-o", traindirec])
        temppve = readpve(traindirec)
        temppve_se = readpve_se(traindirec)
        tempvariants = readvariants(traindirec)

    elif data == "genotype_covariate":
        subprocess.run([
            './gemma',
            "-beta", filedirec + os.sep + filedirec + ".gemma.txt",
            "-bfile", BFILE,
            "-vc", models,
            '-c', traindirec + os.sep + trainfilename + ".covgemma",
            '-o', traindirec
        ])
        temppve = readpve(traindirec)
        temppve_se = readpve_se(traindirec)
        tempvariants = readvariants(traindirec)

    elif data == "genotype_covariate_pca":
        subprocess.run([
            './gemma',
            "-beta", filedirec + os.sep + filedirec + ".gemma.txt",
            "-bfile", BFILE,
            "-vc", models,
            '-c', traindirec + os.sep + trainfilename + ".COV_PCAgemma",
            '-o', traindirec
        ])
        temppve = readpve(traindirec)
        temppve_se = readpve_se(traindirec)
        tempvariants = readvariants(traindirec)

    print("Heritability:", temppve, "SE:", temppve_se, "Variants:", tempvariants)

    global prs_result
    new_row = pd.DataFrame([{
        "clump_p1":       c1_val,
        "clump_r2":       c2_val,
        "clump_kb":       c3_val,
        "p_window_size":  p1_val,
        "p_slide_size":   p2_val,
        "p_LD_threshold": p3_val,
        "numberofpca":    p,
        "h2":             temppve,
        "h2_se":          temppve_se,
        "h2model":        h2modelname + "_" + data,
        "clumpprune":     clumpprune,
        "numberofvariants": tempvariants
    }])
    prs_result = pd.concat([prs_result, new_row], ignore_index=True)

    outdir = traindirec + os.sep + Name
    create_directory(outdir)
    prs_result.to_csv(outdir + os.sep + "Results.csv", index=False)
    return


relatedmatrixs = ["1"]
h2models = ["1"]
datas = ["genotype", "genotype_covariate", "genotype_covariate_pca"]
clumpprunes = ["yes", "no"]

create_directory(folddirec + os.sep + "Gemma2-Heritability")

for p1_val in p_window_size:
    for p2_val in p_slide_size:
        for p3_val in p_LD_threshold:
            for c1_val in clump_p1:
                for c2_val in clump_r2:
                    for c3_val in clump_kb:
                        for p in numberofpca:
                            for model in h2models:
                                for clumpprune in clumpprunes:
                                    for relatedmatrix in relatedmatrixs:
                                        for data in datas:
                                            transform_gemma_data(
                                                folddirec, newtrainfilename, model, p, clumpprune,
                                                relatedmatrix, data,
                                                str(p1_val), str(p2_val), str(p3_val),
                                                str(c1_val), str(c2_val), str(c3_val),
                                                "Gemma2-Heritability"
                                            )

exit(0)