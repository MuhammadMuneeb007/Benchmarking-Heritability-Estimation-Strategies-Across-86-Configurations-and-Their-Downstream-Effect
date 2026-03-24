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
else:
    df["BETA"] = np.log(df["OR"])
    df = df[['CHR', 'BP', 'SNP', 'A1', 'A2', 'N', 'SE', 'P', 'BETA', 'INFO', 'MAF']]
    df['Z'] = df['BETA'] / df['SE']

df.to_csv(filedirec + os.sep + filedirec + ".ldsc.txt", sep="\t", index=False)


def create_directory(directory):
    if not os.path.exists(directory):
        os.makedirs(directory)
    return directory


def readpve(file_path):
    with open(file_path, 'r') as file:
        for line in file:
            if "Total Observed scale h2" in line:
                pve_estimate = line.split(":")[-1].split("(")[0].strip()
                print("PVE Estimate:", pve_estimate)
                return pve_estimate
    return np.nan


def readpve_se(file_path):
    with open(file_path, 'r') as file:
        for line in file:
            if "Total Observed scale h2" in line:
                if "(" in line and ")" in line:
                    pve_se = line.split("(")[-1].replace(")", "").strip()
                    print("PVE SE:", pve_se)
                    return pve_se
    return np.nan


def readvariants(file_path):
    with open(file_path, 'r') as file:
        for line in file:
            if "After merging with regression SNP LD" in line:
                pve_estimate = line.split(",")[-1].strip().split(" ")[0]
                print("Variants:", pve_estimate)
                return pve_estimate
    return np.nan


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


def transform_gemma_data(traindirec, newtrainfilename, p, clumpprune,
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

    if clumpprune == "yes":
        BFILE = traindirec + os.sep + newtrainfilename + ".clumped.pruned"
    else:
        BFILE = traindirec + os.sep + newtrainfilename

    try:
        os.mkdir("LDSCFILES/" + filedirec)
    except:
        pass

    try:
        os.mkdir("LDSCFILES/" + traindirec)
    except:
        pass

    chromosomes = range(1, 23)
    for chr_num in chromosomes:
        command = [
            "./plink",
            "--bfile", BFILE,
            "--chr", str(chr_num),
            "--make-bed",
            "--out", BFILE + "_" + str(chr_num)
        ]
        subprocess.call(' '.join(command), shell=True)

        command = [
            "python", "ldsc.py",
            "--bfile", BFILE + "_" + str(chr_num),
            "--l2",
            "--yes-really",
            "--ld-wind-cm", "1",
            "--out", "LDSCFILES" + os.sep + traindirec + os.sep + str(chr_num)
        ]
        subprocess.call(' '.join(command), shell=True)

    command = [
        "python ldsc.py",
        "--h2", filedirec + os.sep + filedirec + ".ldsc.txt",
        "--ref-ld-chr", "LDSCFILES" + os.sep + traindirec + "/",
        "--w-ld-chr", "LDSCFILES" + os.sep + traindirec + "/",
        "--out", "LDSC" + os.sep + filedirec
    ]
    subprocess.call(" ".join(command), shell=True)

    log_file = "LDSC" + os.sep + filedirec + ".log"
    h2 = readpve(log_file)
    h2_se = readpve_se(log_file)
    variants = readvariants(log_file)

    print("Heritability:", h2, "SE:", h2_se, "Variants:", variants)

    global prs_result
    new_row = pd.DataFrame([{
        "clump_p1":       c1_val,
        "clump_r2":       c2_val,
        "clump_kb":       c3_val,
        "p_window_size":  p1_val,
        "p_slide_size":   p2_val,
        "p_LD_threshold": p3_val,
        "numberofpca":    "-",
        "h2":             h2,
        "h2_se":          h2_se,
        "h2model":        "LDSC_Genotype_Reference",
        "clumpprune":     clumpprune,
        "numberofvariants": variants
    }])
    prs_result = pd.concat([prs_result, new_row], ignore_index=True)

    outdir = traindirec + os.sep + Name
    create_directory(outdir)
    prs_result.to_csv(outdir + os.sep + "Results.csv", index=False)
    return


clumpprunes = ["yes"]

create_directory(folddirec + os.sep + "LDSC-Genotype-Heritability")

for p1_val in p_window_size:
    for p2_val in p_slide_size:
        for p3_val in p_LD_threshold:
            for c1_val in clump_p1:
                for c2_val in clump_r2:
                    for c3_val in clump_kb:
                        for p in numberofpca:
                            for clumpprune in clumpprunes:
                                transform_gemma_data(
                                    folddirec, newtrainfilename, p, clumpprune,
                                    str(p1_val), str(p2_val), str(p3_val),
                                    str(c1_val), str(c2_val), str(c3_val),
                                    "LDSC-Genotype-Heritability"
                                )

exit(0)