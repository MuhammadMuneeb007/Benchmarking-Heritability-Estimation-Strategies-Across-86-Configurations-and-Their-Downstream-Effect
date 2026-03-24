import os
import sys
import pandas as pd
import numpy as np
import subprocess

filedirec = sys.argv[1]

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


prs_result = pd.DataFrame(columns=[
    "h2model", "h2", "h2_se",
    "numberofvariants", "clumpprune", "numberofpca"
])

try:
    os.mkdir("LDSC")
except:
    pass

try:
    os.remove("LDSC" + os.sep + filedirec)
except:
    pass

direc = "LDSC-Heritability"
folds = ["Fold_0", "Fold_1", "Fold_2", "Fold_3", "Fold_4"]
for fold in folds:
    create_directory(filedirec + os.sep + fold + os.sep + direc)
    results_file = filedirec + os.sep + fold + os.sep + direc + os.sep + "Results.csv"
    if os.path.exists(results_file):
        os.remove(results_file)

refs = ["eur_ref_ld_chr", "eur_w_ld_chr"]
weights = ["eur_ref_ld_chr", "eur_w_ld_chr"]

for ref in refs:
    for weight in weights:
        command = [
            "python ldsc.py",
            "--h2", filedirec + os.sep + filedirec + ".ldsc.txt",
            "--ref-ld-chr", "LDSCFILES" + os.sep + ref + "/",
            "--w-ld-chr", "LDSCFILES" + os.sep + weight + "/",
            "--out", "LDSC" + os.sep + filedirec
        ]
        subprocess.call(" ".join(command), shell=True)

        log_file = "LDSC" + os.sep + filedirec + ".log"
        h2 = readpve(log_file)
        h2_se = readpve_se(log_file)
        variants = readvariants(log_file)

        print("Heritability:", h2, "SE:", h2_se, "Variants:", variants)

        new_row = pd.DataFrame([{
            "h2":             h2,
            "h2_se":          h2_se,
            "h2model":        "LDSC" + "_" + ref + "_" + weight,
            "clumpprune":     "no",
            "numberofvariants": variants,
            "numberofpca":    "-"
        }])
        prs_result = pd.concat([prs_result, new_row], ignore_index=True)

        for fold in folds:
            prs_result.to_csv(
                filedirec + os.sep + fold + os.sep + direc + os.sep + "Results.csv",
                index=False
            )

exit(0)