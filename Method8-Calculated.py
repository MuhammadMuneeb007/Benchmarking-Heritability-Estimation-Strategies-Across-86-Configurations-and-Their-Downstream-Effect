import os
import sys
import pandas as pd
import numpy as np
import subprocess
import shutil
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

df_transformed = pd.DataFrame({
    'Predictor': df['SNP'],
    'A1': df['A1'],
    'A2': df['A2'],
    'n': df['N'],
    'Z': df['BETA'] / df['SE'],
    'SNP': df['SNP']
})
print(df_transformed.head())

columns_to_check = ['A1', 'A2']

def specified_single_characters(row, cols):
    return all(len(str(row[col])) == 1 for col in cols)

df_transformed = df_transformed[
    df_transformed.apply(lambda row: specified_single_characters(row, columns_to_check), axis=1)
]
df_transformed.to_csv(filedirec + os.sep + filedirec + ".ldak", sep="\t", index=False)


def create_directory(directory):
    if not os.path.exists(directory):
        os.makedirs(directory)
    return directory


def readpve(direc):
    file_path = direc + os.sep + "snpher1.hers"
    with open(file_path, 'r') as file:
        for line in file:
            if "Her_All" in line:
                pve_estimate = line.split(" ")[1].strip()
                print("Heritability:", pve_estimate)
                return pve_estimate
    return np.nan


def readpve_se(direc):
    file_path = direc + os.sep + "snpher1.hers"
    with open(file_path, 'r') as file:
        for line in file:
            if "Her_All" in line:
                parts = line.split()
                if len(parts) >= 3:
                    pve_se = parts[2].strip()
                    print("Heritability SE:", pve_se)
                    return pve_se
    return np.nan


def readvariants(direc):
    file_path = direc + os.sep + "snpher1.overlap"
    with open(file_path, 'r') as file:
        for line in file:
            if "Summary_Statistic_Predictors" in line:
                pve_estimate = line.split(" ")[1].strip()
                print("Variants:", pve_estimate)
                return pve_estimate
    return np.nan


def delete_previous_results(traindirec):
    for f in ["snpher1.hers", "snpher1.overlap"]:
        fp = traindirec + os.sep + f
        if os.path.exists(fp):
            os.remove(fp)


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


def transform_gemma_data(traindirec, newtrainfilename, p, clumpprune, data,
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
    subprocess.run(command)

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
    subprocess.run(command)

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
    subprocess.run(command)

    command = [
        "./plink",
        "--bfile", traindirec + os.sep + newtrainfilename + ".clumped.pruned",
        "--extract", traindirec + os.sep + trainfilename + ".valid.snp",
        "--pca", p,
        "--out", traindirec + os.sep + trainfilename
    ]
    subprocess.run(command)

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
        names=["FID", "IID"] + [f"PC{str(i)}" for i in range(1, int(p) + 1)]
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

    global prs_result

    if data == "human":
        delete_previous_results(traindirec)
        subprocess.run([
            "./ldak",
            "--calc-tagging", traindirec + os.sep + "HumDef",
            "--bfile", BFILE,
            "--power", "-.25"
        ])
        subprocess.run([
            "./ldak",
            "--sum-hers", traindirec + os.sep + "snpher1",
            "--summary", filedirec + os.sep + filedirec + ".ldak",
            "--tagfile", traindirec + os.sep + "HumDef.tagging",
            "--check-sums", "NO"
        ])
        pve_estimate = readpve(traindirec)
        pve_se = readpve_se(traindirec)
        variants = readvariants(traindirec)
        new_row = pd.DataFrame([{
            "clump_p1": c1_val, "clump_r2": c2_val, "clump_kb": c3_val,
            "p_window_size": p1_val, "p_slide_size": p2_val, "p_LD_threshold": p3_val,
            "numberofpca": p, "h2": pve_estimate, "h2_se": pve_se,
            "h2model": "LDAK_human_" + clumpprune,
            "clumpprune": clumpprune, "numberofvariants": variants
        }])
        prs_result = pd.concat([prs_result, new_row], ignore_index=True)

    if data == "GCTA":
        delete_previous_results(traindirec)
        subprocess.run([
            "./ldak",
            "--calc-tagging", traindirec + os.sep + "HumDef",
            "--bfile", BFILE,
            "--power", "-1"
        ])
        subprocess.run([
            "./ldak",
            "--sum-hers", traindirec + os.sep + "snpher1",
            "--summary", filedirec + os.sep + filedirec + ".ldak",
            "--tagfile", traindirec + os.sep + "HumDef.tagging",
            "--check-sums", "NO"
        ])
        pve_estimate = readpve(traindirec)
        pve_se = readpve_se(traindirec)
        variants = readvariants(traindirec)
        new_row = pd.DataFrame([{
            "clump_p1": c1_val, "clump_r2": c2_val, "clump_kb": c3_val,
            "p_window_size": p1_val, "p_slide_size": p2_val, "p_LD_threshold": p3_val,
            "numberofpca": p, "h2": pve_estimate, "h2_se": pve_se,
            "h2model": "LDAK_GCTA_" + clumpprune,
            "clumpprune": clumpprune, "numberofvariants": variants
        }])
        prs_result = pd.concat([prs_result, new_row], ignore_index=True)

    if data == "BLD-LDAK":
        delete_previous_results(traindirec)
        subprocess.run([
            './ldak', '--cut-weights', traindirec + os.sep + 'sections',
            '--bfile', BFILE
        ])
        subprocess.run([
            './ldak', '--calc-weights-all', traindirec + os.sep + 'sections',
            '--bfile', BFILE
        ])
        shutil.move(
            traindirec + os.sep + 'sections' + os.sep + 'weights.short',
            "LDAKFILES" + os.sep + 'bld65'
        )
        subprocess.run([
            './ldak', '--calc-tagging', traindirec + os.sep + "HumDef",
            '--bfile', BFILE, '--power', '-.25',
            '--annotation-number', '65',
            '--annotation-prefix', "LDAKFILES" + os.sep + 'bld'
        ])
        subprocess.run([
            "./ldak",
            "--sum-hers", traindirec + os.sep + "snpher1",
            "--summary", filedirec + os.sep + filedirec + ".ldak",
            "--tagfile", traindirec + os.sep + "HumDef.tagging",
            "--check-sums", "NO"
        ])
        pve_estimate = readpve(traindirec)
        pve_se = readpve_se(traindirec)
        variants = readvariants(traindirec)
        new_row = pd.DataFrame([{
            "clump_p1": c1_val, "clump_r2": c2_val, "clump_kb": c3_val,
            "p_window_size": p1_val, "p_slide_size": p2_val, "p_LD_threshold": p3_val,
            "numberofpca": p, "h2": pve_estimate, "h2_se": pve_se,
            "h2model": "LDAK_BLD-LDAK_" + clumpprune,
            "clumpprune": clumpprune, "numberofvariants": variants
        }])
        prs_result = pd.concat([prs_result, new_row], ignore_index=True)

    if data == "alpha":
        for j in range(-2, 2):
            alpha = -j / 2.0
            print(alpha)
            delete_previous_results(traindirec)
            subprocess.run([
                './ldak', '--calc-tagging', traindirec + os.sep + 'Alpha' + str(j),
                '--bfile', BFILE, '--power', str(alpha)
            ])
            subprocess.run([
                "./ldak",
                "--sum-hers", traindirec + os.sep + "snpher1",
                "--summary", filedirec + os.sep + filedirec + ".ldak",
                "--tagfile", traindirec + os.sep + 'Alpha' + str(j) + ".tagging",
                "--check-sums", "NO"
            ])
            pve_estimate = readpve(traindirec)
            pve_se = readpve_se(traindirec)
            variants = readvariants(traindirec)
            new_row = pd.DataFrame([{
                "clump_p1": c1_val, "clump_r2": c2_val, "clump_kb": c3_val,
                "p_window_size": p1_val, "p_slide_size": p2_val, "p_LD_threshold": p3_val,
                "numberofpca": p, "h2": pve_estimate, "h2_se": pve_se,
                "h2model": "LDAK_alpha_" + str(alpha) + "_" + clumpprune,
                "clumpprune": clumpprune, "numberofvariants": variants
            }])
            prs_result = pd.concat([prs_result, new_row], ignore_index=True)

    print(prs_result)
    outdir = traindirec + os.sep + Name
    create_directory(outdir)
    prs_result.to_csv(outdir + os.sep + "Results.csv", index=False)
    return


datas = ["human", "GCTA", "BLD-LDAK", "alpha"]
clumpprunes = ["yes"]

create_directory(folddirec + os.sep + "LDAK-Calculated-Heritability")

for p1_val in p_window_size:
    for p2_val in p_slide_size:
        for p3_val in p_LD_threshold:
            for c1_val in clump_p1:
                for c2_val in clump_r2:
                    for c3_val in clump_kb:
                        for p in numberofpca:
                            for clumpprune in clumpprunes:
                                for data in datas:
                                    transform_gemma_data(
                                        folddirec, newtrainfilename, p, clumpprune, data,
                                        str(p1_val), str(p2_val), str(p3_val),
                                        str(c1_val), str(c2_val), str(c3_val),
                                        "LDAK-Calculated-Heritability"
                                    )

exit(0)