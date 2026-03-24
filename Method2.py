import os
import pandas as pd
import numpy as np
import sys
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

output_file = filedirec + os.sep + "SBLUPGWAS.ma"
df_transformed.to_csv(output_file, sep="\t", index=False)


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


def read_hsq(filepath):
    h2 = np.nan
    h2_se = np.nan
    try:
        tempdata = pd.read_csv(filepath, sep="\t")
        row = tempdata[tempdata["Source"] == "V(G)/Vp"]
        if len(row) > 0:
            h2 = row["Variance"].values[0]
            h2_se = row["SE"].values[0]
    except Exception as e:
        print("Could not read hsq file:", filepath, str(e))
    return h2, h2_se


def transform_gcta_data(traindirec, newtrainfilename, models, p, clumpprune,
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
    covariate_train.fillna(0, inplace=True)
    covariate_train.to_csv(traindirec + os.sep + trainfilename + ".cov",
                           sep="\t", index=False)

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

    if clumpprune == "yes":
        BFILE = traindirec + os.sep + newtrainfilename + ".clumped.pruned"
    else:
        BFILE = traindirec + os.sep + newtrainfilename

    variants = len(pd.read_csv(BFILE + ".bim", sep=r"\s+", header=None))

    tempdata = np.nan
    tempdata_se = np.nan

    hsq_file = traindirec + os.sep + "train_data.hsq"
    hsq_file_alt = traindirec + os.sep + "train_dataa.hsq"

    if models == "GCTA_genotype":
        if os.path.exists(hsq_file):
            os.remove(hsq_file)
        subprocess.run([
            './gcta', '--bfile', BFILE,
            '--make-grm', '--out', traindirec + os.sep + trainfilename
        ])
        subprocess.run([
            './gcta',
            '--grm', traindirec + os.sep + trainfilename,
            '--pheno', traindirec + os.sep + trainfilename + ".PHENO",
            '--reml',
            '--out', traindirec + os.sep + trainfilename
        ])
        if os.path.exists(hsq_file):
            tempdata, tempdata_se = read_hsq(hsq_file)
        else:
            print("Heritibility not found using '--reml'")
            subprocess.run([
                './gcta',
                '--grm', traindirec + os.sep + trainfilename,
                '--pheno', traindirec + os.sep + trainfilename + ".PHENO",
                '--reml',
                '--grm-cutoff', str(0.01),
                '--out', traindirec + os.sep + trainfilename + "a"
            ])
            if os.path.exists(hsq_file_alt):
                tempdata, tempdata_se = read_hsq(hsq_file_alt)
        print(tempdata)

    if models == "GCTA_genotype_covariate":
        if os.path.exists(hsq_file):
            os.remove(hsq_file)
        subprocess.run([
            './gcta', '--bfile', BFILE,
            '--make-grm', '--out', traindirec + os.sep + trainfilename
        ])
        subprocess.run([
            './gcta',
            '--grm', traindirec + os.sep + trainfilename,
            '--pheno', traindirec + os.sep + trainfilename + ".PHENO",
            '--qcovar', traindirec + os.sep + trainfilename + ".cov",
            '--reml',
            '--out', traindirec + os.sep + trainfilename
        ])
        if os.path.exists(hsq_file):
            tempdata, tempdata_se = read_hsq(hsq_file)
        else:
            print("Heritibility not found using '--reml'")
            subprocess.run([
                './gcta',
                '--grm', traindirec + os.sep + trainfilename,
                '--pheno', traindirec + os.sep + trainfilename + ".PHENO",
                '--qcovar', traindirec + os.sep + trainfilename + ".cov",
                '--grm-cutoff', str(0.01),
                '--reml',
                '--out', traindirec + os.sep + trainfilename + "a"
            ])
            if os.path.exists(hsq_file_alt):
                tempdata, tempdata_se = read_hsq(hsq_file_alt)

    if models == "GCTA_genotype_covariate_pca":
        if os.path.exists(hsq_file):
            os.remove(hsq_file)
        subprocess.run([
            './gcta', '--bfile', BFILE,
            '--make-grm', '--out', traindirec + os.sep + trainfilename
        ])
        subprocess.run([
            './gcta',
            '--grm', traindirec + os.sep + trainfilename,
            '--pheno', traindirec + os.sep + trainfilename + ".PHENO",
            '--qcovar', traindirec + os.sep + trainfilename + ".COV_PCA",
            '--reml',
            '--out', traindirec + os.sep + trainfilename
        ])
        if os.path.exists(hsq_file):
            tempdata, tempdata_se = read_hsq(hsq_file)
        else:
            print("Heritibility not found using '--reml'")
            subprocess.run([
                './gcta',
                '--grm', traindirec + os.sep + trainfilename,
                '--pheno', traindirec + os.sep + trainfilename + ".PHENO",
                '--qcovar', traindirec + os.sep + trainfilename + ".COV_PCA",
                '--reml',
                '--grm-cutoff', str(0.01),
                '--out', traindirec + os.sep + trainfilename + "a"
            ])
            if os.path.exists(hsq_file_alt):
                tempdata, tempdata_se = read_hsq(hsq_file_alt)

    print("Heritability:", tempdata, "SE:", tempdata_se, "Variants:", variants)

    global prs_result
    new_row = pd.DataFrame([{
        "clump_p1":       c1_val,
        "clump_r2":       c2_val,
        "clump_kb":       c3_val,
        "p_window_size":  p1_val,
        "p_slide_size":   p2_val,
        "p_LD_threshold": p3_val,
        "numberofpca":    p,
        "h2model":        models,
        "h2":             tempdata,
        "h2_se":          tempdata_se,
        "numberofvariants": variants,
        "clumpprune":     clumpprune
    }])
    prs_result = pd.concat([prs_result, new_row], ignore_index=True)

    outdir = traindirec + os.sep + Name
    create_directory(outdir)
    prs_result.to_csv(outdir + os.sep + "Results.csv", index=False)


h2models = ["GCTA_genotype", "GCTA_genotype_covariate", "GCTA_genotype_covariate_pca"]
clumpprunes = ["yes", "no"]

create_directory(folddirec + os.sep + "GCTA-Heritability")

for p1_val in p_window_size:
    for p2_val in p_slide_size:
        for p3_val in p_LD_threshold:
            for c1_val in clump_p1:
                for c2_val in clump_r2:
                    for c3_val in clump_kb:
                        for p in numberofpca:
                            for model in h2models:
                                for clumpprune in clumpprunes:
                                    transform_gcta_data(
                                        folddirec, newtrainfilename, model, p, clumpprune,
                                        str(p1_val), str(p2_val), str(p3_val),
                                        str(c1_val), str(c2_val), str(c3_val),
                                        "GCTA-Heritability"
                                    )

exit(0)