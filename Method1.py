import os
import sys
import pandas as pd
import numpy as np
import subprocess

filedirec = sys.argv[1]
foldnumber = sys.argv[2]

def create_directory(directory):
    if not os.path.exists(directory):
        os.makedirs(directory)
    return directory

GWAS = filedirec + os.sep + filedirec + ".gz"
df = pd.read_csv(GWAS, compression="gzip", sep="\s+")

if "BETA" in df.columns.to_list():
    df = df[['CHR', 'BP', 'SNP', 'A1', 'A2', 'N', 'SE', 'P', 'BETA', 'INFO', 'MAF']]
else:
    df["BETA"] = np.log(df["OR"])
    df = df[['CHR', 'BP', 'SNP', 'A1', 'A2', 'N', 'SE', 'P', 'BETA', 'INFO', 'MAF']]

df.to_csv(filedirec + os.sep + filedirec + ".txt", sep="\t", index=False)

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


def transform_ldpred2_data(traindirec, newtrainfilename, models, p, clumpprune,
                            p1_val, p2_val, p3_val, c1_val, c2_val, c3_val, Name):

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
        sep="\s+", header=None
    )
    phenotype = tempphenotype_train[[0, 1, 5]]
    phenotype.to_csv(traindirec + os.sep + trainfilename + ".PHENO",
                     sep="\t", header=['FID', 'IID', 'PHENO'], index=False)

    pcs_train = pd.read_table(
        traindirec + os.sep + trainfilename + ".eigenvec",
        sep="\s+", header=None,
        names=["FID", "IID"] + [f"PC{str(i)}" for i in range(1, int(p) + 1)]
    )
    covariate_train = pd.read_table(traindirec + os.sep + trainfilename + ".cov", sep="\s+")
    covariate_train.fillna(0, inplace=True)
    covariate_train = covariate_train[
        covariate_train["FID"].isin(pcs_train["FID"].values) &
        covariate_train["IID"].isin(pcs_train["IID"].values)
    ]
    covariate_train['FID'] = covariate_train['FID'].astype(str)
    pcs_train['FID'] = pcs_train['FID'].astype(str)
    covariate_train['IID'] = covariate_train['IID'].astype(str)
    pcs_train['IID'] = pcs_train['IID'].astype(str)
    covandpcs_train = pd.merge(covariate_train, pcs_train, on=["FID", "IID"])
    covandpcs_train.to_csv(traindirec + os.sep + trainfilename + ".COV_PCA",
                           sep="\t", index=False)

    tempdata1 = np.nan
    tempdata1_se = np.nan
    tempdata2 = np.nan

    if models == "LDpred-2_full":
        if clumpprune == "yes":
            os.system("Rscript LDpred-2-Heritability.R " + os.path.join(filedirec) +
                      "  " + traindirec + " " + trainfilename + " " +
                      newtrainfilename + ".clumped.pruned" +
                      " " + "3" + " " + c3_val + " " + c1_val + " " + c2_val + " " + p)
        else:
            os.system("Rscript LDpred-2-Heritability.R " + os.path.join(filedirec) +
                      "  " + traindirec + " " + trainfilename + " " + newtrainfilename +
                      " " + "3" + " " + c3_val + " " + c1_val + " " + c2_val + " " + p)

        h2_file    = traindirec + os.sep + "ldpred_h2_full.txt"
        h2_se_file = traindirec + os.sep + "ldpred_h2_full_se.txt"
        var_file   = traindirec + os.sep + "ldpred_h2_full_variants.txt"

        if os.path.exists(h2_file):
            tempdata1    = pd.read_csv(h2_file, sep="\s+", header=None)[1].values[0]
        if os.path.exists(h2_se_file):
            tempdata1_se = pd.read_csv(h2_se_file, sep="\s+", header=None)[1].values[0]
        if os.path.exists(var_file):
            tempdata2    = pd.read_csv(var_file, sep="\s+", header=None)[1].values[0]

    if models == "LDpred-2_hapmap":
        if clumpprune == "yes":
            os.system("Rscript LDpred-2-Heritability.R " + os.path.join(filedirec) +
                      "  " + traindirec + " " + trainfilename + " " +
                      newtrainfilename + ".clumped.pruned" +
                      " " + "2" + " " + c3_val + " " + c1_val + " " + c2_val + " " + p)
        else:
            os.system("Rscript LDpred-2-Heritability.R " + os.path.join(filedirec) +
                      "  " + traindirec + " " + trainfilename + " " + newtrainfilename +
                      " " + "2" + " " + c3_val + " " + c1_val + " " + c2_val + " " + p)

        h2_file    = traindirec + os.sep + "ldpred_h2_hapmap.txt"
        h2_se_file = traindirec + os.sep + "ldpred_h2_hapmap_se.txt"
        var_file   = traindirec + os.sep + "ldpred_h2_hapmap_variants.txt"

        if os.path.exists(h2_file):
            tempdata1    = pd.read_csv(h2_file, sep="\s+", header=None)[1].values[0]
        if os.path.exists(h2_se_file):
            tempdata1_se = pd.read_csv(h2_se_file, sep="\s+", header=None)[1].values[0]
        if os.path.exists(var_file):
            tempdata2    = pd.read_csv(var_file, sep="\s+", header=None)[1].values[0]

    global prs_result
    new_row = pd.DataFrame([{
        "clump_p1":        c1_val,
        "clump_r2":        c2_val,
        "clump_kb":        c3_val,
        "p_window_size":   p1_val,
        "p_slide_size":    p2_val,
        "p_LD_threshold":  p3_val,
        "numberofpca":     p,
        "h2":              tempdata1,
        "h2_se":           tempdata1_se,
        "h2model":         models,
        "clumpprune":      clumpprune,
        "numberofvariants": tempdata2,
    }])
    prs_result = pd.concat([prs_result, new_row], ignore_index=True)

    outdir = traindirec + os.sep + Name
    create_directory(outdir)
    prs_result.to_csv(outdir + os.sep + "Results.csv", index=False)

    print("Heritability Value:", tempdata1)
    print("Heritability SE:", tempdata1_se)
    print("Heritability Variants:", tempdata2)
    return


h2models = ["LDpred-2_full", "LDpred-2_hapmap"]
clumpprunes = ["yes", "no"]

create_directory(folddirec + os.sep + "LDpred-2-Heritability")

for p1_val in p_window_size:
    for p2_val in p_slide_size:
        for p3_val in p_LD_threshold:
            for c1_val in clump_p1:
                for c2_val in clump_r2:
                    for c3_val in clump_kb:
                        for p in numberofpca:
                            for model in h2models:
                                for clumpprune in clumpprunes:
                                    transform_ldpred2_data(
                                        folddirec, newtrainfilename, model, p, clumpprune,
                                        str(p1_val), str(p2_val), str(p3_val),
                                        str(c1_val), str(c2_val), str(c3_val),
                                        "LDpred-2-Heritability"
                                    )

exit(0)