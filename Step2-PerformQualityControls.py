import os  # Importing the os module for operating system functionalities
import pandas as pd  # Importing pandas for data manipulation and analysis
import subprocess  # Importing subprocess module to spawn new processes
import sys  # Importing sys module for system-specific parameters and functions
import os
import pandas as pd
from sklearn.model_selection import StratifiedKFold
import subprocess
from sklearn.model_selection import KFold, cross_val_score

# Set the directory where the files are located and process one phenotype
filedirec = sys.argv[1]

# Define file paths for different data files
BED = filedirec + os.sep + filedirec
BIM = filedirec + os.sep + filedirec + ".bim"
FAM = filedirec + os.sep + filedirec + ".fam"
COV = filedirec + os.sep + filedirec + ".cov"
Height = filedirec + os.sep + filedirec + ".height"
GWAS = filedirec + os.sep + filedirec + ".gz"

# Read the BIM file and print its columns and first 10 rows
bimfile = pd.read_csv(BIM, sep="\s+", header=None)
print("Columns of BIM file:")
print(bimfile.columns)
print("First 10 rows of BIM file:")
print(bimfile.head(10))

# Create a 'match' column in the BIM file by concatenating certain columns
bimfile["match"] = bimfile[0].astype(str) + "_" + bimfile[3].astype(str) + "_" + bimfile[4].astype(str) + "_" + bimfile[5].astype(str)

# Read the first 10 rows of the FAM file and print its columns and first 10 rows
temp = pd.read_csv(FAM, sep="\s+", header=None, nrows=10)
print("Columns of FAM file:")
print(temp.columns)
print("First 10 rows of FAM file:")
print(temp.head(10))

# Read the first 10 rows of the COV file and print its columns and first 10 rows
temp = pd.read_csv(COV, sep="\s+", nrows=10)
print("Columns of COV file:")
print(temp.columns)
print("First 10 rows of COV file:")
print(temp.head(10))

# Read the first 10 rows of the Height file and print its columns and first 10 rows
temp = pd.read_csv(Height, sep="\s+", nrows=10)
print("Columns of Height file:")
print(temp.columns)
print("First 10 rows of Height file:")
print(temp.head(10))

# Read the first 10 rows of the GWAS file and print its columns and first 10 rows
temp = pd.read_csv(GWAS, sep="\s+", nrows=10)
print("Columns of GWAS file:")
print(temp.columns)
print("First 10 rows of GWAS file:")
print(temp.head(10))

# Read GWAS data from a compressed file using pandas and print the first 10 rows
temp = pd.read_csv(GWAS, compression="gzip", sep="\s+", nrows=10)
print(temp.head(10))

# Read GWAS data from a compressed file using pandas
df = pd.read_csv(GWAS, compression="gzip", sep="\s+", on_bad_lines='skip')

# Display the initial number of rows in the dataframe
print("Initial number of SNPs:", len(df))

# Apply quality control steps: Filter SNPs based on Minor Allele Frequency (MAF) and Imputation Information Score (INFO)
#df = df.loc[(df['MAF'] > 0.01) & (df['INFO'] > 0.8)]

# Display the number of rows after applying the filters
print("Number of SNPs after quality control:", len(df))

# Remove duplicate SNPs based on the 'SNP' column
# Do not use the following command. If the column contains "X" for all SNPs, all the SNPs will be removed.
#df = df.drop_duplicates(subset='SNP')

# Display the number of rows after removing duplicate SNPs
#print("SNPs in GWAS after removing duplicate SNPs:", len(df))

# Remove ambiguous SNPs with complementary alleles (C/G or A/T) to avoid potential errors
df = df[~((df['A1'] == 'A') & (df['A2'] == 'T') |
          (df['A1'] == 'T') & (df['A2'] == 'A') |
          (df['A1'] == 'G') & (df['A2'] == 'C') |
          (df['A1'] == 'C') & (df['A2'] == 'G'))]

# Display the final number of SNPs after removing ambiguous SNPs
print("Final number of SNPs after removing ambiguous SNPs:", len(df))

# Create a 'match' column in the dataframe by concatenating certain columns and remove duplicate rows based on this column
df["match"] = df["CHR"].astype(str) + "_" + df["BP"].astype(str) + "_" + df["A2"].astype(str) + "_" + df["A1"].astype(str)
df.drop_duplicates(subset='match', inplace=True)
del df["match"]

# Print the number of rows after removing SNPs for which any row does not contain the required value
print("Removing SNPs for which even a single row does not contain the required value:", len(df))

# Check if all RSIDs are missing and handle accordingly
if (df['SNP'] == 'X').all():
    print("RSIDs are missing!")
    bimfile = pd.read_csv(filedirec + os.sep + filedirec + ".bim", sep="\s+", header=None)
    bimfile["match"] = bimfile[0].astype(str) + "_" + bimfile[3].astype(str) + "_" + bimfile[4].astype(str) + "_" + bimfile[5].astype(str)
    df["match"] = df["CHR"].astype(str) + "_" + df["BP"].astype(str) + "_" + df["A2"].astype(str) + "_" + df["A1"].astype(str)
    
    df.drop_duplicates(subset='match', inplace=True)
    bimfile.drop_duplicates(subset='match', inplace=True)

    df = df[df['match'].isin(bimfile['match'].values)]
    bimfile = bimfile[bimfile['match'].isin(df['match'].values)]
    df = df[df['match'].isin(bimfile['match'].values)]
    bimfile = bimfile[bimfile['match'].isin(df['match'].values)]

    df.drop_duplicates(subset='match', inplace=True)
    bimfile.drop_duplicates(subset='match', inplace=True)    

    df["SNP"] = bimfile[1].values
    print("match", len(df))
    del df["match"]
    df.to_csv(GWAS, compression="gzip", sep="\t", index=None)   

    pass
else:
    df.drop_duplicates(subset='SNP', inplace=True)
    df.to_csv(GWAS, compression="gzip", sep="\t", index=None)
    print("RSID is present!")
    pass



# ## Individual genotype data (Target Data) Processing
# 
# 
# Ensure that the phenotype file, FAM file, and covariate file contain an identical number of samples. Remove any missing samples based on your data. Note that the extent of missingness in phenotypes and covariates may vary.
# 
# 
# **Note:** Plink needs to be installed or placed in the same directory as this notebook.
# 
# [Download Plink](https://www.cog-genomics.org/plink/)
# 
# We recommend using Linux. In cases where Windows is required due to package installation issues on Linux, we provide the following guidance:
# 
# 1. For Windows, use `plink`.
# 2. For Linux, use `./plink`.
# 

# In[16]:

# New files to be saved with QC suffix
newfilename = filedirec + "_QC"

# Read information from FAM file
f = pd.read_csv(FAM, header=None, sep="\s+", names=["FID", "IID", "Father", "Mother", "Sex", "Phenotype"])
print("FAM file contents:")
print(f.head())
print("Total number of people in FAM file:", len(f))

# Read height information from a separate file
h = pd.read_csv(Height, sep="\t")
# Drop rows with missing height values
h = h.dropna(subset=["Height"])
print("Phenotype information is available for:", len(h), "people")
print(len(h))
# Merge FAM and height data on common columns FID and IID
result = pd.merge(f, h, on=['FID', 'IID'])

# Replace 'Phenotype' column with 'Height' and save to a new PeopleWithPhenotype.txt file
result["Phenotype"] = result["Height"].values
del result["Height"]
print(result)
result.to_csv(filedirec + os.sep + "PeopleWithPhenotype.txt", index=False, header=False, sep="\t")

# Use PLINK to keep only the people with phenotype present
plink_command = [
    './plink',
    '--bfile', filedirec + os.sep + filedirec,
    '--keep', filedirec + os.sep + "PeopleWithPhenotype.txt",
    '--make-bed',
    '--out', filedirec + os.sep + newfilename
]
subprocess.run(plink_command)

# Update the phenotype information in the new FAM file
f = pd.read_csv(filedirec + os.sep + newfilename + ".fam", header=None, sep="\s+",
                names=["FID", "IID", "Father", "Mother", "Sex", "Phenotype"])
f["Phenotype"] = result["Phenotype"].values
f.to_csv(filedirec + os.sep + newfilename + ".fam", index=False, header=False, sep="\t")

# Update the covariate file as well
covfile = filedirec + os.sep + filedirec + '.cov'
covfile = pd.read_csv(covfile, sep="\s+")

print("Covariate file contents:")
print(covfile.head())
print("Total number of people in Covariate file:", len(covfile))

# Match the FID and IID from covariate and height file
covfile = covfile[covfile['FID'].isin(f["FID"].values) & covfile['IID'].isin(f["IID"].values)]
print("Covariate file contents after matching with FAM file:")
print(covfile.head())
print("Total number of people in Covariate file after matching:", len(covfile))
covfile.to_csv(filedirec + os.sep + newfilename + ".cov", index=None, sep="\t")

# Read the Fam file
input_file_path = filedirec + os.sep + newfilename + '.fam'
df = pd.read_csv(input_file_path, sep="\s+", header=None)

# Create 5 directories for storing fold information
output_directory_base = filedirec
os.makedirs(output_directory_base, exist_ok=True)

# Split the data into 5 folds using Stratified K-Folds

# The following code is for continous phenotype.
fold_column = 5  # fifth column contains phenotypes
kf = KFold(n_splits=5, shuffle=True, random_state=42)
phenotype_col = 5

# The following code is for binary phenotype.
skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
column_values = df[phenotype_col].unique()

# Check if the phenotype is binary or continoous.

if set(column_values) == {1, 2}:
    print("The column contains only 0 and 1.")
    #exit(0)
    for fold_id, (train_index, test_index) in enumerate(skf.split(df, df[phenotype_col])):
        fold_directory = os.path.join(output_directory_base, f'Fold_{fold_id}')

        train_file_name = "train_data"
        test_file_name = "test_data"

        new_train_file_name = "train_data.QC"
        new_test_file_name = "test_data.QC"

        os.makedirs(fold_directory, exist_ok=True)
        #"""
        # Save train and test data to separate CSV files
        train_data = df.iloc[train_index]
        test_data = df.iloc[test_index]
        
        train_data.to_csv(os.path.join(fold_directory, 'train_data.fam'),sep="\t",header=False,index=False)
        test_data.to_csv(os.path.join(fold_directory, 'test_data.fam'),sep="\t",header=False,index=False)
        #print(train_data)
        #exit(0)
        
        #exit(0)
        # Step 4: Use PLINK to extract test and train samples for each fold
        plink_train_command = [
            './plink',
            '--bfile', filedirec+os.sep+newfilename,
            '--keep', os.path.join(fold_directory, train_file_name+'.fam'),
            '--make-bed',
            '--out', os.path.join(fold_directory, train_file_name)
        ]

        plink_test_command = [
            './plink',
            '--bfile', filedirec+os.sep+newfilename,
            '--keep', os.path.join(fold_directory, test_file_name+'.fam'),
            '--make-bed',
            '--out', os.path.join(fold_directory, test_file_name)
        ]

        subprocess.run(plink_train_command)
        subprocess.run(plink_test_command)
        
        covfile = filedirec+os.sep+newfilename+'.cov'
        covfile = pd.read_csv(covfile)

        cov_train_data = covfile.iloc[train_index]
        cov_test_data = covfile.iloc[test_index]

        cov_train_data.to_csv(os.path.join(fold_directory, train_file_name+'.cov'),sep=",",index=False)
        cov_test_data.to_csv(os.path.join(fold_directory, test_file_name+'.cov'),sep=",",index=False)

        #exit(0)
        ### perform Quality controls on the training data only.
        print(os.path.join(fold_directory, train_file_name))
        plink_command_1 = [
            './plink',
            '--bfile', os.path.join(fold_directory, train_file_name),
            '--maf', '0.01',
            '--hwe', '1e-6',
            '--geno', '0.1',
            '--mind', '0.1',
            '--write-snplist',
            '--make-just-fam',
            '--out', os.path.join(fold_directory, new_train_file_name)
        ]

        subprocess.run(plink_command_1)
        
        # Command 2
        plink_command_2 = [
            './plink',
            '--bfile', os.path.join(fold_directory, train_file_name),
            '--keep', os.path.join(fold_directory, new_train_file_name+'.fam'),
            '--extract', os.path.join(fold_directory, new_train_file_name+'.snplist'),
            '--indep-pairwise', '200', '50', '0.25',
            '--out', os.path.join(fold_directory, new_train_file_name)
        ]

        #subprocess.run(plink_command_2)
        
        # Command 3
        plink_command_3 = [
            './plink',
            '--bfile', os.path.join(fold_directory, train_file_name),
            #'--extract', os.path.join(fold_directory, new_train_file_name+'.prune.in'),
            '--keep', os.path.join(fold_directory, new_train_file_name+'.fam'),
            #'--extract', os.path.join(fold_directory, new_train_file_name+'.snplist'),
            #'--keep', os.path.join(fold_directory, new_train_file_name+'.fam'),
            '--het',
            '--out', os.path.join(fold_directory, new_train_file_name)
        ]

        subprocess.run(plink_command_3)
        # Invoked R functions.
        
        os.system("Rscript Module1.R "+os.path.join(fold_directory)+"  "+train_file_name+" "+new_train_file_name+ " "+"1")
        
        
        #exit(0)
        
        # Code for sex check: The sample data have 22 chromosome so this operation is skipped.
        #"""
        plink_command = [
        './plink',
        '--bfile', os.path.join(fold_directory, train_file_name),
        #'--extract', os.path.join(fold_directory, 'train_data.QC.prune.in'),
        '--extract', os.path.join(fold_directory, new_train_file_name+'.snplist'),
        '--keep', os.path.join(fold_directory, train_file_name+'.valid.sample'),
        '--check-sex',
        '--out', os.path.join(fold_directory, new_train_file_name)
        ]
        
        # Invoke the PLINK command using subprocess for sex check
        #subprocess.run(plink_command)
        #os.system("Rscript Module1.R "+os.path.join(fold_directory)+"  "+train_file_name+" "+new_train_file_name+ " "+"2")
        
        #"""

        plink_command_1 = [
        './plink',
        '--bfile', os.path.join(fold_directory, train_file_name),
        #'--extract', os.path.join(fold_directory, new_train_file_name+'.prune.in'),
        # If you use sex check then use the following line. Otherwise it is not required.
        '--extract', os.path.join(fold_directory, new_train_file_name+'.snplist'),
        #'--keep', os.path.join(fold_directory, 'train_data.QC.valid'),
        '--rel-cutoff', '0.125',
        '--out', os.path.join(fold_directory, new_train_file_name)
        ]
        
        subprocess.run(plink_command_1)
   
        plink_command_2 = [
            './plink',
            '--bfile', os.path.join(fold_directory, train_file_name),
            '--make-bed',
            '--keep', os.path.join(fold_directory, new_train_file_name+'.rel.id'),
            '--out', os.path.join(fold_directory, new_train_file_name),
            #'--extract', os.path.join(fold_directory, new_train_file_name+'.snplist'),
            # Use the following command for Predictor-based matching! We relied on the heritability tools for matching the GWAS and genotype.
            #'--exclude', os.path.join(fold_directory, train_file_name+'.mismatch'),
            '--a1-allele', os.path.join(fold_directory, train_file_name+'.a1')
        ]

        subprocess.run(plink_command_2)
        #exit(0)
else:
    print("The column does not contain only 0 and 1.")
    for fold_id, (train_index, test_index) in enumerate(kf.split(df, df[phenotype_col])):
        fold_directory = os.path.join(output_directory_base, f'Fold_{fold_id}')

        train_file_name = "train_data"
        test_file_name = "test_data"

        new_train_file_name = "train_data.QC"
        new_test_file_name = "test_data.QC"

        os.makedirs(fold_directory, exist_ok=True)
        #"""
        # Save train and test data to separate CSV files
        train_data = df.iloc[train_index]
        test_data = df.iloc[test_index]
        
        train_data.to_csv(os.path.join(fold_directory, 'train_data.fam'),sep="\t",header=False,index=False)
        test_data.to_csv(os.path.join(fold_directory, 'test_data.fam'),sep="\t",header=False,index=False)
        #print(train_data)
        #exit(0)
        
        #exit(0)
        # Step 4: Use PLINK to extract test and train samples for each fold
        plink_train_command = [
            './plink',
            '--bfile', filedirec+os.sep+newfilename,
            '--keep', os.path.join(fold_directory, train_file_name+'.fam'),
            '--make-bed',
            '--out', os.path.join(fold_directory, train_file_name)
        ]

        plink_test_command = [
            './plink',
            '--bfile', filedirec+os.sep+newfilename,
            '--keep', os.path.join(fold_directory, test_file_name+'.fam'),
            '--make-bed',
            '--out', os.path.join(fold_directory, test_file_name)
        ]

        subprocess.run(plink_train_command)
        subprocess.run(plink_test_command)
        
        covfile = filedirec+os.sep+newfilename+'.cov'
        covfile = pd.read_csv(covfile)

        cov_train_data = covfile.iloc[train_index]
        cov_test_data = covfile.iloc[test_index]

        cov_train_data.to_csv(os.path.join(fold_directory, train_file_name+'.cov'),sep=",",index=False)
        cov_test_data.to_csv(os.path.join(fold_directory, test_file_name+'.cov'),sep=",",index=False)

        #exit(0)
        ### perform Quality controls on the training data only.
        print(os.path.join(fold_directory, train_file_name))
        plink_command_1 = [
            './plink',
            '--bfile', os.path.join(fold_directory, train_file_name),
            '--maf', '0.01',
            '--hwe', '1e-6',
            '--geno', '0.1',
            '--mind', '0.1',
            '--write-snplist',
            '--make-just-fam',
            '--out', os.path.join(fold_directory, new_train_file_name)
        ]

        subprocess.run(plink_command_1)
        
        # Command 2
        plink_command_2 = [
            './plink',
            '--bfile', os.path.join(fold_directory, train_file_name),
            '--keep', os.path.join(fold_directory, new_train_file_name+'.fam'),
            '--extract', os.path.join(fold_directory, new_train_file_name+'.snplist'),
            '--indep-pairwise', '200', '50', '0.25',
            '--out', os.path.join(fold_directory, new_train_file_name)
        ]

        #subprocess.run(plink_command_2)
        
        # Command 3
        plink_command_3 = [
            './plink',
            '--bfile', os.path.join(fold_directory, train_file_name),
            #'--extract', os.path.join(fold_directory, new_train_file_name+'.prune.in'),
            '--keep', os.path.join(fold_directory, new_train_file_name+'.fam'),
            #'--extract', os.path.join(fold_directory, new_train_file_name+'.snplist'),
            #'--keep', os.path.join(fold_directory, new_train_file_name+'.fam'),
            '--het',
            '--out', os.path.join(fold_directory, new_train_file_name)
        ]

        subprocess.run(plink_command_3)
        # Invoked R functions.
        
        os.system("Rscript Module1.R "+os.path.join(fold_directory)+"  "+train_file_name+" "+new_train_file_name+ " "+"1")
        
        
        #exit(0)
        
        # Code for sex check: The sample data have 22 chromosome so this operation is skipped.
        #"""
        plink_command = [
        './plink',
        '--bfile', os.path.join(fold_directory, train_file_name),
        #'--extract', os.path.join(fold_directory, 'train_data.QC.prune.in'),
        '--extract', os.path.join(fold_directory, new_train_file_name+'.snplist'),
        '--keep', os.path.join(fold_directory, train_file_name+'.valid.sample'),
        '--check-sex',
        '--out', os.path.join(fold_directory, new_train_file_name)
        ]
        
        # Invoke the PLINK command using subprocess for sex check
        #subprocess.run(plink_command)
        #os.system("Rscript Module1.R "+os.path.join(fold_directory)+"  "+train_file_name+" "+new_train_file_name+ " "+"2")
        
        #"""

        plink_command_1 = [
        './plink',
        '--bfile', os.path.join(fold_directory, train_file_name),
        #'--extract', os.path.join(fold_directory, new_train_file_name+'.prune.in'),
        # If you use sex check then use the following line. Otherwise it is not required.
        '--extract', os.path.join(fold_directory, new_train_file_name+'.snplist'),
        #'--keep', os.path.join(fold_directory, 'train_data.QC.valid'),
        '--rel-cutoff', '0.125',
        '--out', os.path.join(fold_directory, new_train_file_name)
        ]
        
        subprocess.run(plink_command_1)
   
        plink_command_2 = [
            './plink',
            '--bfile', os.path.join(fold_directory, train_file_name),
            '--make-bed',
            '--keep', os.path.join(fold_directory, new_train_file_name+'.rel.id'),
            '--out', os.path.join(fold_directory, new_train_file_name),
            #'--extract', os.path.join(fold_directory, new_train_file_name+'.snplist'),
            # Use the following command for Predictor-based matching! We relied on the heritability tools for matching the GWAS and genotype.
            #'--exclude', os.path.join(fold_directory, train_file_name+'.mismatch'),
            '--a1-allele', os.path.join(fold_directory, train_file_name+'.a1')
        ]

        subprocess.run(plink_command_2)
        #exit(0)
 

 

import os

# List of file names to check for existence
f = [
    "train_data.QC.bed",
    "train_data.QC.bim",
    "train_data.QC.fam",
    "train_data.cov",
    
    "test_data.bed",
    "test_data.bim",
    "test_data.fam",
    "test_data.cov"
]

# Specify the fold number
fold_number = 1

# Loop through each file name in the list
for loop in f:
    # Check if the file exists in the specified directory for the given fold
    if os.path.exists(filedirec + os.sep + "Fold_" + str(fold_number) + os.sep + loop):
        # Print a message indicating that the file exists
        print(loop, "Yes, the file exists.")
    else:
        # Print a message indicating that the file does not exist
        print(loop, "No, the file does not exist.")

 