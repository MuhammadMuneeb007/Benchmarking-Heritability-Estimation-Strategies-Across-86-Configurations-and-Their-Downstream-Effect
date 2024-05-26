## Heritability tools

Multiple tools calculate heritability, as shown in the table below. Each tool uses different statistical methods to estimate heritability, and the datasets they use also vary. Some tools use the GWAS summary statistic file, some use genotype data, covariates, and PCA, and some use reference panels and SNP tagging.

| Tool        | URL                                                        |
|-------------|------------------------------------------------------------|
| GEMMA       | [https://github.com/genetics-statistics/GEMMA](https://github.com/genetics-statistics/GEMMA) |
| GCTA  | [http://cnsgenomics.com/software/gcta/#Overview](http://cnsgenomics.com/software/gcta/#Overview) |
| LDAK        | [http://dougspeed.com/ldak/](http://dougspeed.com/ldak/)   |
| DPR         | [https://github.com/biostatpzeng/DPR](https://github.com/biostatpzeng/DPR) |
| LDSC        | [https://github.com/bulik/ldsc](https://github.com/bulik/ldsc) |
| SumHer      | [http://dougspeed.com/sumher/](http://dougspeed.com/sumher/) |



## Purpose of this documentation

In this research, we used various heritability tools and created multiple variants of each method to calculate heritability for 11 phenotypes. Two polygenic risk scores (PRS) tools, LDpred-2 and GCTA, rely on heritability estimates for PRS calculation. We investigated whether the method used to calculate heritability impacts the performance of the PRS tools. Benchmarking all these tools is essential to identify the best method for heritability calculation that optimizes PRS calculation.


## Helper tools
 

| Tool | Description | Link |
|------|-------------|------|
| GWASPokerforPRS | A tool for downloading GWAS data from GWAS Catalog| [https://github.com/MuhammadMuneeb007/GWASPokerforPRS](https://github.com/MuhammadMuneeb007/GWASPokerforPRS) |
| Detect genomic build | Detect the genomic build of a dataset | [https://www.biostars.org/p/9495682/#9595219](https://www.biostars.org/p/9495682/#9595219) |
| pyliftover | A Python package for genomic coordinate conversion | [https://pypi.org/project/pyliftover/](https://pypi.org/project/pyliftover/) |

## Dataset

We analyzed 14 phenotypes from the UK Biobank and downloaded the corresponding GWAS files from the GWAS catalog (https://www.ebi.ac.uk/gwas/). After converting the genotype data to match the GWAS file's genotype build, we calculated the number of common variants between the GWAS files and the genotype data. Three phenotypes were removed from further analysis due to a limited number of variants.

## GWAS Data Processing

### Steps for Processing GWAS Data

1. **Download the GWAS file for a specific phenotype.**

2. **Transform the GWAS file to a specific format accepted by most PRS tools.** 
   - The sample transformation code for one phenotype (asthma) is shown below.

#### Original GWAS file for asthma:
| chromosome | base_pair_location | effect_allele | other_allele | effect_allele_frequency | beta   | standard_error | p_value | variant_id |
|------------|--------------------|---------------|--------------|-------------------------|--------|----------------|---------|------------|
| 1          | 100000012          | T             | G            | 0.28029314035229        | 0.0055 | 0.0077         | 0.4756  | NA         |
| 1          | 10000006           | A             | G            | 0.00547967251882486     | -0.1081| 0.0560         | 0.05366 | NA         |
| 1          | 100000135          | T             | A            | 0.00172932855852879     | 0.0255 | 0.1026         | 0.8037  | NA         |
| 1          | 100000374          | C             | G            | 0.000580831145737331    | 0.0204 | 0.2956         | 0.9451  | NA         |
| 1          | 100000827          | T             | C            | 0.329735132840309       | 0.0016 | 0.0074         | 0.8236  | NA         |

#### Transformation Code:
```python
import pandas as pd

df = pd.read_csv('gwas.csv.modified')
print(df.head())
print(df.columns)

column_map = {
    'chromosome': 'CHR', 
    'variant_id': 'SNP', 
    'base_pair_location': 'BP',
    'effect_allele': 'A1', 
    'other_allele': 'A2', 
    'p_value': 'P',
    'effect_allele_frequency': 'MAF', 
    'beta': 'BETA', 
    'standard_error': 'SE'
}

df = df.rename(columns=column_map)

df["SNP"] = "X"
df["INFO"] = 1
df["N"] = 449500

df = df[['CHR', 'BP', 'SNP', 'A1', 'A2', 'N', 'SE', 'P', 'BETA', 'INFO', 'MAF']]

df.to_csv('asthma.gz', compression='gzip', index=False, sep="\t")
```
| CHR | BP         | SNP | A1 | A2 | N     | SE    | P      | BETA   | INFO | MAF               |
|-----|------------|-----|----|----|-------|-------|--------|--------|------|-------------------|
| 1   | 100000012  | X   | T  | G  | 449500| 0.0077| 0.4756 | 0.0055 | 1    | 0.28029314035229  |
| 1   | 10000006   | X   | A  | G  | 449500| 0.056 | 0.05366| -0.1081| 1    | 0.0054796725188248|
| 1   | 100000135  | X   | T  | A  | 449500| 0.1026| 0.8037 | 0.0255 | 1    | 0.0017293285585287|

3.  Detect the genomic build of the GWAS.

4.  Handle missing RSID or SNP information:

- Generate two copies of the genotype data based on the GWAS build (hg19 or hg38).
- Find the common predictors and use the SNPs in the genotype .bim file to fill the SNP column in the GWAS.

| CHR | BP        | SNP         | A1 | A2 | N     | SE     | P      | BETA   | INFO | MAF               |
|-----|-----------|-------------|----|----|-------|--------|--------|--------|------|-------------------|
| 1   | 100008709 | rs55727773  | T  | G  | 449500| 0.0148 | 0.4387 | -0.0115| 1    | 0.881099991887257 |
| 1   | 100011667 | rs113462541 | C  | T  | 449500| 0.0096 | 0.4591 | 0.0071 | 1    | 0.835926208119123 |
| 1   | 100015703 | rs12028261  | T  | C  | 449500| 0.0071 | 0.6983 | 0.0027 | 1    | 0.505985657483862 |


5. Save the final file based on the following format:

```bash
# Define file paths for different data files
BED = asthma/asthma.bed
BIM = asthma/asthma.bim
FAM = asthma/asthma.fam
COV = asthma/asthma.COV
Height = asthma/asthma.PHENO
GWAS = asthma/asthma.gz
```


Below is the diagram showcasing data processing:

<img src="Data.png" width="50%">



## Processed Data

| Phenotype                                | PMID     | File Name                           | SNPs in GWAS (G) | SNPs in Genotype data (GE) | Common in G and GE | DOI                              | Cite            |
|------------------------------------------|----------|-------------------------------------|------------------|----------------------------|--------------------|----------------------------------|-----------------|
| asthma                                   | 34594039 | GCST90018795_buildGRCh37.tsv.gz     | 25837674         | 619653                     | 2867               | 10.1038/s41588-021-00931-x       | Sakaue2021      |
| blood_pressure_medication                | 34662886 | GCST90081464_buildGRCh38.tsv.gz     | 447993           | 619653                     | 56                 | 10.1038/s41586-021-04103-z       | Backman2021     |
| body_mass_index_bmi                      | 34594039 | GCST90018947_buildGRCh37.tsv.gz     | 20538803         | 619653                     | 2866               | 10.1038/s41588-021-00931-x       | Sakaue2021      |
| cholesterol_lowering_medication          | 34662886 | GCST90079486_buildGRCh38.tsv.gz     | 133986           | 619653                     | 56                 | 10.1038/s41586-021-04103-z       | Backman2021     |
| depression                               | 29662059 | UKBiobank_broad_12Jan18.txt         | 7641987          | 619653                     | 545218             | 10.1038/s41467-018-03819-3       | Howard2018      |
| gastro_oesophageal_reflux_gord_gastric_reflux | 34594039 | GCST90018848_buildGRCh37.tsv.gz     | 25843128         | 619653                     | 2906               | 10.1038/s41588-021-00931-x       | Sakaue2021      |
| hayfever_allergic_rhinitis               | 34662886 | GCST90077815_buildGRCh38.tsv.gz     | 458874           | 619653                     | 56                 | 10.1038/s41586-021-04103-z       | Backman2021     |
| high_cholesterol                         | 29892013 | GCST90029021_buildGRCh37.tsv        | 12007882         | 619653                     | 562302             | 10.1038/s41588-018-0144-6        | Loh2018         |
| hypertension                             | 33893285 | GCST90086092_buildGRCh37.tsv        | 15650645         | 619653                     | 619653             | 10.1038/s41467-021-21952-4       | GuindoMartinez2021 |
| hypothyroidism_myxoedema                 | 34594039 | GCST90018862_buildGRCh37.tsv.gz     | 25801084         | 619653                     | 2908               | 10.1038/s41588-021-00931-x       | Sakaue2021      |
| irritable_bowel_syndrome                 | 34741163 | GCST90016564_buildGRCh37.tsv        | 9885499          | 619653                     | 568535             | 10.1038/s41588-021-00950-8       | Eijsbouts2021   |
| migraine                                 | 34737426 | GCST90043745_buildGRCh37.tsv.gz     | 11831933         | 619653                     | 548955             | 10.1038/s41588-021-00954-4       | Jiang2021       |
| osteoarthritis                           | 36411363 | GCST90134279_buildGRCh37.csv        | 24880768         | 619653                     | 619653             | 10.1038/s41588-022-01221-w       | McDonald2022    |


## Perform quality controls on the data.
The initial quality controls varies from dataset ot dataset and the way dataset is being stored processed and handled.
after that execute the 
```
.
├── body_mass_index_bmi
│   ├── GCST90018947_buildGRCh37.tsv.gz
│   ├── gwas.csv.modified
│   └── Output.py
├── genotypes.bed
├── genotypes.bim19
├── genotypes.bim38
├── genotypes.fam
├── phenotype_file.txt
└── Step1-MakeDirectories.py
```



