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


Below is the diagram showcasing data processing:

<img src="Data.png" width="50%">


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





