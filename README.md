# Admixr+ User Guide

*By Pavel Salazar-Fernandez (epsalazarf@gmail.com)*

*Human Population and Evolutionary Genomics Lab | LANGEBIO*

## About
*Source: [admixr [GitHub]](https://github.com/martinsikora/admixr)*

This document explains how to use the R package `admixr` developed by Martin Sikora, including generating the input file and running the analyses, with brief explanations about each module.

## Index
* Generating the Input Files
    * PLINK Binary Files
    * Cleaning
    * FAM File
    * BIM File
    * RAW File

* ADMIXTURE
    * Running

* ADMIXR+
    * Requirements
    * Installing admixr
    * Pipeline Overview
    * Procedure
    * Module Description
        * [1. Preparations]
        * [2. PCA]
        * [3. Allele Frequencies]
        * [4. Fst]
        * [5. f3 Statistics]
        * [6. D Statistics]
        * [7. ADMIXTURE]

* References 

## Generating the Input Files

### PLINK Binary Fileset

`admixr` uses files from the binary PLINK set (BED,BIM,FAM) as input. If your data is in PED/MAP format, convert them with PLINK.

* Program: PLINK
* Input:
	- PED file
	- MAP file

* Command:
`plink --file [FILE] --make-bed --out [FILE]`

* Output:
	- FILE.BED
	- FILE.BIM
	- FILE.FAM

### Linkage Disequilibrium
Due to linkage disequilibrium, not all SNPs are informative. To thin the marker set for linkage disequilibrium:

* Program: PLINK

* Input:
	- FILE.bed
	- FILE.bim
	- FILE.fam

* Command:
`plink --bfile [FILE] --indep-pairwise 50 10 0.1`

* Output:
	- plink.prune.in (SNPs targeted for inclusion)
	- plink.prune.out (SNPs targeted for exclusion)

To apply the filter to the dataset and generate a new cleaned dataset, run now the following command:

`plink --noweb --bfile [FILE] --mind 0.1 --geno 0.05 --extract plink.prune.in --make-bed --out [FILE].clean`

The options used are:
* `--mind 0.1`: filter samples with >10% missing SNPs. Ancient DNA samples may be excluded due to this filter, so this option may be omitted to prevent that.
* `--geno 0.05`: filters SNPs with >5% missing samples.
* `--extract plink.prune.in`: filters SNPs with >0.1 LD in 50-SNP windows.
* `--out [FILE].clean`: The output file name. A name different from the original file is recommended to prevent overwriting (a tag can be added as in the example).


### FAM File

`admixr` uses population codes contained in the first column of the FAM file (usually reserved for FAMID). If your FAM does not contain the population tags in the first column, modify it (with a tool like `awk`) to look like this:

|||||||
|---|---|---|---|---|---|
|**POP1**|pop1.1|0|0|0|1|
|**POP1**|pop1.2|0|0|0|1|
|**POP1**|pop1.3|0|0|0|1|
|**POP2**|pop2.1|0|0|0|1|
|**POP2**|pop2.2|0|0|0|1|
|**POP2**|pop2.3|0|0|0|1|
|**POP3**|pop3.1|0|0|0|1|
|...|...|...|...|...|...|

### BIM File

The BIM file needs no further modification.

### RAW File
*Documentation: [PLINK: Recode](http://pngu.mgh.harvard.edu/~purcell/plink/dataman.shtml#recode)*

`admixr` uses an additive genotype matrix, where 0 represents a reference allele homozygous genotype, 1 an heterozygous genotype and 2 an alternative allele homozygous. The RAW output file is the closest to the format required. While `admixr` does not natively supports reading a RAW file, extra lines of code in the `admixr+.R` will convert it to the appropiate format.

To generate a RAW file:

* Program: PLINK

* Input: PLINK file set (normal or binary)

* Command:
`plink --file (or --bfile) [FILE].clean --recode A --out [FILE].clean`

* Output: [FILE].raw

The resulting RAW file should look like this:

|FID|ID|PAT|MAT|SEX|PHENO|SNP1_A|SNP2_C|SNP3_G|...|
|---|---|:--:|:--:|:--:|:--:|:--:|:--:|:--:|---|
|POP1|pop1.1|0|0|0|1|0|1|1|...|
|POP1|pop1.2|0|0|0|1|0|1|2|...|
|POP1|pop1.3|0|0|0|1|0|1|1|...|
|POP2|pop2.1|0|0|0|1|1|1|0|...|
|POP2|pop2.2|0|0|0|1|1|2|0|...|
|POP2|pop2.3|0|0|0|1|1|1|1|...|
|POP3|pop3.1|0|0|0|1|0|0|2|...|
|...|...|...|...|...|...|...|...|...|...|

## ADMIXTURE
*Documentation: [ADMIXTURE (UCLA)](https://www.genetics.ucla.edu/software/admixture/)*

### Running ADMIXTURE

* Program: ADMIXTURE

* Input:
	- [FILE].clean.bed
	- [FILE].clean.bim
	- [FILE].clean.fam

* Command:

`admixture --cv [FILE].clean.bed [Number of Ks] > [FILE].clean.[Number of Ks].log &`

* Output: ADMIXTURE will produce a set of files for each analysis:
	- .Q file: ancestry fractions for each individual.
	- .P file: allele frequencies in the inferred populations.
	- .log file: Details on parameters used and other values.

Each analysis for a given number of Ks must be run independently. Recommended values are all numbers between 3 and 10, but depends on the dataset. Keep in mind that higher Ks will require greater computing time.

The `--cv` flag estimates the cross-validation error for each run at the bottom of each log file. The lower the CV error, the more accurate the model involving that *K* number of clusters. To quickly view all CV errors, run:

`grep -h CV *.K*.log`

## ADMIXR+

`admixr+.R` is an R script modified from the [example.R](https://github.com/martinsikora/admixr/blob/master/doc/example.R) script given by the developer. Some modifications include: more code commenting and structure, an input section for all analyses in the code and automation tasks for a complete run or only selected modules. Feel free to modify the code to suit your needs.

### Requirements
* Software: R 3.2+
	- R Packages: **admixr**,  ggplot2, data.table, ape. Optional: foreach, doMC.

* Input files:
	- Genotypes: RAW file
	- SNP info: BIM file
	- Populations info: FAM file
	- Admixture: Q file. Optional: P file.

### Installing *admixr*

**Command:**
```R
install.packages("devtools")
library("devtools")
install_github("martinsikora/admixr")
```

### Pipeline Overview
1. Reads input data and transforms it.
2. Generates PCA plots.
3. Gets allele frequencies.
4. Calculates Fst for population pairs.
5. Performs outgroup f3 statistics calculations.
6. Obtains D statistics.
7. Makes Admixture plots and projections.

### Procedure
1. Set the working and output directory.
2. Select the analysis modules to be run.
3. Declare the names of the input files.
4. Adjust the parameters for each analysis. Skipped modules do not need any configuration.
5. Source the script. Messages will inform the progress.
6. Once finished, all the output files will be in the working directory.

### Modules Description
Text in quote block are citations from the Suplementary Material from Allentot *et al.* 2015 [2,3].

#### [1. Preparations]
Receives the input files, analyses parameters and loads the required enviroment for running the script, thus is mandatory to run it. If you need to re-run some analyses with the same input data, you may skip this module.

**Sub-modules**
* **INPUT**: This section includes all the settings for each analysis that must be modified by the user. This is the only section that requires modification, the rest of the code can be run with the default settings. For any special modifications, they must be done inside each module code.
* **PREPARATIONS**: Loads the required libraries, reads the input files and transforms the data as required by the package.

**Output**
* Console message notifying that the ADMIXR+ will begin the analysis. At the end of the run, another message will indicate its conclusion.

#### [2. PCA]
Performs a principal components analysis and plots the results.

**Sub-modules**
* **[2.1] FULL PCA**: Runs a conventional PCA.
* **[2.2] FAST PCA**: A faster approximate PCA.

**Output**
* `pca.plot.pdf`: contains paired PCA plots for the first 10 components.

#### [3. Allele Frequencies]
Calculates allele frequencies for each SNP in all populations. This module is required for Module [6].

**Output**
* No file created, but data generated is required by Module [4,6].

#### [4. Fst]
Calculates Fst for all possible pairs in data and generates plots to help visualize the results. Fst calculations are performed both with diploid genotypes (Weir & Cockerham, 1984) and with allele frequencies (Weir & Hill 2002). The sample-size corrected moment estimator of Weir and Hill is used, restricting the analysis to SNPs with minimum two alleles observed in each population of the pair.


**Output**
* `fst.comparison.pdf`: scatter plot comparing both Fst calculation methodologies. Ideally, all points must be aligned to the solid line representing concordance.
* `fst.nj.plot.pdf`: A dendrogram clustering populations according to their Fst scores. Negative Fst scores will generate branches oriented to the left.
* `fst.matrix.plot.pdf`: Heat map containing Fst values for each paired population. Population order is the same as given by the FAM file, starting from the bottom-left corner.

#### [5. f3 Statistics] (WIP)
Performs different kinds of f3 calculations and plot in designed populations against all other.

> **f3(Outgroup;Population1,Population2)** [Modules 5.1 and 5.3]
This “outgroup”-f3 statistic is expected to be proportional to the amount of shared genetic drift between Population1 and Population2 in their common ancestral population until their divergence. Unlike methods based on pairwise distances such as FST, genetic drift specific to Population1 or Population2 does not affect this statistic.

> **f3(PopulationTest;Population1,Population2)** [Module 5.2]
This is the original formulation of the f3 statistic as a statistical test for admixture. A significantly negative value of this statistic is evidence for a history of admixture in the test population, related to a pair of source populations Population1 and Population2. A positive value on the other hand does not exclude the possibility of admixture, as drift specific to the test population post-admixture enters adds a positive term to the statistic and can therefore obscure a real historical admixture signal.

**Sub-modules**
* **[5.1] Outgroup f3 Statistics**: Calculates f3 statistics of a designated population (`idxX`) against all other populations except the outgroup (`idxO`).
* **[5.2] Paired f3 Statistics**: Calculates and plots pairwise f3 for `popd1` and `popd2`.
* **[5.3] Grouped f3 Statistics**: Calculates f3 statistics by grouping populations as "outgroup" or "derived", which must be designated explicitly in the **INPUT** section.

**Output**
* **[5.1]:** `f3.outgroup.plot.pdf`: Ordered f3 values comparing the title population with all other. Higher values suggest closeness.
* **[5.2]:** `f3.outgroup.pairs.plot.pdf`: Compares f3 values of two populations. Alignment of points to the line indicate same degree of closeness for both populations compred to the third.
* **[5.3]:** `f3.outgroup.grouped.plot.pdf`: A plot with the populations divided in their two designated sections: outgroup and derived.

#### [6. D Statistics]

>**D(Outgroup,PopulationTest)(Population1,Population2)**
This D-statistic measures whether the data is consistent with a four-population tree in which Population1 and Population2 form a clade with each other, to the exclusion of the test population and the outgroup. The expected value in case of consistency with the proposed tree is zero. Significant deviations from zero reject the proposed tree, with negative values indicating that the test population is closer to Population1, and positive values indicating that the test population is closer to Population2.

**Output**
* `D.plot.pdf`: Plot with populations ordered as input, and a line denoting the expected value of 0.

#### [7. ADMIXTURE]
Module used to produce ADMIXTURE plots, one at a time, with automatic colors.

**Sub-modules**
* **[7.1] Admixture Plot**: Generates a simple admixture with populations grouped together.
* **[7.2] Admixture Projection (WIP)**: Used to plot projections made with ADMIXTURE.

**Output**
* `admixture.full.plot.pdf`
* `admixture.projected.plot.pdf`

## References
1. [martinsikora: admixr [GitHub]](https://github.com/martinsikora/admixr)
2. [Population genomics of Bronze Age Eurasia [Nature]](http://www.nature.com/nature/journal/v522/n7555/full/nature14507.html)
3. [Population genomics of Bronze Age Eurasia (Supplementary)[Nature]](http://www.nature.com/nature/journal/v522/n7555/full/nature14507.html#supplementary-information)
4. [PLINK Documentation](http://pngu.mgh.harvard.edu/~purcell/plink/index.shtml)