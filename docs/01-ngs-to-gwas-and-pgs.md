# From NGS to GWAS and Polygenic Scores: A Beginner-Friendly Practical Guide

**A step-by-step introduction to joint-called VCFs, phenotype files, PLINK2, GWAS beta values, covariates, PCA, linkage disequilibrium (LD), and polygenic score (PGS) calculation.**

> **Who this is for:** Students and researchers who understand basic sequencing/variant calling but are new to statistical genetics. No previous GWAS experience is assumed.
>
> **What this is not:** A clinical genetic test, a validated athletic-performance predictor, or a complete production pipeline. **All included participants, variants, effects, and measurements are synthetic.** Do not interpret scores from this repository as health risks.

## Contents

1. [The big picture](#1-the-big-picture)
2. [The five key ideas](#2-the-five-key-ideas)
3. [From FASTQ to a joint-called VCF](#3-from-fastq-to-a-joint-called-vcf)
4. [How a VCF connects to a phenotype form](#4-how-a-vcf-connects-to-a-phenotype-form)
5. [Where a GWAS beta comes from](#5-where-a-gwas-beta-comes-from)
6. [Why we adjust for age, sex, and other covariates](#6-why-we-adjust-for-age-sex-and-other-covariates)
7. [What PC1, PC2 and genetic ancestry mean](#7-what-pc1-pc2-and-genetic-ancestry-mean)
8. [What LD is and why it matters](#8-what-ld-is-and-why-it-matters)
9. [From GWAS results to a PGS](#9-from-gwas-results-to-a-pgs)
10. [Run the synthetic example](#10-run-the-synthetic-example)
11. [Practical PLINK2 workflow](#11-practical-plink2-workflow)
12. [How to evaluate a score without fooling yourself](#12-how-to-evaluate-a-score-without-fooling-yourself)
13. [What changes for a real large cohort](#13-what-changes-for-a-real-large-cohort)
14. [Common mistakes and troubleshooting](#14-common-mistakes-and-troubleshooting)
15. [Glossary](#15-glossary)
16. [Learning resources and references](#16-learning-resources-and-references)

---

## 1. The big picture

Suppose we sequence many volunteers and measure a **quantitative phenotype**, such as *maximal oxygen uptake* (VO₂max, in mL/kg/min). We want to ask:

> **Are differences in genetic variation associated with differences in measured VO₂max across people?**

There are **two different scientific tasks**:

- **GWAS (discovery):** Estimate the association between **each** variant and the trait across many people.
- **PGS (scoring):** Combine previously estimated variant weights with an individual's genotypes to obtain a single summary score.

A third task, **validation**, asks whether the score adds useful information about the trait in **independent** people.

```text
                GENOTYPING SIDE                         MEASUREMENT SIDE

Sequencer: R1/R2 FASTQs                           Assessment forms / lab tests
           |                                                 |
Alignment: BAM / CRAM                                      Athlete ID
           |                                           Age, sex, VO2max, ...
Variant calling: per-sample gVCF                            |
           |                                                 |
Joint genotyping: MULTI-SAMPLE VCF                         TSV / CSV
           |                                                 |
Variant/sample QC                                           |
           |                                                 |
PLINK2: .pgen + .pvar + .psam                               |
           |                                                 |
           +----------------- MATCH SAMPLE IDs --------------+
                                     |
                             Genome-wide analysis
                                     |
             Test SNP 1 across people -> beta 1, SE, P
             Test SNP 2 across people -> beta 2, SE, P
             Test SNP 3 across people -> beta 3, SE, P
                           ... millions of times
                                     |
                             GWAS summary statistics
                                     |
                   Construct LD-aware / validated weights
                                     |
                      Apply weights to target genotypes
                                     |
                          PGS per person + validation
```

**The crucial insight:** One person's VCF and one VO₂max measurement are **not enough to estimate a SNP's effect**. In a GWAS, the same SNP is compared *across many people with measurements*. A published external GWAS can instead supply weights for scoring a single person's genotypes, although a single-person score cannot validate a predictive model.

### A small table of inputs and outputs

| Stage | Input | What is learned | Output |
|---|---|---|---|
| Sequencing | DNA | Sequenced reads | FASTQ |
| Alignment | Reads + reference | Where reads align | BAM/CRAM |
| Variant calling | Alignments + reference | Candidate genotypes | gVCF / VCF |
| Joint genotyping | Many gVCFs | Consistent cohort genotypes | Multi-sample VCF |
| PLINK conversion | Quality-controlled VCF | Efficient genotype representation | PGEN/PVAR/PSAM |
| Phenotype matching | Sample IDs + measurements | Who has which trait value | Phenotype/covariate tables |
| GWAS | Genotypes + trait + covariates | Estimated association **per variant** | BETA, SE, P, effect allele |
| PGS development | GWAS weights + LD information | A set of scoring weights | Scoring file |
| PGS calculation | Scoring file + target genotypes | Score **per person** | `.sscore` |
| Evaluation | Score + independently measured outcome | Predictive value and limitations | R², calibration, etc. |

---

## 2. The five key ideas

### Idea A — A genotype is an observed genetic state

At a SNP with reference allele **A** and alternate allele **G**:

| VCF `GT` | Alleles | G-allele copies |
|---|---|---:|
| `0/0` | A/A | 0 |
| `0/1` | A/G | 1 |
| `1/1` | G/G | 2 |

For a biallelic diploid site, this **0, 1, 2** is often called *dosage*. In real files, unphased genotypes are often separated by `/`, phased by `|`; missing calls, non-diploid chromosomes and multiallelic sites require separate handling.

### Idea B — A phenotype is a measurement

Examples: VO₂max, height, blood pressure, LDL cholesterol, grip strength, diagnosis status. **Phenotypes do not come from the VCF.** They come from measurements, questionnaires, clinical records, or other validated sources.

### Idea C — GWAS estimates beta across people

A GWAS tests each variant in a dataset and estimates the average change in the measured phenotype associated with an additional copy of its tested allele. **Beta is a statistical estimate, not an attribute of the DNA sequence.**

### Idea D — PGS uses many weights to score one person

A polygenic score is a weighted sum:

\[
PGS_i = \sum_{j=1}^{M} w_j G_{ij}
\]

`G` = effect-allele count for the person, `w` = the published or trained weight for the variant. **PGS weights may differ from raw GWAS betas**, because score-building methods adjust them for LD and statistical noise.

### Idea E — A significant association is not necessarily predictive or causal

A SNP can be correlated with VO₂max without directly causing a change in VO₂max. Even a statistically associated SNP, or a whole PGS, may explain very little of the trait's variation.

---

## 3. From FASTQ to a joint-called VCF

This section supplies **an overview of real NGS preprocessing**, not a copy/paste-ready production workflow. Production commands require validated references, read groups, duplicate handling, BQSR or pipeline-equivalent calibration, interval design, caller-specific filtering, and QC.

### 3.1 Sequencing and mapping

After WGS, paired-end reads usually arrive as:

```text
ATH001_R1.fastq.gz
ATH001_R2.fastq.gz
ATH002_R1.fastq.gz
ATH002_R2.fastq.gz
...
```

A common high-level workflow is:

```text
FASTQ -> FastQC/MultiQC -> alignment -> sorted BAM/CRAM -> QC -> variant calling
```

Example **illustrative alignment command** (not a complete production alignment pipeline):

```bash
bwa-mem2 mem -t 16 GRCh38.fa ATH001_R1.fastq.gz ATH001_R2.fastq.gz \
  | samtools sort -@ 8 -o ATH001.sorted.bam
samtools index ATH001.sorted.bam
```

| Part | Meaning |
|---|---|
| `bwa-mem2 mem` | Align sequencing reads to a reference |
| `-t 16` | Use 16 alignment threads |
| `GRCh38.fa` | The reference genome (index files must already exist) |
| `R1` / `R2` | Paired FASTQ reads |
| `samtools sort` | Sort alignments by genomic coordinate |
| `samtools index` | Create an index for querying regions |

**Do not use this bare alignment command as a final clinical/research pipeline.** In particular, add validated `@RG` read groups and all required processing before HaplotypeCaller.

### 3.2 Single-sample variant calling in gVCF mode

After alignment QC, a typical GATK call is:

```bash
gatk HaplotypeCaller \
  -R GRCh38.fa \
  -I ATH001.processed.bam \
  -O ATH001.g.vcf.gz \
  -ERC GVCF
```

`-R` selects the reference; `-I` is the processed alignment; `-O` is the output; `-ERC GVCF` requests the per-sample genomic VCF representation used for later cohort genotyping.

**Why gVCF?** A variant-only VCF often lists only variant positions. A gVCF additionally represents evidence about reference positions/regions. That helps a joint caller tell a confident homozygous-reference genotype from a site where there is insufficient evidence.

### 3.3 Joint genotyping

```text
ATH001.g.vcf.gz ----\
ATH002.g.vcf.gz -----+--> GenomicsDBImport --> GenotypeGVCFs --> cohort.raw.vcf.gz
ATH003.g.vcf.gz ----/
```

For a few samples and a specific interval, GATK syntax resembles:

```bash
gatk GenomicsDBImport \
  -V ATH001.g.vcf.gz \
  -V ATH002.g.vcf.gz \
  -V ATH003.g.vcf.gz \
  -L chr1:1-1000000 \
  --genomicsdb-workspace-path chr1_db

gatk GenotypeGVCFs \
  -R GRCh38.fa \
  -V gendb://chr1_db \
  -O cohort.chr1.raw.vcf.gz
```

For thousands of genomes, shard by intervals, use a validated scalable genotyping strategy and suitable compute/storage resources. If using DRAGEN, follow the documentation for the **specific installed DRAGEN version**; the concepts remain the same but commands and features differ. GATK calls still need appropriate downstream variant filtering before analysis.

**Important:** Joint genotyping makes the genotype representation consistent across people; it does **not** estimate genetic effects on VO₂max.

---

## 4. How a VCF connects to a phenotype form

### 4.1 Understand rows versus columns

Here is a simplified multi-sample VCF excerpt:

```text
#CHROM POS  ID    REF ALT ... FORMAT ATH001 ATH002 ATH003 ATH004
1      1000 rs101 A   G   ... GT     0/0    0/1    1/1    0/1
1      2000 rs102 C   T   ... GT     0/1    0/0    1/1    0/1
2      3000 rs103 G   A   ... GT     1/1    0/1    0/0    0/1
```

- **Each row is a genomic position (variant).**
- **Each sample column is one person.**
- The intersecting cell is that person's genotype at that position.
- This is a simplified display; a legal VCF also has `QUAL`, `FILTER`, `INFO` and other required fields. The downloadable example in this repository is a proper small VCF.

### 4.2 Your participant/assessment sheet is separate

```text
FID  IID     VO2max
0    ATH001  40
0    ATH002  45
0    ATH003  50
0    ATH004  47
```

Here `IID` is the **individual ID** and `FID` is the family ID; `0` denotes an unspecified family in this example. The actual record might also include age, sex, sport, measurement method, date and consent status.

Now view **one SNP across all people**:

| ID | Genotype at rs101 | Copies of G | VO₂max |
|---|---|---:|---:|
| ATH001 | A/A | 0 | 40 |
| ATH002 | A/G | 1 | 45 |
| ATH003 | G/G | 2 | 50 |
| ATH004 | A/G | 1 | 47 |

The common `IID` joins the VCF genotype column to the phenotype table row. We do **not** create one VO₂max value for every SNP; each person has one recorded VO₂max per assessment timepoint, and that same value is considered in the test of each SNP.

> **Timepoint warning:** If VO₂max is measured repeatedly, do not arbitrarily duplicate the same person into ordinary GWAS as if they were independent participants. Define an appropriate baseline/outcome or use a longitudinal model.

### 4.3 Why the match must be audited

If `ATH001` in the genome file is really `ATH002` in the assessment sheet, results can be meaningless. In real studies, check sample identity, duplicate records, ID merges, sex concordance, relatedness and missing measurements **before** GWAS.

A quick sample-ID check, if `bcftools` is installed:

```bash
bcftools query -l cohort.filtered.vcf.gz | sort > genotype_iids.txt
# For an explicitly FID/IID/VO2max, headered TSV:
tail -n +2 phenotypes.tsv | cut -f2 | sort > phenotype_iids.txt
comm -3 genotype_iids.txt phenotype_iids.txt
```

No output from `comm -3` means the two sorted lists match exactly. In real projects this is a **first check**, not a substitute for confirming correct identity or handling valid missing phenotypes.

---

## 5. Where a GWAS beta comes from

### 5.1 GWAS is one SNP compared across people, repeated millions of times

Imagine a toy dataset at **rs101**:

| Person | G copies | VO₂max |
|---|---:|---:|
| P1 | 0 | 40 |
| P2 | 0 | 42 |
| P3 | 1 | 45 |
| P4 | 1 | 47 |
| P5 | 2 | 50 |

A simple linear regression fits:

\[
\text{VO₂max}_i = \alpha + \beta_{\mathrm{rs101}}G_i + \varepsilon_i
\]

- `Gᵢ`: effect/tested allele count (0/1/2) in person `i`.
- `α`: expected value when `G = 0` under the model.
- `β`: fitted *slope*: average trait difference per additional G allele.
- `ε`: residual: the part of the measurement the simple line does not explain.

In this **five-person fictional example**, simple regression yields **β ≈ +4.69 mL/kg/min per G allele**. This is just a classroom calculation; five people cannot support a trustworthy genome-wide finding.

Where does +4.69 come from? It is the best-fit slope of the genotype/phenotype points—not something stored in the VCF, not a functional annotation, and not manually assigned by the analyst.

For a single SNP without covariates, it can be calculated as:

\[
\widehat\beta =
\frac{\sum_i (G_i-\bar G)(Y_i-\bar Y)}
     {\sum_i (G_i-\bar G)^2}
\]

`Y` is the measured trait. The numerator captures whether G and Y vary together; the denominator captures how much G varies across people.

### 5.2 Then the computer moves to the next SNP

```text
GWAS test 1: rs101 genotype for all participants + all their VO2max values
GWAS test 2: rs102 genotype for all participants + all their VO2max values
GWAS test 3: rs103 genotype for all participants + all their VO2max values
... and so on
```

Each SNP can have a different beta, standard error and P-value. No single person's VO₂max generates millions of betas. **An adequately sized set of individuals with both genotypes and phenotypes is needed.**

### 5.3 What the output means

| Field | Plain-English interpretation |
|---|---|
| `ID` | Variant identifier |
| `A1` | Allele whose effect is being estimated in that association result |
| `BETA` | Estimated phenotype change per additional copy of `A1` |
| `SE` | Standard error: estimated uncertainty of the beta |
| `P` | Statistical evidence against the null hypothesis of zero association |
| `OBS_CT` | Number of observations included in that association test |

**A1 is not guaranteed to be ALT.** PLINK2 may select the minor allele as `A1`. Always preserve the effect-allele label alongside beta when moving from GWAS to PGS.

For a continuous phenotype, beta is generally in the phenotype's entered units (unless the phenotype was transformed/standardized). For a binary disease GWAS with logistic regression, the beta is typically a **log-odds coefficient**, not a change in VO₂max.

A tiny P-value does **not** prove causality or good prediction. The conventional threshold `5 × 10⁻⁸` accounts approximately for genome-wide multiple testing of common variants; it is not a universal criterion for all designs.

---

## 6. Why we adjust for age, sex, and other covariates

A **covariate** is a measured characteristic that we include in the statistical model so we compare people more appropriately.

Imagine carriers of G at rs101 appear to have higher VO₂max, but those carriers are also younger and more likely to be endurance athletes. A naive model may incorrectly attribute part of these differences to G.

Without covariates:

\[
Y = \alpha + \beta_G G + \varepsilon
\]

With selected covariates:

\[
Y = \alpha + \beta_G G + \beta_A Age + \beta_S Sex +
    \beta_{PC1} PC1 + \beta_{PC2} PC2 + \varepsilon
\]

Now the reported `β_G` is the SNP association **conditional on the variables in the model**. PLINK estimates the covariate coefficients too, but the SNP beta is usually the primary GWAS result.

| Covariate | Why include it? | Where does it come from? |
|---|---|---|
| Age | Traits can differ by age | Participant record |
| Sex / sex-related variable | Traits can differ across physiological strata | Study metadata, interpreted appropriately |
| Recruitment site | Site-specific assessment or recruitment patterns | Study records |
| Genotyping / sequencing batch | Technical artifacts can differ by batch | Laboratory records |
| PC1–PCn | Broad genetic population structure may confound associations | Calculated from genotype data |
| Sport / training | Sometimes essential for the study question; sometimes overadjustment | Training logs / protocol |

**More adjustment is not automatically better.** If a variable lies on the causal pathway from genotype to phenotype (e.g. a biological mediator), adjusting for it may remove a real component of association. Decide on your scientific question before choosing covariates; do not blindly add every recorded field.

### Confounding versus causation

```text
Genetic ancestry patterns ------> genotype frequencies
             |
             +--> recruitment / environment / sport mixture --> VO2max

If ancestry-related structure is ignored, a genotype may look associated
with VO2max even when it is mostly tagging the cohort's composition.
```

Statistical adjustment reduces some confounding but does not by itself prove that a genetic variant causes a difference in fitness.

---

## 7. What PC1, PC2 and genetic ancestry mean

**PCA = principal component analysis.** It compresses patterns of variation across many genetic markers into a few numerical axes.

Suppose your cohort includes individuals whose allele frequencies differ somewhat because of population history. PCA can reveal these broad differences even when they are spread across thousands of SNPs.

- **PC1** captures the largest axis of variation *in the marker dataset used for PCA*.
- **PC2** captures another independent axis of variation.
- **PC3, PC4, ...** capture further axes.

They are **not genes, variant effects, fitness scores or clinical measurements**. They are coordinates summarizing genomic variation. PCA axes can also reflect relatedness, sequencing batch artifacts or unusual samples: interpret plots before treating PCs as ancestry.

### A typical PCA preparation procedure

1. Keep high-quality common autosomal SNPs and suitable samples.
2. Remove or account for closely related people and outliers as required.
3. **LD-prune** markers so one highly correlated region does not dominate PCA.
4. Run PCA, inspect a PC1-vs-PC2 scatterplot, and use an appropriate number of PCs as covariates.

For **real dense cohort data** (not the 12-SNP toy example):

```bash
# Toy thresholds are common teaching examples, not universal QC requirements.
plink2 --pfile cohort.qc --autosome --maf 0.05 \
  --indep-pairwise 200 50 0.2 --out results/pca_prune

plink2 --pfile cohort.qc --extract results/pca_prune.prune.in \
  --pca 10 approx --out results/cohort_pca
```

**Command meanings**

- `--autosome`: analyze autosomal chromosomes.
- `--maf 0.05`: keep variants whose minor-allele frequency is at least 5% in this step.
- `--indep-pairwise 200 50 0.2`: examine windows of 200 variants, advance 50 variants, and prune strongly correlated pairs above the chosen LD threshold.
- `.prune.in`: selected set of less-correlated variant IDs.
- `--extract`: keep exactly those variants for PCA.
- `--pca 10 approx`: calculate 10 approximate principal components; suitable for large cohorts after QC.
- `.eigenvec`: output with sample IDs and PC values; `.eigenval` contains corresponding eigenvalues.

**Caution:** These commands are schematic. In real cohorts, exclude long-range LD regions when appropriate, verify sample/variant QC, handle ancestry heterogeneity and relatedness, and evaluate whether the PCA captures batch effects. A small example of 12 SNPs cannot generate scientifically meaningful ancestry PCs.

---

## 8. What LD is and why it matters

**LD = linkage disequilibrium**, the non-random association of alleles at different genetic positions. Nearby variants often travel together across generations, although physical distance is not the only determinant of LD.

Imagine three nearby SNPs:

```text
chr1:       1000      1050      1100
            SNP A     SNP B     SNP C
Person 1      0         0         0
Person 2      1         1         1
Person 3      2         2         2
Person 4      0         0         0
```

If SNP A and SNP B always have the same genotype pattern, their genotype correlation has **r² = 1**. With no linear correlation, **r² = 0**. `r²` is often used as a practical measure of LD for biallelic SNPs.

### Why is LD relevant to GWAS?

If one causal variant is correlated with three nearby markers, all three may show an association, even though only one biological mechanism underlies the signal. **Associated SNP ≠ independently causal SNP.** Fine-mapping and functional data are needed for stronger causal claims.

### Why is LD relevant to PGS?

Suppose SNP A, B and C tag mostly the same underlying signal. Naively adding all three unadjusted marginal betas may overcount overlapping information. PGS construction methods therefore account for LD.

| Procedure | What it does | Typical use |
|---|---|---|
| **LD pruning** | Selects a more independent marker set without using GWAS P-values | PCA, other genotype analyses |
| **LD clumping** | Uses association P-values and an LD threshold to retain representative SNPs | Clumping + thresholding PGS |
| **LD-aware shrinkage/modeling** | Adjusts many variant weights while modeling correlation | LDpred2, PRS-CS, etc. |

**Important distinction:** We LD-prune before PCA; for PGS construction we might LD-clump, or instead use a method that models LD. **Do not LD-prune the full dataset indiscriminately and assume the result is the best possible PGS.**

### A practical example

The included synthetic VCF intentionally makes `toy_snp01` and `toy_snp02` strongly correlated. Run `scripts/inspect_toy_data.py` to see their observed genotype-correlation `r²`. This illustrates LD; it does not imply one variant truly causes a trait.

---

## 9. From GWAS results to a PGS

### 9.1 Where do the weights come from?

**Route A: Apply an external published score**

1. Find a relevant validated score in the [PGS Catalog](https://www.pgscatalog.org/).
2. Check the phenotype, training ancestry, validation populations, and target genome build.
3. Download the **scoring file** (variant ID/location, effect allele, effect weight).
4. Harmonize it against your target genotypes, including strand/build/ref-alt handling.
5. Calculate scores and assess variant match rate and independent predictive performance.

No new GWAS is needed to calculate the published score. However, your target cohort may not have the disease events or phenotype data needed to **validate** it. Also confirm that the score's discovery data do not overlap the evaluation cohort.

**Route B: Develop a new score**

1. Define the phenotype and planned covariates.
2. Partition discovery, tuning and test participants **before** model training. Keep relatives together; independent external testing is best.
3. Run GWAS in discovery data to estimate associations.
4. Build PGS weights using methods suitable for LD, population structure and sample size.
5. Choose tuning parameters **without using test-set outcomes**.
6. Apply weights to independent target genotypes and report predictive performance.

A 10,000-person cohort is a valuable sequencing resource, but complex traits often have many small effects; **10,000 does not guarantee enough power to build a useful new PGS**. If people are selected by elite-athlete status, sample selection can also distort transferability to the wider population.

### 9.2 Manual PGS calculation

These effects are **fictional** and pre-specified:

| SNP | Effect allele | Weight |
|---|---|---:|
| SNP1 | G | +2.0 |
| SNP2 | T | +1.5 |
| SNP3 | C | −1.0 |

For an individual with one G at SNP1, two T at SNP2 and one C at SNP3:

\[
PGS = (1\times 2.0) + (2\times 1.5) + (1\times -1.0) = 4.0
\]

This score is **not VO₂max = 4**, **not 4% athletic ability**, and **not a probability of a disease**. It is a weighted sum requiring statistical validation and calibration before meaningful interpretation.

### 9.3 Allele harmonization is essential

If the model weight is for allele **G** but the genotype tool accidentally counts **A**, the contribution points in the wrong direction. Careful handling is needed for:

- Genome build GRCh37 vs GRCh38; coordinate liftover when required.
- REF/ALT and effect allele; `A1` is not always ALT.
- Strand flips and palindromic A/T or C/G variants.
- Multiallelic variants, indels and duplicate variant IDs.
- Missing genotypes, genotype accuracy, dosage/imputation and allele frequencies.
- Score coverage: how many **expected scoring variants** were actually matched/used.

The [PGS Catalog Calculator](https://pgsc-calc.readthedocs.io/en/latest/) automates many steps for supported inputs, but its documentation currently notes **limitations for direct WGS inputs**. WGS-derived genotypes need to be prepared appropriately; do not assume an unmodified WGS VCF will work or that a missing scoring site equals homozygous reference.

---

## 10. Run the synthetic example

This repository includes a **small, reproducible educational example** with 80 fictional participants and 12 fictional SNPs. There are no real patient or athlete data.

### 10.1 Repository layout

```text
genomics-analysis-handbook/
├── docs/01-ngs-to-gwas-and-pgs.md
├── requirements.txt
├── .gitignore
├── data/
│   ├── toy_joint_called.vcf      # VCF: SNP rows, sample columns
│   ├── phenotypes.tsv            # IID -> measured VO2max
│   ├── covariates.tsv            # IID -> age, sex, group, sport
│   └── illustrative_weights.tsv  # pre-specified fictional PGS weights
├── scripts/
│   ├── generate_toy_data.py      # recreates all four example input files
│   └── inspect_toy_data.py       # shows beta, adjusted beta, LD and scores
└── results/                      # generated locally; ignored by Git
```

### 10.2 Install minimal Python packages and run

From the folder containing this README:

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt

python scripts/generate_toy_data.py
python scripts/inspect_toy_data.py
```

Expected explanation of the output (values from the deterministic fictional dataset):

```text
SNP = toy_snp01; number of synthetic athletes = 80
Unadjusted beta = 3.284
Adjusted beta   = 2.717
Genotype correlation r², toy_snp01 vs toy_snp02 = 0.787
ATH001 toy PGS = 3.1
```

The first beta uses only `G copies` to explain VO₂max. The second adjusts for `Age`, `SexCode`, `ToyGroup` and `EnduranceSport`. `ToyGroup` is an **invented categorical variable**, not actual genetic ancestry or PC1. The script calculates coefficients by ordinary least squares to show the mathematics—**not** a publication-grade GWAS.

After running, inspect:

```bash
head -n 5 data/phenotypes.tsv
head -n 5 data/covariates.tsv
head -n 5 results/toy_gwas_betas.tsv
head -n 5 results/toy_pgs_scores.tsv
```

Try opening `data/toy_joint_called.vcf` in a text editor and locate `toy_snp01` and `ATH001`. Then find `ATH001` in the phenotype file. **The identifier is the join key**; the VCF does not contain VO₂max.

### 10.3 What the Python example deliberately does *not* do

- It does not model actual human ancestry or a real genomic build.
- It does not estimate reliable genome-wide P-values or discover true loci.
- It does not perform PCA on 12 SNPs: that would be uninformative.
- It does not construct LD-adjusted PGS weights from a discovery GWAS.
- Its PGS weights are **artificial, pre-specified** for practice; they are not external research estimates.

---

## 11. Practical PLINK2 workflow

> **Prerequisite:** Install [PLINK2](https://www.cog-genomics.org/plink/2.0/) and ensure `plink2` is on your `PATH`. Commands in this section are documented examples; the supplied Python exercises were run and checked, but PLINK2 was **not installed in the environment used to package this guide**. Apply full QC before treating any results as research output.

### 11.1 Import the included synthetic VCF

```bash
mkdir -p results
plink2 --vcf data/toy_joint_called.vcf --make-pgen --out results/toy
```

Meaning:

- `--vcf`: read VCF genotypes.
- `--make-pgen`: save compact PLINK2 binary genotype files.
- `--out results/toy`: prefix the output files with `results/toy`.

Expected files: `toy.pgen` (genotypes), `toy.pvar` (variant metadata) and `toy.psam` (participant IDs). Import does **not** compute beta.

### 11.2 Run a one-variant-at-a-time association analysis

```bash
# Unadjusted teaching example (no additional covariates)
plink2 --pfile results/toy \
  --pheno data/phenotypes.tsv --pheno-name VO2max \
  --glm allow-no-covars hide-covar \
  --out results/toy_unadjusted

# Adjusted example: variables recorded in the synthetic assessment table
plink2 --pfile results/toy \
  --pheno data/phenotypes.tsv --pheno-name VO2max \
  --covar data/covariates.tsv \
  --covar-name Age,SexCode,ToyGroup,EnduranceSport \
  --glm hide-covar \
  --out results/toy_adjusted
```

| Argument | Meaning |
|---|---|
| `--pfile results/toy` | Read `toy.pgen/.pvar/.psam` |
| `--pheno` | Load measured values matched by participant ID |
| `--pheno-name VO2max` | Analyze the column named `VO2max` |
| `--covar` | Load additional columns for adjustment |
| `--covar-name ...` | Use only the listed covariates |
| `--glm` | Fit one regression for each eligible variant |
| `allow-no-covars` | Explicitly permit the unadjusted exercise |
| `hide-covar` | Keep covariate coefficient rows out of the association report |
| `--out` | Write results with the chosen prefix |

Expected main association file names end in `.VO2max.glm.linear`. The exact default output columns are version-dependent. Check the result's **A1, BETA, SE and P** columns before interpreting it. The number of eligible variants and samples may change if you add filters.

The toy `SexCode` is numeric **0/1** in the covariate file; do not confuse it with PLINK's `SEX` field convention (`1` male, `2` female) when importing real metadata.

### 11.3 Calculate a PGS from the fictional weights

The included `illustrative_weights.tsv` contains:

```text
ID             A1  BETA
toy_snp01      G   2.3
toy_snp03      A  -0.4
...
```

Apply these **pre-specified toy weights**:

```bash
plink2 --pfile results/toy \
  --score data/illustrative_weights.tsv 1 2 3 header-read cols=+scoresums \
  --out results/toy_pgs
```

- `1 2 3` = the columns containing variant ID, effect allele and weight, respectively.
- `header-read` = first line is a header containing score names.
- `cols=+scoresums` = include the weighted **sum**, not merely the default allele-normalized score average.
- `toy_pgs.sscore` = per-person scoring output; inspect the `*_SUM` column.

**Audit the number of variants used**. By default, PLINK2 can mean-impute missing dosages when scoring; this affects interpretation. For real projects specify the missingness policy, matching rate, and external score harmonization method rather than simply trusting a numeric `.sscore` output.

### 11.4 What a production analysis might look like

Once you have a **fully QC'd** large dataset `cohort.qc.pgen/.pvar/.psam`, and quality-controlled `phenotypes.tsv` and `covariates_with_pcs.tsv` files, a *schematic* association command is:

```bash
plink2 --pfile cohort.qc \
  --pheno phenotypes.tsv --pheno-name VO2max \
  --covar covariates_with_pcs.tsv \
  --covar-name Age,SexCode,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10 \
  --autosome --mac 20 \
  --glm hide-covar \
  --out results/VO2max_GWAS
```

**What changed compared with the toy?** The `PC` columns now come from *genomic PCA*, and `--mac 20` requires at least 20 observed copies of the less frequent allele for a conventional common-variant analysis. The number of PCs and filters must be justified for the cohort. For related people or complex sample structure, mixed-model association software (e.g. REGENIE or SAIGE, according to trait/study design) may be preferable to ordinary PLINK2 `--glm`.

### 11.5 Typical real-genome preprocessing and QC checks

```bash
# Overview of variants / samples (for a quality-controlled input)
plink2 --pfile cohort.qc --freq --out results/cohort_frequency
plink2 --pfile cohort.qc --missing --out results/cohort_missingness
plink2 --pfile cohort.qc --hardy --out results/cohort_hwe
```

These report allele frequencies, missing calls and Hardy–Weinberg statistics. They are **diagnostics**, not a complete QC pipeline. A research analysis also needs sample contamination/identity checks, relatedness, batch effects, ancestry evaluation, sequencing/genotype quality, sex-chromosome handling, variant normalization and allele harmonization. HWE filtering should not be applied indiscriminately across ancestry groups or to clinically relevant rare variants.

---

## 12. How to evaluate a score without fooling yourself

### GWAS significance is not PGS performance

Suppose a SNP has `P = 10⁻¹⁰`. That suggests an association under the GWAS model, but does not tell you that the SNP—alone or combined with others—predicts VO₂max accurately.

For a quantitative score, compare:

```text
Baseline model:  VO2max ~ age + sex + ancestry PCs + justified covariates
PGS model:       VO2max ~ age + sex + ancestry PCs + justified covariates + PGS
```

In **independent test participants**, report the increase in out-of-sample **R²**, confidence intervals and model performance across relevant groups. Correlations between PGS and VO₂max alone can overstate usefulness if confounding variables are ignored.

For a disease PRS, evaluate discrimination (e.g., AUROC), calibration and meaningful absolute risk estimates **only** when outcome and follow-up data support those claims. Raw PRS is not a disease probability.

### Three major threats

1. **Data leakage:** If you estimate SNP effects using *all* participants, then report predictive accuracy in those same participants, results are optimistic. Separate discovery, tuning and evaluation by person and by family; external replication is ideal.
2. **Ancestry portability:** Score performance may change across populations, owing to LD, allele frequency and environmental differences. Comparing raw PGS means across ancestry groups is not a valid claim about innate ability.
3. **Selection bias:** A cohort selected for elite sports may not represent the general population. Predicting selected athlete status or extrapolating a fitness score to all people requires additional design and controls.

### A publication-oriented reporting checklist

- [ ] Defined primary phenotype and protocol, units and assessment dates.
- [ ] Documented genome reference build, alignment and variant-calling versions.
- [ ] Reported number of individuals and variants before/after QC.
- [ ] Audited sample ID matching, duplicate/related samples and missing phenotypes.
- [ ] Documented ancestry inference, population structure, chosen PCs and covariates.
- [ ] Reported GWAS model, effect allele, effect estimate, SE, P and multiple-testing methods.
- [ ] Explained LD reference/method and scoring weight construction.
- [ ] Reported score IDs, version, harmonization, missingness and matching rates.
- [ ] Prevented overlap/leakage between discovery, tuning and evaluation data.
- [ ] Reported independent predictive performance, uncertainty and limitations.
- [ ] Documented consent, governance, participant privacy and nondiscrimination.

---

## 13. What changes for a real large cohort

The principles remain identical when moving from 80 synthetic individuals to 10,000 or 100,000 human genomes. The engineering and study-design demands increase sharply.

### Data engineering

- Large numbers of BAM/CRAM, gVCF and cohort genotype files require robust storage and workflow orchestration.
- Perform joint genotyping/chromosome or interval sharding with traceable reference and software versions.
- Maintain sample identifiers and secure phenotype joins; use access-controlled, pseudonymized research data.
- Separate **raw variant calls**, **QC'd analysis genotypes**, **GWAS results**, **PGS scoring files** and **evaluation data**.

### Statistical design

- Standardize VO₂max units and testing protocols. Where relevant, distinguish absolute VO₂ (L/min) from body-mass-normalized VO₂max (mL/kg/min).
- Consider age, sex, training exposure, ancestry, relatedness, recruitment site and batch effects.
- Quantify statistical power per trait and subgroup **before** promising new discoveries.
- Validate independently; consider whether external population GWAS weights transfer to the cohort.
- Be careful about using genotype-based scores for individual medical or sporting decisions: available evidence often does not support them.

**Especially for sport-related traits:** environment, training history, selection and measurement variation can be dominant. Genomics is one layer of information, not a deterministic ranking system.

---

## 14. Common mistakes and troubleshooting

| Mistake / question | Correct interpretation or fix |
|---|---|
| “I have one individual's VCF and VO₂max, why no beta?” | A beta is estimated by comparing genotype–phenotype pairs **across people**. Apply external weights if appropriate; one person cannot support GWAS. |
| “PLINK conversion calculated beta, right?” | No. Conversion creates genotype files; `--glm` performs association analysis. |
| “Every SNP has a biological effect because GWAS gives a beta.” | No. Regression returns an estimate even for noise. Look at SE, P, replication and model diagnostics. |
| “The ALT allele always has the GWAS beta.” | No. Check `A1` / tested allele and harmonize against scoring alleles. |
| “A missing VCF site is 0/0.” | Not necessarily. Missing coverage or a variant-only file can hide uncertainty; confirm callable reference positions or imputed genotypes. |
| “PC1 means first/strongest performance gene.” | No. PC1 is an axis of variation across many genomic markers. |
| “LD pruning and LD clumping are the same.” | No. Pruning commonly prepares independent markers for PCA; clumping also uses GWAS strength to select markers for some PGS methods. |
| “Use raw GWAS betas from all SNPs as the score.” | Usually not appropriate: LD, noise, overfitting and portability matter. |
| “Higher PGS means better athlete.” | No. Predictive meaning depends on the trait, weights and independent validation. |
| “Good performance on the GWAS cohort proves accuracy.” | No. Evaluate using truly held-out or independent individuals. |
| “Can I upload real human VCFs to a public GitHub repo?” | Do **not** do this without explicit approvals, consent, governance and lawful data-sharing authorization. Genomes are highly identifying. |

---

## 15. Glossary

| Term | Beginner-friendly meaning |
|---|---|
| **NGS / WGS** | Sequencing DNA; WGS aims to survey the entire genome. |
| **FASTQ** | Raw reads and base-quality scores from sequencing. |
| **BAM / CRAM** | Sequencing reads aligned to a reference genome. |
| **gVCF** | A sample-level variant representation with information about reference regions for joint genotyping. |
| **Joint genotyping** | Calling genotypes consistently at cohort sites across multiple samples. |
| **VCF** | Tabular format holding genomic variants and often per-sample genotypes. |
| **REF / ALT** | Reference and alternate allele at a site. |
| **Genotype** | The alleles a person carries at a position. |
| **Dosage** | Expected or observed number of copies of a specified allele, often 0–2. |
| **Phenotype** | An observed/measured trait, such as VO₂max. |
| **Covariate** | An additional variable adjusted for in statistical analysis. |
| **GWAS** | Tests association between each genetic variant and a trait across people. |
| **Beta / effect size** | Estimated change in phenotype per additional tested allele (for a linear GWAS). |
| **Standard error (SE)** | Estimated uncertainty of an effect estimate. |
| **P-value** | Measure of statistical evidence against a specified null model. |
| **PCA / PC1 / PC2** | Method/coordinates summarizing broad patterns of genomic variation. |
| **Population stratification** | Population structure that can produce confounded genotype–trait associations. |
| **LD** | Statistical non-independence of alleles at genetic positions. |
| **r²** | Squared correlation often used to summarize pairwise LD for biallelic variants. |
| **LD pruning** | Selection of less-correlated markers, often before PCA. |
| **LD clumping** | P-value-guided selection of representative markers from associated LD regions. |
| **PGS** | Weighted sum of many genetic variants for a trait. |
| **PRS** | PGS focused on disease liability/risk; not automatically an absolute risk. |
| **Score harmonization** | Making scoring alleles and genotype alleles refer to exactly the same variants/directions. |
| **Overfitting** | Fitting noise in training data that does not generalize. |
| **Independent validation** | Evaluating a model on people who did not contribute to developing/tuning it. |
| **Incremental R²** | Increase in explained variation attributable to adding PGS to a baseline model, ideally out of sample. |
| **Cross-ancestry portability** | How well a score developed in one population works in another. |

---

## 16. Learning resources and references

**Software and technical documentation**

- [PLINK2: Standard data input](https://www.cog-genomics.org/plink/2.0/input) — VCF and phenotype/covariate files.
- [PLINK2: Association analysis (`--glm`)](https://www.cog-genomics.org/plink/2.0/assoc) — how variant effects are estimated.
- [PLINK2: Population stratification / PCA](https://www.cog-genomics.org/plink/2.0/strat) — ancestry PCs and cautions.
- [PLINK2: Linear scoring (`--score`)](https://www.cog-genomics.org/plink/2.0/score) — calculating a PGS and interpreting sums versus averages.
- [GATK: Joint genotyping using HaplotypeCaller gVCFs](https://gatk.broadinstitute.org/hc/en-us/articles/360035890411-Calling-variants-on-cohorts-of-samples-using-the-HaplotypeCaller-in-GVCF-mode) — from per-sample gVCF to cohort calls.
- [PGS Catalog](https://www.pgscatalog.org/) and [scoring file documentation](https://www.pgscatalog.org/downloads/) — published scores and effect alleles.
- [PGS Catalog Calculator](https://pgsc-calc.readthedocs.io/en/latest/) — supported inputs, allele matching and scoring reports.

**Scientific methods**

- Choi SW, Mak TSH, O'Reilly PF. **Tutorial: a guide to performing polygenic risk score analyses.** *Nature Protocols* (2020). DOI: [10.1038/s41596-020-0353-1](https://doi.org/10.1038/s41596-020-0353-1).
- Price AL et al. **Principal components analysis corrects for stratification in genome-wide association studies.** *Nature Genetics* (2006). DOI: [10.1038/ng1847](https://doi.org/10.1038/ng1847).

### Key takeaway

> **VCF tells us which alleles each person carries. Phenotype records tell us what was measured. GWAS estimates per-variant associations across people. Covariates and ancestry PCs help address confounding. LD describes variant correlation. PGS combines carefully derived variant weights, and independent evaluation tells us whether the score is useful.**

---

**Data and ethics note:** The repository contains only generated synthetic data and must remain separate from any actual individual-level human genome or phenotype records. Publicly sharing research-ready summary statistics or model metadata may still require institutional approval and appropriate disclosure controls.
