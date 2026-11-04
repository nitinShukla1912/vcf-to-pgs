# vcf-to-pgs

A growing, **beginner-friendly and practical** resource for learning genomic data analysis—from sequencing reads and variant calls to statistical genetics, reproducible workflows, and interpretation.

**Start with the illustrated, step-by-step tutorial:**

## [From NGS to GWAS and Polygenic Scores](docs/01-ngs-to-gwas-and-pgs.md)

This guide explains:

- How FASTQ → BAM/CRAM → joint-called VCF → PLINK files are connected.
- How **genotypes** (from a VCF) are matched with **phenotypes** (recorded measurements) by sample ID.
- How a GWAS estimates a **beta** for each SNP by comparing genotypes *across people*.
- Why **age, sex, ancestry PCs and other covariates** matter in association analysis.
- What **linkage disequilibrium (LD)** means, and the difference between LD pruning and LD-aware polygenic scoring.
- How to calculate a **polygenic score (PGS)** and, importantly, how to validate it without data leakage.
- What common PLINK2 commands do and where they fit into a real workflow.

## Overview

```text
Sequencing data            Study measurements
FASTQ files                VO2max, age, sex, etc.
    |                             |
BAM / CRAM                        |
    |                             |
Per-sample gVCF                   |
    |                             |
Joint-called cohort VCF           |
    |                             |
Genotype QC -> PLINK2             |
    |                             |
    +--- match samples by ID -----+
                |
      GWAS: beta/SE/P per SNP
                |
      LD-aware model development
                |
      PGS per person + validation
```

## Try the synthetic demonstration

The repository contains only **fictional** data: 80 simulated participants, 12 simulated SNPs, made-up VO₂max values, and illustrative effect weights. No real genomes or personal data are included.

```bash
python -m pip install -r requirements.txt
python scripts/generate_toy_data.py
python scripts/inspect_toy_data.py
```

**Outputs:** `results/toy_gwas_betas.tsv` and `results/toy_pgs_scores.tsv` (ignored by Git). The Python demonstration illustrates regression coefficients, covariate adjustment, genotype correlation (LD), and score calculation. Read the [complete tutorial](docs/01-ngs-to-gwas-and-pgs.md#10-run-the-synthetic-example) before interpreting them.

## Repository organization

```text
genomics-analysis-handbook/
├── README.md
├── docs/
│   └── 01-ngs-to-gwas-and-pgs.md
├── data/
│   ├── toy_joint_called.vcf
│   ├── phenotypes.tsv
│   ├── covariates.tsv
│   └── illustrative_weights.tsv
├── scripts/
│   ├── generate_toy_data.py
│   └── inspect_toy_data.py
├── requirements.txt
├── ROADMAP.md
└── CONTRIBUTING.md
```

## Future topics

See [ROADMAP.md](ROADMAP.md) for areas to expand: variant annotation, ancestry/relatedness, GWAS QC, fine-mapping, population genetics, rare variants, replication, and reproducible notebooks. Add new guides under `docs/` without disrupting the existing tutorial.

## Scientific limitations and privacy

This is an **educational repository**, not a clinically validated genomic interpretation tool. A statistically associated SNP is not necessarily causal, a PGS is not a diagnosis or fitness percentage, and prediction performance requires independent validation. Never commit raw human sequencing data, identifiable phenotype records, consent documents or sensitive metadata to a public repository.

## Contribute

Suggestions, corrections, and new tutorials are welcome. See [CONTRIBUTING.md](CONTRIBUTING.md). A reuse license has not been selected yet; a public GitHub repository alone does not grant an open-source license.
