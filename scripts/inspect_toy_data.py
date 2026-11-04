#!/usr/bin/env python3
"""Pure teaching example: compute GWAS coefficients, LD r2 and synthetic PGS.

Uses a textbook linear model. This is not a replacement for a validated
GWAS/PGS bioinformatics pipeline.
"""
from __future__ import annotations
import csv
from pathlib import Path
import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parents[1]
DATA=ROOT/'data'
OUT=ROOT/'results'
OUT.mkdir(parents=True,exist_ok=True)

with (DATA/'toy_joint_called.vcf').open() as fh:
    lines=[l.strip().split('\t') for l in fh if not l.startswith('##')]
header=lines[0]
ids=header[9:]
rows=lines[1:]
variant_ids=[r[2] for r in rows]
# G: rows are people, columns are SNP dosages
G=np.array([[int(v.split('/')[0])+int(v.split('/')[1]) for v in r[9:]]
            for r in rows], dtype=float).T
ph=pd.read_csv(DATA/'phenotypes.tsv',sep='\t').set_index('IID').loc[ids]
cov=pd.read_csv(DATA/'covariates.tsv',sep='\t').set_index('IID').loc[ids]
y=ph['VO2max'].to_numpy(dtype=float)

def coefficient(x:np.ndarray, y:np.ndarray, covariates:np.ndarray|None=None):
    """OLS beta of SNP column, with optional covariates (and intercept)."""
    n=len(y)
    x=np.asarray(x).reshape(-1,1)
    C=np.ones((n,1)) if covariates is None else np.column_stack([np.ones(n),covariates])
    X=np.column_stack([C,x])
    if np.linalg.matrix_rank(X) < X.shape[1]:
        return float('nan')
    return float(np.linalg.lstsq(X,y,rcond=None)[0][-1])

all_results=[]
C=cov[['Age','SexCode','ToyGroup','EnduranceSport']].to_numpy(dtype=float)
for j,vid in enumerate(variant_ids):
    all_results.append({'ID':vid, 'ALT':rows[j][4],
                        'BETA_unadjusted':coefficient(G[:,j],y),
                        'BETA_adjusted':coefficient(G[:,j],y,C)})
pd.DataFrame(all_results).to_csv(OUT/'toy_gwas_betas.tsv',sep='\t',index=False,float_format='%.4f')

r2=float(np.corrcoef(G[:,0],G[:,1])[0,1]**2)
print('--- One SNP, all athletes ---')
print('SNP = toy_snp01; number of synthetic athletes =',len(ids))
print('Unadjusted beta =',round(all_results[0]['BETA_unadjusted'],3))
print('Adjusted beta   =',round(all_results[0]['BETA_adjusted'],3))
print('--- LD ---')
print('Genotype correlation r², toy_snp01 vs toy_snp02 =',round(r2,3))

weights=pd.read_csv(DATA/'illustrative_weights.tsv',sep='\t')
indices={vid:j for j,vid in enumerate(variant_ids)}
score=np.zeros(len(ids))
for r in weights.itertuples(index=False):
    j=indices[r.ID]
    # This fixture's effect alleles equal the ALT allele; ALWAYS verify this in real data.
    assert r.A1==rows[j][4]
    score += G[:,j]*r.BETA
pd.DataFrame({'IID':ids,'toy_PGS_sum':score}).to_csv(OUT/'toy_pgs_scores.tsv',sep='\t',index=False,float_format='%.4f')
print('--- PGS ---')
print('ATH001 toy PGS =',round(float(score[0]),3))
print('Wrote results/toy_gwas_betas.tsv and results/toy_pgs_scores.tsv')
