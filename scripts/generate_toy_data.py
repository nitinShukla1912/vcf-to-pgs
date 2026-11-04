#!/usr/bin/env python3
"""Generate synthetic, nonhuman-research example files for teaching VCF/GWAS/PGS.

All variants, measurements, effect sizes and participants are fictional. SNP
positions are illustrative coordinates and NOT asserted as GRCh38 variants.
"""
from __future__ import annotations
import csv
import random
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / 'data'
DATA.mkdir(parents=True, exist_ok=True)
RNG = random.Random(20261104)
N = 80
VARIANTS = [
    ('1', 1000, 'toy_snp01', 'A', 'G'),
    ('1', 1080, 'toy_snp02', 'C', 'T'),  # purposely correlated with toy_snp01
    ('1', 2000, 'toy_snp03', 'G', 'A'),
    ('1', 2500, 'toy_snp04', 'T', 'C'),
    ('1', 5000, 'toy_snp05', 'A', 'C'),
    ('1', 6000, 'toy_snp06', 'T', 'G'),
    ('2', 1200, 'toy_snp07', 'C', 'G'),
    ('2', 1900, 'toy_snp08', 'G', 'T'),
    ('2', 7000, 'toy_snp09', 'A', 'T'),
    ('3', 1500, 'toy_snp10', 'C', 'A'),
    ('3', 3900, 'toy_snp11', 'A', 'G'),
    ('3', 8000, 'toy_snp12', 'G', 'C'),
]

def dosage(p: float) -> int:
    return int(RNG.random() < p) + int(RNG.random() < p)

persons = []
for i in range(1, N+1):
    sid = f'ATH{i:03d}'
    # Two arbitrary simulated demographic groups, NOT genuine ancestry estimates.
    grp = int(RNG.random() < 0.5)
    age = RNG.randint(20, 41) + (3 if grp else 0)
    sex_code = int(RNG.random() < 0.5)
    sport_endurance = int(RNG.random() < (0.73 if grp else 0.37))
    g1 = dosage(0.62 if grp else 0.22)
    # When high correlation is desired, copy SNP 1 genotype most of the time.
    g2 = g1 if RNG.random() < 0.90 else dosage(0.58 if grp else 0.25)
    gs = [g1, g2]
    gs.extend(dosage((0.17 + 0.08 * (j % 5)) + (0.04 if grp else 0)) for j in range(2, len(VARIANTS)))
    # Known data generating rule: phenotype is associated with g1, age, sex,
    # group and training category, plus a relatively large random component.
    vo2 = (49 + 2.3*g1 - 0.22*(age-25) + 2.1*sex_code +
           4.2*sport_endurance + 2.3*grp + RNG.gauss(0, 4.0))
    persons.append({'id':sid, 'group':grp, 'age':age, 'sex':sex_code,
                    'sport':sport_endurance, 'vo2':round(vo2,2), 'gs':gs})

with (DATA/'toy_joint_called.vcf').open('w') as fh:
    fh.write('##fileformat=VCFv4.2\n')
    fh.write('##source=SYNTHETIC_EDUCATIONAL_DATA_NOT_REAL_GENOMES\n')
    fh.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
    fh.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t' +
             '\t'.join(p['id'] for p in persons) + '\n')
    call={0:'0/0',1:'0/1',2:'1/1'}
    for j, (chrom,pos,vid,ref,alt) in enumerate(VARIANTS):
        fh.write('\t'.join([chrom,str(pos),vid,ref,alt,'.','PASS','.','GT']+
                            [call[p['gs'][j]] for p in persons])+'\n')

with (DATA/'phenotypes.tsv').open('w',newline='') as fh:
    out=csv.writer(fh, delimiter='\t')
    out.writerow(['FID','IID','VO2max'])
    for p in persons:
        out.writerow(['0',p['id'],p['vo2']])

with (DATA/'covariates.tsv').open('w',newline='') as fh:
    out=csv.writer(fh,delimiter='\t')
    out.writerow(['FID','IID','Age','SexCode','ToyGroup','EnduranceSport'])
    for p in persons:
        out.writerow(['0',p['id'],p['age'],p['sex'],p['group'],p['sport']])

# Fictional pre-specified weights: they are NOT fitted from these 80 athletes.
# Deliberately leave SNP 2 out to avoid counting the near-duplicate LD signal.
weights=[('toy_snp01','G',2.3),('toy_snp03','A',-0.4),
         ('toy_snp05','C',0.8),('toy_snp07','G',0.25),
         ('toy_snp10','A',-0.3)]
with (DATA/'illustrative_weights.tsv').open('w',newline='') as fh:
    out=csv.writer(fh,delimiter='\t')
    out.writerow(['ID','A1','BETA'])
    out.writerows(weights)

print(f'Created {N} synthetic sample columns, {len(VARIANTS)} SNP rows.')
print('Files: data/toy_joint_called.vcf, data/phenotypes.tsv,')
print('       data/covariates.tsv, data/illustrative_weights.tsv')
