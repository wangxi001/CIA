{\rtf1\ansi\ansicpg1252\cocoartf2708
\cocoatextscaling0\cocoaplatform0{\fonttbl\f0\fswiss\fcharset0 Helvetica;}
{\colortbl;\red255\green255\blue255;}
{\*\expandedcolortbl;;}
\margl1440\margr1440\vieww24220\viewh8400\viewkind0
\pard\tx720\tx1440\tx2160\tx2880\tx3600\tx4320\tx5040\tx5760\tx6480\tx7200\tx7920\tx8640\pardirnatural\partightenfactor0

\f0\fs24 \cf0 import pandas as pd\
import numpy as np\
import matplotlib.pyplot as plt\
import scipy\
from matplotlib.pyplot import figure\
from pybedtools import BedTool\
from sklearn.metrics import average_precision_score, precision_recall_curve\
\
def load_domain(tad_name, ccd_name):\
    #HiC TAD\
    colnames = ['chr','start','end']\
    df_hic_tad = pd.read_csv(tad_name, sep='\\t',header = None,names = colnames)\
    df_hic_tad = df_hic_tad.sort_values(['chr','start','end'])\
    df_hic_tad['value'] = 1\
    pybed_hic_tad = BedTool.from_dataframe(df_hic_tad)\
\
    #CCD\
    colnames = ['chr','start','end','PET']\
    df_ccd = pd.read_csv(ccd_name, sep='\\t',header = None,names = colnames)\
    df_ccd = df_ccd.sort_values(['chr','start','end'])\
    pybed_ccd = BedTool.from_dataframe(df_ccd)\
    \
    return pybed_hic_tad, pybed_ccd\
\
def load_ctcf_loop(ctcf_loop_name, dis_cutoff = 500000, pet_cutoff = 0):\
    #CTCF loop list\
    colnames = ['chr','start','end','PET']\
    df_ctcf_loop = pd.read_csv(ctcf_loop_name,sep='\\t',header = None,names = colnames)\
    df_ctcf_loop['length'] = df_ctcf_loop['end'] - df_ctcf_loop['start']\
    df_ctcf_loop = df_ctcf_loop.sort_values(['chr','start','end'])\
    pybed_ctcf_loop = BedTool.from_dataframe(df_ctcf_loop[(df_ctcf_loop['length'] < dis_cutoff) & (df_ctcf_loop['PET'] > pet_cutoff)])\
    \
    return pybed_ctcf_loop\
\
def load_hic_loop(hic_loop_name):\
    colnames = ['chr','start','end','PET']\
    df_hic_loop = pd.read_csv(hic_loop_name,sep='\\t',header = None,names = colnames)\
    df_hic_loop['length'] = df_hic_loop['end'] - df_hic_loop['start']\
    df_hic_loop = df_hic_loop.sort_values(['chr','start','end'])\
    pybed_hic_loop = BedTool.from_dataframe(df_hic_loop)\
    \
    return pybed_hic_loop\
\
def load_ctcf_binding(ctcf_binding_name):\
    colnames = ['TargetGene','l1','l2','s1','s2','s3']\
    df_ctcf_binding = pd.read_csv(ctcf_binding_name,sep='\\t',header = None,names = colnames)\
    \
    return df_ctcf_binding\
\
def load_promoter(prom_name):\
    colnames = ['chr_gene','start_gene','end_gene','TargetGene','N','strand']\
    df_prom = pd.read_csv(prom_name,sep='\\t',header = None,names = colnames)\
    df_prom = df_prom.sort_values(['chr_gene','start_gene','end_gene']).reset_index(drop=True)\
    pybed_prom = BedTool.from_dataframe(df_prom[['chr_gene','start_gene','end_gene']])\
    \
    return df_prom, pybed_prom\
\
def load_all_pair(all_pair_name):\
    colnames = ['chr','start','end','TargetGene','activity','hic_contact','ABC_Score']\
    df_all_pair = pd.read_csv(all_pair_name,sep='\\t',header = None,names = colnames)\
\
    return df_all_pair \
\
def load_enhancer(enh_name):   \
    colnames = ['chr','start','end','id']\
    df_enh = pd.read_csv(enh_name,sep='\\t',header = None,names = colnames)\
    \
    return df_enh\
\
def compute_3d_features(df):\
    pybed_ep = BedTool.from_dataframe(df2[['chr','left','right']])\
    df['pet_outside'] = pybed_ep.map(pybed_ctcf_loop,c=4,o='sum',f=0.999999999).to_dataframe()['name'].replace('\\.','0', regex=True).astype(int)\
    df['pet_cross'] = pybed_ep.map(pybed_ctcf_loop,c=4,o='sum',f=0.1).to_dataframe().iloc[:,3].replace('\\.','0', regex=True).astype(int) - df2['pet_outside']\
\
    df['hicloop_outside'] = pybed_ep.map(pybed_hic_loop,c=4,o='sum',f=0.999999999).to_dataframe()['name'].replace('\\.','0', regex=True).astype(int)\
    df['hicloop_cross'] = pybed_ep.map(pybed_hic_loop,c=4,o='sum',f=0.1).to_dataframe().iloc[:,3].replace('\\.','0', regex=True).astype(int) - df2['hicloop_outside']\
\
    df['inTAD'] = pybed_ep.map(pybed_hic_tad,c=4,o='sum',f=0.999999999).to_dataframe()['name'].replace('\\.','0', regex=True).astype(int)\
    df['inCCD'] = pybed_ep.map(pybed_ccd,c=4,o='sum',f=0.999999999).to_dataframe()['name'].replace('\\.','0', regex=True).astype(int)\
\
    df['prom_CTCF'] = df[['TargetGene']].merge(df_p_ctcf_binding,how='left')['s3']\
    df['enh_CTCF'] = df[['chr','start']].merge(df_e_ctcf_binding,how='left')['s3']\
    \
    return df\
\
cell='GM12878'\
\
tad_name = '/users/wxi1/Desktop/ABCmodel/FeatureSet/'+cell+'_intactHiC_TAD_hg38.bed'\
ccd_name = '/users/wxi1/Desktop/ABCmodel/FeatureSet/'+cell+'_CTCF_r1_hg38_1kb_coverage_CCD_merged.bed'\
\
if(cell == 'K562'):\
    ctcf_loop_name = '/users/wxi1/Desktop/ENCODE_ChIAPET/CTCF_ChIAPET/K562_CTCF_merged_loop_hg38_loop.bed'\
elif(cell == 'GM12878'):\
    ctcf_loop_name = '/users/wxi1/Desktop/ABCmodel/FeatureSet/GM12878_CTCF_PET_clusters_hg38.bed'\
\
hic_loop_name = '/users/wxi1/Desktop/ABCmodel/FeatureSet/'+cell+'_intactHiC_loops_hg38.bed'\
\
pybed_hic_tad, pybed_ccd = load_domain(tad_name, ccd_name)\
pybed_ctcf_loop = load_ctcf_loop(ctcf_loop_name)\
pybed_hic_loop = load_hic_loop(hic_loop_name)\
\
p_ctcf_binding_name = '/users/wxi1/Desktop/ABCmodel/FeatureSet/'+cell+'_EnhancerPredictionsAllPutative_short_promoters_bw.txt'\
e_ctcf_binding_name = '/users/wxi1/Desktop/ABCmodel/FeatureSet/'+cell+'_EnhancerPredictionsAllPutative_short_enhancers_bw_all.txt'\
\
df_p_ctcf_binding = load_ctcf_binding(p_ctcf_binding_name)\
df_e_ctcf_binding = load_ctcf_binding(e_ctcf_binding_name)\
\
prom_name = '/users/wxi1/Desktop/ABCmodel/RefSeqCurated.170308.bed.CollapsedGeneBounds.TSS500bp.hg38.bed'\
enh_name = '/users/wxi1/Desktop/ABCmodel/FeatureSet/'+cell+'_EnhancerPredictionsAllPutative_short_enhancers.bed'\
all_pair_name = '/users/wxi1/Desktop/ABCmodel/FeatureSet/'+cell+'_EnhancerPredictionsAllPutative_short.txt'\
\
df_prom, pybed_prom = load_promoter(prom_name)\
df_enh = load_enhancer(enh_name)\
df2 = load_all_pair(all_pair_name)\
\
df_prom['prom_PET'] = pybed_prom.map(pybed_ctcf_loop,c=4,o='sum',f=0.9999).to_dataframe().iloc[:,3].replace('\\.','0', regex=True).astype(int)\
df_prom['prom_hicloop'] = pybed_prom.map(pybed_hic_loop,c=4,o='sum',f=0.9999).to_dataframe().iloc[:,3].replace('\\.','0', regex=True).astype(int)\
\
df2[['chr_gene','start_gene','end_gene','prom_PET','prom_hicloop']] = df2[['TargetGene']].merge(df_prom,how='left').iloc[:,[1,2,3,6,7]]\
\
df2 = df2.dropna()\
df2[['start_gene','end_gene','prom_PET','prom_hicloop']] = df2[['start_gene','end_gene','prom_PET','prom_hicloop']].astype(int)\
\
df2['left'] = df2[['start','start_gene']].min(axis=1)\
df2['right'] = df2[['end','end_gene']].max(axis=1)\
\
df2['distance'] = df2['right'] - df2['left']\
\
df2 = df2.sort_values(['chr','left','right']).reset_index().drop('index',1)\
\
df2 = compute_3d_features(df2)\
\
col_out = ['chr','start','end','TargetGene','activity','hic_contact','ABC_Score','pet_outside','pet_cross','prom_PET','inCCD','hicloop_outside','hicloop_cross','prom_hicloop','inTAD']\
df2[col_out].to_csv('/users/wxi1/Desktop/ABCmodel/FeatureSet/'+cell+'_EnhancerPredictionsAllPutative_CTCF_HiC.txt',sep='\\t',index=False)\
\
cell = 'GM12878'\
df2 = pd.read_csv('/users/wxi1/Desktop/ABCmodel/FeatureSet/'+cell+'_EnhancerPredictionsAllPutative_CTCF_HiC.txt',sep='\\t')\
\
def compute_CIA(df):\
    df['CIA_Score'] = 100 * df['activity'] * (df['pet_outside']+0.1) / (df['pet_cross']+1)\
    df_CIA_sum = df.groupby(['TargetGene'])[['CIA_Score']].sum()\
    df['CIA_Score_sum'] = df[['TargetGene']].merge(df_CIA_sum.reset_index(),how='left').iloc[:,1]\
    df['CIA_Score_norm'] = round(df['CIA_Score']/df['CIA_Score_sum'], 5)\
    df['CIA_Score'] = round(df['CIA_Score'],3)\
    \
    return df\
\
df2 = compute_CIA(df2)\
\
col_out = ['chr','start','end','TargetGene','activity','hic_contact','ABC_Score','pet_outside','pet_cross','prom_PET','inCCD','hicloop_outside','hicloop_cross','prom_hicloop','inTAD','CIA_Score','CIA_Score_norm']\
df2[col_out].to_csv('/users/wxi1/Desktop/ABCmodel/FeatureSet/'+cell+'_EnhancerPredictionsAllPutative_CTCF_HiC_CIA.txt',sep='\\t',index=False)\
}