
from pathlib import Path
import pandas as pd
import scanpy as sc
work=Path('/mnt/NexData_file/home/lyh/sc_RNAseq_aftersale/SC2606034175_B/20260622_cell_anno')
evd=work/'annotation_evidence'
markers = {
    'Pan_neuron': ['Snap25','Syt1','Rbfox3','Tubb3','Map2'],
    'Excitatory_neuron': ['Slc17a7','Camk2a','Satb2','Neurod6','Tbr1'],
    'Inhibitory_neuron': ['Gad1','Gad2','Slc32a1','Dlx1','Dlx2'],
    'Astrocyte': ['Aqp4','Aldh1l1','Slc1a3','Gfap','S100b','Aldoc'],
    'Oligodendrocyte': ['Mog','Mbp','Plp1','Mobp','Mag'],
    'OPC': ['Pdgfra','Cspg4','Vcan','Olig1','Olig2'],
    'Microglia': ['P2ry12','Tmem119','Cx3cr1','Aif1','Csf1r','Hexb'],
    'Macrophage_Monocyte': ['Lyz2','C1qa','C1qb','Tyrobp','Itgam','Mrc1'],
    'Endothelial': ['Pecam1','Cldn5','Kdr','Flt1','Slco1a4'],
    'Pericyte_VSMC': ['Pdgfrb','Rgs5','Kcnj8','Acta2','Tagln','Myl9'],
    'Ependymal': ['Foxj1','Tmem212','Dynlrb2','Pifo','Ccdc153'],
    'VLMC_Fibroblast': ['Col1a1','Col1a2','Dcn','Lum','Pdgfra'],
    'Choroid_plexus_epithelial': ['Ttr','Kcne2','Enpp2','Folr1','Aqp1'],
    'T_NK_cell': ['Cd3d','Cd3e','Trac','Nkg7','Klrb1c'],
    'B_Plasma_cell': ['Ms4a1','Cd79a','Cd74','Jchain','Mzb1'],
    'Erythroid': ['Hba-a1','Hba-a2','Hbb-bs','Hbb-bt','Alas2'],
    'Cycling': ['Mki67','Top2a','Cenpf','Ube2c','Anln'],
}
wide=pd.DataFrame({k: pd.Series(v) for k,v in markers.items()})
for name in ['Markers.csv','Marker.csv','Markersplot.csv']:
    wide.to_csv(work/name, index=False)

mapping = {
0:'Microglia', 1:'Astrocyte', 2:'Astrocyte', 3:'Activated_Microglia',
4:'Doublet_or_low_quality', 5:'Microglia', 6:'Activated_Microglia', 7:'Oligodendrocyte',
8:'Endothelial', 9:'Astrocyte', 10:'Oligodendrocyte', 11:'Astrocyte',
12:'Pericyte_VSMC', 13:'Endothelial', 14:'Immature_neuron', 15:'Macrophage_Monocyte',
16:'Endothelial', 17:'Ependymal', 18:'Choroid_plexus_epithelial', 19:'Astrocyte',
20:'Ependymal', 21:'Immune_low_quality', 22:'Doublet_or_low_quality', 23:'Doublet_or_low_quality',
24:'Astrocyte', 25:'Erythroid', 26:'Astrocyte', 27:'VLMC_Fibroblast', 28:'VLMC_Fibroblast'
}
pd.DataFrame({'cluster': list(mapping.keys()), 'celltype': list(mapping.values())}).to_csv(work/'Celltype.csv', index=False)
adata=sc.read_h5ad(work/'cluster/adata_cluster.h5ad')
rg=pd.read_csv(work/'marker/rank_gene_allcluster.csv')
rg['group']=rg['group'].astype(str)
scores=pd.read_csv(evd/'marker_panel_scores.csv')
scores['cluster']=scores['cluster'].astype(str)
counts=adata.obs['leiden'].astype(str).value_counts().to_dict()
sample_ct=pd.crosstab(adata.obs['leiden'].astype(str), adata.obs['sample'].astype(str))
rationales={
'Microglia':'High Hexb/P2ry12/Cx3cr1/C1q microglial program with strongest Microglia marker-panel score.',
'Activated_Microglia':'Microglial Cx3cr1/C1q/Ctss background plus immediate early/inflammatory genes such as Nfkbia, Cd83, Jun/Fos/Egr1.',
'Astrocyte':'Astrocyte-associated Slc1a2/Aldoc/Apoe/Aqp4/Slc4a4/Gria/Grid signatures and highest Astrocyte panel score.',
'Oligodendrocyte':'Mature oligodendrocyte myelin genes Mal/Mbp/Plp1/Mobp/Mag dominate.',
'Endothelial':'Vascular endothelial genes Flt1/Cldn5/Ptprb/Slco1a4/Abcb1a dominate.',
'Pericyte_VSMC':'Mural-cell markers Rgs5/Pdgfrb/Myl9/Vtn/Notch3 dominate.',
'Macrophage_Monocyte':'Border/perivascular macrophage markers Lyz2/Mrc1/Pf4/F13a1/Ms4a7 and high macrophage panel score.',
'Ependymal':'Ciliated ependymal genes Dnah/Cfap/Foxj1-related panel or Ccdc153/Dbi/Nnat/Vim pattern.',
'Choroid_plexus_epithelial':'Choroid epithelial genes Ttr/Enpp2/Folr1/Kcnj13/Prlr dominate.',
'Immature_neuron':'Neurodevelopmental genes Stmn1/Tubb5/Sox11/Sox4/Tuba1a with neuron panel support.',
'Immune_low_quality':'Immune-associated Cd52/Arhgap15 signal but dominated by ribosomal genes; retained as low-confidence immune-like cluster.',
'Doublet_or_low_quality':'Mixed or low-quality signature, including hemoglobin/Ttr/myelin/mitochondrial combinations and/or strong single-sample enrichment.',
'Erythroid':'Hemoglobin and erythroid genes Hba/Hbb/Alas2/Bpgm dominate.',
'VLMC_Fibroblast':'Meningeal/perivascular fibroblast genes Dcn/Col1a2/Colec12/Lum/Pdgfra or Mgp/Rspo3/Igfbp matrix program dominate.'
}
rows=[]
for cl in sorted(mapping, key=int):
    top=rg[rg['group']==str(cl)].head(12)['names'].astype(str).tolist()
    scrow=scores[scores['cluster']==str(cl)]
    top_scores=''
    if not scrow.empty:
        vals=scrow.iloc[0].drop('cluster').sort_values(ascending=False).head(4)
        top_scores='; '.join(f'{k}:{v:.3f}' for k,v in vals.items())
    sample_summary='; '.join(f'{k}:{int(v)}' for k,v in sample_ct.loc[str(cl)].sort_values(ascending=False).head(3).items())
    rows.append({
        'cluster': cl,
        'celltype': mapping[cl],
        'n_cells': int(counts.get(str(cl),0)),
        'top_rank_markers': ','.join(top),
        'marker_support_summary': top_scores,
        'sample_tissue_clues': sample_summary + '; Tissue=mouse brain',
        'rationale': rationales[mapping[cl]],
    })
pd.DataFrame(rows).to_csv(evd/'cluster_annotation_evidence.csv', index=False)
print('wrote Celltype.csv and cluster_annotation_evidence.csv')
