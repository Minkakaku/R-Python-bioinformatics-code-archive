sc.pl.stacked_violin(adata, var_names=marker_genes, groupby=groupby, show=False, swap_axes=True, dendrogram=False)
plt.savefig(os.path.join(output_path, f'{prefix}_violin.pdf'),bbox_inches='tight')
plt.savefig(os.path.join(output_path, f'{prefix}_violin.png'),bbox_inches='tight')

sc.pl.dotplot(adata_raw, var_names=marker_genes, groupby=['cell_type','group'], dendrogram=False, show=False,use_raw=True)
plt.savefig(os.path.join(output_path, f'{prefix}_bubble.pdf'), facecolor='white', bbox_inches="tight", dpi=300)
plt.savefig(os.path.join(output_path, f'{prefix}_bubble.png'), facecolor='white', bbox_inches="tight", dpi=300)

sc.pl.umap(adata, color=marker_genes, show=False, cmap='coolwarm')
save_fig(plt,file_name= f'{prefix}_umap',out_dir=output_path,fig_size={'w':7,'h':7})

count_tbl = pd.crosstab(adata.obs['group'], adata.obs['cell_type'])
count_tbl.to_csv(os.path.join(output_path, "counts.xls"), sep='\t')

prop_tbl = count_tbl.div(count_tbl.sum(axis=1), axis=0).round(2)
prop_tbl.to_csv(os.path.join(output_path, "prop.xls"), sep='\t')

adata=adata[~adata.obs['cell_type'].isin(['Neutrophils']), :]


adata =  adata_raw.concatenate(bdata_map)
adata.raw = adata
adata.raw is None
TRUE
adata = adata.raw.to_adata()
# sc.pp.normalize_total(adata, target_sum=1e4)
# sc.pp.log1p(adata)
adata.uns['log1p']['base'] = 2
sc.pp.highly_variable_genes(self.adata,n_top_genes=2000)
adata.raw = adata