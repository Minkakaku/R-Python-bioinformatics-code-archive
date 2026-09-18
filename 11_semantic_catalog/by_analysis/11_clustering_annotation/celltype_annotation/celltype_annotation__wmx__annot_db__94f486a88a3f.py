# https://ftp.ncbi.nlm.nih.gov/gene/DATA/

import os
import requests
import pandas as pd


def go_annot(gene2gofile, geneinfofile, orgnism):
    gene2go = pd.read_csv(gene2gofile, sep='\t', low_memory=False)
    gene_info = pd.read_csv(geneinfofile, sep='\t', usecols=['GeneID', 'Symbol', 'description', 'type_of_gene'], low_memory=False)

    if orgnism == 'hsa':
        gene2go = gene2go[gene2go['#tax_id']==9606]
    elif orgnism == 'mmu':
        gene2go = gene2go[gene2go['#tax_id']==10090]
    else:
        print('Only human and mouse support now.')

    gene2go = gene2go.groupby('GeneID')['GO_ID', 'GO_term', 'Category'].agg(lambda x: ';'.join(x)).reset_index()
    gene2go['Category'] = gene2go['Category'].str.replace('Process', 'BP').str.replace('Component', 'CC').str.replace('Function', 'MF')

    gene2go = pd.merge(gene_info, gene2go, on='GeneID', how='outer')

    return gene2go
    
def kegg_annot(orgnism):
    def kegg_api(url, col_split=None):
        response = requests.get(url)
        tmp = []
        if response.status_code == 200:
            data = response.text.strip().split('\n')
            for line in data:
                line = line.split('\t')
                if col_split == 1:
                    tmp.append([line[0].split(':')[1], line[1]])
                elif col_split == 2:
                    tmp.append([line[0], line[1].split(':')[1]])
                else:    
                    tmp.append([line[0], line[1]])
        else:
            print("API request failed.")      

        return pd.DataFrame(tmp)

    url = f'https://rest.kegg.jp/conv/{orgnism}/ncbi-geneid'
    geneid2keggid = kegg_api(url, 1)
    geneid2keggid.columns = ['GeneID', 'KEGG id']

    url = f'https://rest.kegg.jp/link/pathway/{orgnism}'
    geneid2pathway = kegg_api(url, 2)
    geneid2pathway.columns = ['KEGG id', 'path']

    url = f'https://rest.kegg.jp/list/pathway/{orgnism}'
    pathway2desc = kegg_api(url)
    pathway2desc.columns = ['path', 'path description']

    url = f'https://rest.kegg.jp/link/ko/{orgnism}'
    keggid2KO = kegg_api(url, 2)
    keggid2KO.columns = ['KEGG id', 'KO']

    tmp = pd.merge(geneid2keggid, keggid2KO, on='KEGG id')
    gene2kegg = pd.merge(pd.merge(tmp, geneid2pathway, on='KEGG id'), pathway2desc, on='path')

    gene2kegg = gene2kegg.groupby('GeneID')['KO', 'path', 'path description'].agg(lambda x: ';'.join(set(x))).reset_index()

    gene2kegg['GeneID'] = gene2kegg['GeneID'].astype(int)

    return gene2kegg

def annot_combine(gene2go, gene2kegg, out_put):
    annotation = pd.merge(gene2go, gene2kegg, on='GeneID', how='outer')
    annotation = annotation[(annotation['type_of_gene']!='biological-region')]
    
    # # Proteomics use below
    # annotation = annotation[(annotation['type_of_gene']=='protein-coding')]

    annotation.to_csv(out_put, sep='\t', index=None)


if __name__ == '__main__':
    # Human
    gene2go = go_annot('./NCBI/ncbi_gene2go', './NCBI/Homo_sapiens.gene_info', 'hsa')
    gene2kegg = kegg_annot('hsa')
    annot_combine(gene2go, gene2kegg, '9606_human_go_kegg_annotation.xls')
    
    # Mouse
    gene2go = go_annot('./NCBI/ncbi_gene2go', './NCBI/Mus_musculus.gene_info', 'mmu')
    gene2kegg = kegg_annot('mmu')
    annot_combine(gene2go, gene2kegg, '10090_mouse_go_kegg_annotation.xls')
