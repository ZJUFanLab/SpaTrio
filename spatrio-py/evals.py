import pandas as pd
import numpy as np
from anndata import AnnData
from .helper import intersect, find_marker, extract_exp, process_anndata
from scipy.stats import pearsonr


def pearson_cal(spax_map: pd.DataFrame,
    ann1: AnnData,
    ann2: AnnData,
    use_marker: bool = True,
    top_n: int = 100,
    run_marker_by: str = "type"
    ):
    adata1 = ann1.copy()
    adata2 = ann2.copy()
    common_genes = intersect(adata1.var.index, adata2.var.index)
    adata1 = adata1[:, common_genes]
    adata2 = adata2[:, common_genes]
    
    adata1 = process_anndata(adata1,scale=True,pca=False)
    adata2 = process_anndata(adata2,scale=True,pca=False)
    
    map = pd.DataFrame(spax_map[['spot','cell']])
    data1 = extract_exp(adata1)
    data2 = extract_exp(adata2)
    data2['cell'] = data2.index

    data2 = pd.merge(map, data2, on='cell',how="left")
    grouped = data2.groupby(['spot'])

    mapped_data = grouped[data2.columns.tolist()[2:]].sum()

    common_genes = intersect(data1.columns, mapped_data.columns)
    data1 = data1[common_genes]
    mapped_data = mapped_data[common_genes]

    if use_marker:
        #sc.pp.filter_cells(adata2_copy, min_genes=200)
        #sc.pp.filter_genes(adata2_copy, min_cells=3)
        marker1 = find_marker(ann1,top_n,maker_by=run_marker_by)
        marker2 = find_marker(ann2,top_n,maker_by=run_marker_by)
        marker1.extend(marker2)
        gene_list = np.unique(marker1).tolist()
        adata1 = adata1[:, gene_list]
        adata2 = adata2[:, gene_list]
        data1 = data1[gene_list]
        mapped_data = mapped_data[gene_list]

    common_index = np.intersect1d(data1.index, mapped_data.index, assume_unique=False, return_indices=False)
    data1 = data1.loc[common_index]
    mapped_data = mapped_data.loc[common_index]
    

    pearson_s = []
    for spot in data1.index.tolist():
        x=data1.loc[spot].values
        y=mapped_data.loc[spot].values
        pc = pearsonr(x,y)
        pearson_s.append(pc[0])
    pearson = np.mean(pearson_s)
    return pearson

def accuracy_cal(map_data: pd.DataFrame) :
    accuracy_data = map_data.copy()
    accuracy = (accuracy_data['spot_type']==accuracy_data['cell_type']).sum()/accuracy_data.shape[0]
    return accuracy
