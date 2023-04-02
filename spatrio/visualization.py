import pandas as pd
import scanpy as sc
import numpy as np

def spatial_plot(adata,info,color_by="type",size=None,mode =1):
    """
    Plots original/mapped spatial coordinates.
    
    Args:
        adata: Spatial dataset to be plotted.
        info: A DataFrame containing spot/cell id and meta information.
        color_by: Which obseravation in adata is used to determine the color of each point in the plot.
        size: The size of each point in the plot.
        mode: Select the type of plot data. If mode=1, plot spot dataset. If mode=2, plot single-cell dataset.

    """
    
    if mode ==1 :
        coor = info.copy()
        coor['id'] = coor.index.tolist()
        coor.columns = ['x','y','id']
    if mode ==2 :
        coor = info.copy()
        coor = coor[['Cell_xcoord','Cell_ycoord','cell']]
        coor.columns = ['x','y','id']

    adata = adata[adata.obs.index.isin(coor.id)]
    left_data = pd.DataFrame(adata.obs.index)
    left_data.columns = ['id']
    right_data = coor
    coordiante = pd.merge(left_data,right_data,on='id',how="left")
    coordiante.drop(['id'],axis=1,inplace=True)
    adata.obsm['spatial'] = np.array(coordiante)
    sc.pl.embedding(adata, basis="spatial", color=color_by,size=size)