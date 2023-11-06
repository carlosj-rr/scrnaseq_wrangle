import scanpy as sp
import pandas as pd
import scipy as sc
import anndata as ad
import numpy as np

def col_matcher(query,reference):
    for i in range(0,len(reference)):
        if query in reference[i]:
            ans=i
    return(ans)

def percent(counts_array: int,num_cells: int) -> float:
    return(100*(counts_array/num_cells))

def col_nc_representation(column: pd.core.series.Series, nc_idcs: list) -> np.ndarray:
    col_ON=np.where(column > 0)
    return(np.intersect1d(col_ON,nc_idcs))

def check_inputs(ad_filename: str, ncgenes_filename: str):
    if ".loom" in ad_filename:
        loom=True
        h5ad=False
    elif ".h5ad" in ad_filename:
        loom=False
        h5ad=True
    else:
        raise TypeError(f"Anndata {ad_filename} filetype not recognized.\n\nPlease make sure it is either '.loom' or '.h5ad'.")
    # check if the list of nc gene markers is in CSV:
    if ".csv" in ncgenes_filename:
        csv=True
        tsv=False
    elif ".tsv" in ncgenes_filename or ".tab" in ncgenes_filename:
        csv=False
        tsv=True
    else:
        raise TypeError(f"List of nerve cell gene markers in an unknown format.\n\nPlease make sure it is in either '.csv' or '.tsv. format.")

def data_importer(ad_filename: str, ncgenes_filename: str) -> (ad._core.anndata.AnnData, np.ndarray):
    if ".loom" in ad_filename:
        ds_ad=sp.read_loom(ad_filename)
    elif ".h5ad" in ad_filename:
        ds_ad=sp.read_h5ad(ad_filename)
    if ".csv" in ncgenes_filename:
        nc_markers=np.array(pd.read_csv(ncgenes_filename, sep=",")).flatten()
    elif ".tsv" in ncgenes_filename:
        nc_markers=np.array(pd.read_csv(ncgenes_filename, sep="\t")).flatten()
    return(ds_ad,nc_markers)

def counts_table_extractor(ad_obj: ad._core.anndata.AnnData) -> pd.DataFrame:
    if type(ad_obj.X) == sc.sparse._csr.csr_matrix:
        out_tab=ad_obj.X.todense()
    else:
        out_tab=ad_obj.X
    colnames=ad_obj.var.iloc[:,0]
    rownames=ad_obj.obs.index
    out_tab=pd.DataFrame(out_tab,columns=colnames,index=rownames)
    return(out_tab)

def count_nc_genes(ad_filename: str, ncgenes_filename: str) -> pd.DataFrame:
    # check if anndata object is in .loom or in .h5ad format:
    check_inputs(ad_filename, ncgenes_filename)

    # import the input data:
    ds_ad,nc_markers=data_importer(ad_filename, ncgenes_filename) #function ensures right datatype output.

    # extract the tables - make sure the output table is a dense matrix, turn it into a pandas DF with the correct column and row names.
    ds_data=counts_table_extractor(ds_ad)

    l=[]
    for k in nc_markers:
        l.append(col_matcher(k,ds_data.columns))
    l.sort()

    list1=[]
    tot_cells=len(ds_data.columns)
    for col in range(0,tot_cells):
        val=len(col_nc_representation(ds_data.iloc[:,col],l))
        list1.append(val)
        print(f"Finished column {col} of {tot_cells}")
       
    gt0=pd.DataFrame(list(zip(ds_data.columns,list1)),columns=("barcode","counts"))
    
    return(gt0,ds_data)

def nc_gene_set_and_counts(cell_by_nc_marker_counts: pd.DataFrame, ds_data: pd.DataFrame) -> pd.DataFrame:
    gt0=cell_by_nc_marker_counts
    unique, count=np.unique(gt0["counts"], return_counts=True)
    max_num_ncmarkers=unique[-1]
    num_cells_max_exp=count[-1]
    nc_cells=np.array(gt0["barcode"][np.where(gt0["counts"] == max_num_ncmarkers)[0]])
    out_arr=np.array([],dtype=str)
    for barcode in nc_cells:
        genes_on_in_cell=ds_data.columns[np.where(ds_data.loc[[barcode]] > 0)[1]]
        out_arr=np.append(out_arr,genes_on_in_cell)

    unique,counts=np.unique(out_arr,return_counts=True)
    perc_occup=percent(counts,num_cells_max_exp)
    out_df=pd.DataFrame(list(zip(unique,counts,perc_occup)),columns=["ID","Counts","Percent"])
    return(out_df)