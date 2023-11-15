import scanpy as sp
import pandas as pd
import scipy as sc
import anndata as ad
import numpy as np
import argparse

#needed
parser = argparse.ArgumentParser()
parser.add_argument("-sc_data","--scrnaseq_input_dataset", help="AnnData or CSV/TSV UMI counts table")
parser.add_argument("-ns_markers","--nervous_system_markers", help="CSV file with all nerve cell markers for species")
parser.add_argument("-o","--output_filename", help="Filename for a table of gene IDs from nervous system cells")

args=parser.parse_args()

in_data_file=args.scrnaseq_input_dataset
in_ns_markers_file=args.nervous_system_markers
out_table=args.output_filename

#needed
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
    print("Input data in a readable format (Loom/H5AD/CSV/TSV)")
#needed
def correct_axes_orientation(df: pd.DataFrame, nc_markers: np.ndarray) -> pd.DataFrame:
    if np.any(np.in1d(nc_markers,df.columns)):
        df=df.transpose()
    return(df)

#needed
def data_importer(ad_filename: str, ncgenes_filename: str):
    if ".loom" in ad_filename:
        ds_ad=sp.read_loom(ad_filename)
    elif ".h5ad" in ad_filename:
        ds_ad=sp.read_h5ad(ad_filename)
    if ".csv" in ncgenes_filename:
        nc_markers=np.array(pd.read_csv(ncgenes_filename, sep=",")).flatten()
    elif ".tsv" in ncgenes_filename:
        nc_markers=np.array(pd.read_csv(ncgenes_filename, sep="\t")).flatten()
    return(ds_ad,nc_markers)

#needed
def counts_table_extractor(ad_obj: ad._core.anndata.AnnData) -> pd.DataFrame:
    if type(ad_obj.X) == sc.sparse._csr.csr_matrix:
        out_tab=ad_obj.X.todense()
    else:
        out_tab=ad_obj.X
    colnames=ad_obj.var.iloc[:,0]
    rownames=ad_obj.obs.index
    out_tab=pd.DataFrame(out_tab,columns=colnames,index=rownames)
    return(out_tab)

#OF number 1
def count_nc_genes(ad_filename: str, ncgenes_filename: str) -> (pd.Series,pd.DataFrame):
    # check if anndata object is in .loom or in .h5ad format:
    check_inputs(ad_filename, ncgenes_filename)

    # import the input data:
    ds_ad,nc_markers=data_importer(ad_filename, ncgenes_filename) #function ensures right datatype output.
    print("Expression data and nervous system markers imported successfully.")
    #Check whether the imported table is an anndata object or not if it is, run 'counts_table_extractor'
    if type(ds_ad) == pd.DataFrame:
        ds_data=correct_axes_orientation(ds_ad,nc_markers)
    elif type(ds_ad) == ad.AnnData:
        # extract the tables - make sure the output table is a dense matrix, turn it into a pandas DF with the correct column and row names.
        ds_data=correct_axes_orientation(counts_table_extractor(ds_ad),nc_markers)
    print("Expression data has the correct axis orientation.")
    nc_gene_ss=ds_data.loc[nc_markers,:]
    gt0=np.sum((nc_gene_ss > 0).astype(int),axis=0)
    
    return(gt0,ds_data)

#OF 2
def top_tenperc_nc_gene_ids_extractor(nc_markers_representation: pd.Series,ds_data: pd.DataFrame, outfile: str) -> None:
    unique3, counts3 = np.unique(nc_markers_representation,return_counts=True)
    tot_cells=nc_markers_representation.shape[0]
    ten_perc=int(np.round(tot_cells*0.1))
    inverted_cumsum=np.cumsum(np.flip(counts3))
    smallest_diff_idx=np.where(abs(inverted_cumsum-ten_perc) == abs(inverted_cumsum-ten_perc).min())[0]
    ten_perc_nc_cell_threshold=np.flip(unique3)[smallest_diff_idx][0]
    nc_occupancy=unique3[np.where(unique3 == ten_perc_nc_cell_threshold)[0][0]:len(unique3)]
    print(f"Top 10% most nerve cell-like cells express {nc_occupancy.size} of a maximum of {unique3[-1]} nervous system markers expressed.")
    out_arr=np.array([],dtype=str)
    for occup in nc_occupancy:
        new=np.array(nc_markers_representation.index[np.where(nc_markers_representation == occup)[0]])
        out_arr=np.append(out_arr,new)

    ten_perc_nc_cells_df=ds_data.loc[:,out_arr]
    print(f"A total of {len(ten_perc_nc_cells_df.columns)} nerve cells were identified.")
    nc_genes=ten_perc_nc_cells_df.index[np.where(ten_perc_nc_cells_df.sum(axis=1) > 0)[0]]
    np.savetxt(outfile,nc_genes,delimiter=",",fmt='%s')
    print(f"Full list of nervous system genes saved to {outfile}")
    return

def main():
    gt0,ds_data=count_nc_genes(in_data_file,in_ns_markers_file)
    top_tenperc_nc_gene_ids_extractor(gt0,ds_data,out_table)
    return

main()