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
parser.add_argument("-o","--output_prefix", help="Prefix to be attached to output files")
parser.add_argument("-ann","--emapper_annot",help="Filename to the *.annotations file from emapper.py.")

args=parser.parse_args()

in_data_file=args.scrnaseq_input_dataset
in_ns_markers_file=args.nervous_system_markers
prefix=args.output_prefix
annot_file=args.emapper_annot

#needed
def check_inputs(ad_filename: str, ncgenes_filename: str):
    if ".loom" in ad_filename or ".h5ad" in ad_filename or ".csv" in ad_filename:
        readable=True
    else:
        raise TypeError(f"Anndata {ad_filename} filetype not recognized.\n\nPlease make sure it is either '.loom' or '.h5ad'.")
    # check if the list of nc gene markers is in CSV:
    if ".csv" in ncgenes_filename or ".tsv" in ncgenes_filename or ".tab" in ncgenes_filename:
        readable=True
    else:
        raise TypeError(f"List of nerve cell gene markers in an unknown format.\n\nPlease make sure it is in either '.csv' or '.tsv. format.")
    print("Input data in a readable format (Loom/H5AD/CSV/TSV)")

#needed
def correct_axes_orientation(df: pd.DataFrame, nc_markers: np.ndarray) -> (pd.DataFrame: np.ndarray):
    if type(df.index) == pd.RangeIndex:
        print("imported indices are not gene IDs - assuming first column is the indices...")
        rownames=df.iloc[:,0]
        df=df.iloc[:,1:len(df.columns)-1]
        df.index=rownames
    orig_nc_marker_num=nc_markers.size
    nc_markers=np.intersect1d(nc_markers,df.index)
    new_num_nc_markers=nc_markers.size
    diff=orig_nc_marker_num-new_num_nc_markers
    if diff != 0:
        print(f"NC markers list had {diff} markers which were not in the scRNA-seq dataset. NC marker list corrected.")
    if np.any(np.in1d(nc_markers,df.columns)):
        df=df.transpose()
    return(df,nc_markers)

#needed
def data_importer(ad_filename: str, ncgenes_filename: str):
    if ".loom" in ad_filename:
        ds_ad=sp.read_loom(ad_filename)
    elif ".h5ad" in ad_filename:
        ds_ad=sp.read_h5ad(ad_filename)
    elif ".csv" in ad_filename:
        ds_ad=pd.read_csv(ad_filename,sep=",", index_col=0,header=0)
        print("scRNA-seq data read directly as a data frame, with the first column as indices and first row as column headers\nplease make sure this is correct: {ds_ad}")
    if ".csv" in ncgenes_filename:
        nc_markers=np.array(pd.read_csv(ncgenes_filename, sep=",")).flatten()
    elif ".tsv" in ncgenes_filename:
        nc_markers=np.array(pd.read_csv(ncgenes_filename, sep="\t")).flatten()
    return(ds_ad,nc_markers)

#needed
def counts_table_extractor(ad_obj: ad._core.anndata.AnnData) -> pd.DataFrame:
    if type(ad_obj) == pd.DataFrame:
        out_tab=ad_obj
        return(out_tab) #early exit if it's already a pandas DF
    if type(ad_obj.X) == sc.sparse._csr.csr_matrix:
        out_tab=ad_obj.X.todense()
    else:
        out_tab=ad_obj.X
    colnames=ad_obj.var.iloc[:,0]
    rownames=ad_obj.obs.index
    out_tab=pd.DataFrame(out_tab,columns=colnames,index=rownames)
    return(out_tab)

#needed
def enog_fetcher(enog_list,taxlev):
        l=[]
        for i in enog_list:
                l.append(taxlev in i)
        try:
                ans=enog_list[np.where(l)[0][0]]
        except IndexError:
                ans=None
        if ans:
                ans=ans.split("@"+taxlev+"|")[0]
        return(ans)


#OF number 1
def count_nc_genes(ad_filename: str, ncgenes_filename: str) -> (pd.Series,pd.DataFrame):
    # check if anndata object is in .loom or in .h5ad format:
    check_inputs(ad_filename, ncgenes_filename)

    # import the input data:
    ds_ad,nc_markers=data_importer(ad_filename, ncgenes_filename) #function ensures right datatype output.
    print("Expression data and nervous system markers imported successfully.")
    #Check whether the imported table is an anndata object or not if it is, run 'counts_table_extractor'
    if type(ds_ad) == pd.DataFrame:
        ds_data, nc_markers=correct_axes_orientation(ds_ad,nc_markers)
    elif type(ds_ad) == ad.AnnData:
        # extract the tables - make sure the output table is a dense matrix, turn it into a pandas DF with the correct column and row names.
        ds_data, nc_markers=correct_axes_orientation(counts_table_extractor(ds_ad),nc_markers)
    print("Expression data has the correct axis orientation.")
    nc_gene_ss=ds_data.loc[nc_markers,:]
    gt0=np.sum((nc_gene_ss > 0).astype(int),axis=0)
    
    return(gt0,ds_data)

#OF 2
def top_perc_nc_gene_ids_extractor(nc_markers_representation: pd.Series,ds_data: pd.DataFrame, outfile: str, perc: float = 0.1) -> np.ndarray:
    unique3, counts3 = np.unique(nc_markers_representation,return_counts=True)
    tot_cells=nc_markers_representation.shape[0]
    perc_red=int(np.round(tot_cells*perc))
    inverted_cumsum=np.cumsum(np.flip(counts3))
    smallest_diff_idx=np.where(abs(inverted_cumsum-perc_red) == abs(inverted_cumsum-perc_red).min())[0]
    ten_perc_nc_cell_threshold=np.flip(unique3)[smallest_diff_idx][0]
    nc_occupancy=unique3[np.where(unique3 == ten_perc_nc_cell_threshold)[0][0]:len(unique3)]
    print(f"Top 10% most nerve cell-like cells express {nc_occupancy} of a maximum of {unique3[-1]} nervous system markers expressed.")
    out_arr=np.array([],dtype=str)
    for occup in nc_occupancy:
        new=np.array(nc_markers_representation.index[np.where(nc_markers_representation == occup)[0]])
        out_arr=np.append(out_arr,new)

    ten_perc_nc_cells_df=ds_data.loc[:,out_arr]
    print(f"A total of {len(ten_perc_nc_cells_df.columns)} nerve cells were identified.")
    nc_genes=ten_perc_nc_cells_df.index[np.where(ten_perc_nc_cells_df.sum(axis=1) > 0)[0]]
    outfile=outfile+"_inferred_ncGenes.list"
    np.savetxt(outfile,nc_genes,delimiter=",",fmt='%s')
    print(f"Full list of {nc_genes.size} nervous system genes saved to {outfile}")
    return(nc_genes)

#OF 3
def gene_IDs_to_enog(annot_file: str,gene_ID_file: str or np.ndarray, prefix: str) -> np.ndarray:
        emap_table=pd.read_csv(annot_file,sep="\t")
        if type(gene_ID_file) == str:
            inferred_ncgenes=np.genfromtxt(gene_ID_file,dtype=str)
        else:
            inferred_ncgenes=gene_ID_file    
        in_common=np.intersect1d(inferred_ncgenes,emap_table.loc[:,"#query"])
        emap_table.index=emap_table.loc[:,"#query"]
        query_genes_enogs=np.array(emap_table.loc[in_common,"eggNOG_OGs"])
        out_arr=np.array([])
        for i in query_genes_enogs:
                out_arr=np.append(out_arr,enog_fetcher(i.split(","),"33154"))
        enogs=np.unique(out_arr[out_arr != None])
        outfile=prefix+"_inferred_ncENOGs.list"
        np.savetxt(outfile,enogs,delimiter=",",fmt="%s")
        return(enogs)

def main():
    gt0,ds_data=count_nc_genes(in_data_file,in_ns_markers_file)
    nc_genes=top_tenperc_nc_gene_ids_extractor(gt0,ds_data,prefix)
    nc_genes_enogs=gene_IDs_to_enog(annot_file,nc_genes,prefix)
    #np.savetxt("NC_ENOGs.list",nc_genes_enogs,delimiter=",",fmt="%s")
    return

main()