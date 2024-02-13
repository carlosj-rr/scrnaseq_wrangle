import pandas as pd
import numpy as np
import math
import scipy as sc
import argparse
import scanpy as sp
import os

parser = argparse.ArgumentParser()
parser.add_argument("-sc_data","--scrnaseq_input_dataset", help="AnnData (loom/h5ad) or CSV/TSV/HDF UMI counts table.")
parser.add_argument("-markers_list","--nervous_system_markers", help="CSV file with all nerve cell markers for species.")
parser.add_argument("-o","--output_prefix", help="Prefix to be attached to output files.")
parser.add_argument("-skimp","--skip_import", default=False ,help="Skip the import of the original data file and go straight to the 'transf' table.")
parser.add_argument("-nperms","--number_of_permutations", default=100,help="Number of permutations for permutation test - defaults to 100.")
parser.add_argument("-reg_bound", "--down_upregulation_bound",default=5, help="Percentile under/above which a gene's UMI counts are considered as down/upregulation - defaults to 5.")
parser.add_argument("-spop_ptile", "--subpop_percentile", default=95, help="Cells expressing marker genes above this percentile are selected as cells of interest. Defaults to 95.")
parser.add_argument("-perms_bound", "--permutation_test_bound", default=5, help="Genes differentially expressed in the top/bottom bound of the permutation's distribution are considered enriched in cell subpopulation. Defaults to 5.")

args=parser.parse_args()

in_data_fname = args.scrnaseq_input_dataset
markers_file = args.nervous_system_markers
output_prefix = args.output_prefix
skip_import = args.skip_import
reg_bound = args.down_upregulation_bound
ptile_markers = args.subpop_percentile
permutations = args.number_of_permutations
perms_bound = args.permutation_test_bound

def umitable_importer(in_table_file, markers):
        if ".h5ad" in in_table_file or ".loom" in in_table_file:
                print("input dataset is a compressed AnnData object.")
                if ".h5ad" in in_table_file:
                        print("AnnData object in H5AD format")
                        adata = sp.read_h5ad(in_table_file)
                if ".loom" in in_table_file:
                        print("AnnData object in LOOM format")
                        adata = sp.read_loom(in_table_file)
                axisA = adata.obs.index
                axisB = adata.var.loc[:,adata.var.columns[0]]
                if type(adata.X) == sc.sparse._csr.csr_matrix:
                        print("Dataset is in sparse matrix format - converting to a dense matrix")
                        dat = pd.DataFrame(data=adata.X.to_dense().astype(int))
                else:
                        dat = pd.DataFrame(data=adata.X.astype(int))
        if ".csv" in in_table_file:
                print("Input data is a CSV table.")
                dat = pd.read_csv(in_table_file, sep=",", header=0, index_col = 0)
                axisA = dat.index
                axisB = dat.columns
        if ".tsv" in in_table_file:
                print("Input data is a TSV table.")
                dat = pd.read_csv(in_table_file,sep="\t",header=0, index_col = 0)
                axisA = dat.index
                axisB = dat.columns
        if ".h5" in in_table_file:
                print("Input data is a table in HD5 format.")
                dat = pd.read_hdf(in_table_file, header=0, index_col=0)
                axisA = dat.index
                axisB = dat.columns
        axes = (axisA, axisB)
        
        gene_ids = [ axis for axis in axes if np.intersect1d(axis, markers).size > 0 ]
        cell_ids = [ axis for axis in axes if np.intersect1d(axis, markers).size == 0 ]
        
        if (len(gene_ids), len(cell_ids)) == dat.shape:
                print("Rows and columns match, ensuring correct row- and column-names.")
                dat.index = gene_ids
                dat.columns = cell_ids
        elif (len(cell_ids), len(gene_ids)) == dat.shape:
                print("Rows and columns don't match - assigning row- and column-names, and transposing table.")
                dat.index = cell_ids
                dat.columns = gene_ids
                dat = dat.T
        return(dat)      
           

def umi_counts_transformer(in_data: pd.DataFrame, bound = 5):
        # Assume a UMI counts table as a pandas Data Frame.
        # check if entire dataset is of integers (even if they're formatted as float) - and make sure they're encoded as int
        lower_bound = bound
        upper_bound = 100-bound
        if np.any(in_data.dtypes != int):
                in_data = in_data.astype(int)
        informative_rows = []
        for i in range(in_data.shape[0]):
                nonzero_idcs = in_data.iloc[i].to_numpy().nonzero()
                nonzero_vals = in_data.iloc[i].iloc[nonzero_idcs]
                # Test if counts are uninformative - if the distribution of counts resembles a uniform distribution, they're not informative.
                dummy_uniform_dist = np.repeat(np.median(nonzero_vals),len(nonzero_vals))
                if len(nonzero_vals) == 0 or sc.stats.kstest(nonzero_vals, dummy_uniform_dist).pvalue > 0.05:
                        print(f"Row {i} is not informative...skipping it")
                else:
                        print(f"Row {i} *is* informative!, recoding values.")
                        informative_rows.append(i)
                        percentile_5 = np.percentile(nonzero_vals,lower_bound)
                        percentile_95 = np.percentile(nonzero_vals,upper_bound)
                        lower_5 = np.where(in_data.iloc[i] < percentile_5)[0]
                        # Set lower 5% to 1
                        in_data.iloc[i,lower_5] = 1
                        upper_5 = np.where(in_data.iloc[i] >= percentile_95)[0]
                        # Set top 5% to 3
                        in_data.iloc[i,upper_5] = 3
                        extremes = np.hstack((lower_5,upper_5))
                        normies = np.setdiff1d(nonzero_idcs,extremes)
                        # Set all other non-zero values to 2
                        in_data.iloc[i,normies] = 2
                out_data=in_data.iloc[informative_rows]
        return out_data
    
def subpop_cell_selector(transformed_dataset: pd.DataFrame, markers: np.ndarray, percentile = 95):
        transf = transformed_dataset
        # Check which of the subpop markers are in the 'informative markers' dataset - exit if none
        markers_present = np.intersect1d(transf.index, markers)
        if len(markers_present) == 0:
                print("No markers are informative in this dataset - :(.")
                exit
        else:
                # Get the row index of each of the subpopulation markers, and sort them
                subpop_markers=list(map((lambda x: np.where(x == transf.index)[0][0]),markers_present))
                subpop_markers.sort()
                # transform the dataset into a boolean data frame where 'True' is values which were in 1 (downreg) OR 3 (upreg), and otherwise False
                bool_dset = np.bitwise_or((transf == 1),(transf == 3))
                # calculate the number of times markers appear on each cell
                subpop_marker_count = np.sum(bool_dset.iloc[:,subpop_markers],axis=0)
                # get the column index for all cells which are on the 95th percentile of expression of markers
                subpop_indices = np.where(subpop_marker_count >= np.percentile(subpop_marker_count,percentile))[0]
        return(bool_dset,subpop_indices)

def diff_test(bool_data: pd.DataFrame, subpop_indices, perms = 100, bound = 5):
        lower_bound = bound
        upper_bound = 100-bound
        # calculate the total of columns in the full dataset
        tot_cols = bool_data.shape[1]
        # identify the indices of the 'background' columns (those not in the subpop of interest)
        background_indices = np.setdiff1d(range(tot_cols),subpop_indices)
        # count the total of subpopulation cells for calculating proportions later
        tot_subpop = subpop_indices.size
        # count the total of background cells for calculating proportions later
        tot_background = background_indices.size
        # initialize the lists where some of the data will end up in
        percentiles = []
        obs_values = []
        for i in range(bool_data.shape[0]):
                # set up random selections of cells for permutation test - defaults to 100 permutations
                subpop_sample = np.random.choice(range(tot_cols),size=(perms,tot_subpop),replace=False) # 'subpopulation permutes'
                backg_sample = np.array(list(map((lambda x: np.setdiff1d(np.arange(tot_cols), x)), subpop_sample))) # 'background permutes'
                coords = tuple(map((lambda x,y: (x,y)),subpop_sample, backg_sample)) # tuple where [0] has the random selection that will be the 'subpopulation', and [1] has the background.
                # append the *actual* ratio of subpopulation expression for gene i
                obs_values.append((sum(bool_data.iloc[i,subpop_indices])/tot_subpop)/(sum(bool_data.iloc[i,background_indices])/tot_background))
                # Calculate the 5th and 95th percentiles for each permutation
                for_ptest = []
                for j in range(len(coords)):
                        val = expression_ratios(bool_data.iloc[i,:],coords[j])
                        for_ptest.append(val)
                percentiles.append((np.percentile(for_ptest,lower_bound),np.percentile(for_ptest,upper_bound)))
        l = list(zip(obs_values,percentiles))
        dexp_idcs = np.where(list(map((lambda x: np.bitwise_or(x[0] < x[1][0], x[0] > x[1][1])),l)))[0]
        dexp_ids = bool_data.index[dexp_idcs]
        return(dexp_ids)

def expression_ratios(series,coords):
        qty1 = len(coords[0])
        qty2 = len(coords[1])
        return (sum(series.iloc[coords[0]]/qty1))/(sum(series.iloc[coords[1]]/qty2))

if not skip_import:
        print("Beginning from a raw CSV UMI counts table...")
        markers = np.genfromtxt(markers_file,delimiter=",",dtype=str)
        in_data = umitable_importer(in_data_fname, markers)
        # Recode dataset into 0 (not expressed), 1 (downregulated), 2 (within 90% of all values), and 3 (upregulated), and remove uninformative genes
        transf = umi_counts_transformer(in_data, reg_bound)
        print(f"Saving transf data as an intermediate step...the transf table has {len(transf.index)} genes.")
        transf.to_csv(output_prefix+"_transf_table.csv",sep=",")
else:
        transf_file = [file for file in os.listdir() if "transf_table.csv" in file][0]
        print(f"Skipped the import of initial table, beginning from the 'transf' table {transf_file}.")
        transf = pd.read_csv(transf_file,sep=",",header=0,index_col=0)
        markers = np.genfromtxt(markers_file,delimiter=",",dtype=str)

        # Select NS cells, and booleanize the dataset into False (0/2) and True (1/3)
transf, subpop_indices = subpop_cell_selector(transf, markers, ptile_markers)

        #FINALLY, output the list of *other* genes differentially expressed in the subpopulation of cells (uses permutation test)
diff_exp_genes = diff_test(transf, subpop_indices, permutations, perms_bound)

        # Save the list to the disk
np.savetxt(output_prefix+"_inferred_neural_markers.csv",diff_exp_genes,fmt="%s",delimiter=",")