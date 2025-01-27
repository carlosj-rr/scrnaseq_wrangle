import pandas as pd
import numpy as np
import scipy as sc
import scanpy as sp

#Reads in a UMI-counts table in H5AD/LOOM/CSV/TSV/H5 and outputs a pandas DF, making sure the column and index headers are correct.
def umitable_importer(in_table_file, markers):
        # H5AD/LOOM import
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
                #Checks if data is in sparse matrix fmt and converts it to a dense matrix
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
        # Extracts the gene and cell ids, to make sure they are in the right place (for this code, cell IDs are columns and gene ids are indices).
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
           
""" Reads the UMI counts table as a pandas DF, removes uninformative rows by:
simulating a monotone distribution of the same size, based on the *median* gene expression value, and removing the whole row if
its own expression does not differ significantly from the monotone distribution, following a Smirnov-Kolmogorov test.
Then, it recodes all others as
1 - downregulated
2 - normally expressed, and
3 - upregulated,
based on percentiles defined by the 'bound' variable.
Ex. bound = 5: any cell expressing the gene at a level < 5th percentile is recoded to 1
any cell expressing the gene between the 5th percentile and the 95th percentile is recoded to 2
any cell expressing the gene => 95th percentile is recoded to 3

*Outputs the recoded table without the uninformative rows.
"""
def umi_counts_transformer(in_data: pd.DataFrame, bound = 5):
        # Assume a UMI counts table as a pandas Data Frame.
        # define upper and lower percentile bounds
        lower_bound = bound
        upper_bound = 100-bound
        # Make sure all data are in integer format (this is necessary for operations down the line)
        if np.any(in_data.dtypes != int):
                in_data = in_data.astype(int)
        # initialise the list where the indices of informative rows will be.
        informative_rows = []
        for i in range(in_data.shape[0]):
                # check for rows that have no expression data at all.
                nonzero_idcs = in_data.iloc[i].to_numpy().nonzero()
                nonzero_vals = in_data.iloc[i].iloc[nonzero_idcs]
                # Create monotone dataset of the same size as the number of cells, but it's all the median value of the gene's expression
                dummy_uniform_dist = np.repeat(np.median(nonzero_vals),len(nonzero_vals))
                #S-K test to see if observed data are distinguishable from a non-informative monotone dataset
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
                subpop_sample = np.array([ np.random.choice(range(tot_cols), tot_subpop, replace=False) for _ in range(perms)]) # 'subpopulation permutes'
                backg_sample = np.array([list(map((lambda x: np.setdiff1d(np.arange(tot_cols), x)), subpop_sample))])[0] # 'background permutes'
                coords = tuple(map((lambda x,y: (x,y)),subpop_sample, backg_sample)) # tuple where [0] has the random selection that will be the 'subpopulation', and [1] has the background.
                # append the *actual* ratio of subpopulation expression for gene i
                obs_values.append(expression_ratios(bool_data.iloc[i,:],(subpop_indices, background_indices)))
                # Calculate the 5th and 95th percentiles for each permutation
                for_ptest = []
                for j in range(len(coords)):
                        val = expression_ratios(bool_data.iloc[i,:],coords[j])
                        for_ptest.append(val)
                percentiles.append((np.nanpercentile(for_ptest,lower_bound),np.nanpercentile(for_ptest,upper_bound)))
        l = list(zip(obs_values,percentiles))
        dexp_idcs = np.where(list(map((lambda x: np.bitwise_or(x[0] < x[1][0], x[0] > x[1][1])),l)))[0]
        dexp_ids = bool_data.index[dexp_idcs]
        return(dexp_ids)

def expression_ratios(series,coords):
        qty1 = len(coords[0])
        qty2 = len(coords[1])
        ratio_on = sum(series.iloc[coords[0]]/qty1)
        ratio_off = sum(series.iloc[coords[1]]/qty2)
        val = ratio_on - ratio_off
        return (val)
