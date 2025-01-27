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
           

def umi_counts_transformer(in_data: pd.DataFrame, bound = 5):
        """
        FUNCTION TO RECODE UMI COUNTS TABLE BASED ON GLOBAL DIFFERENTIAL EXPRESSION AND REMOVE UNINFORMATIVE GENES
        *INPUT: 
        1. UMI counts table as a pandas DF,
        2. 'bound' variable for differential expression detection stringency
        
        First, this function removes uninformative rows by:
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
        
        *OUTPUT: A recoded table without the uninformative rows.
        """
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
        """ 
        FUNCTION TO SELECT A CELL SUBPOPULATION BASED ON GENETIC MARKERS KNOWN A PRIORI
        INPUTS:
        1. the recoded UMI counts pandas DF (after using 'umi_counts_transformer()' above).
        2. a list of known cell markers (eg. genes differentially expressed in our cell population of interest).
        3. an percentile threshold.
        
        This function FIRST checks whether the input DF has some of the informative markers.
        If they have been removed because they are not informative, the function returns a (sad) message and exits.
        This probably means the dataset is not great.
        
        If there ARE markers present, this function turns the entire (recoded) DF into Boolean values:
        * 'True' for any value which is differentially expressed (ie. having a value of '1' - downregulated, or '3' - upregulated: see lines 57-64 above)
        * 'False' for any value which is NOT differentially expressed (ie. having a '2')
        
        ***NOTE that this operation effectively treats up- and downregulation as equally meaningful signs of differential expression.***
        
        THEN, this function counts, for each cell(=column), how many 'True' values it has *in the cell markers of interest*.
        This returns a distribution ranging from the cells which had the least representation of the cell markers to the ones which had the most.
        Then, the subpopulation marked by the marker of interest is chosen as the cells above the 'percentile' (default 95th) percentile\
        of cells differentially expressing the markers of interest.
        
        Example
        100 markers of interest across 1000 cells:
        Most cells will have 0 of those expressed - but taking pleiotropy in mind, probably some markers will be expressed non-specifically.
        Cells from the population will rarely have all markers expressed (their expression will depend on sub-lineage, cell cycle stage, transcriptional noise, etc.),\
        but a well-identified cell will probably express a high proportion of those 100 markers.
        
        Hence, there will be a (potentially left-biased) distribution that will have on the left tail the cells differentially expressing the least of the markers\
        and on the right tail the ones differentially expressing the most markers of interest.
        
        Only the top 5% (by default) are selected as the putative cell subpopulation.
        
        OUTPUTS:
        1. Boolean table as a pandas DF
        2. *Column indices* of the inferred members of the cell subpopulation as a numpy array
        
        """
        transf = transformed_dataset
        # Check which of the subpop markers are in the 'informative markers' dataset - exit if none
        markers_present = np.intersect1d(transf.index, markers)
        if len(markers_present) == 0:
                print("There are no informative markers in this dataset - :(.")
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
        """
        PERMUTATION TEST-TURBOCHARGED FUNCTION TO IDENTIFY NEW DIFFERENTIALLY EXPRESSED GENES
        INPUTS:
        1. 'Booleanised' UMI dataframe
        2. numpy array of cell subpopulation column indices
        3. number of permutations for permutation test
        4. percentile bound for accepting a new gene as differentially expressed

        The general spirit of this function is to identify genes which are differentially expressed in cells of interest\
        by seeing how much their differential expression 'prefers' the cells of interest.

        This preference is measured by how the sets of cells in which they are differentially expressed includes more or less\
        of the cells of interest in comparison to random choices.

        Ex.
        Gene A is differentially expressed (has a 'True' value) in 100 out of 200 cells.
        We have 100 cells of interest. Is Gene A *preferentially* differentially expressed in our cells of interest?
        Case X: If all 100 cells in which Gene A is differentially expressed are the 100 cells of interest, we would say yes.
        Case Y: On the other hand, if none of 100 cells in which Gene A is differentially expressed are the 100 cells of interest,\
        we would say no. This function measures differential expression preference with a proportional expression preference (PEP) metric.
        This metric is the proportion of cells of interest (foreground cell population) in which a cell is differentially expressed *minus*\
        the proportion of the other cells (background cell population) in which a gene is differentially expressed, and it ranges from -1 to 1,\
        with 0 being no preference whatsoever.
        Case X:
        PEP = 1 - all cells of interest are differentially expressed, and none of the other ones are
        Case Y:
        PEP = -1 - no cell of interest is differentially expressed, AND all of the other ones ARE.

        But what about all the cases in between?

        In order to solve this problem, this function selects 100 random sets of cells *with a size equal to the number of cells of interest*.
        Then, for each gene, it calculates the PEP metric, creating a distribution of PEP values.
        
        """
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
