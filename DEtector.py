import percentile_sorter as ps
import numpy as np
import pandas as pd
import argparse
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


if not skip_import:
        print("Beginning from a raw CSV UMI counts table...")
        markers = np.genfromtxt(markers_file,delimiter=",",dtype=str)
        in_data = ps.umitable_importer(in_data_fname, markers)
        # Recode dataset into 0 (not expressed), 1 (downregulated), 2 (within 90% of all values), and 3 (upregulated), and remove uninformative genes
        transf = ps.umi_counts_transformer(in_data, reg_bound)
        print(f"Saving transf data as an intermediate step...the transf table has {len(transf.index)} genes.")
        transf.to_csv(output_prefix+"_transf_table.csv",sep=",")
else:
        transf_file = [file for file in os.listdir() if "transf_table.csv" in file][0]
        print(f"Skipped the import of initial table, beginning from the 'transf' table {transf_file}.")
        transf = pd.read_csv(transf_file,sep=",",header=0,index_col=0)
        markers = np.genfromtxt(markers_file,delimiter=",",dtype=str)

        # Select NS cells, and booleanize the dataset into False (0/2) and True (1/3)
transf, subpop_indices = ps.subpop_cell_selector(transf, markers, ptile_markers)

        #FINALLY, output the list of *other* genes differentially expressed in the subpopulation of cells (uses permutation test)
diff_exp_genes = ps.diff_test(transf, subpop_indices, permutations, perms_bound)

        # Save the list to the disk
np.savetxt(output_prefix+"_inferred_markers.csv",diff_exp_genes,fmt="%s",delimiter=",")