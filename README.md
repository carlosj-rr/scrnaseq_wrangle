# README  
This repository has all the functions used in _Molecular characterisation of neurons across animals identify conserved gene expression_ (Rivera-Rivera CJ, Feuda R, Pisani D) to:   

**I.** parse scRNA-seq datasets to detect genes that are differentially expressed in nerve cells - these functions can be generalised for any other cell population. 

**II.** do a permutation test to see if sets of orthogroups between two cell populations differ as expected by chance or not.

The source datasets - both scRNA-seq and CDS - are referenced in the original paper and are publicly available.

***
### DEtector.py and percentile_sorter.py (item I above)

`DEtector.py` is the wrapper for the functions contained in `percentile_sorter.py`, and it uses them in the correct order. 
The code in `percentile_sorter.py` is annotated with the full explanation of each function, and `DEtector.py` doesn't introduce any new functions.

Here is each of percentile_sorter's functions explained (replicated from source code - but source has more detailed explanations). 

#### ```umitable_importer():```
_Imports UMI scRNA-seq data and ensures it is in the format expected by all other functions_

**INPUTS**
 1. `[string]` Filename of scRNA-seq dataset.
 2. `[string]` Filename of a *CSV* file with a list of gene markers for a cell population of interest.

**OUTPUT** 

`[pandas.DataFrame]` The scRNA-seq data table with genes as rows and cells as columns.

#### ```umi_counts_transformer()```
_Recodes UMI counts tables based on global differential expression and removes uninformative genes_

**INPUTS**  
1. `[pandas.DataFrame]` UMI counts table.
2. `[int]` 'bound' variable for differential expression detection stringency

**OUTPUT**

`[pandas.DataFrame]` The recoded table without the uninformative rows.

#### `subpop_cell_selector()` 
_Selects subpopulation of cells based on genetic markers for said subpopulation_ 

**INPUTS** 
1. `[pandas.DataFrame]` the recoded UMI counts table (after using `umi_counts_transformer()` above).
2. `[numpy.array]` a list of known cell markers (eg. genes differentially expressed in our cell population of interest).
3. `[int]` a percentile threshold.

**OUTPUTS** 
1. `[pandas.DataFrame]` 'Booleanised' table
2. `[numpy.array]` **Column indices** of the inferred members of the cell subpopulation

#### `diff_test()`
_The real workhorse: a permutation test-guided function to identify new differentially expressed genes in the UMI counts tables_ 

**INPUTS**
1. `[pandas.DataFrame]` 'Booleanised' UMI dataframe (output #1 from `subpop_cell_selector()`).
2. `[numpy.array]` cell subpopulation column indices (output #2 from `subpop_cell_selector()`).
3. `[int]` number of permutations for permutation test (default is 100).
4. `[int]` percentile bound for accepting a new gene as differentially expressed.

**OUTPUT** 

`[pandas.Series]` Gene IDs of genes showing differential expression in the cell population of interest.


#### `pep()` 
_Calculates the proportional expression preference value_

**INPUTS**
1. `[pandas.Series]` a single gene's 'Booleanised' row (True = differentially expressed, False = expressed normally).
2. `[tuple]` with two `numpy.array` objects: at index `[0]`, one can access the column indices of the 'foreground' cells (i.e. cells of interest), and at index `[1]`, the column indices of the 'background' cells (all the rest).

**OUTPUT** 

`[float]` the input gene's PEP value with the following range: [-1,1]
***
### expression_overlap_functions.py (item II above)
[...]
