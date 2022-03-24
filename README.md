# scrnaseq_wrangle
The following guide explains step-by-step how to obtain, from the source datasets in [Tarashansky et al. (2021)](https://elifesciences.org/articles/66747), the results from our work in [REF AND LINK].

### Dependencies
For most R scripts, the **APE package** is needed. Versions 3.0 and under will NOT have support for the ```write.FASTA()``` function, so ideally install a version above that, for the following scripts we used _v5.6-1_.
The R package **"data.table"** (_v1.14.2_ for the commands shown here) is also needed in order to read gzipped tables without the need of decompressing the file before.
For mouse, the Bioconductor package **"MouseGastrulationData"** (_v1.4.0_ for the commands shown here) is needed, this project has both a [GitHub repository](https://github.com/MarioniLab/EmbryoTimecourse2018), and a [Bioconductor page](https://bioconductor.org/packages/release/data/experiment/html/MouseGastrulationData.html). Once installed, follow the commands detailed in the next section. We used the data from the Bioconductor package.
***
## Important note
The commands detailed below will work directly for the majority of the datasets, but the nature of some of the original datasets (listed in the reference above, within their "Data availability" statement) required some additional work. The main particular cases are:
- _Danio rerio_: The cells identified as nervous system cells span normalized expression tables from four different time points: 10hpf, 14hpf, 18hpf, and 24hpf. In addition, the reference transcriptome was made available to us via the author of the original dataset, and we would be willing to share it upon request.
- _Xenopus tropicalis_: The column headers of this dataset had more information than the rest, so it had to be pre-porcessed in bash with ```zcat GSE113074_Corrected_combined.annotated_counts.tsv.gz | tail -n+10 | tr \[:blank:] "," > GSE113074_Corr_littleHead.csv```, which extracts the top 9 lines, and translates all the tab characters into commas.
- _Mus musculus_: The raw expression tables for mouse are in a Bioconductor package (see installation instructions above in **Dependencies**), and in order to get it, install and load the package in R (MouseGastrulationData), then run: ```> inData<-MouseGastrulationData::EmbryoAtlasData()``` **RAM will be essential for this command, otherwise the large size of the table can induce a crash**.
***
### Getting the raw data
1. Download the Supplementary File 1 from the reference above, called "elife_66747-suppl1-v2.xslx", and titled ["Cell Atlas metadata and cell annotations"](https://cdn.elifesciences.org/articles/66747/elife-66747-supp1-v2.xlsx). This Excel table has a sheet for each species, and shows to which cluster each cell (by barcode) was assigned.

* Go through each individual sheet by hand, and export each as a CSV file, for the example here, use short species names, followed by "\_clusters.csv". Ex. _Hydra vulgaris_:    Hv_clusters.csv, _Spongilla lacustris_: Sl_clusters.csv, etc. (feel free to check example files in this repo, if you want to be sure it's all looking the same). Make sure the shortened names match the following list:
- _Hydra vulgaris_ - **Hv**
- _Spongilla lacustris_ - **Sl**
- _Schistosoma mansoni_ - **Sa**
- _Schmidtea mediterranea_ - **Sm**
- _Danio rerio_ - **Dr**
- _Mus musculus_ - **Mm**
- _Xenopus tropicalis_ - **Xt**
2. Download the Supplementary File 7 as well, called "elife_66747-suppl7-v2.xslx", titled ["Cell types in the cell type families shown in Figure 5C"](https://cdn.elifesciences.org/articles/66747/elife-66747-supp7-v2.xlsx). This Excel table identifies which clusters belonged to neural cells (NCs), and which belonged to muscle cells, according to their SAMap analysis.
  a. Also export as a CSV file, named "Neural_and_Muscle_clusterNames.csv".
***
### Parsing the raw data
1. From the _"Neural_and_Muscle_clusterNames.csv"_, extract, for each species, the cluster _name_ of all the NC clusters. The script used is called ```cluster_id_getter.sh```, which outputs a CSV file called __"spp_clusterName.csv"__ in which the first column has the 2-letter species name, and the second column has the NC cluster names. Make sure the 2-letter species names of the "\_clusters.csv" files produced on step 1 above match the 2-letter codes that appear in the "spp_clusterName.csv" file. If you followed the codes in point 1 above, you should be fine.
2. Now, we can use the "spp_clusterName.csv" to subsample each of the "\_clusters.csv" files, choosing only the lines that contain the NC names in it, which will also give us the cell barcodes and some other info we will have to get rid of later on. Use ```neural_barcode_getter.sh```.
   * This script will produce a file called "[sppcode]\_nc\_barcodes.csv.tmp" for each species, containing only the rows that matched the NC cluster names from the "Neural_and_Muscle_clusterNames.csv" file.
   * These files have to be checked by eye, so superfluous information can be removed. Once this information has been removed, the files can be renamed, excluding the ".tmp" extension. The exact commands we used are in ```species_specific_NBarcodeExtractor.sh```. Notice that for all species but _Xenopus tropicalis_, it's the first column that has the barcode information. In _Xenopus_, the first three columns have information regarding the barcode, so we're keeping it all just in case.
### Processing each species' gene expression tables
[Repeat for each species]
 1. ```ExprTable_subsampler.R```: Reduce a complete scRNA-seq UMI counts/normalized/TMP/expression table to only the NCs, and only the genes which are ON in these cells. Assumes noise has been removed from each table.
 2. ```geneCountsTabulato.R```: Count on how many cells each gene is expressed, and calculate the percent of the total cells that represents.
 3. ```NC_Trome_make.R```: Using the reference transcriptomes from each study, extract the full length sequence of the canoical (=longest) RNA transcript of all the genes expressed in the NCs.
***
The next steps are detailed in the Methods section of our [paper](empty link).