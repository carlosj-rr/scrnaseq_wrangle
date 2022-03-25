#### __Hydra vulgaris__-specific hacks
In _Hydra_'s gene expression tables, the barcodes are prepended by an X, like so: "X01.D1_ATGCATACACTA" in the gene expression table is written as "01-D1_ATGCATACACTA" in the "Hv_nc_barcodes.csv" file. Notice also that the period (".") in between the first pair of numbers arnd the barcode is swapped for a hyphen ("-") in the lookup table. This can be corrected with bash, but also with R, once the gene expression table was imported.
Ex:
````
x<-"X01.D1_ATGCATACACTA"
> sub("\\.","-",sub("X","",x)) # This nested command first removes the "X", and then swaps the period "." for a hyphen "-". The period character has to be escaped because otherwise it means "any character".
[1] "01-D1_ATGCATACACTA"
````
Ideally though, we would run this command for every single column name, then re-set the column names with the resulting (corrected) vector of column names. This function can be used as a vectorized operation, and we can change all the column names in one full sweep.

````
col_name_cleaner<-function(incolname) {
    corr<-sub("\\.","-",sub("X","",incolname))
    return(corr)
}
to_be_fixed<-colnames(inData)
fixed<-col_name_cleaner(to_be_fixed)
colnames(inData)<-fixed
````
You can run ```colnames(inData)``` to confirm that the results are what you want.

11 mitochondrial genes come up as genes expressed in neural cells, and are not in the reference transcriptome. However, they're expressed in a negligible amount of cells (0.02%-0.33% of all cells), so I'll disregard these. Later on, I will remove the genes expressed in less than 50% of all cells.

````
last_nuc_gene_idx<-which(rownames(NC_cells_NC_genes) == "MT_ATP8")-1 # MT_ATP8 is the first mitochondrial gene of the list, and they're the last 11 row names, so the index before that should be the last gene we get from the reference transcriptome
tromeSubsampleIdcs<-geneIdxFinder(rownames(NC_cells_NC_genes)[1:last_nuc_idx]) # geneIdxFinder()'s call can be modified.
````
