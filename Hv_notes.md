#### _Hydra vulgaris_-specific hacks
In _Hydra_'s gene expression tables, a simple ```fread("Hv_full_expressionTable.csv")``` will import the gene names as the first column, and fill the rownames with numbers. Because 'fread' doesn't implement an option to specify row names from a column's values, we need to do this 'by hand':
````
inData<-fread("Hv_full_expressionTable.csv")
rownames(inData)<-inData$V1
````
You can run ```colnames(inData)``` to confirm that the results are what you want.

When creating the 'outDF' table, the command that produces the 'namevect' has to be modified because the expression data table has a bit of annotation for each gene. The necessary changes are detailed in the "ExprTable_subsampler.R" script.