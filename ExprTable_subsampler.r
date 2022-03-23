library(data.table)

inData<-fread("<filename>.tar.gz") #import the *compressed* cell/gene expression table. Ita can be UMI counts, corrected/normalized counts, etc.
colnames(inData) # Confirm the column names correspond to the barcodes out lists. They should for all except Xenopus
nc_barcodes<-read.table("<filename>",sep=",") # Read in the "[speciescode]_nc_barcodes.csv" tables you produced earlier
colnames(nc_barcodes)<-"barcodes"

# Get only the neural cells from the total data table
nc_column_idcs<-which(colnames(inData) %in% nc_barcodes$barcode) # get the column indices that match with the neural barcodes
NC_cells<-inData[,..nc_column_idcs] # Extract those columns only from the whole scRNA-seq table
# Add the gene names as row names to your table. Usually these will be either as rownames already in "inData", or as the first column name.
gene_names<-inData[,1] # for example, although often it needed a bit more wrangling than that, especially if the gene IDs had a description attached.
colnames(NC_cells)<-gene_names # Because NC_cells has the same number and order of rows as inData, we can assign the row numbers like this.

# Remove genes that are NOT expressed in neural cells (i.e., that their expression value is '0' through all neural cells)
row_tots<-rowSums(NC_cells) # sums row-wise
genes_on_idcs<-which(row_tots > 0) # get row number of genes that have counts/expression above 0.
NC_cells_NC_genes<-NC_cells[genes_on_idcs,]

# Save table to file!
write.table(NC_cells_NC_genes,"<yourfilenamehere>",sep=",", quote=F, col.names=T,row.names=T) # Add any additional formatting flags
# This table will allow you to do two main things:
    # 1: Make the pan neural cell transcriptome by getting the matching sequences from the reference transcriptomes
    # 2. Tabulate on how many cells each gene is expressed.
# Both scripts to follow.