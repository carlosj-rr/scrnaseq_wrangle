library(data.table)

#import the *compressed* cell/gene expression table. Ita can be UMI counts, corrected/normalized counts, etc.
inData<-fread(
    "<filename>.tar.gz"
    ) 

# Confirm the column names correspond to the barcodes out lists. They should, but otherwise make it so (there are species-specific notes where some modifications were needed).
colnames(
    inData
    ) 

# Read in the "[speciescode]_nc_barcodes.csv" tables you produced earlier
nc_barcodes<-read.table(
    "<filename>",
    sep=","
    ) 

colnames(
    nc_barcodes
    )<-"barcodes"

# Get only the neural cells from the total data table
nc_column_idcs<-which(
    # get the column indices that match with the neural barcodes
    colnames(inData) %in% nc_barcodes$barcode
    ) 

# Extract those columns only from the whole scRNA-seq table
NC_cells<-inData[,..nc_column_idcs] 

# OPTIONAL: Add the gene names as row names to your table. Usually these will be either as rownames already in "inData", or as the first column name.
# for example, although often it needed a bit more wrangling than that, especially if the gene IDs had a description attached.
gene_names<-inData[,1] 
rownames(
    # Because NC_cells has the same number and order of rows as inData, we can assign the row numbers like this.
    NC_cells
)<-gene_names 

# FOR THE SPONGE - SEE BELOW
# IF the gene name has a description attached, it will probably have the ID, then a space, then a description. In these cases, do:
short_names<-c() #inintialize a vector to store all the ids withOUT their description
for (idx in 1:length(gene_names))  {
        # Use strsplit to separate the gene names from the description by splitting the string at each space (" "). If the character that separates them is different, modify the "split" parameter. Access the splitted string like a list
        short_bit<-strsplit(
            gene_names[idx],
            split=" "
        )[[1]][1]

        # Add the shortened name into the vecor initialized before
        short_names<-c(
            short_names,
            short_bit
        )
}
# That's it! Now 'short_names' has all the gene IDs, withOUT their description. Now, we can make it the rownames of our NC_cells table:
rownames(
    NC_cells
)<-short_names

### Remove genes that are NOT expressed in neural cells (i.e., that their expression value is '0' through all neural cells)
row_tots<-rowSums(
    NC_cells
) # sums row-wise

genes_on_idcs<-which(
    row_tots > 0
) # get row number of genes that have counts/expression above 0.

NC_cells_NC_genes<-NC_cells[genes_on_idcs,]

# Save table to file!
write.table(NC_cells_NC_genes,
    "<yourfilenamehere>",
    sep=",",
    quote=F,
    col.names=T,
    row.names=T
) # Add any additional formatting flags
# This table will allow you to do two main things:
    # 1: Make the pan neural cell transcriptome by getting the matching sequences from the reference transcriptomes
    # 2. Tabulate on how many cells each gene is expressed - "geneCountsTabulato.R".