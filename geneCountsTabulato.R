library(data.table)

# Give a name to your output file
outFilename<-"<Your output filename>"

# Initialize empty vectors for what will become the 3 columns of the output table:
namevect<-c() # Gene ID
countvect<-c() # in how many cells said gene was ON
percvect<-c() # in what percent of cells said gene was ON

# Iterate through the rows (=genes)...
for (row in 1:nrow(NC_cells_NC_genes)) {
    # Inform what's happening - also useful for debugging
    print(
        paste(
            row,
            "of",
            nrow(
                NC_cells_NC_genes
                )
            )
        )

    # Extract the gene's ID: ASSUMES THE GENE ID IS THE ROWNAME - WILL NOT WORK IF THE GENE ID IS A COLUMN
    gene_ID<-rownames(
        NC_cells_NC_genes[row,]
        )
    # Count how often the gene is expressed in its row (how often the value in the row is nonzero).
    count<-length(
            which(
                NC_cells_NC_genes[row,] != 0
                )
            )
    # Calculate the percent that represents.
    percent<-(
        count/ncol(
            NC_cells_NC_genes
            )
            )*100

    # Add the current gene ID to the gene ID vector.
    namevect<-c(
        namevect,
        gene_ID
        )

    # Add the current counts to the counts vector.
    countvect<-c(
        countvect,
        count
        )

    # Add the current percent to the percents vector.
    percvect<-c(
        percvect,
        percent
        )
} #Iteration finishes

# Create a new dataframe that has the vectors created in the loop above as the columns
outDF<-data.frame(
    gene.ID=namevect,
    numCellsOn=countvect,
    percentCellsOn<-percvect
    )

# Save table to disk
write.table(outDF,
    outFilename,
    quote=F,
    row.names=T
    )