library(data.table)

# Give a name to your output file
outFilename<-"<Your output filename>"

namevect<-rownames(NC_cells_NC_genes)

# Create a boolean version of the table that has TRUE for any value that's nonzero (these tables shouldn't have negative values)
boolData<-NC_cells_NC_genes != 0

# because TRUE is also a 1, we can add the totals of each row, and get the counts of how many TRUE values each row had.
countvect<-rowSums(boolData)

# to get the percentage, just divide by the total number of cells
percvect<-100*(countvect/ncol(NC_cells_NC_genes))

# Create a new dataframe that has the vectors created above as columns.
outDF<-data.frame(
    gene.ID=namevect,
    numCellsOn=countvect,
    percentCellsOn=percvect
)

# Sort data frame by decreasing percent (most common genes at the top)
outDF<-outDF[order(-outDF$percentCellsOn),]

# Save table to disk
write.table(outDF,
    outFilename,
    sep=",",
    quote=F,
    row.names=F
)