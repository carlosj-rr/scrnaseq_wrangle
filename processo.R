library(data.table)

#import the *compressed* cell/gene expression table. Ita can be UMI counts, corrected/normalized counts, etc.
inData<-fread(
    "<filename>.tar.gz"
    ) 

# Confirm the column names correspond to the barcodes out lists. They should, but otherwise make it so (there are species-specific notes where some modifications were needed).
colnames(
    inData
    )

# Confirm the row names correspond to the gene IDs. They should, but some species need a bit of formatting here, ****see Spolac's species-specific notes****
rownames(
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

# Extract those columns only from the whole scRNA-seq table.
NC_cells<-inData[,..nc_column_idcs]

# Re-set the column names of NC_cells.
rownames(NC_cells)<-rownames(inData)

### Remove genes that are NOT expressed in neural cells (i.e., that their expression value is '0' through all neural cells)
row_tots<-rowSums(
    NC_cells
) # sums row-wise

genes_on_idcs<-which(
    row_tots > 0
) # get row number of genes that have counts/expression above 0.

# Subsample the rows, and re-set the rownames again
NC_cells_NC_genes<-NC_cells[genes_on_idcs,]
rownames(NC_cells_NC_genes)<-rownames(NC_cells)[genes_on_idcs]

# Save table to file! (confirm the following things before: 1. row names are gene IDs, 2. column names are cell barcodes, and 3. output file is CSV)
write.table(NC_cells_NC_genes,
    "<yourfilenamehere>",
    sep=",",
    quote=F,
    col.names=T,
    row.names=T
) # Add any additional formatting flags
# This table will allow you to do two main things:
    # 1: Make the pan-neural cell transcriptome by getting the matching sequences from the reference transcriptomes
    # 2. Tabulate on how many cells each gene is expressed

# Tabulate the on how many cells each gene is expressed:

namevect<-rownames(NC_cells_NC_genes) # Hydra's command is different, see species-specific notes for Hydra

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

# Sort data frame by decreasing percent (most common genes at the top), and select only the genes which are expressed in 50% of cells or more
outDF<-outDF[order(-outDF$percentCellsOn),]
outDF<-outDF[outDF$percentCellsOn >= 50,]

# Save table to disk
write.table(outDF,
    "outFile.csv",
    sep=",",
    quote=F,
    row.names=F
)

# Saving neural cell-only transcriptome

library(ape)

Trome<-read.FASTA("<Your reference transcriptome file>")
tromeLabels<-labels(Trome) # Make a vector with only the labels of the transcriptome, so that a simpler data structure can be queried later.

# Define a function that extracts the length of a sequence, given an index number, and a reference transcriptome
seq_len<-function(idx,Trome) {
    size<-as.integer(
        summary(
            Trome[idx]
                )[1]
            )
    return(
        size
        )
}

# Define a function that uses 'seq_len()' to get, from a vector of sequence indices in a transcriptome, the longest sequence. This is useful for cases in which the same gene ID has several transcripts, and you want only the longest - canonical one.
getLongestSeqIdx<-function(inp_vect) {
            topsize<-0
            for (i in inp_vect) {
                    len<-seq_len(
                        i,
                        Trome
                        )
                if (len > topsize) {
                        topsize<-len
                        outidx<-i
                    }
                }
    
    return(outidx)
}


# From a list of gene names, get the index of the longest transcript present in the reference transcriptome. Also output a table that lists all the IDs that were not found.
geneIdxFinder<-function(geneids) {
        idcs_vect<-c()
        orphans<-c()
            for (
                id in geneids
                )
                {
                    print(paste("Getting the sequence index of",id))
                    # 'grep' can be tricky sometimes, because many names can be subsets of other names, this is why a query input must be built with paste.
                    # In this case, it will look for the sequence ID, and en EOL character afterwards.
                    query<-paste(id,"$",sep="")
                    tidx<-grep(query,tromeLabels)
        
                    if (
                        length(
                            tidx) == 1 
                            )
                            {
                                idcs_vect<-c(
                                    idcs_vect,
                                    tidx
                                    )
                            }
                            else if (length(tidx) > 1)
                            {
                                print(tidx)
                                tidx<-getLongestSeqIdx(
                                    tidx
                                    )
                                idcs_vect<-c(
                                    idcs_vect,
                                    tidx
                                    )
                            } else if (length(tidx) == 0) {
                                print(paste(id,"not found"))
                                orphans<-c(orphans,id)
                            }
                }
                    idcs_vect<-sort(idcs_vect)
    return(list(idcs_vect,orphans))
}

# Get the transcriptome indices of all the genes in the table. Because the table is a neural cell- and neural gene-only subselection, we can fetch them all to produce a neural transcriptome.
x<-geneIdxFinder(
    outDF$gene.ID
    )
tromeSubsampleIdcs<-x[[1]]
orphans<-x[[2]]

# Make the neural cell transcriptome by subsampling the reference transcriptome, store orphan hits, if any, and if so, also remove them from the 'outDF' table.
if (length(orphans) == 0 ){
    NC_Trome<-Trome[tromeSubsampleIdcs]
} else {
    NC_Trome<-Trome[tromeSubsampleIdcs]
    print("There is a nonzero amount of orphans, please see list in the file <orphans.csv>")
    write.table(orphans,"orphans.csv",sep=",",quote=F,col.names=F,row.names=F)
    print("In addition, these IDs will be removed from the <outDF> table with the most common genes in NC cells")
    outDF<-outDF[which(!outDF$gene.ID %in% orphans),]
    print("outDF table of NC gene counts has been updated - I kindly suggest that the table is re-saved to disk (see 'write.table()' command above)")
}

# Save it as a FASTA file
write.FASTA(
    NC_Trome,
    "<filename.fas>"
)
