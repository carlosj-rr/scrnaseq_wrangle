library(ape)
# Import the reduced table with only the neural cells, and only the genes that were on. It assumes that the gene names are as the row name.
NC_cells_NC_genes<-read.table("<filename>") # Add any necessary flags for 'read.table()'

# Read in the reference transcriptome in FASTA format, as DNAbin object. Sequence names must match completely or in parts the sequence names in the rows of the 'NC_cells_NC_genes' table
Trome<-read.FASTA("<Your reference transcriptome file>")
tromeLabels<-labels(Trome) # Make a vector with only the labels of the transcriptome, so that a simpler data structure can be queried later.

# Define a function that extracts the length of a sequence, given an index number, and a reference transcriptome
seq_len<-function(
    idx,
    Trome
    )
    {
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
getLongestSeqIdx<-function(
    inp_vect
    )
        {
            topsize<-0

            for (
                i in inp_vect
                )
                {

                    len<-seq_len(
                        i,
                        Trome
                        )

                if (
                    len > topsize
                    )
                    {

                        topsize<-len
                        outidx<-i

                    }
                
                }
    
    return(
        outidx
        )

}

# From a list of gene names, get the 
geneIdxFinder<-function(
    idcs_on
    )
    {
        idcs_vect<-c()
            for (
                id in idcs_on
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
                            else
                            {
                                tidx<-getLongestSeqIdx(
                                    tidx
                                    )
                                idcs_vect<-c(
                                    idcs_vect,
                                    tidx
                                    )
                            }
                }
                    idcs_vect<-sort(idcs_vect)
    return(idcs_vect)
}

# Get the transcriptome indices of all the genes in the table. Because the table is a neural cell- and neural gene-only subselection, we can fetch them all to produce a neural transcriptome.
tromeSubsampleIdcs<-geneIdxFinder(
    rownames(
        NC_cells_NC_genes
        )
    )

# Make the neural cell transcriptome by subsampling the reference transcriptome
NC_Trome<-Trome[tromeSubsampleIdcs]

# Save it as a FASTA file
write.FASTA(NC_Trome,"<filename.fasta>")