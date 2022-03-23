library(data.table)

outFilename<-"<Your output filename>"
namevect<-c()
countvect<-c()
percvect<-c()
for (row in 1:nrow(NC_cells_NC_genes)) {
    print(
        paste(
            row,
            "of",
            nrow(
                NC_cells_NC_genes
                )
            )
        )

    gene_ID<-rownames(
        NC_cells_NC_genes[row,]
        )

    count<-length(
            which(
                NC_cells_NC_genes[row,] != 0
                )
            )

    percent<-(
        count/ncol(
            NC_cells_NC_genes
            )
            )*100

    namevect<-c(
        namevect,
        gene_ID
        )

    countvect<-c(
        countvect,
        count
        )

    percvect<-c(
        percvect,
        percent
        )
    }

outDF<-data.frame(
    gene.ID=namevect,
    numCellsOn=countvect,
    percentCellsOn<-percvect
    )

write.table(outDF,
    outFilename,
    quote=F,
    row.names=T
    )