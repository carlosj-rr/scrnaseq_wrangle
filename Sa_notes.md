### _Schistosoma mansoni_-specic hacks
Schman's sequence IDs in the reference transcriptome have a space in between, so, in the geneIdxFinder() function, a change has to be done.

In the line where the 'query' variable is set, swap the original ```query<-paste(id,"$",sep="")```, which searches for a sequence name, and then an End Of Line (EOL) character, switch to ```query<-paste(id," ",sep="")```, which matches to the sequence ID, then a _space_. Adding the character after the sequence name avoids picking sequences whose id's are subsets of other sequences. For example, querying "gene123" will pick up "gene123", but also "gene123a", "gene123a-like", "gene123b", "gene123c", etc., but this is avoided if we query for the name _and_ the character right after. It's usually an EOL, a space, or a pipe symbol "|".

The gene expression table of Schman includes 8 sequences used to standardize read counting, from the External RNA Controls Consortium (ERCC). These are IDs in between ERCC-[00000-99999]. These entries can be removed, as they do not represent genetic information from the sample, and thus will not be present in the reference transcriptome.

````
# Beginning with 'outDF' as the table of NC genes expressed in 50% or more cells
# Get the IDs of the offending sequences
erccs<-outDF$gene.ID[grep("ERCC",outDF$gene.ID)]

# Only keep the rows in which the offending sequence names are NOT present
outDF<-outDF[which(!outDF$gene.ID %in% erccs),]
````
Voila! The outDF table now only has gene IDs from the actual organisms. Save to disk with ```write.table(outDF,"Sa_geneCounts.csv",sep=",",col.names=T,row.names=F,quote=F)```.