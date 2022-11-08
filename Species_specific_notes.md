## Species-specific hacks (where necessary) for reading scRNA-seq data and obtaining the reference transcriptomes, in alphabetical order:
***
### _Danio rerio_-specific hacks
***
Dr's full expression gene dataset has to be extracted from four different timepoints, 10, 14, 18, and 24hpf. I inspected the files and all of the NC barcodes are included in those four tables. Note that these are the wild-type embryos, which don't have any annotation in the study GSE112294 (i.e. GSM3067192-5).

First, fetch all four datasets, and rename if necessary. For the purposes of this guide, I will keep their original names (e.g. GSM306719#_##hpf_nm.csv.gz)
Then, in __R__:
````
# Read in each of the datasets, and join them all into one large table, and set the rownames using the "Row" column header. Don't worry, I've made sure that the rows are in the same order.
inData<-fread("GSM3067192_10hpf_nm.csv.gz")
inData<-cbind(inData,fread("GSM3067193_14hpf_nm.csv.gz")))
inData<-cbind(inData,fread("GSM3067194_18hpf_nm.csv.gz"))
inData<-cbind(inData,fread("GSM3067195_24hpf_nm.csv.gz"))
rownames(inData)<-inData$Row
````
***
### _Hydra vulgaris_-specific hacks
***
In _Hydra_'s gene expression tables, a simple ```fread("Hv_full_expressionTable.csv")``` will import the gene names as the first column, and fill the rownames with numbers. Because 'fread' doesn't implement an option to specify row names from a column's values, we need to do this 'by hand':
````
inData<-fread("Hv_full_expressionTable.csv")
rownames(inData)<-inData$V1
````
You can run ```rownames(inData)``` to confirm that the results are what you want.

When creating the 'outDF' table, the command that produces the 'namevect' has to be modified because the expression data table has a bit of annotation for each gene:

````
splitter<-function(longname) {
    shortname<-strsplit(longname,split="\\|")[[1]][1];
    return(shortname)
    }

outDF$gene.ID<-as.vector(sapply(outDF$gene.ID,splitter))
namevect<-outDF$gene.ID
````
***
### _Mus musculus_-specific hacks
***
The mouse datasets have to be pre-processed a bit before reaching the point in which we can get to the R scripts.

First, as mentioned in the README, the data has to be imported from a package in Bioconductor, and this is a process of 2 steps (assuming you already imported the library):
````
inStudy<-MouseGastrulationData::EmbryoAtlasData()
inData<-assay(inStudy,"counts")
````
The dataset imported is a data.frame, not a data.table, so the subsampling of the NC columns must be modified from the 'data.table' 'double period' notation: ```NC_cells<-inData[,..nc_column_idcs]```, into the regular subsampling notation of data.frames: ```NC_cells<-inData[,nc_column_idcs]```.

For the reference transcriptome, I extracted the list of gene ids expressed in 50% or more cells from the "Mm_nc_expressionTable.csv" file, and used it as a BioMart query: ```tail -n+2 Mm_nc_expressionTable.csv | cut -d"," -f1 > BM_query.list # From line 2 onwards to avoid column headers.``` This file can now be given to [BioMart](https://www.ensembl.org/biomart/martview), which will fetch all the hits in the current database (GRCm39 as of March 2022), including all the transcripts. Request the cDNA sequences, and only the Gene ID in the FASTA header (not the version, not the transcript ID - this is important later on).

Once you've imported the BioMart results using ```Trome<-read.FASTA("filename")```, you can proceed normally with the analysis. You will notice that there are **6 'orphans'** - IDs with no hits on the references. These are outdated gene IDs that were eliminated when going from GRCm38 to GRCm39, so we can forget about them for the moment. This also means, however, that we have to remove them from the table of common neural genes ("outDF"). I've included an if block that does this, and explains that the object in the session was modified, suggests the older version to be overwritten, and mentions that the list of orphans was saved into a file called "orphans.csv". This block only activates if there were gene IDs that were not found in the reference transcriptome.
***
### _Schistosoma mansoni_-specic hacks
***
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
***
### _Schmidtea mediterranea_-specific notes
***
First, it is not clear from the databases in the [project's SRA site](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE111764), which of the tables is the correct one. My main worry was that whichever dataset I used should have all of the cell barcodes that I want to find, so I got hold of all the matrices, and checked their column names to see if one of them had all the NC barcodes I wanted, or whether I had to combine datasets. I was lucky, and there's one that has them all: __"GSE111764_PrincipalClusteringDigitalExpressionMatrix.dge.txt.gz"__.

Then, the latest reference transcriptome available to date (Spring 2022) is the 'v6' one, so all transcript names have the following nomenclature: "dd_Smed_v6_[Numbers]", however, the gene expression table mapped reads onto an earier version, 'v4', so all of _those_ gene IDs are like so: "dd_Smed_v4_[Numbers]". For this work, I "brute forced" the IDs into the same names by changing the v6 into v4 in the transcriptome (the transcriptome is a smaller file, so changing the names there is faster than in the gene expression table). I assumed here that the numbers remained the same (i.e. that "dd_Smed_v4_1432_0_5" and "dd_Smed_v6_1432_0_5" would be used for the same sequence, on different database versions - a quick cross-check on the different planarian databases confirmed this).
In __bash__:
```sed -i 's/_v6_/_v4_/g' dd_Smed_v6.pcf.contigs.fasta```


In order to keep things relatively consistent, I did not change the file name, and kept the 'v6' there.
***
### _Spongilla lacustris_-specific notes
***
In Spongilla's table of gene expression, each gene has a description attached after the Gene ID, these have to be removed to the querying of the gene IDs in the transcriptome can proceed normally. Once the gene expression table is 

````
gene_names<-inData$V1
short_names<-c() #initialize a vector to store all the ids withOUT their description
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

rownames(
    inData
)<-short_names
````

When extracting the NC transcripts, the end character of the "query" assignment has to be removed. I.e. change ```query<-paste(id,"$",sep="")``` OR ```query<-paste(id," ",sep="")``` for ```query<-paste(id,sep="")```. This change could also be simply ```query<-id``` as well, but keeping the paste command helps to keep the code more or less stable for the analysis of the other species.
***
### _Xenopus tropicalis_-specific notes
***
First off, the list of NC cell barcodes from _Xenopus_ has duplicates, and has to be deduplicated. This is odd, because the first field we extracted, the one that starts with a C, and then a number of 5 OR 6 digits _is_ unique to each row of the table, but the barcodes are not, as are not the "bc[A-Z]{3}" ids. The "C numbers" are not in the gene expression tables, so we don't have an option but to use the "bc codes" and the actual barcode sequences to select the columns, which have to be deduplicated. In one command, I'll deduplicate the list, and create what will be the proper column IDs, which will combine the "bc code" and the barcode sequence with an underscore. Also, R does not permit hyphens "-" in its column names, and automatically changes them to periods ".", so we have to make that change as well so the barcodes match up. For example "C45874,bc0005,TTCGGCCT-AGCGAAGT" will become "bc0005_TTCGGCCT.AGCGAAGT". To do this, run in __bash__:
````
mv Xt_nc_barcodes.csv Xt_nc_barcodes.csv.tmp
cut -d"," -f2,3 Xt_nc_barcodes.csv.tmp | sed 's/,/_/g' | sed 's/-/./g' | sort | uniq > Xt_nc_barcodes.csv
````

This also means that we have to pre-process the gene expression table, so that the sequence headers match the barcodes we constructed for "Xt_nc_barcodes.csv".
1. The actual expression data begin in line 10 of the "GSE113074_Corrected_combined.annotated_counts.tsv.gz" file. Almost everything before (except lines 5 and 6, see below) is information that we cannot use for the moment. We would like to extract line 5 (barcode name, the "bc code") and line 6 (barcode sequence), and join them with an underscore ("_"), to create a new vector of column names. We can do this in two steps, one in bash, and one in R:
__In bash__:
````
zcat GSE113074_Corrected_combined.annotated_counts.tsv.gz | head -n5 | tail -n1 | sed "s/\t/\n/g" > bcode_names.list
zcat GSE113074_Corrected_combined.annotated_counts.tsv.gz | head -n6 | tail -n1 | sed "s/\t/\n/g" > bcode_seqs.list
````
__In R__:
````
library(data.table)
bnames<-as.vector(read.csv("bcode_names.list",header=T)$Barcode_name) # get the barcode names, and import it as a vector
bseqs<-as.vector(read.csv("bcode_seqs.list",header=T)$Barcode_sequence) # import the barcode sequences as a vector
join<-function(a,b) { result<-paste(a,"_",b,sep=""); return(result) } # create a function that joins two strings with an underscore
new_header<-join(bnames,bseqs) # use function in a vectorized manner with both character vectors as input. This is the new header
manyheaders<-fread("GSE113074_Corrected_combined.annotated_counts.tsv.gz") # Read in the whole table, with all the extra headers
headless<-manyheaders[10:nrow(manyheaders),2:ncol(manyheaders)] # extract only the expression data (the rows removed have the extra headers, and the column removed has the gene ids)
colnames(headless)<-new_header # add the column names we constructed earlier
inData<-as.data.frame(lapply(headless,as.numeric)) # set "inData", and switch to numeric (the import is done as characters)
rownames(inData)<-manyheaders$InDrops_version[10:nrow(manyheaders)] # rename the rows so they have the gene IDs


````
When subsampling for the NC cell columns, keep in mind that inData is now a data.frame instead of a data.table. This means that subsampling it with the "double dot data.table notation" command used for most other species: ```NC_cells<-inData[,..nc_column_idcs]
``` will result in an error. Use instead the data.frame syntax (without the "double dots"): ```NC_cells<-inData[,nc_column_idcs]```.

The reference transcriptome for Xt is [9.0](http://ftp.xenbase.org/pub/Genomics/JGI/Xentr9.0/Xtropicalisv9.0.Named.primaryTrs.fa.gz) (link will download the transcriptome), **not** the latest, because the names have different gene IDs, and if we use the latest, we won't be able to automatically find all transcripts. Because we're running analyses on the sequences, this is not a big deal. Besides, this was the exact reference transcriptome used for the [original scRNA-seq study](https://www.science.org/doi/10.1126/science.aar5780), see methods in their [supplementary materials](https://www.science.org/action/downloadSupplement?doi=10.1126%2Fscience.aar5780&file=aar5780_briggs_sm.pdf).
