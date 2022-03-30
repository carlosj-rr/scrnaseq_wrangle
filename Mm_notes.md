### _Mus musculus_-specific hacks
The mouse datasets have to be pre-processed a bit before reaching the point in which we can get to the R scripts.

First, as mentioned in the README, the data has to be imported from a package in Bioconductor, and this is a process of 2 steps (assuming you already imported the library):
````
inStudy<-MouseGastrulationData::EmbryoAtlasData()
inData<-assay(inStudy,"counts")
````
The dataset imported is a data.frame, not a data.table, so the subsampling of the NC columns must be modified from the 'data.table' 'double period' notation: ```NC_cells<-inData[,..nc_column_idcs]```, into the regular subsampling notation of data.frames: ```NC_cells<-inData[,nc_column_idcs]```.

For the reference transcriptome, I extracted the list of gene ids expressed in 50% or more cells from the "Mm_nc_expressionTable.csv" file, and used it as a BioMart query: ```tail -n+2 Mm_nc_expressionTable.csv | cut -d"," -f1 > BM_query.list # From line 2 onwards to avoid column headers.``` This file can now be given to [BioMart](https://www.ensembl.org/biomart/martview), which will fetch all the hits in the current database (GRCm39 as of March 2022), including all the transcripts. Request the cDNA sequences, and only the Gene ID in the FASTA header (not the version, not the transcript ID - this is important later on).

Once you've imported the BioMart results using ```Trome<-read.FASTA("filename")```, you can proceed normally with the analysis. You will notice that there are **6 'orphans'** - IDs with no hits on the references. These are outdated gene IDs that were eliminated when going from GRCm38 to GRCm39, so we can forget about them for the moment. This also means, however, that we have to remove them from the table of common neural genes ("outDF"). I've included an if block that does this, and explains that the object in the session was modified, suggests the older version to be overwritten, and mentions that the list of orphans was saved into a file called "orphans.csv". This block only activates if there were gene IDs that were not found in the reference transcriptome.