### scRNA-seq wrangle README
This repository has all the information needed to produce species-specific lists of transcripts from neural cells based on the SAMap clustering published by [Tarashansky et al. (2021)](https://elifesciences.org/articles/66747). The important files are listed below:
***
1. General_guide
  * Explains the general order of steps that are to be taken to precisely replicate the results in our (yet unpublished) work, and which results from Tarashansky et al. (2021) were needed for this process.
***
2. Dataset_links.md
  * Explicit listing of all the (raw data) databases accessed for each species
    * This covers the reference transcriptomes and the scRNA-seq expression tables
***
3. Species_specific_notes.md
  * Any and all species-specific modifications to the general steps detailed in the "General_guide" file in point #1
    * Note: some `bash` will be needed.
***
4. processo.R
  * A set of commands to be run (by the user, *not* as an automated script - it may need some species-specific modifications, explained in the notes above) in R that take as input a scRNA-seq expression levels table and will result in three files:
    1. A reduced scRNA-seq expression levels table that *only* has neural cells and, only the genes which were expressed in these.
    2. A table of gene counts for each of the genes expressed in 50% or more neural cells.
    3. The list of the canonical transcript of each of the genes in the list in the point (ii) above.
