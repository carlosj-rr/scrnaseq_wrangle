import numpy as np
import pandas as pd

#read in the emapper annotations file
ann_dataset=pd.read_csv("Dr_refGenome-EGNG.emapper.annotations",sep="\t")
#read in the list of inferred nervous system genes
nc_genes=np.genfromtxt("../../Dr_ncGenes.csv",delimiter=",",dtype=str)

#see which one of the inferred NS genes have an annotation entry in the loaded dataset
ng_annotated=np.intersect1d(nc_genes,ann_dataset.iloc[:,0])

# make a list of all the row numbers in which a NS gene has its annotations
l=[]
for i in ng_annotated:
        exp=np.where(i == ann_dataset.iloc[:,0])[0]
        l.append(exp)

#Tidy list and make pretty:
l.sort()
l=np.array(l).flatten()

#BEAUTY!

np.any(np.in1d(nc_markers,ds_data.index))
np.any(np.in1d(nc_markers,ds_data.columns))

row_wise_sum=np.sum((dm_data > 0).astype(int),axis=1)
np.where(row_wise_sum >= dm_data.shape[1])[0].size


import nov_parser as novparser
np=novparser.np
pd=novparser.pd

ds_data=pd.read_csv("Mm_scRNAseq.csv")
nc_markers=np.array(pd.read_csv("NeuralMarkers_mouse.csv", sep=",")).flatten()

ds_data=novparser.correct_axes_orientation(ds_data,nc_markers)
nc_gene_ss=ds_data.loc[nc_markers,:]
gt0=np.sum((nc_gene_ss > 0).astype(int),axis=0)
novparser.top_tenperc_nc_gene_ids_extractor(gt0, ds_data, "Mm_ncGenes.csv")


import scipy as sc
import pandas as pd

ds_data=sc.io.mmread("matrix.mtx")
ds_data=ds_data.todense()
colnames=np.genfromtxt("barcodes.tsv",delimiter="\t",dtype=str)
for_rows=pd.read_csv("features.tsv",sep="\t",header=None)
rownames=np.array(for_rows.iloc[:,0])
cg_data=pd.DataFrame(ds_data,index=rownames,columns=colnames)

np.savetxt("Cg_bioMart_request.list",rownames,fmt='%s', delimiter=",")