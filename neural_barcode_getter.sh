for sp in $(cut -d"," -f1 spp_clusterName.csv | sort | uniq)
do
    queryfile=$(echo ${sp}_clusters.csv)
    outfile=$(echo ${sp}_nc_barcodes.csv)
    for annotation in $(grep ${sp} spp_clusterName.csv | cut -d"," -f2 | sed 's/ /+/g')
    do
        annotation=$(echo $annotation | sed 's/+/ /g')
        grep "$annotation" ${queryfile} >> ${outfile}.tmp
    done
done