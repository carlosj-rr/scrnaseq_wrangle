for cluster in $(cut -d"," -f2 Neural_and_Muscle_clusterNames.csv | sed 's/ /+/g' | tail -n+2 )
do
    spp=$(echo $cluster | cut -d"_" -f1)
    cluster=$(echo $cluster | sed 's/+/ /g' | cut -d"_" -f2-)
    echo $spp,$cluster >> spp_clusterName.csv
done