cut -d"," -f1 Dr_nc_barcodes.csv.tmp > Dr_nc_barcodes.csv
cut -d"," -f1 Hv_nc_barcodes.csv.tmp > Hv_nc_barcodes.csv
cut -d"," -f1 Mm_nc_barcodes.csv.tmp > Mm_nc_barcodes.csv
cut -d"," -f1 Sa_nc_barcodes.csv.tmp > Sa_nc_barcodes.csv
cut -d"," -f1 Sl_nc_barcodes.csv.tmp > Sl_nc_barcodes.csv
cut -d"," -f1 Sm_nc_barcodes.csv.tmp > Sm_nc_barcodes.csv
cut -d"," -f1-3 Xt_nc_barcodes.csv.tmp > Xt_nc_barcodes.csv # Xenopus tropicalis has more barcode information in the other columns, so best to keep a maximum.