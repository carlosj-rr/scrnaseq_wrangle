def col_nc_representation(column: pd.core.series.Series, nc_idcs: list) -> np.ndarray:
    col_ON=np.where(column > 0)
    return(np.intersect1d(col_ON,nc_idcs))

list1=[]
tot_cols=len(dm_data.columns)
for col in range(0,tot_cols):
    val=len(col_nc_representation(dm_data.iloc[:,col],l))
    list1.append(val)
    print(f"Finished column {col} of {tot_cols}")
    