import pandas as pd


for pop in pop2ko2life.keys():
    print(pop)
    pop_list = list(pop2spes[pop])
    pop_list_ = []
    for spe in pop_list:
        if spe in spe2id.keys():
            spe = spe2id[spe]
        pop_list_.append(spe)
    kegg_df_pop = kegg_df[pop_list_]
    kegg_df_pop = kegg_df_pop.T
    kegg_df_pop.fillna(0, inplace=True)
    correlation_matrix = kegg_df_pop.corr()
    nifdf = pd.DataFrame(columns = ['ko','pathway','genename','corr','nif_copy','non-nif_copy','nif_ratio','non-nif_ratio','description'])
    nondf = pd.DataFrame(columns = ['ko','pathway','genename','corr','nif_copy','non-nif_copy','nif_ratio','non-nif_ratio','description'])
    key_ko = 'K02588'
    for ko in correlation_matrix.index:
        corr = correlation_matrix.loc[ko,key_ko]
        if corr >= 0.6:
            nifdf.loc[ko,'ko'] = ko
            nifdf.loc[ko,'genename'] = ko2gene[ko]
            nifdf.loc[ko,'corr'] = round(corr,2)
            nifdf.loc[ko,'nif_copy'] = round(pop2ko2life[pop][ko]['nif']['meancopy'],2)
            nifdf.loc[ko,'non-nif_copy'] = round(pop2ko2life[pop][ko]['non-nif']['meancopy'],2)
            nifdf.loc[ko,'nif_ratio'] = round(pop2ko2life[pop][ko]['nif']['ratio'],2)
            nifdf.loc[ko,'non-nif_ratio'] = round(pop2ko2life[pop][ko]['non-nif']['ratio'],2)
            nifdf.loc[ko,'description'] = ko2descrip[ko]
            sorted_nifdf = nifdf.sort_values(by=["corr", "ko"], ascending=[False, True])
        elif corr <= -0.6:
            nondf.loc[ko,'ko'] = ko
            nondf.loc[ko,'genename'] = ko2gene[ko]
            nondf.loc[ko,'corr'] = round(corr,2)
            nondf.loc[ko,'nif_copy'] = round(pop2ko2life[pop][ko]['nif']['meancopy'],2)
            nondf.loc[ko,'non-nif_copy'] = round(pop2ko2life[pop][ko]['non-nif']['meancopy'],2)
            nondf.loc[ko,'nif_ratio'] = round(pop2ko2life[pop][ko]['nif']['ratio'],2)
            nondf.loc[ko,'non-nif_ratio'] = round(pop2ko2life[pop][ko]['non-nif']['ratio'],2)
            nondf.loc[ko,'description'] = ko2descrip[ko]
            sorted_nondf = nondf.sort_values(by=["corr", "ko"], ascending=[False, True])
    with pd.ExcelWriter(f'{pop}.xlsx',engine='openpyxl') as writer:
        sorted_nifdf.style.apply(color_rows, axis=1).to_excel(writer,sheet_name='nif',index=False)
        sorted_nondf.to_excel(writer,sheet_name='non-nif',index=False)