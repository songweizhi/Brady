import pandas as pd

# https://pandas.pydata.org/docs/user_guide/merging.html



def rm_single_value_cols(df_in):

    # Check if there is only one unique value in each column
    unique_values = df_in.nunique()
    same_value_cols = unique_values[unique_values == 1].index

    # Drop columns with the same value
    df_out = df_in.drop(same_value_cols, axis=1)

    return df_out


# df_dict = {'A': [1,2,3,4], 'B': [5,5,5,5], 'C': [6,7,8,9]}
# df      = pd.DataFrame(df_dict)
# df_out  = rm_single_value_cols(df)


# read in dataframe
df_file                           = '/Users/songweizhi/Desktop/Pathways_PA_carbo.txt'
df_file_without_single_value_cols = '/Users/songweizhi/Desktop/Pathways_PA_carbo_without_single_value_cols.txt'

df = pd.read_csv(df_file, sep='\t', header=0, index_col=0)
df_without_single_value_cols = rm_single_value_cols(df)
df_without_single_value_cols.to_csv(df_file_without_single_value_cols, sep='\t')


# print('df')
# print(df)
# print()
# print('df_out')
# print(df_out)
# print()



df1 = pd.DataFrame(
    {   "A": ["A0", "A1", "A2", "A3"],
        "B": ["B0", "B1", "B2", "B3"],
        "C": ["C0", "C1", "C2", "C3"],
        "D": ["D0", "D1", "D2", "D3"]},
    index=['0', '1', '2', '3'])

df2 = pd.DataFrame(
    {   "A": ["A4", "A5", "A6", "A7"],
        "B": ["B4", "B5", "B6", "B7"],
        "C": ["C4", "C5", "C6", "C7"],
        "D": ["D4", "D5", "D6", "D7"]},
    index=['4', '5', '6', '7'])

df3 = pd.DataFrame(
    {   "A": ["A8", "A9", "A10", "A11"],
        "B": ["B8", "B9", "B10", "B11"],
        "C": ["C8", "C9", "C10", "C11"],
        "D": ["D8", "D9", "D10", "D11"]},
    index=['8', '9', '10', '11'])

#
# print('df1')
# print(df1)
# print()
# print('df2')
# print(df2)
# print()
# print('df3')
# print(df3)
# print()
