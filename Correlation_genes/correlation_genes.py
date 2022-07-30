import pandas as pd
import numpy as np
import seaborn as sns
from scipy.stats import pearsonr
from scipy.stats import spearmanr
import matplotlib.pyplot  as plt
import os
import streamlit as st
import altair as alt
#NOTE: Data must to have genes name in columns and sample in the rows (index)
#this function plot one genes x to all genes in the dataset called data


######################
# Input Text Box
######################

#st.sidebar.header('Enter DNA sequence')
st.header('Upload csv data')



uploaded_file = st.file_uploader("Choose a file format csv")
# Check if file was uploaded
if uploaded_file is not None:
    st.header('Dataset')
    data = pd.read_csv(uploaded_file,index_col=0)
    st.write(data)

    #Transposition data
    Transpose_data = st.text_input("Do you want to transpose data? (Yes/No)")
    NOTE= 'if the genes are in the rows, you must to transpose the dataset'
    NOTE

    if Transpose_data=="Yes":
        data = data.T
        st.success("Data Transposed")
        st.write(data)
    elif Transpose_data=="No":
        st.success("Data not Transposed")
    else:
        st.warning("Please, add answer Yes or No ")

    #OPTIONS
    methods = st.selectbox(
     'What kind of method do you want to use to correlate genes?',
     ('pearsonr', 'spearmanr'))

    if methods=="pearsonr":
        pass
    elif methods=="No":
        st.success("data not Transposed")

    if methods=="spearmanr":
        pass
    elif methods=="No":
        st.success("data not Transposed")

    xgene = st.text_input("what kind of gene do you want to correlate with all genes?")

    if xgene is not None:
        if xgene in list(data.columns):
            st.success("Gene selected")

            #Principal Function
            def significative_correlation(data, xgene=xgene, method=methods):
                dicts={}
                values=data.columns
                for i in data.columns:
                    x= xgene
                    y= i
                    # ax.annotate(stats.pearsonr)
                    if method =="spearmanr":
                        r, p = spearmanr(data[x],data[y])
                        dicts[i] = r,p
                    elif method == "pearsonr":
                        r, p = pearsonr(data[x],data[y])
                        dicts[i] = r,p

                df_r_pvalue=pd.DataFrame(dicts).T.rename(columns={0:"r",1:"p_value"})

                #color bar of dataframe df_r_pvalue
                def color_direct_indirect_corr(val):
                    if val < 0:
                        color ='red'
                    elif val == 0:
                        color = 'white'
                    else:
                        color = 'blue'
                    return f'background-color: {color}'

                #Download Button csv
                @st.cache
                def convert_df(df):
                     # IMPORTANT: Cache the conversion to prevent computation on every rerun
                     return df.to_csv().encode('utf-8')
                csv = convert_df(df_r_pvalue)

                st.download_button(
                     label="Download data as CSV",
                     data=csv,
                     file_name='large_df.csv',
                     mime='text/csv',
                 )

                st.dataframe(df_r_pvalue.style.applymap(color_direct_indirect_corr, subset=['r']))

                #Create a dotplot
                Dotplot = st.text_input("Do you want to generate a Dotplot? (Yes/No)")

                if Dotplot=="Yes":
                    df_r_pvalue["group"]="Gene_dataset"
                    df_r_pvalue=df_r_pvalue[df_r_pvalue["p_value"]<0.05]
                    c=alt.Chart(df_r_pvalue.reset_index()).mark_circle().encode(
                        x=alt.X('group', axis=alt.Axis(grid=True)),
                        y=alt.Y('index', axis=alt.Axis(grid=True)),
                        size=alt.Size("p_value",scale=alt.Scale(domain=[0.05,0])),
                        color=alt.Color("r",scale=alt.Scale(domain=[1,-1],scheme="blueorange")))
                    st.altair_chart(c, use_container_width=False)
                elif Dotplot=="No":
                    st.success("Data not Transposed")
                else:
                    st.warning("Please, add answer Yes or No ")



            #apply the function
            significative_correlation(data)



        elif xgene is not list(data.columns):
            st.warning("gene is not in the data")
        else:
            st.warning("Please add gene name")
