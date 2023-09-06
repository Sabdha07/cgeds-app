import streamlit as st

import requests 
import io
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import pearsonr

### MAIN ###

# ---- PAGE CONFIG ----
st.set_page_config(page_title='CGEDS - Cancer Gene Expression and Drug Sensitivity', layout='wide')

# ---- HIDE STREAMLIT STYLE ----
hide_st_style = """
            <style>
            #MainMenu {visibility: hidden;}
            footer {visibility: hidden;}
            header {visibility: hidden;}
            </style>
            """
st.markdown(hide_st_style, unsafe_allow_html=True)

# Define the app version
app_version = "1.0.0"

#app
def main():
    # ---- HEADER SECTION ----
    with st.container():
        st.title('cGEDs - Cancer Gene Expression and Drug Sensitivity')
        st.write('An application for finding drug effectivity biomarkers in different cancer types')
        st.write('Version:', app_version)

    # ---- SIDEBAR ----
    st.sidebar.title("cGEDs")
    navigation_option = st.sidebar.radio(label= 'navigate',label_visibility='hidden', options = ["What is the app about?","cGEDs: Pick the genes and find the best drug", "User Guide"])    
    
    st.divider()
    # Add a contact section in the sidebar
    st.sidebar.write("## Contact")
    st.sidebar.write("Have questions, comments, or found a bug?")
    st.sidebar.write("Get in touch:")
    st.sidebar.write("Email: [sabwor07@gmail.com](mailto:sabwor07@gmail.com)")
    st.sidebar.write("-----")
    

	    
    if navigation_option == "User Guide":
        st.title("App Usage: User Guide")

        st.header("Dataset Selection")
        st.write("Choose between the GDSC1 and GDSC2 datasets. These datasets contain drug sensitivity and gene expression information for various cancer types and cell lines.")

        st.header("Cancer Type Selection")
        st.write("Select a specific cancer type from the available options. This will narrow down the analysis to the chosen cancer type.")

        st.header("Gene Selection")
        st.write("Pick genes of interest from the available list. The app will analyze the correlation between the selected genes' expression and drug sensitivity.")

        st.header("Correlation Calculation")
        st.write("Set the following parameters for correlation calculation:")
        st.write("- FDR Threshold: Choose the maximum False Discovery Rate (FDR) value for gene-drug pairs to be considered significant.")
        st.write("- Positive Correlation Threshold: Choose the minimum correlation value for positive correlations to be considered significant.")
        st.write("- Negative Correlation Threshold: Choose the maximum correlation value for negative correlations to be considered significant.")

        st.header("Data Visualization")
        st.write("You can choose to display calculated data and visualize it using the following options:")
        st.write("- Display calculated data: Check this option to display the calculated correlations between genes and drugs.")
        st.write("- Display heat map: Check this option to visualize the calculated correlations using a heatmap.")
    

    if navigation_option == "What is the app about?":
        
        st.markdown("## Advancing Cancer Insights for Enhaced Drug Selection")
        st.write("Cancer is a global health issue and a leading cause of death, responsible for approximately 9.6 million deaths in 2018. Traditional treatments face challenges due to the complexity of diverse cancer types and the unique genetic profiles of patients. The effectiveness of a drug depends on the specific interaction between genes and cancer cells. So, tailored treatment solutions become crucial as drugs may not perform uniformly across different cancer types and patients.")
        st.write("To assist with some of the challenges in cancer treatment, CGEDS has been developed with the motive to serve as a guide built upon existing data. CGEDS is designed to aid researchers and healthcare professionals in making more informed decisions regarding drug selection.")

        st.markdown("## How does cGEDs work?")
        st.write("CGEDS merges gene expression and drug sensitivity data(IC50 values) for the selected cancer type and using this combined information, it calculates the correlations between gene expressions and drug sensitivities. These correlations offer valuable insights, enabling researchers and healthcare professionals to analyze drug-gene interactions. This data-driven approach enhances our understanding of how drugs may interact with individual patient profiles, empowering more informed treatment decisions. The data used is sourced from publicly accessible repositories, ensuring comprehensive and credible input for the analysis.")


    elif navigation_option == "cGEDs: Pick the genes and find the best drug":
        #cGEDs app functions

        # ---- DATA ----
        # Load the first dataset
        @st.cache_data
        def load_dataset1():
            dataset = pd.read_csv('data/Drug-sensitivity-data-GDSC1.csv')
            #debug phase1
            dataset = dataset.rename(columns={"Cell line ": "cell_line", "Cell line": "cell_line"})
            return pd.DataFrame(dataset)

        # Load the second dataset
        @st.cache_data
        def load_dataset2():
            dataset = pd.read_csv('data/Drug-sensitivity-data-GDSC2.csv')
            #debug phase1
            dataset = dataset.rename(columns={"Cell line ": "cell_line", "Cell line": "cell_line"})
            return pd.DataFrame(dataset)

        #Load gene expression data
        @st.cache_data
        def load_geneexp():
            url = 'https://media.githubusercontent.com/media/Sabdha07/cgeds-app/main/data/Gene-expression-data-GDSC.csv'
            s=requests.get(url).content
            ex =pd.read_csv(io.StringIO(s.decode('utf-8')))
            
            #debug phase1
            ex = ex.rename(columns={"Cell line ": "cell_line", "Cell line": "cell_line"})
            return ex

        ## FUNCTIONS ##

        #heatmap
        def heatmap_plot(df):
            num_genes = len(df['Gene'].unique())
            num_drugs = len(df['Drug'].unique())
            figsize = max((5,3),(num_genes * 0.2, num_drugs * 0.2))
            
            plt.figure(figsize=figsize)

            # Pivot the DataFrame to get genes as rows and drugs as columns
            heatmap_data = df.pivot(index='Drug', columns='Gene', values='Correlation')
            
            ax = sns.heatmap(heatmap_data, center=0, annot=False, linewidth = .5)
            ax.set(xlabel= "Gene", ylabel="Drug")
            plt.xticks(fontsize=6) 
            plt.yticks(fontsize=max(4,num_drugs*0.05)) 
            #ax.xaxis.tick_top()
            plt.title("Gene-Drug Correlations Heatmap")
            #plt.tight_layout()
            st.pyplot(plt, use_container_width=False)

        #convert df to csv for download
        @st.cache_resource
        def convert_df(df):
	        return df.to_csv().encode('utf-8')

        #correlation calculation
        @st.cache_resource
        def calculate_correlations(df, fdr_thr, pos_cor_thr, neg_cor_thr, show_data, visualize):
            # Provide a vector of drug names
            # Remove drugs with less than two cell lines
            drug_counts = df.groupby('Drug.name')['cell_line'].nunique()
            drugs_to_retain = drug_counts[drug_counts > 2].index
            filtered_df = df[df['Drug.name'].isin(drugs_to_retain)]
            #initialize corrs dataframe
            corrs = pd.DataFrame(columns=['Gene','Drug','Correlation', 'FDR']) 

            #update progress 0 
            empty_placeholder = st.empty()
            t=0
            if t==0:
                empty_placeholder.write("Calculating, please wait..")
            
            # Calculate the correlations and FDRs for each drug separately    
            for drug in drugs_to_retain:
                drug_df = filtered_df[filtered_df['Drug.name'] == drug]
                gene_exp = drug_df.iloc[:,3:]
                ic50 = drug_df.iloc[:,2]
                corr_values = []
                p_values = []

                #update progress 1
                t=1
                if t==1:
                    empty_placeholder.write("Calculating, please wait...")

                #iterate over each gene
                for gene_column in gene_exp.columns:
                    gene_values = gene_exp[gene_column]
                    # Check if the gene expression values are constant
                    if np.all(gene_values == gene_values.iloc[0]):
                        corr = 0.0  # You can set any default value here
                        p_value = 1.0  # You can set any default value here
                    else:
                        #st.write(gene_values, ic50)
                        corr, p_value = pearsonr(gene_values, ic50)
                        
                    corr_values.append(corr)
                    p_values.append(p_value)
                    
                    #progress update 2        
                    t=2
                    if t==2:
                        empty_placeholder.write("Calculating, please wait......")

                # After calculating correlations for all genes within the current drug,
                # append the data to the 'corrs' DataFrame
                for i in range(len(gene_exp.columns)):
                    gene_column = gene_exp.columns[i]
                    corr = corr_values[i]
                    p_value = p_values[i]

                    if not ((corrs['Drug'] == drug) & (corrs['Gene'] == gene_column)).any():
                        corrs = pd.concat([corrs, pd.DataFrame({'Gene': [gene_column], 'Drug': [drug],
                                                                'Correlation': [corr],
                                                                'FDR': [p_value]})],
                                                                ignore_index=True)
                        #progress update 2        
                        t=2
                        if t==2:
                            empty_placeholder.write("Calculating, please wait..")

            #filter by thresholds
            sigcors1 = corrs[(corrs['FDR'] < fdr_thr) & (corrs['Correlation'] > pos_cor_thr)]
            sigcors2 = corrs[(corrs['FDR'] < fdr_thr) & (corrs['Correlation'] < neg_cor_thr)]
            sigcors = pd.concat([sigcors1, sigcors2])
            
            #final progress update
            t = 3
            if t==3:
                empty_placeholder.write("Calculation completed!")
            
            #final result
            if show_data:
                if visualize:
                    return sigcors, True
                else:
                    return sigcors, False
            else:
                if visualize:
                    return None, True
                else:
                    return None, False

        ## ---- PIPELINE ---- ##

        #SELECT DATASET
        col1, col2 = st.columns(2)
        dataset_selected = col1.selectbox("Select a drug sensitivity and gene expression dataset", ('Select dataset',"GDSC1", "GDSC2"))


        col2.markdown(
            """
            <details>
                <summary><strong>Info</strong></summary>
                <p>You can choose the drug sensitivity and gene expression dataset among these publicly available datasets:</p>
                <p>GDSC1: 970 Cell lines and 403 Compounds</p>
                <p>GDSC2: 969 Cell lines and 297 Compounds</p>
            </details>
            """,
            unsafe_allow_html=True
        )
        if dataset_selected == "GDSC1":
            selected_dataset = load_dataset1()
        else:
            selected_dataset = load_dataset2()


        #SELECT CANCER TYPE
        options = np.insert(selected_dataset["Cancer-Type"].unique(),[0],'Select cancer type')
        cancer = st.selectbox(
            "Select a cancer type", options)

        ds = pd.DataFrame(selected_dataset[selected_dataset["Cancer-Type"] == cancer])
        ds = ds.drop(columns=["Cancer-Type"])


        #SELECT GENES
        ex = load_geneexp()
        #st.write(ex.shape)
        genes = st.multiselect("Select your desired genes", options=(ex.columns[1:]))
        st.markdown(
            """
            <details>
                <summary><strong>Note:</strong></summary>
                <p>For better visualization in the heat map, it is recommended to select at least three genes. Consider this as an optional suggestion, not a requirement.</p>
            </details>
            """,
            unsafe_allow_html=True
        )

        #filtering off the unselected genes
        ex_filtered = pd.DataFrame(columns=['cell_line'])
        if genes:
                    ex_filtered = pd.DataFrame(pd.concat([ex["cell_line"], ex[genes]], axis=1))
                    #st.write("Filtered DataFrame:")
                    #st.dataframe(ex_filtered)
        else:
                    st.warning("No genes selected.")

        #merge drugs and gene expression
        df = pd.merge(ds, ex_filtered, on="cell_line")
        if not df.empty:
            st.write('Raw data:', df)
            #has columns drug, cell_line, ic50, gene expression

        with st.container():
            fdr_thr = str(st.number_input("Choose Gene-drug pairs with FDRs less than:", value=0.05))
            pos_cor_thr = str(st.slider("Choose Gene-drug pairs with correlations more than:", 0.0, 1.0, 0.7, step=0.1))
            neg_cor_thr = st.slider("Choose Gene-drug pairs with correlations less than:", -1.0, 0.0, -0.7, step=0.1)

        with st.container():
            show_data = st.checkbox("Display calculated data")
            visualize = st.checkbox('Display heat map')

        #Calculate correlations
        if st.button('Calculate Correlations'):  
            fdr = float(fdr_thr)  
            pos = float(pos_cor_thr)  
            neg = float(neg_cor_thr)
            showdata = show_data
            vis = visualize
            calculated_corrs, show_visualization = calculate_correlations(df, fdr, pos, neg,showdata,vis)
            if not calculated_corrs.empty:
                if show_data:
                    st.write("Calculated and Filtered Data:")
                    st.write(calculated_corrs)
    
                    if show_visualization:
                        st.write("Heat Map for all genes and drug pairs:")
                        heatmap_plot(calculated_corrs)
                            
                else:
                    if show_visualization:
                        st.write("Heat Map for all genes and drug pairs:")
                        heatmap_plot(calculated_corrs)
                

		        
            else:
                st.write("No gene-drug pairs with the given threshold values were found.")

        #download button
        if st.button("Download data as CSV"):
            fdr = float(fdr_thr)  
            pos = float(pos_cor_thr)  
            neg = float(neg_cor_thr)
            showdata = show_data
            vis = visualize
            calculated_corrs, show_visualization = calculate_correlations(df, fdr, pos, neg,showdata,vis)
            download_csv = convert_df(calculated_corrs)

            st.download_button(label="Download", data=download_csv, file_name='calculated_correlations.csv', mime='text/csv')


    
if __name__ == "__main__":
    main()

