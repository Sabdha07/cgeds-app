import streamlit as st

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


def main():
    st.sidebar.title("cGEDs - Battling Cancer, One Gene at a Time!")
    navigation_option = st.sidebar.radio(label= 'navigate',label_visibility='hidden', options = ["What is the app about?","cGEDs: Pick the genes and find the best drug", "User Guide"])    
    
    st.divider()
    # Add a contact section in the sidebar
    st.sidebar.write("## Contact")
    st.sidebar.write("Have questions, comments, or found a bug?")
    st.sidebar.write("Get in touch:")
    st.sidebar.write("Email: [sabwor07@gmail.com](mailto:sabwor07@gmail.com)")
    

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
        st.write("- Display heat map and cluster bar plot: Check this option to visualize the calculated correlations using a heatmap and a clustered bar plot.")



    if navigation_option == "What is the app about?":
        st.title("About cGEDs - Cancer Gene Expression and Drug Sensitivity")
        st.markdown("🔬 **Unleashing Precision Medicine: CGEDS - The Cancer Genomics Guide for Drug Selection**")
        st.write("Battling Cancer, One Gene at a Time")
        
        st.markdown("🌍 **Global Impact of Cancer**")
        st.write("- Cancer, a relentless adversary, affects millions of lives globally.")
        st.write("- Its indiscriminate nature and the dire need for personalized treatments have intensified the quest for groundbreaking solutions.")
        st.write("- According to the World Health Organization (WHO), cancer is a leading cause of death worldwide, with an estimated 9.6 million deaths in 2018.")
        st.write("- Conventional treatments often fall short due to the complex interplay between diverse cancer cell lines and the intricate genetic makeup of individual patients.")
        
        st.markdown("## The Challenge: A Symphony of Variability")
        st.write("🎭 **Diverse and Complex Cancers**")
        st.write("- Cancers are as diverse as the patients they afflict, with varying gene expressions that drive their behaviors.")
        st.write("- The efficacy of a drug depends on this complex interaction between genes and cancer cell lines.")
        st.write("- The same drug might hold a silver bullet for one type of cancer but fall short for another.")
        st.write("- Traditional treatments, often generalized, miss the mark due to this intricate variability.")
        
        st.markdown("## CGEDS - The Personalized Cancer Guide for Drug Selection")
        st.write("🚀 **Revolutionizing Cancer Treatment**")
        st.write("- CGEDS is designed to revolutionize the landscape of cancer treatment.")
        st.write("- The app is a manifestation of the fusion of gene expression data and drug sensitivity profiles, all united with the shared goal of creating a better future for cancer patients.")
        
        st.markdown("### Empowering Precision Medicine:")
        st.write("💡 **Predicting Tumor Response**")
        st.write("- At the heart of CGEDS lies the ability to predict tumor response to drugs based on gene-expression biomarkers.")
        st.write("- This unprecedented power stems from the integration of diverse datasets from cancer cell lines, providing a rich source of information about drug efficacy.")
        
        st.markdown("## How It Works")
        st.write("🚀 **Empowering Clinicians and Researchers**")
        st.write("- CGEDS integrates gene expression and drug sensitivity profiles from an extensive range of cancer cell line datasets.")
        st.write("- By doing so, it calculates correlations between gene expressions and drug sensitivities across multiple cancer types, transcending the limitations of conventional analyses.")
        
        st.markdown("## Why CGEDS is the Game-Changer")
        st.write("1. 🎯 **Personalization in Action**")
        st.write("   The app transforms a one-size-fits-all approach into a personalized strategy, offering insights into the most effective treatments based on individual cancer profiles.")
        st.write("2. 🧪 **Research Meets Reality**")
        st.write("   The app's predictive models bridge the gap between laboratory discoveries and clinical applications.")
        st.write("   CGEDS translates intricate gene-expression interactions into tangible treatment recommendations, emboldening healthcare professionals in their pursuit of optimal patient outcomes.")
        st.write("3. 🌌 **A Universe of Insights**")
        st.write("   With the power to analyze multiple cancer types, CGEDS reveals a lot of novel insights that were previously hidden in the labyrinth of data.")
        st.write("   By bringing together data from different sources, the app paves the way for unprecedented discoveries.")
        st.write("4. 🎯 **Targeted Medicine, Redefined**")
        st.write("   CGEDS harnesses the potential of targeted medicine by identifying drugs that hold promise against specific cancer types.")
        
    

    elif navigation_option == "cGEDs: Pick the genes and find the best drug":
        #cGEDs app functions

        # ---- DATA ----
        # Load the first dataset
        @st.cache_data
        def load_dataset1():
            dataset = pd.read_csv('data/Drug-sensitivity-data-GDSC1.csv')
            return pd.DataFrame(dataset)

        # Load the second dataset
        @st.cache_data
        def load_dataset2():
            dataset = pd.read_csv('data/Drug-sensitivity-data-GDSC2.csv')
            return pd.DataFrame(dataset)

        #Load gene expression data
        @st.cache_data
        def load_geneexp():
            ex = pd.read_csv('data/Gene-expression-data-GDSC.csv')
            return pd.DataFrame(ex)

        ## FUNCTIONS ##

        #clustered bar plot
        def create_clustered_bar_plot(df):
            num_genes = len(df['Gene'].unique())
            num_drugs = len(df['Drug'].unique())
            figsize = (num_genes * 3, 10)  # Adjust the scaling factor as needed
            
            plt.figure(figsize=figsize)
            sns.set_palette('bright')
            ax = sns.barplot(data=df, x='Gene', y='Correlation', hue='Drug')
            
            # Set labels and title
            plt.xlabel('Gene')
            plt.ylabel('Correlation')
            plt.title('Clustered Bar Plot of Gene-Drug Correlations')
            
            plt.xticks(rotation=45)
            plt.legend(title='Drug', bbox_to_anchor=(1.05, 1), loc='upper left')
            
            # Annotate the bars with their correlation values
            for p in ax.patches:
                height = p.get_height()
                annotation_text = format(height, ".2f")
                text_y = height + 0.05 if height >= 0 else height - 0.05
                text_color = 'black' 

                ax.annotate(annotation_text, 
                            (p.get_x() + p.get_width() / 2., height), 
                            ha='center', va='center', 
                            xytext=(0, 10) if height >= 0 else (0, -10), 
                            textcoords='offset points',
                            rotation=90,
                            color=text_color)
            
            plt.tight_layout()
            
            st.pyplot(plt)

        #heatmap
        def heatmap_plot(df):
            num_genes = len(df['Gene'].unique())
            num_drugs = len(df['Drug'].unique())
            figsize = (num_genes * 1.3, num_drugs * 1.3)  # Adjust the scaling factor as needed
            
            plt.figure(figsize=figsize)

            # Pivot the DataFrame to get genes as rows and drugs as columns
            heatmap_data = df.pivot(index='Drug', columns='Gene', values='Correlation')
            
            sns.heatmap(heatmap_data, cmap='Spectral', center=0, annot=True, fmt=".2f")
            plt.title("Gene-Drug Correlations Heatmap")
            plt.tight_layout()
            st.pyplot(plt)


        #correlation calculation
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


        # ---- HEADER SECTION ----
        with st.container():
            st.title('cGEDs - Cancer Gene Expression and Drug Sensitivity')
            st.write('An application for finding drug effectivity biomarkers in different cancer types')
            
        st.divider()

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

        #filtering off the unselected genes
        ex_filtered = pd.DataFrame(pd.concat([ex["Cell line"], ex[genes]], axis=1))

        #debug phase1
        ds = ds.rename(columns={"Cell line ": "cell_line", "Cell line": "cell_line"})
        ex_filtered = ex_filtered.rename(columns={"Cell line": "cell_line"})

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
            visualize = st.checkbox('Display heat map and cluster bar plot')

        #Calculate correlations
        if st.button('Calculate Correlations'):  
            fdr = float(fdr_thr)  
            pos = float(pos_cor_thr)  
            neg = float(neg_cor_thr)
            showdata = show_data
            vis = visualize
            calculated_corrs, show_visualization = calculate_correlations(df, fdr, pos, neg,showdata,vis)
            
            if show_data:
                st.write("Calculated Data:")
                st.write(calculated_corrs)

                if show_visualization:
                    st.write("Heat Map for all genes and drug pairs:")
                    heatmap_plot(calculated_corrs)
                    st.divider()
                    st.write("Clustered bar plot for all genes and drug pairs:")
                    create_clustered_bar_plot(calculated_corrs)

            else:
                if show_visualization:
                    st.write("Heat Map for all genes and drug pairs:")
                    heatmap_plot(calculated_corrs)
                    st.divider()
                    st.write("Clustered bar plot for all genes and drug pairs:")
                    create_clustered_bar_plot(calculated_corrs)



    
if __name__ == "__main__":
    main()

