import pandas as pd
import gseapy as gp

def run_gsva_and_update_adata(adata, gene_sets, min_size=1, max_size=500, outdir='gsva_results'):
    """
    Run GSVA on AnnData object with user-defined gene sets and update the AnnData object with GSVA scores.

    Parameters:
    adata (AnnData): The AnnData object containing single-cell RNA-seq data.
    gene_sets (dict): A dictionary of gene sets. Each key is a gene set name and the value is a list of genes.
    min_size (int): Minimum size of gene sets to be considered.
    max_size (int): Maximum size of gene sets to be considered.
    outdir (str): Directory to store GSVA results.
    
    Returns:
    AnnData: Updated AnnData object with GSVA scores in adata.obs
    """
    
    # Convert AnnData to DataFrame and transpose
    expr_data = adata.to_df().T
    
    # Run GSVA
    gsva_results = gp.gsva(expr_data, gene_sets, method='gsva', outdir=outdir, 
                            min_size=min_size, max_size=max_size)

    # Extract GSVA scores
    gsva_score = gsva_results.res2d.pivot(index='Term', columns='Name', values='ES').T

    # Update AnnData with GSVA scores
    adata.obs = pd.concat([adata.obs, gsva_score], join="inner", axis=1)

    return adata
