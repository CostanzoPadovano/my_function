import scipy.stats as stats
import pandas as pd
import numpy as np
import seaborn as sns
from scipy.stats import pearsonr
from scipy.stats import spearmanr
import matplotlib.pyplot  as plt
import os
#NOTE: Data must to have genes name in columns and sample in the rows (index)
#this function plot one genes x to all genes in the dataset called data

def significative_correlation(data, xgene="", method="pearsonr", correlation="direct", savefig=False):
    try:
        if savefig==True:
            # define the name of the directory to be created
            try:
                npath=0
                path = os.getcwd()+"/"+f"{methods}_Pictures_genes_plot_{npath}"
                os.makedirs(path)
            except:
                npath=+1
                path=os.getcwd()+"/"+f"{methods}_Pictures_genes_plot_{npath}"
                os.makedirs(path)
                print ("Successfully created the directory %s" % path)
        else:
            pass        
        if methods=="pearsonr" and correlation=="direct":
            for i in list(data.columns):
                x=xgene
                y= i
                df=data
                # ax.annotate(stats.pearsonr)
                r, p = stats.pearsonr(df[x],df[y])
                if p<0.05 and r > 0:
                    g = sns.jointplot(data=df,x=x, y=y, kind='reg', color='royalblue')
                    g.ax_joint.annotate(f'$\\rho = {r:.3f}, p = {p:.3f}$',
                                        xy=(0.1, 0.9), xycoords='axes fraction',
                                        ha='left', va='center',
                                        bbox={'boxstyle': 'round', 'fc': 'powderblue', 'ec': 'navy'})
                    plt.tight_layout()
                    if savefig==True:
                        plt.savefig(f"{os.getcwd()}/Pictures_genes_plot_{npath}/{i}vs{xgene}.pdf")
                    else:
                        pass
                    plt.show()
                else:
                    pass
        elif methods=="spearmanr" and correlation=="direct": 
            for i in list(data.columns):
                x=xgene
                y= i
                df=data
                # ax.annotate(stats.pearsonr)
                r, p = stats.spearmanr(df[x],df[y])
                if p<0.05 and r > 0:
                    g = sns.jointplot(data=df,x=x, y=y, kind='reg', color='royalblue')
                    g.ax_joint.annotate(f'$\\rho = {r:.3f}, p = {p:.3f}$',
                                         xy=(0.1, 0.9), xycoords='axes fraction',
                                         ha='left', va='center',
                                         bbox={'boxstyle': 'round', 'fc': 'powderblue', 'ec': 'navy'})
                    plt.tight_layout()
                    if savefig==True:
                        plt.savefig(f"{os.getcwd()}/Pictures_genes_plot_{npath}/{i}vs{xgene}.pdf")
                    else:
                        pass
                    plt.show()
        elif methods=="pearsonr" and correlation=="indirect": 
            for i in list(data.columns):
                x=xgene
                y= i
                df=data
                # ax.annotate(stats.pearsonr)
                r, p = stats.spearmanr(df[x],df[y])
                if p<0.05 and r < 0:
                    g = sns.jointplot(data=df,x=x, y=y, kind='reg', color='royalblue')
                    g.ax_joint.annotate(f'$\\rho = {r:.3f}, p = {p:.3f}$',
                                         xy=(0.1, 0.9), xycoords='axes fraction',
                                         ha='left', va='center',
                                         bbox={'boxstyle': 'round', 'fc': 'powderblue', 'ec': 'navy'})
                    plt.tight_layout()
                    if savefig==True:
                         plt.savefig(f"{os.getcwd()}/Pictures_genes_plot_{npath}/{i}vs{xgene}.pdf")
                    else:
                        pass
                    plt.show()
        elif methods=="spearmanr" and correlation=="indirect" : 
            for i in list(data.columns):
                x=xgene
                y= i
                df=data
                # ax.annotate(stats.pearsonr)
                r, p = stats.spearmanr(df[x],df[y])
                if p<0.05 and r < 0:
                    g = sns.jointplot(data=df,x=x, y=y, kind='reg', color='royalblue')
                    g.ax_joint.annotate(f'$\\rho = {r:.3f}, p = {p:.3f}$',
                                         xy=(0.1, 0.9), xycoords='axes fraction',
                                         ha='left', va='center',
                                         bbox={'boxstyle': 'round', 'fc': 'powderblue', 'ec': 'navy'})
                    plt.tight_layout()
                    if savefig==True:
                        plt.savefig(f"{os.getcwd()}/Pictures_genes_plot_{npath}/{i}vs{xgene}.pdf")
                    else:
                        pass
                    plt.show()
        elif methods=="spearmanr" and correlation=="all" : 
            for i in list(data.columns):
                x=xgene
                y= i
                df=data
                # ax.annotate(stats.pearsonr)
                r, p = stats.spearmanr(df[x],df[y])
                g = sns.jointplot(data=df,x=x, y=y, kind='reg', color='royalblue')
                g.ax_joint.annotate(f'$\\rho = {r:.3f}, p = {p:.3f}$',
                                         xy=(0.1, 0.9), xycoords='axes fraction',
                                         ha='left', va='center',
                                         bbox={'boxstyle': 'round', 'fc': 'powderblue', 'ec': 'navy'})
                plt.tight_layout()
                if savefig==True:
                    plt.savefig(f"{os.getcwd()}/Pictures_genes_plot_{npath}/{i}vs{xgene}.pdf")
                else:
                    pass
                plt.show()        
    except KeyError:
        print("You must add xgene='name of gene' and correlation= 'direct' or 'indirect'or 'all'")
        

        
