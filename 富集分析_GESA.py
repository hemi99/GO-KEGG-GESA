import pandas as pd
import gseapy as gp
import requests
import sleep
import matplotlib.pyplot as plt
from gseapy.parser import Biomart
from gseapy.plot import barplot, dotplot


gene_list = pd.read_csv(r"D:\graduate\gene_list.txt",header=None, sep="\t")
gene_list1 = pd.read_csv(r"D:\graduate\non_geneID.csv")
gene_list1.head()

glist = gene_list1.squeeze().str.strip().tolist()
names = gp.get_library_name() # default: Human

s = requests.session()
s.keep_alive = False



enr = gp.enrichr(gene_list=r"D:/graduate/gene_list.txt",
     # or gene_list=glist
     description='',
     gene_sets=['KEGG_2019_Human'],
     outdir='test/enrichr_kegg',
     cutoff=0.5 # test dataset, use lower value from range(0,1)
    )



barplot(enr.res2d,title='KEGG_2019_Human',)
dotplot(enr.res2d, title='KEGG_2019_Human',)
