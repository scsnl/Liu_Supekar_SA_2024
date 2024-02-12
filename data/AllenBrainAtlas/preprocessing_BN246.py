
# please install abagen in terminal first
# by running the following:
# "pip install abagen
# 

import abagen
import numpy
import pandas as pd
from abagen import images

files = abagen.fetch_microarray(donors = 'all', resume=True, data_dir = '~/data/') # downloding data
# atlas = abagen.fetch_desikan_killiany() 
# set parcellation
atlas = {'image':'C:\\Users\\jinliu5\\AppData\\Local\\Programs\\Python\\Python310\\Lib\\site-packages\\abagen\\data\\Reslice_BN_Atlas_246_1mm.nii.gz'} 
# preprocessing of gene expression with default setting
expression = abagen.get_expression_data(atlas['image'], donors='all', tolerance=0)
# save the preprocessed gene expression data
print(type(expression))
expression.to_csv("BN246_geneexpression.csv")

# report the setting
from abagen import reporting
generator = reporting.Report(atlas['image'], tolerance=0)
report = generator.gen_report().encode("utf8")
with  open("preprocessing.txt", "wb") as file:
	file.write(report)
	file.close()
