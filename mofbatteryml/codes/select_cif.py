import pandas as pd
import glob
import shutil

df =pd.read_csv('../../mofbatteryml/data/csv/porosity_data.csv')

selected_items = df[df['PLD'] > 2.5]['Recode'].to_list()

print (len(selected_items))
all_cifs = glob.glob('../../../../MOF_structures/data/Valid/*cif')

# for cif in all_cifs:
#     basename = cif.split('/')[-1].split('.')[0]
#     if basename in selected_items:
#         shutil.copy (cif, '../data/selected_cifs/')


