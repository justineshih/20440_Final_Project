#!/usr/bin/env python
# coding: utf-8

# In[1]:


## Import libraries ##
import pandas as pd

## Create dataframe for RNA-seq data ##
df = pd.read_csv('GSE102152_FPKM_normalized_expression.txt.gz', sep= '\t') # Read csv to dataframe

# Rename columns to conditions
df = df.rename(columns={'G1':'Common FW1', 'G2':'Common FW2', 'G3':'Common FW3', 
                   'G4':'SR86 FW1', 'G5':'SR86 FW2', 'G6':'SR86 FW3', 'G7':'SR86 SW1', 
                   'G8':'SR86 SW2', 'G9':'SR86 SW3'})

# Drop initial non-gene elements
df.drop(df.index[0:150], inplace=True)
df = df.reset_index()
df = df.drop(['index'], axis=1)

# Drop final non-gene elements
df.drop(df.index[57624:57781], inplace=True)
df = df.reset_index()
df = df.drop(['index'], axis=1)

# Remove and store locus information from dataframe
df = df.sort_values(by=['gene_id'])
df = df.reset_index()
df = df.drop(['index'], axis=1)
gene_id = df.pop('gene_id')
gene_id = gene_id.values.flatten().tolist()


## Open locus annotations and convert to dataframe ##
loc_file = open("all.locus_brief_info.7.0.txt", "r")

# Locus file parameters
Lines = loc_file.readlines()
loc_names = []
gene_names = []
count = 0

# Read each line annotation and save to lists
for line in Lines:
    count += 1
    if (count != 1) and (line.split()[1][0] == 'L'):
        loc_names.append(line.split()[1])
        gene_line = line.split()[9:]
        gene_names.append(' '.join(gene_line))

# Create dataframe of locus annotations
loc_convert_df = pd.DataFrame()
loc_convert_df['loc'] = loc_names
loc_convert_df['gene_name'] = gene_names

# Remove duplicate locus annotations and sort
loc_convert_df = loc_convert_df.drop_duplicates(subset=['loc'])
loc_convert_df = loc_convert_df.sort_values(by=['loc'])
loc_convert_df = loc_convert_df.reset_index()
loc_convert_df = loc_convert_df.drop(['index'], axis=1)
loc_names = loc_convert_df['loc'].values.flatten().tolist()


## Find unnannotated gene loci and save to dataframe ##
drop_ids = list(set(gene_id) - set(loc_names))
drop_ids_df = pd.DataFrame()
drop_ids_df['loc'] = drop_ids
drop_ids_df = drop_ids_df.sort_values(by=['loc'])
drop_ids_df = drop_ids_df.reset_index()
drop_ids_df = drop_ids_df.drop(['index'], axis=1)


## Split RNA-seq dataframe by annotation ##

# Initial dataframes
drop_df = df.copy() # Unannotated dropped
unan_df = df.copy() # Unannoteated retained
gene_id_df = pd.DataFrame() 
gene_id_df['id'] = gene_id

# Save unannotated loci and annotated loci to dataframes
for i in range(len(drop_df)):
    if gene_id[i] not in drop_ids:
        unan_df.drop([i], inplace=True)
    else:
        drop_df.drop([i], inplace=True)
        gene_id_df.drop([i], inplace=True)

# Add loci back to annotated dataframe
drop_df['loc_o'] = gene_id_df['id']
drop_df = drop_df.reset_index()
drop_df = drop_df.drop(['index'], axis=1)

# Add annotation loci and gene name to annotated dataframe
drop_df['loc'] = loc_convert_df['loc']
drop_df['gene'] = loc_convert_df['gene_name']
print_csv = drop_df['loc_o'].equals(drop_df['loc']) # Ensure loci match annotation

# Reformat annotated loci dataframe
drop_df_s = drop_df.sort_values(by=['gene'])
drop_df_s = drop_df_s.drop(['loc_o'], axis=1)
drop_df_s = drop_df_s.reset_index()
drop_df_s = drop_df_s.drop(['index'], axis=1)

# Reformat unattotated loci dataframe and add annotation
unan_df = unan_df.reset_index()
unan_df = unan_df.drop(['index'], axis=1)
unan_df['loc'] = drop_ids_df['loc']
unan_df['gene'] = 'Unannotated'


## Recombine annotated and unannotated dataframes
df_con = pd.concat([drop_df_s, unan_df], ignore_index=True)

# Edit dataframe for downstream analysis
loc = df_con.pop('loc') # Pop loci values
gene = df_con.pop('gene') # Pop gene annotations
df = df_con # Save FPKM values


## Remove non-expressing genes from dataframe ##
df_cm = df.copy()
for i in range(len(df)):
    row = df.iloc[[i]].values.flatten().tolist()
    if max(row) == 0:
        df_cm.drop([i], inplace=True)
        loc.drop([i], inplace=True)
        gene.drop([i], inplace=True)

# Reset indices
df_cm = df_cm.reset_index()
df_cm = df_cm.drop(['index'], axis=1)
loc = loc.reset_index()
loc = loc.drop(['index'], axis=1)
gene = gene.reset_index()
gene = gene.drop(['index'], axis=1)


## Save annotated non-zero gene dataframe to csv ##
df_csv = df_cm.copy()
df_csv['loc'] = loc
df_csv['gene'] = gene
if print_csv:
    df_csv.to_csv('Rice_RNA_FPKM_Data_Annotated.csv')


# In[ ]:




