#!/usr/bin/env python
# coding: utf-8

# In[1]:


## Import libraries ##
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
import plotly.express as px
import scipy as sci
from scipy.cluster.hierarchy import fcluster
from sklearn.cluster import KMeans
import sklearn.metrics as sm

## Load annotated RNA csv to dataframe ##
df_cm = pd.read_csv('Rice_RNA_FPKM_Data_Annotated.csv')
df_cm = df_cm.drop(['Unnamed: 0'], axis=1)

## Filter data for highest expressing genes and cluster ##
# Generate dataframe for filtering
df_cm_top = df_cm.copy()
loc_top = df_cm_top.pop('loc')
gene_top = df_cm_top.pop('gene')
df_cm_2 = df_cm_top.copy()
scaler = StandardScaler() # Create standard scaler

# Filter for genes expressing at > 2000 FPKM
for i in range(len(df_cm_2)):
    row = df_cm_2.iloc[[i]].values.flatten().tolist()
    if max(row) < 2000:
        df_cm_top.drop([i], inplace=True)
        loc_top.drop([i], inplace=True)
        gene_top.drop([i], inplace=True)
        
# Scale dataframe
df_cm_top_s = scaler.fit_transform(df_cm_top.T)
df_cm_top_s = df_cm_top_s.T

# Cluster by flat cluster distance = 2
dist_mtx = sci.spatial.distance.pdist(df_cm_top_s)
link_mtx = sci.cluster.hierarchy.linkage(df_cm_top_s)
flat_clust = fcluster(link_mtx, t=2.0, criterion='distance').tolist() # Flat cluster with distance

# Set parameters for clustermap
x_labels = ['C_FW1','C_FW2','C_FW3','S_FW1','S_FW2','S_FW3','S_SW1','S_SW2','S_SW3']
y_labels = gene_top.values.flatten().tolist()
lut = dict(zip(set(flat_clust), sns.hls_palette(len(set(flat_clust)), l=0.6, s=0.8))) # Set group colors
row_colors = pd.DataFrame(flat_clust)[0].map(lut) # Map row colors to group

# Plot clustermap of top expressing proteins
cm = sns.clustermap(df_cm_top_s, xticklabels=x_labels, yticklabels=y_labels, 
                    cmap = "vlag", col_cluster=False, row_colors=row_colors.to_numpy())

# Resize clustermap height and width
hm = cm.ax_heatmap.get_position()
plt.setp(cm.ax_heatmap.yaxis.get_majorticklabels(), fontsize=12)
cm.ax_heatmap.set_position([hm.x0, hm.y0, hm.width*2, hm.height])
col = cm.ax_col_dendrogram.get_position()
cm.ax_col_dendrogram.set_position([col.x0, col.y0, col.width*2, col.height*2])

# Save clustermap
plt.savefig('20440_Proj_Top_FPKM_Heatmap.png', bbox_inches='tight')

# Calculate flat cluster silhouette score
sq_dist_mtx = sci.spatial.distance.squareform(dist_mtx)
scc_o = sm.silhouette_score(sq_dist_mtx,flat_clust)
print('Original silhouette score (fclust, distance): ' + str(scc_o))

## Run PCA on data and print results ##
matplotlib.rcdefaults()

# Prepare dataframes
df_cm_t = df_cm.copy()
loc = df_cm_t.pop('loc')
gene = df_cm_t.pop('gene')
df_cm_t = df_cm_t.T
df_cm_t_s = df_cm_t.copy()
df_cm_t['Condition'] = ['Common FW','Common FW','Common FW','SR86 FW','SR86 FW',
                        'SR86 FW','SR86 SW','SR86 SW','SR86 SW']

# Scale dataframe and perform PCA
df_cm_t_s = scaler.fit_transform(df_cm_t_s)
pca = PCA(n_components=9)
pca.fit(df_cm_t_s)
comps = [1,2,3,4,5,6,7,8,9] # List first 9 components
var = pca.explained_variance_ratio_ # Generate variance explained ratios
print('Variance explained by PC1+2: ' + str(var[0]+var[1])) # Add ratios for PC1+PC2

# Plot variance across principal components
fig = plt.figure()
plt.plot(comps,[i*100 for i in var[0:10]]) # PLot variance of first 9 components
plt.xlabel('PC (#)')
plt.ylabel('Variance Explained (%)')
plt.title('Variance Explained by Components')
plt.savefig('20440_Proj_Variance_Explained_by_Components_Plot.png')

# Plot clustering by PCA
pca = PCA()
components = pca.fit_transform(df_cm_t_s)
labels = {str(i): f"PC {i+1} ({var:.1f}%)"
    for i, var in enumerate(pca.explained_variance_ratio_ * 100)}
fig = px.scatter_matrix(components,labels=labels,dimensions=range(6),
    color=df_cm_t['Condition'], width = 1500, height = 1000)
fig.update_traces(diagonal_visible=False)
fig.show()

## Calculate loadings of each gene on PCs ##

# Create loadings dataframe
df_cm_d = df_cm_t.drop(['Condition'], axis=1)
loadings = pd.DataFrame(pca.components_, columns=df_cm_d.columns) # Create loadings dataframe
labels = list(loadings.columns) # Create label list
loadings_t = loadings.copy().T
loadings_t.columns = ['PC1','PC2','PC3','PC4','PC5','PC6','PC7','PC8','PC9']

# Find indices of maximal loading
ind_PC1 = loadings_t['PC1'].idxmax()
ind_PC2 = loadings_t['PC2'].idxmax()

# Plot top gene for PC1 vs PC2
df_cm_t_sc = df_cm_t_s.copy()
df_cm_t_sc = df_cm_t_sc.T
fig = plt.figure()
sns.scatterplot(x=df_cm_t_sc[ind_PC1], y=df_cm_t_sc[ind_PC2], hue=df_cm_t['Condition'])
print('PC1 Top Feature Locus: ' + str(loc.iloc[ind_PC1]))
print('PC2 Top Feature Locus: ' + str(loc.iloc[ind_PC2]))
plt.xlabel(gene.iloc[ind_PC1])
plt.ylabel(gene.iloc[ind_PC2])
plt.title('Top Feature for PC1 and PC2')
plt.savefig('20440_Proj_Top_Feature_for_PC1_and_PC2_Plot.png')
plt.show()

## Generate and analyze dataframes for gene families ##
df_cm['loc'] = loc
df_cm['gene'] = gene
gene_grouped = df_cm.groupby(['gene'])

# Generate dataframe for gene family names
gene_nodup = gene.copy()
gene_nodup = gene_nodup.drop_duplicates()
gene_nodup = gene_nodup.reset_index()
gene_nodup = gene_nodup.drop(['index'], axis=1)

# Generate list of gene family dataframes
gene_dfs = []
for i in range(len(gene_nodup)):
    gene_dfs.append(gene_grouped.get_group(gene_nodup.iloc[[i]].values.flatten()[0]))

## Cluster gene families with flat clustering ##
sscs = []
gene_dfs_i = []
gene_dfs_n = []
dup_ind = [4373,4384,5993] # Families with all duplicate reads
for i in range(len(gene_dfs)):
    # If families are large enough to cluster with no duplicates
    if len(gene_dfs[i]) > 3 and (i not in dup_ind):
        # Generate dataframe for clustering
        gene_ssc = gene_dfs[i].copy()
        gene_dfs_n.append(gene_ssc.iloc[0]['gene']) # Append gene family name
        gene_dfs_i.append(i) # Append index
        gene_ssc = gene_ssc.drop(['loc','gene'], axis=1)
        # Scale dataframe
        gene_ssc = scaler.fit_transform(gene_ssc.T)
        # Calculate distance and linkage matrices
        dist_mtx = sci.spatial.distance.pdist(gene_ssc.T)
        link_mtx = sci.cluster.hierarchy.linkage(gene_ssc.T)
        sq_dist_mtx = sci.spatial.distance.squareform(dist_mtx)
        # Flate cluster with maxclust = 3
        flat_clust = fcluster(link_mtx, t=3, criterion='maxclust').tolist()
        # Calculate silhouette scores for flat clustering and append
        ssc = sm.silhouette_score(sq_dist_mtx,flat_clust)
        sscs.append(ssc)

# Create silhouette score dataframe
ssc_df = pd.DataFrame()
ssc_df['gene_group'] = gene_dfs_i
ssc_df['gene_family'] = gene_dfs_n
ssc_df['silhouette_score'] = sscs
ssc_df = ssc_df.sort_values(by=['silhouette_score'],ascending=False)


## Cluster gene families by K-means ##
kmeans_complete = []
kmeans_sil = []
labels_true = [0,0,0,1,1,1,2,2,2] # True cluster labels for completeness scores
for i in range(len(ssc_df)):
    # Generate dataframe for clustering
    g = int(ssc_df.iloc[i]['gene_group'])
    gene_test = gene_dfs[g].copy()
    gene_test = gene_test.drop(['loc','gene'], axis=1)
    # Scale dataframe
    gene_test_s = scaler.fit_transform(gene_test.T)
    # Cluster by k-means with 3 clusters
    pca = PCA(2) # Define PCA by number of components
    pca_data = pd.DataFrame(pca.fit_transform(gene_test_s),columns=['PC1','PC2']) # Create PCA dataframe
    kmeans = KMeans(n_clusters=3,random_state=42).fit(gene_test_s) # Identify clusters by k-means
    kmeans_pred = KMeans(n_clusters=3,random_state=42)
    kmeans_pred.fit_predict(gene_test_s)
    pca_data['cluster'] = pd.Categorical(kmeans.labels_) # Generate cluster labels
    # Calculate completeness and silhouette scores and append
    kmeans_complete.append(sm.completeness_score(labels_true, kmeans.labels_))
    kmeans_sil.append(sm.silhouette_score(gene_test_s, kmeans_pred.labels_, metric='euclidean'))

# Add k-means scores to silhouette score dataframe
ssc_df['kmeans_complete'] = kmeans_complete
ssc_df['kmeans_sil'] = kmeans_sil

# Calculate average silhouette score between flat and k-means clustering
ssc_df['Avg'] = ssc_df[['silhouette_score', 'kmeans_sil']].mean(axis=1)

## Plot gene family heatmaps with high k-means scores ##

# Set plotting parameters
x_labels = ['C_FW1','C_FW2','C_FW3','S_FW1','S_FW2','S_FW3','S_SW1','S_SW2','S_SW3']
vmin = -2.1
vmax = 2.1

# Plot gene families with high silhouette scores
for i in range(len(ssc_df)):
    # Read silhouette and completeness scores
    s = ssc_df.iloc[i]['kmeans_sil']
    k = ssc_df.iloc[i]['kmeans_complete']
    if  (s >= 0.6) and (k >= 1): # If silhouette and completeness are above standards
        # Print gene group and silhouette score
        g = ssc_df.iloc[i]['gene_group']
        f = ssc_df.iloc[i]['gene_family']
        caption = 'Gene Family_' + str(f) + '_silhouette score_=_' + str(s)
        print(caption)
        # Generate dataframe for clustermap
        gene_test = gene_dfs[g].copy()
        y_labels = gene_test['loc'].values.flatten().tolist()
        gene_test = gene_test.drop(['loc','gene'], axis=1)
        # Scale dataframe
        gene_test_s = scaler.fit_transform(gene_test.T)
        gene_test_s = gene_test_s.T
        # Print clustermap
        cm = sns.clustermap(gene_test_s, xticklabels=x_labels, yticklabels=y_labels,
                            cmap = "vlag", col_cluster=False, vmin=vmin, vmax=vmax)
        plt.setp(cm.ax_heatmap.yaxis.get_majorticklabels(), rotation=0)
        plt.savefig('20440_Proj_' + caption + '_Heatmap.png',bbox_inches='tight')

# Get gene loci for highest expressing proteins
gene_top_df = pd.DataFrame()
gene_top_df['loc'] = loc_top.values.flatten().tolist()
gene_top_df['gene'] = gene_top.values.flatten().tolist()
display(gene_top_df)

# Get gene loci for top ten loadings PC1+2
pc1_ind = list(loadings_t[['PC1']].sort_values('PC1').tail(10).index)
pc2_ind = list(loadings_t[['PC2']].sort_values('PC2').tail(10).index)

pc1_top_genes = [loc.iloc[i] for i in pc1_ind]
pc2_top_genes = [loc.iloc[i] for i in pc2_ind]
pc1_top_genes_fam = [gene.iloc[i] for i in pc1_ind]
pc2_top_genes_fam = [gene.iloc[i] for i in pc2_ind]

print(pc1_top_genes)
print(pc1_top_genes)

