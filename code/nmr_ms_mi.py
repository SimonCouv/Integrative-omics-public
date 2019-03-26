
# setup ----------------------------------------------------------------------------------------------------------------
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
from sklearn import preprocessing
from matplotlib_venn import venn3_unweighted
import missingno as msno
from scipy.spatial import distance
from scipy.cluster import hierarchy
import knnimpute

# prepare data ---------------------------------------------------------------------------------------------------------

prot_ms = pd.read_excel(
    io="/home/simon/OneDrive/KCL/Mayr_proteomics/multi-omics-project/input_data/bruneck/Bruneck 2000 Agilent.xlsx",
    index_col=0, header=0).transpose()
lipid_ms = pd.read_excel(
    io="/home/simon/OneDrive/KCL/Mayr_proteomics/multi-omics-project/input_data/bruneck/Lipidomics + case-control status_withID.xlsx",
    index_col=0, header=0)
cvd = pd.concat([lipid_ms.pop(x) for x in ['CVD0010', 'CVD0010W']], 1)
lipid_nmr = (pd.read_excel(
    io="/home/simon/OneDrive/KCL/Mayr_proteomics/multi-omics-project/input_data/bruneck/Bruneck2000NMRData_Cleaned.xlsx",
    na_values=' ',
    index_col=0, header=0, sheet_name='data')
             .drop(columns=['Isopropyl_', 'glutamin', 'low_protein', 'ethanol','A2MG', 'LCN2', 'SELP', 'PLF4#k', 'PPBP'])
             .drop(index='F3183')) # patient has no NMR data

lipid_nmr_m3 = lipid_nmr.drop(columns=['Pyr', 'Gln', 'Glol'])

lipid_ms.shape
lipid_nmr.shape
prot_ms.shape

# OLD innner join approach

# merged = lipid_nmr.join(lipid_ms, how='inner').join(prot_ms, how='inner')
# merged_complete_cases = merged.dropna(axis='index', how='any')
# merged_complete_feat = merged.dropna(axis='columns', how='any')
#
# merged_m3 = lipid_nmr_m3.join(lipid_ms, how='inner').join(prot_ms, how='inner')
# merged_complete_cases_m3 = merged_m3.dropna(axis='index', how='any')
# merged_complete_feat_m3 = merged_m3.dropna(axis='columns', how='any')
#
# msno.matrix(lipid_nmr)
# msno.matrix(merged)
# msno.matrix(lipid_nmr_m3)
# msno.matrix(merged_m3)
#
#
# merged_complete_cases.shape
# merged_complete_feat.shape
# merged_complete_cases_m3.shape
# merged_complete_feat_m3.shape

# Venn diagrams for correspondance of patient IDs ----------------------------------------------------------------------

plt.figure
venn3_unweighted([set(prot_ms.index.values), set(lipid_nmr.index.values), set(lipid_ms.index.values)],
                 set_labels=('prot_ms', 'lipid_nmr', 'lipid_ms'))
plt.show()

# left-join + impute approach (~ Bhawana)  ----------------------------------------------------------------------------
# merged = lipid_ms.join(lipid_nmr, how="left").join(prot_ms, how="left")
merged = lipid_nmr.join(prot_ms, how="left")  # exclude lipid_ms from analysis
merged.shape
msno.matrix(merged)
prop_missing = merged.isna().mean()

plt.figure()
ax = prop_missing.plot.box(title='feature-wise proportion missing')
ax.hlines(y=0.2,xmin=0,xmax=2, color='r')
# plt.hlines(y=0.2)
plt.show()

# # ECDF
# plt.figure()
# plt.plot(np.sort(prop_missing.values), np.linspace(0, 1, len(prop_missing.values), endpoint=False))
# plt.show()

prop_missing[prop_missing.values > 0.05]
merged = merged[prop_missing.index[prop_missing.values < 0.2]]
# removes Pyr and Glol (Gln remains in at 17.95% missing)
msno.matrix(merged)
missing_gln = merged.Gln.index[merged.Gln.isna()]


# min-max scaling
min_max_scaler = preprocessing.MinMaxScaler()
# test = min_max_scaler.fit_transform(merged)
# test[:, merged.columns.get_loc('Gln')] #na conserved in scaling
merged_minmax = pd.DataFrame(min_max_scaler.fit_transform(merged), columns=merged.columns, index=merged.index)

merged_minmax.GlnW

# imputation
knnimpute.knn_impute_optimistic(merged_minmax.values, k=20, missing_mask=np.isnan(merged_minmax.values))  #in-place!

# export for WGCNA in R
with open('../../data/bruneck_merged_minmax_imp.tsv', mode='w') as f:
    merged_minmax.to_csv(f, sep='\t')
with open('../../data/bruneck_outcomes.tsv', mode='w') as f:
    cvd.to_csv(f, sep='\t')

# correlations ---------------------------------------------------------------------------------------------------------

merged_complete_cases_corr = merged_complete_cases.corr(method='spearman')
merged_complete_cases_m3_corr = merged_complete_cases_m3.corr(method='spearman')
merged_complete_feat_corr = merged_complete_feat.corr(method='spearman')


# Create clustered corplots -------------------------------------------------------------------------------------------

# sns.set_context("notebook", font_scale=0.2)
sns.set_context("notebook", font_scale=1)

# https://seaborn.pydata.org/examples/many_pairwise_correlations.html
# Generate a mask for the upper triangle
# mask = np.zeros_like(prot_corr, dtype=np.bool)
# mask[np.triu_indices_from(mask)] = True

# Generate a custom diverging colormap
cmap = sns.diverging_palette(220, 10, as_cmap=True)

# Generate color code for data source
colors = np.array([('r', 'w', 'w'),
                   ('w', 'g', 'w'),
                   ('w', 'w', 'b')])

colors = np.array(
    [[['r'] * 10, ['w'] * 10, ['w'] * 10],
     (['w'] * 10, ['g'] * 10, ['w'] * 10),
     (['w'] * 10, ['w'] * 10, ['b'] * 10)
     ])

col_width = 3
a=[['r'] * col_width, ['w'] * col_width, ['w'] * col_width]
b=[['w'] * col_width, ['b'] * col_width, ['w'] * col_width]
c=[['w'] * col_width, ['w'] * col_width, ['g'] * col_width]

test = np.concatenate((a,b,c),axis=1)


d = []
for feat in merged.columns.values:
    feat_origin = np.array([feat in lipid_ms.columns.values,
                            feat in lipid_nmr.columns.values,
                            feat in prot_ms.columns.values])
    d.append(test[np.argwhere(feat_origin).tolist()[0][0],])

feat_colors = pd.DataFrame(d, index=merged.columns.values)
feat_colors.columns = [''] * len(feat_colors.columns)
merged_complete_cases_corr_colors = feat_colors.loc[merged_complete_cases_corr.index.values]
merged_complete_cases_m3_corr_colors = feat_colors.loc[merged_complete_cases_m3_corr.index.values]
merged_complete_feat_corr_colors = feat_colors.loc[merged_complete_feat_corr.index.values]

# generate plots
p_merged_complete_features = (sns.clustermap(merged_complete_feat_corr, method="ward", cmap=cmap,center=0,
                                             square=True, cbar_kws={"shrink": .3},
                                             col_colors=merged_complete_feat_corr_colors,
                                             row_colors=merged_complete_feat_corr_colors)
                              .fig.suptitle('Complete features (438 features on 634 cases)'))
p_merged_complete_features.figure.savefig('Complete features (438 features on 634 cases)')
p_merged_complete_cases = (sns.clustermap(merged_complete_cases_corr, method="ward", cmap=cmap,center=0,
                                          square=True, cbar_kws={"shrink": .3},
                                          col_colors=merged_complete_feat_corr_colors,
                                          row_colors=merged_complete_cases_corr_colors)
                           .fig.suptitle('Complete cases (469 features on 245 cases)'))
p_merged_complete_cases.figure.savefig("Complete cases (469 features on 245 cases).png")
p_merged_complete_cases_m3 = (sns.clustermap(merged_complete_cases_m3_corr, method="ward", cmap=cmap,center=0,
                                             square=True, cbar_kws={"shrink": .3},
                                             col_colors=merged_complete_feat_corr_colors,
                                             row_colors=merged_complete_cases_m3_corr_colors)
                              .fig.suptitle('Complete cases, excl Pyr-Gln-Glol (466 features on 613 cases)'))
p_merged_complete_cases_m3.figure.savefig("Complete cases, excl Pyr-Gln-Glol (466 features on 613 cases).png")



# using manually calculated dissimilarity matrix

# complete cases min 3

merged_complete_cases_m3_diss = 1 - merged_complete_cases_m3_corr
merged_complete_cases_m3_condens_diss = distance.squareform(merged_complete_cases_m3_diss)
merged_complete_cases_m3_UPGMAlink = hierarchy.linkage(merged_complete_cases_m3_condens_diss, method='average', optimal_ordering=True)
merged_complete_cases_m3_wardlink = hierarchy.linkage(merged_complete_cases_m3_condens_diss, method='ward', optimal_ordering=True)

p_merged_complete_cases_m3_UPGMAlink = (sns.clustermap(merged_complete_cases_m3_corr, cmap=cmap, center=0,
                                                 square=True, cbar_kws={"shrink": .3},
                                                 row_colors=merged_complete_cases_m3_corr_colors,
                                                 col_colors=merged_complete_cases_m3_corr_colors,
                                                 row_linkage=merged_complete_cases_m3_UPGMAlink,
                                                 col_linkage=merged_complete_cases_m3_UPGMAlink)
                                  .fig.suptitle('Complete cases, excl Pyr-Gln-Glol (466 features on 613 cases) - man UPGMA'))
p_merged_complete_cases_m3_UPGMAlink.figure.savefig('Complete cases, excl Pyr-Gln-Glol (466 features on 613 cases) - man UPGMA')

p_merged_complete_cases_m3_wardlink = (sns.clustermap(merged_complete_cases_m3_corr, cmap=cmap, center=0,
                                                 square=True, cbar_kws={"shrink": .3},
                                                 row_colors=merged_complete_cases_m3_corr_colors,
                                                 col_colors=merged_complete_cases_m3_corr_colors,
                                                 row_linkage=merged_complete_cases_m3_wardlink,
                                                 col_linkage=merged_complete_cases_m3_wardlink)
                                  .fig.suptitle('Complete cases, excl Pyr-Gln-Glol (466 features on 613 cases) - man ward'))
p_merged_complete_cases_m3_wardlink.figure.savefig('Complete cases, excl Pyr-Gln-Glol (466 features on 613 cases) - man ward')


# complete cases

merged_complete_cases_diss = 1 - merged_complete_cases_corr
merged_complete_cases_condens_diss = distance.squareform(merged_complete_cases_diss)
merged_complete_cases_UPGMAlink = hierarchy.linkage(merged_complete_cases_condens_diss, method='average', optimal_ordering=True)
merged_complete_cases_wardlink = hierarchy.linkage(merged_complete_cases_condens_diss, method='ward', optimal_ordering=True)

p_merged_complete_cases_UPGMAlink = (sns.clustermap(merged_complete_cases_corr, cmap=cmap, center=0,
                                              square=True, cbar_kws={"shrink": .3},
                                              row_colors=merged_complete_cases_corr_colors,
                                              col_colors=merged_complete_cases_corr_colors,
                                              row_linkage=merged_complete_cases_UPGMAlink,
                                              col_linkage=merged_complete_cases_UPGMAlink)
                               .fig.suptitle('Complete cases (469 features on 245 cases)- man UPGMA'))
p_merged_complete_cases_UPGMAlink.figure.savefig('Complete cases (469 features on 245 cases)- man UPGMA')

p_merged_complete_cases_wardlink = (sns.clustermap(merged_complete_cases_corr, cmap=cmap, center=0,
                                              square=True, cbar_kws={"shrink": .3},
                                              row_colors=merged_complete_cases_corr_colors,
                                              col_colors=merged_complete_cases_corr_colors,
                                              row_linkage=merged_complete_cases_wardlink,
                                              col_linkage=merged_complete_cases_wardlink)
                               .fig.suptitle('Complete cases (469 features on 245 cases)- man ward'))
p_merged_complete_cases_wardlink.figure.savefig('Complete cases (469 features on 245 cases)- man ward')
#complete feat


merged_complete_feat_diss = 1 - merged_complete_feat_corr
merged_complete_feat_condens_diss = distance.squareform(merged_complete_feat_diss)
merged_complete_feat_UPGMAlink = hierarchy.linkage(merged_complete_feat_condens_diss, method='average', optimal_ordering=True)
merged_complete_feat_wardlink = hierarchy.linkage(merged_complete_feat_condens_diss, method='ward', optimal_ordering=True)

p_merged_complete_feat_UPGMAlink = (sns.clustermap(merged_complete_feat_corr, cmap=cmap, center=0,
                                              square=True, cbar_kws={"shrink": .3},
                                              row_colors=merged_complete_feat_corr_colors,
                                              col_colors=merged_complete_feat_corr_colors,
                                              row_linkage=merged_complete_feat_UPGMAlink,
                                              col_linkage=merged_complete_feat_UPGMAlink)
                               .fig.suptitle('Complete feat (438 features on 634 cases)- man UPGMA'))
p_merged_complete_feat_UPGMAlink.figure.savefig('Complete feat (438 features on 634 cases)- man UPGMA')

p_merged_complete_feat_wardlink = (sns.clustermap(merged_complete_feat_corr, cmap=cmap, center=0,
                                              square=True, cbar_kws={"shrink": .3},
                                              row_colors=merged_complete_feat_corr_colors,
                                              col_colors=merged_complete_feat_corr_colors,
                                              row_linkage=merged_complete_feat_wardlink,
                                              col_linkage=merged_complete_feat_wardlink)
                               .fig.suptitle('Complete feat (438 features on 634 cases)- man ward'))
p_merged_complete_feat_wardlink.figure.savefig('Complete feat (438 features on 634 cases)- man ward')


#alternative linkage methods
# p_prot = sns.clustermap(prot_corr, method="single", cmap=cmap,center=0,
#             square=True, cbar_kws={"shrink": .3})
# p_prot = sns.clustermap(prot_corr, method="weighted", cmap=cmap,center=0,
#             square=True, cbar_kws={"shrink": .3})
# p_prot = sns.clustermap(prot_corr, method="average", cmap=cmap,center=0,
#             square=True, cbar_kws={"shrink": .3})


