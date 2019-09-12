
import pickle
from sklearn import cluster
from sklearn import metrics
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import glob
import itertools
from heapq import nsmallest, nlargest
from collections import Counter
import seaborn as sns


vjs = pickle.load(open('155p_feats_of_int.p', 'rb'))


b = np.load('/home/ligia/Desktop/files/TCRBV08-02TCRBJ01-03_index.npy') # 11 seqs
c = np.load('/home/ligia/Desktop/files/TCRBV08-02TCRBJ01-05_index.npy') #5 seqs

d = np.load('/home/ligia/Desktop/files/TCRBV05-03unresolved_index.npy') #38 seqs

e = np.load('/home/ligia/Desktop/files/TCRBV29-or09_02TCRBJ01-01_index.npy') #1

f = np.load('/home/ligia/Desktop/files/TCRBV05-02TCRBJ02-07_index.npy') #2

g= np.load('/home/ligia/Desktop/files/TCRBV22-01TCRBJ02-06_index.npy') #2


bib = np.concatenate([b ,c, d, e, f, g])


filt = vjs.loc[bib]

filt['splt']= filt.sample_tags.apply(lambda row: row.split(','))

filt['splt_hla']= filt.splt.apply(lambda row: list(filter(lambda i: 'MHC class I' in i, row)))


filt['A']= filt.splt_hla.apply(lambda row: tuple(row[0:2]))# if 'HLA-A' in row[1:1] else row[0:1])


#gives u the count of unique HLA (not combo)
cnt = []
for i in splt:
     for s in i:
             if "MHC class I" in s:
                    cnt.append(s)

Counter(cnt)

def plot_confusion_matrix(cm, cmap='magma', labels=None):
    fig = plt.figure(figsize=(20,15))
    ax = fig.add_subplot(111)
    # I also added cmap=cmap here, to make use of the
    # colormap you specify in the function call
    cax = ax.matshow(cm,cmap=cmap)
    fig.colorbar(cax)
    if labels:
        plt.xticks(range(len(labels)), labels)
        plt.yticks(range(len(labels)), labels)
    plt.xticks(rotation=90)
    plt.tight_layout()
    plt.savefig('tri.png')
    plt.close()


plot_confusion_matrix(squareform(cos), cmap='magma', labels=flat_lst)



res = np.flatnonzero(np.isin(name,ix))
ve = np.delete(v, res, axis =0) #array filtered for
                        # the res (array of indexs that you dont want in your final array)






''' get a column with which you can plot all epitopes related to
a pathogen while everything else is tagged as 'other' '''

#this function only works for epitope the one below is more flexible and works for any
# feature you want to plot individually per epitope species

# def get_epitope_per_pathogen(d, pathogen):
#     epitope = d[d['Epitope species'] == str(pathogen)]
#     epi_path = epitope.Epitope.unique().tolist()
#     def ep_val(row):
#         for i in epi_path:
#             if row == i:
#                 return i
#         if row not in epi_path:
#             return 0
#     d['{}_epitope'.format(str(pathogen))] = d['Epitope'].apply(lambda row: ep_val(row))
#
#
#
# for i in d['Epitope species'].unique().tolist():
#     get_epitope_per_pathogen(d, i)


def get_feature_per_pathogen(d, pathogen, feature):
    ep_spec = d[d['Epitope species'] == str(pathogen)]
    feat_path = ep_spec[str(feature)].unique().tolist()
    def feat_val(row):
        for i in feat_path:
            if row[str(feature)] == i and row['Epitope species'] == str(pathogen):
                return i
        else:
            return 0
    d['{}_'.format(str(pathogen))+str(feature)] = d.apply(lambda row: feat_val(row), axis =1)



for i in d['Epitope species'].unique().tolist():
    get_feature_per_pathogen(d, i, 'hla_ag')
    print(i)
    print(d[d['Epitope species'] == str(i)]['aggr_v'].value_counts())
    print(d.groupby(['Epitope species', 'aggr_v'])['aggr_v'].counts())




r[r['aggr_j'] == 1].Epitope.value_counts()

#a list of lists with average cosine sim per j gene per pathogen
F = pd.DataFrame(allmeans:).T
F= F[:-2]
#change INDEX of file to the labels you want to add in "hue"
n = d['Epitope species'].unique().tolist()
n = n[:-2]
F.index = n
avgs = [0.8380632049397447, 0.8421457412706613, 0.8407460008086931, 0.8353110338216968, 0.8411001305665917, 0.8274343351660187, 0.8399568245221785]

DF = pd.concat([F,pd.DataFrame(avgs).T], axis=0)

DF.rename(index={0:'J-Gene Average'},inplace=True)
DF = DF.fillna(DF.mean())


df = pd.melt(DF.reset_index(),id_vars='index', var_name='J-Gene', value_name='Cosine Similarity')
df['Epitope Species'] = df['index']

df['J-Gene'] = df['J-Gene'].apply(lambda i: i+1)



ax = sns.swarmplot(x='J-Gene', y='Cosine Similarity', data=df, hue='Epitope Species')
plt.legend(bbox_to_anchor=(1.04,0), loc= 'lower left')
leg = ax.get_legend()
hl_dict = {handle.get_label(): handle for handle in leg.legendHandles}
hl_dict['J-Gene Average'].set_color('r')
plt.savefig('viol.png', bbox_inches='tight')#, dpi=fig.dpi)
plt.close()








#taken from annot seqs file
''' features'''

d['p']= d['cmv_newlabel'].apply(lambda i:1 if i == "Positive" else 0)


d['len'] = d['amino_acid'].apply(lambda i: len(i))
d[d['Epitope species'] == 'HomoSapiens'].Epitope.value_counts().count()

d[d['Epitope species'] == 'HomoSapiens'].Epitope.value_counts().count()
d[d['aggr_j'] == 1].Epitope.value_counts()

eeee = r[r['Epitope species'] == 'HIV-1'][r[r['Epitope species'] == 'HIV-1']['aggr_j'] == 7][['Epitope', 0, 1]].sort_values(1)
# v = np.load('/home/ligia/Desktop/Emerson2017/vec_weights100D_allP_vdjdb.npy')
# d = pd.read_csv('/home/ligia/Desktop/Emerson2017/vdjdb_notsorted.csv')


# generate column of disease type features from additional dataset on CDR3 annotations
# def gen_table(self):
#
#     temp = pd.DataFrame()
#     temp[['amino_acid', 'type', 'V', 'J','MHCA', 'MHCB']] = data[['CDR3', 'Epitope species', 'V', 'J', 'MHC A', 'MHC B']]
#
# #    temp['type'] = data['Epitope species']
#
#     c = temp['amino_acid'].isin(df_features['amino_acid'].astype(str))
#
#     print(c.value_counts()) #we see that most of these CDRs are not in our sample (only 9 of them /125)
#
#
#     #select the annot seqs that are in sample
#
#     pres_seqs = temp[c]
#
#     return pres_seqs
# z = gen_table(seqs)
# pres_seqs.drop_duplicates('amino_acid', keep='first')
# %%
# #generate stats
# duplicates = z.duplicated(subset=['amino_acid'], keep='first')
# z[duplicates].sort_values(by='amino_acid')
#
# z = z[~duplicates]
# #%%
# df_features = df_features.merge(z, how='outer')
#
# #%%
#
# df_features['type'] = df_features['type'].replace('nan', np.nan).fillna(0)
#
# #%%
#
# #check if these seqs are public or private
# mask = np.in1d(df_features['type'].values,['InfluenzaA', 'CMV', 'YellowFeverVirus', 'HIV-1', 'HomoSapiens',
#        'HCV', 'EBV', 'DENV1', 'DENV3/4', 'RSV', 'LCMV', 'MCMV', 'HTLV-1',
#        'SIV', 'DENV2', 'HSV-2', 'GallusGallus', 'TriticumAestivum'])
#
# table = df_features[['sample_name', 'type', 'amino_acid', 'J', 'V']][mask]


# %%
embedding_matrix = np.load('/home/ligia/Desktop/Emerson2017/vec_weights100D_vdjdb.npy')
d = pd.DataFrame(embedding_matrix)
j = pd.concat([x['single'], d], axis=1)

b2 = j[j['single'] == 'B2']
b2 = np.array(b2.drop('single', axis=1))
x_b2 = x[x['single'] == 'B2']
d = x_b2.reset_index(drop=True)
# %%

d[d['Epitope species'] == 'EBV']['J'].value_counts()


# %%
# aggregate J genes for VDJDB
def aggr_j(row):
    if "B" in row:
        return "B{0}".format(row[6])
    if "A" in row:
        return row[0:6]
    else:
        return "other"


#
#
beta['aggr_j'] = beta['J'].apply(lambda row: aggr_j(row))
beta['single'] = beta['aggr_j'].apply(lambda row: spec_val(row))
#
x_b1 = d[d['single'] == 'B1']
# %%

''' SAVE NPY FILE OF DATA FILTERED BY GENE OF INTEREST '''


def ge_file(gene, file, single=True, dbl=False):
    embedding_matrix = np.load(
        'vec_weights100D_from300p_vdjdb.npy')  # these vec weights have been reconstructed from OG sample to VDJDB
    e = pd.DataFrame(embedding_matrix)
    d = beta[beta['single'].str.contains(str(gene))].reset_index(drop=True)
    j = pd.concat([beta['single'], e],
                  axis=1)  # x already contains the column single which the B genes have been aggregated (both family 1 and 2)
    if single:
        j_genes = j[j['single'] == str(gene)]
    if dbl:
        j_genes = j[j['single'].str.contains(str(gene))]  # format of gene 'd|f'    #when selecting 2 J genes to view
    j_genes = np.array(j_genes.drop('single', axis=1))
    return d, np.save(str(file), j_genes)


d = ge_file('B1|B7', 'b1b7_300p', single=False, dbl=True)


# %%

def ge_file(antigen, file):
    embedding_matrix = np.load('vec_weights100D_from300p_vdjdb.npy')
    e = pd.DataFrame(embedding_matrix)
    tempbeta = beta.sort_values(['aggr_j', 'V'])
    d = tempbeta[tempbeta['Epitope species'].str.contains(str(antigen))].reset_index(drop=True)
    j = pd.concat([beta[['Epitope species', 'aggr_j', 'V']], e],
                  axis=1)  # x already contains the column single which the B genes have been aggregated (both family 1 and 2)
    j = j.sort_values(['J', 'V'])
    ant = j[j['Epitope species'] == str(antigen)]
    ant = np.array(ant.drop(['Epitope species', 'J', 'V'], axis=1))
    return d, tempbeta, np.save(str(file), ant)


d = ge_file('EBV', 'EBV_300p')

# %%

'''plot cosine similarity '''

import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd

plt.figure(figsize=(15, 15))
plt.matshow(np.load('cos.npy'), fignum=1)
plt.colorbar()
plt.savefig('d.png')
# %%

'''PLOT ONLY V GENES ASSOCIATED WITH SPECIFC ANTIGEN // what v gene corresponds wit a specif antigen sequence'''


def get_v(colname, antigen):
    d[str(colname)] = d.loc[antigen == 1, 'V']
    d[str(colname)] = d[str(colname)].replace(np.nan, 0)
    print(d[str(colname)].value_counts())
    return d[str(colname)]


get_v('vHTLV1', d.HTLV1)

# %%
''' PLOT A SINGLE VALUE ON TSNE '''


def spec_val(row):
    if row == 'HTLV-1':  # or row == 'TRAJ11':
        return 1
    else:
        return 0


d['HTLV1'] = d['Epitope species'].apply(lambda row: spec_val(row))


# %%
def spec_val(row):
    if row == 'B1':
        return 'B1'
    if row == 'B2':
        return 'B2'
    if row == 'B3':
        return 'B3'
    if row == 'B4':
        return 'B4'
    if row == 'B5':
        return 'B5'
    if row == 'B6':
        return 'B6'
    if row == 'B7':
        return 'B7'
    if row == 'other':
        return 'NA'
    else:
        return "Alpha"







'''FROM LEV_DIS_COMPUTE // it works'''
arr_labels = d.sort_values('aggr_j')['aggr_j']

lcm = dict(zip(set(arr_labels), sns.color_palette('PuOr',
                                                      len(set(arr_labels)))))
rc = pd.DataFrame(arr_labels)['aggr_j'].map(lcm)
#sqmat =squareform(np.load('distancematrix_vdjdb.npy'))
g = sns.clustermap(squareform(distance_matrix), row_cluster=False,
               col_cluster=False, linewidth=0,cmap='magma',
               row_colors=[rc], col_colors=[rc])
for al in set(arr_labels):
    g.ax_col_dendrogram.bar(0,0,color=lcm[al], label=al, linewidth=0)

g.ax_col_dendrogram.legend(loc='center', ncol=7)
g.cax.set_position([.15,.2,.03,.45])
# plt.matshow(squareform(distance_matrix), fignum=1, cmap = 'magma')
# plt.xticks(rotation=90)
# plt.colorbar()
#plt.matshow(squareform(distance_matrix), fignum=1, cmap = 'magma')
#ax = g.ax_heatmap
plt.axis('off')
plt.savefig('vdjaggr_jo_mat.png', bbox_inches='tight')
plt.close()

# ss_epi = ss.sort_values('hla_ag')




''' FROM FISH ANALYSIS'''

# # pval_eth = []
# # cn_eth = []
# xt = []
# # for ss in ze:
# #     xt = []
# #     p_values = []
# #     sig = []
# #     combi_name = []
# #     fail=[]
# for unq_j in bin_j:
#     for vars in varlist2:
#         for v in vars:
#             # try:
#             xtab = pd.crosstab(d[unq_j], d[v], margins = False)
#             if xtab.shape!= (2,2):
#                  continue
#             xt.append(xtab)
#             combi_name.append(tuple([unq_j, e]))
#             r, p = stats.fisher_exact(xtab, alternative='two-sided')
#             if p <= .05:
#                 sig.append(xtab)
#             p_values.append(p)
# pval_eth.append(p_values)
# cn_eth.append(combi_name)
# xt_eth.append(xt)
# #
# # res_per_ss = []
# #
# # # ze = dfs[:2]
# # for ss in dfs:
# #     res_per_ss.append(get_fish(ss, unique_cols, vs))



#Testing the code
pval_eth = []
cn_eth = []
xt_eth = []
counting = 0
fail_eth = []
for ss in ze:
    xt = []
    p_values = []
    sig = []
    combi_name = []
    fail = []
    counting += 1
    for v1 in bin_v:
        for var2 in varlist:
            for v2 in var2:
                xtab = pd.crosstab(ss[v1], ss[v2], margins=False)
                if xtab.shape != (2, 2):
                    fail.append(tuple([str(counting) + v1, v2]))
                    continue
                xt.append(xtab)
                combi_name.append(tuple([str(counting) + v1, v2]))
                r, p = stats.fisher_exact(xtab, alternative='two-sided')
                if p <= .05:
                    sig.append(xtab)
                p_values.append(p)
    pval_eth.append(p_values)
    cn_eth.append(combi_name)
    xt_eth.append(xt)
    fail_eth.append(fail)
    # except:
    #     fail.append(combi_name)

pvalues = []
for v3 in bin_j:
    for vars in varlist2:
        for v in vars:
            # for v2 in var2:
            #     # try
            xtab = pd.crosstab(d[v3], d[v], margins=False)
            if xtab.shape != (2, 2):
                fail_eth.append(tuple([v3, v]))
                continue
            xt_eth.append(xtab)
            cn_eth.append(tuple([v3, v]))
            r, p = stats.fisher_exact(xtab, alternative='two-sided')
            #     if p <= .05:
            #         sig.append(xtab)
            #     p_values.append(p)
            pvalues.append(p)

pval_eth.append(pvalues)





#
# def get_fish(d, epi, js):
#     xt = []
#     p_values = []
#     sig = []
#     combi_name = []
#     fail=[]
#     for unq_j in js:
#         for e in epi:
#             # try:
#             xtab = pd.crosstab(d[unq_j], d[e],margins = False)
#             if xtab.shape!= (2,2):
#                 fail.append(tuple([unq_j,e]))
#                 continue
#             xt.append(xtab)
#             combi_name.append(tuple([unq_j,e]))
#             r, p = stats.fisher_exact(xtab, alternative='two-sided')
#             if p <= .05:
#                 sig.append(xtab)
#             p_values.append(p)
#             # except:
#             #     fail.append(combi_name)
#     r, corr, sidak, bonfer = multi.multipletests(p_values,.05, method='fdr_bh')
#     rejected_xtabs_ix = np.where(r == True)
#     rej_xtabs = [xt[i] for i in rejected_xtabs_ix[0]]
#     combi_name =[combi_name[i] for i in rejected_xtabs_ix[0]]
#     #rej[0][0][0] -->first column top element / switch last 0 to 1 you get bottom
#     #rej[0][1][0] --> second column top element / switch last 0 to 1 you get bottom
#     enriched_count=0
#     nrch = []
#     for rej, cn in zip(rej_xtabs, combi_name):
#         try:
#             ath = rej[1][0]/ rej[0][0]
#             j = np.float(rej[1][1]) / np.float(rej[0][1])
#             print('Other:{}'.format(ath) , 'J:{}'.format(j))
#             if j>ath:
#                 print('ENRICHED!')
#                 enriched_count += 1
#                 nrch.append(cn)
#             else:
#                 print('NOT ENRICHED')
#         except:
#             print(cn)
#     return rej_xtabs, p_values, r, nrch, enriched_count, combi_name, fail
#
#


#WORKS // old way of doing the analysis (was only J and VJ now its J and J + all other feature values)
def get_fish(dfs, varlist, varlist2, var1, var3):
    pval_eth = []
    cn_eth = []
    xt_eth = []
    counting=0
    fail_eth=[]
    for ss in dfs:
        xt = []
        p_values = []
        sig = []
        combi_name = []
        fail = []
        counting+=1
        for v1 in var1:
            for var2 in varlist:
                for v2 in var2:
                    xtab = pd.crosstab(ss[v1], ss[v2], margins=False)
                    if xtab.shape != (2, 2):
                        fail.append(tuple([str(counting) + v1, v2]))
                        continue
                    xt.append(xtab)
                    combi_name.append(tuple([str(counting) + v1, v2]))
                    r, p = stats.fisher_exact(xtab, alternative='two-sided')
                    if p <= .05:
                        sig.append(xtab)
                    p_values.append(p)
        pval_eth.append(p_values)
        cn_eth.append(combi_name)
        xt_eth.append(xt)
        fail_eth.append(fail)
            # except:
            #     fail.append(combi_name)
    j_pvalues = []
    j_combiname = []
    j_fail = []
    j_xtabs = []
    for v3 in var3:
        for vars in varlist2:
            for v in vars:
            # for v2 in var2:
            #     # try
                xtab = pd.crosstab(d[v3], d[v], margins = False)
                if xtab.shape!= (2,2):
                    j_fail.append(tuple([v3, v]))
                    continue
                j_xtabs.append(xtab)
                j_combiname.append(tuple([v3, v]))
                r, p = stats.fisher_exact(xtab, alternative='two-sided')
            #     if p <= .05:
            #         sig.append(xtab)
            #     p_values.append(p)
                j_pvalues.append(p)
    pval_eth.append(j_pvalues)
    cn_eth.append(j_combiname)
    xt_eth.append(j_xtabs)
    fail_eth.append(j_fail)
    pval_eth= list(itertools.chain(*pval_eth))
    cn_eth = list(itertools.chain(*cn_eth))
    xt_eth = list(itertools.chain(*xt_eth))
    fail_eth = list(itertools.chain(*fail_eth))
    r, corr, sidak, bonfer = multi.multipletests(pval_eth,.05, method='fdr_bh')
    rejected_xtabs_ix = np.where(r == True)
    rej_xtabs = [xt_eth[i] for i in rejected_xtabs_ix[0]]
    combi_name =[cn_eth[i] for i in rejected_xtabs_ix[0]]
    #rej[0][0][0] -->first column top element / switch last 0 to 1 you get bottom
    #rej[0][1][0] --> second column top element / switch last 0 to 1 you get bottom
    enriched_count=0
    nrch = []
    for rej, cn in zip(rej_xtabs, combi_name):
        try:
            ath = rej[1][0]/ rej[0][0]
            j = rej[1][1] / rej[0][1]
            print('Other:{}'.format(ath) , 'J:{}'.format(j))
            if j>ath:
                print('ENRICHED!')
                enriched_count += 1
                nrch.append(cn)
            else:
                print('NOT ENRICHED')
        except:
            print(cn)
    return rej_xtabs, pval_eth, r, nrch, enriched_count, combi_name, fail_eth, sidak, bonfer




'''runnin_cossimv2'''

#   clrrect cmv
ext_cmv = np.load('164cmvgenerator.npy')
#v = np.load('vec_weights100D_155p_tr330_unique.npy')
ext_cmv = [vec.reshape(1, -1) for vec in ext_cmv]

ext_cmv = np.load('150inflgenerator.npy')
#v = np.load('vec_weights100D_155p_tr330_unique.npy')
ext_cmv = [vec.reshape(1, -1) for vec in ext_cmv]


## commented is to take a subsample to test code
#ext_cmv_index = np.load('164cmvgenerator_index.npy')
#
#shuffle(ext_cmv_index)

#ix = ext_cmv_index[:10]
#np.save('ix_164_10samples', ix)

#ix = np.load('ix_164_10samples.npy')
#ix = ix[:2]
#print(ix)

ext_cmv = ve[np.random.choice(ve.shape[0], size= 164), :]
ext_cmv = [vec.reshape(1, -1) for vec in ext_cmv]

ix = np.load('164cmvgenerator_index.npy')
#ix = np.load('150inflgenerator_index.npy')


#ext_cmv = []
#for i in ix:
#    ext_cmv.append(vecs[i])
#
#ext_cmv = [vec.reshape(1, -1) for vec in ext_cmv]





import numpy as np
from sklearn.metrics.pairwise import cosine_similarity


class top_k_analyses():
    def __init__(self, ext_cmv, all_vj_files):
        self.maxlen = 500
        self.ext_cmv = ext_cmv
        self.all_vj_files = all_vj_files
    def generator(self):
        '''
        :return: a generator for each embedding we want to compute
                cosine similarity across all 844 combination
                36 cores and 164 embeddings takes 24h (bottleneck is insort algorithm)
        '''
        for i in self.ext_cmv:
            yield i
    def insort(self, a, x, lo=0, hi=None):
        '''
        disclaimer: this algorithm is taken from Python's standard library
        the function is "bisect" but we modify to include a maximum len
        :param a: an array (ex. topk)
        :param x: the cosine similarity you want to put in the array
        :return: insert cosine similarity into a sorted array
                maintain the length of the array at 500 cosine similarities
        '''
        if lo < 0:
            raise ValueError('lo must be non-negative')
        hi = len(a)
        while lo < hi:
            mid = (lo + hi) // 2
            if x < a[mid]:
                hi = mid
            else:
                lo = mid + 1
        if len(a) < self.maxlen:
            a.insert(lo, x)
    def get_topk(self, seq_of_int):
        '''
        :param seq_of_int: each i from the generator
        :return: for all seq_of_int a list of the average of the top 500 cosine
                similarities per VJ combination
        '''
        all_files = []
        for vj_file in self.all_vj_files:
            topk = []
            for i in vj_file:
                cos = cosine_similarity(seq_of_int, i)
                # topk.append(cos)
                # bisect.insort(topk, cos)
                self.insort(topk, cos)
                #vec_index.append(i)
                #cmv_seqs.append(seq_of_int)
                #print(seq_of_int)
            #topk=topk
            topk = np.mean(topk)
            all_files.append(topk)
            #cs.append(cmv_seqs)
        return all_files




'''COSINE ANALYSIS FILE'''




'''ALTERNATIVE OF PLOTTING COSINE FOR SEQUENCE ACROSS VDJDB
 WHERE U ONLY NEED TO PLOT 3 AT THE SAME TIME AND NOT ALL OF THEM'''
def singlePlot(dict):
    fig, ax = plt.subplots(3, 1, figsize=(10,8))
    y1=dict[1]
    n1 = len(y1)
    x1 = np.linspace(1,n1,n1)
    y2=dict[20]
    n2 = len(y2)
    x2 = np.linspace(1,n2,n2)
    y3=dict[100]
    n3 = len(y3)
    x3 = np.linspace(1,n3,n3)
    ax[0].plot(x1, y1,color='#221631')
    # ax[0].set_xlabel('
    ax[1].plot(x2,y2, color='#221631')
    ax[2].plot(x3,y3, color ='#221631')
    plt.savefig('cosrandom_vj.png')
    plt.close()


b = np.load('/home/ligia/Desktop/files/TCRBV08-02TCRBJ01-03_index.npy') # 11 seqs
c = np.load('/home/ligia/Desktop/files/TCRBV08-02TCRBJ01-05_index.npy') #5 seqs

d = np.load('/home/ligia/Desktop/files/TCRBV05-03unresolved_index.npy') #38 seqs

e = np.load('/home/ligia/Desktop/files/TCRBV29-or09_02TCRBJ01-01_index.npy') #1

f = np.load('/home/ligia/Desktop/files/TCRBV05-02TCRBJ02-07_index.npy') #2

g= np.load('/home/ligia/Desktop/files/TCRBV22-01TCRBJ02-06_index.npy') #2


bib = np.concatenate([b ,c, d, e, f, g])


filt = vjs.loc[bib]

filt['sample_name'].value_counts()

filt['sta'].value_counts()


flat_rand = [val for inner in rand for val in inner]

np.save('flat_rand', flat_rand)

t.test(cos, rand, paired = TRUE, alternative = "greater")





# #%%flat_lst = [val for inner in lwstn_all for val in inner]
#
#
# # #OLD WAY OF DOING IT, NEW WAY IS BELOW
# # def get_count(lst, lrgst=False, lwst=False):
# #     flat_lst = [val for inner in lst for val in inner]
# #     counter = Counter(flat_lst)
# #     most_occur_ix = counter.most_common()
# #     most_occur_vj = []
# #     for i in most_occur_ix[:100]:
# #         most_occur_vj.append(seq_idx_allvjs[i[0]])
# #     splt = [vj.split('TCR') for vj in most_occur_vj]
# #     for i in splt:
# #         i.remove('')
# #     unresolved = []
# #     v = []
# #     j = []
# #     for i in splt:
# #         if len(i) == 2:
# #             v.append(i[0][2:4])
# #             j.append(i[1][-2:])
# #         if len(i) == 1:
# #             unresolved.append(i)
# #     return flat_lst, most_occur_ix, most_occur_vj, splt, v, j, unresolved
# #
# #
# # flat_lst, most_occur_ix, most_occur_vj, splt, v, j, unresolved = get_count(lrgst_all, lrgst=False, lwst=False)
#
# # flat_lst_rand = [val for inner in lw_cos_all for val in inner]
#
# # #manual dicitonary i created to retrieve the n lowest VJ combs for plotting
# # dic = {398:'NA-V09',370:'J03-V08', 486:'J05-V08', 402:'NA-V05', 371:'J01-V29orV09', 239:'j07-V05', 607:'J06-V22'}
# #
# # flat_lst = [dic[i] for i in flat_lst]
#
# # most_occur_vj = []
# # for i in flat_lst:
# #     name = seq_idx_allvjs[i][38:-4]
# #     name.split()
#
# # vj_splt =  [vj.split('TCR') for vj in most_occur_vj]
# #PLOTTING N LOWEST VJ COMBS
# fig, ax = plt.subplots()
# sns.countplot(pd.DataFrame(flat_lst)[0], color='#221631',
#               order=pd.DataFrame(flat_lst)[0].value_counts().sort_values().index)
# ax.set(ylabel="Count", xlabel= 'J V Combination')
# plt.xticks(rotation=45)
# # a=ax.get_xticks().tolist()
# # a[1]='change'
# # ax.set_xticklabels(a)
# plt.tight_layout()
# plt.savefig('164_lowestvjcombs.png')
# plt.close()



#for the largest VJ comb cosines
flat_lst = [val for inner in lrgst_all for val in inner]
counter = Counter(flat_lst)
most_occur_ix = counter.most_common()
most_occur_vj = []
for i in most_occur_ix[:100]:
    most_occur_vj.append(seq_idx_allvjs[i[0]][38:-4])




# cmv['v'] = cmv['V gene'].apply(lambda row: row[5:7])

# pd.crosstab(df.v, df.j)
# df.pivot_table(index='j', columns='v', values='cnt'.sum())
#analysis of each Vj on its own
# tab = df.groupby(['j', 'v'])['cnt'].sum().unstack().fillna(0)
tab = df.groupby(['j', 'v'])['cnt'].sum()
tab_percent = tab.groupby(level = 0).apply(lambda x: round(100*x / float(x.sum()))).unstack().fillna(0)
cmap = sns.cubehelix_palette(light=1, as_cmap=True)
plt.subplots(figsize=(14,5))
ax_hmj = sns.heatmap(tab_percent, annot=True, fmt="g", cmap= cmap, xticklabels=True, yticklabels=True)
ax_hmj.set(ylabel="J-Gene", xlabel= 'V-Family')
ax_hmj.tick_params(left=False, bottom=False)
plt.savefig('jgene_vfam_heatmap_rand_lrgst_perc.png', bbox_inches='tight')
plt.close()


# all_split = [[i.split('unresolved') for i in vj] for vj in splt]



# ''' this is the old way of seprating v n js from a tuple'''
# sep = []
# n = 12
# for vj in most_occur_vj:
#     sep.append(tuple(vj[i:i + n] for i in range(0, len(vj), n)))
#
# def get_count(list):
#     counter = Counter(list)
#     most_occur = counter.most_common()
#     return most_occur
#
# v_count_influenza = get_count(v)
# j_count_influenza = get_count(j)
#
#
# v_count_cmv, j_count_cmv  = get_count(v), get_count(j)
#



'''DATA VISUALIZATION'''

'''PLOT DENDROGRAM W/ CLUSTERMAP'''


label_1 = np.array(d['len'])
label_colormap_v = dict(zip(set(label_1), sns.color_palette('PuOr', len(set(label_1)))))
row_colors_v = pd.DataFrame(label_1)[0].map(label_colormap_v)



g=sns.clustermap(squareform(cos), row_cluster=False, col_cluster=False, linewidths=0, cmap='magma',
                 row_colors=[row_colors_v], col_colors=[row_colors_v])


for l_v in set(label_1):
    g.ax_col_dendrogram.bar(0, 0, color=label_colormap_v[l_v], label=l_v, linewidth=0)

g.ax_col_dendrogram.legend(loc="center", ncol=9)

g.cax.set_position([.15,.2, .03, .45])
ax = g.ax_heatmap
ax.axis('off')
plt.savefig('clust_vdjdb_sortbylen.png')
plt.close()


plt.matshow(squareform(cos), fignum=1, cmap = 'magma')
plt.axis('off')
plt.savefig('vdjdbbylen.png')
plt.close()

''' WHEN YOU WANT TO PLOT 2 LABELS AT ONCE'''
label_1 = np.array(d['len'])
label_colormap_v = dict(zip(set(label_1), sns.color_palette('PuOr', len(set(label_1)))))
row_colors_v = pd.DataFrame(label_1)[0].map(label_colormap_v)

#Create additional row_colors here
label_j = np.array(j1['Epitope species'])
label_colormap_j = dict(zip(set(label_j), sns.color_palette('tab20c', len(set(label_j)))))
row_colors_j = pd.DataFrame(label_j)[0].map(label_colormap_j)


# row_linkage, col_linkage = (hc.linkage(sp.distance.pdist(v, 'cosine'), method='ward')
#   for x in (v, v.T))

# linkage = hierarchy.linkage(vr, method='average')
# g=sns.clustermap(cos, row_cluster=False, col_cluster=False, linewidths=0, cmap='magma',
#                  row_colors=[row_colors_j, row_colors_v], col_colors=[row_colors_j, row_colors_v])

##old version before 14aug when you were trying out the method with CMV sequences
# g=sns.clustermap(cos, row_cluster=False, col_cluster=False, linewidths=0, cmap='magma',
#                  row_colors=[row_colors_j, row_colors_v], col_colors=[row_colors_j, row_colors_v])
#

g=sns.clustermap(squareform(cos), linewidths=0, cmap='magma',
                 row_colors=[row_colors_j, row_colors_v], col_colors=[row_colors_j, row_colors_v])


for l_v, l_j in zip(set(label_1), set(label_j)):
    g.ax_col_dendrogram.bar(0, 0, color=label_colormap_j[l_j], label=l_j, linewidth=0)
    g.ax_col_dendrogram.bar(0, 0, color=label_colormap_v[l_v], label=l_v, linewidth=0)

g.ax_col_dendrogram.legend(loc="center", ncol=9)


g.cax.set_position([.15,.2, .03, .45])
plt.savefig('hclust_j7.png')



all_zips = os.listdir("/home/ligia/Desktop/Emerson2017/zip_trial2/")

main = ef.create_main_file(all_zips)


#count sequences with stop codons
main_150['amino_acid'].str.contains('\*').value_counts()[True]


#drop sequences with stop codons
main_150 = main_150[~main_150.amino_acid.str.contains('\*')]

'''GENERATOR'''

#   Create dicitionary of sequence name and index

x = np.load('330_files_unique_155seqname.npy')

arr_fl = x.ravel()

seq_idx_dic = {i: arr_fl[i] for i in range(0, len(arr_fl))}

#   Load vectors & reshape them to compute top k cosine similarity

v = np.load('vec_weights100D_155p_tr330_unique.npy')

vecs = [vec.reshape(1, -1) for vec in v]

#   Create generator for sequences of interest
def gen_prep(cmv=True, antigen=True):#, hla=True, age=True):
    feature_unique_all = []
    if cmv:
        ext_seqs = pd.read_csv('164cmv.csv')
        temp = pd.DataFrame()
        temp['amino_acid'] = ext_seqs['CDR3']
        c = temp['amino_acid'].isin(arr_fl.astype(str))
        kn_seqs = temp[c]
        feature_unique_all.append(np.unique(kn_seqs['amino_acid'].values))
    if antigen:
#   find which sequence match sequences from VDJDB with known antigen association --> below for ALL antigens in VDJDB
        ext_seqs = pd.read_csv('SearchTable-2019-05-28 21_49_18.859.tsv', sep='\t')
        ext_seqs = ext_seqs[ext_seqs['Species'] == 'HomoSapiens']
        temp = pd.DataFrame()
        temp[['amino_acid', 'type']] = ext_seqs[['CDR3', 'Epitope species']]
        c = temp['amino_acid'].isin(arr_fl.astype(str))
        kn_seqs = temp[c]
        feature = kn_seqs[kn_seqs['type'] == 'CMV']
        feature_unique_all.append(np.unique(feature['amino_acid'].values))
    # if hla:
    #
    # if age:
    index_list = []
    for feat in feature_unique_all:
        temp1 = []
        for seq in feat:
            temp1.append(list(seq_idx_dic.keys())[list(seq_idx_dic.values()).index(str(seq))])
        index_list.append(temp1)
    gen_seqs = []
    for feat in index_list:
        temp2 = []
        for i in feat:
            gen_seqs.append(vecs[i])
        gen_seqs.append(temp2)
    return gen_seqs, index_list


gen_seqs, index_list = gen_prep(cmv=True, antigen=True)

# s = gen_seqs[:10]
# np.save('164cmvgenerator', gen_seqs)
#
# np.save('164cmvgenerator_index', index_list)

# def generator():
#     for i in s:
#         yield i
#
#
# maxlen = 100
#
#
# def insort(a, x, lo=0, hi=None):
#     if lo < 0:
#         raise ValueError('lo must be non-negative')
#
#     hi = len(a)
#     while lo < hi:
#         mid = (lo + hi) // 2
#         if x < a[mid]:
#             hi = mid
#         else:
#             lo = mid + 1
#     if len(a) < maxlen:
#         a.insert(lo, x)
#
#
# def get_topk(seq_of_int):
#     topk = []
#     vec_index = []
#     for i in vecs:
#         cos = cosine_similarity(seq_of_int, i)
#         # topk.append(cos)
#
#         # bisect.insort(topk, cos)
#         insort(topk, cos)
#         vec_index.append(i)
#
#     return topk, vec_index
#
#
# pool = mp.Pool(processes=10)
#
# print('generating top_k cosine similarity')
# results = pool.map(get_topk, generator())
#
# np.save('cosresults', results)
# print('Done.')

#   find which sequence match sequences from VDJDB with known antigen association --> below for ALL antigens in VDJDB

data = pd.read_csv('SearchTable-2019-05-28 21_49_18.859.tsv', sep='\t')
data = ext_seqs[ext_seqs['Species'] == 'HomoSapiens']
temp = pd.DataFrame()
temp[['amino_acid', 'type']] = ext_seqs[['CDR3', 'Epitope species']]
c = data['CDR3'].isin(arr_fl.astype(str))
kn_seqs = data[c]
feature = kn_seqs[kn_seqs['type'] == 'InfluenzaA']
feature = np.unique(feature['amino_acid'].values)

from random import shuffle

shuffle(feature)
feat = feature[:150]

index_list = []
for seq in feat:
    index_list.append(list(seq_idx_dic.keys())[list(seq_idx_dic.values()).index(str(seq))])

gen_seqs = []
for ix in index_list:
    gen_seqs.append(vecs[ix])
np.save('150inflgenerator', gen_seqs)

np.save('150inflgenerator_index', index_list)




'''VJ MATRIX FILES'''


all_zips = os.listdir("/home/ligia/Desktop/Emerson2017/temp/")

main = create_main_file(all_zips)

print('prepro done')

# np.save('vj_comb', main)

filename = 'vj_comb.npy'
vec_w, seq_idx_dic, vj_as_int, vj_unique_str = get_dicts(filename)


# def manual_func(data):
#     lookup = {}
#     result = []
#     #z = []
#     for ix, element in enumerate(data):
#         if element not in lookup:
#             target = lookup[element] = [element]
#             result.append(target)
#             #print(ix)
#         else:
#             lookup[element].append(element) #this works
#             #z.append(ix)
#     return result


# use this with cat_as_int
def get_ixdict(data):
    vj_ix_dic = {}
    result = []
    #z = []
    for ix, element in enumerate(data):
        if element not in vj_ix_dic:
            target = vj_ix_dic[element] = [ix]
            result.append(target)
            #print(ix)
        else:
            vj_ix_dic[element].append(ix) #this works
            #z.append(ix)
    val_only = [v for v in vj_ix_dic.values()]
    key_only = [k for k, v in vj_ix_dic.items()]

    return vj_ix_dic, val_only, key_only


vj_ix_dic, val_only, key_only = get_ixdict(vj_as_int)



#   creating a list of lists where each inner list corresponds to embeddings grouped by vj combo
#res_list = list(map(v.__getitem__, only_values))
grp_by_vj = [vec_w[i] for i in val_only]


for arr, vj in zip(grp_by_vj, key_only):
    np.save('/home/ligia/Desktop/Emerson2017/trial/{0}.npy'.format(str(vj_unique_str[vj])), arr)


# for k, v in dic.items():
#     for ar in s:
#          #np.save(str(cat_as_str[k])+'.npy', ar)
#          #print(a
#         np.save('/home/ligia/Desktop/Emerson2017/trial/{0}.npy'.format(str(b)),f)
#         #print(ar.shape)


#   sanity check
for i in glob.glob('/home/ligia/Desktop/Emerson2017/trial/*.npy'):
    s = np.load(i)
    print(s.shape)

def create_main_file(main_file):
    main_file = main_file[~main_file.CDR3.str.contains('\*')]
    main_file['vj'] = main_file[['V', 'J']].apply(lambda x: ''.join(x), axis=1)
    main_arr_temp = np.array(main_file['vj'].values.tolist())
    main_arr = main_arr_temp.reshape(main_arr_temp.shape[0], -1)
    return main_arr

x = pd.read_csv('SearchTable-2019-05-28 21_49_18.859.tsv', delimiter ='\t')

x = x[x['Species'] == 'HomoSapiens']

x = x.reset_index(drop=True)
x = x.replace(np.nan, 'no')



tcr = []
for i in x:
    for j in i:
        if "B" in j:
            tcr.append(j)

tcrb = np.array(tcr)

vj = tcrb.reshape(tcrb.shape[0],-1)

np.save('vdjdb_vjcomb', vj)






''' for single plot of cosine'''
cos_neg = cosine_similarity(np.load('50000_main_pos_vj_only.npy'))
np.save('50000_main_pos_vj_only.npy', cos)
#%%
cos = np.load('50000_main_pos_vj_only.npy')

plt.close()
plt.figure(figsize=(15,15))
plt.matshow(cos_neg, fignum=1, cmap = 'magma')
plt.colorbar()
plt.axis('off')
plt.savefig('/home/ligia/Desktop/Emerson2017/cos_antigen/cos_neggonly_vjonly_vers2.png')

plt.close()
plt.figure(figsize=(15,15))
plt.matshow(cos, fignum=1, cmap = 'magma')
plt.colorbar()
plt.axis('off')
plt.savefig('/home/ligia/Desktop/Emerson2017/cos_antigen/cos_posonly_vjonly_vers2.png')


#%%
plt.savefig('/home/ligia/Desktop/Emerson2017/cos_antigen/{0}_org.png'.format(str(i[44:-4])))
plt.close()




my_palette = dict(zip(df.cyl.unique(), ["orange", "yellow", "brown"]))
row_colors = df.cyl.map(my_palette)

# plot
sns.clustermap(df, metric="correlation", method="single", cmap="Blues", standard_scale=1, row_colors=row_colors)



# labels = np.random.random_integers(0,5, size=50)
# lut = dict(zip(set(labels), sns.hls_palette(len(set(labels)), l=0.5, s=0.8)))
# row_colors = pd.DataFrame(labels)[0].map(lut)
#
# #Create additional row_colors here
# labels2 = np.random.random_integers(0,1, size=50)
# lut2 = dict(zip(set(labels2), sns.hls_palette(len(set(labels2)), l=0.5, s=0.8)))
# row_colors2 = pd.DataFrame(labels2)[0].map(lut2)

g=sns.clustermap(matrix, col_cluster=False, linewidths=0.1, cmap='coolwarm', row_colors=[row_colors, row_colors2])
plt.show()



''' NEWEST WAY OF PLOTTING COSINE (INCLUDE ROW COLOR CODE ON HEATMAP'''

cmv = pd.read_csv('164cmv.csv')
cmv = cmv.sort_values(['aggr_j', 'aggr_v'])
#
#
# labels = [np.array(cmv['aggr_v']), np.array(cmv['aggr_j'])]
#
# def get_labels(labels):
#     rc = []
#     lcm = []
#     for i in labels:
#         label_colormap = dict(zip(set(i), sns.color_palette('PuOr', len(set(i)))))
#         row_colors = pd.DataFrame(i)[0].map(label_colormap)
#         rc.append(row_colors)
#         lcm.append(label_colormap)
#     return rc, lcm
#
#
# rc, lcm = get_labels(labels)
#
#
# g=sns.clustermap(cos, row_cluster=False, col_cluster=False, linewidths=0, cmap='magma',
#                  row_colors=[rc[1], rc[0]], col_colors=[rc[1], rc[0]])
#
# for label in set(labels[1]):
#     g.ax_col_dendrogram.bar(0, 0, color=lcm[1][label], label =label, linewidth=0)
#
# g.ax_col_dendrogram.legend(loc="center", ncol=1)
#
# for label in set(labels[0]):
#     g.ax_col_dendrogram.bar(0, 0, color=lcm[0][label], label =label, linewidth=0)
#
# g.ax_col_dendrogram.legend(loc="center", ncol=8)
#
# g.cax.set_position([.15,.2, .03, .45])
# plt.savefig('tri.png')





v = np.load('j7_allP.npy')
j1= pd.read_csv('j7_allP.csv')
j1_l = j1['CDR3'].tolist()
squareform(cos)
# labels = [np.array(cmv['aggr_v']), np.array(cmv['aggr_j'])]
labels = [np.array(j1['aggr_v']), np.array(j1['Epitope species'])]


label_v = np.array(d['len'])
label_colormap_v = dict(zip(set(label_v),  sns.color_palette('PuOr', len(set(label_v)))))
row_colors_v = pd.DataFrame(label_v)[0].map(label_colormap_v)

#Create additional row_colors here
label_j = np.array(j1['Epitope species'])
label_colormap_j = dict(zip(set(label_j), sns.color_palette('tab20c', len(set(label_j)))))
row_colors_j = pd.DataFrame(label_j)[0].map(label_colormap_j)


# row_linkage, col_linkage = (hc.linkage(sp.distance.pdist(v, 'cosine'), method='ward')
#   for x in (v, v.T))

# linkage = hierarchy.linkage(vr, method='average')
# g=sns.clustermap(cos, row_cluster=False, col_cluster=False, linewidths=0, cmap='magma',
#                  row_colors=[row_colors_j, row_colors_v], col_colors=[row_colors_j, row_colors_v])

##old version before 14aug when you were trying out the method with CMV sequences
# g=sns.clustermap(cos, row_cluster=False, col_cluster=False, linewidths=0, cmap='magma',
#                  row_colors=[row_colors_j, row_colors_v], col_colors=[row_colors_j, row_colors_v])
#

g=sns.clustermap(squareform(cos), linewidths=0, cmap='magma',
                 row_colors=[row_colors_j, row_colors_v], col_colors=[row_colors_j, row_colors_v])


for l_v, l_j in zip(set(label_v), set(label_j)):
    g.ax_col_dendrogram.bar(0, 0, color=label_colormap_j[l_j], label=l_j, linewidth=0)
    g.ax_col_dendrogram.bar(0, 0, color=label_colormap_v[l_v], label=l_v, linewidth=0)

g.ax_col_dendrogram.legend(loc="center", ncol=9)


g.cax.set_position([.15,.2, .03, .45])
plt.savefig('hclust_j7.png')



#beta = open_csv('cmvpos_5000.csv', 'cmvneg_5000.csv', concat=True)
#
#
#
#
#beta = open_csv('/home/ligia/Desktop/Emerson2017/subsamples_cosine/5000seqs_nonimmun_sample.csv', '', ss=True)

#%%


'''COS_ANTIGEN_VDJDB file left overs'''


#for main
#aggr_cols(df_neg, 'v', 'j', main=True)

#   for vdjdb
aggr_cols(d6, 'v_family', 'j_gene', main = True)

#
# v = np.load('/home/ligia/Desktop/Emerson2017/vfam_155p.npy')
# j = np.load('/home/ligia/Desktop/Emerson2017/jgene_155p.npy')
#
# rand_pos = np.random.choice(pos_ix[0], 50000, replace=False)
#
# rand_neg = np.random.choice(neg_ix[0], 50000, replace=False)
#
# df_pos = vj.loc[rand_pos]
# df_neg = vj.loc[rand_neg]
#
# aggr_cols(df_pos, 'v', 'j', main=True)
# aggr_cols(df_neg, 'v', 'j', main=True)
#
# df_pos = df_pos.sort_values(['aggr_j', 'aggr_v'])
#
# df_neg = df_neg.sort_values(['aggr_j', 'aggr_v'])
#
# dfn = df_neg.reset_index()
# dfp = df_pos.reset_index()
#
# v_pos = [v[ix] for ix in np.array(dfp['index'])]
# v_neg = [v[ix] for ix in np.array(dfn['index'])]
#
# cos = cosine_similarity(v_pos)
# cos_neg = cosine_similarity(v_neg)




#v = np.load('vec_weights100D_155p_tr330_unique.npy') #    all vectors


# v= np.load('vec_weights100D_from300p_vdjdb.npy')

v = np.load('vec_weights100D_allP_vdjdb.npy')


#ss_noim_ix = np.load('/home/ligia/Desktop/Emerson2017/subsamples_cosine/5000seqs_noimm_ix.npy')
#all164cmv_v = np.load('/home/ligia/Desktop/Emerson2017/subsamples_cosine/all_164cmv_main_vecs.npy')
#pos_ix =  np.load('/home/ligia/Desktop/Emerson2017/subsamples_cosine/cmvpos_5000_ix.npy')
#neg_ix =  np.load('/home/ligia/Desktop/Emerson2017/subsamples_cosine/cmvneg_5000_ix.npy')
sngl = np.array(beta['index'])#np.load('j2_ix.npy') #np.load('/home/ligia/Desktop/Emerson2017/subsamples_cosine/200000_main_ix.npy')
d = get_file(beta,'vdjdb_vecs_nodups', single=True)
#check if it worked
o = np.load('/home/ligia/Desktop/Emerson2017/vdjdb_5_vecs.npy')




#   this one until u_antigen is for creating concat of vecs + seqs of antigen (vdj db) when you wanna plot multiple
#def ge_file(antigen, file):
#    embedding_matrix = np.load('vec_weights100D_from300p_vdjdb.npy')
#    e = pd.DataFrame(embedding_matrix)
#    tempbeta = beta.sort_values(['aggr_j', 'aggr_v'])
#    d = tempbeta[tempbeta['Epitope species'].str.contains(str(antigen))].reset_index(drop=True)
#    j = pd.concat([beta[['Epitope species','aggr_j', 'aggr_v', 'CDR3']], e], axis =1)#x already contains the column single which the B genes have been aggregated (both family 1 and 2)
#    j = j.sort_values(['aggr_j', 'aggr_v'])
#    j = j.drop_duplicates('CDR3') #LATEST EDIT
#    ant = j[j['Epitope species'] == str(antigen)]
#    ant = np.array(ant.drop(['Epitope species','aggr_j', 'aggr_v', 'CDR3'], axis =1))
#    return d,tempbeta, np.save(str(file), ant)
#
#beta['Epitope species'] = [i.replace('/','') for i in beta['Epitope species']]
#u_antigen = beta['Epitope species'].unique()
#for i in u_antigen:
#    ge_file(str(i), '/home/ligia/Desktop/Emerson2017/cos_antigen/{0}_300p_V2'.format(str(i)))



#%%

#cos_vdjdb_300p.npy

'''KMEANS FILE'''

# file = np.load("/home/ligia/Desktop/Emerson2017/temp_kmean/TCRBV08-02TCRBJ01-03.npy")
# silhouette_values = []
# k_values = range(2, 6)
# for k in k_values:
#     kmeans = cluster.KMeans(n_clusters=k)
#     kmeans.fit(file)
#     labels = kmeans.labels_
#     # print(labels)
#     centroids = kmeans.cluster_centers_
#     # print(centroids)
#     silhouette_score = metrics.silhouette_score(file, labels, metric='euclidean')
#     silhouette_score = silhouette_score.mean()
#     #print(silhouette_score)
#     silhouette_values.append(silhouette_score)
#
# plot.plot(k_values, silhouette_values, 'b-')
# plot.xlabel('k')
# plot.ylabel('silhouette score')
# plot.ylim([0, 1])
# plot.show()
# plot.savefig('tri2.png')
# plot.close()
#     #n =+1
#     #plot.savefig('/home/ligia/Desktop/Emerson2017/kmeans_sil/{0}.png'.format(str(all_sil.index(v))))
#     #plot.close()
#
#
#
# #all files



#   get glob index (bc the previous files were all over the place)
# glob_index = []
# for i in glob.glob('/home/ligia/Desktop/Emerson2017/trial/*.npy'):
#     glob_index.append(i)
#
# seq_idx_dic = {i: glob_index[i] for i in range(0, len(glob_index))}
#



glob_index = []
for i in sorted(glob.glob('/home/ligia/Desktop/Emerson2017/temp_kmean/*.npy')):
    glob_index.append(i)

seq_idx_dic = {i: glob_index[i] for i in range(0, len(glob_index))}
# all_files = []
# for i in glob.glob('/home/ligia/Desktop/Emerson2017/trial/*.npy'):
#     file = np.load(i)
#     all_files.append(file)
#
# def get_silh(all_files):
#     k_values = range(2, 10)
#     all_sil = []
#     for i in all_files:
#         silhouette_values = []
#         k_values = range(2, 10)
#     for k in k_values:
#         # print('k = ', k)
#         kmeans = cluster.KMeans(n_clusters=k)
#         kmeans.fit(i)
#         labels = kmeans.labels_
#         # print(labels)
#         centroids = kmeans.cluster_centers_
#         # print(centroids)
#         silhouette_score = metrics.silhouette_score(i, labels, metric='euclidean')
#         silhouette_score = silhouette_score.mean()
#         # dataset_cluster = clusteringNH.k_means_over_instances(copy.deepcopy(dataset), ['acc_phone_x', 'acc_phone_y', 'acc_phone_z'], k, 'default', 20, 10)
#         # silhouette_score = dataset_cluster['silhouette'].mean()
#         print('silhouette = ', silhouette_score)
#         silhouette_values.append(silhouette_score)
#     all_sil.append(silhouette_values)
#     for v in all_sil:
#         plot.plot(k_values, v, 'b-')
#         plot.xlabel('k')
#         plot.ylabel('silhouette score')
#         plot.ylim([0,1])
#         # plot.show()
#         plot.savefig('/home/ligia/Desktop/Emerson2017/kmeans_sil/{0}.png'.format(str(v)))
#
#
#
# def get_silh():
#     k_values = range(2, 10)
#     # all_sil = []
#     for i in glob.glob('/home/ligia/Desktop/Emerson2017/trial/*.npy'):
#         file = np.load(i)
#     # for i in all_files:
#         silhouette_values = []
#         # k_values = range(2, 10)
#         for k in k_values:
#             # print('k = ', k)
#             kmeans = cluster.KMeans(n_clusters=k)
#             kmeans.fit(file)
#             labels = kmeans.labels_
#             # print(labels)
#             centroids = kmeans.cluster_centers_
#             # print(centroids)
#             silhouette_score = metrics.silhouette_score(file, labels, metric='euclidean')
#             silhouette_score = silhouette_score.mean()
#             # dataset_cluster = clusteringNH.k_means_over_instances(copy.deepcopy(dataset), ['acc_phone_x', 'acc_phone_y', 'acc_phone_z'], k, 'default', 20, 10)
#             # silhouette_score = dataset_cluster['silhouette'].mean()
#             print('silhouette = ', silhouette_score)
#             silhouette_values.append(silhouette_score)
#         # all_sil.append(silhouette_values)
#         plot.plot(k_values, silhouette_values, 'b-')
#         plot.xlabel('k')
#         plot.ylabel('silhouette score')
#         plot.ylim([0, 1])
#         # plot.show()
#         s = plot.savefig('/home/ligia/Desktop/Emerson2017/kmeans_sil/{0}.png'.format(str()))
#     return s
#
# for i in glob.glob('/home/ligia/Desktop/Emerson2017/trial/*.npy'):
#
#
# import itertools
# cos = np.load('cosresultvj.npy')
# cos = np.load('cosresultvj_infl150.npy')
# #%%
# cos_list = cos.tolist()
# dict={}
# for i in cos_list:
#     dict[cos_list.index(i)]=i[:]
#
# def chunked(it, size):
#     it = iter(it)
#     while True:
#         p = tuple(itertools.islice(it, size))
#         if not p:
#             break
#         yield p
#
# for chunk in chunked(dict.items(), 20):
#     fig, ax = plt.subplots(20, 1, figsize=(20,18))
#     z += 1
#     for a,key, in zip(ax,dict.keys()):
#         y=dict[key]
#         print(key)
#         n=len(y)
#         x = np.linspace(1,n,n)
#         a.plot(x,y)
#         plt.savefig('/home/ligia/Desktop/Emerson2017/cosinesim_pics/150influenza{0}.png'.format(str(z)))
#
# #   compute n lowest similarities for each cmv sequences across all vjs
# lwstn_all = []
# lrgst_all =[]
# for seq in cos_list:
#     nsmall = nsmallest(3, seq)
#     nhigh = nlargest(10, seq)
#     lwstn = []
#     lrgstn = []
#     for n in nsmall:
#         lwstn.append(seq.index(n))
#     for n in nhigh:
#         lrgstn.append(seq.index(n))
#     lwstn_all.append(lwstn)
#     lrgst_all.append(lrgstn)
#
# #%%
#
# flat_lrgst = [val for inner in lrgst_all for val in inner]
# counter = Counter(flat_lrgst)
# most_occur_ix = counter.most_common()
#
# most_occur_vj = []
# for i in most_occur_ix[:100]:
#     most_occur_vj.append(seq_idx_allvjs[i[0]][38:-4])
#
# #   all cosine computations that used glob were done on the unsorted glob
# #   (so order files are stored ≠ from order it appears in folder) / use below dic
#
# glob_index_allvjs = []
# for i in glob.glob('/home/ligia/Desktop/Emerson2017/trial/*.npy'):
#     glob_index_allvjs.append(i)
#
#
# seq_idx_allvjs = {i: glob_index_allvjs[i] for i in range(0, len(glob_index_allvjs))}
# seq_idx_allvjs[653][38:-4] #    checks which are file names/vj combinations
#
#
# sep = []
# n = 12
# for vj in most_occur_vj:
#     sep.append(tuple(vj[i:i+n] for i in range(0, len(vj), n)))
#
# v = []
# j = []
# for i in sep:
#     v.append(i[0])
#     j.append(i[1])
#
# def get_count(list):
#     counter = Counter(list)
#     most_occur = counter.most_common()
#     return most_occur
#
# v_count_influenza = get_count(v)
# j_count_influenza = get_count(j)
