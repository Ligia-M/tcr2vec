
import glob
import pickle
from heapq import nsmallest
from itertools import groupby
from operator import itemgetter

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import ward
from scipy.spatial.distance import squareform, pdist

v = np.load('j7_allP.npy')
j1= pd.read_csv('j7_allP.csv')
j1_l = j1['CDR3'].tolist()
cos = pdist(v, 'cosine')
linkage_matrix = ward(cos)

linkage = hc.linkage(sp.distance.squareform(cos), method='average')

#C = hierarchy.fcluster(linkage_matrix, threshold, criterion="distance")

fig, ax = plt.subplots(figsize=(15, 20)) # set size
g = sns.clustermap(squareform(cos))
#ax = dendrogram(C, orientation="right", labels=j1_l)

plt.tick_params(\
    axis= 'x',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    bottom='off',      # ticks along the bottom edge are off
    top='off',         # ticks along the top edge are off
    labelbottom='off')

plt.tight_layout() #show plot with tight layout

#uncomment below to save figure
plt.savefig('ward_clusters_j1_vdjdb.png', dpi=200)

#%%

# v = np.load('j1_allP.npy')
# vecs = np.load('vec_weights100D_allP_vdjdb.npy')
#
# tsne = np.load('tsne_all1.npy')
#
# d j1= pd.read_csv('j1_allP.csv')
# cc = pd.concat([d, pd.DataFrame(tsne)], axis=1)


def collist(jgene):
    d = pd.read_csv('j{0}_allP.csv'.format(jgene))
    d['length'] = d['CDR3'].apply(lambda row: len(row))
    tsne = np.load('tsne_all{}.npy'.format(str(jgene)))
    dt = pd.concat([d, pd.DataFrame(tsne, columns=['x', 'y'])], axis=1)
    cols = dt.columns.tolist()
    return cols

cols = collist(3)


def concs(jgene):
    d = pd.read_csv('j{0}_allP.csv'.format(jgene))
    v = pd.DataFrame(np.load('j{}_allP.npy'.format(jgene)))
    d['length'] = d['CDR3'].apply(lambda row: len(row))
    tsne = np.load('tsne_all{}.npy'.format(str(jgene)))
    dt = pd.concat([d, pd.DataFrame(tsne, columns=['x', 'y'])], axis=1)
    cols = dt.columns.tolist()
    cc = pd.concat([d, pd.DataFrame(tsne, columns=['x', 'y']), v], axis=1)
    return cc



def get_s(cc, var1, val1, var2, val2, var3, val3, filtBy2=False, filtBy3=False):
    ''' cc[(cc['Epitope species'] == 'EBV') & (r['aggr_v']== 20)
    & (r['Epitope']== 'GLCTLVAML')]'''
    if filtBy2:
        filt = cc[(cc[str(var1)] == str(val1)) & (cc[str(var2)]== val2)]
    if filtBy3:
        filt = cc[(cc[str(var1)] == str(val1)) & (cc[str(var2)]== val2)
                 & (cc[str(var3)]== str(val3))]
    s = filt.sort_values(['x', 'y'])
    return s

#
# for i in nrch:
#     if 'bin_typ' in i[1]:
#         print(i)

#only focus on vj combinations
nrch = pickle.load(open('/home/ligia/Desktop/Emerson2017/fish_nrched_hlaag.p', 'rb'))
vj_of_int = []
j_of_int = []
for i in nrch:
    if i[0][0].isdigit():
        vj_of_int.append(i)
    else:
        j_of_int.append(i)

#create a dictionary per J per V and then export to a CSV if at least 8 rows have == type n v
sorter = sorted(vj_of_int, key=itemgetter(0))
s_j = sorted(j_of_int, key=itemgetter(0))
grouper = groupby(sorter, key=itemgetter(0))
g_j = groupby(s_j, key=itemgetter(0))

res = {k: list(map(itemgetter(1), v)) for k, v in grouper}
dic_j = {k: list(map(itemgetter(1), v)) for k, v in g_j}

for key in res:
    v_fam = int(''.join(list(filter(str.isdigit, key[-3:]))))
    jgene = int(''.join(list(filter(str.isdigit, key[:3]))))
    type = []
    # epi = []
    for i in res[key]:
        if 'bin_typ' in i:
            type.append(i[7:])
        # if 'bin_ep' in i:
        #     epi.append(i[6:])
    type = ['DENV34' if w == 'DENV3/4' else w for w in type]
    for i in type:
        s = get_s(concs(jgene), 'Epitope species', str(i), 'aggr_v',
                  v_fam, 'Epitope', 'IVTDFSVIK', filtBy2=True)
        if s.shape[0] >=8:
            s.to_csv('/home/ligia/Desktop/trial_removed/'+str(jgene)+str(i)+'_{}.csv'.format(str(v_fam)), index = False)



df_name = []
dfs = []
vecs_mean = []
avg = []
smol = []
nsmal_locs = []
for i in sorted(glob.glob('/home/ligia/Desktop/trial_removed/*.csv')):
    df_name.append(i)
    dfs.append(pd.read_csv(i))
    df = pd.read_csv(i)
    vec_arr = np.array(df.drop(cols, axis=1))
    vecs_mean.append(np.mean(vec_arr, axis=0))
    cos = pdist(v, 'cosine')
    cos = 1 - cos
    avg.append(cos.mean())
    nsmal = nsmallest(10, cos)
    smol.append(nsmal)
    ri = df.reset_index()
    locs = []
    for n in nsmal:
        ix_smol = np.where(squareform(cos) == 0.00986668)
        for idx in ix_smol:
            locs.append([ri.loc[idx[0]], ri.loc[idx[1]]])
    nsmal_locs.append(locs)




np.where(np.asarray(df_name) == '/home/ligia/Desktop/trial_removed/3HIV-1_2.csv')

dfs[32][['Epitope','x', 'y', 'length', 'MHC B', 'MHC A', 'aggr_v', 'aggr_j']]


dfs[42][['Epitope', 'x', 'y', 'length', 'CDR3', 'Epitope species']]


bad_df = dfs[41].index.isin([12, 9, 13, 16])

dfs[41][~bad_df]
# np.where(np.asarray(df_name) == '/home/ligia/Desktop/trial_removed/2EBV_2.csv')


dfs[43][16:38][~bad_df]['Epitope'].value_counts()



#
def get_xy(row):
    if -20 <= row['x'] <= 0 and -40<= row ['y'] <= -20: #in np.arange(-10.0, 0.0):# and row[1] in np.arange(-20, -10):
        return 'c1'
    if 0 <= row['x'] <= 20 and -50 <= row['y'] <= -35:
        return 'c2'
    # if 10 <= row[0] <= 20 and -20 <= row[1] <= -9:
    #     return 'c3'
    # if 25 <= row[0] <= 35 and 25 <= row[1] <= 39:
    #     return 'c4'
    else:
        return 0
#
dfs[20]['clusts'] = dfs[20].apply(lambda row: get_xy(row), axis=1)
# #
# # dfs[20]['clusts'].value_counts()
#
# # avg = []
# # smol = []
# # nsmal_locs = [ ]
# # for i in s['j1hivv2'].unique():
# #     c1 = s[s['j1hivv2']== i]
# #     ix_ = c1.index.tolist()
# #
# #     embedding_matrix = [v[ix] for ix in ix_]
# #     vecs_mean = np.mean(embedding_matrix, axis=0)
# #
# #     cos = pdist(embedding_matrix, 'cosine')
# #     cos = 1 - cos
# #     avg.append(cos.mean())
# #     nsmal = nsmallest(10, cos)
# #     smol.append(nsmal)
# #     ri = c1.reset_index()
# #     locs = []
# #     for n in nsmal:
# #         ix_smol = np.where(squareform(cos) == n)
# #         for idx in ix_smol:
# #             locs.append([ri.loc[idx[0]], ri.loc[idx[1]]])
# #     nsmal_locs.append(locs)
#
#
#
#
#
#
# # # vecs = np.load('ebv_cluster_j2_vecs.npy')
# # '''THIS ONE WILL REMOVE THE 1s/THE MIRRORED VALUES '''
# # cos = pdist(embedding_matrix,'cosine')
# # cos = 1-cos
# # cos_mean = cos.mean()
# #
# # nsmallest(10, cos)
# #
# # np.where(squareform(cos) == 0.9945136417521261)
# #
# # vecs_mean = np.mean(embedding_matrix, axis=0)
# #
# # c1.loc[[3733,3767]]
# #
# # vecs_mean = np.mean(embedding_matrix, axis=0)
# # np.save('1hivv2_mean',vecs_mean )
#
#
# # vec_arr = np.array(e.drop(cols, axis=1))
#
dfs[41][~bad_df][['Epitope', 'length', 'MHC B', 'MHC A', 'aggr_v']]
#
#
# for i in sorted(glob.glob('/home/ligia/Desktop/trial_removed/*.csv')):
#     dfs.append(pd.read_csv(i))
#
#
ep = dfs[20][dfs[20]['clusts'] == 'c1']
ep2 = dfs[20][dfs[20]['clusts'] == 'c2']
#
ep2 = ep2.drop('clusts', axis =1 )
#
#




vec_arr = np.array(dfs[41][~bad_df].drop(cols, axis=1))

hiv3 = np.mean(vec_arr, axis=0)
HIV3_filt= np.mean(vec_arr, axis=0)
HIVfiv = np.mean(vec_arr, axis=0)

cos = 1 - pdist([HIV5_filt, hiv3], 'cosine')
cos = 1 - cos
cos.mean()


mm = np.mean(vec_arr, axis=0).shape
np.save('/home/ligia/Desktop/Emerson2017/vec_means/5HIV-1_9_min12_9_13_16ix_MEAN', HIV5_filt)


chungus = []
name = []
for i in sorted(glob.glob('/home/ligia/Desktop/Emerson2017/vec_means/*.npy')):
    # name.append(i.lower())
    name.append(i)
    chun = np.load(i)
    chungus.append(chun)
cos = 1 - pdist(chungus,'cosine')



cmv = []
denv = []
htlv = []
ebv = []
hcv = []
hiv = []

for n, c in zip(name,chungus):
    if 'cmv' in n:
        cmv.append(c)
    if 'hiv' in n:
        hiv.append(c)
    if 'denv' in n:
        denv.append(c)
    if 'htlv' in n:
        htlv.append(c)
    if 'ebv' in n:
        ebv.append(c)
    if 'hcv' in n:
        hcv.append(c)




cos = 1 - pdist(hiv,'cosine')





np.where(np.asarray(name) == '/home/ligia/Desktop/Emerson2017/vec_means/')

vecs = np.load('/home/ligia/Desktop/Emerson2017/vec_means/flu_c2_vecmean.npy')
cos = 1 - pdist(vec_arr, 'cosine')
np.where(squareform(cos) == 0.9820428418637471)
cos = 1 - cos
cos.mean()

hiv3 = np.load('/home/ligia/Desktop/Emerson2017/vec_means/hivj3_v2_103116_mean.npy')


n = []
for i in name:
    n.append(i[42:-9])


fig = plt.figure()
ax = fig.add_subplot(111)
cax = ax.matshow(sqmat)
fig.colorbar(cax)
# ax.xticks(rotation=90)
ax.set_xticks(np.arange(len(n)))
ax.set_xticklabels(['']+n, rotation=90)
ax.set_yticks(np.arange(len(n)))
ax.set_yticklabels(['']+n)
# plt.colorbar()
plt.savefig('vecmeans.png', bbox_inches='tight')
plt.show()


# v_rand_1 = []
# v_rand_2 = []
# v_rand_3 = []
# v_rand_4 = []
# v_rand_5 = []
# v_rand_6 = []
# v_rand_7 = []





df = pd.read_csv('/home/ligia/Desktop/Emerson2017/vdjdb_notsorted.csv')
vecs = np.load('vec_weights100D_allP_vdjdb.npy')


epitopes

nmb = range(1, 8)
times = 10000
all_rand = []
for i in nmb:
    #idx = df[df['aggr_j'] == i].index.tolist()
    v_randpath = []
    for p in d['Epitope species'].unique().tolist():
        if p == 'HSV-2' == 'TriticumAestivum':
            continue
        idx_path = d[(d['aggr_j'] == i) & (d['Epitope species'] == p)].\
            index.tolist()
        v_rand = []
        for t in range(times):
            try:
                v = vecs[idx_path][np.random.choice(vecs[idx_path].shape[0], 30), :]
                v_rand.append(v)
            except:
                print(idx_path)
        v_randpath.append(v_rand)
    all_rand.append(v_randpath)


allmeans = []
for i in all_rand:
    cos_per_j = []
    for pa in i:
        cosMean = []
        for ii in pa:
            cos = 1 - pdist(ii, 'cosine')
            cosMean.append(cos.mean())
        cos_per_j.append(np.mean(cosMean))
    allmeans.append(cos_per_j)


nmb = range(1, 8)
times = 10000
all_rand = []
for i in nmb:
    idx = d[d['aggr_j'] == i].index.tolist()
    v_rand = []
    for t in range(times):
        v = vecs[idx][np.random.choice(vecs[idx].shape[0], 30, replace=False), :]
        v_rand.append(v)
    all_rand.append(v_rand)

allmeans = []
for i in all_rand:
    cosMean = []
    for ii in i:
        cos = 1 - pdist(ii, 'cosine')
        cosMean.append(cos.mean())
    allmeans.append(np.mean(cosMean))

result from above^
[0.877107366701917, 0.8729527318927919, 0.9106933820594567, 0.8748666813518962, 0.8630085152914778, 0.8521075608786568, 0.9267567278891188]


times = 10000
all_rand = []
for t in range(times):
    v = vecs[np.random.choice(vecs.shape[0], 30, replace=False), :]
    all_rand.append(v)


allmeans = []
for i in all_rand:

cosMean = []
for ii in all_rand:
    cos = 1 - pdist(ii, 'cosine')
    cosMean.append(cos.mean())

np.mean(cosMean)


nmb = range(1, 8)
type_all = []
epi_all = []
hla_all = []
len_all = []
for n in nmb:
    newDict = { key:value for (key,value) in res.items() if int(''.join(list(filter(str.isdigit, key[:3])))) == n}
    type = []
    epi = []
    hla = []
    len = []
    for key in newDict:
        # type = []
        # epi = []
        # hla = []
        # len = []
        for i in newDict[key]:
            if 'bin_typ' in i:
                type.append(i[7:])
            if 'bin_ep' in i:
                epi.append(i[6:])
            if 'bin_len' in i:
                len.append(i[7:])
            if 'bin_mhc' in i:
                hla.append(i[7:])
    type_all.append(type)
    epi_all.append(epi)
    hla_all.append(hla)
    len_all.append(len)





import numpy as np
import matplotlib.pyplot as plt
x = map(lambda x: x[0], type_all)
y = map(lambda x: x[1], type_all)
plt.bar(x,y)
plt.savefig('bar.png')
