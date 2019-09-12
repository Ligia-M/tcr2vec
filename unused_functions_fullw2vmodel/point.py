from math import sqrt
from scipy.spatial.distance import squareform, pdist

# vecs = np.load('ebv_cluster_j2_vecs.npy')

cos = pdist(embedding_matrix,'cosine')
cos = 1-cos
cos_mean = cos.mean()

def nnd(points, density):
    cos = pdist(vecs, 'cosine')
    cos = 1 - cos
    R0 = cos.mean()
    Re = 1.0 / (2.0 * sqrt(density))
    R = R0 / Re
    Var = (4 - pi) / (4 * density * pi)
    z = (R - 1) / sqrt(Var)
    return R0, Re, R, Var, z



import pandas as pd


j2 = pd.read_csv('/home/ligia/Desktop/Emerson2017/j2_vdjdb.csv')

j2 = read.csv('/home/ligia/Desktop/Emerson2017/j2_vdjdb.csv')

def get_bin(row, nmbr, ep, v=False, epitope=False):
    if v:
        if row == nmbr:
            return 1
        else:
            return 0
    if epitope:
        if row == str(ep):
            return 1
        else:
            return 0


for i in j2.aggr_v.unique():
    j2['bin_v{}'.format(str(i))] = j2['aggr_v'].apply(lambda row: get_bin(row, i, '',  v=True))



j2['bin_v'] = j2['aggr_v'].apply(lambda row: get_bin(row, 20, '',  v=True))
j2['bin_epitope'] = j2['Epitope species'].apply(lambda row: get_bin(row,0,'EBV', epitope=True))


j2['bin_cmv'] = j2['Epitope species'].apply(lambda row: get_bin(row,0,'CMV', epitope=True))
j2['bin_hiv'] = j2['Epitope species'].apply(lambda row: get_bin(row,0,'HIV-1', epitope=True))
j2['bin_im'] = j2['Epitope species'].apply(lambda row: get_bin(row,0,'HomoSapiens', epitope=True))
j2['bin_flu'] = j2['Epitope species'].apply(lambda row: get_bin(row,0,'InfluenzaA', epitope=True))
j2['bin_len'] = j2['length'].apply(lambda row: get_bin(row,13, '', v=True ))


j2['length'] = j2['CDR3'].apply(lambda row: len(row))


j2.to_csv('j2_vdjdb.csv')


xtab = pd.crosstab(y['bin_len'], y['bin_epitope'],margins = False)

scipy.stats.fisher_exact(xtab, alternative='greater')


ss = j2[j2$aggr_v == 19,]

lss = ss[ss$length == 15,]
dim(ss)
#R code
xtab = table(ss$MHC.A, ss$bin_epitope)
xtab
fisher.test(xtab, workspace = 200000, hybrid = FALSE,
            hybridPars = c(expect = 5, percent = 80, Emin = 1),
            control = list(), or = 1, alternative = "two.sided",
            conf.int = TRUE, conf.level = 0.95,
            simulate.p.value = TRUE, B = 10000000)


j2.groupby(['aggr_v'])['Epitope species'].count()

j2.groupby(['Epitope species'])['aggr_v'].count()

j2[j2['aggr_v']==2]['Epitope species'].value_counts()


library(reticulate)
np = import('numpy')
v = np$load('ebv_cluster_j2_vecs.npy')





e = np.load('/home/ligia/Desktop/Emerson2017/tsne_100_vdjdb_j2_ppl30_LR500_iter1000.npy')
d = pd.read_csv('j2_vdjdb.csv')

cc = pd.concat([d, pd.DataFrame(tsne)], axis=1)


def get_s(cc, var1, val1, var2, val2, var3, val3, filtBy2=False, filtBy3=False):
    ''' cc[(cc['Epitope species'] == 'EBV') & (r['aggr_v']== 20)
    & (r['Epitope']== 'GLCTLVAML')]'''
    if filtBy2:
        ebv = cc[(cc[str(var1)] == str(val1)) & (cc[str(var2)]== val2)]
    if filtBy3:
        ebv = cc[(cc[str(var1)] == str(val1)) & (cc[str(var2)]== val2)
                 & (cc[str(var3)]== str(val3))]
    s = ebv.sort_values(0)
    return s

s = get_s(cc,'Epitope species', 'InfluenzaA', 'aggr_v',
          19, 'Epitope', 'GLCTLVAML', filtBy2=True)

s[4:]






cc[(cc['Epitope species'] == 'EBV') & (r['aggr_v']== 20) & (r['Epitope']== 'GLCTLVAML')]
flu = cc[cc['Epitope species'] == 'InfuenzaA']

#1st one worked for ebv
def get_xy(row):
    if -10<=row[0]<=0 and -20<=row[1]<=-10: #in np.arange(-10.0, 0.0):# and row[1] in np.arange(-20, -10):
        return 'c1'
    if -10 <= row[0] <= 10 and 0 <= row[1] <= 20:
        return 'c2'
    if 10 <= row[0] <= 20 and -22 <= row[1] <= -10:
        return 'c3'
    if 25 <= row[0] <= 35 and 25 <= row[1] <= 39:
        return 'c4'



# c2
# y = 0, 20
# x = -10, 10
#
# c3
# y = -10, -22
# x = 10, 20
#
# c4
# y = 25, 39
# x= 25, 35


s['c1'] = s.apply(lambda row: get_xy(row), axis=1)
s['c1'].value_counts()

c1= s[s['c1']== 'c1']
c2= s[s['c1']== 'c2']
c3= s[s['c1']== 'c3']
c4= s[s['c1']== 'c4']


c2.to_csv('flu_c2.csv')
c3.to_csv('flu_c3.csv')
c4.to_csv('flu_c4.csv')

l = [c1, c2,c3,c4]

for i in l:

c1 = c1.drop([0,1], axis =1)
c2 = c2.drop([0,1], axis =1)
c3 = c3.drop([0,1], axis =1)
c4 = c4.drop([0,1], axis =1)






fl1 = np.load('flu_cluster20_subc1_j2_vecs.npy')
fl1m = np.mean(fl1, axis=0)


np.save('flu_c1_vecmean', fl1m)
np.save('flu_c2_vecmean', fl2m)
np.save('flu_c3_vecmean', fl3m)
np.save('flu_c4_vecmean', fl4m)












'''filter data by J, sort values by ascending for 0 and 1
compare both columsn 0 and 1 --> per row, select all similar epitope species, check 0 and 1 if at least 5 sequences 
have 0 and 1 fall between
 arange of 10 for both append values to list named epitopespecies_c1 etc for how many clust3ers you have. IF value
 outside of 10 range for 0 and 1 by a margin of +/- 10 --> take this value as new center and repeat the above
 repeat above for every epitope species with every new "center" u call this c1 c2 c3 etc
   
   
       controls: 
    for above retrieve all Vs and Js that have clusters on, retrieve the # of sequences per cluster 
    make a cluster of == size, within == vj 
    make a new cluster of a specific epitope == J but the XY distance is huge (so != a cluster) 
    select random sequences (≠ vjs) n create cluters 
    
    
   
   try to create a new column per epitope spcies w/ values c1 c2 c3 etc and NOT A VARIABLE PER THING 
   once this loop is done --> drop columns 0 and 1 
   per column for each unique value of the column (i.e c1 c2 etc) 
   filter the data with it, concatenate the main VECS w. it, take the mean and save each mean
    
    
    
    compute similarity within == J if ≠ clusters of == epitope 
    compute between J's if == epitope --> all possible combinations
    compute within J btw ≠ disease if both high is it bc they share epitope? save what variables they are == on
    compute btw Js for ≠ diseases
    
    then repeat the above but w/ the controls 
    what do you get 

     
     '''



def get_s(cc, var1, val1, var2, val2, var3, val3, filtBy2=False, filtBy3=False):
    ''' cc[(cc['Epitope species'] == 'EBV') & (r['aggr_v']== 20)
    & (r['Epitope']== 'GLCTLVAML')]'''
    if filtBy2:
        ebv = cc[(cc[str(var1)] == str(val1)) & (cc[str(var2)]== val2)]
    if filtBy3:
        ebv = cc[(cc[str(var1)] == str(val1)) & (cc[str(var2)]== val2)
                 & (cc[str(var3)]== str(val3))]
    s = ebv.sort_values(0)
    return s

s = get_s(cc,'Epitope species', 'InfluenzaA', 'aggr_v',
          19, 'Epitope', 'GLCTLVAML', filtBy2=True)




vecs = np.load('vec_weights100D_allP_vdjdb.npy')

tsne = np.load('tsne_100_vdjdb_vecsALLP_ppl40_LR500_iter3000.npy')

d = pd.read_csv('vdjdb_notsorted.csv')
cc = pd.concat([d, pd.DataFrame(tsne)], axis=1)


j1 = cc[cc['aggr_j']==1]
j2 = cc[cc['aggr_j']==2]
j3 = cc[cc['aggr_j']==3]
j4 = cc[cc['aggr_j']==4]
j5 = cc[cc['aggr_j']==5]
j6 = cc[cc['aggr_j']==6]
j7 = cc[cc['aggr_j']==7]


l = [j1, j2, j3, j4, j5, j6, j7]
ln = ['j1', 'j2', 'j3', 'j4', 'j5', 'j6', 'j7']
cnames

for i, j in zip(l, ln):
    np.save('{}_allP'.format(j), np.array(i.drop(columns = cnames, axis =1)))
    i[cnames].to_csv('{}_allP.csv'.format(j), index=False)




ref0 = 0
ref1 = 0
cnt = 0
out = []
for r0, r1 in zip(sort['0'],sort['1']):
    if r0-ref0 > 10 or r0-ref0 <-10 and r1-ref1 > 10 or r1-ref1 <-10:
        cnt += 1
        ref0 = r0
        ref1 = r1
    out.append(cnt)

sort['trial'] = out


# sort['trial'] = sort.apply(lambda row: clusters(row), axis =1)



def get_clusters():
    for j in sort['aggr_j'].unique():
        temp_jdf = sort[sort['aggr_j'] == j]
        for i in temp_jdf['Epitope species'].unique():
            sort[str(j)+'{}_clusters'.format(i)] = cc.apply(lambda row: clusters(row, i), axis =1  )


# def get_bin(row, nmbr, ep, number_val=False, epitope=False):
#     if number_val:
#         if row == nmbr:
#             return 1
#         else:
#             return 0
#     if epitope:
#         if row == str(ep):
#             return 1
#         else:
#             return 0


