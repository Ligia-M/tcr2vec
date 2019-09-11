#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 11 16:50:16 2019

@author: ligia
"""

import glob
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

sns.set()

sns.set_style("white")
from sklearn.metrics.pairwise import cosine_similarity


class CosinePathogenAnalysis:
    def __init__(self, dirpathplot, v):
        self.dirpathplot = dirpathplot
        self.v = np.load(v) # ex)  'vec_weights100D_allP_vdjdb.npy'

    def get_beta(self, row):
        if "B" in row:
            return 1
        else:
            return 0

    def aggr_j(self, row, main=False, vdjdb=False):
        if vdjdb:
            if "B" in row:
                return "B{0}".format(row[6])
            if "A" in row:
                return row[0:6]
            else:
                return "other"
        if main:
            if "B" in row:
                return row[9]
            if "unresolved" in row:
                return 0

    def aggr_v(self, row, vdjdb=False, main=False):
        if vdjdb:
            return row[4:6]
        if main:
            return row[-2:]

    def open_csv(self, file, file2, vdjdb=False): #, ss=False, concat=False):
        '''
        :param file: str, file name
        :param file2: str, file name 2
        :param vdjdb: for vdjdb dataset
        :return: returns file with only beta sequences
        '''
        if vdjdb:
            x = pd.read_csv(str(file), delimiter ='\t')
            x = x[x['Species'] == 'HomoSapiens']
            x = x.reset_index(drop=True)
            x = x.replace(np.nan, 'no')
            x['beta'] = x['J'].apply(lambda row: self.get_beta(row))
            beta = x[x['beta'] == 1]
        # if ss:
        #     x = pd.read_csv(str(file))
        #     beta = x.replace(np.nan, 'no')
        # if concat:
        #     x, x2 = pd.read_csv(str(file)), pd.read_csv(str(file2))
        #     x['sample'], x2['sample'] = 1, 2
        #     cc = pd.concat([x, x2], axis = 0)
        #     beta = cc.replace(np.nan, 'no')
        #     beta = beta.reset_index(drop = True)
        return beta

    def aggr_cols(self, beta, colv, colj, vdjdb=False, main=False):
        '''
        :param beta: dataframe
        :param colv: str, name of v family column
        :param colj: str, name of j-gene column
        :param vdjdb: for vdjdb dataset
        :param main: for emerson dataset
        :return:
        '''
        if vdjdb:
            beta['aggr_v'] = beta[str(colv)].apply(lambda row: self.aggr_v(row, vdjdb=True))
            beta['aggr_j'] = beta[str(colj)].apply(lambda row: self.aggr_j(row, vdjdb=True))
            beta['aggr_v'] = [i.replace('*','') for i in beta['aggr_v']]
            beta['aggr_v'] = [i.replace('-','') for i in beta['aggr_v']]
            beta['aggr_j'] = [i.replace('B','') for i in beta['aggr_j']]
        if main:
            beta['aggr_v'] = beta[str(colv)].apply(lambda row: self.aggr_v(row, main=True))
            beta['aggr_j'] = beta[str(colj)].apply(lambda row: self.aggr_j(row, main=True))
            beta['aggr_v'] = [i.replace('na', '0') for i in beta['aggr_v']]
            beta['aggr_v'] = [i.replace('VA', '0') for i in beta['aggr_v']]
        beta['aggr_j'], beta['aggr_v'] = beta['aggr_j'].astype(int), beta['aggr_v'].astype(int)
        return beta['aggr_j'], beta['aggr_v']

    ''' SORT AND REMOVE DUPLCIATES IN ALL VDJDB SEQUENCES '''
    def get_file(self, data, file, single=True):#, vdjdb=False, main=False, concat=False):
        if single:
            data = data.reset_index()
            embedding_matrix = [self.v[ix] for ix in np.array(data['index'])]
        # if vdjdb:
        #     embedding_matrix = np.load('/home/ligia/Desktop/Emerson2017/subsamples_cosine/all_164cmv_main_vecs.npy')
        # if concat:
        #     #        ss_noim_v = [v[ix] for ix in ss_noim_ix]
        #     pos_v = [self.v[ix] for ix in pos_ix]
        #     neg_v = [self.v[ix] for ix in neg_ix]
        #     embedding_matrix = np.concatenate([pos_v, neg_v], axis=0)
        # if main:
        #     pos = [self.v[ix] for ix in pos_ix]
        #     neg = [self.v[ix] for ix in neg_ix]
        #     embedding_matrix = np.concatenate([pos, neg], axis=0)
        e = pd.DataFrame(embedding_matrix)
        j = pd.concat([data.reset_index(drop=True), e],
                      axis=1)  # x already contains the column single which the
                                # B genes have been aggregated (both family 1 and 2)
        # j = j[j['status'] == 'Cytomegalovirus +']
        # j = j.sample(50000)
        j = j.sort_values(['aggr_j', 'aggr_v']).drop_duplicates(
            'CDR3')  # LATEST EDIT #latestlatest edit: CDR3 when using vdjdb
        #    j = j.sort_values(['status','aggr_j', 'aggr_v']).drop_duplicates('amino_acid')
        #    j = j.sort_values('length').drop_duplicates('amino_acid')
        #    j = j.drop_duplicates('amino_acid') #LATEST EDIT #latestlatest edit: CDR3 when using vdjdb
        # ant = j[j['Epitope species'] == str(antigen)]
        columns_beta = data.columns.tolist()
        d = j[columns_beta]
        ant = np.array(j.drop(columns=columns_beta, axis=1))
        return d, np.save(str(file), ant)

    '''GET COSINE NPY FILES FOR MULTIPLE FILES (ANTIGENS) PLOT EACH AND SAVE FIGURE'''
    def get_pathogen_distance_plots(self):
        fail = []
        for i in glob.glob('{}*.npy'.format(self.dirpathplot)):
            try:
                cos = cosine_similarity(np.load(i))
                np.save('cos_{0}'.format(str(i[44:-4])), cos)
                plt.figure(figsize=(15, 15))
                plt.matshow(cos, fignum=1)
                plt.colorbar()
                plt.savefig('{0}_org.png'.format(str(i[44:-4])))
                plt.close()
            except:
                fail.append(i)
        np.save('fail', fail)

    '''LOOPS FOR PASSING THE ABOVE FUNCTION WHICH TAKES THE CORRECT VECTORS per J-gene
    to save them & saving the corresponding csv file w/ og values'''
    def get_vec_csvfiles(self):
        beta = self.open_csv('SearchTable-2019-05-28 21_49_18.859.tsv', '', vdjdb=True)
        beta = self.aggr_cols(beta, 'v', 'j', vdjdb = True)
        j3 = beta[beta['aggr_j'] == 3].reset_index()
        j1 = beta[beta['aggr_j'] == 1].reset_index()
        j7 = beta[beta['aggr_j'] == 7].reset_index()
        j4 = beta[beta['aggr_j'] == 4].reset_index()
        j6 = beta[beta['aggr_j'] == 6].reset_index()
        j5 = beta[beta['aggr_j'] == 5].reset_index()
        l_name = ["j1", "j3", "j4", "j5", "j6", "j7"]
        l = [j1, j3, j4, j5, j6, j7]
        j_clean = []
        z=0
        for i in l:
            z += 1
            # ix_arr = np.array(i['index'])
            # columns_beta = i.columns.tolist()
            d = self.get_file(i, 'vdjdb_{0}_vecs'.format(str(z)), single=True)
            j_clean.append(d[0])
        for d, n in zip(j_clean, l_name):
            d.to_csv('vdjdb_{0}_cleaned.csv'.format(n))


#ss_noim_ix = np.load('/home/ligia/Desktop/Emerson2017/subsamples_cosine/5000seqs_noimm_ix.npy')
#all164cmv_v = np.load('/home/ligia/Desktop/Emerson2017/subsamples_cosine/all_164cmv_main_vecs.npy')
#pos_ix =  np.load('/home/ligia/Desktop/Emerson2017/subsamples_cosine/cmvpos_5000_ix.npy')
#neg_ix =  np.load('/home/ligia/Desktop/Emerson2017/subsamples_cosine/cmvneg_5000_ix.npy')
# sngl = np.array(beta['index'])#np.load('j2_ix.npy') #np.load('/home/ligia/Desktop/Emerson2017/subsamples_cosine/200000_main_ix.npy')
# d = get_file(beta,'vdjdb_vecs_nodups', single=True)
# #check if it worked
# o = np.load('/home/ligia/Desktop/Emerson2017/vdjdb_5_vecs.npy')
#

