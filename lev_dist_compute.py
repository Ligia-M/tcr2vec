#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 24 10:29:11 2019

@author: ligia
"""

import numpy as np    
from Levenshtein import distance
from scipy.spatial.distance import pdist, squareform
import matplotlib.pyplot as plt
# my list of strings
import seaborn as sns
import pickle

''' compute levenshtein distance '''
class lev_analyses():
    def __init__(self, d):
        self.d = d
    def get_lev(self, feat, fname, cdr3col, col1, col_val1, col2, col_val2,
                farr, outname, strings_given=False, strings_not=False,
                per_1=False, per_2=False, per_all=False, plot_col_colors=False):
        '''d: pandas dataframe
        feat = str, a feature to sort the data (hla, epitope, length,  etc)
        fname = npy file that contains the strings if boolean strings_given = True
        farr = npy file with corresponding feature values if boolean strings_given = True
        strings = variable in string format (sequence)
        cdr3col = str name of cdr3 column in d
        col1, col2 = str, column name ex) 'aggr_j'
        col_val1, colval2 = column value for filtering dataframe ex) 1 or 'HLA-DRA*01'
        outname = str, name of outbound file
        if strings_not = you want to extract yourself the strings + feature to sort strings with
            if per_1 = later visualize for a single characteristic (ex. J-gene)
            if per_2 = visualize data filtered for 2 characteristics (ex. J-gene / V-family)
            if per_all = visualize entire data (d) sorted by feature
                        (directs visualization to exhibit certain patterns)
        if strings_given = strings already extracted in npy file format '''
        if strings_not:
            if per_1:
                d = self.d[self.d[str(col1)] == col_val1].sort_values(str(feat))
            if per_2:
                d = self.d[(self.d[str(col1)] == col_val1) & (d[str(col2)] == col_val2)].sort_values(str(feat))
            if per_all:
                d = self.d.sort_values(str(feat))
            arr_labels = np.array(d[str(feat)])
            strings = np.array(d[str(cdr3col)])
        if strings_given:
            strings = np.load(str(fname))
            arr_labels = np.load(farr)
        # prepare 2 dimensional array
        transformed_strings = np.array(strings).reshape(-1,1)
        print('Calculating distance matrix ...')
        # calculate condensed distance matrix by wrapping the Levenshtein distance function
        distance_matrix = pdist(transformed_strings,lambda x,y: distance(x[0],y[0]))
        np.save('levmatrix_{}'.format(str(outname)), distance_matrix)
        print('model saved')
        # plot
        if plot_col_colors:
            lcm = dict(zip(set(arr_labels), sns.color_palette('PuOr',
                                                              len(set(arr_labels)))))
            rc = pd.DataFrame(arr_labels)[0].map(lcm)
            #sqmat =squareform(np.load('distancematrix_vdjdb.npy'))
            print('Generating plot...')
            g = sns.clustermap(squareform(disstance_matrix), row_cluster=False,
                               col_cluster=False, linewidth=0,cmap='magma',
                               row_colors=[rc], col_colors=[rc])
            for al in set(arr_labels):
                g.ax_col_dendrogram.bar(0,0,color=lcm[al], label=al, linewidth=0)
            g.ax_col_dendrogram.legend(loc='center', ncol=7)
            g.cax.set_position([.15,.2,.03,.45])
            # plt.matshow(sqmat, fignum=1)
            # plt.xticks(rotation=90)
            # plt.colorbar()
            ax = g.ax_heatmap
            ax.axis('off')
            plt.savefig('lev_{}.png'.format((str(outname))))
            plt.close()
    '''For any pairwise matrix (M x M) if you want to compute mean and SD'''
    def matr_stats(self, fname):
        distance_matrix = np.load(fname)
        print('Average:{}'.format(np.mean(distance_matrix)))
        print('Standard deviation:{}'.format(np.std(distance_matrix)))
    '''Generate random levensthein matrix'''
    def rndm_lev_mat(self, perm_n, n, col1, col2, col_val1, col_val2, feat_vals, feat_col, filt_1=False, filt_2=False, filt_more=False):
        ''' d: pandas dataframe
            perm_n: int, permuation number, ex) 10 000
            n = int, size of random sample
            col1, col2 = str, column name ex) 'aggr_j' or 'aggr_v'
            col_val1, colval2 = int or str, column value for filtering dataframe ex) 1 or 'HLA-DRA*01'
            feat_vals: a list of classes of feature of interest, ex) ['HLA-B*07', 'HLA-DRA*01', 'HLA-A*02','HLA-A*24']
            feat_col: str name of feature column feat_vals come from, ex) 'hla_ag'
            '''
        avg =[]
        sd = []
        for i in range(1, perm_n):
            if filt_1:
                rand = self.d[self.d[col1] == col_val1].sample(n)
            if filt_2:
                rand = self.d[(self.d[col1]==col_val1)& (self.d[col2]== col_val2)].sample(n)
            if filt_more:
                rand = self.d[self.d[col1] == col_val1]
                rand = rand[self.d[feat_col].isin(feat_vals)].sample(n)
            strings = np.asarray(rand['CDR3'])
            transformed_strings = np.array(strings).reshape(-1,1)
            distance_matrix = pdist(transformed_strings,lambda x,y: distance(x[0],y[0]))
            avg.append(np.mean(distance_matrix))
            sd.append(np.std(distance_matrix))
        print('Average for {} permutations'.format(perm_n), np.mean(avg))
        print('Standard deviation for {} permutations'.format(perm_n), np.std(sd))
    ''' Generate permuation test between two distance matric, ex) test vs random'''
    def exact_mc_perm_test(self, mat1, mat2, perm_n):
        n, k = len(mat1), 0
        diff = np.abs(np.mean(mat1) - np.mean(mat2))
        zs = np.concatenate([mat1, mat2])
        for j in range(perm_n):
            np.random.shuffle(zs)
            k += diff <= np.abs(np.mean(zs[:n]) - np.mean(zs[n:]))
        return k / perm_n
