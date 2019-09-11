#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 14 13:04:31 2019

@author: ligia
"""

# from multiprocessing import Pool, Process
import numpy as np
import pandas as pd
from gensim.models.word2vec import Word2Vec
from generate_embeddings import W2vProcessing


class GetGenerators:
    def __init__(self, x, n):
        '''
        :param x: npy file with all cdr3 sequence names from emerson dataset
        :param n: int, number of random sequences to create generator
        '''
        self.arr_fl = x.ravel()
        self.w2v = Word2Vec.load('w2vmodel_100D_allP.bin')
        self.n = n

    def gen_prep(self, filename, epispe, cmv=False, pathogen=False):
        '''
        :param filename: str,  name of generators
        :param epispe: str, if pathogen=True, which pathogen to extract sequences (ex. 'InfluenzaA')
        :param cmv: bool, true = 164 cmv sequences
        :param pathogen: bool, true = pathogen of choice
        :return: reconstructs vectors from w2v model, used in topk analysis generator function
        '''
        if cmv:
            ext_seqs = pd.read_csv('164cmv.csv')
        if pathogen:
            data = pd.read_csv('SearchTable-2019-05-28 21_49_18.859.tsv', sep='\t')
            data = data[data['Species'] == 'HomoSapiens']
            c = data['CDR3'].isin(self.arr_fl.astype(str))
            kn_seqs = data[c]
            ext_seqs = kn_seqs[kn_seqs['Epitope species'] == epispe].\
                drop_duplicates('CDR3').sample(self.n)
            ext_seqs.to_csv('{}_cdr3name.csv'.format(epispe))
        ge = W2vProcessing('/home/ligia/Desktop/Emerson2017/')
        tokens_sum = ge.processing_token_sum(ext_seqs, 'CDR3')
        ext_seqs['tokens_sum'] = [inner_list for inner_list in tokens_sum]
        ext_seqs['vec_sequence'] = list(map(lambda exp_seq:
                                        ge.get_sequence_vectors(self.w2v, exp_seq),
                                        ext_seqs.tokens_sum))
        vec_weights = np.array(list(map(np.array, ext_seqs.vec_sequence)))
        print(vec_weights)
        np.save('generatorseqs_{}'.format(filename), vec_weights)


