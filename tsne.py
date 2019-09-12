#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May  6 16:20:20 2019

@author: ligia
"""
import numpy as np
import pylab as plt
import seaborn as sns; sns.set()
import pickle
import sys; sys.path.append('/home/ligia/FIt-SNE')
from fast_tsne import fast_tsne  #package needs to be downloaded from Github page


class Tsne:
    def fit_sne(self, vec_weights, filename):
        '''
        :param vec_weights: npy file with vector weights
        :param filename: outbound file name (both plot and embedding vecs)
        :return:
        '''
        embedding_matrix = vec_weights
        model =  fast_tsne(embedding_matrix, perplexity=30, max_iter = 1000, learning_rate = 500 , seed=30)
        print('Saving model...')
        np.save('{}'.format(filename), model)
        print('Done.')
        fig = plt.figure(figsize=(15,15))
        plt.scatter(model[:,0], model[:,1], cmap='hsv', marker= ',', s = 1)
        plt.tight_layout()
        fig.savefig('{}.png'.format(filename), dpi=fig.dpi)
        return model

    def DS_init(self, vec_weights, filename):
        embedding_matrix = np.load(vec_weights)
        np.random.seed(42)
        pcainit = embedding_matrix[:, :2] / np.std(embedding_matrix[:, 0]) * 0.0001
        ind25k = np.random.choice(embedding_matrix.shape[0], 25000, replace=False)
        Z25k = fast_tsne(embedding_matrix[ind25k, :], perplexity=500, initialization=pcainit[ind25k, :], max_iter=5000,
                         seed=42)
        print('Phase 1 downsampling: Done.')
        fig = plt.figure(figsize=(15, 15))
        plt.scatter(Z25k[:, 0], Z25k[:, 1], cmap='hsv', marker=',', s=1)
        fig.savefig('/DS_phase1_{}.png'.format(filename), dpi=fig.dpi)
        pickle.dump([Z25k, []], open("downsampling_100d_unique.pickle", "wb"))
        def pdist2(A, B):
            return np.sum(A ** 2, axis=1)[:, None] + np.sum(B ** 2, axis=1)[None, :] - 2 * A @ B.T
        batchsize = 1000
        steps = int(np.floor(embedding_matrix.shape[0] / batchsize) + 1)
        position_id = np.zeros(embedding_matrix.shape[0], dtype=int)
        for i in range(steps):
            print('.', end='', flush=True)
            if i > 1 and i % 100 == 0:
                print('', flush=True)
            endind = np.min(((i + 1) * batchsize, embedding_matrix.shape[0]))
            D = pdist2(embedding_matrix[ind25k, :], embedding_matrix[i * batchsize:endind, :])
            m = np.argmin(D, axis=0)
            position_id[i * batchsize:endind] = m
        print('', flush=True)
        pickle.dump([Z25k, position_id], open("downsampling_100d_unique.pickle", "wb"))
        X = np.load(vec_weights)
        print('Phase 2 downsampling...')
        Z25k, position_id = pickle.load(open("downsampling_100d_unique.pickle", "rb"))
        init25k = Z25k[position_id, :] / np.std(Z25k[position_id, 0]) * 0.0001
        Z = fast_tsne(X, perplexity=30, initialization=init25k, max_iter=10000,
                      learning_rate=1000, seed=42, load_affinities='save')
        np.save('tsne_{}'.format(filename), Z)
