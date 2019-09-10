#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 21 15:32:17 2019

@author: ligia
"""

import numpy as np
import pandas as pd
import os
import matplotlib.pyplot as plt
import seaborn as sns
sns.set()
from scipy.spatial.distance import squareform, pdist
import annot_seqs as anseq


class DataVis:
    def __init__(self, emerson, vdjdb, dist_val):
        self.emerson = emerson
        self.vdjdb = vdjdb
        self.clr = "#221631"
        self.pal1 = 'tab20c'
        self.pal2 = 'PuOr'
        self.dist = dist_val
        self.feat_list = ['aggr_v', 'aggr_j', 'len', 'Epitopespecies', 'Epitope', 'hla_ag']
        self.feat_emerson = ['v_family', 'j_gene', 'sample_name', 'Age', 'Biological Sex', 'lengh']

    ''' feature counts '''
    def barplot_percentages(self, filename, v1, v1name, exp1=False , exp2=False):
        '''
        :param filename: str, name of barplot
        :param exp1: true: for emerson dataset show count frequency as percentages
        :param exp2: true: for vdjdb dataset
        :param v1: str, data attribute you want to look at frequency percentage (ex. 'j-gene')
        :param v1name: str, name on plot axis of the attribute 'J-Gene'
        :return: a bar plot with y axis = percentages, x axis = classes of attribute
        '''
        if exp1:
            size = range(len(self.emerson))
            temp = pd.DataFrame({'size': size, 'var':self.emerson[str(v1)]})
        if exp2:
            size = range(len(self.vdjdb))
            temp = pd.DataFrame({'size': size, 'var':self.vdjdb[str(v1)]})
        sns.set_style("white")
        ax = sns.barplot(x="var", y="size", data = temp, estimator=lambda x: len(x) / len(temp) * 100,
                         color=self.clr, order = temp['var'].value_counts().sort_values().index)
        ax.set(ylabel="Percentage", xlabel= str(v1name))
        plt.xticks(rotation=90)
        plt.savefig('{}.png'.format(str(filename)), bbox_inches='tight')
        plt.close()

    ''' Heatmap J family by V family '''
    def heatmap_per_exp(self,v1, v1_name, v2, v2_name,filename, exp1=False, exp2=False):
        '''
        :param v1: str, variable on y axis (ex. 'aggr_j')
        :param v1_name: str, name on the yaxis (ex. J-Gene)
        :param v2: str, variable on x axis (ex. 'aggr_v' or any other data attribute: age,
                    sample name, epitope species, epitope etc)
        :param v2_name: str,name on the xaxis (ex. V-Family)
        :param filename: str, name of plot for saving
        :param exp1: bool, emerson dataset (cdr3 column has different name)
        :param exp2: bool, vdjdb dataset (cdr3 column has different name)
        :return: a heatmap to show enrichment of attribute per another attribute
        '''
        if exp2:
            plot = self.vdjdb.groupby([v1, v2])['CDR3'].count().sort_values(
                ascending = True).unstack()
        if exp1:
            plot = self.emerson.groupby([v1, v2])['amino_acid'].count().sort_values(
                ascending=True).unstack()
        cmap = sns.cubehelix_palette(light=1, as_cmap=True)
        ax_hmj = sns.heatmap(plot, annot=False, fmt="g", cmap= cmap, xticklabels=True, yticklabels=True)
        ax_hmj.set(ylabel=v1_name, xlabel= v2_name)
        #ax_hmj.set_xticklabels(ax_hmj.get_yticklabels(), rotation = 45)
        plt.savefig('{}.png'.format(str(filename)), bbox_inches='tight')
        plt.close()

    ''' WHEN YOU WANT TO PLOT 2 LABELS AT ONCE'''
    def mat_rowcolors(self, d, v1, v2, mat, filename):
        '''
        :param d: dataframe a matrix is generated from
        :param v1: str, 1st variable to specify matrix ordering
                (ex. 10 seqs are J-1, next 20 are J-2 etc)
        :param v2: str, 2st variable to specify matrix ordering
        :param mat: a matrix, can be levenshtein, or cosine similarity
                    (need to have column and row ordering)
        :param filename: name of file when saving plot
        :return: a matrix plot with color code specifying 2 labels (ex. J-gene and V-family)
        '''
        label_1 = np.array(d[v1])
        label_colormap_1 = dict(zip(set(label_1), sns.color_palette(self.pal2, len(set(label_1)))))
        row_colors_1 = pd.DataFrame(label_1)[0].map(label_colormap_1)
        label_2 = np.array(d[v2])
        label_colormap_2 = dict(zip(set(label_2), sns.color_palette(self.pal1, len(set(label_2)))))
        row_colors_2 = pd.DataFrame(label_2)[0].map(label_colormap_2)
        g=sns.clustermap(squareform(mat), linewidths=0, cmap='magma',
                         row_colors=[row_colors_2, row_colors_1], col_colors=
                         [row_colors_2, row_colors_1])
        for l_1, l_2 in zip(set(label_1), set(label_2)):
            g.ax_col_dendrogram.bar(0, 0, color=label_colormap_2[l_2], label=l_2, linewidth=0)
            g.ax_col_dendrogram.bar(0, 0, color=label_colormap_1[l_1], label=l_1, linewidth=0)
        g.ax_col_dendrogram.legend(loc="center", ncol=9)
        g.cax.set_position([.15,.2, .03, .45])
        plt.savefig('{}.png'.format(filename))
        plt.close()

    def projection_plots(self, d, e, feature, fname, um=False, fsne=False):
        '''
        :param d: dataframe
        :param e: npy file with umap or tsne output
        :param feature: attribute that is projected on plot
        :param fname: file name to save graph
        :param um: true = projection on umap
        :param fsne: true = projection on FIt-SNE
        :return: plot for a given attribute
        '''
        cm = plt.get_cmap(self.pal1)  # tab20c
        # NUM_COLORS = len(x['single'].unique())
        NUM_COLORS = len(d[str(feature)].unique())
        # y = x['single']
        y = d[str(feature)]
        fig, ax = plt.subplots(figsize=(8, 8))
        ax.set_prop_cycle('color', [cm(1. * i / NUM_COLORS) for i in range(NUM_COLORS)])
        groups = pd.DataFrame(e, columns=['e', 'y']).assign(category=y).groupby('category')
        for name, points in groups:
            ax.scatter(points.e, points.y, marker=',', s=1, label=name)
        # fig.patch.set_variable(False)
        ax.axis('off')
        # ax.set_title('VDJDB: J-2 Epitope Species', weight= 'bold')
        ax.legend(bbox_to_anchor=(1.04, 0), loc='lower left',
                  markerscale=4)  # bbox_to_anchor=(1.04,0),#the bbanchor & loc is only needed for when lgend = too many
        if um:
            fig.savefig(str(fname) + '_umap_{}.png'.format((str(feature))), bbox_inches='tight',
                        dpi=fig.dpi)  # , dpi=fig.dpi)
        if fsne:
            fig.savefig(str(fname) + '_tsne_{}_.png'.format((str(feature))), bbox_inches='tight', dpi=fig.dpi)
        plt.close()  # fig.savefi

    def plots_per_j(self, um=False, fsne=False):
        '''
        :param um: true = umap
        :param fsne: true = fit-sne
        :return: a plot per j-gene for each of the attributes on VDJdb
                both require dataset to be divided into subfiles (npy and csv) filtered per j gene
        '''
        for n in range(1, 8):
            if fsne:
                e = np.load('tsne_all{}.npy'.format(str(n)))
                d = pd.read_csv('j{}_allP.csv'.format(str(n)))
                d['len'] = d['CDR3'].apply(lambda i: len(i))
                # d['hs'] = d['Epitope species'].apply(lambda i: 1 if i == 'HomoSapiens' else 0)
            if um:
                v = np.load('j{}_allP.npy'.format(str(n)))
                d = pd.read_csv('j{}_allP.csv'.format(str(n)))
                d['len'] = d['CDR3'].apply(lambda i: len(i))
                d['hla_ag'] = d['MHC A'].apply(lambda row: anseq.aggr_hla(row))
                e = umap.UMAP(min_dist=self.dist, metric='cosine').fit_transform(v)
            cols = ['hla_ag', 'Epitope species', 'Epitope', 'aggr_v', 'len']
            for i in cols:
                self.projection_plots(d, e, i, '{}'.format(str(n)), um=True)

    def plots_per_epispe(self, e, ump=False, ft=False):
        d = self.vdjdb.rename({'DENV3/4_aggr_v': 'DENV34_aggr_v'}, axis=1)
        for fts in self.feat_list:
            for i in d.columns.tolist():
                if ft:
                    if '_{}'.format(fts) in i:
                        self.projection_plots(d, e, i, 'vdjdb', fsne=True)
                if ump:
                    if '_{}'.format(fts) in i:
                        self.projection_plots(d, e, i, 'vdjdb', um=True)

    def plots_emerson(self, e):
        for i in self.feat_emerson:
            self.projection_plots(self.emerson, e, i, 'emerson', fsne=True)




