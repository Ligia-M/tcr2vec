#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 28 20:28:05 2019

@author: ligia
"""


import pandas as pd
from cos_antigen_vdjdb import CosinePathogenAnalysis

#%%
# seqs = pd.read_csv('SearchTable-2019-05-28 21_49_18.859.tsv', sep ='\t')
# seqs = seqs[seqs['Species']=='HomoSapiens']


class AnnotateSeqs:
    def __init__(self, vdjdb, emerson):
        self.vdjdb = vdjdb
        self.emerson = emerson
        self.feat_list = ['aggr_v', 'aggr_j', 'len', 'Epitope species', 'Epitope', 'hla_ag']

    # self.df_features
    #generate a column that tells me which amino acid is private vs shared
    #true means private false is shared
    def fill_rows(self, row):
        if row['dup'] == False and row['Virus Diseases'] == 'Cytomegalovirus -':
            return 'Negative status'
        if row['dup'] == False and row['Virus Diseases'] == 'Cytomegalovirus +':
            return 'Positive status'
        if row['dup'] == True:
            return 'Shared'

    def cmv_cat(self):
        '''
        iterate through data

        if amino acid has no duplicate & cmv == '- cmv '--> third column = 0

        ifelse amino acid no dupli & positive --> 1

        if duplicate --> 2 '''
        duplicates = self.emerson.duplicated(subset=['amino_acid'], keep='first')
        self.emerson['dup'] = duplicates
        self.emerson['cmv_newlabel'] = self.emerson.apply(lambda row: self.fill_rows(row), axis=1)
        return self.emerson

    # df_features['hla'] = df_features['HLA MHC class I'].str[4:8]
    #
    ##st = st.set_index('HLA')
    # st['hla'] = st['HLA'].str[0:4]
    #
    # df_features.join(st.set_index('hla'), on='hla')

    def priv_share(self):
        df2 = self.emerson.groupby('amino_acid').apply(lambda x: x['sample_name'].unique())
        df2 = df2.reset_index().rename({0: 'subjs'}, axis=1)
        merged = self.emerson.merge(df2, how='outer')

        def all_same(row):
            return all(x == row['subjs'][0] for x in row['subjs'])
        merged['sh_pr'] = merged.apply(lambda row: all_same(row), axis=1)
        return merged

    ''' extract features from description "sample_tags"
        returns a dataframe of all features of interest '''
    def get_emerson_features(self, cmv_new=True, priv_share=True, len=True):
        # specify features of interest with separate variable for features in textual description
        # feats_og = data[['v_family', 'amino_acid', 'v_gene', 'j_family','j_gene', 'cdr3_length', 'sample_name']]
        feats_og = self.emerson[['v_family', 'amino_acid', 'j_gene', 'sample_name']]
        feats_description = self.emerson['sample_tags']
        # extracting features from description
        split_comma = feats_description.apply(lambda x: x.split(','))
        split_colon = [[tuple(pair.split(':')) for pair in row] for row in split_comma]
        df_features = pd.DataFrame([{x[0]: x[1] for x in row} for row in split_colon])
        # concatenate all features into single dataframe
        df_features = pd.concat([feats_og.reset_index(drop=True), df_features.reset_index(drop=True)], axis=1)
        st = pd.read_csv('supertype2.csv')
        st = st.dropna()
        df_features['hla'] = df_features['HLA MHC class I'].str[4:8]
        # st = st.set_index('HLA')
        st['hla'] = st['HLA'].str[0:4]
        st = st.drop_duplicates('hla', keep='first')
        st = st.drop('HLA', axis=1)
        df_features = df_features.merge(st, how='outer')
        df_features = df_features[0:1215354]
        # adding feature
        if cmv_new:
            df_features = self.cmv_cat()
        if priv_share:
            df_features = self.priv_share()
        if len:
            self.get_len(vdj=False, emrsn=True)
        cpa = CosinePathogenAnalysis('/home/ligia/Desktop/Emerson2017/','vec_weights100D_allP_vdjdb.npy')
        cpa.aggr_cols(df_features, 'v_family', 'j_gene', main=True)
        #    df_features = df_features.drop('amino_acid', axis = 1)
        return df_features

    def get_len(self, vdj=True, emrsn=True):
        if vdj:
            self.vdjdb['len'] = self.vdjdb['CDR3'].apply(lambda i: len(i))
        if emrsn:
            self.emerson['len'] = self.emerson['amino_acid'].apply(lambda i: len(i))

    ''' get a column with which you can plot all epitopes related to
    a pathogen while everything else is tagged as 'other' '''
    def get_feature_per_pathogen(self, pathogen, feature):
        ep_spec = self.vdjdb[self.vdjdb['Epitope species'] == str(pathogen)]
        feat_path = ep_spec[str(feature)].unique().tolist()

        def feat_val(row):
            for i in feat_path:
                if row[str(feature)] == i and row['Epitope species'] == str(pathogen):
                    return i
            else:
                return 0
        self.vdjdb['{}_'.format(str(pathogen))+str(feature)] = self.vdjdb.apply(lambda row: feat_val(row), axis =1)

    ''' Aggregate HLA for VDJdb'''
    def aggr_hla(self, row):
        if 'HLA-A' in row:
            return row[0:8]
        if 'HLA-B' in row:
            return row[0:8]
        if 'HLA-DRA' in row:
            return row[0:10]
        if 'HLA-DQA1' in row:
            return row[0:11]
        if 'HLA-DQB1' in row:
            return row[0:11]
        if 'HLA-DRB' in row:
            return row[0:11]
        if 'B2M' in row:
            return row

    def generate_pathfeat_cols(self):
        cpa = CosinePathogenAnalysis('/home/ligia/Desktop/Emerson2017/','vec_weights100D_allP_vdjdb.npy')
        cpa.aggr_cols(self.vdjdb, 'V', 'J', vdjdb=True)
        self.vdjdb['hla_ag'] = self.vdjdb['MHC A'].apply(lambda row: self.aggr_hla(row))
        self.vdjdb = self.get_len(vdjdb=True)
        for fts in self.feat_list:
            for i in self.vdjdb['Epitope species'].unique().tolist():
                self.get_feature_per_pathogen(i, str(fts))
            # print(i)
            # print(d[d['Epitope species'] == str(i)]['aggr_v'].value_counts())
            # print(d.groupby(['Epitope species', 'aggr_v'])['aggr_v'].counts())
        return self.vdjdb

    def generate_labels(self, df_cols):
        encoded_features = pd.DataFrame({col: df_cols[col].astype('category').
                                        cat.codes for col in df_cols}, index=df_cols.index)
        # generate dictionary of mapping between label and feature
        feature_label_dic = {col: {n: cat for n, cat in enumerate(df_cols[col].astype('category').cat.categories)}
                             for col in df_cols}
        # retrieve array of labels per sequence
        feature_arrays = encoded_features.values
        return feature_arrays, feature_label_dic


