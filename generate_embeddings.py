#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 17 10:59:44 2019

@author: ligia
"""

import os
import pickle
import zipfile

import numpy as np
import pandas as pd
from gensim.models.word2vec import Word2Vec


class W2vProcessing:
    def __init__(self, dir_loc):
        self.fcols = ['sample_name', 'amino_acid', 'sample_tags', 'v_family', 'j_gene']
        self.tcols = ['sample_name', 'amino_acid', 'sample_tags']
        self.all_zips = os.listdir(dir_loc) #folder with zipfiles ex."/home/ligia/Desktop/Emerson2017/temp/"
        self.v_vals = ['TCRBV20', 'TCRBV12', 'TCRBV25', 'TCRBV13']

    def create_main_file(self, dirloc, train_sA_w2v=False, train_all_w2v=True, feat_only=False):
        '''
        :param dirloc: folder with all zip files
        :param train_sA_w2v: if sample A dataset (filter for V-families)
        :param train_all_w2v: if train a W2V model with all dataset
        :param feat_only: extract attributes
        :return: a dataframe without duplicates, stop codons and missing CDR3 values
        '''
        main_file_temp = []
        for zip_file in self.all_zips:
            zip = zipfile.ZipFile(str(dirloc) + str(zip_file))
            files_in_zip = zip.infolist()
            if train_sA_w2v:
                for file in files_in_zip:
                   retrieved_data = pd.read_csv(zip.open(file), sep = '\t', usecols=self.fcols)
                   mask = np.in1d(retrieved_data['v_family'].values, self.v_vals)
                   retrieved_data = retrieved_data[mask]
                   main_file_temp.append(retrieved_data[self.fcols])
                   main_file = pd.DataFrame()
                   main_file = pd.concat(main_file_temp, axis = 0, ignore_index = True)
            if train_all_w2v:
                for file in files_in_zip:
                    #            retrieved_data = pd.DataFrame.from_csv(zip.open(file), sep = '\t')
                    retrieved_data = pd.read_csv(zip.open(file), sep='\t',
                                                 usecols=self.tcols)
                    main_file_temp.append(retrieved_data[self.tcols])
                    main_file = pd.DataFrame()
                    main_file = pd.concat(main_file_temp, axis=0, ignore_index=True)
            if feat_only:
                for file in files_in_zip:
        #            retrieved_data = pd.DataFrame.from_csv(zip.open(file), sep = '\t')
                    retrieved_data = pd.read_csv(zip.open(file), sep = '\t', usecols=self.fcols)
                    main_file_temp.append(retrieved_data[self.fcols])
                    main_file = pd.DataFrame()
                    main_file = pd.concat(main_file_temp, axis = 0, ignore_index = True)
                main_file['v_family'] = [i.replace(str(np.nan),'na') for i in main_file['v_family'].astype(str)]
                main_file['j_gene'] = [i.replace(str(np.nan),'na') for i in main_file['j_gene'].astype(str)]
            main_file = main_file.dropna().reset_index(drop = True)
            main_file = main_file[~main_file.duplicated(subset = ['amino_acid', 'sample_name'], keep = 'first')]
        return main_file

    def processing_training_tokens(self, main_data, cdr3):
        '''
        :param main_data: dataframe
        :param cdr3: str, either 'amino_acid' for emerson and 'CDR3' for vdjdb
        :return: list of lists of sequence split into 3 non-overlapping reading frames
        '''
        processed_4_w2v_training = []
        for x in main_data[cdr3]:
            inner1 = []
            inner2 = []
            inner3 = []
            for i in range(0, len(x), 3):
               inner1.append( "".join(x[i:i + 3]))
            for i in range(1, len(x), 3):
                inner2.append( "".join(x[i:i + 3]))
            for i in range(2, len(x), 3):
                inner3.append( "".join(x[i:i + 3]))
            processed_4_w2v_training.extend([inner1,inner2, inner3])
        return processed_4_w2v_training

    def processing_token_sum(self, main_data, cdr3):
        '''
        :param main_data: dataframe
        :param cdr3: str, either 'amino_acid' for emerson and 'CDR3' for vdjdb
        :return: flattened list of sequence split into 3 non-overlapping reading frames
        '''
        processed_4_w2v_sum = []
        for x in main_data[cdr3]:
            inner1 = []
            inner2 = []
            inner3 = []
            for i in range(0, len(x), 3):
               inner1.append( "".join(x[i:i + 3]))
            for i in range(1, len(x), 3):
                inner2.append( "".join(x[i:i + 3]))
            for i in range(2, len(x), 3):
                inner3.append( "".join(x[i:i + 3]))
            processed_4_w2v_sum.append([inner1,inner2, inner3])
        return processed_4_w2v_sum

    def load_model(self, model_name):
        '''
        :param model_name: bin file where model is saved in
        :return: opens the w2v model file
        '''
        return Word2Vec.load(model_name)

    def save_file(self, filename):
        with open(filename, 'wb') as fp:
            pickle.dump(filename, fp)

    def open_file(self, filename):
        file = open(filename, 'rb')
        opened_file = pickle.load(file)
        return opened_file

    def init_model(self, tokens_training, filename):
        '''
        :param tokens_training: variable from processing_training_tokens function
        :return: initializes, train and save model
        '''
        num_features = 100  # Word vector dimensionality
        min_word_count = 3  # Minimum word count
        num_workers = 8  # Number of threads to run in parallel
        context = 2  # Context window size
        downsampling = 1e-3  # Downsample setting for frequent words
        W2Vmodel = Word2Vec(sentences=tokens_training,
                            sg=1,
                            hs=0,
                            workers=num_workers,
                            size=num_features,
                            min_count=min_word_count,
                            window=context,
                            sample=downsampling,
                            negative=5,
                            iter=6)
        # train model
        W2Vmodel.train(tokens_training, total_examples=len(tokens_training), epochs=10)
        W2Vmodel.save('w2vmodel_{}.bin'.format(filename))
    
    ''' Sum vectors across sequence '''
    def get_sequence_vectors(self, vec_model, expanded_sequences):
        k_mers = np.concatenate(expanded_sequences)  # kmer in text
        kmer_keys = set(vec_model.wv.vocab.keys())  # known kmer
        dim_4_row = np.zeros(vec_model.vector_size, dtype="float32")
        kmer_count = 0
        for kmer in k_mers:
            if kmer in kmer_keys:
                dim_4_row = np.add(dim_4_row, vec_model[kmer])
                kmer_count += 1
        return dim_4_row

    def call_vecs_fncts(self, dirloc, main_data, cdr3, filename, train=True, exp1=False, exp2=False):
        '''
        :param dirloc: folder with all zip files
        :param main_data: dataframe
        :param cdr3: str, either 'amino_acid' for emerson and 'CDR3' for vdjdb
        :param filename: str, name of the outbound vector file and w2v model
        :param train: True = will train w2v from emerson or vdjdb
        :param exp1: True = emerson processing and training
        :param exp2: True = vdjdb processing and training
        :return: call all functions to train W2V model or reconstruct embeddings for sequences
        '''
        if train:
            if exp1:
                print('Extracting files...')
                main = self.create_main_file(dirloc, train_all_w2v=True)
                print('All files read')
                print(main['amino_acid'].str.contains('\*').value_counts()[True])
                # drop sequences with stop codons
                main_data = main[~main.amino_acid.str.contains('\*')]
                print('preprocessing done')
            if exp2:
                main_data = main_data.drop_duplicates('CDR3')
            print(' Begin tokenization...')
            tokens_training = self.processing_training_tokens(main_data, cdr3)
            print('training model ...')
            self.init_model(tokens_training, filename)
        tokens_sum = self.processing_token_sum(main_data, cdr3)
        main_data['tokens_sum'] = [inner_list for inner_list in tokens_sum]
        loaded_model = self.load_model('w2vmodel_100D_allP.bin')
        print('Generating embedding per sequence...')
        main_data['vec_sequence'] = list(map(lambda exp_seq:
                                        self.get_sequence_vectors(loaded_model, exp_seq),
                                        main_data.tokens_sum))
        vec_weights = np.array(list(map(np.array, main_data.vec_sequence)))
        np.save(filename, vec_weights)
        print('Done.')


