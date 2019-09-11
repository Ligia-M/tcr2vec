#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 17 10:59:44 2019

@author: ligia
"""

import os
import zipfile
import glob
import numpy as np
import pandas as pd

from code.generate_embeddings import W2vProcessing


class VJmatrix:
    def __init__(self, vjfilename, vecweightfile, all_zips):
        '''
        :param vjfilename: str, outbound file name (ex. 'vj_combinations')
        :param vecweightfile: npy file, the vector file for all sequences generated in "generate_embeddings.py"
        :param all_zips: location of folder with all data zip files
        '''
        self.vjfilename = vjfilename
        self.vecweightfile = vecweightfile
        self.all_zips = all_zips

    def create_main_file(self):
        '''
        :return: returns an array with all VJ combinations
        '''
        main_file_temp = []
        for zip_file in os.listdir(self.all_zips):
            zip = zipfile.ZipFile(str(self.all_zips) + str(zip_file))
            files_in_zip = zip.infolist()
            for file in files_in_zip:
                #            retrieved_data = pd.DataFrame.from_csv(zip.open(file), sep = '\t')
                retrieved_data = pd.read_csv(zip.open(file), sep='\t', usecols=['sample_name',
                                                                                'amino_acid', 'v_gene', 'j_gene'])
                main_file_temp.append(retrieved_data[['sample_name', 'amino_acid', 'v_gene', 'j_gene']])
                main_file = pd.DataFrame()
                main_file = pd.concat(main_file_temp, axis=0, ignore_index=True)
                main_file = main_file.dropna().reset_index(drop=True)
                main_file = main_file[~main_file.duplicated(subset=['amino_acid', 'sample_name'], keep='first')]
                main_file = main_file[~main_file.amino_acid.str.contains('\*')]
                main_file = main_file.drop(['amino_acid', 'sample_name'], axis=1)
                main_file['vj'] = main_file[['v_gene', 'j_gene']].apply(lambda x: ''.join(x), axis=1)
                main_arr_temp = np.array(main_file['vj'].values.tolist())
                main_arr = main_arr_temp.reshape(main_arr_temp.shape[0], -1)
                np.save(self.vjfilename, main_arr)
            return main_arr

    def get_dicts(self):
        '''
        :return: generate dictionaries for later functions
        '''
        vj_array = self.create_main_file()
        vj_unique_str, vj_as_int = np.unique(vj_array, return_inverse=True)
        #   to view OG mapping of number 2 class --> cat_as_str[cat_as_int] or cat_as_str[# of specific index]
        #   transform cat_as_int into a index dictionary to map to OG vec weights
        seq_idx_dic = {i: vj_as_int[i] for i in range(0, len(vj_as_int))}
        vec_w = np.load(self.vecweightfile) # ex. 'vec_weights100D_155p_tr330_unique.npy'
        return vec_w, seq_idx_dic, vj_as_int, vj_unique_str

    def get_ixdict(self, data):
        '''
        :return: turns vj_as_int into a dictionary with indexes
        '''
        vj_ix_dic = {}
        result = []
        # z = []
        for ix, element in enumerate(data):
            if element not in vj_ix_dic:
                target = vj_ix_dic[element] = [ix]
                result.append(target)
                # print(ix)
            else:
                vj_ix_dic[element].append(ix)  # this works
                # z.append(ix)
        val_only = [v for v in vj_ix_dic.values()]
        key_only = [k for k, v in vj_ix_dic.items()]
        return vj_ix_dic, val_only, key_only

    def generate_files(self, vjfilesloc):
        '''
        :param vjfilesloc: str, location of where you want to save all VJ combination files
        :return: groups sequences in the embedding space per VJ combination --> each file represents 1 combination
        '''
        vec_w, seq_idx_dic, vj_as_int, vj_unique_str = self.get_dicts()
        vj_ix_dic, val_only, key_only = self.get_ixdict(vj_as_int)
        grp_by_vj = [vec_w[i] for i in val_only]
        for arr, vj in zip(grp_by_vj, key_only):
            np.save(str(vjfilesloc)+'{0}.npy'.format(str(vj_unique_str[vj])), arr)
        print('Done.')
        print('Sanity check for data shape')
        for i in glob.glob('{}*.npy'.format(vjfilesloc)):
            s = np.load(i)
            print(s.shape)

