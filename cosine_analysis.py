import itertools
from collections import Counter
from heapq import nsmallest, nlargest

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


class OutputTopkAnalysis:
    def __init__(self, glob_index_allvjs):
        # self.cos = cos
        self.glob_index_allvjs = glob_index_allvjs

    def get_dic(self, cos):
        '''
        :param cos: npy file with results from topkanalysis loaded
        :return: returns a list of lists (164 x 844) and a dictionary of same size
                used for later analysis
        '''
        cos_list = cos.tolist()
        dict={}
        for i in cos_list:
            dict[cos_list.index(i)]=i[:]
        return cos_list, dict
    # cos_list, dict = get_dic(cos)
    # cos_list, dict = get_dic(rand)
    #   all cosine computations that used glob were done on the unsorted glob
    #   (so order files are stored ≠ from order it appears in folder) / use below dic

    #plot the cosine similarity across VJ
    def chunked(self, it, size):
        '''
        :param it: dictionary we iterate over
        :param size: int, number of sequences to divide the dictionary in
        :return: skeleton for dividing data to plot subplots into multiple files
        '''
        it = iter(it)
        while True:
            p = tuple(itertools.islice(it, size))
            if not p:
                break
            yield p

    def get_topkcosim_plot(self, cos, filename, n):
        '''
        :param cos: npy file
        :param filename: str, directory + name of file to save
        :param n: int, plot number per page (164 plots with 20 per page was most legible)
        :return: call function chunked and will plot dictionary (in order of presentation)
        '''
        cos_list, dict = self.get_dic(cos)
        z = 0
        for chunk in self.chunked(dict.items(), n):
            fig, ax = plt.subplots(20, 1, figsize=(20,18))
            z += 1
            for a,key, in zip(ax,dict.keys()):
                y=dict[key]
                print(key)
                n=len(y)
                x = np.linspace(1,n,n)
                a.plot(x,y, color='#221631')
                plt.savefig(str(filename)+'_sortedbyvj{0}.png'.format(str(z)))

    def cos_stats(self, cos, nlw, nlr):
        '''
        :param cos: npy file
        :param nlw: int, across all sequences of interest (164 in this case), find n lowest
        :param nlr: int, across all sequences of interest (164 in this case), find n largest
        :return: list of lists: lwst and lrgst indexes (to know name of VJ combination)
                + ist of lists: the actual cosine results
                + lrgst default is 30 (164 x 30) and lwst default is 3 (164 x 3)
        '''
        cos_list, dict = self.get_dic(cos)
        lwst_ix = []
        lwst_cos = []
        lrgst_ix =[]
        lrgst_cos = []
        for seq in cos_list:
            nsmall = nsmallest(nlw, seq)
            nhigh = nlargest(nlr, seq)
            lwstn = []
            lwst_cos = []
            lrgstn = []
            lrgst_cos = []
            for n in nsmall:
                lwstn.append(seq.index(n))
                lwst_cos.append(n)
            for n in nhigh:
                lrgstn.append(seq.index(n))
                lrgst_cos.append(n)
            lwst_ix.append(lwstn)
            lwst_cos.append(lwst_cos)
            lrgst_ix.append(lrgstn)
            lrgst_cos.append(lrgst_cos)
        flat_lst_lr = [val for inner in lrgst_cos for val in inner]
        mean_lr, sd_lr = np.mean(flat_lst_lr), np.std(flat_lst_lr)
        flat_lst_lw = [val for inner in lwst_cos for val in inner]
        mean_lw, sd_lw = np.mean(flat_lst_lw), np.std(flat_lst_lw)
        return lwst_ix, lwst_cos, lrgst_ix, lrgst_cos, mean_lr, sd_lr, mean_lw, sd_lw

    # function to add a space between V and J
    def insert_space(self, string, integer):
        '''
        :param string: str
        :param integer: int
        :return: after n integer in string add space
                use fnction inside get_mostoccur()
        '''
        return string[0:integer] + ' ' + string[integer:]

    def get_mostoccur(self, cos, filename, lrgst=False, lwst=False):
        '''
        :param cos: npy file
        :param filename: str, directory + name of file to save
        :param lrgst: return function for n largest combination
        :param lwst: return function for n lowest combination
        :return: vcounter: V-family only frequency count
                jcounter: J-gene only frequency count
                most_occur_vj_name: combination of VJs frequency count
                plots a heatmap of VJ occurences
        '''
        seq_idx_allvjs = {i: self.glob_index_allvjs[i] for i in range(0, len(self.glob_index_allvjs))}
        lwst_ix, lwst_cos, lrgst_ix, lrgst_cos = self.cos_stats(cos)
        if lrgst:
            flat_lst = [val for inner in lrgst_ix for val in inner]
        if lwst:
            flat_lst = [val for inner in lwst_ix for val in inner]
        add_space = []
        splt_str = []  # a flat list of all most occuring VJ comb using the name
        for i in flat_lst:
            vj_name = seq_idx_allvjs[i][38:-4]
            if 'or' in vj_name:
                add_space.append(self.insert_space(vj_name, 15))
                splt_str.append(self.insert_space(vj_name, 15).split(' '))
            else:
                add_space.append(self.insert_space(vj_name, 10))
                splt_str.append(self.insert_space(vj_name, 10).split(' '))
            splt_str.append(self.insert_space(vj_name, 10).split(' '))
        counter = Counter(add_space)
        most_occur_vj_name = counter.most_common()
        v = []
        j = []
        for i in splt_str:
            if len(i) == 2:
                v.append(i[0][5:7])
                j.append(i[1][-2:])
            # if len(i) == 1:
            #     unresolved.append(i)
        vcounter, jcounter = Counter(v), Counter(j)
        # most_occur_vj_name = counter.most_common()
        # generate a dataframe with the most occuring vj combs + the count for plotting
        df = pd.DataFrame(list(counter.items()), columns=['vj', 'cnt'])
        df[['v', 'j']] = df.vj.str.split(' ', expand=True)
        df['v'] = df['v'].apply(lambda row: row[5:7])
        df['v'] = df['v'].replace('ol', '0')
        df['v'] = df['v'].replace('A-', '0')
        df['j'] = df['j'].apply(lambda row: row[-2:])
        df['j'] = df['j'].replace('ed', '0')
        tab = df.groupby(['j', 'v'])['cnt'].sum()
        tab_percent = tab.groupby(level=0).apply(lambda x: round(100 * x / float(x.sum()))).unstack().fillna(0)
        cmap = sns.cubehelix_palette(light=1, as_cmap=True)
        plt.subplots(figsize=(14, 5))
        ax_hmj = sns.heatmap(tab_percent, annot=True, fmt="g", cmap=cmap, xticklabels=True, yticklabels=True)
        ax_hmj.set(ylabel="J-Gene", xlabel='V-Family')
        ax_hmj.tick_params(left=False, bottom=False)
        plt.savefig('{}.png'.format(str(filename)), bbox_inches='tight')
        plt.close()
        return vcounter, jcounter, most_occur_vj_name


