
import pandas as pd
from scipy.stats import stats
import statsmodels.stats.multitest as multi
import numpy as np
import itertools


class FishAnalyses:
    def __init__(self, d):
        self.d = d

    def get_bin(self, row, feat_val):#, ep):# number_val=False): #, epitope=False):
        '''
        :param row: current row function is iterating over
        :param feat_val: int or str, value of a feature ex) in aggr_j it will be int 1 - 7
        :return: encode each feature value as binary (1 = present, 0 = anything else)
        '''
        if row == feat_val:
            return 1
        else:
            return 0

    def gener_bincols(self, colname, newname):
        '''
        :param colname: column of interest present in d dataframe
        :param d: the dataframe you want to get binary columns
        :param newname: new name you give to the column (ex. bin_v or bin_len)
        :return: create a new binary column for each unique value of a column of interest
                the name will be ex. bin_v_5 --> all instances of a V family == 5 will be 1, else 0
        '''
        for i in self.d[str(colname)].unique():
            self.d[str(newname)+'{}'.format(str(i))] = self.d[str(colname)].apply\
                (lambda row: self.get_bin(row, i))

    def get_cols(self, colname):
        '''
        :param colname: each column with 'bin'
        :return: return a list of all unique columns that are values of a feature
            ex. bin_j =  ['bin_j7', 'bin_j1', 'bin_j2', 'bin_j5', 'bin_j3', 'bin_j4', 'bin_j6']
        '''
        listname = []
        for i in self.d:
            if str(colname) in i:
                listname.append(i)
        return listname

    '''function to filter data'''
    def filter_by(self, var1, var2, filtby2=False, filtby1=False):
        '''
        :param var1: str, feature name #1 to filter data with (ex. 'aggr_j')
        :param var2: str,  feature #2 to filter data with
        :param filtby2: subsample d by both vars
        :param filtby1: subsample d by 1 var --> the default that I used
        :return: list of lists where each list is a dataframe subsampled for a specific value
                ex. [[1], [2], [etc]] --> 1 is all cells in df that have aggr_j==1
        '''
        dfs = []
        if filtby2:
            for v1 in self.d[str(var1)].unique():
                for v2 in self.d[str(var2)].unique():
                    dfs.append(self.d[(self.d[str(var1)] == v1) & (self.d[str(var2)] == v2)])
        if filtby1:
            for v1 in self.d[str(var1)].unique():
                dfs.append(self.d[self.d[str(var1)] == v1])
        return dfs

    '''Compute Fisher analysis'''
    def get_fish(self, varlist, var3):
        '''
        :param varlist: list with variables outputted from get_cols()
        :param var3: additional variable for computing fisher's exact test (FET)
                    aggr_j is default --> will compute most global (FET)
        :return: computes both fine-grain FET (for a given J, for a given feature value in varlist) &
                coarse-grain FET (present/absent value of a feature in a given J
                Compute global FDR correction after all FETs
                return variables of interest:
                    fail_eth = feature values without enough observations to compute xtabs
                    pval_eth = total p_values computed from all FET
                    nrch = enriched combinations after FDR correction
                    enriched_count = count of these nrch combinations
        '''
        dfs = self.filter_by('aggr_j', 'aggr_v', filtby1=True)
        pval_eth = []
        cn_eth = []
        xt_eth = []
        counting=0
        fail_eth=[]
        for ss in dfs:
            xt = []
            p_values = []
            sig = []
            combi_name = []
            fail = []
            counting+=1
            for allvars in varlist:
                for v1 in allvars:
                    for var2 in varlist:
                        for v2 in var2:
                            if allvars == var2:
                                continue
                            xtab = pd.crosstab(ss[v1], ss[v2], margins=False)
                            if xtab.shape != (2, 2):
                                fail.append(tuple([str(counting) + v1, v2]))
                                continue
                            xt.append(xtab)
                            combi_name.append(tuple([str(counting) + v1, v2]))
                            r, p = stats.fisher_exact(xtab, alternative='two-sided')
                            if p <= .05:
                                sig.append(xtab)
                            p_values.append(p)
            pval_eth.append(p_values)
            cn_eth.append(combi_name)
            xt_eth.append(xt)
            fail_eth.append(fail)
        j_pvalues = []
        j_combiname = []
        j_fail = []
        j_xtabs = []
        for v3 in var3:
            for vars in varlist:
                for v in vars:
                # for v2 in var2:
                #     # try
                    xtab = pd.crosstab(d[v3], d[v], margins = False)
                    if xtab.shape!= (2,2):
                        j_fail.append(tuple([v3, v]))
                        continue
                    j_xtabs.append(xtab)
                    j_combiname.append(tuple([v3, v]))
                    r, p = stats.fisher_exact(xtab, alternative='two-sided')
                #     if p <= .05:
                #         sig.append(xtab)
                #     p_values.append(p)
                    j_pvalues.append(p)
        pval_eth.append(j_pvalues)
        cn_eth.append(j_combiname)
        xt_eth.append(j_xtabs)
        fail_eth.append(j_fail)
        pval_eth= list(itertools.chain(*pval_eth))
        cn_eth = list(itertools.chain(*cn_eth))
        xt_eth = list(itertools.chain(*xt_eth))
        fail_eth = list(itertools.chain(*fail_eth))
        r, corr, sidak, bonfer = multi.multipletests(pval_eth,.05, method='fdr_bh')
        rejected_xtabs_ix = np.where(r == True)
        rej_xtabs = [xt_eth[i] for i in rejected_xtabs_ix[0]]
        combi_name =[cn_eth[i] for i in rejected_xtabs_ix[0]]
        enriched_count=0
        nrch = []
        for rej, cn in zip(rej_xtabs, combi_name):
            try:
                ath = rej[1][0]/ rej[0][0]
                j = rej[1][1] / rej[0][1]
                print('Other:{}'.format(ath) , 'J:{}'.format(j))
                if j>ath:
                    print('ENRICHED!')
                    enriched_count += 1
                    nrch.append(cn)
                else:
                    print('NOT ENRICHED')
            except:
                print(cn)
        return rej_xtabs, pval_eth, r, nrch, enriched_count, combi_name, fail_eth, sidak, bonfer
