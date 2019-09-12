from annot_seqs import AnnotateSeqs
from cos_antigen_vdjdb import CosinePathogenAnalysis
from lev_dist_compute import LevAnalyses
from fisher_analysis import FishAnalyses
from runnin_cossimv2 import top_k_analyses #the one that doesnt work --> use non class file instead
from cosine_analysis import OutputTopkAnalysis
from data_visualization import DataVis
from getting_generator_sequences import GetGenerators
from generate_embeddings import W2vProcessing
from vj_matrix import VJmatrix
from kmedoids_per_vj import Kmeans_VJ
from tsne import Tsne
import multiprocessing as mp
import numpy as np
import glob


# open_csv has booleans vdjdb = False, ss= False, concat = False
#first argument is file name
#function drops NaNs, reset index
# -> SS can be for taking a subsamble of the emerson dataset as CSV
# -> concat is for the concatenation between 2 CSVs (ex. immunized vs non-immunized)
#second string space is for second file in concat, otherwise leave blank


'''Generate embeddings'''
w2v = W2vProcessing('/home/ligia/Desktop/Emerson2017/')
cpa = CosinePathogenAnalysis('/home/ligia/Desktop/Emerson2017/','vec_weights100D_allP_vdjdb.npy')
beta = cpa.open_csv('SearchTable-2019-05-28 21_49_18.859.tsv', '', vdjdb=True)
vdjdb = w2v.call_vecs_fncts(beta, 'CDR3', 'vecweights_emerson', train=True, exp2=True)
emerson = w2v.call_vecs_fncts(beta, 'amino_acid', 'vecweights_emerson', train=True, exp1=True)



''' Feature generation'''
#initiate annotation class
annotate = AnnotateSeqs(vdjdb, emerson)

#private versus shared columns (private = True / Shared = False)
emerson = annotate.get_emerson_features()
vdjdb = annotate.generate_pathfeat_cols()



''' Data Visualization'''
t = Tsne()
e = t.fit_sne('vec_weights100D_allP_vdjdb.npy', 'vdjdb')
dv = DataVis(emerson, vdjdb, .1) #  3rd argument is minimum distance for clustering (1-cosine similarity of choice)

#   Barplot with frequency as percentages, exp1=emerson dataset / exp2 = vdjdb dataset
dv.barplot_percentages('vdjdb_barplot', 'aggr_j', 'J-Gene', exp1=False , exp2=True)

#heatmap with example of variables
dv.heatmap_per_exp('aggr_j', 'J-Gene', 'aggr_v', 'V-Family','vdjdb_percentages', exp1=False, exp2=True)

dv.plots_per_j(um=True, fsne=False) #set either to True or false for umap or fit-sne

#   returns plot per epitope species for all features in feat_list (feat for pathogen = 1 everything else = 0 )
dv.plots_per_epispe(e, ft=True, ump=False) #set either to True or false for umap or fit-sne

#   make FIt-sne plots emerson dataset
dv.plots_emerson(e)



'''Top-K cosine analysis '''
#   Generate all 844 combinations
vjm = VJmatrix('vj_comb', 'vec_weights100D_155p_tr330_unique.npy', '/home/ligia/Desktop/Emerson2017/temp/')
vjm.generate_files('/home/ligia/Desktop/Emerson2017/trial/') # location for saving the files

#   get vectors & reshape to compare against 844 combinations
x = np.load('330_files_unique_155seqname.npy')
gg = GetGenerators(x, 150)
ext = gg.gen_prep('150flu', 'InfluenzaA', cmv=True, pathogen=False) #if cmv=true only change first argument
# ext = np.load('164cmvgenerator.npy')
ext = [vec.reshape(1, -1) for vec in ext]

#   initiate class
all_vj_files = []
for i in sorted(glob.glob('/home/ligia/Desktop/Emerson2017/trial/*.npy')):  # location of vj embeddings are
    vector = np.load(i)
    vct = [vec.reshape(1, -1) for vec in vector]
    all_vj_files.append(vct)

topk = top_k_analyses(ext, all_vj_files) #, 'home/ligia/Desktop/Emerson2017/trial')

#   initiate multi-processing
##   processes = n where n is number of available cores for parallelization
pool = mp.Pool(processes= 15)

#   begin analysis
results = pool.map(topk.get_topk, topk.generator())
np.save('cosresultvj_', results) #save results for later analysis

#   output analysis
cos = np.load('cosresultvj_.npy')

glob_index_allvjs = []
for i in glob.glob('/home/ligia/Desktop/Emerson2017/trial/*.npy'):
    glob_index_allvjs.append(i)

ca = OutputTopkAnalysis(glob_index_allvjs)

#   Plot sequence of interest across 844 VJ combinations
ca.get_topkcosim_plot(cos, 'filename')

#   n largest / n lowest cosine similarities across all analyses
lwst_ix, lwst_cos, lrgst_ix, lrgst_cos, \
mean_lr, sd_lr, mean_lw, sd_lw = ca.cos_stats(cos, 3, 30) # n lowest for each 1 x 844
                                                        #   n largest for each 1 x 844
#   returns counter + plot of the counters
vcounter, jcounter, most_occur_vj_name = ca.get_mostoccur\
    (cos, 'filename', lrgst=True, lwst=False)   #Set lrgst or lwst to true (not both)



'''K-means analysis'''
#   Initiate class
kvj = Kmeans_VJ('/home/ligia/Desktop/Emerson2017/trial/')
twoclst, threeclst, fourclst, fiveclst, sixclst, sevenclst,\
eightclst, nineclst, check = kvj.run_all()



'''Levensthein analysis'''
#   initiate class
lev = LevAnalyses(vdjdb)

#  call functions
#random control
lev.rndm_lev_mat(10, 20, 'aggr_j', 'aggr_v', 1, 4, '', '',filt_2=True)

#   compute levensthein matrix
lev.get_lev('vdjdb_nodups.csv', 'amino_acid', 2, 'farr', 'J2_A', strings_not=True)  #for main
lev.get_lev('vdjdb_nodups.csv', 'CDR3', 2, 'farr', 'J2_vdj', strings_not=True)  #for vdjdb

#   permutation test for random vs target matrices
lev.exact_mc_perm_test(distance_matrix, distance_matrix_rand, 50000)



'''Fisher's exact test analysis'''
#   initiate class
fish = FishAnalyses(vdjdb)

#   initiate global variables
curr_cols = ['aggr_v', 'hla_ag', 'aggr_j', 'Epitope species', 'Epitope']
cnames = ['bin_v', 'bin_mhc', 'bin_j', 'bin_len', 'bin_typ', 'bin_ep']

for col, cnam in zip(curr_cols,cnames):
    fish.gener_bincols(col, vdjdb, cnam)

bin_j = fish.get_cols('bin_j')
bin_v = fish.get_cols('bin_v')
bin_len = fish.get_cols('bin_len')
bin_typ = fish.get_cols('bin_typ')
bin_ep = fish.get_cols('bin_ep')
bin_mhc = fish.get_cols('bin_mhc')

varlist = [bin_v, bin_len, bin_typ, bin_ep, bin_mhc]

#   generate fisher's exact test coarse and fine grain effects w/ global FDR correction
rej_xtabs, pval_eth, r, nrch, enriched_count, combi_name, fail_eth, sidak, bonfer = \
    fish.get_fish(varlist, bin_j)
#uncomment to save in pickle format:
# pickle.dump(pval_eth, open('fish_pvaleth_hlaag_allfeats.p','wb'))