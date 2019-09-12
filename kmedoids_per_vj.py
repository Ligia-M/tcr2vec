
import glob
import matplotlib.pyplot as plt
import numpy as np
from sklearn import cluster
from sklearn import metrics


# vj_file = np.load("/home/ligia/Desktop/Emerson2017/trial/TCRBV01-01TCRBJ01-01.npy")
# vj_file_ix = np.load("/home/ligia/Desktop/files/TCRBV01-01TCRBJ01-01_index.npy")

class Kmeans_VJ:
    def __init__(self, pathdirvjcombs):
        '''
        :param pathdirvjcombs: directory specified in VJ matrix class (outbound loc for files)
        '''
        self.pathdirvjcombs = pathdirvjcombs

    def getsilscores(self):
        '''
        :return: generate kmeans clustering and find optimal cluster number with silhouette score
        '''
        all_sil = []
        all_files =[]
        # for i in glob.glob('/home/ligia/Desktop/Emerson2017/kmeans_sil/*npy'):
        for i in sorted(glob.glob('{}*.npy'.format(self.pathdirvjcombs))):
            file = np.load(i)
            all_files.append(i)
            notwork=[]
            check=[]
        # for i in f:
            silhouette_values = []
            k_values = range(2, 10)
            for k in k_values:
                try:
                    kmeans = cluster.KMeans(n_clusters=k)
                    kmeans.fit(file)
                    labels = kmeans.labels_
                    # print(labels)
                    centroids = kmeans.cluster_centers_
                    # print(centroids)
                    silhouette_score = metrics.silhouette_score(file, labels, metric='euclidean')
                    silhouette_score = silhouette_score.mean()
                    #print(silhouette_score)
                    silhouette_values.append(silhouette_score)
                except:
                    notwork.append(file)
            all_sil.append(silhouette_values)
        return all_sil, all_files, k_values

    def plot_silscores(self):
        '''
        :return: plots silhouette score per vj combination
        '''
        all_sil, all_files, k_values = self.getsilscores()
        check=[]
        ix_notwork = []
        for v, f in zip(all_sil, all_files):
            try:
                plt.plot(k_values, v, 'b-')
                plt.xlabel('k')
                plt.ylabel('silhouette score')
                plt.ylim([0, 1])
                plt.show()
                #n =+1
                plt.savefig('{0}.png'.format(str(f[43:-4])))
                plt.close()
            except:
                ix_notwork.append(all_sil.index(v))
                check.append(f[43:-4])

    def get_clust_ix(self):
        '''
        :return: indexes for which combination has n clusters
        '''
        all_sil, all_files, k_values = self.getsilscores()
        twoclst = []
        threeclst = []
        fourclst = []
        fiveclst = []
        sixclst = []
        sevenclst = []
        eightclst = []
        nineclst = []
        check=[]
        for i in all_sil:
            try:
                if i[0] == max(i):
                    twoclst.append(all_sil.index(i))
                if i[1] == max(i):
                    threeclst.append(all_sil.index(i))
                if i[2] == max(i):
                    fourclst.append(all_sil.index(i))
                if i[4] == max(i):
                    sixclst.append(all_sil.index(i))
                if i[3] == max(i):
                    fiveclst.append(all_sil.index(i))
                if i[5] == max(i):
                    sevenclst.append(all_sil.index(i))
                if i[6] == max(i):
                    eightclst.append(all_sil.index(i))
                if i[7] == max(i):
                    nineclst.append(all_sil.index(i))
            except:
                check.append(i)
        return twoclst, threeclst, fourclst, fiveclst, sixclst, sevenclst, eightclst, nineclst, check

    def run_all(self):
        print('Generating silhouette scores...')
        all_sil, all_files, k_values = self.getsilscores()
        print('Generating plots...')
        self.plot_silscores()
        print('Generating indexes per cluster type...')
        self.get_clust_ix()
        print('Done.')


