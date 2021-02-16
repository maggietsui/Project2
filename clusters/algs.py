import pandas as pd
import numpy as np
import random
import copy
from scipy.special import comb

"""
algs.py
====================================
Classes and methods for project 2
"""
class Ligand():
    """Class that stores data about ligands including LigandID, score, SMILES, on bits"""
    def __init__(self, ID, score, SMILES, onbits):
        """
        Initialize ligand object
        Parameters
        ---------
        ID
            Path to csv with ligand info
        score
        SMILES
        onbits
        """
        self._ID = ID
        self._score = score
        self._SMILES = SMILES
        self._onbits = onbits

    @property
    def ID(self):
        return self._ID
    
    @ID.setter
    def ID(self, ID):
        self._ID = ID

    @property
    def score(self):
        return self._score
    
    @score.setter
    def score(self, score):
        self._score = score

    @property
    def SMILES(self):
        return self._SMILES
    
    @SMILES.setter
    def SMILES(self, SMILES):
        self._SMILES = SMILES

    @property
    def onbits(self):
        return self._onbits
    
    @onbits.setter
    def onbits(self, onbits):
        self._onbits = onbits

   
class Clustering():
    """Base class for clustering"""
    def __init__(self):
        """
        Blah blah blah.
        Parameters
        ---------
        name
            A string to assign to the `name` instance attribute.
        """

    def get_ligands(self, n):
        """
        Gets the first n ligands from the csv
        Parameters
        ---------
        n
            Number of ligands to return

        Returns: dictionary where keys are ligand IDs
        and values are ligand objects
        """
        table = pd.read_csv("./ligand_information.csv")
        ligands = {}
        for i in range(n):
            onbits = list(map(int, table.iloc[i]['OnBits'].split(",")))
            score = table.iloc[i]['Score']
            smiles = table.iloc[i]['SMILES']
            lig = Ligand(i,score,smiles,onbits)
            ligands[i] = lig

        return ligands

    def calculate_distance(self, A, B):
        """
        The jaccard index takes the intersect of onbits over
        the union of onbits. Distance here is 1 - jaccard index 
        Parameters
        ---------
        A
            list of onbits for first ligand
        B
            list of onbits for second ligand

        Returns: distance between two ligands
        """
        return 1 - len(set(A) & set(B))/len(set(A+B))
    
    def distance_matrix(self, ligands):
        """
        Computes distance matrix for a set of ligands
        Parameters
        ---------
        ligands
            A dict where keys are ligandIDs and values
            are ligands
        Returns: distance matrix
        """
        # Initialize distance matrix
        dist = np.full((len(ligands),len(ligands)), float("inf"))

        # loop through matrix and fill in distances
        for i in range(len(ligands)):
            for j in range(len(ligands)):
                if i == j:
                    continue
                dist[i,j] = self.calculate_distance(ligands[i].onbits,ligands[j].onbits)
        
        return dist
    
    def silhouette_score(self, ligands, clusters):
        """
        Calculates average silhouette score for a set of clustered
        ligands to assess quality. 
        Parameters
        ---------
        clusters
            A list of clusters to calculate the score for
        Returns: Average silhouette score
        """
        # Initialize distance matrix
        dist_mat = self.distance_matrix(ligands)
        scores = []
        
        # loop through each ligand in each cluster
        for curr_cluster in range(len(clusters)):
            for ligand in clusters[curr_cluster]:
                intra_dists = []
                mean_inter_dists = []
                
                # mean distance betw ligand and all members of same cluster
                for member in clusters[curr_cluster]:
                    if member == ligand:
                        continue
                    intra_dists.append(dist_mat[ligand,member])
                avg_intra = mean(intra_dists)
                
                # mean distance betw ligand and members of the nearest cluster
                # loop through all other clusters
                for other_cluster in range(len(clusters)):
                    if other_cluster == curr_cluster: # skip current cluster
                        continue
                    inter_dists = []
                    # calculate dist betw ligand and other cluster members
                    for member in clusters[other_cluster]:
                        inter_dists.append(dist_mat[ligand,member])
                    
                    mean_inter_dists.append(mean(inter_dists))
                    
                avg_inter = min(inter_dists) # get mean dist for nearest cluster
                score = (avg_inter - avg_intra) / max(avg_inter, avg_intra)
                scores.append(score)
        
        return mean(scores)
    
    def rand_index(self, ligands,clustering1, clustering2):
        """
        Rand index measures similarity between two sets of clusters.
        The Rand index has a value between 0 and 1, with 0 indicating
        that the two data clusterings do not agree on any pair of points
        and 1 indicating that the data clusterings are exactly the same.
        Parameters
        ---------
        ligands
            dict of ligands
        clustering1
            First set of clusters
        clustering2
            Second set of clusters
        Returns: Rand index for the two sets of clusters
        """
        # Convert list of clusters to list of counts
        c1 = [0] * len(ligands)
        for i in range(len(clustering1)):
            for j in range(len(clustering1[i])):
                idx = clustering1[i][j]
                c1[idx] = i

        c2 = [0] * len(ligands)
        for i in range(len(clustering2)):
            for j in range(len(clustering2[i])):
                idx = clustering2[i][j]
                c2[idx] = i

        # number of different pairs for each clustering
        combos_c1 = comb(np.bincount(c1),2).sum()
        combos_c2 = comb(np.bincount(c2),2).sum()

        # combine counts into one mat
        mat = np.c_[(c1, c2)]

        # go through mat and find matching pairs (same
        # in both clusterings)
        F_11 = 0
        for i in range(len(clustering1)):
            F_11 += comb(np.bincount(mat[mat[:, 0] == i, 1]), 2).sum()
        F_01 = combos_c1 - F_11
        F_10 = combos_c2 - F_11
        F_00 = comb(len(mat), 2) - F_11 - F_10 - F_01
        return (F_11 + F_00)/(F_11 + F_01 + F_10 + F_00)
    
class HierarchicalClustering(Clustering):
    """Implementation of HC using single linkage"""
    def __init__(self):
        """
        Class that implements hierarchical clustering using single
        linkage
        Parameters
        ---------
        n_clusters
            Default 1. Stops clustering when algorithm reaches
            n_clusters
        """
        super().__init__()

    def cluster(self, ligands, k=1):
        """
        Method that takes a set of ligands and clusters them
        Parameters
        ---------
        ligands
            dict where keys are ligand IDs and values are ligands
        k
            number of clusters to stop at
        Returns: Dictionary of clusters with assigned ligands
        """
        # Initialize cluster list
        clusters = []
        for key in ligands.keys():
            clusters.append([key]) # put all ligands in their own cluster

        # Initialize distance matrix
        dist = self.distance_matrix(ligands)
        
        while len(clusters) > k:
            # find min 
            min_i, min_j = np.unravel_index(np.argmin(dist), dist.shape)

            # find the two ligands in clusters
            found_i = -1
            found_j = -1
            for idx in range(len(clusters)):
                #if found_i >= 0 and if found_j >= 0:
                 #   break
                if found_i < 0 and min_i in clusters[idx]:
                    found_i = idx
                if found_j < 0 and min_j in clusters[idx]:
                    found_j = idx

            # merge the clusters
            clusters[found_i]+= clusters[found_j]
            del clusters[found_j]


            # update distance matrix using single linkage
            for idx in range(len(ligands)):
                if dist[min_i, idx] == np.inf: #
                    continue
                minimum = min(dist[min_i, idx], dist[min_j,idx])
                dist[min_i, idx] = minimum
                dist[idx, min_i] = minimum

            dist[min_j, :] = np.inf
            dist[:, min_j] = np.inf
        
        return clusters


class PartitionClustering(Clustering):
    """Implementation of partition clustering (kmeans)"""
    def __init__(self):
        """
        Blah blah blah.
        Parameters
        ---------
        name
            A string to assign to the `name` instance attribute.
        """
        super().__init__()

    def cluster(self, ligands, k):
        """
        Method that takes a set of ligands and clusters them
        Parameters
        ---------
        ligands
            dict where keys are ligand IDs and values are ligands
        k
            number of clusters to initialize
        Returns: List clusters with assigned ligands
        """
        # choose k random ligands as centroids
        clusters = random.sample(range(len(ligands)), k)
        centroids = []
        for i in range(k):
            centroids.append(ligands[clusters[i]].onbits)
            clusters[i] = [clusters[i]]
        
        #same_clusters = False
        #same_iterations = 0
        old_clusters = copy.deepcopy(clusters)
        iters = 0

        # Recompute centroids until they clusters haven't changed for 2 iterations
        #while same_clusters == False:
        while iters < 300:
            new_clusters = copy.deepcopy(clusters)
            new_centroids = copy.deepcopy(centroids)
            for ligand in ligands.keys(): # Loop through ligands
                distances = []
                # Compute distance to each centroid
                for centroid in centroids:
                    dist = pc.calculate_distance(ligands[ligand].onbits, centroid)
                    distances.append(dist)

                # assign ligand to cluster with min distance (tie breaking?)
                closest = distances.index(min(distances))
                new_clusters[closest].append(ligand)
                new_clusters[closest] = list(set(new_clusters[closest]))

                # Add onbits to new centroids
                new_centroids[closest].extend(ligands[ligand].onbits)
                new_centroids[closest] = list(set(new_centroids[closest]))

            # recalculate centroids 
            centroids = copy.deepcopy(new_centroids)

            # Keep track of how many times the clusters have not changed
            #if new_clusters == old_clusters:
            #    same_iterations += 1
            #if same_iterations == 10:
            #    same_clusters = True
            iters += 1
            old_clusters = copy.deepcopy(new_clusters)

        return new_clusters
        