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
        pass

    def get_ligands(self, n):
        """
        Gets the first n ligands from the csv

        Parameters
        ---------
        n
            Number of ligands to return
        
        Returns
        ---------
        dictionary where keys are ligand IDs and values are ligand objects
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

        Returns
        ---------
        distance between two ligands
        """
        return 1 - len(set(A) & set(B))/len(set(A+B))
    
    def distance_matrix(self, ligands):
        """
        Computes distance matrix for a set of ligands

        Parameters
        ---------
        ligands
            A dict where keys are ligandIDs and values are ligands

        Returns
        ---------
        distance matrix
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
    
    def silhouette_score(self, ligands, clustering):
        """
        Calculates average silhouette score for a set of clustered
        ligands to assess quality. 

        Parameters
        ---------
        clusters
            A list of clusters to calculate the score for

        Returns
        ---------
        Average silhouette score
        """
        # Initialize distance matrix
        dist_mat = self.distance_matrix(ligands)
        scores = []
        
        clusters = []
        for i in set(clustering):
            clusters.append([j for j in range(len(clustering)) if clustering[j] == i])

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
                avg_intra = np.mean(intra_dists)
                
                # mean distance betw ligand and members of the nearest cluster
                # loop through all other clusters
                for other_cluster in range(len(clusters)):
                    if other_cluster == curr_cluster: # skip current cluster
                        continue
                    inter_dists = []
                    # calculate dist betw ligand and other cluster members
                    for member in clusters[other_cluster]:
                        inter_dists.append(dist_mat[ligand,member])
                    
                    mean_inter_dists.append(np.mean(inter_dists))
                    
                avg_inter = min(inter_dists) # get mean dist for nearest cluster
                score = (avg_inter - avg_intra) / max(avg_inter, avg_intra)
                scores.append(score)
        
        return np.mean(scores)
    
    def rand_index(self, clustering1, clustering2):
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

        Returns
        ---------
        Rand index for the two sets of clusters
        """
        # number of different pairs for each clustering
        combos_c1 = comb(np.bincount(clustering1),2).sum()
        combos_c2 = comb(np.bincount(clustering2),2).sum()

        # combine counts into one mat
        mat = np.c_[(clustering1, clustering2)]

        # go through mat and find matching pairs (same
        # in both clusterings)
        F_11 = 0
        for i in range(len(clustering1)):
            F_11 += comb(np.bincount(mat[mat[:, 0] == i, 1]), 2).sum()
        F_01 = combos_c1 - F_11 # pairs that are together in c1
        F_10 = combos_c2 - F_11 # pairs that are together in c2
        # F_00 is the rest of the pairs that are different in both
        F_00 = comb(len(mat), 2) - F_11 - F_10 - F_01
        return (F_11 + F_00)/(F_11 + F_01 + F_10 + F_00)
    
class HierarchicalClustering(Clustering):
    """Implementation of HC using single linkage"""
    def __init__(self):
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

        Returns
        ---------
        Dictionary of clusters with assigned ligands
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
        
        c = [0] * len(ligands)
        for i in range(len(clusters)):
            for j in range(len(clusters[i])):
                idx = clusters[i][j]
                c[idx] = i
                
        return c


class PartitionClustering(Clustering):
    """Implementation of partition clustering (k++ means)"""
    def __init__(self):
        super().__init__()

    def cluster(self, ligands, k):
        """
        Method that takes a set of ligands and clusters them
        using k++ means
        Parameters
        ---------
        ligands
            dict where keys are ligand IDs and values are ligands
        k
            number of clusters to initialize

        Returns
        ---------
        List clusters with assigned ligands
        """
        # initialize centroids with k++ method
        d = self.distance_matrix(ligands)
        centroids = []
        centroids = []
        clusters = [-1] * len(ligands)
        chosen = []
        # start with equal weights for all ligands
        weights = np.ones(len(ligands))

        for i in range(k):
            lig = random.choices(range(len(ligands)), weights)[0]
            # prevent it from choosing the same centroids again
            while(lig in chosen == True):
                lig = random.choices(range(len(ligands)), weights)[0]
            chosen.append(lig)
            clusters[lig] = i
            centroids.append(ligands[lig].onbits)
            # update weights to the distance between lig and all other ligs
            weights = copy.deepcopy(d[lig,:])
            weights[np.where(weights == float("inf"))] = -10 # downweight the same ligand

        print(clusters)
        iters = 0

        # Recompute centroids until a certain number of iterations is reached
        while iters < 100:
            if iters > 0:
                clusters = [0] * len(ligands)
            for ligand in ligands.keys(): # Loop through ligands               
                distances = []
                # Compute distance to each centroid
                for centroid in centroids:
                    dist = self.calculate_distance(ligands[ligand].onbits, centroid)
                    distances.append(dist)

                # assign ligand to cluster with min distance (tie breaking?)
                closest = distances.index(min(distances))
                clusters[ligand] = closest
                
                # Add onbits to new centroids
                centroids[closest].extend(ligands[ligand].onbits)

            # recalculate centroids 
            # compute which onbits are on for at least half of the ligands
            # assign those onbits as new centroids
            new_centroids = []
            for i in range(len(centroids)):
                bits = list(set(centroids[i]))
                new_bits = []
                for bit in bits:
                    if (centroids[i].count(bit)) / clusters.count(i) >= 0.5:
                        new_bits.append(bit)
                new_centroids.append(new_bits)
                        
            
            centroids = new_centroids
            iters += 1
                
        # Check for empty clusters
        for i in range(k):
            # if empty, fill in a random ligand
            if i not in clusters:
                lig = random.randint(0,len(ligands))
                clusters[lig] = i
        return clusters
        