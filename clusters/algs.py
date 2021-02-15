import pandas as pd
import numpy as np
import random
import copy

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
        dist = np.full((len(ligands),len(ligands)), float("inf"))

        # loop through bottom half of matrix and fill in distances
        for i in range(len(ligands)):
            for j in range(len(ligands)):
                if i == j:
                    continue
                dist[i,j] = self.calculate_distance(ligands[i].onbits,ligands[j].onbits)

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
        
        same_clusters = False
        same_iterations = 0
        old_clusters = copy.deepcopy(clusters)
        

        # Recompute centroids until they clusters haven't changed for 2 iterations
        while same_clusters == False:
            new_clusters = copy.deepcopy(clusters)
            new_centroids = copy.deepcopy(centroids)
            print(new_clusters)
            print(new_centroids)
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
            if new_clusters == old_clusters:
                same_iterations += 1
            if same_iterations == 4:
                same_clusters = True

            old_clusters = copy.deepcopy(new_clusters)
            
            print(new_clusters)
            print(new_centroids)
        return new_clusters
   