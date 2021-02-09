import pandas as pd

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
	def __init__(self):
		"""
	    Blah blah blah.
	    Parameters
	    ---------
	    name
	        A string to assign to the `name` instance attribute.
	    """
	
	def get_ligands(n):
		"""
	    Gets the first n ligands from the csv
	    Parameters
	    ---------
	    n
	        Number of ligands to return

	    Returns: dictionary where keys are ligand IDs
	    and values are ligand objects
	    """
	    table = pd.read_csv("../ligand_information.csv")
	    ligands = {}
	    for i in range(n):
	        onbits = table.iloc[i]['OnBits'].split(",")
	        score = table.iloc[i]['Score']
	        smiles = table.iloc[i]['SMILES']
	        lig = Ligand(i,score,smiles,onbits)
	        ligands[i] = lig

	    return ligands

	def calculate_distance(self, lig1, lig2):
		"""
	    Distance function
	    Parameters
	    ---------
	    lig1
	        First ligand
	    lig2
	    	Second ligand

	    Returns: distance between two ligands
	    """
	    pass

class HierarchicalClustering(Clustering):
	"""Implementation of HC using """
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
		# Initialize cluster dictionary
		clusters = []
		for key in ligands.keys():
    		clusters.append([key]) # put all ligands in their own cluster
			
		# Initialize distance matrix
		dist = np.zeroes(shape=(len(ligands),len(ligands)))

		# loop through bottom half of matrix and fill in distances
		for i in range(len(ligands)):
			for j in in range(len(ligands)):
				if i == j:
					break
				dist[i,j] = calculate_distance(i,j)


		while len(set(clusters.values())) > k
			# find min 
			min_i, min_j = np.unravel_index(np.argmin(dist), dist.shape)

			# find the two ligands in clusters
			found_i = -1
			found_j = -1
			for idx in range(len(clusters)):
				if found_i >= 0 and if found_j >= 0:
					break
    			if found_i < 0 and min_i in clusters[idx]:
    				found_i = idx
    			if found_j < 0 and min_j in clusters[idx]:
    				found_j = idx

    		# merge the clusters
    		clusters[found_i]+= clusters[found_j]
    		del clusters[found_j]

      
			
			# update distance matrix using single linkage
			

		


class PartitionClustering(Clustering):
	"""An example docstring for a class definition."""
	def __init__(self):
		"""
	    Blah blah blah.
	    Parameters
	    ---------
	    name
	        A string to assign to the `name` instance attribute.
	    """
	    super().__init__()
		self.name = name

	def cluster(self, ligands):
		"""
		Return information about an instance created from ExampleClass.
		"""
		pass