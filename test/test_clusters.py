
import pytest
import pandas as pd
from clusters.algs import *

@pytest.fixture

def test_get_ligands():
	pc = PartitionClustering()
	ligands = pc.get_ligands(3)
	assert ligands.length == 3
	assert ligands == {1:[360,489,915], 2:[53,623,650], 3:[332,342,650]}

def test_calculate_distance():
	pc = PartitionClustering()
	assert pc.calculate_distance([10], [10]) == 0
	assert pc.calculate_distance([4], [10]) == 1
	assert pc.calculate_distance([4, 10], [10]) == 0.5

def test_distance_matrix():
	pc = PartitionClustering()
	ligands = pc.get_ligands(4)
	d = pc.distance_matrix(ligands)
	assert d.shape == (4, 4)
	assert d[0,0] == float("inf")
	assert d[1,2] == 1/5

def test_silhouette_score():
	pc = PartitionClustering()
	ligands = pc.get_ligands(4)
	assert pc.silhouette_score(ligands, [0,0,1,1]) == -0.16249999999999998

def test_rand_index():
	pc = PartitionClustering()
	# see if my rand index matches the example in lecture slides
	clustering1 = [0, 0, 1, 0, 1, 2]
	clustering2 = [0, 0, 1, 1, 0, 2]
	assert pc.rand_index(clustering1, clustering2) == 0.6

def test_partitioning():
	pc = PartitionClustering()
	ligands = pc.get_ligands(4)
	clusters = pc.cluster(ligands, k = 2)
	assert set(clusters) == {0, 1}
	assert len(clusters) == 4

def test_hierarchical():
	hc = HierarchicalClustering()
	ligands = hc.get_ligands(4)
	clusters = hc.cluster(ligands, k = 1)
	assert set(clusters) == {0}
	assert len(clusters) == 4
