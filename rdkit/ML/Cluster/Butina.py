# Copyright (C) 2007-2024 Greg Landrum and other RDKit contributors
#   All Rights Reserved
#
""" Implementation of the clustering algorithm published in:
  Butina JCICS 39 747-750 (1999)

"""
import numpy as np

def EuclideanDist(pi, pj):
  """Calculate the Euclidean distance between two points."""
  pi, pj = np.asarray(pi), np.asarray(pj)
  return np.sqrt(np.sum((pi - pj) ** 2))

def compute_distance_matrix(data, n_pts, dist_func):
  """Compute the distance matrix for the given data."""
  dist_matrix = np.zeros((n_pts, n_pts))
  for i in range(n_pts):
    for j in range(i):
      dist_matrix[i, j] = dist_matrix[j, i] = dist_func(data[i], data[j])
  return dist_matrix

def ClusterData(data, nPts, distThresh, isDistData=False, distFunc=EuclideanDist, reordering=False):
  """  clusters the data points passed in and returns the list of clusters

    **Arguments**

      - data: a list, tuple, or numpy array of items with the input data
        (see discussion of _isDistData_ argument for the exception)

      - nPts: the number of points to be used

      - distThresh: elements within this range of each other are considered
        to be neighbors

      - isDistData: set this toggle when the data passed in is a
          distance matrix.  The distance matrix should be stored 
          in one of two formats: as an nxn NumPy array, or as a 
          symmetrically stored list or 1D array generated using a 
          similar process to the example below:

            dists = []
            for i in range(nPts):
              for j in range(i):
                dists.append( distfunc(i,j) )

      - distFunc: a function to calculate distances between points.
           Receives 2 points as arguments, should return a float

      - reordering: if this toggle is set, the number of neighbors is updated
           for the unassigned molecules after a new cluster is created such
           that always the molecule with the largest number of unassigned
           neighbors is selected as the next cluster center.

    **Returns**

      - a tuple of tuples containing information about the clusters:
         ( (cluster1_elem1, cluster1_elem2, ...),
           (cluster2_elem1, cluster2_elem2, ...),
           ...
         )
         The first element for each cluster is its centroid.

  """
  if isDistData:
    # Check if data is a supported type
    if not isinstance(data, (list, tuple, np.ndarray)):
      raise TypeError(f"Unsupported type for data, {type(data)}")
    
    # Check if data is a 1D array or list
    if isinstance(data, (list, tuple)) or (isinstance(data, np.ndarray) and data.ndim == 1):
      # Check if data length matches the required number of points
      if len(data) != (nPts * (nPts - 1)) // 2:
        raise ValueError("Mismatched input data dimension and nPts")
        
      # Create a distance matrix from the 1D data
      dist_matrix = np.zeros((nPts, nPts))
      idx = np.tril_indices(nPts, -1)
      dist_matrix[idx] = data
      dist_matrix += dist_matrix.T
    else:
      # Check if data is a matrix of the correct shape and use it as distance matrix
      if data.shape != (nPts, nPts):
          raise ValueError(f"Input data with shape {data.shape} is not a matrix of the required shape {(nPts, nPts)}")
      dist_matrix = data
  else:
    # Compute distance matrix from the data points
    dist_matrix = compute_distance_matrix(data, nPts, distFunc)

  # Initialize neighbor lists
  neighbor_lists = [np.where(dist_matrix[i] <= distThresh)[0].tolist() for i in range(nPts)]

  # Sort points by the number of neighbors in descending order
  sorted_indices = [(len(neighbors), idx) for idx, neighbors in enumerate(neighbor_lists)]
  sorted_indices.sort(reverse=True)

  # Initialize clusters and a seen array to keep track of processed points
  clusters = []
  seen = np.zeros(nPts, dtype=bool)

  # Process all candidate clusters that have at least two members
  while sorted_indices and sorted_indices[0][0] > 1:
    _, idx = sorted_indices.pop(0)
    if seen[idx]:
      continue

    # Create a new cluster and mark points as seen
    cluster = [idx]
    seen[idx] = True
    for neighbor in neighbor_lists[idx]:
      if not seen[neighbor]:
        cluster.append(neighbor)
        seen[neighbor] = True
    
    clusters.append(tuple(cluster))

    # Update the number of neighbors:
    # remove all members of the new cluster from the list of
    # neighbors and reorder the sorted_indices
    if reordering:
      # Get the set of unassigned and affected molecules, i.e. all unseen molecules
      # which have at least one of the members of the new cluster
      # as a neighbor
      affected = set(neighbor for point in cluster for neighbor in neighbor_lists[point] if not seen[neighbor])

      # Loop over all remaining molecules in sorted_indices but only
      # consider unassigned and affected compounds
      for ii, element in enumerate(sorted_indices):
        affected_point = element[1]
        if affected_point in affected:
          # Update the number of neighbors
          new_neighbors = [nbr for nbr in neighbor_lists[affected_point] if not seen[nbr]]
          neighbor_lists[affected_point] = new_neighbors
          sorted_indices[ii] = (len(new_neighbors), affected_point)
      # Reorder the list
      sorted_indices.sort(reverse=True)
  
  # Process any remaining single-point clusters
  while sorted_indices:
    _, idx = sorted_indices.pop(0)
    if seen[idx]:
      continue
    clusters.append(tuple([idx]))
  return tuple(clusters)
