def louvain(i, j, val, dim, partition_method, initial_membership, weights, resolution, node_sizes, seed, verbose):
  import louvain
  import igraph as ig
  import numpy
  from scipy.sparse import csc_matrix
  data = csc_matrix((val, (i, j)), shape = dim)
  # vcount = max(data.shape)
  sources, targets = data.nonzero()
  edgelist = zip(sources.tolist(), targets.tolist())
  G = ig.Graph(edges = list(edgelist))

  # G = ig.Graph.Adjacency(data.tolist())
  
  if partition_method == 'ModularityVertexPartition':
    partition = louvain.CPMVertexPartition(G, initial_membership = initial_membership, weights = weights)
  elif partition_method == 'RBConfigurationVertexPartition':
    partition = louvain.CPMVertexPartition(G, initial_membership = initial_membership, weights = weights, resolution_parameter = resolution)
  elif partition_method == 'RBERVertexPartition':
    partition = louvain.CPMVertexPartition(G, initial_membership = initial_membership, weights = weights, node_sizes = node_sizes, resolution_parameter = resolution)
  elif partition_method == 'CPMVertexPartition':
    partition = louvain.CPMVertexPartition(G, initial_membership = initial_membership, weights = weights, node_sizes = node_sizes, resolution_parameter = resolution)
  elif partition_method == 'SignificanceVertexPartition':
    partition = louvain.CPMVertexPartition(G, initial_membership = initial_membership, node_sizes = node_sizes)
  elif partition_method == 'SurpriseVertexPartition':
    partition = louvain.CPMVertexPartition(G, initial_membership = initial_membership, weights = weights, node_sizes = node_sizes)
  else: 
    raise ValueError('partition_method ' + partition_method + ' is NOT supported.')
  
  if seed != None: 
    louvain.set_rng_seed(seed)

  optimiser = louvain.Optimiser()
  diff = optimiser.optimise_partition(partition)
  
  # ig.plot(partition)
  return partition
