def umap(i, j, val, dim, n, n_c, metric, n_epochs, negative_sample_rate, learning_rate, init, mdist, spread, set_op_mix_ratio, local_connectivity, repulsion_strength, a, b, random_state, metric_kwds, angular_rp_forest, target_n_neighbors, target_metric, target_metric_kwds, target_weight, transform_seed, verbose):
  import umap
  import numpy
  from scipy.sparse import csc_matrix
  data = csc_matrix((val, (i, j)), shape = dim)
  res = umap.UMAP(n_neighbors = n, 
				  	n_components = n_c, 
				  	metric = metric, 
				  	n_epochs = n_epochs, 
				  	negative_sample_rate = negative_sample_rate,
				  	learning_rate = learning_rate,
				  	init = init,
				  	min_dist = mdist,
				  	spread = spread,
				  	set_op_mix_ratio = set_op_mix_ratio,
				  	local_connectivity = local_connectivity,
				  	repulsion_strength = repulsion_strength,
				  	a = a, 
				  	b = b, 
				  	random_state = random_state,
				  	metric_kwds = metric_kwds, 
				  	angular_rp_forest = angular_rp_forest,
				  	target_n_neighbors = target_n_neighbors,
				  	target_metric = target_metric, 
				  	target_metric_kwds = target_metric_kwds, 
				  	target_weight = target_weight, 
				  	transform_seed = transform_seed, 
				  	verbose = verbose).fit(data.toarray()) 
  return res
