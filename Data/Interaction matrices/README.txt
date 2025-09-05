Index ordering of the tipping elements:
['AMAZ', 'AMOC', 'ASI', 'ENSO', 'GIS', 'PERM', 'REEF', 'WAIS', 'WAM']

Adjacency matrices here is ordered as A[source,target]

## Sampling:
seed = 45643564852856852352111
sampler = qmc.LatinHypercube(d=len(df_W24),rng=Generator(PCG64(seed)))
sample = sampler.random(n=500) # still in a [0,1]^18 hypercube 500 samples!
