import scipy as sp
from scipy.stats import qmc
from numpy.random import Generator, PCG64

no_samples = 500
no_edges = 17
seed = 45643564852856852352111
sampler = qmc.LatinHypercube(d = no_edges,rng=Generator(PCG64(seed))) # 17 edges
sample = sampler.random(n=no_samples) # still in a [0,1]^17 hypercube 500 samples!
TEs = ['AMAZ', 'AMOC', 'ASI', 'ENSO', 'GIS', 'PERM', 'REEF', 'WAIS', 'WAM'] # see README.txt
TE_to_ind = {te:i for i,te in enumerate(TEs)}

folder = '../Data/'
df_interactions = pd.read_excel(folder+'Interactions_TableA2.xlsx',usecols='A:F',header=[0])
# df_interactions.head()

a = df_interactions['min'].to_numpy()
b = df_interactions['max'].to_numpy()
A = []
for r in range(no_samples):
    weights = a+sample[r,:]*(b-a)
    mat = np.zeros([len(TEs)]*2)
    for index, row in df_interactions.iterrows():
        i = TE_to_ind[row['Source']]
        j = TE_to_ind[row['Target']]
        mat[i,j] = weights[index]
    A.append(mat)
    # np.savetxt("../Data/Interaction matrices/A_%s.csv"%r,A[r])
