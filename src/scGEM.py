import scot as sc


# load data
X = np.genfromtxt("../data/scGEM/Cheow_expression.csv", delimiter=',') #Row-formatted data; shape (#DataPoints x #Features)
y = np.genfromtxt("../data/scGEM/Cheow_methylation.csv", delimiter=',') #Row-formatted data; shape (#DataPoints x #Features)

# set parameters
k = 20
eps =  .005

# perform SCOT alignment
X_transported, y_normalized = sc.scot(X, y, k, eps)

# evaluate results
fracs, xs = em.get_fracs(X_transported, y_normalized)
print("Average fraction of samples closer than true match is ", np.mean(fracs))




