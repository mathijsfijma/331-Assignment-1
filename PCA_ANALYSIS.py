from statistics import variance
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

from isOrthogonal import isOrthogonal

with open("PCA_EIGENVALUES.csv") as PCA_EIGENVALUES:
    eigenvalues = np.loadtxt(PCA_EIGENVALUES, delimiter=",", skiprows=1)

eigenvectors = np.loadtxt(open("PCA_EIGENVECTOR.csv", "rb"), delimiter=",", skiprows=1)
g = [[sum(a*b for a,b in zip(X_row,Y_col)) for Y_col in zip(*eigenvectors.transpose())] for X_row in eigenvectors]
#basically orthogonal
print(np.size(eigenvectors))
print(isOrthogonal(eigenvectors, 56, 56))
#assert eigenvectors*eigenvectors.transpose() == np.eye(len(eigenvectors)), f"not orthogonal"
print(eigenvalues)
varianceRatios = []
for i in range(len(eigenvalues)):
    varianceRatios.append(sum(eigenvalues[0:i])/sum(eigenvalues))

print(varianceRatios)
plt.plot(np.arange(len(eigenvalues)), varianceRatios)
plt.show()
