from audioop import avg
from cProfile import label
from cmath import log
from statistics import variance
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd


#Opening Eigenpair data
with open("PCA_EIGENVALUES.csv") as PCA_EIGENVALUES:
    eigenvalues = np.loadtxt(PCA_EIGENVALUES, delimiter=",", skiprows=1)

eigenvectors = np.loadtxt(open("PCA_EIGENVECTOR.csv", "rb"), delimiter=",", skiprows=1)

print(np.size(eigenvectors))

print(eigenvalues)

#Calculating explained variance and plotting to find number of principle components
varianceRatios = []
for i in range(len(eigenvalues)):
    #varianceRatios.append((eigenvalues[i])/sum(eigenvalues))
    varianceRatios.append(sum(eigenvalues[0:i])/sum(eigenvalues))

print(varianceRatios)
plt.plot(varianceRatios)
#plt.yscale('log')
plt.plot([0, 6], [0.9923, 0.9923], color='r', linestyle='-', linewidth=1)
plt.plot([6, 6], [0, 0.9923], color='r', linestyle='-', linewidth=1)
plt.xlabel("No. of Principal Components")
plt.ylabel("Explained Variance (Cumulative)")
plt.show()
#Ajusted explained variance calculated (Var(1)>>Var(2:N))
# varianceRatios2 = []
# for i in range(len(eigenvalues)):
#     varianceRatios2.append(sum(eigenvalues[1:i])/sum(eigenvalues))

# print(varianceRatios2)

spectra = np.loadtxt(open("DS19hH2_dk0_FTIR_Spectra_instant_coffee.csv", "rb"), delimiter=",")
mu = spectra.mean(axis=0)
plt.plot(mu)
plt.xlabel("Wavelength (nm)")
plt.ylabel("Reflectance")
plt.show()


#Projecting spectra onto most important principal components
principalComponents = np.arange(0,5)

for i in range(len(principalComponents)):
    principalComponents[i] = (np.dot(spectra[i,:],eigenvectors[i,:]))
print(principalComponents)

principalComponents = np.arange(0,6)
for i in range(len(principalComponents)):

    for j in range(286):

            principalComponents[i] += (np.dot(spectra[i,:],eigenvectors[j,:]))

truePrincipalComponents = np.arange(0,56)
for i in range(len(truePrincipalComponents)):

    for j in range(286):

            truePrincipalComponents[i] += (np.dot(spectra[i,:],eigenvectors[j,:]))
    #truePrincipalComponents[i] = truePrincipalComponents[i]/56

proper_projections = np.empty((56,286))
for i in range(56):
    row_mean = np.mean(spectra[i,:])
    for j in range(286):
        proper_projections[i,j] = row_mean
        for k in range(len(principalComponents)):
            proper_projections[i,j] += ((principalComponents[k])*eigenvectors[k,j])

true_proper_projections = np.empty((56,286))
for i in range(56):
    row_mean = np.mean(spectra[i,:])
    for j in range(286):
        true_proper_projections[i,j] = row_mean
        for k in range(len(truePrincipalComponents)):
            true_proper_projections[i,j] += ((truePrincipalComponents[k])*eigenvectors[k,j])


mu = proper_projections.mean(axis=0)
true_mu = true_proper_projections.mean(axis=0)
plt.plot(true_mu,label = 'original')
plt.plot(mu,label='reconstruction')
plt.legend()
plt.xlabel("Wavelength (nm)")
plt.ylabel("Reflectance")
plt.show()
#print(proper_projections)

#Explained variance plot

# plt.semilogy(varianceRatios, 'o')
# plt.axvline(max(varianceRatios)-min(varianceRatios))
# plt.show()

final_projections = np.empty((56,6))
for i in range(6):
    for j in range(29):
        final_projections[j,i] = np.dot(spectra[j,:], eigenvectors[i,:])
    for k in range(27):
        final_projections[29+k,i] = np.dot(spectra[29+k,:], eigenvectors[i,:])
chur =2
target1 = np.arange(0,29)
target2 = np.arange(29,56) #0,2;#1,2;#3,2
# plt.scatter(final_projections[target1,5], final_projections[target1,2])
# plt.scatter(final_projections[target2,5], final_projections[target2,2])
# plt.show()
# plt.scatter(final_projections[target1,1], final_projections[target1,2])
# plt.scatter(final_projections[target2,1], final_projections[target2,2])
# plt.show()
# plt.scatter(final_projections[target1,0], final_projections[target1,2])
# plt.scatter(final_projections[target2,0], final_projections[target2,2])
# plt.show()
# plt.scatter(final_projections[target1,3], final_projections[target1,2])
# plt.scatter(final_projections[target2,3], final_projections[target2,2])
# plt.show()
# plt.scatter(final_projections[:,4], final_projections[:,2])
# plt.scatter(final_projections[:,2], final_projections[:,4])
# plt.show()

# arabica = np.empty((10,286))
# robusta = np.empty((10,286))
# for i in range(10):
#     row_mean = np.mean(spectra[i,:])
#     for j in range(286):
#         arabica[i,j] = row_mean
#         for k in range(len(principalComponents)):
#             arabica[i,j] += (principalComponents[k])*eigenvectors[k,j]

# for i in range(10):
#     row_mean = np.mean(spectra[i+29,:])
#     for j in range(286):
#         robusta[i,j] = row_mean
#         for k in range(len(principalComponents)):
#             robusta[i,j] += (principalComponents[k])*eigenvectors[k,j]

# arabica_mu = arabica.mean(axis=0)
# robusta_mu = robusta.mean(axis=0)
# plt.plot(arabica_mu,label= 'arabica')
# plt.plot(robusta_mu,label = 'robusta')
# plt.legend()
# plt.show()


arabica = np.empty((1,286))
robusta = np.empty((1,286))

i = 0
row_mean = np.mean(spectra[i,:])
for j in range(286):
    arabica[0,j] = row_mean
    for k in range(len(principalComponents)):
        arabica[0,j] += ((principalComponents[k])*eigenvectors[k,j])

row_mean = np.mean(spectra[i+29,:])
for j in range(286):
    robusta[0,j] = row_mean
    for k in range(len(principalComponents)):
        robusta[0,j] += ((principalComponents[k])*eigenvectors[k,j])

arabica_mu = arabica.mean(axis=0)
robusta_mu = robusta.mean(axis=0)
plt.plot(arabica_mu,label= 'arabica')
plt.plot(robusta_mu,label = 'robusta')
plt.xlabel("Wavelength (nm)")
plt.ylabel("Reflectance")
plt.legend()
plt.show()
