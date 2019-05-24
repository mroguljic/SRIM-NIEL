import ROOT as r
import matplotlib.pyplot as plt
import numpy as np
import scipy.optimize

def targetFunction(x,Rmin,Rmax):
    #everything in microns
    if(x<Rmin):
        return(8.0)
    elif(x>Rmax):
        return(0.0)
    else:
        return(20.0)

def calculateSpectrum(fileNames,Rmin,Rmax,plot=True):

    vacancies       = []
    depths          = []

    for fileName in fileNames:

        f = open(fileName,'r')
        print('Opened '+fileName)
        tempVacancies   = []
        tempDepths      = []

        for line in f.readlines():
            splitLine = line.split()
            if(len(splitLine)==3):
                #print(line)
                try:
                    tempDepth   = float(splitLine[0])/10000. #in microns
                    tempVacancy = (float(splitLine[1])+float(splitLine[2]))*100000 #arbitrary units
                    #print(tempDepth, tempVacancy)
                    tempVacancies.append(tempVacancy)
                    tempDepths.append(tempDepth)
                except:
                    continue

        f.close()
        if(len(depths)==0):
            depths.append(tempDepths)
        else:
            if(depths[0]!=tempDepths):
                print(fileName+' ---Depths not matching!---')

        vacancies.append(tempVacancies)

    A = np.array(vacancies)
    A = np.transpose(A)
    #print(np.shape(A))
    B = [] #sampling fitting function

    for x in depths[0]:
        B.append(targetFunction(x,Rmin,Rmax))

    B = np.array(B)
    #print(np.shape(B))

    fitResults = scipy.optimize.nnls(A, B)
    spectrum = np.array(fitResults[0])
    norm= np.linalg.norm(spectrum, ord=1)
    if norm==0.:
        print("Norm is zero!")

    print('Calculated spectrum\n')
    for i,fileName in enumerate(fileNames):
        energyString = fileName.replace("VACANCY_","").replace(".txt","")
        weight       = spectrum[i]/norm
        print(energyString+" --- ",weight)

    if(plot):
        y_fit = np.matmul(A,spectrum)

        plt.plot(depths[0], y_fit, 'r--', depths[0], B, 'bs')
        plt.show()




fileNames = ['VACANCY_3.3MeV.txt','VACANCY_3.7MeV.txt','VACANCY_4.1MeV.txt','VACANCY_4.5MeV.txt','VACANCY_4.9MeV.txt','VACANCY_5.2MeV.txt','VACANCY_5.6MeV.txt',
'VACANCY_3.4MeV.txt','VACANCY_3.8MeV.txt','VACANCY_4.2MeV.txt','VACANCY_4.6MeV.txt','VACANCY_5.3MeV.txt','VACANCY_5.7MeV.txt','VACANCY_6.0MeV.txt',
'VACANCY_3.5MeV.txt','VACANCY_3.9MeV.txt','VACANCY_4.3MeV.txt','VACANCY_4.7MeV.txt','VACANCY_5.0MeV.txt','VACANCY_5.4MeV.txt','VACANCY_5.8MeV.txt','VACANCY_6.1MeV.txt',
'VACANCY_3.6MeV.txt','VACANCY_4.0MeV.txt','VACANCY_4.4MeV.txt','VACANCY_4.8MeV.txt','VACANCY_5.1MeV.txt','VACANCY_5.5MeV.txt','VACANCY_5.9MeV.txt']
Rmin            = 100. #microns
Rmax            = 300.
calculateSpectrum(fileNames,Rmin,Rmax,plot=True)