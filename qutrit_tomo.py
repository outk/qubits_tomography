'''
This is a implementation of Quantum State Tomography for Qutrits,
using techniques of following papars.

'Iterative algorithm for reconstruction of entangled states(10.1103/PhysRevA.63.040303)'
'Diluted maximum-likelihood algorithm for quantum tomography(10.1103/PhysRevA.75.042108)'
'Qudit Quantum State Tomography(10.1103/PhysRevA.66.012303)'
'On-chip generation of high-dimensional entangled quantum states and their coherent control(Nature volume 546, pages622-626(2017))'

'''

import numpy as np
from numpy import array, kron, trace, identity, sqrt, zeros, exp, pi, conjugate, random
from scipy.linalg import sqrtm
from datetime import datetime
from concurrent import futures
import os
from pathlib import Path
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import pickle


""" 
Definition of Three Frequency Bases: 

    fb1 = array([1, 0, 0])
    fb2 = array([0, 1, 0])
    fb3 = array([0, 0, 1])

"""

zero_base_array1 = zeros((1,3))
zero_base_array1[0][0] = 1
fb1 = zero_base_array1

zero_base_array2 = zeros((1,3))
zero_base_array2[0][1] = 1
fb2 = zero_base_array2

zero_base_array3 = zeros((1,3))
zero_base_array3[0][2] = 1
fb3 = zero_base_array3



""" Make Measurement Bases """

mb1 = (conjugate((fb1 + fb2).T) @ (fb1 + fb2)) / 2
mb2 = (conjugate((fb1 + fb3).T) @ (fb1 + fb3)) / 2
mb3 = (conjugate((fb2 + fb3).T) @ (fb2 + fb3)) / 2
mb4 = (conjugate((exp( 2*pi*1j/3) * fb1 + (exp(-2*pi*1j/3)) * fb2).T) @ (exp( 2*pi*1j/3) * fb1 + (exp(-2*pi*1j/3) * fb2))) / 2
mb5 = (conjugate((exp(-2*pi*1j/3) * fb1 + (exp( 2*pi*1j/3)) * fb2).T) @ (exp(-2*pi*1j/3) * fb1 + (exp( 2*pi*1j/3) * fb2))) / 2
mb6 = (conjugate((exp( 2*pi*1j/3) * fb1 + (exp(-2*pi*1j/3)) * fb3).T) @ (exp( 2*pi*1j/3) * fb1 + (exp(-2*pi*1j/3) * fb3))) / 2
mb7 = (conjugate((exp(-2*pi*1j/3) * fb1 + (exp( 2*pi*1j/3)) * fb3).T) @ (exp(-2*pi*1j/3) * fb1 + (exp( 2*pi*1j/3) * fb3))) / 2
mb8 = (conjugate((exp( 2*pi*1j/3) * fb2 + (exp(-2*pi*1j/3)) * fb3).T) @ (exp( 2*pi*1j/3) * fb2 + (exp(-2*pi*1j/3) * fb3))) / 2
mb9 = (conjugate((exp(-2*pi*1j/3) * fb2 + (exp( 2*pi*1j/3)) * fb3).T) @ (exp(-2*pi*1j/3) * fb2 + (exp( 2*pi*1j/3) * fb3))) / 2

bases = array([mb1, mb2, mb3, mb4, mb5, mb6, mb7, mb8, mb9])

def makeBases(numberOfQutrits, bases):
    for _ in range(numberOfQutrits-1):
        fixedBases = []
        for base1 in bases:
            fixedBases.extend([kron(base1, base2) for base2 in bases])
        bases = fixedBases.copy()

    return array(bases)



""" Get Experimental Data """

def getExperimentalData(pathOfExperimentalData):
    """
    getExperimentalData(pathOfExperimentalData(string)):

        This function is getting experimental data from file at "pathOfExperimentalData",
        
        ----------------------------------------------------------------------------------------
        return:

            np.array of them.

    """

    with open(pathOfExperimentalData) as f:
        experimentalData = []
        for s in f.readlines():
            experimentalData.extend(map(int, s.strip().split()))

    return array(experimentalData)



""" Iterative Algorithm for Quantum Tomography """

def doIterativeAlgorithm(numberOfQutrits, bases, listOfExperimentalData):
    """
    doIterativeAlgorithm():

        This function is to do iterative algorithm(10.1103/PhysRevA.63.040303) and diluted MLE algorithm(10.1103/PhysRevA.75.042108) to a set of data given from a experiment.
        This recieve four variables (numberOfQutrits, bases, maxNumberOfIteration, listAsExperimentalData),
        and return most likely estimated density matrix (np.array) and total time of calculation(datetime.timedelta).

        First quantum state is an Identity: Maximum Mixed State.

        --------------------------------------------------------------------------------------------------------------
        Return:

            most likely estimated density matrix(np.array)

    """

    """ Setting initial parameters """
    iter = 0
    epsilon = 1000
    endDiff = 10e-10
    diff = 100
    # TolFun = 10e-11
    # traceDistance = 100
    maxNumberOfIteration = 100000

    dataList = listOfExperimentalData
    totalCountOfData = sum(dataList)
    nDataList = dataList / totalCountOfData # nDataList is a list of normarized data
    densityMatrix = identity(3 ** numberOfQutrits) # Initial State: Maximun Mixed State
    
    """ Start iteration """
    # while traceDistance > TolFun and iter <= maxNumberOfIteration:
    while diff > endDiff and iter <= maxNumberOfIteration:

        probList = [trace(bases[i] @ densityMatrix) for i in range(len(bases))]
        nProbList = probList / sum(probList)
        rotationMatrix = sum([(nDataList[i] / probList[i]) * bases[i] for i in range(9 ** numberOfQutrits)])

        """ Normalization of Matrices for Measurement Bases """
        U = np.linalg.inv(sum(bases)) / sum(probList)   
        rotationMatrixLeft = (identity(3 ** numberOfQutrits) + epsilon * U @ rotationMatrix) / (1 + epsilon)
        rotationMatrixRight = (identity(3 ** numberOfQutrits) + epsilon * rotationMatrix @ U) / (1 + epsilon)

        """ Calculation of updated density matrix """
        modifiedDensityMatrix = rotationMatrixLeft @ densityMatrix @ rotationMatrixRight / trace(rotationMatrixLeft @ densityMatrix @ rotationMatrixRight)
        eigValueArray, eigVectors = np.linalg.eig(densityMatrix - modifiedDensityMatrix)
        traceDistance = sum(np.absolute(eigValueArray)) / 2

        """ Update Likelihood Function, and Compared with older one """
        LikelihoodFunction = sum([nDataList[i] * np.log(nProbList[i]) for i in range(9 ** numberOfQutrits)])
        probList = [trace(bases[i] @ modifiedDensityMatrix) for i in range(len(bases))]
        nProbList = probList / sum(probList)
        modifiedLikelihoodFunction = sum([nDataList[i] * np.log(nProbList[i]) for i in range(9 ** numberOfQutrits)]) 
        diff = modifiedLikelihoodFunction - LikelihoodFunction

        """ Show Progress of Calculation """
        progress = 100 * iter / maxNumberOfIteration
        if progress % 5 == 0:
            msg = "Progress of calculation: " + str(int(progress)) + "%"
            print(msg)

        """ Increment """
        iter += 1

        """ Check Increasing of Likelihood Function  """
        if diff < 0:
            epsilon = epsilon * 0.1
            continue

        """ Update Density Matrix """
        densityMatrix = modifiedDensityMatrix.copy()

    
    """ Check That Max Iteration Number was appropriate """
    if iter >= maxNumberOfIteration:
        print("----------------------------------------------")
        print("Iteration time reached max iteration number.")
        print("The number of max iteration times is too small.")
        print("----------------------------------------------")
    
    """ Show the total number of iteration """
    endIterationTimes = iter
    emsg = "Iteration was '" + str(endIterationTimes) + "' times."
    print(emsg)

    return modifiedDensityMatrix



""" Calculate Fidelity """

def calculateFidelity(idealDensityMatrix, estimatedDensityMatrix):
    """
    calculateFidelity(idealDensityMatrix, estimatedDensityMatrix):

    """

    fidelity = np.real(trace(sqrtm(sqrtm(idealDensityMatrix) @ estimatedDensityMatrix @ sqrtm(idealDensityMatrix))))

    return fidelity



""" Iterative Simulation """

def doIterativeSimulation(numberOfQutrits, bases, pathOfExperimentalData, idealDensityMatrix, resultDirectoryName):
    """
    doIterativeSimulation(numberOfQutrits, bases, pathOfExperimentalData, idealDensityMatrix, resultDirectoryName)


    """

    """ Get Experimental Data"""
    listOfExperimentalData = getExperimentalData(pathOfExperimentalData)

    """ Calculate """
    estimatedDensityMatrix = doIterativeAlgorithm(numberOfQutrits, bases, listOfExperimentalData)
    fidelity = calculateFidelity(idealDensityMatrix, estimatedDensityMatrix)

    """ Make File Name of result """
    l = 0
    r = len(pathOfExperimentalData)-1
    for i in range(len(pathOfExperimentalData)):
        if pathOfExperimentalData[len(pathOfExperimentalData)-1-i] == ".":
            r = len(pathOfExperimentalData)-1-i
        if pathOfExperimentalData[len(pathOfExperimentalData)-1-i] == "/" or pathOfExperimentalData[len(pathOfExperimentalData)-1-i] == "\\":
            l = len(pathOfExperimentalData)-i
            break
    resultFileName = pathOfExperimentalData[l:r]
    resultFilePath = '.\\result\\qutrit\\iterative\\' + resultDirectoryName + '\\' + resultFileName + '_result' + '.txt'

    """ Save Result """
    with open(resultFilePath, mode='a') as f:
        f.writelines(str(fidelity) + '\n')

    """ Make 3D Plot """
    plotResult(numberOfQutrits, estimatedDensityMatrix)



""" Poisson Distributed Simulation """

def doPoissonDistributedSimulation(numberOfQutrits, bases, pathOfExperimentalData, idealDensityMatrix, resultDirectoryName):
    """
    doPoissonDistributedSimulation(numberOfQutrits, bases, pathOfExperimentalData, idealDensityMatrix, resultDirectoryName)


    """

    """ Get Experimental Data"""
    listOfExperimentalData = getExperimentalData(pathOfExperimentalData)

    """ Calculate """
    estimatedDensityMatrix = doIterativeAlgorithm(numberOfQutrits, bases, random.poisson(listOfExperimentalData))
    fidelity = calculateFidelity(idealDensityMatrix, estimatedDensityMatrix)

    """ Make File Name of result """
    l = 0
    r = len(pathOfExperimentalData)-1
    for i in range(len(pathOfExperimentalData)):
        if pathOfExperimentalData[len(pathOfExperimentalData)-1-i] == ".":
            r = len(pathOfExperimentalData)-1-i
        if pathOfExperimentalData[len(pathOfExperimentalData)-1-i] == "/" or pathOfExperimentalData[len(pathOfExperimentalData)-1-i] == "\\":
            l = len(pathOfExperimentalData)-i
            break
    resultFileName = pathOfExperimentalData[l:r]
    resultFilePath = '.\\result\\qutrit\\iterative\\' + resultDirectoryName + "\\" + resultFileName + '_result' + '.txt'

    """ Save Result """
    with open(resultFilePath, mode='a') as f:
        f.write(str(fidelity) + '\n')


def plotResult(numberOfQutrits, densityMatrix):
    """
    plotResult(numberOfQutrits, densityMatrix)
    
    """

    baseNames = np.array([])

    """ Plot Setting """
    xedges, yedges = np.arange(3**numberOfQutrits), np.arange(3**numberOfQutrits)
 
    xpos, ypos = np.meshgrid(xedges, yedges) # x,y座標を3D用の形式に変換（その１）
    zpos = 0 # zは常に0を始点にする
    
    dx = 1 # x座標の幅を設定
    dy = 1 # y座標の幅を設定
    dz = densityMatrix.ravel() # z座標の幅は棒の長さに相当
    
    xpos = xpos.ravel() # x座標を3D用の形式に変換（その２）
    ypos = ypos.ravel() # y座標を3D用の形式に変換（その２)

    fig = plt.figure() # 描画領域を作成
    ax1 = fig.add_subplot(121, projection="3d") # 3Dの軸を作成    
    ax1.bar3d(xpos,ypos,zpos,dx,dy,np.real(dz), edgecolor='black') # ヒストグラムを3D空間に表示
    plt.title("Real Part") # タイトル表示
    # plt.xlabel("X") # x軸の内容表示
    plt.xticks(np.arange(0, 3**numberOfQutrits, 1), labels=baseNames)
    # plt.ylabel("Y") # y軸の内容表示
    plt.yticks(np.arange(0, 3**numberOfQutrits, 1), labels=baseNames)
    # ax1.set_zlabel("Z") # z軸の内容表示
    ax1.set_zlim(-0.1, 0.5)
    
    ax2 = fig.add_subplot(122, projection="3d") # 3Dの軸を作成
    ax2.bar3d(xpos,ypos,zpos,dx,dy,np.imag(dz), edgecolor='black') # ヒストグラムを3D空間に表示
    plt.title("Imaginary Part") # タイトル表示
    # plt.xlabel("X") # x軸の内容表示
    plt.xticks(np.arange(0, 3**numberOfQutrits, 1), labels=baseNames)
    # plt.ylabel("Y") # y軸の内容表示
    plt.yticks(np.arange(0, 3**numberOfQutrits, 1), labels=baseNames)
    # ax2.set_zlabel("Z") # z軸の内容表示
    ax2.set_zlim(-0.1, 0.5)

    plt.show()

    with open('qutritplottest'+'_plot.pkl', mode='wb') as f:
        pickle.dump(fig, f)



""" Get Number of Qutrits """

def getNumberOfQutrits():
    """
    getNumberOfQutrits()

    """

    print("------------------------------------------------------------")
    print("PLEASE ENTER NUMBER OF QUTRITS")
    print("------------------------------------------------------------")
    print(">>")

    numberOfQutrits = int(input())

    return numberOfQutrits



""" Get Path of Experimental Data Directory """

def getExperimentalDataDirectoryPath():
    """
    getExperimentalDataDirectoryPath()

    """

    print("------------------------------------------------------------")
    print("PLEASE ENTER PATH OF EXPERIMENTAL DATA DIRECTORY")
    print("")
    print("LIKE THIS >> .\\datadirectory")
    print("------------------------------------------------------------")
    print(">>")

    return Path(input())



""" Get Name of Result Directory AND FILE """

def getNameOfResultDirectory():
    """
    getNameOfResultDirectory()

    """

    print("------------------------------------------------------------")
    print("PLEASE ENTER NAME OF RESULT DIRECTORY ")
    print("")
    print("THE RESULT DATA WILL SAVED AT ")
    print("'.\\result\\qutrit\\iterative(or poisson)\\{ YOUR ENTED DIRECTORY NAME }\\{ EXPERIMENTAL DATA FILE NAME }_result.txt'")
    print("")
    print("IF EMPTY, THE NAME OF RESULT DIRECTORY IS 'default'")
    print("------------------------------------------------------------")
    print(">>")

    nameOfResultDirectory = input()
    if nameOfResultDirectory == "":
        nameOfResultDirectory = "default"

    return nameOfResultDirectory



""" Whether Do Poisson Distributed Simulation """

def checkPoisson():
    """
    checkPoisson()

    """

    print("------------------------------------------------------------")
    print("PLEASE ENTER ANSWER WHETHER DO POISSON DISTRIBUTED SIMULATION")
    print("IF YOU DO, PLEASE ENTER 'yes'")
    print("IF YOU ENTER ANOTHER WORD OR EMPTY, YOUR ANSWER IS REGARED AS 'no'")
    print("------------------------------------------------------------")
    print(">>")

    answer = input()
    if answer == "yes" or answer == "Yes" or answer == "YES":
        print("YOUR ANSWER IS: 'yes'")
        poissonPaths = getExperimentalDataPaths()
        eachIterationTime = getEachIterationTime()
        return True, poissonPaths*eachIterationTime
    else:
        print("YOUR ANSWER IS: 'no'")
        return False, []



""" Get Each Iteration Time """

def getEachIterationTime():
    """
    getEachIterationTime()

    """

    print("------------------------------------------------------------")
    print("PLEASE ENTER ITERATION TIME OF EACH POISSON SIMULATION")
    print("------------------------------------------------------------")
    print(">>")
    
    eachIterationTime = input()
    if eachIterationTime == "":
        eachIterationTime = 0
    else:
        eachIterationTime = int(eachIterationTime)

    return eachIterationTime



""" Get Number of Parallel Comuting """

def getNumberOfParallelComputing():
    """
    getNumberOfParallelComputing()

    """

    print("------------------------------------------------------------")
    print("HOW MANY TIMES DO YOU WANT TO PARALLELIZE?")
    print("IF THE NUMBER IS TOO LARGE, THE PARFORMANCE OF SIMULATION BECOME LOWER.")
    print("THE NUMBER OF LOGICAL PROCESSOR OF YOUR COMPUTER IS >>")
    print(os.cpu_count())
    print("RECOMENDED NUMBER IS LESS THAN THE ABOVE NUMBER.")
    print("------------------------------------------------------------")
    print(">>")
    
    n = input()
    if n != '':
        numberOfParallelComputing = int(n)
    else:
        numberOfParallelComputing = 1

    return numberOfParallelComputing



if __name__ == "__main__":

    """ Get Number of Qutrits """
    numberOfQutrits = getNumberOfQutrits()
    
    """ Get Paths of Experimental Data Directory """
    directoryPath = getExperimentalDataDirectoryPath()
    paths = list(directoryPath.glob("*.txt"))

    """ Get Name of Result Directory """
    resultDirectoryName = getNameOfResultDirectory()

    """ Check Poisson Distributed Simulation """
    check, poissonPaths = checkPoisson()

    """ Get Number of Parallel Computing """
    numberOfParallelComputing = getNumberOfParallelComputing()

    """ Make Bases """
    basesOfQutrits = makeBases(numberOfQutrits, bases)

    """ Make Ideal Density Matrix """
    baseVecter = np.zeros([1, 3**numberOfQutrits])
    baseVecter[0][0] = 1 / sqrt(3)
    baseVecter[0][4] = 1 / sqrt(3)
    baseVecter[0][3**numberOfQutrits-1] = 1 / sqrt(3)
    idealDensityMatrix = baseVecter.T @ baseVecter

    """ Make Result Directory """
    if not os.path.exists('.\\result\\qutrit\\iterative\\' + resultDirectoryName):
        os.makedirs('.\\result\\qutrit\\iterative\\' + resultDirectoryName)

    """ Start Tomography """
    with futures.ProcessPoolExecutor(max_workers=numberOfParallelComputing) as executor:
        for path in paths:
            executor.submit(fn=doIterativeSimulation, numberOfQutrits=numberOfQutrits, bases=basesOfQutrits, pathOfExperimentalData=str(path), idealDensityMatrix=idealDensityMatrix, resultDirectoryName=resultDirectoryName)

    """ Start Poisson Distributed Simulation """
    if check:
        """ Make Result Directory for Poisson Distributed Simulation """
        if not os.path.exists('.\\result\\qutrit\\posson\\' + resultDirectoryName):
            os.makedirs('.\\result\\qutrit\\poisson\\' + resultDirectoryName)

        with futures.ProcessPoolExecutor(max_workers=numberOfParallelComputing) as executor:
            for poissonPath in poissonPaths:
                executor.submit(fn=doPoissonDistributedSimulation, numberOfQutrits=numberOfQutrits, bases=basesOfQutrits, pathOfExperimentalData=poissonPath, idealDensityMatrix=idealDensityMatrix, resultDirectoryName=resultDirectoryName)

    


