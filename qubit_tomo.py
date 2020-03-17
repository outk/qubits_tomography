'''
This is a implementation of Quantum State Tomography for Qubits,
using techniques of following papars.

'Iterative algorithm for reconstruction of entangled states(10.1103/PhysRevA.63.040303)'
'Diluted maximum-likelihood algorithm for quantum tomography(10.1103/PhysRevA.75.042108)'
'Qudit Quantum State Tomography(10.1103/PhysRevA.66.012303)'

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


su2b = array([
    [[   1,  0], [   0,  1]],
    [[   0,  1], [   1,  0]],
    [[   0,-1j], [  1j,  0]],
    [[   1,  0], [   0, -1]]
]
)

su2Bases = []
newbases = su2b.copy()

def makeSU2Bases(numberOfQubits):
    global newbases, su2Bases
    for _ in range(numberOfQubits-1):
        for i in range(len(newbases)):
            su2Bases.extend([kron(newbases[i], su2b[j]) for j in range(4)])
        newbases = su2Bases.copy()
        su2Bases = []
        
    su2Bases = array(newbases) / (2**numberOfQubits)


bH = array([[1,0],[0,0]])
bV = array([[0,0],[0,1]])
bD = array([[1/2,1/2],[1/2,1/2]])
bR = array([[1/2,1j/2],[-1j/2,1/2]])
bL = array([[1/2,-1j/2],[1j/2,1/2]])

initialBases = array([bH, bV, bR, bD])

cycleBases1 = array([bH, bV, bR, bD])
cycleBases2 = array([bD, bR, bV, bH])


def makeBases(numberOfQubits):
    beforeBases = initialBases
    for _ in range(numberOfQubits - 1):
        afterBases = []
        for i in range(len(beforeBases)):
            if i % 2 == 0:
                afterBases.extend([kron(beforeBases[i], cycleBase) for cycleBase in cycleBases1])
            else:
                afterBases.extend([kron(beforeBases[i], cycleBase) for cycleBase in cycleBases2])
        beforeBases = afterBases
    
    baseName = []

    for i in range(2**numberOfQubits):
        if i%4 == 0:
            baseName.append("|"+ str(np.binary_repr(i, width=4)) +">")

    return array(afterBases), baseName



def makeBMatrix(numberOfQubits, bases):
    global su2Bases
    B = np.zeros((4**numberOfQubits, 4**numberOfQubits))

    for i in range(4**numberOfQubits):
        for j in range(4**numberOfQubits):
            B[i][j] = np.trace(bases[i] @ su2Bases[j])

    return B


def makeMMatrix(numberOfQubits, bases):
    global su2Bases
    B = makeBMatrix(numberOfQubits, bases)

    BInverse = np.linalg.inv(B)

    M = []

    for i in range(4**numberOfQubits):
        M.append(sum([BInverse[j][i] * su2Bases[j] for j in range(4**numberOfQubits)]))
    
    return array(M)


def makeInitialDensityMatrix(numberOfQubits, dataList, bases, M):
    # M = makeMMatrix(numberOfQubits, bases)

    N = sum([np.trace(M[i]) * dataList[i] for i in range(4**numberOfQubits)])

    densityMatrix = sum([dataList[i] * M[i] for i in range(4**numberOfQubits)]) / N

    initialDensityMatrix = choleskyDecomposition(numberOfQubits, densityMatrix)

    return initialDensityMatrix


""" cholesky decomposition """
def choleskyDecomposition(numberOfQubits, matrix):

    L = np.zeros([2**numberOfQubits, 2**numberOfQubits], dtype=np.complex)

    for i in range(2**numberOfQubits-1, -1, -1):
        s = matrix[i][i]
        for k in range(2**numberOfQubits-1, i, -1):
            s -= np.conjugate(L[k][i]) * L[k][i]
        if s >= 0:
            L[i][i] = np.sqrt(s)
        else:
            L[i][i] = np.sqrt(s)
        for j in range(i):
            t = matrix[i][j]
            for k in range(2**numberOfQubits-1, i, -1):
                t -= (np.conjugate(L[k][i]) * L[k][j])
            if L[i][i] != 0:
                L[i][j] = t / np.conjugate(L[i][i])
            else:
                L[i][j] = t / 1e-9

    for i in range(2**numberOfQubits):
        L[i][i] = np.real(L[i][i])

    return (np.conjugate(L).T @ L) / np.trace(np.conjugate(L).T @ L)


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



def doIterativeAlgorithm(numberOfQubits, bases, listOfExperimentalDatas, MMatrix):
    """
    doIterativeAlgorithm():

        This function is to do iterative algorithm(10.1103/PhysRevA.63.040303) and diluted MLE algorithm(10.1103/PhysRevA.75.042108) to a set of datas given from a experiment.
        This recieve four variables (numberOfQubits, bases, listAsExperimentalDatas),
        and return most likely estimated density matrix (np.array) and total time of calculation(datetime.timedelta).

        First quantum state matrix for this algorithm is a identity matrix.


        --------------------------------------------------------------------------------------------------------------
        Return:

            most likely estimated density matrix(np.array), 

    """

    """ Setting initial parameters """
    iter = 0
    epsilon = 1000
    endDiff = 1e-11
    diff = 100
    # TolFun = 10e-11
    # traceDistance = 100
    maxNumberOfIteration = 100000

    dataList = listOfExperimentalDatas
    totalCountOfData = sum(dataList)
    nDataList = dataList / totalCountOfData # nDataList is a list of normarized datas
    densityMatrix = makeInitialDensityMatrix(numberOfQubits, dataList, bases, MMatrix)
    # densityMatrix = identity(2 ** numberOfQubits)

    """ Start iteration """
    # while traceDistance > TolFun and iter <= maxNumberOfIteration:
    while diff > endDiff and iter <= maxNumberOfIteration and epsilon > 1e-6:

        probList = [trace(base @ densityMatrix) for base in bases]
        nProbList = probList / sum(probList)
        rotationMatrix = sum([(nDataList[i] / probList[i]) * bases[i] for i in range(4 ** numberOfQubits) if probList[i] != 0])

        """ Normalization of Matrices for Measurement Bases """
        U = np.linalg.inv(sum(bases)) / sum(probList)   
        rotationMatrixLeft = (identity(2 ** numberOfQubits) + epsilon * U @ rotationMatrix) / (1 + epsilon)
        rotationMatrixRight = (identity(2 ** numberOfQubits) + epsilon * rotationMatrix @ U) / (1 + epsilon)

        """ Calculation of updated density matrix """
        modifiedDensityMatrix = rotationMatrixLeft @ densityMatrix @ rotationMatrixRight / trace(rotationMatrixLeft @ densityMatrix @ rotationMatrixRight)
        eigValueArray, eigVectors = np.linalg.eig(densityMatrix - modifiedDensityMatrix)
        traceDistance = sum(np.absolute(eigValueArray)) / 2

        """ Update Likelihood Function, and Compared with older one """
        LikelihoodFunction = sum([nDataList[i] * np.log(nProbList[i]) for i in range(4 ** numberOfQubits) if nProbList[i] != 0])
        probList = [trace(base @ modifiedDensityMatrix) for base in bases]
        nProbList = probList / sum(probList)
        modifiedLikelihoodFunction = sum([nDataList[i] * np.log(nProbList[i]) for i in range(4 ** numberOfQubits) if nProbList[i] != 0])
        nowdiff = np.real(modifiedLikelihoodFunction - LikelihoodFunction)

        """ Show Progress of Calculation """
        progress = 100 * iter / maxNumberOfIteration
        if progress % 5 == 0:
            msg = "Progress of calculation: " + str(int(progress)) + "%"
            print(msg)

        """ Increment """
        iter += 1

        """ Check Increasing of Likelihood Function  """
        if nowdiff < 0:
            epsilon = epsilon * 0.1
            continue
        else:
            diff = nowdiff
        
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
    print("Iteration was '" + str(endIterationTimes) + "' times.")

    return modifiedDensityMatrix, endIterationTimes



""" Calculate Fidelity """

def calculateFidelity(idealDensityMatrix, estimatedDensityMatrix):
    """
    calculateFidelity(idealDensityMatrix, estimatedDensityMatrix):

    """

    fidelity = np.real(trace(sqrtm(sqrtm(idealDensityMatrix) @ estimatedDensityMatrix @ sqrtm(idealDensityMatrix)))) ** 2

    return fidelity



""" Iterative Simulation """

def doIterativeSimulation(numberOfQubits, bases, pathOfExperimentalData, idealDensityMatrix, resultDirectoryName, MMatrix, baseNames):
    """
    doIterativeSimulation(numberOfQubits, bases, pathOfExperimentalData, idealDensityMatrix, resultDirectoryName, MMatrix)


    """

    """ Get Experimental Data"""
    listOfExperimentalData = getExperimentalData(pathOfExperimentalData)

    """ Calculate """
    estimatedDensityMatrix, totalIterationTime = doIterativeAlgorithm(numberOfQubits, bases, listOfExperimentalData, MMatrix)
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
    resultFilePath = '.\\result\\qubit\\iterative\\' + resultDirectoryName + '\\' + 'result' + '.txt'
    resultIterationTimeFilePath = '.\\result\\qubit\\iterative\\' + resultDirectoryName + '\\' + 'resultIterationTime' + '.txt'

    """ Save Result """
    with open(resultFilePath, mode='a') as f:
        f.writelines(str(fidelity) + '\n')
    with open(resultIterationTimeFilePath, mode='a') as f:
        f.writelines(str(totalIterationTime) + '\n')

    """ Make 3D Plot """
    plotResult(numberOfQubits, estimatedDensityMatrix, baseNames)



""" Poisson Distributed Simulation """

def doPoissonDistributedSimulation(numberOfQubits, bases, pathOfExperimentalData, idealDensityMatrix, resultDirectoryName, MMatrix):
    """
    doPoissonDistributedSimulation(numberOfQubits, bases, pathOfExperimentalData, idealDensityMatrix, resultDirectoryName, MMatrix)


    """

    """ Get Experimental Data"""
    listOfExperimentalData = getExperimentalData(pathOfExperimentalData)

    """ Calculate """
    estimatedDensityMatrix, totalIterationTime = doIterativeAlgorithm(numberOfQubits, bases, random.poisson(listOfExperimentalData), MMatrix)
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
    resultFilePath = '.\\result\\qubit\\poisson\\' + resultDirectoryName + "\\" + resultFileName + '_result' + '.txt'

    """ Save Result """
    with open(resultFilePath, mode='a') as f:
        f.write(str(fidelity) + '\n')



def plotResult(numberOfQubits, densityMatrix, baseNames):
    """
    plotResult(densityMatrix)
    
    """

    """ Plot Setting """
    xedges, yedges = np.arange(2**numberOfQubits), np.arange(2**numberOfQubits)
 
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
    plt.xticks(np.arange(0, 2**numberOfQubits, 4), labels=baseNames)
    # plt.ylabel("Y") # y軸の内容表示
    plt.yticks(np.arange(0, 2**numberOfQubits, 4), labels=baseNames)
    # ax1.set_zlabel("Z") # z軸の内容表示
    ax1.set_zlim(-0.1, 0.6)
    
    ax2 = fig.add_subplot(122, projection="3d") # 3Dの軸を作成
    ax2.bar3d(xpos,ypos,zpos,dx,dy,np.imag(dz), edgecolor='black') # ヒストグラムを3D空間に表示
    plt.title("Imaginary Part") # タイトル表示
    # plt.xlabel("X") # x軸の内容表示
    plt.xticks(np.arange(0, 2**numberOfQubits, 4), labels=baseNames)
    # plt.ylabel("Y") # y軸の内容表示
    plt.yticks(np.arange(0, 2**numberOfQubits, 4), labels=baseNames)
    # ax2.set_zlabel("Z") # z軸の内容表示
    ax2.set_zlim(-0.1, 0.6)
    
    plt.show()

    with open('firstplottest'+'_plot.pkl', mode='wb') as f:
        pickle.dump(fig, f)


""" Get Number of Qubits """

def getNumberOfQubits():
    """
    getNumberOfQubits()

    """

    print("------------------------------------------------------------")
    print("PLEASE ENTER NUMBER OF QUBITS")
    print("------------------------------------------------------------")
    print(">>")

    numberOfQubits = int(input())

    return numberOfQubits



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



""" Get Paths of Experimental Data """

def getExperimentalDataPaths():
    """
    getExperimentalDataPaths()

    """

    print("------------------------------------------------------------")
    print("PLEASE ENTER PATHS OF EXPERIMENTAL DATA")
    print("")
    print("IF THERE ARE MULTIPLE DATA FILE YOU WANT TO TOMOGRAPHY,")
    print("ENTER ALL PATHS SEPARATED WITH SPACE.")
    print("LIKE THIS >> .\\datadirectory\\ex1.txt .\\datadirectory\\ex2.txt ...")
    print("------------------------------------------------------------")
    print(">>")

    paths = list(input().split())

    return paths



""" Get Name of Result Directory AND FILE """

def getNameOfResultDirectory():
    """
    getNameOfResultDirectory()

    """

    print("------------------------------------------------------------")
    print("PLEASE ENTER NAME OF RESULT DIRECTORY ")
    print("")
    print("THE RESULT DATA WILL SAVED AT ")
    print("'.\\result\\qubit\\iterative(or poisson)\\{ YOUR ENTED DIRECTORY NAME }\\{ EXPERIMENTAL DATA FILE NAME }_result.txt'")
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

    """ Get Number of Qubits """
    numberOfQubits = getNumberOfQubits()

    """ Make SU2 Bases """
    makeSU2Bases(numberOfQubits)   
    
    """ Get Path of Experimental Data Directory """
    directoryPath = getExperimentalDataDirectoryPath()
    paths = list(directoryPath.glob("*.txt"))

    """ Get Name of Result Directory """
    resultDirectoryName = getNameOfResultDirectory()

    """ Check Poisson Distributed Simulation """
    check, poissonPaths = checkPoisson()

    """ Get Number of Parallel Computing """
    numberOfParallelComputing = getNumberOfParallelComputing()

    """ Make Bases """
    basesOfQubits, baseNames = makeBases(numberOfQubits)

    """ Make M Matrix """
    M = makeMMatrix(numberOfQubits, basesOfQubits)

    """ Make Ideal Density Matrix """
    baseVecter = np.zeros([1, 2**numberOfQubits])
    baseVecter[0][0] = 1 / sqrt(2)
    baseVecter[0][2**numberOfQubits-1] = 1 / sqrt(2)
    idealDensityMatrix = baseVecter.T @ baseVecter

    # baseVecter[0][1] = 1 / 2
    # baseVecter[0][2] = 1 / 2
    # baseVecter[0][4] = 1 / 2
    # baseVecter[0][8] = 1 / 2
    # baseVecter = np.full([1, 2**numberOfQubits], 1/np.sqrt(2**numberOfQubits), dtype=np.complex)
    # idealDensityMatrix = baseVecter.T @ baseVecter

    # matrix = np.zeros([2**numberOfQubits, 2**numberOfQubits]) # (|0000>+|1111>)(<0000|+<1111|) + |0001><0001| + |0010><0010| + |0100><0100| + |1000><1000|
    # baseVecter = np.zeros([1, 2**numberOfQubits])
    # baseVecter[0][1] = 1
    # matrix += baseVecter.T @ baseVecter
    # baseVecter = np.zeros([1, 2**numberOfQubits])
    # baseVecter[0][2] = 1
    # matrix += baseVecter.T @ baseVecter
    # baseVecter = np.zeros([1, 2**numberOfQubits])
    # baseVecter[0][4] = 1
    # matrix += baseVecter.T @ baseVecter
    # baseVecter = np.zeros([1, 2**numberOfQubits])
    # baseVecter[0][8] = 1
    # matrix += baseVecter.T @ baseVecter
    # idealDensityMatrix = baseVecter.T @ baseVecter

    # baseVecter = np.zeros([1, 2**numberOfQubits])
    # baseVecter[0][0] = 1
    # baseVecter[0][2**numberOfQubits-1] = 1
    # matrix += baseVecter.T @ baseVecter
    # matrix = matrix/np.trace(matrix)
    # idealDensityMatrix = matrix

    """ Make Result Directory """
    if not os.path.exists('.\\result\\qubit\\iterative\\' + resultDirectoryName):
        os.makedirs('.\\result\\qubit\\iterative\\' + resultDirectoryName)

    """ Start Tomography """
    with futures.ProcessPoolExecutor(max_workers=numberOfParallelComputing) as executor:
        for path in paths:
            executor.submit(fn=doIterativeSimulation, numberOfQubits=numberOfQubits, bases=basesOfQubits, pathOfExperimentalData=str(path), idealDensityMatrix=idealDensityMatrix, resultDirectoryName=resultDirectoryName, MMatrix=M, baseNames=baseNames)

    """ Start Poisson Distributed Simulation """
    if check:
        """ Make Result Directory for Poisson Distributed Simulation """
        if not os.path.exists('.\\result\\qubit\\poisson\\' + resultDirectoryName):
            os.makedirs('.\\result\\qubit\\poisson\\' + resultDirectoryName)

        with futures.ProcessPoolExecutor(max_workers=numberOfParallelComputing) as executor:
            for poissonPath in poissonPaths:
                executor.submit(fn=doPoissonDistributedSimulation, numberOfQubits=numberOfQubits, bases=basesOfQubits, pathOfExperimentalData=poissonPath, idealDensityMatrix=idealDensityMatrix, resultDirectoryName=resultDirectoryName, MMatrix=M)


    # if not os.path.exists('.\\result\\4qubit\\poisson\\benchmark'):
    #         os.makedirs('.\\result\\4qubit\\poisson\\benchmark')

    # with open('benchmark'+str(numberOfQubits)+'qubits.txt', mode='a') as f:
    #     f.write("total time: " + str(end_time - start_time) + "\n")

