#!/usr/bin/env python3 
import math,pickle,os,sys,json,time,glob,string,subprocess
from sys import platform as sys_pf
if sys_pf == 'darwin':
    import matplotlib
    matplotlib.use("TkAgg")
import tkinter as tk
import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns
import pandas as pd
from sklearn.preprocessing import StandardScaler
from modifyDataFrames import returnModifiedDf
from miscFunctions import reindexDataFrame,printProgressBar,extractValues,reorderDfByInputOrder,returnSpecificExtensionFiles
idx = pd.IndexSlice

plateRowLetters = list(string.ascii_uppercase)[:16]
plateColumnNumbers = list(range(1,25))

def produceSingleCellHeaders(cellTypes):
    newMultiIndexList = []
    for cellType in cellTypes:
        newMultiIndexList.append([cellType])
    return newMultiIndexList

def grabCellTypeList(experimentParameters):
    path = 'inputData/singleCellCSVFiles/'
    if experimentParameters['format'] == 'plate':
        nrows = experimentParameters['overallPlateDimensions'][0]
        ncols = experimentParameters['overallPlateDimensions'][1]
        if 'unpackingDict' in experimentParameters.keys():
            orderingList = [str(x).zfill(3) for x in range(1,385)]
        else:
            orderingList = [str(x).zfill(3) for x in range(1,nrows*ncols+1)]

        #Plate case
        #Walk through a plate, grab each unique cell type
        cellTypeList = []
        folder = ''
        for fileName in os.listdir(path):
            if '.DS' not in fileName:
                folder = fileName
        
        for fileName in os.listdir(path+folder+'/'):
            if '.DS' not in fileName:
                splitFile = fileName.split('_')
                parsingVal = ''
                for split in splitFile:
                    if split in orderingList:
                        parsingVal = split
                parsingPose = fileName.rindex(parsingVal)+4
                cellType = fileName[parsingPose:].split('.')[0].replace('/','\/')
                if cellType not in cellTypeList:
                    cellTypeList.append(cellType)
    else:
        cellTypeList = []
        sampleNameDf = pd.read_pickle('misc/tubeLayout-'+os.getcwd().split('/')[-1]+'-cell.pkl')
        fullSampleFileName = sampleNameDf.iloc[0,0]
        dotIndex = fullSampleFileName.rfind('.')
        sampleFileName = fullSampleFileName[:dotIndex]
        for fileName in os.listdir(path):
            if '.DS' not in fileName and sampleFileName in fileName:
                splitFile = fileName.split(sampleFileName[1:]+'_')
                cellType = splitFile[1].split('.')[0]
                if cellType not in cellTypeList:
                    cellTypeList.append(cellType)

    return cellTypeList

def createTubeSingleCellDataFrame(folderName,experimentParameters,fileNameDf):
    
    path = 'inputData/singleCellCSVFiles/'
    
    fileNameDf = fileNameDf.stack().to_frame('fileName')

    #Grab common prefix for files
    fullFileName = ''
    for fileName in os.listdir(path):
        if '.DS' not in fileName:
            print(fileName)
            fullFileName = fileName
            break
    prefix = fullFileName.split('_')[0]

    #Remove extraneous decorations in marker names
    newColumns = []
    fcsDf = pd.read_csv(path+fullFileName,header=0)
    for column in fcsDf.columns:
        if '::' in column:
            column = column.split(' :: ')[1]
        #CyTOF
        if column[0].isdigit():
            column = column.split('_')[1]
        newColumns.append(column)

    cellTypeList = grabCellTypeList(experimentParameters)
    
    completeDfList = []
    for cellType in cellTypeList:
        for row in range(fileNameDf.shape[0]):
            fullFileName = fileNameDf.iloc[row,0]
            dotIndex = fullFileName.rfind('.')
            fileName = fullFileName[:dotIndex]
            trueFileName = prefix+'_'+fileName+'_'+cellType+'.csv'
            levelValues = list(fileNameDf.iloc[row,:].name) 
            fcsDf = pd.read_csv(path+trueFileName,header=0)
            eventNumber = fcsDf.shape[0]
            eventList = range(1,eventNumber+1)
            allLevelValues = []
            for event in eventList:
                allLevelValues.append([cellType]+levelValues+[event])
            newMultiIndex = pd.MultiIndex.from_tuples(allLevelValues,names=['CellType']+list(fileNameDf.index.names)+['Event'])
            newDf = pd.DataFrame(fcsDf.values,index=newMultiIndex,columns=newColumns)
            completeDfList.append(newDf)
            printProgressBar(row + 1, fileNameDf.shape[0], prefix = ' Concatenating samples:', suffix = 'Complete', length = 50)

    completeDataFrame = pd.concat(completeDfList)
    completeDataFrame.columns.name = 'Marker'

    #Remove extraneous markers (namely -h parameters)
    columnsToKeep = []
    for col,column in enumerate(completeDataFrame.columns):
        if '-H' not in column and 'Time' not in column and '-W' not in column:
            columnsToKeep.append(col)
    completeDataFrame = completeDataFrame.iloc[:,columnsToKeep]

    completeDataFrame.to_hdf('outputData/pickleFiles/initialSingleCellDf-channel-'+folderName+'.h5', key='df', mode='w')
    print(completeDataFrame)

#If multiplexing option chosen
def demultiplexSingleCellData(experimentParameters):
    #"multiplexingOption": "96->384 well", "unpackingDict": {"A1-2_B1-2": ["A1", "A2", "B2", "B1"]
    unpackingDict = experimentParameters['unpackingDict']
    unpackingPositionDict = {(0,0):0,(0,1):1,(1,0):2,(1,1):3}
    cellTypeList = grabCellTypeList(experimentParameters)

    #Currently in format A1; A2, B1; B2
    #Need to change to format A1: B1, A2, B2
    plateRowLetters = string.ascii_uppercase[:16]
    plateColumnNumbers = list(range(1,25))

    wellPlateRowLetters = string.ascii_uppercase[:8]
    wellPlateColumnNumbers = list(range(1,13))

    #Create appropriate folders in each population
    for combinedPlateName in list(unpackingDict.keys()):
        for unpackedPlateName in unpackingDict[combinedPlateName]:
            if unpackedPlateName != '':
                if unpackedPlateName not in returnSpecificExtensionFiles('inputData/singleCellCSVFiles/','',False):
                    subprocess.run(['mkdir','inputData/singleCellCSVFiles/'+unpackedPlateName])
    fileNameDict = {}
    combinedPlateNames = list(unpackingDict.keys())
    for combinedPlateName in combinedPlateNames:
        for populationName in cellTypeList:
            #scale_Specimen_001_P9_P09_369_TCells.csv
            allFileNames = returnSpecificExtensionFiles('inputData/singleCellCSVFiles/'+combinedPlateName,'',False)
            unpackedPlateNames = unpackingDict[combinedPlateName]
            for k,fileName in enumerate(allFileNames):
                sampleID = fileName.split('_')[3]
                currentRowLetter = sampleID[0]
                currentColumnNumber = int(sampleID[1:])
                #Get index of current row and column position
                currentRowLetterIndex = plateRowLetters.index(currentRowLetter)
                currentColumnNumberIndex = plateColumnNumbers.index(currentColumnNumber)
                #Demultiplex sample ids 
                wellPlateRowLetter = wellPlateRowLetters[int(currentRowLetterIndex/2)]
                wellPlateColumnNumber = wellPlateColumnNumbers[int(currentColumnNumberIndex/2)]
                newSampleID = str(wellPlateRowLetter)+str(wellPlateColumnNumber)
                newSampleID2 = str(wellPlateRowLetter)+str(wellPlateColumnNumber).zfill(2)
                newFileName = '_'.join(['_'.join(fileName.split('_')[:3]),newSampleID,newSampleID2,'_'.join(fileName.split('_')[-2:])])
                unpackedPlateIndex = unpackingPositionDict[(currentRowLetterIndex%2,currentColumnNumberIndex%2)]
                unpackedFolder = unpackedPlateNames[unpackedPlateIndex]
                trueFileName = newFileName
                #print(fileName+'->'+newFileName)
                fileNameDict['_'.join(trueFileName.split('_')[1:-1])] = '_'.join(fileName.split('_')[1:-1])
                completeNewFileName = unpackedFolder+'/'+trueFileName
                subprocess.run(['cp','inputData/singleCellCSVFiles/'+combinedPlateName+'/'+fileName,'inputData/singleCellCSVFiles/'+completeNewFileName])
                printProgressBar(k + 1, len(allFileNames), prefix = ' Demultiplexing '+combinedPlateName+','+populationName+':', suffix = 'Complete', length = 50)
    with open('misc/fileNameDict.pkl','wb') as f:
        pickle.dump(fileNameDict,f)

def createPlateSingleCellDataFrame(folderName,experimentParameters,levelLayout):
    
    path = 'inputData/singleCellCSVFiles/'
    if 'unpackingDict' in experimentParameters:
        demultiplexSingleCellData(experimentParameters)
    
    nrows = experimentParameters['overallPlateDimensions'][0]
    ncols = experimentParameters['overallPlateDimensions'][1]

    if 'unpackingDict' in experimentParameters.keys():
        orderingList = [str(x).zfill(3) for x in range(1,385)]
    else:
        orderingList = [str(x).zfill(3) for x in range(1,nrows*ncols+1)]

    completeKeyMatrix = np.dstack(list(levelLayout['keys'].values()))
    unraveledKeyMatrix = np.reshape(completeKeyMatrix,(completeKeyMatrix.shape[0]*completeKeyMatrix.shape[1],completeKeyMatrix.shape[2]))
    unraveledBlankMatrix = levelLayout['blank'].ravel()
    #print(unraveledBlankMatrix)

    sampleIndex = pd.MultiIndex.from_arrays([levelLayout['plateID'].ravel(),levelLayout['wellID'].ravel()],names=['Plate','Well'])
    sampleKeyDf = pd.DataFrame(unraveledKeyMatrix,index=sampleIndex,columns=list(experimentParameters['levelLabelDict'].keys()))
    sampleDf = sampleKeyDf.copy()
    rowsToKeep = []
    
    #print(sampleKeyDf.iloc[:,3].values.ravel())
    for row in range(sampleDf.shape[0]):
        print('row')
        print(row)
        for col in range(sampleDf.shape[1]):
            level = list(experimentParameters['levelLabelDict'].keys())[col]
            print(level)
            levelValueIndex = sampleKeyDf.iloc[row,col]
            if unraveledBlankMatrix[row] == -1:
                print(experimentParameters['levelLabelDict'][level])
                levelValue = experimentParameters['levelLabelDict'][level][levelValueIndex]
                sampleDf.iloc[row,col] = levelValue
            else:
                sampleDf.iloc[row,col] = 'Blank'
    #Drop blanks
    sampleDf = sampleDf.query("Time != 'Blank'")
    sampleDf.to_excel('outputData/excelFiles/fcsLabelingKey.xlsx')
    
    cellTypeList = grabCellTypeList(experimentParameters)

    newSampleDfList = []
    for cellType in cellTypeList:
        #Create dict that relates sample file locations to plate/wellIDs
        orderValDict = {}
        for plateName in list(np.unique(levelLayout['plateID'])):
            orderValDict2 = {}
            if 'DS' not in plateName:
                for fileName in os.listdir(path+plateName+'/'):
                    if 'DS' not in fileName:
                        splitFile = fileName.split('_')
                        parsingVal = ''
                        for split in splitFile:
                            if split in orderingList:
                                parsingVal = split
                        parsingPose = fileName.rindex(parsingVal)
                        currentCellType = fileName[parsingPose+4:].split('.')[0]
                        if currentCellType == cellType:
                            orderVal = fileName[parsingPose:parsingPose+3]
                            wellID = fileName[:parsingPose-1].split('_')[-2]
                            orderValDict2[wellID] = path+plateName+'/'+fileName
            orderValDict[plateName] = orderValDict2
        sampleTupleList,sampleList = [],[]
        for row in range(sampleDf.shape[0]):
            sampleID = list(sampleDf.iloc[row,:].name)
            plateID = sampleID[0]
            wellID = sampleID[1]
            sampleFileName = orderValDict[plateID][wellID]
            sampleList.append(sampleFileName)
            if type(cellType) != list:
                if type(cellType) == tuple:
                    cellType = list(cellType)
                else:
                    cellType = [cellType]
            sampleTuple = cellType+sampleDf.loc[idx[plateID,wellID],:].values.tolist()
            sampleTupleList.append(sampleTuple)
        sampleMI = pd.MultiIndex.from_tuples(sampleTupleList,names=['CellType']+list(sampleDf.columns))
        newSampleDf = pd.DataFrame(sampleList,index=sampleMI,columns=['fileName'])
        newSampleDfList.append(newSampleDf)
    
    fileNameDf = pd.concat(newSampleDfList)

    bothBool = False
    fullFileName = fileNameDf.iloc[0,0]
    fcsDf = pd.read_csv(fullFileName,header=0)
    for column in fcsDf.columns:
        if '::' in column:
            bothBool = True
    
    completeDfList = []
    for row in range(fileNameDf.shape[0]):
        fullFileName = fileNameDf.iloc[row,0]
        levelValues = list(fileNameDf.iloc[row,:].name) 
        fcsDf = pd.read_csv(fullFileName,header=0)
        eventNumber = fcsDf.shape[0]
        if eventNumber == 0:
            tk.messagebox.showerror("Error", "Filename:\n"+fullFileName+"\nhas no events. Please re-export this file and try again.")
            sys.exit(0)
        eventList = range(1,eventNumber+1)
        allLevelValues = []
        for event in eventList:
            allLevelValues.append(levelValues+[event])
        newMultiIndex = pd.MultiIndex.from_tuples(allLevelValues,names=list(fileNameDf.index.names)+['Event'])
        if bothBool:
            newColumns = [x.split(' :: ')[1] if '::' in x else x for x in fcsDf.columns]
        else:
            newColumns = fcsDf.columns
        newDf = pd.DataFrame(fcsDf.values,index=newMultiIndex,columns=newColumns)
        completeDfList.append(newDf)
        printProgressBar(row + 1, fileNameDf.shape[0], prefix = ' Concatenating samples:', suffix = 'Complete', length = 50)

    completeDataFrame = pd.concat(completeDfList)
    completeDataFrame.columns.name = 'Marker'

    #Remove extraneous markers (namely -h parameters)
    columnsToKeep = []
    for col,column in enumerate(completeDataFrame.columns):
        if '-H' not in column and 'Time' not in column and '-W' not in column:
            columnsToKeep.append(col)
    completeDataFrame = completeDataFrame.iloc[:,columnsToKeep]

    completeDataFrame.to_hdf('outputData/pickleFiles/initialSingleCellDf-channel-'+folderName+'.h5', key='df', mode='w')
    print(completeDataFrame)
