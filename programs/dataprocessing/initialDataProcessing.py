#!/usr/bin/env python3 
# -*- coding: utf-8 -*-
"""
created on sat jul 21 13:12:56 2018

@author: acharsr
"""
import json,pickle,math,matplotlib,sys,os,string,subprocess
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit
from itertools import groupby
from miscFunctions import Hill,InverseHill,r_squared,cleanUpFlowjoCSV,extractValues,reorderDfByInputOrder,returnSpecificExtensionFiles
from modifyDataFrames import returnModifiedDf
import cytokineDataProcessing,singleCellDataProcessing,cellDataProcessing

dataTypeLevelNames = {'cyt':['Cytokine'],'cell':['CellType','Marker','Statistic'],'prolif':['Statistic'],'singlecell':['CellType']}
dataTypeDataFrameFileNames = {'cyt':'cytokineConcentrationPickleFile','cell':'cellStatisticPickleFile','prolif':'proliferationStatisticPickleFile','singlecell':'initialSingleCellPickleFile'}
plateRowLetters = list(string.ascii_uppercase)[:16]
plateColumnNumbers = list(range(1,25))

def returnMultiIndex(sortedData,sortedFiles,dataType,folderName):
    if(dataType == 'cyt'):
        newMultiIndex = cytokineDataProcessing.parseCytokineCSVHeaders(pd.read_csv('inputData/bulkCSVFiles/A1_'+dataType+'.csv').columns)
    elif(dataType == 'cell'):
        if 'antibodyPanel-'+folderName+'.csv' in os.listdir('misc'):
            panelData = pd.read_csv('misc/antibodyPanel-'+folderName+'.csv',)
            newMultiIndex = cellDataProcessing.parseCellCSVHeaders(pd.read_csv('inputData/bulkCSVFiles/A1_'+dataType+'.csv').columns,panelData=panelData)
        else:
            newMultiIndex = cellDataProcessing.parseCellCSVHeaders(pd.read_csv('inputData/bulkCSVFiles/A1_'+dataType+'.csv').columns)
    elif(dataType == 'singlecell'):
        #Grabs a file from samples to read marker names off of
        cellTypeList = []
        for fileName in os.listdir('inputData/singleCellCSVFiles/A1/'):
            if 'DS' not in fileName:
                cellTypeList.append(fileName)
        newMultiIndex = singleCellDataProcessing.produceSingleCellHeaders(cellTypeList)
    elif(dataType == 'cytcorr'):
        newMultiIndex = []
    if dataType != 'singlecell':
        return sortedData,newMultiIndex
    else:
        return sortedFiles,newMultiIndex

def unpackMultiplexedPlates(experimentParameters,folderName,dataType):
    #A1->A1,A2->A2,A3- >A1,B1->A4,B2->A3
    #1.1,1.3,1.5,3.1,3.3,3.5 -> A1; 1.2,1.4,1.6,3.2,3.4,3.6->A2; 2.1,2.3,3.5,4.1,4.3,4.5->A4; 2.2,2.4,2.6,4.2,4.4,4.6->A3
    #1-24
    #A-P
    wellIDConversionDict = {}
    for i,plateLetter in enumerate(plateRowLetters):
        for j,plateNumber in enumerate(plateColumnNumbers):
            wellIDConversionDict[plateLetter+str(plateNumber)] = plateRowLetters[int(i/2)]+str(plateColumnNumbers[int(j/2)])
    plateIDConversionDict = {}
    for multiplexedPlateName in experimentParameters['unpackingDict']:
        for i,multiplexedWellPos in enumerate(experimentParameters['unpackingDict'][multiplexedPlateName]):
            wellIDs = []
            if multiplexedWellPos != '':
                if i < 2:
                    j = 0
                else:
                    j = 1
                for rowIndex in range(j,len(plateRowLetters),2):
                    for colIndex in range(i%2,len(plateColumnNumbers),2):
                        wellID = plateRowLetters[rowIndex]+str(plateColumnNumbers[colIndex])
                        wellIDs.append(wellID)
            plateIDConversionDict[multiplexedWellPos] = wellIDs
    #Specimen_001_A1_A01_001.fcs
    if experimentParameters['multiplexingOption'] == '96->384 well':
        multiplexedPlateNames = experimentParameters['unpackingDict'].keys()
        #sortedMultiplexedData,sortedMultiplexedFiles = cleanUpFlowjoCSV(multiplexedPlateNames,folderName,dataType)
        for multiplexedPlateName in multiplexedPlateNames:
            multiplexedWellPoses = experimentParameters['unpackingDict'][multiplexedPlateName]
            print(multiplexedWellPoses)
            with open('inputData/bulkCSVFiles/'+multiplexedPlateName+'_'+dataType+'.csv', 'r') as f:
                multiplexedCSVLines = f.readlines()
            for multiplexedWellPos in multiplexedWellPoses:
                if multiplexedWellPos != '':
                    wellIDsInThisPos = plateIDConversionDict[multiplexedWellPos]
                    linesToMove = []
                    newCSVLines = []
                    for lineNum,line in enumerate(multiplexedCSVLines):
                        if lineNum in [0,len(multiplexedCSVLines)-2,len(multiplexedCSVLines)-1]:
                            newCSVLines.append(line)
                        else:
                            fileName = line.split(',')[0]
                            wellID = fileName.split('_')[2]
                            #Well ID is in pos or first line or last two lines
                            if wellID in wellIDsInThisPos:
                                #underscorePoses = [pos for pos, char in enumerate(line) if char == '_']
                                newWellID = wellIDConversionDict[wellID]
                                newLine = line.replace(wellID,newWellID)
                                newCSVLines.append(newLine)
                    print(multiplexedWellPos)
                    with open('inputData/bulkCSVFiles/'+multiplexedWellPos+'_'+dataType+'.csv', 'w') as f:
                        for item in newCSVLines:
                            f.write("%s" % item)

#If multiplexing option chosen
def demultiplexSingleCellData(experimentParameters):
    #"multiplexingOption": "96->384 well", "unpackingDict": {"A1-2_B1-2": ["A1", "A2", "B2", "B1"]
    unpackingDict = experimentParameters['unpackingDict']
    unpackingPositionDict = {(0,0):0,(0,1):1,(1,0):2,(1,1):3}

    #Currently in format A1; A2, B1; B2
    #Need to change to format A1: B1, A2, B2
    plateRowLetters = string.ascii_uppercase[:16]
    plateColumnNumbers = list(range(1,25))

    wellPlateRowLetters = string.ascii_uppercase[:8]
    wellPlateColumnNumbers = list(range(1,13))

    #Create appropriate folders in each population
    for combinedPlateName in list(unpackingDict.keys()):
        for unpackedPlateName in unpackingDict[combinedPlateName]:
            for populationName in returnSpecificExtensionFiles('inputData/singleCellCSVFiles/'+combinedPlateName,'',False):
                if unpackedPlateName not in returnSpecificExtensionFiles('inputData/singleCellCSVFiles/','',False):
                    subprocess.run(['mkdir','inputData/singleCellCSVFiles/'+unpackedPlateName])
                    if populationName not in returnSpecificExtensionFiles('inputData/singleCellCSVFiles/'+unpackedPlateName,'',False):
                        subprocess.run(['mkdir','inputData/singleCellCSVFiles/'+unpackedPlateName+'/'+populationName])
                else:
                    if populationName not in returnSpecificExtensionFiles('inputData/singleCellCSVFiles/'+unpackedPlateName,'',False):
                        subprocess.run(['mkdir','inputData/singleCellCSVFiles/'+unpackedPlateName+'/'+populationName])
    fileNameDict = {}
    for combinedPlateName in list(unpackingDict.keys()):
        for populationName in returnSpecificExtensionFiles('inputData/singleCellCSVFiles/'+combinedPlateName,'',False):
            #scale_Specimen_001_P9_P09_369_TCells.csv
            allFileNames = returnSpecificExtensionFiles('inputData/singleCellCSVFIles/'+combinedPlateName+'/'+populationName,'',False)
            unpackedPlateNames = unpackingDict[combinedPlateName]
            for fileName in allFileNames:
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
                newFileName = '_'.join(['_'.join(fileName.split('_')[:3]),newSampleID,'_'.join(fileName.split('_')[-3:])])
                unpackedPlateIndex = unpackingPositionDict[(currentRowLetterIndex%2,currentColumnNumberIndex%2)]
                unpackedFolder = unpackedPlateNames[unpackedPlateIndex]
                
                realFileNameCSV = pd.read_csv('inputData/bulkCSVFiles/'+unpackedFolder+'_cell.csv',header=0)
                fileNameColumn = realFileNameCSV.iloc[:-2,0]
                indexList = []
                for row in range(fileNameColumn.shape[0]):
                    indexList.append(fileNameColumn.iloc[row].split('_')[2])
                fileIndex = indexList.index(newSampleID)
                trueFileName = newFileName.split('_')[0]+'_'+fileNameColumn.iloc[fileIndex][:-4]+'_'+populationName+'.csv'
                print(newFileName+'->'+trueFileName)
                fileNameDict['_'.join(trueFileName.split('_')[1:-1])] = '_'.join(fileName.split('_')[1:-1])
                completeNewFileName = unpackedFolder+'/'+populationName+'/'+trueFileName
                subprocess.run(['cp','inputData/singleCellCSVFiles/'+combinedPlateName+'/'+populationName+'/'+fileName,'inputData/singleCellCSVFiles/'+completeNewFileName])
    with open('misc/fileNameDict.pkl','wb') as f:
        pickle.dump(fileNameDict,f)

def createBaseDataFrame(experimentParameters,folderName,experimentNumber,dataType,layoutDict):
    if experimentParameters['format'] == 'tube':
        fullFormatDf = pickle.load(open('misc/tubeLayout-'+folderName+'-cell.pkl','rb'))
        dfList = []
        for fileName in os.listdir('inputData/bulkCSVFiles/'):
            if '.csv' in fileName:
                bulkTubeCSVFileName = fileName
                columnMultiIndexTuples = cellDataProcessing.parseCellCSVHeaders(pd.read_csv('inputData/bulkCSVFiles/'+bulkTubeCSVFileName).columns)
                columnMultiIndex = pd.MultiIndex.from_tuples(columnMultiIndexTuples,names=['CellType','Marker','Statistic'])

                fullData = pd.read_csv('inputData/bulkCSVFiles/'+bulkTubeCSVFileName,header=0)
                if 'Unnamed' in fullData.columns[-1]:
                    data = fullData.iloc[:-2,1:-1].values
                else:
                    data = fullData.iloc[:-2,1:].values
                sampleNames = fullData.iloc[:-2,0].values.ravel()
                sampleIndexStart = fullFormatDf.values.ravel().tolist().index(sampleNames[0])
                sampleIndexEnd = fullFormatDf.values.ravel().tolist().index(sampleNames[-1])
                rowMultiIndex = fullFormatDf.iloc[sampleIndexStart:sampleIndexEnd+1,:].index
                
                timeDataList = []
                timeSubsets = []
                times = experimentParameters['columnLevelValues']

                #Can use sample name file to assign time values
                if 'sampleNameFile.xlsx' in os.listdir('misc') or 'sampleNameFile.csv' in os.listdir('misc'):
                    if 'sampleNameFile.xlsx' in os.listdir('misc'): 
                        sampleNameDf = pd.read_excel('misc/sampleNameFile.xlsx')
                    else:
                        sampleNameDf = pd.read_csv('misc/sampleNameFile.csv')
                    if 'Time' in sampleNameDf.columns:
                        for time in times:
                            timeIndices = []
                            for row in range(sampleNameDf.shape[0]):
                                if sampleNameDf[experimentParameters['columnVariableName']].values[row] == time:
                                    timeIndices.append(row)
                            timeSubsets.append(timeIndices)
                    #Otherwise just assume 1 timepoint (HACK NEED TO FIX EVENTUALLY)
                    else:
                        timeSubsets.append(list(range(data.shape[0])))
                #Otherwise just assume 1 timepoint (HACK NEED TO FIX EVENTUALLY)
                else:
                    timeSubsets.append(list(range(data.shape[0])))

                for timeSubset in timeSubsets: 
                    dataList = []
                    columnTupleList = []
                    for i,columnTuple in enumerate(columnMultiIndexTuples):
                        ser = pd.Series(data[timeSubset,i],index=rowMultiIndex)
                        dataList.append(ser)
                        columnTupleList.append(tuple(columnTuple))
                    fullExperimentDf = pd.concat(dataList,keys=columnTupleList,names=['CellType','Marker','Statistic'])
                    timeDataList.append(fullExperimentDf)
                
                k = pd.concat(timeDataList,keys=times,names=[experimentParameters['columnVariableName']])
                repeatList = []
                for name in k.index:
                    if name not in repeatList:
                        repeatList.append(name)
                    else:
                        print('Repeated:')
                        print(name)
                partialExperimentDf = pd.concat(timeDataList,keys=times,names=[experimentParameters['columnVariableName']]).unstack(experimentParameters['columnVariableName'])
                dfList.append(partialExperimentDf)

        fullExperimentDf = pd.concat(dfList)
    else:
        if dataType == 'singlecell':
            realDataType = 'singlecell'
            dataType = 'cell'
            print(experimentParameters)
            if 'multiplexingOption' in experimentParameters:
                if experimentParameters['multiplexingOption'] != 'None':
                    demultiplexSingleCellData(experimentParameters)
        else:
            realDataType = dataType
        
        if 'multiplexingOption' in experimentParameters.keys():
            if experimentParameters['multiplexingOption'] != 'None':
                unpackMultiplexedPlates(experimentParameters,folderName,dataType)

        #Legacy experiment parameter files compatibility
        if 'paired' in experimentParameters.keys():
            if experimentParameters['paired']:
                numRowPlates = 2
            else:
                numRowPlates = 1
            numColumnPlates = int(experimentParameters['numPlates']/numRowPlates)
        else:
            numRowPlates = experimentParameters['numRowPlates'] 
            numColumnPlates = experimentParameters['numColumnPlates']

        #Combine plate and well IDs into a single ID val for every single sample in the experiment
        identificationMatrix = np.empty(layoutDict['plateID'].shape,dtype=object)
        for row in range(identificationMatrix.shape[0]):
            for col in range(identificationMatrix.shape[1]):
                wellID = layoutDict['wellID'][row,col]
                plateID = layoutDict['plateID'][row,col]
                fullID = plateID+'-'+wellID
                identificationMatrix[row,col] = fullID

        plateNames = np.unique(layoutDict['plateID'])

        plateDimensions = experimentParameters['overallPlateDimensions'] 
        levels = experimentParameters['allLevelNames']
        conditionLevels = experimentParameters['conditionLevelNames']
        conditionLevelValues = experimentParameters['conditionLevelValues']
        allLevelValues = experimentParameters['allLevelValues']

        sortedData,sortedFiles = cleanUpFlowjoCSV(plateNames,folderName,dataType,experimentParameters)
        allRawData,newLevelList = returnMultiIndex(sortedData,sortedFiles,realDataType,folderName)
        dfList = []
        for rawData,plateID in zip(allRawData,plateNames):
            fullTupleList = []
            index = 0
            for row in range(rawData.shape[0]):
                sampleID = int(rawData.iloc[row,0])-1
                wellID = plateRowLetters[int(sampleID/plateDimensions[1])]+str(plateColumnNumbers[sampleID % plateDimensions[1]])
                fullID = plateID+'-'+wellID
                sampleLocation = np.argwhere(identificationMatrix == fullID)[0]
                column = []
                tupleList = []
                for levelID in layoutDict['keys']:
                    levelValueID = layoutDict['keys'][levelID][sampleLocation[0],sampleLocation[1]]
                    if levelValueID != 'blank':
                        level = levels[levelID]
                        levelValue = allLevelValues[level][levelValueID]
                        tupleList.append(levelValue)
                    else:
                        print('wat2')
                if len(tupleList) != 0:
                    fullTupleList.append(tupleList)
                index+=1
            mi = pd.MultiIndex.from_tuples(fullTupleList,names=levels)
            columnSeriesList = []
            columnTupleList = []
            #single cell dataframe only has 2 columns
            if realDataType == 'singlecell':
                for columnTuple in newLevelList:
                    columnSeries = pd.Series(rawData.values[:,1],index=mi)
                    columnSeriesList.append(columnSeries)
                    columnTupleList.append(tuple(columnTuple))
            else:
                for column,columnTuple in enumerate(newLevelList):
                    columnSeries = pd.Series(rawData.values[:,column+1],index=mi)
                    columnSeriesList.append(columnSeries)
                    columnTupleList.append(tuple(columnTuple))
            platedf = pd.concat(columnSeriesList,axis=0,keys=columnTupleList,names=dataTypeLevelNames[realDataType])
            dfList.append(platedf)

        idx=pd.IndexSlice 
        fullExperimentDf = pd.concat(dfList)
        """
        fullExperimentDf = pd.Series(np.nan_to_num(fullExperimentDf.values),index=fullExperimentDf.index)
        dflist = []
        for i in pd.unique(fullExperimentDf.index.get_level_values('Time')):
            if i != 77.0:
                dflist.append(fullExperimentDf.loc[idx[:,:,:,i]]+i)
            else:
                dflist.append(fullExperimentDf.loc[idx[:,:,:,i]])
        fullExperimentDf = pd.concat(dflist,keys=list(pd.unique(fullExperimentDf.index.get_level_values('Time'))),names=['Time'],axis=1)
        print(fullExperimentDf)
        print(fullExperimentDf[fullExperimentDf.index.duplicated()])
        print('wat') 
        """
        #dfl = [fullExperimentDf.xs([12.0],level=['Time']),fullExperimentDf.xs([60.0],level=['Time']),fullExperimentDf.xs([96.0],level=['Time']),fullExperimentDf.xs([156.0],level=['Time'])]
        tempdf = fullExperimentDf.to_frame('temp')
        temp = []
        for row in range(fullExperimentDf.shape[0]):
            name = list(tempdf.iloc[row,:].name)
            if name in temp:
                print(name)
                print(row)
            else:
                temp.append(name)
        #Remove blanks
        for i,level in enumerate(fullExperimentDf.index.names):
            tempLevelValues = pd.unique(fullExperimentDf.index.get_level_values(level))
            if 'Blank' in tempLevelValues:
                fullExperimentDf = fullExperimentDf.drop('Blank',level=i)
        temp = []
        temp2 = []
        tempdf = fullExperimentDf.to_frame('wat')
        for row in range(fullExperimentDf.shape[0]):
            name = list(tempdf.iloc[row,:].name)
            if name not in temp:
                temp.append(name)
            else:
                temp2.append(name)
        fullExperimentDf = fullExperimentDf.unstack(experimentParameters['columnVariableName'])
    
    fullExperimentDf = reorderDfByInputOrder(experimentParameters,fullExperimentDf)
    return fullExperimentDf

def convertDataFramesToExcel(folderName,secondPath,dataType,df,useModifiedDf):
    if useModifiedDf:
        modifiedString = '-modified'
    else:
        modifiedString = ''
    writer = pd.ExcelWriter('outputData/excelFiles/excelFile-'+folderName+'-'+dataType+modifiedString+'.xlsx')
    if dataType == 'cyt':
        dfg = pickle.load(open('outputData/pickleFiles/cytokineGFIPickleFile-'+folderName+'.pkl','rb'))
        dfc = pickle.load(open('outputData/pickleFiles/'+dataTypeDataFrameFileNames[dataType]+'-'+folderName+modifiedString+'.pkl','rb'))
        dfg.to_excel(writer,'GFI')
        dfc.to_excel(writer,'Concentration')
    else:
        for statistic in list(pd.unique(df.index.get_level_values('Statistic'))):
            statisticDf = df.xs(statistic,level='Statistic')
            statisticDf.to_excel(writer,statistic)
    writer.save()
    print(dataType[0].upper()+dataType[1:]+' Excel file Saved')

def saveFinalDataFrames(folderName,secondPath,experimentNumber,dataType,fullExperimentDf,excel_data):
    with open('outputData/pickleFiles/'+dataTypeDataFrameFileNames[dataType]+'-'+folderName+'.pkl', "wb") as f:
        pickle.dump(fullExperimentDf, f)
    convertDataFramesToExcel(folderName,secondPath,dataType,fullExperimentDf,False)
    print(fullExperimentDf)
    
    modifiedFullExperimentDf = returnModifiedDf(experimentNumber,fullExperimentDf,dataType,excel_data)
    with open('outputData/pickleFiles/'+dataTypeDataFrameFileNames[dataType]+'-'+folderName+'-modified.pkl', "wb") as f:
        pickle.dump(modifiedFullExperimentDf, f)
    convertDataFramesToExcel(folderName,secondPath,dataType,modifiedFullExperimentDf,True)
    print(modifiedFullExperimentDf)
