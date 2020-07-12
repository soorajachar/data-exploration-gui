#!/usr/bin/env python3 
import math,pickle,os,sys,json,time,glob
from sys import platform as sys_pf
if sys_pf == 'darwin':
    import matplotlib
    matplotlib.use("TkAgg")
import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns
import pandas as pd
from sklearn.preprocessing import StandardScaler
from modifyDataFrames import returnModifiedDf
from miscFunctions import reindexDataFrame

def produceSingleCellHeaders(cellTypes):
    newMultiIndexList = []
    for cellType in cellTypes:
        newMultiIndexList.append([cellType])
    return newMultiIndexList

def createTubeSingleCellDataFrame(folderName,fileNameDf):
    fileNameDf = fileNameDf.stack()
    dflist = []
    #Multiple subfolders for each cell type
    k = os.listdir('inputData/singleCellCSVFiles/')
    dirList = [x[0].split('/')[2] for x in os.walk('inputData/singleCellCSVFiles/')]
    if len(dirList)-1 > 0:
        #Find markers common to all fcs files
        allMarkers = []
        for cellType in dirList[1:]:
            scalingTypeList = os.listdir('inputData/singleCellCSVFiles/'+cellType+'/')
            if '.csv' not in scalingTypeList[0]:
                scalingType = scalingTypeList[1].split('_')[0]
            else:
                scalingType = scalingTypeList[0].split('_')[0]
            for row in range(fileNameDf.shape[0]):
                fileName = fileNameDf.iloc[row].split('.')[0]
                possibleScale = fileName.split('_')
                if possibleScale[0] == scalingType:
                    tempScalingType = ''
                else:
                    tempScalingType = scalingType+'_'
                csv = pd.read_csv('inputData/singleCellCSVFiles/'+cellType+'/'+tempScalingType+fileName+'_'+cellType+'.csv')
                newcolumns = []
                for i,column in enumerate(csv.columns):
                    if 'FSC' in column:
                        newcolumns.append('Size')
                    elif 'SSC' in column:
                        newcolumns.append('Granularity')
                    else:
                        if ' :: ' in column:
                            newcolumns.append(column.split(' :: ')[1])
                        else:
                            newcolumns.append(column)
                allMarkers.append(newcolumns)
        commonMarkers = list(set.intersection(*map(set,allMarkers)))

        for cellType in dirList[1:]:
            scalingTypeList = os.listdir('inputData/singleCellCSVFiles/'+cellType+'/')
            if '.csv' not in scalingTypeList[0]:
                scalingType = scalingTypeList[1].split('_')[0]
            else:
                scalingType = scalingTypeList[0].split('_')[0]
            for row in range(fileNameDf.shape[0]):
                fileName = fileNameDf.iloc[row].split('.')[0]
                possibleScale = fileName.split('_')
                if possibleScale[0] == scalingType:
                    tempScalingType = ''
                else:
                    tempScalingType = scalingType+'_'
                csv = pd.read_csv('inputData/singleCellCSVFiles/'+cellType+'/'+tempScalingType+fileName+'_'+cellType+'.csv')
                newcolumns = []
                for i,column in enumerate(csv.columns):
                    if 'FSC' in column:
                        newcolumns.append('Size')
                    elif 'SSC' in column:
                        newcolumns.append('Granularity')
                    else:
                        if ' :: ' in column:
                            newcolumns.append(column.split(' :: ')[1])
                        else:
                            newcolumns.append(column)
                indexTuple = []
                for row2 in range(csv.shape[0]):
                    indexTuple.append([cellType]+list(fileNameDf.index.tolist()[row])+[row2+1])
                mi = pd.MultiIndex.from_tuples(indexTuple,names=['CellType']+list(fileNameDf.index.names)+['Event'])
                newdf = pd.DataFrame(csv.values,index=mi,columns=newcolumns)
                #Subset by common markers
                commonMarkerDf = newdf.loc[:,commonMarkers]
                dflist.append(commonMarkerDf)
    #Only one celltype/celltype doesn't matter
    else:
        #Find markers common to all fcs files
        allMarkers = []
        cellType = 'allCells'
        scalingTypeList = os.listdir('inputData/singleCellCSVFiles/')
        if '.csv' not in scalingTypeList[0]:
            scalingType = scalingTypeList[1].split('_')[0]
        else:
            scalingType = scalingTypeList[0].split('_')[0]
        for row in range(fileNameDf.shape[0]):
            fileName = fileNameDf.iloc[row].split('.')[0]
            possibleScale = fileName.split('_')
            if possibleScale[0] == scalingType:
                tempScalingType = ''
            else:
                tempScalingType = scalingType+'_'
            csv = pd.read_csv('inputData/singleCellCSVFiles/'+tempScalingType+fileName+'.csv')
            newcolumns = []
            for i,column in enumerate(csv.columns):
                if 'FSC' in column:
                    newcolumns.append('Size')
                elif 'SSC' in column:
                    newcolumns.append('Granularity')
                else:
                    if ' :: ' in column:
                        newcolumns.append(column.split(' :: ')[1])
                    else:
                        newcolumns.append(column)
            allMarkers.append(newcolumns)
        commonMarkers = list(set.intersection(*map(set,allMarkers)))

        for row in range(fileNameDf.shape[0]):
            fileName = fileNameDf.iloc[row].split('.')[0]
            possibleScale = fileName.split('_')
            if possibleScale[0] == scalingType:
                tempScalingType = ''
            else:
                tempScalingType = scalingType+'_'
            csv = pd.read_csv('inputData/singleCellCSVFiles/'+tempScalingType+fileName+'.csv')
            newcolumns = []
            for i,column in enumerate(csv.columns):
                if 'FSC' in column:
                    newcolumns.append('Size')
                elif 'SSC' in column:
                    newcolumns.append('Granularity')
                else:
                    if ' :: ' in column:
                        newcolumns.append(column.split(' :: ')[1])
                    else:
                        newcolumns.append(column)
            indexTuple = []
            for row2 in range(csv.shape[0]):
                indexTuple.append([cellType]+list(fileNameDf.index.tolist()[row])+[row2+1])
            mi = pd.MultiIndex.from_tuples(indexTuple,names=['CellType']+list(fileNameDf.index.names)+['Event'])
            newdf = pd.DataFrame(csv.values,index=mi,columns=newcolumns)
            commonMarkerDf = newdf.loc[:,commonMarkers]
            dflist.append(commonMarkerDf)
    
    completeDataFrame = pd.concat(dflist,axis=0)
    completeDataFrame.columns.name = 'Marker'
    print(completeDataFrame)
    completeDataFrame.to_hdf('outputData/pickleFiles/initialSingleCellDf-'+scalingType+'-'+folderName+'.h5', key='df', mode='w')
    #with open('outputData/pickleFiles/initialSingleCellDf-'+scalingType+'-'+folderName+'.pkl','wb') as f:
    #    pickle.dump(completeDataFrame,f)

def createInitialSingleCellDataFrame(folderName,experimentNumber,fileNameDataFrame):
    
    #Grabs a file from samples to read marker names off of
    for fileName in os.listdir('inputData/singleCellCSVFiles/A1/'):
        if 'DS' not in fileName:
            cellType = fileName
    tempFilePath = 'inputData/singleCellCSVFiles/A1/'+cellType+'/' 
    fileExtension = '.csv'
    tempFileName = glob.glob(tempFilePath+'*'+fileExtension)[0]
    experimentalChannelDf = pd.read_csv(tempFileName, header=0)
    experimentalChannelNames = experimentalChannelDf.columns.tolist()
    experimentalMarkerNames = []
    if 'gateVals' in os.listdir('misc'):
        gatingMarkers = pickle.load(open('inputFiles/gateVals.pkl','rb'))
    else: 
        gatingMarkers = {}
    if 'TCell_Gate' not in gatingMarkers.keys():
        gatingMarkers['TCell_Gate'] = 'none'
    if 'APC_Gate' not in gatingMarkers.keys():
        gatingMarkers['APC_Gate'] = 'none'
    #Creates column headings for all measured parameters
    for i in range(len(experimentalChannelNames)):
        #"Comp-APC-A :: CTFR"
        if('::' in experimentalChannelNames[i]):
            experimentalMarkerName = experimentalChannelNames[i].split(' :: ')[1]
        else:
            experimentalMarkerName = experimentalChannelNames[i]
        if len(gatingMarkers) > 0:
            if gatingMarkers['TCell_Gate'] == experimentalMarkerName:
                experimentalMarkerNames.append('TCell_Gate')
            elif gatingMarkers['APC_Gate'] == experimentalMarkerName:
                experimentalMarkerNames.append('APC_Gate')
            else:
                if 'FSC-A' in experimentalChannelNames[i]:
                    experimentalMarkerNames.append('Size')
                elif 'SSC-A' in experimentalChannelNames[i]:
                    experimentalMarkerNames.append('Granularity')
                else:
                    experimentalMarkerNames.append(experimentalMarkerName)
        else:
            if 'FSC-A' in experimentalChannelNames[i]:
                experimentalMarkerNames.append('Size')
            elif 'SSC-A' in experimentalChannelNames[i]:
                experimentalMarkerNames.append('Granularity')
            else:
                experimentalMarkerNames.append(experimentalMarkerName)
    stackedFileFrame = fileNameDataFrame.stack()
    levelNames = list(stackedFileFrame.index.names)
    singleCellLevelNames = levelNames+['Event']
    channelBool = False
    scaleBool = False
    scalingLevels = [] 
    for fileName in glob.glob(tempFilePath+'*'+fileExtension):
        if 'channel' in fileName:
            channelBool = True
        if 'scale' in fileName:
            scaleBool = True
    if not channelBool and not scaleBool:
        print('No csv files with channel or scale value suffixes found. Please try again.')
    else:
        if channelBool:
            scalingLevels.append('channel')
        if scaleBool:
            scalingLevels.append('scale')
        for scalingType in scalingLevels:
            fullFileFrameTemp = stackedFileFrame.copy().to_frame('Temp')
            completeDataFrameList = []
            for row in range(stackedFileFrame.shape[0]):
                levelValues = fullFileFrameTemp.iloc[row].name
                cellType = levelValues[0]
                fileIndex = stackedFileFrame.iloc[row].rfind('/')
                beforeCellType = stackedFileFrame.iloc[row][:fileIndex]
                afterCellType = scalingType+'_'+stackedFileFrame.iloc[row][fileIndex+1:]+'_'+cellType+fileExtension
                fullFileName = beforeCellType+'/'+cellType+'/'+afterCellType
                
                fcsDf = pd.read_csv(fullFileName,header=0)
                eventNumber = fcsDf.shape[0]
                eventList = range(1,eventNumber+1)
                allLevelValues = []
                for event in eventList:
                    allLevelValues.append(list(levelValues)+[event])
                newMultiIndex = pd.MultiIndex.from_tuples(allLevelValues,names=singleCellLevelNames)
                """
                temp = []
                if 'CD45R' in fcsDf.columns or (cols[2] == 'CD25' and cols[3] == 'CD4'):
                    tempDf1 = fcsDf.copy().iloc[:,2]
                    tempDf2 = fcsDf.copy().iloc[:,3]
                    fcsDf.iloc[:,2] = tempDf2
                    fcsDf.iloc[:,3] = tempDf1
                cols = list(fcsDf.columns)
                for i,col in enumerate(fcsDf.columns):
                    if 'CD45R' not in col:
                        temp.append(i)
                fcsDf = fcsDf.iloc[:,temp]
                """
                newDf = pd.DataFrame(fcsDf.values,index=newMultiIndex,columns=experimentalMarkerNames)
                #newDf = pd.DataFrame(fcsDf.values,index=newMultiIndex,columns=fcsDf.columns)
                completeDataFrameList.append(newDf)
            
            completeDataFrame = pd.concat(completeDataFrameList)
            completeDataFrame.columns.name = 'Marker'
            
            #Remove extraneous markers (namely -h parameters)
            columnsToKeep = []
            for col,column in enumerate(completeDataFrame.columns):
                if '-H' not in column and 'Time' not in column:
                    columnsToKeep.append(col)
            completeDataFrame = completeDataFrame.iloc[:,columnsToKeep]

            completeDataFrame.to_hdf('outputData/pickleFiles/initialSingleCellDf-'+scalingType+'-'+folderName+'.h5', key='df', mode='w')
            #with open('outputData/pickleFiles/initialSingleCellDf-'+scalingType+'-'+folderName+'.pkl','wb') as f:
            #    pickle.dump(completeDataFrame,f)
            #Legacy proliferation code
            if 'TCell_Gate' in completeDataFrame.columns:
                tcellGateName = 'TCell_Gate'
            else:
                if 'CTV' in completeDataFrame.columns:
                    for column in completeDataFrame.columns:
                        if 'CTV' in column:
                            tcellGateName = column
                            break
                else:
                    tcellGateName = 'None'
            
            if scalingType == 'channel':
                if 'TCells' in pd.unique(completeDataFrame.index.get_level_values('CellType')) and tcellGateName != 'None':
                    logicleDataProliferation = completeDataFrame[tcellGateName].xs(['TCells'],level=['CellType'])
                    with open('misc/logicleProliferationDf.pkl','wb') as f:
                        pickle.dump(logicleDataProliferation,f)
            else:
                if 'TCells' in pd.unique(completeDataFrame.index.get_level_values('CellType')) and tcellGateName != 'None':
                    rawDataProliferation = completeDataFrame[tcellGateName].xs(['TCells'],level=['CellType'])
                    with open('misc/rawProliferationDf.pkl','wb') as f:
                        pickle.dump(rawDataProliferation,f)

def createCompleteSingleCellDf(folderName):
   
    idx = pd.IndexSlice
    
    initialSingleCellDf = pickle.load(open('outputData/pickleFiles/initialSingleCellDf-channel-'+folderName+'.pkl','rb'))
    initialSingleCellDf = initialSingleCellDf.drop(['APC_Gate','TCell_Gate'],axis=1)
    bulkCytokineConcentrationDf = pickle.load(open('semiProcessedData/cytokineConcentrationPickleFile-'+folderName+'-modified.pkl','rb'))
    proliferationDf = pickle.load(open('semiProcessedData/singleCellDataFrame-proliferation-'+folderName+'.pkl','rb'))

    bulkCellStatisticDf = pickle.load(open('semiProcessedData/cellStatisticPickleFile-'+folderName+'-modified.pkl','rb'))

    cytokineLevelNamesUM = []
    cytokineLevelNames = []
    cellLevelNames = []

    tempCellLevels = list(bulkCellStatisticDf.iloc[0,:].name)[:3]
    tempCellDf = bulkCellStatisticDf.loc[tempCellLevels[0]].loc[tempCellLevels[1]].loc[tempCellLevels[2]]

    tempCytokineLevels = list(bulkCytokineConcentrationDf.iloc[0,:].name)[:1]
    tempCytokineDf = bulkCytokineConcentrationDf.loc[tempCytokineLevels[0]]
    
    unmodifiedCytokineConcentrationDf = pickle.load(open('semiProcessedData/cytokineConcentrationPickleFile-'+folderName+'.pkl','rb')).loc[tempCytokineLevels[0]]
    
    for row in range(tempCellDf.shape[0]):
        levelNames = list(tempCellDf.iloc[row,:].name)
        for column in tempCellDf.columns:
            cellLevelNames.append(tuple(levelNames+[column]))
    for row in range(tempCytokineDf.shape[0]):
        levelNames = list(tempCytokineDf.iloc[row,:].name)
        for column in tempCytokineDf.columns:
            cytokineLevelNames.append(tuple(levelNames+[column]))
    for row in range(unmodifiedCytokineConcentrationDf.shape[0]):
        levelNames = list(unmodifiedCytokineConcentrationDf.iloc[row,:].name)
        for column in unmodifiedCytokineConcentrationDf.columns:
            cytokineLevelNamesUM.append(tuple(levelNames+[column]))
    
    differences = list((set(tuple(cytokineLevelNamesUM)) | set(tuple(cytokineLevelNames)) | set(tuple(cellLevelNames))) - (set(tuple(cytokineLevelNamesUM)) & set(tuple(cytokineLevelNames)) & set(tuple(cellLevelNames))))
    levelsToKeep = [0]*initialSingleCellDf.shape[0]
    k=0
    row = 0
    while row < initialSingleCellDf.shape[0]:
        levelNames = tuple(initialSingleCellDf.iloc[row,:].name)[:-1]
        stackedLevelNames = tuple(list(levelNames)+[slice(None)])
        stackedLength = initialSingleCellDf.loc[stackedLevelNames,:].shape[0]
        #Check logic here later; should do for now
        if (tuple(levelNames) in differences):
            print(levelNames)
            pass
        else:
            rowVals = range(row,row+stackedLength)
            levelsToKeep[k:k+stackedLength] = rowVals
            k+=stackedLength
        row+=stackedLength
    initialSingleCellDf = initialSingleCellDf.iloc[levelsToKeep[:k],:]
    proliferationDf = proliferationDf.iloc[levelsToKeep[:k],:]
    indexList = []
    numEventsList = []
    for elem in pd.unique(initialSingleCellDf.index):
        indexList.append(elem[:-1])
    indexList = pd.unique(indexList)
    for index in indexList:
        numEventsList.append(initialSingleCellDf.loc[idx[index],:].shape[0])
    completeSingleCellCytokineValues = []
    for cytokine in pd.unique(bulkCytokineConcentrationDf.index.get_level_values(0)):
        individualCytokineSingleCellValues = []
        for index,numEvents in zip(indexList,numEventsList):
            bulkIndex = tuple([cytokine]+list(index[:-1]))
            bulkCytokineValue = bulkCytokineConcentrationDf.loc[idx[bulkIndex],index[-1]]
            singleCellCytokineValues = np.repeat(bulkCytokineValue,numEvents)
            individualCytokineSingleCellValues.append(singleCellCytokineValues)
        completeSingleCellCytokineValues.append(np.concatenate(individualCytokineSingleCellValues))
    singleCellCytokineMatrix = np.stack(completeSingleCellCytokineValues,axis=1)
    singleCellCytokineDf = pd.DataFrame(singleCellCytokineMatrix,index=initialSingleCellDf.index,columns=pd.unique(bulkCytokineConcentrationDf.index.get_level_values(0)))
     
    completeSingleCellDf = pd.concat([initialSingleCellDf,singleCellCytokineDf,proliferationDf],keys=['Markers','Cytokines','Proliferation'],names=['DataType','Parameter'],axis=1)
    print(completeSingleCellDf)
    with open('outputData/pickleFiles/singleCellDataFrame-complete-'+folderName+'.pkl','wb') as f:
        pickle.dump(completeSingleCellDf,f)
