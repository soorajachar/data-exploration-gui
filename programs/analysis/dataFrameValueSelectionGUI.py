#!/usr/bin/env python3 
import json,pickle,math,matplotlib,sys,os,string,subprocess
from sys import platform as sys_pf
if sys_pf == 'darwin':
    import matplotlib
    matplotlib.use("TkAgg")
import numpy as np
import pandas as pd
import seaborn as sns
import tkinter as tk
import itertools
from matplotlib import pyplot as plt
from scipy import stats
from sklearn.manifold import TSNE
from sklearn.preprocessing import MinMaxScaler
sys.path.insert(0, '../dataprocessing/')
from miscFunctions import sortSINumerically,reindexDataFrame,setMaxWidth
sys.path.insert(0, '../plotting/')
from plottingGUI import checkUncheckAllButton,selectLevelValuesPage
import facetPlotLibrary as fpl
from umap import UMAP
from sklearn.manifold import Isomap
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans, DBSCAN, AgglomerativeClustering
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
#from reduceDimensions import reduceDimensions 
import operateOnDataSelection as ods

idx = pd.IndexSlice
dataTypeObservableRangeDict = {'cyt':1,'cell':3,'prolif':1,'singlecell':1}
realDataTypeNameDict = {'cyt':'Supernatant','cell':'Surface/Intracellular Marker','prolif':'Proliferation','singlecell':'Single Cell'}

#Get level names and values into an easily accessible dictionary
def createLabelDict(df,levelRange):
    fulldf = df.stack()
    labelDict = {}
    for i in range(levelRange[0],levelRange[1]):
        levelName = fulldf.index.levels[i].name
        if levelName not in ['Event','event']:
            labelDict[levelName] = list(pd.unique(fulldf.index.get_level_values(levelName)))
    return labelDict

def sampleDataFrame(df,sampleType,sampleSubset,fraction='',nmax=''):
    #whole dataframe
    if sampleSubset == 'all':
        if sampleType == 'fraction':
            if float(fraction) == 1.0:
                sampledDf = df.copy()
            else:
                sampledDf = df.sample(frac=float(fraction))
        else:
            sampledDf = df.sample(n=int(nmax))
    #per condition
    else:
        grouped = df.groupby(list(df.index.names)[:-1])
        if sampleType == 'fraction':
            if float(fraction) == 1.0:
                sampledDf = df.copy()
            else:
                sampledDf = grouped.apply(lambda x: x.sample(frac=float(fraction)))
                sampledDf = sampledDf.droplevel(list(range(len(df.index.names)-1)),axis=0)
        else:
            sampledDf = grouped.apply(lambda x: x.sample(int(nmax)) if len(x) > int(nmax) else x)
            sampledDf = sampledDf.droplevel(list(range(len(df.index.names)-1)),axis=0)
    return sampledDf

#(DataSelectionHomePage,AnalysisStartPage,folderName,AnalysisStartPage,'cluster',mainhomepage)
class DataSelectionHomePage(tk.Frame):
    def __init__(self, master,fsp,fName,bp,pt,shp):
        
        global finalSwitchPage,folderName,backpage,secondaryhomepage
        finalSwitchPage = fsp
        folderName = fName
        backpage = bp
        processType = pt
        secondaryhomepage = shp

        dataTypeFileDict = {'cyt':'cytokineConcentrationPickleFile-'+folderName+'-modified.pkl','cell':'cellStatisticPickleFile-'+folderName+'-modified.pkl','prolif':'proliferationStatisticPickleFile-'+folderName+'-modified.pkl','singlecell':'initialSingleCellDf-channel-'+folderName+'.pkl'}
        tk.Frame.__init__(self, master)
        
        mainWindow = tk.Frame(self)
        mainWindow.pack(side=tk.TOP,padx=10,pady=10)
        
        l2 = tk.Label(mainWindow, text="DataTypes to postprocess:").grid(row=0,column=1)
        dataTypeCheckButtons = []
        checkButtonVariableList = []
        for row,dataType in enumerate(dataTypeFileDict):
            dataTypeBool = tk.BooleanVar()
            cb = tk.Checkbutton(mainWindow, text=realDataTypeNameDict[dataType],padx = 20, variable=dataTypeBool)
            cb.grid(row=row+1,column=1,sticky=tk.W)
            if dataTypeFileDict[dataType] not in os.listdir('outputData/pickleFiles'):
                print(os.listdir('outputData/pickleFiles'))
                print(dataTypeFileDict[dataType].split('.'))
                if dataTypeFileDict[dataType].split('.')[0]+'.h5' not in os.listdir('outputData/pickleFiles'):
                    cb.config(state=tk.DISABLED)
            else:
                if dataType != 'singlecell':
                    cb.select()
            dataTypeCheckButtons.append(cb)
            checkButtonVariableList.append(dataTypeBool)
        
        tk.Label(mainWindow,text='Downsample by:').grid(row=len(dataTypeFileDict.keys()),column=2)
        sampleMethodVar = tk.StringVar(value='fraction')
        sampleRbT1 = tk.Radiobutton(mainWindow,text='fraction=',value='fraction',variable=sampleMethodVar)
        sampleRbT1.grid(row=len(dataTypeFileDict.keys()),column=3,sticky=tk.W)
        sampleRbT2 = tk.Radiobutton(mainWindow,text='count=',value='count',variable=sampleMethodVar)
        sampleRbT2.grid(row=len(dataTypeFileDict.keys())+1,column=3,sticky=tk.W)
        fractionEntry = tk.Entry(mainWindow)
        fractionEntry.grid(row=len(dataTypeFileDict.keys()),column=4,sticky=tk.W)
        fractionEntry.insert(tk.END, '0.1')
        countEntry = tk.Entry(mainWindow)
        countEntry.grid(row=len(dataTypeFileDict.keys())+1,column=4,sticky=tk.W)
        countEntry.insert(tk.END, '1000')
        
        tk.Label(mainWindow,text='across:').grid(row=len(dataTypeFileDict.keys()),column=5)
        sampleTypeVar = tk.StringVar(value='all')
        sampleRb1 = tk.Radiobutton(mainWindow,text='all',value='all',variable=sampleTypeVar)
        sampleRb1.grid(row=len(dataTypeFileDict.keys()),column=6,sticky=tk.W)
        sampleRb2 = tk.Radiobutton(mainWindow,text='perSample',value='perSample',variable=sampleTypeVar)
        sampleRb2.grid(row=len(dataTypeFileDict.keys())+1,column=6,sticky=tk.W)

        def collectInputs():
            global dataTypeDfDict
            dataTypeDfDict = {}
            if checkButtonVariableList[-1].get():
                h5FileBool = False
                for fileName in os.listdir('outputData/pickleFiles'):
                    if '.h5' in fileName:
                        h5FileBool = True
                        break
                if h5FileBool:
                    ogname = 'outputData/pickleFiles/'+dataTypeFileDict['singlecell']
                    newname = ogname.split('.')[0]+'.h5'
                    newdf = pd.read_hdf(newname, 'df')
                    dataTypeDfDict['singlecell'] = newdf
                else:
                    dataTypeDfDict['singlecell'] = pickle.load(open('outputData/pickleFiles/'+dataTypeFileDict['singlecell'],'rb'))
                dataTypeDfDict['singlecell'] = sampleDataFrame(dataTypeDfDict['singlecell'],sampleMethodVar.get(),sampleTypeVar.get(),fraction=fractionEntry.get(),nmax=countEntry.get())
            else:
                for checkButtonVariable,dataType in zip(checkButtonVariableList,dataTypeFileDict):
                    if checkButtonVariable.get():
                        dataTypeDfDict[dataType] = pickle.load(open('outputData/pickleFiles/'+dataTypeFileDict[dataType],'rb'))
                #print(dataTypeDfDict['cell'].loc[idx[:,:,'MFI'],:])
            
            global trueLabelDict
            sampleDataType = next(iter(dataTypeDfDict))
            sampleDf = dataTypeDfDict[sampleDataType]
            #dataTypeDf = [dataTypeObservableRangeDict[dataType],]
            trueLabelDict = createLabelDict(sampleDf,[dataTypeObservableRangeDict[sampleDataType],sampleDf.index.nlevels+1])
            #switchPage,trueLabelDict,backPage,fsp,fName,shp,pt
            print('wat3')
            print(secondaryhomepage)
            with open('misc/switchingParameters.pkl','wb') as f:
                pickle.dump([finalSwitchPage,folderName,backpage,backpage,processType,secondaryhomepage,[sampleMethodVar.get(),sampleTypeVar.get(),fractionEntry.get(),countEntry.get()]],f)
            master.switch_frame(selectLevelValuesPage,SelectDimensionsPage,trueLabelDict,DataSelectionHomePage,backpage,finalSwitchPage,folderName,secondaryhomepage,processType)
            #master.switch_frame(SelectDimensionsPage,folderName,dataTypeDfDict)

        buttonWindow = tk.Frame(self)
        buttonWindow.pack(side=tk.TOP,pady=10)

        tk.Button(buttonWindow, text="OK",command=lambda: collectInputs()).grid(row=5,column=0)
        #if backpage.num_args == 1:
        #    tk.Button(buttonWindow, text="Back",command=lambda: master.switch_frame(backpage,folderName,secondaryhomepage)).grid(row=5,column=1)
        #else:
        #    tk.Button(buttonWindow, text="Back",command=lambda: master.switch_frame(backpage,folderName,secondaryhomepage)).grid(row=5,column=1)
        print('wat')
        print(folderName)
        print(backpage)
        print(secondaryhomepage)
        tk.Button(buttonWindow, text="Back",command=lambda: master.switch_frame(backpage,folderName,secondaryhomepage)).grid(row=5,column=1)
        tk.Button(buttonWindow, text="Quit",command=quit).grid(row=5,column=2)

class SelectDimensionsPage(tk.Frame):
    def __init__(self, master,temp):
        tk.Frame.__init__(self, master)
        switchingParameters = pickle.load(open('misc/switchingParameters.pkl','rb'))
        if len(switchingParameters) == 6:
            finalSwitchPage,folderName,backpage,secondaryBackPage,processType,secondaryhomepage = switchingParameters
        else:
            finalSwitchPage,folderName,backpage,secondaryBackPage,processType,secondaryhomepage,downsampleFraction = switchingParameters
        print('wat2')
        print(folderName)
        print(backpage)
        print(secondaryhomepage)
        global includeLevelValues2
        includeLevelValues2 = temp
        labelWindow = tk.Frame(self)
        labelWindow.pack(side=tk.TOP,padx=10,fill=tk.X,expand=True)
        l1 = tk.Label(labelWindow, text='Which specific measurables do you want to include in the dimension reduction?',pady=10).pack()
        dataTypeLevelCheckButtonList = []
        dataTypeLevelCheckButtonVariableList = []
        for dataType in dataTypeDfDict:
            dataTypeDf = dataTypeDfDict[dataType]
            observableLevelDict = createLabelDict(dataTypeDf,[0,dataTypeObservableRangeDict[dataType]])
            levelValueCheckButtonList = []
            overallCheckButtonVariableList = []
            checkAllButtonList = []
            uncheckAllButtonList = []
            dataTypeWindow = tk.Frame(labelWindow,borderwidth=1,relief='groove')
            dataTypeLabel = tk.Label(dataTypeWindow,text=realDataTypeNameDict[dataType]+':',font="-weight bold").grid(row=0,column=0,columnspan=len(observableLevelDict)*6)
            i=0
            maxNumLevelValues = 0
            for levelName in observableLevelDict:
                j=0
                levelCheckButtonList = []
                levelCheckButtonVariableList = []
                levelLabel = tk.Label(dataTypeWindow, text=levelName+':')
                levelLabel.grid(row=1,column = i*6,sticky=tk.N,columnspan=5)
                for levelValue in observableLevelDict[levelName]:
                    includeLevelValueBool = tk.BooleanVar()
                    cb = tk.Checkbutton(dataTypeWindow, text=levelValue, variable=includeLevelValueBool)
                    cb.grid(row=j+4,column=i*6+2,columnspan=2,sticky=tk.W)
                    dataTypeWindow.grid_columnconfigure(i*6+3,weight=1)
                    cb.select()
                    levelCheckButtonList.append(cb)
                    levelCheckButtonVariableList.append(includeLevelValueBool)
                    j+=1
                
                checkAllButton1 = checkUncheckAllButton(dataTypeWindow,levelCheckButtonList, text='Check All')
                checkAllButton1.configure(command=checkAllButton1.checkAll)
                checkAllButton1.grid(row=2,column=i*6,sticky=tk.N,columnspan=3)
                checkAllButtonList.append(checkAllButton1)
                
                uncheckAllButton1 = checkUncheckAllButton(dataTypeWindow,levelCheckButtonList, text='Uncheck All')
                uncheckAllButton1.configure(command=checkAllButton1.uncheckAll)
                uncheckAllButton1.grid(row=2,column=i*6+3,sticky=tk.N,columnspan=3)
                uncheckAllButtonList.append(checkAllButton1)

                levelValueCheckButtonList.append(levelCheckButtonList)
                overallCheckButtonVariableList.append(levelCheckButtonVariableList)
                if len(observableLevelDict[levelName]) > maxNumLevelValues:
                    maxNumLevelValues = len(observableLevelDict[levelName])
                i+=1
            dataTypeLevelCheckButtonList.append(levelValueCheckButtonList)
            dataTypeLevelCheckButtonVariableList.append(overallCheckButtonVariableList)
            dataTypeWindow.pack(side=tk.LEFT,fill=tk.Y)
        
        def collectInputs():
            global dimensionDict
            dimensionDict = {}
            for obl,dataType in zip(dataTypeLevelCheckButtonVariableList,dataTypeDfDict):
                includeLevelValueList = []
                for checkButtonVariableList in obl:
                    tempLevelValueList = []
                    for checkButtonVariable in checkButtonVariableList:
                        tempLevelValueList.append(checkButtonVariable.get())
                    includeLevelValueList.append(tempLevelValueList)
                dimensionDict[dataType] = includeLevelValueList
            print(dimensionDict)
            master.switch_frame(FinalDimensionSelectionPage)

        buttonWindow = tk.Frame(self)
        buttonWindow.pack(side=tk.BOTTOM,pady=10)
        tk.Button(buttonWindow, text="OK",command=lambda: collectInputs()).grid(row=maxNumLevelValues+4,column=0)
        tk.Button(buttonWindow, text="Back",command=lambda: master.switch_frame(selectLevelValuesPage,SelectDimensionsPage,trueLabelDict,DataSelectionHomePage,secondaryBackPage,finalSwitchPage,folderName,secondaryhomepage,processType)).grid(row=maxNumLevelValues+4,column=1)
        tk.Button(buttonWindow, text="Quit",command=lambda: quit()).grid(row=maxNumLevelValues+4,column=2)

class FinalDimensionSelectionPage(tk.Frame):
    def __init__(self, master):
        tk.Frame.__init__(self, master)
        switchingParameters = pickle.load(open('misc/switchingParameters.pkl','rb'))
        if len(switchingParameters) == 6:
            finalSwitchPage,folderName,backpage,secondaryBackPage,processType,secondaryhomepage = switchingParameters
        else:
            finalSwitchPage,folderName,backpage,secondaryBackPage,processType,secondaryhomepage,downsampleFraction = switchingParameters

        mainWindow = tk.Frame(self)
        mainWindow.pack(side=tk.TOP,padx=10,fill=tk.X,expand=True)
         
        everyLevelCombo = []
        for dataType in dataTypeDfDict:
            rowList = []
            if dataType == 'cell':
                #Iterate through each row in datatype df; grab dimension names
                for row in range(dataTypeDfDict[dataType].shape[0]):
                    names = dataTypeDfDict[dataType].iloc[row,:].name
                    dimensionNames = names[:dataTypeObservableRangeDict[dataType]]
                    rowList.append(dimensionNames)
            else:
                for dimension in pd.unique(dataTypeDfDict[dataType].index.get_level_values(0)):
                    rowList.append(tuple([dimension]))
            print(rowList)
            #Go thorugh each level that was selected, add level values of each level
            selectedLevelList = []
            for i,level in enumerate(dataTypeDfDict[dataType].index.names[:len(dimensionDict[dataType])]):
                levelList = []
                for levelValue,includeLevelValue in zip(list(pd.unique(dataTypeDfDict[dataType].index.get_level_values(level))),dimensionDict[dataType][i]):
                    if includeLevelValue:
                        levelList.append(levelValue)
                selectedLevelList.append(levelList)
        
            #Get all possible combinations of level values from each dimension
            allPossibleSelectedLevelCombinations = itertools.product(*selectedLevelList)
            print(allPossibleSelectedLevelCombinations)
            for levelCombination in allPossibleSelectedLevelCombinations:
                print(levelCombination)
                if levelCombination in rowList:
                    everyLevelCombo.append(levelCombination)

        checkButtonWindow = tk.Frame(mainWindow)
        checkButtonWindow.grid(row=0,column=0)
        cblist = []
        cbvarlist = []
        print('wat')
        print(everyLevelCombo)
        tk.Label(checkButtonWindow,text='Select final dimensions to include in data subset:').grid(row=0,column=0)
        for i,levelcombo in enumerate(everyLevelCombo):
            includeLevelValueBool = tk.BooleanVar()
            cb = tk.Checkbutton(checkButtonWindow, text=','.join(levelcombo), variable=includeLevelValueBool)
            cb.grid(row=i+1,column=0,sticky=tk.W)
            cb.select()
            cblist.append(cb)
            cbvarlist.append(includeLevelValueBool)

        def collectInputs():
            #Create data subset data frame, with features in columns and sample names in rows
            dataSelectionTitle = titleEntry.get()
            if dataSelectionTitle == '':
                dataSelectionTitle = 'defaultTitle'
            dataSelectionTitle = dataSelectionTitle.replace('-','_')
            finalDimensions = []
            for cbvar,levelcombo in zip(cbvarlist,everyLevelCombo):
                if cbvar.get():
                    finalDimensions.append(levelcombo)
            print(finalDimensions)
            dataSelectionDf = createSubsettedDataFrame(dimensionDict,dataTypeDfDict,finalDimensions)
            with open('outputData/analysisFiles/subsettedData/'+dataSelectionTitle+'.pkl','wb') as f:
                pickle.dump(dataSelectionDf,f)
            #Use data subset to create preprocessed dataframes for each scaling option selected
            for checkbutton in cblist:
                if 'singlecell' in dataTypeDfDict.keys():
                    df = ods.operateOnData(dataSelectionDf,dataSelectionTitle,'scale',{'scalingMethod':checkbutton.var.get(),'downsampleFraction':downsampleFraction})
                else:
                    df = ods.operateOnData(dataSelectionDf,dataSelectionTitle,'scale',{'scalingMethod':checkbutton.var.get()})
            master.switch_frame(finalSwitchPage,folderName,secondaryhomepage)

        titleWindow = tk.Frame(self)
        titleWindow.pack(side=tk.TOP,pady=10)
        titleprefix = 'Data Selection'
        tk.Label(titleWindow,text=titleprefix+' Title: ').grid(row=0,column=0)
        titleEntry = tk.Entry(titleWindow)
        titleEntry.grid(row=0,column=1)
        
        scalingWindow = tk.Frame(self)
        scalingWindow.pack(side=tk.TOP,pady=10)
        tk.Label(scalingWindow,text='Data Scaling Options').grid(row=0,column=0,columnspan=3)
        cblist = []
        for i,scalingOption in enumerate(['minmax','maxabs','standard']):
            scalingVar = tk.StringVar()
            cb = tk.Checkbutton(scalingWindow,text=scalingOption,variable=scalingVar,onvalue=scalingOption,offvalue='none')
            cb.grid(row=1,column=i)
            cb.var = scalingVar
            if scalingOption == 'minmax':
                cb.select()
            else:
                cb.deselect()
            cblist.append(cb)
        
        buttonWindow = tk.Frame(self)
        buttonWindow.pack(side=tk.BOTTOM,pady=10)
        tk.Button(buttonWindow, text="OK",command=lambda: collectInputs()).grid(row=4,column=0)
        tk.Button(buttonWindow, text="Back",command=lambda: master.switch_frame(selectLevelValuesPage,SelectDimensionsPage,trueLabelDict,DataSelectionHomePage,secondaryBackPage,finalSwitchPage,folderName,secondaryhomepage,processType)).grid(row=4,column=1)
        tk.Button(buttonWindow, text="Quit",command=lambda: quit()).grid(row=4,column=2)

def createSubsettedDataFrame(dimensionDict,dataTypeDfDict,finalDimensions):
    dimensionReductionMatrixList = []
    
    #Subset index (samples)
    for dataType in dataTypeDfDict:
        tempdf = dataTypeDfDict[dataType].stack().to_frame('temp')
        k = 0
        subsettingList = [] 
        for i in range(len(tempdf.index.names)):
            if i < dataTypeObservableRangeDict[dataType] or tempdf.index.names[i] in ['Event','event']:
                subsettingList.append(slice(None))
            else:
                levelName = tempdf.index.names[i]
                levelValues = list(pd.unique(tempdf.index.get_level_values(levelName)))
                levelValueBooleanList = includeLevelValues2[k]
                levelValueList = []
                print(levelName)
                for j in range(len(levelValues)):
                    if levelValueBooleanList[j]:
                        levelValueList.append(levelValues[j])
                subsettingList.append(levelValueList)
                k+=1
        
        print('watwat')
        newdf = tempdf.loc[tuple(subsettingList),:]
        #Comment here
        #Undoes the stacking performed at begining of loop
        if dataType != 'singlecell':
            stackingVariable = 'Time'
        else:
            stackingVariable = 'Marker'
        newdf2 = newdf.unstack(stackingVariable).droplevel(axis=1,level=0)
        print(newdf2)
        wat = newdf.groupby(newdf.index.names[:-1],sort=False).first()
        watdf = pd.DataFrame(newdf2.values,index=wat.index,columns=newdf2.columns)
        watdf.columns.name = stackingVariable
        if dataType != 'singlecell':
            newdf3 = reindexDataFrame(newdf2,watdf,False,sortDataTypeLevels=False)
        else:
            newdf3 = newdf2.copy()
        print(newdf3)
        dataTypeDfDict[dataType] = newdf3
        #dataTypeDfDict[dataType] = newdf
    
    if dataType != 'singlecell':
        #Subset columns (measurables)
        postProcessingMatrices = []
        postProcessingFeatures = []
        allIndexList = []
        allIndexDictList = []
        #Iterate through datatypes
        for dataType in dataTypeDfDict:
            rowList = []
            #Iterate through each row in datatype df; grab dimension names
            for row in range(dataTypeDfDict[dataType].shape[0]):
                names = dataTypeDfDict[dataType].iloc[row,:].name
                dimensionNames = names[:dataTypeObservableRangeDict[dataType]]
                rowList.append(dimensionNames)
            #Go thorugh each level that was selected, add level values of each level
            selectedLevelList = []
            for i,level in enumerate(dataTypeDfDict[dataType].index.names[:len(dimensionDict[dataType])]):
                levelList = []
                for levelValue,includeLevelValue in zip(list(pd.unique(dataTypeDfDict[dataType].index.get_level_values(level))),dimensionDict[dataType][i]):
                    if includeLevelValue:
                        levelList.append(levelValue)
                selectedLevelList.append(levelList)
            
            #Get all possible combinations of level values from each dimension
            print(selectedLevelList)
            allPossibleSelectedLevelCombinations = itertools.product(*selectedLevelList)
            print(allPossibleSelectedLevelCombinations)
            rowindexlist = []
            #From original dataframe; select all rows that appear in the all possible combination list
            selectedDimensions = []
            for levelCombination in allPossibleSelectedLevelCombinations:
                if levelCombination in rowList and levelCombination in finalDimensions:
                    indices = [i for i, x in enumerate(rowList) if x == levelCombination]
                    rowindexlist+=indices
            subsettedDf = dataTypeDfDict[dataType].iloc[rowindexlist,:]
            #Move measuarables to columns
            postProcessingDf = subsettedDf.stack().unstack(dataTypeDfDict[dataType].index.names[:dataTypeObservableRangeDict[dataType]])
            indexList = []
            indexDict = {}
            for row in range(postProcessingDf.shape[0]):
                key = ','.join(list(map(str,list(postProcessingDf.iloc[row,:].name))))
                indexDict[key] = row
                indexList.append(key)
            allIndexList.append(indexList)
            allIndexDictList.append(indexDict)
            postProcessingCommonIndex = postProcessingDf.index
            print(postProcessingDf.index)
            postProcessingFeatures.append(list(postProcessingDf.columns))
            postProcessingMatrices.append(postProcessingDf.values)
        if len(dataTypeDfDict.keys()) > 1:
            result = list(set(allIndexList[0]).intersection(*allIndexList[1:]))
            reorderedResult = []
            for value in allIndexList[0]:
                if value in result:
                    reorderedResult.append(value)
            result = reorderedResult
        else:
            result = allIndexList[0]
        for i,postProcessingMatrix,indexDict in zip(range(len(postProcessingMatrices)),postProcessingMatrices,allIndexDictList):
            rows = []
            for key in result:
                rows.append(indexDict[key])
            postProcessingMatrices[i] = postProcessingMatrix[rows,:]
        fullPostProcessingMatrix = np.hstack(postProcessingMatrices)
        commonFeatures = [item for sublist in postProcessingFeatures for item in sublist]
        fullPostProcessingDf = pd.DataFrame(fullPostProcessingMatrix,index=postProcessingCommonIndex,columns=commonFeatures)
    else:
        subsettingDimensions = []
        for dim in finalDimensions:
            subsettingDimensions.append(dim[0])
        fullPostProcessingDf = newdf3.loc[subsettingDimensions]
    return fullPostProcessingDf
