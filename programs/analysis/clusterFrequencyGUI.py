#! /usr/bin/env python3
import json,pickle,math,matplotlib,sys,os,string,subprocess
from sys import platform as sys_pf
if sys_pf == 'darwin':
    import matplotlib
    matplotlib.use("TkAgg")
import numpy as np
import pandas as pd
import tkinter as tk
import tkinter.ttk
sys.path.insert(0, '../dataprocessing/')
from miscFunctions import setMaxWidth
sys.path.insert(0, '../plotting/')
from plottingGUI import PlotTypePage 

class ClusterFrequencyHomePage(tk.Frame):
    num_args = 2
    def __init__(self, master,fName,bp,shp):
        global folderName,secondaryhomepage,backpage
        folderName = fName
        secondaryhomepage = shp
        backpage = bp
        tk.Frame.__init__(self, master)
        
        mainWindow = tk.Frame(self)
        mainWindow.pack(side=tk.TOP,padx=10,pady=10)
        
        tk.Label(mainWindow,text='Select a single cell clustered dataframe: ').pack()
        clusteredDatasetVar = tk.StringVar()
        clusteredDatasets = [x for x in os.listdir('outputData/analysisFiles/clusteredData') if 'cluster' in x]
        clusteredDatasets = ['new'] + clusteredDatasets
        clusteredDatasetVarMenu = tk.OptionMenu(mainWindow,clusteredDatasetVar,*clusteredDatasets)
        clusteredDatasetVarMenu.pack()
        clusteredDatasetVar.set(clusteredDatasets[0])
        setMaxWidth(clusteredDatasets,clusteredDatasetVarMenu)

        def collectInputs(action):
            dataset = clusteredDatasetVar.get()
            if dataset == 'new':
                pass
            else:
                datasetName = dataset.split('.')[0]
                if action == 'create':
                    df = pd.read_pickle('outputData/analysisFiles/clusteredData/'+dataset)
                    createClusterFrequencyDataframe(datasetName,df)
                    tk.messagebox.showinfo("Info", "Cluster Frequency Dataframe for dataset: "+datasetName+ " completed!")
                else:
                    df = pd.read_pickle('outputData/analysisFiles/clusterFrequencyData/clusterFrequency-'+datasetName+'.pkl')
                    with open('misc/plottingParams.pkl','wb') as f:
                        pickle.dump({'df':df,'folderName':folderName,'homepage':ClusterFrequencyHomePage,'bp':backpage,'shp':secondaryhomepage},f)
                    with open('misc/normalPlottingBool.pkl','wb') as f:
                        pickle.dump(False,f)
                    master.switch_frame(PlotTypePage) 

        actionWindow = tk.Frame(self)
        actionWindow.pack()
        tk.Label(actionWindow,text='Select action:').grid(row=0,column=0,columnspan=2)
        tk.Button(actionWindow, text="Create dataframe",command=lambda: collectInputs('create')).grid(row=1,column=0)
        tk.Button(actionWindow, text="Plot dataframe",command=lambda: collectInputs('plot')).grid(row=1,column=1)
        
        buttonWindow = tk.Frame(self)
        buttonWindow.pack(side=tk.TOP,pady=10)

        tk.Button(buttonWindow, text="Back",command=lambda: master.switch_frame(backpage,folderName,secondaryhomepage)).grid(row=5,column=1)
        tk.Button(buttonWindow, text="Quit",command=quit).grid(row=5,column=2)

def createClusterFrequencyDataframe(datasetName,df):
    
    clusterdf = pd.DataFrame({'Cluster':list(df.index.get_level_values('Cluster'))},index=df.index).droplevel('Cluster')
    indexingDf = clusterdf.groupby(list(clusterdf.index.names)[:-1]).first()
    newdfTuples = []
    newDataMatrixList = []
    newDataMatrixList2 = []
    clusterOrder = []
    clusterDictKeys = list(pd.unique(df.index.get_level_values('Cluster')))
    clusterDictKeys.sort(key=int)
    emptyClusterDict = {}
    
    for clusterKey in clusterDictKeys:
        emptyClusterDict[clusterKey] = 0
    for row in range(indexingDf.shape[0]):
        index = indexingDf.iloc[row,:].name
        #Percent
        sampleDf = clusterdf.loc[index].squeeze().value_counts(normalize=True).mul(100).T
        #Count
        sampleDf2 = clusterdf.loc[index].squeeze().value_counts().T
        clusterDict = emptyClusterDict.copy()
        clusterDict2 = emptyClusterDict.copy()
        cvals = []
        cvals2 = []
        for clusterKey in sampleDf.index:
            clusterDict[clusterKey] = sampleDf.loc[clusterKey]
        for clusterKey2 in sampleDf2.index:
            clusterDict2[clusterKey2] = sampleDf2.loc[clusterKey2]
        newDataMatrixList.append(list(clusterDict.values()))
        newDataMatrixList2.append(list(clusterDict2.values()))
        newdfTuples.append(list(index))
    
    newDataMatrix = np.vstack(newDataMatrixList)
    newDataMatrix2 = np.vstack(newDataMatrixList2)
    mi = pd.MultiIndex.from_tuples(newdfTuples,names=indexingDf.index.names)
    percentdf = pd.DataFrame(newDataMatrix,index=mi,columns=list(clusterDict.keys()))
    countdf = pd.DataFrame(newDataMatrix2,index=mi,columns=list(clusterDict.keys()))
    
    frequencydf = pd.concat([percentdf,countdf],keys=['percent','count'],names=['Statistic']).swaplevel(0,1)
    frequencydf.columns.name = 'Cluster'
    if 'clusterFrequencyData' not in os.listdir('outputData/analysisFiles/'):
        subprocess.run(['mkdir','outputData/analysisFiles/clusterFrequencyData'])
    
    frequencydf.to_pickle('outputData/analysisFiles/clusterFrequencyData/clusterFrequency-'+datasetName+'.pkl')
