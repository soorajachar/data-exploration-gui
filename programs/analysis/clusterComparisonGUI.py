#!/usr/bin/env python3 
import json,pickle,math,matplotlib,sys,os,string
from sys import platform as sys_pf
if sys_pf == 'darwin':
    import matplotlib
    matplotlib.use("TkAgg")
import numpy as np
import pandas as pd
import seaborn as sns
import tkinter as tk
import tkinter.ttk
import itertools
from matplotlib import pyplot as plt
from scipy import stats
from sklearn.manifold import TSNE
from sklearn.preprocessing import MinMaxScaler
sys.path.insert(0, '../dataprocessing/')
from miscFunctions import sortSINumerically,reindexDataFrame,setMaxWidth,returnSpecificExtensionFiles,returnTicks,get_cluster_centroids,rainbow_text
sys.path.insert(0, '../plotting/')
from plottingGUI import checkUncheckAllButton,selectLevelValuesPage,createLabelDict
import facetPlotLibrary as fpl
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.widgets import RectangleSelector
from umap import UMAP
from sklearn.manifold import Isomap
from sklearn.decomposition import PCA
#import clusterPlottingLibrary as cpl
import itertools
from dataFrameValueSelectionGUI import DataSelectionHomePage
import operateOnDataSelection as ods
import clusterPlottingLibrary as cpl 
from dimensionReductionGUI import DimensionReductionHomePage 
sys.path.insert(0, '../plotting/')
import interactiveGUIElements as ipe

idx = pd.IndexSlice
letters = string.ascii_uppercase

class ClusterComparisonHomePage(tk.Frame):
    num_args = 2
    def __init__(self, master,fName,bp,shp):
        global folderName,secondaryhomepage,backpage
        folderName = fName
        secondaryhomepage = shp
        backpage = bp
        tk.Frame.__init__(self, master)
        
        mainWindow = tk.Frame(self)
        mainWindow.pack(side=tk.TOP,padx=10,pady=10)
        
        self.preprocessedDict = {}
        self.preprocessedDict2 = {}
        def getUpdateData(event):
            self.DimRedCombo['values'] = self.preprocessedDict[self.PreprocessedCombo.get()]
            if len(self.DimRedCombo['values']) == 1:
                self.DimRedCombo.set(self.DimRedCombo['values'][0])
            elif len(self.DimRedCombo['values']) == 2:
                self.DimRedCombo.set(self.DimRedCombo['values'][1])
            self.ClusterCombo['values'] = self.preprocessedDict2[self.PreprocessedCombo.get()]
            if len(self.ClusterCombo['values']) == 1:
                self.ClusterCombo.set(self.ClusterCombo['values'][0])
            elif len(self.ClusterCombo['values']) == 2:
                self.ClusterCombo.set(self.ClusterCombo['values'][1])

        maxFnLen = 0
        maxDfLen = 0
        maxDfLen2 = 0
        for scaledData in os.listdir('outputData/analysisFiles/scaledData'):
            dimreds = []
            clusters = []
            if '.DS' not in scaledData:
                if len(scaledData) > maxFnLen:
                    maxFnLen = len(scaledData)
                for dimRed in os.listdir('outputData/analysisFiles/reducedData'):
                    if scaledData.split('-')[0] == dimRed.split('-')[0]:
                        dimreds.append(dimRed)
                    if len(dimRed) > maxDfLen:
                        maxDfLen = len(dimRed)
                if not isinstance(dimreds,list):
                    dimreds = [dimreds]
                self.preprocessedDict[scaledData] = ['new']+dimreds
                for cluster in os.listdir('outputData/analysisFiles/clusteredData'):
                    if scaledData.split('-')[0] == cluster.split('-')[0]:
                        clusters.append(cluster)
                    if len(cluster) > maxDfLen2:
                        maxDfLen2 = len(cluster)
                if not isinstance(clusters,list):
                    clusters = [clusters]
                self.preprocessedDict2[scaledData] = ['none']+clusters

        l1 = tk.Label(mainWindow, text="""Select Preprocessed Subset: """).grid(row=0,column=0,sticky=tk.W)
        self.PreprocessedCombo = tkinter.ttk.Combobox(mainWindow,values = list(self.preprocessedDict.keys()))
        self.PreprocessedCombo['width'] = max([10,maxFnLen])
        self.PreprocessedCombo.bind('<<ComboboxSelected>>', getUpdateData)
        self.PreprocessedCombo.grid(row = 0,column = 1,sticky=tk.W)
        
        l2 = tk.Label(mainWindow, text="""Dimensional Reduction (for visualization): """).grid(row=1,column=0,sticky=tk.W)
        self.DimRedCombo = tkinter.ttk.Combobox(mainWindow,state='readonly')
        self.DimRedCombo['width'] = max([10,maxDfLen])
        self.DimRedCombo.grid(row = 1,column = 1)
        
        l3 = tk.Label(mainWindow, text="""Select clustered data subset: """).grid(row=2,column=0,sticky=tk.W)
        self.ClusterCombo = tkinter.ttk.Combobox(mainWindow,state='readonly')
        self.ClusterCombo['width'] = max([5,maxDfLen2])
        self.ClusterCombo.grid(row = 2,column = 1,sticky=tk.W)
        if len(self.PreprocessedCombo['values']) == 1:
            self.PreprocessedCombo.set(self.PreprocessedCombo['values'][0])
            self.DimRedCombo['values'] = self.preprocessedDict[self.PreprocessedCombo.get()]
            if len(self.DimRedCombo['values']) == 1:
                self.DimRedCombo.set(self.DimRedCombo['values'][0])
            elif len(self.DimRedCombo['values']) == 2:
                self.DimRedCombo.set(self.DimRedCombo['values'][1])
            self.ClusterCombo['values'] = self.preprocessedDict2[self.PreprocessedCombo.get()]
            self.ClusterCombo.set(self.ClusterCombo['values'][0])
        
        supervisedBool = tk.BooleanVar(value=False)
        supervisedCheckbox = tk.Checkbutton(mainWindow,text='Use clusters to supervise dimensional reduction (only for umap)',variable=supervisedBool)
        supervisedCheckbox.grid(row=2,column=2,sticky=tk.W)

        def collectInputs():
            dataSelectionFileName = self.PreprocessedCombo.get()
            dataSubsetTitle = dataSelectionFileName.split('-scaledBy')[0]
            scaledData = pickle.load(open('outputData/analysisFiles/scaledData/'+dataSelectionFileName,'rb'))
            #If we are not using cluster labels
            if self.ClusterCombo.get() == 'none':
                #If we are using a pre existing dimensional reduction, we load it and use it 
                if self.DimRedCombo.get() != 'new':
                    reductionFileName = self.DimRedCombo.get()
                    reducedData = pickle.load(open('outputData/analysisFiles/reducedData/'+reductionFileName,'rb'))
                    master.switch_frame(ClusterComparisonPage,scaledData,reducedData)
                #If we are not using a pre existing dimensional reduction, we make one
                else:
                    master.switch_frame(DimensionReductionHomePage,folderName,secondaryhomepage)
            #If we are using cluster labels
            else:
                clusterSelectionFileName = self.ClusterCombo.get()
                #new Justin: default is pkl; if generated from biowulf, you need to merge npy with scaled data
                if clusterSelectionFileName.endswith(".pkl"):
                    clusteredData = pickle.load(open('outputData/analysisFiles/clusteredData/'+clusterSelectionFileName,'rb'))
                elif clusterSelectionFileName.endswith(".npy"):
                    clusteredData = np.load(open('outputData/analysisFiles/clusteredData/'+clusterSelectionFileName,'rb'))
                    clusteredData = scaledData.assign(Cluster=list(map(str,clusteredData))).set_index('Cluster', append=True)
                #If we use pre-existing cluster labels to create a supervised umap, we will automatically be creating a new dimensional reduction; there will be no user input 
                if supervisedBool.get():
                    cluster_labels = list(clusteredData.index.get_level_values('Cluster'))
                    reductionDict = ods.parseFileName(self.DimRedCombo.get())
                    reducedData = ods.operateOnData(scaledData,dataSubsetTitle,'reduce',reductionDict,cluster_labels=cluster_labels)
                    master.switch_frame(ClusterComparisonPage,scaledData,reducedData,clusteredData)
                else:
                    reductionFileName = self.DimRedCombo.get()
                  #new Justin: default is a pkl.  If you genereated files from the biowulf, you would have an HDF, which needs minor renaming
                    if reductionFileName.endswith(".pkl"):
                        reducedData = pickle.load(open('outputData/analysisFiles/reducedData/'+reductionFileName,'rb'))
                    elif reductionFileName.endswith(".hdf"):
                        reducedData = pd.read_hdf('outputData/analysisFiles/reducedData/'+reductionFileName, key='df')
                        reducedData.rename(columns={'UMAP 1': 'Dimension 1', 'UMAP 2': 'Dimension 2'}, inplace=True)
                    master.switch_frame(ClusterComparisonPage,scaledData,reducedData,clusteredData)

        buttonWindow = tk.Frame(self)
        buttonWindow.pack(side=tk.TOP,pady=10)

        tk.Button(buttonWindow, text="OK",command=lambda: collectInputs()).grid(row=5,column=0)
        tk.Button(buttonWindow, text="Back",command=lambda: master.switch_frame(backpage,folderName,secondaryhomepage)).grid(row=5,column=1)
        tk.Button(buttonWindow, text="Quit",command=quit).grid(row=5,column=2)
        
class ClusterComparisonPage(tk.Frame):
    def __init__(self, master,scaledData,reducedData,clusteredData = []):
        #Initialize class
        self.root = master.root
        tk.Frame.__init__(self, master)
         
        #Initialize default variables
        clusterComparisonList = []
        clusterComparisonDict = {}
        clusterCentroidList2 = [] 
        clusterCentroidList = [] 
        clusterCentroidDict = {}
        clusterCentroidDict2 = {}
        palette = sns.color_palette().as_hex()
        self.paletteIndex = 0
        self.GUI_Selection_Index=0
        #Use order of variables in original dataframe
        experimentParametersBool = False
        for fn in os.listdir('misc'):
            for dataType in ['cell','cyt']:
                if 'experimentParameters' in fn:
                    if dataType in fn:
                        experimentParametersBool = True
                        experimentParameters = json.load(open('misc/experimentParameters-'+folderName+'-'+dataType+'.json','r'))
                        break
        if experimentParametersBool:
            trueLabelDict = experimentParameters['levelLabelDict'] 
        else:
            if 'Event' in scaledData.index.names or 'event' in scaledData.index.names:
                h5FileBool = False
                for fileName in os.listdir('outputData/pickleFiles'):
                    if '.h5' in fileName:
                        h5FileBool = True
                        break
                if h5FileBool:
                    ogname = 'outputData/pickleFiles/'+'initialSingleCellDf-channel-'+folderName+'.pkl'
                    newname = ogname.split('.')[0]+'.h5'
                    newdf = pd.read_hdf(newname, 'df')
                    originalDf = newdf
                else:
                    originalDf = pickle.load(open('outputData/pickleFiles/'+'initialSingleCellDf-channel-'+folderName+'.pkl','rb'))
            else:
                originalDf = scaledData.copy()
            trueLabelDict = createLabelDict(originalDf)
            with open('misc/experimentParameters-'+folderName+'-cell.json', 'w') as fp:
                json.dump({'levelLabelDict':trueLabelDict}, fp)

        figheight = 5 
        figwidth = 5
        widthscale = 1.4
        heightscale = 1.2
        plotFrameVisual = tk.Frame(self,borderwidth=1,relief='groove')
        plotFrameVisual.grid(row=0,column=0,sticky=tk.W+tk.E)
        figVisual = plt.figure(figsize=(figwidth*(widthscale+0.1), figheight*heightscale))
        figVisual.subplots_adjust(bottom=heightscale-1,left=widthscale-1)
        gsVisual = figVisual.add_gridspec(1, 1)
        maxFeatureLen = 0
        for feature in scaledData.columns:
            for splitfeature in feature.split(','):
                if len(splitfeature) > maxFeatureLen:
                    maxFeatureLen = len(splitfeature)
                    print(splitfeature)
        featureScaleFactor = 5 
        #featureScaleFactor = 5*2.5
        self.canvasVisual = FigureCanvasTkAgg(figVisual,master=plotFrameVisual)
        self.canvasVisual.draw()
        self.canvasVisual.get_tk_widget().pack()
        
        #Cluster visualization
        defaultGroupCompDict  = {}
        levelPlotAxis = figVisual.add_subplot(gsVisual[0])
        levelPlotAxis.set_title('Data Visualization')
        if isinstance(clusteredData,list):
            visualPlottingDfForLegend = reducedData.reset_index()
            visualPlottingDf = pd.concat([reducedData,scaledData],axis=1).reset_index()
            levelList = list(scaledData.index.names)+['None']
            featureList = list(scaledData.columns)
            defaultClusterBool = False
        else:
            visualPlottingDfForLegend = reducedData.reset_index()
            visualPlottingDf = pd.concat([reducedData,scaledData],axis=1).reset_index()
            visualPlottingDf['Cluster'] = list(clusteredData.index.get_level_values('Cluster'))
            levelList = list(clusteredData.index.names)+['None']
            featureList = list(scaledData.columns)
            defaultClusterBool = True
        defaultplotkwargs,defaultDict = ipe.getDefaultKwargs(visualPlottingDfForLegend)
        if defaultClusterBool:
            defaultDict['hue'] = 'Cluster'
        #Widgets
        #Frame
        visualizationParameterWindow = tk.Frame(self)
        visualizationParameterWindow.grid(row=1,column=0)
        levelPlotWindow = tk.Frame(visualizationParameterWindow)
        levelPlotWindow.grid(row=0,column=0)
        #Dropdowns
        levelParameterList = ['hue','style','size']
        levelParameterValueDict = {}
        for level in levelParameterList:
            if level == 'hue' or level == 'size':
                levelParameterValueDict[level] = levelList.copy()+featureList.copy()
            else:
                levelParameterValueDict[level] = levelList.copy()
        dropdownList,dropdownVarsDict,levelDropdown,levelValueDropdown = ipe.createParameterSelectionDropdownsWithIndividualLevels(levelPlotWindow,levelParameterList,levelParameterValueDict,defaultDict,visualPlottingDf,experimentParameters)
        #Plot
        sizeParam = {}
        if 'Event' in visualPlottingDf.columns or 'event' in visualPlottingDf.columns:
            sizeParam['s'] = 3 
        ipe.updateDropdownControlledPlot(self.canvasVisual,levelPlotAxis,visualPlottingDf,dropdownVarsDict,'Dimension 1','Dimension 2',alpha=0.4,legendoffset=-0.2,trueLabelDict = trueLabelDict,levelDropdown = levelDropdown,levelValueDropdown=levelValueDropdown)
        levelPlotAxis.set_title('Data Visualization')
        self.originalxlims = levelPlotAxis.get_xlim()
        self.originalylims = levelPlotAxis.get_ylim()
        self.currentxlims = levelPlotAxis.get_xlim()
        self.currentylims = levelPlotAxis.get_ylim()
        #Button
        def zoomIn():
            clusterSelectionBox = toggle_selector2.RS.corners
            ll = np.array([clusterSelectionBox[0][0], clusterSelectionBox[1][0]])  # lower-left
            ur = np.array([clusterSelectionBox[0][2], clusterSelectionBox[1][2]])  # upper-right
            
            inidx = np.all(np.logical_and(ll <= reducedData.values[:,:2], reducedData.values[:,:2] <= ur), axis=1)
            inbox = reducedData.loc[inidx]
            bufferval = 0.1
            xlims = [min(inbox['Dimension 1'])-bufferval,max(inbox['Dimension 1'])+bufferval]
            ylims = [min(inbox['Dimension 2'])-bufferval,max(inbox['Dimension 2'])+bufferval]
            self.currentxlims = xlims
            self.currentylims = ylims
            levelPlotAxis.set_xlim(self.currentxlims)
            levelPlotAxis.set_ylim(self.currentylims)
            clusterSelectionAxis.set_xlim(self.currentxlims)
            clusterSelectionAxis.set_ylim(self.currentylims)
            self.canvasVisual.draw()
            self.canvasSelection.draw()
        def zoomOut():
            levelPlotAxis.set_xlim(self.originalxlims)
            levelPlotAxis.set_ylim(self.originalylims)
            clusterSelectionAxis.set_xlim(self.originalxlims)
            clusterSelectionAxis.set_ylim(self.originalylims)
            self.currentxlims = self.originalxlims 
            self.currentylims = self.originalylims
            self.canvasVisual.draw()
            self.canvasSelection.draw()
        def update():
            ipe.updateDropdownControlledPlot(self.canvasVisual,levelPlotAxis,visualPlottingDf,dropdownVarsDict,'Dimension 1','Dimension 2',alpha=0.4,legendoffset=-0.2,trueLabelDict = trueLabelDict,levelDropdown=levelDropdown,levelValueDropdown=levelValueDropdown,axisLimits = [levelPlotAxis.get_xlim(),levelPlotAxis.get_ylim()])
            levelPlotAxis.set_title('Data Visualization')
            levelPlotAxis.set_xlim(self.currentxlims)
            levelPlotAxis.set_ylim(self.currentylims)
            clusterSelectionAxis.set_xlim(self.currentxlims)
            clusterSelectionAxis.set_ylim(self.currentylims)
            self.canvasVisual.draw()
            self.canvasSelection.draw()
        #Click and drag widget
        def line_select_callback2(eclick, erelease):
            'eclick and erelease are the press and release events'
            x1, y1 = eclick.xdata, eclick.ydata
            x2, y2 = erelease.xdata, erelease.ydata
        def toggle_selector2(event):
            if event.key in ['Q', 'q'] and toggle_selector2.RS.active:
                toggle_selector2.RS.set_active(False)
            if event.key in ['A', 'a'] and not toggle_selector2.RS.active:
                toggle_selector2.RS.set_active(True)
        rectpropsdict2 = {'facecolor':'grey','alpha':0.2,'edgecolor':'grey'}
        toggle_selector2.RS = RectangleSelector(levelPlotAxis, line_select_callback2,drawtype='box', useblit=True,button=[1, 3], minspanx=5, minspany=5,spancoords='pixels',interactive=True,rectprops=rectpropsdict2)
        self.ts2 = toggle_selector2.RS
        self.canvasVisual.mpl_connect('key_press_event', toggle_selector2)
        
        zoomWindow = tk.Frame(visualizationParameterWindow)
        zoomWindow.grid(row=0,column=1)
        tk.Button(zoomWindow, text="Zoom in",command=lambda: zoomIn()).grid(row=0,column=0,sticky=tk.W)
        tk.Button(zoomWindow, text="Zoom out",command=lambda: zoomOut()).grid(row=1,column=0,sticky=tk.W)
        
        tk.Button(visualizationParameterWindow, text="Update visualization plot",command=lambda: update()).grid(row=1,column=0,columnspan=2)
        tk.Button(visualizationParameterWindow, text="Save visualization plot",command=lambda: exportFigure('visualization')).grid(row=2,column=0,columnspan=2)

        #Cluster selection
        plotFrameSelection = tk.Frame(self,borderwidth=1,relief='groove')
        plotFrameSelection.grid(row=0,column=1,sticky=tk.W+tk.E)
        figSelection = plt.figure(figsize=(figwidth, figheight*heightscale))
        figSelection.subplots_adjust(bottom=heightscale-1)#,left=widthscale-1)
        gsSelection = figSelection.add_gridspec(1, 1)
        self.canvasSelection = FigureCanvasTkAgg(figSelection,master=plotFrameSelection)
        self.canvasSelection.draw()
        self.canvasSelection.get_tk_widget().pack()
        #Plot
        clusterSelectionAxis = figSelection.add_subplot(gsSelection[0])
        clusterSelectionAxis.set_title('Group Selection')
        reducedPlottingDf = reducedData.reset_index()
        g1 = sns.scatterplot(data=reducedPlottingDf,x='Dimension 1',y='Dimension 2',ax=clusterSelectionAxis,alpha=0.7,color='#808080',**sizeParam)
        if clusterSelectionAxis.legend_ is not None:
            clusterSelectionAxis.legend_.remove()
        #Click and drag widget
        def line_select_callback(eclick, erelease):
            'eclick and erelease are the press and release events'
            x1, y1 = eclick.xdata, eclick.ydata
            x2, y2 = erelease.xdata, erelease.ydata

        def toggle_selector(event):
            if event.key in ['Q', 'q'] and toggle_selector.RS.active:
                toggle_selector.RS.set_active(False)
            if event.key in ['A', 'a'] and not toggle_selector.RS.active:
                toggle_selector.RS.set_active(True)
        rectpropsdict = {'facecolor':palette[self.paletteIndex],'alpha':0.2,'edgecolor':palette[self.paletteIndex]}
        toggle_selector.RS = RectangleSelector(clusterSelectionAxis, line_select_callback,drawtype='box', useblit=True,button=[1, 3], minspanx=5, minspany=5,spancoords='pixels',interactive=True,rectprops=rectpropsdict)
        self.ts = toggle_selector.RS
        self.canvasSelection.mpl_connect('key_press_event', toggle_selector)
        
        #Functions
        def updateSelectionPlot():
            featureDfList = []
            selectionPlottingDf = reducedPlottingDf.copy()
            selectionPlottingDf['Cluster'] = ['0']*reducedData.shape[0]
            for i in clusterComparisonDict:
                featureDfList.append(scaledData.iloc[clusterComparisonDict[i],:])
                selectionPlottingDf.loc[:,'Cluster'].iloc[clusterComparisonDict[i]] = str(i+1)
            clusterSelectionAxis.clear()
            #selectionPlottingDf['Cluster'] = clusterLabelList
            newpalette = ['#808080']+palette[:len(clusterComparisonDict.keys())]
            if len(selectionPlottingDf[selectionPlottingDf['Cluster'] == '0']) == 0:
                newpalette = newpalette[1:]
            g1 = sns.scatterplot(data=selectionPlottingDf,x='Dimension 1',y='Dimension 2',hue='Cluster',ax=clusterSelectionAxis,alpha=0.7,palette=newpalette,**sizeParam)
            if clusterSelectionAxis.legend_ is not None:
                clusterSelectionAxis.legend_.remove()
            for clusterCentroid in clusterCentroidList:
                g1.annotate(clusterCentroid[0],xy=clusterCentroid[1])
            clusterSelectionAxis.set_xlim(self.currentxlims)
            clusterSelectionAxis.set_ylim(self.currentylims)
            #rainbow_text(0.2, 1.05, "Group1;vs.;4,5;vs.;6".split(';'),[palette[0],'black', palette[1],'black', palette[2]],size=18)
            #returnGroupTitle(clusterSelectionAxis)
            clusterSelectionAxis.set_title('Group Selection')
            self.canvasSelection.draw()
        def addPredeterminedToCluster():
            predeterminedCluster = clusterVar.get()

            tempDf = pd.concat([reducedData,scaledData],axis=1).reset_index()
            tempDf['Cluster'] = list(clusteredData.index.get_level_values('Cluster'))
            #tempDf = visualPlottingDf.copy()
            tempDf['index'] = range(tempDf.shape[0])
            #inidx = np.all(tempDf['Cluster'] == predeterminedCluster)
            if predeterminedCluster == 'All Not Selected':
                tempSelDf = reducedPlottingDf.copy() 
                tempSelDf['Cluster'] = ['0']*reducedData.shape[0]
                for clusterIndex in clusterComparisonDict:
                    tempSelDf.loc[:,'Cluster'].iloc[clusterComparisonDict[clusterIndex]] = '-1' 
                inbox = tempDf.loc[tempSelDf['Cluster'] != '-1']
            else:
                inbox = tempDf.loc[tempDf['Cluster'] == predeterminedCluster]
            predetCentroid = get_cluster_centroids(inbox,singleCluster=True)
            tupleList = []
            if self.paletteIndex not in clusterComparisonDict.keys():
                clusterComparisonDict[self.paletteIndex] = list(map(int,list(inbox.values[:,-1])))
            else:
                clusterComparisonDict[self.paletteIndex] += list(map(int,list(inbox.values[:,-1])))
            clusterCentroidList.append(predetCentroid[0])
            clusterCentroidList2.append(predetCentroid[0])
            updateSelectionPlot()
        def addToCluster():
            clusterSelectionBox = toggle_selector.RS.corners
            ll = np.array([clusterSelectionBox[0][0], clusterSelectionBox[1][0]])  # lower-left
            ur = np.array([clusterSelectionBox[0][2], clusterSelectionBox[1][2]])  # upper-right
            
            tempDf = reducedData.copy()
            tempDf['index'] = range(tempDf.shape[0])

            inidx = np.all(np.logical_and(ll <= tempDf.values[:,:2], tempDf.values[:,:2] <= ur), axis=1)
            inbox = tempDf.loc[inidx]
            selectionCentroid = get_cluster_centroids(inbox,singleCluster=True)
            selectionCentroid[0][0] = letters[self.GUI_Selection_Index] 
            self.GUI_Selection_Index+=1
            tupleList = []
            if self.paletteIndex not in clusterComparisonDict.keys():
                clusterComparisonDict[self.paletteIndex] = list(map(int,list(inbox.values[:,-1])))
            else:
                clusterComparisonDict[self.paletteIndex] += list(map(int,list(inbox.values[:,-1])))
            clusterCentroidList.append(selectionCentroid[0])
            clusterCentroidList2.append(selectionCentroid[0])
            updateSelectionPlot()
        def returnGroupTitle(axis):
            #rainbow_text(0.2, 1.05, "1,2,3;vs.;4,5;vs.;6".split(';'),[palette[0],'black', palette[1],'black', palette[2]],size=18)
            titleColorList = []
            fullTitleList = []
            palette = sns.color_palette(sns.color_palette(),len(clusterCentroidDict))
            for i,group in enumerate(clusterCentroidDict):
                titleList = []
                for clusterCentroid in clusterCentroidDict[group]:
                    fullName = clusterCentroid[0]
                    titleList.append(fullName)
                selectionsInGroup = ','.join(titleList)
                fullTitleList.append(selectionsInGroup)
                titleColorList+=[palette[i],'black']
            fullTitleString = ';vs.;'.join(fullTitleList)
            titleColorList = titleColorList[:-1]
            axis.set_title('')
            xpos = ((max(visualPlottingDf['Dimension 1'])-min(visualPlottingDf['Dimension 1']))/2)+min(visualPlottingDf['Dimension 1'])
            ypos = max(visualPlottingDf['Dimension 2'])+2
            rainbow_text(xpos, ypos, fullTitleString.split(';'),titleColorList,size=10,ax=axis)

        def grabSamplesInCluster(): 
            clusterComparisonList.append(clusterComparisonDict[self.paletteIndex])
            clusterCentroidDict[str(self.paletteIndex)] = clusterCentroidList.copy()
            clusterCentroidDict2[str(self.paletteIndex)] = clusterCentroidList2.copy()
            del clusterCentroidList2[:]
            self.paletteIndex+=1
            #returnGroupTitle(clusterSelectionAxis)
            rectpropsdict = {'facecolor':palette[self.paletteIndex],'alpha':0.2,'edgecolor':palette[self.paletteIndex]}
            toggle_selector.RS = RectangleSelector(clusterSelectionAxis, line_select_callback,drawtype='box', useblit=True,button=[1, 3], minspanx=5, minspany=5,spancoords='pixels',interactive=True,rectprops=rectpropsdict)
            self.ts = toggle_selector.RS
            self.canvasSelection.draw()
        def clearClusters():
            clusterComparisonAxis.clear()
            clusterSelectionAxis.clear()
            groupCompositionAxis.clear()
            clusterSelectionAxis.set_title('Group Selection')
            clusterComparisonAxis.set_title('Group Comparison')
            groupCompositionAxis.set_title('Group Composition')
            g1 = sns.scatterplot(data=reducedPlottingDf,x='Dimension 1',y='Dimension 2',ax=clusterSelectionAxis,alpha=0.7,color='#808080',**sizeParam)
            if clusterSelectionAxis.legend_ is not None:
                clusterSelectionAxis.legend_.remove()
            del clusterComparisonList[:]
            clusterComparisonDict.clear()
            del clusterCentroidList[:]
            clusterCentroidDict.clear()
            del clusterCentroidList2[:]
            clusterCentroidDict2.clear()
            self.paletteIndex = 0
            self.GUI_Selection_Index=0
            rectpropsdict = {'facecolor':palette[self.paletteIndex],'alpha':0.2,'edgecolor':palette[self.paletteIndex]}
            toggle_selector.RS = RectangleSelector(clusterSelectionAxis, line_select_callback,drawtype='box', useblit=True,button=[1, 3], minspanx=5, minspany=5,spancoords='pixels',interactive=True,rectprops=rectpropsdict)
            self.ts = toggle_selector.RS
            clusterSelectionAxis.set_xlim(self.currentxlims)
            clusterSelectionAxis.set_ylim(self.currentylims)
            self.canvasSelection.draw()
        def compareClusters():
            updateSelectionPlot()
            updateComparisonPlot(sliderList,radiobuttonVars,radiobuttonVars2)
        #Widgets
        #Frame
        selectionParameterWindow = tk.Frame(self)
        selectionParameterWindow.grid(row=1,column=1)
        plotButtonWindow = tk.Frame(selectionParameterWindow)
        plotButtonWindow.grid(row=0,column=0)
        plotButtonWindow.root = master.root
        #Buttons
        #Interactive selection
        tk.Label(plotButtonWindow,text='GUI Selection').grid(row=0,column=0)
        tk.Button(plotButtonWindow, text="Add To Current Group",command=lambda: addToCluster()).grid(row=1,column=0,sticky=tk.W)
        
        #Cluster dropdown selection
        preDetWindow = tk.Frame(plotButtonWindow)
        preDetWindow.grid(row=0,column=1)
        if not isinstance(clusteredData,list):
            tk.Label(preDetWindow,text='Cluster: ').grid(row=0,column=0,sticky=tk.E)
            clusters = list(pd.unique(clusteredData.index.get_level_values('Cluster')))
            sortedClusters = list(map(str,sorted(list(map(int,clusters)))))+['All Not Selected']
            clusterVar = tk.StringVar()
            clusterMenu = tk.OptionMenu(preDetWindow,clusterVar,*sortedClusters)
            clusterMenu.grid(row=0,column=1)
            setMaxWidth(sortedClusters,clusterMenu)
            tk.Button(plotButtonWindow, text="Add To Current Group",command=lambda: addPredeterminedToCluster()).grid(row=1,column=1,sticky=tk.E)

        #Common to both
        tk.Button(plotButtonWindow, text="Store Current Group",command=lambda: grabSamplesInCluster()).grid(row=2,column=0,columnspan=2)
        tk.Button(plotButtonWindow, text="Clear Groups",command=lambda: clearClusters()).grid(row=3,column=0,columnspan=2)
        tk.Button(selectionParameterWindow, text="Compare Groups",command=lambda: compareClusters()).grid(row=1,column=0)
        tk.Button(selectionParameterWindow, text="Save selection plot",command=lambda: exportFigure('selection')).grid(row=2,column=0)
        
        #Group Comparison
        plotFrameComparison = tk.Frame(self,borderwidth=1,relief='groove')
        plotFrameComparison.grid(row=0,column=2,sticky=tk.W+tk.E)
        figComparison = plt.figure(figsize=(2*figwidth, figheight*heightscale))
        figComparison.subplots_adjust(bottom=heightscale-1)#,left=widthscale-1)
        gsComparison = figComparison.add_gridspec(1, 2)
        self.canvasComparison = FigureCanvasTkAgg(figComparison,master=plotFrameComparison)
        self.canvasComparison.draw()
        self.canvasComparison.get_tk_widget().pack()
        #Plot
        clusterComparisonAxis = figComparison.add_subplot(gsComparison[0:2])
        clusterComparisonAxis.set_title('Group Comparison')
        #Group composition
        self.comparisonPlottingDf = []
        #Functions
        def updateComparisonPlot(sliderList,rbVarList,rbVarList2):
            clusterComparisonAxis.clear()
            featureDfList = []
            featureDfWithClusters = scaledData.copy()
            featureDfWithClusters['Cluster'] = ['0']*scaledData.shape[0]
            for i,clusterComparison in enumerate(clusterComparisonList):
                featureDfList.append(scaledData.iloc[clusterComparison,:])
                featureDfWithClusters.loc[:,'Cluster'].iloc[clusterComparison] = str(i+1)
            clusterLabelList = featureDfWithClusters['Cluster'] 
            comparisonDf = pd.concat(featureDfList,keys=list(map(str,list(range(1,len(clusterComparisonList)+1)))),names=['Cluster'])
            comparisonDf.columns.name = 'Feature'
            comparisonPlottingDf = comparisonDf.stack().to_frame('Metric').reset_index()
            self.comparisonPlottingDf = comparisonPlottingDf.copy()
            #featureDfWithClusters['Cluster'] = clusterLabelList
            #Get unique cluster labels for hue order and pairwise significance comparisons
            clusterLabelListNoZero = list(filter(lambda a: a != '0', clusterLabelList))
            uniqueClusterLabelListNoZero = list(pd.unique(clusterLabelListNoZero))
            orderedClusterLabelListNoZero = []
            for uniqueClusterLabel in uniqueClusterLabelListNoZero:
                intval = int(uniqueClusterLabel)
                orderedClusterLabelListNoZero.append(intval)
            orderedClusterLabelListNoZero = sorted(orderedClusterLabelListNoZero)
            for i,orderedClusterLabel in enumerate(orderedClusterLabelListNoZero):
                orderedClusterLabelListNoZero[i] = str(orderedClusterLabel)
            #Get all pairwise combinations of clusters
            pairwiseClusterCombinations = list(itertools.combinations(orderedClusterLabelListNoZero,2))
            #Conduct significance testing to get features that define cluster differences
            comparisonKwargs = {'confidence':sliderList[0].get(),'foldchange':sliderList[1].get(),'responseCutoff':sliderList[2].get(),'errorCorrection':rbVarList['error correction method'].get()}
            self.comparisonKwargs = comparisonKwargs.copy()
            featureDfWithClusters.columns.name = 'Feature'
            significanceArray = []
            if len(clusterComparisonList) > 1:
                dataMatrix,significanceArray = cpl.significanceTesting(featureDfWithClusters,pairwiseClusterCombinations,**comparisonKwargs)
            else:
                significanceArray = list(scaledData.columns)
                print(significanceArray)
            if len(significanceArray) != 0:
                comparisonPlottingDf = comparisonPlottingDf[comparisonPlottingDf['Feature'].isin(significanceArray)]
                plotType = rbVarList2['comparison plot type'].get()
                self.comparisonKwargs['plotType'] = plotType
                comparisonPlottingDfTemp = comparisonPlottingDf.copy()
                newcols = []
                for col in comparisonPlottingDf['Feature']:
                    newcols.append(col.replace(',','\n'))
                comparisonPlottingDfTemp['Feature'] = newcols
                comparisonPlottingDfTemp['Metric'] = comparisonPlottingDfTemp['Metric'].astype(float)
                if plotType == 'violin':
                    if len(clusterComparisonList) == 1:
                        g2 = sns.violinplot(data=comparisonPlottingDfTemp,x='Feature',y='Metric',ax=clusterComparisonAxis,scale='width',inner='quartile',legend=False)
                    elif len(clusterComparisonList) == 2:
                        g2 = sns.violinplot(data=comparisonPlottingDfTemp,hue='Cluster',hue_order=orderedClusterLabelListNoZero,x='Feature',y='Metric',ax=clusterComparisonAxis,split=True,scale='width',inner='quartile',legend=False)
                    else:
                        g2 = sns.violinplot(data=comparisonPlottingDfTemp,hue='Cluster',hue_order=orderedClusterLabelListNoZero,x='Feature',y='Metric',ax=clusterComparisonAxis,dodge=True,scale='width',inner='quartile',legend=False)
                elif plotType == 'box':
                    g2 = sns.boxplot(data=comparisonPlottingDfTemp,hue='Cluster',hue_order=orderedClusterLabelListNoZero,x='Feature',y='Metric',ax=clusterComparisonAxis,showfliers=False)
                elif plotType == 'bar':
                    g2 = sns.barplot(data=comparisonPlottingDfTemp,hue='Cluster',hue_order=orderedClusterLabelListNoZero,x='Feature',y='Metric',ax=clusterComparisonAxis,ci='sd',errwidth=1,capsize=0.05,edgecolor='k')
                plt.setp(clusterComparisonAxis.xaxis.get_majorticklabels(), rotation=45)
                if len(clusterComparisonList) > 1:
                    clusterComparisonAxis.legend_.remove()
                if max(comparisonPlottingDfTemp['Metric']) > 100:
                    yticksToUse = [-1000,100,1000,10000,100000]
                    ytickValues,ytickLabels = returnTicks(yticksToUse)
                    #clusterComparisonAxis.set_ylim([0,1000])
                    clusterComparisonAxis.set_yticks(ytickValues)
                    clusterComparisonAxis.set_yticklabels(ytickLabels)
                    #clusterComparisonAxis.set_ylim([0,1000])
            else:
                clusterComparisonAxis.text(0.5, 0.5,'No significant differences',horizontalalignment='center',verticalalignment='center',transform = clusterComparisonAxis.transAxes)
            clusterComparisonAxis.set_title('Group Comparison')
            self.canvasComparison.draw()
        #Widgets
        #Frame
        clusterComparisonParameterWindow = tk.Frame(self)
        clusterComparisonParameterWindow.grid(row=1,column=2)
        sliderWindow = tk.Frame(clusterComparisonParameterWindow)
        sliderWindow.grid(row=0,column=0)
        sliderWindow.root = master.root
        #Sliders
        sliderList = ipe.createParameterAdjustmentSliders(sliderWindow,ods.clusterComparisonParameterList,ods.clusterComparisonParameterBoundsDict)
        #sliderList = ipe.createParameterAdjustmentSliders(clusterComparisonParameterWindow,ods.clusterComparisonParameterList,ods.clusterComparisonParameterBoundsDict)
        #Radio buttons (for plot type)
        radioWindow = tk.Frame(clusterComparisonParameterWindow)
        radioWindow.grid(row=0,column=1)
        radioWindow.root = master.root
        radiobuttonList2,radiobuttonVars2 = ipe.createParameterSelectionRadiobuttons(radioWindow,ods.clusterComparisonParameterList3,ods.clusterComparisonParameterValueDict3)
        
        #Radio buttons (for error correction)
        radioWindow2 = tk.Frame(clusterComparisonParameterWindow)
        radioWindow2.grid(row=0,column=2)
        radioWindow2.root = master.root
        radiobuttonList,radiobuttonVars = ipe.createParameterSelectionRadiobuttons(radioWindow2,ods.clusterComparisonParameterList2,ods.clusterComparisonParameterValueDict2)
        #Buttons
        tk.Button(clusterComparisonParameterWindow, text="Update comparison plot",command=lambda: updateComparisonPlot(sliderList,radiobuttonVars,radiobuttonVars2)).grid(row=1,column=0,columnspan=3)
        tk.Button(clusterComparisonParameterWindow, text="Save comparison plot",command=lambda: exportFigure('comparison')).grid(row=2,column=0,columnspan=3)
        
        #Group Composition
        plotFrameComposition = tk.Frame(self,borderwidth=1,relief='groove')
        plotFrameComposition.grid(row=0,column=3,sticky=tk.W+tk.E)
        figComposition = plt.figure(figsize=(figwidth*widthscale, figheight*heightscale))
        figComposition.subplots_adjust(bottom=heightscale-1,right=2-widthscale)
        gsComposition = figComposition.add_gridspec(1, 1)
        self.canvasComposition = FigureCanvasTkAgg(figComposition,master=plotFrameComposition)
        self.canvasComposition.draw()
        self.canvasComposition.get_tk_widget().pack()
        #Widgets
        groupCompositionAxis = figComposition.add_subplot(gsComposition[0])
        groupCompositionAxis.set_title('Group Composition')
        #Frame
        groupCompositionWindow = tk.Frame(self)
        groupCompositionWindow.grid(row=1,column=3)
        #Dropdowns
        groupCompParameterList = ['x','y','hue']
        groupCompParameterValueDict = {}
        #levelList2 = ['Cluster']+list(scaledData.index.names)+['None']
        levelList2 = ['Group']+list(scaledData.index.names)+['None']
        for level in groupCompParameterList:
            if level == 'x':
                groupCompParameterValueDict[level] = levelList2.copy()+list(scaledData.columns)
            elif level == 'y':
                groupCompParameterValueDict[level] = ['frequency','log-frequency','percent']+list(scaledData.columns)+['None']
            else:
                groupCompParameterValueDict[level] = levelList2.copy()+list(scaledData.columns)
        compositionDropdownWindow = tk.Frame(groupCompositionWindow)
        compositionDropdownWindow.root = master.root
        compositionDropdownWindow.grid(row=0,column=0)
        dropdownListComposition,dropdownVarsDictComposition = ipe.createParameterSelectionDropdowns(compositionDropdownWindow,groupCompParameterList,groupCompParameterValueDict,{'x':'Cluster','y':'frequency','hue':list(scaledData.index.names)[0]})
        
        def updateComposition():
            if isinstance(self.comparisonPlottingDf,list):
                featureDfList = []
                featureDfWithClusters = scaledData.copy()
                featureDfWithClusters['Cluster'] = ['0']*scaledData.shape[0]
                for i,clusterComparison in enumerate(clusterComparisonList):
                    featureDfList.append(scaledData.iloc[clusterComparison,:])
                    featureDfWithClusters.loc[:,'Cluster'].iloc[clusterComparison] = str(i+1)
                clusterLabelList = featureDfWithClusters['Cluster'] 
                comparisonDf = pd.concat(featureDfList,keys=list(map(str,list(range(1,len(clusterComparisonList)+1)))),names=['Cluster'])
                comparisonDf.columns.name = 'Feature'
                comparisonPlottingDf = comparisonDf.stack().to_frame('Metric').reset_index()
                self.comparisonPlottingDf = comparisonPlottingDf.copy()
            #compositionParameter = self.compositionRbs['Cluster composition y axis'].get()
            #ipe.updatedropdowncontrolledcompositionplot(self.canvascomposition,groupcompositionaxis,self.comparisonplottingdf,truelabeldict,dropdownvarsdictcomposition,compositionparameter)
            ipe.updateDropdownControlledCompositionPlot(self.canvasComposition,groupCompositionAxis,self.comparisonPlottingDf,trueLabelDict,dropdownVarsDictComposition)
            print(ipe.getDropdownValues(dropdownVarsDictComposition))
            print(clusterCentroidDict2)
            groupCompositionAxis.set_title('Group Composition')
            self.canvasComposition.draw()
        
        #Radio buttons (for composition y axis)
        #compositionRadioButtonWindow = tk.Frame(groupCompositionWindow)
        #compositionRadioButtonWindow.root = master.root
        #compositionRadioButtonWindow.grid(row=0,column=1)
        #compositionParameters = ['Cluster composition y axis'] 
        #compositionParameterValues = {'Cluster composition y axis':['frequency','percent']}
        #radiobuttonListComp,self.compositionRbs = ipe.createParameterSelectionRadiobuttons(compositionRadioButtonWindow,compositionParameters,compositionParameterValues)
        #Plot
        tk.Button(groupCompositionWindow, text="Update composition plot",command=lambda: updateComposition()).grid(row=1,column=0)
        tk.Button(groupCompositionWindow, text="Save composition plot",command=lambda: exportFigure('composition')).grid(row=2,column=0)

        #Figure exporting frame
        titleFrame = tk.Frame(self)
        titleFrame.grid(row=3,column=0,columnspan=4)
        titleEntryFrame = tk.Frame(titleFrame)
        titleEntryFrame.grid(row=0,column=0)
        titleEntryFrame.root = master.root
        tk.Label(titleEntryFrame,text='Title:').grid(row=0,column=0,sticky=tk.W)
        titleEntry = tk.Entry(titleEntryFrame)
        titleEntry.grid(row=0,column=1,sticky=tk.W)
        titleEntry.insert(tk.END,'group1vs2')
        
        def createDescriptiveTitle(figureName):
            #Group title segment
            groupList = []
            for group in clusterCentroidDict2:
                memberlist = []
                centroids = clusterCentroidDict2[group]
                for centroid in centroids:
                    memberlist.append(centroid[0])
                groupList.append(','.join(memberlist))
            #Figure parameter title segment
            if figureName == 'visualization':
                parameters1 = ipe.getDropdownValues(dropdownVarsDict)
                if levelValueDropdown.get() == 'All':
                    parameters2 = {'subset':levelValueDropdown.get()} 
                else:
                    parameters2 = {'subset':levelDropdown.get()+'-'+levelValueDropdown.get()} 
                parameters = {**parameters1,**parameters2}
            elif figureName == 'comparison':
                parameters = self.comparisonKwargs.copy()
            elif figureName == 'composition':
                parameters = ipe.getDropdownValues(dropdownVarsDictComposition)
            else:
                parameters = {}
            groupTitle = 'vs'.join(groupList)
            parameterTitle = ods.returnParameterString(parameters)
            return groupTitle,parameterTitle

        def exportFigure(figureName):
            figureDict = {'visualization':figVisual,'selection':figSelection,'comparison':figComparison,'composition':figComposition}
            customTitle = titleEntry.get()
            groupTitle,parameterTitle = createDescriptiveTitle(figureName)
            if figureName != 'all':
                fullFigureTitle = 'plots/'+'-'.join([customTitle,groupTitle,figureName,parameterTitle])+'.'+self.plotFormatRbs['Figure Format'].get()
                figureDict[figureName].savefig(fullFigureTitle,bbox_inches='tight')
                print(figureName+' Figure Saved')
            else:
                for figureNm in list(figureDict.keys()):
                    exportFigure(figureNm)
        
        plotFormatRadioButtonWindow = tk.Frame(titleFrame)
        plotFormatRadioButtonWindow.root = master.root
        plotFormatRadioButtonWindow.grid(row=0,column=1)
        plotFormatParameters = ['Figure Format'] 
        plotFormatParameterValues = {'Figure Format':['png','pdf']}
        radiobuttonListPlotFormat,self.plotFormatRbs = ipe.createParameterSelectionRadiobuttons(plotFormatRadioButtonWindow,plotFormatParameters,plotFormatParameterValues)
        tk.Button(titleFrame, text="Export all figures",command=lambda: exportFigure('all')).grid(row=1,column=0,columnspan=2)
        
        #Default save and quit buttons
        buttonWindow = tk.Frame(self)
        buttonWindow.grid(row=4,column=0,columnspan=4)
        tk.Button(buttonWindow, text="OK",command=lambda: master.switch_frame(backpage,folderName,secondaryhomepage)).grid(row=0,column=0)
        tk.Button(buttonWindow, text="Back",command=lambda: master.switch_frame(ClusterComparisonHomePage,folderName,backpage,secondaryhomepage)).grid(row=0,column=1)
        tk.Button(buttonWindow, text="Quit",command=lambda: quit()).grid(row=0,column=2)
