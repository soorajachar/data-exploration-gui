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
from miscFunctions import sortSINumerically,reindexDataFrame,setMaxWidth,returnSpecificExtensionFiles,returnTicks
sys.path.insert(0, '../plotting/')
from plottingGUI import checkUncheckAllButton,selectLevelValuesPage
import facetPlotLibrary as fpl
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.widgets import RectangleSelector
from umap import UMAP
from sklearn.manifold import Isomap
from sklearn.decomposition import PCA
import clusterPlottingLibrary as cpl
import itertools
import operateOnDataSelection as ods
import interactiveGUIElements as ipe
from dimensionReductionGUI import DimensionReductionHomePage 
idx = pd.IndexSlice

class ClusteringHomePage(tk.Frame):
    num_args = 2
    def __init__(self, master,fName,bp,shp):
        global folderName,backpage,secondaryhomepage
        folderName = fName
        secondaryhomepage = shp
        backpage = bp
        tk.Frame.__init__(self, master)
        
        mainWindow = tk.Frame(self)
        mainWindow.pack(side=tk.TOP,padx=10,pady=10)
        
        self.preprocessedDict = {}
        def getUpdateData(event):
            self.DimRedCombo['values'] = self.preprocessedDict[self.PreprocessedCombo.get()]
            if len(self.DimRedCombo['values']) == 1:
                self.DimRedCombo.set(self.DimRedCombo['values'][0])
            elif len(self.DimRedCombo['values']) == 2:
                self.DimRedCombo.set(self.DimRedCombo['values'][1])

        maxFnLen = 0
        maxDfLen = 0
        for scaledData in os.listdir('outputData/analysisFiles/scaledData'):
            dimreds = []
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

        l1 = tk.Label(mainWindow, text="""Select Preprocessed Subset: """).grid(row=0,column=0,sticky=tk.W)
        self.PreprocessedCombo = tkinter.ttk.Combobox(mainWindow,values = list(self.preprocessedDict.keys()))
        self.PreprocessedCombo['width'] = maxFnLen 
        self.PreprocessedCombo.bind('<<ComboboxSelected>>', getUpdateData)
        self.PreprocessedCombo.grid(row = 0,column = 1,sticky=tk.W)
        
        l2 = tk.Label(mainWindow, text="""Dimensional Reduction (for visualization): """).grid(row=1,column=0,sticky=tk.W)
        self.DimRedCombo = tkinter.ttk.Combobox(mainWindow,state='readonly')
        self.DimRedCombo['width'] = maxDfLen
        self.DimRedCombo.grid(row = 1,column = 1)
        if len(self.PreprocessedCombo['values']) == 1:
            self.PreprocessedCombo.set(self.PreprocessedCombo['values'][0])
            self.DimRedCombo['values'] = self.preprocessedDict[self.PreprocessedCombo.get()]
            if len(self.DimRedCombo['values']) == 1:
                self.DimRedCombo.set(self.DimRedCombo['values'][0])
            elif len(self.DimRedCombo['values']) == 2:
                self.DimRedCombo.set(self.DimRedCombo['values'][1])

        l3 = tk.Label(mainWindow, text="""Clustering Method: """).grid(row=2,column=0,sticky=tk.W)
        v3 = tk.StringVar()
        v3.set(list(ods.clusteringFunctionDict.keys())[0])
        clusterRbList = []
        for i,clusteringFunc in enumerate(ods.clusteringFunctionDict):
            rb = tk.Radiobutton(mainWindow,text=clusteringFunc,padx = 20, variable=v3, value=clusteringFunc)
            rb.grid(row=i+2,column=1,sticky=tk.W)
            clusterRbList.append(rb)

        def collectInputs():
            dataSelectionFileName = self.PreprocessedCombo.get()
            scaledData = pickle.load(open('outputData/analysisFiles/scaledData/'+dataSelectionFileName,'rb'))
            dataSubsetTitle = dataSelectionFileName.split('-scaledBy')[0]
            if self.DimRedCombo.get() != 'new':
                reductionFileName = self.DimRedCombo.get()
                reducedData = pickle.load(open('outputData/analysisFiles/reducedData/'+reductionFileName,'rb'))
                master.switch_frame(InteractiveClusteringPage,scaledData,reducedData,dataSubsetTitle,v3.get())
            else:
                master.switch_frame(DimensionReductionHomePage,folderName,backpage,secondaryhomepage)

        buttonWindow = tk.Frame(self)
        buttonWindow.pack(side=tk.TOP,pady=10)

        tk.Button(buttonWindow, text="OK",command=lambda: collectInputs()).grid(row=5,column=0)
        tk.Button(buttonWindow, text="Back",command=lambda: master.switch_frame(backpage,folderName,secondaryhomepage)).grid(row=5,column=1)
        tk.Button(buttonWindow, text="Quit",command=quit).grid(row=5,column=2)

class InteractiveClusteringPage(tk.Frame):
    def __init__(self, master,scaledData,reducedData,dataSubsetTitle,clusteringMethod):
        tk.Frame.__init__(self, master)
        
        loff = -0.2

        #Initialize 2x1 canvas for interactive plots
        plotFrame = tk.Frame(self)
        plotFrame.grid(row=0,column=0,columnspan=2)
        fig = plt.figure(figsize=(15, 6))
        gs = fig.add_gridspec(1, 2)
        gs.update(wspace=0.3)
        fig.subplots_adjust(left=0.2)
        levelPlotAxis = fig.add_subplot(gs[0])
        clusterPlotAxis = fig.add_subplot(gs[1])
        self.canvas = FigureCanvasTkAgg(fig,master=plotFrame)
        self.canvas.draw()
        self.canvas.get_tk_widget().pack()
        
        #Dimensional reduction plot (can be colored/resized/restyled by different level values in dropdowns)
        levelPlotWindow = tk.Frame(self)
        levelPlotWindow.grid(row=1,column=0,sticky=tk.N)
        levelParameterList = ['hue','style','size']
        levelParameterValueDict = {}
        levelList = list(scaledData.index.names)+['None']
        featureList = list(scaledData.columns)
        for level in levelParameterList:
            if level == 'hue' or level == 'size':
                levelParameterValueDict[level] = levelList.copy()+featureList.copy()
            else:
                levelParameterValueDict[level] = levelList.copy()
        plottingDfReducedForLegend = reducedData.reset_index()
        kwargs,defaultDict = ipe.getDefaultKwargs(plottingDfReducedForLegend)
        plottingDfReduced = pd.concat([reducedData,scaledData],axis=1).reset_index()
        dropdownList,dropdownVarsDict = ipe.createParameterSelectionDropdowns(levelPlotWindow,levelParameterList,levelParameterValueDict,defaultDict)
        ipe.updateDropdownControlledPlot(self.canvas,levelPlotAxis,plottingDfReduced,dropdownVarsDict,'Dimension 1','Dimension 2',legendoffset=loff)
        
        self.originalxlims = levelPlotAxis.get_xlim()
        self.originalylims = levelPlotAxis.get_ylim()
        self.currentxlims = levelPlotAxis.get_xlim()
        self.currentylims = levelPlotAxis.get_ylim()
        
        #Click and drag widget
        def line_select_callback(eclick, erelease):
            'eclick and erelease are the press and release events'
            x1, y1 = eclick.xdata, eclick.ydata
            x2, y2 = erelease.xdata, erelease.ydata

        def toggle_selector(event):
            print(' Key pressed.')
            if event.key in ['Q', 'q'] and toggle_selector.RS.active:
                print(' RectangleSelector deactivated.')
                toggle_selector.RS.set_active(False)
            if event.key in ['A', 'a'] and not toggle_selector.RS.active:
                print(' RectangleSelector activated.')
                toggle_selector.RS.set_active(True)
        rectpropsdict = {'facecolor':'grey','alpha':0.2,'edgecolor':'grey'}
        toggle_selector.RS = RectangleSelector(levelPlotAxis, line_select_callback,drawtype='box', useblit=True,button=[1, 3], minspanx=5, minspany=5,spancoords='pixels',interactive=True,rectprops=rectpropsdict)
        self.ts = toggle_selector.RS
        self.canvas.mpl_connect('key_press_event', toggle_selector)

        def zoomIn():
            clusterSelectionBox = toggle_selector.RS.corners
            ll = np.array([clusterSelectionBox[0][0], clusterSelectionBox[1][0]])  # lower-left
            ur = np.array([clusterSelectionBox[0][2], clusterSelectionBox[1][2]])  # upper-right
            
            inidx = np.all(np.logical_and(ll <= reducedData.values, reducedData.values <= ur), axis=1)
            inbox = reducedData.loc[inidx]
            bufferval = 0.1
            xlims = [min(inbox['Dimension 1'])-bufferval,max(inbox['Dimension 1'])+bufferval]
            ylims = [min(inbox['Dimension 2'])-bufferval,max(inbox['Dimension 2'])+bufferval]
            self.currentxlims = xlims
            self.currentylims = ylims
            levelPlotAxis.set_xlim(xlims)
            levelPlotAxis.set_ylim(ylims)
            clusterPlotAxis.set_xlim(self.currentxlims)
            clusterPlotAxis.set_ylim(self.currentylims)
            self.canvas.draw()
        
        def zoomOut():
            levelPlotAxis.set_xlim(self.originalxlims)
            levelPlotAxis.set_ylim(self.originalylims)
            clusterPlotAxis.set_xlim(self.originalxlims)
            clusterPlotAxis.set_ylim(self.originalylims)
            self.currentxlims = self.originalxlims 
            self.currentylims = self.originalylims
            self.canvas.draw()
        
        def update():
            ipe.updateDropdownControlledPlot(self.canvas,levelPlotAxis,plottingDfReduced,dropdownVarsDict,'Dimension 1','Dimension 2',legendoffset=loff)
            levelPlotAxis.set_xlim(self.currentxlims)
            levelPlotAxis.set_ylim(self.currentylims)
            clusterPlotAxis.set_xlim(self.currentxlims)
            clusterPlotAxis.set_ylim(self.currentylims)
            self.canvas.draw()

        tk.Button(levelPlotWindow, text="Update level plot",command=lambda: update()).grid(row=3,column=0,columnspan = 2)
        tk.Button(levelPlotWindow, text="Zoom in",command=lambda: zoomIn()).grid(row=4,column=0,columnspan = 2)
        tk.Button(levelPlotWindow, text="Zoom out",command=lambda: zoomOut()).grid(row=5,column=0,columnspan = 2)
        
        #Clustering plot (dimensional reduction is recolored based on clustering parameters selected from sliders)
        clusterParameterWindow = tk.Frame(self)
        clusterParameterWindow.grid(row=1,column=1,sticky=tk.N)
        sliderList = ipe.createParameterAdjustmentSliders(clusterParameterWindow,ods.clusterParameterDict[clusteringMethod],ods.clusterParameterBoundsDict)
        def updateClusterPlot(sliders):
            clusterPlotAxis.clear()
            parametersForClusteringFunction = ipe.getSliderValues(sliderList,ods.clusterParameterDict[clusteringMethod])
            self.currentClusteringParameters = parametersForClusteringFunction.copy()
            self.clusterdf = ods.clusterData(scaledData,clusteringMethod,parametersForClusteringFunction) 
            self.clusterdf.columns.name = 'Feature'
            reducedDataWithClusters = reducedData.copy()
            reducedDataWithClusters['Cluster'] = list(self.clusterdf.index.get_level_values('Cluster'))
            plottingDfClustered = reducedDataWithClusters.reset_index()
            clusterPalette = sns.color_palette(sns.color_palette(),len(pd.unique(reducedDataWithClusters['Cluster'])))
            g1 = sns.scatterplot(data=plottingDfClustered,x='Dimension 1',y='Dimension 2',s=3,ax=clusterPlotAxis,alpha=0.7,hue='Cluster',palette=clusterPalette)
            clusterPlotAxis.legend_.remove()
            clusterPlotAxis.set_xlim(self.currentxlims)
            clusterPlotAxis.set_ylim(self.currentylims)
            self.canvas.draw()
        updateClusterPlot(sliderList)
        tk.Button(clusterParameterWindow, text="Update cluster plot",command=lambda: updateClusterPlot(sliderList)).grid(row=2,column=0)
        def exportDataFrames():
            ods.savePostProcessedFile(self.clusterdf,dataSubsetTitle,'cluster',clusteringMethod,self.currentClusteringParameters)
            df2 = self.clusterdf.copy() 
            df3 = df2.groupby(list(df2.index.names)[-1]).mean()
            print(df2)
            if self.clusterdf.copy().values.max() > 100:
                mfiTicks = [-1000,100,1000,10000,100000]
                mfiTickValues,mfiTickLabels = returnTicks(mfiTicks)
                cg = sns.clustermap(df3.T,cbar_kws={'label':'MFI','ticks':mfiTickValues})
                cg.cax.set_yticklabels(mfiTickLabels)
            else:
                cg = sns.clustermap(df3.T,cbar_kws={'label':'Metric'})
            plt.setp(cg.ax_heatmap.xaxis.get_majorticklabels(), rotation=0)
            plt.setp(cg.ax_heatmap.yaxis.get_majorticklabels(), rotation=0)
            clustermapName = ods.getFileName(dataSubsetTitle,'cluster',clusteringMethod,self.currentClusteringParameters,fileExtension = '.png')[0]
            plt.savefig('plots/clustermap-'+clustermapName,bbox_inches='tight')
            plt.clf()
            print('Clustered Data Frame And Phenotype Plot Saved')
        tk.Button(clusterParameterWindow, text="Save Cluster",command=lambda: exportDataFrames()).grid(row=2,column=1)

        def okCommand():
            exportDataFrames()
            master.switch_frame(backpage,folderName,secondaryhomepage)

        #Default save and quit buttons
        buttonWindow = tk.Frame(self)
        buttonWindow.grid(row=2,column=0,columnspan=2)
        tk.Button(buttonWindow, text="OK",command=lambda: okCommand()).grid(row=0,column=0)
        tk.Button(buttonWindow, text="Back",command=lambda: master.switch_frame(ClusteringHomePage,folderName,backpage,secondaryhomepage)).grid(row=0,column=1)
        tk.Button(buttonWindow, text="Quit",command=lambda: quit()).grid(row=0,column=2)
