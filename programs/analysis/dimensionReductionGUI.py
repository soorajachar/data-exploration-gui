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
#import clusterPlottingLibrary as cpl
import itertools
from dataFrameValueSelectionGUI import DataSelectionHomePage
import operateOnDataSelection as ods
import interactiveGUIElements as ipe

idx = pd.IndexSlice

class DimensionReductionHomePage(tk.Frame):
    num_args = 2
    def __init__(self, master,fName,bp,shp):
        global folderName,backpage,secondaryhomepage
        folderName = fName
        backpage = bp
        secondaryhomepage = shp
        tk.Frame.__init__(self, master)
        
        mainWindow = tk.Frame(self)
        mainWindow.pack(side=tk.TOP,padx=10,pady=10)
        
        self.preprocessedDict = {}
        def getUpdateData(event):
            self.DimRedCombo['values'] = self.preprocessedDict[self.PreprocessedCombo.get()]

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
        if len(self.PreprocessedCombo['values']) == 1:
            self.PreprocessedCombo.set(self.PreprocessedCombo['values'][0])
    
        l2 = tk.Label(mainWindow, text="Dimensional Reduction Type: ").grid(row=1,column=0,sticky=tk.W,pady=(0,10))
        v2 = tk.StringVar()
        v2.set('umap')
        rb2a = tk.Radiobutton(mainWindow,text="umap",padx = 20, variable=v2, value='umap')
        rb2b = tk.Radiobutton(mainWindow,text="tsne",padx = 20, variable=v2, value='tsne')
        rb2d = tk.Radiobutton(mainWindow,text="isomap",padx = 20, variable=v2, value='isomap')
        rb2e = tk.Radiobutton(mainWindow,text="pca",padx = 20, variable=v2, value='pca')
        rb2a.grid(row=1,column=1,sticky=tk.W)
        rb2b.grid(row=2,column=1,sticky=tk.W)
        rb2d.grid(row=3,column=1,sticky=tk.W)
        rb2e.grid(row=4,column=1,sticky=tk.W)
        
        """
        l1 = tk.Label(mainWindow, text="Select Preprocessed Data Subset: ").grid(row=0,column=0,sticky=tk.W)
        
        dropVar = tk.StringVar()
        fileList = returnSpecificExtensionFiles('postProcessedData/scaledData','.pkl',True)
        dropMenu = tk.OptionMenu(mainWindow,dropVar,*fileList)
        dropVar.set(fileList[0])
        setMaxWidth(fileList,dropMenu)
        dropMenu.grid(row=1,column=0,sticky=tk.W)
        
        """

        l3 = tk.Label(mainWindow, text="""Action: """).grid(row=6,column=0,sticky=tk.W,pady=(0,10))
        v3 = tk.StringVar()
        v3.set('create')
        rb3a = tk.Radiobutton(mainWindow,text="create",padx = 20, variable=v3, value='create')
        rb3b = tk.Radiobutton(mainWindow,text="plot",padx = 20, variable=v3, value='plot')
        rb3a.grid(row=6,column=1,sticky=tk.W)
        rb3b.grid(row=7,column=1,sticky=tk.W)
        
        self.DimRedCombo = tkinter.ttk.Combobox(mainWindow,state='readonly')
        self.DimRedCombo['width'] = maxDfLen
        self.DimRedCombo.grid(row = 7,column = 2)
        
        def collectInputs():
            dataSelectionFileName = self.PreprocessedCombo.get()
            print('outputData/analysisFiles/scaledData/'+dataSelectionFileName)
            scaledData = pickle.load(open('outputData/analysisFiles/scaledData/'+dataSelectionFileName,'rb'))
            dataSubsetTitle = dataSelectionFileName.split('-scaledBy')[0]
            if v3.get() == 'create':
                master.switch_frame(InteractiveDimensionReductionPage,scaledData,dataSubsetTitle,v2.get(),[])
            else:
                reductionFileName = self.DimRedCombo.get()
                reducedData = pickle.load(open('outputData/analysisFiles/reducedData/'+reductionFileName,'rb'))
                master.switch_frame(InteractiveDimensionReductionPage,scaledData,dataSubsetTitle,v2.get(),reducedData)

        buttonWindow = tk.Frame(self)
        buttonWindow.pack(side=tk.TOP,pady=10)

        tk.Button(buttonWindow, text="OK",command=lambda: collectInputs()).grid(row=5,column=0)
        tk.Button(buttonWindow, text="Back",command=lambda: master.switch_frame(backpage,folderName,secondaryhomepage)).grid(row=5,column=1)
        tk.Button(buttonWindow, text="Quit",command=quit).grid(row=5,column=2)

class InteractiveDimensionReductionPage(tk.Frame):
    def __init__(self, master,scaledData,dataSubsetTitle,dimRedType,plottingReduction):
        tk.Frame.__init__(self, master)
        
        #Initialize 1x1 canvas for interactive plots
        plotFrame = tk.Frame(self)
        plotFrame.grid(row=0,column=0,columnspan=3)
        fig = plt.figure(figsize=(10, 7))
        gs = fig.add_gridspec(1, 1)
        fig.subplots_adjust(left=0.25)
        offset = -0.1
        levelPlotAxis = fig.add_subplot(gs[0])
        self.canvas = FigureCanvasTkAgg(fig,master=plotFrame)
        self.canvas.draw()
        self.canvas.get_tk_widget().pack()
        
        #Dimensional reduction plot (can be colored/resized/restyled by different level values in dropdowns)
        #Level parameter controls
        levelPlotWindow = tk.Frame(self)
        levelPlotWindow.grid(row=1,column=0,sticky=tk.N)
        levelParameterList = ['hue','style','size']
        levelParameterValueDict = {}
        #Dimensional reduction parameter controls
        dimRedNumericParameterWindow = tk.Frame(self)
        dimRedNumericParameterWindow.grid(row=1,column=1,sticky=tk.N)
        sliderList = ipe.createParameterAdjustmentSliders(dimRedNumericParameterWindow,ods.dimReductionNumericParameterDict[dimRedType],ods.dimReductionNumericParameterBounds)
        dimRedQualitativeParameterWindow = tk.Frame(self)
        dimRedQualitativeParameterWindow.grid(row=1,column=2,sticky=tk.N)
        radiobuttonList,radiobuttonVars = ipe.createParameterSelectionRadiobuttons(dimRedQualitativeParameterWindow,ods.dimReductionQualitativeParameterDict[dimRedType],ods.dimReductionQualitativeParameterValues)
        #Plot
        numericParametersForDimensionReduction = ipe.getSliderValues(sliderList,ods.dimReductionNumericParameterDict[dimRedType])
        qualitativeParametersForDimensionReduction = ipe.getRadiobuttonValues(radiobuttonVars)
        self.currentReductionParameters = numericParametersForDimensionReduction.copy()
        self.currentReductionParameters.update(qualitativeParametersForDimensionReduction.copy())
        if len(plottingReduction) == 0:
            self.reducedData = ods.reduceDimensions(scaledData,dimRedType,[],self.currentReductionParameters) 
        else:
            self.reducedData = plottingReduction.copy()
        self.reducedDataWithFeatures = pd.concat([self.reducedData,scaledData],axis=1)
        levelList = list(self.reducedData.index.names)+['None']
        featureList = list(scaledData.columns)
        for level in levelParameterList:
            if level == 'hue' or level == 'size':
                levelParameterValueDict[level] = levelList.copy()+featureList.copy()
            else:
                levelParameterValueDict[level] = levelList.copy()
        plottingDfReduced = self.reducedData.reset_index()
        plottingDfReducedWithFeatures = self.reducedDataWithFeatures.reset_index()
        kwargs,defaultDict = ipe.getDefaultKwargs(plottingDfReduced)
        dropdownList,dropdownVarsDict = ipe.createParameterSelectionDropdowns(levelPlotWindow,levelParameterList,levelParameterValueDict,defaultDict)
        ipe.updateDropdownControlledPlot(self.canvas,levelPlotAxis,plottingDfReducedWithFeatures,dropdownVarsDict,'Dimension 1','Dimension 2',legendoffset=offset)
        
        self.originalxlims = levelPlotAxis.get_xlim()
        self.originalylims = levelPlotAxis.get_ylim()
        self.currentxlims = levelPlotAxis.get_xlim()
        self.currentylims = levelPlotAxis.get_ylim()

        def updateDimRedPlot(sliderList,radiobuttonVars):
            levelPlotAxis.clear()
            numericParametersForDimensionReduction = ipe.getSliderValues(sliderList,ods.dimReductionNumericParameterDict[dimRedType])
            qualitativeParametersForDimensionReduction = ipe.getRadiobuttonValues(radiobuttonVars)
            self.currentReductionParameters = numericParametersForDimensionReduction.copy()
            self.currentReductionParameters.update(qualitativeParametersForDimensionReduction.copy())
            self.reducedData = ods.reduceDimensions(scaledData,dimRedType,[],self.currentReductionParameters)
            plottingDfReduced = self.reducedData.reset_index()
            self.reducedDataWithFeatures = pd.concat([self.reducedData,scaledData],axis=1)
            plottingDfReducedWithFeatures = self.reducedDataWithFeatures.reset_index()
            ipe.updateDropdownControlledPlot(self.canvas,levelPlotAxis,plottingDfReducedWithFeatures,dropdownVarsDict,'Dimension 1','Dimension 2',legendoffset=offset)
            levelPlotAxis.set_xlim(self.currentxlims)
            levelPlotAxis.set_ylim(self.currentylims)
            self.canvas.draw()
        
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
            
            inidx = np.all(np.logical_and(ll <= self.reducedData.values, self.reducedData.values <= ur), axis=1)
            inbox = self.reducedData.loc[inidx]
            bufferval = 0.1
            xlims = [min(inbox['Dimension 1'])-bufferval,max(inbox['Dimension 1'])+bufferval]
            ylims = [min(inbox['Dimension 2'])-bufferval,max(inbox['Dimension 2'])+bufferval]
            self.currentxlims = xlims
            self.currentylims = ylims
            levelPlotAxis.set_xlim(xlims)
            levelPlotAxis.set_ylim(ylims)
            self.canvas.draw()
        
        def zoomOut():
            levelPlotAxis.set_xlim(self.originalxlims)
            levelPlotAxis.set_ylim(self.originalylims)
            self.currentxlims = self.originalxlims 
            self.currentylims = self.originalylims
            self.canvas.draw()
        
        def update():
            ipe.updateDropdownControlledPlot(self.canvas,levelPlotAxis,pd.concat([self.reducedData,scaledData],axis=1).reset_index(),dropdownVarsDict,'Dimension 1','Dimension 2') 
            levelPlotAxis.set_xlim(self.currentxlims)
            levelPlotAxis.set_ylim(self.currentylims)
            self.canvas.draw()

        def exportDimRed():
            ods.savePostProcessedFile(self.reducedData,dataSubsetTitle,'reduce',dimRedType,self.currentReductionParameters)
            print('Dimensional Reduction Saved!')
            if plotAllDimReds.get():
                if 'Event' in plottingDfReducedWithFeatures.columns or 'event' in plottingDfReducedWithFeatures.columns:
                    sizeParam = {'s':5}
                else:
                    sizeParam = {}
                for feature in pd.unique(scaledData.columns):
                    g = sns.relplot(data=plottingDfReducedWithFeatures,x='Dimension 1',y='Dimension 2',hue=feature,palette='coolwarm',alpha=0.7,**sizeParam)
                    leg = g._legend
                    if max(plottingDfReducedWithFeatures[feature]) > 100 and len(sizeParam) > 0:
                        a,b = returnTicks([-1000,100,10000,100000])
                        for t, l in zip(leg.texts[1:],(b)):
                            t.set_text(l)
                    else:
                        for t in leg.texts:
                            # truncate label text to 4 characters
                            if t.get_text() == '1.2000000000000002':
                                t.set_text('1.0')
                            else:
                                if '.' in t.get_text():
                                    if t.get_text().replace('.','',1).isdigit():
                                        decimalIndex = t.get_text().find('.')
                                        t.set_text(round(float(t.get_text()),2))
                    reducedName = ods.getFileName(dataSubsetTitle,'reduce',dimRedType,self.currentReductionParameters,fileExtension = '.png')[0]
                    subprocess.run(['mkdir','plots/'+reducedName[:-4]+'-featureColoredPlots'])
                    plt.savefig('plots/'+reducedName[:-4]+'-featureColoredPlots/'+str(feature)+'.png',bbox_inches='tight')
                    plt.clf()
                    print(feature+' plot saved')

        tk.Button(self, text="Update plot styling",command=lambda: update()).grid(row=2,column=0)
        tk.Button(self, text="Zoom in",command=lambda: zoomIn()).grid(row=4,column=0)
        tk.Button(self, text="Zoom out",command=lambda: zoomOut()).grid(row=5,column=0)
        tk.Button(self, text="Update hyperparameters",command=lambda: updateDimRedPlot(sliderList,radiobuttonVars)).grid(row=2,column=1,columnspan=2)

        tk.Button(self,text='Save plot',command=lambda: fig.savefig('plots/'+ods.getFileName(dataSubsetTitle,'reduce',dimRedType,self.currentReductionParameters,plotParameterDict=ipe.getDropdownValues(dropdownVarsDict),fileExtension='.png')[0],bbox_inches='tight')).grid(row=3,column=0)
        tk.Button(self,text='Save dimensional reduction',command=lambda: exportDimRed(),font='Helvetica 14 bold').grid(row=3,column=1,columnspan=2)
        plotAllDimReds = tk.BooleanVar()
        plotAllDimRedsButton = tk.Checkbutton(self,text='Save all feature-colored dimension reduction plots?',variable=plotAllDimReds)
        plotAllDimRedsButton.select()
        plotAllDimRedsButton.grid(row=4,column=1,columnspan=2)
    
        def okCommand():
            exportDimRed()
            master.switch_frame(backpage,folderName,secondaryhomepage)

        #Default save and quit buttons
        buttonWindow = tk.Frame(self)
        buttonWindow.grid(row=6,column=0,columnspan=2)
        try:
            k=backpage
        except NameError: 
            backpage,folderName,secondaryhomepage = pickle.load(open('misc/dimRedPlottingParamList.pkl','rb'))
            
        tk.Button(buttonWindow, text="OK",command=lambda: master.switch_frame(backpage,folderName,secondaryhomepage)).grid(row=0,column=0)
        tk.Button(buttonWindow, text="Back",command=lambda: master.switch_frame(DimensionReductionHomePage,folderName,backpage,secondaryhomepage)).grid(row=0,column=1)
        tk.Button(buttonWindow, text="Quit",command=lambda: quit()).grid(row=0,column=2)

