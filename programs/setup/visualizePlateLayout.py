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
import itertools
from matplotlib import pyplot as plt
from scipy import stats
from sklearn.manifold import TSNE
from sklearn.preprocessing import MinMaxScaler
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.widgets import RectangleSelector
from createPlateLayout import returnBaseLayout  
import itertools
sys.path.insert(0, '../dataprocessing/')
from miscFunctions import sortSINumerically,reindexDataFrame,setMaxWidth,returnSpecificExtensionFiles,returnTicks,get_cluster_centroids,rainbow_text

idx = pd.IndexSlice
dataTypeObservableRangeDict = {'cyt':1,'cell':3,'prolif':1}
realDataTypeNameDict = {'cyt':'Supernatant','cell':'Surface/Intracellular Marker','prolif':'Proliferation'}
plateRowLetters = string.ascii_uppercase[:16]
plateColumnNumbers = list(range(1,25))

#For macbook by itself
#figLengthScaling = 0.5
#For work monitor
figLengthScaling = 1
#For home monitor
#figLengthScaling = 0.75

def createLayoutVisual(baseLayoutDf,currentLayout,levelIndex,currentLevel,levelValues,plateDimensions,numRowPlates,numColumnPlates):
    
        fig = plt.figure(figsize=(3*numColumnPlates*(plateDimensions[0]/8), 2.5*numRowPlates*(plateDimensions[1]/12)*figLengthScaling),tight_layout=True)
        #fig_ax1 = fig.add_subplot(111)
        gs = fig.add_gridspec(1, 1)
        fig_ax1 = fig.add_subplot(gs[0])
        
        #Remove all axis elements
        fig_ax1.set_xlim((0, max(baseLayoutDf['x'])+1))
        fig_ax1.set_ylim((0, max(baseLayoutDf['y'])+1))
        fig_ax1.set_xticks([], [])
        fig_ax1.set_yticks([], [])
        fig_ax1.set_xlabel('')
        fig_ax1.set_ylabel('')
        #Add plate dividing lines
        for vlinex in vlinelist:
            fig_ax1.axvline(vlinex,linestyle=':',color='k')
        for hliney in hlinelist:
            fig_ax1.axhline(hliney,linestyle=':',color='k')
        #Add well IDs
        for row in range(baseLayoutDf.shape[0]):
            fig_ax1.annotate(infoDf.values[row][0],(baseLayoutDf.iloc[row,0],baseLayoutDf.iloc[row,1]),ha='center',va='center',size=7)
        #Add plateIDs
        for row in range(baseLayoutDf.shape[0]):
            if infoDf.values[row][1] != 'DoNotLabel':
                fig_ax1.annotate(infoDf.values[row][1],(baseLayoutDf.iloc[row,0]-0.6,baseLayoutDf.iloc[row,1]+0.5),size=7,ha='center',va='center')
        
        currentpalette = ['#808080']+sns.color_palette('husl',len(levelValues[levelIndex])).as_hex()
        newLayoutDf = baseLayoutDf.copy()
        newLayoutDf['key'] = currentLayout.flatten()
            
        g1 = sns.scatterplot(data=baseLayoutDf,x='x',y='y',ax=fig_ax1,color='#ffffff',s=200,marker='o')
        if -1 in list(newLayoutDf['key']):
            hueorder = [-1]
            modifiedPalette = ['#808080']
        else:
            hueorder = []
            modifiedPalette = []
        for i in range(len(levelValues[levelIndex])):
            if i in list(pd.unique(newLayoutDf['key'])):
                hueorder.append(i)
                modifiedPalette.append(currentpalette[i+1])
        g1 = sns.scatterplot(data=newLayoutDf,x='x',y='y',ax=fig_ax1,hue='key',hue_order=hueorder,palette=modifiedPalette,s=200,markers=['X','o'],alpha=0.5,style='blank',style_order=[-1,0])
        fig_ax1.legend_.remove()
        titlelabels = [currentLevel+': ']
        titlecolors = ['black']
        print(levelValues)
        for i,levelValue in enumerate(levelValues[levelIndex]):
            titlelabels.append(str(levelValue))
            titlecolors.append(modifiedPalette[i])
        rainbow_text(2,max(baseLayoutDf['y'])+2, titlelabels,titlecolors,size=14,ax=fig_ax1)
        #plt.savefig('plateLayout-'+currentLevel+'.png',bbox_inches='tight')
        subprocess.run(['mkdir','plots/plateLayouts'])
        plt.savefig('plots/plateLayouts/plateLayout-'+currentLevel+'.png',bbox_inches='tight')
        plt.clf()

root = '../../../allExperiments/QualityQuantityDeconvolution/20200120-TumorTimeseries_2/misc/'
folderName = '20200120-TumorTimeseries_2'

experimentParameters = json.load(open(root+'experimentParameters-'+folderName+'-cyt.json','r'))
numRowPlates = experimentParameters['numRowPlates']
numColumnPlates = experimentParameters['numColumnPlates']
plateDimensions = experimentParameters['overallPlateDimensions']
levelValues = []
for level in experimentParameters['allLevelValues']:
    levelValues.append(experimentParameters['allLevelValues'][level])
newLevelValues = [levelValues[-1]]
levelValues = newLevelValues+levelValues[:-1]

layoutDict = pickle.load(open(root+'layoutDict-'+folderName+'-cyt.pkl','rb'))
baseLayoutDf,infoDf,vlinelist,hlinelist = returnBaseLayout(plateDimensions,numRowPlates,numColumnPlates)
for i in layoutDict['keys']:
    level = experimentParameters['allLevelNames'][i]
    createLayoutVisual(baseLayoutDf,layoutDict['keys'][i],i,level,levelValues,plateDimensions,numRowPlates,numColumnPlates)
