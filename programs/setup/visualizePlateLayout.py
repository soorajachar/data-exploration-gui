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

colwrap = 6
#For macbook by itself
#figLengthScaling = 0.5
#For work monitor
figLengthScaling = 1
#For home monitor
#figLengthScaling = 0.75

def createLayoutVisual(baseLayoutDf,currentLayout,levelIndex,currentLevel,levelValues,plateDimensions,numRowPlates,numColumnPlates,dt,infoDf,vlinelist,hlinelist,root=''):
    
    maxTextLen = len(str(currentLevel))
    for levelVal in levelValues[levelIndex]:
        if len(str(levelVal)) > maxTextLen:
            maxTextLen = len(str(levelVal))
    
    if numRowPlates == 1 and numColumnPlates > colwrap:
        colwrapBool = True
    else:
        colwrapBool = False

    if colwrapBool:
        visualNumRowPlates = math.ceil(numColumnPlates/colwrap)
        visualNumColumnPlates = colwrap
    else:
        visualNumRowPlates = numRowPlates
        visualNumColumnPlates = numColumnPlates 
    fig = plt.figure(figsize=(3*visualNumColumnPlates*(plateDimensions[0]/8)+1+(0.125*maxTextLen), 2.5*visualNumRowPlates*(plateDimensions[1]/12)*figLengthScaling),tight_layout=True)
    #fig_ax1 = fig.add_subplot(111)
    gs = fig.add_gridspec(1, 1)
    fig_ax1 = fig.add_subplot(gs[0])
    
    #Remove all axis elements
    fig_ax1.set_xlim((0, max(baseLayoutDf['x'])+1))
    fig_ax1.set_ylim((0, max(baseLayoutDf['y'])+1))
    fig_ax1.set_xticks([])
    fig_ax1.set_yticks([])
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
    print(newLayoutDf.key.unique())
    if -1 in list(newLayoutDf['key']) or 'blank' in list(newLayoutDf['key']):
        hueorder = [-1]
        modifiedPalette = ['#808080']
        ogpalette = modifiedPalette.copy()
    else:
        hueorder = []
        modifiedPalette = []
        ogpalette = modifiedPalette.copy()
    for i in range(len(levelValues[levelIndex])):
        if i in list(pd.unique(newLayoutDf['key'])):
            hueorder.append(i)
            modifiedPalette.append(currentpalette[i+1])
    blanksExist = False 
    for i,val in enumerate(newLayoutDf['key']):
        if val == 'blank':
            newLayoutDf['key'].iloc[i] = -1
            newLayoutDf['blank'].iloc[i] = -1
            blanksExist = True
        else:
            newLayoutDf['blank'].iloc[i] = 0 
    
    g1 = sns.scatterplot(data=newLayoutDf,x='x',y='y',ax=fig_ax1,hue='key',hue_order=hueorder,palette=modifiedPalette,s=200,markers=['X','o'],alpha=0.5,style='blank',style_order=[-1,0])
    titlelabels = [currentLevel+': ']
    titlecolors = ['black']
    
    print(modifiedPalette)
    print(ogpalette)
    print(levelValues[levelIndex])
    for i,levelValue in enumerate(levelValues[levelIndex]):
        titlelabels.append(str(levelValue))
        titlecolors.append(modifiedPalette[i+len(ogpalette)])
    #plt.savefig('plateLayout-'+currentLevel+'.png',bbox_inches='tight')
    if 'plateLayouts' not in os.listdir(root+'plots'):
        subprocess.run(['mkdir',root+'plots/plateLayouts'])
    legendHandlesLabels = fig_ax1.get_legend_handles_labels()
    i=0
    newLegendLabels = []
    newLegendHandles = []
    legendHandles = legendHandlesLabels[0]
    currentLevelValues = levelValues[levelIndex]
    if blanksExist:
        offset = 1
    else:
        offset = 0
    for legendHandle in legendHandles:
        #Skips the style legend handles
        if i < len(currentLevelValues)+1:
            if i == 0:
                modifiedLevel = currentLevel.translate(str.maketrans({"-":  r"\-","]":  r"\]","\\": r"\\","^":  r"\^","$":  r"\$","*":  r"\*",".":  r"\.","_":  r"\_","(":  r"\(",")":  r"\)","[":  r"\[","%": r"\%"}))
                newLegendLabels.append('$\\bf'+modifiedLevel+'$')
                newLegendHandles.append(legendHandle)
            else:
                newLegendLabels.append(currentLevelValues[i-1])
                newLegendHandles.append(legendHandles[i+offset])
        i+=1
    fig_ax1.legend(bbox_to_anchor=(1, 1),frameon=False,handles=newLegendHandles, labels=newLegendLabels)
    plt.savefig(root+'plots/plateLayouts/plateLayout-'+currentLevel+'-'+dt+'.png',bbox_inches='tight')
    plt.clf()

dt = 'cell'
root = '/Users/acharsr/Documents/allExperiments/Naomi-CAR-Collaboration/20200129-Human-CAR-T_3/'
folderName = root.split('/')[-2]
experimentParameters = json.load(open(root+'misc/experimentParameters-'+folderName+'-'+dt+'.json','r'))
numRowPlates = 1
numColumnPlates = 1
plateDimensions = experimentParameters['overallPlateDimensions']
levelValues = []
for level in experimentParameters['levelLabelDict']:
    levelValues.append(experimentParameters['levelLabelDict'][level])

baseLayoutDf,infoDf,vlinelist,hlinelist = returnBaseLayout(plateDimensions,numRowPlates,numColumnPlates)
levels = list(experimentParameters['levelLabelDict'].keys())
layoutDict = pickle.load(open(root+'misc/layoutDict-'+folderName+'-'+dt+'.pkl','rb'))
for i in layoutDict['keys']:
    #print(layoutDict['keys'][i][:,:8])
    #print(layoutDict['keys'][i][:,3:-2])
    level = levels[i]
    createLayoutVisual(baseLayoutDf,layoutDict['keys'][i],i,level,levelValues,plateDimensions,numRowPlates,numColumnPlates,dt,infoDf,vlinelist,hlinelist,root=root)
