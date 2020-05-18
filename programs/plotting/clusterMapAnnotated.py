#!/usr/bin/env python3  
# -*- coding: utf-8 -*-
"""
created on sat jul 21 13:12:56 2018

@author: acharsr
"""
import pickle
import os
import sys,subprocess
from sklearn import preprocessing
import matplotlib
import numpy as np
import seaborn as sns
import pandas as pd
from matplotlib import pyplot as plt
from itertools import groupby
import sklearn.metrics as skm
import scipy.cluster.hierarchy as shc
from operator import itemgetter
import matplotlib.gridspec as gridspec
import scipy
import collections 
from matplotlib.patches import Rectangle
from sklearn.preprocessing import MinMaxScaler,MaxAbsScaler,StandardScaler

#Adds vertical lines starting at xpos,ypos of length length; positive length values cause the line to go up from starting point while negatives cause it to go down
def add_vline(ax, xpos, ypos,length):
    line = plt.Line2D([xpos, xpos], [ypos, ypos+length],transform=ax.transAxes, color='black',linewidth=0.5)
    line.set_clip_on(False)
    ax.add_line(line)

#Adds horizontal lines starting at xpos,ypos of length length; positive length values cause the line to go to the right from starting point while negatives cause it to go to the left
def add_hline(ax, xpos, ypos,length):
    line = plt.Line2D([xpos, xpos+length], [ypos, ypos],transform=ax.transAxes, color='black',linewidth=0.5)
    line.set_clip_on(False)
    ax.add_line(line)

#Adds borders around heatmap
def add_borders(ax):
    add_hline(ax,0,0,1)
    add_hline(ax,0,1,1)
    add_vline(ax,0,0,1)
    add_vline(ax,1,0,1)

def label_len(my_index,level):
    labels = my_index.get_level_values(level)
    return [(k, sum(1 for i in g)) for k,g in groupby(labels)]

#Creates row labels (not row levels), starting at position given by xstart
def label_rows(ax, df,xstart,labelboxwidth,fontsize=10,original_df=[]):
    xpos = xstart
    scale = 1./df.index.size
    for level in range(len(list(df.index.names)))[::-1]:
        levelName = list(df.index.names)[level]
        pos = len(df.index)
        palette = sns.color_palette(sns.color_palette(),len(pd.unique(df.index.get_level_values(level))))
        lut={}
        if not isinstance(original_df,list):
            for i,label in enumerate(pd.unique(original_df.index.get_level_values(levelName))):
                lut[label] = palette[i]
        else:
            for i,label in enumerate(pd.unique(df.index.get_level_values(level))):
                lut[label] = palette[i]
        for label, rpos in label_len(df.index,level):
            lypos = (pos - .6 * rpos)*scale
            t = ax.text(xpos+labelboxwidth/2, lypos, label, ha='center', va='center',transform=ax.transAxes,size=fontsize)
            add_hline(ax, xpos, pos*scale,labelboxwidth)
            ax.add_patch(Rectangle((xpos, pos*scale), labelboxwidth, -2*(pos*scale-(pos-0.5*rpos)*scale), facecolor=lut[label],transform=ax.transAxes,clip_on=False,alpha=0.7))
            #ax.add_patch(Rectangle((xpos, pos*scale), labelboxwidth, (lypos-pos*scale)*2, facecolor=lut[label],transform=ax.transAxes,clip_on=False,alpha=0.7))
            pos -= rpos
        add_hline(ax, xpos, pos*scale,labelboxwidth)
        add_vline(ax,xpos+labelboxwidth,0,1)
        add_vline(ax,xstart,0,1)
        xpos += labelboxwidth
    add_hline(ax,xstart,1,df.index.nlevels*labelboxwidth)

#Creates column labels (not column levels) starting at position given by ystart
def label_columns(ax,df,ystart,labelboxheight,fontsize=10):
    ypos = ystart
    scale = 1./df.columns.size
    for level in range(len(list(df.columns.names)))[::-1]:
        pos = 0
        for label, rpos in label_len(df.columns,level):
            lxpos = (pos + .5 * rpos)*scale
            t = ax.text(lxpos, ypos+labelboxheight/2, label, ha='center', va='center',transform=ax.transAxes,size=fontsize)
            add_vline(ax, pos*scale, ypos, labelboxheight)
            pos += rpos
        add_vline(ax,pos*scale,ypos,labelboxheight)
        add_hline(ax,0,ypos+labelboxheight,1)
        ypos += labelboxheight
    #Adds beginning vline
    add_vline(ax,0,ystart,df.columns.nlevels*labelboxheight)

#Creates row level labels starting at xstart,ystart
def label_row_headers(ax,df,xstart,ystart,labelboxwidth,labelboxheight,fontsize=10):
    scale = len(df.index.names)
    xpos = xstart
    for name in df.index.names[::-1]:
        add_hline(ax,xpos,ystart-labelboxheight,labelboxwidth)
        add_vline(ax,xpos,ystart-labelboxheight,labelboxheight)
        ax.text(xpos+labelboxwidth/2,ystart-labelboxheight/2,name,ha='center',va='center',transform=ax.transAxes,fontweight='bold',size=fontsize)
        xpos+=labelboxwidth
    add_vline(ax,xpos,ystart-labelboxheight,labelboxheight)

#Creates column level labels starting at xstart,ystart
def label_column_headers(ax,df,xstart,ystart,labelboxwidth,labelboxheight,fontsize=10):
    scale = len(df.columns.names)
    ypos = ystart
    for name in df.columns.names[::-1]:
        add_hline(ax,ystart-labelboxwidth/2,ypos,labelboxwidth/2)
        add_vline(ax,ystart-labelboxwidth/2,ypos,labelboxheight)
        ax.text(ystart-labelboxwidth/4,ypos+labelboxheight/2,name,ha='center',va='center',transform=ax.transAxes,fontweight='bold',size=fontsize)
        ypos+=labelboxheight
    add_hline(ax,ystart-labelboxwidth/2,ypos,labelboxwidth)
    
#WILL MAKE FOLDER CALLED CLUSTERMAPS IN YOUR CURRENT DIRECTORY, WHERE THE CLUSTERMAPS YOU GENERATE WILL BE SAVED
def prepareClusterMaps(df=[],outputFileName='default',binSize=1,outputPath='.',fontsize=5,labelboxwidth=0.09,boxheightfactor=4,scalingMethod = 'minmax'):
    if 'clustermaps' not in os.listdir(outputPath):
        subprocess.run(['mkdir','clustermaps'])
    #Bin timepoints together to avoid overlapping column labels
    binLength = binSize
    i=-1*binLength
    timepoints = list(pd.unique(df.columns.get_level_values('Time')))
    newcolumns = []
    for col in range(df.shape[1]):
        if col%binLength == 0:
            i+=binLength
        timeRange = str(str(int(df.columns[i]))).zfill(2)+'-'+str(int(df.columns[i+binLength-1]))
        newcolumns.append(timeRange)
    df.columns = newcolumns
    df.columns.name = 'Time'

    if 'Cytokine' in df.index.names or 'Cytokine' in df.columns.names:
        unstackingList = ['Cytokine']
    else:
        unstackingList = ['CellType','Marker','Statistic']

    df2 = df.unstack(unstackingList).stack('Time')
    
    scalingFunctionDict = {'minmax':MinMaxScaler,'maxabs':MaxAbsScaler,'standard':StandardScaler}
    #Scale with chosen method
    if scalingMethod.lower() != 'none':
        x = df2.values #returns a numpy array
        scaler = scalingFunctionDict[scalingMethod.lower()]  
        x_scaled = scaler().fit_transform(x)
        df2 = pd.DataFrame(x_scaled,index=df2.index,columns=df2.columns)
    
    #Normalized across all times
    df2 = df2.unstack('Time')
    
    linkage = shc.linkage(df2,method='ward')
    if '/' != outputPath[-1]:
        outputPath+='/'
    
    original_df = df.copy()
    df = df2.copy()
    linkage = linkage
    folderName = outputPath+'clustermaps/'+outputFileName+'-annotatedClusterMap.pdf'

    h = 20
    clax = sns.clustermap(df,xticklabels=False,yticklabels=False,row_linkage=linkage,col_cluster=False,cbar_pos=None,dendrogram_ratio=(0.2,0))
    fig1=plt.gcf()
    fig1.dpi=500
    hax = clax.ax_heatmap
    hax.set_xlabel('')
    hax.set_ylabel('')
    
    #CHANGE BOX WIDTHS AND TEXT SIZES HERE
    #For dendrogram with row labels at right of heatmap
    #Change this to change width of boxes used for row labels
    #Change this to change height of boxes used for column labels
    xstart = 1
    ystart = 0

    labelboxheight = -1./df.index.size/boxheightfactor
    
    newindex = clax.dendrogram_row.reordered_ind 
    newdf = pd.DataFrame(df.values,index=df.index[newindex],columns=df.columns)
    label_rows(hax,newdf,xstart,labelboxwidth,fontsize=fontsize,original_df=original_df)
    label_columns(hax,newdf,ystart,labelboxheight,fontsize=fontsize)
    label_row_headers(hax,newdf,xstart,1,labelboxwidth,labelboxheight,fontsize=fontsize)
    label_column_headers(hax,newdf,0,ystart,labelboxwidth*df.index.nlevels,labelboxheight,fontsize=fontsize)
    add_borders(hax)
    fig1.savefig(folderName,bbox_inches='tight') 
    plt.clf()
