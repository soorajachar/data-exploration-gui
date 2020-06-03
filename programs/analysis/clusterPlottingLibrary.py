#!/usr/bin/env python3  
# -*- coding: utf-8 -*-
"""
created on sat jul 21 13:12:56 2018

@author: acharsr
"""
import pickle
import os
import sys
from sys import platform as sys_pf
if sys_pf == 'darwin':
    import matplotlib
    matplotlib.use("TkAgg")
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
from statistics import mean,median
from scipy.stats import shapiro
from statsmodels.multivariate.manova import MANOVA
from statsmodels.stats.multicomp import pairwise_tukeyhsd,MultiComparison
sys.path.insert(0, '../dataprocessing/')
from miscFunctions import sortSINumerically,reindexDataFrame,setMaxWidth

def significanceTesting(featureDf2,pairwiseClustersToCompare,confidence=0.05,foldchange=2,responseCutoff=0.1,errorCorrection='bonferroni'):
    n = len(featureDf2.columns)-1
    if errorCorrection == 'bonferroni':
        alpha = confidence/n
    else:
        alpha = confidence
    uniqueClusters = [list(x) for x in set(tuple(x) for x in pairwiseClustersToCompare)]
    
    #Kruskal Wallis is unecessary; one way anova seems to be relatively robust to non-normality: http://www.biostathandbook.com/kruskalwallis.html
    endog = featureDf2.iloc[:,:-1]
    exog = featureDf2.iloc[:,-1]
    modelFormula = " + ".join("Q(\'"+featureDf2.columns[:-1]+"\')")+" ~ Cluster"
    print(featureDf2)
    sys.exit(0)
    manova = MANOVA.from_formula(modelFormula,data=featureDf2)
    #Pillai's trace is most robust against deviations from assumptions of manova
    manovapval = manova.mv_test().results['Cluster']['stat'].iloc[1,4]
    print(manovapval)
    #Need to think about how to handle multiple clusters; for now just iterate through all pairs
    if manovapval < confidence:
        allDataMatrices = []
        allSignificantDifferences = []
        for clustersToCompare in pairwiseClustersToCompare:
            comp1 = clustersToCompare[0] 
            comp2 = clustersToCompare[1]
            group1 = featureDf2[featureDf2['Cluster'] == str(comp1)].iloc[:,:-1]
            group2 = featureDf2[featureDf2['Cluster'] == str(comp2)].iloc[:,:-1]
            anova = scipy.stats.kruskal(group1,group2)
            pval2 = anova[1]
            stat = anova[0]
            if pval2 < 0.01:
                print('Different')
            significantArray = []
            allBoxPairs = []
            pvalList = []
            meanFoldChangeList = []
            medianFoldChangeList = []
            foldChangeList = []
            normalityList = []
            tempnormalityList = []

            for col in range(featureDf2.shape[1]-1):
                group1 = featureDf2[featureDf2['Cluster'] == str(comp1)].iloc[:,col]
                group2 = featureDf2[featureDf2['Cluster'] == str(comp2)].iloc[:,col]
                normalitypval = shapiro(group1)[1]
                normalitypval2 = shapiro(group2)[1]
                normalityCondition = False
                if normalitypval < 0.05 and normalitypval2 < 0.05:
                    normalityCondition = True
                    try:
                        pval = scipy.stats.ttest_ind(group1,group2)[1]
                    except:
                        pval = 0.5
                else:
                    try:
                        pval = scipy.stats.mannwhitneyu(group1,group2)[1]
                    except:
                        pval = 0.5
                pvalList.append(pval)
                tempnormalityList.append(normalityCondition)
            
            #For holm bonferroni
            ordered_pval_list = sorted(pvalList)

            for col in range(featureDf2.shape[1]-1):
                pvalCondition = False
                foldChangeCondition = False
                group1 = featureDf2[featureDf2['Cluster'] == str(comp1)].iloc[:,col]
                group2 = featureDf2[featureDf2['Cluster'] == str(comp2)].iloc[:,col]
                
                pval = pvalList[col]
                if errorCorrection != 'holm-bonferroni':
                    if pval < alpha:
                        pvalCondition = True
                else:
                    rank = ordered_pval_list.index(pval)+1
                    modifiedAlpha = alpha / (n - rank + 1)
                    if pval < modifiedAlpha:
                        pvalCondition = True

                normalityCondition = tempnormalityList[col]
                if normalityCondition:
                    if np.nanmean(group1) < responseCutoff:
                        if np.nanmean(group2) >= responseCutoff:
                            meanFoldChangeList.append(4)
                            foldChangeList.append(4)
                        else:
                            meanFoldChangeList.append(0.0001)
                            foldChangeList.append(0.0001)
                    else:
                        if np.nanmean(group2) < responseCutoff:
                            meanFoldChangeList.append(4)
                            foldChangeList.append(4)
                        else:
                            meanFoldChangeList.append(np.nanmean(group1)/np.nanmean(group2))
                            foldChangeList.append(np.nanmean(group1)/np.nanmean(group2))
                else:
                    if np.nanmedian(group1) < responseCutoff:
                        if np.nanmedian(group2) >= responseCutoff:
                            medianFoldChangeList.append(4)
                            foldChangeList.append(4)
                        else:
                            medianFoldChangeList.append(0.0001)
                            foldChangeList.append(0.0001)
                    else:
                        if np.nanmedian(group2) < responseCutoff:
                            medianFoldChangeList.append(4)
                            foldChangeList.append(4)
                        else:
                            medianFoldChangeList.append(np.nanmedian(group1)/np.nanmedian(group2))
                            foldChangeList.append(np.nanmedian(group1)/np.nanmedian(group2))
                if pvalCondition:
                    if abs(np.log2(foldChangeList[-1])) >= np.log2(foldchange):
                        significantArray.append(featureDf2.columns.get_level_values('Feature')[col])
                        allBoxPairs.append(((featureDf2.columns.get_level_values('Feature')[col],str(comp1)),(featureDf2.columns.get_level_values('Feature')[col],str(comp2))))        
                        normalityList.append(normalityCondition)
            
            foldChangeArray = np.log2(np.array(foldChangeList))
            pvalArray = -np.log10(np.array(pvalList))
            
            dataMatrix = np.vstack([foldChangeArray,pvalArray])
            allSignificantDifferences.append(significantArray)
            allDataMatrices.append(dataMatrix)
        significantArray = list(set().union(*allSignificantDifferences))
        dataMatrix = np.vstack(allDataMatrices)
    else:
        significantArray = []
        dataMatrix = []
    
    print(significantArray)
    return dataMatrix,significantArray
    
#Volcano plot
def volcanoPlot(featureDf,comparisonDictEntry,dataMatrix,ax1a,xaxistitle,yaxistitle):
    comp1 = comparisonDictEntry['clustersToCompare'][0] 
    comp2 = comparisonDictEntry['clustersToCompare'][1]
    
    volcanoPlotDf = pd.DataFrame(dataMatrix.T,index=featureDf.columns,columns=[xaxistitle,yaxistitle])
    plottingDfVolcano = volcanoPlotDf.reset_index()
    for k, v in volcanoPlotDf.iterrows():
        if v[1] > -1*np.log10(0.01/len(featureDf.columns)) and abs(v[0]) >= 1:
            ax1a.annotate(k,v,xytext=(10,-5),textcoords='offset points',family='sans-serif')
    
    largerClusterList = []
    for i in range(plottingDfVolcano.shape[0]):
        if plottingDfVolcano[xaxistitle][i] >= 1:
            largerClusterList.append(str(comp1))
        elif plottingDfVolcano[xaxistitle][i] <= -1:
            largerClusterList.append(str(comp2))
        else:
            largerClusterList.append('neither')

    plottingDfVolcano['LargerCluster'] = largerClusterList
    return plottingDfVolcano     

def addFeatureMultiIndex(df):
    newcolumnTuples = []
    for column in df.columns:
        if not isinstance(column,tuple):
            fullColumnTuple = ['Cytokines','allCells','Concentration']+[column]
        else:
            column2 = list(column)
            column3 = column2.copy()
            column3[1] = column2[2]
            column3[2] = column2[1]
            if column3[-2] == 'Negative GFI':
                column3[-2] = 'MFI'
            fullColumnTuple = ['Surface Markers']+column3
        newcolumnTuples.append(fullColumnTuple)
    newcolumns = pd.MultiIndex.from_tuples(newcolumnTuples,names=['DataType','Population','Statistic','Feature'])
    newdf = pd.DataFrame(df.values,index=df.index,columns=newcolumns)
    return newdf

def mergeTimepointLabels(df,reindexingDf,timepointLabelMergingList):

    dfToReindex = df.unstack('Time')
    reindexedDf = reindexDataFrame(dfToReindex,reindexingDf,False)
    
    newcolumnsTuples = []
    for col in range(reindexedDf.shape[1]):
        names = list(reindexedDf.iloc[:,col].name)
        for i,timeRange in enumerate(timepointLabelMergingList):
            lowerTimebound = float(timeRange.split('-')[0])
            upperTimebound = float(timeRange.split('-')[1])
            if float(names[-1]) > lowerTimebound and float(names[-1]) <= upperTimebound:
                names[-1] = timeRange
                break
        newcolumnsTuples.append(names)
    newcolumns = pd.MultiIndex.from_tuples(newcolumnsTuples,names=['DataType','Population','Statistic','Feature','Time'])
    newdf = pd.DataFrame(reindexedDf.values,index=reindexedDf.index,columns=newcolumns)
    return newdf

def reorderDataframeByLinkage(df,linkage):
    tempdendrogram = augmented_dendrogram(1,linkage, color_threshold=None,orientation='left',no_labels=False,count_sort=True, no_plot=True) 
    ll = tempdendrogram ['leaves'][::-1]
    temp = df.reset_index().reindex(ll)
    tuplesForReindex = []
    for row in range(temp.shape[0]):
        tuplesForReindex.append(tuple(temp.iloc[row,:len(df.index.names)]))
    newmi = pd.MultiIndex.from_tuples(tuplesForReindex,names=df.index.names)
    newdf = pd.DataFrame(temp.values[:,len(df.index.names):].astype(float),index=newmi,columns=df.columns)
    return newdf

def shrinkFeatureMultiIndex(df,startIndex=1):
    newcolList = []
    if 'Time' not in df.index.names:
        newdf = df.stack('Time')
    #Do not include levels of cell features that only have a single value in the feature titles
    cellFeatureComponentList = [[],[],[]]
    for col in df.columns:
        if isinstance(col,tuple):
            for i,cellFeatureComponent in enumerate(col):
                cellFeatureComponentList[i].append(col[i])
    uniqueCellFeatures = []
    for featureComponentList in cellFeatureComponentList:
        uniqueCellFeatures.append(list(pd.unique(featureComponentList)))
    if 'NotApplicable' in uniqueCellFeatures[1]:
        uniqueCellFeatures[1].remove('NotApplicable')
    cellComponentsToInclude = []
    for i in range(len(uniqueCellFeatures)):
        if len(uniqueCellFeatures[i]) > 1:
            cellComponentsToInclude.append(i)

    for col in df.columns:
        if isinstance(col,tuple):
            includedElements = itemgetter(*cellComponentsToInclude)(col)
            if not isinstance(includedElements,tuple) and not isinstance(includedElements,list):
                includedElements = [includedElements]
            newcol = ','.join(includedElements)
        else:
            newcol = col
        newcolList.append(newcol)
    newdf = df.copy()
    newdf.columns = newcolList
    newdf.columns.name='Feature'
    return newdf

def augmented_dendrogram(breakpointDistance,*args, **kwargs):
    #Draw in dendrogram with annotations
    ddata = shc.dendrogram(*args, **kwargs)
    if not kwargs.get('no_plot', False):
        matrixlist = []
        for i,d in zip(ddata['icoord'],ddata['dcoord']):
            x = 0.5 * sum(i[1:3])
            y = d[1]
            matrixlist.append([x,y])
        matrix = np.matrix(matrixlist)
        tempdf = pd.DataFrame(matrix,columns=['y','x'])
        sortedDf = tempdf.sort_values('x',ascending=True)
        for ind in range(sortedDf.shape[0]):
            x = sortedDf.iloc[ind,1] 
            y = sortedDf.iloc[ind,0]
            if x == breakpointDistance: 
                kwargs['ax'].plot(x, y, 'k',marker='o')
                #ax2.annotate(labels[ind], (x, y), xytext=(0,-5),textcoords='offset points',va='top', ha='center',color=sns.color_palette()[2])
            #else:
            #    ax2.plot(x, y, color='k',marker='o')
                #ax2.annotate(labels[ind], (x, y), xytext=(0,-5),textcoords='offset points',va='top', ha='center')
    return ddata

def returnColoredDendrogramLinkingFunction(linkage,comparisonDictEntry):
    #Draw in colored branches on dendrogram
    pal = sns.color_palette()
    palcolorshex = pal.as_hex()
    
    ordered_cluster_compare_labels = list(map(str,sorted(comparisonDictEntry['clustersToCompare'])))
    tempdendrogram = augmented_dendrogram(1,linkage, color_threshold=None,orientation='left',no_labels=False,count_sort=True, no_plot=True) 
    ll = tempdendrogram ['leaves']
    clusters = list(map(str,comparisonDictEntry['clusterLabels']))
     
    clusterPalette = sns.color_palette(sns.color_palette(),len(comparisonDictEntry['clustersToCompare']))+[(0.5,0.5,0.5)]
    notInClusterColor = "#808080"   # Unclustered gray
    colorDict = {ordered_cluster_compare_labels[0]:palcolorshex[0],ordered_cluster_compare_labels[1]:palcolorshex[1]}
    for cluster in pd.unique(comparisonDictEntry['clusterLabels']):
        if cluster not in comparisonDictEntry['clustersToCompare']:
            colorDict[str(cluster)] = notInClusterColor
    
    D_leaf_colors = {}
    for i,l in enumerate(ll):
        D_leaf_colors[str(l)] = colorDict[str(clusters[l])]  
    link_cols = {}
    for i, i12 in enumerate(linkage[:,:2].astype(int)):
        c1, c2 = (link_cols[x] if x > len(linkage) else D_leaf_colors["%d"%x] for x in i12)
        link_cols[i+1+len(linkage)] = c1 if c1 == c2 else notInClusterColor 
    return link_cols    

def createClusterComparisonDictionary(linkage,bpd=[0]):
    
    if len(bpd) == 1:
        breakPointDistances = []
        for i in range(len(linkage)):
            breakPointDistances.append(linkage[i][2])
        breakpointDistances = sorted(breakPointDistances,reverse=True)
    else:
        breakpointDistances = bpd.copy()
    comparisonDict = {}
    increment = 0.00001
    for i in range(len(breakpointDistances)):
        clusters = list(shc.fcluster(linkage,breakpointDistances[i]+increment,'distance'))
        ival=10
        if i == ival:
            print(clusters)
        if i != 0:
            prevClusters = list(shc.fcluster(linkage,breakpointDistances[i-1]+increment,'distance'))
            diffClusters = list(pd.unique(prevClusters))
            prevClusterIndexDict = [] 
            for diffCluster in diffClusters:
                temp = []
                for j in range(len(prevClusters)):
                    if diffCluster == prevClusters[j]:
                        temp.append(j)
                prevClusterIndexDict.append(temp)
            if i == ival:
                print(prevClusters)
                print(diffClusters)
                print(prevClusterIndexDict)
            for prevClusterIndices in prevClusterIndexDict:
                brokenPrevCluster = itemgetter(*prevClusterIndices)(prevClusters)
                brokenCurrentClusters = itemgetter(*prevClusterIndices)(clusters)
                if not isinstance(brokenPrevCluster, collections.Sized):
                    brokenPrevCluster = [brokenPrevCluster]
                if not isinstance(brokenCurrentClusters, collections.Sized):
                    brokenCurrentClusters = [brokenCurrentClusters]
                if i == ival:
                    print(brokenPrevCluster)
                    print(brokenCurrentClusters)
                if len(pd.unique(brokenPrevCluster)) != len(pd.unique(brokenCurrentClusters)):
                    comparisonDict[i] = {'clusterSubsetIndices':prevClusterIndices,'clusterLabels':clusters,'clustersToCompare':list(pd.unique(brokenCurrentClusters)),'breakpointDistance':breakpointDistances[i-1]}
                    break
        else:
            comparisonDict[i] = {'clusterSubsetIndices':range(len(clusters)),'clusterLabels':clusters,'clustersToCompare':list(pd.unique(clusters))}
    return comparisonDict

def createClusterPlottingDf(df,linkage,comparisonDictEntry):
    clusters = list(map(str,comparisonDictEntry['clusterLabels']))
    clusterSubsetIndices = comparisonDictEntry['clusterSubsetIndices']
    
    numTimepoints = len(pd.unique(df.index.get_level_values('Time')))
    numFeatures = len(df.columns)

    stackedClusterSubsetIndices = []
    for clusterSubsetIndex in clusterSubsetIndices:
        startIndex = clusterSubsetIndex*numTimepoints*numFeatures
        k = 0
        for i in range(numTimepoints):
            for j in range(numFeatures):
                stackedClusterSubsetIndices.append(startIndex+k)
                k+=1
    newclusters = []
    for cluster in clusters:
        for i in range(numTimepoints):
            for j in range(numFeatures):
                newclusters.append(cluster)
    
    fullPlottingDf = df.stack().to_frame('Metric').reset_index()
    fullPlottingDf['Cluster'] = newclusters
    return fullPlottingDf

def createClusterDf(df,linkage,comparisonDictEntry):
    clusters = list(map(str,comparisonDictEntry['clusterLabels']))
    clusterSubsetIndices = comparisonDictEntry['clusterSubsetIndices']
    
    numTimepoints = len(pd.unique(df.index.get_level_values('Time')))

    stackedClusterSubsetIndices = []
    for clusterSubsetIndex in clusterSubsetIndices:
        startIndex = clusterSubsetIndex*numTimepoints
        k = 0
        for i in range(numTimepoints):
            stackedClusterSubsetIndices.append(startIndex+k)
            k+=1
    newclusters = []
    for cluster in clusters:
        for i in range(numTimepoints):
            newclusters.append(cluster)
    
    df['Cluster'] = newclusters
    return df

def returnMasks(df,comparisonDictEntry,clusters):
    
    clusters = list(map(str,list(map(int,clusters))))
    ordered_cluster_compare_labels = list(map(str,sorted(comparisonDictEntry['clustersToCompare'])))
    
    print(clusters)
    print(ordered_cluster_compare_labels)
    bluemaskrows = []
    orangemaskrows = []
    for row in range(df.shape[0]):
        if clusters[row] == ordered_cluster_compare_labels[0]:
            bluemaskrows.append(row)
        elif clusters[row] == ordered_cluster_compare_labels[1]:
            orangemaskrows.append(row)
    
    bluemask = np.ones(df.shape)
    orangemask = np.ones(df.shape)

    for bluerow in bluemaskrows:
        for col in range(df.shape[1]):
            bluemask[bluerow,col] = False
    for orangerow in orangemaskrows:
        for col in range(df.shape[1]):
            orangemask[orangerow,col] = False

    return bluemask,orangemask
