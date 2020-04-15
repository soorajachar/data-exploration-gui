#!/usr/bin/env python3 
import json,pickle,math,matplotlib,sys,os,string,subprocess
from sys import platform as sys_pf
if sys_pf == 'darwin':
    import matplotlib
    matplotlib.use("TkAgg")
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.preprocessing import MinMaxScaler,MaxAbsScaler,StandardScaler
from sklearn.decomposition import PCA
from sklearn.manifold import Isomap
from sklearn.manifold import TSNE
from fitsne import FItSNE
from umap import UMAP
from sklearn.cluster import KMeans, DBSCAN, AgglomerativeClustering
import phenograph
import clusterPlottingLibrary as cpl

#Postprocessing attributes currently implemented
scalingFunctionDict = {'minmax':MinMaxScaler,'maxabs':MaxAbsScaler,'standard':StandardScaler}

dimReductionFunctionDict = {'umap':UMAP,'tsne':TSNE,'FItSNE':FItSNE,'isomap':Isomap,'pca':PCA}
dimReductionNumericParameterDict = {'umap':['n_neighbors','min_dist'],'tsne':['perplexity'],'FItSNE':['perplexity'],'isomap':['n_neighbors'],'pca':[]}
dimReductionNumericParameterBounds = {'n_neighbors':[2,100,1,15],'min_dist':[0,1,0.05,0.1],'perplexity':[5,50,1,30]}
dimReductionQualitativeParameterDict = {'umap':['metric'],'tsne':['method'],'FItSNE':[],'isomap':[]}
dimReductionQualitativeParameterValues = {'method':['barnes_hut','exact'],'metric':['euclidean','manhattan','chebyshev','minkowski']}

clusteringFunctionDict = {'k-means':KMeans,'db-scan':DBSCAN,'hierarchical':AgglomerativeClustering,'phenograph':phenograph.cluster}
clusterParameterDict = {'k-means':['n_clusters'],'db-scan':['eps'],'hierarchical':['n_clusters','distance_threshold'],'phenograph':['k']}
clusterParameterBoundsDict = {'n_clusters':[0,20,1,2],'eps':[0,5,0.05,0.5],'distance_threshold':[0,20,0.2,0],'k':[2,100,2,30]}
        
clusterComparisonParameterList = ['confidence','fold change','response cutoff (%)']
clusterComparisonParameterBoundsDict = {'confidence':[0.1,0.01,0.01,0.05],'fold change':[1,5,0.2,2],'response cutoff (%)':[0,0.5,0.05,0.1]}

clusterComparisonParameterList2 = ['error correction method']
clusterComparisonParameterValueDict2 = {'error correction method':['holm-bonferroni','holm','none']}

clusterComparisonParameterList3 = ['comparison plot type']
clusterComparisonParameterValueDict3 = {'comparison plot type':['violin','box','bar']}

def operateOnData(df,dftitle,operationClass,operationParameters,cluster_labels=[]):
    if len(cluster_labels) > 0:
        supervised = True
    else:
        supervised = False
    
    if operationClass == 'scale':
        scalingMethod = operationParameters['scalingMethod']
        if 'downsampleFraction' in operationParameters:
            downsampleParameters = operationParameters['downsampleFraction']
            hyperparameterDict = {}
            if downsampleParameters[0] == 'fraction':
                hyperparameterDict['fraction'] = downsampleParameters[2]
            else:
                hyperparameterDict['count'] = downsampleParameters[3]
            hyperparameterDict['sampledAcross'] = downsampleParameters[1]
        else:
            hyperparameterDict = {}
        operatedDf = preprocessData(df,scalingMethod)
        operationName = scalingMethod
    elif operationClass == 'reduce':
        dimRedType = operationParameters['dimRedType']
        reductionParameters = operationParameters['reductionParameters']
        operatedDf = reduceDimensions(df,dimRedType,cluster_labels,reductionParameters)
        operationName = dimRedType
        hyperparameterDict = reductionParameters
    elif operationClass == 'cluster':
        clusteringMethod = operationParameters['clusteringMethod']
        clusteringParameters = operationParameters['clusteringParameters']
        operatedDf = clusterData(df,clusteringMethod,clusteringParameters)
        operationName = clusteringMethod
        hyperparameterDict = clusteringParameters
        
    savePostProcessedFile(operatedDf,dftitle,operationClass,operationName,hyperparameterDict,supervised=supervised)
    return operatedDf

def isFloat(string):
    try:
        float(string)
        return True
    except ValueError:
        return False

def parseFileName(fileName):
    parameterDict = {}
    if 'reducedBy' in fileName:
        root = fileName.split('reducedBy-')[1]
        parameterDict['dimRedType'] = root.split('-')[0]
        tempdict = {}
        parameterList = root.split('-')[1].split('.p')[0].split(',')
        for parameterString in parameterList:
            parameter = parameterString.split('=')[0]
            parameterVal = parameterString.split('=')[1]
            if parameterVal.isnumeric():
                parameterVal = int(parameterVal)
            elif isFloat(parameterVal):
                parameterVal = float(parameterVal)
            tempdict[parameter] = parameterVal
        parameterDict['reductionParameters'] = tempdict
    return parameterDict

def returnParameterString(pDict):
    parameterStringList = []
    for parameter in pDict:
        parameterString = parameter+'='+str(pDict[parameter])
        parameterStringList.append(parameterString)
    return ','.join(parameterStringList)

def getFileName(dataSelectionTitle,operationClass,operationName,hyperparameterDict,plotParameterDict = [],fileExtension = '.pkl'):
    if operationClass[-1] == 'e':
        operationClass = operationClass[:] + 'd'
    else:
        operationClass = operationClass[:] + 'ed'
    fileFolderName = 'outputData/analysisFiles/'+operationClass+'Data/'
    fileNamePrefix = dataSelectionTitle+'-'+operationClass+'By-'+operationName

    if len(hyperparameterDict) > 0:
        hyperParameterSegment = '-'+returnParameterString(hyperparameterDict)
    else:
        hyperParameterSegment = ''
    if len(plotParameterDict) > 0:
        print('wat')
        print(plotParameterDict)
        plotParameterSegment = '-'+returnParameterString(plotParameterDict)
        print(plotParameterSegment)
    else:
        plotParameterSegment = ''
    
    fileName = fileNamePrefix+hyperParameterSegment+plotParameterSegment+fileExtension
    
    return fileName,fileFolderName

def savePostProcessedFile(df,dataSelectionTitle,operationClass,operationName,hyperparameterDict,supervised=False):
    fileName,fileFolderName = getFileName(dataSelectionTitle,operationClass,operationName,hyperparameterDict)
    if supervised:
        fileName = 'supervised-'+fileName
    if operationClass+'Data' not in os.listdir('outputData/analysisFiles'):
        print(os.listdir('outputData/analysisFiles'))
        subprocess.run(['mkdir',fileFolderName])
    with open(fileFolderName+fileName,'wb') as f:
        pickle.dump(df,f)

def preprocessData(dataSelectionDf,scalingMethod):
    
    #Convert log-normal features to approximately gaussian distributions
    #Take log of supernatant features, remove negative mfi and take log of icell mfi features, take log of cell frequency features
    dataSelectionMatrix = dataSelectionDf.values
    for col,feature in enumerate(dataSelectionDf.columns):
        featureValues = dataSelectionMatrix[:,col]
        if 'Supernatant' in feature:
            featureValues[np.isnan(featureValues)]=1
            featureValues = np.log10(featureValues)
        elif 'Cells' in feature:
            if 'MFI' in feature or 'GFI' in feature:
                minMFI = 1
                #Move nans to 1 to allow for appropriate log scaling
                featureValues[np.isnan(featureValues)]=minMFI
                #If negative terms exist in feature, add absolute value of most negative term to allow for future log scaling
                constantToAdd = np.amin(featuerValues)
                if constantToAdd <= 0:
                    featureValues =  featureValues + minMFI + abs(constantToAdd)
                featureValues = np.log10(featureValues)
            #Linear; nans can be zero
            else:
                featureValues = np.nan_to_num(featureValues)
        #Linear; nans can be zero
        else:
            featureValues = np.nan_to_num(featureValues)
        dataSelectionMatrix[:,col] = featureValues
    
    if scalingMethod != 'none':
        scalingFunc = scalingFunctionDict[scalingMethod]
        scaledMatrix = scalingFunc().fit_transform(dataSelectionMatrix)
    else:
        scaledMatrix = dataSelectionMatrix.copy()
    preprocessedDf = pd.DataFrame(scaledMatrix,index=dataSelectionDf.index,columns=dataSelectionDf.columns)
    if 'Time' in preprocessedDf.columns.names or 'Time' in preprocessedDf.index.names:
        preprocessedDf = cpl.shrinkFeatureMultiIndex(preprocessedDf)
    else:
        preprocessedDf.columns.name = 'Feature'
    return preprocessedDf

def reduceDimensions(scaledDf,dimRedType,target,allParameters):
    dimensionReductionFunc = dimReductionFunctionDict[dimRedType]
    if dimRedType == 'pca':
        allParameters = {'n_components':2}
    if len(target) > 0:
        target = list(map(int,target))
        dimRedMatrix = dimensionReductionFunc(**allParameters).fit_transform(scaledDf.values,y=target)
    else:
        if dimRedType == 'FItSNE':
            dimRedMatrix = dimensionReductionFunc(scaledDf.values.astype(float).copy(order='C'),**allParameters)
        else:
            dimRedMatrix = dimensionReductionFunc(**allParameters).fit_transform(scaledDf.values)
    dimRedDf = pd.DataFrame(dimRedMatrix,index=scaledDf.index,columns=['Dimension 1','Dimension 2'])
    
    return dimRedDf

def clusterData(scaledDf,clusteringMethod,clusteringParameters):
    
    clusterFunc = clusteringFunctionDict[clusteringMethod]
    if clusteringMethod == 'phenograph':
        cluster_labels,graph,Q = clusterFunc(scaledDf,**clusteringParameters)
    else:
        if 'n_clusters' in clusteringParameters.keys() and 'distance_threshold' in clusteringParameters.keys():
            clusteringParameters['n_clusters'] = None
        cluster_labels = clusterFunc(**clusteringParameters).fit_predict(scaledDf)
    
    tupleList = []
    for row in range(scaledDf.shape[0]):
        names = list(scaledDf.iloc[row,:].name)
        tupleList.append(names+[str(cluster_labels[row])])
    clusteredMultiIndex = pd.MultiIndex.from_tuples(tupleList,names=list(scaledDf.index.names)+['Cluster'])
    clusteredDf = pd.DataFrame(scaledDf.values,index=clusteredMultiIndex,columns=scaledDf.columns) 
    
    return clusteredDf
