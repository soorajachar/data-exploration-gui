#!/usr/bin/env python3 
# -*- coding: utf-8 -*-
"""
created on sat jul 21 13:12:56 2018

@author: acharsr
"""
import json,pickle,math,matplotlib,sys,os,string
from sys import platform as sys_pf
if sys_pf == 'darwin':
    import matplotlib
    matplotlib.use("TkAgg")
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
from itertools import groupby
from miscFunctions import Hill,InverseHill,r_squared,cleanUpFlowjoCSV,extractValues  
from modifyDataFrames import returnModifiedDf
import initialDataProcessing as idp
import cytokineDataProcessing as cydp
import cellDataProcessing as cdp
import proliferationDataProcessing as pdp
import singleCellDataProcessing as scdp
#import automatedCBAProcessingGUI as autoCBA
import tkinter as tk


pathToExperimentSpreadsheet = '../../experiments/'
secondPath = '../../outputData'

concUnit = 1e9
unitPrefixDictionary = {1e12:'pM',1e9:'nM',1e6:'uM',1e3:'mM',1e0:'M'}
concUnitPrefix = unitPrefixDictionary[concUnit]

FACS_Detector_Names = ['BV421-A','BV510-A','BV605-A','BV605-A','BV650-A','BV711-A','BV786-A','BUV396-A','DAPI-A','BUV737-A','APC-A','Alexa Fluor 700-A','APC-Cy7-A','FITC-A','PerCP-Cy5-5-A','PE ( 561 )-A','PE-CF594-A','PE-Cy5-A','PE-Cy7-A']
CyTOF_Detector_Names = []
Full_Detector_Dict = {'FACS':FACS_Detector_Names,'CyTOF':CyTOF_Detector_Names}

class DataProcessingStartPage(tk.Frame):
    def __init__(self, master,folderName,expNum,ex_data,bPage):
        tk.Frame.__init__(self, master)
            
        #os.chdir(master.homedirectory+'/'+folderName)
        backPage = bPage

        experimentNameWindow = tk.Frame(self)
        experimentNameWindow.pack(side=tk.TOP,padx=10,pady=10)
        experimentNameLabel = tk.Label(experimentNameWindow,text=folderName+':').pack()

        mainWindow = tk.Frame(self)
        mainWindow.pack(side=tk.TOP,padx=10,pady=10)
        
        l2 = tk.Label(mainWindow, text="""Datatype: """,padx = 20).grid(row=0,column=0,sticky=tk.W)
        l2a = tk.Label(mainWindow, text="Cytokine:",padx = 20).grid(row=1,column=0,sticky=tk.W)
        l2b = tk.Label(mainWindow,text="Cell:",padx = 20).grid(row=2,column=0,sticky=tk.W)
        l2c = tk.Label(mainWindow,text="Proliferation:",padx = 20).grid(row=3,column=0,sticky=tk.W)
        l2d = tk.Label(mainWindow,text="Single Cell:",padx = 20).grid(row=4,column=0,sticky=tk.W)
        
        def createDataFrame(dataType):
            dataProcessingMaster(folderName,expNum,dataType,ex_data,v3.get())
        
        l3 = tk.Label(mainWindow, text="""Action: """).grid(row=0,column=1,sticky=tk.W)
        
        cytCalibrationParametersButton = tk.Button(mainWindow,text='Enter CBA bead calibration parameters',command=lambda: master.switch_frame(cydp.CalibrationParameterPage,folderName,expNum,ex_data,DataProcessingStartPage,bPage))
        cytCalibrationParametersButton.grid(row=1,column=1,sticky=tk.W)
        #cytCalibrationParametersButton = tk.Button(mainWindow,text='Create CBA gates',command=lambda: master.switch_frame(autoCBA.AutomateCBAStartPage,folderName,expNum,ex_data,DataProcessingStartPage,bPage))
        #cytCalibrationParametersButton.grid(row=1,column=2,sticky=tk.W)
        cytDfButton = tk.Button(mainWindow,text='Create dataframe',command=lambda: createDataFrame('cyt'))
        cytDfButton.grid(row=1,column=2,sticky=tk.W)

        cellDfButton = tk.Button(mainWindow,text='Create dataframe',command=lambda: createDataFrame('cell'))
        cellDfButton.grid(row=2,column=2,sticky=tk.W)
        #cellAbPanelButton = tk.Button(mainWindow,text='Edit antibody panel',command=lambda: master.switch_frame(cdp.MarkerNumberPage,folderName,expNum,ex_data,DataProcessingStartPage,bPage))
        #cellAbPanelButton.grid(row=2,column=2,sticky=tk.W)
        
        prolifGenerationGatesButton = tk.Button(mainWindow,text='Edit generation gates',command=lambda: master.switch_frame(pdp.GatingPage,folderName,expNum,ex_data,DataProcessingStartPage,bPage))
        prolifGenerationGatesButton.grid(row=3,column=1,sticky=tk.W)
        prolifDfButton = tk.Button(mainWindow,text='Create dataframe',command=lambda: createDataFrame('prolif'))
        prolifDfButton.grid(row=3,column=2,sticky=tk.W)

        #completeSingleCellDfButton = tk.Button(mainWindow,text='Create complete dataframes')
        #completeSingleCellDfButton.grid(row=4,column=2,sticky=tk.W)
        singleCellDfButton = tk.Button(mainWindow,text='Create dataframe',command=lambda: createDataFrame('singlecell'))
        singleCellDfButton.grid(row=4,column=2,sticky=tk.W)
        cbWindow = tk.Frame(mainWindow)
        cbWindow.grid(row=4,column=1,sticky=tk.W)
        l3 = tk.Label(cbWindow,text='Use empty wells?').grid(row=0,column=0,sticky=tk.W)
        v3 = tk.BooleanVar(value=False)
        cb = tk.Checkbutton(cbWindow, variable=v3)
        cb.grid(row=0,column=1,sticky=tk.W)

        for i,button in enumerate([cytDfButton,cellDfButton,prolifDfButton,prolifGenerationGatesButton,singleCellDfButton]):
            if i == 0:
                requiredFiles = ['CBAcalibrationParameters-'+folderName+'.json']
            elif i == 1:
                requiredFiles = []
            elif i == 2:
                requiredFiles = ['singleCellDataFrame-proliferation-'+folderName+'.pkl']
            elif i == 3:
                requiredFiles = ['logicleProliferationDf.pkl','rawProliferationDf.pkl']
            elif i ==4:
                requiredFiles = []
                #requiredFiles = ['A1_cell.csv']
            else:
                requiredFiles = ['initialSingleCellDf-channel-'+folderName+'.pkl']
            for requiredFile in requiredFiles:
                if requiredFile not in os.listdir('misc')+os.listdir('inputData/bulkCSVFiles')+os.listdir('outputData/pickleFiles'):
                    button.config(state=tk.DISABLED)
        
        buttonWindow = tk.Frame(self)
        buttonWindow.pack(side=tk.TOP,pady=10)

        #tk.Button(buttonWindow, text="OK",command=lambda: collectInputs()).grid(row=5,column=0)
        tk.Button(buttonWindow, text="Back",command=lambda: master.switch_frame(backPage,folderName)).grid(row=5,column=1)
        tk.Button(buttonWindow, text="Quit",command=quit).grid(row=5,column=2)

def dataProcessingMaster(folderName,expNum,dataType,ex_data,useBlankWells):
    print('Creating Dataframes for: '+str(folderName))
    if dataType == 'singlecell' or dataType == 'prolif':
        parameterExtension = 'cell'
    else:
        parameterExtension = dataType
    experimentParameters = json.load(open('misc/experimentParameters-'+folderName+'-'+parameterExtension+'.json','r'))
    if experimentParameters['format'] == 'plate':
        experimentFormat = 'plate'
        experimentLevelLayoutDict = pickle.load(open('misc/layoutDict-'+folderName+'-'+parameterExtension+'.pkl','rb'))
    else:
        experimentFormat = 'tube'
        experimentLevelLayoutDict = pickle.load(open('misc/tubeLayout-'+folderName+'-'+parameterExtension+'.pkl','rb'))
    #experimentLevelLayoutDict = idp.tilePlateLayouts(experimentParameters,levelLayouts)
    if(dataType == 'cyt'):
        calibrationParameters = json.load(open('misc/CBAcalibrationParameters-'+folderName+'.json','r'))
        numberOfCalibrationSamples = calibrationParameters['Number']
        initialStandardVolume = calibrationParameters['Volume']
        cydp.calibrateExperiment(folderName,secondPath,concUnit,concUnitPrefix,numberOfCalibrationSamples,initialStandardVolume)
        basecytdf = idp.createBaseDataFrame(experimentParameters,folderName,expNum,dataType,experimentLevelLayoutDict)
        cytdf = cydp.createCytokineDataFrame(folderName,basecytdf,concUnitPrefix)
        idp.saveFinalDataFrames(folderName,secondPath,expNum,dataType,cytdf,ex_data) 
    elif(dataType == 'cell'):
        celldf = idp.createBaseDataFrame(experimentParameters,folderName,expNum,dataType,experimentLevelLayoutDict)
        idp.saveFinalDataFrames(folderName,secondPath,expNum,dataType,celldf,ex_data) 
    elif(dataType == 'prolif'):
        prolifdf = pdp.generateBulkProliferationStatistics(folderName,expNum)
        idp.saveFinalDataFrames(folderName,secondPath,expNum,dataType,prolifdf,ex_data) 
    elif(dataType == 'singlecell'):
        dataType = 'singlecell'
        if experimentFormat == 'plate':
            scdp.createPlateSingleCellDataFrame(folderName,experimentParameters,experimentLevelLayoutDict,useBlankWells)
        else:
            scdf = scdp.createTubeSingleCellDataFrame(folderName,experimentParameters,experimentLevelLayoutDict)
    print(dataType+' dataframe created!')
