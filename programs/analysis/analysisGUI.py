#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thursday July 5 13:52:27 2018

@author: acharsr
"""
import pickle,sys,os,json,math,subprocess,string
from sys import platform as sys_pf
if sys_pf == 'darwin':
    import matplotlib
    matplotlib.use("TkAgg")
from pathlib import Path
import pickle
from itertools import product,combinations
import numpy as np
import tkinter as tk
import pandas as pd
import facetPlotLibrary as fpl
sys.path.insert(0, '../dataProcessing/')
from miscFunctions import sortSINumerically,reindexDataFrame
from dimensionReductionGUI import DimensionReductionHomePage
from clusteringGUI import ClusteringHomePage
from clusterComparisonGUI import ClusterComparisonHomePage
from dataFrameValueSelectionGUI import DataSelectionHomePage
from clusterFrequencyGUI import ClusterFrequencyHomePage

class AnalysisStartPage(tk.Frame):
    num_args = 1
    def __init__(self, master,folderName,mainhomepage):
        tk.Frame.__init__(self, master)
        
        experimentNameWindow = tk.Frame(self)
        experimentNameWindow.pack(side=tk.TOP,padx=10,pady=10)
        experimentNameLabel = tk.Label(experimentNameWindow,text=folderName+':').pack()

        mainWindow = tk.Frame(self)
        mainWindow.pack(side=tk.TOP,padx=10,pady=10)
        
        l1 = tk.Label(mainWindow, text="""Postprocessing Type:""").grid(row=0,column=0)
        postProcessTypeVar = tk.StringVar(value='subset')
        rb1b = tk.Radiobutton(mainWindow,text='Select Data Subet',variable=postProcessTypeVar,value='subset')
        rb1b.grid(row=1,column=0,sticky=tk.W)
        rb1a = tk.Radiobutton(mainWindow,text='Dimensionality Reduction',variable=postProcessTypeVar,value='dr')
        rb1a.grid(row=2,column=0,sticky=tk.W)
        rb2a = tk.Radiobutton(mainWindow,text='Clustering',variable=postProcessTypeVar,value='clustering')
        rb2a.grid(row=3,column=0,sticky=tk.W)
        rb3a = tk.Radiobutton(mainWindow,text='Cluster Comparison',variable=postProcessTypeVar,value='comparison')
        rb3a.grid(row=4,column=0,sticky=tk.W)
        rb4a = tk.Radiobutton(mainWindow,text='Cluster Frequency Analysis',variable=postProcessTypeVar,value='frequency')
        rb4a.grid(row=5,column=0,sticky=tk.W)

        def collectInputs():
            postProcessType = postProcessTypeVar.get()
            if postProcessType == 'subset':
                master.switch_frame(DataSelectionHomePage,AnalysisStartPage,folderName,AnalysisStartPage,'cluster',mainhomepage)
            elif postProcessType == 'dr':
                master.switch_frame(DimensionReductionHomePage,folderName,AnalysisStartPage,mainhomepage)
            elif postProcessType == 'clustering':
                master.switch_frame(ClusteringHomePage,folderName,AnalysisStartPage,mainhomepage)
            elif postProcessType == 'comparison':
                master.switch_frame(ClusterComparisonHomePage,folderName,AnalysisStartPage,mainhomepage)
            elif postProcessType == 'frequency':
                master.switch_frame(ClusterFrequencyHomePage,folderName,AnalysisStartPage,mainhomepage)

        buttonWindow = tk.Frame(self)
        buttonWindow.pack(side=tk.TOP,pady=10)
        tk.Button(buttonWindow, text="OK",command=lambda: collectInputs()).pack(in_=buttonWindow,side=tk.LEFT)
        tk.Button(buttonWindow, text="Back",command=lambda: master.switch_frame(mainhomepage,folderName)).pack(in_=buttonWindow,side=tk.LEFT)
        tk.Button(buttonWindow, text="Quit",command=lambda: quit()).pack(in_=buttonWindow,side=tk.LEFT)
