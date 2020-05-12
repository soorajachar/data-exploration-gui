#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thursday July 5 13:52:27 2018

@author: acharsr
"""
from pathlib import Path
import pickle
from itertools import product,combinations
import numpy as np
import pickle,sys,os,json,math,subprocess,string
from tkinter import *
import tkinter as tk
from createTubeLayout import TubeLayoutPage 
from createPlateLayout import BlankSelectionPage
import pandas as pd

experimentParameters = {}
parametersUpdatedByGridGUI = {}

class ExperimentSetupStartPage(tk.Frame):
    def __init__(self, master,fName,bPage):
        global folderName,backPage
        folderName = fName
        backPage = bPage
        tk.Frame.__init__(self, master)
        experimentNameWindow = tk.Frame(self)
        experimentNameWindow.pack(side=tk.TOP,padx=10,pady=10)
        experimentNameLabel = tk.Label(experimentNameWindow,text=folderName+':').pack()
        
        mainWindow = tk.Frame(self)
        mainWindow.pack(side=tk.TOP,padx=10,pady=10)
        
        v2 = tk.StringVar(value='inpt')
        rb2a = tk.Radiobutton(mainWindow, text="Enter experiment labels ",padx = 20, variable=v2, value='inpt')
        rb2b = tk.Radiobutton(mainWindow,text="Create experiment layout ",padx = 20, variable=v2, value='pl')
        rb2a.grid(row=0,column=0,sticky=tk.W)
        rb2b.grid(row=1,column=0,sticky=tk.W)
        
        v3 = tk.StringVar(value='both')
        l3 = tk.Label(mainWindow, text="for:     ")
        l3.grid(row=0,column=1,rowspan=2,sticky=tk.EW)
        rb3a = tk.Radiobutton(mainWindow, text="Cells",padx = 20, variable=v3, value='cell')
        rb3b = tk.Radiobutton(mainWindow,text="Cytokines",padx = 20, variable=v3, value='cyt')
        rb3c = tk.Radiobutton(mainWindow,text="Cells and Cytokines",padx = 20, variable=v3, value='both')
        rb3a.grid(row=0,column=2,sticky=tk.W)
        rb3b.grid(row=1,column=2,sticky=tk.W)
        rb3c.grid(row=2,column=2,sticky=tk.W)
        
        def experimentLayout():
            if v3.get() == 'both':
                if 'experimentParameters-'+folderName+'-cyt.json' in os.listdir('misc'):
                    experimentParameters = json.load(open('misc/experimentParameters-'+folderName+'-cyt.json','r'))
                elif 'experimentParameters-'+folderName+'-cell.json' in os.listdir('misc'):
                    experimentParameters = json.load(open('misc/experimentParameters-'+folderName+'-cell.json','r'))
                else:
                    experimentParameters = json.load(open('misc/experimentParameters-'+folderName+'.json','r'))
            else:
                experimentParameters = json.load(open('misc/experimentParameters-'+folderName+'-'+v3.get()+'.json','r'))
            print(experimentParameters)
            if 'format' not in experimentParameters.keys():
                experimentParameters['format'] = 'plate'
            if experimentParameters['format'] == 'plate':
                if 'paired' in experimentParameters.keys():
                    if experimentParameters['paired']:
                        numRowPlates = 2
                    else:
                        numRowPlates = 1
                    numColumnPlates = int(experimentParameters['numPlates']/numRowPlates)
                else:
                    numRowPlates = experimentParameters['numRowPlates'] 
                    numColumnPlates = experimentParameters['numColumnPlates']
                levels = experimentParameters['allLevelNames']
                conditionLevelValues = experimentParameters['conditionLevelValues']
                plateDimensions = experimentParameters['overallPlateDimensions']
                levelValues = []
                for level in levels:
                    levelValues.append(experimentParameters['allLevelValues'][level])
                maxNumLevelValues = len(max(levelValues,key=len))

                master.switch_frame(BlankSelectionPage,folderName,levels,levelValues,maxNumLevelValues,numRowPlates,numColumnPlates,plateDimensions,v3.get(),ExperimentSetupStartPage,bPage)
            #Tube mode
            else:
                master.switch_frame(TubeLayoutPage,folderName,experimentParameters['conditionLevelValues'],experimentParameters['columnLevelValues'],experimentParameters['numSamples'],experimentParameters['numericLevels'],experimentParameters['allLevelNames'],v3.get(),ExperimentSetupStartPage,bPage)
        
        def collectInput():
            global dataType
            dataType = v3.get()
            if v2.get() == 'inpt':
                master.switch_frame(ExperimentFormatPage,folderName)
            elif v2.get() == 'pl':
                experimentLayout()

        buttonWindow = tk.Frame(self)
        buttonWindow.pack(side=TOP,pady=10)

        tk.Button(buttonWindow, text="OK",command=lambda: collectInput()).grid(row=5,column=0)
        tk.Button(buttonWindow, text="Back",command=lambda: master.switch_frame(backPage,folderName)).grid(row=5,column=1)
        tk.Button(buttonWindow, text="Quit",command=quit).grid(row=5,column=2)

class ExperimentFormatPage(tk.Frame):
    def __init__(self, master,folderName):
        tk.Frame.__init__(self, master)
        
        mainWindow = tk.Frame(self)
        mainWindow.pack(side=tk.TOP,padx=10,pady=10)
        
        v1 = tk.StringVar()
        l1 = tk.Label(mainWindow,text='Collection Format: ')
        rb2a = tk.Radiobutton(mainWindow, text="Tube",padx = 20, variable=v1, value='tube')
        rb2b = tk.Radiobutton(mainWindow,text="Plate",padx = 20, variable=v1, value='plate')
        v1.set('tube')
        
        l1.grid(row=0,column=0)
        rb2a.grid(row=1,column=0)
        rb2b.grid(row=2,column=0)
        
        def collectInputs():
            if v1.get() == 'tube':
                experimentParameters['format'] = 'tube'
                master.switch_frame(TubeExperimentParameterPage,folderName)
            else:
                experimentParameters['format'] = 'plate'
                master.switch_frame(PlateExperimentParameterPage,folderName)

        buttonWindow = tk.Frame(self)
        buttonWindow.pack(side=TOP,pady=10)

        tk.Button(buttonWindow, text="OK",command=lambda: collectInputs()).grid(row=5,column=0)
        tk.Button(buttonWindow, text="Back",command=lambda: master.switch_frame(ExperimentSetupStartPage,folderName,backPage)).grid(row=5,column=1)
        tk.Button(buttonWindow, text="Quit",command=quit).grid(row=5,column=2)

class TubeExperimentParameterPage(tk.Frame):
    def __init__(self, master,folderName):
        tk.Frame.__init__(self, master)
        
        mainWindow = tk.Frame(self)
        mainWindow.pack(side=tk.TOP,padx=10,pady=10)
        

        l4 = tk.Label(mainWindow, text="How many samples are there in this experiment?")
        e0 = tk.Entry(mainWindow)
        l4.grid(row=0,column=0)
        e0.grid(row=0,column=1)
        
        l5 = tk.Label(mainWindow, text="How many unique levels are there in this experiment (excluding Time)?")
        e1 = tk.Entry(mainWindow)
        l5.grid(row=2,column=0)
        e1.grid(row=2,column=1)
        
        if 'sampleNameFile.xlsx' in os.listdir('misc'):
            df = pd.read_excel('misc/sampleNameFile.xlsx')
            columns = []
            for column in df:
                if column not in ['','FileName','Time']:
                    columns.append(column)
            e0.insert(tk.END,str(df.shape[0]))
            e1.insert(tk.END,str(len(columns)))

        def collectInputs():
            experimentParameters['numSamples'] = int(e0.get())
            experimentParameters['numAllLevels'] = int(e1.get())+1
            master.switch_frame(allLevelNamePage,folderName)

        buttonWindow = tk.Frame(self)
        buttonWindow.pack(side=TOP,pady=10)

        tk.Button(buttonWindow, text="OK",command=lambda: collectInputs()).grid(row=5,column=0)
        tk.Button(buttonWindow, text="Back",command=lambda: master.switch_frame(ExperimentFormatPage,folderName)).grid(row=5,column=1)
        tk.Button(buttonWindow, text="Quit",command=quit).grid(row=5,column=2)

class PlateExperimentParameterPage(tk.Frame):
    def __init__(self, master,folderName):
        tk.Frame.__init__(self, master)
        
        Multiplexing_Options = {'96->384 Well Plate','Barcoding'}

        mainWindow = tk.Frame(self)
        mainWindow.pack(side=tk.TOP,padx=10,pady=10)
        
        v = tk.IntVar()
        v2 = tk.IntVar()
        v3 = tk.IntVar()

        l3 = tk.Label(mainWindow, text="""What format were the samples collected in?:""")
        rb3a = tk.Radiobutton(mainWindow, text="96 well plate",padx = 20, variable=v3, value=96)
        rb3b = tk.Radiobutton(mainWindow,text="384 well plate",padx = 20, variable=v3, value=384)
        l3.grid(row=0,column=0)
        rb3a.grid(row=0,column=1)
        rb3b.grid(row=0,column=2)
        v3.set(96)

        l4 = tk.Label(mainWindow, text="How many plates was a single timepoint's conditions spread over (rowPlates)?")
        e0 = tk.Entry(mainWindow)
        l4.grid(row=1,column=0)
        e0.grid(row=1,column=1)
        
        l5 = tk.Label(mainWindow, text="How many plates was a single condition's timepoints spread over (columnPlates)?")
        e1 = tk.Entry(mainWindow)
        l5.grid(row=2,column=0)
        e1.grid(row=2,column=1)

        l6 = tk.Label(mainWindow, text="Enter the number of condition levels (excluding Time): ")
        e2 = tk.Entry(mainWindow)
        l6.grid(row=3,column=0)
        e2.grid(row=3,column=1)
        
        l7 = tk.Label(mainWindow, text="Multiplexing Options:")
        multiplex = tk.StringVar(value='None')
        rb1 = tk.Radiobutton(mainWindow,variable=multiplex,value='None',text='None')
        rb2 = tk.Radiobutton(mainWindow,variable=multiplex,value='96->384 well',text='96->384 well')
        rb3 = tk.Radiobutton(mainWindow,variable=multiplex,value='Barcoding',text='Barcoding')
        rb4 = tk.Radiobutton(mainWindow,variable=multiplex,value='96->384 well + Barcoding',text='96->384 well + Barcoding')

        l7.grid(row=4,column=0)
        rb1.grid(row=4,column=1)
        rb2.grid(row=4,column=2)
        rb3.grid(row=4,column=3)
        rb4.grid(row=4,column=4)
         
        def collectInputs():
            print(e0.get())
            print(e1.get())
            experimentParameters['numPlates'] = int(e0.get())*int(e1.get())
            experimentParameters['numRowPlates'] = int(e0.get())
            experimentParameters['numColumnPlates'] = int(e1.get())
            experimentParameters['numAllLevels'] = int(e2.get())+1
            if v3.get() == 384:
                experimentParameters['overallPlateDimensions'] = [16,24]
                parametersUpdatedByGridGUI['currentPlateDimensions'] = [16,24]
            else:
                experimentParameters['overallPlateDimensions'] = [8,12]
                parametersUpdatedByGridGUI['currentPlateDimensions'] = [8,12]
            
            multiplexingOption = multiplex.get()
            experimentParameters['multiplexingOption'] = multiplexingOption
            if multiplexingOption == 'None':
                master.switch_frame(allLevelNamePage,folderName)
            else:
                master.switch_frame(multiplexingPage,multiplexingOption)

        buttonWindow = tk.Frame(self)
        buttonWindow.pack(side=TOP,pady=10)

        tk.Button(buttonWindow, text="OK",command=lambda: collectInputs()).grid(row=5,column=0)
        tk.Button(buttonWindow, text="Back",command=lambda: master.switch_frame(ExperimentFormatPage,folderName)).grid(row=5,column=1)
        tk.Button(buttonWindow, text="Quit",command=quit).grid(row=5,column=2)
    
class multiplexingPage(tk.Frame):
    def __init__(self, master,mpo):
        tk.Frame.__init__(self, master)
        multiplexingOption = mpo
        Well384ConversionWindow = tk.Frame(self)
        Well384ConversionWindow.pack()
        if multiplexingOption == '96->384 well':
            maxNum384Plates = math.ceil(experimentParameters['numPlates'] / 4)
            Unpacking384PlateNameList = []
            Unpacking384List = []
            l1 = tk.Label(Well384ConversionWindow,text='Combined Plate Names: ')
            l1.grid(row=0,column=0)
            for plateNum in range(maxNum384Plates):
                e1 = tk.Entry(Well384ConversionWindow,width=5)
                e1.grid(row=0,column=plateNum*2+1)
                Unpacking384PlateNameList.append(e1)
                UnpackingWellPosList = []
                for wellUnpackingRow in range(2):
                    for wellUnpackingCol in range(2):
                        e2 = tk.Entry(Well384ConversionWindow,width=5)
                        if wellUnpackingCol == 1:
                            padvar = 5
                        else:
                            padvar = 0
                        e2.grid(row=wellUnpackingRow+1,column=wellUnpackingCol+1+plateNum*2,padx=padvar)
                        UnpackingWellPosList.append(e2)
                Unpacking384List.append(UnpackingWellPosList)

        elif multiplexingOption == 'Barcoding':
            print('wat2')
            pass
        else:
            print('wat3')
            pass
        
        def collectInputs():
            unpackingDict = {}
            for unpackingPlateNameEntry,unpackingWellPosList in zip(Unpacking384PlateNameList,Unpacking384List):
                tempList = []
                for unpackingEntry in unpackingWellPosList:
                    tempList.append(unpackingEntry.get())
                unpackingDict[unpackingPlateNameEntry.get()] = tempList
            experimentParameters['unpackingDict'] = unpackingDict
            master.switch_frame(allLevelNamePage,folderName)

        buttonWindow = tk.Frame(self)
        buttonWindow.pack(side=TOP,pady=10)

        tk.Button(buttonWindow, text="OK",command=lambda: collectInputs()).grid(row=5,column=0)
        tk.Button(buttonWindow, text="Back",command=lambda: master.switch_frame(PlateExperimentParameterPage,folderName)).grid(row=5,column=1)
        tk.Button(buttonWindow, text="Quit",command=quit).grid(row=5,column=2)
        
class allLevelNamePage(tk.Frame):
    def __init__(self, master,folderName):
        tk.Frame.__init__(self, master)
        numAllLevels = experimentParameters['numAllLevels']
            
        entryList1 = []
        entryList2 = []
        
        numericCheckBoxes = []
        numericCheckBoxVars = []

        mainWindow = Frame(self)
        mainWindow.pack(side=TOP,padx=10,pady=10)
        
        lt1 = tk.Label(mainWindow, text="Level Name").grid(row=0,column=1)
        lt2 = tk.Label(mainWindow, text="Number of Level Values").grid(row=0,column=2)
        lt3 = tk.Label(mainWindow, text="Numeric?").grid(row=0,column=3)
        
        if 'sampleNameFile.xlsx' in os.listdir('misc'):
            df = pd.read_excel('misc/sampleNameFile.xlsx')
            levels = []
            levelValueNums = []
            for level in df.columns:
                if level not in ['','FileName']:
                    levels.append(level)
                    levelValueNums.append(len(pd.unique(df[level])))
            j=0
        
        for conditionLevelNumber in range(1,numAllLevels+1):
            l1 = tk.Label(mainWindow, text="Condition "+str(conditionLevelNumber))
            e1 = tk.Entry(mainWindow)
            v = tk.IntVar()
            if conditionLevelNumber == 1:
                #e1 = tk.Entry(mainWindow,state='readonly')
                e1.insert(tk.END, 'Time')
                e2 = tk.Entry(mainWindow)
                if 'sampleNameFile.xlsx' in os.listdir('misc'):
                    if 'Time' in df.columns:
                        e2.insert(tk.END,str(len(pd.unique(df['Time']))))
                    else: 
                        e2.insert(tk.END, '1')
                else:
                    e2.insert(tk.END, '1')
                v.set(1)
            else:
                e1 = tk.Entry(mainWindow)
                e2 = tk.Entry(mainWindow)
                v.set(0)
                if 'sampleNameFile.xlsx' in os.listdir('misc'):
                    e1.insert(tk.END, levels[j])
                    e2.insert(tk.END, str(levelValueNums[j]))
                    j+=1
                                    
            
            cb1 = tk.Checkbutton(mainWindow, text="",variable=v,onvalue=1,offvalue=0)
             
            l1.grid(row=conditionLevelNumber+1,column=0)
            e1.grid(row=conditionLevelNumber+1,column=1)
            e2.grid(row=conditionLevelNumber+1,column=2)
            entryList1.append(e1)
            entryList2.append(e2)
            cb1.grid(row=conditionLevelNumber+1,column=3)
            numericCheckBoxes.append(cb1)
            numericCheckBoxVars.append(v)
        
        def collectInputs():
            conditionNames = []
            numConditionLevelValues = []
            tiledLevels = []
            experimentParameters['columnVariableName'] = 'Time' 
            timebool=False
            numericlevels = []
            for allLevelNumber in range(numAllLevels):
                #Remove column variable from condition name list
                if entryList1[allLevelNumber].get() == 'Time':
                    timebool=True
                    experimentParameters['numColumnLevelValues'] = int(entryList2[allLevelNumber].get())
                    experimentParameters['numColumnLevelValues'] = int(entryList2[allLevelNumber].get())
                else:
                    conditionNames.append(str(entryList1[allLevelNumber].get()))
                    numConditionLevelValues.append(int(entryList2[allLevelNumber].get()))
                numericlevels.append(numericCheckBoxVars[allLevelNumber].get() == 1)

            experimentParameters['numConditionLevels'] = numAllLevels - 1
            experimentParameters['conditionLevelNames'] = conditionNames
            experimentParameters['allLevelNames'] = [experimentParameters['columnVariableName']]+conditionNames
            experimentParameters['numConditionLevelValues'] = numConditionLevelValues
            experimentParameters['numericLevels'] = numericlevels
            parametersUpdatedByGridGUI['numLevelsUnparsed'] = numAllLevels
            experimentParameters[''] = tiledLevels
            with open('misc/gui-parametersUpdatedByGridGUI.pkl','wb') as f:
                pickle.dump(parametersUpdatedByGridGUI,f)
            master.switch_frame(columnLevelValuesPage,folderName)
        
        def backCommand():
            if experimentParameters['format'] == 'tube':
                master.switch_frame(TubeExperimentParameterPage,folderName)
            else:
                master.switch_frame(PlateExperimentParameterPage,folderName)
        
        buttonWindow = Frame(self)
        buttonWindow.pack(side=TOP,pady=10)
        
        tk.Button(buttonWindow, text="OK",command=lambda: collectInputs()).grid(row=numAllLevels+1,column=0)
        tk.Button(buttonWindow, text="Back",command=lambda: backCommand()).grid(row=numAllLevels+1,column=1)
        tk.Button(buttonWindow, text="Quit",command=quit).grid(row=numAllLevels+1,column=2)

class columnLevelValuesPage(tk.Frame):
    def __init__(self, master,folderName):
        numColumnLevelValues = experimentParameters['numColumnLevelValues']
        
        tk.Frame.__init__(self, master)
        
        mainWindow = Frame(self)
        mainWindow.pack(side=TOP,padx=10,pady=10)
        
        lt = tk.Label(mainWindow,text=experimentParameters['columnVariableName']+':').grid(row=0,column=0)
        col_wrap = 12
        for col in range(1,numColumnLevelValues+1):
            lt1 = tk.Label(mainWindow, text='Level Value '+str(col),width=10).grid(row=int((col-1)/col_wrap)*2,column=((col-1)%col_wrap)+1)
        entryList = []
        if 'sampleNameFile.xlsx' in os.listdir('misc'):
            df = pd.read_excel('misc/sampleNameFile.xlsx')
            if 'Time' in df.columns:
                times = list(map(float,list(pd.unique(df['Time']))))
                sortedTimes = list(map(str,sorted(times)))
                j=0

        for columnLevelValueNumber in range(numColumnLevelValues):
            e1 = tk.Entry(mainWindow,width=10)
            e1.grid(row=2*int(columnLevelValueNumber/col_wrap)+1,column=(columnLevelValueNumber%col_wrap)+1)
            if 'sampleNameFile.xlsx' in os.listdir('misc'):
                if 'Time' in df.columns:
                    e1.insert(tk.END,sortedTimes[j]) 
                    j+=1
            entryList.append(e1)

        def collectInputs():
            columnLevelValues = []
            for entry in entryList:
                columnLevelValues.append(float(entry.get()))
            experimentParameters['columnLevelValues'] = columnLevelValues
            master.switch_frame(conditionLevelValuesPage,folderName)
        
        buttonWindow = Frame(self)
        buttonWindow.pack(side=TOP,pady=10)
        
        tk.Button(buttonWindow, text="OK",command=lambda: collectInputs()).grid(row=2*int(numColumnLevelValues/col_wrap)+2,column=5)
        tk.Button(buttonWindow, text="Back",command=lambda: master.switch_frame(allLevelNamePage,folderName)).grid(row=2*int(numColumnLevelValues/col_wrap)+2,column=6)
        tk.Button(buttonWindow, text="Quit",command=quit).grid(row=2*int(numColumnLevelValues/col_wrap)+2,column=7)

class LinkedEntryButton(tk.Button):
    def __init__(self,parent,linkedentries,**kwargs):
        tk.Button.__init__(self,parent,**kwargs)
        self.parent = parent
        self.linkedentries = linkedentries

    def fillNumerically(self):
        for i,entry in enumerate(self.linkedentries):
            entry.insert(tk.END,str(i+1))

class conditionLevelValuesPage(tk.Frame):
    def __init__(self, master,folderName):
        numConditionLevels = experimentParameters['numConditionLevels']
        maxLevelValues = max(experimentParameters['numConditionLevelValues'])
        
        tk.Frame.__init__(self, master)
        
        mainWindow = Frame(self)
        mainWindow.pack(side=TOP,padx=10,pady=10)
        
        entryWrap = 16
        fullEntryList = []
        blist = []

        labelList = experimentParameters['conditionLevelNames']
        rowNum = 0
        if 'sampleNameFile.xlsx' in os.listdir('misc'):
            df = pd.read_excel('misc/sampleNameFile.xlsx')
        for conditionLevelNumber in range(numConditionLevels):
            if 'sampleNameFile.xlsx' in os.listdir('misc'):
                j = 0
            currentLabel = labelList[conditionLevelNumber]
            initialRowNum = rowNum
            l1 = tk.Label(mainWindow, text='Level values for \"'+currentLabel+'\":').grid(row=rowNum,column=1,sticky=tk.W) 
            levelEntryList = []
            rowNum-=1
            for col in range(1,maxLevelValues+1):
                if 'sampleNameFile.xlsx' in os.listdir('misc'):
                    for column in df.columns:
                        if column == currentLabel:
                            filledInLevelValues = list(pd.unique(df[column]))
                if (col-1)%entryWrap == 0:
                    rowNum+=1
                if col < experimentParameters['numConditionLevelValues'][conditionLevelNumber]+1:
                    e1 = tk.Entry(mainWindow,width=8)
                    e1.grid(row=rowNum,column=(col-1)%entryWrap+2,sticky=tk.W)
                    if 'sampleNameFile.xlsx' in os.listdir('misc'):
                        e1.insert(tk.END,filledInLevelValues[j])
                        j+=1
                    levelEntryList.append(e1)
            b = LinkedEntryButton(mainWindow,levelEntryList,text='Fill Values Numerically')
            b.configure(command=b.fillNumerically)
            b.grid(row=initialRowNum,column=0,sticky=tk.W)
            blist.append(b)
            fullEntryList.append(levelEntryList)
            rowNum+=1

        def collectInputs():
            conditionLevels = {}
            for lvlentrylist,i in zip(fullEntryList,range(numConditionLevels)):
                tempLevels = []
                for entry in lvlentrylist:
                    tempLevels.append(str(entry.get()))
                conditionLevels[experimentParameters['conditionLevelNames'][i]] = tempLevels
            experimentParameters['conditionLevelValues'] = conditionLevels.copy()
            experimentParameters['allLevelValues'] = conditionLevels.copy()
            experimentParameters['allLevelValues'][experimentParameters['columnVariableName']] = experimentParameters['columnLevelValues']
            if dataType == 'both':
                with open('misc/experimentParameters-'+folderName+'-cell.json', 'w') as fp:
                    json.dump(experimentParameters, fp)
                with open('misc/experimentParameters-'+folderName+'-cyt.json', 'w') as fp:
                    json.dump(experimentParameters, fp)
            else:
                with open('misc/experimentParameters-'+folderName+'-'+dataType+'.json', 'w') as fp:
                    json.dump(experimentParameters, fp)
            master.switch_frame(ExperimentSetupStartPage,folderName,backPage)

        buttonWindow = Frame(self)
        buttonWindow.pack(side=TOP,pady=10)
        
        tk.Button(buttonWindow, text="OK",command=lambda: collectInputs()).grid(row=numConditionLevels+1,column=int(maxLevelValues/2))
        tk.Button(buttonWindow, text="Back",command=lambda: master.switch_frame(columnLevelValuesPage,folderName)).grid(row=numConditionLevels+1,column=int(maxLevelValues/2)+1)
        tk.Button(buttonWindow, text="Quit",command=quit).grid(row=numConditionLevels+1,column=int(maxLevelValues/2)+2)
