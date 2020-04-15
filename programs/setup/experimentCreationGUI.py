#!/usr/bin/env python3
import sys,os,subprocess,pickle
from sys import platform as sys_pf
if sys_pf == 'darwin':
    import matplotlib
    matplotlib.use("TkAgg")
import tkinter as tk
import tkinter.ttk

class NewProjectWindow(tk.Frame):
    def __init__(self,master,backpage):
        tk.Frame.__init__(self, master)
        mainWindow = tk.Frame(self)
        mainWindow.pack(side=tk.TOP,padx=10)
        
        l1 = tk.Label(mainWindow, text="Researcher/Project Name: ")
        e1 = tk.Entry(mainWindow)

        l2 = tk.Label(mainWindow, text="Path to Researcher/Project folder (must start with a ~ or /): ")
        e2 = tk.Entry(mainWindow)
        
        def createProject():
            projectName = e1.get()
            rawPathName = e2.get()
            if rawPathName[-1] != '/':
                rawPathName+='/'
            if rawPathName.split('/')[-2] == projectName:
                rawPathName = rawPathName[:rawPathName[:-1].rfind('/')]+'/'
            if projectName not in os.listdir(rawPathName):
                subprocess.run(['mkdir',rawPathName+projectName])
            if 'pathDict.pkl' in os.listdir('misc'):
                pathDict = pickle.load(open('misc/pathDict.pkl','rb'))
            else:
                pathDict = {}
            if projectName in pathDict.keys():
                del pathDict[projectName]
            pathDict[projectName] = rawPathName
            with open('misc/pathDict.pkl','wb') as f:
                pickle.dump(pathDict,f)
            print('Researcher/Project created!')

        b = tk.Button(mainWindow,text='Create new researcher/project',command=lambda:createProject())
        
        l1.grid(row=0,column=0,sticky=tk.W)
        e1.grid(row=0,column=1,sticky=tk.W)
        l2.grid(row=1,column=0,sticky=tk.W)
        e2.grid(row=1,column=1,sticky=tk.W)
        b.grid(row=2,column=0,columnspan=2)

        buttonWindow = tk.Frame(self)
        buttonWindow.pack(side=tk.TOP)
        tk.Button(buttonWindow, text="Back",command=lambda: master.switch_frame(backpage)).pack(side=tk.LEFT)
        tk.Button(buttonWindow, text="Quit",command=lambda: quit()).pack(side=tk.LEFT)

class NewExperimentWindow(tk.Frame):
    def __init__(self,master,backpage):
        tk.Frame.__init__(self, master)
        mainWindow = tk.Frame(self)
        mainWindow.pack(side=tk.TOP,padx=10)
         
        pathDict = pickle.load(open('misc/pathDict.pkl','rb'))
        projects = list(pathDict.keys())
        projectTitle = tk.Label(mainWindow,text='Researcher/Project name: ')
        projectMenu = tkinter.ttk.Combobox(mainWindow,values = projects)
        projectMenu.width = len(max(projects,key=len))
        if len(projects) == 1:
            projectMenu.set(projectMenu['values'][0])

        l1 = tk.Label(mainWindow, text="Experiment Date (YYYYMMDD): ")
        e1 = tk.Entry(mainWindow)

        l2 = tk.Label(mainWindow, text="Experiment Name: ")
        e2 = tk.Entry(mainWindow)
        
        def createExperiment():
            experimentName = e1.get()+'-'+e2.get()
            projectName = projectMenu.get()
            pathName = pathDict[projectName]
            subprocess.run(['mkdir',pathName+projectName+'/'+experimentName])
            subfolders = ['inputData','outputData','plots','misc']
            subsubfoldersDict = {'inputData':['fcsFiles','singleCellCSVFiles','bulkCSVFiles'],'outputData':['excelFiles','pickleFiles','analysisFiles']}
            subsubsubfoldersDict = {'analysisFiles':['scaledData','reducedData','clusteredData','subsettedData']}
            for subfolder in subfolders:
                subprocess.run(['mkdir',pathName+projectName+'/'+experimentName+'/'+subfolder])
                if subfolder in subsubfoldersDict.keys():
                    subsubfolders = subsubfoldersDict[subfolder]
                    for subsubfolder in subsubfolders:
                        subprocess.run(['mkdir',pathName+projectName+'/'+experimentName+'/'+subfolder+'/'+subsubfolder])
                        if subsubfolder in subsubsubfoldersDict.keys():
                            subsubsubfolders = subsubsubfoldersDict[subsubfolder]
                            for subsubsubfolder in subsubsubfolders:
                                subprocess.run(['mkdir',pathName+projectName+'/'+experimentName+'/'+subfolder+'/'+subsubfolder+'/'+subsubsubfolder])
            
            print('Experiment created!')

        b = tk.Button(mainWindow,text='Create experiment',command=lambda:createExperiment())
        
        projectTitle.grid(row=0,column=0)
        projectMenu.grid(row=0,column=1)
        l1.grid(row=1,column=0,sticky=tk.W)
        e1.grid(row=1,column=1,sticky=tk.W)
        l2.grid(row=2,column=0,sticky=tk.W)
        e2.grid(row=2,column=1,sticky=tk.W)
        b.grid(row=3,column=0,columnspan=2)

        buttonWindow = tk.Frame(self)
        buttonWindow.pack(side=tk.TOP)
        tk.Button(buttonWindow, text="Back",command=lambda: master.switch_frame(backpage)).pack(side=tk.LEFT)
        tk.Button(buttonWindow, text="Quit",command=lambda: quit()).pack(side=tk.LEFT)
