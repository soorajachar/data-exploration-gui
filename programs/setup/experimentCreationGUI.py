#!/usr/bin/env python3
import sys,os,subprocess,pickle,math
from sys import platform as sys_pf
if sys_pf == 'darwin':
    import matplotlib
    matplotlib.use("TkAgg")
import tkinter as tk
import tkinter.ttk
import tkinter.font as tkfont

def setMaxWidth(stringList, element):
    f = tkfont.nametofont(element.cget("font"))
    zerowidth=f.measure("M")
    w=max([f.measure(i) for i in stringList])/zerowidth
    element.config(width=max([6,math.ceil(1.5*w)]))

class NewProjectWindow(tk.Frame):
    def __init__(self,master,backpage):
        tk.Frame.__init__(self, master)
        mainWindow = tk.Frame(self)
        mainWindow.pack(side=tk.TOP,padx=10)
        
        l1 = tk.Label(mainWindow, text="Project Name: ")
        e1 = tk.Entry(mainWindow)

        l2 = tk.Label(mainWindow, text="Path to Project folder (must start with a ~ or /)")
        e2 = tk.Entry(mainWindow)
        
        l3 = tk.Label(mainWindow, text="For reference, path to current folder is: \n"+os.getcwd()+'/\n')
        
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
            print('Project created!')

        b = tk.Button(mainWindow,text='Create new project',command=lambda:createProject())
        
        l1.grid(row=0,column=0,sticky=tk.W)
        e1.grid(row=0,column=1,sticky=tk.W)
        l2.grid(row=1,column=0,sticky=tk.W)
        e2.grid(row=1,column=1,sticky=tk.W)
        l3.grid(row=3,column=0,columnspan=2)
        b.grid(row=2,column=0,columnspan=2)

        buttonWindow = tk.Frame(self)
        buttonWindow.pack(side=tk.TOP)
        tk.Button(buttonWindow, text="Back",command=lambda: master.switch_frame(backpage)).pack(side=tk.LEFT)
        tk.Button(buttonWindow, text="Quit",command=lambda: quit()).pack(side=tk.LEFT)

class RemoveProjectWindow(tk.Frame):
    def __init__(self,master,backpage):
        tk.Frame.__init__(self, master)
        mainWindow = tk.Frame(self)
        mainWindow.pack(side=tk.TOP,padx=10)
        
        l1 = tk.Label(mainWindow, text="Project name to remove (will not remove files, only path): ")
        pathDict = pickle.load(open('misc/pathDict.pkl','rb'))
        projects = list(pathDict.keys())
        projectVar = tk.StringVar()
        projectMenu = tk.OptionMenu(mainWindow,projectVar,*projects)
        setMaxWidth(projects,projectMenu)

        def removeProject():
            projectName = projectVar.get()
            del pathDict[projectName]
            with open('misc/pathDict.pkl','wb') as f:
                pickle.dump(pathDict,f)
            print('Project removed!')

        b = tk.Button(mainWindow,text='REMOVE project',command=lambda:removeProject())
        
        l1.grid(row=0,column=0,sticky=tk.W)
        projectMenu.grid(row=0,column=1)
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
        projectTitle = tk.Label(mainWindow,text='Project name: ')
        projectMenu = tkinter.ttk.Combobox(mainWindow,values = projects)
        projectMenu.width = len(max(projects,key=len))
        if len(projects) == 1:
            projectMenu.set(projectMenu['values'][0])

        l1 = tk.Label(mainWindow, text="Experiment Date (YYYYMMDD): ")
        e1 = tk.Entry(mainWindow)

        l2 = tk.Label(mainWindow, text="Experiment Name: ")
        e2 = tk.Entry(mainWindow)
        
        def createExperiment():
            experimentName = e2.get()
            amendedExperimentName = experimentName.replace('-','_')
            amendedExperimentName = amendedExperimentName.replace('/','_')
            amendedExperimentName = amendedExperimentName.replace(' ','_')
            experimentName = e1.get()+'-'+amendedExperimentName
            projectName = projectMenu.get()
            pathName = pathDict[projectName]
            subprocess.run(['mkdir',pathName+projectName+'/'+experimentName])
            subfolders = ['inputData','outputData','plots','misc']
            subsubfoldersDict = {'inputData':['fcsFiles','singleCellCSVFiles','bulkCSVFiles'],'outputData':['excelFiles','pickleFiles','analysisFiles']}
            subsubsubfoldersDict = {'analysisFiles':['scaledData','reducedData','clusteredData','subsettedData','clusterFrequencyData']}
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
