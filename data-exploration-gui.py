#!/usr/bin/env python3
import pickle,sys,os,json,math,subprocess,string
from sys import platform as sys_pf
if sys_pf == 'darwin':
    import matplotlib
    matplotlib.use("TkAgg")
import tkinter as tk
import tkinter.ttk
import pandas as pd
sys.path.insert(0, 'programs/setup')
from experimentCreationGUI import NewExperimentWindow,NewProjectWindow,RemoveProjectWindow
from experimentSetupGUI import ExperimentSetupStartPage
sys.path.insert(0, 'programs/dataprocessing')
from miscFunctions import setMaxWidth
from dataProcessingGUI import DataProcessingStartPage
sys.path.insert(0, 'programs/plotting')
from plottingGUI import PlotExperimentWindow 

#Root class; handles frame switching in gui
class MainApp(tk.Tk):
    def __init__(self):
        self.root = tk.Tk.__init__(self)
        self._frame = None
        self.homedirectory = os.getcwd()
        if self.homedirectory[-1] != '/':
            self.homedirectory+='/'
        self.switch_frame(ExperimentSelectionPage)

    def switch_frame(self, frame_class,*args):
        """Destroys current frame and replaces it with a new one."""
        new_frame = frame_class(self,*args)
        if self._frame is not None:
            self._frame.destroy()
        self._frame = new_frame
        self._frame.pack()

#Top level actions for experiments
class ExperimentSelectionPage(tk.Frame):
    def __init__(self,master):
        tk.Frame.__init__(self, master)
        mainWindow = tk.Frame(self)
        mainWindow.pack(side=tk.TOP,padx=10)
        
        l = tk.Label(mainWindow,text='Data Exploration GUI', font='Helvetica 18 bold')
        l.grid(row=0,column=0,columnspan=3)
        v = tk.StringVar(value='slt')
        rb1a = tk.Button(mainWindow, text="Create new project",padx = 20, command=lambda:master.switch_frame(NewProjectWindow,ExperimentSelectionPage))
        rb1b = tk.Button(mainWindow, text="Remove project",padx = 20, command=lambda:master.switch_frame(RemoveProjectWindow,ExperimentSelectionPage))
        rb1c = tk.Button(mainWindow, text="Associate new experiment with project",padx = 20, command=lambda:master.switch_frame(NewExperimentWindow,ExperimentSelectionPage))
        rb1a.grid(row=1,column=0,sticky=tk.EW,pady=(20,1))
        rb1b.grid(row=2,column=0,sticky=tk.EW,pady=1)
        rb1c.grid(row=3,column=0,sticky=tk.EW,pady=1)
        
        expSelectionFrame = tk.Frame(mainWindow,borderwidth=1,relief=tk.SOLID)
        expSelectionFrame.grid(row=4,column=0,pady=(10,0))
        rb1d = tk.Button(expSelectionFrame,text="Work on existing experiment ",padx = 20, command=lambda: collectInput())
        rb1d.grid(row=8,column=0,pady=(5,0),sticky=tk.EW)
        
        def getUpdateData(event):
            projectName = self.projectMenu.get()
            pathName = self.pathDict[projectName]
            experiments = []
            for experimentName in os.listdir(pathName+projectName):
                if '.DS' not in experimentName:
                    experiments.append(experimentName)
            experiments = sorted(experiments)[::-1]
            self.experimentMenu['values'] = experiments
            if len(experiments) == 1:
                self.experimentMenu.set(self.experimentMenu['values'][0])
            self.experimentMenu['width'] = len(max(experiments,key=len))
        
        if 'misc' not in os.listdir():
            subprocess.run(['mkdir','misc'])
        if 'pathDict.pkl' not in os.listdir(master.homedirectory+'misc'):
            self.pathDict = {}
        else:
            self.pathDict = pickle.load(open(master.homedirectory+'misc/pathDict.pkl','rb'))
        projects = list(self.pathDict.keys())
        self.projectMenu = tkinter.ttk.Combobox(expSelectionFrame,values = projects)
        if len(self.pathDict) > 0:
            self.projectMenu['width'] = len(max(projects,key=len))
        tk.Label(expSelectionFrame,text='Project name: ').grid(row=4,column=0)
        self.projectMenu.grid(row=5,column=0)
        #tk.Label(mainWindow,text='Project name: ').grid(row=4,column=1)
        #self.projectMenu.grid(row=4,column=2,sticky=tk.W)
        self.projectMenu.bind('<<ComboboxSelected>>', getUpdateData)

        self.experimentMenu = tkinter.ttk.Combobox(expSelectionFrame)
        tk.Label(expSelectionFrame,text='Experiment name: ').grid(row=6,column=0)
        self.experimentMenu.grid(row=7,column=0)
        #tk.Label(mainWindow,text='Experiment name: ').grid(row=5,column=1,sticky=tk.W)
        #self.experimentMenu.grid(row=5,column=2)

        def collectInput():
            projectName = self.projectMenu.get()
            pathName = self.pathDict[projectName]
            selectedExperiment = self.experimentMenu.get()
            os.chdir(pathName+projectName+'/'+selectedExperiment)
            master.switch_frame(ExperimentActionWindow,selectedExperiment)

        buttonWindow = tk.Frame(self)
        buttonWindow.pack(side=tk.TOP,padx=10,pady=(50,10))
        #tk.Button(buttonWindow, text="OK",command=lambda: collectInput()).pack(side=tk.LEFT)
        tk.Button(buttonWindow, text="Quit",command=lambda: quit()).pack(side=tk.LEFT)

#Top level actions for experiments
class ExperimentSelectionPage2(tk.Frame):
    def __init__(self,master):
        tk.Frame.__init__(self, master)
        mainWindow = tk.Frame(self)
        mainWindow.pack(side=tk.TOP,padx=10)
        
        l = tk.Label(mainWindow,text='Data Exploration GUI', font='Helvetica 18 bold')
        l.grid(row=0,column=0,columnspan=3)
        v = tk.StringVar(value='slt')
        rb1a = tk.Radiobutton(mainWindow, text="Create new project",padx = 20, variable=v, value='rp')
        rb1b = tk.Radiobutton(mainWindow, text="Remove project",padx = 20, variable=v, value='rmv')
        rb1c = tk.Radiobutton(mainWindow, text="Associate new experiment with project",padx = 20, variable=v, value='ce')
        rb1d = tk.Radiobutton(mainWindow,text="Work on existing experiment: ",padx = 20, variable=v, value='slt')
        rb1a.grid(row=1,column=0,sticky=tk.W)
        rb1b.grid(row=2,column=0,sticky=tk.W)
        rb1c.grid(row=3,column=0,sticky=tk.W)
        rb1d.grid(row=4,column=0,sticky=tk.W)
        
        def getUpdateData(event):
            projectName = self.projectMenu.get()
            pathName = self.pathDict[projectName]
            experiments = []
            for experimentName in os.listdir(pathName+projectName):
                if '.DS' not in experimentName:
                    experiments.append(experimentName)
            experiments = sorted(experiments)[::-1]
            self.experimentMenu['values'] = experiments
            if len(experiments) == 1:
                self.experimentMenu.set(self.experimentMenu['values'][0])
            self.experimentMenu['width'] = len(max(experiments,key=len))
        
        if 'misc' not in os.listdir():
            subprocess.run(['mkdir','misc'])
        if 'pathDict.pkl' not in os.listdir(master.homedirectory+'misc'):
            self.pathDict = {}
        else:
            self.pathDict = pickle.load(open(master.homedirectory+'misc/pathDict.pkl','rb'))
        projects = list(self.pathDict.keys())
        self.projectMenu = tkinter.ttk.Combobox(mainWindow,values = projects)
        if len(self.pathDict) > 0:
            self.projectMenu['width'] = len(max(projects,key=len))
        tk.Label(mainWindow,text='Project name: ').grid(row=5,column=0)
        self.projectMenu.grid(row=6,column=0)
        #tk.Label(mainWindow,text='Project name: ').grid(row=4,column=1)
        #self.projectMenu.grid(row=4,column=2,sticky=tk.W)
        self.projectMenu.bind('<<ComboboxSelected>>', getUpdateData)

        self.experimentMenu = tkinter.ttk.Combobox(mainWindow)
        tk.Label(mainWindow,text='Experiment name: ').grid(row=7,column=0)
        self.experimentMenu.grid(row=8,column=0)
        #tk.Label(mainWindow,text='Experiment name: ').grid(row=5,column=1,sticky=tk.W)
        #self.experimentMenu.grid(row=5,column=2)

        def collectInput():
            action = v.get()
            if action == 'rp':
                master.switch_frame(NewProjectWindow,ExperimentSelectionPage)
            elif action == 'rmv':
                master.switch_frame(RemoveProjectWindow,ExperimentSelectionPage)
            elif action == 'ce':
                master.switch_frame(NewExperimentWindow,ExperimentSelectionPage)
            elif action == 'slt':
                projectName = self.projectMenu.get()
                pathName = self.pathDict[projectName]
                selectedExperiment = self.experimentMenu.get()
                os.chdir(pathName+projectName+'/'+selectedExperiment)
                master.switch_frame(ExperimentActionWindow,selectedExperiment)

        buttonWindow = tk.Frame(self)
        buttonWindow.pack(side=tk.TOP,padx=10,pady=(50,10))
        tk.Button(buttonWindow, text="OK",command=lambda: collectInput()).pack(side=tk.LEFT)
        tk.Button(buttonWindow, text="Quit",command=lambda: quit()).pack(side=tk.LEFT)

class ExperimentActionWindow(tk.Frame):
    def __init__(self,master,selectedExperiment):
        tk.Frame.__init__(self, master)
        mainWindow = tk.Frame(self)
        mainWindow.pack(side=tk.TOP,padx=10)
        
        l1 = tk.Label(mainWindow,text='Choose an action:')
        v = tk.StringVar(value='plt')
        rb1a = tk.Radiobutton(mainWindow, text="Setup experiment",padx = 20, variable=v, value='se')
        rb1b = tk.Radiobutton(mainWindow,text="Process experiment",padx = 20, variable=v, value='pd')
        rb1c = tk.Radiobutton(mainWindow,text="Plot experiment",padx = 20, variable=v, value='plt')
        #rb1d = tk.Radiobutton(mainWindow,text="Analyze experiment",padx = 20, variable=v, value='ae')
        l1.grid(row=0,column=0)
        rb1a.grid(row=1,column=0,sticky=tk.W)
        rb1b.grid(row=2,column=0,sticky=tk.W)
        rb1c.grid(row=3,column=0,sticky=tk.W)
        #rb1d.grid(row=4,column=0,sticky=tk.W)
        
        def collectInput():
            action = v.get()
            if action == 'se':
                master.switch_frame(ExperimentSetupStartPage,selectedExperiment,ExperimentActionWindow)
            elif action == 'pd':
                master.switch_frame(DataProcessingStartPage,selectedExperiment,32,[],ExperimentActionWindow)
            elif action == 'plt':
                master.switch_frame(PlotExperimentWindow,selectedExperiment,ExperimentActionWindow)
        
        def backCommand():
            os.chdir(master.homedirectory)
            master.switch_frame(ExperimentSelectionPage)

        buttonWindow = tk.Frame(self)
        buttonWindow.pack(side=tk.TOP,padx=10,pady=10)
        tk.Button(buttonWindow, text="OK",command=lambda: collectInput()).pack(side=tk.LEFT)
        tk.Button(buttonWindow, text="Back",command=lambda: backCommand()).pack(side=tk.LEFT)
        tk.Button(buttonWindow, text="Quit",command=lambda: quit()).pack(side=tk.LEFT)

if __name__== "__main__":
    app = MainApp()
    app.mainloop()
