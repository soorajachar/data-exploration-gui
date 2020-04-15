#!/usr/bin/env python3
import sys,os,subprocess
from sys import platform as sys_pf
if sys_pf == 'darwin':
    import matplotlib
    matplotlib.use("TkAgg")
import tkinter as tk
from tkinter import ttk
from scrollable import Scrollable

sys.path.insert(0, 'programs/setup')
sys.path.insert(0, 'programs/dataprocessing')
sys.path.insert(0, 'programs/plotting')

root = tk.Tk()

header = ttk.Frame(root)
body = ttk.Frame(root)
footer = ttk.Frame(root)
header.pack()
#body.pack()
body.pack(fill=tk.BOTH, expand=True)
footer.pack()

ttk.Label(header, text="The header").pack()
ttk.Label(footer, text="The Footer").pack()


scrollable_body = Scrollable(body, width=32)

for i in range(30):
    ttk.Button(scrollable_body, text="I'm a button in the scrollable frame").grid(row=i,column=0)

scrollable_body.update()

root.mainloop()
