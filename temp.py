#! /usr/bin/env python3
import pickle

path = '/Users/acharsr/Documents/allExperiments/singleCellDataGUI_Testing/20200822-BALL_RUN2/misc/'

df = pickle.load(open(path+'tubeLayout-20200822-BALL_RUN2-cell.pkl','rb'))
print(df)
