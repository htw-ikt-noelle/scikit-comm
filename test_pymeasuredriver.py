# -*- coding: utf-8 -*-
"""
Created on Tue Jun 23 16:09:58 2020

@author: UET-Labor
"""
import numpy as np
import matplotlib.pyplot as plt

import sys
sys.path.append("D:/Python/pymeasure")       #from lab
#sys.path.append("C:/PythonSignal/pymeasure") #from home
import pymeasure
from pymeasure.instruments.agilent import Agilent33500
from pymeasure.adapters import VISAAdapter  # importing Visaadapters from pymeasure

print(pymeasure.__version__)

adapter = VISAAdapter('TCPIP::192.168.1.44::INSTR') # interface for the adapter must be defined
awg = Agilent33500(adapter) 

awg.id
"""awg.id
Out[11]: ['Agilent Technologies', '33522A', 'MY50002322', '5.02-1.19-2.00-58-00']"""

# Value of SQUARE is not in the discrete set ['SIN', 'SQU', 
#'TRI', 'RAMP', 'PULS', 'PRBS', 'NOIS', 'ARB', 'DC']

#awg.shape = "SQU" # or 'SQU'
# #awg.beep() # gives beep output
#awg.output = True # selects the channel 1 output True = on, False = off
awg.frequency = 100
#awg.amplitude = 10
awg.amplitude_unit = 'VPP' #VRMS, DBM
#awg.offset =2

awg.shape = (input("Which function would you like to enter?  :  "))

try: 
    functionstr = str(awg.shape)
    print("The function you chose is:", awg.shape)
    
except ValueError as e:
    e

else :
    print ( " Try Again ")
    awg.shape = (input("Which function would you like to enter?  :  "))


    

    
    
   

