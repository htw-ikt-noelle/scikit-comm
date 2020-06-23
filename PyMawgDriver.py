# -*- coding: utf-8 -*-
"""
Created on Mon Jun 22 20:21:08 2020

@author: Prabesh
"""
import numpy as np
import matplotlib.pyplot as plt

import sys
sys.path.append("D:/Python/pymeasure")       #from lab
#sys.path.append("C:/PythonSignal/pymeasure") #from home
import pymeasure
from pymeasure.instruments.agilent import Agilent33521A
from pymeasure.adapters import VISAAdapter  # importing Visaadapters from pymeasure

print(pymeasure.__version__)

adapter = VISAAdapter('TCPIP::192.168.1.44::INSTR') # interface for the adapter must be defined
awg = Agilent33521A(adapter) 
awg.id # returns the result of *IDN? command #awg.ask("*IDN?")
awg.write(':OUTPut1:STATe off')
#awg.write('APPLy:SQUare')
print(awg.id , ' The device Agilent33522 has been identified')

