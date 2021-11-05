# -*- coding: utf-8 -*-
"""
Created on Fri Nov  5 15:52:07 2021

@author: UET-Labor
"""


import sys
import os
if not any(os.path.abspath('..') == p for p in sys.path): 
    print('adding comm module to path...')
    sys.path.insert(0, os.path.abspath('..'))
import time
import numpy as np
import scipy.signal as ssignal
import matplotlib.pyplot as plt
import scipy.interpolate as sinterp
import comm as comm
import copy




#Test 1 
#Testing 
test_1 = comm.instrument_control.get_samples_HP_8153A_lightwave_multimeter()
print('<<-----TEST 1----->>')
print(test_1)
print('\n\n')

#Test 2
#Testing if user dont use the list for channel 
# test_2 = comm.instrument_control.get_samples_HP_8153A_lightwave_multimeter(channels = ('1','2'), GPIB_address='22',power_units = ['DBM','DBM'], wavelengths=[1500.0,1550.0] , log_mode = False)
# print('<<-----TEST 2----->>')
# print(test_2)
# print('\n\n')


#Test 3
#Testing if user dont use the string type to input the channel
# test_3 = comm.instrument_control.get_samples_HP_8153A_lightwave_multimeter(channels = [1,2], GPIB_address='22',power_units = ['DBM','DBM'], wavelengths=[1500.0,1550.0] , log_mode = False)
# print('<<-----TEST 3----->>')
# print(test_3)
# print('\n\n')



#Test 4
#Testing if user dont use GBIB number as string number 
# test_4 = comm.instrument_control.get_samples_HP_8153A_lightwave_multimeter(channels = ['1','2'], GPIB_address=22 ,power_units = ['DBM','DBM'], wavelengths=[1500.0,1550.0] , log_mode = False)
# print('<<-----TEST 4----->>')
# print(test_4)
# print('\n\n')



#Test 5
#Testing if user dont use the list for power_units 
# test_5 = comm.instrument_control.get_samples_HP_8153A_lightwave_multimeter(channels = ['1','2'], GPIB_address='22',power_units = ('DBM','DBM'), wavelengths=[1500.0,1550.0] , log_mode = False)
# print('<<-----TEST 5----->>')
# print(test_5)
# print('\n\n')



#Test 6
#Testing if user  use the wrong type in power unit 
# test_6 = comm.instrument_control.get_samples_HP_8153A_lightwave_multimeter(channels = ['1','2'], GPIB_address='22',power_units = ['HTZ','DBM'], wavelengths=[1500.0,1550.0] , log_mode = False)
# print('<<-----TEST 6----->>')
# print(test_6)
# print('\n\n')


#Test 7
#Testing if user dont use the float number for wavelengths
# test_7 = comm.instrument_control.get_samples_HP_8153A_lightwave_multimeter(channels = ['1','2'], GPIB_address='22',power_units = ['DBM','DBM'], wavelengths=[1500,1550] , log_mode = False)
# print('<<-----TEST 7----->>')
# print(test_7)
# print('\n\n')


#Test 8
#Testing if user dont use list for wavelengths
# test_8 = comm.instrument_control.get_samples_HP_8153A_lightwave_multimeter(channels = ['1','2'], GPIB_address='22',power_units = ['DBM','DBM'], wavelengths=(1500.0,1550.0) , log_mode = False)
# print('<<-----TEST 8----->>')
# print(test_8)
# print('\n\n')


#Test 9
#Testing if user dont use the boolean type for log_mode
# test_9 = comm.instrument_control.get_samples_HP_8153A_lightwave_multimeter(channels = ['1','2'], GPIB_address='22',power_units = ['DBM','DBM'], wavelengths=[1500.0,1550.0] , log_mode = 20)
# print('<<-----TEST 9----->>')
# print(test_9)
# print('\n\n')



#Test 10
#Testing if user put more channel 
# test_10 = comm.instrument_control.get_samples_HP_8153A_lightwave_multimeter(channels = ['1','2','3'], GPIB_address='22',power_units = ['DBM','DBM'], wavelengths=[1500.0,1550.0] , log_mode = False)
# print('<<-----TEST 10----->>')
# print(test_10)
# print('\n\n')
 

# #Test 11
# #Testing if user use the wrong name of channel
# test_11 = comm.instrument_control.get_samples_HP_8153A_lightwave_multimeter(channels = ['one','two'], GPIB_address='22',power_units = ['DBM','DBM'], wavelengths=[1500.0,1550.0] , log_mode = False)
# print('<<-----TEST 11----->>')
# print(test_11)
# print('\n\n')


#Test 12
#Testing if user use the wrong range of wavelengths
# test_12 = comm.instrument_control.get_samples_HP_8153A_lightwave_multimeter(channels = ['1','2'], GPIB_address='22',power_units = ['DBM','DBM'], wavelengths=[500.0,1850.0] , log_mode = False)
# print('<<-----TEST 12----->>')
# print(test_12)
# print('\n\n')



#Test 13
#Testing if user use the wrong range of wavelengths
test_13 = comm.instrument_control.get_samples_HP_8153A_lightwave_multimeter(channels = ['1','2'], GPIB_address='22',power_units = ['DBM','DBM'], wavelengths=[850.0,1900.0] , log_mode = False)
print('<<-----TEST 13----->>')
print(test_13)
print('\n\n')