# Used modules
import sys
import os
if not any(os.path.abspath('..') == p for p in sys.path): 
    print('adding comm module to path...')
    sys.path.insert(0, os.path.abspath('.'))
    print(os.path.abspath('.'))
import comm as comm

# Test 1
# Checking the default values
# Expectetd result: The data from cassette 1 should be readed.
test_1 = comm.instrument_control.set_attenuation_MTA_150()
print('<---Test 1--->')
print(test_1)
print('\n\n')

# Test 2
# Checking with two cassets
test_2 = comm.instrument_control.set_attenuation_MTA_150(cassetts=['1','2'],attenuations=[0,0], offsets=[0,0], wavelengths=[1500,1500])
print('<---Test 2--->')
print(test_2)
print('\n\n')

# Test 3
# Checking negativ offset
test_3 = comm.instrument_control.set_attenuation_MTA_150(cassetts=['1','2'],attenuations=[50,50], offsets=[-10,-10], wavelengths=[1500,1500])
print('<---Test 3--->')
print(test_3)
print('\n\n')

# Test 4
# Checking maximal attenuation
test_3 = comm.instrument_control.set_attenuation_MTA_150(cassetts=['1','2'],attenuations=[60,60], offsets=[60,60], wavelengths=[1500,1500])
print('<---Test 4--->')
print(test_3)
print('\n\n')