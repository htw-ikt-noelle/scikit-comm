import sys
import os
if not any(os.path.abspath('..') == p for p in sys.path): 
    print('adding comm module to path...')
    sys.path.insert(0, os.path.abspath('.'))
    print(os.path.abspath('.'))
import numpy as np  # Version 1.20.1
import matplotlib.pyplot as plt # Version 3.3.4
import comm as comm

##### Test with one trace
print('-----Test with one trace-----')
test_one_trace = comm.instrument_control.get_samples_HP_71450B_OSA(traces = ['A'], GPIB_address='4', log_mode=False)
# Get trace data
test_one_trace_data = test_one_trace['A']['trace_data']
print('Size of trace data: '+ str(test_one_trace_data.size))
print('Data type of trace data: '+ str(type(test_one_trace_data)))
# Get Unit
test_one_trace_unit = test_one_trace['A']['Unit']
print('Unit of trace data: ' + test_one_trace_unit)
# Get start WL
test_one_trace_Start_WL = test_one_trace['A']['Start_WL']
print('Start wavelength: ' + str(test_one_trace_Start_WL)) 
# Get stop WL
test_one_trace_Stop_WL = test_one_trace['A']['Start_WL']
print('Stop wavelength: ' + str(test_one_trace_Stop_WL)) 
# Get wavelength vector
test_one_trace_WL_vector = test_one_trace['A']['WL_Vector']
print('Size of wavelength vector: '+ str(test_one_trace_WL_vector.size))
print('Data type of wavelength vector: '+ str(type(test_one_trace_WL_vector)))
# Plot received spectrum
plt.figure(1)
plt.title('Test with one trace')
plt.xlabel('Wavelength in nm')
plt.ylabel('Amplitude in ' + test_one_trace_unit)
plt.plot(test_one_trace_WL_vector,test_one_trace_data)

##### Test with two traces
test_two_traces = comm.instrument_control.get_samples_HP_71450B_OSA(traces = ['A','B'], GPIB_address='4', log_mode=False)


##### Test with three traces
test_three_traces = comm.instrument_control.get_samples_HP_71450B_OSA(traces = ['A','B','C'], GPIB_address='4', log_mode=False)
