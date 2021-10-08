
import sys
import os
if not any(os.path.abspath('..') == p for p in sys.path): 
    print('adding comm module to path...')
    sys.path.insert(0, os.path.abspath('..'))
    print(os.path.abspath('..'))
import numpy as np  # Version 1.20.1
import matplotlib.pyplot as plt # Version 3.3.4
import comm as comm

wanted_trace = 'A'

##### Test with one trace
print('-----Test with one trace-----')
test_one_trace = comm.instrument_control.get_samples_HP_71450B_OSA(traces = [wanted_trace], GPIB_address='13', log_mode=False, single_sweep = False)
# Get trace data
test_one_trace_data = test_one_trace[wanted_trace]['Trace_data']
print('Size of trace data: '+ str(test_one_trace_data.size))
print('Data type of trace data: '+ str(type(test_one_trace_data)))
# Get Unit
test_one_trace_unit = test_one_trace[wanted_trace]['Unit']
print('Unit of trace data: ' + test_one_trace_unit)
# Get start WL
test_one_trace_Start_WL = test_one_trace[wanted_trace]['Start_WL']
print('Start wavelength: ' + str(test_one_trace_Start_WL)) 
# Get stop WL
test_one_trace_Stop_WL = test_one_trace[wanted_trace]['Stop_WL']
print('Stop wavelength: ' + str(test_one_trace_Stop_WL)) 
# Get wavelength vector
test_one_trace_WL_vector = test_one_trace[wanted_trace]['WL_Vector']
print('Size of wavelength vector: '+ str(test_one_trace_WL_vector.size))
print('Data type of wavelength vector: '+ str(type(test_one_trace_WL_vector)))
# Plot received spectrum
plt.figure(1)
plt.title('Test with one trace')
plt.xlabel('Wavelength in nm')
plt.ylabel('Amplitude in ' + test_one_trace_unit)
plt.plot(test_one_trace_WL_vector,test_one_trace_data)
plt.show()

##### Test with two traces
test_two_traces = comm.instrument_control.get_samples_HP_71450B_OSA(traces = ['A','B'], GPIB_address='13', log_mode=False)
print('-----Test with two traces-----')
for trace in test_two_traces:
    # Get trace data
    test_two_traces_data = test_two_traces[trace]['Trace_data']
    print('Size of trace data: '+ str(test_two_traces_data.size))
    print('Data type of trace data: '+ str(type(test_two_traces_data)))
    # Get Unit
    test_two_traces_unit = test_two_traces[trace]['Unit']
    print('Unit of trace data: ' + test_two_traces_unit)
    # Get start WL
    test_two_traces_Start_WL = test_two_traces[trace]['Start_WL']
    print('Start wavelength: ' + str(test_two_traces_Start_WL)) 
    # Get stop WL
    test_two_traces_Stop_WL = test_two_traces[trace]['Stop_WL']
    print('Stop wavelength: ' + str(test_two_traces_Stop_WL)) 
    # Get wavelength vector
    test_two_traces_WL_vector = test_two_traces[trace]['WL_Vector']
    print('Size of wavelength vector: '+ str(test_two_traces_WL_vector.size))
    print('Data type of wavelength vector: '+ str(type(test_two_traces_WL_vector)))
    # Plot received spectrum
    plt.figure(2)
    plt.title('Test with two traces')
    plt.xlabel('Wavelength in nm')
    plt.ylabel('Amplitude in ' + test_two_traces_unit)
    plt.plot(test_one_trace_WL_vector,test_two_traces_data)
plt.show()


##### Test with three traces
test_three_traces = comm.instrument_control.get_samples_HP_71450B_OSA(traces = ['A','B','C'], GPIB_address='13', log_mode=False)
print('-----Test with three traces-----')
for trace in test_three_traces:

    # Get trace data
    test_three_traces_data = test_three_traces[trace]['Trace_data']
    print('Size of trace data: '+ str(test_three_traces_data.size))
    print('Data type of trace data: '+ str(type(test_three_traces_data)))
    # Get Unit
    test_three_traces_unit = test_three_traces[trace]['Unit']
    print('Unit of trace data: ' + test_three_traces_unit)
    # Get start WL
    test_three_traces_Start_WL = test_three_traces[trace]['Start_WL']
    print('Start wavelength: ' + str(test_one_trace_Start_WL)) 
    # Get stop WL
    test_three_traces_Stop_WL = test_three_traces[trace]['Stop_WL']
    print('Stop wavelength: ' + str(test_one_trace_Stop_WL)) 
    # Get wavelength vector
    test_three_traces_WL_vector = test_three_traces[trace]['WL_Vector']
    print('Size of wavelength vector: '+ str(test_three_traces_WL_vector.size))
    print('Data type of wavelength vector: '+ str(type(test_three_traces_WL_vector)))
    # Plot received spectrum
    plt.figure(2)
    plt.title('Test with three traces')
    plt.xlabel('Wavelength in nm')
    plt.ylabel('Amplitude in ' + test_three_traces_unit)
    plt.plot(test_three_traces_WL_vector,test_three_traces_data)
plt.show()

# #### Test with to many traces
# print('Test with too many traces')
# test_to_many_traces = comm.instrument_control.get_samples_HP_71450B_OSA(traces = ['A','A','B','C'], GPIB_address='4', log_mode=False)


# #### Test with to less traces
# print('Test with too less traces')
# test_to_less_traces = comm.instrument_control.get_samples_HP_71450B_OSA(traces = [], GPIB_address='4', log_mode=False)


# #### Test with wrong labeling
# print('Test with wrong labeling 1')
# test_wrong_labeling_1 = comm.instrument_control.get_samples_HP_71450B_OSA(traces = ['D','B','C'], GPIB_address='4', log_mode=False)


# print('Test with wrong labeling 2')
# test_wrong_labeling_2 = comm.instrument_control.get_samples_HP_71450B_OSA(traces = ['AA','B','C'], GPIB_address='4', log_mode=False)


# #### Test worng GPIB address
# try:
#     test_three_traces = comm.instrument_control.get_samples_HP_71450B_OSA(traces = ['A','B','C'], GPIB_address='100', log_mode=False)
# except Exception as e:
#     print(e)