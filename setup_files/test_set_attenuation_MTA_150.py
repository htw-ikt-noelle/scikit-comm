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
# Expectetd result: The data from cassette 1 should be readed. The properties of the device should not be changed.
test_1 = comm.instrument_control.set_attenuation_MTA_150()
print('<---Test 1--->')
for cassette in test_1:
    print (cassette)
    for cassette_item in test_1[cassette]:
        print (cassette_item,':',test_1[cassette][cassette_item])
print('\n\n')

# Test 2
# Checking with two cassets
# Expecated result: Get the data from the cassette 1 and 2. The properties of the device should not be changed.
test_2 = comm.instrument_control.set_attenuation_MTA_150(cassetts=['1','2'])
print('<---Test 2--->')
for cassette in test_2:
    print (cassette)
    for cassette_item in test_2[cassette]:
        print (cassette_item,':',test_2[cassette][cassette_item])
print('\n\n')
3
# Test 3
# Change values of attenuation, offset and wavelength for cassette 2.
test_3 = comm.instrument_control.set_attenuation_MTA_150(cassetts=['2'],attenuations=[10], offsets=[5], wavelengths=[1300])
print('<---Test 3--->')
for cassette in test_3:
    print (cassette)
    for cassette_item in test_3[cassette]:
        print (cassette_item,':',test_3[cassette][cassette_item])
print('\n\n')

# Test 4
# Change only attenuation of cassette 2.
test_4 = comm.instrument_control.set_attenuation_MTA_150(cassetts=['2'],attenuations=[0])
print('<---Test 4--->')
for cassette in test_4:
    print (cassette)
    for cassette_item in test_4[cassette]:
        print (cassette_item,':',test_4[cassette][cassette_item])
print('\n\n')

# Test 5
# Change wavelength attenuation of cassette 2.
test_5 = comm.instrument_control.set_attenuation_MTA_150(cassetts=['2'],wavelength=[1350])
print('<---Test 5--->')
for cassette in test_5:
    print (cassette)
    for cassette_item in test_5[cassette]:
        print (cassette_item,':',test_5[cassette][cassette_item])
print('\n\n')

# Test 6
# Change only offset of cassette 2.
test_6 = comm.instrument_control.set_attenuation_MTA_150(cassetts=['2'],offset=[0])
print('<---Test 6--->')
for cassette in test_6:
    print (cassette)
    for cassette_item in test_6[cassette]:
        print (cassette_item,':',test_6[cassette][cassette_item])
print('\n\n')