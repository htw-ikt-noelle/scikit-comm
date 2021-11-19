# Used modules
import sys
import os
if not any(os.path.abspath('.') == p for p in sys.path): 
    print('adding comm module to path...')
    sys.path.insert(0, os.path.abspath('.'))
    print(os.path.abspath('.'))
import comm as comm

# # Functional control

# # Test 1
# # Checking the default values
# # Expectetd result: The data from cassette 1 should be readed. The properties of the device should not be changed.
# test_1 = comm.instrument_control.set_attenuation_MTA_150()
# print('<---Test 1--->')
# for cassette in test_1:
#     print (cassette)
#     for cassette_item in test_1[cassette]:
#         print (cassette_item,':',test_1[cassette][cassette_item])
# print('\n\n')

# # Test 2
# # Checking with two cassets
# # Expecated result: Get the data from the cassette 1 and 2. The properties of the device should not be changed.
# test_2 = comm.instrument_control.set_attenuation_MTA_150(cassetts=['1','2'])
# print('<---Test 2--->')
# for cassette in test_2:
#     print (cassette)
#     for cassette_item in test_2[cassette]:
#         print (cassette_item,':',test_2[cassette][cassette_item])
# print('\n\n')
# 3
# # Test 3
# # Change values of attenuation, offset and wavelength for cassette 2.
# test_3 = comm.instrument_control.set_attenuation_MTA_150(cassetts=['2'],attenuations=[5.0], offsets=[5.0], wavelengths=[1300.5])
# print('<---Test 3--->')
# for cassette in test_3:
#     print (cassette)
#     for cassette_item in test_3[cassette]:
#         print (cassette_item,':',test_3[cassette][cassette_item])
# print('\n\n')

# # Test 4
# # Change only attenuation of cassette 2.
# test_4 = comm.instrument_control.set_attenuation_MTA_150(cassetts=['2'],attenuations=[8.0])
# print('<---Test 4--->')
# for cassette in test_4:
#     print (cassette)
#     for cassette_item in test_4[cassette]:
#         print (cassette_item,':',test_4[cassette][cassette_item])
# print('\n\n')

# # Test 5
# # Change wavelength of cassette 2.
# test_5 = comm.instrument_control.set_attenuation_MTA_150(cassetts=['2'],wavelengths=[1400.0])
# print('<---Test 5--->')
# for cassette in test_5:
#     print (cassette)
#     for cassette_item in test_5[cassette]:
#         print (cassette_item,':',test_5[cassette][cassette_item])
# print('\n\n')

# # Test 6
# # Change only offset of cassette 2.
# test_6 = comm.instrument_control.set_attenuation_MTA_150(cassetts=['2'],offsets=[0.0])
# print('<---Test 6--->')
# for cassette in test_6:
#     print (cassette)
#     for cassette_item in test_6[cassette]:
#         print (cassette_item,':',test_6[cassette][cassette_item])
# print('\n\n')

# # Test 7
# # Change attenuation of cassette 1 and wavelength of cassette 2
# test_7 = comm.instrument_control.set_attenuation_MTA_150(cassetts=['1','2'],attenuations=[5.0,None],wavelengths=[None,1550.0])
# print('<---Test 7--->')
# for cassette in test_7:
#     print (cassette)
#     for cassette_item in test_7[cassette]:
#         print (cassette_item,':',test_7[cassette][cassette_item])
# print('\n\n')

# # Test 7
# # Change order of cassetts
# test_8 = comm.instrument_control.set_attenuation_MTA_150(cassetts=['2','1'],attenuations=[None,0.0],wavelengths=[1500.0,None])
# print('<---Test 8--->')
# for cassette in test_8:
#     print (cassette)
#     for cassette_item in test_8[cassette]:
#         print (cassette_item,':',test_8[cassette][cassette_item])
# print('\n\n')

# # Test 9
# # Getting values from cassette 1 and cassette 2
# test_9 = comm.instrument_control.set_attenuation_MTA_150(cassetts=['1','2'])
# print('<---Test 9--->')
# for cassette in test_9:
#     print (cassette)
#     for cassette_item in test_9[cassette]:
#         print (cassette_item,':',test_9[cassette][cassette_item])
# print('\n\n')


# Correct input control

# # Test 10
# # Casettes is not list
# comm.instrument_control.set_attenuation_MTA_150(cassetts='1')

# # Test 11
# # Item of cassette is not a string
# comm.instrument_control.set_attenuation_MTA_150(cassetts=[1,2])

# # Test 12
# # To less items for casetts
# comm.instrument_control.set_attenuation_MTA_150(cassetts=[])

# # Test 13
# # To many items for casetts
# comm.instrument_control.set_attenuation_MTA_150(cassetts=['1','2','3','4','5','6','7','8','9'])

# # Test 14 
# # Attenuations is not list
# comm.instrument_control.set_attenuation_MTA_150(attenuations=None)

# # Test 15
# # Items of attenuations are not float or None
# comm.instrument_control.set_attenuation_MTA_150(attenuations=[12])

# # Test 16
# # An attenuations value is to high
# comm.instrument_control.set_attenuation_MTA_150(attenuations=[65.0])

# # Test 17
# # Attenuations value is to low
# comm.instrument_control.set_attenuation_MTA_150(attenuations=[-5.0])

# # Test 18
# # Offsets is not list
# comm.instrument_control.set_attenuation_MTA_150(offsets=40)

# # Test 19
# # Items of Offsets are not float
# comm.instrument_control.set_attenuation_MTA_150(offsets=[10])

# # Test 20
# # An Offsets value is to high
# comm.instrument_control.set_attenuation_MTA_150(offsets=[70.0])

# # Test 21
# # An Offsets value is to low
# comm.instrument_control.set_attenuation_MTA_150(offsets=[-70.0])

# # Test 22
# # Wavelengths is not list
# comm.instrument_control.set_attenuation_MTA_150(wavelengths=1400.0)

# # Test 23
# # Items of Wavelengths are not float
# comm.instrument_control.set_attenuation_MTA_150(wavelengths=[1400])

# # Test 24
# # An Wavelengths value is to high
# comm.instrument_control.set_attenuation_MTA_150(wavelengths=[1800.0])

# # Test 25
# # An Wavelengths value is to low
# comm.instrument_control.set_attenuation_MTA_150(wavelengths=[1100.0])

# # Test 26
# # GPIB address is not a string
# comm.instrument_control.set_attenuation_MTA_150(GPIB_address=12)

# # Test 27
# # Wrong cassete number
# comm.instrument_control.set_attenuation_MTA_150(cassetts=['9'])

# Test 28
# Wrong cassete number
comm.instrument_control.set_attenuation_MTA_150(cassetts=['0'])

# Test 29
# Wrong GPIB address
comm.instrument_control.set_attenuation_MTA_150(GPIB_address='22')
