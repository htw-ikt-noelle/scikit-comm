# -*- coding: utf-8 -*-
"""
Created on Tue Jun 23 14:48:32 2020

@author: UET-Labor
"""


import visa
rm = visa.ResourceManager('@py')
address = '192.168.1.44'
awg = rm.open_resource('TCPIP::' + address + '::INSTR')
print(awg.query('*IDN?'))

awg.write('OUTPUT OFF')
awg.write('SOUR:APPL:PRBS')
# APPLy:ARBitrary  APPLy:DC
# APPLy:NOISe 

# APPLy:PRBS 

# APPLy:PULSe 

# APPLy:RAMP 

# APPLy:SINusoid 

# APPLy:SQUare 

# VOLTage:UNIT 

awg.write('SOUR:FREQ 10') # source1/2 [SOURce[1|2]:]FREQuency {<frequency>|MINimum|MAXimum}

awg.write('SOUR:VOLT 2') # [SOURce[1|2]:]VOLTage {<amplitude>|MINimum|MAXimum}
# ampuntit VPP by default

del rm
del awg

