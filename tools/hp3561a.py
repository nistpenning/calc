import struct
import serial
import sys
import time
import os.path

import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime as dt

"""
Created on Wed Jun 11 16:38:51 2014

@author: Justin Bohnet

hp3561a.py is a script designed to pull off the current MAGNITUDE trace
off the hp signal analyzer and process the data for analysis and display
in other programs.

This script requires the use of the Prologix GPIB to USB converter
to handle GPIB communication protocol.

Inputs, given as argument to running the script:
1: COM port
2: GPIB address of hp3561a
3: 'fileID'

Outputs:
1: raw binary and ASCII data files at the location of the script
2: a processed data file 'fileID_data_XX' where XX is current time minutes
    this data file appears in the folder C:\data\.  The script prompts you to
    create this folder if it does not exist
"""


class HP3561A:
    """This class is an interface with the HP3561A dynamic signal analyzer. It makes
    use of the Prologix USB-GPIB adapter. This adapter provides a simple serial (USB)
    interface to GPIB.

    Example Use
        del(sa)  # in case it already exists
        sa = HP3561A(comport=4, gpibAddr=11, fileName='out.dat')
        (traceData, dispData) = sa.getTrace()
    """

    def __init__(self, comport=4, gpibAddr=11, fileName='out_hp3561a.dat'):

        # some private members
        self.__comport = comport
        self.__gpibAddr = gpibAddr
        self.__fileID = fileName
        self.__gpibRwDelay = 0.25  # is time between reads and writes

        # establish serial connection with the Prologix controller
        # note that pySerial starts counting serial ports at 0
        self.ser = serial.Serial( self.__comport-1, 9600, timeout=0.5 )

        # check that serial communication with Prologix USB-GPIB
        # adapter is working OK
        self.prologix_check_version()
        self.gpib_device_identity_check()

    def __del__(self):
        self.ser.close()

    def prologix_config(self):
        '''This function estabilishes GPIB with the Prologix controller
            all commands proceeded by ++ talk to the controller'''
        self.ser.flush()
        self.ser.write('++mode 1\n')
        self.ser.write('++auto 1\n')
        s = '+addr %d \n' % self.__gpibAddr
        self.ser.write(s)

    def prologix_check_version(self):
        '''Communicate with prologix USB-GPIB adapter. Get hardware version number.'''
        expectedResponse = 'Prologix GPIB-USB Controller version 6.101\r\n'
        r = self.gpib_send('++ver')
        if r != expectedResponse:
            print('ERROR :: hp3561a.prologixCheckVersion() did not get expected response. Is this really a Prologix GPIB-USB adapter?')
            print('ERROR :: hp3561a.prologixCheckVersion() response to ++ver was {}'.format(r))
            sys.exit()
        else:
            print(r)
        return

    def gpib_device_identity_check(self):
        '''Get GPIB device identity string.
        This makes sure we're actually talking to the HP3561. '''
        self.prologix_config()

        expected_response = 'HP3561A\r\n'
        r = self.gpib_send('ID?')
        if r != expected_response:
            print('ERROR :: hp3561a.gpibDeviceIdentityCheck() '
                  'did not get expected response. Is this really a HP 3561A?')
            print('ERROR :: hp3561a.gpibDeviceIdentityCheck() '
                  'response to ++ver was {}'.format(r))
            sys.exit()
        else:
            print(r)

    def config_hp3561a(self, mode='psd'):
        '''This is the default configuration for the HP3561A which
        is suitable for making particular types of measurements. Options
        are:
            mode = 'psd' :: power spectral density
        '''
        if mode == 'psd':
            # Display :: Format :: Single trace display
            self.gpib_send('SNGL;')

    def gpib_send(self,myStr):
        '''Routine for interacting with the GPIB device. It sends strings
        and listens for a response terminated by a carriage return.

        Example:
        gpibSend('SP10KHZ;') sets the span to 10 kHz
        '''
        s = myStr + '\n'
        try:
            time.sleep(self.__gpibRwDelay)
            self.ser.flush()
            self.ser.write(s)
            time.sleep(self.__gpibRwDelay)
            result = self.ser.read(1028*10)
        except serial.SerialException:
            print("serial.SerialException")
        except KeyboardInterrupt:
            print("KeyboardInterrupt")
            self.ser.close()
        return result

    def get_trace(self):
        '''(traceData,displayData) = getTrace() \n
        Get trace data from the HP 3561A. This is modeled on Example
        4-11 in the User Manual.'''

        # tell 3561A to dump active trace and header data
        # trace files are always 1028 bytes
        traceDataRaw = self.gpib_send('DSTB;')
        print('read trace date of length {}'.formatlen(traceDataRaw))

        # tell 3561A to dump current trace display info
        displayInfoRaw = self.gpib_send('DDSA;')
        print('display info raw length {}'.format(len(displayInfoRaw)))
        return (traceDataRaw, displayInfoRaw)


    def process_trace_data(self, traceData):

        #process the trace data which has the following format
        # byte 1: #
        # byte 2: A
        # byte 3: LB1 is MSB of total number of bytes transmitted
        # byte 4: LB2 is LSB of total number of bytes transmitted
        # byte 5...5+N: data bytes
        # byte 5+N+1: EOI
        # NOTE that the data is 2 bytes per integer; big-endian order

        #first two bytes are expected to be #A
        expectedValue = '#A'
        if traceData[0:2] != expectedValue :
            print('ERROR :: hp3561a.processTraceData() :: unexpected preamble')
        print('hp3561a.processTraceData() :: read {} bytes'.format(len(traceData)))

        #read 2 byte words in big-endian order
        pattern = '>'+'h'*400
        data = struct.unpack(pattern, traceData[4:4+2*400])
        data = np.array(data) #unpack returns a tuple; this converts to a list
        magdata = 0.005*data
        return magdata


    def process_display_data(self, line):
        # process the screen data to get frequency inf
            line = line[4:-1]  # strip off meaningless characters at beginning

            aline_full = line.split('   ')   # create array of info for access
            line = line.replace('   ','')  # make line more readable for footer

            # remove empty elements
            aline = [x.strip() for x in aline_full if x.strip() != '']
            aline = [x.replace(' ','') for x in aline] #remove spaces

            # if 'OVLD' is present, it messes up the elements for readout

            if aline[5] == 'OVLD':
                aline = np.delete(aline,5)
                line = line + ' OVLD' # note an overload condition for footer

            start = float(''.join([x for x in aline[-5] if x.isdigit()]))
            stop = float(''.join([x for x in aline[-3] if x.isdigit()]))

    def alltherest(self):
        # process the screen data to get frequency info
        with open("plot_settings.txt", "r") as f:
            line = f.readline()
        f.close()
        line = line[4:-1]  # strip off meaningless characters at beginning

        aline_full = line.split('   ')   # create array of info for access
        line = line.replace('   ', '')  # make line more readable for footer

        # remove empty elements
        aline = [x.strip() for x in aline_full if x.strip() != '']
        aline = [x.replace(' ','') for x in aline] #remove spaces

        # if 'OVLD' is present, it messes up the elements for readout

        if aline[5] == 'OVLD':
            aline = np.delete(aline,5)
            line += ' OVLD'  # note an overload condition for footer

        start = float(''.join([x for x in aline[-5] if x.isdigit()]))
        stop = float(''.join([x for x in aline[-3] if x.isdigit()]))
        numpts = magdata.size
        freq = np.linspace(start,stop,num=numpts)

        dataout = np.array([freq,magdata]).transpose()

        # get time for fileID
        minute = str(dt.now().minute)

        # note that BW gives the effective noise bandwidth ENBW used to convert
        # to a power spectral density if desired

        if os.path.exists("C:\\data\\"):
            np.savetxt('C:\\data\\' + fileID +'_data_' + minute + '.txt', dataout, delimiter = ',',
                       header = 'Freq, '+ aline[6], footer = line )
        else:
            print("Need to create C:\\data\\")

    @staticmethod
    def dbv_to_vsd(dbv, gain, bw):
        """ Regardless UNITS mode on HP3146A, the output to file is always
        in dBv/Hz.

        :param dbv: 'dbv' value from display of HP3146A.
        :param gain: gain from other instruments
        :param bw: bandwidth of measurement (see last line of data file)
        :return: dbV (V^2/Hz)
        """
        return (np.sqrt(2)*10**(dbv/20)/np.sqrt(bw)/gain)**2

    @staticmethod
    def dbv_to_v(dbv):
        """ UNITS mode on HP3146A is set to 'VOLT(dbV)'.

        using: dB_V = 20 Log_10(V_rms/1Vrms)

        :param dbv: amplitude indicated on display of HP3146A
        :return: volts amplitude
        """
        return np.sqrt(2)*10**(dbv/20)



