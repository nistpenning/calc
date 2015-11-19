# -*- coding: utf-8 -*-
"""
Created on Wed May 06 14:31:01 2015

@author: jgb

Collection of functions for getting data from the files generated by hfGUI3
"""
import os, csv, itertools
import numpy as np

"""
Here I save a few long names from the ion properities for easy access later
"""
brightMean = "detection%det_brightMean"
darkMean = "detection%det_darkMean"
det_t = "detection%det_t"

def n_slice(name):  
    #function for getting a short relevant part of filename, since I can't remember the numbers
    return name[-20:-8]

def get_immediate_subdirectories(a_dir):
    return [name for name in os.listdir(a_dir)
            if os.path.isdir(os.path.join(a_dir, name))]

def get_gen_csv(first_name, skip_header=False):
    #Get the data file name
    file_name = False
    for file in os.listdir(os.getcwd()):
        if file.startswith(first_name) and file.endswith(".csv"):
            file_name = file
    if file_name is False:
        print("Did not find file")
        return 0
    else:
        #Get data from the file
        if skip_header is False:
            data = np.genfromtxt(file_name, delimiter=",", names=True, dtype=None)
        else: data = np.genfromtxt(file_name, delimiter=",",
                                   skip_header=skip_header,
                                   names=None,
                                   dtype=None,
                                   comments="#")
        return file_name, data

def get_ionProp_value(prop, where='props'):
    if where is 'props':
        os.chdir('props')
    file_name = False
    for file in os.listdir(os.getcwd()):
        if file.startswith('IonProperties') and file.endswith(".txt"):
            file_name = file
    if file_name is False:
        print("Did not find file")
        if where is 'props':
            os.chdir('..')
        return 0
    else:
        ionprop = file_name
        data = np.genfromtxt(ionprop, dtype='str', delimiter=' = ')
        ind = np.where(data == '{'+prop+'}')
        val_str = data[ind[0].astype(int)[0]][1]
        if where is 'props':
            os.chdir('..')
        return np.float(val_str[1:-1])

def get_ionProp_dict(keys, where='props'):
    if where is 'props':
        os.chdir('props')
    file_name = False
    for file in os.listdir(os.getcwd()):
        if file.startswith('IonProperties') and file.endswith(".txt"):
            file_name = file
    if file_name is False:
        print("Did not find file")
        if where is 'props':
            os.chdir('..')
        return 0
    else:
        ionprop = file_name
        data = np.genfromtxt(ionprop, dtype='str', delimiter=' = ')
        vals =[]
        for prop in keys:
            ind = np.where(data == '{'+prop+'}')
            val_str = data[ind[0].astype(int)[0]][1]
            vals.append(np.float(val_str[1:-1]))
        val_dict = dict(zip(keys,vals))
        if where is 'props':
            os.chdir('..')
        return val_dict

def parse_raw_counts(array):
    bad = 0
    verbose = False
    for x in np.nditer(array, op_flags=['readwrite']):
        if x == -1:
            print('Found bad data point')
            bad += 1
            x[...] = -1
        elif np.isnan(x) == True:
            print('Found bad data point Nan')
            bad += 1
            x[...] = -1
        else:
            x[...] = int(x) & 0x1fff
    if bad > 0:
        if verbose is True:
            print("# of bad points: {}".format(bad))
            print("removing all bad points in return value")
        array = array[array!=-1]
    return array

def parse_raw_counts_data(file_name):
    #Get data from the file
    with open(file_name) as f:
        reader = csv.reader(f)
        colnames = next(reader)

    colnames = [i for i in colnames if i] # gets rid of column names that are empty strings
    num_reg_columns = len(colnames) - 1
    data = np.genfromtxt(file_name, delimiter=",", names=True, dtype=None, usecols=range(num_reg_columns+1))
    non_hist_cols = num_reg_columns  # defined by HFGUI expt type
    
    counts_data = np.genfromtxt(file_name, delimiter=",", skip_header=1, dtype='float')
    if len(np.shape(counts_data))>1:
        num_cols_total = len(counts_data[0])
    else:
        num_cols_total = len(counts_data)
    counts_cols = range(num_reg_columns, num_cols_total)
    counts_data = np.genfromtxt(file_name, delimiter=",", skip_header=1, dtype='int', usecols=counts_cols)

    def parse_raw_counts(array):
        bad = 0
        for x in np.nditer(array, op_flags=['readwrite']):
            if x == -1:
                print('Found bad data point')
                bad += 1
                x[...] = -1
            else:
                x[...] = int(x) & 0x1fff
        if bad > 0:
            print("# of bad points: {}".format(bad))

    def make_detect_array(array):
        det_array = np.copy(array)
        for x in np.nditer(det_array, op_flags=['readwrite']):
            x[...] = ((int(x) & (0xf<<13))>>13)-1
        return det_array

    parse_raw_counts(counts_data)
    det_array = make_detect_array(counts_data)

    return colnames, data, counts_data, det_array

def get_raw_counts():
    file_name = False
    for file in os.listdir(os.getcwd()):
        if file.startswith("rawData.") and file.endswith(".csv"):
            file_name = file
    if file_name is False:
        print("Did not find file")
        return 0
    else:
        colnames, data, counts_data, det_array = parse_raw_counts_data(file_name)

        x = data[colnames[1]]
        scandata = data[colnames[2]]

        num_scans = len(scandata[scandata == scandata[0]])
        points_in_scan = np.size(scandata)/num_scans

        trials = np.size(counts_data,axis=1)*1.0

        # Make mask for bad points
        #avg_pmt_counts = np.mean(counts_data[counts_data!=0], axis=1)
        avg_pmt_counts = np.array([np.mean(i[i!=-1]) for i in counts_data ])
        #pmterr = np.std(counts_data[counts_data!=0], axis=1)
        pmterr = np.array([np.std(i[i!=-1]) for i in counts_data ])

        return file_name, scandata, avg_pmt_counts, pmterr, trials, data

def get_raw_counts_hist():
    file_name = False
    for file in os.listdir(os.getcwd()):
        if file.startswith("rawData.") and file.endswith(".csv"):
            file_name = file
    if file_name is False:
        print("Did not find file")
        return 0
    else:
        colnames, data, counts_data, det_array = parse_raw_counts_data(file_name)
        
        x = data[colnames[1]]
        scandata = data[colnames[2]]

        if len(np.shape(data))>1:

            num_scans = len(scandata[scandata == scandata[0]])
            points_in_scan = np.size(scandata)/num_scans
    
            trials = np.size(counts_data,axis=1)*1.0

        # Make mask for bad points
        #avg_pmt_counts = np.mean(counts_data[counts_data!=0], axis=1)
        avg_pmt_counts = np.array([np.mean(i[i!=-1]) for i in counts_data ])
        #pmterr = np.std(counts_data[counts_data!=0], axis=1)
        pmterr = np.array([np.std(i[i!=-1]) for i in counts_data ])

        return file_name, scandata, counts_data, data

def bright_fraction(counts, err=False):
    #Convert counts into fraction of ions bright, from IonProp Brightmean
    max_count = get_ionProp_value("detection%det_brightMean")
    min_count = get_ionProp_value("detection%det_darkMean")
    b_prob = (counts - min_count)/(float(max_count - min_count))
    if err is False:
        return b_prob
    else:
        b_prob_err = err/(float(max_count - min_count))
        return b_prob, b_prob_err

def get_histData(max_count=100, min_count=0):
    #Get the data file name
    file_name = False
    for file in os.listdir(os.getcwd()):
        if file.startswith("histData.") and file.endswith(".csv"):
            file_name = file
    if file_name is False:
        print("Did not find file")
        return 0
    else:
        #Get data from the file

        data = np.genfromtxt(file_name, delimiter=",", names=True, dtype=None)
        columns = data.dtype.names
        non_hist_cols = columns.index('hist00')  # defined by HFGUI expt type
        avg_pmt_col = columns.index('ave0')
        x, x_val, avg_pmt_counts = np.genfromtxt(file_name,
                                                 unpack=True,
                                                 skip_header=1,
                                                 dtype=None,
                                                 usecols=(1, 2, avg_pmt_col),
                                                 delimiter=',')
        xname = data.dtype.names[2]
        scandata = x_val

        num_scans = len(scandata[scandata == scandata[0]])
        points_in_scan = np.size(scandata)/num_scans

        if num_scans > 1:
            #have to average the data together
            scandata = scandata.reshape((num_scans, points_in_scan))
            scandata = np.mean(scandata, axis=0)
            avg_pmt_counts = avg_pmt_counts.reshape((num_scans, points_in_scan))
            avg_pmt_counts = np.mean(avg_pmt_counts, axis = 0)

    #get just the histogram data
        num_cols_total = len(data[0])
        num_hist_cols = num_cols_total - non_hist_cols-1
        hist_data_cols = np.arange(non_hist_cols+1, non_hist_cols+num_hist_cols+1)
        histdata = np.genfromtxt(file_name,delimiter=",",names=None,
                                 skip_header = 1,dtype=None,usecols=hist_data_cols)
        if num_scans > 1:
            #accumulate in bins from all data sets
            accu_hist_data = np.array([histdata[i:i+points_in_scan] for i in range(num_scans)])
            histdata = np.sum(accu_hist_data, axis=0)

    #calc error bars from histogram data
        counts_m = np.arange(num_hist_cols)
        counts = np.arange(1, num_hist_cols + 1)
        mom1 = (counts+counts_m)/2.0
        mom2 = (counts**2+counts_m*counts+counts_m**2)/3.0
        trials = np.sum(histdata[0], dtype=float)
        prob = histdata/trials

        hist_mean = np.array([np.dot(p, mom1) for p in prob])
        hist_2mean = np.array([np.dot(p, mom2) for p in prob])

        var = hist_2mean - hist_mean**2

        pmterr = np.sqrt(var)

    #scale avg_pmt_counts to make a bright state probability
        b_prob = (avg_pmt_counts - min_count)/(float(max_count - min_count))
        histextent = [-min_count/float(max_count-min_count),
                      (num_hist_cols-min_count)/float(max_count-min_count)]
        pmterr = (pmterr)/(float(max_count - min_count))
        return file_name, scandata, b_prob, pmterr, histdata, histextent, xname

def get_scanData(max_count=100, min_count=0):
    file_name, scandata, b_prob, pmterr, histdata, histextent, xname = get_histData(max_count=max_count, min_count=min_count)
    return file_name, scandata, b_prob, pmterr

def get_squeezeData(max_count=100, min_count=0):
    #Get the data file name
    file_name = False
    for file in os.listdir(os.getcwd()):
        if file.startswith("histData.") and file.endswith(".csv"):
            file_name = file
    if file_name is False:
        print("Did not find file")
        return 0
    else:
        #Get data from the file

        data = np.genfromtxt(file_name, delimiter=",", names=True, dtype=None)
        columns = data.dtype.names
        non_hist_cols = columns.index('hist00')  # defined by HFGUI expt type
        avg_pmt_col = columns.index('ave0')
        x, x_val, avg_pmt_counts = np.genfromtxt(file_name,
                                                 unpack=True,
                                                 skip_header=1,
                                                 dtype=None,
                                                 usecols=(1, 2, avg_pmt_col),
                                                 delimiter=',')
        xname = data.dtype.names[2]
        scandata = x_val

        num_scans = len(scandata[scandata == scandata[0]])
        points_in_scan = np.size(scandata)/num_scans

        if num_scans > 1:
            #have to average the data together
            scandata = scandata.reshape((num_scans, points_in_scan))
            scandata = np.mean(scandata, axis=0)
            avg_pmt_counts = avg_pmt_counts.reshape((num_scans, points_in_scan))
            avg_pmt_counts = np.mean(avg_pmt_counts, axis = 0)

    #get just the histogram data
        num_cols_total = len(data[0])
        num_hist_cols = num_cols_total - non_hist_cols-1
        hist_data_cols = np.arange(non_hist_cols+1, non_hist_cols+num_hist_cols+1)
        histdata = np.genfromtxt(file_name,delimiter=",",names=None,
                                 skip_header = 1,dtype=None,usecols=hist_data_cols)
        if num_scans > 1:
            #accumulate in bins from all data sets
            accu_hist_data = np.array([histdata[i:i+points_in_scan] for i in range(num_scans)])
            histdata = np.sum(accu_hist_data, axis=0)

    #calc error bars from histogram data
        counts_m = np.arange(num_hist_cols)
        counts = np.arange(1, num_hist_cols + 1)
        mom1 = (counts+counts_m)/2.0
        mom2 = (counts**2+counts_m*counts+counts_m**2)/3.0
        trials = np.sum(histdata[0], dtype=float)
        prob = histdata/trials

        hist_mean = np.array([np.dot(p, mom1) for p in prob])
        hist_2mean = np.array([np.dot(p, mom2) for p in prob])

        var = hist_2mean - hist_mean**2

        pmterr = np.sqrt(var)

    #scale avg_pmt_counts to make a bright state probability
        b_prob = (avg_pmt_counts - min_count)/(float(max_count - min_count))
        histextent = [-min_count/float(max_count-min_count),
                      (num_hist_cols-min_count)/float(max_count-min_count)]
        pmterr = (pmterr)/(float(max_count - min_count))
        return file_name, scandata, b_prob, pmterr, avg_pmt_counts, trials

# do a new defintion so I can just work on counts, not jz
def get_countsData():
    #Get the data file name
    file_name = False
    for file in os.listdir(os.getcwd()):
        if file.startswith("histData.") and file.endswith(".csv"):
            file_name = file
    if file_name is False:
        print("Did not find file")
        return 0
    else:
        #Get data from the file

        data = np.genfromtxt(file_name, delimiter=",", names=True, dtype=None)
        columns = data.dtype.names
        non_hist_cols = columns.index('hist00')  # defined by HFGUI expt type
        avg_pmt_col = columns.index('ave0')
        x, x_val, avg_pmt_counts = np.genfromtxt(file_name,
                                                 unpack=True,
                                                 skip_header=1,
                                                 dtype=None,
                                                 usecols=(1, 2, avg_pmt_col),
                                                 delimiter=',')

        scandata = x_val

        num_scans = len(scandata[scandata == scandata[0]])
        points_in_scan = np.size(scandata)/num_scans

        if num_scans > 1:
            #have to average the data together
            scandata = scandata.reshape((num_scans, points_in_scan))
            scandata = np.mean(scandata, axis=0)
            avg_pmt_counts = avg_pmt_counts.reshape((num_scans, points_in_scan))
            avg_pmt_counts = np.mean(avg_pmt_counts, axis = 0)

    #get just the histogram data
        num_cols_total = len(data[0])
        num_hist_cols = num_cols_total - non_hist_cols-1
        hist_data_cols = np.arange(non_hist_cols+1, non_hist_cols+num_hist_cols+1)
        histdata = np.genfromtxt(file_name,delimiter=",",names=None,
                                 skip_header = 1,dtype=None,usecols=hist_data_cols)
        if num_scans > 1:
            #accumulate in bins from all data sets
            accu_hist_data = np.array([histdata[i:i+points_in_scan] for i in range(num_scans)])
            histdata = np.sum(accu_hist_data, axis=0)

    #calc error bars from histogram data
        counts_m = np.arange(num_hist_cols)
        counts = np.arange(1, num_hist_cols + 1)
        mom1 = (counts+counts_m)/2.0
        mom2 = (counts**2+counts_m*counts+counts_m**2)/3.0
        trials = np.sum(histdata[0], dtype=float)
        prob = histdata/trials

        hist_mean = np.array([np.dot(p, mom1) for p in prob])
        hist_2mean = np.array([np.dot(p, mom2) for p in prob])

        var = hist_2mean - hist_mean**2

        pmterr = np.sqrt(var)

        return file_name, scandata, avg_pmt_counts, pmterr, trials

