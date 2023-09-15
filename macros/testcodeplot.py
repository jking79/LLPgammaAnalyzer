import pandas as pd
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

# Read the text file into a Pandas DataFrame
#inputfile = 'time_codeing_test.log'
#inputfile = 'time_codeing_test_ampcut75.log'
#inputfile = 'time_codeing_test_ecut2.log'
#inputfile = 'time_codeing_test_resids_ecut2.log'
#inputfile = 'time_codeing_test_resids_ecut05.log'
#inputfile = 'time_codeing_test_reds_ampcut05.log'
#inputfile = 'time_codeing_test_reds_noampcut.log'
#inputfile = 'time_codeing_test_reds_ampcut0.log'
#inputfile = 'time_codeing_ns_test_reds_ampcut0.log'
#inputfile = 'time_codeing_v2_ns_test_reds_ampcut0.log'
#inputfile = 'time_codeing_ns_test_reds_ampcut10.log'
#inputfile = 'time_codeing_ns_test_reds_ampcut25.log'
#inputfile = 'time_codeing_ns_test_reds_ampcut75.log'
#inputfile = 'time_codeing_ns_test_reds_evts1k_ampcut75.log'
#inputfile = 'time_mycodeing_v2_ns_test_reds_ampcut0.log'
#inputfile = 'time_mycodeing_v2_ns_test_reds_ampcut10.log'
#inputfile = 'time_mycodeing_v2_ns_test_dt_ampcut10.log'
#inputfile = 'time_mycodeing_v2_ns_test_dt_ampcut0.log'
#inputfile = 'time_codeing_ns_test_difVencdif_ampcut0.log'
#inputfile = 'time_codeing_ns_test_difVencdif_ampcut10.log'
#inputfile = 'time_codeing_ns_test_resdif_ampcut0.log'
#inputfile = 'time_codeing_ns_test_resdif_ampcut10.log'
#inputfile = 'time_mycodeing_v2_ns_test_tvt_ampcut0.log'
#inputfile = 'time_mycodeing_v2_ns_test_tvt_ampcut10.log'
#inputfile = 'time_mycodeing_v2_ns_test_adjdif1_ampcut0.log'
#inputfile = 'time_mycodeing_v3_ns_test_adjtvt110_ampcut0.log'
#inputfile = 'time_mycodeing_v2_ns_test_adjdif2_ampcut0.log'
#inputfile = 'time_mycodeing_v3_ns_test_adjdif2_ampcut0.log'
#inputfile = 'time_mycodeing_v3_ns_test_adjtvt210_ampcut0.log'
#inputfile = 'time_mycodeing_v2_ns_test_adjdif3_ampcut0.log'
#inputfile = 'time_mycodeing_v3_ns_test_adjtvt2_ampcut0.log'

#inputfile = 'time_mycodeing_v3_ns_nic_resd_ampcut0.log'
#inputfile = 'time_mycodeing_v3_ns_cornic_resd_ampcut0.log'
#inputfile = 'time_mycodeing_v3_ns_jwknic_resd_ampcut0.log'
#inputfile = 'time_mycodeing_v3_ns_nicfixed_resd_ampcut0.log'

#inputfile = 'time_mycodeing_v3_ns_ampvct_ampcut0.log'
#inputfile = 'time_mycodeing_v3_ns_jwkamp_resd_ampcut0.log'#range 8 ns
#inputfile = 'time_mycodeing_v3_ns_jwk9amp_resd_ampcut0.log'
#inputfile = 'time_mycodeing_v3_ns_jwk10amp_resd_ampcut0.log'
#inputfile = 'time_mycodeing_v3_ns_jwk12amp_resd_ampcut0.log'

#inputfile = 'time_mycodeing_v3_ns_test_difdifvut_ampcut0.log'
#inputfile = 'time_mycodeing_v3_ns_test_difdifvut_ampcut1.log'
#inputfile = 'time_mycodeing_v3_ns_13_utvadjdif_ampcut2.log'
#inputfile = 'time_mycodeing_v3_ns_14_utvadjdif_ampcut2.log'
#inputfile = 'time_mycodeing_v3_ns_14m125_utvadjdif_ampcut5.log'

#inputfile = 'time_mycodeing_v3_ns_14e175_utvadjdif_ampcut0.log'
#inputfile = 'time_mycodeing_v3_ns_14e05_utvadjdif_ampcut0.log'
#inputfile = 'time_mycodeing_v3_ns_14e025_utvadjdif_ampcut0.log'
#inputfile = 'time_mycodeing_v3_ns_14e01_utvadjdif_ampcut0.log'
#inputfile = 'time_mycodeing_v3_ns_14e005_utvadjdif_ampcut0.log'
#inputfile = 'time_mycodeing_v3_ns_16e005_utvadjdif_ampcut0.log'
#inputfile = 'time_mycodeing_v3_ns_20e005_utvadjdif_ampcut0.log'
#inputfile = 'time_mycodeing_v3_ns_20e025_utvadjdif_ampcut0.log'
#inputfile = 'time_mycodeing_v3_ns_18e025_utvadjdif_ampcut0.log'
#inputfile = 'time_mycodeing_v3_ns_18e020_utvadjdif_ampcut0.log'
#inputfile = 'time_mycodeing_v3_ns_18e020b20_utvadjdif_ampcut2.log'
#inputfile = 'time_mycodeing_v3_ns_18e020b20_utvadjdif_ampcut4.log'
#inputfile = 'time_mycodeing_v3_ns_jwkslopejitter2_adjdiff_ampcut0.log'

inputfile = 'koot_encoding_log.txt'

#data = pd.read_csv( inputfile, header=None, names=['Values'] )
#data = pd.read_csv( inputfile, sep=' ', header=None, names=['x', 'y']) #, low_memory=False)
data = pd.read_csv( inputfile, sep=' ', header=None, names=['jit', 'ncjit', 'atime', 'koot', 'itc', 'osc' ])

#m, b = np.polyfit(data['x'], data['y'], 1)
#print( 'Fit m: ', m, ' b: ', b )

#max_val = 12.0
#min_val = -12.0
#max_val = 0.1
#min_val = 0.0
#mask = ( data['y'] < min_val ) | ( data['y'] > max_val )
#base_cnt = len(data)
#mask_cnt = len(data[mask])

#print( "Acceptance : ", mask_cnt, base_cnt, 1-float(mask_cnt)/base_cnt )

# Set the range for the x-axis
#x_min = -0.045
#x_max = 0.045
#range_vals = (x_min, x_max)

# define the bin range for the x and y dimensions
#x_range = [-26, 26 ]
#x_bin = 520
#x_range = [-15, 15 ]
#x_bin = 300
#x_range = [-10, 10 ]
#x_bin = 200
#x_label = 'nocorrtime [ns]'
#x_label = 'corrtime [ns]'
#x_label = 'dt_(time-nocorrtime) [ns]'
#y_range = [-0.04, 0.04 ]
#y_bin = 160
#y_range = [0.05,0.1 ]
#y_bin = 500
#y_range = [-0.05, 0.05 ]
#y_bin = 100
#y_range = [-0.1, 0.1 ]
#y_bin = 200
#y_range = [-0.25, 0.25 ]
#y_bin = 500
#y_range = [-1, 1 ]
#y_bin = 200
#y_range = [-10, 10 ]
#y_bin = 200
#y_range = [-0.20, 0.20 ]
#y_bin = 40
#y_range = [-15, 15 ]
#y_bin = 300
#y_range = [-46, 46 ]
#y_bin = 920
#y_range = [-10, 10 ]
#y_bin = 200
#y_range = [-7.5, 7.5 ]
#y_bin = 150
#y_range = [-5, 5 ]
#y_bin = 100
#y_range = [-26, 26 ]
#y_bin = 520
#y_range = [0, 200 ]
#y_bin = 200
#y_range = [0, 25 ]
#y_bin = 25
#y_label = 'residuals [ns]'
#y_label = 'nocorrtime [ns]'
#y_label = 'adjnocorrtime [ns]'
#y_label = 'dt_(time-nocorrtime) [ns]'
#y_label = 'dt_(adjtime-nocorrtime) [ns]'
#y_label = 'dt_(time-encnocorrtime) [ns]'
#y_label = 'dt residuals [ns]'
#y_label = 'Amplitude'

x_range = [-26, 26 ]
x_bin = 520
x_label = 'nocorredtime [ns]'
y_range = [-26, 26 ]
y_bin = 520
y_label = 'dt_(adjtime-nocorrtime) [ns]'

# Plot a histogram of the data using Pandas
#data.hist(bins=1000)
#data.plot.hist(bins=225, range=range_vals, log=True)
#histogram = data.plot(kind='scatter', x='x', y='y', bins=(255,100), range=[x_range, y_range], cmap='Blues')
#hist, xedges, yedges = 
#plt.hist2d(data['x'], data['y'], bins=[x_bin,y_bin], range=[x_range, y_range], cmap='winter', norm=matplotlib.colors.LogNorm() )
#plt.hist2d(data['x'], data['y'], bins=[x_bin,y_bin], range=[x_range, y_range], cmap='winter')

kootmask = data['koot'] > 0
data_masked = data[kootmask]
plt.hist( data_masked['atime'], bins=x_bin, range=x_range )

# set the axis labels and title
#histogram.set_xlabel('X-axis')
#histogram.set_ylabel('Y-axis')
#histogram.set_title('2D Histogram')
plt.grid(True)
#plt.colorbar()
plt.xlabel(x_label)
#plt.ylabel(y_label)
#plt.title('1.8xE0.20-2.0 RH Amp > 2')
plt.title('kOOT True')
#plt.title('1.4x RH Amp > 2')
#plt.title('ct v amplitude')
#plt.title('jwkamp r12 tuc v residual')
#plt.title('nicola tuc v residual')

# Add labels and title to the plot
#plt.xlabel('$\delta_(t-uncorrt)$')
#plt.ylabel('Frequency')
#plt.title('Diff in encoded time v orignal time')

# Display the plot
plt.show()
