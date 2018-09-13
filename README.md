# IWF-DNN
Python and Matlab code to produce MSc final project  ----- Deep neural network for MIMO BC wireless resource allocation.

All codes have been tested successfully on Matlab 2016b and Python 3.6.0 with Jupyter notebook and some package support, such as TensorFlow 1.8.0, Keras 2.2.0, h5py2.8.0, Numpy 1.14.5, pandas. 


Folder matlab_GenerateData is the Matlab files to generate dataset for training DNN model by using iterative water filling algorithm. 
In the folder: 
	matlab_GenerateData/generateMIMOBC.m: the function of generate samples, 
									    input: (number of transmit antenna, number of receive antenna, users, samples ),
									    output: generated dataset.
	matlab_GenerateData/changeComplex.m: transfer the complex number of the generated dataset to real part and imaginary part.
	matlab_GenerateData/testperformance2D.m: compare the sum rate performance between DNN model and IWF algorithm. 
	matlab_GenerateData/iterative_waterfill.m :implement iterative water-filling algorithm, which is also called on generateMIMOBC.m.  [1]
	matlab_GenerateData/waterfill_block.m: implement water-filling algorithm, which is callded on iterative_waterfill.m[1]
	matlab_GenerateData/MAC_to_BC.m: transform uplink covariance matrix to downlink covariance matrix.[1]

MIMO100000_3.mat: DNN training and validation dataset have 100000 samples (3 users, 2 transmit antennas, 2 receive antennas)
MIMO5000_3.mat: DNN test dataset have 5000 samples (3 users, 2 transmit antennas, 2 receive antennas)

Final2.ipynb: the ipython code for DNN model, and some functions about sum rate capacity evaluation method, sum power value, etc.
		     this file also show part of results of one DNN model.

model5layerslr001.json: store model structure
model5layerslr001.h5:    store model parameters 
				     these two files can be called on Matlab. 

errors.npz: training errors and validation errors in the training process. 


Reference
[1] N Jindal, W Rhee, S Vishwanath, S Jafar, A Goldsmith, ‘Sum Power Iterative Waterfilling Code’, (June 2005), http://people.ece.umn.edu/~nihar/iterative_wf_code.html 

version 1.0——August 2018.

Written by Ruoqi ZHAO (k1757932@kcl.ac.uk)

