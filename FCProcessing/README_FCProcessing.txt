4.9.2018, CG

This code was modified from original FCPROCESS_MSC.m code associated with Gordon et al., 2017, Neuron, to work for task data.

Different versions:
*task: task analysis - points to residual files
*filtered: filters motion parameters before calculating FD. 
	   Necessary for MSC03 and MSC10, who have high frequency noise in the motion parameters. Shouldn't hurt to apply to others.
* preGLM2: task analysis on NON-residual data
