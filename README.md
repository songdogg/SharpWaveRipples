# SharpWaveRipples
Source code for custom MATLAB scripts used in the undergraduate dissertation project **"The Properties and Behavioural Correlates of Human Hippocampal Sharp-Wave Ripples During an Associative Memory Task"**.

# SharpWaveRipples
Source code for custom MATLAB scripts used in the the undergraduate dissertation project **"The Properties and Behavioural Correlates of Human Hippocampal Sharp-Wave Ripples During an Associative Memory Task"**.

In addition to the following custom MATLAB scripts, the NPMK and shadedErrorBar toolboxes, and a custom script for Morlet wavelet convolution were also used.

1. ripple_analysis.m              
--- Function to import raw .ns6 intracranial EEG data into MATLAB and identify sharp-wave ripples (SWRs). (See Section 2.4)
2. wrapper_ripple_analysis.m      
--- Wrapper script for ripple_analysis.m

3. spike_analysis.m               
--- Function to analyse single neuron spiking data. (See Section 2.5)
4. task_analysis.m                
--- Function to analyse SWRs identified by ripple_analysis.m based on the portion of the associative memory task they occur in. (See section 2.6)
5. wrapper_task_spike_analysis.m  
--- Wrapper script for task_analysis.m and spike_analysis.m

6. plot_ripples_LFP_method.m      
--- Script for plotting figures for various stages of iEEG processing and SWR identification (Figure 3).
7. plot_ripples.m                 
--- Script for plotting SWR characteristics (Figures 4-7).
8. stats_analysis.m               
--- Script for statistical analysis and plotting of SWR characteristics in different task phases (Figure 8).
9. x_corr_calculate.m             
--- Script for calculating cross correlation statistics (Section 3.2).

