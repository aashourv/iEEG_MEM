# iEEG_MEM
  Code to reproduce geometric structrual connectivity null model in ashourvan et al. 2020 
  ("Pairwise maximum entropy model explains the role of white matter structure in shaping 
  emergent co-activation states and functional connectivity in intracranial EEG").
  
  The Distance_SC_Null.m generates structural connectivity nulls for a 
  subset of anatomical regions of interests (ROIs). In this demo script the 
  empirical subset of ROIs correspond to ROIs closest to a sample patinet's 
  intra-cranial EEG (iEEG) electrodes. This algorithm creates structural 
  nulls by iteratively resampling ROIs aimed at keeping the pairwise distance
  between the null ROI pairs closest to the empirical distance profiles. 
  "NULL_ROIS" and "Null_Dist" variables contain the null ROI labels 
  (i.e., numbers) and the distance between null ROI pairs, respectively.
  
  The average time to create a single null model in the sample script is around 100 seconds
  using Macbook pro OS X 10.13.6 2.8 GHz Intel Core i7. To increase the speed, we included 
  Distance_SC_Null_ShellScript.text script to allow the Distance_SC_Null.m script to be 
  submitted to a computing cluster using a Sun Grid Engine job scheduler (qsub)
  
  Arian Ashourvan, University of Pennsylvania, july 2020 
  

  Please contact Arian Ashourvan (ashourv@seas.upenn.edu) with any questions regarding this code.
  ========================================================================
