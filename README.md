# iEEG_MEM
  The Distance_SC_Null.m generates structural connectivity nulls for a 
  subset of anatomical regions of interests (ROIs). In this demo script the 
  empirical subset of ROIs correspond to ROIs closest to a sample patinet's 
  intra-cranial EEG (iEEG) electrodes. This algorithm creates structural 
  nulls by iteratively resampling ROIs aimed at keeping the pairwise distance
  between the null ROI pairs closest to the empirical distance profiles. 
  "NULL_ROIS" and "Null_Dist" variables contain the null ROI labels 
  (i.e., numbers) and the distance between null ROI pairs, respectively.
 
  
  Arian Ashourvan, University of Pennsylvania, july 2020 
  

  Please contact Arian Ashourvan (ashourv@seas.upenn.edu) with any questions regarding this code.
  ========================================================================
