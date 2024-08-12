# MH_Paramater_Estimation-3MPEA_Crystallization

This code is used to generate the results presented in 

Moving Horizon Parameter Estimation for an Enzyme Catalyzed Transamination Reaction with Integrated Product Removal - 
M.Sc. Jessica Behrens, M.Sc. Sven Tiedemann, M.Sc. Tom Kunde, Prof. Dr. Jan von Langermann  and Prof. Dr.-Ing. Achim Kienle

For more information please refer to the paper [link will follow]().
The code is structured as follows: 

**General Functions:**
- folder with subfolders containing all necessary files and functions to make some general definitions, store parameters etc. for all further calculations
* folder needs to be incoperated into every 'main_[].m' by adding the file path - The location where to add the path is indicated in the 'main_[].m'-files

**MHE_batch:** 
- moving horizon parameter estimation for the batch scenario

**MHE_rep_batch:** 
- moving horizon parameter estimation for the repetitive batch scenario 

**Profile_Likelihood:** 
- profile likelihood calculations
  - The confidence intervals can be tested by your own generation of test samples or by loading '5000_para_est.mat'
    where 5000 test samples were generated and saved (see respective main file)

  - subfolders:
    - PL_in_MHE_batch - calculation of profile likelihoods in the moving horizon window for the batch scenario
    - PL_in_MHE_rep_batch - calculation of profile likelihoods in the moving horizon window for the repetitive batch scenario
     
    - _Note:_ both calculations need additional inputs generated by the MHE_batch/MHE_rep_batch calculation. So first run them and afterwards the 
	      PL_in_MHE_... calculations!

The main file, where the calculations are performed is always named 'main_[].m' - Run those files. All other m-files simply generate some nice figures. Contact me in case some things do not work for you or in case of additional questions. Have fun!
