Test Data, updated MHA 18.3.2019 - tested on 87fa4eecc852e3933047909178e8f90a36757839
Should run on all branches stable/beta/HLRS/skyregions

0a_basic_data.zip - contains all 'elementary' data files required to run a sample fixed-tilt solar project from scratch,
	including generation of module-model interpolant structure, and horizon-profiles, which take a few minutes.

0a_sample_project.ssp - self-contained project file, ready to run shading analysis (~5000 light hours in one full year)

Copy either file (unzip) to a new directory (do not mess up the GIT folder!), set as MATLAB working directory and run splitGUI() 
or splitscript(). Using SplitGUI, follow all simulation steps in order. Using splitscript() grab a coffee and cross your fingers.

case_generator.zip - contains the folder ./common, which along testbase.m can generate a full set (all tracker types) of tests. 

Unzip to a new directory, copy testbase.m along with the extracted ./common folder. Set as working directory and run 
testbase('solve',true) [see TESTBASE doc. for additional options]. 
File ./common/AddSimOptions.m contains by default the option meteofilter = '100rnd', which means a reduced, 100-random-sample of
hours will be used for each simulation. Option 'solve' = true handles all cases consecutively in splitscript() mode.
