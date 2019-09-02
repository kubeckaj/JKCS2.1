# JKCS (=JKCS2.0) version 2.00
Jammy Key for Configurational Sampling

  HOW TO SETUP:

    cd ${YOUR_APPS_DIR}                               #enter to the directory with your scripts
    git clone https://github.com/kubeckaj/JKCS2.0.git #copy JKCS (JKCS folder is created)
    cd JKCS2.0
    ls
    sh setup.sh                                       #run basic "setup" (installation)
    source ~/.bashrc                                  #required after the first "setup"
   
    # adjust ~/.JKCSusersetup.txt:
    # ABCluster: rigidmoloptimizer path
    # XTB: xtb path, XTBHOME path
    # PYTHON: please setup how to use python
    # If neaded load modules: Mpython, Mgaussian ...
   
  Simple TEST (on local computer):
   
    cd $WRKDIR                 #go to your working directory
    mkdir TESTING_JKCS         #create a new forder 
    cd TESTING_JKCS            #and enter it
    JKCS0_copy --help          #CHECKING POINT:This fail if instal. wasn't succesful or you didn't source ~/.bashrc  
    JKCS0_copy SA A            #create input.txt for sulphuric acid and ammonia cluster
    ls
    #Here you can adjust the input.txt file [vim input.txt]
    JKCS1_prepare --help       #explains what is doing this script
    JKCS1_prepare              #prepare folders and some files based on input.txt
    JKCS2_runABC --help        #explains what is doing this script
    JKCS2_runABC -gen 10 -lm 10 -pop 10  #CHECK POINT:This fail if you didn't setup rigidmoloptimizer
                                         #run 10 bee colony for 10 generations and save 10 best struc.
                                         #enter each folder and perform ABCluster exploration
    #optional: you can perform [JKCS4_collect ABC] also here
    JKCS3_runXTB --help        
    JKCS3_runXTB               #CHECKING POINT:This fail if you didn't setup xtb paths in programs.txt properly
    JKCS4_collect XTB          #CHECKING POINT:This fail if you didn't setup python properly
                               #this script collect semi-empirically(xtb) optimized structures
    cd SYS_1SA_1A              #Enter folder with cluster 1sa1a
    molden movieXTB.xyz or vim resultsXTB.dat [structure_file gyration_radius energy] to see results.
     
    #######################
    # If everything worked, continue, otherwise, contact Jakub Kubecka (jakub.kubecka@helsinki.fi).
    # adjust the rest of file ~/.JKCSusersetup.txt:
    # QUEING, QUANTUM CHEMISTRY PROGRAMS etc.
    
  and the rest you have to consult with manual 
    
    
   
 
