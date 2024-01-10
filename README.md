&theme=cobalt
# JKCS (=JKCS2.1) 
Jammy Key for Configurational Sampling

  HOW TO SETUP:
  
  (INSTALLATION: https://youtu.be/9mw45tbj1G4)

    git clone https://github.com/kubeckaj/JKCS2.1.git #copy JKCS (JKCS folder is created)
    cd JKCS2.1
    sh setup.sh                                       #run basic "setup" (installation)
    source ~/.bashrc                                  #required after the first "setup"
    
 ADJUST USER SETUP:
 
 (H.Vehkamäki group: Puhti users do not have to adjust it)
 
    vim/nano/emacs ~/.JKCSusersetup.txt               #adjust user setup
    #     CHANGE:
    # ABCluster: rigidmoloptimizer path
    # XTB: xtb path
    # PYTHON: please setup how to use python2.0 (python 3.0 is just for GoodVibes in JKCS4)
    # Gaussian: G16 path and module
    # Orca: Orca path and modules
    
    sh test.sh                                        #test that all your paths#       are set correctly
   
  Simple TEST (on local computer):
  
  (VERY FAST SHOW OF JKCS: https://youtu.be/xKKWZrO-EmU)
  
  (USING JKCS: https://youtu.be/C4dAkhU7O8E)
   
    cd $WRKDIR                 #go to your working directory
    mkdir TESTING_JKCS         #create a new forder 
    cd TESTING_JKCS            #and enter it
    JKCS0_copy --help          #CHECKING POINT:This fail if instal. wasn't succesful or you didn't source ~/.bashrc  
    JKCS0_copy SA A            #create input.txt for sulphuric acid and ammonia cluster
    ls
    #Here you can adjust the input.txt file [vim input.txt]
    JKCS1_prepare --help       #explains what is doing this script
    JKCS1_prepare              #prepare folders and some files based on input.txt
    ls                         #list files and folders
    JKCS2_runABC --help        #explains what is doing this script
    JKCS2_runABC -gen 10 -lm 10 -pop 10  -loc
                               #CHECK POINT:This fail if you didn't setup rigidmoloptimizer
                               #run 10 bee colony for 10 generations and save 10 best struc.
                               #enter each folder and perform ABCluster exploration
    JKcheck                    #check that all calculations were finished
    #optional: you can perform [JKCS4_collect ABC -loc] to collect ABC results
    JKCS3_runXTB --help        #shows help
    JKCS3_runXTB -loc          #CHECKING POINT:This fail if you didn't setup xtb paths in programs.txt properly
    JKCS4_collect XTB -loc     #CHECKING POINT:This fail if you didn't setup python properly
                               #this script collect semi-empirically(xtb) optimized structures
    cd SYS_1SA_1A              #Enter folder with cluster 1sa1a
    cat resultsXTB.dat         #see results [structure | gyration rad. | energy]
    molden movieXTB.xyz        #visualize molecules
     
    #######################
    # If everything worked, continue, otherwise, contact Jakub Kubecka (ja-kub-ecka@chem.au.dk).
    # adjust the rest of file ~/.JKCSusersetup.txt:
    # QUEING, QUANTUM CHEMISTRY PROGRAMS etc.
    
  and the rest, you have to consult with manual 
    
    
   
 
