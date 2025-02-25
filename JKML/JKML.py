####################################################################################################
####################################################################################################
#                                                                                                  #
#    JJJJ K  K L    M   M                                                                          #
#       J K K  L    MM MM                                                                          #
#       J KK   L    M M M                                                                          #
#    J  J KKK  L    M   M                                                                          #
#     JJ  K  K LLLL M   M                                                                          #
#                                                                                                  #
####################################################################################################
####################################################################################################

print(
    """
   JJJJJJ KK   KK MMM    MMM LL
   JJ  JJ KK  KK  MMMM  MMMM LL
       JJ KK KK   MM MMMM MM LL
   JJ  JJ KKKK    MM  MM  MM LL
   JJ  JJ KKKKK   MM      MM LL
   JJJJJJ KK KKK  MM      MM LLLLLL
    JJJJ  KK  KKK MM      MM LLLLLLL
"""
)

print("JKML has started", flush=True)

# TREATING ARGUMENTS
from src.arguments import arguments
from sys import argv

locals().update(arguments(argv[1:]))

# LOADING DATABASES
from src.load_databases import load_databases

(
    train_high_database,
    train_low_database,
    test_high_database,
    test_low_database,
    monomers_high_database,
    monomers_low_database,
) = load_databases(
    Qtrain,
    Qeval,
    Qopt,
    Qmonomers,
    method,
    VARS_PKL,
    TRAIN_HIGH,
    TRAIN_LOW,
    TEST_HIGH,
    TEST_LOW,
    MONOMERS_HIGH,
    MONOMERS_LOW,
)

####################################################################################################
####################################################################################################

###################################
###### SAMPLEEACH/SIMILARITY ######
###################################

if Qtrain == 0:
    print("JKML: training option is missing [EXITING]", flush=True)
    exit()

# SAMPLE EACH SECTION:
if Qsampleeach > 0:
    from src.sampleeach import sampleeach_mbtr

    sampleeach_all, sample_train, sample_test = sampleeach_mbtr(
        train_high_database, test_high_database
    )
elif Qsampleeach < 0:
    from src.sampleeach import sampleeach_fchl

    sampleeach_all, sample_train, sample_test = sampleeach_fchl(
        train_high_database, test_high_database, krr_cutoff
    )
else:
    sampleeach_all = ["once"]


####################################################################################################
####################################################################################################

### Doing the whole process process for all or sampleeach/similarity
for sampleeach_i in sampleeach_all:
    ### SIMILARITY AND SAMPLE EACH
    if Qsampleeach > 0:
        from src.sampleeach import sampledist_mbtr

        sampledist = sampledist_mbtr(
            sample_train, sample_test, sampleeach_i, Qsampleeach
        )
    elif Qsampleeach < 0:
        from src.sampleeach import sampledist_fchl

        sampledist = sampledist_fchl(
            sample_train,
            sample_test,
            sampleeach_i,
            Qsampleeach,
            Qkernel,
            train_high_database,
        )
    else:
        sampledist = None

    ####################################################################################################
    ####################################################################################################

    #####################################
    ### TRAINING ########################
    #####################################

    if Qtrain == 1:

        ### DATABASE LOADING ###
        print("JKML is preparing things for training.")
        from src.data import prepare_data_for_training as prepare_data

        # returns: strs, Y_train, F_train, Qforces, Q_charge, Q_charges, Qcharge, D_dipole, Qdipole, size
        locals().update(
            prepare_data(
                train_high_database,
                monomers_high_database,
                train_low_database,
                monomers_low_database,
                seed,
                size,
                method,
                column_name_1,
                column_name_2,
                Qmin,
                Qifforces,
                Qifcharges,
                Qifdipole,
                Qsampleeach,
                Qforcemonomers,
                sampledist,
                Qmonomers,
            )
        )

        #####################################
        print("JKML is starting the training.")
        if Qmethod == "krr":
            from src.QML import training

            locals().update(
                training(
                    Qrepresentation,
                    Qkernel,
                    Qsplit,
                    strs,
                    Y_train,
                    krr_cutoff,
                    lambdas,
                    sigmas,
                    varsoutfile,
                    Qsplit_i,
                    Qsplit_j,
                )
            )
        elif Qmethod == "knn":
            from JKML.src.QKNN import training

            locals().update(
                training(
                    Qrepresentation,
                    strs,
                    Y_train,
                    krr_cutoff,
                    varsoutfile,
                )
            )
        #####################################
        elif Qmethod == "nn":
            from src.SchNetPack import training

            training(
                Qforces,
                Y_train,
                F_train,
                Qenergytradoff,
                strs,
                nn_tvv,
                nn_cutoff,
                nw,
                nn_rbf,
                Qrepresentation,
                nn_atom_basis,
                nn_interactions,
                Qbatch_size,
                Qlearningrate,
                parentdir,
                seed,
                varsoutfile,
                Qearlystop,
                nn_epochs,
                Qcheckpoint,
                Qtime,
                Qifcharges,
                Q_charges,
                Qifdipole,
                D_dipole,
            )
        ###################################
        elif Qmethod == "physnet":
            from src.PhysNetInterface import training

            training(
                Qforces,
                Y_train,
                F_train,
                D_dipole,
                Q_charge,
                Q_charges,
                strs,
                nn_atom_basis,
                nn_rbf,
                nn_interactions,
                nn_cutoff,
                nn_tvv,
                Qbatch_size,
                seed,
                nn_epochs,
                Qlearningrate,
                varsoutfile,
            )
        ###################################
        else:
            print(
                "JKML: Wrong method ("
                + Qmethod
                + ") or representation ("
                + Qrepresentation
                + ") chosen."
            )
            exit()
        ###################################

    ####################################################################################################
    ####################################################################################################

    ######################
    ### LOAD TRAINING ####
    ######################

    # TODO collect splitting: /home/kubeckaj/ML_SA_B/ML/TRAIN/TEST/SEPARATE/cho_solve.pkl
    if Qtrain == 2:
        import os

        if not os.path.exists(VARS_PKL):
            print("JKML: Error reading trained model. VARS_PKL = " + VARS_PKL)
            exit()
        if Qmethod == "krr":
            f = open(VARS_PKL, "rb")
            import pickle

            if Qrepresentation == "fchl":
                X_train, sigmas, alpha = pickle.load(f)
            elif Qrepresentation == "mbdf":
                X_train, X_atoms, sigmas, alpha = pickle.load(f)
            if len(alpha) != 1:
                alpha = [alpha]
            f.close()
            print("JKML: Trained model loaded.", flush=True)
        elif Qmethod == "knn":
            import pickle

            with open(VARS_PKL, "wb") as f:
                X_train, X_atoms, A, knn = pickle.load(f)
        elif Qmethod == "nn" or Qmethod == "physnet":
            varsoutfile = VARS_PKL
            print("JKML: Trained model found.")
        else:
            print("JKML: Wrong method or representation chosen.")
            exit()

    ####################################################################################################
    ####################################################################################################

    ######################
    ## EVALUATE/TESTING ##
    ######################

    if Qeval == 1 or Qeval == 2:

        ### DATABASE LOADING ###
        print("JKML is preparing things for testing/evaluation.")
        from src.data import prepare_data_for_testing as prepare_data

        # returns: Qforces, Q_charges, Qcharge
        locals().update(
            prepare_data(
                test_high_database,
                test_low_database,
                monomers_high_database,
                monomers_low_database,
                Qsampleeach,
                sampleeach_i,
                method,
                size,
                seed,
                column_name_1,
                column_name_2,
                Qeval,
                Qifforces,
                Qmonomers,
                Qmin,
                Qprintforces,
                Qifcharges,
            )
        )

        #####################################
        print("JKML is starting the testing/evaluation.")
        if Qmethod == "krr":
            from src.QML import evaluate

            Y_predicted = evaluate(
                Qrepresentation, krr_cutoff, X_train, sigmas, alpha, strs, Qkernel
            )
            Qforces = 0
            F_predicted = None
            Qa_predicted = None
        #####################################
        elif Qmethod == "knn":
            from src.QKNN import evaluate

            Y_predicted = evaluate(Qrepresentation, krr_cutoff, X_train, strs, knn)
        elif Qmethod == "nn":
            from src.SchNetPack import evaluate

            Y_predicted, F_predicted, Qa_predicted = evaluate(
                Qforces, varsoutfile, nn_cutoff, clusters_df, method, Qmin, Qifcharges
            )
        #####################################
        elif Qmethod == "physnet":
            from src.PhysNetInterface import evaluate

            Y_predicted, F_predicted = evaluate(varsoutfile, clusters_df, method, Qmin)
            Qa_predicted = None
        #####################################
        else:
            print("JKML: Wrong method or representation chosen.")
            exit()
        #####################################

        ### PRINTING THE RESULTS
        from src.print_output import print_results

        print_results(
            clusters_df,
            Y_predicted,
            Y_validation,
            F_predicted,
            F_test,
            ens,
            ens_correction,
            form_ens,
            ens2,
            ens2_correction,
            form_ens2,
            method,
            Qeval,
            Qforces,
            Qwolfram,
            Qprintforces,
            outfile,
            sampleeach_i,
            Qsampleeach,
            column_name_1,
            column_name_2,
            clustersout_df,
            sigmas,
            Qcharge,
            Q_charges,
            Qa_predicted,
        )

    ####################################################################################################
    ####################################################################################################

    ######################
    ###### OPTIMIZE ######
    ######################
    # TODO OPTIMIZE SECTION IS NOT READY YET

    if Qopt > 0:
        print("JKML: CAREFULLY: OPTIMIZATION IS IN PROGRESS.")
        ### DATABASE LOADING ###
        ## The high level of theory
        clusters_df = test_high_database
        strs = clusters_df["xyz"]["structure"]

        # TODO: This part should be unified!
        #####################################
        if Qmethod == "krr" and Qrepresentation == "fchl":
            from src.QML import optimize

            optimize(
                strs,
                Qrepresentation,
                krr_cutoff,
                varsoutfile,
                Qsplit,
                nn_cutoff,
                Qkernel,
                Qopt,
                opt_maxstep,
                opt_dump,
                opt_steps,
                md_temperature,
                md_thermostatfriction,
                md_dump,
                md_steps,
            )
        #####################################
        elif Qmethod == "nn":
            from src.SchNetPack import optimize

            optimize(
                test_high_database,
                varsoutfile,
                nn_cutoff,
                Qopt,
                opt_maxstep,
                opt_dump,
                opt_steps,
                md_temperature,
                Qmd_timestep,
                md_thermostatfriction,
                md_dump,
                md_steps,
            )
        #####################################
        elif Qmethod == "physnet":
            print("JKML: Sorry this part was not implemented yet.")
            exit()
        #####################################
        else:
            raise Exception("JKML: Wrong method or representation chosen.")
        #####################################

        # TODO unify printing the results
        ### PRINTING THE RESULTS
        # print("Ypredicted = {" + ",".join([str(i) for i in save_energy])+"};", flush = True)

        ### PRINTING THE QML PICKLES
        # clustersout_df = clusters_df.copy()
        # for i in range(len(clustersout_df)):
        #  clustersout_df.loc[clustersout_df.iloc[i].name,(column_name_1,column_name_2)] = Y_predicted[0][i]
        # clustersout_df.to_pickle(outfile)
    ########

####################################################################################################
####################################################################################################
print("\nJKML has finished", flush=True)
