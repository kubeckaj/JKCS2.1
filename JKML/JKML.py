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

# timing
from timeit import default_timer as timer

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
    Qhyper,
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

if Qtrain == 0 and not Qhyper:
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
        train_start = timer()
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
        #####################################
        elif Qmethod == "knn":
            from src.QKNN import training

            locals().update(
                training(
                    Qrepresentation,
                    strs,
                    Y_train,
                    varsoutfile,
                    no_metric=no_metric,
                    hyper_cache=hyper_cache,
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
        elif Qmethod == "aimnet":
            from src.AIMNetInterface import training
            
            training(
                Y_train,
                F_train,
                Z_atoms,
                N_atoms,
                Qenergytradoff,
                strs,nn_tvv,
                Qbatch_size,
                Qlearningrate,
                parentdir,
                varsoutfile,
                Qearlystop,
                nn_epochs,
                Q_charge,
                Q_charges,
                file_basenames,
                seed,
                Qmonomers
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

        train_time = timer() - train_start
        print(
            f"JKML training done. Took {train_time:.3f} s ({train_time / X_train.shape[0]:.4f} per sample)."
        )
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
                try:
                  X_train, sigmas, alpha, train_metadata = pickle.load(f)
                except:
                  f.close()
                  f = open(VARS_PKL, "rb")
                  train_metadata = {}
                  X_atoms = None
                  X_train, sigmas, alpha = pickle.load(f)
            elif Qrepresentation == "fchl19":
                X_train, X_atoms, sigmas, alpha, train_metadata = pickle.load(f)
            elif Qrepresentation == "mbdf":
                X_train, X_atoms, sigmas, alpha, train_metadata = pickle.load(f)
            if len(alpha) != 1:
                alpha = [alpha]
            f.close()
            print("JKML: Trained model loaded.", flush=True)
            # store the training metadata to locals
            locals().update(train_metadata)
        elif Qmethod == "knn":
            import pickle
            from sklearn.neighbors import KNeighborsRegressor

            with open(VARS_PKL, "rb") as f:
                if no_metric:
                    X_train, Y_train, X_atoms, knn_params, train_metadata = pickle.load(
                        f
                    )
                elif Qrepresentation == "fchl-kernel":
                    (
                        X_train,
                        Y_train,
                        X_atoms,
                        knn_params,
                        vp_params,
                        train_metadata,
                    ) = pickle.load(f)
                else:
                    (
                        X_train,
                        Y_train,
                        X_atoms,
                        A,
                        mlkr,
                        knn_params,
                        train_metadata,
                    ) = pickle.load(f)

            # need to recreate the model due to not being able to pickle the custom metric
            if not no_metric and Qrepresentation != "fchl-kernel":
                knn_params["metric"] = mlkr.get_metric()
            if Qrepresentation == "fchl-kernel":
                from src.QKNN import load_vp_knn

                knn = load_vp_knn(X_train, Y_train, vp_params, **knn_params)
            else:
                knn = KNeighborsRegressor(**knn_params)
                knn.fit(X_train, Y_train)
            # store the training metadata to locals
            locals().update(train_metadata)
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

        eval_start = timer()
        #####################################
        print("JKML is starting the testing/evaluation.")
        if Qmethod == "krr":
            from src.QML import evaluate

            Y_predicted, repr_test_wall, repr_test_cpu, test_wall, test_cpu, d_test = (
                evaluate(
                    Qrepresentation, krr_cutoff, X_train, X_atoms, sigmas, alpha, strs, Qkernel
                )
            )
            Qforces = 0
            F_predicted = None
            Qa_predicted = None
        #####################################
        elif Qmethod == "knn":
            from src.QKNN import evaluate

            Y_predicted, repr_test_wall, repr_test_cpu, test_wall, test_cpu, d_test = (
                evaluate(
                    Qrepresentation,
                    X_train,
                    strs,
                    knn,
                    hyper_cache=hyper_cache,
                )
            )
            Qforces = 0
            F_predicted = None
            Qa_predicted = None
        #####################################
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
        elif Qmethod == "aimnet":
            from src.AIMNetInterface import evaluate
            
            Y_predicted, F_predicted = evaluate(varsoutfile, clusters_df, method, Qmin, parentdir, Qmonomers, Z_atoms)
            Qa_predicted = None
        #####################################
        else:
            print("JKML: Wrong method or representation chosen.")
            exit()
        #####################################

        eval_time = timer() - eval_start
        time_per_sample = eval_time / X_train.shape[0]
        print(
            f"JKML evaluation done. Took {eval_time:.3f} s ({time_per_sample:.4f} per sample)."
        )
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
            repr_train_wall,
            repr_train_cpu,
            repr_test_wall,
            repr_test_cpu,
            train_wall,
            train_cpu,
            test_wall,
            test_cpu,
            n_train,
            d_train,
            d_test,
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

    #########################
    ###### HYPERPARAMS ######
    #########################

    if Qhyper:
        ### DATABASE LOADING ###
        print("JKML is preparing things for training.", flush=True)
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
        print("JKML is starting hyperparameter optimisation.", flush=True)
        if Qmethod == "knn":
            from src.QKNN import hyperopt

            params = hyperopt(
                Qrepresentation,
                strs,
                Y_train,
                varsoutfile,
                no_metric,
                hyper_cache=hyper_cache,
            )
        else:
            raise ValueError(
                f"Hyperparameter optimisation not implemented for {Qmethod}!"
            )
        print(f"JKML: optimal parameters:")
        if "representation" in params:
            print("-------Representation-------")
            for k, v in params["representation"].items():
                print(f"{k}: {v}")
        if "knn" in params:
            print("-------KNN-------")
            for k, v in params["knn"].items():
                print(f"{k}: {v}")
        print("", flush=True)

####################################################################################################
####################################################################################################
print("\nJKML has finished", flush=True)
