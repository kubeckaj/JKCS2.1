def flatten(matrix):
    from numpy import array

    return array([item for row in matrix for item in row])


def print_results(
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
):

    from numpy import array, isnan, mean, std, abs, sqrt, median, ndarray

    if Qwolfram == 0:
        lb = "["
        rb = "]"
    else:
        lb = "{"
        rb = "}"

    print("\n########   RESULTS   #########\n")
    #print("Y_predicted = ", Y_predicted[0], flush=True)
    #print("ens = ", ens, flush=True)
    #print("ens_correction = ", ens_correction, flush=True)
    ### PRINT THE PROPERTIES ###
    if Qeval == 2:
        print("Ytest = " + lb + ",".join([str(i) for i in ens]) + rb + ";", flush=True)
    if method != "delta":
        print(
            "Ypredicted = "
            + lb
            + ",".join([str(i) for i in Y_predicted[0] + ens_correction])
            + rb
            + ";",
            flush=True,
        )
    else:
        print(
            "Ypredicted = "
            + lb
            + ",".join([str(i) for i in Y_predicted[0] + form_ens2 + ens_correction])
            + rb
            + ";",
            flush=True,
        )

    ### PRINT THE FORCES
    if Qforces == 1:
        if Qprintforces == 0:
            if Qeval == 2:
                print(
                    "Ftest = I will not print the forces becuase it would be too many numbers. (Use -printforces)",
                    flush=True,
                )
            print(
                "Fpredicted = I will not print the forces becuase it would be too many numbers. (Use -printforces)",
                flush=True,
            )
        else:
            if Qeval == 2:
                print(
                    "Fpredicted = "
                    + lb
                    + ",".join(
                        [
                            lb
                            + ",".join(
                                [lb + ",".join([str(k) for k in j]) + rb for j in i]
                            )
                            + rb
                            for i in F_predicted[0]
                        ]
                    )
                    + rb
                    + ";",
                    flush=True,
                )
            print(
                "Ftest = "
                + lb
                + ",".join(
                    [
                        lb
                        + ",".join([lb + ",".join([str(k) for k in j]) + rb for j in i])
                        + rb
                        for i in F_test
                    ]
                )
                + rb
                + ";",
                flush=True,
            )

    ### PRINT THE STATISTICS
    # If possible print MAE and RMSE of energies and forces
    if Qeval == 2:
        ### Kick-out nans
        err = len(Y_predicted[0][isnan(Y_predicted[0])])
        if err > 0:
            print(
                "\nERROR: "
                + str(err)
                + " NaNs in the predicted set.\nThose have been removed from statistics.\n",
                flush=True,
            )
            valid_indices = ~isnan(Y_predicted[0])
            Y_validation_rem = Y_validation
            Y_validation = Y_validation[valid_indices]
            Y_predicted_rem = Y_predicted
            Y_predicted = [Y_predicted[0][valid_indices]]
            if Qforces == 1:
                F_predicted_rem = F_predicted
                F_predicted = [list(array(F_predicted[0])[valid_indices])]
                F_test = list(array(F_test)[valid_indices])
        ### Calculate mean-absolute-error (MAE):
        if column_name_1 == "log" and column_name_2 == "electronic_energy":
            multiplier = 627.503
            units = " kcal/mol"
        else:
            multiplier = 1.0
            units = " [?]"
        print("", flush=True)
        print("Results:", flush=True)
        mae = [
            multiplier * mean(abs(Y_predicted[i] - Y_validation))
            for i in range(len(sigmas))
        ]
        sd = [
            multiplier * std(abs(Y_predicted[i] - Y_validation))
            for i in range(len(sigmas))
        ]
        print("Mean Absolute Error:", flush=True)
        print(
            "mae = "
            + ",".join([str(i) for i in mae])
            + units
            + "(+- "
            + str(sd[0])
            + " )",
            flush=True,
        )
        ### Calculate root-mean-squared-error (RMSE):
        rmse = [
            multiplier * sqrt(mean(abs(Y_predicted[i] - Y_validation) ** 2))
            for i in range(len(sigmas))
        ]
        print("Root Mean Squared Error", flush=True)
        print("rmse = " + ",".join([str(i) for i in rmse]) + units, flush=True)
        ### Calculate mean-absolute-relative-error (MArE):
        diff = [mean(Y_predicted[i]) - mean(Y_validation) for i in range(len(sigmas))]
        mare = [
            multiplier * mean(abs(Y_predicted[i] - Y_validation - diff[i]))
            for i in range(len(sigmas))
        ]
        print("Mean Absolute (mean-)Relative Error:", flush=True)
        print("mare = " + ",".join([str(i) for i in mare]) + units, flush=True)
        ### Calculate root-mean-squared-relative-error (RMSrE):
        rmsre = [
            multiplier * sqrt(mean(abs(Y_predicted[i] - Y_validation - diff[i]) ** 2))
            for i in range(len(sigmas))
        ]
        print("Root Mean Squared (mean-)Relative Error", flush=True)
        print("rmsre = " + ",".join([str(i) for i in rmsre]) + units, flush=True)
        diff = [
            median(Y_predicted[i]) - median(Y_validation) for i in range(len(sigmas))
        ]
        mare = [
            multiplier * mean(abs(Y_predicted[i] - Y_validation - diff[i]))
            for i in range(len(sigmas))
        ]
        print("Mean Absolute (median-)Relative Error:", flush=True)
        print("mare = " + ",".join([str(i) for i in mare]) + units, flush=True)
        ### Calculate root-mean-squared-relative-error (RMSrE):
        rmsre = [
            multiplier * sqrt(mean(abs(Y_predicted[i] - Y_validation - diff[i]) ** 2))
            for i in range(len(sigmas))
        ]
        print("Root Mean Squared (median-)Relative Error", flush=True)
        print("rmsre = " + ",".join([str(i) for i in rmsre]) + units, flush=True)

        if Qforces == 1:
            multiplier = 627.503
            units = " [kcal/mol/Angstrom]"
            print("", flush=True)
            print("Results for forces:", flush=True)
            mae = [
                multiplier * mean(abs(flatten(F_predicted[i]) - flatten(F_test)))
                for i in range(len(sigmas))
            ]
            std = [
                multiplier * std(abs(flatten(F_predicted[i]) - flatten(F_test)))
                for i in range(len(sigmas))
            ]
            print(
                "MAE = "
                + ",".join([str(i) for i in mae])
                + units
                + "(+- "
                + str(std[0])
                + " )",
                flush=True,
            )
            rmse = [
                multiplier
                * sqrt(mean(abs(flatten(F_predicted[i]) - flatten(F_test)) ** 2))
                for i in range(len(sigmas))
            ]
            print("RMSE = " + ",".join([str(i) for i in rmse]) + units, flush=True)

        if Qcharge == 1 and Qa_predicted is not None:
            print("", flush=True)
            print("Results for charges:", flush=True)
            Qa_predicted = flatten(array([flatten(array(i)) for i in Qa_predicted]))
            Q_charges = flatten(array([array(i) for i in Q_charges]))
            mae = [mean(abs(Qa_predicted - Q_charges))]
            print("MAE = " + ",".join([str(i) for i in mae]) + " e", flush=True)
            rmse = [sqrt(mean(abs(Qa_predicted - Q_charges) ** 2))]
            print("RMSE = " + ",".join([str(i) for i in rmse]) + " e", flush=True)

        if err > 0:
            Y_predicted = Y_predicted_rem
            Y_validation = Y_validation_rem
            if Qforces == 1:
                F_predicted = F_predicted_rem

    ### PRINTING THE QML PICKLES
    if Qforces == 1:
      if ("extra", "forces") not in clustersout_df.columns:
          # clustersout_df.loc[clustersout_df.iloc[i].name,("extra","forces")] = [array(force) for force in F_predicted[0][i]]
          clustersout_df.loc[:, ("extra", "forces")] = None
      if ("extra", "forces_true") not in clustersout_df.columns:
          clustersout_df.loc[:, ("extra", "forces_true")] = None
    for i in range(len(clustersout_df)):
        if method != "delta":
            clustersout_df.loc[
                clustersout_df.iloc[i].name, (column_name_1, column_name_2)
            ] = (Y_predicted[0][i] + ens_correction[i])
        else:
            clustersout_df.loc[
                clustersout_df.iloc[i].name, (column_name_1, column_name_2)
            ] = (Y_predicted[0][i] + form_ens2[i] + ens_correction[i])
        if Qeval == 2:
            clustersout_df.loc[clustersout_df.iloc[i].name, ("extra", "error")] = (
                multiplier * abs(Y_predicted[0][i] - Y_validation[i])
            )
            clustersout_df.loc[clustersout_df.iloc[i].name, ("extra", "Y_true")] = (
                multiplier * Y_validation[i]
            )
        if Qforces == 1:
            clustersout_df.at[clustersout_df.iloc[i].name, ("extra", "forces_true")] = [
                array(force) for force in F_test[i]
            ]
            clustersout_df.at[clustersout_df.iloc[i].name, ("extra", "forces")] = [
                array(force) for force in F_predicted[0][i]
            ]
            
    ### SAVE TIME AND SHAPE INFORMATION
    clustersout_df.loc[:, ("extra", "repr_train_wall")] = repr_train_wall
    clustersout_df.loc[:, ("extra", "repr_train_cpu")] = repr_train_cpu
    clustersout_df.loc[:, ("extra", "repr_test_wall")] = repr_test_wall
    clustersout_df.loc[:, ("extra", "repr_test_cpu")] = repr_test_cpu
    clustersout_df.loc[:, ("extra", "train_wall")] = train_wall
    clustersout_df.loc[:, ("extra", "train_cpu")] = train_cpu
    clustersout_df.loc[:, ("extra", "test_wall")] = test_wall
    clustersout_df.loc[:, ("extra", "test_cpu")] = test_cpu
    clustersout_df.loc[:, ("extra", "n_train")] = n_train
    clustersout_df.loc[:, ("extra", "d_train")] = d_train
    if isinstance(Y_predicted, list):
        clustersout_df.loc[:, ("extra", "n_test")] = array(Y_predicted).shape[1]
    elif isinstance(Y_predicted, ndarray):
        clustersout_df.loc[:, ("extra", "n_test")] = Y_predicted.shape[1]
    clustersout_df.loc[:, ("extra", "d_test")] = d_test

    clustersout_df.to_pickle(outfile)
    print(f"Saved results to {outfile}", flush=True)
    if Qsampleeach > 0:
        import os
        if sampleeach_i == 0:
            os.system("JKQC " + outfile + " -out predicted_QML_FULL.pkl -noex")
        else:
            os.system(
                "JKQC "
                + outfile
                + " predicted_QML_FULL.pkl -out predicted_QML_FULL.pkl -noex"
            )

    return
