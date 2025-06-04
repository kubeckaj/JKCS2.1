# PREDEFINED ARGUMENTS
import os.path


def print_help():
    print(
        "#################################################################", flush=True
    )
    word = "JKML"
    colors = ["\033[34m", "\033[31m", "\033[32m", "\033[33m"]

    # 31red, 32green, 33yellow, 34blue, 35magenta, 36cyan
    def pJKML():
        for i, letter in enumerate(word):
            color = colors[i % len(colors)]  # cycle through colors
            print(f"{color}{letter}\033[0m", end="", flush=True)

    pJKML()
    print(" HELP:", flush=True)
    print("This script interfaces ML model and pickled structures.", flush=True)
    print("", flush=True)
    print(f"  ", end="")
    pJKML()
    print(" [OPTION(s)]", flush=True)
    print("", flush=True)
    print("  HELP OPTIONS:", flush=True)
    print("    -help                            prints basic help", flush=True)
    print(
        "    -help_nn                         prints help for neural network methods (e.g. PaiNN,SchNet,PhysNet)",
        flush=True,
    )
    print(
        "    -help_krr                        prints help for kernel ridge regression methods (e.g. FCHL)",
        flush=True,
    )
    print(
        "    -help_adv                        print some advanced features (e.g., OPT, MD, seed, splits)"
    )
    print("", flush=True)
    print("  OPTIONS:", flush=True)
    print(
        "    -qml                             use KRR with FCHL19 [set by default]",
        flush=True,
    )
    print("    -mbdf                            use KRR with MBDF", flush=True)
    print(
        "    -nn,-painn                       switch to NN = neural network with PaiNN",
        flush=True,
    )
    print(
        "    -schnet                          switch to NN = neural network with SchNet",
        flush=True,
    )
    print(
        "    -physnet                         switch to NN = neural network with PhysNet",
        flush=True,
    )
    print(
        "    -aimnet                          switch to NN = neural network with AIMNet2",
        flush=True,
    )
    print("", flush=True)
    print("  INPUT FILES:", flush=True)
    print(
        "    -train <HIGH.pkl> [<LOW.pkl>]    train on given pikled files", flush=True
    )
    print("    -trained <model.pkl>             take pre-trained ML model", flush=True)
    print(
        "    -monomers <HIGH.pkl> [<LOW.pkl>] properties with respect to monomer(s) in pickled file(s)",
        flush=True,
    )
    print(
        "    -test <HIGH.pkl> [<LOW.pkl>]     predict energies of structures in pickled file(s) with known truth",
        flush=True,
    )
    print(
        "    -eval <STRS.pkl> [<LOW.pkl>]     predict energies of structures in pickled file(s) with uknown truth (no verification)",
        flush=True,
    )
    print("", flush=True)
    print("  OUTPUT FILES:", flush=True)
    print(
        "    -out <file.pkl>                  name of file with predicted values [def = predicted.pkl]",
        flush=True,
    )
    print(
        "    -modelout <file.pkl>             name of file with saved model [def = model.pkl]",
        flush=True,
    )
    print("", flush=True)
    print("  EXAMPLES:", flush=True)
    print(f"    ", end="")
    pJKML()
    print(
        " -loc -train collected_structures.pkl -modelout trained_example.pkl",
        flush=True,
    )
    print(f"    ", end="")
    pJKML()
    print(
        " -par q64 -cpu 64 -train tr_high.pkl tr_low.pkl -test te_high.pkl te_low.pkl",
        flush=True,
    )
    print(f"    ", end="")
    pJKML()
    print(
        " -par qgpu -cpu 2 -nn -epochs 10000 -train x.pkl -eval str.pkl -monomers monomers.pkl",
        flush=True,
    )
    print(f"    ", end="")
    pJKML()
    print(
        " -loc -trained model_for_atomization_en.pkl -eval diff_molecules.pkl -monomers atoms.pkl ",
        flush=True,
    )
    print("    ", end="")
    pJKML()
    print(
        " -loc -train train.pkl -physnet -nn_features 128 -nn_basis 64 -nn_blocks 5 -epochs 1000",
        flush=True,
    )
    print("", flush=True)


def help_adv():
    print("  ADVANCED OPTIONS:", flush=True)
    print(
        "    -column <str> <str> selects a different column from the pickled files (e.g. log entropy)",
        flush=True,
    )
    print(
        "    -noforces           forces are not trained/tested even if it is possible"
    )
    print(
        "    -size <int>         randomly selects only portion of the trainin database (e.g. 200)",
        flush=True,
    )
    print("    -seed <int>         seed for random number generators [def = 42]")
    print(
        "    -categorize <int>   selects structures with similar MBTR (our similarity definition)",
        flush=True,
    )
    print(
        "    -similarity <int>   selects structures with similar FCHL (uses kernel)",
        flush=True,
    )
    print(
        "    -forcemonomers      adds (extra) monomers to sampleeach/selection",
        flush=True,
    )
    print(
        "    -printforces        print out all forces (this might be a lot of numbers)",
        flush=True,
    )
    print("   OPT:", flush=True)
    print("    -opt <STRS.pkl>     optimize structure based on model [NN]", flush=True)
    print("    -opt_maxs <float>   max step in Angstrom in optimization [def = 0.02]")
    print(
        "    -opt_steps <int>    maximal number of the optimization steps [def = 100]",
        flush=True,
    )
    print(
        "    -opt_dump <int>                 dump structure every n-th step [def = 1]",
        flush=True,
    )
    print("   MD:", flush=True)
    print(
        "    -md <STRS.pkl>                 run md starting from provided structure(s) based on model [NN]",
        flush=True,
    )
    print(
        "    -md_temperature <float>        temperature of the simulation [def = 300.0]",
        flush=True,
    )
    print(
        "    -md_timestep <int>             integration time step [def = 0.1]",
        flush=True,
    )
    print(
        "    -md_thermostatfriction <float> strength of the themrostat [def = 0.002]",
        flush=True,
    )
    print("    -md_steps <int>                number of steps [def = 2000]", flush=True)
    print(
        "    -md_dump <int>                 dump structure every n-th step [def = 1]",
        flush=True,
    )
    print(
        "    -nn_cutoff <float>             cutoff function (Angstrom) [def = 5.0]",
        flush=True,
    )
    print("   SPKMD:", flush=True)
    print("    -spkmd                         use spkmd script", flush=True)
    print("    -langevin                      use langevin thermostat", flush=True)
    print(
        "    -rpmd <int>                    run ring polymer MD with n simulations",
        flush=True,
    )
    print(
        "    -md <STRS.xyz>                 run md starting from provided structure(s) based on model [NN]",
        flush=True,
    )
    print(
        "    -md_timestep <int>             integration time step [def = 0.1]",
        flush=True,
    )
    print(
        "    -md_thermostatconstant <float> frequency of the themrostat [def = 100]",
        flush=True,
    )
    print("    -md_steps <int>                number of steps [def = 2000]", flush=True)
    print(
        "    -spkmd_extra <string>          extra arguments based on spkmd", flush=True
    )
    print(
        "    -nn_cutoff <float>             cutoff function (Angstrom) [def = 5.0]",
        flush=True,
    )
    print("", flush=True)
    print("  EXTRA ADVANCED OPTIONS:", flush=True)
    print(
        "    -so3net             switch to NN = neural network with SO3net (from SchNetPack)",
        flush=True,
    )
    print(
        "    -split <int>        only with -krr/-fchl how many splits of KRR matrix do you do",
        flush=True,
    )
    print(
        "    -startsplit <int>   the same like above but only construct kernels",
        flush=True,
    )
    print(
        "    -finishsplit <int>  (see -split) combines the splitted kernel and creates model",
        flush=True,
    )
    print("    -wolfram            prints {} instead of []", flush=True)
    print("", flush=True)


def help_krr():
    print("  OPTIONS FOR KERNEL RIDGE REGRESSION:", flush=True)
    print(
        "    -sigma <int>        Gaussian width hyperparameter [def = 1.0]", flush=True
    )
    print(
        "    -lambda <int>       numerical stability (for matrix inv.) hyperparameter [def = 1e-4]"
    )
    print("    -laplacian          switch to Laplacian kernel (FCHL)")
    print("    -krr_cutoff <float> cutoff function (Angstrom) [def = 10.0]", flush=True)
    print("", flush=True)


def help_nn():
    print("  OPTIONS FOR NEURAL NETWORKS:", flush=True)
    print("    -epochs <int>              number of epochs [def = 1000]", flush=True)
    print(
        "    -batch_size,-bs <int>      batch size [def = 16], the same size is used for validation",
        flush=True,
    )
    print(
        "    -nn_train <float>          portion of training data (exlc. validation) [def = 0.9]",
        flush=True,
    )
    print(
        "    -nn_ESpatience <int>       Early-Stop patience of epochs for no improvement [def = 200]",
        flush=True,
    )
    print(
        "    -nn_energytradeoff <float> trade-off [energy, force] = [<float>, 1] [def = 0.01]",
        flush=True,
    )
    print("    -nn_lr                     learning rate [def = 1e-4]", flush=True)
    print(
        "    -nw                        number of workers for database manipulation [def = 1] {SchNetPack}",
        flush=True,
    )
    print(
        "    -ckpt,-chkp <file>         resume from last check point [def = None] {SchNetPack}",
        flush=True,
    )
    print("", flush=True)
    print("  OPTIONS FOR REPRESENTATION:", flush=True)
    print(
        "    -nn_ab, -nn_features <int>    number of atom basis/features/size of embeding vector [def = 256]",
        flush=True,
    )
    print(
        "    -nn_int, -nn_blocks <int>     number of (interaction) blocks [def = 5]",
        flush=True,
    )
    print(
        "    -nn_rb, -nn_basis <int>       number of (radial) basis [def = 20]",
        flush=True,
    )
    print(
        "    -nn_cutoff <float>            cutoff function (Angstrom) [def = 5.0]",
        flush=True,
    )
    print("", flush=True)


def arguments(argument_list=[]):
    # Predefined arguments
    method = "direct"  # direct/delta/min
    Qmin = 0
    size = "full"
    seed = 42
    TRAIN_HIGH = ""
    TRAIN_LOW = ""
    TEST_DATABASE = ""
    TEST_LOW = ""
    TEST_HIGH = ""
    MONOMERS_LOW = ""
    MONOMERS_HIGH = ""
    VARS_PKL = ""
    Qtrain = 0  # 0=nothing (fails), 1=train, 2=trained
    Qsplit = 1  # 1=no split, how many splits to do; ONLY FOR TRAINING
    Qsplit_i = 1
    Qsplit_j = 1
    Qeval = 0  # 0=nothing (possible), 1=validate, 2=eval
    Qopt = 0  # 0=nothing (possible), 1=optimize
    Qwolfram = 0  # 1 means that {} will be printed insead of []
    Qprintforces = 0  # print forces?
    Qmonomers = 2  # 0=monomers taken from database, 1=monomers in separate files, 2=no monomer subtraction
    Qsampleeach = 0
    Qforcemonomers = 0
    column_name_1 = "log"
    column_name_2 = "electronic_energy"
    Qifforces = 1  # IF forces exist use them in calculations
    Qifcharges = 0 
    Qifdipole = 0
    Qifeldisp = 0
    Qifforcedispi = 0

    Qmethod = "krr"
    Qrepresentation = "fchl"

    # knn default
    no_metric = False

    # Predefined QML
    Qkernel = "Gaussian"
    sigmas = [1.0]
    lambdas = [1e-4] * len(sigmas)
    krr_cutoff = 10.0

    # Predefined NN - PaiNN
    nn_rbf = 20
    nn_tvv = 0.9
    nn_cutoff = 5.0
    nn_atom_basis = 256
    nn_interactions = 5
    nn_epochs = 1000
    Qlearningrate = 1e-4
    Qearlystop = 200
    Qenergytradoff = 0.01  # if forces are trained on: [energy, force] = [X, 1]
    nw = 1
    Qbatch_size = 16
    Qcheckpoint = None
    Qtime = None
    parentdir = "./"

    # OPT/MD
    opt_maxstep = 0.02
    opt_steps = 100
    opt_dump = 1
    md_temperature = 300.0
    md_timestep = 0.1
    md_thermostatfriction = 0.002
    md_steps = 2000
    md_dump = 1

    # hyperparams
    Qhyper = False
    hyper_cache = None

    # Output files
    outfile = "predicted.pkl"
    varsoutfile = "model.pkl"

    ################################################################################################

    last = ""
    for arg in argument_list:
        # Help
        help_list = ["-help", "-help_nn", "-help_adv", "-help_krr"]
        if arg in help_list:
            print_help()
            if arg == "-help_nn":
                help_nn()
            if arg == "-help_adv":
                help_adv()
            if arg == "-help_krr":
                help_krr()
            print(
                "#################################################################",
                flush=True,
            )
            exit()

        # Vars out file
        var_out_file_list = ["-varsout", "-modelout"]
        if arg in var_out_file_list:
            last = "-varsout"
            continue

        if last == "-varsout":
            varsoutfile = arg
            varsoutfile = varsoutfile.split(".pkl")[0] + ".pkl"
            last = ""
            continue

        # Out file
        if arg == "-out":
            last = "-out"
            continue

        if last == "-out":
            outfile = arg
            outfile = outfile.split(".pkl")[0] + ".pkl"
            last = ""
            continue

        # Method
        if arg == "-method":
            last = "-method"
            continue

        if arg == "direct" or last == "-method" or arg == "delta":
            last = ""
            if arg == "direct":
                method = "direct"
            elif arg == "delta":
                method = "delta"
            else:
                raise Exception(f"I do not understand your input:\n{arg}")
            continue

        # Force monomers
        if arg == "-forcemonomers":
            Qforcemonomers = 1
            continue

        # Time
        if arg == "-time":
            last = "-time"
            continue
        if last == "-time":
            last = ""
            from datetime import timedelta

            def parse_duration(duration):
                try:
                    if "-" in duration:
                        days, time_str = duration.split("-")
                    else:
                        days = 0
                        time_str = duration
                    time_str_split = time_str.split(":")
                    days = int(days)
                    if len(time_str_split) == 4:
                        days, hours, minutes, seconds = map(int, time_str_split)
                    elif len(time_str_split) == 3:
                        hours, minutes, seconds = map(int, time_str_split)
                    elif len(time_str_split) == 1:
                        raise Exception("It is too short time")
                    else:
                        hours = 0
                        minutes, seconds = map(int, time_str_split)
                except Exception:
                    raise Exception(f"I do not understand your input:\n{arg}")

                return timedelta(
                    days=days, hours=hours, minutes=minutes, seconds=seconds
                ) - timedelta(days=0, hours=0, minutes=20, seconds=0)

            if arg is not None:
                Qtime = parse_duration(arg)
            continue

        # Seed
        if arg == "-seed":
            last = "-seed"
            continue
        if last == "-seed":
            last = ""
            try:
                seed = int(arg)
            except ValueError:
                raise Exception(f"Seed must be an integer. Got {arg}. [EXITING]")
            continue

        # Parent dir
        if arg == "-dir":
            last = "-dir"
            continue
        if last == "-dir":
            last = ""
            try:
                parentdir = str(arg)
            except Exception:
                raise Exception(f"I do not understand your input:\n{arg}")
            continue

        # Wolfram
        if arg == "-wolfram":
            Qwolfram = 1
            continue

        # Print forces
        if arg == "-printforces":
            Qprintforces = 1
            continue

        # Column
        if arg == "-column":
            last = "-column"
            continue
        if last == "-column":
            column_name_1 = arg
            last = "-column2"
            continue
        if last == "-column2":
            column_name_2 = arg
            last = ""
            continue

        # Hyperparameters
        if arg == "-sigma":
            last = "-sigma"
            continue
        if last == "-sigma":
            last = ""
            try:
                sigmas = [float(arg)]
            except ValueError:
                raise Exception(f"Sigma must be a number. Got {arg}. [EXITING]")
            continue
        if arg == "-lambda":
            last = "-lambda"
            continue
        if last == "-lambda":
            last = ""
            lambdas = [float(arg)] * len(sigmas)
            continue

        # Training size
        if arg == "-size":
            last = "-size"
            continue
        if last == "-size":
            last = ""
            try:
                if arg != "full":
                    size = int(arg)
            except ValueError:
                raise Exception(f"Size must be an integer. Got {arg}. [EXITING]")
            except Exception as e:
                print(e)
            continue

        # Databases
        if arg == "-train":
            last = "-train"
            continue
        if arg == "-test":
            last = "-test"
            continue
        if arg == "-eval":
            last = "-eval"
            continue
        opt_arg_list = ["-opt", "-optimize"]
        if arg in opt_arg_list:
            Qopt = 1
            last = "-opt"
            continue
        if arg == "-md":
            Qopt = 2
            last = "-opt"
            continue
        monomer_arg_list = ["-monomers", "-mon", "-atoms"]
        if arg in monomer_arg_list:
            last = "-monomers"
            continue
        if arg == "-trained":
            if Qtrain == 0:
                Qtrain = 2
            else:
                raise Exception("Cannot take trained if training [EXITING]")
            last = "-trained"
            continue

        # Min
        if arg == "-min":
            last = "-min"
            continue
        if last == "-min":
            last = ""
            try:
                Qmin = float(arg)
            except ValueError:
                raise Exception(f"Min must be a number. Got {arg}. [EXITING]")
            method = "min"
            continue

        # Split
        if arg == "-finishsplit":
            Qsplit = -1
        split_arg_list = ["-split", "-startsplit", "-finishsplit"]
        if arg in split_arg_list:
            last = "-split"
            continue
        if last == "-split":
            last = "-splitt"
            try:
                Qsplit = int(arg)
            except ValueError:
                raise Exception(f"Split must be an integer. Got {arg}. [EXITING]")
            continue
        if last == "-splitt":
            last = "-splittt"
            Qsplit_i = int(arg) - 1
            continue
        if last == "-splittt":
            last = ""
            Qsplit_j = int(arg) - 1
            continue

        # Sample each
        sample_each_arg_list = ["-sampleeach", "-se", "-categorize"]
        if arg in sample_each_arg_list:
            last = "-sampleeach"
            continue
        if last == "-sampleeach":
            last = ""
            if Qtrain == 0:
                Qtrain = 1
            if Qeval == 0:
                Qeval = 1
            try:
                Qsampleeach = int(arg)
            except ValueError:
                raise Exception(f"Sample each must be an integer. Got {arg}. [EXITING]")
            continue

        # No forces
        if arg == "-noforces":
            Qifforces = 0
            continue

        # Similarity
        similarity_arg_list = ["-similarity", "-sim"]
        if arg in similarity_arg_list:
            last = "-similarity"
            continue
        if last == "-similarity":
            last = ""
            if Qtrain == 0:
                Qtrain = 1
            if Qeval == 0:
                Qeval = 1
            try:
                Qsampleeval = -int(arg)
            except ValueError:
                raise Exception(f"Similarity must be an integer. Got {arg}. [EXITING]")
            continue

        # Laplacian
        if arg == "-laplacian":
            Qkernel = "Laplacian"
            continue

        # Train database(s)
        if last == "-trained":
            last = ""
            VARS_PKL = arg
            continue
        if last == "-train":
            TRAIN_HIGH = arg
            if Qtrain == 0 or Qsampleeach != 0:
                Qtrain = 1
            else:
                raise Exception("Cannot train if taking trained [EXITING]")
            last = "-train2"
            continue
        if last == "-train2":
            if os.path.exists(arg):
                method = "delta"
                TRAIN_LOW = arg
                last = ""
                continue

        # Test/Eval/Opt database(s)
        if last == "-eval":
            TEST_HIGH = arg
            Qeval = 1
            last = "-test2"
            continue
        if last == "-test":
            TEST_HIGH = arg
            Qeval = 2
            last = "-test2"
            continue
        if last == "-opt":
            TEST_HIGH = arg
            last = "-test2"
            continue
        if last == "-test2":
            if os.path.exists(arg):
                method = "delta"
                TEST_LOW = arg
                last = ""
                continue

        # Monomer database(s)
        if last == "-monomers":
            no_monomers_arg_list = ["0", "none", "no"]
            if arg in no_monomers_arg_list:
                Qmonomers = 2
                last = ""
            else:
                MONOMERS_HIGH = arg
                Qmonomers = 1
            last = "-monomers2"
            continue
        if last == "-monomers2":
            if os.path.exists(arg):
                method = "delta"
                MONOMERS_LOW = arg
                last = ""
                continue

        # Models and representations
        if arg == "-painn" or arg == "-nn":
            Qmethod = "nn"
            Qrepresentation = "painn"
            continue
        if arg == "-painn_EFQ":
            Qmethod = "nn"
            Qrepresentation = "painn"
            Qifcharges = 1
            continue
        if arg == "-painn_EFD":
            Qmethod = "nn"
            Qrepresentation = "painn"
            Qifdipole = 1
            Qifcharges = 1
            continue
        if arg == "-schnet":
            Qmethod = "nn"
            Qrepresentation = "schnet"
            continue
        if arg == "-physnet":
            Qmethod = "physnet"
            Qrepresentation = "physnet"
            Qifcharges = 1
            Qifdipole = 1
            continue
        if arg == "-so3net":
            Qmethod = "nn"
            Qrepresentation = "so3net"
            continue
        if arg == "-knn":
            Qmethod = "knn"
            continue
        if arg == "-aimnet":
            Qmethod = "aimnet"
            Qrepresentation = "aimnet"
            Qifeldisp = 1
            Qifforcedisp = 1
            Qifcharges = 1

        # override representation
        if arg == "-repr":
            last = "-repr"
            continue
        if last == "-repr":
            Qrepresentation = arg
            last = ""
            continue

        # turn off metric learning for k-NN
        if arg == "-nometric":
            no_metric = True
            continue

        # Epochs
        if arg == "-nn_epochs" or arg == "-epochs":
            last = "-nn_epochs"
            continue
        if last == "-nn_epochs":
            last = ""
            try:
                nn_epochs = int(arg)
            except ValueError:
                raise Exception(
                    f"Number of epochs must be an integer. Got {arg}. [EXITING]"
                )
            continue

        # Epochs
        if arg == "-nn_tvv" or arg == "-nn_train":
            last = "-nn_tvv"
            continue
        if last == "-nn_tvv":
            last = ""
            try:
                nn_tvv = float(arg)  # TODO Exceptions up from here
            except ValueError:
                raise Exception(
                    f"Number of epochs must be an integer. Got {arg}. [EXITING]"
                )
            continue

        # Radial basis
        if arg == "-nn_rbf" or arg == "-nn_rb" or arg == "-nn_basis":
            last = "-nn_rbf"
            continue
        if last == "-nn_rbf":
            last = ""
            try:
                nn_rbf = int(arg)
            except ValueError:
                raise Exception(
                    f"Number of radial basis must be an integer. Got {arg}. [EXITING]"
                )
            continue

        # NN cutoff
        cutoff_arg_list = ["-nn_cutoff", "-krr_cutoff", "-cutoff"]
        if arg in cutoff_arg_list:
            last = "-cutoff"
            continue
        if last == "-cutoff":
            last = ""
            try:
                nn_cutoff = float(arg)
                krr_cutoff = float(arg)
            except ValueError:
                raise Exception(f"Cutoff must be a number. Got {arg}. [EXITING]")
            continue

        # Opt max step
        max_step_arg_list = ["-opt_maxstep", "-opt_maxs"]
        if arg in max_step_arg_list:
            last = "-opt_maxstep"
            continue
        if last == "-opt_maxstep":
            last = ""
            try:
                opt_maxstep = float(arg)
            except ValueError:
                raise Exception(f"Opt max step must be a number. Got {arg}. [EXITING]")
            continue

        # Opt dump
        if arg == "-opt_dump":
            last = "-opt_dump"
            continue
        if last == "-opt_dump":
            last = ""
            try:
                opt_dump = int(arg)
            except ValueError:
                raise Exception(f"Opt dump must be an integer. Got {arg}. [EXITING]")
            continue

        # Opt steps
        if arg == "-opt_steps":
            last = "-opt_steps"
            continue
        if last == "-opt_steps":
            last = ""
            try:
                opt_steps = int(arg)
            except ValueError:
                raise Exception(
                    f"Number of opt steps must be an integer. Got {arg} [EXITING]"
                )
            continue

        # MD time step
        if arg == "-md_timestep":
            last = "-md_timestep"
            continue
        if last == "-md_timestep":
            last = ""
            try:
                md_timestep = float(arg)
            except ValueError:
                raise Exception(f"Timestep must be a number. Got {arg} [EXITING]")
            continue

        # MD thermostat friction
        if arg == "-md_thermostatfriction":
            last = "-md_thermostatfriction"
            continue
        if last == "-md_thermostatfriction":
            last = ""
            try:
                md_thermostatfriction = float(arg)
            except ValueError:
                raise Exception(
                    f"Thermostatfriction must be an integer. Got {arg}. [EXITING]"
                )
            continue

        # MD steps
        if arg == "-md_steps":
            last = "-md_steps"
            continue
        if last == "-md_steps":
            last = ""
            try:
                md_steps = int(arg)
            except ValueError:
                raise Exception(
                    f"Number of MD step must be an integer. Got {arg}. [EXITING]"
                )
            continue

        # MD dump
        if arg == "-md_dump":
            last = "-md_dump"
            continue
        if last == "-md_dump":
            last = ""
            try:
                md_dump = int(arg)
            except ValueError:
                raise Exception(f"MD dump must be an integer. Got {arg}. [EXITING]")
            continue

        # Atom basis
        atom_basis_arg_list = ["-nn_ab", "-nn_atom_basis", "-nn_features"]
        if arg in atom_basis_arg_list:
            last = "-nn_ab"
            continue
        if last == "-nn_ab":
            last = ""
            try:
                nn_atom_basis = int(arg)
            except ValueError:
                raise Exception(f"Atom basis must be an integer. Got {arg}. [EXITING]")
            continue

        # NN interactions
        nn_int_arg_list = ["-nn_int", "-nn_interctions", "-nn_blocks"]
        if arg in nn_int_arg_list:
            last = "-nn_int"
            continue
        if last == "-nn_int":
            last = ""
            try:
                nn_interactions = int(arg)
            except ValueError:
                raise Exception(
                    f"NN interactions argument must be an integer. Got {arg}. [EXITING]"
                )
            continue

        # NN interactions
        if arg == "-nw":
            last = "-nw"
            continue
        if last == "-nw":
            last = ""
            try:
                nw = int(arg)
            except ValueError:
                raise Exception(
                    f"NN interactions must be an integer. Got {arg}. [EXITING]"
                )
            continue

        # Energy tradeoff
        if arg == "-nn_energytradeoff":
            last = "-nn_energytradeoff"
            continue
        if last == "-nn_energytradeoff":
            last = ""
            try:
                nn_energytradeoff = float(arg)
            except ValueError:
                raise Exception(
                    f"Energy tradeoff must be a number. Got {arg}. [EXITING]"
                )
            continue

        # Checkpoint
        if arg == "-ckpt" or arg == "-chkp":
            last = "-ckpt"
            continue
        if last == "-ckpt":
            last = ""
            if os.path.isfile(arg):
                Qcheckpoint = arg
            else:
                raise Exception(
                    f"Path {arg} to the checkpoint does not exist [EXITING]"
                )
            continue

        # Early stop
        if arg == "-nn_ESpatience" or arg == "-nn_espatience":
            last = "-nn_espatience"
            continue
        if last == "-nn_espatience":
            last = ""
            try:
                Qearlystop = int(arg)
            except ValueError:
                raise Exception(f"Early stop must be an integer. Got {arg}. [EXITING]")
            continue

        # Learning rate
        if arg == "-nn_lr":
            last = "-nn_lr"
            continue
        if last == "-nn_lr":
            last = ""
            try:
                Qlearningrate = float(arg)
            except ValueError:
                raise Exception(f"Learning rate must be a number. Got {arg}. [EXITING]")
            continue

        # Batch size
        if arg == "-batch_size" or arg == "-bs":
            last = "-bs"
            continue
        if last == "-bs":
            last = ""
            try:
                Qbatch_size = int(arg)  # float 16.0 can be int 16
            except ValueError:
                raise Exception(f"Batch size must be an integer. Got {arg}. [EXITING]")
            continue

        # MD temperature
        if arg == "-md_temperature" or arg == "-temperature":
            last = "-md_temperature"
            continue
        if last == "-md_temperature":
            last = ""
            try:
                md_temperature = float(arg)
            except ValueError:
                raise Exception(f"Temperature must be number. Got {arg}. [EXITING]")
            continue

        # KRR (by default)
        krr_arg_list = ["-krr", "-fchl", "-qml"]
        if arg in krr_arg_list:
            Qmethod = "krr"
            Qrepresentation = "fchl"
            continue
        if arg == "-fchl19":
            Qmethod = "krr"
            Qrepresentation = "fchl19"
            continue
        if arg == "-mbdf":
            Qmethod = "krr"
            Qrepresentation = "mbdf"
            continue

        # Hyperparameter optimisation
        if arg == "-hyper":
            Qhyper = True
            Qtrain = 0  # set Qtrain to zero, as the "training" be part of the tuning
            Qeval = 0  # similarly as above
            continue

        # Load hyperparameter cache
        if last == "-hyper-cache":
            hyper_cache = arg
            last = ""
            continue
        if arg == "-hyper-cache":
            last = arg
            continue

        # Unknown argument
        raise Exception(f"Sorry cannot understand this argument: {arg} [EXITING]")

    if last != "":
        print(last)
        raise Exception(f"Hey looser, the last argument is incomplete. [EXITING]")

    return locals()
