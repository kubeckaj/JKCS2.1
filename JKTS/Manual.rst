============================================
About & Usage
============================================

.. contents:: Table of Contents
   :depth: 2

Introduction
============

``JKTS`` is a tool designed to facilitate a streamlined and automatic process of generating and monitoring files integral for transition state search in atmospheric chemistry, through the use of the quantum chemical programs ORCA and Gaussian16.
It implements multiconformer transition state theory (MC-TST) to achieve a realistic understanding and prediction of the kinetics of atmospheric chemical reactions between organic species and oxidants. ``JKTS`` supports hydrogen abstraction by the OH radical (``-OH``), the Cl radical (``-Cl``), and the NO3 radical (``-NO3``, nighttime chemistry).
Tunneling effects are accounted for with the Eckart tunneling correction factor, by fitting an unsymmetrical Eckart potential to the reactant, transition state, and product energies. For reactions involving light atoms, such as hydrogen, the tunneling correction has been shown to be crucial.

.. note::
   JKTS is ready to go once JKCS has been installed.


Workflow Overview
=================

Each working directory runs one of three predefined step sequences:

Transition State Workflow
-------------------------

#. Constrained preoptimization of the initial TS guess (active site frozen).
#. Transition state optimization towards a first-order saddle point, validated by normal-mode analysis.
#. Conformer sampling of the validated TS with CREST (active site constrained).
#. Constrained optimization of the sampled conformers.
#. Transition state optimization of the conformers.
#. Final energy refinement with DLPNO-CCSD(T) single points.

Reactants and Products Workflow
-------------------------------

#. Conformer sampling with CREST.
#. Geometry optimization of the conformers.
#. Final energy refinement with DLPNO-CCSD(T) single points.

The small partner molecules (OH/H2O, Cl/HCl, NO3/HNO3) skip conformer sampling and run only geometry optimization followed by DLPNO-CCSD(T); they are submitted automatically when the reactant or product batch finishes.

Using JKTS
==========

JKTS can be configured with various options for different computational scenarios. As an example, to calculate the reaction dynamics of hydrogen abstraction from methane by the OH radical, use the following command:

.. code-block:: bash

    JKTS CH4.xyz -OH

The input can also be given as a SMILES string with ``-smiles`` (e.g. ``JKTS -smiles C -OH`` for methane). This command creates three directories: **reactants**, **products**, and **CH4_H1**.

- The **reactants** directory contains the methane molecule (the OH radical is added automatically), following the workflow for reactants.
- The **products** directory contains the alkyl radical products of hydrogen abstraction (H2O is added automatically). For methane there is only one product since all hydrogens are equivalent.
- One **{molecule}_H{n}** directory is created per *unique* abstraction site. Chemically equivalent hydrogens are detected with RDKit symmetry-equivalence classes: only one hydrogen per class generates a transition state, and the class size is recorded as the reaction path degeneracy :math:`\sigma` of that channel (methane :math:`\sigma = 4`, ethane one TS with :math:`\sigma = 6`, a CH2 group :math:`\sigma = 2`). For chiral molecules, where the symmetry analysis is skipped, JKTS falls back to treating only methyl and formaldehyde hydrogens as equivalent.

Each directory is self-contained and holds:

- ``log`` — one timestamped log file recording state changes (submissions, running transitions, convergences, errors, step transitions).
- ``.metadata`` — a hidden JSON file with the run settings: reaction type, method/basis set/backend, SLURM resources, the original command line, and (for TS directories) the active-site indexes and :math:`\sigma`.
- ``{dir}_checkpoint.pkl`` — a crash-safe checkpoint of all molecules in the directory, rewritten atomically at every state change (a ``.bak`` of the previous good copy is kept).
- ``input_files/``, ``log_files/``, ``failed_logs/``, ``slurm_output/`` — job artifacts sorted as the workflow progresses.

Within each directory the workflow steps run as SLURM jobs. Conformer sampling uses the xtb program CREST [1]_ and results in a file named ``collection{molecule name}.pkl`` containing the sampled conformers; JKTS reads this file once it is created and generates input files for the next step. Geometry relaxations use the DFT method and basis set given by ``-method`` and ``-basis_set``. For transition state structures, every TS optimization is validated: the imaginary frequency must lie below ``-freq_cutoff`` (default -100 cm⁻¹) and its normal mode must correspond to hydrogen transfer at a geometrically sane active site; structures failing validation are corrected and resubmitted once before being dropped. Duplicate conformers are filtered out automatically with the ArbAlign program [4]_ after the constrained conformer optimization and after the DLPNO single points.

The final step is an energy correction with DLPNO-CCSD(T) [2]_ [3]_ (or CCSD(T)-F12 with ``-F12``) on the optimized geometries. Using the thermochemical contributions of the optimized geometries and the DLPNO single-point energies, the rate constant of each abstraction channel is calculated with the MC-TST equation:

.. math::
   k = \sigma \, \kappa \, \frac{k_b T}{h} \left( \frac{\sum_{i}^{\text{TS conf}} \exp\left(-\frac{\Delta E_i}{k_b T}\right) Q_{\text{TS},i}}{\sum_{j}^{\text{reac conf}} \exp\left(-\frac{\Delta E_j}{k_b T}\right) Q_{\text{reac},j}} \right) \exp\left(-\frac{E_{\text{TS},0} - E_{\text{R},0}}{k_b T}\right)

where :math:`\sigma` is the reaction path degeneracy of the channel and :math:`\kappa` the Eckart tunneling coefficient. The results are written to ``molecules.txt`` (per-molecule summary), ``Final_reactants_*.pkl`` / ``Final_products_*.pkl`` / ``Final_TS_*.pkl`` (result pickles), and ``Rate_constants.txt`` (per-channel table with the total rate constant).

.. _CREST: https://crest-lab.github.io/crest-docs/


Stopping, Restarting and Monitoring Time
----------------------------------------

For smaller molecules where the computational task is not as intensive, such as methane, the monitoring duration can be modified with the ``-time`` argument:

.. code-block:: bash

    JKTS CH4.xyz -OH -time 1:00:00

If the monitoring ends prematurely — the wall time ran out, the node died, or the run was stopped deliberately with ``-stop`` — the workflow can be resumed from the directory's checkpoint. ``cd`` into the affected directory (e.g. **CH4_H1**) and run:

.. code-block:: bash

    JKTS -restart

This reads ``{dir}_checkpoint.pkl`` (falling back to its ``.bak`` copy if the main file is corrupt), restores the original run settings from ``.metadata`` (reaction type, method, basis set, backend, SLURM resources — any flag given explicitly on the restart command line wins), and reconciles every molecule against the SLURM queue and the on-disk log files: finished jobs are recognized and advanced, still-queued jobs are reattached to, and jobs that vanished without finishing are resubmitted — queued jobs are never duplicated. A specific pickle or set of log files can also be given explicitly, e.g. ``JKTS CH4_H1_TS_opt.pkl -restart`` or ``JKTS *.log -restart``. A restart refuses to start while a previous JKTS monitor is still alive in the directory (tracked via ``monitor_pid`` in ``.metadata``).

The same resilience applies during a run: a job that disappears from the queue without a termination string in its log (node failure) is automatically resubmitted after a grace period, with its own per-molecule budget separate from the error-termination budget.

.. tip::
   If the user does not wish the program to automatically continue to the next step in the workflow, the ``-stop`` option is available: ``JKTS CH4.xyz -OH -stop``. The run then halts after the current step completes and can be resumed later with ``JKTS -restart``.


Advanced Usage
--------------

To run JKTS with specific settings, like a custom level of theory:

.. code-block:: bash

    JKTS yourfile.xyz -OH -method "B3LYP" -basis_set "6-31+g(d)"

Keep in mind the natural limitation of ORCA and Gaussian16 in relation to which basis sets and methods have been implemented in the respective programs. For methods that deploy their own basis set, such as B97-3c, r2scan-3c, and PM7, specifying a basis set is not needed. With ``-F12`` the final single points use CCSD(T)-F12/cc-pVTZ-F12 instead of DLPNO-CCSD(T).

Monitoring of log files
~~~~~~~~~~~~~~~~~~~~~~~~

JKTS monitors the log files at certain intervals to avoid overwhelming communication between computers. By default the program checks up to 100 times, with a time interval between checks calculated from the size of the input molecule and the type of job running. The maximum number of attempts, the interval, and the initial delay can be set manually:

.. code-block:: bash

    JKTS yourfile.xyz -OH -interval 500 -attempts 200 -initial_delay 2000

Resulting in an initial delay of 2000 seconds before the log files are checked, with 500 seconds between checks, performed up to 200 times.


Command Line Arguments
======================

``JKTS`` accepts various arguments to control its behavior:

.. list-table::
   :widths: 35 65
   :header-rows: 1

   * - Input Commands
     - Description
   * - ``-h``, ``--help``
     - Print help page
   * - ``-OH``
     - H abstraction by OH radical
   * - ``-Cl``
     - H abstraction by Cl radical
   * - ``-NO3``
     - H abstraction by NO3 radical (nighttime chemistry)
   * - ``-G16``
     - Use Gaussian16 for QC calculations (default)
   * - ``-ORCA``
     - Use ORCA for QC calculations
   * - ``-reactants``
     - Skip the reactants workflow (bare flag or explicit false) [def = run it]
   * - ``-products``
     - Skip the products workflow (bare flag or explicit false) [def = run it]
   * - ``-stop``
     - Stop after the current workflow step completes instead of continuing automatically (resume with ``-restart``)
   * - ``-restart``
     - Resume the workflow after a crash or ``-stop``: reattaches to queued jobs and resubmits lost ones. Reads the given .pkl/.log/.out files, or the ``{dir}_checkpoint.pkl`` in the current directory if no file is given
   * - ``-k``
     - Calculate the MC-TST rate constant at the end [def = True]
   * - ``-method``
     - QC method for optimization and TS search [def = wB97XD]
   * - ``-basis_set``
     - Basis set for the QC method [def = 6-31++g(d,p)]
   * - ``-F12``
     - Use CCSD(T)-F12/cc-pVTZ-F12 instead of DLPNO for single points
   * - ``--gfn``
     - GFN version for CREST (1 or 2) [def = 2]
   * - ``-ewin``
     - Energy window for CREST conformer sampling in kcal/mol [def = 5]
   * - ``-energy_cutoff``
     - Drop conformers this much above the lowest after preoptimization, in kcal/mol [def = 5]
   * - ``-max_conformers``
     - Max number of conformers taken from CREST [def = 1000]
   * - ``-freq_cutoff``
     - Cutoff for the TS imaginary frequency in cm^-1 [def = -100]
   * - ``-cpu``
     - Number of CPUs per job [def = 4]
   * - ``-mem``
     - Memory per job in MB [def = 8000]
   * - ``-par``
     - SLURM partition to use [def = q64]
   * - ``-time``
     - SLURM wall time for the monitoring job [def = 240:00:00]
   * - ``-interval``
     - Time interval between log file checks in seconds [def = based on molecule size]
   * - ``-initial_delay``
     - Initial delay before checking log files in seconds [def = based on molecule size]
   * - ``-attempts``
     - How many times the log files are checked [def = 100]
   * - ``-smiles``
     - Input molecule as a SMILES string instead of an .xyz file
   * - ``-info``
     - Print information of molecules in log files or .pkl files
   * - ``-movie``
     - Write movie.xyz of the given structures for viewing
   * - ``-collect``
     - Collect DFT thermochemistry and DLPNO single points into a .pkl file
   * - ``-pickle``
     - Store the given log files into a .pkl file
   * - ``-CHO``
     - Atom indexes of the active site (C H O), 1-indexed


Code Overview
=============

This section is the developer reference: the role of each file in ``src/`` and its most important functions.

The engine is a circular cycle in ``monitoring.py``: ``handle_termination`` prepares input files for the current workflow step and calls ``submit_and_monitor``, which submits a SLURM job (or array) and spawns a polling thread (``check_convergence`` or ``check_crest``); when all jobs of a step have converged, the poller calls ``handle_termination`` again for the next step. After the final DLPNO step, ``main.py`` assembles the MC-TST rate constant from the reactant, product and TS results.

src/main.py
-----------

Entry point: argument parsing and mode dispatch.

- ``build_parser()`` — the argparse CLI (the Command Line Arguments table above mirrors it).
- ``main()`` — mode dispatch: ``-init`` (build TS guesses, create directories, print the banner), ``-info``, ``-movie``; otherwise start or resume the workflow of the current directory, wait for all monitor threads, then write ``molecules.txt``, the ``Final_*.pkl`` result pickle, and — for TS directories with ``-k`` — compute and record the rate constant.
- ``read_input()`` — load the given .pkl/.log/.out/.com/.inp/.xyz files into ``Molecule`` objects, sorted by conformer number and capped at ``-max_conformers``.
- ``restore_settings_from_metadata()`` — on ``-restart``/``-rerun``, restore settings from ``.metadata`` for every flag the user left at its default; explicit flags always win.
- ``guard_and_write_pid_file()`` — record the monitor pid in ``.metadata`` and refuse to restart while a previous monitor is still alive.

src/classes.py
--------------

The data model. ``Step`` is a string-compatible enum of the workflow steps, defined next to the three workflow tuples (``TS_workflow``, ``reactant_product_workflow``, ``OH_H2O_workflow``). ``Vector`` provides static geometry helpers (distances, angles, Rodrigues rotation). ``Molecule`` holds one conformer/structure and its full workflow state: geometry, energies, frequencies, partition function, workflow position, SLURM job id and error counters, and the TS active site (``constrained_indexes``, 1-indexed dict with keys ``C``/``H``/``X`` and — for OH — ``XH``; ``reaction_path_degeneracy`` is :math:`\sigma`). Most important members:

- ``H_abstraction()`` — the TS-guess generator: for the representative H of each symmetry class, place the abstractor (OH, Cl, or NO3) at TS distances/angle (special parameters for aldehyde hydrogens), record the active site and :math:`\sigma`; optionally build the product radicals.
- ``equivalent_hydrogen_classes()`` — RDKit canonical ranks group abstractable hydrogens into symmetry classes; returns ``{representative_H: class_size}`` or ``None`` for chiral molecules (caller falls back to methyl/formaldehyde-only equivalence).
- ``find_active_site()`` / ``set_active_site()`` / ``perturb_active_site()`` — locate the C/H/X indexes (from ``-CHO``, ``.metadata``, or a geometric search) and repair/displace the active-site geometry before resubmissions.
- ``update_energy()`` — parse a G16/ORCA log for energies, frequencies, rotational data and dipole; ``partition_function()`` computes Q = q_vib·q_rot·q_trans·q_elec (with spin-orbit-corrected q_elec for OH and Cl).
- ``log2xyz()`` / ``log2program()`` / ``log2method()`` — extract geometry, program and method from log files; ``determine_workflow()`` / ``set_current_step()`` / ``update_step()`` manage the workflow position.
- ``create_OH/H2O/Cl/HCl/NO3/HNO3()`` — classmethods building the small partner molecules.
- ``save_to_pickle()`` / ``molecules_to_pickle()`` / ``load_from_pickle()`` — atomic pickling with backup fallback (via ``checkpoint``); ``move_files()`` / ``move_converged()`` / ``move_inputfile()`` / ``move_failed()`` tidy job artifacts; ``print_items()`` is the ``-info`` report.

src/monitoring.py
-----------------

Job monitoring and workflow engine. Module constants: ``FILTER_STEPS`` (steps after which conformers are deduplicated), ``VANISHED_GRACE_POLLS`` and ``MAX_NODE_FAILURES`` (node-failure handling).

- ``check_convergence()`` — the polling thread for QC steps: sweep the logs, classify each finished job (converged / error / TS-validation failure / vanished), resubmit or drop failures (two-error budget per molecule), checkpoint every state change, and advance the workflow when all jobs converged.
- ``check_crest()`` — the polling thread for CREST jobs: wait for the ``collection*.pkl`` files, resubmitting jobs that vanish without producing them.
- ``handle_termination()`` — step-transition hub: advance converged molecules, run ArbAlign filtering after ``FILTER_STEPS``, prepare inputs via ``STEP_HANDLERS``, checkpoint and submit.
- ``submit_and_monitor()`` — submit one job or an array, checkpoint the job id (so a crash can reattach), start the matching polling thread.
- ``termination_status()`` / ``handle_error_termination()`` / ``resubmit_job()`` — log classification and error triage (G16 convergence errors resubmit from the last geometry, intervention errors are dropped for manual inspection, SCF failures at a good active site get a geometry repair first).
- ``handle_input_molecules()`` / ``reconcile_group()`` — the restart entry point: group loaded molecules by step and reconcile each group against the SLURM queue and the on-disk logs (finished → advance, queued → reattach, lost → resubmit, never duplicating).
- ``submit_missing_partner()`` — submit the small partner molecule (OH/H2O etc.) after a reactant/product batch finishes DLPNO, if it is not already known.

src/slurm_submit.py
-------------------

HPC integration: writes per-job sbatch wrapper scripts (``{program}_submit.sh``, or ``{step}_submit.sh`` + ``array.txt`` for arrays) in the molecule directory and polls job status. Set ``JKTS_DRYRUN=1`` to print submit commands instead of submitting.

- ``submit_job()`` / ``submit_array_job()`` — write and run the CREST/G16/ORCA wrapper script; wall time is derived from the polling interval, and DLPNO memory escalates with the error count.
- ``update_molecules_status()`` / ``parse_job_statuses()`` — query ``squeue`` per unique job id and set each molecule's status (``running`` / ``pending`` / ``completed or not found`` / ``unknown`` when squeue itself fails — never treated as job loss).
- ``get_interval_seconds()`` — heuristic polling interval from heavy-atom count and step type.
- ``_run_submit_script()`` — runs sbatch with retries on transient QOS/submit-limit rejections.

src/checkpoint.py
-----------------

Crash-safe checkpointing; all pickle writes in the tree go through this module.

- ``atomic_pickle_dump()`` — temp file + fsync + ``os.replace``, previous copy rotated to ``.bak`` first, so a crash always leaves one readable copy.
- ``load_pickle_with_fallback()`` — load a pickle, falling back to ``.bak`` if the main file is missing or corrupt.
- ``save_checkpoint()`` / ``load_checkpoint()`` — the canonical per-directory ``{dir}_checkpoint.pkl``: the active batch merged with the already-finished molecules, written under a lock and never raising.

src/metadata.py
---------------

The single per-directory ``.metadata`` JSON file (see Using JKTS for the fields).

- ``load_metadata()`` — return the metadata dict, falling back read-only to the legacy ``.method``/``.constrain``/``.symmetry`` dotfiles of directories created by older versions.
- ``update_metadata()`` — atomic, lock-protected read-modify-write; never raises. The first call migrates legacy dotfile values into ``.metadata``.

src/qc_input.py
---------------

QC input file generation.

- ``QC_input()`` — write the .com/.inp input for the current step: G16 or ORCA route per step type (constrained opt / TS opt / plain opt / DLPNO single point), frozen active-site coordinates when constrained, TS-mode and Hessian settings (tighter on retries), DLPNO memory escalation, F12 variant with ``-F12``.
- ``crest_constrain()`` — write ``constrain.inp`` for a TS CREST run (C–H and H–X distances and the C–H–X angle fixed).
- ``mkdir()`` — create a working directory with its subfolders and record the run settings in ``.metadata``.

src/ts_validation.py
--------------------

TS geometry/frequency validation after a TS optimization converges.

- ``check_transition_state()`` — the full verdict: imaginary frequency below ``-freq_cutoff``, then the relative C–H/H–X bond-length changes along the imaginary mode must exceed a threshold (looser for aldehyde TSs) at a geometrically sane active site. Returns ``(passed, message)``.
- ``good_active_site()`` — distance/angle sanity check of the active site (wider thresholds for Cl, looser angles for aldehydes).
- ``extract_normal_coordinates()`` — imaginary-mode displacement vectors from a G16 log; ``is_aldehyde()`` classifies the abstraction site.

src/conformer_tools.py
----------------------

Conformer filtering and loading. ArbAlign is not a JKTS file: the JKQC installer copies ``TOOLS/SCRIPTS/modifiedArbAlign.py`` into the JKQC venv's site-packages as ``ArbAlign.py``, and JKTS imports it with ``from ArbAlign import compare``.

- ``filter_molecules()`` — RMSD/energy/dipole-based duplicate removal (keeping the lower-energy conformer of each match), applied automatically after the steps in ``monitoring.FILTER_STEPS``.
- ``energy_cutoff()`` — drop conformers more than ``-energy_cutoff`` kcal/mol above the lowest.
- ``initiate_conformers()`` — expand a CREST ``collection*.pkl`` (a JKQC pandas DataFrame) into ``Molecule`` conformers sorted by energy.
- ``collect_DFT_and_DLPNO()`` — pair DFT-optimized structures with their DLPNO single points and merge the energies (the ``-collect`` utility).

src/rate_constant.py
--------------------

MC-TST rate constants with tunneling.

- ``rate_constant()`` — split off the radical and small product by name, Boltzmann-sum the conformer partition functions, compute the tunneling coefficient, and return a ``RateResult`` with k scaled by :math:`\sigma`.
- ``eckart()`` — Eckart tunneling coefficient :math:`\kappa` from the forward/reverse barriers and the imaginary frequency; returns 1 (with a warning) on failure.

src/results.py
--------------

Result formatting: ``RateResult`` (namedtuple with k, :math:`\kappa`, Ea, partition functions, :math:`\sigma`, T), ``format_rate()`` (the single canonical rate-constant string used in the log, on stdout and in ``Rate_constants.txt``), ``record_rate()`` (maintains the machine-readable ``.rates.tsv`` and rewrites the ``Rate_constants.txt`` summary), and ``write_molecule_summary()`` (``molecules.txt``).

src/output.py
-------------

Logging. One ``Logger`` per working directory writes timestamped, severity-tagged lines to the directory's ``log`` file and echoes state changes to the terminal. Severities: ``info`` (file-only detail), ``event``/``success`` (file + stdout), ``warning``/``error`` (prefixed, errors to stderr), ``results`` (framed between ``=`` rules). Thread-safe (per-logger file lock + process-wide terminal lock). ``console`` is a module-level terminal-only logger for code without directory context; ``banner()`` prints the startup summary.

src/runtime.py
--------------

Shared mutable state, set once by ``main()`` after argument parsing: ``args``, ``QC_program``, the per-program termination/error strings, ``global_molecules`` (+ lock) holding the molecules that finished the workflow, and ``start_dir``.



Citations
=========

  .. [1] https://pubs.rsc.org/en/content/articlelanding/2020/CP/C9CP06869D
  .. [2] C. Riplinger and F. Neese, "An efficient and near linear scaling pair natural orbital based local coupled cluster method," J. Chem. Phys., vol. 138, p. 034106, 2013
  .. [3] C. Riplinger, B. Sandhoefer, A. Hansen, and F. Neese, "Natural triple excitations in local coupled cluster calculations with pair natural orbitals," J. Chem. Phys., vol. 139, p. 134101, 2013.
  .. [4] B. Temelso, J. M. Mabey, T. Kubota, N. Appiah-Padi, and G. C. Shields, "Arbalign: A tool for optimal alignment of arbitrarily ordered isomers using the kuhn-munkres algorithm," Journal of Chemical Information and Modeling, vol. 57, no. 5, pp. 1045-1054, 2017.
