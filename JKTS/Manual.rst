============================================
About & Usage
============================================

.. contents:: Table of Contents
   :depth: 2

Introduction
============

``JKTS`` is a tool designed to facilitate a streamlined and automatic process of generating and monitoring files integral for transition state search in atmospheric chemistry, through the use of quantum chemical programs ORCA and Gaussian16.
It implements the multiconformer transition state theory (MC-TST) to achieve realistic understanding and prediction of the kinetics in atmospheric chemical reactions between organic species and oxidants in the atmosphere. Currently, ``JKTS`` supports the  reaction of hydrogen abstraction by OH radical and OH radical addition to carbon-carbon double bonds. 
Furthermore, tunneling effects are accounted for using the Eckart tunneling correction factor, by fitting a unsymmetrical Eckart potential to the reactant complex, transition state, and product complex energies. In the case of reactions involving light atoms, such as hydrogen, the tunneling correction has shown to be crucial.

.. note::
   JKTS would be ready to go once JKCS has been installed
   

Workflow Overview
=================

The JKTS tool processes transition state molecules and reactants/products using distinct workflows:

Transition State Molecules Workflow
-----------------------------------

#. Conformer sampling with CREST.
#. Constrained optimization of the conformers.
#. Transition state optimization.
#. Final energy refinement with DLPNO-CCSD(T) calculations.

Reactants and Products Workflow
-------------------------------

#. Conformer sampling with CREST.
#. Geometry optimization of the conformers.
#. Final energy refinement with DLPNO-CCSD(T) calculations.

.. note::
   Insert a workflow diagram here.

Using JKTS
==========

JKTS can be configured with various options for different computational scenarios. As an example, to calculate the reaction dynamics of hydrogen abstraction from methane with the OH radical, use the following command:

.. code-block:: bash

    JKTS CH4.xyz -OH

This command will create three directories: **products**, **reactants**, and **CH4_H1**.

- The **products** directory includes the H2O molecule and configurations of products from hydrogen abstraction. For methane, there is only one product type since all hydrogens are equivalent.
- The **CH4_H1** directory is for the case of single hydrogen abstraction from methane. THe JKTS program tries to treat chemically equivalent hydrogens. A rather simple but efficient approach is used by storing the indexes of carbons with three (or more) hydrogens bonded to it during the generation of the initial guess for the TS state. These are assumed to correspond to the molecules methyl groups and that these hydrogens are chemically equivalent. This assumption is justified based on the initial use of the CREST program for the conformer sampling and therefore the subsequent sampling of the methyl groups.
- The **reactants** directory contains the methane molecule and the OH radical, following the workflow for reactants.

Within each directory, the workflow begins with conformer sampling using the xtb program, CREST [1]_. This process results in a file named `collection{molecule name}.pkl`, which contains the sampled conformers. JKTS reads this file once its created and generates input files for the next step of the workflow. All steps and progress inside each directory is logged to a file ``"log"``.

Post conformer sampling, the default subsequent step is geometry relaxation, performed with the DFT method and basis set given by ``-method`` and ``-basis_set``. For transition state (TS) structures, a constrained preoptimization of the active site precedes the optimization towards a first-order saddle point.

The final step involves an energy correction using DLPNO-CCSD(T) [#]_ [#]_, applied to the geometries derived from the geometry relaxations of product and reactant structures, as well as TS structures from first order saddle point optimization. Utilizing the optimized geometries' thermochemical contributions and single-point energy from DLPNO calculations, the rate constant for specific abstraction sites is calculated using the MC-TST equation:

.. math::
   k = \kappa \frac{k_b T}{h} \left( \frac{\sum_{i}^{\text{TS conf}} \exp\left(-\frac{\Delta_{E_i}}{k_b T}\right) Q_{\text{TS},i}}{\sum_{j}^{\text{reac conf}} \exp\left(-\frac{\Delta_{E_j}}{k_b T}\right) Q_{\text{reac},j}} \right) \exp\left(-\frac{E_{\text{TS},0} - E_{\text{R},0}}{k_b T}\right)
   

.. _CREST: https://crest-lab.github.io/crest-docs/


Adjusting Monitoring Time and Restarting Jobs
---------------------------------------------

For smaller molecules where the computational task is not as intensive, such as for methane, the monitoring duration can be modified with the ``-time`` argument. To set the monitoring time to five hours, the following can be specified:

.. code-block:: bash

    JKTS CH4.xyz -OH -time 5:00:00

However, imagine this wasn't quite enough time and the monitoring ended prematurely, for instance during the optimization towards a first order saddle point. For methane we would ``cd`` into the **CH4_H1** directory. In here a number of log files will exist for the current step in the workflow that is proceeding. That is if the monitoring ended before the program counted all the log files to have converged and therefor wasn't able to have move them into their respective folders. However, we can simply restart the calculations from their last set of geometries with the command:

.. code-block:: bash

    JKTS *.log

The wildcard symbol (*) matches all `.log` files in the directory. The JKTS program will go through all the passed log files and access which have terminated correctly or with and error and also those who perhaps didnt finish within their allowed wall-time.
Log files with error termination of who didnt finish will restart from the last geometry from the log file and log files deemed to have converged correctly will be waiting for the non-converged log files to finish. The workflow will subsequently resume from there on.


.. tip::
	If the user does not wish the program to automatically continue to the next step in the workflow for the submitted job, the ``-stop`` option is available:        ``JKTS *.log -stop``
       
    
Advanced Usage
--------------

To run JKTS with specific settings, like a custom level of theory:

.. code-block:: bash

    JKTS yourfile.xyz -OH -method "B3LYP" -basis_set "6-31+g(d)"
    
Keep in mind the natural limitation of ORCA and Gaussian16 in relation to which basis sets and methods have been implemented into the respective programs. For the case of methods who utilize a self deployed basis set, such as B97-3c, r2scan-3c, and PM7, the need for specifying basis set is not needed.

Monitoring of log files
~~~~~~~~~~~~~~~~~~~~~~~~    

JKTS monitors the log file with certain intervals to avoid overwhelming communication between computers. By default the program allows this communication a `100` times with a certain time interval between each check determined by ``interval``. By default the time between checks is calculated based on the size of the input molecule and the current job running. However, the maximum number of attempts to check the log files and the interval between them can be manually set with command line arguments:

.. code-block:: bash

    JKTS yourfile.xyz -OH -interval 500 -attempts 200 -initial_delay 2000
    
Resulting in an initial delay of 2000 seconds before the log files are checked with 500 seconds interval between each check and this check is performed up to 200 times.


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
     - Perform H abstraction with OH radical
   * - ``-Cl``
     - Perform H abstraction with Cl radical
   * - ``-NO3``
     - Perform H abstraction with NO3 radical (nighttime chemistry)
   * - ``-G16``
     - Gaussian16 is used for QC calculations (default)
   * - ``-ORCA``
     - ORCA is used for QC calculations
   * - ``-reactants``
     - Skip the reactants workflow (bare flag or explicit false) [def = run it]
   * - ``-products``
     - Skip the products workflow (bare flag or explicit false) [def = run it]
   * - ``-stop``
     - Stop after the current workflow step completes instead of continuing automatically (resume with ``-restart``)
   * - ``-restart``
     - Resume the workflow from the current step of the given .pkl/.log/.out files
   * - ``-k``
     - Calculate Multiconformer Transition State rate constant [def = True]
   * - ``-method``
     - QC method for optimization and TS search [def = wB97XD]
   * - ``-basis_set``
     - Basis set for the QC method [def = 6-31++g(d,p)]
   * - ``-cpu``
     - Number of CPUs [def = 4]
   * - ``-mem``
     - Amount of memory allocated for job [def = 8000mb]
   * - ``-par``
     - Partition to use [def = q64]
   * - ``-time``
     - SLURM wall time for the monitoring job [def = 240:00:00]
   * - ``-interval``
     - Set time interval between checks of log files [def = based on molecule size]
   * - ``-initial_delay``
     - Set an initial delay before checking log files [def = based on molecule size]
   * - ``-attempts``
     - Set how many times a log files should be checked [def = 100]
   * - ``-max_conformers``
     - Set max number of conformers from CREST [def = 1000]
   * - ``-freq_cutoff``
     - Set cutoff for TS imaginary frequency to [int] cm^-1 [def = -100]
   * - ``-ewin``
     - Set energy threshold to [int] kcal/mol for CREST conformer sampling [def = 5]
   * - ``-info``
     - Print information of molecules in log files or .pkl file

Duplicate conformers are filtered out automatically with the ArbAlign program [#]_ after the constrained conformer optimization and after the DLPNO single points.


Citations
=========
  
  .. [#] https://pubs.rsc.org/en/content/articlelanding/2020/CP/C9CP06869D
  .. [#] C. Riplinger and F. Neese, “An efficient and near linear scaling pair natural orbital based local coupled cluster method,” J. Chem. Phys., vol. 138, p. 034106, 2013
  .. [#] C. Riplinger, B. Sandhoefer, A. Hansen, and F. Neese, “Natural triple excitations in local coupled cluster calculations with pair natural orbitals,” J. Chem. Phys., vol. 139, p. 134101,2013.
  .. [#] B. Temelso, J. M. Mabey, T. Kubota, N. Appiah-Padi, and G. C. Shields, “Arbalign: A tool for optimal alignment of arbitrarily ordered isomers using the kuhn–munkres algorithm,” Journal of Chemical Information and Modeling, vol. 57, no. 5, pp. 1045–1054, 2017.
 
 


