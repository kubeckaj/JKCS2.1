#! /usr/bin/perl -w
#-------------------------------------------------

################################################################################
##
##                       ACDC
##          Atmospheric Cluster Dynamics Code
##   by Oona Kupiainen-Määttä and Tinja Olenius 2012-
##
##
##   Table of contents
##   #################
##
##   1. Processing command-line input
##
##   2. Reading the input file
##      2.1 Reading the info on the molecule types
##      2.2 Reading the nucleation rules
##      2.3 Reading the coagulation sink cluster
##      2.4 Reading the cluster types
##      2.5 Adding molecule types for the proton and removed proton if needed
##      2.6 Determining the system size for the loop option
##      2.7 Adding cluster numbers for generic ions if they are used
##      2.8 Adding cluster numbers for other fluxes
##      2.9 Checking which molecule types are used in the system
##
##   3. Dealing with forbidden reactions 
##      and reactions with nonstandard products
##
##   4. Charge balance for generic positive or negative ions
##
##   5. Reading the Delta G file
##      5.1 Reading Delta G (or H and S) for each cluster and possibly hydrates
##      5.2 Calculating hydrate distributions
##
##   6. Making the equations
##      6.1 Looping over the clusters and going through all possible collisions
##      6.2 Save the collision term if the collision can happen
##      6.3 Save the corresponding evaporation if it can happen
##      6.4 External losses
##      6.5 Finding out the number of MCMC coefficients
##      6.6 Time dependent source terms
##      6.7 Reactions inside clusters
##
##   7. Reading the cluster radius file
##
##   8. Starting to print the output files
##      8.1 Cluster properties for Matlab driver and Fortran system file
##      8.2 dC/dt for Fortran
##      8.3 Determining if parameters for ions are used in the Matlab files
##      8.4 dC/dt for Matlab
##      8.5 Fluxes for Matlab
##      8.6 Formation rate for Fortran
##      8.7 Formation rate for Matlab
##      8.8 Arrays for collisions, evaporations and losses for Fortran
##      8.9 Masses, radii and molecular content for the loop version
##
##   9. Loss coefficients
##      9.1 Coagulation sink
##         9.11 Coagulation sink parameterizations
##            9.111 Size-dependent sink of the form A*(d/dref)^B,
##                  where d is the diameter
##            9.112 Size-dependent sink based on a constant concentration
##                  and size of background particles
##         9.12 Coagulation sink given by one number
##         9.13 Printing out size-dependent sink in dry case
##         9.14 Hydrate average for size-dependent sink
##      9.2 Wall losses
##         9.21 Reading in wall losses for all the clusters from a given file
##         9.22 Wall loss parameterizations
##            9.221 Cloud 4 wall loss parametrization from Jasper's exel file
##            9.222 Cloud 4 wall loss parametrization from Andreas K\"urten
##            9.223 Cloud 4 wall loss, simplified parametrization
##            9.224 Cloud 3 wall loss parametrization
##            9.225 Diffusion wall loss for a flow tube (originally for the tube of Hanson and Eisele)
##         9.23 Wall losses given by one number
##            9.231 Ift flow tube wall loss
##            9.232 Other size-independent wall loss
##            9.233 Wall losses for specific clusters from the command line
##         9.24 Printing out size-dependent wall losses in dry case
##         9.25 Hydrate average for size-dependent wall losses
##
##   10. Collision rates
##      10.1 Reading in rates for specific collisions
##      10.2 Reading in sticking factors
##      10.3 Reading in dipole moments and polarizabilities
##           for collision enhancement factors
##      10.4 Calculating the collision rates
##         10.41 Determining if this collision has been given
##               a fixed rate value in the collision rate file
##         10.42 Printing a function for calculating
##               the collision rate when using loops
##         10.43 Calculating the hard sphere collision rate
##      10.5 Calculating the ion collision enhancement factors
##         10.51 Su and Bowers, 1973
##         10.52 Su and Chesnavich, 1982
##         10.53 Constant enhancement for all collisions
##      10.6 Hydrate averaged collision rate including ion enhancement
##      10.7 Figuring out sticking factors for each collision
##      10.8 Checking if this is an MCMC coefficient
##      10.9 Printing out the collision rates
##
##   11. Evaporation rates
##      11.1 Function for computing evaporation rates when using loops
##			11.11 Evaporation rates from DeltaGs
##			11.12 Evaporation rates from the Kelvin formula
##      11.2 Reading in rates for specific evaporations
##      11.3 Reading in scaling factors for the evaporation rates
##      11.4 Calculating the evaporation rates
##         11.41 Determining if this evaporation has been given
##               a fixed rate value in the evaporation rate file
##         11.42 Calculating the free energy -based evaporation rate
##            11.421 Figuring out scaling factors for each evaporation
##      11.5 Hydrate averaged evaporation rate
##      11.6 Temperature dependent evaporation rate
##      11.7 Printing out the evaporation rates
##
##   12. Printing the ion loss enhancement factors for Fortran
##
##   13. "True" coefficient values for MCMC
##
##   14. Finishing the Matlab driver file
##
##   15. Subroutines
##
################################################################################

# Modules
use strict;  # makes us declare all variables
use Getopt::Long;  # command line input
#use 5.010;   # use a newer Perl version
no if $] >= 5.018, warnings =>  "experimental::smartmatch";
use feature "switch";
use List::Util qw(max min sum);  # use some basic functions

my(@comments,@label,%num_from_label,%mol_num_from_name);
my($max_cluster,$string,$iclus);
my($ncols,$this_clus,$icol,$i,$j,$k,$l,$n,$ii,$jj,$kk,$ll,$jc,$kc,$lc,$combined);
my($line,$max_terms_per_line);
my($recomb_coeff,$jclus,$mass1,$mass2,$vol1,$vol2,$lmonomer1,$lmonomer2);
my(@temp_array,@temp_array_2);
my($valid,$radius,$string1,$string2,$string3,$string4);
my(@columns);
my($temperature,$imol,$jmol,$kmol,$nmol,$nmol_types);
my($iline,@input_file,$number_of_lines,$nend,$iend,@mass,@density,@psat,@surface_tension,$sigma);
my($temp_val,$temp_val2,$temp_val3,$size_array,@delta_g,$boltz,$Na,$kcal_per_mol_to_J,@evap_coef,$coll_coef_temp);
my($pi,@coll_coef,@sticking_coef,$volume,$mass_conv,$mass_nion,$dens_nion,$mass_pion,$dens_pion,$evap_coef_temp);
my(@lmonomer);
my(@evap_data,$pressure);
my($driver_file_name,$driver_function_name,$input_file_name,$output_file_name,$output_function_name,$ldo_table_out,$include_negative_ion_terms,$include_positive_ion_terms,$disable_coag_sinks,$coag_terms,$coag_terms_in,$disable_wall_terms,$wall_terms,$wall_terms_in,$disable_evap,@disable_evap_only,$free_energy_file_name,$disable_nonmonomers,$disable_nonmonomer_evaps,$enable_rec_nonmonomers,@nonmonomer_collisions_only,@no_derivative_clus,@no_derivative);
my($lfound);
my($disable_excluded,$disable_excluded_ionization);
my($cs_term_temp,$wl_term_temp,$found_negative,$found_positive,$diameter,$ion_coll_method);
my($kclus);
my($disable_dilution,$l_use_get_losses);
my($iline2,$iline0);
my(@molecule_name,@neutral_mols,$n_neutral_mols,@charge,@corr_neutral_mol,@corr_negative_mol,@corr_positive_mol,@n_corr_neutral_mol, @n_corr_negative_mol, @n_corr_positive_mol,@lnegative_cluster,@lpositive_cluster,@lnegative, @lpositive,@lneutral_cluster,@lneutral,$proton,$missing_proton,$generic_negative_ions,$generic_positive_ions);
my($corr_mol);
my($max_cluster_number);
my($wanted_mol,$valid_evap,$valid_coll);
my($flux_file_name,$flux_function_name,$j_file_name,$j_function_name,$input_for_eq,$input_for_eq_C_fun,$input_for_flux,$input_for_flux_temp,$n_generic_neg,$n_generic_pos,$n_flux_out,$n_flux_out_neg,$n_flux_out_pos,$n_flux_out_max,$n_flux_coag,$n_flux_wall,$n_flux_dil,$n_flux_bound,$n_flux_source,$n_flux_rec,$max_n_flux,$l_save_outgoing,$size_C_out);
my($l_all_output,$l_increase_time,$help);
my($first_neg,$first_pos,@n_max_type,@n_max_type_neg,@n_max_type_pos);
my($temperature_inp,$l_change_temperature,$element_name,$lout);
my($l_charge_balance,$generic_positive_ions_corr,$generic_negative_ions_corr);
my(%lines);
my($temp_label,$temp_label2,$label,@temp_composition,$lsame);
my(%nmonomers,$monomer_clus,@n_coll_products);
my(@nucl_rules,@nucl_rules_neg,@nucl_rules_pos,$n_nucl_rules,$n_nucl_rules_neg,$n_nucl_rules_pos,$irule);
my($dip_file_name,@dip_data,$c_mon,$c_clust,@dip,@pol,$nclus,@fcr);
my($nonstandard_reaction_file_name,$l_nonstandard_allowed,$l_nonstand,@nonstandard_reaction_data,@nonstandard_reactions);
my($collision_coef_file_name,@collision_coef_in,@lcoll_coef_in,$evaporation_coef_file_name,@evaporation_coef_in,@levap_coef_in,$wl_coef_file_name);
my($sticking_factor_file_name,$l_sticking_factor,@sticking_factor,$wanted_clus);
my($scale_evap_file_name, $l_scale_evap, @scale_evap_factor,@evap_coef_scale);
my($l_bulk_density,$radius_file_name,@radius_data,@radius_in,$r1,$r2);
my($lnitrate);
my($air_viscosity,$viscosity_N2,$radius_N2,$mass_N2,$D_acid_N2,$D0_factor,$D_factor,$D_to_wl,$tube_radius,$wl_term_0,$tube_radius_in,$tube_pressure_in);
my($exp_loss_coefficient,$exp_loss_exponent,$exp_loss_ref_cluster,$default_monomer,$lfound_several,$lfound_default);
my($bg_concentration,$bg_diameter,$bg_density,$Diff1,$Diff2,$c1,$c2,$Kn);
my($lfortran,$output_file_name_fortran,$system_file_name_fortran,%coef_quad, %coef_quad_form,%coef_lin);
my(%coef_reac,@ind_reac_form,@ind_reac_loss);
my(@ind_quad_loss,@ind_quad_form,@ind_quad_loss_extra,@ind_quad_form_extra, @ind_lin_loss,@ind_lin_form, %l_quad_bound,$l_coll_factor_mcmc_all,$l_coll_factor_mcmc_n_n,$l_coll_factor_mcmc_n_i,$l_coll_factor_mcmc_n_g,$l_coll_factor_mcmc_rec,$l_coll_factor_mcmc_any,$l_q_neg_mcmc,$l_q_pos_mcmc,$l_wl_mcmc,@wl_mcmc_only,@wl_only,@wl_only_hydr,$l_mcmc_any,$n_mcmc_coefs, $missing_energies,$missing_fcr,$l_mcmc_coef,@l_special_mcmc_wl,$icoef_wl);
my($icoef_coll,@l_coll_not_used,$nmolecules,@molecules);
my($lhydrate_average,$rh,$psw,$name_water,$mass_water,$dens_water,$volume_water,$iwater,$nwaters,$lfound_hydr,$ihydr,$jhydr,$khydr,$nhydr,@lhas_hydr,@delta_g_hydr,@delta_g_temp, @distr,@distr_temp,@distr_temp2,@distr_temp3,$sum,$sum2);
my(@coll_coef_hydr,@evap_coef_hydr,@fcr_hydr,@dip_hydr,@pol_hydr,@coll_coef_eff,@evap_coef_eff,%radius_in_hydr,@cs_coef_hydr,@wl_coef_hydr,@hydrates_str);
my($l_evap_factor_mcmc,$icoef_evap,$ncoef_evap,$ncoef_coll,$ncoef_coll_n_n,$ncoef_coll_n_i,$ncoef_coll_n_g,$ncoef_coll_rec,$ncoef_wl);
my($dil_value,$dil_str,$use_coag_sinks,$use_wall_terms,$use_dilution,$include_generic_ions,$no_generic_ions,$include_generic_neg,$include_generic_pos,$no_generic_neg,$no_generic_pos);
my($variable_temp,$variable_cs,$cs_value,$cs_value_in,$variable_ion_source);
my($append_to_file_names,$append_to_subroutine_names,$label_generic_neg, $label_generic_pos,$label_proton,$label_missing_proton,@coef_true);
my(@negative_clusters,@positive_clusters);
my($stopping);
my($iclus_label,$jclus_label,$kclus_label,$disable_useless_collisions,$keep_useless_collisions,$coef_temp,$one_cs_for_all,$one_wl_for_all,$one_wl_for_almost_all,$include_ion_terms,$line_temp,$fcr_temp);
my(%l_quad_nst,@clust_num_of_monomer,$fcs,$fwl,$evap_comment);
my($parameters_in_get_fcr_1,$parameters_in_get_fcr_2, $parameters_in_get_fcr_3,$parameters_in_get_fcr_4);
my(@base_strength,@acid_strength,$l_strength_neutral,$acid_max,$base_max,$l_print_boundary);
my($str_temp,$str_temp_form,$str_temp2,$str_temp3,$value_string,$coll_string,$hydr_string,$evap_string,$fixed_string,$fcr_string);
my(@cs_cluster,@cs_cluster_neg,@cs_cluster_pos,$use_coag_to_outgoing, $n_coag_on_formed_clusters,$l_coag,$cs_cluster_label,$cs_cluster_neg_label,$cs_cluster_pos_label);
my($keep_boundary_clusters,$max_n_flux_inc_boundary_clusters,@label_inc_boundary_clusters, $kclus_bound, %clus_bound_ind,$lfound_bound,$print_3d_flux);
my($l_old_output,$ref,$comp_ok);
my(@clust_comp,@clust_mass,@clust_radius,@clust_volume);
my(@l_can_be_lost_bound,@l_can_evaporate);
my($clust_mass_max,$diameter_max,$mob_diameter_max);
my($lloop,$lloop_diag,$lloop_max_ratio,$lloop_cs,$nmol_types_used,@nmol_ranges_array,$nmol_ranges_list,$nmol_ranges_ind,$cluster_from_indices,@cluster_from_indices_array,$monomer_indices,@source_function,@source_function_clust,@mol_types_used,@mol_names_used,$loop_coll_coef,$loop_evap_coef,$rlim_no_evap,$rlim_no_evap_str,$lloop_j_frac,$lloop_j_tests,$lgrowth_for_j_tests,$c_growth_compound,$c_growth_compound_str,$rlim_growth_compound_str);
my($max_ratio_basis,$lfree_mon,$nmol_inp_max_ratio_basis,$nmol_inp_not_basis,$nmol_max_ratio_basis,$nmol_not_basis,$max_wrt_mon);
my($l_jlim,$nlim,$jlim_sizes_list,$l_j_in,$j_in_function);
my($reactions_in_clusters_file,@cluster_reactions,$n_react_in_clust,$l_reac_in_clus);
my($n_react,$n_prod,@comp_react,@comp_prod,$rate_fwd,$rate_bck,@needed_mols,@forbidden_mols);
my($n_reacting_W,$n_formed_W,$n_catalyzing_W,$l_forbidden_W);

################################################################################
################################################################################
##   1. Processing command-line input                                         ##
################################################################################
################################################################################

# All variables beginning with l are treated as Boolean logical variables: 1 is true, 0 is false

# These should not be changed
$label_generic_neg='neg';
$label_generic_pos='pos';

##### Initialize input variables #####

# Output of this Perl script (formats etc.)
$lfortran = 0;
$l_all_output = 1;
$l_old_output = 0;
$ldo_table_out = 1;
$disable_nonmonomers = 0;
@nonmonomer_collisions_only = ();
$enable_rec_nonmonomers = 0;
$disable_evap = 0;
$disable_useless_collisions = 1;
$keep_useless_collisions = 0;
$l_print_boundary = 0;
$keep_boundary_clusters = 0;
$print_3d_flux = 0;
@no_derivative_clus = ();
$l_increase_time = 0;

$lloop = 0;
$lloop_diag = 0;
$lloop_max_ratio = 0;
$max_ratio_basis = '';
$lfree_mon = 0;
$lloop_cs = 0;
$loop_coll_coef = 'hard_spheres';
$loop_evap_coef = 'DeltaG';
$rlim_no_evap = '';
$lloop_j_frac = 0;
$l_jlim = 0;
$jlim_sizes_list = '';
$l_j_in = 0;
$j_in_function = '';
$lloop_j_tests = 0;
$lgrowth_for_j_tests = 0;
$c_growth_compound = 1e7; # cm^-3

$driver_file_name = 'driver_acdc.m';
$output_file_name = 'equations_acdc.m';
$output_file_name_fortran = 'acdc_equations.f90';
$system_file_name_fortran = 'acdc_system.f90';
$flux_file_name = 'dofluxes.m';
$j_file_name = 'formationrate.m';
$append_to_file_names = '';
$append_to_subroutine_names = '';
$max_terms_per_line = 10;

$l_coll_factor_mcmc_any = 0;
$l_coll_factor_mcmc_all = 0;
$l_coll_factor_mcmc_n_n = 0;
$l_coll_factor_mcmc_n_i = 0;
$l_coll_factor_mcmc_n_g = 0;
$l_coll_factor_mcmc_rec = 0;
$l_q_neg_mcmc = 0;
$l_q_pos_mcmc = 0;
$l_wl_mcmc = 0;
@wl_mcmc_only = ();
@wl_only = ();
$l_evap_factor_mcmc = 0;

$l_change_temperature = 0;

@source_function = ();

# Default input files or variables of this Perl script
$input_file_name = 'input_AD.inp';
$free_energy_file_name = 'G278K.txt';
$dip_file_name = 'dip_pol_278K.txt';

$rh = 0;

$radius_file_name = '';
$nonstandard_reaction_file_name = '';
$reactions_in_clusters_file = '';
$l_reac_in_clus = 0;

$collision_coef_file_name = '';
$evaporation_coef_file_name = '';
$sticking_factor_file_name = '';
$sticking_factor[0][0] = 1;
$sticking_factor[1][0] = 1;
$scale_evap_file_name = '';
$scale_evap_factor[0][0] = 0;

# Variables related to external sinks
$disable_coag_sinks = 0;
$disable_wall_terms = 0;
$disable_dilution = 0;
$use_coag_sinks = 0;
$use_coag_to_outgoing = 0;
$use_wall_terms = 0;
$use_dilution = 0;

$coag_terms_in = '';                # Coagulation losses
$cs_value_in = '';
$cs_value = "2.6e-3";
$coag_terms = 'exp_loss';			# Options: exp_loss, bg_loss, number in 1/s
$exp_loss_coefficient = $cs_value;
$exp_loss_exponent = -1.6;
$exp_loss_ref_cluster = '';
$default_monomer = 'A';
$bg_concentration = 1e3; # cm^-3
$bg_diameter = 100; # nm
$bg_density = 1000; # kg/m^3

$wall_terms_in = '';                # Wall losses
$wl_coef_file_name = '';
$wall_terms = 'CLOUD4';				# Options: cloud4, cloud4_exel_file, cloud4_Andreas, cloud3, ift, diffusion, number in 1/s
$tube_radius_in = '';
$tube_pressure_in = '';

$dil_value = 9.6e-5;				# Dilution from CLOUD: dil = flow / vol / 60000

# Variables related to ions
$include_generic_ions = 1;
$no_generic_ions = 0;
$include_generic_neg = 0;
$include_generic_pos = 0;
$no_generic_neg = 0;
$no_generic_pos = 0;
$l_charge_balance = 0;				# 1 (-1) -> set [pos] ([neg]) to match overall negative (positive) charge, 0 -> use source terms for both
$ion_coll_method = 'Su82';			# Options: Su82, Su73, constant
$fcs = 1;
$fwl = 3.3;
$recomb_coeff = 1.6e-12;

# Variables related to parameters that can be set to have non-fixed values
$variable_temp = 0;
$variable_cs = 0;
$variable_ion_source = 0;


##### Read in the command line input #####

#getopts('d:o:e:f:wit:xcm', \%opts); # for Getopt::Std
GetOptions('help'=>\$help, 'fortran'=>\$lfortran, 'all_output|a!'=>\$l_all_output, 'old_output_order|old_output'=>\$l_old_output, 
'disable_nonmonomers|no_nonmonomers|m'=>\$disable_nonmonomers, 'nonmonomer_only|cluster_collisions_only=s' => \@nonmonomer_collisions_only, 
'enable_rec_nonmonomers'=>\$enable_rec_nonmonomers, 'disable_evap|no_evap'=>\$disable_evap, 
'disable_nonmonomer_evaps|no_nonmonomer_evaps|me'=>\$disable_nonmonomer_evaps, 'all_collisions|keep_useless_collisions'=>\$keep_useless_collisions, 
'print_boundary'=>\$l_print_boundary, 'keep_boundary_clusters'=>\$keep_boundary_clusters, 'print_3d_flux|flux_3d' => \$print_3d_flux, 
'disable_excluded'=>\$disable_excluded, 'disable_excluded_ionization'=>\$disable_excluded_ionization, 'save_outgoing'=>\$l_save_outgoing, 'no_derivative_clus|no_eq=s' => \@no_derivative_clus, 'increase_time'=>\$l_increase_time, 
'loop' => \$lloop, 'loop_diag' => \$lloop_diag, 'loop_max_ratio=s' => \$max_ratio_basis, 'free_mon' => \$lfree_mon, 
'loop_cs' => \$lloop_cs, 'loop_coll_coef=s' => \$loop_coll_coef, 'loop_evap_coef=s' => \$loop_evap_coef, 'rlim_no_evap|evaplim=f' => \$rlim_no_evap, 
'loop_j_frac|j_frac' => \$lloop_j_frac, 'jlim=s' => \$jlim_sizes_list, 'j_in' => \$l_j_in, 'j_in_function=s' => \$j_in_function, 
'loop_j_tests' => \$lloop_j_tests, 'growth_for_j_tests' => \$lgrowth_for_j_tests, 'c_growth_compound=f' => \$c_growth_compound, 
'driver_file_name|d=s' => \$driver_file_name, 'output_file_name|o=s' => \$output_file_name,
'system_file_name|syst=s' => \$system_file_name_fortran, 'flux_file_name|f=s' => \$flux_file_name,
'append_to_file_names|append=s'=>\$append_to_file_names, 'append_to_subroutine_names|append_to_subroutines=s'=>\$append_to_subroutine_names, 
'coll_coef_mcmc|coll_coef_mcmc_all' => \$l_coll_factor_mcmc_all, 
'coll_coef_mcmc_neutral|coll_coef_mcmc_neutral_neutral|coll_coef_mcmc_n_n' => \$l_coll_factor_mcmc_n_n, 
'coll_coef_mcmc_ion|coll_coef_mcmc_neutral_ion|coll_coef_mcmc_n_i' => \$l_coll_factor_mcmc_n_i, 
'coll_coef_mcmc_charger|coll_coef_mcmc_neutral_charger|coll_coef_mcmc_n_g' => \$l_coll_factor_mcmc_n_g, 
'coll_coef_mcmc_rec' => \$l_coll_factor_mcmc_rec, 'q_neg_mcmc' => \$l_q_neg_mcmc, 'q_pos_mcmc' => \$l_q_pos_mcmc,
'wl_mcmc' => \$l_wl_mcmc, 'wl_mcmc_only=s' => \@wl_mcmc_only, 'wl_only=s' => \@wl_only, 'evap_coef_mcmc' => \$l_evap_factor_mcmc, 
'temperature|t=s' => \$temperature, 'source_function=s' => \@source_function,
'cluster_set_file_name|cluster_set|c|input_file_name|i=s' => \$input_file_name, 'free_energy_file_name|evap_file_name|e=s' => \$free_energy_file_name, 
'dip_file_name|dip=s' => \$dip_file_name, 'rh=f' => \$rh, 'radius_file_name|rad=s' => \$radius_file_name, 
'nonstandard_reaction_file|nst=s' => \$nonstandard_reaction_file_name, 
'reactions_in_clusters_file|cluster_reaction_file|clust_reac=s'  => \$reactions_in_clusters_file, 
'collision_coef_file_name=s' => \$collision_coef_file_name, 'evaporation_coef_file_name=s' => \$evaporation_coef_file_name,
'sticking_factor_file_name|enhancement_factor_file_name=s' => \$sticking_factor_file_name, 
'sticking_factor|enhancement_factor=f' => \$sticking_factor[0][0], 'sticking_factor_ion_neutral=f' => \$sticking_factor[1][0],
'scale_evap_file_name=s' => \$scale_evap_file_name, 'scale_evap_factor=f' => \$scale_evap_factor[0][0], 
'disable_coag_sinks|disable_coag'=>\$disable_coag_sinks, 'disable_wall_terms|disable_wall_losses|disable_wl'=>\$disable_wall_terms, 
'disable_dilution'=>\$disable_dilution, 'use_coag_sinks|use_coag|use_cs'=>\$use_coag_sinks, 'use_wall_terms|use_wall_losses|use_wl'=>\$use_wall_terms, 
'use_dilution|use_dil'=>\$use_dilution, 'coag_terms|cs=s'=>\$coag_terms_in, 'cs_value=f'=>\$cs_value_in, 
'exp_loss_coefficient=f' => \$exp_loss_coefficient, 'exp_loss_exponent=f' => \$exp_loss_exponent, 'exp_loss_ref_cluster=s' => \$exp_loss_ref_cluster, 
'bg_concentration|bg_c=f' => \$bg_concentration, 'bg_diameter|bg_d=f' => \$bg_diameter, 'bg_density|bg_rho=f' => \$bg_density, 
'wall_terms|wl=s'=>\$wall_terms_in, 'wl_coef_file_name|wlfile=s'=>\$wl_coef_file_name, 
'flowtube_radius=f' => \$tube_radius_in, 'flowtube_pressure=f' => \$tube_pressure_in, 'dil_value=f' => \$dil_value, 
'no_generic_ions|no_gen'=>\$no_generic_ions, 'no_generic_neg|no_neg_charger'=>\$no_generic_neg, 'no_generic_pos|no_pos_charger'=>\$no_generic_pos, 
'nitrate'=>\$lnitrate, 'charge_balance|cb=i' => \$l_charge_balance, 'ion_coll_method=s'=> \$ion_coll_method, 'fcs=f' => \$fcs, 'fwl=f' => \$fwl, 
'variable_temp'=>\$variable_temp, 'variable_cs'=>\$variable_cs, 'variable_ion_source'=>\$variable_ion_source);

if(!$l_all_output){
	die "Option noall_output has been disabled since no-one was using it.\n";
}

if($help){
		die "Usage: \n\t\tacdc.pl [--options]\n
		Input and output options:\n
		\t--print_boundary
		\t\tprint out how each cluster that has grown too far in a wrong direction is brought back to the boundary\n
		\t--keep_boundary_clusters
		\t\tsave the fluxes to and from each boundary cluster in the Matlab version (to track what exactly is going to the boundary and how it comes back)\n
		\t--no_eq clust
		\t\tdon't print out the derivative/arrays for it for cluster 'clust', and set the concentration to be constant\n
		\t--fortran
		\t\tprint the code in fortran instead of matlab\n
		\t--driver_file_name name
		\t--d
		\t\tmatlab function for running the simulation\n
		\t--output_file_name name
		\t--o
		\t\tmatlab/fortran file containing the birth-death equations\n
		\t--flux_file_name name
		\t--f
		\t\tmatlab/fortran file containing the flux equations\n
		\t--system_file_name name
		\t--syst
		\t\tfortran file specifying the system\n
		\t--append_to_file_names ending
		\t--append
		\t\tappend the same ending to all output file names, including get_coll etc.\n
		\t--free_energy_file_name name
		\t--evap_file_name
		\t--e
		\t\tfile containing cluster free energies\n
		\t--input_file_name name
		\t--i
		\t\tfile specifying the system\n
		\t--dip_file_name name
		\t--dip
		\t\tfile containing dipole moments (D) and polarizabilities (Å^3) of neutrals\n
		\t--radius_file_name name
		\t--rad
		\t\tfile containing radii of molecules and clusters\n
		\t--old_output_order
		\t--old_output
		\t\told output order: converged is the last output argument (this was the default up to version 2014-11-03)\n
		\t--save_outgoing
		\t\tconcentration output includes concentration of clusters lost from the system by nucleation, coagulation and wall-losses\n
		\t--all_collisions
		\t--keep_useless_collisions
		\t\tprint equations also for reactions of the form x + y -> x + y; by default these are left out\n
		\t--loop
		\t\tuse a loop in the Fortran codes instead of writing out all equations and rate constants
		\t\tNote: use this option only for a big box-shaped system - no fancy boundary conditions are possible, except the loop_max_ratio option\n
		\t--loop_diag
		\t\tsame as loop, but for a 2-component system with growth only along the diagonal\n
		\t--loop_max_ratio name
		\t\tsame as loop, but for a 2-component system in which the principal building block is the molecule with the given name,\n
		\t\tand the maximum number of molecules of the other type is according to what is given in the input file
		\t\t(e.g. largest cluster = 10A20N -> max. N:A ratio 2:1, when A is the building block)
		\t--free_mon
		\t\tdisable the formation of clusters containing one building block molecule and one or more ligand molecules when using loop_max_ratio,
		\t\ti.e. the clustering process only starts with the formation of the dimer consisting of two building block molecules
		\t--loop_cs
		\t\tsame as loop, but outgoing collision products stay at the boundary instead of going out,
		\t\tcan also be used together with --loop_diag or --loop_max_ratio\n\n
		\t--loop_j_frac or --j_frac
		\t\trecord the fractions of different growth routes over the given size limits jlim
		\t--jlim list
		\t\tlist of cluster sizes (radius in nm, separated by commas) over which the net formation flux will be recorded, e.g. 0.6,1.5,3
		\t\tcan be used only in the loop mode
		\t--j_in
		\t--j_in_function function_name
		\t\tuse a source (formation) rate given by function_name (optional) for some specific size, and omit the distribution below this size
		\t\t(i.e. don't consider it in the birth-death equations), with the exception of monomers
		\t\tcan be used only in the loop mode
		\t--loop_coll_coef
		\t\tfunction for calculating collision coefficients in the loop mode
		\t\tavailable functions/parameterizations: hard_spheres, Dahneke
		\t--loop_evap_coef
		\t\tfunction for calculating evaporation coefficients in the loop mode
		\t\tavailable functions/parameterizations: DeltaG, Kelvin
		\t--rlim_no_evap | evaplim
		\t\tcluster size (radius in nm) at and beyond which evaporation is set to zero (works only for the loop mode)
		Losses:\n
		\t--disable_coag_sinks
		\t\tdisables coagulation sinks\n
		\t--use_coag_sinks
		\t\tuses coagulation sinks\n
		\t--disable_wall_terms
		\t\tdisables wall losses\n
		\t--use_wall_terms
		\t\tuses wall losses\n
		\t--disable_dilution
		\t\tdisables dilution losses\n
		\t--use_dilution
		\t\tuses dilution losses\n
		\tNote: if some losses are disabled, all others are used, and if some losses are used, all others are disabled. You should either specify which losses to use or which losses to disable, but not both.\n
		\t--coag_terms
		\t--cs
		\t\tavailable coagulation loss parameterizations: exp_loss, bg_loss, value in 1/s\n
		\t--exp_loss_coefficient
		\t\tconstant pre-factor A for wall loss option 'exp_loss', which is of the form A*(d/dref)^B, where d is the diameter\n
		\t--exp_loss_exponent
		\t\texponent B for wall loss option 'exp_loss', which is of the form A*(d/dref)^B, where d is the diameter\n
		\t--exp_loss_ref_cluster
		\t\treference cluster for dref for wall loss option 'exp_loss', which is of the form A*(d/dref)^B, where d is the diameter\n
		\t--bg_concentration
		\t--bg_c
		\t\tconcentration of background particles in cm^-3\n
		\t--bg_diameter
		\t--bg_d
		\t\tdiameter of background particles in nm\n
		\t--bg_density
		\t--bg_rho
		\t\tdensity of background particles in kg m^-3\n
		\t--wall_terms
		\t--wl
		\t\tavailable wall loss parameterizations: cloud4, cloud4_exel_file, cloud4_Andreas, cloud3, ift, diffusion, value in 1/s\n
		\t--flowtube_radius
		\t\tflowtube radius in cm for wall loss option 'diffusion'\n
		\t--flowtube_pressure
		\t\tpressure inside flowtube in Pa for wall loss option 'diffusion'\n\n
		\t--wl_coef_file_name name
		\t--wlfile name
		\t\tfile containing wall loss rates for all the clusters in s^-1 (if no built-in parameterization is used)
		\t\tformat: <loss rate> <cluster>\n
		\t--wl_only clust,value
		\t\tuse for cluster 'clust' the wall loss 'value' instead of the parameterization or value used for other clusters\n\n
		\t--dil_value
		\t\tdilution loss value in 1/s (default: 9.6e-5)\n
		MCMC options:\n
		\t--evap_coef_mcmc
		\t\tmcmc coefficients for evaporation rates\n
		\t--coll_coef_mcmc_neutral
		\t\tmcmc coefficients for neutral-neutral collision rates\n
		\t--coll_coef_mcmc_ion
		\t\tmcmc coefficients for neutral-ion collision rates (including charger ions)\n
		\t--coll_coef_mcmc_rec
		\t\tmcmc coefficients for ion-ion recombination rates\n
		\t--coll_coef_mcmc_charger
		\t\tmcmc coefficients for neutral-charger ion collision rates\n
		\t--coll_coef_mcmc
		\t\tmcmc coefficients for all collision rates\n
		\t--l_q_neg_mcmc
		\t\tmcmc coefficients for negative charger ion source term\n
		\t--l_q_pos_mcmc
		\t\tmcmc coefficients for positive charger ion source term\n
		\t--l_wl_mcmc
		\t\tmcmc coefficients for wall loss term (same for all clusters, no charge enhancement)\n
		\t--wl_mcmc_only clust
		\t\tuse mcmc wall losses only for cluster 'clust'\n\n
		Other simulation physics options:\n
		\t--source_function clust[,function name]
		\t\tuse for cluster 'clust' a time-dependent source function; the function name is optional and if it's not given, a default name will be used\n
		\t--disable_evap or --no_evap
		\t\tdisables all evaporations
		\t--no_generic_ions
		\t--no_gen
		\t\tdisables generic ions\n
		\t--no_generic_neg
		\t--no_neg_charger
		\t\tdisables generic negative ions\n
		\t--no_generic_pos
		\t--no_pos_charger
		\t\tdisables generic positive ions\n
		\t--charge_balance
		\t--cb
		\t\t1 (-1) -> set positive (negative) charger ion concentration to match overall negative (positive) charge\n
		\t--disable_boundary
		\t\tdisables collisions and ionization/recombination resulting in clusters outside of what we are simulating unless they lead to 'nucleation'\n
		\t--disable_boundary_ionization
		\t\tdisables ionization/recombination resulting in clusters outside of what we are simulating unless they lead to 'nucleation'\n
		\t--disable_excluded
		\t\tdisables collisions and ionization/recombination resulting in clusters outside of what we are simulating\n
		\t--disable_excluded_ionization
		\t\tdisables ionization/recombination if we don't have the correponding ionic/neutral cluster in the system\n
		\t--increase_time
		\t\tincrease the length of the simulation if it doesn't reach a steady state (default: run the simulation for a fixed time)\n
		\t\tthis can also be accomplished with the keyword 'repeat' when running the Matlab driver\n
		\t--ion_coll_method
		\t\tdictates how we enhance ion-neutral collisions (Su73, Su82, constant, small)\n
		\t--fwl
		\t\tion enhancement for wall losses (default = 3.3)\n
		\t--fcs
		\t\tion enhancement for coagulation sink (default = 1)\n
		\t--nitrate
		\t\tuse nitrate dimers mass and density for negative charger ions\n
		\t--disable_nonmonomer_evaps
		\t--no_nonmonomer_evaps
		\t--me
		\t\tonly monomers are allowed to evaporate\n
		\t--disable_nonmonomers
		\t--no_nonmonomers
		\t--m
		\t\tonly monomer collisions and evaporations\n
		\t--nonmonomer_only clust
		\t--cluster_collisions_only clust
		\t\texclude cluster-cluster collisions, except when one of the partners in 'clust' (you can define several of these)\n
		\t--enable_rec_nonmonomers
		\t\tallow recombination of clusters when otherwise only monomer processes are considered
		\t--nonstandard_reaction_file name
		\t--nst
		\t\tcollisions with non-standard products (eg. breaking due to energy nonaccommodation)
		\t\tformat e.g.: 4A4N 1A .3 out_neu .7 4A4N .7 1A\n
		\t--rh number
		\t\tcalculates the equilibrium hydrate distribution for all the clusters for which hydrate energies are available at the given rh (in percent), and uses the average hydrate collision and evaporation rates\n
		\t--collision_coef_file_name name
		\t\tfile containing rates for specific collisions in m^3 s^-1; these are used as they are, i.e. sticking factors, fcr's etc. are not added
		\t\tformat: <collision rate> <collider 1> <collider 2>\n
		\t--evaporation_coef_file_name name
		\t\tfile containing rates for specific evaporations in s^-1; these are used as they are, i.e. scaling factors, fcr's etc. are not added
		\t\tformat: <evaporation rate> <mother cluster> <evaporator (optional, otherwise the rate applies to all evaporations of the mother cluster; which is useful e.g. if the rate is set to zero)>\n		
		\t--sticking_factor number
		\t\tadds a sticking coefficient to neutral-neutral collisions\n	
		\t--sticking_factor_ion_neutral number
		\t\tadds a sticking coefficient to ion-neutral collisions\n
		\t--sticking_factor_file_name name
		\t\tfile containing scaling for collision rates, given as sticking probabilities
		\t\tformat: <sticking probability> <optional (otherwise the probability is added to all neutral-neutral collisions): molecule name label (then the probability is added to all collisions involving the type) or colliding cluster> <optional: if the previous column is a cluster name, the other colliding cluster can be given here if the probability is only for this specific collision>\n
		\t--scale_evap_factor number
		\t\tscaling for evaporation rates, given as a correction to the reaction free energy change\n
		\t--scale_evap_file_name name
		\t\tfile containing scaling for evaporation rates, given as a correction to the reaction free energy change
		\t\tformat: <change to Delta G> <evaporator> <mother cluster (optional, otherwise evaporations from all clusters)>\n
		\t--temperature
		\t--t
		\t\ttemperature (mainly if you are using H and S at a different temperature to get G, or don't have evaporations at all)\n
		\t--variable_temp
		\t\tprints the collision and evaporation rates and charge enhancement factors as functions taking in the temperature, so simulations can be done at different temperatures without rerunning the perl script\n\n";
}

##### Update the input variables after reading the command line #####

### Determining which losses are used ###
if($cs_value_in ne ''){
	$cs_value=$cs_value_in;
	if($coag_terms_in ne ''){
		die "\nGive coagulation losses either as the name of the parameterization (or a single number) or as the variable cs_value, but not both!\n\n";
	}else{
		$coag_terms_in=$cs_value_in;
	}
}
if($wl_coef_file_name ne ''){
	if($wall_terms_in ne ''){
		die "\nGive wall losses either as the name of the parameterization (or a single number) or as the name of the file containing the losses, but not both!\n\n";
	}else{
		$wall_terms_in='wlfile';
	}
}
if(($use_coag_sinks+$use_wall_terms+$use_dilution+$variable_cs+$l_wl_mcmc)>0 || $coag_terms_in ne '' || $wall_terms_in ne '' || $#wl_only>=0 || $#wl_mcmc_only>=0){
	if(($disable_coag_sinks+$disable_wall_terms+$disable_dilution)>0){
		die "\nSpecify which losses to use or which losses to disable, but not both!\n\n";
	}
	if($use_coag_sinks==0 && $variable_cs==0 && $coag_terms_in eq ''){
		$disable_coag_sinks = 1;
		
	}else{
		print "Using coagulation sinks.\n";
	}
	if($use_wall_terms==0 && $l_wl_mcmc==0 && $#wl_mcmc_only<0 && $wall_terms_in eq '' && $#wl_only<0){
		$disable_wall_terms = 1;
		
	}elsif(($wall_terms_in ne '') && $l_wl_mcmc==1){
		die "\nSpecify your wall loss term or use an mcmc coefficient for it, but not both!\n\n";
	}elsif(($wall_terms_in eq '') && $#wl_only>=0){
		$wall_terms_in = '0.0';
		print "Using wall losses (but only for some clusters).\n";
	}elsif($#wl_only>=0 && ($l_wl_mcmc || $#wl_mcmc_only>=0)){
		die "wl_only cannot be used together with mcmc wall losses";
	}else{
		print "Using wall losses.\n";
	}
	if($use_dilution==0){
		$disable_dilution = 1;
		
	}else{
		print "Using dilution.\n";
	}
}elsif(($disable_coag_sinks+$disable_wall_terms+$disable_dilution)==0){
	$disable_coag_sinks = 1;
	$disable_wall_terms = 1;
	$disable_dilution = 1;
}
if($coag_terms_in ne ''){
	$coag_terms = $coag_terms_in;
}
if($wall_terms_in ne ''){
	$wall_terms = $wall_terms_in;
}
if(&is_pos_number($coag_terms)){
	$one_cs_for_all = 1;
}else{
	$one_cs_for_all = 0;
}
if($l_wl_mcmc==1 || &is_pos_number($wall_terms) || $wall_terms =~ m/ift/i){
	if($#wl_mcmc_only>=0){
		$one_wl_for_all = 0;
		$one_wl_for_almost_all = 1;
	}elsif($#wl_only>=0){
		$one_wl_for_all = 0;
		$one_wl_for_almost_all = 1;
	}else{
		$one_wl_for_all = 1;
		$one_wl_for_almost_all = 0;
	}
}else{
	$one_wl_for_all = 0;
	$one_wl_for_almost_all = 0;
}
# Determining if we are using some losses that go into get_losses
if((!$disable_coag_sinks && (!$one_cs_for_all || !$variable_cs)) || !$disable_wall_terms || ($lfortran && !$disable_dilution)){
	$l_use_get_losses = 1;
}else{
	$l_use_get_losses = 0;
}
# Unit conversions etc. related to some of the loss options
if(!$disable_coag_sinks && $coag_terms =~ m/^bg_loss$/i){
	$bg_concentration *= 1e6; # cm^-3 -> m^-3
	$bg_diameter *= 1e-9; # nm -> m
}
if(!$disable_dilution){
	$dil_str = sprintf('%.14e',$dil_value);
	if($lfortran){
		$dil_str =~ s/e/d/;
	}
}

### Determining if all generic ions are disabled ###
if($no_generic_ions){
	$include_generic_ions = 0;
	$include_generic_neg = 0;
	$include_generic_pos = 0;
	$no_generic_neg = 1;
	$no_generic_pos = 1;
	print "Not using generic ions.\n";
}

### Settings related to disabled cluster-cluster processes
if($disable_nonmonomers && !$disable_nonmonomer_evaps){
	$disable_nonmonomer_evaps = 1;
}
# Determining if cluster-cluster processes are allowed only when recombination is involved ###
if($enable_rec_nonmonomers && !$disable_nonmonomers){
	$disable_nonmonomers = 1;
	print "Disabling cluster-cluster processes except for collisions involving recombination.\n";
}

### Checking the option for disabling the printing out of some derivatives ###
if($#no_derivative_clus>=0 && !$lfortran){
	die "The option to not print out derivatives is only for the Fortran version, for now.";
}

### Determining if the equations are written as loops ###
if($max_ratio_basis ne ''){
	$lloop_max_ratio = 1;
}
if($lloop_diag || $lloop_max_ratio || $lloop_cs){
	$lloop = 1;
}
if($lfree_mon && !$lloop_max_ratio){
	die "\nOption --free_mon works only with --loop_max_ratio.\n\n";
}
# Disabling cluster fissions if the Kelvin evaporation is used in the loop mode
if($lloop && $loop_evap_coef =~ m/^Kelvin$/i && !$disable_nonmonomer_evaps){
	$disable_nonmonomer_evaps = 1;
	print "Disabling cluster fissions.\n";
}

### Determining if reactions in clusters are used ###
if($reactions_in_clusters_file ne ''){
	$l_reac_in_clus = 1;
	if($lfortran){
		die "\nThe reactions in clusters option hasn't yet been implemented to the fortran version.\n\n";
	}
	if($variable_temp){
		die "\nThe reactions in clusters option cannot be used with a variable temperature.\n\n";
	}
}

### Determinig if the net formation fluxes over some specific size limits are recorded ###
# (works only for the loop mode)
if($jlim_sizes_list ne ''){
	if(!$lloop){
		die "\nSaving the fluxes over the given size limits is a feature of the loop mode\n";
	}elsif($lloop_diag){
		die "\nSaving the fluxes over the given size limits is not set to work for loop_diag right now,\nsince I have no idea what's going on with that option\n";
	}else{
		$l_jlim = 1;
		@temp_array=split(',',$jlim_sizes_list);
		$nlim = $#temp_array + 1;
		if($lfortran){
			for(@temp_array){
				if($_ =~ m/\./){
					$_ .= "d0";
				}else{
					$_ .= ".d0";
				}
			}
			$jlim_sizes_list = join(',',@temp_array);
		}
	}
}

### Determinig if the distribution is simulated only above some "formation size", for which a source function is used ###
# (works only for the one-component loop mode)
if($j_in_function ne ''){
	$l_j_in = 1;
}
if($l_j_in){
	if(!$lloop){
		die "\nUsing an input formation rate for some specific size is a feature of the loop mode\n";
	}
	# Add the "formation size" to the list of source functions
	$str_temp = "n_j_in";
	if($j_in_function ne ''){
		$str_temp .= ",$j_in_function";
	}
	push(@source_function,$str_temp);
}
# This is only used in the J validation study (on Triolith/Fafner), and it
# (1) prints the output as in the tests for J(t) vs. J(steady-state)/J(extrapolated steady-state) etc.,
# [(2) uses the specific modules etc. to get the bg concentration,] -> this is now a default for the loop mode
# (3) considers evaporation only up to a specific size,
# (4) adds the irreversible condensation of a so-called implicit growth compound
#     to particles larger than a specific threshold size
if($lloop_j_tests){
	print "\nPrinting the output for the J validation tests!\n";
	$lgrowth_for_j_tests = 1;
}
if($lgrowth_for_j_tests){
	$rlim_no_evap = 0.75; # nm

	$rlim_growth_compound_str = '1.5d0*1.d-9'; # nm -> m
	
	$c_growth_compound *= 1e6; # cm^-3 -> m^-3
	$c_growth_compound *= 0.5; # assuming two different masses (each have the same number concentration so that Ctot=$c_growth_compound)
	$c_growth_compound_str = sprintf('%.14e',$c_growth_compound);
	$c_growth_compound_str =~ s/e/d/;
	#$c_growth_compound_str = '5.d6*1.d6'; # cm^-3 -> m^-3
	
	print "Using implicit 'growth compounds' with C = $c_growth_compound_str m^-3,\n";
	print "that condense onto particles r >= $rlim_growth_compound_str m,\n";
	print "also restricting evaporation to sizes of r < $rlim_no_evap nm\n";
}

### Determining if a limit size for evaporation is used
if($rlim_no_evap !~ /^\s*$/){
	if(!$lloop){
		die "\nOption --rlim_no_evap works only for the loop mode.\n\n";
	}
	$rlim_no_evap_str = sprintf('%.14e',$rlim_no_evap*1e-9); # nm -> m
	$rlim_no_evap_str =~ s/e/d/;
}

### Determining if useless collisions are kept in the output files ###
if($keep_useless_collisions == 1){
	$disable_useless_collisions = 0;
}

### Determining if clusters outside the system boundaries are tracked ###
# (works only for the Matlab version)
if($keep_boundary_clusters && $lfortran){
	die "\nSaving the boundary fluxes is a feature of the Matlab version\n";
}

### Determining if bulk densities are used ###
if($radius_file_name eq ''){
	$l_bulk_density = 1;
}else{
	$l_bulk_density = 0;
}

### Names of the output files ###
if($lfortran){
	if($output_file_name ne 'equations_acdc.m'){
		$output_file_name_fortran = $output_file_name;
	}
}
# Appending an ending to all output file names
if($append_to_file_names ne ''){
	if($lfortran){
		$output_file_name_fortran =~ s/\.f90$/$append_to_file_names\.f90/;
		$system_file_name_fortran =~ s/\.f90$/$append_to_file_names\.f90/;
	}else{
		$driver_file_name =~ s/\.m$/$append_to_file_names\.m/;
		$output_file_name =~ s/\.m$/$append_to_file_names\.m/;
		$flux_file_name =~ s/\.m$/$append_to_file_names\.m/;
		$j_file_name =~ s/\.m$/$append_to_file_names\.m/;
	}
}
# Extracting the Matlab function names from the files names (removing the .m ending)
if(!$lfortran){
	$driver_function_name = $driver_file_name;
	$driver_function_name =~ s/\.\S*$//;
	
	$output_function_name = $output_file_name;
	$output_function_name =~ s/\.\S*$//;
	
	$flux_function_name = $flux_file_name;
	$flux_function_name =~ s/\.\S*$//;
	
	$j_function_name = $j_file_name;
	$j_function_name =~ s/\.\S*$//;
}

### MCMC stuff ###
if($l_coll_factor_mcmc_all){
	$l_coll_factor_mcmc_n_n = 1;
	$l_coll_factor_mcmc_n_i = 1;
	$l_coll_factor_mcmc_rec = 1;
}
if($l_coll_factor_mcmc_n_i){
	$l_coll_factor_mcmc_n_g = 1;
}

if($l_coll_factor_mcmc_n_n || $l_coll_factor_mcmc_n_i || $l_coll_factor_mcmc_rec || $l_coll_factor_mcmc_n_g){
	$l_coll_factor_mcmc_any = 1;
}
if($l_coll_factor_mcmc_any || $l_evap_factor_mcmc || $l_q_neg_mcmc || $l_q_pos_mcmc || $l_wl_mcmc || $#wl_mcmc_only>=0){
	$l_mcmc_any = 1;
	$ncoef_coll = 0;
	$ncoef_coll_n_n = 0;
	$ncoef_coll_n_i = 0;
	$ncoef_coll_n_g = 0;
	$ncoef_coll_rec = 0;
	$ncoef_evap = 0;
	$ncoef_wl = 0;
	if($lfortran == 0){
		die "\nThe MCMC coefficients work only with the fortran version!\n\n";
	}
	if($variable_temp){
		die "\nThe MCMC coefficients cannot be combined with a varying temperature!\n\n";
	}
	if(($sticking_factor[0][0] != 1) || ($sticking_factor[1][0] != 1) || ($sticking_factor_file_name ne '')){
		die "\nThe MCMC coefficients cannot be combined with other sticking factors!\n\n";
	}
}

### Parameters with non-fixed values ###
if($variable_cs){
	if(!$lfortran || $lloop){
		die "\nThe --variable_cs option works only with the fortran version and the non-loop format.\n\n";
	}
}
if($variable_ion_source ){
	if($lfortran == 0){
		die "\nThe --variable_ion_source option works only with the fortran version.\n\n";
	}
	if($no_generic_neg || $no_generic_pos){
		die "\nThe --variable_ion_source option works only if both charger ions are enabled.\n\n";
	}
}

### Constants and properties of generic compounds (charger ions, water...) ###

# Values for constants
$boltz=1.3806504e-23; # J/K
$pi=3.14159265359;
$Na=6.02214179e23;
$mass_conv=1.0/$Na/1000.0; # convert g/mol to kg
$kcal_per_mol_to_J=4.184*1000.0/$Na;

# Properties of charger ions
if($lnitrate){
	# nitrate ion dimer
	$mass_nion=125*$mass_conv; # kg
	$dens_nion=1512.9; # liquid density in kg/m**3
}else{
	# the generic negative ion now has the mass and density of an O2 molecule
	$mass_nion=32.00*$mass_conv; # kg
	$dens_nion=1141.0; # liquid density in kg/m**3
}
# the generic positive ion now has the mass and density of a (charged) water molecule
$mass_pion=19.02*$mass_conv; # kg
$dens_pion=997.0; # kg/m**3

# Use equilibrium hydrate distributions?
if($rh > 0){
	$lhydrate_average=1;	
	if($variable_temp==1){
		die "\nThe option --variable_temp is not compatible with hydrate distributions!\n\n";
	}
}else{
	$lhydrate_average=0;
}

# Determine if evaporations are disabled, and read in the temperature (use it also for the hydrate distribution calculations)
if(($disable_evap && !$lhydrate_average) || $lloop){
	if(!$temperature){
		die "\nYou must give the temperature (for the collision coefficients) if evaporations are disabled.\n\n";
	}
}elsif( -f $free_energy_file_name){
	
	# read in the temperature
	open(GFILE, "<$free_energy_file_name") || die "\nCould not open the delta G file $free_energy_file_name\n\n";
	chomp(@evap_data=<GFILE>);
	# the first line has the pressure in Pa, and the second one has the temperature
	if(&is_pos_number($evap_data[1])){
		$temperature_inp=$evap_data[1];
		if ($temperature && ($temperature != $temperature_inp)){
			print "Converting thermodynamics data from temperature $temperature_inp K to $temperature K.\n";
			$l_change_temperature = 1;
		}else{
			$temperature = $temperature_inp;
		}
	}else{
		die "\nSecond line of $free_energy_file_name has to be the reference temperature, in K.\n\n";
	}
	close(GFILE);
}

# Properties of water for the hydrate distribution calculations
$name_water='W';
$mass_water=18.02*$mass_conv; # kg
$dens_water=997.0; # kg/m**3
$volume_water=$mass_water/$dens_water;
# Saturation vapor pressure from Arnold Wexler: 'Vapor Pressure Formulation for Water in Range 0 to 100 °C. A Revision'. JOURNAL OF RESEARCH of the National Bureau of Standards, 80A, 5-6, 775-785, 1976. 
$psw = exp(-2991.2729*$temperature**(-2)-6017.0128/$temperature+18.87643854-0.028354721*$temperature+0.17838301e-4*$temperature**2-0.84150417e-9*$temperature**3+0.44412543e-12*$temperature**4+2.858487*log($temperature));
# Previously used formula from Hanna Vehkamäki's CNT book
#$psw = exp(77.344913-7235.4247/$temperature-8.2*log($temperature)+0.0057113*$temperature);


################################################################################
################################################################################
##   2. Reading the input file                                                ##
################################################################################
################################################################################


######################################################################
### This reads in the input file, which has a listing of all the
### different cluster compositions, using a compact notation...for example,
### s	b   a   d
### 0-5  0   1   0
### is a series of six clusters with one ammonia, 0 ions, 0 dimethylamines, and 0-5 SA molecules.
### This neatly divides the output file into rows, which makes it more readable.
### A row is a grouping of clusters where the number of ions, ammonias, and DMAs
### are all constant.
### Any lines beginning with # are ignored.
######################################################################



$max_cluster=0;  #this is the number of the cluster
# maximum amounts of different kinds of molecules (for any cluster), separately for neutral, neg. and pos. clusters 
@n_max_type = (0);
@n_max_type_pos = (0);
@n_max_type_neg = (0);
@clust_comp = ([0]);
@temp_array = ();
# molecule numbers of neutral molecules
@neutral_mols = ();
$n_neutral_mols = 0;
# logical telling if acid and/or base strengths are defined for neutral molecules
$l_strength_neutral = 0;

# Arrays to save nucleation criteria and their sizes
$n_nucl_rules = 0;
$n_nucl_rules_neg = 0;
$n_nucl_rules_pos = 0;

open(INP, "<$input_file_name") || die "\nCould not open input file $input_file_name\n\n";
chomp(@input_file=<INP>); # find out how many molecules we have
close(INP);

$number_of_lines=$#input_file;
@columns=split(' ',$input_file[0]); # first line should be e.g. '#   sulphuric_acid   bisulphate ...', basically the names of the molecules 
if($columns[0] eq '#'){
	shift(@columns);
}
$nmol_types=$#columns+1;

for($imol=0;$imol<$nmol_types;$imol++){
	push(@n_max_type,0);
	push(@n_max_type_neg,0);
	push(@n_max_type_pos,0);
	$l_can_evaporate[$imol] = 1;
	$l_can_be_lost_bound[$imol] = 1;
}
$found_negative=0; # changes to true if we find a cluster with a negative ion in it
$found_positive=0; # changes to true if we find a cluster with a positive ion in it

# Arrays for the indices of negative and positive clusters
@negative_clusters=();
@positive_clusters=();

# This is the molecule number for proton (positive hydrogen ion)
$proton=-1;
# This is the molecule number for missing_proton (removed from an acid to charge it negatively)
$missing_proton=-1;
$nmolecules = 0;

# Keep track of the highest acid and base strength, so we can give reasonable values for the proton and missing proton.
$acid_max = 0;
$base_max = 0;

for($iline=1;$iline <= $number_of_lines;$iline++){
	$line=$input_file[$iline];
	next if($line =~ /^\s*$/);
	@columns=split(' ',$line);
	next if($columns[0] =~ /^\#/);

################################################################################
################################################################################
##      2.1 Reading the info on the molecule types                            ##
################################################################################
################################################################################

	if($#columns < $nmol_types-1){
		die("Can't understand line ".($iline+1)." ($line) of input file - too few columns.\n");
	}else{
		
		# the name labels of all the molecules
		if($input_file[$iline] =~ /name/i){

			$nend=$#columns;
			$iend=1;
			while($iend <= $nmol_types){
				$molecule_name[$nmol_types-$iend]=$columns[$nend+1-$iend];
				$mol_num_from_name{$columns[$nend+1-$iend]}=$nmol_types-$iend;
				$iend++;
			}

		# the charges of all the molecules
		}elsif($input_file[$iline] =~ /charge/i){

			$nend=$#columns;
			$iend=1;
			while($iend <= $nmol_types){
				$charge[$nmol_types-$iend]=$columns[$nend+1-$iend];
				$iend++;
			}
			for($imol=0;$imol<$nmol_types;$imol++){
				die "\nCharge should be 0, -1 or 1!\n\n" unless($charge[$imol] == -1 || $charge[$imol] == 0 || $charge[$imol] == 1);
				
				# If this is the proton, let's save the molecule number for it
				if($charge[$imol] == 1){
					$proton=$imol;
				}elsif($charge[$imol] == 0){
					push(@neutral_mols,$imol);
					$n_neutral_mols++;
				}
			}
		
		# the corresponding neutral molecules for the ions
		}elsif($input_file[$iline] =~ /corresponding neutral molecule/i){

			$nend=$#columns;
			$iend=1;
			while($iend <= $nmol_types){
				$imol = -1;
				if($columns[$nend+1-$iend] eq '-'){
					$corr_neutral_mol[$nmol_types-$iend]='';
				}else{
					$corr_neutral_mol[$nmol_types-$iend]=$columns[$nend+1-$iend];
					if(exists $mol_num_from_name{$corr_neutral_mol[$nmol_types-$iend]}){
						$imol=$mol_num_from_name{$corr_neutral_mol[$nmol_types-$iend]};
					}
				}
				$n_corr_neutral_mol[$nmol_types-$iend] = $imol;
				$iend++;
			}


		# the corresponding positive ions for bases
		}elsif($input_file[$iline] =~ /corresponding positive ion/i){

			$nend=$#columns;
			$iend=1;
			while($iend <= $nmol_types){
				$imol = -1;
				if($columns[$nend+1-$iend] eq '-'){
					$corr_positive_mol[$nmol_types-$iend]='';
				}else{
					$corr_positive_mol[$nmol_types-$iend]=$columns[$nend+1-$iend];
					if(exists $mol_num_from_name{$corr_positive_mol[$nmol_types-$iend]}){
						$imol=$mol_num_from_name{$corr_positive_mol[$nmol_types-$iend]};
					}
				}
				$n_corr_positive_mol[$nmol_types-$iend] = $imol;
				$iend++;
			}

			# If we want to do the table output, we need to know if the first molecule has a positive form.
			if($ldo_table_out){ 
				$first_pos = -1;
				if(length($corr_positive_mol[0])>0){
					if(exists $mol_num_from_name{$corr_positive_mol[0]}){
						$first_pos = $mol_num_from_name{$corr_positive_mol[0]};
					}
				}
			}

		# the corresponding negative ions for acids
		}elsif($input_file[$iline] =~ /corresponding negative ion/i){

			$nend=$#columns;
			$iend=1;
			while($iend <= $nmol_types){
				$imol = -1;
				if($columns[$nend+1-$iend] eq '-'){
					$corr_negative_mol[$nmol_types-$iend]='';
				}else{
					$corr_negative_mol[$nmol_types-$iend]=$columns[$nend+1-$iend];
					if(exists $mol_num_from_name{$corr_negative_mol[$nmol_types-$iend]}){
						$imol=$mol_num_from_name{$corr_negative_mol[$nmol_types-$iend]};
					}
				}
				$n_corr_negative_mol[$nmol_types-$iend] = $imol;
				$iend++;
			}

			# If we want to do the table output, we need to know if the first molecule has a negative form.
			if($ldo_table_out){ 
				$first_neg = -1;
				if(length($corr_negative_mol[0])>0){
					if(exists $mol_num_from_name{$corr_negative_mol[0]}){
						$first_neg = $mol_num_from_name{$corr_negative_mol[0]};
					}
				}
			}
		
		# the masses of all the molecules
		}elsif($input_file[$iline] =~ /mass \[\S+\]/i){

			$nend=$#columns;
			$iend=1;
			while($iend <= $nmol_types){
				$mass[$nmol_types-$iend]=$columns[$nend+1-$iend]*$mass_conv; # convert to kg
				if($mass[$nmol_types-$iend] < 0){
					$missing_proton=$nmol_types-$iend;
				}
				$iend++;
			}

		# the densities of all the molecules
		}elsif($input_file[$iline] =~ /density \[\S+\]/i){

			$nend=$#columns;
			$iend=1;
			while($iend <= $nmol_types){
				$density[$nmol_types-$iend]=$columns[$nend+1-$iend];
				$iend++;
			}
			
		# the saturation vapor pressures of all the pure compounds
		# (optional)
		}elsif($input_file[$iline] =~ /saturation vapor pressure \[\S+\]/i){

			$nend=$#columns;
			$iend=1;
			while($iend <= $nmol_types){
				$psat[$nmol_types-$iend]=$columns[$nend+1-$iend];
				$iend++;
			}
			
		# the surface tensions of all the pure compounds
		# (optional)
		}elsif($input_file[$iline] =~ /surface tension \[\S+\]/i){

			$nend=$#columns;
			$iend=1;
			while($iend <= $nmol_types){
				$surface_tension[$nmol_types-$iend]=$columns[$nend+1-$iend];
				$iend++;
			}

		# the acidities of all the molecules
		# (-1 = base, 0 = neither acid or base, 1 = weak acid, 2 = stronger acid, 3 = even stronger acid, ...)
		}elsif($input_file[$iline] =~ /acid/i){

			$nend=$#columns;
			$iend=1;
			while($iend <= $nmol_types){
				$acid_strength[$nmol_types-$iend]=$columns[$nend+1-$iend];
				$acid_max = max($acid_max,$acid_strength[$nmol_types-$iend]);
				$iend++;
			}

		# the basicities of all the molecules
		# (-1 = acid, 0 = neither acid or base, 1 = weak base, 2 = stronger base, 3 = even stronger base, ...)
		}elsif($input_file[$iline] =~ /base/i){

			$nend=$#columns;
			$iend=1;
			while($iend <= $nmol_types){
				$base_strength[$nmol_types-$iend]=$columns[$nend+1-$iend];
				$base_max = max($base_max,$base_strength[$nmol_types-$iend]);
				$iend++;
			}
		
		
		# can this molecule evaporate ...
		}elsif($input_file[$iline] =~ /can be lost/i || $input_file[$iline] =~ /can evaporate/i){
			
			# ... at all ...
			if($input_file[$iline] !~ /bound/i){
				$nend=$#columns;
				$iend=1;
				while($iend <= $nmol_types){
					$l_can_evaporate[$nmol_types-$iend]=$columns[$nend+1-$iend];
					$iend++;
				}
			}else{
			# ... and when bringing clusters back to the boundary
				$nend=$#columns;
				$iend=1;
				while($iend <= $nmol_types){
					$l_can_be_lost_bound[$nmol_types-$iend]=$columns[$nend+1-$iend];
					$iend++;
				}
			}

################################################################################
################################################################################
##      2.2 Reading the nucleation rules                                      ##
################################################################################
################################################################################

		# the "nucleation criteria" for neutrals
		}elsif($input_file[$iline] =~ /out neutral/i){

			$nend=$#columns;
			$iend=1;
			while($iend <= $nmol_types){
				$nucl_rules[$n_nucl_rules][$nmol_types-$iend]=$columns[$nend+1-$iend];
				$iend++;
			}
			$n_nucl_rules++;

		# the "nucleation criteria" for negatives
		}elsif($input_file[$iline] =~ /out negative/i){

			$nend=$#columns;
			$iend=1;
			while($iend <= $nmol_types){
				$nucl_rules_neg[$n_nucl_rules_neg][$nmol_types-$iend]=$columns[$nend+1-$iend];
				$iend++;
			}
			$n_nucl_rules_neg++;

		# the "nucleation criteria" for positives
		}elsif($input_file[$iline] =~ /out positive/i){

			$nend=$#columns;
			$iend=1;
			while($iend <= $nmol_types){
				$nucl_rules_pos[$n_nucl_rules_pos][$nmol_types-$iend]=$columns[$nend+1-$iend];
				$iend++;
			}
			# keeping track of the maximum number that the 1st molecule can have in negative clusters of the system
			$n_nucl_rules_pos++;

################################################################################
################################################################################
##      2.3 Reading the coagulation sink cluster                              ##
################################################################################
################################################################################

		}elsif($input_file[$iline] =~ /coag/i && $input_file[$iline] =~ /neg/i){
			$use_coag_to_outgoing = 1;
			$nend=$#columns;
			$iend=1;
			if(defined $cs_cluster_neg[0]){
				die "You can only have one negative coagulation sink cluster.\n";
			}
			@cs_cluster_neg = ();
			while($iend <= $nmol_types){
				$cs_cluster_neg[$nmol_types-$iend]=$columns[$nend+1-$iend];
				$iend++;
			}
		}elsif($input_file[$iline] =~ /coag/i && $input_file[$iline] =~ /pos/i){
			$use_coag_to_outgoing = 1;
			$nend=$#columns;
			$iend=1;
			if(defined $cs_cluster_pos[0]){
				die "You can only have one positive coagulation sink cluster.\n";
			}
			@cs_cluster_pos = ();
			while($iend <= $nmol_types){
				$cs_cluster_pos[$nmol_types-$iend]=$columns[$nend+1-$iend];
				$iend++;
			}
		}elsif($input_file[$iline] =~ /coag/i && $input_file[$iline] =~ /neutr/i){
			$use_coag_to_outgoing = 1;
			$nend=$#columns;
			$iend=1;
			if(defined $cs_cluster[0]){
				die "You can only have one neutral coagulation sink cluster.\n";
			}
			@cs_cluster = ();
			while($iend <= $nmol_types){
				$cs_cluster[$nmol_types-$iend]=$columns[$nend+1-$iend];
				$iend++;
			}
		}elsif($input_file[$iline] =~ /coag/i){
			$use_coag_to_outgoing = 1;
			$nend=$#columns;
			$iend=1;
			@temp_array = ();
			$i = 0;
			while($iend <= $nmol_types){
				$temp_array[$nmol_types-$iend]=$columns[$nend+1-$iend];
				if($temp_array[$nmol_types-$iend]>0 && $charge[$nmol_types-$iend] < 0){
					$i -= $temp_array[$nmol_types-$iend];
				}elsif($temp_array[$nmol_types-$iend]>0 && $charge[$nmol_types-$iend] > 0){
					$i += $temp_array[$nmol_types-$iend];
				}
				$iend++;
			}
			if($i==0){
				if(defined $cs_cluster[0]){
					die "You can only have one neutral coagulation sink cluster.\n";
				}
				@cs_cluster = @temp_array;
			}elsif($i<0){
				if(defined $cs_cluster_neg[0]){
					die "You can only have one negative coagulation sink cluster.\n";
				}
				@cs_cluster_neg = @temp_array;
			}else{
				if(defined $cs_cluster_pos[0]){
					die "You can only have one positive coagulation sink cluster.\n";
				}
				@cs_cluster_pos = @temp_array;
			}


	
################################################################################
################################################################################
##      2.4 Reading the cluster types                                         ##
################################################################################
################################################################################
		
		}elsif($#columns == $nmol_types-1 || $input_file[$iline] =~ /no evap/i){
		
			# Removing the columns not related to molecule numbers from the lines starting with "no evap"
			while($nmol_types-1<$#columns){
				shift(@columns);
			}
		
			# Checking which column (if any) contains a range of molecule numbers
			for($icol=0;$icol<$nmol_types;$icol++){
				$lfound = 0;
				if($columns[$icol] =~ /(\d+)\-(\d+)/){
					$lfound = 1;
					last;
				}
			}
			if(!$lfound){
#				print "\nIf you only want to specify one cluster by the line\n'$input_file[$iline]',\nyou should use $columns[0]-$columns[0] for the first molecule type. Trying to proceed anyway, but fix your input file to be safe.\n\n";
				$columns[0] = "$columns[0]-$columns[0]";
				$icol = 0;
			}
		
			# Creating cluster labels and recording their compositions
			$columns[$icol] =~ /(\d+)\-(\d+)/;
			for($jclus=$1;$jclus<=$2;$jclus++){
				$max_cluster++;
				$columns[$icol]=$jclus;
				$nmol = 0;
				for($imol=0;$imol<$nmol_types;$imol++){
					$nmol += $columns[$imol];
				}
				if($nmol < 1){
					die "\nYou have a cluster with no molecules in the input file!\n\n";
				}
				push(@clust_comp,[@columns]);
			
				# Recording whether this cluster can evaporate
				if($input_file[$iline] =~ /no evap/i){
					$disable_evap_only[$max_cluster]=1;
				}else{
					$disable_evap_only[$max_cluster]=0;
				}
			
				# these logicals tell if the cluster is negative, positive or neutral (1 is true, 0 is false)
				$lnegative_cluster[$max_cluster]=0;
				$lnegative[$max_cluster]=0;
				$lpositive_cluster[$max_cluster]=0;
				$lpositive[$max_cluster]=0;
				$lneutral_cluster[$max_cluster]=1;
				$lneutral[$max_cluster]=1;
				for($imol=0;$imol<$nmol_types;$imol++){
					if($clust_comp[$max_cluster][$imol] > 0){
						if($charge[$imol] < 0){
							$found_negative =1;
							$lnegative_cluster[$max_cluster]=1;
							$lnegative[$max_cluster]=1;
							$lneutral_cluster[$max_cluster]=0;
							$lneutral[$max_cluster]=0;
							push(@negative_clusters,$max_cluster);
							if($clust_comp[$max_cluster][$imol] > 1){
								die "\nCan't have more than one charged molecule per cluster!\n\n";
							}
						}elsif($charge[$imol] > 0){
							$found_positive =1;
							$lpositive_cluster[$max_cluster]=1;
							$lpositive[$max_cluster]=1;
							$lneutral_cluster[$max_cluster]=0;
							$lneutral[$max_cluster]=0;
							push(@positive_clusters,$max_cluster);
							if($clust_comp[$max_cluster][$imol] > 1){
								die "\nCan't have more than one charged molecule per cluster!\n\n";
							}
						}
					}
				}
				if($lnegative_cluster[$max_cluster] && $lpositive_cluster[$max_cluster]){
					die "\nCan't have both negative and positive molecules in the same cluster!\n\n";
				}
				if($lnegative_cluster[$max_cluster]){
					for($imol=0;$imol<$nmol_types;$imol++){
							$nmol = $clust_comp[$max_cluster][$imol];
							if($n_max_type_neg[$imol] < $nmol){
								$n_max_type_neg[$imol] = $nmol;
							}
#							if($charge[$imol] == 0 && $n_corr_negative_mol[$imol] >= 0){
#								$nmol += $clust_comp[$max_cluster][$n_corr_negative_mol[$imol]];
#							}
#							if($charge[$imol] < 0 && $n_corr_neutral_mol[$imol] >= 0){
#								$nmol -= $clust_comp[$max_cluster][$n_corr_neutral_mol[$imol]];
#								$nmol = max($nmol,0);
#							}
#							if($n_max_type_neg_0[$imol] < $nmol){
#								$n_max_type_neg_0[$imol] = $nmol;
#							}
					}
				}elsif($lpositive_cluster[$max_cluster]){
					for($imol=0;$imol<$nmol_types;$imol++){
							$nmol = $clust_comp[$max_cluster][$imol];
							if($n_max_type_pos[$imol] < $nmol){
								$n_max_type_pos[$imol] = $nmol;
							}
#							if($charge[$imol] == 0 && $n_corr_positive_mol[$imol] >= 0){
#								$nmol = $nmol+$clust_comp[$max_cluster][$n_corr_positive_mol[$imol]];
#							}
#							if($charge[$imol] > 0 && $n_corr_neutral_mol[$imol] >= 0){
#								$nmol -= $clust_comp[$max_cluster][$n_corr_neutral_mol[$imol]];
#								$nmol = max($nmol,0);
#							}
#							if($n_max_type_pos[$imol] < $nmol){
#								$n_max_type_pos[$imol] = $nmol;
#							}
					}
				}else{
					for($imol=0;$imol<$nmol_types;$imol++){
							if($n_max_type[$imol] < $clust_comp[$max_cluster][$imol]){
								$n_max_type[$imol] = $clust_comp[$max_cluster][$imol];
							}
					}
				}
			
				# Make the cluster label
				$label[$max_cluster]=&create_label(\@columns);
				if($l_print_boundary){
					print "Cluster $max_cluster: $label[$max_cluster]\n";
				}
				$comments[$max_cluster]="% Cluster $max_cluster: $label[$max_cluster]";
				$num_from_label{"$label[$max_cluster]"} = $max_cluster;

				# check if this is a monomer
				$lmonomer[$max_cluster]=&check_monomer(\@columns);
				if ($lmonomer[$max_cluster]){
					$nmolecules += 1;
					$molecules[$nmolecules] = $max_cluster;
				}
			
				# get the mass and radius of the cluster
				$clust_mass[$max_cluster][0] = &calculate_mass(\@columns);
				if($l_bulk_density){
					if($lhydrate_average){
						($clust_volume[$max_cluster],$clust_radius[$max_cluster][0]) = &calculate_volume_and_radius(\@columns);
					}else{
						$clust_radius[$max_cluster][0] = &calculate_radius(\@columns);
					}
				}
				#print "$label[$max_cluster]: mass = $clust_mass[$max_cluster][0], radius =  $clust_radius[$max_cluster][0]\n";
			}
		}else{
			print "Cannot understand this line of the input file. Line number: ".($iline+1)."\nData on line: $input_file[$iline]\ncolumns:$#columns\nnmol_types-1:".($nmol_types-1)."\n";
			die;
		}
	}
}
## Done reading the input file

if($lloop && $max_cluster>1){
	die "Too many clusters in the input file! In loop mode, you should only give the composition of the biggest cluster.\n";
}
if($lloop && ($n_nucl_rules>0 || $n_nucl_rules_neg>0 || $n_nucl_rules_pos>0)){
	die "You should not define nucleation rules in loop mode - all outgoing collisions are allowed.\n";
}

# Adding a cluster number for the formed particles if they act as a coag. sink
if($use_coag_to_outgoing){
	if($n_nucl_rules>0){
		if(!defined($cs_cluster[0])){
			die "Define a coagulation sink cluster also for neutrals, or none at all.";
		}
		$max_cluster++;
		$n_flux_out=$max_cluster;
		$n_flux_out_max=$n_flux_out;
		push(@clust_comp,[@cs_cluster]);
		# Make the cluster label
		$cs_cluster_label=&create_label(\@cs_cluster);
		$label[$max_cluster]='cs_cluster';
		if($l_print_boundary){
			print "Cluster $max_cluster: $label[$max_cluster] (=formed neutral particles)\n";
		}
		$comments[$max_cluster]="% Cluster $max_cluster: $label[$max_cluster] (=formed particles)";
		$lnegative_cluster[$max_cluster]=0;
		$lnegative[$max_cluster]=0;
		$lpositive_cluster[$max_cluster]=0;
		$lpositive[$max_cluster]=0;
		$lneutral_cluster[$max_cluster]=1;
		$lneutral[$max_cluster]=1;
		# get the mass and radius of the cluster
		$clust_mass[$max_cluster][0] = &calculate_mass(\@cs_cluster);
		if($l_bulk_density){
			$clust_radius[$max_cluster][0] = &calculate_radius(\@cs_cluster);
		}
	}
	if($n_nucl_rules_neg>0){
		if(!defined($cs_cluster_neg[0])){
			die "Define a coagulation sink cluster also for negatives, or none at all.";
		}
		$max_cluster++;
		$n_flux_out_neg=$max_cluster;
		$n_flux_out_max=$n_flux_out_neg;
		push(@clust_comp,[@cs_cluster_neg]);
		# Make the cluster label
		$cs_cluster_neg_label=&create_label(\@cs_cluster_neg);
		$label[$max_cluster]='cs_cluster_neg';
		if($l_print_boundary){
			print "Cluster $max_cluster: $label[$max_cluster] (=formed negative particles)\n";
		}
		$comments[$max_cluster]="% Cluster $max_cluster: $label[$max_cluster] (=formed negative particles)";
		$lnegative_cluster[$max_cluster]=1;
		$lnegative[$max_cluster]=1;
		$lpositive_cluster[$max_cluster]=0;
		$lpositive[$max_cluster]=0;
		$lneutral_cluster[$max_cluster]=0;
		$lneutral[$max_cluster]=0;
		# get the mass and radius of the cluster
		$clust_mass[$max_cluster][0] = &calculate_mass(\@cs_cluster_neg);
		if($l_bulk_density){
			$clust_radius[$max_cluster][0] = &calculate_radius(\@cs_cluster_neg);
		}
	}
	if($n_nucl_rules_pos>0){
		if(!defined($cs_cluster_pos[0])){
			die "Define a coagulation sink cluster also for positives, or none at all.";
		}
		$max_cluster++;
		$n_flux_out_pos=$max_cluster;
		$n_flux_out_max=$n_flux_out_pos;
		push(@clust_comp,[@cs_cluster_pos]);
		# Make the cluster label
		$cs_cluster_pos_label=&create_label(\@cs_cluster_pos);
		$label[$max_cluster]='cs_cluster_pos';
		if($l_print_boundary){
			print "Cluster $max_cluster: $label[$max_cluster] (=formed positive particles)\n";
		}
		$comments[$max_cluster]="% Cluster $max_cluster: $label[$max_cluster] (=formed positive particles)";
		$lnegative_cluster[$max_cluster]=0;
		$lnegative[$max_cluster]=0;
		$lpositive_cluster[$max_cluster]=1;
		$lpositive[$max_cluster]=1;
		$lneutral_cluster[$max_cluster]=0;
		$lneutral[$max_cluster]=0;
		# get the mass and radius of the cluster
		$clust_mass[$max_cluster][0] = &calculate_mass(\@cs_cluster_pos);
		if($l_bulk_density){
			$clust_radius[$max_cluster][0] = &calculate_radius(\@cs_cluster_pos);
		}
	}
}


for($imol=0;$imol<$nmol_types;$imol++){
	# Array for cluster numbers corresponding to monomers
	$clust_num_of_monomer[$imol] = &get_cluster_number("1$molecule_name[$imol]");
	# If the acid and base strengths were not given, give the same value for all molecules
	if(!(defined $acid_strength[$imol])){
		$acid_strength[$imol] = 0;
	}
	if(!(defined $base_strength[$imol])){
		$base_strength[$imol] = 0;
	}
	# See if any acid and/or base strengths are given for neutral molecules
	if($charge[$imol] == 0 && ($acid_strength[$imol] > 0 || $base_strength[$imol] > 0)){
		$l_strength_neutral = 1;
	}
}

################################################################################
################################################################################
##      2.5 Adding molecule types for the proton and removed proton if needed ##
################################################################################
################################################################################

$include_negative_ion_terms = 0;
$include_positive_ion_terms = 0;
$include_ion_terms = 0;

if($found_negative){
	$include_negative_ion_terms = 1;
	$include_ion_terms = 1;
}
if($found_positive){
	$include_positive_ion_terms = 1;
	$include_ion_terms = 1;
}

if(!$include_ion_terms){
	# if we don't have any ions, no need to include generic charger ions
	$include_generic_ions = 0;
}

if($include_generic_ions){
	if($no_generic_neg){
		$include_generic_neg = 0;
	}else{
		$include_generic_neg = 1;
	}
	if($no_generic_pos){
		$include_generic_pos = 0;
	}else{
		$include_generic_pos = 1;
	}
}

# if we have charger ions or several molecule types that can get charged to the same polarity, 
# we need molecule types for the proton and/or removed proton

# first checking the proton ...
if($proton==-1){
	$label_proton = '';
	$lfound = -1;
	for($imol=0;$imol<$nmol_types;$imol++){
		if($charge[$imol] == 1 && $n_corr_neutral_mol[$imol] >= 0){
			$lfound += 1;
			if($label_proton eq ''){
				$label_proton = "$molecule_name[$imol]-1$molecule_name[$n_corr_neutral_mol[$imol]]";
			}else{
				print "Warning: the positive charge might not be able to move from one molecule type to another.\n";
			}
		}
	}
	if($include_generic_pos){
		$lfound = 1;
	}
	if($lfound>=1 && $label_proton eq ''){
		$label_proton = 'P';
		$proton = $nmol_types;
		$nmol_types += 1;
		$molecule_name[$proton] = 'P';
		$mol_num_from_name{'P'}=$proton;
		$charge[$proton] = 1;
		$corr_neutral_mol[$proton] = '';
		$corr_negative_mol[$proton] = '';
		$corr_positive_mol[$proton] = '';
		$n_corr_neutral_mol[$proton] = -1;
		$n_corr_negative_mol[$proton] = -1;
		$n_corr_positive_mol[$proton] = -1;
		$mass[$proton] = 1.*$mass_conv;
		$density[$proton] = '';
		$acid_max += 1;
		$acid_strength[$proton] = $acid_max;
		# Add the new molecule type also in the cluster composition table ...
		for($i=1;$i<=$max_cluster;$i++){
			$clust_comp[$i][$proton] = 0;
		}
		# ... and in the nucleation rules		
		for($i=0;$i<$n_nucl_rules;$i++){
			$nucl_rules[$i][$proton] = 0;
		}
		for($i=0;$i<$n_nucl_rules_neg;$i++){
			$nucl_rules_neg[$i][$proton] = 0;
		}
		for($i=0;$i<$n_nucl_rules_pos;$i++){
			$nucl_rules_pos[$i][$proton] = 0;
		}
	}
}else{
	$label_proton = $molecule_name[$proton];
}

# ... then the missing proton
if($missing_proton==-1){
	$label_missing_proton = '';
	$lfound = -1;
	for($imol=0;$imol<$nmol_types;$imol++){
		if($charge[$imol] == -1 && $n_corr_neutral_mol[$imol] >= 0){
			$lfound += 1;
			if($label_missing_proton eq ''){
				$label_missing_proton = "$molecule_name[$imol]-1$molecule_name[$n_corr_neutral_mol[$imol]]";
			}else{
				print "Warning: the negative charge might not be able to move from one molecule type to another.\n";
			}
		}
	}
	if($include_generic_neg){
		$lfound = 1;
	}
	if($lfound>=1 && $label_missing_proton eq ''){
		$label_missing_proton = 'RP';
		$missing_proton = $nmol_types;
		$nmol_types += 1;
		$molecule_name[$missing_proton] = 'RP';
		$mol_num_from_name{'RP'}=$missing_proton;
		$charge[$missing_proton] = -1;
		$corr_neutral_mol[$missing_proton] = '';
		$corr_negative_mol[$missing_proton] = '';
		$corr_positive_mol[$missing_proton] = '';
		$n_corr_neutral_mol[$missing_proton] = -1;
		$n_corr_negative_mol[$missing_proton] = -1;
		$n_corr_positive_mol[$missing_proton] = -1;
		$mass[$missing_proton] = -1.*$mass_conv;
		$density[$missing_proton] = '';	
		$base_max += 1;
		$base_strength[$missing_proton] = $base_max;
		# Add the new molecule type also in the cluster composition table ...
		for($i=1;$i<=$max_cluster;$i++){
			$clust_comp[$i][$missing_proton] = 0;
		}
		# ... and in the nucleation rules		
		for($i=0;$i<$n_nucl_rules;$i++){
			$nucl_rules[$i][$missing_proton] = 0;
		}
		for($i=0;$i<$n_nucl_rules_neg;$i++){
			$nucl_rules_neg[$i][$missing_proton] = 0;
		}
		for($i=0;$i<$n_nucl_rules_pos;$i++){
			$nucl_rules_pos[$i][$missing_proton] = 0;
		}	
	}
}else{
	$label_missing_proton = $molecule_name[$missing_proton];
}

################################################################################
################################################################################
##      2.6 Determining the system size for the loop option                   ##
################################################################################
################################################################################
if($lloop){
	if(!$lloop_max_ratio){
		$max_cluster = 1;
		for($imol=0;$imol<$nmol_types;$imol++){
			if($n_max_type[$imol]>0){
				$max_cluster *= ($n_max_type[$imol]+1); # the "+1" corresponds to zero molecules
			}
		}
		$max_cluster -= 1; # the total number of all molecules can't be zero
	}else{
		# Find the number of the molecule type max_ratio_basis, and the max. number of other molecules with respect to a monomer
		$lfound = 0;
		for($imol=0;$imol<$nmol_types;$imol++){
			if($molecule_name[$imol] eq $max_ratio_basis){
				$lfound = 1;
				$nmol_inp_max_ratio_basis = $imol;
				last;
			}
		}
		die "Did not find the given basis molecule from the input file\n" if(!$lfound);
		$lfound = 0;
		for($imol=0;$imol<$nmol_types;$imol++){
			if($n_max_type[$imol]>0 && $imol != $nmol_inp_max_ratio_basis){
				die "Too many molecule types in the input file; you should have two types for the max_ratio_basis option\n" if($lfound);
				$lfound = 1;
				$nmol_inp_not_basis = $imol;
				$max_wrt_mon = $n_max_type[$nmol_inp_not_basis]/$n_max_type[$nmol_inp_max_ratio_basis];
				print "\nUsing a max. $molecule_name[$nmol_inp_not_basis]:$molecule_name[$nmol_inp_max_ratio_basis] ratio of $max_wrt_mon\n";
				if($max_wrt_mon !~ /^\d+$/){
					die "Max. number of molecules $molecule_name[$nmol_inp_not_basis] per one molecule $molecule_name[$nmol_inp_max_ratio_basis] is not an integer\n";
				}
			}
		}
		die "Did not find other molecule types in addition to the given basis molecule from the input file\n" if(!$lfound);
		$max_cluster = 0;
		for($i=0;$i<=$n_max_type[$nmol_inp_max_ratio_basis];$i++){ # $i=0 corresponds to the monomer of the non-basis type
			$max_cluster += $i*$max_wrt_mon+1; # the "+1" corresponds to the cluster consisting solely of the basis type molecules when $i > 0
		}
	}
	# Determine the maximum numbers for fluxes and clusters
	if($l_save_outgoing){
		$max_n_flux = $max_cluster+1;
	}else{
		$max_n_flux = $max_cluster;
	}
	$max_cluster_number = $max_cluster;
	# Determine the index label for the outgoing flux
	if($lloop_cs){
		$n_flux_out = 'n_cs';
	}else{
		$n_flux_out = 'nclust+1';
	}
}else{
################################################################################
################################################################################
##      2.7 Adding cluster numbers for generic ions if they are used          ##
################################################################################
################################################################################

	# $max_cluster is the biggest number for a proper cluster (possibly including the formed particle)
	# $max_cluster_number is the biggest number for a cluster or generic ion
	$max_cluster_number=$max_cluster;
	$n_generic_neg=-1;
	$n_generic_pos=-1;

	if($include_generic_neg){
		$n_generic_neg=$max_cluster_number+1;
		$lnegative[$n_generic_neg] = 1;
		$lpositive[$n_generic_neg] = 0;
		$lneutral_cluster[$n_generic_neg] = 0;
		$label[$n_generic_neg]=$label_generic_neg;
		$comments[$n_generic_neg]="% Cluster $n_generic_neg: generic negative ion ";
		$max_cluster_number=$n_generic_neg; 
	
		# get the "composition" of the cluster	
		($ref,$comp_ok)=&determine_cluster_composition($label_generic_neg);
		if($comp_ok ne ''){
			die "Couldn't determine the composition of the negative charger ion.\n";
		}else{
			$clust_comp[$n_generic_neg] = [@$ref];
		}
	
		# get the mass and radius of the cluster
		$clust_mass[$max_cluster_number][0] = $mass_nion;
		if($l_bulk_density){
			$clust_radius[$max_cluster_number][0] = ($mass_nion/$dens_nion/4.0*3.0/$pi)**(1.0/3.0);
		}
	}
	if($include_generic_pos){
		$n_generic_pos=$max_cluster_number+1;
		$lnegative[$n_generic_pos] = 0;
		$lpositive[$n_generic_pos] = 1;
		$lneutral_cluster[$n_generic_pos] = 0;
		$label[$n_generic_pos]=$label_generic_pos;
		$comments[$n_generic_pos]="% Cluster $n_generic_pos: generic positive ion ";
		$max_cluster_number=$n_generic_pos;
	
		# get the "composition" of the cluster	
		($ref,$comp_ok)=&determine_cluster_composition($label_generic_pos);
		if($comp_ok ne ''){
			die "Couldn't determine the composition of the positive charger ion.\n";
		}else{
			$clust_comp[$n_generic_pos] = [@$ref];
		}
	
		# get the mass and radius of the cluster
		$clust_mass[$max_cluster_number][0] = $mass_pion;
		if($l_bulk_density){
			$clust_radius[$max_cluster_number][0] = ($mass_pion/$dens_pion/4.0*3.0/$pi)**(1.0/3.0);
		}
	}

################################################################################
################################################################################
##      2.8 Adding cluster numbers for other fluxes                           ##
################################################################################
################################################################################

	$i = 0;
	$i = $i + 1;
	$n_flux_source=$max_cluster_number+$i;
	$i = $i + 1;
	$n_flux_coag=$max_cluster_number+$i;
	$i = $i + 1;
	$n_flux_wall=$max_cluster_number+$i;
	$i = $i + 1;
	$n_flux_dil=$max_cluster_number+$i;
	$i = $i + 1;
	$n_flux_rec=$max_cluster_number+$i;
	$label[$n_flux_source]="source";
	$label[$n_flux_coag]="coag";
	$label[$n_flux_wall]="wall";
	$label[$n_flux_dil]="dil";
	$label[$n_flux_rec]="rec";
	$comments[$n_flux_source]="% $n_flux_source is for source fluxes";
	$comments[$n_flux_coag]="% $n_flux_coag is for coagulation losses ";
	$comments[$n_flux_wall]="% $n_flux_wall is for wall losses ";
	$comments[$n_flux_dil]="% $n_flux_dil is for dilution losses ";
	$comments[$n_flux_rec]="% $n_flux_rec is for recombination of positive and negative charger ions with each other";
	if($use_coag_to_outgoing){
		$i = $i + 1;
		$n_coag_on_formed_clusters = $max_cluster_number+$i;
		$label[$n_coag_on_formed_clusters] = "coag2";
		$comments[$n_coag_on_formed_clusters] = "% $n_coag_on_formed_clusters is for coagulation losses on formed clusters";
	}
	if(!defined($n_flux_out)){
		$i = $i + 1;
		$n_flux_out=$max_cluster_number+$i;
		$n_flux_out_max=$n_flux_out;
		$label[$n_flux_out]="out_neu";
		$comments[$n_flux_out]="% $n_flux_out is for collisions that lead succesfully out of the system as neutrals";
		$lneutral[$n_flux_out] = 1;
		$lnegative[$n_flux_out] = 0;
		$lpositive[$n_flux_out] = 0;
	}
	if($include_negative_ion_terms && !defined($n_flux_out_neg)){
		$i = $i + 1;
		$n_flux_out_neg=$max_cluster_number+$i;
		$n_flux_out_max=$n_flux_out_neg;
		$label[$n_flux_out_neg]="out_neg";
		$comments[$n_flux_out_neg]="% $n_flux_out_neg is for collisions that lead succesfully out of the system as negatives";
		$lneutral[$n_flux_out_neg] = 0;
		$lnegative[$n_flux_out_neg] = 1;
		$lpositive[$n_flux_out_neg] = 0;
	}elsif(!$include_negative_ion_terms){
		$n_flux_out_neg=-1
	}
	if($include_positive_ion_terms && !defined($n_flux_out_pos)){
		$i = $i + 1;
		$n_flux_out_pos=$max_cluster_number+$i;
		$n_flux_out_max=$n_flux_out_pos;
		$label[$n_flux_out_pos]="out_pos";
		$comments[$n_flux_out_pos]="% $n_flux_out_pos is for collisions that lead succesfully out of the system as positives";
		$lneutral[$n_flux_out_pos] = 0;
		$lnegative[$n_flux_out_pos] = 0;
		$lpositive[$n_flux_out_pos] = 1;
	}elsif(!$include_positive_ion_terms){
		$n_flux_out_pos=-1
	}
	if($l_reac_in_clus){
		$i = $i + 1;
		$n_react_in_clust = $max_cluster_number+$i;
		$label[$n_react_in_clust] = "react";
		$comments[$n_react_in_clust] = "% $n_react_in_clust is for chemical reactions inside clusters";
	}

	if(!$keep_boundary_clusters){
		$i = $i + 1;
		$n_flux_bound=$max_cluster_number+$i;
		$label[$n_flux_bound]="bound";
		$comments[$n_flux_bound]="% $n_flux_bound is for collisions that lead out of the system, but the resulting cluster is brought back";
	}

	# This is the size of the flux matrix, used to loop over all fluxes EXCEPT the explicit boundary clusters, if we're keeping them
	$max_n_flux=$max_cluster_number+$i;
	# Size of the flux matrix including the boundary clusters
	# This will be updated when we find out the boundary clusters while looping through the collisions
	$max_n_flux_inc_boundary_clusters=$max_n_flux;
	@label_inc_boundary_clusters=@label;

	if($l_print_boundary){
		for($i=$max_cluster+1;$i<=$max_cluster_number;$i++){
			$str_temp = $comments[$i];
			$str_temp =~ s/% //;
			print "$str_temp\n";
		}
		for($i=$max_cluster_number+1;$i<=$max_n_flux;$i++){
			$str_temp = $comments[$i];
			$str_temp =~ s/%/Flux/;
			print "$str_temp\n";
		}
	}

}

if($l_save_outgoing){
	$size_C_out = $max_n_flux;
}else{
	$size_C_out = $max_cluster_number;
}

################################################################################
################################################################################
##      2.9 Checking which molecule types are used in the system              ##
################################################################################
################################################################################

$nmol_types_used = 0;
@mol_types_used = ();
@mol_names_used = ();
# Determine the numbers of used molecule types and their indices, as well as monomer indices
if($lloop){
	@nmol_ranges_array = ();
	$nmol_ranges_list = '';
	$nmol_ranges_ind = '';
	for($imol=0;$imol<$nmol_types;$imol++){
		if($n_max_type[$imol]>0){
			$nmol_types_used += 1;
			$mol_types_used[$nmol_types_used] = $imol;
			$mol_names_used[$nmol_types_used] = $molecule_name[$imol];
			$nmol_ranges_array[$nmol_types_used] = $n_max_type[$imol];
			$nmol_ranges_list .= "$n_max_type[$imol],";
			$nmol_ranges_ind .= "i".($nmol_types_used).", ";
			if($lloop_max_ratio){
				if($imol==$nmol_inp_max_ratio_basis){
					$nmol_max_ratio_basis = $nmol_types_used;
				}elsif($imol==$nmol_inp_not_basis){
					$nmol_not_basis = $nmol_types_used;
				}
			}
		}
	}
	$nmol_ranges_list =~ s/,$//;
	$nmol_ranges_ind =~ s/, $//;
	if(($lloop_diag || $lloop_max_ratio) && $nmol_types_used != 2){
		die "Found $nmol_types_used molecule types in the input file; it should be two for the options you try to use\n";
	}elsif($l_j_in && $nmol_types_used > 1){
		die "Found $nmol_types_used molecule types in the input file; option j_in only works with one\n";
	}
	@cluster_from_indices_array = ();
	if(!$lloop_max_ratio){
		$i = 1;
		$string1 = "n1$mol_names_used[$nmol_types_used] = $i";
		$string2 = "$i";
		$cluster_from_indices_array[$nmol_types_used] = $i;
		$string3 = "ij_ind($nmol_types_used)";
		$string4 = "";
		for($imol=$nmol_types_used-1;$imol>=1;$imol--){
			$i *= ($nmol_ranges_array[$imol+1]+1);
			$cluster_from_indices_array[$imol] = $i;
			$string1 = "n1$mol_names_used[$imol] = $i, $string1";
			$string2 = "$i, $string2";
			$string4 = "*".($nmol_ranges_array[$imol+1]+1)."$string4";
			$string3 = "ij_ind($imol)$string4+$string3";
		}
		if($nmol_types_used>1){
			# Cluster index label corresponding to composition ij_ind
			$cluster_from_indices=$string3;
		}
	}else{
		# In the max. ratio mode, the non-basis and basis molecule types are the 1st and 2nd clusters, respectively
		$string1 = '';
		$string2 = '';
		for($imol=$nmol_types_used;$imol>=1;$imol--){
			if($imol==$nmol_max_ratio_basis){
				$i = 2;
			}elsif($imol==$nmol_not_basis){
				$i = 1;
			}
			$string1 = "n1$mol_names_used[$imol] = $i, $string1";
			$cluster_from_indices_array[$imol] = $i;
			$string2 = "$i, $string2";
		}
		$string1 =~ s/, $//;
		$string2 =~ s/, $//;
		
		$cluster_from_indices="n_pure(ij_ind($nmol_max_ratio_basis))+ij_ind($nmol_not_basis)";
	}
	$lines{'monomers'} = $string1;
	$monomer_indices = $string2;
}else{
	for($imol=0;$imol<$nmol_types;$imol++){
		if($n_max_type[$imol]>0 || $n_max_type_neg[$imol]>0 || $n_max_type_pos[$imol]>0){
			$nmol_types_used += 1;
			$mol_types_used[$nmol_types_used] = $imol;
			$mol_names_used[$nmol_types_used] = $molecule_name[$imol];
		}
	}
	$string2 = '';
	for($imol=1;$imol <= $nmolecules;$imol++){
		$string2 .= "$molecules[$imol], ";
	}
	$string2 =~ s/, $//;
	$monomer_indices = $string2;
}

################################################################################
################################################################################
##   3. Dealing with forbidden reactions                                      ##
##      and reactions with nonstandard products                               ##
################################################################################
################################################################################

# If some cluster-cluster collisions are allowed, determining the corresponding cluster numbers
if($#nonmonomer_collisions_only>=0){
	@temp_array = @nonmonomer_collisions_only;
	@nonmonomer_collisions_only = ();
	for($i=0;$i<=$#temp_array;$i++){
		# Checking if the cluster definition contains a range of molecule numbers
		if($temp_array[$i] =~ /(\d+)\-(\d+)/){
			@temp_array_2 = split("$1-$2",$temp_array[$i]);
			for($j=$1;$j<=$2;$j++){
				$iclus=&get_cluster_number("$temp_array_2[0]$j$temp_array_2[1]");
				if($iclus eq ''){
					print "\nCould not find cluster $temp_array_2[0]$j$temp_array_2[1] in the system -> not using it to define allowed cluster-cluster collisions.\n";
				}else{
					push(@nonmonomer_collisions_only,$iclus);
				}
			}
		}else{
			$iclus=&get_cluster_number($temp_array[$i]);
			if($iclus eq ''){
				print "\nCould not find cluster $temp_array[$i] in the system -> not using it to define allowed cluster-cluster collisions.\n";
			}else{
				push(@nonmonomer_collisions_only,$iclus);
			}
		}
	}
}

# Now reading the file for reactions with nonstandard products
if($nonstandard_reaction_file_name eq ''){
	$l_nonstandard_allowed = 0;
}else{
	$l_nonstandard_allowed = 1;

	# array for collisions that are in principle valid but are not used
	for($iclus=1;$iclus<=$max_cluster;$iclus++){
		for($jclus=1;$jclus<=$max_cluster;$jclus++){
			$l_coll_not_used[$iclus][$jclus] = 0;
		}
	}
	# array for other non-standard reactions
	for($iclus=1;$iclus<=$max_cluster_number;$iclus++){
		for($jclus=1;$jclus<=$max_cluster_number;$jclus++){
			$nonstandard_reactions[$iclus][$jclus][0] = ();
		}
	}
	
	open(NONACC, "<$nonstandard_reaction_file_name") || die "\nCould not open energy non-accommodation file $nonstandard_reaction_file_name\n\n";
	chomp(@nonstandard_reaction_data=<NONACC>);
	close(NONACC);
	$number_of_lines=$#nonstandard_reaction_data+1;
	for($iline=0;$iline<$number_of_lines;$iline++){
		next if (($nonstandard_reaction_data[$iline] =~ /^\#/) || ($nonstandard_reaction_data[$iline] =~ /^\s*$/));
		@columns=split(' ',$nonstandard_reaction_data[$iline]);
		# the line must have an even number of entries (colliding molecules + products preceded by a coefficient)
		die "\nLines in nonaccommodation file $nonstandard_reaction_file_name have to contain colliding molecules + products preceded by coefficients.\n\n" if ($#columns % 2 == 0);
		#first two entries on the line should be the colliding molecules
		$iclus=&get_cluster_number($columns[0]);
		if($iclus eq ''){
			print "\nCould not find cluster $columns[0] in the system -> not using non-standard reaction '$nonstandard_reaction_data[$iline]'.\n";
			next;
		}
		$lmonomer1=$lmonomer[$iclus];
		$temp_label=$columns[1];
		# Now check if this line is for forbidding cluster-cluster collisions
		if(!$lmonomer1 && $#columns==1 && ($temp_label =~ /clusters/i)){
			for($jclus=1;$jclus<=$max_cluster;$jclus++){
				$lmonomer2=$lmonomer[$jclus];
				if(!$lmonomer2 && $jclus!=$n_flux_out && $jclus!=$n_flux_out_neg && $jclus!=$n_flux_out_pos){
					$l_coll_not_used[$iclus][$jclus] = 1;
					$l_coll_not_used[$jclus][$iclus] = 1;
				}
			}
			next;
		}
		$jclus=&get_cluster_number($temp_label);
		if($jclus eq ''){
			print "\nCould not find cluster $temp_label in the system -> not using non-standard reaction '$nonstandard_reaction_data[$iline]'.\n";
			next;
		}
		$lmonomer2=$lmonomer[$jclus];
		
		# If the line only contains the colliding clusters, this is interpreted as meaning that the collision is forbidden.
		if($#columns==1){
			$columns[2] = 1;
			$columns[3] = $columns[0];
			$columns[4] = 1;
			$columns[5] = $columns[1];
		}
		
		# if this line defines a forbidden collision, no point in printing equations of the type x + y -> x + y
		if($disable_useless_collisions){
			if($#columns==1 || (($#columns==5)&&(&compare_clusters($columns[3],$label[$iclus]) || &compare_clusters($columns[5],$label[$iclus]))&&($columns[2]==1.00)) || (($#columns==3)&&($lmonomer1)&&($columns[2]==2.00))){
				if(($columns[3] eq $label[$iclus]) || ($columns[3] eq $label[$jclus])){
					$l_coll_not_used[$iclus][$jclus] = 1;
					$l_coll_not_used[$jclus][$iclus] = 1;
				}
			}
		}
		# going through the rest of the line
		for($icol=2;$icol<$#columns;$icol+=2){
			# every other entry should be a positive number (coefficient for the reaction product)
			die "\n$columns[$icol] on line $iline in file $nonstandard_reaction_file_name should be a positive nuber\n\n" if (!&is_pos_number($columns[$icol]));
			if($lfortran && $columns[$icol] =~ m/\./){
				die "Non-standard reactions with non-integer coefficients don't work in the Fortran version.\n";
			}
			$coef_temp = sprintf("%.2f",$columns[$icol]);
			if($lfortran){
				$coef_temp =~ s/e/d/;
			}
			push(@{ $nonstandard_reactions[$iclus][$jclus] },$columns[$icol]);
			if($iclus != $jclus){
				push(@{ $nonstandard_reactions[$jclus][$iclus] },$columns[$icol]);
			}
			# the next entry should be the corresponding reaction product
			$temp_label=$columns[$icol+1];
			if($temp_label =~ /^out_neu$/i || $temp_label =~ /^out$/i){
				$kclus = $n_flux_out;
				$lfound = 1;
			}elsif($temp_label =~ /^out_neg$/i){
				$kclus = $n_flux_out_neg;
				$lfound = 1;
			}elsif($temp_label =~ /^out_pos$/i){
				$kclus = $n_flux_out_pos;
				$lfound = 1;
			}else{
				$kclus=&get_cluster_number($temp_label);
			}
			die "\nReaction product $temp_label on line $iline in file $nonstandard_reaction_file_name is not a valid cluster\n\n" if ($kclus eq '');
			push(@{ $nonstandard_reactions[$iclus][$jclus] },$kclus);
			if($iclus != $jclus){
				push(@{ $nonstandard_reactions[$jclus][$iclus] },$kclus);
			}
		}
	}
}

################################################################################
################################################################################
##   4. Charge balance for generic positive or negative ions                  ##
################################################################################
################################################################################

if($include_generic_neg && $include_generic_pos){
	$generic_positive_ions = '';
	$generic_negative_ions = '';
	$generic_positive_ions_corr = '';
	$generic_negative_ions_corr = '';
	
	if($l_charge_balance!=0){
		
		# Make strings for printing out the arrays of indices for negative and positive clusters
		if($include_negative_ion_terms){
			if($negative_clusters[$#negative_clusters]-$negative_clusters[0]==$#negative_clusters){
				$string1 = "$negative_clusters[0]:$negative_clusters[$#negative_clusters]";
			}else{
				$string1 = $negative_clusters[0];
				for($iclus=1;$iclus<=$#negative_clusters;$iclus++){
					# Making sure the line won't be too long for the fortran compiler
					if($lfortran && $iclus%20==0){
						$string1 = "$string1, &\n\t\t\t\t\t\t&$negative_clusters[$iclus]";
					}else{
						$string1 = "$string1, $negative_clusters[$iclus]";
					}
				}
				if($lfortran){
					$string1 = "(/$string1/)";
				}
			}
		}
	
		if($include_positive_ion_terms){
			if($positive_clusters[$#positive_clusters]-$positive_clusters[0]==$#positive_clusters){
				$string2 = "$positive_clusters[0]:$positive_clusters[$#positive_clusters]";
			}else{
				$string2 = $positive_clusters[0];
				for($iclus=1;$iclus<=$#positive_clusters;$iclus++){
					# Making sure the line won't be too long for the fortran compiler
					if($lfortran && $iclus%20==0){
						$string2 = "$string2, &\n\t\t\t\t\t\t&$positive_clusters[$iclus]";
					}else{
						$string2 = "$string2, $positive_clusters[$iclus]";
					}
				}
			}
		}
	
		# First the case where we have just negative clusters
		if($include_negative_ion_terms && !$include_positive_ion_terms){

			if($l_charge_balance>0){
				$generic_negative_ions = '';
				$generic_positive_ions = "c($n_generic_neg)+sum(c($string1))";
				if(defined($cs_cluster_neg[0])){
					$generic_positive_ions .= "+c($n_flux_out_neg)";
				}
			}elsif($l_charge_balance<0){
				die "\nNo point in balancing generic negative ions with positive clusters, if there are no positive clusters!\n\n";
			}

		# Then the same in the case that we have just positive clusters
		}elsif(!$include_negative_ion_terms && $include_positive_ion_terms){

			if($l_charge_balance<0){
				$generic_positive_ions = '';
				$generic_negative_ions = "c($n_generic_pos)+sum(c($string2))";
				if(defined($cs_cluster_pos[0])){
					$generic_negative_ions .= "+c($n_flux_out_pos)";
				}
			}elsif($l_charge_balance>0){
				die "\nNo point in balancing generic positive ions with negative clusters, if there are no negative clusters!\n\n";
			}

		# If we have them both
		}elsif($include_negative_ion_terms && $include_positive_ion_terms){

			# This time we want the difference between positive and negative clusters, not including the chargers
			# The concentrations of charger ions of both signs will be set in the Matlab or Fortran equation files
			if($l_charge_balance>0){
				# Sum of negative clusters minus the sum of positive clusters
				if($lfortran && index($string1,',')>=0 && index($string2,',')>=0){
					# Making sure the line won't be too long for the fortran compiler
					$generic_positive_ions_corr = "sum(c($string1))&\n\t\t\t\t\t\t&-sum(c($string2))";
				}else{
					$generic_positive_ions_corr = "sum(c($string1))-sum(c($string2))";
				}
				if(defined($cs_cluster_neg[0])){
					$generic_positive_ions_corr .= "+c($n_flux_out_neg)";
				}
				if(defined($cs_cluster_pos[0])){
					$generic_positive_ions_corr .= "-c($n_flux_out_pos)";
				}
			}elsif($l_charge_balance<0){
				# Sum of positive clusters minus the sum of negative clusters
				if($lfortran && index($string1,',')>=0 && index($string2,',')>=0){
					# Making sure the line won't be too long for the fortran compiler
					$generic_negative_ions_corr = "sum(c($string2))&\n\t\t\t\t\t\t&-sum(c($string1))";
				}else{
					$generic_negative_ions_corr = "sum(c($string2))-sum(c($string1))";
				}
				if(defined($cs_cluster_pos[0])){
					$generic_negative_ions_corr .= "+c($n_flux_out_pos)";
				}
				if(defined($cs_cluster_neg[0])){
					$generic_negative_ions_corr .= "-c($n_flux_out_neg)";
				}
			}
		}
	}
}

################################################################################
################################################################################
##   5. Reading the Delta G file                                              ##
################################################################################
################################################################################

### Read the DeltaG file for the evaporation coefficients, and find out which clusters
### have hydrates if we're using the hydrate distributions

# Keep track whether there are missing or multiply defined energies, and if so, stop after going through all clusters.
$stopping=0;

if($lhydrate_average){
	# Logical that tells if the cluster has hydrates
	for($iclus=1;$iclus <= $max_cluster_number;$iclus++){
		$lhas_hydr[$iclus] = 0;
	}
}

# Energies are needed for evaporations and hydrate distributions, but not necessarily for evaporation rate mcmc.
$evap_comment = '%Note: no energy file was given/read';
if((!$disable_evap || $lhydrate_average) && !($l_evap_factor_mcmc && !$free_energy_file_name) && !$lloop){
	
	if(!$free_energy_file_name){
		die "\nYou didn't give cluster energies!\n\n";
	}
	# read in the energy values, in kcal/mol
	open(GFILE, "<$free_energy_file_name") || die "\nCould not open the delta G file $free_energy_file_name\n\n";
	chomp(@evap_data=<GFILE>);
	close(GFILE);
	$number_of_lines=$#evap_data+1;
	# first line has the pressure in Pa
	if(&is_pos_number($evap_data[0])){
		$pressure=$evap_data[0];
	}else{
		die "\nFirst line of $free_energy_file_name has to be the reference pressure, in Pa.\n\n";
	}
	
################################################################################
################################################################################
##      5.1 Reading Delta G (or H and S) for each cluster and possibly hydrates#
################################################################################
################################################################################

	if ($l_change_temperature){
		if($lfortran){
			$evap_comment = "\t!  Converting thermodynamics data from temperature $temperature_inp K to $temperature K.\n";
		}else{
			$evap_comment = "% Converting thermodynamics data from temperature $temperature_inp K to $temperature K.\n";
		}
	}else{
		$evap_comment = '';
	}
	$temp_val=$pressure/$boltz/$temperature;

	# Going through the clusters in the energy file
	for($iline=2;$iline<$number_of_lines;$iline++){
		next if ($evap_data[$iline] =~ /^\#/);
		next if ($evap_data[$iline] =~ /^\s*$/);
		$evap_data[$iline] =~ s/^\s*//;
		@columns=split(/\b\s+/,$evap_data[$iline]);
#		# if there's white space on the line before the cluster name, remove the corresponding empty column
#		if($columns[0] eq ''){
#			shift(@columns);
#		}
		if ($#columns != 1 && $#columns != 2){
			die "Error on line " . ($iline+1) . " of the energy file: each line must have a cluster name, one or two numbers (G or H and S, respectively) and nothing else.\n";
		}

		$temp_label=$columns[0];
		
		# Checking if we have the cluster (or the corresponding dry one 
		# if we are using hydrates) in the system
		($iclus,$nwaters)=&get_corresponding_dry($temp_label);

		# Save the energy value if we have the cluster (or its dry counterpart) in the system
		if($iclus ne ''){
			# Each cluster must appear only once in the energy file, so the energy should not be defined yet
			if(exists $delta_g_hydr[$iclus][$nwaters]){
				$stopping = 1;
				print "Energy of cluster $temp_label is multiply defined in the energy file $free_energy_file_name.\n";
			}
			# If there are two values, they should be the S and H values
			if($#columns == 2){
				if($variable_temp){
					if($lfortran){
						if($columns[1] =~ m/\./){
							$temp_val = "($columns[1]d0/temperature";
						}else{
							$temp_val = "($columns[1].d0/temperature";
						}
						if($columns[2] =~ m/\./){
							$temp_val = "$temp_val-($columns[2]d0)/1.d3)";
						}else{
							$temp_val = "$temp_val-($columns[2].d0)/1.d3)";
						}
					}else{
						$temp_val = "($columns[1]/temperature-($columns[2])/1000.0)";
					}
				}else{
					$temp_val = $columns[1]-$temperature*$columns[2]/1000.0; # H (and G) in kcal/mol,S in cal/mol/K
					# convert to J from kcal/mol
					$temp_val = $temp_val*$kcal_per_mol_to_J;
				}

			# If there is just one value, it should be the deltaG value in kcal/mol
			}else{
				if($variable_temp){
					die "\nThe option --variable_temp can only be used with H and S values!\n\n";
				}elsif($l_change_temperature){
					die "\nEnergies can only be converted to a different temperature using H and S values!\n\n";
				}
				$temp_val = $columns[1];
				# convert to J from kcal/mol
				$temp_val = $temp_val*$kcal_per_mol_to_J;
			}
			if($nwaters==0){
				$delta_g[$iclus] = $temp_val;
				$delta_g_hydr[$iclus][0] = $temp_val;
			}else{
				$delta_g_hydr[$iclus][$nwaters] = $temp_val;
				$lhas_hydr[$iclus]=1;
				$clust_mass[$iclus][$nwaters] = $clust_mass[$iclus][0]+$nwaters*$mass_water;
				if($l_bulk_density){
					$volume = $clust_volume[$iclus]+$nwaters*$volume_water;
#print "$clust_volume[$iclus]+$nwaters*$volume_water=".($clust_volume[$iclus]+$nwaters*$volume_water)."\n";
					$clust_radius[$iclus][$nwaters] = ($volume/4.0*3.0/$pi)**(1.0/3.0);
#print "$temp_label: mass = $clust_mass[$iclus][$nwaters], radius =  $clust_radius[$iclus][$nwaters]\n";
				}
			}
		}
	}
	# checking that all energies were given and setting undefined monomer energies to zero
	for($iclus=1;$iclus<=$max_cluster;$iclus++){
		next if($iclus==$n_flux_out || $iclus==$n_flux_out_neg || $iclus==$n_flux_out_pos);
		if(!(exists $delta_g[$iclus])){
			if($lmonomer[$iclus]){
				if($variable_temp){
					if($lfortran){
						$delta_g[$iclus]='0.d0';
					}else{
						$delta_g[$iclus]='0';
					}
				}else{
					$delta_g[$iclus]=0.0;
					$delta_g_hydr[$iclus][0]=$delta_g[$iclus];
				}
			}elsif($disable_evap_only[$iclus]==0){
				$stopping = 1;
				print "Energy of cluster $label[$iclus] was not given in the energy file $free_energy_file_name.\n";
			}
		}
	}
}

if($stopping){
	die;
}

################################################################################
################################################################################
##      5.2 Calculating hydrate distributions                                 ##
################################################################################
################################################################################

for($iclus=1;$iclus<=$max_cluster_number;$iclus++){
	# Initializing the distribution to the dry cluster for all clusters, also when hydrates are not used
	$distr[$iclus]=[(1)];
	if($lhydrate_average){
		# See if we have the hydrates of this cluster
		$hydrates_str[$iclus]='';
		if($lhas_hydr[$iclus]==1){
			@delta_g_temp = @{$delta_g_hydr[$iclus]};
			
			# Unnormalized hydrate distribution
			$sum=0.0;
			@distr_temp=();
			for($ihydr=0;$ihydr<=$#delta_g_temp;$ihydr++){
				if(defined $delta_g_temp[$ihydr]){
					if($ihydr==0){
						# The first element is zero (dry cluster)
						$delta_g_temp[$ihydr] = 0.0;
					}else{
						#Substracting the energy of the dry cluster
						$delta_g_temp[$ihydr] -= $delta_g[$iclus];
						$hydrates_str[$iclus]="$hydrates_str[$iclus]+ $label[$iclus]$ihydr$name_water\t";
					}
					$distr_temp[$ihydr] = ($rh/100*$psw/$pressure)**$ihydr*exp(-$delta_g_temp[$ihydr]/$boltz/$temperature);
					$sum += $distr_temp[$ihydr];
				}else{
					print"DeltaG for the hydrate of $label[$iclus] with $ihydr waters is not given, although the cluster has hydrates with more than $ihydr waters.\n";
				}
			}
			# Normalizing the hydrate distribution
			$sum2=0.0;
			for($ihydr=0;$ihydr<=$#distr_temp;$ihydr++){
				if(defined $distr_temp[$ihydr]){
					$distr_temp[$ihydr] /= $sum;
					$sum2 += $distr_temp[$ihydr];
				}
			}
			
			# If the sum of the normalized distribution is not one, don't use the hydrates for this cluster
			if($sum2 > 0.99){
				$distr[$iclus]=[@distr_temp];
			}else{
				print"Found hydrates for $label[$iclus], but the sum of the normalized distribution is $sum2\n";
				print"\t--> Not using the hydrate distribution for $label[$iclus]\n";
			}			
		}
	}	
}

################################################################################
################################################################################
##   6. Making the equations                                                  ##
################################################################################
################################################################################

if(!$lloop){

	# Arrays for indices of collision, evaporation and cluster reaction rates
	for($iclus=1;$iclus <= $max_n_flux;$iclus++){
		$ind_quad_loss[$iclus][0] = 0;
		$ind_quad_form[$iclus][0] = 0;
		$ind_quad_loss_extra[$iclus][0] = 0;
		$ind_quad_form_extra[$iclus][0] = 0;
		$ind_lin_loss[$iclus][0] = 0;
		$ind_lin_form[$iclus][0] = 0;
		$ind_reac_loss[$iclus][0] = 0;
		$ind_reac_form[$iclus][0] = 0;
	}

	# Keep track whether there are missing collision products, and if so, stop at the end of the loops.
	$stopping = 0;

################################################################################
################################################################################
##      6.1 Looping over the clusters and going through all possible collisions#
################################################################################
################################################################################

	EQUATION: for($iclus=1;$iclus<=$max_cluster_number;$iclus++){

		$iclus_label=$label[$iclus];

########################################################################
# Collisions that remove this cluster and evaporations that create it
########################################################################

		for($jclus=$iclus;$jclus<=$max_cluster_number;$jclus++){
			$jclus_label=$label[$jclus];

			# First find out if the collision and evaporation can happen
			$lout=0;
			$l_nonstand=0;
			$l_coag=0;
			%nmonomers = ();
			$coll_coef_temp='';
			$evap_coef_temp='';
			$str_temp = '';
		
			# If the collision is forbidden, forget about it.
			if($l_nonstandard_allowed && $l_coll_not_used[$iclus][$jclus]){
				$valid_coll=0;
				$valid_evap=0;
				if($l_print_boundary){
					print "Not using collision $iclus_label + $jclus_label\n";
				}
				
			# Next option: do we have some special rule for this collision?
			}elsif($l_nonstandard_allowed && $#{ $nonstandard_reactions[$iclus][$jclus] }>0){
				$valid_coll = 1;
				$valid_evap = 0;
				$l_nonstand = 1;
			
			# Next option: these are the generic ions and nothing is formed from their recombination
			}elsif(($iclus==$n_generic_neg && $jclus==$n_generic_pos) || ($iclus==$n_generic_pos && $jclus==$n_generic_neg)){
				$valid_coll = 1;
				$valid_evap = 0;
				$kclus = $n_flux_rec;
			
			# Next option: one or both of the clusters correspond to the pool of formed clusters
			}elsif($jclus==$n_flux_out || $jclus==$n_flux_out_neg || $jclus==$n_flux_out_pos){
				if(($lnegative[$iclus] && $lnegative[$jclus]) || ($lpositive[$iclus] && $lpositive[$jclus])){
					$valid_coll = 0;
					next;
				}else{
					$valid_coll = 1;
					$valid_evap = 0;
					$kclus = $n_coag_on_formed_clusters;
					if($l_print_boundary){
						if($iclus==$n_flux_out || $iclus==$n_flux_out_neg || $iclus==$n_flux_out_pos){
							print "Self coagulation of formed particles: $label[$iclus] + $label[$jclus] -> ";
						}else{
							print "Coagulation on formed particles: $label[$iclus] + $label[$jclus] -> ";
						}
					}
					if($lnegative[$iclus]+$lnegative[$jclus]-$lpositive[$iclus]-$lpositive[$jclus]==0){
						if($l_print_boundary){
							print "cs_cluster\n";
						}
						$l_coag = 1;
					}elsif($lnegative[$iclus]+$lnegative[$jclus]-$lpositive[$iclus]-$lpositive[$jclus]>0){
						if($l_print_boundary){
							print "cs_cluster_neg\n";
						}
						$l_coag = 2;
					}else{
						if($l_print_boundary){
							print "cs_cluster_pos\n";
						}
						$l_coag = 3;
					}
				}
			# Otherwise just collide the clusters and see what we get:
			}else{
				($combined,$valid_coll,$valid_evap,$valid) = &combine_labels($iclus,$jclus);
			
				# If the collision is not possible (e.g. two negative ions), go on to the next one
				if(!$valid_coll){
					next;
				}
			
				if($combined eq ''){
					die "\nCan't figure out what comes from $label[$iclus] + $label[$jclus]\n\n";
				}
			
				# If not, does it satisfy the nucleation criteria or should it be brought back to the system?
				if($valid_coll && !$valid){
					$str_temp = $combined;
					($combined,$lout,%nmonomers)=&check_boundary($combined);
					if($combined eq ''){
						print "Problem with: $label[$iclus] + $label[$jclus] -> $str_temp -> ?\n\n";
						$stopping = 1;
						next;
					}
					$valid_evap=0;
					# Probably no need to keep this collision, if the products are the same as the reactants
					if($disable_useless_collisions){
						$lmonomer1=$lmonomer[$iclus];
						$lmonomer2=$lmonomer[$jclus];
						if(($lmonomer1 && &compare_clusters($jclus_label,$combined)) ||($lmonomer2 && &compare_clusters($iclus_label,$combined))){
							$valid_coll=0;
							if($l_print_boundary){
								print "Not using boundary collision $iclus_label + $jclus_label -> $iclus_label + $jclus_label\n";
							}
						}
					}
					if($lout>0){
						if($lout==1){
							$temp_val = '_neu';
							$kclus = $n_flux_out;
						}elsif($lout==2){
							$temp_val = '_neg';
							$kclus = $n_flux_out_neg;
						}else{
							$temp_val = '_pos';
							$kclus = $n_flux_out_pos;
						}
						if($l_print_boundary){
							print "Growing succesfully out: $iclus_label + $jclus_label -> $str_temp -> out$temp_val";
							if($use_coag_to_outgoing){
								print " -> cs_cluster$temp_val\n";
							}else{
								print "\n";
							}
						}
					}elsif($valid_coll && $l_print_boundary){
						print "Brought back to boundary: $iclus_label + $jclus_label -> $str_temp -> $combined";
						for($imol=0;$imol<$nmol_types;$imol++){
							if($nmonomers{$molecule_name[$imol]}>0){
								print " + $nmonomers{$molecule_name[$imol]} $molecule_name[$imol]";
							}
						}
						print "\n";
					}
				}

				# Find the number of the (main) collision product
				if($lout==0){
					$kclus=&get_cluster_number($combined);
					if($kclus eq ''){
						die "Can't find cluster number for $label[$iclus] + $label[$jclus] -> $combined\n";
					}
				}
			}

################################################################################
################################################################################
##      6.2 Save the collision term if the collision can happen               ##
################################################################################
################################################################################
		
			if($valid_coll){
		
				# Collision coefficient for two different clusters ...
				if($iclus != $jclus){
					$coll_coef_temp="K($iclus,$jclus)";
				# ... and for two identical clusters.
				}else{
					if($lfortran){
						$coll_coef_temp="0.5d0*K($iclus,$jclus)";
					}else{
						$coll_coef_temp="0.5*K($iclus,$jclus)";
					}
				}
			
				# Keep track of the number of collision processes (needed for mcmc)
				$ncoef_coll += 1;
				if($lneutral_cluster[$iclus] && $lneutral_cluster[$jclus]){
					$ncoef_coll_n_n += 1;
				}elsif($lneutral_cluster[$iclus] || $lneutral_cluster[$jclus]){
					if($iclus>$max_cluster || $jclus>$max_cluster){
						$ncoef_coll_n_g += 1;
					}else{
						$ncoef_coll_n_i += 1;
					}
				}else{
					$ncoef_coll_rec += 1;
				}
				
			
				# Add this collision into the flux array
				if(!$l_nonstand && $l_coag==0){
					#reaction coefficients
					$coef_quad{"$iclus,$jclus,$kclus"} = "$coll_coef_temp";
					$coef_quad{"$jclus,$iclus,$kclus"} = "$coll_coef_temp";
					#keeping track of indices
					$ind_quad_loss[$iclus][0] += 1;
					push(@{ $ind_quad_loss[$iclus] },$jclus);
					push(@{ $ind_quad_loss[$iclus] },$kclus);
					if(($iclus!=$jclus) || $lfortran){
						$ind_quad_loss[$jclus][0] += 1;
						push(@{ $ind_quad_loss[$jclus] },$iclus);
						push(@{ $ind_quad_loss[$jclus] },$kclus);
					}
					$ind_quad_form[$kclus][0] += 1;
					push(@{ $ind_quad_form[$kclus] },$iclus);
					push(@{ $ind_quad_form[$kclus] },$jclus);
					if(!$valid && $lout==0 && $kclus!=$n_flux_rec){
						# Brought back from the boundary
						$l_quad_bound{"$iclus,$jclus,$kclus"} = 1;
						$l_quad_bound{"$jclus,$iclus,$kclus"} = 1;
						# Keeping track of the number of products
						$n_coll_products[$iclus][$jclus] = 1;
						# Saving the boundary cluster, if needed
						if($keep_boundary_clusters){
							# See if this boundary cluster is already listed
							# If not, create a new element
							if($str_temp ~~ @label_inc_boundary_clusters){
								$lfound_bound = 0;
								for($i=$max_n_flux;$i<=$max_n_flux_inc_boundary_clusters;$i++){
									if($str_temp eq $label_inc_boundary_clusters[$i]){
										$kclus_bound = $i;
										$lfound_bound = 1;
										#print "Found $str_temp\n";
										last;
									}
								}
								if(!$lfound_bound){
									die "Can't find boundary cluster $str_temp\n";
								}
							}else{
								#print "Adding $str_temp\n";
								$kclus_bound = $max_n_flux_inc_boundary_clusters+1;
								$max_n_flux_inc_boundary_clusters = $kclus_bound;
								$label_inc_boundary_clusters[$kclus_bound] = $str_temp;
								$comments[$kclus_bound]="% Cluster $kclus_bound: $label_inc_boundary_clusters[$kclus_bound], a boundary cluster that is brought back into the system";
							}						
							$clus_bound_ind{"$iclus,$jclus,$kclus"} = $kclus_bound;
							$clus_bound_ind{"$jclus,$iclus,$kclus"} = $kclus_bound;					
						}
						# Add the fluxes to monomers into a separate array
						for($imol=0;$imol<$nmol_types;$imol++){
							if($nmonomers{$molecule_name[$imol]} > 0){
								$monomer_clus=$clust_num_of_monomer[$imol];
								if($monomer_clus eq ''){
									die "How can $molecule_name[$imol] evaporate, since it is not included in the system?\n";
								}
								#reaction coefficients
								$coef_quad{"$iclus,$jclus,$monomer_clus"} = $coll_coef_temp;
								$coef_quad{"$jclus,$iclus,$monomer_clus"} = $coll_coef_temp;
								#keeping track of indices
								$ind_quad_loss_extra[$iclus][0] += 1;
								push(@{ $ind_quad_loss_extra[$iclus] },$jclus);
								push(@{ $ind_quad_loss_extra[$iclus] },$monomer_clus);
								push(@{ $ind_quad_loss_extra[$iclus] },$nmonomers{$molecule_name[$imol]});
								$ind_quad_loss_extra[$jclus][0] += 1;
								push(@{ $ind_quad_loss_extra[$jclus] },$iclus);
								push(@{ $ind_quad_loss_extra[$jclus] },$monomer_clus);
								push(@{ $ind_quad_loss_extra[$jclus] },$nmonomers{$molecule_name[$imol]});
								$ind_quad_form_extra[$monomer_clus][0] += 1;
								push(@{ $ind_quad_form_extra[$monomer_clus] },$iclus);
								push(@{ $ind_quad_form_extra[$monomer_clus] },$jclus);
								push(@{ $ind_quad_form_extra[$monomer_clus] },$nmonomers{$molecule_name[$imol]});
								$l_quad_bound{"$iclus,$jclus,$monomer_clus"} = 1;
								$l_quad_bound{"$jclus,$iclus,$monomer_clus"} = 1;
								# Keeping track of the number of products
								$n_coll_products[$iclus][$jclus] += $nmonomers{$molecule_name[$imol]};
								# Marking the boundary cluster also for the monomers
								if($keep_boundary_clusters){
									$clus_bound_ind{"$iclus,$jclus,$monomer_clus"} = $kclus_bound;
									$clus_bound_ind{"$jclus,$iclus,$monomer_clus"} = $kclus_bound;
								}
							}
						}
					}
				}elsif($l_nonstand){
					# Now let's deal with the non-standard reactions.
					# The first product is considered the main one and included in the main flux array.
				
					$coef_temp = $nonstandard_reactions[$iclus][$jclus][1];
					$n_coll_products[$iclus][$jclus] = $coef_temp;
					$kclus=$nonstandard_reactions[$iclus][$jclus][2];
					if($l_print_boundary){
						print "Nonstandard collision: $label[$iclus] + $label[$jclus] -> $coef_temp $label[$kclus]";
					}
				
					#reaction coefficients
					$coef_quad{"$iclus,$jclus,$kclus"} = $coll_coef_temp;
					$coef_quad{"$jclus,$iclus,$kclus"} = $coef_quad{"$iclus,$jclus,$kclus"};
					if($coef_temp == 1){
						$coef_quad_form{"$iclus,$jclus,$kclus"} = $coll_coef_temp;
					}else{
						if($lfortran){
							if($2 =~ m/\./){
								$coef_temp = $coef_temp . "d0";
							}else{
								$coef_temp = $coef_temp . ".d0";
							}
						}
						$coef_quad_form{"$iclus,$jclus,$kclus"} = "$coef_temp*$coll_coef_temp";
					}
					$coef_quad_form{"$jclus,$iclus,$kclus"} = $coef_quad_form{"$iclus,$jclus,$kclus"};
					# Keeping track of the number of products
				
					#keeping track of indices
					$ind_quad_loss[$iclus][0] += 1;
					push(@{ $ind_quad_loss[$iclus] },$jclus);
					push(@{ $ind_quad_loss[$iclus] },$kclus);
					if(($iclus!=$jclus) || $lfortran){
						$ind_quad_loss[$jclus][0] += 1;
						push(@{ $ind_quad_loss[$jclus] },$iclus);
						push(@{ $ind_quad_loss[$jclus] },$kclus);
					}
					$ind_quad_form[$kclus][0] += 1;
					push(@{ $ind_quad_form[$kclus] },$iclus);
					push(@{ $ind_quad_form[$kclus] },$jclus);
				
					# Mark this as a non-standard reaction
					$l_quad_nst{"$iclus,$jclus,$kclus"} = 1;
				
					# Then the other products, these go to the same array as molecules coming from the boundary
					for($icol=3;$icol<$#{ $nonstandard_reactions[$iclus][$jclus] }; $icol+=2){
					
						$coef_temp = $nonstandard_reactions[$iclus][$jclus][$icol];
						# Keeping track of the number of products
						$n_coll_products[$iclus][$jclus] += $coef_temp;
						$kclus=$nonstandard_reactions[$iclus][$jclus][$icol+1];
						if($l_print_boundary){
							print " + $coef_temp $label[$kclus]";
						}
					
						$coef_quad{"$iclus,$jclus,$kclus"} = $coll_coef_temp;
						$coef_quad{"$jclus,$iclus,$kclus"} = $coef_quad{"$iclus,$jclus,$kclus"};
						# If the coefficient is an integer, put it in the collision arrays like for boundary collisions
						if($coef_temp =~ /^\d+$/){
							$coef_quad_form{"$iclus,$jclus,$kclus"} = $coll_coef_temp;
						# otherwise print it with the collision coefficient
						}else{
							if($lfortran){
								$coef_temp = $coef_temp . "d0";
							}
							$coef_quad_form{"$iclus,$jclus,$kclus"} = "$coef_temp*$coll_coef_temp";
							$coef_temp = 1;
						}
						$coef_quad_form{"$jclus,$iclus,$kclus"} = $coef_quad_form{"$iclus,$jclus,$kclus"};
					
						#keeping track of indices
						$ind_quad_loss_extra[$iclus][0] += 1;
						push(@{ $ind_quad_loss_extra[$iclus] },$jclus);
						push(@{ $ind_quad_loss_extra[$iclus] },$kclus);
						push(@{ $ind_quad_loss_extra[$iclus] },$coef_temp);
						$ind_quad_loss_extra[$jclus][0] += 1;
						push(@{ $ind_quad_loss_extra[$jclus] },$iclus);
						push(@{ $ind_quad_loss_extra[$jclus] },$kclus);
						push(@{ $ind_quad_loss_extra[$jclus] },$coef_temp);
						$ind_quad_form_extra[$kclus][0] += 1;
						push(@{ $ind_quad_form_extra[$kclus] },$iclus);
						push(@{ $ind_quad_form_extra[$kclus] },$jclus);
						push(@{ $ind_quad_form_extra[$kclus] },$coef_temp);
					}
					if($l_print_boundary){
						print "\n";
					}
				}else{
					# Now let's deal with the coagulation to formed clusters.
				
					#reaction coefficients
					$coef_quad{"$iclus,$jclus,$kclus"} = "$coll_coef_temp";
					$coef_quad{"$jclus,$iclus,$kclus"} = "$coll_coef_temp";
				
					# $iclus is lost, nothing happens to $jclus or it changes charging state, no new clusters are formed
					#keeping track of indices
					$ind_quad_loss[$iclus][0] += 1;
					push(@{ $ind_quad_loss[$iclus] },$jclus);
					push(@{ $ind_quad_loss[$iclus] },$kclus);
					$ind_quad_form[$kclus][0] += 1;
					push(@{ $ind_quad_form[$kclus] },$iclus);
					push(@{ $ind_quad_form[$kclus] },$jclus);
					if(!$lneutral[$iclus]){
						if($l_coag == 1){
							$kclus = $n_flux_out;
						}elsif($l_coag == 2){
							$kclus = $n_flux_out_neg;
						}else{
							$kclus = $n_flux_out_pos;
						}
						$ind_quad_loss[$jclus][0] += 1;
						push(@{ $ind_quad_loss[$jclus] },$iclus);
						push(@{ $ind_quad_loss[$jclus] },$kclus);
						$ind_quad_form[$kclus][0] += 1;
						push(@{ $ind_quad_form[$kclus] },$iclus);
						push(@{ $ind_quad_form[$kclus] },$jclus);
					}
				}
					
			}	

################################################################################
################################################################################
##      6.3 Save the corresponding evaporation if it can happen               ##
################################################################################
################################################################################

			if($valid_evap){

				# Evaporation coefficient
				$evap_coef_temp= "E($iclus,$jclus)";
			
				#reaction coefficients
				$coef_lin{"$iclus,$jclus,$kclus"} = "$evap_coef_temp";
				$coef_lin{"$jclus,$iclus,$kclus"} = "$evap_coef_temp";
			
				#keeping track of indices
				$ind_lin_form[$iclus][0] += 1;
				push(@{ $ind_lin_form[$iclus] },$jclus);
				push(@{ $ind_lin_form[$iclus] },$kclus);
				if(($iclus!=$jclus) || $lfortran){
					$ind_lin_form[$jclus][0] += 1;
					push(@{ $ind_lin_form[$jclus] },$iclus);
					push(@{ $ind_lin_form[$jclus] },$kclus);
				}
				$ind_lin_loss[$kclus][0] += 1;
				push(@{ $ind_lin_loss[$kclus] },$iclus);
				push(@{ $ind_lin_loss[$kclus] },$jclus);
			
				# Keep track of the number of evaporation processes (needed for mcmc)
				$ncoef_evap += 1;
			}
		} # End of $jclus loop
	
################################################################################
################################################################################
##      6.4 External losses                                                   ##
################################################################################
################################################################################

###########################
## Coagulation sink terms
###########################

		if(!$disable_coag_sinks){
			if(! $lneutral_cluster[$iclus]){
				$cs_term_temp="fcs*";
			}else{
				$cs_term_temp="";
			}
			if($one_cs_for_all){
				$coef_lin{"$n_flux_coag,$n_flux_coag,$iclus"} = $cs_term_temp . "CS";
			}else{
				$coef_lin{"$n_flux_coag,$n_flux_coag,$iclus"} = $cs_term_temp . "CS($iclus)";
			}
			if($lfortran){
				$coef_lin{"$n_flux_coag,$n_flux_coag,$iclus"} =~ s/CS/cs/;
			}
			$ind_lin_loss[$iclus][0] += 1;
			push(@{ $ind_lin_loss[$iclus] },$n_flux_coag);
			push(@{ $ind_lin_loss[$iclus] },$n_flux_coag);
			$ind_lin_form[$n_flux_coag][0] += 1;
			push(@{ $ind_lin_form[$n_flux_coag] },$n_flux_coag);
			push(@{ $ind_lin_form[$n_flux_coag] },$iclus);
		}

###########################
## Wall loss terms
###########################

		if(!$disable_wall_terms){
			if(! $lneutral_cluster[$iclus] && !$l_wl_mcmc){
				$wl_term_temp="fwl*";
			}else{
				$wl_term_temp="";
			}
			$l_special_mcmc_wl[$iclus] = 0;
			if($one_wl_for_all){
				$coef_lin{"$n_flux_wall,$n_flux_wall,$iclus"} = $wl_term_temp . "WL";
			}else{
				$lfound = 0;
				if($#wl_mcmc_only>=0){
					for($icol=0;$icol<=$#wl_mcmc_only;$icol++){
						if(&compare_clusters($label[$iclus],$wl_mcmc_only[$icol])){
							$coef_lin{"$n_flux_wall,$n_flux_wall,$iclus"} = "WL($iclus)";
							$ncoef_wl++;
							$lfound = 1;
							$l_special_mcmc_wl[$iclus] = 1;
							last;
						}
					}
				}
				if(!$lfound){
					$coef_lin{"$n_flux_wall,$n_flux_wall,$iclus"} = $wl_term_temp . "WL($iclus)";
				}
			}
			if($lfortran){
				$coef_lin{"$n_flux_wall,$n_flux_wall,$iclus"} =~ s/WL/wl/;
			}
			$ind_lin_loss[$iclus][0] += 1;
			push(@{ $ind_lin_loss[$iclus] },$n_flux_wall);
			push(@{ $ind_lin_loss[$iclus] },$n_flux_wall);
			$ind_lin_form[$n_flux_wall][0] += 1;
			push(@{ $ind_lin_form[$n_flux_wall] },$n_flux_wall);
			push(@{ $ind_lin_form[$n_flux_wall] },$iclus);
		}

##################
## Dilution loss
##################

		if(!$disable_dilution){
			$coef_lin{"$n_flux_dil,$n_flux_dil,$iclus"} = "dil";
			$ind_lin_loss[$iclus][0] += 1;
			push(@{ $ind_lin_loss[$iclus] },$n_flux_dil);
			push(@{ $ind_lin_loss[$iclus] },$n_flux_dil);
			$ind_lin_form[$n_flux_dil][0] += 1;
			push(@{ $ind_lin_form[$n_flux_dil] },$n_flux_dil);
			push(@{ $ind_lin_form[$n_flux_dil] },$iclus);
		}


	} # End of $iclus loop
	if($stopping == 1){
		die "Use --nonstandard_reaction_file to specify the outcome of the collisions that haven't been brought back succesfully from the boundary\n";
	}

################################################################################
################################################################################
##      6.5 Finding out the number of MCMC coefficients                       ##
################################################################################
################################################################################

	# number of mcmc collision coefficients
	if($l_mcmc_any){
		$ncoef_coll = 0;
		if($l_coll_factor_mcmc_n_n){
			$ncoef_coll = $ncoef_coll+$ncoef_coll_n_n;
		}
		if($l_coll_factor_mcmc_n_i){
			$ncoef_coll = $ncoef_coll+$ncoef_coll_n_i;
		}
		if($l_coll_factor_mcmc_rec){
			$ncoef_coll = $ncoef_coll+$ncoef_coll_rec;
		}
		if($l_coll_factor_mcmc_n_g){
			$ncoef_coll = $ncoef_coll+$ncoef_coll_n_g;
		}
		if($l_wl_mcmc){
			$ncoef_wl++;
		}
		$n_mcmc_coefs = $ncoef_coll+$ncoef_evap+$l_q_neg_mcmc+$l_q_pos_mcmc+$ncoef_wl;
		if(!$l_use_get_losses && $ncoef_wl>0){
			$l_use_get_losses = 1;
		}
	}

	###		Done with the equations
}

################################################################################
################################################################################
##      6.6 Time dependent source terms                                       ##
################################################################################
################################################################################

$lines{'source'} = '';
$lines{'source_def'} = '';
if($#source_function>=0){
	$lines{'source'} = "\n\t% time dependent source terms defined by functions\n";
	if($lfortran){
		$lines{'source'} =~ s/%/!/;
		$lines{'source_def'} = "\n\treal(kind(1.d0)), external :: ";
	}
	$ii=1;
	for($i=0;$i<=$#source_function;$i++){	
		@columns=split(',',$source_function[$i]);
		$source_function_clust[$i] = $columns[0];
		if ($#columns==1){
			$source_function[$i] = $columns[1];
		}elsif($#columns==0){
			$source_function[$i] = "source_$source_function_clust[$i](t)";
		}else{
			die "Can't understand the source function definition $source_function[$i] - the format should be either 'clust,function name' or 'clust'.\n";
		}
		if($source_function_clust[$i] !~ /^n_j_in$/){
			if(!$lloop){
				$valid=&check_validity($source_function_clust[$i]);
			}else{
				$valid = 1;
				($string1,$string2) = &determine_cluster_composition($source_function_clust[$i]);
				@temp_array=@$string1;
				if($string2 ne '' || sum(@temp_array)>sum(@temp_array[@mol_types_used[1..$nmol_types_used]])){
					# cluster contains molecule types not included in the system
					$valid = 0;
				}elsif($lloop_max_ratio && (($temp_array[$nmol_inp_max_ratio_basis]==0 && $temp_array[$nmol_inp_not_basis]>1) || ($temp_array[$nmol_inp_max_ratio_basis]>0 && $temp_array[$nmol_inp_not_basis]>$max_wrt_mon*$temp_array[$nmol_inp_max_ratio_basis]))){
					# wrong composition
					$valid = 0;
				}else{
					for($imol=0;$imol<$nmol_types;$imol++){
						if($temp_array[$imol]>$n_max_type[$imol]){
							# cluster outside the system
							$valid = 0;
							last;
						}
					}
				}
			}
			if(!$valid){
				die "Source term cluster $source_function_clust[$i] is not a valid cluster\n";
			}
			if(!$lloop){
				$iclus=&get_cluster_number($source_function_clust[$i]);
			}elsif($nmol_types_used>1){
				if(!$lloop_max_ratio){
					$iclus = 0;
					for($imol=1;$imol<=$nmol_types_used;$imol++){
						$iclus = $iclus+$temp_array[$mol_types_used[$imol]]*$cluster_from_indices_array[$imol];
					}
				}else{
					$iclus = 1;
					if($temp_array[$nmol_inp_max_ratio_basis] > 0){
						for($j=1;$j<=$temp_array[$nmol_inp_max_ratio_basis];$j++){
							$iclus = $iclus+$max_wrt_mon*($j-1)+1;
						}
						$iclus = $iclus+$temp_array[$nmol_inp_not_basis];
					}
				}
			}else{
				$imol=1;
				$iclus = $temp_array[$mol_types_used[$imol]];
			}
		}else{
			$iclus = $source_function_clust[$i];
			if($lloop_j_tests){
				$lines{'source'} .= "\tif (use_jtable) then\n";
				if($temperature =~ m/\./){
					$str_temp = $temperature . "d0";
				}else{
					$str_temp = $temperature . ".d0";
				}
				if($rh =~ m/\./){
					$str_temp2 = $rh . "d0/100.d0";
				}else{
					$str_temp2 = $rh . ".d0/100.d0";
				}
				$lines{'source'} .= "\t\tcs_for_table = loss(1)+sum(K(1,n_j_in:nclust)*c(n_j_in:nclust))\n";
				$lines{'source'} .= "\t\tcall interp_J_table($str_temp,$str_temp2,cs_for_table,c(1)*1.d-6,0.d0,j_temp)\n";
				$lines{'source'} .= "\t\tj_temp = j_temp*1.d6\n";
				$lines{'source'} .= "\t\tif (use_KK) then\n";
				$lines{'source'} .= "\t\t\tif (r_j_init*1.d-9 .ge. $rlim_growth_compound_str) then\n";
				$lines{'source'} .= "\t\t\t\twrite(*,*) 'Add the growth compound to the KK extrapolation!'\n";
				$lines{'source'} .= "\t\t\t\tstop\n";
				$lines{'source'} .= "\t\t\tend if\n";
				$temp_val = sprintf('%.14e',$mass[0]/$density[0]);
				$temp_val =~ s/e/d/g;
				$lines{'source'} .= "\t\t\tgr = (6.d0*$temp_val/pi)**(1.d0/3.d0)/3.d0/n_j_init**(2.d0/3.d0)*K(1,n_j_init)*c(1)\n";
				$lines{'source'} .= "\t\t\tcoags = loss(n_j_init)+sum(K(n_j_init,n_j_in:nclust)*c(n_j_in:nclust))\n";
				$lines{'source'} .= "\t\t\tsource(n_j_in) = extrapolate_J(r_j_init*1.d-9,r_j_in*1.d-9,j_temp,gr,coags)\n";
				$lines{'source'} .= "\t\telse\n";
				$lines{'source'} .= "\t\t\tsource(n_j_in) = j_temp\n";
				$lines{'source'} .= "\t\tend if\n";
				$lines{'source'} .= "\telse\n\t";
			}
		}
		print "\n---------------------------------------------------------\n";
		print "\tYou will need the function $source_function[$i]\n\tin order to compile and run the code.\n";
		print "---------------------------------------------------------\n\n";
		$lines{'source'} .= "\tsource($iclus) = $source_function[$i];\n";
		if($lloop_j_tests && $iclus =~ /^$source_function_clust[$i]$/){
			$lines{'source'} .= "\tend if\n";
		}
		if($lfortran){
			if($ii==4){
				$lines{'source_def'} =~ s/, $//;
				$lines{'source_def'} .= "\n\treal(kind(1.d0)), external :: ";
				$ii = 1;
			}
			$str_temp = $source_function[$i];
			$str_temp =~ s/\(\S+\)$//;
			$lines{'source_def'} .= "$str_temp, ";
			$ii++;
		}
	}
	$lines{'source'} .= "\n";
	if($lfortran){
		$lines{'source'} =~ s/;//g;
		$lines{'source_def'} =~ s/, $/\n/;
		if($lloop_j_tests && $l_j_in){
			$lines{'source_def'} .= "\texternal :: interp_J_table\n";
			$lines{'source_def'} .= "\treal(kind(1.d0)), external :: extrapolate_J\n";
		}
	}else{
		$lines{'source'} =~ s/\t//g;
	}
}

################################################################################
################################################################################
##      6.7 Reactions inside clusters                                         ##
################################################################################
################################################################################

if($l_reac_in_clus){
	open(CRFILE, "<$reactions_in_clusters_file") || die "\nCould not open the cluster reaction file $reactions_in_clusters_file\n\n";
	chomp(@cluster_reactions=<CRFILE>);
	close(GFILE);
	$number_of_lines=$#cluster_reactions+1;
	REACTIONLINES: for($iline=0;$iline<$number_of_lines;$iline++){
		# Reading lines from the cluster reaction file
		@columns = split(' ',$cluster_reactions[$iline]);
		# First two entries: number of reactants and number of products
		$n_react = $columns[0];
		$n_prod = $columns[1];
		# Numbers of water molecules participating in the reaction as reactants, products or catalysts
		$n_reacting_W = 0;
		$n_formed_W = 0;
		$n_catalyzing_W = 0;
		$l_forbidden_W = 0;
		# empty cluster composition array corresponding to the set of reacting molecules and the set of product molecules
		for($imol=0;$imol<$nmol_types;$imol++){
			$comp_react[$imol] = 0;
			$comp_prod[$imol] = 0;
		}
		# Checking that we have the reactant molecules in the system ...
		for($icol=2;$icol<=$n_react+1;$icol++){
			if($lhydrate_average && $columns[$icol] eq 'W'){
				$n_reacting_W += 1;
				next;
			}elsif(! exists $mol_num_from_name{$columns[$icol]}){
				print "The reactant $columns[$icol] on line ".($iline+1)." of the cluster reaction file is not a valid molecule. -> Ignoring the reaction.\n\n";
				next REACTIONLINES;
			}
			# ... and adding them to the cluster composition array corresponding to the reacting molecules unless they are hydration waters
			$comp_react[$mol_num_from_name{$columns[$icol]}] += 1;
		}
		# Checking that we have the product molecules in the system ...
		for($icol=2+$n_react;$icol<=$n_react+$n_prod+1;$icol++){
			if($lhydrate_average && $columns[$icol] eq 'W'){
				$n_formed_W += 1;
				next;
			}elsif(! exists $mol_num_from_name{$columns[$icol]}){
				print "The product $columns[$icol] on line ".($iline+1)." of the cluster reaction file is not a valid molecule. -> Ignoring the reaction.\n\n";
				next REACTIONLINES;
			}
			# ... and adding them to the cluster composition array corresponding to the produced molecules unless they are hydration waters
			$comp_prod[$mol_num_from_name{$columns[$icol]}] += 1;
		}
		$rate_fwd = $columns[$n_react+$n_prod+2];
		$rate_bck = $columns[$n_react+$n_prod+3];
		@needed_mols = ();
		@forbidden_mols = ();
		for($icol=$n_react+$n_prod+4;$icol<=$#columns;$icol++){
			if($columns[$icol] =~ /^!/){
				$temp_val = $columns[$icol];
				$temp_val =~ s/!//;
				if(exists $mol_num_from_name{$temp_val}){
					push(@forbidden_mols,$mol_num_from_name{$temp_val});
				# Water is treated separately when the $lhydrate_average option is used
				}elsif(($columns[$icol] eq 'W') && $lhydrate_average && ($rh>0)){
					$l_forbidden_W = 1;
				}
			}else{
				if(exists $mol_num_from_name{$columns[$icol]}){
					push(@needed_mols,$mol_num_from_name{$columns[$icol]});
				# Water is treated separately when the $lhydrate_average option is used
				}elsif(($columns[$icol] eq 'W') && $lhydrate_average && ($rh>0)){
					$n_catalyzing_W += 1;
				}else{
					print "The catalyzing molecule $columns[$icol] on line ".($iline+1)." of the cluster reaction file is not a valid molecule. -> Ignoring the reaction.\n\n";
					next REACTIONLINES;
				}
			}
		}
		# Going through all clusters
		for($iclus=1;$iclus<=$max_cluster;$iclus++){
			# checking both directions of the reaction
			DIRECTIONS: for($i=1;$i<=2;$i++){
				@temp_array = @{ $clust_comp[$iclus] };
				if($i==2){
					@temp_array_2 = @comp_react;
					@comp_react = @comp_prod;
					@comp_prod = @temp_array_2;
					$temp_val = $rate_fwd;
					$rate_fwd = $rate_bck;
					$rate_bck = $temp_val;
					$temp_val = $n_reacting_W;
					$n_reacting_W = $n_formed_W;
					$n_formed_W = $temp_val;
				}
				# checking if the cluster contains the reactants of this reaction (i.e. products of the initial reaction for $i=2)
				for($imol=0;$imol<$nmol_types;$imol++){
					$temp_array[$imol] -= $comp_react[$imol];
					if($temp_array[$imol]<0){
						next DIRECTIONS;
					}
				}
				
				
				# checking if it in addition contains the other needed (catalyzing) molecules
				@temp_array_2 = @temp_array;
				for($icol=0;$icol<=$#needed_mols;$icol++){
					$temp_array_2[$needed_mols[$icol]] -= 1;
					if($temp_array_2[$needed_mols[$icol]]<0){
						next DIRECTIONS;
					}
				}
				# if there are hydrate-averaged waters as catalyzing molecules, check that the cluster has enough hydrate
				if($lhydrate_average && $#{ $distr[$iclus] }<$n_catalyzing_W){
					next DIRECTIONS;
				}
				# checking that it doesn't contain the forbidden molecules (other than possibly hydration W)
				for($icol=0;$icol<=$#forbidden_mols;$icol++){
					if($temp_array_2[$forbidden_mols[$icol]]>0){
						next DIRECTIONS;
					}
				}
				
				# determining the new (dry) cluster and checking if it is a valid one
				for($imol=0;$imol<$nmol_types;$imol++){
					$temp_array[$imol] += $comp_prod[$imol];
				}
				$label = &create_label(\@temp_array);
				$jclus=&get_cluster_number($label);
				if($jclus eq ''){
					next DIRECTIONS;
				}
				
				# now we have a product, last thing to save it, unless we already have a higher rate for the process
				if(!exists $coef_reac{"$iclus,$jclus"}){
					$coef_reac{"$iclus,$jclus"} = 0;
					if($lhydrate_average){
						for($ihydr=0;$ihydr<=$#{ $distr[$iclus] };$ihydr++){
							$coef_reac{"$iclus,$jclus,$ihydr"} = 0;
						}
					}
					#keeping track of indices if this is a new reaction
					$ind_reac_form[$jclus][0] += 1;
					push(@{ $ind_reac_form[$jclus] },$n_react_in_clust);
					push(@{ $ind_reac_form[$jclus] },$iclus);
					$ind_reac_loss[$iclus][0] += 1;
					push(@{ $ind_reac_loss[$iclus] },$jclus);
					push(@{ $ind_reac_loss[$iclus] },$n_react_in_clust);
				}
				if(!$lhydrate_average){
					if($rate_fwd>$coef_reac{"$iclus,$jclus"}){
						# save reaction coefficient
						$coef_reac{"$iclus,$jclus"} = "$rate_fwd";
					}
				}else{
					# when hydrates are used, keep track of the minimum/maximum numger of waters in the initial cluster
					for($ihydr=$n_reacting_W;$ihydr<=$#{ $distr[$iclus] };$ihydr++){
						if(($ihydr>=$n_reacting_W + $n_catalyzing_W) && !$l_forbidden_W && $rate_fwd>$coef_reac{"$iclus,$jclus,$ihydr"}){
							$coef_reac{"$iclus,$jclus,$ihydr"} = $rate_fwd;
						}
					}
				}
			}
		}
	}
	if(!$lfortran){
		open(REAC, ">get_reac$append_to_file_names.m");
#		if($variable_temp){
#			print REAC "function R = get_reac(temperature)\n\n";
#		}
		print REAC "% Reaction coefficients for chemical reactions inside clusters\n\n";
		print REAC "R = nan($max_cluster_number,$max_cluster_number);\n\n";
		for($iclus=1;$iclus<=$max_cluster;$iclus++){
			for($jclus=$iclus+1;$jclus<=$max_cluster;$jclus++){
				print "$iclus, $jclus, ".$coef_reac{"$iclus,$jclus"}."\n";
				if(exists $coef_reac{"$iclus,$jclus"}){
					$string1 = '';
					$string2 = '';
					# calculating the hydrate averaged rates
					if($lhydrate_average){
						$string1 .= ', hydrate average: ';
						$string2 .= ', hydrate average: ';
						for($ihydr=0;$ihydr<=$#{ $distr[$iclus] };$ihydr++){
							$coef_reac{"$iclus,$jclus"} += $distr[$iclus][$ihydr]*$coef_reac{"$iclus,$jclus,$ihydr"};
							$string1 .= sprintf("%.2f",$distr[$iclus][$ihydr])."*".$coef_reac{"$iclus,$jclus,$ihydr"}."+";
						}
						$string1 =~ s/\+$//;
						for($jhydr=0;$jhydr<=$#{ $distr[$jclus] };$jhydr++){
							$coef_reac{"$jclus,$iclus"} += $distr[$jclus][$jhydr]*$coef_reac{"$jclus,$iclus,$jhydr"};
							$string2 .= sprintf("%.2f",$distr[$jclus][$jhydr])."*".$coef_reac{"$jclus,$iclus,$jhydr"}."+";
						}
						$string2 =~ s/\+$//;
					}
					print REAC "R($iclus,$jclus) = ".$coef_reac{"$iclus,$jclus"}.";\t% $label[$iclus] -> $label[$jclus]$string1\n";
					print REAC "R($jclus,$iclus) = ".$coef_reac{"$jclus,$iclus"}.";\t% $label[$jclus] -> $label[$iclus]$string2\n\n";
				}
			}
		}
#		if($variable_temp){
#			print REAC "\nend";
#		}
		close(REAC);
	}
}
		

################################################################################
################################################################################
##   7. Reading the cluster radius file                                       ##
################################################################################
################################################################################

### Reading in the radii, if we are not using bulk densities

if(!$l_bulk_density){
	open(RADII, "<$radius_file_name") || die "\nCould not open radius file $radius_file_name\n\n";
	chomp(@radius_data=<RADII>);
	close(RADII);
	$number_of_lines=$#radius_data+1;
	for($iline=0;$iline<$number_of_lines;$iline++){
		next if ($radius_data[$iline] =~ /^\#/);
		if($radius_data[$iline] =~ /^\s*(\S+)\s+(\S+)/){

			$temp_label=$1;

			# Checking if we have the cluster (or the corresponding dry one 
			# if we are using hydrates) in the system
			($iclus,$nwaters)=&get_corresponding_dry($temp_label);
			
			# save the radius value if we have the cluster in the system
			if($iclus ne ''){
				if(!&is_pos_number($2)){
					die "The radius given for cluster $label[$iclus] is not a number!\n";
				};
				# Each cluster must appear only once in the radius file, so the radius should not be defined yet
				if(exists $clust_radius[$iclus][$nwaters]){
					die "Radius of cluster $temp_label is multiply defined in the radius file $radius_file_name.\n";
				}
				$clust_radius[$iclus][$nwaters] = $2*1.0e-9;	# convert from nm to m
			}
		}
	}
	# Now all the cluster should have a radius value
	for($iclus=1;$iclus<=$max_cluster_number;$iclus++){
		for($ihydr=0;$ihydr<=$#{ $distr[$iclus] };$ihydr++){
			if(!(exists $clust_radius[$iclus][$ihydr])){
				if($ihydr==0){
					die "\nRadius not found in $radius_file_name for $label[$iclus]!\n\n";
				}else{
					die "\nRadius not found in $radius_file_name for $label[$iclus]$ihydr$name_water!\n\n";
				}
			}
		}
	}
}

################################################################################
################################################################################
##   8. Starting to print the output files                                    ##
################################################################################
################################################################################

################################################################################
################################################################################
##      8.1 Cluster properties for Matlab driver and Fortran system file      ##
################################################################################
################################################################################

if($lfortran){
	open(SYST, ">$system_file_name_fortran") || die "\nCould not open system file $system_file_name_fortran\n\n";
	print SYST "module acdc_system\n\nimplicit none\n\n";
	print SYST "integer, parameter :: nclust = $max_cluster_number\t\t\t\t\t\t! number of clusters, molecules and ions\n";
	print SYST "integer, parameter :: neq = $size_C_out\t\t\t\t\t\t\t! number of equations\n";
	if($l_jlim){
		print SYST "integer, parameter :: nlim = $nlim\t\t\t\t\t\t\t! number of limits over which the net flux is recorded\n";
		print SYST "real(kind(1.d0)), parameter :: r_lim(nlim) = (/$jlim_sizes_list/)*1.d-9\n";
	}
	print SYST "\n";
	if($l_mcmc_any){
		print SYST "integer, parameter :: ncoefs = " . ($n_mcmc_coefs) . "\t\t\t\t\t\t! number of mcmc coefficients\n";
		print SYST "integer, parameter :: ncoefs_coll = " . ($ncoef_coll) . "\t\t\t\t! number of sticking coefficients\n";
		print SYST "integer, parameter :: ncoefs_evap = " . ($ncoef_evap) . "\t\t\t\t! number of evaporation coefficients\n";
		if($l_q_neg_mcmc){
			print SYST "integer, parameter :: n_q_neg_coef = " . ($ncoef_coll+$ncoef_evap+1) . "\t\t\t\t! mcmc coefficient corresponding to Q_neg\n";
		}
		if($l_q_pos_mcmc){
			print SYST "integer, parameter :: n_q_pos_coef = " . ($ncoef_coll+$ncoef_evap+$l_q_neg_mcmc+1) . "\t\t\t\t! mcmc coefficient corresponding to Q_pos\n";
		}
		if($ncoef_wl>0){
			print SYST "integer, parameter :: n_wl_coefs = $ncoef_wl\t\t\t\t\t! mcmc coefficient corresponding to the wall loss term\n";
		}
		print SYST "\n";
	}
	if($temperature =~ m/\./){
		print SYST "real(kind(1.d0)), parameter :: temp = $temperature" . "d0\t\t\t\t! temperature in K\n";
	}else{
		print SYST "real(kind(1.d0)), parameter :: temp = $temperature" . ".d0\t\t\t\t! temperature in K\n";
	}
	if($rh =~ m/\./){
		print SYST "real(kind(1.d0)), parameter :: rh = $rh" . "d0\t\t\t\t\t! RH in %\n";
	}else{
		print SYST "real(kind(1.d0)), parameter :: rh = $rh" . ".d0\t\t\t\t\t! RH in %\n";
	}
	print SYST "\n";
	if(!$disable_coag_sinks){
		if($coag_terms =~ m/^exp_loss$/i){
			$str_temp = sprintf('%.14e',$exp_loss_exponent);
			$str_temp =~ s/e/d/;
			print SYST "real(kind(1.d0)), parameter :: cs_exponent_default = $str_temp\t\t\t! parameters for the exponential loss\n";
			$str_temp = sprintf('%.14e',$exp_loss_coefficient);
			$str_temp =~ s/e/d/;
			print SYST "real(kind(1.d0)), parameter :: cs_coefficient_default = $str_temp\n\n";
		}elsif($coag_terms =~ m/^bg_loss$/i){
			$str_temp = sprintf('%.14e',$bg_concentration);
			$str_temp =~ s/e/d/;
			print SYST "real(kind(1.d0)), parameter :: cbg_default = $str_temp\t\t\t! default properties of the scavengers (in SI units)\n";
			$str_temp = sprintf('%.14e',$bg_diameter);
			$str_temp =~ s/e/d/;
			print SYST "real(kind(1.d0)), parameter :: dpbg_default = $str_temp\n";
			$str_temp = sprintf('%.14e',$bg_density);
			$str_temp =~ s/e/d/;
			print SYST "real(kind(1.d0)), parameter :: rhobg_default = $str_temp\n\n";
		}
	}
	if(!$lloop){
		print SYST "integer, parameter :: ";
		for($imol=1;$imol <= $nmolecules;$imol++){
			if($imol>1){
				print SYST ", ";
			}
			print SYST "n$label[$molecules[$imol]] = $molecules[$imol]";
		}
		if($include_generic_neg){
			print SYST ", n$label[$n_generic_neg] = $n_generic_neg";
		}
		if($include_generic_pos){
			print SYST ", n$label[$n_generic_pos] = $n_generic_pos";
		}
		print SYST "\t\t\t! cluster indices for monomers and ions";
		if($l_save_outgoing){
			$str_temp =  "\ninteger, parameter :: ";
			for($k=$n_flux_out; $k<=$n_flux_out_max; $k++){
				$str_temp .= "n$label[$k] = $k, ";
			}
			$str_temp =~ s/, $//;
			print SYST "$str_temp\n";
		}
	}else{
		print SYST "integer, parameter :: $lines{'monomers'}\n";
	}
	print SYST "\n";
	print SYST "integer, parameter :: n_mol_types = $nmol_types_used\n";
	print SYST "integer, parameter :: ";
	for($imol=1;$imol<=$nmol_types_used;$imol++){
		if($imol>1){
			print SYST ", ";
		}
		print SYST "nmol$mol_names_used[$imol] = $imol";
	}
	print SYST "\t\t\t! molecule indices for the used species\n\n";
	
}else{
	
	open(DRIVER, ">$driver_file_name") || die "\nCould not open driver $driver_file_name.\n\n";
	
	# print some explanations in the beginning
	if($l_old_output){
		print DRIVER "function [C, T, clust, Cf, J_out, flux, out, outflux_matrix, sources, form, converged] = $driver_function_name(Tmax, varargin)\n";
		print DRIVER "%function [C, T, clust, Cf, J_out, flux, out, outflux_matrix, sources, form, converged] = $driver_function_name(Tmax, {C0, {T0}}, {'Option', {Value}, {...}}))\n";
	}else{
		print DRIVER "function [C, T, converged, clust, Cf, varargout] = $driver_function_name(Tmax, varargin)\n";
		print DRIVER "%function [C, T, converged, clust, Cf, J_out, flux, out, outflux_matrix, sources, form] = $driver_function_name(Tmax, {C0, {T0}}, {'File', Filename, {...}}))\n";
		print DRIVER "%function [C, T, converged, clust, Cf, J_out] = $driver_function_name(Tmax, {C0, {T0}}, 'no_dofluxes',  {'Option', {Value}, {...}}))\n";
		print DRIVER "%function [C, T, converged, clust, Cf] = $driver_function_name(Tmax, {C0, {T0}}, 'no_fluxes',  {'Option', {Value}, {...}}))\n";
	}
		
	print DRIVER "% \n";
	print DRIVER "% Input parameters:\n";
	print DRIVER "% Tmax:		     duration of the simulation in seconds (added to T0), or vector containing output time points (not added to T0).\n";
	print DRIVER "% C0:			 an optional row vector (or matrix) containing the initial concentrations (on the last row).\n";
	print DRIVER "% T0:			 if C0 contains concentrations as a function of time, the corresponding times should be given in the column vector T0.\n";
	print DRIVER "% 'No_dofluxes':  the driver does not call $flux_file_name and no related output is printed\n";
	print DRIVER "% 'No_fluxes':    the driver does not call $flux_file_name or $j_file_name and no related output is printed\n";
	print DRIVER "% 'Sources_in', filename:    input file containing source terms, concentrations to be kept constant and initial concentrations\n";
	print DRIVER "% 'Sources_out', filename:   output file for monomer source terms, only used if steady-state is reached\n";
	print DRIVER "% 'Constants_out', filename: output file for constant monomer concentrations, only used if steady-state is reached\n";
	print DRIVER "% 'C_neutr', filename:	    output file for concentrations of neutrals\n";
	print DRIVER "% 'C_neg', filename:		    output file for concentrations of negatives\n";
	print DRIVER "% 'C_pos', filename:		    output file for concentrations of positives\n";
	print DRIVER "% 'Fluxes', filename:	    output file for fluxes\n";
	print DRIVER "% 'Outmat', filename:	    output file for matrix containing collisions going out from the system\n";
	print DRIVER "% 'Cluster_data', filename:	    output mat-file for cluster names, masses, diameters etc.\n";
	print DRIVER "% 'Cfun', {cluster function variable1 variable2 ...}: function giving the concentration of some molecule or cluster as a function of time or other concentrations, example: 'Cfun', {'1A' 'acidfunction' 'T'} (acidfunction.m gives acid concentration as a function of time) or 'Cfun', {'1D' 'dma_from_ammonia' '1N' 'T'} (dma_from_ammonia.m gives DMA concentration as a function of ammonia concentration and time)\n";
	if($variable_temp){
		print DRIVER "% 'Temperature', number:	   temperature if not $temperature K\n";
	}
	print DRIVER "% \n";
	print DRIVER "% Output parameters:\n";
	print DRIVER "% C:			  cluster distribution in cm^(-3) as a function of time T\n";
	print DRIVER "% T:			  time points (seconds) corresponding to C\n";
	if(!$l_old_output){
		print DRIVER "% converged:	  1 if the simulation has converged, and 0 if not, and -1 if there are negative concentrations lower than the threshold value -1e-12 cm^(-3)\n";
	}
	print DRIVER "% clust:		  cluster names\n";
	print DRIVER "% Cf:          final cluster distribution in cm^(-3)\n";
	print DRIVER "% J_out:       particle formation rate (rate out of the system) in cm^(-3) s^(-1) as a function of time\n";
	print DRIVER "% flux:        net fluxes from i to j (cm^(-3) s^(-1)) are given in the matrix elements flux(i,j)\n";
	print DRIVER "%              in collisions/fissions involving two similar clusters, the flux is from the point of view of the\n 
";
	print DRIVER "%              SMALLER cluster, i.e. must be divided by 2 when considering the point of view of the larger cluster\n";
	print DRIVER "% out:         fluxes from i out to a cluster with j mol1 (cm^(-3) s^(-1)) are in out(i,j+1,k) with\n";
	print DRIVER "%     		  k=1: neutral, k=2: neg, k=3: pos, k=4: recombination\n";
	print DRIVER "%     		  in collisions involving two similar clusters, the flux is from the point of view of \n";
	print DRIVER "%     		  the outgrowing clusters (i.e. does NOT need to be divided by 2)\n";
	print DRIVER "% outflux_matrix:	  fluxes going out of the system by collision of clusters i and j (cm^(-3) s^(-1)) \n";
	print DRIVER "%     		  are given in the matrix elements outflux_matrix(i,j)\n";
	print DRIVER "%     		  for i=j, the flux is from the point of view of the outgrowing clusters (i.e. does NOT need to be divided by 2)\n";
	print DRIVER "% sources:	  sources of each species (cm^(-3) s^(-1))\n";
	print DRIVER "% form:	      formation rates of clusters with j mol1 in the system (cm^(-3) s^(-1)); \n";
	print DRIVER "%     		  DO NOT trust this 100% because it hasn't really been used or tested\n";
	if($l_old_output){
		print DRIVER "% converged:	  1 if the simulation has converged, 0 if not, and -1 if there are negative concentrations lower than the threshold value -1e-12 cm^(-3)\n\n";
	}

}

if(!$lloop){
	$lines{'nof1'} = "";
	$lines{'1of1'} = "";
	$lines{'0'} = "";
	$lines{'n'} = "";
	$lines{'p'} = "";
	$lines{'0c'} = "";
	$lines{'nc'} = "";
	$lines{'pc'} = "";
	$lines{'m'} = "";
	$lines{'d'} = "";
	$lines{'dm'} = "";
	$lines{'cl'} = "";
	$lines{'mon'} = "";
	$lines{'nmon'} = "";
	$lines{'t0'} = "";
	$lines{'tn'} = "";
	$lines{'tp'} = "";
	$lines{'out0'} = "";
	$lines{'outn'} = "";
	$lines{'outp'} = "";
	$i = 0;
	$j = 0;
	$k = 0;
	$l = 0;
	$ii = 0;
	$jj = 0;
	$kk = 0;
	$ll = 0;
	$jc = 0;
	$kc = 0;
	$lc = 0;
	$jclus = 0;
	$clust_mass_max = 0;
	$diameter_max = 0;
	$mob_diameter_max = 0;
	for($iclus=1;$iclus<=$max_cluster_number;$iclus++){
		$lines{'cl'} .= "'" . $label[$iclus] . "' ";
		if($lmonomer[$iclus]){
			$lines{'mon'} .= "$iclus ";
		}else{
			$lines{'nmon'} .= "$iclus ";
		}
		if($jclus==10){
			$lines{'m'} .= "&\n\t\t&";
			$lines{'d'} .= "&\n\t\t&";
			$lines{'dm'} .= "&\n\t\t&";
			$jclus = 0;
		}
		$jclus += 1;
		$mass1 = $clust_mass[$iclus][0];
		$mass1 = sprintf("%.2f",$mass1/$mass_conv);
		$clust_mass_max = max($clust_mass_max,$mass1);
		$mass1 =~ s/e/d/;
		$lines{'m'} .= "$mass1, ";
		$radius = $clust_radius[$iclus][0];
		$diameter = sprintf("%.2f",$radius*2e9);
		$diameter_max = max($diameter_max,$diameter);
		$diameter =~ s/e/d/;
		$lines{'d'} .= "$diameter, ";
		$diameter = sprintf("%.2f",(2.0e9*$radius+0.3)*sqrt(1+28.8*$mass_conv/$mass1));
		$mob_diameter_max = max($mob_diameter_max,$diameter);
		$diameter =~ s/e/d/;
		$lines{'dm'} .= "$diameter, ";
		if($lneutral_cluster[$iclus]){
			if($clust_comp[$iclus][0] == 1){
				if($ii==15){
					$lines{'1of1'} .= "&\n\t\t&";
					$ii = 0;
				}
				$lines{'1of1'} .= "$iclus, ";
				$i += 1;
				$ii += 1;
			}
			if($jj==20){
				$lines{'0'} .= "&\n\t\t&";
				$jj = 0;
				$lines{'0c'} .= "&\n\t\t&";
				$lines{'nof1'} .= "&\n\t\t&";
			}
			$lines{'0'} .= "$iclus, ";
			if(!$lmonomer[$iclus]){
				$lines{'0c'} .= "$iclus, ";
				$jc += 1;
			}
			$j += 1;
			$jj += 1;
			$jmol = 1;
			$nmol = 1;
			for($n=1;$n<$n_neutral_mols;$n++){
				$imol = $neutral_mols[$n];
				$jmol = $jmol+$nmol*$clust_comp[$iclus][$imol];
				$nmol = $nmol*($n_max_type[$imol]+1);
			}
			$lines{'nof1'} .= "$clust_comp[$iclus][0], ";
			$lines{'t0'} .= "\[" . ($clust_comp[$iclus][0]+1) . " $jmol\] ";
		}elsif($lnegative_cluster[$iclus]){
			if($kk==15){
				$lines{'n'} .= "&\n\t\t&";
				$lines{'nc'} .= "&\n\t\t&";
				$kk = 0;
			}
			$lines{'n'} .= "$iclus, ";
			if(!$lmonomer[$iclus]){
				$lines{'nc'} .= "$iclus, ";
				$kc += 1;
			}
			$k += 1;
			$kk += 1;
			$jmol = 1;
			$nmol = 1;
			for($n=1;$n<$n_neutral_mols;$n++){
				$imol = $neutral_mols[$n];
				$jmol = $jmol+$nmol*$clust_comp[$iclus][$imol];
				if($n_corr_negative_mol[$imol] > 0){
						$jmol = $jmol+$nmol*$clust_comp[$iclus][$n_corr_negative_mol[$imol]];
				}
				$nmol = $nmol*($n_max_type_neg[$imol]+1);
			}
			if($first_neg>0){
				$lines{'nof1'} .= ($clust_comp[$iclus][0]+$clust_comp[$iclus][$first_neg]) . ", ";
				$lines{'tn'} .= "\[" . ($clust_comp[$iclus][0]+$clust_comp[$iclus][$first_neg]+1) . " $jmol\] ";
			}else{
				$lines{'nof1'} .= "$clust_comp[$iclus][0], ";
				$lines{'tn'} .= "\[" . ($clust_comp[$iclus][0]+1) . " $jmol\] ";
			}
		}elsif($lpositive_cluster[$iclus]){ 
			if($ll==15){
				$lines{'p'} .= "&\n\t\t&";
				$lines{'pc'} .= "&\n\t\t&";
				$ll = 0;
			}
			$lines{'p'} .= "$iclus, ";
			if(!$lmonomer[$iclus]){
				$lines{'pc'} .= "$iclus, ";
				$lc += 1;
			}
			$l += 1;
			$ll += 1;
			$jmol = 1;
			$nmol = 1;
			for($n=1;$n<$n_neutral_mols;$n++){
				$imol = $neutral_mols[$n];
				$jmol = $jmol+$nmol*$clust_comp[$iclus][$imol];
				if($n_corr_positive_mol[$imol] > 0){
						$jmol = $jmol+$nmol*$clust_comp[$iclus][$n_corr_positive_mol[$imol]];
				}
				$nmol = $nmol*($n_max_type_pos[$imol]+1);
			}
			if($first_pos > 0){
				$lines{'nof1'} .= ($clust_comp[$iclus][0]+$clust_comp[$iclus][$first_pos]) . ", ";
				$lines{'tp'} .= "\[" . ($clust_comp[$iclus][0]+$clust_comp[$iclus][$first_pos]+1) . " $jmol\] ";
			}else{
				$lines{'nof1'} .= "$clust_comp[$iclus][0], ";
				$lines{'tp'} .= "\[" . ($clust_comp[$iclus][0]+1) . " $jmol\] ";
			}
		}
	}
	for($iclus=$max_cluster_number+1;$iclus<=$size_C_out;$iclus++){
		$lines{'cl'} .= "'" . $label[$iclus] . "' ";
	}
	if($include_generic_neg){
		$lines{'n'} .= "$n_generic_neg ";
		$lines{'nof1'} .= "0, ";
	}
	if($include_generic_pos){
		$lines{'p'} .= "$n_generic_pos ";
		$lines{'nof1'} .= "0, ";
	}
	for($imol=1;$imol<=$nmol_types_used;$imol++){
		for($irule=0;$irule<$n_nucl_rules;$irule++){
			$lines{'out0'} .= "$nucl_rules[$irule][$mol_types_used[$imol]], ";
		}
		for($irule=0;$irule<$n_nucl_rules_neg;$irule++){
			$lines{'outn'} .= "$nucl_rules_neg[$irule][$mol_types_used[$imol]], ";
		}
		for($irule=0;$irule<$n_nucl_rules_pos;$irule++){
			$lines{'outp'} .= "$nucl_rules_pos[$irule][$mol_types_used[$imol]], ";
		}
	}

	$lines{'1of1'} =~ s/, $//;
	$lines{'nof1'} =~ s/, $//;
	$lines{'0'} =~ s/, $//;
	$lines{'n'} =~ s/, $//;
	$lines{'p'} =~ s/, $//;
	$lines{'0c'} =~ s/, $//;
	$lines{'nc'} =~ s/, $//;
	$lines{'pc'} =~ s/, $//;
	$lines{'m'} =~ s/, $//;
	$lines{'d'} =~ s/, $//;
	$lines{'dm'} =~ s/, $//;
	$lines{'mon'} =~ s/ $//;
	$lines{'nmon'} =~ s/ $//;
	$lines{'cl'} =~ s/ $//;
	$lines{'t0'} =~ s/ $//;
	$lines{'tn'} =~ s/ $//;
	$lines{'tp'} =~ s/ $//;
	$lines{'out0'} =~ s/, $//;
	$lines{'outn'} =~ s/, $//;
	$lines{'outp'} =~ s/, $//;
	
	if($lines{'out0'} !~ /^\s*$/){
		$lines{'out0'} = "reshape((/$lines{'out0'}/),(/$n_nucl_rules, $nmol_types_used/))";
	}
	if($lines{'outn'} !~ /^\s*$/){
		$lines{'outn'} = "reshape((/$lines{'outn'}/),(/$n_nucl_rules_neg, $nmol_types_used/))";
	}
	if($lines{'outp'} !~ /^\s*$/){
		$lines{'outp'} = "reshape((/$lines{'outp'}/),(/$n_nucl_rules_pos, $nmol_types_used/))";
	}
	
	if($lfortran){
		print SYST "integer, parameter :: n_1$molecule_name[0]_clusters = $i\t\t\t\t! number molecules and clusters containing 1 $molecule_name[0] molecule\n";
		print SYST "\ninteger, parameter :: n_neutrals = $j\t\t\t! number of neutral molecules and clusters\n";
		print SYST "integer, parameter :: n_negatives = $k\t\t\t! number of negative molecules and clusters\n";
		print SYST "integer, parameter :: n_positives = $l\t\t\t! number of positive molecules and clusters\n";
		print SYST "\ninteger, parameter :: n_neutral_clusters = $jc\t\t\t! number of neutral clusters\n";
		print SYST "integer, parameter :: n_negative_clusters = $kc\t\t\t! number of negative clusters\n";
		print SYST "integer, parameter :: n_positive_clusters = $lc\t\t\t! number of positive clusters\n\n";
	
		$clust_mass_max = sprintf("%.2f",$clust_mass_max);
		$clust_mass_max =~ s/e/d/;
		print SYST "real(kind(1.d0)), parameter :: mass_max = $clust_mass_max\n";
		$diameter_max = sprintf("%.2f",$diameter_max);
		$diameter_max =~ s/e/d/;
		print SYST "real(kind(1.d0)), parameter :: diameter_max = $diameter_max\n";
		$mob_diameter_max = sprintf("%.2f",$mob_diameter_max);
		$mob_diameter_max =~ s/e/d/;
		print SYST "real(kind(1.d0)), parameter :: mob_diameter_max = $mob_diameter_max\n";
		
		if($lines{'out0'} !~ /^\s*$/){
			print SYST "\ninteger, parameter :: nmols_out_neutral($n_nucl_rules, $nmol_types_used) = $lines{'out0'}\t\t\t! criteria for outgrowing neutrals\n";
		}
		if($lines{'outn'} !~ /^\s*$/){
			print SYST "integer, parameter :: nmols_out_negative($n_nucl_rules_neg, $nmol_types_used) = $lines{'outn'}\t\t\t! criteria for outgrowing negatives\n";
		}
		if($lines{'outp'} !~ /^\s*$/){
			print SYST "integer, parameter :: nmols_out_positive($n_nucl_rules_pos, $nmol_types_used) = $lines{'outp'}\t\t\t! criteria for outgrowing positives\n";
		}
	
		$lines{'c_inp'} = '';
		$lines{'c_size'} = '';
		$lines{'c_def'} = '';
		$lines{'a_inp'} = '';
		$lines{'a_size'} = '';
		$lines{'a_def'} = '';
		print SYST "\n\ncontains\n\nsubroutine n_$molecule_name[0]_in_clusters(n_$molecule_name[0])\n\timplicit none\n";
		print SYST "\tinteger :: n_$molecule_name[0]($max_cluster_number)\n\n";
		print SYST "\tn_$molecule_name[0] = (/$lines{'nof1'}/)\n\n";
		print SYST "end subroutine n_$molecule_name[0]_in_clusters\n\n";
		print SYST "subroutine clusters_with_1_$molecule_name[0](cluster_numbers)\n\timplicit none\n";
		print SYST "\tinteger :: cluster_numbers($i)\n\n\tcluster_numbers = (/$lines{'1of1'}/)\n\n";
		print SYST "end subroutine clusters_with_1_$molecule_name[0]\n\n";
		if($j>0 && ($k+$l)>0 && $lines{'0'} ne ''){
			$lines{'a_inp'} .= "neutrals, ";
			$lines{'a_size'} .= "neutrals($j), ";
			$lines{'a_def'} .= "\tneutrals = (/$lines{'0'}/)\n";
		}
		if($lines{'0c'} ne ''){
			$lines{'c_inp'} .= "neutral_clusters, ";
			$lines{'c_size'} .= "neutral_clusters($jc), ";
			$lines{'c_def'} .= "\tneutral_clusters = (/$lines{'0c'}/)\n";
		}
	
		if($k>0){
			if(($i+$l)>0 && $lines{'n'} ne ''){
				$lines{'a_inp'} .= "negatives, ";
				$lines{'a_size'} .= "negatives($k), ";
				$lines{'a_def'} .= "\tnegatives = (/$lines{'n'}/)\n";
			}
			if($lines{'nc'} ne ''){
				$lines{'c_inp'} .= "negative_clusters, ";
				$lines{'c_size'} .= "negative_clusters($kc), ";
				$lines{'c_def'} .= "\tnegative_clusters = (/$lines{'nc'}/)\n";
			}
		}
		if($l>0){
			if(($i+$k)>0 && $lines{'p'} ne ''){
				$lines{'a_inp'} .= "positives, ";
				$lines{'a_size'} .= "positives($l), ";
				$lines{'a_def'} .= "\tpositives = (/$lines{'p'}/)\n";
			}
			if($lines{'pc'} ne ''){
				$lines{'c_inp'} .= "positive_clusters, ";
				$lines{'c_size'} .= "positive_clusters($lc), ";
				$lines{'c_def'} .= "\tpositive_clusters = (/$lines{'pc'}/)\n";
			}
		}
		$lines{'a_inp'} =~ s/, $//;
		if($lines{'a_inp'} ne ''){
			print SYST "subroutine arrays$append_to_subroutine_names($lines{'a_inp'})\n\timplicit none\n";
			$lines{'a_size'} =~ s/, $//;
			print SYST "\tinteger :: $lines{'a_size'}\n\n";
			print SYST $lines{'a_def'};
			print SYST "\nend subroutine arrays$append_to_subroutine_names\n\n";
		}
	
		$lines{'c_inp'} =~ s/, $//;
		print SYST "subroutine cluster_arrays$append_to_subroutine_names($lines{'c_inp'})\n\timplicit none\n";
		$lines{'c_size'} =~ s/, $//;
		print SYST "\tinteger :: $lines{'c_size'}\n\n";
		print SYST $lines{'c_def'};
		print SYST "\nend subroutine cluster_arrays$append_to_subroutine_names\n\n";
	
	
		print SYST "subroutine get_mass$append_to_subroutine_names(mass)\n\timplicit none\n";
		print SYST "\treal(kind(1.d0)) :: mass($max_cluster_number)\n\n\tmass = (/$lines{'m'}/)\n";
		print SYST "\nend subroutine get_mass$append_to_subroutine_names\n\n";
	
	
		print SYST "subroutine get_diameter$append_to_subroutine_names(diameter)\n\timplicit none\n\n";
		print SYST "\treal(kind(1.d0)) :: diameter($max_cluster_number)\n\n\t diameter = (/$lines{'d'}/)\t! dry value\n";
		print SYST "\nend subroutine get_diameter$append_to_subroutine_names\n\n";
	
	
		print SYST "subroutine get_mob_diameter$append_to_subroutine_names(mob_diameter)\n\timplicit none\n\n";
		print SYST "\treal(kind(1.d0)) :: mob_diameter($max_cluster_number)\n\n\t mob_diameter = (/$lines{'dm'}/)\t! dry value\n";
	
		print SYST "\nend subroutine get_mob_diameter$append_to_subroutine_names\n\n";
	
		print SYST "subroutine cluster_names$append_to_subroutine_names(clust)\n\timplicit none\n\tcharacter(len=11), dimension($size_C_out) :: clust\n\n";
		for($iclus=1;$iclus <= $size_C_out;$iclus++){
			if ($lhydrate_average && $lhas_hydr[$iclus]==1){
				$line = "\t! $hydrates_str[$iclus]";
			}else{
				$line = '';
			}
			print SYST "\tclust($iclus)(:) = '$label[$iclus]'$line\n";
		}
		print SYST "\nend subroutine cluster_names$append_to_subroutine_names\n\n";

	}else{

		# print the cluster names and some info on the clusters

		print DRIVER "\n\n% Cluster definitions etc.\n\n";

		$lines{'nof1'} =~ s/&\n\t\t&/...\n\t\t/g;
		$lines{'0'} =~ s/&\n\t\t&/...\n\t\t/g;
		$lines{'n'} =~ s/&\n\t\t&/...\n\t\t/g;
		$lines{'p'} =~ s/&\n\t\t&/...\n\t\t/g;
		$lines{'0c'} =~ s/&\n\t\t&/...\n\t\t/g;
		$lines{'nc'} =~ s/&\n\t\t&/...\n\t\t/g;
		$lines{'pc'} =~ s/&\n\t\t&/...\n\t\t/g;
		$lines{'m'} =~ s/&\n\t\t&/...\n\t\t/g;
		$lines{'d'} =~ s/&\n\t\t&/...\n\t\t/g;
		$lines{'dm'} =~ s/&\n\t\t&/...\n\t\t/g;
	
		$lines{'nof1'} =~ s/, / /g;
		$lines{'0'} =~ s/, / /g;
		$lines{'n'} =~ s/, / /g;
		$lines{'p'} =~ s/, / /g;
		$lines{'m'} =~ s/, / /g;
		$lines{'d'} =~ s/, / /g;
		$lines{'dm'} =~ s/, / /g;
		$lines{'0c'} =~ s/, / /g;
		$lines{'nc'} =~ s/, / /g;
		$lines{'pc'} =~ s/, / /g;
		$lines{'out0'} =~ s/, / /g;
		$lines{'outn'} =~ s/, / /g;
		$lines{'outp'} =~ s/, / /g;
		
		$line = "\t% molecule names\n\tlabels = {";
		for($imol=1;$imol<=$nmol_types_used;$imol++){
			$line = $line . "'" . $mol_names_used[$imol] . "' ";
		}
		$line =~ s/ $//;
		$line = "$line};\n";
		print DRIVER $line;
		
		print DRIVER "\t% cluster names\n\tclust = {$lines{'cl'}};\n";
		# print out also the cluster names plus the names of the other fluxes (sources, losses etc.)
		$line = "\t% name vector containing also the names of fluxes other than clusters (sources, losses etc.)\n\tclust_flux = [clust ";
		for($iclus=$max_cluster_number+1;$iclus<=$max_n_flux_inc_boundary_clusters;$iclus++){
			$line = $line . "'" . $label_inc_boundary_clusters[$iclus] . "' ";
		}
		$line =~ s/ $//;
		$line = "$line];\n";
		print DRIVER $line;
		
		print DRIVER "\t% dry values\n\tdiameters = [$lines{'d'}];\n";
		print DRIVER "\t% dry values\n\tmobility_diameters = [$lines{'dm'}]; % mobility diameter (2*r+0.3nm)*sqrt(1+28.8u/m)\n";
		print DRIVER "\tmasses = [$lines{'m'}]; % mass in amus\n";
		print DRIVER "\tmonomers = [$lines{'mon'}];\n";
		print DRIVER "\tnonmonomers = [$lines{'nmon'}];\n";
		print DRIVER "\tneutrals = [$lines{'0'}];\n";
		print DRIVER "\tnegatives = [$lines{'n'}];\n";
		print DRIVER "\tpositives = [$lines{'p'}];\n";
		print DRIVER "\tneutral_clusters = [$lines{'0c'}];\n";
		print DRIVER "\tnegative_clusters = [$lines{'nc'}];\n";
		print DRIVER "\tpositive_clusters = [$lines{'pc'}];\n";
		print DRIVER "\t" . $molecule_name[0] . "_in_clust = [$lines{'nof1'}];\n";
		print DRIVER "\ttable_0 = {$lines{'t0'}};\n";
		print DRIVER "\ttable_neg = {$lines{'tn'}};\n";
		print DRIVER "\ttable_pos = {$lines{'tp'}};\n";
		
		if($lines{'out0'} !~ /^\s*$/){
			$lines{'out0'} =~ s/\(\//\[/g;
			$lines{'out0'} =~ s/\/\)/\]/g;
			print DRIVER "\t% criteria for outgrowing neutrals\n\tnmols_out_neutral = $lines{'out0'};\n";
		}
		if($lines{'outn'} !~ /^\s*$/){
			$lines{'outn'} =~ s/\(\//\[/g;
			$lines{'outn'} =~ s/\/\)/\]/g;
			print DRIVER "\t% criteria for outgrowing negatives\n\tnmols_out_negative = $lines{'outn'};\n";
		}
		if($lines{'outp'} !~ /^\s*$/){
			$lines{'outp'} =~ s/\(\//\[/g;
			$lines{'outp'} =~ s/\/\)/\]/g;
			print DRIVER "\t% criteria for outgrowing positives\n\tnmols_out_positive = $lines{'outp'};\n";
		}
		
	}
}else{
	$str_temp = '';
	for($imol=1;$imol<=$nmol_types_used;$imol++){
		$str_temp .= "$n_max_type[$mol_types_used[$imol]]$mol_names_used[$imol]";
	}
	($string1,$string2) = &determine_cluster_composition($str_temp);
	@temp_array=@$string1;
	$diameter_max = &calculate_radius(\@temp_array)*2e9;
	$diameter_max = sprintf("%.2f",$diameter_max);
	$diameter_max =~ s/e/d/;
	print SYST "real(kind(1.d0)), parameter :: diameter_max = $diameter_max\n";

	print SYST "\ncontains\n\n";
}

if($lfortran){
	print SYST "subroutine molecule_names$append_to_subroutine_names(labels)\n\timplicit none\n";
	print SYST "\tcharacter(len=11), dimension($nmol_types_used) :: labels\n\n";
	for($imol=1;$imol<=$nmol_types_used;$imol++){
		print SYST "\tlabels($imol)(:) = '$mol_names_used[$imol]'\n";
	}
	print SYST "\nend subroutine molecule_names$append_to_subroutine_names\n\n";
	
	print SYST "subroutine monomer_indices$append_to_subroutine_names(n_monomers)\n\timplicit none\n";
	print SYST "\tinteger :: n_monomers($nmol_types_used)\n\n";
	print SYST "\tn_monomers = (/$monomer_indices/)\n\n";
	print SYST "end subroutine monomer_indices$append_to_subroutine_names\n\n";
	print SYST "\nend module acdc_system\n\n";
	close(SYST);
}


################################################################################
################################################################################
##      8.2 dC/dt for Fortran                                                 ##
################################################################################
################################################################################

if($lfortran){
	# Opening the Fortran output file
	open(OUT, ">$output_file_name_fortran") || die "\nCould not open output file $output_file_name_fortran\n\n";
	
	if($lloop){
		print OUT "!!!!!!!!!!!!!!!!!!!!!!!! KEY !!!!!!!!!!!!!!!!!!!!!!!!\n!\n";
		if($nmol_types==1){
			print OUT "! $max_cluster clusters with 1-$n_max_type[0] $molecule_name[0] molecules\n";
		}else{
			print OUT "! $max_cluster clusters with\n";
			for($imol=0;$imol<$nmol_types;$imol++){
				if($n_max_type[$imol]>0){
					print OUT "!   0-$n_max_type[$imol] $molecule_name[$imol] molecules\n";
				}
			}
			if($lloop_max_ratio){
				print OUT "! so that the maximum number of molecules $molecule_name[$nmol_inp_not_basis] per each molecule $molecule_name[$nmol_inp_max_ratio_basis] is $max_wrt_mon";
				if($lfree_mon){
					print OUT ",\n! except for the monomer 1$molecule_name[$nmol_inp_max_ratio_basis] which is free of $molecule_name[$nmol_inp_not_basis]";
				}
				print OUT "\n";
			}
			if($l_save_outgoing){
				print OUT "! and element $max_n_flux is the outgoing flux\n";
			}
		}
		if($l_j_in){
			print OUT "!\n! using a source for cluster n_j_in and omitting the distribution below it\n";
		}
		print OUT "!\n!!!!!!!!!!!!!!!!!!!!!! END KEY !!!!!!!!!!!!!!!!!!!!!!\n\n";
	}else{
		print OUT "!!!!!!!!!!!!!!!!!!!!!!!! KEY !!!!!!!!!!!!!!!!!!!!!!!!\n";
		for($iclus=1;$iclus<=$max_n_flux;$iclus++){
			$line = $comments[$iclus];
			$line =~ s/%/!/;
			print OUT "$line\n";
		}
		print OUT "!!!!!!!!!!!!!!!!!!!!!! END KEY !!!!!!!!!!!!!!!!!!!!!!\n\n";
		
		#maximum number of columns in the index tables
		$ind_quad_loss[0] = 0;
		$ind_quad_form[0] = 0;
		$ind_quad_loss_extra[0] = 0;
		$ind_quad_form_extra[0] = 0;
		$ind_lin_loss[0] = 0;
		$ind_lin_form[0] = 0;
		$ind_reac_loss[0] = 0;
		$ind_reac_form[0] = 0;
		for($iclus=1;$iclus <= $max_n_flux;$iclus++){
			if($ind_quad_loss[$iclus][0]*2>$ind_quad_loss[0]){
				$ind_quad_loss[0]=$ind_quad_loss[$iclus][0]*2;
			}
			if($ind_quad_form[$iclus][0]*2>$ind_quad_form[0]){
				$ind_quad_form[0]=$ind_quad_form[$iclus][0]*2;
			}
			if($ind_quad_loss_extra[$iclus][0]*3>$ind_quad_loss_extra[0]){
				$ind_quad_loss_extra[0]=$ind_quad_loss_extra[$iclus][0]*3;
			}
			if($ind_quad_form_extra[$iclus][0]*3>$ind_quad_form_extra[0]){
				$ind_quad_form_extra[0]=$ind_quad_form_extra[$iclus][0]*3;
			}
			if($ind_lin_loss[$iclus][0]*2>$ind_lin_loss[0]){
				$ind_lin_loss[0]=$ind_lin_loss[$iclus][0]*2;
			}
			if($ind_lin_form[$iclus][0]*2>$ind_lin_form[0]){
				$ind_lin_form[0]=$ind_lin_form[$iclus][0]*2;
			}
			if($ind_reac_loss[$iclus][0]*2>$ind_reac_loss[0]){
				$ind_reac_loss[0]=$ind_reac_loss[$iclus][0]*2;
			}
			if($ind_reac_form[$iclus][0]*2>$ind_reac_form[0]){
				$ind_reac_form[0]=$ind_reac_form[$iclus][0]*2;
			}
		}
	}
	
	print OUT "! differential equations: f = dc/dt\n";
	if($lloop && !$lloop_cs && !$l_save_outgoing){
		print OUT "subroutine feval$append_to_subroutine_names(neqn,t,c,f_out,coef,ipar)\n";
	}else{
		print OUT "subroutine feval$append_to_subroutine_names(neqn,t,c,f,coef,ipar)\n";
	}
	if($lloop){
		print OUT "use monomer_settings, only : sources_and_constants\n";
		if($l_j_in){
			print OUT "use shared_input, only : r_j_in => rlim_for_J";
			print "\nUsing module shared_input in the feval subroutine to get r_j_in\n";
			if($lloop_j_tests){
				print OUT ", r_j_init => rlim_for_J_init\n";
				print OUT "use shared_input, only : use_jtable, use_KK";
			}
			print OUT "\n";
		}
	}
	print OUT "\timplicit none\n";
	if($lloop_cs && $nmol_types_used==1){
		print OUT "\tinteger, parameter :: nclust = $max_cluster_number, n_cs = $max_cluster_number\n";
	}elsif($lloop && $nmol_types_used>1){
		print OUT "\tinteger, parameter :: nclust = $max_cluster_number, ij_ind_max($nmol_types_used) = (/$nmol_ranges_list/)";
		if($disable_nonmonomers){
			print OUT ", n_monomers($nmol_types_used) = (/$monomer_indices/)";
		}
		print OUT "\n";
	}else{
		print OUT "\tinteger, parameter :: nclust = $max_cluster_number\n";
	}
	if(!$lloop){
		print OUT "\tinteger :: neqn, ipar(4), i, j, k, n\n";
	}else{
		print OUT "\tinteger :: neqn, ipar(4), i, j, ij, n";
		if($nmol_types_used>1){
			print OUT ", ij_ind($nmol_types_used), imol";
		}
		print OUT "\n";
	}
	if($l_j_in){
		print OUT "\tinteger, save :: n_j_in = 0";
		if($lloop_j_tests){
			print OUT ", n_j_init = 0";
		}
		print OUT "\n";
	}
	if($lloop && !$lloop_cs && !$l_save_outgoing){
		print OUT "\treal(kind(1.d0)) :: t,c(neqn),f_out(neqn)";
	}else{
		print OUT "\treal(kind(1.d0)) :: t,c(neqn),f(neqn)";
	}
	if($l_mcmc_any){
		print OUT ",coef(" . ($n_mcmc_coefs+$variable_cs+2*$variable_ion_source) . ")\n";
	}elsif($variable_temp || $variable_cs || $variable_ion_source){
		print OUT ",coef(" . ($variable_temp+$variable_cs+2*$variable_ion_source) . ")\n";
	}elsif($lloop && $nmol_types_used>1){
		print OUT ",coef($nmol_types_used),pfac\n";
	}elsif($lloop){
		print OUT ",coef,pfac";
		if($lloop_j_tests || $lgrowth_for_j_tests || $l_j_in){
			print OUT ",m,r";
			if($lloop_j_tests && $l_j_in){
				print OUT ",cs_for_table,gr,coags,j_temp";
			}
		}
		print OUT "\n";
	}else{
		print OUT ",coef(:)\n";
	}
	if($lloop && !$lloop_cs && !$l_save_outgoing){
		print OUT "\treal(kind(1.d0)) :: f(neqn+1)\n";
	}
	if($lloop){
		print OUT "\tlogical, save :: isconst($size_C_out) = .false.\n";
		print OUT "\treal(kind(1.d0)), save :: K(nclust,nclust)";
		if(!$disable_evap){
			print OUT ", E(";
			if($disable_nonmonomer_evaps && $nmol_types_used==1){
				print OUT "1";
			}else{
				print OUT "nclust";
			}
			print OUT ",nclust)";
		}
		if($l_use_get_losses){
			print OUT ", loss(nclust)";
		}
		print OUT ", source($size_C_out)\n";
		print OUT "\tinteger, save :: fitted(nclust,0:nclust)=0";
		if(!$lloop_diag && $nmol_types_used>1){
			print OUT ", indices(nclust,$nmol_types_used)";
			if($lloop_max_ratio){
				print OUT ", n_pure($nmol_ranges_array[$nmol_max_ratio_basis])";
			}
		}
	}else{
		print OUT "\tlogical, save :: isconst($size_C_out) = .false.\n";
		print OUT "\treal(kind(1.d0)), save :: coef_quad(nclust,nclust,$max_n_flux)=0.d0,coef_lin($max_n_flux,$max_n_flux,nclust)=0.d0,source($size_C_out)=0.d0\n";
		print OUT "\tinteger, save :: ind_quad_loss($size_C_out,0:$ind_quad_loss[0])=0,ind_quad_form($size_C_out,0:$ind_quad_form[0])=0\n";
		print OUT "\tinteger, save :: ind_quad_loss_extra($size_C_out,0:$ind_quad_loss_extra[0])=0,ind_quad_form_extra($size_C_out,0:$ind_quad_form_extra[0])=0\n";
		print OUT "\tinteger, save :: ind_lin_loss($size_C_out,0:$ind_lin_loss[0])=0,ind_lin_form($size_C_out,0:$ind_lin_form[0])=0,fitted($max_cluster_number,0:$max_cluster_number)=0\n";
		if($l_reac_in_clus){
			print OUT "\tinteger, save :: ind_reac_loss(nclust,0:$ind_reac_loss[0])=0,ind_reac_form(nclust,0:$ind_reac_form[0])=0\n";
		}
		print OUT "\treal(kind(1.d0)), save :: n_quad_form_extra($size_C_out,0:" . ($ind_quad_form_extra[0])/3 . ")=0.d0";
	}
	print OUT $lines{'source_def'};
	print OUT "\n\treal(kind(1.d0)), parameter :: pi=4.d0*atan(1.d0)\n";
	print OUT "\n\n\t! the parameters are read at the very beginning and eg. if source terms or collision rates have changed\n";
	print OUT "\tif (ipar(1) .eq. 0) then\n";
	print OUT "\t\t! after this the same values are used until some other routine tells otherwise\n\t\tipar(1) = 1\n";
	if($lloop){
		print OUT "\t\tcall sources_and_constants(source,isconst,fitted)\n";
		if($l_mcmc_any){
			print OUT "\t\tcall get_coll$append_to_subroutine_names(K,coef)\n";
		}elsif( $variable_temp){
			print OUT "\t\tcall get_coll$append_to_subroutine_names(K,coef(1))\n";
		}else{
			print OUT "\t\tcall get_coll$append_to_subroutine_names(K)\n";
		}
		if(!$disable_evap){
			if($l_mcmc_any){
				print OUT "\t\tcall get_evap$append_to_subroutine_names(E,coef)\n";
			}elsif( $variable_temp){
				print OUT "\t\tcall get_evap$append_to_subroutine_names(E,K,coef(1))\n";
			}else{
				print OUT "\t\tcall get_evap$append_to_subroutine_names(E)\n";
			}
		}
		if($l_use_get_losses){
			print OUT "\t\tcall get_losses$append_to_subroutine_names(loss)\n";
		}
		if(!$lloop_diag && $nmol_types_used>1){
			print OUT "\t\tcall get_molecule_numbers$append_to_subroutine_names(indices)\n";
			if($lloop_max_ratio){
				print OUT "\t\tcall pure_clust_indices$append_to_subroutine_names(n_pure)\n";
			}
		}		
		if($l_j_in){
			print OUT "\t\tdo i=1,nclust\n";
			print OUT "\t\t\tcall get_masses_and_radii$append_to_subroutine_names(m,r,i)\n";
			print OUT "\t\t\tif (r*1.d9 .ge. r_j_in) then\n";
			print OUT "\t\t\t\tn_j_in = i\n";
			print OUT "\t\t\t\twrite(*,*) 'Size for which the input J will be used:'\n";
			print OUT "\t\t\t\twrite(*,*) r_j_in,' nm, ',n_j_in,' molecules'\n";
			print OUT "\t\t\t\texit\n";
			print OUT "\t\t\tend if\n";
			print OUT "\t\tend do\n";
			if($lloop_j_tests){
				print OUT "\t\tif (use_KK) then\n";
				print OUT "\t\t\tdo i=1,nclust\n";
				print OUT "\t\t\t\tcall get_masses_and_radii(m,r,i)\n";
				print OUT "\t\t\t\tif (r*1.d9 .ge. r_j_init) then\n";
				print OUT "\t\t\t\t\tn_j_init = i\n";
				print OUT "\t\t\t\t\twrite(*,*) 'Size from which the input J will be extrapolated:'\n";
				print OUT "\t\t\t\t\twrite(*,*) r_j_init,' nm, ',n_j_init,' molecules'\n";
				print OUT "\t\t\t\t\texit\n";
				print OUT "\t\t\t\tend if\n";
				print OUT "\t\t\tend do\n";
				print OUT "\t\tend if\n";
			}
		}
	}else{
		print OUT "\t\tcall initialize_parameters$append_to_subroutine_names(coef_quad,coef_lin,ind_quad_loss,ind_quad_form,&\n\t\t&\tind_lin_loss,ind_lin_form,ind_quad_loss_extra,ind_quad_form_extra,&\n\t\t&\tsource,isconst,fitted";
		if($l_mcmc_any || $variable_temp || $variable_cs || $variable_ion_source){
			print OUT ",coef";
		}
		print OUT ",ipar)\n";
		if($ind_quad_form_extra[0]>0){
			print OUT "\t\tn_quad_form_extra(:,1:" . ($ind_quad_form_extra[0]/3) . ") = real(ind_quad_form_extra(:,3:$ind_quad_form_extra[0]:3),kind=kind(1.d0))\n";
		}
	}
	print OUT "\tend if\n";
	if($lloop && !$lloop_diag){
		print OUT "\n\t! override the monomer sources given by the subroutine call\n";
		print OUT "\tsource((/$monomer_indices/)) = coef\n";
		print OUT $lines{'source'};
		print OUT "\n\t! override the isconst settings, if needed\n";
		print OUT "\tif (ipar(4) .eq. 1) then\n\t\tisconst = .false.\n\tend if\n\n";
		if($l_use_get_losses){
			print OUT "\tf = (/source-loss*c,0.d0/)\n\n";
		}else{
			print OUT "\tf = (/source,0.d0/)\n\n";
		}
		if($l_j_in){
			# Assume that the formation at size n_j_in directly reduces the vapor concentration,
			# as is usually done in growth models
			print OUT "\tf(1) = f(1)-n_j_in*source(n_j_in)\n\n";
		}
		# For a one-component system with only monomer evaporation, print the evaporation terms here for a neater output;
		# otherwise add them into the loop after collision terms
		if(!$disable_evap && $disable_nonmonomer_evaps && $nmol_types_used == 1){
			if(!$l_j_in){
				print OUT "\tf(1) = f(1)+2.d0*E(1,1)*c(2)\n";
				if($lloop_cs){
					print OUT "\tf(1) = f(1)+sum(E(1,2:nclust-2)*c(3:nclust-1))\n";
					print OUT "\tf(2:nclust-2) = f(2:nclust-2)+E(1,2:nclust-2)*c(3:nclust-1)\n";
					print OUT "\tf(2:nclust-1) = f(2:nclust-1)-E(1,1:nclust-2)*c(2:nclust-1)\n";
				}else{
					print OUT "\tf(1) = f(1)+sum(E(1,2:nclust-1)*c(3:nclust))\n";
					print OUT "\tf(2:nclust-1) = f(2:nclust-1)+E(1,2:nclust-1)*c(3:nclust)\n";
					print OUT "\tf(2:nclust) = f(2:nclust)-E(1,1:nclust-1)*c(2:nclust)\n";
				}
			}else{
				if($lloop_cs){
					print OUT "\tf(1) = f(1)+sum(E(1,n_j_in:nclust-2)*c(n_j_in+1:nclust-1))\n";
					print OUT "\tf(n_j_in:nclust-2) = f(n_j_in:nclust-2)+E(1,n_j_in:nclust-2)*c(n_j_in+1:nclust-1)\n";
					print OUT "\tf(n_j_in+1:nclust-1) = f(n_j_in+1:nclust-1)-E(1,n_j_in:nclust-2)*c(n_j_in+1:nclust-1)\n";
				}else{
					print OUT "\tf(1) = f(1)+sum(E(1,n_j_in:nclust-1)*c(n_j_in+1:nclust))\n";
					print OUT "\tf(n_j_in:nclust-1) = f(n_j_in:nclust-1)+E(1,n_j_in:nclust-1)*c(n_j_in+1:nclust)\n";
					print OUT "\tf(n_j_in+1:nclust) = f(n_j_in+1:nclust)-E(1,n_j_in:nclust-1)*c(n_j_in+1:nclust)\n";
				}
			}
		}
		if(!$l_j_in){
			print OUT "\tdo i=1,nclust\n";
		}else{
			print OUT "\tdo i=n_j_in,nclust\n";
			print OUT "\t\t! collisions among clusters taken into account in the distribution\n";
		}
		if($lloop_max_ratio && $lfree_mon){
			print OUT "\t\tif ((indices(i,$nmol_max_ratio_basis) .eq. 1) .and. (indices(i,$nmol_not_basis) .gt. 0)) then\n";
			print OUT "\t\t\tif (c(i) .gt. 1.d-10) then\n";
			print OUT "\t\t\t\twrite(*,*)'Cluster ',i,' has a non-zero concentration of ',c(i)*1.d-6,' cm^-3?'\n";
			print OUT "\t\t\tend if\n";
			print OUT "\t\t\tf(i) = 0.d0\n";
			print OUT "\t\t\tcycle\n";
			print OUT "\t\tend if\n";
		}
		if($disable_nonmonomers){
			if($nmol_types_used>1){
				print OUT "\t\tdo imol = 1,$nmol_types_used\n";
				print OUT "\t\t\tj = n_monomers(imol)\n";
				# Don't take the same process twice
				print OUT "\t\t\tif ((sum(indices(i,:)).eq.1) .and. (j.lt.i)) then\n";
				print OUT "\t\t\t\tcycle\n";
				print OUT "\t\t\tend if\n";
			}else{
				print OUT "\t\tdo j=1,1\n";
			}
		}else{
			print OUT "\t\tdo j=i,nclust\n";
		}
		if($nmol_types_used>1){
			print OUT "\t\t\tij_ind = indices(i,:)+indices(j,:)\n";
			if($lloop_max_ratio){
				# skip the collisions resulting in a wrong composition (i.e. those of the "non-basis" molecule with clusters that already have the maximum number of these molecules)
				if($lfree_mon){
					print OUT "\t\t\tif (((ij_ind($nmol_max_ratio_basis) .eq. 1) .and. (ij_ind($nmol_not_basis) .gt. 0)) .or. (ij_ind($nmol_not_basis) .gt. $max_wrt_mon*ij_ind($nmol_max_ratio_basis))) then\n";
				}else{
					print OUT "\t\t\tif (ij_ind($nmol_not_basis) .gt. $max_wrt_mon*ij_ind($nmol_max_ratio_basis)) then\n";
				}
				print OUT "\t\t\t\tcycle\n";
				print OUT "\t\t\tend if\n";
			}
			if($lloop_cs){
				print OUT "\t\t\tij_ind = min(ij_ind,ij_ind_max)\n";
				print OUT "\t\t\tij = $cluster_from_indices\n";
			}else{
				print OUT "\t\t\tif (any(ij_ind .gt. ij_ind_max)) then\n";
				print OUT "\t\t\t\tij = neqn+1\n";
				print OUT "\t\t\telse\n";
				print OUT "\t\t\t\tij = $cluster_from_indices\n";
				print OUT "\t\t\tend if\n";
			}
		}else{
			print OUT "\t\t\tij = min($n_flux_out,i+j)\n";
		}
		# Add the factor that deals with the permutations for processes involving two similar clusters
		$lines{'pfac'} = "\tpfac = 1.d0\n";
		$lines{'pfac'} .= "\tif (i .eq. j) then\n";
		$lines{'pfac'} .= "\t\tpfac = 0.5d0\n";
		$lines{'pfac'} .= "\tend if\n";		
		$str_temp = $lines{'pfac'};
		$str_temp =~ s/\t(?!\t)/\t\t\t/g;
		print OUT $str_temp;
		print OUT "\t\t\tf(i) = f(i)-pfac*K(i,j)*c(i)*c(j)\n";
		print OUT "\t\t\tf(j) = f(j)-pfac*K(i,j)*c(i)*c(j)\n";
		print OUT "\t\t\tf(ij) = f(ij)+pfac*K(i,j)*c(i)*c(j)\n";
		# Create a string containing the condition for evaporation of cluster i+j (if (...) then)
		$lines{'l_evap'} = "if (";
		if($nmol_types_used>1){
			if($disable_nonmonomer_evaps){
				$lines{'l_evap'} .= "min(sum(indices(i,:)),sum(indices(j,:))).eq.1 .and. ";					
			}
			if($lloop_cs){
				$lines{'l_evap'} .= "all(ij_ind.lt.ij_ind_max)";
			}else{
				$lines{'l_evap'} .= "all(ij_ind.le.ij_ind_max)";
			}
		}else{
			if($disable_nonmonomer_evaps){
				$lines{'l_evap'} .= "min(i,j).eq.1 .and. ";					
			}
			$lines{'l_evap'} .= "i+j.lt.$n_flux_out";
		}
		$lines{'l_evap'} .= ") then\n";
		if(!$disable_evap && (!$disable_nonmonomer_evaps || $nmol_types_used>1)){
			print OUT "\t\t\t$lines{'l_evap'}";
			print OUT "\t\t\t\tf(i) = f(i)+E(i,j)*c(ij)\n";
			print OUT "\t\t\t\tf(j) = f(j)+E(i,j)*c(ij)\n";
			print OUT "\t\t\t\tf(ij) = f(ij)-E(i,j)*c(ij)\n";
			print OUT "\t\t\tend if\n";
		}
		print OUT "\t\tend do\n";
		if($l_j_in && !$disable_nonmonomers){
			# Deal with the clusters colliding with a monomer here, since the monomer isn't included in the equation loop
			# when all collisions are taken into account
			print OUT "\t\t! monomer collisions\n";
			print OUT "\t\tj = 1\n";
			print OUT "\t\tij = min($n_flux_out,i+j)\n";
			print OUT "\t\tf(i) = f(i)-K(i,j)*c(i)*c(j)\n";
			print OUT "\t\tf(j) = f(j)-K(i,j)*c(i)*c(j)\n";
			print OUT "\t\tf(ij) = f(ij)+K(i,j)*c(i)*c(j)\n";
			if(!$disable_evap && !$disable_nonmonomer_evaps){
				print OUT "\t\t$lines{'l_evap'}";
				print OUT "\t\t\tf(i) = f(i)+E(i,j)*c(ij)\n";
				print OUT "\t\t\tf(j) = f(j)+E(i,j)*c(ij)\n";
				print OUT "\t\t\tf(ij) = f(ij)-E(i,j)*c(ij)\n";
				print OUT "\t\tend if\n";
			}
		}
		if($lgrowth_for_j_tests){
			print OUT "\t\t! collisions with the implicit growth compound\n";
			print OUT "\t\tcall get_masses_and_radii(m,r,i)\n";
			print OUT "\t\tif (r .ge. $rlim_growth_compound_str) then\n";
			print OUT "\t\t\tj = 2\n";
			print OUT "\t\t\tij = min($n_flux_out,i+j)\n";
			print OUT "\t\t\tf(i) = f(i)-K(i,j)*c(i)*$c_growth_compound_str\n";
			print OUT "\t\t\tf(ij) = f(ij)+K(i,j)*c(i)*$c_growth_compound_str\n";
			print OUT "\n";
			print OUT "\t\t\tj = 3\n";
			print OUT "\t\t\tij = min($n_flux_out,i+j)\n";
			print OUT "\t\t\tf(i) = f(i)-K(i,j)*c(i)*$c_growth_compound_str\n";
			print OUT "\t\t\tf(ij) = f(ij)+K(i,j)*c(i)*$c_growth_compound_str\n";
			print OUT "\t\tend if\n";
		}
		print OUT "\tend do\n";
		print OUT "\tdo i=1,nclust\n";
		print OUT "\t\tif (isconst(i)) then\n";
		print OUT "\t\t\tf(i) = 0.d0\n";
		print OUT "\t\tend if\n";
		print OUT "\tend do\n";
		if(!$lloop_cs && !$l_save_outgoing){
			print OUT "\tf_out = f(1:nclust)\n";
		}
	}else{
		print OUT "\n\tf = 0.d0\n";
		print OUT "\tf = f + source ! add source term\n";
		if($include_generic_neg && $include_generic_pos){
			if($generic_positive_ions_corr ne ''){
				print OUT "\n\t! setting positive charger ion concentration to give charge balance\n";
				print OUT "\tif (c($n_generic_neg)+$generic_positive_ions_corr>0) then\n";
				print OUT "\tc($n_generic_pos) = c($n_generic_neg)+$generic_positive_ions_corr\n";
				print OUT "\telse\n";
				print OUT "\t\tc($n_generic_neg) = -($generic_positive_ions_corr)\n";
				print OUT "\t\tc($n_generic_pos) = 0\n";
				print OUT "\tend if\n\n";
			}elsif($generic_negative_ions_corr ne ''){
				print OUT "\n\t! setting negative charger ion concentration to give charge balance\n";
				print OUT "\tif (c($n_generic_pos)+$generic_negative_ions_corr>0) then\n";
				print OUT "\tc($n_generic_neg) = c($n_generic_pos)+$generic_negative_ions_corr\n";
				print OUT "\telse\n";
				print OUT "\t\tc($n_generic_pos) = -($generic_negative_ions_corr)\n";
				print OUT "\t\tc($n_generic_neg) = 0\n";
				print OUT "\tend if\n\n";
			}elsif($generic_positive_ions ne ''){
				print OUT "\n\t! setting positive charger ion concentration to give charge balance\n";
				print OUT "\tc($n_generic_pos) = $generic_positive_ions\n\n";
			}elsif($generic_negative_ions ne ''){
				print OUT "\n\t! setting negative charger ion concentration to give charge balance\n";
				print OUT "\tc($n_generic_neg) = $generic_negative_ions\n\n";
			}
			if($l_q_neg_mcmc){
				print OUT "\n\t! negative charger ion source term from mcmc coefficient\n";
				print OUT "\tsource($n_generic_neg) = coef(" . ($ncoef_coll+$ncoef_evap+1) . ")\t\t\t\t! Q_{neg}\n\n";
			}
			if($l_q_pos_mcmc){
				print OUT "\n\t! positive charger ion source term from mcmc coefficient\n";
				print OUT "\tsource($n_generic_pos) = coef(" . ($ncoef_coll+$ncoef_evap+$l_q_neg_mcmc+1) . ")\t\t\t\t! Q_{pos}\n\n";
			}
			print OUT "\tdo i=1, neqn ! loop over all the clusters + generic charger ions";
			if($l_save_outgoing){
				print OUT " + other fluxes\n";
			}else{
				print OUT "\n";
			}
		}else{
			print OUT "\tdo i=1, neqn ! loop over all the clusters";
			if($l_save_outgoing){
				print OUT " + other fluxes\n";
			}else{
				print OUT "\n";
			}
		}
		print OUT "\t\tif (.not. isconst(i)) then\n";
		print OUT "\t\t\t! first calculate coefficients for all loss terms\n";
		print OUT "\t\t\tdo n=1, 2*ind_quad_loss(i,0)-1, 2 ! loop over all collisions removing this cluster\n";
		print OUT "\t\t\t\tf(i) = f(i)-coef_quad(i,ind_quad_loss(i,n),ind_quad_loss(i,n+1))*c(ind_quad_loss(i,n))\n";
		print OUT "\t\t\tend do\n";
		print OUT "\t\t\tdo n=1, 2*ind_lin_loss(i,0)-1, 2 ! loop over all evaporations + wall losses etc. removing this cluster\n";
		print OUT "\t\t\t\tf(i) = f(i)-coef_lin(ind_lin_loss(i,n),ind_lin_loss(i,n+1),i)\n";
		print OUT "\t\t\tend do\n";
		print OUT "\t\t\tf(i) = f(i)*c(i) ! multiplying loss coefficients with current concentration\n";
		print OUT "\t\t\t! then add all terms that form this cluster\n";
		print OUT "\t\t\tdo n=1, 2*ind_quad_form(i,0)-1, 2 ! loop over all collisions forming this cluster\n";
		print OUT "\t\t\t\tf(i) = f(i)+coef_quad(ind_quad_form(i,n),ind_quad_form(i,n+1),i)*&\n";
		print OUT "\t\t\t\t&      c(ind_quad_form(i,n))*c(ind_quad_form(i,n+1))\n";
		print OUT "\t\t\tend do\n";
		print OUT "\t\t\tdo n=1, 3*ind_quad_form_extra(i,0)-2, 3 ! loop over all collisions forming this cluster as an extra product\n";
		print OUT "\t\t\t\tf(i) = f(i)+coef_quad(ind_quad_form_extra(i,n),ind_quad_form_extra(i,n+1),i)*&\n";
		print OUT "\t\t\t\t&      n_quad_form_extra(i,(n+2)/3)*c(ind_quad_form_extra(i,n))*c(ind_quad_form_extra(i,n+1))\n";
		print OUT "\t\t\tend do\n";
		print OUT "\t\t\tdo n=1, 2*ind_lin_form(i,0)-1, 2 ! loop over all evaporations forming this cluster\n";
		print OUT "\t\t\t\tf(i) = f(i)+coef_lin(i,ind_lin_form(i,n),ind_lin_form(i,n+1))*c(ind_lin_form(i,n+1))\n";
		print OUT "\t\t\tend do\n";
		print OUT "\t\tend if\n";
		print OUT "\tend do\n";
	}
	print OUT "\tdo n=1, fitted(1,0) ! loop over the clusters whose concentrations are fitted (e.g. [1A]+[1A1D]=const.)\n";
	if($l_j_in){
		print OUT "\t\tif ((fitted(n+1,0).ne.1) .and. (fitted(n+1,0).lt.n_j_in)) then\n";
		print OUT "\t\t\tcycle\n";
		print OUT "\t\tend if\n";
	}
	print OUT "\t\tf(fitted(n+1,0)) = 0.d0 ! don't use value from birth-death equations\n";
	print OUT "\t\tdo i=1, fitted(n+1,1) ! loop over the clusters used for the fitting\n";
	print OUT "\t\t\tf(fitted(n+1,0)) = f(fitted(n+1,0))-f(fitted(n+1,i+1)) ! (e.g. d[1A]/dt = -d[1A1D]/dt)\n";
	print OUT "\t\tend do\n";
	print OUT "\tend do\n";
	print OUT "\nend subroutine feval$append_to_subroutine_names\n\n!-----------------------------------------------------------\n\n";
	
	
	
	print OUT "! jacobian of the differential equations: dfdc(i,j) = df(i)/dc(j)\n! not using this, since the solver works faster using a numerical jacobian...\n";
	print OUT "subroutine jeval$append_to_subroutine_names(neqn,t,c,ml,mu,dfdc,ldim,coef,ipar)\n";
	print OUT "\timplicit none\n\tinteger :: ldim,neqn,ierr,ipar(4),i,j,k,n\n";
	print OUT "\treal(kind(1.d0)) :: t,c(neqn),dfdc(ldim,neqn)";
	if($l_mcmc_any){
		print OUT ",coef(" . ($n_mcmc_coefs+$variable_cs+2*$variable_ion_source) . ")";
	}elsif($variable_temp || $variable_cs || $variable_ion_source){
		print OUT ",coef(" . ($variable_temp+$variable_cs+2*$variable_ion_source) . ")";
	}else{
		print OUT ",coef(:)";
	}
	print OUT ",ml,mu\n";
	print OUT "\nend subroutine jeval$append_to_subroutine_names\n\n!-----------------------------------------------------------\n\n";
}

################################################################################
################################################################################
##      8.3 Determining if parameters for ions are used in the Matlab files   ##
################################################################################
################################################################################

# Check if we have ions; if not, don't print stuff related to them to the Matlab files
if(!$lfortran){
	if(!$include_ion_terms){
		# Input for the equation and flux files
		$input_for_eq='t,c,K,E,WL,CS,source,isconst,isfitted,C_fit';
		# this is for the case where a function is given for the concentration of some cluster
		$input_for_eq_C_fun='K,E,WL,CS,source,isconst,isfitted,C_fit';
		$input_for_flux='c,K,E,WL,CS,source,calc_sources';
	}else{
		$input_for_eq='t,c,K,E,WL,CS,source,isconst,isfitted,C_fit,fwl,fcs';
		$input_for_eq_C_fun='K,E,WL,CS,source,isconst,isfitted,C_fit,fwl,fcs';
		$input_for_flux='c,K,E,WL,CS,source,calc_sources,fwl,fcs';
	}
	if($l_reac_in_clus){
		$input_for_eq =~ s/K,E,/K,E,R,/;
		$input_for_eq_C_fun =~ s/K,E,/K,E,R,/;
		$input_for_flux =~ s/K,E,/K,E,R,/;
	}
	if($#source_function>=0){
		$input_for_flux =~ s/^c,/c,t,/;
	}

################################################################################
################################################################################
##      8.4 dC/dt for Matlab                                                  ##
################################################################################
################################################################################

	open(OUT, ">$output_file_name") || die "\nCould not open output file $output_file_name\n\n";
	
	print OUT "function cdot = $output_function_name($input_for_eq)\n\n%%%%%%%%%%%%%%%%%%%%% KEY %%%%%%%%%%%%%%%%\n";
	
	# Printing the key of what equation is what
	for($iclus=1;$iclus<=$max_n_flux;$iclus++){
		print OUT "$comments[$iclus]\n";
	}
	print OUT "%%%%%%%%%%%%%%%%%%%%% END KEY %%%%%%%%%%%%%%%%\n\n";

	print OUT "% Temperature = $temperature K\n\n";
	if(!$disable_dilution){
		print OUT "dil = $dil_str;\n\n";
	}

	print OUT "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n%%% Begin Equations Here\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n\ncdot = zeros(size(c));\n\n";
	
	# Fitting generic positive or negative ions to give charge balance
	if($include_generic_neg && $include_generic_pos){
		if($generic_positive_ions_corr ne ''){
			print OUT "% setting positive charger ion concentration to give charge balance\n";
			print OUT "if c($n_generic_neg)+$generic_positive_ions_corr>0\n";
			print OUT "c($n_generic_pos) = c($n_generic_neg)+$generic_positive_ions_corr;\n\n";
			print OUT "else\n";
			print OUT "\tc($n_generic_neg) = -($generic_positive_ions_corr);\n";
			print OUT "\tc($n_generic_pos) = 0;\n";
			print OUT "end\n\n";
			print OUT "isconst($n_generic_pos) = 1;\n\n";
		}elsif($generic_negative_ions_corr ne ''){
			print OUT "% setting negative charger ion concentration to give charge balance\n";
			print OUT "if c($n_generic_pos)+$generic_negative_ions_corr>0\n";
			print OUT "c($n_generic_neg) = c($n_generic_pos)+$generic_negative_ions_corr;\n\n";
			print OUT "else\n";
			print OUT "\tc($n_generic_pos) = -($generic_negative_ions_corr);\n";
			print OUT "\tc($n_generic_neg) = 0;\n";
			print OUT "end\n\n";
			print OUT "isconst($n_generic_neg) = 1;\n\n";
		}elsif($generic_positive_ions ne ''){
			print OUT "% setting positive charger ion concentration to give charge balance\n";
			print OUT "c($n_generic_pos) = $generic_positive_ions;\n";
			print OUT "isconst($n_generic_pos) = 1;\n\n";
		}elsif($generic_negative_ions ne ''){
			print OUT "% setting negative charger ion concentration to give charge balance\n";
			print OUT "c($n_generic_neg) = $generic_negative_ions;\n";
			print OUT "isconst($n_generic_neg) = 1;\n\n";
		}
	}

	# Equations for keeping a sum of concentrations constant
	print OUT "% Equations for keeping a sum of concentrations constant\n";
	print OUT "if any(isfitted)\n";
	print OUT "\tfor i = 1:$max_cluster\n";
	print OUT "\t\tif isfitted(i)\n";
	print OUT "\t\t\tc(i) = C_fit{i,1}-sum(c(C_fit{i,2}));\n";
	print OUT "\t\t\tif c(i) < 0\n";
	print OUT "\t\t\t\tc(C_fit{i,2}) = c(C_fit{i,2})*C_fit{i,1}/sum(c(C_fit{i,2}));\n";
	print OUT "\t\t\t\tc(i) = 0;\n";
	print OUT "\t\t\tend\n";
	print OUT "\t\tend\n";
	print OUT "\tend\n";
	print OUT "end\n$lines{'source'}";
	

	# Looping through the equations for all clusters and charger ions + the outgoing flux
	for($iclus=1;$iclus<=$max_n_flux;$iclus++){
		if($iclus>$max_cluster_number && (($ind_quad_form[$iclus][0]==0 && $ind_lin_form[$iclus][0]==0) || $iclus==$n_flux_rec)){
			next;
		}
		print OUT "\n$comments[$iclus]\n";
		if($iclus > $max_cluster_number){
			$line = "cdot($iclus) =";
		}else{
			$line = "if isconst($iclus),\n\tcdot($iclus) = 0;\nelse\n\tcdot($iclus) = source($iclus)";
		}
		$line_temp = "";
		# Losses in collisions
		for($i=1;$i<=2*$ind_quad_loss[$iclus][0];$i+=2){
			$jclus = $ind_quad_loss[$iclus][$i];
			$kclus = $ind_quad_loss[$iclus][$i+1];
			if($iclus == $jclus){
				$str_temp = "2*" . $coef_quad{"$iclus,$jclus,$kclus"} . "*c($iclus)*c($jclus)";
			}else{
				$str_temp = $coef_quad{"$iclus,$jclus,$kclus"} . "*c($iclus)*c($jclus)";
			}
			# Adding blank spaces to make all terms the same length -> output easier to read
			$line_temp = "$line_temp- $str_temp " . ' ' x max(30-length($str_temp),0);
			# Only three terms per line -> output easier to read
			if(($i+1)%6==0){
				$line = "$line...\n\t\t$line_temp";
				$line_temp = "";
			}
		}
		if($line_temp ne ''){
			$line = "$line...\n\t\t$line_temp";
			$line_temp = "";
		}
		# Formation in collisions
		for($i=1;$i<=2*$ind_quad_form[$iclus][0];$i+=2){
			$jclus = $ind_quad_form[$iclus][$i];
			$kclus = $ind_quad_form[$iclus][$i+1];
			if(exists $coef_quad_form{"$jclus,$kclus,$iclus"}){
				$str_temp = $coef_quad_form{"$jclus,$kclus,$iclus"} . "*c($jclus)*c($kclus)";
			}else{
				$str_temp = $coef_quad{"$jclus,$kclus,$iclus"} . "*c($jclus)*c($kclus)";
			}
			$line_temp = "$line_temp+ $str_temp " . ' ' x max(30-length($str_temp),0);
			if(($i+1)%6==0){
				$line = "$line...\n\t\t$line_temp";
				$line_temp = "";
			}
		}
		if($line_temp ne ''){
			$line = "$line...\n\t\t$line_temp";
			$line_temp = "";
		}
		# Extra formation in boundary collisions
		for($i=1;$i<=3*$ind_quad_form_extra[$iclus][0];$i+=3){
			$jclus = $ind_quad_form_extra[$iclus][$i];
			$kclus = $ind_quad_form_extra[$iclus][$i+1];
			if(exists $coef_quad_form{"$jclus,$kclus,$iclus"}){
				$str_temp = $coef_quad_form{"$jclus,$kclus,$iclus"} . "*c($jclus)*c($kclus)";
			}elsif($ind_quad_form_extra[$iclus][$i+2] != 1){
				$str_temp = $ind_quad_form_extra[$iclus][$i+2] . "*" . $coef_quad{"$jclus,$kclus,$iclus"} . "*c($jclus)*c($kclus)";
			}else{
				$str_temp = $coef_quad{"$jclus,$kclus,$iclus"} . "*c($jclus)*c($kclus)";
			}
			$line_temp = "$line_temp+ $str_temp " . ' ' x max(30-length($str_temp),0);
			if(($i+2)%9==0){
				$line = "$line...\n\t\t$line_temp";
				$line_temp = "";
			}
		}
		if($line_temp ne ''){
			$line = "$line...\n\t\t$line_temp";
			$line_temp = "";
		}
		# Losses by evaporation or external sinks
		for($i=1;$i<=2*$ind_lin_loss[$iclus][0];$i+=2){
			$jclus = $ind_lin_loss[$iclus][$i];
			$kclus = $ind_lin_loss[$iclus][$i+1];
			$str_temp = $coef_lin{"$jclus,$kclus,$iclus"} . "*c($iclus)";
			$line_temp = "$line_temp- $str_temp " . ' ' x max(30-length($str_temp),0);
			if(($i+1)%6==0){
				$line = "$line...\n\t\t$line_temp";
				$line_temp = "";
			}
		}
		if($line_temp ne ''){
			$line = "$line...\n\t\t$line_temp";
			$line_temp = "";
		}
		# Formation from evaporations
		for($i=1;$i<=2*$ind_lin_form[$iclus][0];$i+=2){
			$jclus = $ind_lin_form[$iclus][$i];
			$kclus = $ind_lin_form[$iclus][$i+1];
			if($iclus == $jclus && $iclus<=$max_cluster_number){
				$str_temp = "2*" . $coef_lin{"$iclus,$jclus,$kclus"} . "*c($kclus)";
			}else{
				$str_temp = $coef_lin{"$iclus,$jclus,$kclus"} . "*c($kclus)";
			}
			$line_temp = "$line_temp+ $str_temp " . ' ' x max(30-length($str_temp),0);
			if(($i+1)%6==0){
				$line = "$line...\n\t\t$line_temp";
				$line_temp = "";
			}
		}
		if($line_temp ne ''){
			$line = "$line...\n\t\t$line_temp";
			$line_temp = "";
		}
		# Losses by reactions inside clusters
		for($i=1;$i<=2*$ind_reac_loss[$iclus][0];$i+=2){
			$jclus = $ind_reac_loss[$iclus][$i];
			$str_temp = "R($iclus,$jclus)*c($iclus)";
			$line_temp = "$line_temp- $str_temp " . ' ' x max(30-length($str_temp),0);
			if(($i+1)%6==0){
				$line = "$line...\n\t\t$line_temp";
				$line_temp = "";
			}
		}
		if($line_temp ne ''){
			$line = "$line...\n\t\t$line_temp";
			$line_temp = "";
		}
		# Formation by reactions inside clusters
		for($i=1;$i<=2*$ind_reac_form[$iclus][0];$i+=2){
			$kclus = $ind_reac_form[$iclus][$i+1];
			$str_temp = "R($kclus,$iclus)*c($kclus)";
			$line_temp = "$line_temp+ $str_temp " . ' ' x max(30-length($str_temp),0);
			if(($i+1)%6==0){
				$line = "$line...\n\t\t$line_temp";
				$line_temp = "";
			}
		}
		if($line_temp ne ''){
			$line = "$line...\n\t\t$line_temp";
			$line_temp = "";
		}
		$line = "$line;\n";
		if($iclus<=$max_cluster_number){
			$line = $line . "end\n";
		}
		print OUT $line;
	}
	
	#Equations for keeping a sum of concentrations constant
	print OUT "\n\n% Equations for keeping a sum of concentrations constant\n";
	print OUT "if any(isfitted)\n";
	print OUT "\tfor i = 1:$max_cluster\n";
	print OUT "\t\tif isfitted(i)\n";
	print OUT "\t\t\tcdot(i) = -sum(cdot(C_fit{i,2}));\n";
	print OUT "\t\tend\n";
	print OUT "\tend\n";
	print OUT "end\n";
	print OUT "\nend\n";
	close(OUT);

################################################################################
################################################################################
##      8.5 Fluxes for Matlab                                                 ##
################################################################################
################################################################################

	open(FLUX, ">$flux_file_name") || die "\nCould not open flux file $flux_file_name\n\n";
	
	print FLUX "function [flux_2d, coll_evap_2d, sources, flux_3d, flux_out, clust_flux] = $flux_function_name($input_for_flux)\n\n%%%%%%%%%%%%%%%%%%%%% KEY %%%%%%%%%%%%%%%%\n";
	
	# Printing the key of what equation is what
	for($iclus=1;$iclus<=$max_n_flux_inc_boundary_clusters;$iclus++){
		print FLUX "$comments[$iclus]\n";
	}
	print FLUX "%%%%%%%%%%%%%%%%%%%%% END KEY %%%%%%%%%%%%%%%%\n\n";
	
	# print out the vector containing the cluster names plus the names of the other fluxes (sources, losses etc.)
	$line = "% name vector containing all the names (clusters, sources, losses etc.) in the flux matrices\nclust_flux = {";
	for($iclus=1;$iclus<=$max_n_flux_inc_boundary_clusters;$iclus++){
		$line = $line . "'" . $label_inc_boundary_clusters[$iclus] . "' ";
	}
	$line = "$line};\n\n";
	print FLUX $line;
	
	if(!$disable_dilution){
		print FLUX "dil = $dil_str;\n\n";
	}
	
	print FLUX "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n%%% Begin Equations Here\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n\n";
	
	if($include_generic_neg && $include_generic_pos){
		if($generic_positive_ions_corr ne ''){
			print FLUX "% setting positive charger ion concentration to give charge balance\n";
			print FLUX "if c($n_generic_neg)+$generic_positive_ions_corr>0\n";
			print FLUX "c($n_generic_pos) = c($n_generic_neg)+$generic_positive_ions_corr;\n\n";
			print FLUX "else\n";
			print FLUX "\tc($n_generic_neg) = -($generic_positive_ions_corr);\n";
			print FLUX "\tc($n_generic_pos) = 0;\n";
			print FLUX "end\n\n";
		}elsif($generic_negative_ions_corr ne ''){
			print FLUX "% setting negative charger ion concentration to give charge balance\n";
			print FLUX "if c($n_generic_pos)+$generic_negative_ions_corr>0\n";
			print FLUX "c($n_generic_neg) = c($n_generic_pos)+$generic_negative_ions_corr;\n\n";
			print FLUX "else\n";
			print FLUX "\tc($n_generic_pos) = -($generic_negative_ions_corr);\n";
			print FLUX "\tc($n_generic_neg) = 0;\n";
			print FLUX "end\n\n";
		}elsif($generic_positive_ions ne ''){
			print FLUX "% setting positive charger ion concentration to give charge balance\n"; 
			print FLUX "c($n_generic_pos) = $generic_positive_ions;\n";
			print FLUX "isconst($n_generic_pos) = 1;\n\n";
		}elsif($generic_negative_ions ne ''){
			print FLUX "% setting negative charger ion concentration to give charge balance\n";
			print FLUX "c($n_generic_neg) = $generic_negative_ions;\n";
			print FLUX "isconst($n_generic_neg) = 1;\n\n";
		}
	}

	print FLUX "flux_2d = zeros($max_n_flux_inc_boundary_clusters,$max_n_flux_inc_boundary_clusters);\n";
	if($print_3d_flux){
		print FLUX "flux_3d = zeros($max_n_flux_inc_boundary_clusters,$max_n_flux_inc_boundary_clusters,$max_n_flux_inc_boundary_clusters);\n";
	}else{
		print FLUX "flux_3d = [];\n";
	}
	print FLUX "coll_evap_2d = zeros($max_n_flux_inc_boundary_clusters,$max_n_flux_inc_boundary_clusters);\n";
	print FLUX "flux_out = zeros($max_cluster_number,$max_cluster_number,4);\n";
	print FLUX "sources = zeros($max_cluster_number,1);\n$lines{'source'}";
	
	if($print_3d_flux){
		print FLUX "% Note on 3D boundary fluxes:\n";
		print FLUX "% if you're keeping the boundary clusters, these are not working!\n";
		print FLUX "% sum(flux_3d(i,j,:)) is correct\n";
		print FLUX "% sum(flux_3d(:,:,k)) is correct but some of it comes from the boundary as flux_3d((index of the boundary cluster),(index of the boundary cluster),k)\n";
		print FLUX "% flux_3d(i,j,k1)/flux_3d(i,j,k2) is correct\n\n";
	}
	
	# Looping through the equations
	for($iclus=1;$iclus<=$max_cluster_number;$iclus++){
		print FLUX "\n$comments[$iclus]\n";
		if($ind_quad_loss[$iclus][0]>0){
			print FLUX "\t% Collisions removing this cluster\n";
		}
		# Losses in collisions + corresponding formation
		for($i=1;$i<=2*$ind_quad_loss[$iclus][0];$i+=2){
			$jclus = $ind_quad_loss[$iclus][$i];
			$kclus = $ind_quad_loss[$iclus][$i+1];
			$str_temp = $coef_quad{"$iclus,$jclus,$kclus"} . "*c($iclus)*c($jclus);\n";
			if(exists $coef_quad_form{"$iclus,$jclus,$kclus"}){
				$str_temp_form = $coef_quad_form{"$iclus,$jclus,$kclus"} . "*c($iclus)*c($jclus);\n";
			}else{
				$str_temp_form = $str_temp;
			}
			if($iclus==$jclus){
				# 2d-fluxes have a factor 2 for collisions of identical clusters except in fluxes from the boundary
				$str_temp2 = "2*$str_temp";
			}else{
				$str_temp2 = $str_temp;
			} 
			if(exists $l_quad_bound{"$iclus,$jclus,$kclus"} || exists $l_quad_nst{"$iclus,$jclus,$kclus"}){
				# 2-d fluxes to boundary and from boundary to the main product
				if(exists $l_quad_bound{"$iclus,$jclus,$kclus"} && $keep_boundary_clusters){
					$kclus_bound = $clus_bound_ind{"$iclus,$jclus,$kclus"};
					if($iclus<=$jclus && $print_3d_flux){
						print FLUX "\tflux_3d($iclus,$jclus,$kclus_bound) = flux_3d($iclus,$jclus,$kclus_bound)+1/$n_coll_products[$iclus][$jclus]*$str_temp_form";
						print FLUX "\tflux_3d($kclus_bound,$kclus_bound,$kclus) = flux_3d($kclus_bound,$kclus_bound,$kclus)+(1-1/$n_coll_products[$iclus][$jclus])*$str_temp_form";
					}
					print FLUX "\tcoll_evap_2d($iclus,$kclus_bound) = coll_evap_2d($iclus,$kclus_bound)+$str_temp2";
					# Print the fluxes from the boundary only once
					if($iclus>=$jclus){
						print FLUX "\tcoll_evap_2d($kclus_bound,$kclus) = coll_evap_2d($kclus_bound,$kclus)+$str_temp_form";
					}
				}else{
					if($iclus<=$jclus && $print_3d_flux){
						print FLUX "\tflux_3d($iclus,$jclus,$kclus) = flux_3d($iclus,$jclus,$kclus)+1/$n_coll_products[$iclus][$jclus]*$str_temp_form";
						print FLUX "\tflux_3d($n_flux_bound,$n_flux_bound,$kclus) = flux_3d($n_flux_bound,$n_flux_bound,$kclus)+(1-1/$n_coll_products[$iclus][$jclus])*$str_temp_form";
					}
					print FLUX "\tcoll_evap_2d($iclus,$n_flux_bound) = coll_evap_2d($iclus,$n_flux_bound)+$str_temp2";
					# Print the fluxes from the boundary only once
					if($iclus>=$jclus){
						print FLUX "\tcoll_evap_2d($n_flux_bound,$kclus) = coll_evap_2d($n_flux_bound,$kclus)+$str_temp_form";
					}
				}
			}else{
				if($iclus<=$jclus && $print_3d_flux){
					print FLUX "\tflux_3d($iclus,$jclus,$kclus) = flux_3d($iclus,$jclus,$kclus)+$str_temp";
				}
				print FLUX "\tcoll_evap_2d($iclus,$kclus) = coll_evap_2d($iclus,$kclus)+$str_temp2";
			}
		}
		# Extra loss pathways in boundary collisions + corresponding formation
		if($ind_quad_loss_extra[$iclus][0]>0){
			print FLUX "\t% More products from the same collisions (monomers from the boundary etc.)\n";
		}
		for($i=1;$i<=3*$ind_quad_loss_extra[$iclus][0];$i+=3){
			$jclus = $ind_quad_loss_extra[$iclus][$i];
			# Print these only once
			if($iclus<=$jclus && !($i>1 && $jclus==$ind_quad_loss_extra[$iclus][$i-3])){
				$kclus = $ind_quad_loss_extra[$iclus][$i+1];
				if(exists $coef_quad_form{"$iclus,$jclus,$kclus"}){
					$str_temp_form = $coef_quad_form{"$iclus,$jclus,$kclus"} . "*c($iclus)*c($jclus);\n";
				}else{
					$str_temp_form = $coef_quad{"$iclus,$jclus,$kclus"} . "*c($iclus)*c($jclus);\n";
				}
				if($print_3d_flux){
					if(exists $l_quad_bound{"$iclus,$jclus,$kclus"} && $keep_boundary_clusters){
						$kclus_bound = $clus_bound_ind{"$iclus,$jclus,$kclus"};
						print FLUX "\tflux_3d($iclus,$jclus,$kclus_bound) = flux_3d($iclus,$jclus,$kclus_bound)+$ind_quad_loss_extra[$iclus][$i+2]/$n_coll_products[$iclus][$jclus]*$str_temp_form";
						print FLUX "\tflux_3d($kclus_bound,$kclus_bound,$kclus) = flux_3d($kclus_bound,$kclus_bound,$kclus)+$ind_quad_loss_extra[$iclus][$i+2]*(1-1/$n_coll_products[$iclus][$jclus])*$str_temp_form";
					}else{
						print FLUX "\tflux_3d($iclus,$jclus,$kclus) = flux_3d($iclus,$jclus,$kclus)+$ind_quad_loss_extra[$iclus][$i+2]/$n_coll_products[$iclus][$jclus]*$str_temp_form";
						print FLUX "\tflux_3d($n_flux_bound,$n_flux_bound,$kclus) = flux_3d($n_flux_bound,$n_flux_bound,$kclus)+$ind_quad_loss_extra[$iclus][$i+2]*(1-1/$n_coll_products[$iclus][$jclus])*$str_temp_form";
					}
				}
				if($ind_quad_loss_extra[$iclus][$i+2] != 1){
					$str_temp_form =  "$ind_quad_loss_extra[$iclus][$i+2]*$str_temp_form";
				}
				if(exists $l_quad_bound{"$iclus,$jclus,$kclus"} && $keep_boundary_clusters){
					$kclus_bound = $clus_bound_ind{"$iclus,$jclus,$kclus"};
					print FLUX "\tcoll_evap_2d($kclus_bound,$kclus) = coll_evap_2d($kclus_bound,$kclus)+$str_temp_form";
				}else{
					print FLUX "\tcoll_evap_2d($n_flux_bound,$kclus) = coll_evap_2d($n_flux_bound,$kclus)+$str_temp_form";
				}
			}
		}
		#  Formation by evaporation + corresponding losses
		if($ind_lin_form[$iclus][0]>0){
			print FLUX "\t% Evaporations resulting in this cluster\n";
		}
		for($i=1;$i<=2*$ind_lin_form[$iclus][0];$i+=2){
			$jclus = $ind_lin_form[$iclus][$i];
			$kclus = $ind_lin_form[$iclus][$i+1];
			$str_temp = $coef_lin{"$iclus,$jclus,$kclus"} . "*c($kclus);\n";
			if($iclus<=$jclus && $print_3d_flux){
				print FLUX "\tflux_3d($iclus,$jclus,$kclus) = flux_3d($iclus,$jclus,$kclus)-$str_temp";
			}
			if($iclus==$jclus){
				# 2d-fluxes have a factor 2 for evaporations producing identical clusters
				$str_temp = "2*$str_temp";
			} 
			print FLUX "\tcoll_evap_2d($kclus,$iclus) = coll_evap_2d($kclus,$iclus)+$str_temp";
		}
		# Losses to external sinks
		if(!$disable_wall_terms || !$disable_dilution || !$disable_coag_sinks){
			print FLUX "\t% External losses\n";
		}
		for($i=1;$i<=2*$ind_lin_loss[$iclus][0];$i+=2){
			$jclus = $ind_lin_loss[$iclus][$i];
			$kclus = $ind_lin_loss[$iclus][$i+1];
			if($jclus>$max_cluster_number){
				$str_temp = $coef_lin{"$jclus,$kclus,$iclus"} . "*c($iclus);\n";
				if($print_3d_flux){
					print FLUX "\tflux_3d($jclus,$kclus,$iclus) = flux_3d($jclus,$kclus,$iclus)-$str_temp";
				}
				print FLUX "\tcoll_evap_2d($iclus,$kclus) = coll_evap_2d($iclus,$kclus)+$str_temp";
			}
		}
		# Losses by reactions inside clusters + corresponding formation
		if($ind_reac_loss[$iclus][0]>0){
			print FLUX "\t% Chemical reactions removing this cluster\n";
		}
		for($i=1;$i<=2*$ind_reac_loss[$iclus][0];$i+=2){
			$jclus = $ind_reac_loss[$iclus][$i];
			$kclus = $ind_reac_loss[$iclus][$i+1];
			$str_temp = "R($iclus,$jclus)*c($iclus);\n";
			print FLUX "\tflux_3d($jclus,$kclus,$iclus) = flux_3d($jclus,$kclus,$iclus)-$str_temp";
			print FLUX "\tcoll_evap_2d($iclus,$jclus) = coll_evap_2d($iclus,$jclus)+$str_temp";
		}
	}

	# Flux matrix for outgoing collisions: flux_out(cluster 1, cluster 2, charge of formed big cluster)
	print FLUX "\n\n% Matrix for collisions that lead succesfully out of the system:\n% flux_out(i,j,k) is the outgoing flux from clusters i and j, where\n% k=1 for neutral-neutral, 2 for neutral-negative, 3 for neutral-positive and 4 for recombination collisions\n";
	for($k=$n_flux_out; $k<=$n_flux_out_max; $k++){
		for($i=1;$i<=2*$ind_quad_form[$k][0];$i+=2){
			$iclus = $ind_quad_form[$k][$i];
			$jclus = $ind_quad_form[$k][$i+1];
			if(exists $coef_quad_form{"$iclus,$jclus,$k"}){
				$line_temp = $coef_quad_form{"$iclus,$jclus,$k"} . "*c($iclus)*c($jclus);\n";
			}else{
				$line_temp = $coef_quad{"$iclus,$jclus,$k"} . "*c($iclus)*c($jclus);\n";
			}
			if($k==$n_flux_out && $lneutral_cluster[$iclus] && $lneutral_cluster[$jclus]){
				$j = 1;
			}elsif($k==$n_flux_out_neg){
				$j = 2;
			}elsif($k==$n_flux_out_pos){
				$j = 3;
			}else{
				$j = 4;
			}
			print FLUX "\tflux_out($iclus,$jclus,$j) = $line_temp";
		}
	}

	# Combine the collision and evaporation fluxes to net fluxes
	print FLUX "\n% Combining collisions and evaporations into net fluxes.\n";
	print FLUX "\tfor i = 1:" . ($max_n_flux_inc_boundary_clusters) . "\n\t\tfor j = (i+1):" . ($max_n_flux_inc_boundary_clusters) . "\n";
	print FLUX "\t\t\tnet = coll_evap_2d(i,j)-coll_evap_2d(j,i);\n\t\t\tif net > 0\n\t\t\t\tflux_2d(i,j) = net;\n\t\t\telse\n";
	print FLUX "\t\t\t\tflux_2d(j,i) = -net;\n\t\t\tend\n\t\tend\n\tend\n\n";
	
	# Adding source terms to the 2d flux matrices
	print FLUX "\n\n% Adding given source terms to the flux matrix and solving source terms of monomers, if needed\n";
	print FLUX "\tflux_2d($n_flux_source,1:$max_cluster_number) = source';\n\n";
	print FLUX "\tif calc_sources == 1\n";
	for($iclus=1;$iclus<=$max_cluster_number;$iclus++){
		if ($lneutral_cluster[$iclus] && $lmonomer[$iclus]){
			print FLUX "\t\tsources($iclus) = sum(coll_evap_2d($iclus,:))-sum(coll_evap_2d([1:$max_n_flux_inc_boundary_clusters],$iclus));\n";
			print FLUX "\t\tif source($iclus) == 0\n\t\t\tflux_2d($n_flux_source,$iclus) = sources($iclus);\n\t\tend\n";
		}elsif($iclus==$n_generic_neg || $iclus==$n_generic_pos){
			print FLUX "\t\tsources($iclus) = sum(coll_evap_2d($iclus,:));\n";
			print FLUX "\t\tif source($iclus) == 0\n\t\t\tflux_2d($n_flux_source,$iclus) = sources($iclus);\n\t\tend\n";
		}
	}
	print FLUX "\tend\n";
	print FLUX "\nend\n";
	close(FLUX);
}

################################################################################
################################################################################
##      8.6 Formation rate for Fortran                                        ##
################################################################################
################################################################################

if($lfortran){
	print OUT "subroutine formation$append_to_subroutine_names(neqn,c";
	if($lloop){
		print OUT ",j_out";
		if($l_jlim){
			print OUT ",j_lim";
			if($lloop_j_frac){
				print OUT ",j_lim_frac";
			}
		}
	}else{
		print OUT ",j_tot,j_by_charge,j_by_cluster,j_all";
	}
	print OUT ",ipar";
	if($l_mcmc_any || $variable_temp || $variable_cs || $variable_ion_source){
		print OUT ",coef";
	}
	print OUT ")\n";
	if($l_j_in){
		print OUT "use shared_input, only : r_j_in => rlim_for_J\n";
		print "\nUsing module shared_input in the formation subroutine to get r_j_in\n";
	}
	print OUT "\timplicit none\n";
	if($lloop_cs && $nmol_types_used==1){
		print OUT "\tinteger, parameter :: nclust = $max_cluster_number, n_cs = $max_cluster_number\n";
	}elsif($lloop && $nmol_types_used>1){
		print OUT "\tinteger, parameter :: nclust = $max_cluster_number, ij_ind_max($nmol_types_used) = (/$nmol_ranges_list/)";
		if($disable_nonmonomers){
			print OUT ", n_monomers($nmol_types_used) = (/$monomer_indices/)";
		}
		print OUT "\n";
	}else{
		print OUT "\tinteger, parameter :: nclust = $max_cluster_number\n";
	}
	print OUT "\tinteger :: neqn, i";
	if($include_negative_ion_terms || $include_positive_ion_terms || $lloop){
		print OUT ", j";
	}
	if($lloop && $nmol_types_used>1){
		print OUT ", ij_ind($nmol_types_used), imol";
	}
	print OUT ", n, ipar(4)\n";
	if($lloop_j_tests || $lgrowth_for_j_tests || $l_j_in){
		print OUT "\treal(kind(1.d0)) :: r, m\n";
		if($l_j_in){
			print OUT "\tinteger, save :: n_j_in = 0\n";
		}
	}
	print OUT "\treal(kind(1.d0)) :: c(neqn)";
	if($lloop){
		print OUT ", j_out";
	}else{
		print OUT ", j_add, j_tot,j_by_charge(4),j_by_cluster(neqn),j_all(neqn,4)";
	}
	if($l_mcmc_any){
		print OUT ",coef(" . ($n_mcmc_coefs+$variable_cs+2*$variable_ion_source) . ")";
	}elsif($variable_temp || $variable_cs || $variable_ion_source){
		print OUT ",coef(" . ($variable_temp+$variable_cs+2*$variable_ion_source) . ")";
	}
	print OUT "\n";
	if($lloop){
		print OUT "\treal(kind(1.d0)) :: pfac\n";
		print OUT "\treal(kind(1.d0)), save :: K(nclust,nclust)\n";
		if($l_jlim && !$disable_evap){
			print OUT "\treal(kind(1.d0)), save :: E(";
			if($disable_nonmonomer_evaps && $nmol_types_used==1){
				print OUT "1";
			}else{
				print OUT "nclust";
			}
			print OUT ",nclust)\n";
		}
		if(!$lloop_diag && $nmol_types_used>1){
			print OUT "\tinteger, save :: indices(nclust,$nmol_types_used)";
			if($l_jlim && $lloop_max_ratio){
				print OUT ", n_pure($nmol_ranges_array[$nmol_max_ratio_basis])";
			}
			print OUT "\n";
		}
	}else{
		print OUT "\treal(kind(1.d0)), save :: coef_quad(nclust,nclust,$max_n_flux)=0.d0,coef_lin($max_n_flux,$max_n_flux,nclust)=0.d0,source($size_C_out)=0.d0\n";
		print OUT "\tinteger, save :: charge(nclust)=0,ind_quad_loss($size_C_out,0:$ind_quad_loss[0])=0,ind_quad_form($size_C_out,0:$ind_quad_form[0])=0\n";
		print OUT "\tinteger, save :: ind_quad_loss_extra($size_C_out,0:$ind_quad_loss_extra[0])=0,ind_quad_form_extra($size_C_out,0:$ind_quad_form_extra[0])=0\n";
		print OUT "\tinteger, save :: ind_lin_loss($size_C_out,0:$ind_lin_loss[0])=0,ind_lin_form($size_C_out,0:$ind_lin_form[0])=0,fitted(nclust,0:nclust)=0\n";
		print OUT "\tlogical, save :: isconst($size_C_out)=.false.\n\n";
	}
	if($l_jlim){
		print OUT "\tinteger, parameter :: nlim = $nlim\n";
		print OUT "\treal(kind(1.d0)), parameter :: r_lim(nlim) = (/$jlim_sizes_list/)*1.d-9\n";
		print OUT "\tinteger :: l, ij, nmin\n";
		print OUT "\treal(kind(1.d0)) :: j_lim(nlim)";
		if($lloop_j_frac){
			print OUT ", j_lim_frac(nlim,neqn)";
		}
		print OUT "\n";
		print OUT "\treal(kind(1.d0)) :: ri, rj, rk, mi, mj, mk\n\n";
	}
	print OUT "\tif (ipar(3) .eq. 0) then\n";
	print OUT "\t\tipar(3) = 1\n";
	if($lloop){
		# The current print-out of the formation rate in the loop mode is a nightmare;
		# it should be re-organized, but meanwhile be careful when you make any modifications
		print OUT "\t\tcall get_coll$append_to_subroutine_names(K)\n";
		if($l_jlim && !$disable_evap){
			print OUT "\t\tcall get_evap$append_to_subroutine_names(E)\n";
		}
		if(!$lloop_diag && $nmol_types_used>1){
			print OUT "\t\tcall get_molecule_numbers$append_to_subroutine_names(indices)\n";
			if($l_jlim && $lloop_max_ratio){
				print OUT "\t\tcall pure_clust_indices$append_to_subroutine_names(n_pure)\n";
			}
		}
		if($l_j_in){
			print OUT "\t\tdo i=1,nclust\n";
			print OUT "\t\t\tcall get_masses_and_radii$append_to_subroutine_names(mi,ri,i)\n";
			print OUT "\t\t\tif (ri*1.d9 .ge. r_j_in) then\n";
			print OUT "\t\t\t\tn_j_in = i\n";
			print OUT "\t\t\t\twrite(*,*) 'Size for which the input J will be used:'\n";
			print OUT "\t\t\t\twrite(*,*) r_j_in,' nm, ',n_j_in,' molecules'\n";
			print OUT "\t\t\t\texit\n";
			print OUT "\t\t\tend if\n";
			print OUT "\t\tend do\n";
		}
		if($l_jlim){
			if(!$lloop_diag && $nmol_types_used>1){
				print OUT "\t\tcall get_masses_and_radii$append_to_subroutine_names(mk,rk,ij_ind_max)\n";
			}elsif(!$lloop_diag && $nmol_types_used==1){
				print OUT "\t\tcall get_masses_and_radii$append_to_subroutine_names(mk,rk,nclust)\n";
			}
			# some kind of check (not all comprehensive) for if the given size limits are reasonable
			print OUT "\t\tif (any(r_lim.gt.rk)) then\n";
			print OUT "\t\t\twrite(*,*)'FYI: at least one of the given size limits is outside the largest system boundary limit r=',rk*1.d9,'nm!'\n";
			if(!$lloop_diag && $nmol_types_used==1){
				# another check (mainly related to generating look-up tables for different sizes)
				print OUT "\t\t\tcall get_masses_and_radii$append_to_subroutine_names(mk,rk,nclust+1)\n";
				print OUT "\t\t\twrite(*,*)'First outgrown cluster: r=',rk*1.d9,'nm'\n";
			}
			print OUT "\t\tend if\n";
		}
		print OUT "\tend if\n\n";
		print OUT "\tj_out = 0.d0\n";
		if($l_jlim){
			print OUT "\tj_lim = 0.d0\n";
			if($lloop_j_frac){
				print OUT "\tj_lim_frac = 0.d0\n";
			}
		}
		if($nmol_types_used==1 && $disable_nonmonomers){
			print OUT "\tj_out = j_out+K(nclust,1)*c(nclust)*c(1)\n";
		}
		if($nmol_types_used>1 || $l_jlim || !$disable_nonmonomers || $lgrowth_for_j_tests){
			# Print a loop over the cluser sizes to go through specific collisions
			if(!$l_j_in){
				print OUT "\tdo i=1,nclust\n";
			}else{
				print OUT "\tdo i=n_j_in,nclust\n";
				print OUT "\t\t! collisions among clusters taken into account in the distribution\n";
			}
			if(!$lloop_diag && $nmol_types_used>1){
				if($lloop_max_ratio){
					if($lfree_mon){
						print OUT "\t\tif ((indices(i,$nmol_max_ratio_basis) .eq. 1) .and. (indices(i,$nmol_not_basis) .gt. 0)) then\n";
						print OUT "\t\t\tcycle\n";
						print OUT "\t\tend if\n";
					}
				}
				if($disable_nonmonomers){
					print OUT "\t\tdo imol = 1,$nmol_types_used\n";
					print OUT "\t\t\tj = n_monomers(imol)\n";
					# Don't take the same process twice
					print OUT "\t\t\tif ((sum(indices(i,:)).eq.1) .and. (j.lt.i)) then\n";
					print OUT "\t\t\t\tcycle\n";
					print OUT "\t\t\tend if\n";
				}else{
					print OUT "\t\tdo j=i,nclust\n";
				}
				print OUT "\t\t\tij_ind = indices(i,:)+indices(j,:)\n";
				if($lloop_max_ratio){
					# skip the collisions resulting in a wrong composition (i.e. those of the "non-basis" molecule with clusters that already have the maximum number of these molecules)
					if($lfree_mon){
						print OUT "\t\t\tif (((ij_ind($nmol_max_ratio_basis) .eq. 1) .and. (ij_ind($nmol_not_basis) .gt. 0)) .or. (ij_ind($nmol_not_basis) .gt. $max_wrt_mon*ij_ind($nmol_max_ratio_basis))) then\n";
					}else{
						print OUT "\t\t\tif (ij_ind($nmol_not_basis) .gt. $max_wrt_mon*ij_ind($nmol_max_ratio_basis)) then\n";
					}
					print OUT "\t\t\t\tcycle\n";
					print OUT "\t\t\tend if\n";
				}
				$str_temp = $lines{'pfac'};
				$str_temp =~ s/\t(?!\t)/\t\t\t/g;
				print OUT $str_temp;
				if($lloop_cs){
					print OUT "\t\t\tif (any(ij_ind.ge.ij_ind_max)) then\n";
				}else{
					print OUT "\t\t\tif (any(ij_ind.gt.ij_ind_max)) then\n";
				}
				print OUT "\t\t\t\tj_out = j_out+pfac*K(i,j)*c(i)*c(j)\n";
				print OUT "\t\t\tend if\n";
			}elsif(!$lloop_diag && $nmol_types_used==1 && !$disable_nonmonomers){
				print OUT "\t\tdo j=max(i,$n_flux_out-i),nclust\n";
				$str_temp = $lines{'pfac'};
				$str_temp =~ s/\t(?!\t)/\t\t\t/g;
				print OUT $str_temp;
				print OUT "\t\t\tj_out = j_out+pfac*K(i,j)*c(i)*c(j)\n";
				print OUT "\t\tend do\n";
			}
			if($l_jlim){
				if($nmol_types_used>1){
					print OUT "\t\t\tcall get_masses_and_radii$append_to_subroutine_names(mi,ri,indices(i,:))\n";
					print OUT "\t\t\tcall get_masses_and_radii$append_to_subroutine_names(mj,rj,indices(j,:))\n";
					print OUT "\t\t\tcall get_masses_and_radii$append_to_subroutine_names(mk,rk,ij_ind)\n";
				}else{
					print OUT "\t\tcall get_masses_and_radii$append_to_subroutine_names(mi,ri,i)\n";
					if($disable_nonmonomers){
						print OUT "\t\tdo j=1,1\n";
					}else{
						print OUT "\t\tdo j=i,nclust\n";
					}
					$str_temp = $lines{'pfac'};
					$str_temp =~ s/\t(?!\t)/\t\t\t/g;
					print OUT $str_temp;
					print OUT "\t\t\tcall get_masses_and_radii$append_to_subroutine_names(mj,rj,j)\n";
					print OUT "\t\t\tcall get_masses_and_radii$append_to_subroutine_names(mk,rk,i+j)\n";
				}
				print OUT "\t\t\tdo l=1,nlim\n";
				print OUT "\t\t\t\tif ((ri.lt.r_lim(l)) .and. (rj.lt.r_lim(l)) .and. (rk.ge.r_lim(l))) then\n";
				print OUT "\t\t\t\t\tj_lim(l) = j_lim(l)+pfac*K(i,j)*c(i)*c(j)\n";
				if($lloop_j_frac){
					if($nmol_types_used>1){
						print OUT "\t\t\t\t\tif (rj .lt. ri) then\n";
						print OUT "\t\t\t\t\t\tnmin = j\n";
						print OUT "\t\t\t\t\telse\n";
						print OUT "\t\t\t\t\t\tnmin = i\n";
						print OUT "\t\t\t\t\tend if\n";
					}else{
						print OUT "\t\t\t\t\tnmin = min(i,j)\n";
					}				
					print OUT "\t\t\t\t\tj_lim_frac(l,nmin) = j_lim_frac(l,nmin)+pfac*K(i,j)*c(i)*c(j)\n";
				}
				if(!$disable_evap){
					print OUT "\t\t\t\t\t$lines{'l_evap'}";
					if($nmol_types_used>1){
						print OUT "\t\t\t\t\t\tij = $cluster_from_indices\n";
					}else{
						print OUT "\t\t\t\t\t\tij = i+j\n";
					}
					print OUT "\t\t\t\t\t\tj_lim(l) = j_lim(l)-E(i,j)*c(ij)\n";
					if($lloop_j_frac){
						print OUT "\t\t\t\t\t\tj_lim_frac(l,nmin) = j_lim_frac(l,nmin)-E(i,j)*c(ij)\n";
					}
					print OUT "\t\t\t\t\tend if\n";
				}
				print OUT "\t\t\t\tend if\n";
				print OUT "\t\t\tend do\n";
			}
			if($nmol_types_used>1 || $l_jlim){
				print OUT "\t\tend do\n";
			}
			if($l_j_in && !$disable_nonmonomers){
				# Deal with the clusters colliding with a monomer here, since the monomer isn't included in the equation loop
				print OUT "\t\t! monomer collisions\n";
				print OUT "\t\tj = 1\n";
				print OUT "\t\tif (i+j.ge.$n_flux_out) then\n";
				print OUT "\t\t\tj_out = j_out+K(i,j)*c(i)*c(j)\n";
				print OUT "\t\tend if\n";
				if($l_jlim){
					print OUT "\t\tcall get_masses_and_radii$append_to_subroutine_names(mj,rj,j)\n";
					print OUT "\t\tcall get_masses_and_radii$append_to_subroutine_names(mk,rk,i+j)\n";
					print OUT "\t\tdo l=1,nlim\n";
					print OUT "\t\t\tif ((ri.lt.r_lim(l)) .and. (rj.lt.r_lim(l)) .and. (rk.ge.r_lim(l))) then\n";
					print OUT "\t\t\t\tj_lim(l) = j_lim(l)+K(i,j)*c(i)*c(j)\n";
					if($lloop_j_frac){
						if($nmol_types_used>1){
							print OUT "\t\t\t\tif (rj .lt. ri) then\n";
							print OUT "\t\t\t\t\tnmin = j\n";
							print OUT "\t\t\t\telse\n";
							print OUT "\t\t\t\t\tnmin = i\n";
							print OUT "\t\t\t\tend if\n";
						}else{
							print OUT "\t\t\t\tnmin = min(i,j)\n";
						}
						print OUT "\t\t\t\tj_lim_frac(l,nmin) = j_lim_frac(l,nmin)+K(i,j)*c(i)*c(j)\n";
					}
					if(!$disable_evap){
						print OUT "\t\t\t\t$lines{'l_evap'}";
						if($nmol_types_used>1){
							print OUT "\t\t\t\t\tij = $cluster_from_indices\n";
						}else{
							print OUT "\t\t\t\t\tij = i+j\n";
						}
						print OUT "\t\t\t\t\tj_lim(l) = j_lim(l)-E(j,i)*c(ij)\n";
						if($lloop_j_frac){
							print OUT "\t\t\t\t\tj_lim_frac(l,nmin) = j_lim_frac(l,nmin)-E(j,i)*c(ij)\n";
						}
						print OUT "\t\t\t\tend if\n";
					}
					print OUT "\t\t\tend if\n";
					print OUT "\t\tend do\n";
				}
			}
			if($lgrowth_for_j_tests){
				print OUT "\t\t! collisions with the implicit growth compound\n";
				print OUT "\t\tcall get_masses_and_radii(m,r,i)\n";
				print OUT "\t\tif (r .ge. $rlim_growth_compound_str) then\n";
				print OUT "\t\t\tj = 2\n";
				print OUT "\t\t\tif (i+j.ge.$n_flux_out) then\n";
				print OUT "\t\t\t\tj_out = j_out+K(i,j)*c(i)*$c_growth_compound_str\n";
				print OUT "\t\t\tend if\n";
				if($l_jlim){
					print OUT "\t\t\tcall get_masses_and_radii$append_to_subroutine_names(mj,rj,j)\n";
					print OUT "\t\t\tcall get_masses_and_radii$append_to_subroutine_names(mk,rk,i+j)\n";
					print OUT "\t\t\tdo l=1,nlim\n";
					print OUT "\t\t\t\tif ((ri.lt.r_lim(l)) .and. (rj.lt.r_lim(l)) .and. (rk.ge.r_lim(l))) then\n";
					print OUT "\t\t\t\t\tj_lim(l) = j_lim(l)+K(i,j)*c(i)*$c_growth_compound_str\n";
					if($lloop_j_frac){
						print OUT "\t\t\t\t\tnmin = min(i,j)\n";
						print OUT "\t\t\t\t\tj_lim_frac(l,nmin) = j_lim_frac(l,nmin)+K(i,j)*c(i)*$c_growth_compound_str\n";
					}
					print OUT "\t\t\t\tend if\n";
					print OUT "\t\t\tend do\n";
				}
				print OUT "\n";
				print OUT "\t\t\tj = 3\n";
				print OUT "\t\t\tif (i+j.ge.$n_flux_out) then\n";
				print OUT "\t\t\t\tj_out = j_out+K(i,j)*c(i)*$c_growth_compound_str\n";
				print OUT "\t\t\tend if\n";
				if($l_jlim){
					print OUT "\t\t\tcall get_masses_and_radii$append_to_subroutine_names(mj,rj,j)\n";
					print OUT "\t\t\tcall get_masses_and_radii$append_to_subroutine_names(mk,rk,i+j)\n";
					print OUT "\t\t\tdo l=1,nlim\n";
					print OUT "\t\t\t\tif ((ri.lt.r_lim(l)) .and. (rj.lt.r_lim(l)) .and. (rk.ge.r_lim(l))) then\n";
					print OUT "\t\t\t\t\tj_lim(l) = j_lim(l)+K(i,j)*c(i)*$c_growth_compound_str\n";
					if($lloop_j_frac){
						print OUT "\t\t\t\t\tnmin = min(i,j)\n";
						print OUT "\t\t\t\t\tj_lim_frac(l,nmin) = j_lim_frac(l,nmin)+K(i,j)*c(i)*$c_growth_compound_str\n";
					}
					print OUT "\t\t\t\tend if\n";
					print OUT "\t\t\tend do\n";
				}
				print OUT "\t\tend if\n";
			}
			print OUT "\tend do\n";
		}
		
	}else{
		print OUT "\t\tcall initialize_parameters$append_to_subroutine_names(coef_quad,coef_lin,ind_quad_loss,ind_quad_form,&\n\t\t\t&ind_lin_loss,ind_lin_form,ind_quad_loss_extra,ind_quad_form_extra,&\n\t\t&\tsource,isconst,fitted";
		if($l_mcmc_any || $variable_temp || $variable_cs || $variable_ion_source){
			print OUT ",coef";
		}
		print OUT ",ipar)\n";
		print OUT "\t\t! cluster charges\n\t\tcharge = (/";
		$line = '';
		$i = 1;
		for($iclus=1;$iclus <= $max_cluster_number;$iclus++){
			if($i==30){
				$line = "$line&\n\t\t\t\t&";
				$i = 0;
			}
			$i += 1;
			if($lneutral_cluster[$iclus]){
				$line = "$line 0,";
			}elsif(($iclus==$n_generic_neg) || $lnegative_cluster[$iclus]){
				$line = "$line-1,";
			}elsif(($iclus==$n_generic_pos) || $lpositive_cluster[$iclus]){
				$line = "$line 1,";
			}else{
				die "\nThe charge of cluster $iclus ($label[$iclus]) is not defined!\n\n";
			}
		}
		$line =~ s/,$//;
		print OUT "$line /)\n";
		print OUT "\tend if\n";
		print OUT "\tj_tot = 0.d0\t\t\t! total formation rate\n";
		print OUT "\tj_by_charge = 0.d0\t\t! contributions of neutrals, negatives, positives and recombinations\n";
		print OUT "\tj_by_cluster = 0.d0\t\t! contributions of different clusters\n";
		print OUT "\tj_all = 0.d0\t\t\t! contributions of different clusters and different charging states\n";
		if($include_generic_neg && $include_generic_pos){
			if($generic_positive_ions_corr ne ''){
				print OUT "\n\t! setting positive charger ion concentration to give charge balance\n";
				print OUT "\tif (c($n_generic_neg)+$generic_positive_ions_corr>0) then\n";
				print OUT "\tc($n_generic_pos) = c($n_generic_neg)+$generic_positive_ions_corr\n";
				print OUT "\telse\n";
				print OUT "\t\tc($n_generic_neg) = -($generic_positive_ions_corr)\n";
				print OUT "\t\tc($n_generic_pos) = 0\n";
				print OUT "\tend if\n\n";
			}elsif($generic_negative_ions_corr ne ''){
				print OUT "\n\t! setting negative charger ion concentration to give charge balance\n";
				print OUT "\tif (c($n_generic_pos)+$generic_negative_ions_corr>0) then\n";
				print OUT "\tc($n_generic_neg) = c($n_generic_pos)+$generic_negative_ions_corr\n";
				print OUT "\telse\n";
				print OUT "\t\tc($n_generic_pos) = -($generic_negative_ions_corr)\n";
				print OUT "\t\tc($n_generic_neg) = 0\n";
				print OUT "\tend if\n\n";
			}elsif($generic_positive_ions ne ''){
				print OUT "\n\t! setting positive charger ion concentration to give charge balance\n";
				print OUT "\tc($n_generic_pos) = $generic_positive_ions\n\n";
			}elsif($generic_negative_ions ne ''){
				print OUT "\n\t! setting negative charger ion concentration to give charge balance\n";
				print OUT "\tc($n_generic_neg) = $generic_negative_ions\n\n";
			}
			if($l_q_neg_mcmc){
				print OUT "\n\t! negative charger ion source term from mcmc coefficient\n";
				print OUT "\tsource($n_generic_neg) = coef(" . ($ncoef_coll+$ncoef_evap+1) . ")\n\n";
			}
			if($l_q_pos_mcmc){
				print OUT "\n\t! positive charger ion source term from mcmc coefficient\n";
				print OUT "\tsource($n_generic_pos) = coef(" . ($ncoef_coll+$ncoef_evap+$l_q_neg_mcmc+1) . ")\n\n";
			}
		}
		if($include_negative_ion_terms || $include_positive_ion_terms){
			print OUT "\tdo j=$n_flux_out, $n_flux_out_max\n";
			print OUT "\t\tdo n=1, 2*ind_quad_form(j,0)-1, 2 ! loop over all collisions leading out of the system\n";
			if($include_negative_ion_terms && $include_positive_ion_terms){
				print OUT "\t\t\tif ((j .eq. $n_flux_out) .and. (charge(ind_quad_form(j,n)) .eq. 0) .and. (charge(ind_quad_form(j,n+1)) .eq. 0)) then\n";
				print OUT "\t\t\t\ti = 1 ! neutral-neutral collision\n";
				print OUT "\t\t\telseif (j .eq. $n_flux_out_neg) then\n";
				print OUT "\t\t\t\ti = 2 ! negative ion-neutral collision\n";
				print OUT "\t\t\telseif (j .eq. $n_flux_out_pos) then\n";
				print OUT "\t\t\t\ti = 3 ! positive ion-neutral collision\n";
				print OUT "\t\t\telse\n";
				print OUT "\t\t\t\ti = 4 ! negative ion-positive ion recombination\n";
			}elsif($include_negative_ion_terms){
				print OUT "\t\t\tif ((j .eq. $n_flux_out)) then\n";
				print OUT "\t\t\t\ti = 1 ! neutral-neutral collision\n";
				print OUT "\t\t\telse\n";
				print OUT "\t\t\t\ti = 2 ! negative ion-neutral collision\n";
			}else{
				print OUT "\t\t\tif ((j .eq. $n_flux_out)) then\n";
				print OUT "\t\t\t\ti = 1 ! neutral-neutral collision\n";
				print OUT "\t\t\telse\n";
				print OUT "\t\t\t\ti = 3 ! positive ion-neutral collision\n";
			}
			print OUT "\t\t\tend if\n";
			print OUT "\t\t\tj_add = coef_quad(ind_quad_form(j,n),ind_quad_form(j,n+1),j)*c(ind_quad_form(j,n))*c(ind_quad_form(j,n+1))\n";
			print OUT "\t\t\tj_tot = j_tot + j_add\n";
			print OUT "\t\t\tj_by_charge(i) = j_by_charge(i) + j_add\n";
			print OUT "\t\t\tj_by_cluster(ind_quad_form(j,n)) = j_by_cluster(ind_quad_form(j,n)) + j_add\n";
			print OUT "\t\t\tj_by_cluster(ind_quad_form(j,n+1)) = j_by_cluster(ind_quad_form(j,n+1)) + j_add\n";
			print OUT "\t\t\tj_all(ind_quad_form(j,n),i) = j_all(ind_quad_form(j,n),i) + j_add\n";
			print OUT "\t\t\tj_all(ind_quad_form(j,n+1),i) = j_all(ind_quad_form(j,n+1),i) + j_add\n";
			print OUT "\t\tend do\n";
			print OUT "\tend do\n";
		}else{
			print OUT "\tdo n=1, 2*ind_quad_form($n_flux_out,0)-1, 2 ! loop over all collisions leading out of the system\n";
			print OUT "\t\tj_add = coef_quad(ind_quad_form($n_flux_out,n),ind_quad_form($n_flux_out,n+1),$n_flux_out)*c(ind_quad_form($n_flux_out,n))*c(ind_quad_form($n_flux_out,n+1))\n";
			print OUT "\t\tj_tot = j_tot + j_add\n";
			print OUT "\t\tj_by_charge(1) = j_by_charge(1) + j_add\n";
			print OUT "\t\tj_by_cluster(ind_quad_form($n_flux_out,n)) = j_by_cluster(ind_quad_form($n_flux_out,n)) + j_add\n";
			print OUT "\t\tj_by_cluster(ind_quad_form($n_flux_out,n+1)) = j_by_cluster(ind_quad_form($n_flux_out,n+1)) + j_add\n";
			print OUT "\t\tj_all(ind_quad_form($n_flux_out,n),1) = j_all(ind_quad_form($n_flux_out,n),1) + j_add\n";
			print OUT "\t\tj_all(ind_quad_form($n_flux_out,n+1),1) = j_all(ind_quad_form($n_flux_out,n+1),1) + j_add\n";
			print OUT "\tend do\n";
		}
	}
	print OUT "\nend subroutine formation$append_to_subroutine_names\n\n!-----------------------------------------------------------\n\n";
}

################################################################################
################################################################################
##      8.7 Formation rate for Matlab                                         ##
################################################################################
################################################################################

if(!$lfortran){
	open(JFILE, ">$j_file_name") || die "\nCould not open formation rate file $j_file_name\n\n";
	
	if(!$variable_temp){
		print JFILE "function J = $j_function_name(c,K)\n\n";
		print JFILE "c = c*1e6; % converting from cm^-3 to m^-3\n\n";
		print JFILE "if nargin==1\n";
		print JFILE "\tget_coll$append_to_file_names; % reading collision rates\n";
		print JFILE "end\n";
	}else{
		print JFILE "function J = $j_function_name(c,temperature,K)\n\n";
		print JFILE "c = c*1e6; % converting from cm^-3 to m^-3\n\n";
		print JFILE "if nargin==2\n";
		print JFILE "\tK = get_coll$append_to_file_names(temperature); % reading collision rates\n";
		print JFILE "end\n";
	}
	
	print JFILE "\nJ = ";
	$line_temp = "";
	$lfound = 0;
	for($j=$n_flux_out;$j<=$n_flux_out_max;$j++){
		for($i=1;$i<=2*$ind_quad_form[$j][0];$i+=2){
			$iclus = $ind_quad_form[$j][$i];
			$jclus = $ind_quad_form[$j][$i+1];
			if(exists $coef_quad_form{"$iclus,$jclus,$j"}){
				$line_temp = "$line_temp+" . $coef_quad_form{"$iclus,$jclus,$j"} . "*c($iclus)*c($jclus) ";
				$lfound = 1;
			}else{
				$line_temp = "$line_temp+" . $coef_quad{"$iclus,$jclus,$j"} . "*c($iclus)*c($jclus) ";
				$lfound = 1;
			}
			if(($i+1)==6){
				print JFILE "$line_temp";
				$line_temp = "";
			}elsif(($i+1)%6==0){
				print JFILE "...\n\t$line_temp";
				$line_temp = "";
			}
		}
	}
	if(!$lfound){
		print JFILE "0;";
	}else{
		print JFILE "$line_temp;";
	}
	
	print JFILE "\n\nJ = J*1e-6; % converting from m^-3s^-1 to cm^-3s^-1\n\nend\n";
	close(JFILE);
}

################################################################################
################################################################################
##      8.8 Arrays for collisions, evaporations and losses for Fortran        ##
################################################################################
################################################################################

# First see if some derivatives are not used, i.e. the arrays don't need to be printed out
for($iclus=1;$iclus <= $size_C_out;$iclus++){
	$no_derivative[$iclus]=0;
}
for($icol=0;$icol<=$#no_derivative_clus;$icol++){
	($iclus,$nwaters)=&get_corresponding_dry($no_derivative_clus[$icol]);
	if($iclus eq ''){
		die "Cluster $no_derivative_clus[$icol] for which the derivative should not be printed out does not exist.\n";
	}
	$no_derivative[$iclus]=1;
}

if($lfortran && !$lloop){
	print OUT "subroutine initialize_parameters$append_to_subroutine_names(coef_quad,coef_lin,ind_quad_loss,ind_quad_form,&\n\t\t\t&ind_lin_loss,ind_lin_form,ind_quad_loss_extra,ind_quad_form_extra,&\n\t\t&\tsource,isconst,fitted";
	if($l_mcmc_any || $variable_temp || $variable_cs || $variable_ion_source){
		print OUT ",coef";
	}
	print OUT ",ipar)\n";
	print OUT "\n\tuse monomer_settings, only : sources_and_constants\n\n\timplicit none\n";
	print OUT "\tinteger, parameter :: neqn=$size_C_out, nclust=$max_cluster_number\n";
	print OUT "\tlogical :: isconst($size_C_out)\n";
	print OUT "\treal(kind(1.d0)) :: coef_quad(nclust,nclust,$max_n_flux),coef_lin($max_n_flux,$max_n_flux,nclust)\n";
	print OUT "\treal(kind(1.d0)) :: source($size_C_out)\n";
	print OUT "\tinteger :: ind_quad_loss($size_C_out,0:$ind_quad_loss[0]),ind_quad_form($size_C_out,0:$ind_quad_form[0]),ind_quad_loss_extra($size_C_out,0:$ind_quad_loss_extra[0]),ind_quad_form_extra($size_C_out,0:$ind_quad_form_extra[0])\n";
	print OUT "\tinteger :: ind_lin_loss($size_C_out,0:$ind_lin_loss[0]),ind_lin_form($size_C_out,0:$ind_lin_form[0])\n";
	print OUT "\tinteger :: fitted(nclust,0:nclust)\n";	
	if($l_mcmc_any){
		print OUT "\treal(kind(1.d0)) :: coef($n_mcmc_coefs)\n";
	}elsif($variable_temp || $variable_cs || $variable_ion_source){
		print OUT "\treal(kind(1.d0)) :: coef(" . ($variable_temp+$variable_cs+2*$variable_ion_source) . ")\n";
	}
	print OUT "\tinteger :: ipar(4)\n";
	
	print OUT "\n\tcall sources_and_constants$append_to_subroutine_names(source,isconst,fitted,ipar)\n";
	if($include_generic_neg && $include_generic_pos){
		if($generic_positive_ions_corr ne ''){
			print OUT "\n\t! setting positive charger ion concentration to give charge balance\n";
			print OUT "\tisconst($n_generic_pos) = .true.\n\n";
		}elsif($generic_negative_ions_corr ne ''){
			print OUT "\n\t! setting negative charger ion concentration to give charge balance\n";
			print OUT "\tisconst($n_generic_neg) = .true.\n\n";
		}elsif($generic_positive_ions ne ''){
			print OUT "\n\t! setting positive charger ion concentration to give charge balance\n";
			print OUT "\tisconst($n_generic_pos) = .true.\n\n";
		}elsif($generic_negative_ions ne ''){
			print OUT "\n\t! setting negative charger ion concentration to give charge balance\n";
			print OUT "\tisconst($n_generic_neg) = .true.\n\n";
		}
	}
	if($variable_ion_source){
		print OUT "\tisconst($n_generic_neg) = .false.\n";
		print OUT "\tsource($n_generic_neg) = coef(" . ($variable_temp+$variable_cs+1) . ")\n";
		print OUT "\tisconst($n_generic_pos) = .false.\n";
		print OUT "\tsource($n_generic_pos) = coef(" . ($variable_temp+$variable_cs+2) . ")\n";
	}
	print OUT "\n\tind_quad_loss = 0\n\tind_quad_form = 0\n\tind_lin_loss = 0\n\tind_lin_form = 0\n\tind_quad_loss_extra = 0\n\tind_quad_form_extra = 0\n\n";
	print OUT "\t! ind_quad_loss(i,0:2n) = (/n,j1,k1,j2,k2,...,jn,kn/) gives the cluster numbers for all the collisions\n";
	print OUT "\t!  i + j1 -> k1 etc. through which cluster i is lost\n\n";
	print OUT "\t! ind_quad_form(k,0:2n) = (/n,i1,j1,i2,j2,...,in,jn/) gives the cluster numbers for all the collisions\n";
	print OUT "\t!  i1 + j1 -> k etc. through which cluster k is formed\n\n";
	print OUT "\t! ind_lin_loss(k,0:2n) = (/n,i1,j1,i2,j2,...,lossn,lossn/) gives the cluster numbers for all the evaporations\n";
	print OUT "\t!  k -> i1 + j1 etc. and other losses k -> wall etc. through which cluster k is lost\n\n";
	print OUT "\t! ind_lin_form(i,0:2n) = (/n,j1,k1,j2,k2,...,jn,kn/) gives the cluster numbers for all the evaporations\n";
	print OUT "\t!  k1 -> i + j1 etc. through which cluster i is formed\n\n";
	print OUT "\t! ind_quad_loss_extra(i,0:3n) = (/n,j1,k1,c1,...,jn,kn,cn/) gives the cluster numbers and coefficients \n";
	print OUT "\t!  i + j1 -> c1*k1 etc. for additional collision products k (e.g. monomers from the boundary)\n";
	print OUT "\t!  in collisions where cluster i is lost\n\n";
	print OUT "\t! ind_quad_form_extra(k,0:2n) = (/n,i1,j1,c1,...,in,jn,cn/) gives the cluster numbers and coefficients \n";
	print OUT "\t!  i1 + j1 -> c1*k etc. for additional ways of forming k (e.g. as a monomer from the boundary)\n\n";
	
	for($iclus=1;$iclus <= $size_C_out;$iclus++){
		if($no_derivative[$iclus]==1){
			print OUT "\n\t! Cluster $iclus: $label[$iclus]\n";
			print OUT "\t! Derivative not used and concentration set to be constant\n";
			print OUT "\tisconst($iclus) = .true.\n";
			next;
		}
		if($ind_quad_loss[$iclus][0]+$ind_quad_form[$iclus][0]+$ind_lin_loss[$iclus][0]+$ind_lin_form[$iclus][0]>0){
			print OUT "\n\t! Cluster $iclus: $label[$iclus]\n";
		}
		if($ind_quad_loss[$iclus][0]>0){
			print OUT "\tind_quad_loss($iclus,0:" . $ind_quad_loss[$iclus][0]*2 . ") = (/ $ind_quad_loss[$iclus][0]";
			$i = 1;
			for($icol=1;$icol <= $ind_quad_loss[$iclus][0];$icol++){
				if($i==$max_terms_per_line){
					print OUT "&\n\t\t\t&";
					$i = 0;
				}
				$i += 1;
				print OUT ", $ind_quad_loss[$iclus][2*$icol-1],$ind_quad_loss[$iclus][2*$icol]";
			}
			print OUT " /)\n";
		}
		if($ind_quad_loss_extra[$iclus][0]>0){
			print OUT "\tind_quad_loss_extra($iclus,0:" . $ind_quad_loss_extra[$iclus][0]*3 . ") = (/ $ind_quad_loss_extra[$iclus][0]";
			$i = 1;
			for($icol=1;$icol <= $ind_quad_loss_extra[$iclus][0];$icol++){
				if($i==$max_terms_per_line){
					print OUT "&\n\t\t\t&";
					$i = 0;
				}
				$i += 1;
				print OUT ", $ind_quad_loss_extra[$iclus][3*$icol-2],$ind_quad_loss_extra[$iclus][3*$icol-1],$ind_quad_loss_extra[$iclus][3*$icol]";
			}
			print OUT " /)\n";
		}
		if($ind_quad_form[$iclus][0]>0){
			print OUT "\tind_quad_form($iclus,0:" . $ind_quad_form[$iclus][0]*2 . ") = (/ $ind_quad_form[$iclus][0]";
			$i = 1;
			for($icol=1;$icol <= $ind_quad_form[$iclus][0];$icol++){
				if($i==$max_terms_per_line){
					print OUT "&\n\t\t\t&";
					$i = 0;
				}
				$i += 1;
				print OUT ", $ind_quad_form[$iclus][2*$icol-1],$ind_quad_form[$iclus][2*$icol]";
			}
			print OUT " /)\n";
		}
		if($ind_quad_form_extra[$iclus][0]>0){
			print OUT "\tind_quad_form_extra($iclus,0:" . $ind_quad_form_extra[$iclus][0]*3 . ") = (/ $ind_quad_form_extra[$iclus][0]";
			$i = 1;
			for($icol=1;$icol <= $ind_quad_form_extra[$iclus][0];$icol++){
				if($i==$max_terms_per_line){
					print OUT "&\n\t\t\t&";
					$i = 0;
				}
				$i += 1;
				print OUT ", $ind_quad_form_extra[$iclus][3*$icol-2],$ind_quad_form_extra[$iclus][3*$icol-1],$ind_quad_form_extra[$iclus][3*$icol]";
			}
			print OUT " /)\n";
		}
		if($ind_lin_loss[$iclus][0]>0){
			print OUT "\tind_lin_loss($iclus,0:" . $ind_lin_loss[$iclus][0]*2 . ") = (/ $ind_lin_loss[$iclus][0]";
			$i = 1;
			for($icol=1;$icol <= $ind_lin_loss[$iclus][0];$icol++){
				if($i==$max_terms_per_line){
					print OUT "&\n\t\t\t&";
					$i = 0;
				}
				$i += 1;
				print OUT ", $ind_lin_loss[$iclus][2*$icol-1],$ind_lin_loss[$iclus][2*$icol]";
			}
			print OUT " /)\n";
		}
		if($ind_lin_form[$iclus][0]>0){
			print OUT "\tind_lin_form($iclus,0:" . $ind_lin_form[$iclus][0]*2 . ") = (/ $ind_lin_form[$iclus][0]";
			$i = 1;
			for($icol=1;$icol <= $ind_lin_form[$iclus][0];$icol++){
				if($i==$max_terms_per_line){
					print OUT "&\n\t\t\t&";
					$i = 0;
				}
				$i += 1;
				print OUT ", $ind_lin_form[$iclus][2*$icol-1],$ind_lin_form[$iclus][2*$icol]";
			}
			print OUT " /)\n";
		}
	}
	for($iclus=$size_C_out+1;$iclus <= $max_n_flux;$iclus++){
		if($ind_quad_form[$iclus][0]>0){
			print OUT "\n\t! $label[$iclus]\n";
			print OUT "\tind_quad_form($iclus,0:" . $ind_quad_form[$iclus][0]*2 . ") = (/ $ind_quad_form[$iclus][0]";
			$i = 1;
			for($icol=1;$icol <= $ind_quad_form[$iclus][0];$icol++){
				if($i==$max_terms_per_line){
					print OUT "&\n\t\t\t&";
					$i = 0;
				}
				$i += 1;
				print OUT ", $ind_quad_form[$iclus][2*$icol-1],$ind_quad_form[$iclus][2*$icol]";
			}
			print OUT " /)\n";
		}
	}
	print OUT "\n\tcall get_rate_coefs$append_to_subroutine_names(coef_quad,coef_lin";
	if($l_mcmc_any || $variable_temp || $variable_cs || $variable_ion_source){
		print OUT ",coef";
	}
	print OUT ")\n";
	print OUT "\nend subroutine initialize_parameters$append_to_subroutine_names\n\n!-----------------------------------------------------------\n\n";
	
	# Figuring out the input/output parameters of get_fcr -- these are used several times
	# These are used in the subroutine call ...
	$parameters_in_get_fcr_1 = '';
	$parameters_in_get_fcr_2 = '';
	# ... and these are for variable declarations
	$parameters_in_get_fcr_3 = '';
	$parameters_in_get_fcr_4 = '';
	if($include_ion_terms){
		if(!$disable_coag_sinks){
			$parameters_in_get_fcr_1 = "fcs";
			$parameters_in_get_fcr_3 = "fcs";
		}
		if(!$disable_wall_terms && !$l_wl_mcmc){
			if($parameters_in_get_fcr_1 ne ''){
				$parameters_in_get_fcr_1 = "$parameters_in_get_fcr_1,fwl";
				$parameters_in_get_fcr_3 = "$parameters_in_get_fcr_3,fwl";
			}else{
				$parameters_in_get_fcr_1 = "fwl";
				$parameters_in_get_fcr_3 = "fwl";
			}
		}
		$parameters_in_get_fcr_2 = $parameters_in_get_fcr_1;
		$parameters_in_get_fcr_4 = $parameters_in_get_fcr_3;
	}

	print OUT "subroutine get_rate_coefs$append_to_subroutine_names(coef_quad,coef_lin";
	if($l_mcmc_any || $variable_temp || $variable_cs || $variable_ion_source){
		print OUT ",coef";
	}
	print OUT ")\n";
	print OUT "\timplicit none\n";
	print OUT "\treal(kind(1.d0)) :: coef_quad($max_cluster_number,$max_cluster_number,$max_n_flux),coef_lin($max_n_flux,$max_n_flux,$max_cluster_number)\n";
	print OUT "\treal(kind(1.d0)) :: K($max_cluster_number,$max_cluster_number),E($max_cluster_number,$max_cluster_number)";
	if(!$disable_coag_sinks){
		if($one_cs_for_all){
			print OUT ",cs";
		}else{
			print OUT ",cs($max_cluster_number)";
		}
	}
	if(!$disable_wall_terms){
		if($one_wl_for_all){
			print OUT ",wl";
		}else{
			print OUT ",wl($max_cluster_number)";
		}
	}
	if(!$disable_dilution){
		print OUT ", dil";
	}
	if($parameters_in_get_fcr_3 ne ''){
		print OUT ",$parameters_in_get_fcr_3";
	}
	print OUT "\n";
	if($l_mcmc_any){
		print OUT "\treal(kind(1.d0)) :: coef($n_mcmc_coefs)\n";
	}elsif($variable_temp || $variable_cs || $variable_ion_source){
		print OUT "\treal(kind(1.d0)) :: coef(" . ($variable_temp+$variable_cs+2*$variable_ion_source) . ")\n";
	}

	print OUT "\n\tcoef_quad = 0.d0\n\tcoef_lin = 0.d0\n";

	if($l_mcmc_any){
		print OUT "\tcall get_coll$append_to_subroutine_names(K,coef)\n";
	}elsif( $variable_temp){
		print OUT "\tcall get_coll$append_to_subroutine_names(K,coef(1))\n";
	}else{
		print OUT "\tcall get_coll$append_to_subroutine_names(K)\n";
	}

	if($parameters_in_get_fcr_1 ne ''){
		print OUT "\tcall get_fcr$append_to_subroutine_names($parameters_in_get_fcr_1)\n";
	}

	if(!$disable_evap){
		if($l_mcmc_any){
			print OUT "\tcall get_evap$append_to_subroutine_names(E,coef)\n";
		}elsif( $variable_temp){
			print OUT "\tcall get_evap$append_to_subroutine_names(E,K,coef(1))\n";
		}else{
			print OUT "\tcall get_evap$append_to_subroutine_names(E)\n";
		}
	}
	if($l_use_get_losses){
		$line = '';
		if(!$disable_coag_sinks && (!$one_cs_for_all || !$variable_cs)){
			$line = 'cs';
		}
		if(!$disable_wall_terms || ($l_mcmc_any && $ncoef_wl>0)){
			if($line ne ''){
				$line = "$line,wl";
			}else{
				$line = 'wl';
			}
		}
		if(!$disable_dilution){
			if($line ne ''){
				$line = "$line,dil";
			}else{
				$line = 'dil';
			}
		}
		print OUT "\tcall get_losses$append_to_subroutine_names($line)\n";
	}
	if($variable_cs){
		if($one_cs_for_all){
			print OUT "\tcs = coef(" . ($variable_temp+1) . ")\n";
		}else{
			print OUT "\tcs = coef(" . ($variable_temp+1) . ")*cs\n";
		}
	}
	if($l_wl_mcmc && $ncoef_wl==1){
		print OUT "\twl = coef($n_mcmc_coefs)\t\t\t\t! wall loss\n";
	}
	print OUT "\n";
	for($iclus=1;$iclus <= $max_n_flux;$iclus++){
		$line = '';
		for($jclus=1;$jclus <= $iclus;$jclus++){
			for($kclus=$max_n_flux;$kclus >=1 ;$kclus-=1){
				if(exists $coef_quad{"$iclus,$jclus,$kclus"}){
					if(exists $l_quad_bound{"$iclus,$jclus,$kclus"}){
						$line = $line . "\tcoef_quad($iclus,$jclus,$kclus) = " . $coef_quad{"$iclus,$jclus,$kclus"} . "\t! $label[$iclus] + $label[$jclus] -> boundary -> $label[$kclus]\n";
					}elsif($use_coag_to_outgoing && $kclus==$n_flux_out){
						$line = $line . "\tcoef_quad($iclus,$jclus,$kclus) = " . $coef_quad{"$iclus,$jclus,$kclus"} . "\t! $label[$iclus] + $label[$jclus] -> out -> cs_cluster\n";
					}elsif($use_coag_to_outgoing && $kclus==$n_flux_out_neg){
						$line = $line . "\tcoef_quad($iclus,$jclus,$kclus) = " . $coef_quad{"$iclus,$jclus,$kclus"} . "\t! $label[$iclus] + $label[$jclus] -> out -> cs_cluster_neg\n";
					}elsif($use_coag_to_outgoing && $kclus==$n_flux_out_pos){
						$line = $line . "\tcoef_quad($iclus,$jclus,$kclus) = " . $coef_quad{"$iclus,$jclus,$kclus"} . "\t! $label[$iclus] + $label[$jclus] -> out -> cs_cluster_pos\n";
					}else{
						$line = $line . "\tcoef_quad($iclus,$jclus,$kclus) = " . $coef_quad{"$iclus,$jclus,$kclus"} . "\t! $label[$iclus] + $label[$jclus] -> $label[$kclus]\n";
					}
					if($iclus != $jclus){
						$line = $line . "\tcoef_quad($jclus,$iclus,$kclus) = " . $coef_quad{"$jclus,$iclus,$kclus"} . "\n";
					}
				}
				if(exists $coef_lin{"$iclus,$jclus,$kclus"}){
					if($iclus > $max_cluster_number){
						$line = $line . "\tcoef_lin($iclus,$jclus,$kclus) = " . $coef_lin{"$iclus,$jclus,$kclus"} . "\t! $label[$kclus] -> $label[$iclus]\n";
					}else{
						$line = $line . "\tcoef_lin($iclus,$jclus,$kclus) = " . $coef_lin{"$iclus,$jclus,$kclus"} . "\t! $label[$kclus] -> $label[$iclus] + $label[$jclus]\n";
					}
					if($iclus != $jclus){
						$line = $line . "\tcoef_lin($jclus,$iclus,$kclus) = " . $coef_lin{"$jclus,$iclus,$kclus"} . "\n";
					}
				}
			}
		}
		if($line ne ''){
			print OUT "\n$line";
		}
	}
	print OUT "\nend subroutine get_rate_coefs$append_to_subroutine_names\n\n!-----------------------------------------------------------\n\n";
}	

################################################################################
################################################################################
##      8.9 Masses, radii and molecular content for the loop version
################################################################################
################################################################################

if($lloop){
	if($nmol_types_used>1){
		print OUT "subroutine get_masses_and_radii$append_to_subroutine_names(m,r,i_ind)\n";
	}else{
		print OUT "subroutine get_masses_and_radii$append_to_subroutine_names(m,r,i)\n";
	}
	print OUT "\timplicit none\n";
	print OUT "\treal(kind(1.d0)), parameter :: u2kg = 1.660538782d-27, kB=1.3806504d-23, pi=4.d0*atan(1.d0)\n";
	
	$str_temp = '';
	$str_temp2 = '';
	$string = '';
	$string1 = '';
	$string2 = '';
	$i = 0;
	for($imol=0;$imol<$nmol_types;$imol++){
		if($n_max_type[$imol]>0){
			$i ++;
			$temp_val = sprintf('%.14e',$mass[$imol]/$mass_conv);
			$temp_val =~ s/e/d/g;
			$str_temp .= "\treal(kind(1.d0)), parameter :: m$i = $temp_val\t\t! in g/mol\n";
			$temp_val = sprintf('%.14e',$mass[$imol]/$density[$imol]);
			$temp_val =~ s/e/d/g;
			$str_temp2 .= "\treal(kind(1.d0)), parameter :: monvol$i = $temp_val\t\t! in m^3\n";
			$string .= "i$i+";
			$string1 .= "i_ind($i)*m$i+";
			$string2 .= "i_ind($i)*monvol$i+";
		}
	}
	print OUT "$str_temp";
	print OUT "$str_temp2";
	print OUT "\treal(kind(1.d0)) :: r, m\n";
	if($nmol_types_used>1){
		print OUT "\tinteger :: i_ind($nmol_types_used)\n\n";
		$string1 =~ s/\+$//;
		$string2 =~ s/\+$//;
		$string1 = "($string1)";
		$string2 = "($string2)";
	}else{
		print OUT "\tinteger :: i\n\n";
		$string1 = "i*m1";
		$string2 = "i*monvol1";
	}
	print OUT "\tm = $string1\n";
	print OUT "\tr = (3.d0/4.d0/pi*$string2)**(1.d0/3.d0)\n\n";
	print OUT "\nend subroutine get_masses_and_radii$append_to_subroutine_names\n\n!-----------------------------------------------------------\n\n";
	if($nmol_types_used>1){
		print OUT "subroutine get_molecule_numbers$append_to_subroutine_names(indices)\n";
		print OUT "\timplicit none\n";
		print OUT "\tinteger, parameter :: nclust = $max_cluster_number\n";
		print OUT "\tinteger :: n, $nmol_ranges_ind, indices(nclust,$nmol_types_used)";
		if($lloop_max_ratio){
			print OUT ", n_pure($nmol_ranges_array[$nmol_max_ratio_basis])\n\n";
			print OUT "\tcall pure_clust_indices$append_to_subroutine_names(n_pure)\n\n";
		}else{
			print OUT "\n\n";
		}
		if(!$lloop_max_ratio){
			print OUT "\tn = 0\n";
			$str_temp = "\t";
			$str_temp2 = '';
			# Printing loops over all molecule numbers
			for($imol=1;$imol<=$nmol_types_used;$imol++){
				# If there are many molecule types, the molecule number loops start from 0 ...
				print OUT $str_temp."do i$imol = 0,$nmol_ranges_array[$imol]\n";
				$str_temp2 = $str_temp."end do\n".$str_temp2;
				$str_temp .= "\t";
			}
			# ... but they can't all be 0 ...
			$string =~ s/\+$//;
			print OUT $str_temp."if ($string>0) then\n";
			$str_temp2 = $str_temp."end if\n".$str_temp2;
			$str_temp .= "\t";
			print OUT $str_temp."n = n+1\n";
			$string =~ s/\+/,/g;
			print OUT $str_temp."indices(n,:) = (/$string/)\n";
			print OUT $str_temp2;
		}else{
			$string =~ s/\+$//;
			$string =~ s/\+/,/g;
			print OUT "\tn = 1\n";
			print OUT "\ti$nmol_max_ratio_basis = 0\n";
			print OUT "\ti$nmol_not_basis = 1\n";
			print OUT "\tindices(n,:) = (/$string/)\n";
			print OUT "\tdo i$nmol_max_ratio_basis = 1,$nmol_ranges_array[$nmol_max_ratio_basis]\n";
			print OUT "\t\tdo i$nmol_not_basis = 0,$max_wrt_mon*i$nmol_max_ratio_basis\n";
			print OUT "\t\t\tn = n+1\n";
			print OUT "\t\t\tindices(n,:) = (/$string/)\n";
			print OUT "\t\tend do\n";
			print OUT "\tend do\n";
		}
		print OUT "\nend subroutine get_molecule_numbers$append_to_subroutine_names\n\n!-----------------------------------------------------------\n\n";
		if($lloop_max_ratio){
			# subroutine for determining the cluster numbers of n-mers of pure "basis" molecule type
			# (for easily finding the cluster numbers in the equation loops)
			print OUT "subroutine pure_clust_indices$append_to_subroutine_names(n_pure)\n\timplicit none\n";
			print OUT "\tinteger :: i, n_pure($nmol_ranges_array[$nmol_max_ratio_basis])\n\n";			
			print OUT "\tn_pure(1) = 2\n";
			print OUT "\tdo i=2,$nmol_ranges_array[$nmol_max_ratio_basis]\n";
			print OUT "\t\tn_pure(i) = n_pure(i-1)+$max_wrt_mon*(i-1)+1\n";
			print OUT "\tend do\n";
			print OUT "\nend subroutine pure_clust_indices$append_to_subroutine_names\n\n!-----------------------------------------------------------\n\n";
		}
	}
}


################################################################################
################################################################################
##   9. Loss coefficients                                                     ##
################################################################################
################################################################################

if($lfortran){
	if($l_use_get_losses){
		$lines{'inp'} = '';
		$lines{'def'} = '';
		if(!$disable_coag_sinks){
			if(!$one_cs_for_all){
				$lines{'inp'} = 'cs';
				$lines{'def'} = "cs($max_cluster_number)";
			}elsif(!$variable_cs){
				$lines{'inp'} = 'cs';
				$lines{'def'} = "cs";
			}
		}
		if(!$disable_wall_terms || ($l_mcmc_any && $ncoef_wl>0)){
			if($lines{'inp'} ne ''){
				$lines{'inp'} .= ",wl";
				if($one_wl_for_all){
					$lines{'def'} .= ",wl";
				}else{
					$lines{'def'} .= ",wl($max_cluster_number)";
				}
			}else{
				$lines{'inp'} = 'wl';
				if(!$one_wl_for_all){
					$lines{'def'} = "wl($max_cluster_number)";
				}else{
					$lines{'def'} = "wl";
				}
			}
		}
		if(!$disable_dilution){
			if($lines{'inp'} ne ''){
				$lines{'inp'} .= ", dil";
				$lines{'def'} .= ", dil";
			}else{
				$lines{'inp'} = ' dil';
				$lines{'def'} = ' dil';
			}
		}
		if(!$lloop){
			print OUT "subroutine get_losses$append_to_subroutine_names($lines{'inp'})\n";
			print OUT "\timplicit none\n\treal(kind(1.d0)) :: $lines{'def'}\n\n";
		}else{
			print OUT "subroutine get_losses$append_to_subroutine_names(loss)\n";
			if(!$disable_coag_sinks && $coag_terms =~ m/^bg_loss$/i){
				# Get cbg from a module by default
				#if($lloop_j_tests){
				print OUT "use shared_input, only : cbg => bg_concentration\n";
				print "\nUsing module shared_input in the get_losses subroutine to get cbg\n";
				#}
			}
			print OUT "\timplicit none\n\treal(kind(1.d0)) :: $lines{'def'}, loss";
			if((!$disable_coag_sinks && !$one_cs_for_all) || ((!$disable_wall_terms || ($l_mcmc_any && $ncoef_wl>0)) && !$one_wl_for_all)){
				print OUT "($max_cluster_number), m, r";
				if(!$disable_wall_terms && ($wall_terms =~ m/^cloud4_exel_file$/i || $wall_terms =~ m/^cloud4_Andreas$/i)){
					print OUT ", d, sc";
				}
				if(!$disable_coag_sinks){
					if($coag_terms =~ m/^exp_loss$/i){
						print OUT ", r_ref";
					}elsif($coag_terms =~ m/^bg_loss$/i){
						#if($lloop_j_tests){
							print OUT ", vth, D, mbg, rbg, vthbg, Dbg, Kn";
						#}else{
							#print OUT ", vth, D, mbg, rbg, vthbg, Dbg, cbg, Kn";
						#}
					}
				}
			}
			print OUT "\n\tinteger, parameter :: nclust = $max_cluster_number\n\tinteger :: i\n";
			if(!$lloop_diag && $nmol_types_used>1){
				print OUT "\tinteger :: indices(nclust,$nmol_types_used)\n\n";
				print OUT "\tcall get_molecule_numbers$append_to_subroutine_names(indices)\n";
			}
			print OUT "\n";
		}
		if(!$disable_dilution){
			print OUT "\tdil = $dil_str\n\n";
		}
	}
}else{
	open(CS, ">get_cs$append_to_file_names.m");
	open(WL, ">get_wl$append_to_file_names.m");
}

################################################################################
################################################################################
##      9.1 Coagulation sink                                                  ##
################################################################################
################################################################################

if(! $lfortran){
	print CS "% Coagulation sink terms\n\n";
	if($include_ion_terms){
		print CS "% enhancement factor for ions\nfcs = $fcs;\n\n";
	}
}
if(!$disable_coag_sinks){
	if($lfortran && $l_use_get_losses){
		print OUT "\t! coagulation sink\n\n";
		if($variable_cs){
			if($one_cs_for_all){
				print OUT "\t! variable cs: the value will be given as input\n\n";
			}else{
				print OUT "\t! variable cs: the following values only give the size dependence of the sink,\n";
				print OUT "\t! and will be multiplied by the size-independent factor given as input\n\n";
			}
		}
	}
	# array for the coagulation loss coefficients of the hydrates
	@cs_coef_hydr=();

################################################################################
################################################################################
##         9.11 Coagulation sink parameterizations                            ##
################################################################################
################################################################################

################################################################################
################################################################################
##            9.111 Size-dependent sink of the form A*(d/dref)^B,             ##
##                  where d is the diameter                                   ##
################################################################################
################################################################################
	
	if($coag_terms =~ m/^exp_loss$/i){
		
		if(!$lloop){
			# Find the reference size
			if($exp_loss_ref_cluster ne ''){
				($iclus,$nwaters)=&get_corresponding_dry($exp_loss_ref_cluster);
				if($iclus ne ''){
					$r1=$clust_radius[$iclus][$nwaters];
				}else{
					die "\nCan't calculate loss coefficients relative to $exp_loss_ref_cluster, if you don't have it in the system\n\n";
				}
			}else{
				# Use the monomer (the default monomer if there are several monomers) as reference cluster
	
				# find the monomer(s)
				$lfound=0;
				$lfound_several=0;
				$lfound_default=0;
				for($iclus=1;$iclus<=$max_cluster_number;$iclus++){
					if($lmonomer[$iclus]){
						if($lfound==1){
							$lfound_several=1;
						}
						$lfound=1;
						$radius=$clust_radius[$iclus][0];
						if($label[$iclus] eq "1$default_monomer"){
							$lfound_default=1;
							$r1=$radius;
						}
					}
				}
				if($lfound==0){
					die "\nCan't calculate loss coefficients relative to that of the monomer, if you don't have monomers in the system\n\n";
				}elsif($lfound_several==0){
					$r1=$radius;
				}elsif($lfound_several==1 && $lfound_default==0){
					die "\nDid not find the default monomer $default_monomer among several monomer species for calculating the loss coefficients -> you need to specify the reference cluster as eg. --exp_loss_ref_cluster $label[1]\n";
				}
			}
	
			for($iclus=1;$iclus<=$max_cluster_number;$iclus++){
				$this_clus=$label[$iclus];						
				# loss coefficients for the dry cluster and possible hydrates
				for($ihydr=0;$ihydr<=$#{ $distr[$iclus] };$ihydr++){
					if(defined $distr[$iclus][$ihydr]){
						if($ihydr==0){
							$temp_label=$this_clus;
						}else{
							$temp_label=$this_clus.$ihydr.$name_water;
						}
						$radius=$clust_radius[$iclus][$ihydr];

						# Coagulation loss of the form A*(d/dref)^B
						$cs_term_temp = ($radius/$r1)**$exp_loss_exponent;
						if(!$variable_cs){
							# If the sink is not varied, include the reference loss rate in it
							$cs_term_temp *= $exp_loss_coefficient;
						}
						$cs_coef_hydr[$iclus][$ihydr] = $cs_term_temp;
					}
				}
			}
		}else{
			$exp_loss_coefficient = sprintf('%.14e',$exp_loss_coefficient);
			$exp_loss_coefficient =~ s/e/d/;
			$exp_loss_exponent = sprintf('%.14e',$exp_loss_exponent);
			if(!&is_pos_number($exp_loss_exponent)){
				$exp_loss_exponent = "($exp_loss_exponent)";
			}
			$exp_loss_exponent =~ s/e/d/;
			$cs_coef_hydr[0][0] = "\t\tcs(i) = $exp_loss_coefficient*(r/r_ref)**$exp_loss_exponent\n";
			
			# Find the reference size
			if($exp_loss_ref_cluster ne ''){
				($string1,$string2) = &determine_cluster_composition($exp_loss_ref_cluster);
				@temp_array=@$string1;
				$valid = 1;
				if($string2 ne '' || sum(@temp_array)>sum(@temp_array[@mol_types_used[1..$nmol_types_used]])){
					# cluster contains molecule types not included in the system
					$valid = 0;
				}elsif($lloop_max_ratio && (($temp_array[$nmol_inp_max_ratio_basis]==0 && $temp_array[$nmol_inp_not_basis]>1) || ($temp_array[$nmol_inp_max_ratio_basis]>0 && $temp_array[$nmol_inp_not_basis]>$max_wrt_mon*$temp_array[$nmol_inp_max_ratio_basis]))){
					# wrong composition
					$valid = 0;
				}else{
					for($imol=0;$imol<$nmol_types;$imol++){
						if($temp_array[$imol]>$n_max_type[$imol]){
							# cluster outside the system
							$valid = 0;
							last;
						}
					}
				}
				if(!$valid){
					die "Loss term reference cluster $exp_loss_ref_cluster is not a valid cluster\n";
				}
				if($nmol_types_used>1){
					if(!$lloop_max_ratio){
						$i = 0;
						for($imol=1;$imol<=$nmol_types_used;$imol++){
							$i = $i+$temp_array[$mol_types_used[$imol]]*$cluster_from_indices_array[$imol];
						}
					}else{
						$i = 1;
						if($temp_array[$nmol_inp_max_ratio_basis] > 0){
							for($j=1;$j<=$temp_array[$nmol_inp_max_ratio_basis];$j++){
								$i = $i+$max_wrt_mon*($j-1)+1;
							}
							$i = $i+$temp_array[$nmol_inp_not_basis];
						}
					}
				}else{
					$imol = 1;
					$i = $temp_array[$mol_types_used[$imol]];
				}
				print OUT "\t! reference radius corresponding to $exp_loss_ref_cluster\n"
			}elsif($nmol_types_used==1){
				# Use the monomer as reference cluster
				$i = 1;
				print OUT "\t! reference radius corresponding to 1$molecule_name[$mol_types_used[1]]\n"
			}else{
				# Use the default monomer as reference cluster
				$lfound_default=0;
				for($imol=1;$imol<=$nmol_types_used;$imol++){
					if($molecule_name[$mol_types_used[$imol]] eq $default_monomer){
						$lfound_default = 1;
						$i = $cluster_from_indices_array[$imol];
					}
				}
				if($lfound_default==0){
					die "\nDid not find the default monomer $default_monomer among several monomer species for calculating the loss coefficients -> you need to specify the reference cluster as eg. --exp_loss_ref_cluster 1$molecule_name[$mol_types_used[1]]\n";
				}
				print OUT "\t! reference radius corresponding to 1$default_monomer\n"
			}
			if($nmol_types_used>1){
				print OUT "\tcall get_masses_and_radii$append_to_subroutine_names(m,r_ref,indices($i,:))\n\n"
			}else{
				print OUT "\tcall get_masses_and_radii$append_to_subroutine_names(m,r_ref,$i)\n\n"
			}
		}
		
################################################################################		
##            9.112 Size-dependent sink based on a constant concentration     ##
##                  and size of background particles                          ##
################################################################################

	}elsif($coag_terms =~ m/^bg_loss$/i){
		# Loss calculated as the frequency at which the clusters collide with
		# background particles that have constant concentration, diameter and density
		
		$air_viscosity = 2.5277e-7*$temperature**0.75302; # from DMAN
		# mean free path of air, assuming standard pressure
		$temp_val = 2*$air_viscosity/(101325*sqrt(8*0.0289/($pi*$Na*$boltz*$temperature)));
		# S&P Eq. 9.73
		$temp_val2 = $boltz*$temperature/(3*$pi*$air_viscosity);
		# factor for calculating particle velocity
		$temp_val3 = 8*$boltz*$temperature/$pi;
		
		# properties of the background particles
		$mass2=$bg_density*$pi/6*$bg_diameter**3;
		$r2=$bg_diameter/2;
		$Diff2 = $temp_val2/(2*$r2)*((5+4*($temp_val/$r2)+6*($temp_val/$r2)**2+18*($temp_val/$r2)**3)/(5-($temp_val/$r2)+(8+$pi)*($temp_val/$r2)**2));
		$c2 = sqrt($temp_val3/$mass2);

		if(!$lloop){
			for($iclus=1;$iclus<=$max_cluster_number;$iclus++){
				$this_clus=$label[$iclus];						
				# loss coefficients for the dry cluster and possible hydrates
				for($ihydr=0;$ihydr<=$#{ $distr[$iclus] };$ihydr++){
					if(defined $distr[$iclus][$ihydr]){
						$mass1=$clust_mass[$iclus][$ihydr];
						$r1=$clust_radius[$iclus][$ihydr];
						# Use the parameterization by Dahneke (1983) for the collision coefficients
						# (see the 'Dahneke' option of the loop mode)
						$Diff1 = $temp_val2/(2*$r1)*((5+4*($temp_val/$r1)+6*($temp_val/$r1)**2+18*($temp_val/$r1)**3)/(5-($temp_val/$r1)+(8+$pi)*($temp_val/$r1)**2));
						$c1 = sqrt($temp_val3/$mass1);
						$Kn = 2*($Diff1+$Diff2)/(sqrt($c1**2+$c2**2)*($r1+$r2));
						$cs_term_temp = 4*$pi*($r1+$r2)*($Diff1+$Diff2)*(1+$Kn)/(1+2*$Kn*(1+$Kn));
						if(!$variable_cs){
							# If the sink is not varied, include the bg concentration in it
							$cs_term_temp *= $bg_concentration;
						}
						$cs_coef_hydr[$iclus][$ihydr] = $cs_term_temp;
					}
				}
			}
		}else{
			# print out the mass and radius of the background particles
			$string1 = sprintf('%.14e',$mass2/$mass_conv); # kg -> g/mol, because the masses in get_masses_and_radii of the loop version are g/mol
			$string1 =~ s/e/d/;
			print OUT "\tmbg = $string1\n";
			$string1 = sprintf('%.14e',$r2);
			$string1 =~ s/e/d/;
			print OUT "\trbg = $string1\n";
			#if(!$lloop_j_tests){
				#$string1 = sprintf('%.14e',$bg_concentration);
				#$string1 =~ s/e/d/;
				#print OUT "\tcbg = $string1\n";
			#}
		
			$string1 = sprintf('%.14e',$temp_val);
			$string1 =~ s/e/d/;
			$string2 = sprintf('%.14e',$temp_val2);
			$string2 =~ s/e/d/;
			$string3 = sprintf('%.14e',$temp_val3/$mass_conv); # g/mol -> kg
			$string3 =~ s/e/d/;
			$string4 = sprintf('%.14e',$pi);
			$string4 =~ s/e/d/;
			print OUT "\tDbg = $string2/(2.d0*rbg)*((5.d0+4.d0*($string1/rbg)+&\n\t\t\t\t& 6.d0*($string1/rbg)**2+18.d0*($string1/rbg)**3)/&\n\t\t\t\t& (5.d0-($string1/rbg)+(8.d0+$string4)*($string1/rbg)**2))\n";
			print OUT "\tvthbg = sqrt($string3/mbg)\n\n";
				
			$cs_coef_hydr[0][0] = "";
				
			# For comparisons with DMAN/TOMAS, coefficients for collisions involving an acid or ammonia monomer are calculated in a different way
			# -> uncomment the following block if you want to use this feature
			# if('A' ~~ @molecule_name[@mol_types_used[1..$nmol_types_used]] || 'N' ~~ @molecule_name[@mol_types_used[1..$nmol_types_used]]){
				# print "Applying the DMAN condensation coefficients to the loss of monomers A and N.\n\n";					
				# if($nmol_types_used>1){
					# die "Add this if you need it; don't have time to fix these stupidities now.\n";
				# }else{
					# $cs_coef_hydr[0][0] .= "\t\tif(i.eq.1) then\n";
					# # gas diffusivity in air from gasdiff.f, assuming standard pressure
					# $str_temp = sprintf('%.14e',1e-7*$temperature**1.75/101325*1e5/(42.88**(1/3)+20.1**(1/3))**2);
					# $str_temp =~ s/e/d/;
					# $cs_coef_hydr[0][0] .= "\t\t\tD = $str_temp*sqrt((m+28.9d0)/(m*28.9d0))\n";
					# $cs_coef_hydr[0][0] .= "\t\t\tKn = $string1/rbg\n";
					# if($molecule_name[$mol_types_used[1]] eq 'A'){
						# # the so-called accommodation coefficient is 0.65 for acid
						# $cs_coef_hydr[0][0] .= "\t\t\tcs(i) = 4.d0*$string4*D*rbg*(1.d0+Kn)/(1.d0+2.d0*Kn*(1.d0+Kn)/0.65)*cbg\n";
					# }else{
						# die "You have a one-component system, but the component is not acid (A)??\n";
					# }
					# $cs_coef_hydr[0][0] .= "\t\telse\n";
				# }
				# $str_temp = "\t";
				# $str_temp2 = "\t\tend if\n";
			# }else{
				# $str_temp = "";
				# $str_temp2 = "";
			# }
			$str_temp = "";
			$str_temp2 = "";
				
			$cs_coef_hydr[0][0] .= "$str_temp\t\tD = $string2/(2.d0*r)*((5.d0+4.d0*($string1/r)+&\n\t\t\t\t& 6.d0*($string1/r)**2+18.d0*($string1/r)**3)/&\n\t\t\t\t& (5.d0-($string1/r)+(8.d0+$string4)*($string1/r)**2))\n";
			$cs_coef_hydr[0][0] .= "$str_temp\t\tvth = sqrt($string3/m)\n";
			$cs_coef_hydr[0][0] .= "$str_temp\t\tKn = 2.d0*(D+Dbg)/(sqrt(vth**2+vthbg**2)*(r+rbg))\n";
            $cs_coef_hydr[0][0] .= "$str_temp\t\tcs(i) = 4.d0*$string4*(r+rbg)*(D+Dbg)*(1.d0+Kn)/(1.d0+2.d0*Kn*(1.d0+Kn))*cbg\n";
				
			$cs_coef_hydr[0][0] .= $str_temp2;
		}

################################################################################
################################################################################
##         9.12 Coagulation sink given by one number                          ##
################################################################################
################################################################################

	}elsif(&is_pos_number($coag_terms)){
		# Same coagulation frequency for all
		$coag_terms =~ s/E/e/;
		if($lfortran && !$variable_cs){
			$str_temp=sprintf('%.14e',$coag_terms);
			$str_temp =~ s/e/d/;
			print OUT "\tcs = $str_temp\n";
		}elsif(!$lfortran){
			print CS "CS = $coag_terms;\n";
		}
	}else{
		die "\nCan't understand coagulation sink definition $coag_terms.\n\n";
	}


################################################################################
################################################################################
##         9.13 Printing out size-dependent sink in dry case                  ##
################################################################################
################################################################################
	
	if(!$lhydrate_average && !$one_cs_for_all){
		if(!$lloop){
			for($iclus=1;$iclus<=$max_cluster_number;$iclus++){
				$cs_term_temp=sprintf('%.14e',$cs_coef_hydr[$iclus][0]);
				if($lfortran){
					$cs_term_temp =~ s/e/d/;
					print OUT "\tcs($iclus) = $cs_term_temp\t! coagulation loss of $label[$iclus]\n";
				}else{
					print CS "CS($iclus) = $cs_term_temp;\t% coagulation loss of $label[$iclus]\n";
				}
			}
		}else{
			print OUT "\tdo i = 1,nclust\n";
			if($nmol_types_used>1){
				print OUT "\t\tcall get_masses_and_radii$append_to_subroutine_names(m,r,indices(i,:))\n"
			}else{
				print OUT "\t\tcall get_masses_and_radii$append_to_subroutine_names(m,r,i)\n"
			}
			print OUT $cs_coef_hydr[0][0];
			print OUT "\tend do\n\n";
		}
	
################################################################################
################################################################################
##         9.14 Hydrate average for size-dependent sink                       ##
################################################################################
################################################################################	
	
	}elsif(!$one_cs_for_all){
		
		if($lfortran){
			print OUT "\t! The coagulation loss coefficients are averaged over all the available hydrates\n\t! Dry values are printed as comments for reference\n\n";
		}else{
			print CS "% The coagulation loss coefficients are averaged over all the available hydrates\n% Dry values are printed as comments for reference\n\n";
		}
		
		for($iclus=1;$iclus<=$max_cluster_number;$iclus++){
			$this_clus=$label[$iclus];
			$cs_term_temp=0.0;
			for($ihydr=0;$ihydr<=$#{ $distr[$iclus] };$ihydr++){
				if(defined $distr[$iclus][$ihydr]){
					$cs_term_temp += $distr[$iclus][$ihydr]*$cs_coef_hydr[$iclus][$ihydr];
				}
			}
			
			$cs_term_temp=sprintf('%.14e',$cs_term_temp);
			# dry value
			$str_temp=sprintf('%.14e',$cs_coef_hydr[$iclus][0]);

			if($lfortran){
				$cs_term_temp =~ s/e/d/;
				$str_temp =~ s/e/d/;
				if($lhas_hydr[$iclus]){
					$string="$cs_term_temp\t! coagulation loss of $label[$iclus], dry value $str_temp";
				}else{
					$string="$str_temp\t! coagulation loss of $label[$iclus], no hydrates";
				}
				print OUT "\tcs($iclus) = $string\n";
			}else{
				if($lhas_hydr[$iclus]){
					$string="CS($iclus) = $cs_term_temp;\t% coagulation loss of $label[$iclus], dry value $str_temp";
				}else{
					$string="CS($iclus) = $str_temp;\t% coagulation loss of $label[$iclus], no hydrates";
				}
				print CS "$string\n";
			}
		}
	}
		
	if(! $lfortran){
		print CS "\n";
		close(CS);
	}

}elsif(! $lfortran){
	print CS "% CS terms disabled\n\n";
	if($include_ion_terms){
			print CS "% enhancement factor for ions\nfcs = nan;\n\n";
	}
	print CS "CS = [];\n";
	close(CS);
}


################################################################################
################################################################################
##      9.2 Wall losses                                                       ##
################################################################################
################################################################################

if(! $lfortran){
	print WL "% Wall losses\n\n";
	if($include_ion_terms){
		print WL "% enhancement factor for ions\nfwl = $fwl;\n\n";
	}
}
if(!$disable_wall_terms || ($l_mcmc_any && $ncoef_wl>0)){
	if($lfortran){
		print OUT "\n\t! wall losses\n\n";
	}
	# array for the wall loss coefficients of the hydrates
	@wl_coef_hydr=();
	
################################################################################
################################################################################	
##         9.21 Reading in wall losses for all the clusters from a given file ##
################################################################################
################################################################################

	if($wall_terms eq 'wlfile'){
		print "\nwall loss rate file: $wl_coef_file_name\n";
		open(WL_IN, "<$wl_coef_file_name") || die "\nCould not open wall loss rate file $wl_coef_file_name\n\n";
		chomp(@temp_array=<WL_IN>);
		close(WL_IN);
		$number_of_lines=$#temp_array;
		for($iline=0;$iline<=$number_of_lines;$iline++){
			next if ($temp_array[$iline] =~ /^\#/);
			if($temp_array[$iline] =~ /(\d+$name_water-?\d+)|(\d+$name_water\s*$)/){
				die "$wl_coef_file_name: Wall loss rate cannot be given for hydrates for now (if hydrate averaging is used, the given rate for dry clusters is taken to be the effective rate)";
			}
			@columns=split(' ',$temp_array[$iline]);			
			die "\n$columns[0] on line $iline in file $wl_coef_file_name should be a number\n\n" if (!&is_pos_number($columns[0]));
			die "\nThere should be two columns on line $iline in file $wl_coef_file_name: the rate and the cluster\n\n" if ($#columns != 1);
		
			$iclus=&get_cluster_number($columns[1]);
			if($iclus ne ''){
				$wl_coef_hydr[$iclus][0]=$columns[0];
			}
		}
		for($iclus=1;$iclus<=$max_cluster_number;$iclus++){
			if(!defined $wl_coef_hydr[$iclus][0]){
				die "\nCould not find loss rate for $label[$iclus] from $wl_coef_file_name!\n\n";
			}
		}

		
################################################################################
################################################################################
##         9.22 Wall loss parameterizations                                   ##
################################################################################
################################################################################

	
################################################################################
################################################################################
##            9.221 Cloud 4 wall loss parametrization from Jasper's exel file ##
################################################################################
################################################################################

	}elsif($wall_terms =~ m/^cloud4_exel_file$/i){
	# CLOUD wall losses (parametrization from experiments)
	
		$air_viscosity = 1.708e-5*($temperature/273.15)**1.5*393.396/($temperature+120.246);
		
		if(!$lloop){
			for($iclus=1;$iclus<=$max_cluster_number;$iclus++){
			
				$this_clus=$label[$iclus];
						
				# loss coefficients for the dry cluster and possible hydrates
				for($ihydr=0;$ihydr<=$#{ $distr[$iclus] };$ihydr++){
					if(defined $distr[$iclus][$ihydr]){
						$mass1=$clust_mass[$iclus][$ihydr];
						$radius=$clust_radius[$iclus][$ihydr];
						# mobility diameter:
						$diameter = (2.0*$radius+0.3e-9)*sqrt(1+28.8*$mass_conv/$mass1);
						# slip correction factor:
						$wl_term_temp = 1+2.e-9/(101300*$diameter*.752e-6)*(6.32+2.01*exp(-0.1095e9*101300*$diameter*.752e-6));
						# diffusion coefficient (m^2/s):
						$wl_term_temp *= $boltz*$temperature/(3.*$pi*$air_viscosity*$diameter);
						# CLOUD4 wall loss (1/s):
						$wl_coef_hydr[$iclus][$ihydr] = .774*sqrt($wl_term_temp);	
					}
				}
			}
		}else{
			$wl_coef_hydr[0][0] = "\t\td = (2.d0*r+0.3d-9)*sqrt(1.d0+28.8d0/m)\n";
			$str_temp = sprintf('%.14e',2.e-9/(101300*.752e-6));
			$str_temp =~ s/e/d/;
			$str_temp2 = sprintf('%.14e',-0.1095e9*101300*.752e-6);
			$str_temp2 =~ s/e/d/;
			$wl_coef_hydr[0][0] .= "\t\tsc = 1.d0+$str_temp*(6.32d0+2.01d0*exp($str_temp2*d))/d\n";
			$str_temp = sprintf('%.14e',$boltz*$temperature/(3.*$pi*$air_viscosity));
			$str_temp =~ s/e/d/;
			$wl_coef_hydr[0][0] .= "\t\twl(i) = .774d0*sqrt(sc*$str_temp/d)\n";
		}

################################################################################
################################################################################
##            9.222 Cloud 4 wall loss parametrization from Andreas K\"urten   ##
################################################################################
################################################################################

	}elsif($wall_terms =~ m/^cloud4_Andreas$/i){
	# CLOUD wall losses (parametrization from experiments)
	
		$air_viscosity = (11.798+0.629976*$temperature-1.81158e-04*$temperature**2)*1e-07;
		
		if(!$lloop){
			for($iclus=1;$iclus<=$max_cluster_number;$iclus++){
			
				# loss coefficients for the dry cluster and possible hydrates
				for($ihydr=0;$ihydr<=$#{ $distr[$iclus] };$ihydr++){
					if(defined $distr[$iclus][$ihydr]){
						$mass1=$clust_mass[$iclus][$ihydr];
						$radius=$clust_radius[$iclus][$ihydr];
						# mobility diameter:
						$diameter = 2.0*$radius;   				##### Note: this is for some reason the geometric diameter.
						# lambda (at pressure 100000 Pa)
						$wl_term_temp = $boltz*$temperature/(sqrt(2)*100000*$pi*(0.37e-09)**2);
						# Knudsen number
						$wl_term_temp /= $diameter;
						# slip correction factor:
						$wl_term_temp = 1+$wl_term_temp*(1.142+0.558*exp(-0.999/$wl_term_temp));
						# diffusion coefficient (m^2/s):
						$wl_term_temp *= $boltz*$temperature/(3.*$pi*$air_viscosity*$diameter);
						# CLOUD4 wall loss (1/s):
						$wl_coef_hydr[$iclus][$ihydr] = .77*sqrt($wl_term_temp);	
					}
				}
			}
		}else{
			$wl_coef_hydr[0][0] = "\t\td = 2.d0*r\n";
			$str_temp = sprintf('%.14e',$boltz*$temperature/(sqrt(2)*100000*$pi*(0.37e-09)**2));
			$str_temp =~ s/e/d/;
			$wl_coef_hydr[0][0] .= "\t\tsc = 1.d0+$str_temp*(1.142d0+0.558d0*exp(-0.999d0/$str_temp*d))/d\n";
			$str_temp = sprintf('%.14e',$boltz*$temperature/(3.*$pi*$air_viscosity));
			$str_temp =~ s/e/d/;
			$wl_coef_hydr[0][0] .= "\t\twl(i) = .77d0*sqrt(sc*$str_temp/d)\n";
		}

################################################################################
################################################################################
##            9.223 Cloud 4 wall loss, simplified parametrization             ##
################################################################################
################################################################################
	
	}elsif($wall_terms =~ m/^cloud4$/i){
	# CLOUD wall losses (parametrization from experiments)
		if(!$lloop){
			for($iclus=1;$iclus<=$max_cluster_number;$iclus++){
				$this_clus=$label[$iclus];
						
				# loss coefficients for the dry cluster and possible hydrates
				for($ihydr=0;$ihydr<=$#{ $distr[$iclus] };$ihydr++){
					if(defined $distr[$iclus][$ihydr]){
						$mass1=$clust_mass[$iclus][$ihydr];
						$radius=$clust_radius[$iclus][$ihydr];
					
						# WL factor proportional to 1/mobility diameter
						$wl_term_temp = 1.0/((2.0*$radius+0.3e-9)*sqrt(1+28.8*$mass_conv/$mass1));
						$wl_term_temp *= 1.66e-12;
						$wl_coef_hydr[$iclus][$ihydr] = $wl_term_temp;
					}
				}
			}
		}else{
			$wl_coef_hydr[0][0] = "\t\twl(i) = 1.66d-12/((2.d0*r+0.3d-9)*sqrt(1.d0+28.8d0/m))\n";
		}


################################################################################
################################################################################
##            9.224 Cloud 3 wall loss parametrization                         ##
################################################################################
################################################################################
	
	}elsif($wall_terms =~ m/^cloud3$/i){
	# CLOUD wall losses (parametrization from experiments)
		if(!$lloop){
			for($iclus=1;$iclus<=$max_cluster_number;$iclus++){
				$this_clus=$label[$iclus];
						
				# loss coefficients for the dry cluster and possible hydrates
				for($ihydr=0;$ihydr<=$#{ $distr[$iclus] };$ihydr++){
					if(defined $distr[$iclus][$ihydr]){
						$mass1=$clust_mass[$iclus][$ihydr];
						$radius=$clust_radius[$iclus][$ihydr];

						# WL factor proportional to 1/mobility diameter
						$wl_term_temp = 1.0/((2.0*$radius+0.3e-9)*sqrt(1+28.8*$mass_conv/$mass1));
						$wl_term_temp *= 1.310e-12;
						$wl_coef_hydr[$iclus][$ihydr] = $wl_term_temp;
					}
				}
			}
		}else{
			$wl_coef_hydr[0][0] = "\t\twl(i) = 1.31d-12/((2.d0*r+0.3d-9)*sqrt(1.d0+28.8d0/m))\n";
		}
		
###############################################################################################################
###############################################################################################################
##            9.225 Diffusion wall loss for a flow tube (originally for the tube of Hanson and Eisele)       ##
###############################################################################################################
###############################################################################################################
	
	}elsif($wall_terms =~ m/^diffusion$/i){
		# Loss terms based on diffusion in a flow tube (originally for modeling the setup of Hanson and Eisele, 2000)
		# Diffusion coefficients are calculated as in the kinetic gas theory, and relative to the coefficient of sulfuric acid
	
		# Conversion from diffusion coefficient to wall loss factor
		if($tube_radius_in ne ''){
			$tube_radius=$tube_radius_in/100;
		}else{
			$tube_radius=0.049/2; # m (Hanson and Eisele, 2000)
		}		
		$D_to_wl=3.65/$tube_radius**2; # for laminar diffusion-limited flow (Brown, R.L.: Tubular flow reactor with first-order kinetics, J. Res. Natl. Bur. Stand. (U.S.), 83, 1, 1-6, 1978)
		
		# Properties of N2, the radius is calculated from viscosity (e.g. eq. 11-67 in Present: Kinetic theory of gases, 1958)
		$mass_N2=28.01*$mass_conv;
		
		
		# Sutherland formula, coefficients for N2:
		# C=111 from Crane (1988) (haven't seen the original paper, but lots of references citing the number)
		# mu0=17.9 at T0=300 from CRC
		
		$viscosity_N2=17.9*1e-6*(300+111)/($temperature+111)*($temperature/300)**1.5;
		$radius_N2=0.5*(5/16/$viscosity_N2)**0.5*($mass_N2*$boltz*$temperature/$pi)**0.25;
	
		# find the acid monomer
		$iclus=&get_cluster_number('1A');
		if($iclus eq ''){
			die "\nCan't calculate diffusion coefficients relative to that of sulfuric acid, if you don't have acid in the system\n\n";
		}else{
			$mass1=$clust_mass[$iclus][0];
			$radius=$clust_radius[$iclus][0];
		}

		# Loss coefficient of the acid monomer
		# Diffusion coefficient of acid in N2 (e.g. eq. 8-87 in Present: Kinetic theory of gases, 1958)
		# (experimental coefficient from Hanson and Eisele is 0.094 +- 0.0012 atm cm^2 s^-1)
		$D_acid_N2=3/8/101325*$boltz*$temperature/($radius+$radius_N2)**2*sqrt($boltz*$temperature/2/$pi*(1/$mass1+1/$mass_N2)); # at 1 atm
		if($tube_pressure_in ne ''){
			$wl_term_0=$D_acid_N2/($tube_pressure_in/101325)*$D_to_wl; # s^-1, converted to the given pressure
		}else{
			$wl_term_0=$D_acid_N2/(133.322/101325*620)*$D_to_wl; # s^-1, converted to 620 Torr (Hanson and Eisele, 2000)
		}

		$D0_factor=1/($radius+$radius_N2)**2*sqrt(1/$mass1+1/$mass_N2);
	
		if(!$lloop){
			for($iclus=1;$iclus<=$max_cluster;$iclus++){
				$this_clus=$label[$iclus];
						
				# loss coefficients for the dry cluster and possible hydrates
				for($ihydr=0;$ihydr<=$#{ $distr[$iclus] };$ihydr++){
					if(defined $distr[$iclus][$ihydr]){
						$mass1=$clust_mass[$iclus][$ihydr];
						$radius=$clust_radius[$iclus][$ihydr];
					
						$D_factor=1/($radius+$radius_N2)**2*sqrt(1/$mass1+1/$mass_N2);
						$wl_coef_hydr[$iclus][$ihydr]=$D_factor/$D0_factor*$wl_term_0;
					}
				}
			}
		}else{
			$radius_N2 = sprintf('%.14e',$radius_N2);
			$radius_N2 =~ s/e/d/;
			$mass_N2 = sprintf('%.14e',$mass_N2);
			$mass_N2 =~ s/e/d/;
			$D0_factor = sprintf('%.14e',$D0_factor);
			$D0_factor =~ s/e/d/;
			$wl_term_0 = sprintf('%.14e',$wl_term_0);
			$wl_term_0 =~ s/e/d/;
			$wl_coef_hydr[0][0] = "\t\twl(i) = 1.d0/(r+$radius_N2)**2*sqrt(1.d0/m+1.d0/$mass_N2)/$D0_factor*$wl_term_0\n";
		}
		
		
################################################################################
################################################################################
##         9.23 Wall losses given by one number                               ##
################################################################################
################################################################################


################################################################################
################################################################################
##            9.231 Ift flow tube wall loss                                   ##
################################################################################
################################################################################
	
	}elsif($wall_terms =~ m/ift/i){
	# Ift-lft wall losses (estimated, Berndt & Richters, 2011, http://dx.doi.org/10.1016/j.atmosenv.2011.10.060)
		
		$wl_term_temp = 2.3e-2;
		
		if($#wl_only>=0){
			if($lfortran){
				$wl_term_temp = sprintf('%.14e',$wl_term_temp);
			}
		}else{
			if($lfortran){
				$str_temp=sprintf('%.14e',$wl_term_temp);
				$str_temp =~ s/e/d/;
				print OUT "\twl = $str_temp\n";
			}else{
				print WL "WL = $wl_term_temp;\n";
			}
		}
	
################################################################################
################################################################################
##            9.232 Other size-independent wall loss                          ##
################################################################################
################################################################################
	
	}elsif(&is_pos_number($wall_terms)){
	# Other wall losses (same number for all)
		$wall_terms =~ s/E/e/;
		if($#wl_only>=0){
			if($lfortran){
				$wl_term_temp = sprintf('%.14e',$wall_terms);
			}else{
				$wl_term_temp = $wall_terms;
			}
		}else{
			if($lfortran){
				$str_temp=sprintf('%.14e',$wall_terms);
				$str_temp =~ s/e/d/;
				print OUT "\twl = $str_temp\n";
			}else{
				print WL "WL = $wall_terms;\n";
			}
		}

	}else{	
		die "\nCan't understand wall loss definition $wall_terms.\n\n";
	}
	
	
############################################################################################
############################################################################################
##            9.233 Wall losses for specific clusters from the command line               ##
############################################################################################
############################################################################################
	
	# If we have special wall losses for some clusters,
	if($#wl_only>=0){
		# go through the input and separate the cluster and the rate
		for($iclus=1;$iclus<=$max_cluster_number;$iclus++){
			$wl_only_hydr[$iclus] = ();
		}
		for($icol=0;$icol<=$#wl_only;$icol++){
			@columns=split(',',$wl_only[$icol]);
			if ($#columns!=1){
				die "Can't understand the special wall loss $wl_only[$icol] - the format should be 'clust,number'.\n";
			}
			($iclus,$nwaters)=&get_corresponding_dry($columns[0]);
			if($iclus eq ''){
				die "Cluster $columns[0] from the special wall loss $wl_only[$icol] does not exist.\n";
			}
			$wl_only_hydr[$iclus][$nwaters] = sprintf('%.14e',$columns[1]);
		}
		
		# Then calculate the hydrate averages if needed or just put the values back in the same array
		@wl_only = ();
		$icol = -1;
		for($iclus=1;$iclus<=$max_cluster_number;$iclus++){
			next if ($#{ $wl_only_hydr[$iclus] }==-1);
			$icol += 1;
			$wl_only[$icol][0] = $iclus;
			$wl_only[$icol][1] = 0;
			for($ihydr=0;$ihydr<=$#{ $distr[$iclus] };$ihydr++){
				if(defined $distr[$iclus][$ihydr]){
					if(defined $wl_only_hydr[$iclus][$ihydr]){
						$wl_only[$icol][1] += $distr[$iclus][$ihydr]*$wl_only_hydr[$iclus][$ihydr];
					}else{
						if($#{ $wl_only_hydr[$iclus] }==0){
							print "You have defined special wall losses for $label[$iclus] but not its hydrates -> using the given value for the hydrate average.\n";
							$wl_only[$icol][1] = $wl_only_hydr[$iclus][0];
							last;
						}else{
							print "You have defined special wall losses for some but not all hydrates of $label[$iclus] -> using the standard wall losses for the others.\n";
							if($one_wl_for_almost_all){
								$wl_only[$icol][1] += $distr[$iclus][$ihydr]*$wl_term_temp;
							}else{
								$wl_only[$icol][1] += $distr[$iclus][$ihydr]*$wl_coef_hydr[$iclus][$ihydr];
							}
						}
					}
				}
			}
		}
	}
	if($one_wl_for_almost_all){
		# If we have a constant wall loss for most clusters and a special value for others,
		# print the coefficients here.
		if($#wl_only>=0){
			for($iclus=1;$iclus<=$max_cluster;$iclus++){
				$str_temp = $wl_term_temp;
				for($icol=0;$icol<=$#wl_only;$icol++){
					if($wl_only[$icol][0]==$iclus){
						$str_temp = $wl_only[$icol][1];
						if(!&is_pos_number($str_temp)){
							die "\nCan't understand wall loss definition $str_temp for cluster $label[$iclus]\n";
						}else{
							print "Using special wall loss $str_temp for cluster $label[$iclus].\n";
						}
						last;
					}
				}
				if($lfortran){
					$str_temp =~ s/e/d/;
					print OUT "\twl($iclus) = $str_temp\t! wall loss of $label[$iclus]\n";
				}else{
					print WL "WL($iclus) = $str_temp;\t% wall loss of $label[$iclus]\n";
				}
			}
		}
		
		
		# If we have MCMC wall losses for some clusters and a size-independent value for others,
		# print the MCMC coefficients here.
		if($l_mcmc_any){
			if($l_wl_mcmc){
				print OUT "\twl = coef($n_mcmc_coefs)\t\t\t\t! wall loss\n";
				$coef_true[$n_mcmc_coefs] = '';
			}
			$icoef_wl = $n_mcmc_coefs-$ncoef_wl;
			for($iclus=1;$iclus<=$max_cluster_number;$iclus++){
				if($l_special_mcmc_wl[$iclus]){
					$icoef_wl++;
					print OUT "\twl($iclus) = coef($icoef_wl)\t\t\t! wall loss of $label[$iclus]\n";
					$coef_true[$icoef_wl] = '';
				}
			}
		}
	
	
################################################################################
################################################################################
##         9.24 Printing out size-dependent wall losses in dry case           ##
################################################################################
################################################################################
	
	}elsif((!$lhydrate_average || $wall_terms eq 'wlfile') && !$one_wl_for_all){
		if($wall_terms eq 'wlfile'){
			$string=", given as input";
		}else{
			$string="";
		}
		if($l_mcmc_any){
			$icoef_wl = $n_mcmc_coefs-$ncoef_wl;
		}
		if(!$lloop){
			for($iclus=1;$iclus<=$max_cluster_number;$iclus++){
				$wl_term_temp=sprintf('%.14e',$wl_coef_hydr[$iclus][0]);
				# If there are some special cases, going through those now
				if($#wl_only>=0){
					for($icol=0;$icol<=$#wl_only;$icol++){
						if($wl_only[$icol][0]==$iclus){
							$wl_term_temp = $wl_only[$icol][1];
							if(!&is_pos_number($wl_term_temp)){
								die "\nCan't understand wall loss definition $wl_term_temp for cluster $label[$iclus]\n";
							}else{
								print "Using special wall loss $wl_term_temp for cluster $label[$iclus].\n";
							}
							last;
						}
					}
				}

				if($lfortran){
					$wl_term_temp =~ s/e/d/;
					if($l_special_mcmc_wl[$iclus]){
						# If this is an MCMC coefficient, save the true value to be printed later ...
						$icoef_wl++;
						$coef_true[$icoef_wl] = "$wl_term_temp\t! wall loss of $label[$iclus]$string";
						# ... and print here just the MCMC coefficient
						print OUT "\twl($iclus) = coef($icoef_wl)\t\t\t\t! wall loss of $label[$iclus]\n";
					}else{
						print OUT "\twl($iclus) = $wl_term_temp\t! wall loss of $label[$iclus]$string\n";
					}
				}else{
					print WL "WL($iclus) = $wl_term_temp;\t% wall loss of $label[$iclus]$string\n";
				}
			}
		}else{
			print OUT "\tdo i = 1,nclust\n";
			if($nmol_types_used>1){
				print OUT "\t\tcall get_masses_and_radii$append_to_subroutine_names(m,r,indices(i,:))\n"
			}else{
				print OUT "\t\tcall get_masses_and_radii$append_to_subroutine_names(m,r,i)\n"
			}
			print OUT $wl_coef_hydr[0][0];
			print OUT "\tend do\n\n";
		}
	
################################################################################
################################################################################
##         9.25 Hydrate average for size-dependent wall losses                ##
################################################################################
################################################################################	
	
	}elsif(!$one_wl_for_all){
	
		if($l_mcmc_any){
			$icoef_wl = $n_mcmc_coefs-$ncoef_wl;
		}
		
		if($lfortran){
			print OUT "\t! The wall loss coefficients are averaged over all the available hydrates\n\t! Dry values are printed as comments for reference\n\n";
		}else{
			print WL "% The wall loss coefficients are averaged over all the available hydrates\n% Dry values are printed as comments for reference\n\n";
		}
		
		for($iclus=1;$iclus<=$max_cluster_number;$iclus++){
			$this_clus=$label[$iclus];
			$wl_term_temp=0.0;
			for($ihydr=0;$ihydr<=$#{ $distr[$iclus] };$ihydr++){
				if(defined $distr[$iclus][$ihydr]){
					$wl_term_temp += $distr[$iclus][$ihydr]*$wl_coef_hydr[$iclus][$ihydr];
				}
			}
			
			$wl_term_temp=sprintf('%.14e',$wl_term_temp);
			# dry value
			$str_temp=sprintf('%.14e',$wl_coef_hydr[$iclus][0]);
			
			# If there are some special cases, going through those now
			if($#wl_only>=0){
				for($icol=0;$icol<=$#wl_only;$icol++){
					if($wl_only[$icol][0]==$iclus){
						$wl_term_temp = $wl_only[$icol][1];
						$str_temp = "$str_temp (using a special wall loss)";
						if(!&is_pos_number($wl_term_temp)){
							die "\nCan't understand wall loss definition $wl_term_temp for cluster $label[$iclus]\n";
						}else{
							print "Using special wall loss $wl_term_temp for cluster $label[$iclus].\n";
						}
						last;
					}
				}
			}

			if($lfortran){
				$wl_term_temp =~ s/e/d/;
				$str_temp =~ s/e/d/;
				if($lhas_hydr[$iclus]){
					$string="$wl_term_temp\t! wall loss of $label[$iclus], dry value $str_temp";
				}else{
					$string="$str_temp\t! wall loss of $label[$iclus], no hydrates";
				}
				if($l_special_mcmc_wl[$iclus]){
					# If this is an MCMC coefficient, save the true value to be printed later ...
					$icoef_wl++;
					$coef_true[$icoef_wl] = "$string";
					# ... and print here just the MCMC coefficient
					print OUT "\twl($iclus) = coef($icoef_wl)\t\t\t\t! wall loss of $label[$iclus]\n";
				}else{
					print OUT "\twl($iclus) = $string\n";
				}
			}else{
				if($lhas_hydr[$iclus]){
					$string="WL($iclus) = $wl_term_temp;\t% wall loss of $label[$iclus], dry value $str_temp";
				}else{
					$string="WL($iclus) = $str_temp;\t% wall loss of $label[$iclus], no hydrates";
				}
				print WL "$string\n";
			}
		}
	}
		
	if(! $lfortran){
		print WL "\n";
		close(WL);
	}
}elsif(! $lfortran){
	print WL "% WL terms disabled\n\n";
	if($include_ion_terms){
		print WL "% enhancement factor for ions\nfwl = nan;\n\n";
	}
	print WL "WL = [];\n";
	close(WL);
}

if($lfortran && $l_use_get_losses){
	if($lloop){
		$str_temp = $lines{'inp'};
		$str_temp =~ s/,/+/;
		print OUT "\tloss = $str_temp\n";
	}
	print OUT "\nend subroutine get_losses$append_to_subroutine_names\n\n!-----------------------------------------------------------\n\n";
}

################################################################################
################################################################################
##   10. Collision rates                                                      ##
################################################################################
################################################################################

if($lfortran){
	if($l_coll_factor_mcmc_any){
		$icoef_coll=0;
	}
	if($l_evap_factor_mcmc){
		$icoef_evap=0;
	}
	print OUT "subroutine get_coll$append_to_subroutine_names(K";
	if($l_mcmc_any){
		print OUT ",coef";
	}elsif($variable_temp){
		print OUT ",temperature";
	}
	print OUT ")\n\timplicit none\n";
	print OUT "\tinteger, parameter :: nclust = $max_cluster_number\n";
	print OUT "\treal(kind(1.d0)) :: K(nclust,nclust)";
	if($l_mcmc_any){
		print OUT ", coef($n_mcmc_coefs)";
	}elsif($variable_temp){
		print OUT ", temperature";
	}
	print OUT "\n";
	if(!$lloop){
		print OUT "\n\t! collision coefficients\n\n";
		print OUT "\tK = 0.d0\n";
	}
}else{
	open(COLL, ">get_coll$append_to_file_names.m");
	if($variable_temp){
		print COLL "function K = get_coll(temperature)\n\n";
	}
	print COLL "% Collision coefficients\n\n";
	print COLL "K = nan($max_cluster_number,$max_cluster_number);\n\n";
}

################################################################################
################################################################################
##      10.1 Reading in rates for specific collisions                         ##
################################################################################
################################################################################

if($collision_coef_file_name ne ''){
	print "\ncollision rate file: $collision_coef_file_name\n";
	open(COLLISION, "<$collision_coef_file_name") || die "\nCould not open collision rate file $collision_coef_file_name\n\n";
	chomp(@temp_array=<COLLISION>);
	close(COLLISION);
	$number_of_lines=$#temp_array;
	$iline2 = 0;
	for($iline=0;$iline<=$number_of_lines;$iline++){
		next if ($temp_array[$iline] =~ /^\#/);
		if($temp_array[$iline] =~ /(\d+$name_water-?\d+)|(\d+$name_water\s*$)/){
			die "$collision_coef_file_name: Collision rate cannot be given for hydrates for now (if hydrate averaging is used, the given rate for dry clusters is taken to be the effective rate)";
		}
		@columns=split(' ',$temp_array[$iline]);			
		die "\n$columns[0] on line $iline in file $collision_coef_file_name should be a number\n\n" if (!&is_pos_number($columns[0]));
		die "\nThere should be three columns on line $iline in file $collision_coef_file_name: the rate and the two collision parties\n\n" if ($#columns != 2);
		$collision_coef_in[$iline2][0] = $columns[0];
		push(@{ $collision_coef_in[$iline2] },$columns[1]);
		push(@{ $collision_coef_in[$iline2] },$columns[2]);
		$iline2++;
	}
}

################################################################################
################################################################################
##      10.2 Reading in sticking factors                                      ##
################################################################################
################################################################################

if(($sticking_factor[0][0] == 1) && ($sticking_factor[1][0] == 1) && ($sticking_factor_file_name eq '')){
	$l_sticking_factor = 0;
}else{
	$l_sticking_factor = 1;
	if($sticking_factor[1][0] != 1){
		if($sticking_factor[0][0] == 1){
			$temp_val = $sticking_factor[1][0];
			@sticking_factor = ();
			$sticking_factor[0][0] = $temp_val;
			$sticking_factor[0][1] = 'ion-neutral';
		}else{
			$sticking_factor[1][0] = $sticking_factor[1][0];
			$sticking_factor[1][1] = 'ion-neutral';
		}
	}elsif($sticking_factor[0][0]){
		$temp_val = $sticking_factor[0][0];
		@sticking_factor = ();
		$sticking_factor[0][0] = $temp_val;
	}else{
		@sticking_factor = ();
	}
	if($lfortran && !$lloop){
		print OUT "\t! Using sticking factors that are visible in front of the collision coefficients\n";
	}elsif(!$lfortran){
		print COLL "% Using sticking factors that are visible in front of the collision coefficients\n";
	}
	if($sticking_factor_file_name ne ''){
		print "sticking factor file: $sticking_factor_file_name\n";
		open(STICKING, "<$sticking_factor_file_name") || die "\nCould not open sticking factor file $sticking_factor_file_name\n\n";
		chomp(@temp_array=<STICKING>);
		close(STICKING);
		$number_of_lines=$#temp_array;
		$iline2 = $#sticking_factor+1;
		for($iline=0;$iline<=$number_of_lines;$iline++){
			next if ($temp_array[$iline] =~ /^\#/);
			@columns=split(' ',$temp_array[$iline]);			
			die "\n$columns[0] on line $iline in file $sticking_factor_file_name should be a number\n\n" if (!&is_pos_number($columns[0]));
			$sticking_factor[$iline2][0] = $columns[0];
			if($#columns>0){
				push(@{ $sticking_factor[$iline2] },$columns[1]);
				if($#columns>1){
					push(@{ $sticking_factor[$iline2] },$columns[2]);
				}
			}
			$iline2++;
		}
	}
}

################################################################################
################################################################################
##      10.3 Reading in dipole moments and polarizabilities                   ##
##           for collision enhancement factors                                ##
################################################################################
################################################################################

if($include_ion_terms){
	# The Su73 and Su82 ion collision parametrizations require dipole moments and polarizabilities
	# but if all ion collisions have mcmc collision rates, the data is not obligatory.
	if($ion_coll_method =~ /Su73/i || $ion_coll_method =~ /Su82/i && !(!$dip_file_name && $l_coll_factor_mcmc_n_i && ($l_coll_factor_mcmc_n_g || (!$include_generic_neg && !$include_generic_pos)))){
		if(!$dip_file_name){
			die "\nYou didn't give dipole data file name!\n\n";
		}
		
		# Checking whether the energy and dipole files seem to belong together if both are used
		
		

		# read in the dipole moments (D) and polarizabilities (Å^3)
		open(DIP, "<$dip_file_name") || die "\nCould not open dipole file $dip_file_name\n\n";
		chomp(@dip_data=<DIP>);
		close(DIP);
		$number_of_lines=$#dip_data+1;
		$iline0 = 0;
		# first two lines have the dipole locking coefficient for monomers and clusters in Su 73 method
		if(&is_pos_number($dip_data[0])){
			$c_mon=$dip_data[0];
			if ($c_mon < 0 || $c_mon > 1){
				die "\nDipole locking coefficients in $dip_file_name have to be in the range [0,1].\n\n";
			}
			$iline0++;
		}elsif($ion_coll_method =~ /Su73/i){
			die "\nFirst line of $dip_file_name has to be the dipole locking coefficient for monomers.\n\n";
		}
		if(&is_pos_number($dip_data[1])){
			$c_clust=$dip_data[1];
			if ($c_clust < 0 || $c_clust > 1){
				die "\nDipole locking coefficients in $dip_file_name have to be in the range [0,1].\n\n";

			}
			$iline0++;
		}elsif($ion_coll_method =~ /Su73/i){
			die "\nSecond line of $dip_file_name has to be the dipole locking coefficient for clusters\n\n";
		}
		for($iline=$iline0;$iline<$number_of_lines;$iline++){
			next if ($dip_data[$iline] =~ /^\#/);

			if($dip_data[$iline] =~ /^\s*(\S+)/){

				$temp_label=$1;
				
				# Checking if we have the cluster (or the corresponding dry one 
				# if we are using hydrates) in the system
				($iclus,$nwaters)=&get_corresponding_dry($temp_label);
				
				# save the dipole data if we have the cluster in the system
				if($iclus ne ''){
					# Each cluster must appear only once in the file, so the values should not be defined yet
					if(exists $dip_hydr[$iclus][$nwaters]){
						die "Diple data of cluster $temp_label is multiply defined in the file $dip_file_name.\n";
					}
					# there should be two values
					if($dip_data[$iline] =~ /^\s*(\S+)\s+(\S+)\s+(\S+)\s*/){
						$temp_val = $2;
						$temp_val2 = $3;
						if($ion_coll_method =~ /Su73/i){
							if($lmonomer[$iclus]){
								$temp_val *= $c_mon;
							}else{
								$temp_val *= $c_clust;
							}
						}
						if($nwaters==0){
							$dip[$iclus] = $temp_val;
							$pol[$iclus] = $temp_val2;
						}
						$dip_hydr[$iclus][$nwaters] = $temp_val;
						$pol_hydr[$iclus][$nwaters] = $temp_val2;
					}else{
						die "\nCouldn't read dipole moment and polarizability for $label in $dip_file_name.\n\n";
					}
				}
			}
		}

		# now all the clusters should have a dipole moment and polarizability values
		for($iclus=1;$iclus <= $max_cluster;$iclus++){
			next if($iclus==$n_flux_out || $iclus==$n_flux_out_neg || $iclus==$n_flux_out_pos);
			if($lneutral_cluster[$iclus]){
				for($ihydr=0;$ihydr<=$#{ $distr[$iclus] };$ihydr++){
					if(!(exists $dip_hydr[$iclus][$ihydr])){
						if($ihydr==0){
							die "\nDipole data not found in $dip_file_name for $label[$iclus].\n\n";
						}else{
							die "\nDipole data not found in $dip_file_name for $label[$iclus]$ihydr$name_water.\n\n";
						}
					}
				}
			}
		}
	}
}

################################################################################
################################################################################
##      10.4 Calculating the collision rates                                  ##
################################################################################
################################################################################

$missing_fcr = 0;

# This is just to go through the collision rate loop once
if($lloop){
	$n_flux_out_max = 1;
	$ind_quad_loss[1][0] = 1;
	$ind_quad_loss[1][1] = 1;
	$ind_quad_loss[1][2] = 1;
}

# Looping through the possible collisions
for($iclus=1;$iclus<=$n_flux_out_max;$iclus++){
	$lsame = 0;
	for($i=1;$i<=2*$ind_quad_loss[$iclus][0];$i+=2){
		$jclus = $ind_quad_loss[$iclus][$i];
		$kclus = $ind_quad_loss[$iclus][$i+1];
		
		# In the fortran version collisions of identical clusters are listed twice, 
		# but we only want to print the collision rates once here
		if($iclus==$jclus){
			if($lsame){
				$lsame = 0;
				next;
			}else{
				$lsame = 1;
			}
		}
		# Also collisions of different clusters need to be printed only once
		if(!$lloop && $jclus!=$n_flux_out && $jclus!=$n_flux_out_neg && $jclus!=$n_flux_out_pos && $iclus<$jclus){
			next;
		}
		
		$iclus_label=$label[$iclus];
		$jclus_label=$label[$jclus];
		
################################################################################
################################################################################
##         10.41 Determining if this collision has been given                 ##
##               a fixed rate value in the collision rate file                ##
################################################################################
################################################################################

		$lcoll_coef_in[$iclus][$jclus] = 0;
		
		if($collision_coef_file_name ne ''){
			for($iline=0;$iline<=$#collision_coef_in;$iline++){
				if((&compare_clusters($iclus_label,$collision_coef_in[$iline][1]) && &compare_clusters($jclus_label,$collision_coef_in[$iline][2])) || (&compare_clusters($jclus_label,$collision_coef_in[$iline][1]) && &compare_clusters($iclus_label,$collision_coef_in[$iline][2]))){
					$lcoll_coef_in[$iclus][$jclus] = 1;
					$coll_coef[$iclus][$jclus] = $collision_coef_in[$iline][0];
					last;
				}
			}
		}
		
		if($lcoll_coef_in[$iclus][$jclus] && !$lloop){
			
			$coll_coef_eff[$iclus][$jclus] = $coll_coef[$iclus][$jclus];
			$coll_coef_hydr[$iclus][$jclus][0][0] = $coll_coef[$iclus][$jclus]; # this is for calculating the evaporation rate, if it's not given as input
			$sticking_coef[$iclus][$jclus]="";
			$sticking_coef[$jclus][$iclus]="";
			
			if($lneutral_cluster[$iclus]+$lneutral_cluster[$jclus]==1){
				$fcr[$iclus][$jclus] = 1;
				$fcr_hydr[$iclus][$jclus][0][0] = 1; # this is for calculating the evaporation rate, if it's not given as input
			}
		
################################################################################
################################################################################
##         10.42 Printing a function for calculating                          ##
##               the collision rate when using loops                          ##
################################################################################
################################################################################

			
		}elsif($lloop){
			print OUT "\treal(kind(1.d0)) :: ri, rj, mi, mj\n";
			if($loop_coll_coef =~ m/^Dahneke$/i){
				print OUT "\treal(kind(1.d0)) :: Di, Dj, vthi, vthj, Kn\n";
			}
			print OUT "\tinteger :: i, j\n\n";
			if(!$lloop_diag && $nmol_types_used>1){
				print OUT "\tinteger :: indices(nclust,$nmol_types_used)\n\n";
				print OUT "\tcall get_molecule_numbers$append_to_subroutine_names(indices)\n\n";
			}
			print OUT "\tK = 0.d0\n\n";
			print OUT "\tdo i = 1,nclust\n";
			if($nmol_types_used>1){
				print OUT "\t\tcall get_masses_and_radii$append_to_subroutine_names(mi,ri,indices(i,:))\n";
				print OUT "\t\tdo j = i,nclust\n";
				print OUT "\t\t\tcall get_masses_and_radii$append_to_subroutine_names(mj,rj,indices(j,:))\n";
			}else{
				print OUT "\t\tcall get_masses_and_radii$append_to_subroutine_names(mi,ri,i)\n";
				print OUT "\t\tdo j = i,nclust\n";
				print OUT "\t\t\tcall get_masses_and_radii$append_to_subroutine_names(mj,rj,j)\n";
			}			
			
			if($loop_coll_coef =~ m/^hard_spheres$/i){
				# Hard-sphere collision rate
				$temp_val = 8*$pi*$boltz*$temperature/$mass_conv;
				$string1 = sprintf('%.14e',$temp_val);
				$string1 =~ s/e/d/;
				print OUT "\t\t\tK(i,j) = sqrt($string1*(1.d0/mi+1.d0/mj))*(ri+rj)**2\n";
			}elsif($loop_coll_coef =~ m/^Dahneke$/i){
				# Collision rate parameterization by Dahneke (1983) with a non-continuum correction
				# Seinfeld and Pandis, 2nd edition, chapters 9, 12 and 13
				$air_viscosity = 2.5277e-7*$temperature**0.75302; # from DMAN
				# mean free path of air from S&P Eq. 9.6, assuming standard pressure
				$string1 = sprintf('%.14e',2*$air_viscosity/(101325*sqrt(8*0.0289/($pi*$Na*$boltz*$temperature))));
				$string1 =~ s/e/d/;
				# diffusivities
				# S&P Eq. 9.73, but don't know where the factor (Cc?) with the Kn term comes from (it is like this in DMAN)
				$string2 = sprintf('%.14e',$boltz*$temperature/(3*$pi*$air_viscosity));
				$string2 =~ s/e/d/;
				# factor for the thermal velocities
				$string3 = sprintf('%.14e',8*$boltz*$temperature/$pi/$mass_conv);
				$string3 =~ s/e/d/;
				$string4 = sprintf('%.14e',$pi);
				$string4 =~ s/e/d/;
				
				# For comparisons with DMAN/TOMAS, coefficients for collisions involving an acid or ammonia monomer are calculated in a different way
				# -> uncomment the following block if you want to use this feature
				# Note that when two monomers collide, one of them is treated as a molecule and one as a particle (for mathematical consistency)
				# if('A' ~~ @molecule_name[@mol_types_used[1..$nmol_types_used]] || 'N' ~~ @molecule_name[@mol_types_used[1..$nmol_types_used]]){
					# print "Applying the DMAN condensation coefficients to collisions involving monomers A and N.\n\n";					
					# if($nmol_types_used>1){
					    # die "Add this if you need it; don't have time to fix these stupidities now.\n";
						# #print OUT "\t\t\tif(sum(indices(i,:)).eq.1 .and. sum(indices(j,:)).eq.1) then\n";
						# #print OUT "\t\t\tend if\n";
					# }else{
						# print OUT "\t\t\tif(i.eq.1) then\n";
						# # gas diffusivity in air from gasdiff.f, assuming standard pressure
						# $str_temp = sprintf('%.14e',1e-7*$temperature**1.75/101325*1e5/(42.88**(1/3)+20.1**(1/3))**2);
						# $str_temp =~ s/e/d/;
						# print OUT "\t\t\t\tDi = $str_temp*sqrt((mi+28.9d0)/(mi*28.9d0))\n";
						# print OUT "\t\t\t\tKn = $string1/rj\n";
						# if($molecule_name[$mol_types_used[1]] eq 'A'){
							# # the so-called accommodation coefficient is 0.65 for acid
							# print OUT "\t\t\t\tK(i,j) = 4.d0*$string4*Di*rj*(1.d0+Kn)/(1.d0+2.d0*Kn*(1.d0+Kn)/0.65)\n";
						# }else{
							# die "You have a one-component system, but the component is not acid (A)??\n";
						# }
						# print OUT "\t\t\telse\n";
					# }
					# $str_temp = "\t";
					# $str_temp2 = "\t\t\tend if\n";
				# }else{
					# $str_temp = "";
					# $str_temp2 = "";
				# }
				$str_temp = "";
				$str_temp2 = "";
				
				print OUT "$str_temp\t\t\tDi = $string2/(2.d0*ri)*((5.d0+4.d0*($string1/ri)+&\n\t\t\t\t& 6.d0*($string1/ri)**2+18.d0*($string1/ri)**3)/&\n\t\t\t\t& (5.d0-($string1/ri)+(8.d0+$string4)*($string1/ri)**2))\n";
				print OUT "$str_temp\t\t\tDj = $string2/(2.d0*rj)*((5.d0+4.d0*($string1/rj)+&\n\t\t\t\t& 6.d0*($string1/rj)**2+18.d0*($string1/rj)**3)/&\n\t\t\t\t& (5.d0-($string1/rj)+(8.d0+$string4)*($string1/rj)**2))\n";
				print OUT "$str_temp\t\t\tvthi = sqrt($string3/mi)\n";
				print OUT "$str_temp\t\t\tvthj = sqrt($string3/mj)\n";
			    # S&P eqn 12.35 with the relative velocity
				print OUT "$str_temp\t\t\tKn = 2.d0*(Di+Dj)/(sqrt(vthi**2+vthj**2)*(ri+rj))\n";
                # S&P Eq. 13.50 with a non-continuum correction beta as in S&P Table 12.1
                print OUT "$str_temp\t\t\tK(i,j) = 4.d0*$string4*(ri+rj)*(Di+Dj)*(1.d0+Kn)/(1.d0+2.d0*Kn*(1.d0+Kn))\n";
				
				print OUT $str_temp2;
				
			}else{
				die "\nCan't understand collision coefficient definition $loop_coll_coef.\n\n";
			}
				
			print OUT "\t\t\tK(j,i) = K(i,j)\n";
			print OUT "\t\tend do\n";
			print OUT "\tend do\n\n";
			if($sticking_factor[0][0]!=1){
				print OUT "\t! sticking/enhancement factor\n";
				print OUT "\tK = $sticking_factor[0][0]*K\n\n";
			}
		
################################################################################
################################################################################
##         10.43 Calculating the hard sphere collision rate                   ##
################################################################################
################################################################################

		}else{
		
			@distr_temp = @{$distr[$iclus]};
			@distr_temp2 = @{$distr[$jclus]};
		
			# This will be the effective collision rate summed over all hydrates
			$coll_coef_eff[$iclus][$jclus]=0.0;
#print "\n";
		
			for($ihydr=0;$ihydr<=$#distr_temp;$ihydr++){
				if($ihydr != 0){
					$temp_label=$iclus_label.$ihydr.$name_water;
				}else{
					$temp_label=$iclus_label;
				}
				for($jhydr=0;$jhydr<=$#distr_temp2;$jhydr++){
					if($jhydr != 0){
						$temp_label2=$jclus_label.$jhydr.$name_water;
					}else{
						$temp_label2=$jclus_label;
					}
				
					if(defined $distr_temp[$ihydr] && defined $distr_temp2[$jhydr]){
					
						# If the clusters are ions with opposite charge, we'll use a constant recombination coefficient
						if(($lnegative[$iclus]&&$lpositive[$jclus]) || ($lpositive[$iclus]&&$lnegative[$jclus])){
							$coll_coef_hydr[$iclus][$jclus][$ihydr][$jhydr] = $recomb_coeff;
							$coll_coef_eff[$iclus][$jclus] = $recomb_coeff;
							next;
					
						# Otherwise calculate the collision coefficient from kinetic gas theory
						# Ion-neutral collision enhancement factors are done later
						}else{
							$mass1=$clust_mass[$iclus][$ihydr];
							$mass2=$clust_mass[$jclus][$jhydr];

							$r1=$clust_radius[$iclus][$ihydr];
							$r2=$clust_radius[$jclus][$jhydr];
#print "mass($label[$iclus]$ihydr" . "W)=$mass1\n";
#print "mass($label[$jclus]$jhydr" . "W)=$mass2\n";
#print "$label[$jclus]$jhydr" . "W: $lhas_hydr[$jclus],$#distr_temp2, @distr_temp2\n\n";
							$coll_coef_temp=(8.0*$pi*$boltz*$temperature)**0.5;
							$coll_coef_temp*=(1.0/$mass1+1.0/$mass2)**0.5;
							$coll_coef_temp*=($r1+$r2)**2;
							
							$coll_coef_hydr[$iclus][$jclus][$ihydr][$jhydr] = $coll_coef_temp;
						}
#print "+$distr_temp[$ihydr]*$distr_temp2[$jhydr]*$coll_coef_temp ";
						if($lneutral_cluster[$iclus] && $lneutral_cluster[$jclus]){
							# Neutral-neutral collision rate, weighted with the hydrate distributions
							$coll_coef_eff[$iclus][$jclus] += $distr_temp[$ihydr]*$distr_temp2[$jhydr]*$coll_coef_temp;
							next;
						}
#						a1=.0757; a3=.0015; b0=.0151;b1=-.186;b3=-.0163;
#					for i=1:length(n), for j=i:length(n), AA=6.4e-20/kB/T*4*r(i)*r(j)/(r(i)+r(j))^2; E(i,j)=1+sqrt(AA/3)/(1+b0*sqrt(AA))+b1*log(1+AA)+b3*(log(1+AA))^3; fprintf(fid,'%.14g\t%8s\t%8s\n',E(i,j),[num2str(i),'Ad'],[num2str(j),'Ad']); end, end
################################################################################
################################################################################
##      10.5 Calculating the ion collision enhancement factors                ##
################################################################################
################################################################################

						# Now we have left only the cases with one ion and one neutral. Which is the neutral one?
						if($lneutral_cluster[$iclus]){
							$nclus = $iclus;
							$nhydr = $ihydr;
						}else{
							$nclus = $jclus;
							$nhydr = $jhydr;
						}
					
################################################################################
################################################################################
##         10.51 Su and Bowers, 1973                                          ##
################################################################################
################################################################################
		
						if($ion_coll_method =~ /Su73/i){
							if(!defined($dip_hydr[$nclus][$nhydr])){
								$fcr_temp = -1;
								$missing_fcr += 1;
							}elsif(!$variable_temp){
								$fcr_temp = 9.5436e-29*sqrt($pol_hydr[$nclus][$nhydr])+6.4805e-27*$dip_hydr[$nclus][$nhydr]/sqrt($temperature);
								$fcr_temp *= sqrt(1./$mass1+1./$mass2);
								$fcr_temp /= $coll_coef_temp;
								if($fcr_temp < 1){
									$fcr_temp = 1;
								}
							}else{
								$fcr_temp = sprintf('%.14e',9.5436e-29*sqrt($pol_hydr[$nclus][$nhydr])*sqrt(1./$mass1+1./$mass2));
								$str_temp = sprintf('%.14e',6.4805e-27*$dip_hydr[$nclus][$nhydr]*sqrt(1./$mass1+1./$mass2));
								if($lfortran){
									$fcr_temp =~ s/e/d/;
									$str_temp =~ s/e/d/;
								}
								$fcr_temp = "($fcr_temp+$str_temp/sqrt(temperature))";
							}
					
################################################################################
################################################################################
##         10.52 Su and Chesnavich, 1982                                      ##
################################################################################
################################################################################
		
						}elsif($ion_coll_method =~ /Su82/i){
							if(!defined($dip_hydr[$nclus][$nhydr])){
								$fcr_temp = -1;
								$missing_fcr += 1;
							}elsif(!$variable_temp){

								$fcr_temp = 100.*$dip_hydr[$nclus][$nhydr]/sqrt(2.*$pol_hydr[$nclus][$nhydr]*$boltz*1.e23*$temperature);
								if($fcr_temp<2){
									$fcr_temp = ($fcr_temp+0.5090)**2/10.526+0.9754;
								}else{
									$fcr_temp = 0.4767*$fcr_temp+0.62;
								}
								$fcr_temp *= 2.*$pi*4.8032e-16*sqrt($pol_hydr[$nclus][$nhydr]*(1./$mass1+1./$mass2)/1.e27);
								$fcr_temp /= $coll_coef_temp;
								if($fcr_temp < 1){
									$fcr_temp = 1;
								}
							}elsif($lfortran){
								$temp_val = sprintf('%.14e',(50.*$dip_hydr[$nclus][$nhydr])**2/(2.*$pol_hydr[$nclus][$nhydr]*$boltz*1.e23));
								$temp_val =~ s/e/d/;
								$fcr_temp = "(.5d0+sign(.5d0,temperature-$temp_val))&\n\t\t&";
								$str_temp = "(.5d0-sign(.5d0,temperature-$temp_val))&\n\t\t&";
								$temp_val = sprintf('%.14e',100.*$dip_hydr[$nclus][$nhydr]/sqrt(2.*$pol_hydr[$nclus][$nhydr]*$boltz*1.e23));
								$temp_val =~ s/e/d/;
								$fcr_temp .= "*(($temp_val/sqrt(temperature)+0.5090d0)**2/10.526d0+0.9754d0)";
								$str_temp .= "*(0.4767d0*$temp_val/sqrt(temperature)+0.62d0)";
								$temp_val = sprintf('%.14e',2.*$pi*4.8032e-16*sqrt($pol_hydr[$nclus][$nhydr]*(1./$mass1+1./$mass2)/1.e27));
								$temp_val =~ s/e/d/;
								$fcr_temp = "(&\n\t\t&$fcr_temp&\n\t\t&+$str_temp)&\n\t\t&*$temp_val";
							}else{
								$temp_val = sprintf('%.14e',(50.*$dip_hydr[$nclus][$nhydr])**2/(2.*$pol_hydr[$nclus][$nhydr]*$boltz*1.e23));
								$fcr_temp = "stepfun(temperature,$temp_val)";
								$str_temp = "stepfun($temp_val,temperature)";
								$temp_val = sprintf('%.14e',100.*$dip_hydr[$nclus][$nhydr]/sqrt(2.*$pol_hydr[$nclus][$nhydr]*$boltz*1.e23));
								$fcr_temp .= "*(($temp_val/sqrt(temperature)+0.5090)^2/10.526+0.9754)";
								$str_temp .= "*(0.4767*$temp_val/sqrt(temperature)+0.62)";
								$temp_val = sprintf('%.14e',2.*$pi*4.8032e-16*sqrt($pol_hydr[$nclus][$nhydr]*(1./$mass1+1./$mass2)/1.e27));
								$fcr_temp = "(...\n\t$fcr_temp...\n\t+$str_temp)...\n\t*$temp_val";
							}

################################################################################
################################################################################
##         10.53 Constant enhancement for all collisions                      ##
################################################################################
################################################################################

						}elsif($ion_coll_method =~ m/constant/i){
							if(!$variable_temp){
								$fcr_temp = 10.0;
							}elsif($lfortran){
								$fcr_temp = "10.d0";
							}else{
								$fcr_temp = "10.0";
							}
						}else{
							die "Error: trying to use undefined ion collision method $ion_coll_method.\n";
						}
					
################################################################################
################################################################################
##      10.6 Hydrate averaged collision rate including ion enhancement        ##
################################################################################
################################################################################

						$fcr_hydr[$iclus][$jclus][$ihydr][$jhydr] = $fcr_temp;
						if(!$lhydrate_average){
							$fcr[$iclus][$jclus] = $fcr_temp;
						}else{
							$fcr[$iclus][$jclus] = '';
							# Collision rate with ion enhancement factor, weighted with the hydrate distributions
							if($coll_coef_eff[$iclus][$jclus]>=0 && $fcr_temp>=0){
								$coll_coef_eff[$iclus][$jclus] += $distr_temp[$ihydr]*$distr_temp2[$jhydr]*$coll_coef_temp*$fcr_temp;
							}else{
							# This is for an MCMC simulation where dipole data is missing for some hydrate
								$coll_coef_eff[$iclus][$jclus] = -1;
							}
						}
					}
				}
			}
			if(!$lhydrate_average){
				$coll_coef[$iclus][$jclus] = $coll_coef_hydr[$iclus][$jclus][0][0];
			}

################################################################################
################################################################################
##      10.7 Figuring out sticking factors for each collision                 ##
################################################################################
################################################################################

			if($l_sticking_factor == 1){
				$temp_val = 1.0;
				for($iline=0;$iline<=$#sticking_factor;$iline++){
					if($#{ $sticking_factor[$iline] } == 0 || $sticking_factor[$iline][1] ~~ @molecule_name){
						# if only one sticking factor (and molecule type) is given, 
						# using it only for neutral-neutral collisions
						if($lneutral_cluster[$iclus] && $lneutral_cluster[$jclus]){
							if($#{ $sticking_factor[$iline] } == 0){
								$temp_val = $sticking_factor[$iline][0];
							}elsif($iclus_label =~ /\d+$sticking_factor[$iline][1]\d+/ || $iclus_label =~ /\d+$sticking_factor[$iline][1]$/ || $jclus_label =~ /\d+$sticking_factor[$iline][1]\d+/ || $jclus_label =~ /\d+$sticking_factor[$iline][1]$/){
								$temp_val = $sticking_factor[$iline][0];
							}
						}
					# next option: this is a sticking factor for all ion-neutral collisions
					}elsif($#{ $sticking_factor[$iline] } == 1 && $sticking_factor[$iline][1] eq 'ion-neutral'){
						if($lneutral_cluster[$iclus] + $lneutral_cluster[$jclus] == 1){
							$temp_val = $sticking_factor[$iline][0];
						}
					# if the sticking factor is defined for specific collisions
					}elsif(&compare_clusters($iclus_label,$sticking_factor[$iline][1]) || &compare_clusters($jclus_label,$sticking_factor[$iline][1])){
						if($#{ $sticking_factor[$iline] } == 1){
							$temp_val = $sticking_factor[$iline][0];
						}else{
							if(&compare_clusters($iclus_label,$sticking_factor[$iline][1])){
								$wanted_clus=$jclus_label;
							}else{
								$wanted_clus=$iclus_label;
							}
							if(&compare_clusters($wanted_clus,$sticking_factor[$iline][2])){
								$temp_val = $sticking_factor[$iline][0];
							}
						}
					}
				}
				if($temp_val != 1.){
					$temp_val=sprintf('%.4e',$temp_val);
					if($lfortran){
						$temp_val =~ s/e/d/;
					}
					$sticking_coef[$iclus][$jclus]="$temp_val*";
				}else{
					$sticking_coef[$iclus][$jclus]="";
				}
				$sticking_coef[$jclus][$iclus]=$sticking_coef[$iclus][$jclus];
			}else{
				$sticking_coef[$iclus][$jclus]="";
				$sticking_coef[$jclus][$iclus]="";
			}
		
################################################################################
################################################################################
##      10.8 Checking if this is an MCMC coefficient                          ##
################################################################################
################################################################################
		
			$l_mcmc_coef = 0;
			if($l_coll_factor_mcmc_all){
				$l_mcmc_coef = 1;
			}elsif($l_coll_factor_mcmc_n_n && $lneutral_cluster[$iclus] && $lneutral_cluster[$jclus]){
				$l_mcmc_coef = 1;
			}elsif($l_coll_factor_mcmc_n_i && ($lneutral_cluster[$iclus]+$lneutral_cluster[$jclus] == 1)){
				$l_mcmc_coef = 1;
			}elsif($l_coll_factor_mcmc_n_g && $lneutral_cluster[$iclus] && $jclus>$max_cluster){
				$l_mcmc_coef = 1;
			}elsif($l_coll_factor_mcmc_n_g && $lneutral_cluster[$jclus] && $iclus>$max_cluster){
				$l_mcmc_coef = 1;
			}elsif($l_coll_factor_mcmc_rec && !$lneutral_cluster[$iclus] && !$lneutral_cluster[$jclus]){
				$l_mcmc_coef = 1;
			}
			if($l_mcmc_coef){
				$icoef_coll+=1;
			}
		
		}
				
################################################################################
################################################################################
##      10.9 Printing out the collision rates                                 ##
################################################################################
################################################################################
		
		$value_string = "";
		$coll_string = "$iclus_label + $jclus_label";
		$fixed_string = "";
		$fcr_string = "";
		$hydr_string = "";
		
		if(!$lloop && $lcoll_coef_in[$iclus][$jclus]){
		
			$fixed_string = ", collision rate given as input";
			$value_string = sprintf('%.14e',$coll_coef[$iclus][$jclus]);
			if($lfortran){
				$value_string =~ s/e/d/;
			}
		
		}elsif(!$lloop){
		
			# When using hydrates, print also the corresponding dry value (FCR & collision coefficient) for comparison.
			if($lhydrate_average){
				if($lhas_hydr[$iclus] || $lhas_hydr[$jclus]){
					if($lneutral_cluster[$iclus]+$lneutral_cluster[$jclus] == 1){
						$hydr_string = ", dry value: " . sprintf('%.14e',$fcr_hydr[$iclus][$jclus][0][0]*$coll_coef_hydr[$iclus][$jclus][0][0]) . ", including ion enhancement " . sprintf('%.3g',$fcr_hydr[$iclus][$jclus][0][0]) . "*" . sprintf('%.2e',$coll_coef_hydr[$iclus][$jclus][0][0]);
					}else{
						$hydr_string = ", dry value: " . sprintf('%.14e',$coll_coef_hydr[$iclus][$jclus][0][0]);
					}
					if($sticking_coef[$iclus][$jclus] ne ''){
						$hydr_string = "$hydr_string (sticking coef.  $sticking_coef[$iclus][$jclus])";
					}
				}else{
					$hydr_string = ", no hydrates";
				}
				$value_string = sprintf('%.14e',$coll_coef_eff[$iclus][$jclus]);
				if($lfortran){
					$value_string =~ s/e/d/;
				}
			}else{
				# When using variable_temp, take the T dependence out of the collision rate (not for recombination)
				if($variable_temp && ($lneutral_cluster[$iclus] || $lneutral_cluster[$jclus])){
				
					$str_temp = sprintf('%.14e',$coll_coef[$iclus][$jclus]/$temperature**0.5);
					if($lfortran){
						$str_temp =~ s/e/d/;
					}
					$str_temp = "$str_temp*sqrt(temperature)"; # this is the bare hard-sphere collision rate
				
					if(($lneutral_cluster[$iclus]+$lneutral_cluster[$jclus] == 1) && ($ion_coll_method !~ m/constant/i)){
						$value_string = "max($fcr[$iclus][$jclus],$str_temp)";
						$fcr_string = ", selecting the higher one of ionic rate and hard-sphere rate";
					}else{
						$value_string = $str_temp;
					}
					
				# Standard case: no hydrates and no visible temperature dependence
				}else{
					if($lneutral_cluster[$iclus]+$lneutral_cluster[$jclus] == 1){
						$value_string = sprintf('%.14e',$fcr[$iclus][$jclus]*$coll_coef[$iclus][$jclus]);
						$fcr_string = ", including ion enhancement " . sprintf('%.3g',$fcr[$iclus][$jclus]) . "*" . sprintf('%.2e',$coll_coef[$iclus][$jclus]);
					}else{
						$value_string = sprintf('%.14e',$coll_coef[$iclus][$jclus]);
					}
					if($lfortran){
						$value_string =~ s/e/d/;
					}
				}
			}
			
			# Add a possible sticking coefficient
			$value_string = "$sticking_coef[$iclus][$jclus]$value_string";
			
		}
		
		if(!$lfortran){
			print COLL "K($iclus,$jclus) = $value_string;\t% $coll_string$fixed_string$fcr_string$hydr_string\n";
			print COLL "K($jclus,$iclus) = K($iclus,$jclus);\n" if($iclus != $jclus);
		}elsif(!$lloop){
			$string = "$value_string\t! $coll_string$fixed_string$fcr_string$hydr_string";
			if($l_mcmc_coef){
				# If this is an MCMC coefficient, save the true value to be printed later ...
				$coef_true[$icoef_coll] = $string;
				# ... and print here just the MCMC coefficient
				print OUT "\tK($iclus,$jclus) = coef($icoef_coll)\t! $coll_string\n";
			}else{
				print OUT "\tK($iclus,$jclus) = $string\n";
			}
			print OUT "\tK($jclus,$iclus) = K($iclus,$jclus)\n" if($iclus != $jclus);
		}
		
	}  # end of the jclus loop
}  # end of the iclus loop

if($lfortran){
	print OUT "\nend subroutine get_coll$append_to_subroutine_names\n\n!-----------------------------------------------------------\n\n";
}else{
	if($variable_temp){
		print COLL "\nend";
	}
	close(COLL);
}

################################################################################
################################################################################
##   11. Evaporation rates                                                    ##
################################################################################
################################################################################

if(!$lfortran){
	open(EVAP, ">get_evap$append_to_file_names.m");
	if($variable_temp){
		print EVAP "function E = get_evap(temperature,K)\n\n";
	}
}
$missing_energies = 0;
if(!$disable_evap){

	if($lfortran){
		print OUT "subroutine get_evap$append_to_subroutine_names(E";
		if($l_mcmc_any){
			print OUT ",coef";
		}elsif($variable_temp){
			print OUT ",K,temperature";
		}
		print OUT ")\n\timplicit none\n";
		if(!$lloop){
			print OUT "\treal(kind(1.d0)) :: E($max_cluster_number,$max_cluster_number)";
			if($variable_temp){
				print OUT ", K($max_cluster_number,$max_cluster_number), temperature";
			}elsif($l_mcmc_any){
				print OUT ", coef(" . ($n_mcmc_coefs+$variable_cs+2*$variable_ion_source) . ")";
			}		
			print OUT "\n\n\t! evaporation coefficients\n";
			print OUT $evap_comment;
			print OUT "\n\tE = 0.d0\n";
		}
	}elsif(!$lfortran){
		print EVAP "% Evaporation coefficients\n";
		print EVAP $evap_comment;
		print EVAP "\nE = nan($max_cluster_number,$max_cluster_number);\n\n";
	}


################################################################################
################################################################################
##      11.1 Function for computing evaporation rates when using loops        ##
################################################################################
################################################################################

	if($lloop){
		print OUT "\tinteger, parameter :: nclust = $max_cluster_number\n";
		print OUT "\treal(kind(1.d0)) :: E(";
		if($disable_nonmonomer_evaps && $nmol_types_used==1){
			print OUT "1";
		}else{
			print OUT "nclust";
		}
		print OUT ",nclust), K(nclust,nclust), temperature\n";
		if($nmol_types_used>1){
			print OUT "\tinteger, parameter :: ij_ind_max($nmol_types_used) = (/$nmol_ranges_list/)\n";
			print OUT "\tinteger :: indices(nclust,$nmol_types_used), ij_ind($nmol_types_used)";
			if($lloop_max_ratio){
				print OUT ", n_pure($nmol_ranges_array[$nmol_max_ratio_basis])";
			}
			print OUT "\n";
		}
		print OUT "\tinteger :: i, j, ij, imol\n";
		print OUT "\treal(kind(1.d0)) :: m, r, xmol\n";
		
		$temp_val = $boltz*$temperature;
		$string1 = sprintf('%.14e',$temp_val);
		$string1 =~ s/e/d/;
		
		$temp_val = $kcal_per_mol_to_J/$boltz/$temperature;
		$string2 = sprintf('%.14e',$temp_val);
		$string2 =~ s/e/d/;

################################################################################
################################################################################
##   11.11 Evaporation rates from DeltaGs                                     ##
################################################################################
################################################################################
		
		if($loop_evap_coef =~ m/^DeltaG$/i){
			# Evaporation coefficients via DeltaGs and detailed balance
			
			print OUT "\treal(kind(1.d0)) :: gi, gj, gij, p\n\n";
			
			print OUT "\tcall get_coll$append_to_subroutine_names(K)\n";
			if($nmol_types_used>1){
				print OUT "\tcall get_molecule_numbers$append_to_subroutine_names(indices)\n";
				if($lloop_max_ratio){
					print OUT "\tcall pure_clust_indices$append_to_subroutine_names(n_pure)\n";
				}
			}
			
			if($disable_nonmonomer_evaps && $nmol_types_used==1){
				print OUT "\tdo i = 1,nclust-1\n";
				if($rlim_no_evap !~ /^\s*$/){
					print OUT "\t\tcall get_masses_and_radii$append_to_subroutine_names(m,r,i+1)\n";
					print OUT "\t\tif(r .ge. $rlim_no_evap_str) then\n";
					print OUT "\t\t\texit\n";
					print OUT "\t\tend if\n";
				}
				print OUT "\t\tif(i .eq. 1) then\n";
				print OUT "\t\t\tcall cluster_energies(gi,p,1,i)\n";
				print OUT "\t\t\tgj = gi\n";
				print OUT "\t\telse\n";
				print OUT "\t\t\tgi = gij\n";
				print OUT "\t\tend if\n";
				print OUT "\t\tcall cluster_energies(gij,p,1,i+1)\n";
				print OUT "\t\tE(1,i) = K(1,i)*p/$string1*exp((gij-gi-gj)*$string2)\n";
				print OUT "\t\tif (i .eq. 1) then\n";
				print OUT "\t\t\tE(1,i) = .5d0*E(1,i)\n";
				print OUT "\t\tend if\n";
				print OUT "\tend do\n\n";
			}else{
				print OUT "\tdo i = 1,nclust\n";
				if($nmol_types_used>1){
					print OUT "\t\tcall cluster_energies(gi,p,$nmol_types_used,indices(i,:))\n";
				}else{
					print OUT "\t\tcall cluster_energies(gi,p,1,i)\n";
				}
				print OUT "\t\tdo j = i,nclust\n";
				if($nmol_types_used>1){
					print OUT "\t\t\tij_ind = indices(i,:)+indices(j,:)\n";
					if($rlim_no_evap !~ /^\s*$/){
						print OUT "\t\t\tcall get_masses_and_radii$append_to_subroutine_names(m,r,ij_ind)\n";
						print OUT "\t\t\tif(r .ge. $rlim_no_evap_str) then\n";
						print OUT "\t\t\t\tcycle\n";
						print OUT "\t\t\tend if\n";
					}
					print OUT "\t\t\tif(";
					if($lloop_max_ratio){
						print OUT "(ij_ind($nmol_not_basis).le.$max_wrt_mon*ij_ind($nmol_max_ratio_basis)) .and. ";
					}
					if($lloop_cs){
						print OUT "all(ij_ind .lt. ij_ind_max)) then\n";
					}else{
						print OUT "all(ij_ind .le. ij_ind_max)) then\n";
					}
					print OUT "\t\t\t\tij = $cluster_from_indices\n";
					print OUT "\t\t\t\tcall cluster_energies(gj,p,$nmol_types_used,indices(j,:))\n";
					print OUT "\t\t\t\tcall cluster_energies(gij,p,$nmol_types_used,ij_ind)\n";
				}else{
					print OUT "\t\t\tij=i+j\n";
					if($rlim_no_evap !~ /^\s*$/){
						print OUT "\t\tcall get_masses_and_radii$append_to_subroutine_names(m,r,ij)\n";
						print OUT "\t\tif(r .ge. $rlim_no_evap_str) then\n";
						print OUT "\t\t\texit\n";
						print OUT "\t\tend if\n";
					}
					print OUT "\t\t\tif(ij .le. nclust) then\n";
					print OUT "\t\t\t\tcall cluster_energies(gj,p,1,j)\n";
					print OUT "\t\t\t\tcall cluster_energies(gij,p,1,ij)\n";
				}
				print OUT "\t\t\t\tE(i,j) = K(i,j)*p/$string1*exp((gij-gi-gj)*$string2)\n";
				print OUT "\t\t\t\tif(i .eq. j) then\n";
				print OUT "\t\t\t\t\tE(i,j) = .5d0*E(i,j)\n";
				print OUT "\t\t\t\tend if\n";
				print OUT "\t\t\t\tE(j,i) = E(i,j)\n";
				print OUT "\t\t\tend if\n";
				print OUT "\t\tend do\n";
				print OUT "\tend do\n\n";
			}
			
################################################################################
################################################################################
##   11.12 Evaporation rates from the Kelvin formula                          ##
################################################################################
################################################################################
			
		}elsif($loop_evap_coef =~ m/^Kelvin$/i){
		    # No fissions, monomer evaporation according to the Kelvin equation
			# = an approximation based on liquid drop DeltaG and non-discrete cluster limit

			# Suggestions for values for parameters psat and sigma that produce evaporation frequencies of realistic orders of magnitude:
			# sigma is assumed to be temperature-independent and to correspond to the sulfuric acid value of 0.05 N/m
			# psat is assumed to be temperature-dependent, but the dependence is a very rough approximation, for now:
			# The parameters are set so that at 280 K, the value is 5e-11 Pa which produces evaporation rates around QC values for H2SO4 and NH3;
			# the functional form is psat=exp(A/T+B) for H2SO4 vapor pressure over solid ammonium sulfate from Marti et al., JGR 102, 3725-3735, 1997, 
			# in which A = -5928 (+-891) and B = -3.77 (+-2.69);
			# however, here B is set to B = -2.55 in order to produce the wanted values at 280 K (an often used temperature e.g. in CLOUD).
			#if($temperature == 280){
			#	$temp_val = 5e-11;
			#}else{
			#	$temp_val = exp(-5928/$temperature-2.55);
			#}
			#$string3 = sprintf('%.14e',$temp_val);
			#$string3 =~ s/e/d/; # psat
			
			print OUT "\tinteger :: n_monomers($nmol_types_used) = (/$monomer_indices/)\n";
			
			$str_temp = "";
			$str_temp2 = "";
			$str_temp3 = "";
			for($imol=1;$imol<=$nmol_types_used;$imol++){
				$str_temp .= sprintf('%.14e',$psat[$mol_types_used[$imol]]).",";
				$str_temp2 .= sprintf('%.14e',$mass[$mol_types_used[$imol]]/$density[$mol_types_used[$imol]]).",";
				if($imol==1){
					$sigma = $surface_tension[$mol_types_used[$imol]];
					$str_temp3 = sprintf('%.14e',$sigma);
				}else{
					if($surface_tension[$mol_types_used[$imol]] != $sigma){
						print "The surface tension is not the same for all compounds\n";
						print "-> Using the one for the first compound (no parameterization for a mixture implemented yet)\n";
					}
				}
			}
			$str_temp =~ s/,$//;
			$str_temp =~ s/e/d/g;
			$str_temp2 =~ s/,$//;
			$str_temp2 =~ s/e/d/g;
			$str_temp3 =~ s/e/d/g;
			if($nmol_types_used==1){
				print OUT "\treal(kind(1.d0)), parameter :: psat = $str_temp\n";
				print OUT "\treal(kind(1.d0)), parameter :: monvol = $str_temp2\n";
			}else{
				print OUT "\treal(kind(1.d0)), parameter :: psat($nmol_types_used) = (/$str_temp/)\n";
				print OUT "\treal(kind(1.d0)), parameter :: monvol($nmol_types_used) = (/$str_temp2/)\n";
			}
			print OUT "\treal(kind(1.d0)), parameter :: sigma = $str_temp3\n\n";

			print OUT "\tE = 0.d0\n\n";
			print OUT "\tcall get_coll$append_to_subroutine_names(K)\n";
			if($nmol_types_used>1){
				print OUT "\tcall get_molecule_numbers$append_to_subroutine_names(indices)\n";
				if($lloop_max_ratio){
					print OUT "\tcall pure_clust_indices$append_to_subroutine_names(n_pure)\n";
				}
			}

			if($nmol_types_used==1){
				print OUT "\tdo i = 1,nclust-1\n";
				print OUT "\t\tcall get_masses_and_radii$append_to_subroutine_names(m,r,i+1)\n";
				if($rlim_no_evap !~ /^\s*$/){
					print OUT "\t\tif(r .ge. $rlim_no_evap_str) then\n";
					print OUT "\t\t\texit\n";
					print OUT "\t\tend if\n";
				}
				print OUT "\t\tE(1,i) = K(1,i)*psat/$string1*exp((2.d0*sigma*monvol)/($string1*r))\n";
				print OUT "\t\tif (i .eq. 1) then\n";
				print OUT "\t\t\tE(1,i) = .5d0*E(1,i)\n";
				print OUT "\t\tend if\n";
				print OUT "\tend do\n\n";
			}else{
				# Loop over the evaporation of each component
				print OUT "\tdo imol = 1,$nmol_types_used\n";
				print OUT "\t\tj = n_monomers(imol)\n";
				print OUT "\t\tdo i = 1,nclust\n";
				print OUT "\t\t\tij_ind = indices(i,:)+indices(j,:)\n";
				print OUT "\t\t\tcall get_masses_and_radii$append_to_subroutine_names(m,r,ij_ind)\n";
				if($rlim_no_evap !~ /^\s*$/){
					print OUT "\t\t\tif(r .ge. $rlim_no_evap_str) then\n";
					print OUT "\t\t\t\tcycle\n";
					print OUT "\t\t\tend if\n";
				}
				print OUT "\t\t\tif(";
				if($lloop_max_ratio){
					print OUT "(ij_ind($nmol_not_basis).le.$max_wrt_mon*ij_ind($nmol_max_ratio_basis)) .and. ";
				}
				if($lloop_cs){
					print OUT "all(ij_ind .lt. ij_ind_max)) then\n";
				}else{
					print OUT "all(ij_ind .le. ij_ind_max)) then\n";
				}
				#print OUT "\t\t\t\tij = $cluster_from_indices\n";
				print OUT "\t\t\t\txmol = real(ij_ind(imol),kind=kind(1.d0))/real(sum(ij_ind),kind=kind(1.d0))\n";
				print OUT "\t\t\t\tE(i,j) = K(i,j)*psat(imol)*xmol/$string1*exp((2.d0*sigma*monvol(imol))/($string1*r))\n";
				print OUT "\t\t\t\tif(i .eq. j) then\n";
				print OUT "\t\t\t\t\tE(i,j) = .5d0*E(i,j)\n";
				# For heterodimers, the Kelvin approximation gives different results depending on which
				# molecule is assumed to be the "particle" -> use the higher value for the evaporation rate
				print OUT "\t\t\t\telse if(E(j,i) .gt. E(i,j)) then\n";
				print OUT "\t\t\t\t\tE(i,j) = E(j,i)\n";
				print OUT "\t\t\t\tend if\n";
				print OUT "\t\t\t\tE(j,i) = E(i,j)\n";
				print OUT "\t\t\tend if\n";
				print OUT "\t\tend do\n";
				print OUT "\tend do\n\n";
			}
			
		}else{
			die "\nCan't understand evaporation coefficient definition $loop_evap_coef.\n\n";
		}
		
	}else{
	
################################################################################
################################################################################
##      11.2 Reading in rates for specific evaporations                       ##
################################################################################
################################################################################

		if($evaporation_coef_file_name ne ''){
			open(EVAPORATION, "<$evaporation_coef_file_name") || die "\nCould not open evaporation rate file $evaporation_coef_file_name\n\n";
			chomp(@temp_array=<EVAPORATION>);
			close(EVAPORATION);
			$number_of_lines=$#temp_array;
			$iline2 = 0;
			for($iline=0;$iline<=$number_of_lines;$iline++){
				next if ($temp_array[$iline] =~ /^\#/);
				if($temp_array[$iline] =~ /(\d+$name_water-?\d+)|(\d+$name_water\s*$)/){
					die "$evaporation_coef_file_name: Evaporation rate cannot be given for hydrates for now (if hydrate averaging is used, the given rate for dry clusters is taken to be the effective rate)";
				}
				@columns=split(' ',$temp_array[$iline]);			
				die "\n$columns[0] on line $iline in file $evaporation_coef_file_name should be a number\n\n" if (!&is_pos_number($columns[0]));
				die "\nThere should be two or three columns on line $iline in file $evaporation_coef_file_name: the rate and the evaporating cluster, or the rate and the mother cluster and the evaporator\n\n" if ($#columns < 1 || $#columns > 2);
				$evaporation_coef_in[$iline2][0] = $columns[0];
				push(@{ $evaporation_coef_in[$iline2] },$columns[1]);
				if($#columns>1){
					push(@{ $evaporation_coef_in[$iline2] },$columns[2]);
				}
				$iline2++;
			}
		}

	################################################################################
	################################################################################
	##      11.3 Reading in scaling factors for the evaporation rates             ##
	################################################################################
	################################################################################

		if(($scale_evap_factor[0][0] == 0) && ($scale_evap_file_name eq '')){
			$l_scale_evap = 0;
		}else{
			$l_scale_evap = 1;
			if($scale_evap_file_name ne ''){
				print "evaporation rate scaling file: $scale_evap_file_name\n";
				open(SCALE, "<$scale_evap_file_name") || die "\nCould not open evaporation rate scaling file $scale_evap_file_name\n\n";
				chomp(@temp_array=<SCALE>);
				close(SCALE);
				$number_of_lines=$#temp_array;
				if($scale_evap_factor[0][0] == 0){
					$iline2 = 0;
				}else{
					$iline2 = 1;
				}
				for($iline=0;$iline<=$number_of_lines;$iline++){
					next if ($temp_array[$iline] =~ /^\#/);
					@columns=split(' ',$temp_array[$iline]);
					die "\n$columns[0] on line $iline in file $scale_evap_file_name should be a number\n\n" if (!&is_pos_number($columns[0]));
					$scale_evap_factor[$iline2][0] = $columns[0];
					if($#columns>0){
						push(@{ $scale_evap_factor[$iline2] },$columns[1]);
						if($#columns>1){
							push(@{ $scale_evap_factor[$iline2] },$columns[2]);
						}
					}
					$iline2++;
				}
			}
		}

	################################################################################
	################################################################################
	##      11.4 Calculating the evaporation rates                                ##
	################################################################################
	################################################################################

		# Looping through the possible evaporations
		for($kclus=1;$kclus<=$max_cluster;$kclus++){
			for($i=1;$i<=2*$ind_lin_loss[$kclus][0];$i+=2){
				$jclus = $ind_lin_loss[$kclus][$i];
				$iclus = $ind_lin_loss[$kclus][$i+1];
			
				# Skipping, if this is a loss term
				if($iclus>$max_cluster || $jclus>$max_cluster){
					next;
				}
			
				$iclus_label=$label[$iclus];
				$jclus_label=$label[$jclus];
				$kclus_label=$label[$kclus];
			
	################################################################################
	################################################################################
	##         11.41 Determining if this evaporation has been given               ##
	##               a fixed rate value in the evaporation rate file              ##
	################################################################################
	################################################################################
			
				$levap_coef_in[$iclus][$jclus] = 0;
		
				if($evaporation_coef_file_name ne ''){
					for($iline=0;$iline<=$#evaporation_coef_in;$iline++){
						if(&compare_clusters($kclus_label,$evaporation_coef_in[$iline][1])){
							if(($#{ $evaporation_coef_in[$iline] } == 1) || (&compare_clusters($iclus_label,$evaporation_coef_in[$iline][2]) || &compare_clusters($jclus_label,$evaporation_coef_in[$iline][2]))){
								$levap_coef_in[$iclus][$jclus] = 1;
								$evap_coef[$iclus][$jclus] = $evaporation_coef_in[$iline][0];
								last;
							}
						}
					}
				}
			
				if($levap_coef_in[$iclus][$jclus]){
			
					$evap_coef_eff[$iclus][$jclus] = $evap_coef[$iclus][$jclus];
			
	################################################################################
	################################################################################
	##         11.42 Calculating the free energy -based evaporation rate          ##
	################################################################################
	################################################################################

				}else{

					# make sure we have the free energies (not necessary when using mcmc for evaporation)
					if(!$l_evap_factor_mcmc){
						die "\nMissing free energy for $iclus_label.\n\n" if($delta_g[$iclus] eq '');
						die "\nMissing free energy for $jclus_label.\n\n" if($delta_g[$jclus] eq '');
						die "\nMissing free energy for $kclus_label.\n\n" if($delta_g[$kclus] eq '');
					}
			
	################################################################################
	################################################################################
	##            11.421 Figuring out scaling factors for each evaporation        ##
	################################################################################
	################################################################################
			
					if($l_scale_evap == 1){
						$temp_val = 0;
						for($iline=0;$iline<=$#scale_evap_factor;$iline++){
							if($#{ $scale_evap_factor[$iline] } == 0){
								$temp_val = $scale_evap_factor[$iline][0];
							}elsif(&compare_clusters($iclus_label,$scale_evap_factor[$iline][1]) || &compare_clusters($jclus_label,$scale_evap_factor[$iline][1])){
								if(($#{ $scale_evap_factor[$iline] } == 1) || &compare_clusters($kclus_label,$scale_evap_factor[$iline][2])){
									$temp_val = $scale_evap_factor[$iline][0];
								}
							}
						}
						if($temp_val != 0){
							$evap_coef_scale[$iclus][$jclus]=$temp_val*$kcal_per_mol_to_J;
							$evap_coef_scale[$jclus][$iclus]=$evap_coef_scale[$iclus][$jclus];
						}
					}

					if($l_evap_factor_mcmc && (!defined $delta_g[$iclus] || !defined $delta_g[$jclus] || !defined $delta_g[$kclus])){
						$evap_coef[$iclus][$jclus] = -1;
						$missing_energies += 1;
					}elsif(!$variable_temp){
						# calculate the evaporation rates of all possible combinations of hydrates of these dry clusters
						@distr_temp = @{$distr[$iclus]};
						@distr_temp2 = @{$distr[$jclus]};
						@distr_temp3 = @{$distr[$kclus]};
				
						# Effective evaporation rate
						$evap_coef_eff[$iclus][$jclus]=0.0;
				
						for($khydr=0;$khydr<=$#distr_temp3;$khydr++){
							$lfound = 0;
							for($ihydr=0;$ihydr<=$#distr_temp;$ihydr++){
								if(!(defined $distr_temp[$ihydr])){
									next;
								}
								$jhydr = $khydr-$ihydr;
								if($jhydr<0 || $jhydr>$#distr_temp2 || !(defined $distr_temp2[$jhydr])){
									next;
								# Don't count the same evaporation twice
								}elsif($iclus==$jclus && $jhydr<$ihydr){
									next;
								}else{
									$lfound = 1;
								}
							
								if($lcoll_coef_in[$iclus][$jclus] && ($ihydr>0 || $jhydr>0)){
									die "Cannot give one collision rate for $iclus_label and $jclus_label if you want to average the evaporation rate over hydrates";
								}
							
								if($ihydr != 0){
									$temp_label=$iclus_label.$ihydr.$name_water;
								}else{
									$temp_label=$iclus_label;
								}
								if($jhydr != 0){
									$temp_label2=$jclus_label.$jhydr.$name_water;
								}else{
									$temp_label2=$jclus_label;
								}
					
								$evap_coef_temp=$delta_g_hydr[$kclus][$khydr]-$delta_g_hydr[$iclus][$ihydr]-$delta_g_hydr[$jclus][$jhydr];

								if(defined $evap_coef_scale[$iclus][$jclus]){
									$evap_coef_temp+=$evap_coef_scale[$iclus][$jclus];
								}

								$evap_coef_temp /= ($boltz*$temperature);
								$evap_coef_temp += log($pressure*1.0/$boltz/$temperature);
								$evap_coef_temp = exp($evap_coef_temp)*$coll_coef_hydr[$iclus][$jclus][$ihydr][$jhydr];
					
								# Divide by two if the evaporators are similar
								if($iclus == $jclus && $ihydr == $jhydr){
									$evap_coef_temp *= 0.5;
								}
								$evap_coef_hydr[$iclus][$jclus][$ihydr][$jhydr] = $evap_coef_temp;
							
	################################################################################
	################################################################################
	##      11.5 Hydrate averaged evaporation rate                                ##
	################################################################################
	################################################################################

						
								# Evaporation rate without ion enhancement, weighted with the hydrate distribution
								$temp_val2 = $distr_temp3[$khydr]*$evap_coef_temp;
								# Ion enhancement for ion-neutral collisions
								if($lneutral_cluster[$iclus]+$lneutral_cluster[$jclus] == 1){
									$temp_val2 *= $fcr_hydr[$iclus][$jclus][$ihydr][$jhydr];
								}
								$evap_coef_eff[$iclus][$jclus] += $temp_val2;
						

							}
							if(!$lfound){
								$string = $iclus_label . "i$name_water + " . $jclus_label . "j$name_water";
								die "No possible evaporation products for $kclus_label$khydr$name_water -> $string.\n";
							}
						}
						if(!$lhydrate_average){
							$evap_coef[$iclus][$jclus] = $evap_coef_hydr[$iclus][$jclus][0][0];
						}else{
							$evap_coef[$iclus][$jclus] = $evap_coef_eff[$iclus][$jclus];
						}
				
	################################################################################
	################################################################################
	##      11.6 Temperature dependent evaporation rate                           ##
	################################################################################
	################################################################################

					# If we are using a variable temperature, we don't have a hydrate distribution.
					}else{
						if($lfortran){
							$string = "&\n\t\t\t &";
						}else{
							$string = " ...\n\t\t\t";
						}
						if($delta_g[$iclus] ne '0'){
							$string1 = "$string-$delta_g[$iclus]";
						}else{
							$string1 = '';
						}
						if($delta_g[$jclus] ne '0'){
							$string2 = "$string-$delta_g[$jclus]";
						}else{
							$string2 = '';
						}
						$evap_coef[$iclus][$jclus]="$delta_g[$kclus]$string1$string2";
						if(defined $evap_coef_scale[$iclus][$jclus]){
							$str_temp = sprintf('%.14e',$evap_coef_scale[$iclus][$jclus]);
							if($lfortran){
								$str_temp =~ s/e/d/;
							}
							$evap_coef[$iclus][$jclus]="$evap_coef[$iclus][$jclus]+$str_temp";
						}
						$temp_val = $kcal_per_mol_to_J/$boltz;
						$str_temp = sprintf('%.14e',$temp_val);
						if($lfortran){
							$str_temp =~ s/e/d/;
						}
						$evap_coef[$iclus][$jclus] = "exp($string($evap_coef[$iclus][$jclus])$string*$str_temp)";
						$temp_val = $pressure/$boltz;
						$str_temp = sprintf('%.14e',$temp_val);
						if($lfortran){
							$str_temp =~ s/e/d/;
						}
						$evap_coef[$iclus][$jclus] = "$str_temp/temperature*$evap_coef[$iclus][$jclus]*K($iclus,$jclus)";
	
						#### Divide by two if the evaporators are similar, since we want to have the evaporation coefficient from the point of view of the evaporator
						if($iclus == $jclus){
							$evap_coef[$iclus][$jclus] = "0.5*$evap_coef[$iclus][$jclus]";
						}else{
							$evap_coef[$jclus][$iclus] = $evap_coef[$iclus][$jclus];
						}
					}
				
				}	
			
	################################################################################
	################################################################################
	##      11.7 Printing out the evaporation rates                               ##
	################################################################################
	################################################################################

				# Make the 0.5 factor visible for clarity, if products are similar (already done for variable_temp)
			
				$value_string = "";
				$evap_string = "$kclus_label -> $iclus_label + $jclus_label";
				$fixed_string = "";
				$fcr_string = "";
				$hydr_string = "";
			
				if($levap_coef_in[$iclus][$jclus]){
			
					$fixed_string = ", evaporation rate given as input";
					$value_string = sprintf('%.14e',$evap_coef[$iclus][$jclus]);
					if($lfortran){
						$value_string =~ s/e/d/;
					}
			
				}else{
			
					if($lcoll_coef_in[$iclus][$jclus]){
						$fixed_string = ", collision rate given as input";
					}
			
					# When using hydrates, print also the dry value (FCR & collision coefficient) for comparison.
					if($lhydrate_average){
						if($lhas_hydr[$kclus]){
							if($iclus == $jclus){
								$hydr_string = ", dry value: 0.5*" . sprintf('%.14e',2*$evap_coef_hydr[$iclus][$jclus][0][0]);
							}else{
								if($lneutral_cluster[$iclus]+$lneutral_cluster[$jclus] == 1){
									$hydr_string = ", dry value: " . sprintf('%.14e',$fcr_hydr[$iclus][$jclus][0][0]*$evap_coef_hydr[$iclus][$jclus][0][0]) . ", including ion enhancement " . sprintf('%.3g',$fcr_hydr[$iclus][$jclus][0][0]) . "*" . sprintf('%.2e',$evap_coef_hydr[$iclus][$jclus][0][0]);
								}else{
									$hydr_string = ", dry value: " . sprintf('%.14e',$evap_coef_hydr[$iclus][$jclus][0][0]);
								}
							}
							if($sticking_coef[$iclus][$jclus] ne ''){
								$hydr_string = "$hydr_string (sticking coef.  $sticking_coef[$iclus][$jclus])";
							}
						}else{
							$hydr_string = ", no hydrates";
						}
						if($iclus == $jclus){
							$value_string = sprintf('%.14e',2*$evap_coef_eff[$iclus][$jclus]);
							if($lfortran){
								$value_string =~ s/e/d/;
								$value_string = ".5d0*$value_string";
							}else{
								$value_string = ".5*$value_string";
							}
						}else{
							$value_string = sprintf('%.14e',$evap_coef_eff[$iclus][$jclus]);
							if($lfortran){
								$value_string =~ s/e/d/;
							}
						}
					}else{					
						if($variable_temp){
							$value_string = $evap_coef[$iclus][$jclus];
							if($lneutral_cluster[$iclus]+$lneutral_cluster[$jclus] == 1){
								$fcr_string = ", including ion enhancement" if(!$lcoll_coef_in[$iclus][$jclus]);
							}					
						}elsif($iclus == $jclus){
							$value_string = sprintf('%.14e',2*$evap_coef[$iclus][$jclus]);
							if($lfortran){
								$value_string =~ s/e/d/;
								$value_string = "0.5d0*$value_string";
							}else{
								$value_string = "0.5*$value_string";
							}
						}else{
							if($lneutral_cluster[$iclus]+$lneutral_cluster[$jclus] == 1){
								$value_string = sprintf('%.14e',$fcr[$iclus][$jclus]*$evap_coef[$iclus][$jclus]);
								if(!$lcoll_coef_in[$iclus][$jclus]){
									$fcr_string = ", including ion enhancement " . sprintf('%.3g',$fcr[$iclus][$jclus]) . "*" . sprintf('%.2e',$evap_coef[$iclus][$jclus]);
								}
							}else{
								$value_string = sprintf('%.14e',$evap_coef[$iclus][$jclus]);
							}
						
							if($lfortran){
								$value_string =~ s/e/d/;
							}
						
						}
					
					}
				
					# Add a possible sticking coefficient
					$value_string = "$sticking_coef[$iclus][$jclus]$value_string" if(!$lcoll_coef_in[$iclus][$jclus]);
			
				}
			
				if($lfortran){
					$string = "$value_string\t! $evap_string$fixed_string$fcr_string$hydr_string";
					if($l_evap_factor_mcmc == 1){
						$icoef_evap += 1;
						# If this is an MCMC coefficient, save the true value to be printed later ...
						if($evap_coef[$iclus][$jclus]<0){
							$coef_true[$ncoef_coll+$icoef_evap] = "-1.d0";
						}else{
							$coef_true[$ncoef_coll+$icoef_evap] = $string;
						}
						# ... and print here just the MCMC coefficient
						print OUT "\tE($iclus,$jclus) = coef(" . ($ncoef_coll+$icoef_evap) . ")\t! $evap_string\n";
					}else{
						print OUT "\tE($iclus,$jclus) = $string\n";
					}
					print OUT "\tE($jclus,$iclus) = E($iclus,$jclus)\n" if($iclus != $jclus);
				}else{
					print EVAP "E($iclus,$jclus) = $value_string;\t% $evap_string$fixed_string$fcr_string$hydr_string\n";
					print EVAP "E($jclus,$iclus) = E($iclus,$jclus);\n" if($iclus != $jclus);
				}
			}
		}
	}
	if($lfortran){
		print OUT "\nend subroutine get_evap$append_to_subroutine_names\n";
	}
}else{
	if(!$lfortran){
		print EVAP "% Evaporation disabled\n\nE = [];\n";
	}
}

if(!$lfortran){
	if($variable_temp){
		print EVAP "\nend";
	}
	close(EVAP);
}

################################################################################
################################################################################
##   12. Printing the ion loss enhancement factors for Fortran                ##
################################################################################
################################################################################

# In the Matlab version, these are in the get_cs and get_wl files

if($include_ion_terms && $lfortran && ($parameters_in_get_fcr_1 ne '')){
	print OUT "\n!-----------------------------------------------------------\n\n";
	print OUT "subroutine get_fcr$append_to_subroutine_names($parameters_in_get_fcr_1)\n\timplicit none\n";
	print OUT "\treal(kind(1.d0)) :: $parameters_in_get_fcr_3\n\n";
	print OUT "\t! enhancement factors for ion losses\n\n";
	
	if(!$disable_coag_sinks){
		if($fcs =~ m/\./){
			print OUT "\tfcs = $fcs"."d0\n";
		}else{
			print OUT "\tfcs = $fcs.d0\n";
		}
	}
	if(!$disable_wall_terms || ($l_mcmc_any && $ncoef_wl>0)){
		if($fwl =~ m/\./){
			print OUT "\tfwl = $fwl"."d0\n";
		}else{
			print OUT "\tfwl = $fwl.d0\n";
		}
	}	
	print OUT "\nend subroutine get_fcr$append_to_subroutine_names\n";
}

################################################################################
################################################################################
##   13. "True" coefficient values for MCMC                                   ##
################################################################################
################################################################################
	
if($lfortran && ($l_coll_factor_mcmc_n_i || $l_coll_factor_mcmc_n_g || ($l_evap_factor_mcmc && $missing_energies<$ncoef_evap))){
	print OUT "\n!-----------------------------------------------------------\n\n";
	print OUT "subroutine coefs_from_data$append_to_subroutine_names(coef)";
	print OUT "\n\timplicit none\n\treal(kind(1.d0)) :: coef(" . ($n_mcmc_coefs+$variable_cs+2*$variable_ion_source) . ")";
	if($variable_cs || $variable_ion_source){
		print OUT "\n\n\tcoef(1:" . ($n_mcmc_coefs) . ") = 1.d0\n";
	}else{
		print OUT "\n\n\tcoef = 0.d0\n";
	}
	if($l_coll_factor_mcmc_any){
		print OUT "\n\t! " . '"true"' . " values of the collision coefficients (based on hard spheres and possibly ion enhancement)\n";
		if($missing_fcr>0){
			print OUT "\n\t! (-1 for coefficients for which dipole data are missing)\n";
		}
	}
	for($i=1;$i<=$ncoef_coll;$i++){
		if($coef_true[$i] ne ''){
			print OUT "\tcoef($i) = $coef_true[$i]\n";
		}
	}
	if($l_evap_factor_mcmc){
		print OUT "\n\t! " . '"true"' . " values of the evaporation coefficients (based on input energies)\n";
		if($missing_energies+$missing_fcr>0){
			print OUT "\n\t! (-1 for coefficients for which energies or dipole data are missing)\n";
		}
	}
	for($i=$ncoef_coll+1;$i<=$ncoef_coll+$ncoef_evap;$i++){
		if($coef_true[$i] ne ''){
			print OUT "\tcoef($i) = $coef_true[$i]\n";
		}
	}
	if($l_q_neg_mcmc){
		print OUT "\n\tcoef(" . ($ncoef_coll+$ncoef_evap+1) . ") = 3.d6\t\t\t\t\t! typical GCR ion production rate in m^-3s^-1";
	}
	if($l_q_pos_mcmc){
		print OUT "\n\tcoef(" . ($ncoef_coll+$ncoef_evap+$l_q_neg_mcmc+1) . ") = 3.d6\t\t\t\t\t! typical GCR ion production rate in m^-3s^-1";
	}
	if($l_wl_mcmc && $ncoef_wl==1){
		print OUT "\n\tcoef($n_mcmc_coefs) = $cs_value\t\t\t\t! example external loss (Hyytiälä coagulation sink)";
	}elsif($l_mcmc_any && $ncoef_wl>0){
		for($i=$n_mcmc_coefs-$ncoef_wl+1;$i<=$n_mcmc_coefs;$i++){
			if($coef_true[$i] ne ''){
				print OUT "\n\tcoef($i) = $coef_true[$i]";
			}else{
				print OUT "\n\tcoef($i) = $cs_value\t\t\t\t! example external loss (Hyytiälä coagulation sink)";
			}
		}
	}
	print OUT "\n\nend subroutine coefs_from_data$append_to_subroutine_names";
}

################################################################################
################################################################################
##   14. Finishing the Matlab driver file                                     ##
################################################################################
################################################################################

if(!$lfortran){
	
	print DRIVER "\n\n% Initializing concentrations etc.\n\n";
	print DRIVER "\tC0 = zeros(1,$max_n_flux);\n\n\tisconst = zeros($max_cluster_number,1);\n\tisfitted = zeros($max_cluster,1);\n\tC_fit = cell($max_cluster,2);\n\tsource = zeros($max_cluster_number,1);\n\n";

	print DRIVER "% Checking what input we have and processing it\n\n";
	print DRIVER "\tif size(Tmax,1)>size(Tmax,2), Tmax = Tmax'; end\n";
	if(!$l_increase_time){
		print DRIVER "\tfixed_time = 1;\n";
	}else{
		print DRIVER "\tfixed_time = 0;\n";
	}
	if(!$l_old_output){
		print DRIVER "\tif nargout<=5\n";
		print DRIVER "\t\tno_fluxes = 1;\n";
		print DRIVER "\t\tno_j = 1;\n";
		print DRIVER "\t\tcalc_sources = 0;\n";
		print DRIVER "\telseif nargout<=6\n";
		print DRIVER "\t\tno_fluxes = 1;\n";
		print DRIVER "\t\tno_j = 0;\n";
		print DRIVER "\t\tcalc_sources = 0;\n";
		print DRIVER "\telse\n";
		print DRIVER "\t\tno_fluxes = 0;\n";
		print DRIVER "\t\tno_j = 0;\n";
		print DRIVER "\t\tcalc_sources = 1;\n";
		print DRIVER "\tend\n";
	}else{
		print DRIVER "\tno_fluxes = 0;\n\tno_j = 0;\n";
	}
	print DRIVER "\tsources_in = 0;\n\tsources_out = 0;\n\tconsts_out = 0;\n\topt_str = '';\n\tc_neutr = 0;\n\tc_neg = 0;\n\tc_pos = 0;\n\tfl_out = 0;\n\toutmat = 0;\n\tcluster_data = 0;\n\tC_fun_str = {};\n\tn_fun = [];\n";
	
	$line = "\tC0_in = 0;\n\tCSfactor = 1;\n\tWLfactor = 1;\n";
	if($variable_temp){
		$line = $line . "\ttemperature = $temperature;\n";
	}
	$line = $line . "\tif ~isempty(varargin)\n";
	$line = $line . "\t\ti = 0;\n";
	$line = $line . "\t\tif isnumeric(varargin{1})\t% Initial concentrations cm^3 -> m^3\n";
	$line = $line . "\t\t\tC00 = varargin{1}*1e6;\n";
	$line = $line . "\t\t\tC0_in = 1;\n";
	$line = $line . "\t\t\ti = 1;\n";
	$line = $line . "\t\t\tif size(C00,2)==$size_C_out\n";
	if(! $l_save_outgoing){
		$line = $line . "\t\t\t\tC00 = [C00 zeros(size(C00,1)," . ($max_n_flux-$size_C_out) . ")];\n";
	}
	$line = $line . "\t\t\t\tC0 = C00(end,:);\n";
	$line = $line . "\t\t\t\tif (length(varargin)>1 && isnumeric(varargin{2}))\t% Time (vector) corresponding to C00\n";
	$line = $line . "\t\t\t\t\tT00 = varargin{2};\n";
	$line = $line . "\t\t\t\t\tif (isvector(T00) && numel(T00)==size(C00,1))\n";
	$line = $line . "\t\t\t\t\t\ti = 2;\n";
	$line = $line . "\t\t\t\t\t\tif size(T00,2)>size(T00,1), T00 = T00'; end\n";
	$line = $line . "\t\t\t\t\t\tif numel(Tmax)==1\n";
	$line = $line . "\t\t\t\t\t\t\tTmax = Tmax+T00(end);\n";
	$line = $line . "\t\t\t\t\t\tend;\n";
	$line = $line . "\t\t\t\t\telse\n";
	$line = $line . "\t\t\t\t\t\terror('Sizes of initial concentration array and corresponding time array don''t match.')\n";
	$line = $line . "\t\t\t\t\tend\n";
	$line = $line . "\t\t\t\telse\n";
	$line = $line . "\t\t\t\t\tT00 = 0;\n";
	$line = $line . "\t\t\t\tend\n";
	$line = $line . "\t\t\telse\n";
	$line = $line . "\t\t\t\terror('Bad initial concentrations.')\n";
	$line = $line . "\t\t\tend\n";
	$line = $line . "\t\telse\n";
	$line = $line . "\t\t\tT00 = 0;\n";
	$line = $line . "\t\t\tC00 = C0;\n";
	$line = $line . "\t\tend\n";
	print DRIVER $line;
	$line = "\t\t% Check for single keywords in the input arguments\n";
	$line = $line . "\t\tif ismember(1,strcmpi('fixed_time',varargin))\n";
	$line = $line . "\t\t\tfixed_time = 1;\n";
	$line = $line . "\t\t\tvarargin = varargin(~strcmpi('fixed_time',varargin));\n";
	$line = $line . "\t\telseif ismember(1,strcmpi('repeat',varargin))\n";
	$line = $line . "\t\t\tfixed_time = 0;\n";
	$line = $line . "\t\t\tvarargin = varargin(~strcmpi('repeat',varargin));\n";
	$line = $line . "\t\telseif ismember(1,strcmpi('no_dofluxes',varargin))\n";
	$line = $line . "\t\t\tno_fluxes = 1;\n";
	$line = $line . "\t\t\tvarargin = varargin(~strcmpi('no_dofluxes',varargin));\n";
	$line = $line . "\t\telseif ismember(1,strcmpi('no_fluxes',varargin))\n";
	$line = $line . "\t\t\tno_fluxes = 1;\n";
	$line = $line . "\t\t\tno_j = 1;\n";
	$line = $line . "\t\t\tvarargin = varargin(~strcmpi('no_fluxes',varargin));\n";
	$line = $line . "\t\tend\n";
	$line = $line . "\t\t% Go through the paired input arguments\n";
	$line = $line . "\t\tfor j = i+1:2:length(varargin)\n";
	$line = $line . "\t\t\tif strcmpi(varargin{j},'Sources_in')\n";
	$line = $line . "\t\t\t\tsources_in = 1;\n";
	$line = $line . "\t\t\t\tfname1 = varargin{j+1};\n";
	$line = $line . "\t\t\telseif strcmpi(varargin{j},'Sources_out')\n";
	$line = $line . "\t\t\t\tsources_out = 1;\n";
	$line = $line . "\t\t\t\tfname2 = varargin{j+1};\n";
	$line = $line . "\t\t\t\tno_fluxes = 0;\n";
	$line = $line . "\t\t\t\tcalc_sources = 1;\n";
	$line = $line . "\t\t\telseif strcmpi(varargin{j},'Constants_out')\n";
	$line = $line . "\t\t\t\tconsts_out = 1;\n";
	$line = $line . "\t\t\t\tfname2c = varargin{j+1};\n";
	$line = $line . "\t\t\telseif strcmpi(varargin{j},'Options')\n";
	$line = $line . "\t\t\t\topt_str = varargin{j+1};\n";
	$line = $line . "\t\t\telseif strcmpi(varargin{j},'C_neutr')\n";
	$line = $line . "\t\t\t\tc_neutr = 1;\n";
	$line = $line . "\t\t\t\tfname3 = varargin{j+1};\n";
	$line = $line . "\t\t\telseif strcmpi(varargin{j},'C_neg')\n";
	$line = $line . "\t\t\t\tc_neg = 1;\n";
	$line = $line . "\t\t\t\tfname4 = varargin{j+1};\n";
	$line = $line . "\t\t\telseif strcmpi(varargin{j},'C_pos')\n";
	$line = $line . "\t\t\t\tc_pos = 1;\n";
	$line = $line . "\t\t\t\tfname5 = varargin{j+1};\n";
	$line = $line . "\t\t\telseif strcmpi(varargin{j},'Fluxes')\n";
	$line = $line . "\t\t\t\tfl_out = 1;\n";
	$line = $line . "\t\t\t\tfname6 = varargin{j+1};\n";
	$line = $line . "\t\t\t\tno_fluxes = 0;\n";
	$line = $line . "\t\t\t\tcalc_sources = 1;\n";
	$line = $line . "\t\t\telseif strcmpi(varargin{j},'CSfactor')\n";
	$line = $line . "\t\t\t\tCSfactor = varargin{j+1};\n";
	$line = $line . "\t\t\telseif strcmpi(varargin{j},'WLfactor')\n";
	$line = $line . "\t\t\t\tWLfactor = varargin{j+1};\n";
	$line = $line . "\t\t\telseif strcmpi(varargin{j},'Outmat')\n";
	$line = $line . "\t\t\t\toutmat = 1;\n";
	$line = $line . "\t\t\t\tfname7 = varargin{j+1};\n";
	$line = $line . "\t\t\t\tno_fluxes = 0;\n";
	$line = $line . "\t\t\telseif strcmpi(varargin{j},'Cluster_data')\n";
	$line = $line . "\t\t\t\tcluster_data = 1;\n";
	$line = $line . "\t\t\t\tfname8 = varargin{j+1};\n";
	$line = $line . "\t\t\t\tno_fluxes = 0;\n";
	$line = $line . "\t\t\t\tcalc_sources = 1;\n";
	$line = $line . "\t\t\telseif strcmpi(varargin{j},'Cfun')\n";
	$line = $line . "\t\t\t\tn_fun = [n_fun strmatch(varargin{j+1}{1},clust,'exact')];\n";
	$line = $line . "\t\t\t\tisconst(n_fun(end)) = 1;\n";
	$line = $line . "\t\t\t\tif length(varargin{j+1})<2 || isempty(varargin{j+1}{2})\n";
	$line = $line . "\t\t\t\t\terror(['The function related to cluster ',varargin{j+1}{1},' was not given.'])\n";
	$line = $line . "\t\t\t\tend\n";
	$line = $line . "\t\t\t\tC_fun_str{n_fun(end)} = ['1e6*', varargin{j+1}{2},'('];\n";
	$line = $line . "\t\t\t\tC_fun_str_2{n_fun(end)} = ['1e6*', varargin{j+1}{2},'('];\n";
	$line = $line . "\t\t\t\tC_fun_str_3{n_fun(end)} = ['1e6*', varargin{j+1}{2},'('];\n";
	$line = $line . "\t\t\t\tfor k = 3:length(varargin{j+1})\n";
	$line = $line . "\t\t\t\t\tif strcmpi(varargin{j+1}{k},'t')\n";
	$line = $line . "\t\t\t\t\t\tC_fun_str{n_fun(end)} = [C_fun_str{n_fun(end)}, 't'];\n";
	$line = $line . "\t\t\t\t\t\tC_fun_str_2{n_fun(end)} = [C_fun_str_2{n_fun(end)}, 'T(j)'];\n";
	$line = $line . "\t\t\t\t\t\tC_fun_str_3{n_fun(end)} = [C_fun_str_3{n_fun(end)}, 'T'];\n";
	$line = $line . "\t\t\t\t\telse\n";
	$line = $line . "\t\t\t\t\t\tC_fun_str{n_fun(end)} = [C_fun_str{n_fun(end)}, '1e-6*c(', num2str(strmatch(varargin{j+1}{k},clust,'exact')), ')'];\n";
	$line = $line . "\t\t\t\t\t\tC_fun_str_2{n_fun(end)} = [C_fun_str_2{n_fun(end)}, '1e-6*C(j,', num2str(strmatch(varargin{j+1}{k},clust,'exact')), ')'];\n";
	$line = $line . "\t\t\t\t\t\tC_fun_str_3{n_fun(end)} = [C_fun_str_3{n_fun(end)}, '1e-6*C(j,', num2str(strmatch(varargin{j+1}{k},clust,'exact')), ')'];\n";
	$line = $line . "\t\t\t\t\tend\n";
	$line = $line . "\t\t\t\t\tif k<length(varargin{j+1})\n";
	$line = $line . "\t\t\t\t\t\tC_fun_str{n_fun(end)} = [C_fun_str{n_fun(end)}, ','];\n";
	$line = $line . "\t\t\t\t\t\tC_fun_str_2{n_fun(end)} = [C_fun_str_2{n_fun(end)}, ','];\n";
	$line = $line . "\t\t\t\t\t\tC_fun_str_3{n_fun(end)} = [C_fun_str_3{n_fun(end)}, ','];\n";
	$line = $line . "\t\t\t\t\tend\n\t\t\t\tend\n";
	$line = $line . "\t\t\t\tC_fun_str{n_fun(end)} = [C_fun_str{n_fun(end)}, ')'];\n";
	$line = $line . "\t\t\t\tC_fun_str_2{n_fun(end)} = [C_fun_str_2{n_fun(end)}, ')'];\n";
	$line = $line . "\t\t\t\tC_fun_str_3{n_fun(end)} = [C_fun_str_3{n_fun(end)}, ')'];\n";
	if($variable_temp){
		$line = $line . "\t\t\telseif strcmpi(varargin{j},'Temperature')\n";
		$line = $line . "\t\t\t\ttemperature = varargin{j+1};\n";
	}

	$line = $line . "\t\t\tend\n";
	$line = $line . "\t\tend\n";
	$line = $line . "\tend\n";
	print DRIVER $line;

	print DRIVER "\tif ~isempty(n_fun)\n";
	print DRIVER "\t\tn_fun = sort(n_fun);\n";
	print DRIVER "\t\tdC_fun_str = 'dcdt = \@($input_for_eq) $output_function_name(t,[c(1:';\n";
	print DRIVER "\t\tfor i = 1:length(n_fun)\n";
	print DRIVER "\t\t\tdC_fun_str = [dC_fun_str, num2str(n_fun(i)-1),');',C_fun_str{n_fun(i)}, ';c(', num2str(n_fun(i)+1), ':'];\n";
	print DRIVER "\t\tend\n";
	print DRIVER "\t\tdC_fun_str = [dC_fun_str, 'end)],$input_for_eq_C_fun);'];\n";
	print DRIVER "\t\teval(dC_fun_str);\n";
	print DRIVER "\telse\n";
	print DRIVER "\t\tdcdt = \@($input_for_eq) $output_function_name($input_for_eq);\n";
	print DRIVER "\tend\n\n";
	
	print DRIVER "\tif numel(Tmax)>1\n";
	print DRIVER "\t\tif Tmax(1)==T00(end)\n";
	print DRIVER "\t\t\tTmax = Tmax(2:end);\n";
	print DRIVER "\t\tend\n";
	print DRIVER "\tend\n\n";
	
	print DRIVER "\tif cluster_data\n";
	print DRIVER "\t\tsave(fname8, 'clust', 'diameters', 'mobility_diameters', 'masses', 'monomers', 'nonmonomers', 'neutrals', 'negatives', 'positives', 'clust_flux','" . $molecule_name[0] . "_in_clust')\n";
	print DRIVER "\tend\n";
	
	print DRIVER "\n% If input source file was given, reading its content and converting cm^3 -> m^3\n\tif sources_in";
	print DRIVER "\n\n\t\tsfid = fopen(fname1, 'r');\n";
	print DRIVER "\t\tif sfid == -1\n\t\terror('Source file not found!')\n\t\tend\n";
	print DRIVER "\t\tA = textscan(sfid, '\%s \%s \%f \%s', 'commentStyle', '\%');\n\t\tfclose(sfid);\n\n";

	print DRIVER "\t\tfor i=1:$max_cluster_number\n";
	print DRIVER "\t\t\tj = strmatch(clust{i},A{2},'exact');\n";
	print DRIVER "\t\t\tdefined = 0;\n";
	print DRIVER "\t\t\tfor k=1:length(j)\n";
	print DRIVER "\t\t\t\tif (~C0_in && ~isempty(strmatch('i',A{1}(j(k)))))\n";
	print DRIVER "\t\t\t\t\tif (defined == 1 || defined == 3), error('Source term multiply defined.'), end\n";
	print DRIVER "\t\t\t\t\tC0(i) = A{3}(j(k))*1e6;\n";
	print DRIVER "\t\t\t\t\tdefined = 1;\n";
	print DRIVER "\t\t\t\telseif ~isempty(strmatch('s',A{1}(j(k))))\n";
	print DRIVER "\t\t\t\t\tif (defined > 1), error('Source term multiply defined.'), end\n";
	print DRIVER "\t\t\t\t\tsource(i) = A{3}(j(k))*1e6;\n";
	print DRIVER "\t\t\t\t\tdefined = 2;\n";
	print DRIVER "\t\t\t\telseif ~isempty(strmatch('c',A{1}(j(k))))\n";
	print DRIVER "\t\t\t\t\tif (defined > 0), error('Source term multiply defined.'), end\n";
	print DRIVER "\t\t\t\t\tif length(A{4})<j(k) || isempty(A{4}{j(k)})\n";
	print DRIVER "\t\t\t\t\t\tisconst(i) = 1;\n";
	print DRIVER "\t\t\t\t\t\tC0(i) = A{3}(j(k))*1e6;\n";
	print DRIVER "\t\t\t\t\telse\n";
	print DRIVER "\t\t\t\t\t\tC_fit{i,1} = A{3}(j(k))*1e6;\n";
	print DRIVER "\t\t\t\t\t\tC_fit{i,2} = [];\n";
	print DRIVER "\t\t\t\t\t\tfor i2=1:$max_cluster_number\n";
	print DRIVER "\t\t\t\t\t\t\tfor l=1:length(strmatch(clust{i2},regexp(A{4}{j(k)},'-','split'),'exact'))\n";
	print DRIVER "\t\t\t\t\t\t\t\tC_fit{i,2} = [C_fit{i,2},i2];\n";
	print DRIVER "\t\t\t\t\t\t\tend\n";
	print DRIVER "\t\t\t\t\t\tend\n";
	print DRIVER "\t\t\t\t\t\tif ~isempty(C_fit{i,2})\n";
	print DRIVER "\t\t\t\t\t\t\tisfitted(i) = 1;\n";
	print DRIVER "\t\t\t\t\t\telse\n";
	print DRIVER "\t\t\t\t\t\t\tisconst(i) = 1;\n";
	print DRIVER "\t\t\t\t\t\tend\n";
	print DRIVER "\t\t\t\t\t\tC0(i) = A{3}(j(k))*1e6-sum(C0(C_fit{i,2}));\n";
	print DRIVER "\t\t\t\t\tend\n";
	print DRIVER "\t\t\t\t\tdefined = 3;\n";
	print DRIVER "\t\t\t\tend\n\t\t\tend\n\t\tend\n\n";
	if($include_generic_neg){
		print DRIVER "\t\tj = strmatch('nion',A{1});\n";
		print DRIVER "\t\tif length(j) > 1\n\t\t\terror('Ion source term multiply defined.')\n\t\tend\n";
		print DRIVER "\t\tif ~isempty(j),\n";
		print DRIVER "\t\t\tif source($n_generic_neg) ~= 0 || isconst($n_generic_neg) == 1\n";
		print DRIVER "\t\t\t\terror('Ion source term multiply defined.')\n\t\t\tend\n";
		print DRIVER "\t\t\tsource($n_generic_neg) = A{3}(j)*1e6;\n";
		print DRIVER "\t\tend\n";
	}
	if($include_generic_pos){
		print DRIVER "\t\tj = strmatch('pion',A{1});\n";
		print DRIVER "\t\tif length(j) > 1\n\t\t\terror('Ion source term multiply defined.')\n\t\tend\n";
		print DRIVER "\t\tif ~isempty(j),\n";
		print DRIVER "\t\t\tif source($n_generic_pos) ~= 0 || isconst($n_generic_pos) == 1\n";
		print DRIVER "\t\t\t\terror('Ion source term multiply defined.')\n\t\t\tend\n";
		print DRIVER "\t\t\tsource($n_generic_pos) = A{3}(j)*1e6;\n";
		print DRIVER "\t\tend\n";
	}
	print DRIVER "\t\tC00(end,:) = C0;\n";
	print DRIVER "\t\tj = strmatch('wall',A{1});\n";
	print DRIVER "\t\tif length(j) > 1\n\t\t\terror('Wall loss enhancement multiply defined.')\n\t\tend\n";
	print DRIVER "\t\tif ~isempty(j),\n";
	print DRIVER "\t\t\tfwl_in = A{3}(j);\n";
	print DRIVER "\t\tend\n";
	print DRIVER "\t\tj = strmatch('coag',A{1});\n";
	print DRIVER "\t\tif length(j) > 1\n\t\t\terror('Coagulation loss enhancement multiply defined.')\n\t\tend\n";
	print DRIVER "\t\tif ~isempty(j),\n";
	print DRIVER "\t\t\tfcs_in = A{3}(j);\n";
	print DRIVER "\t\tend\n";
	print DRIVER "\telseif C0_in == 0\n\t\terror('Input concentrations or source terms must be given either in the vector C0 or the file ''Sources_in''!');\n\tend\n\n";
	
	print DRIVER "% Setting the options for the solver\n";
	print DRIVER "\teval(['options = odeset(',opt_str,');']);\n\n";
	
	if(!$variable_temp){
		if($l_reac_in_clus){
			print DRIVER "% Reading in the reaction, collision and evaporation rates, coagulation losses and wall losses (enhancement factors for ion losses included in the loss rate files)\n";
			print DRIVER "\tget_reac$append_to_file_names;\n";
		}else{
			print DRIVER "% Reading in the collision and evaporation rates, coagulation losses and wall losses (enhancement factors for ion losses included in the loss rate files)\n";
		}
		print DRIVER "\tget_coll$append_to_file_names;\n\tget_evap$append_to_file_names;\n\tget_cs$append_to_file_names;\n\tget_wl$append_to_file_names;\n";
	}else{
		if($l_reac_in_clus){
			print DRIVER "% Reading in the reaction, collision and evaporation rates, coagulation losses and wall losses (enhancement factors for ion losses included in the loss rate files)\n";
			print DRIVER "\tR = get_reac$append_to_file_names(temperature);\n";
		}else{
			print DRIVER "% Reading in the collision and evaporation rates, coagulation losses and wall losses (enhancement factors for ion losses included in the loss rate files)\n";
		}
		print DRIVER "\tK = get_coll$append_to_file_names(temperature);\n\tE = get_evap$append_to_file_names(temperature,K);\n\tget_cs$append_to_file_names;\n\tget_wl$append_to_file_names;\n";
	}
	
	print DRIVER "\n\tCS = CS*CSfactor;\n\tWL = WL*WLfactor;\n";
	
	if($include_ion_terms){
		print DRIVER "\tif exist('fcs_in','var')\n\t\tfcs = fcs_in;\n\tend\n";
		print DRIVER "\tif exist('fwl_in','var')\n\t\tfwl = fwl_in;\n\tend\n";
	}
	
	print DRIVER "\n% Solving the birth-death equations\n";

	print DRIVER "\tconverged = 0;\n";
	print DRIVER "\ti = 0;\n";
	print DRIVER "\tTfin = Tmax;\n";
	print DRIVER "\tif fixed_time\n";
	print DRIVER "\t\timax = 1;\n";
	print DRIVER "\telse\n";
	print DRIVER "\t\timax = 5;\n";
	print DRIVER "\t\tTadd = Tmax-T00(end);\n";
	print DRIVER "\tend\n";
	print DRIVER "\twhile ~converged\n";
	print DRIVER "\t\ti = i+1;\n";
	print DRIVER "% Checking for negative concentrations\n";
	print DRIVER "\t\tif min(min(C00))<-1e-6\n";
	print DRIVER "\t\t\tconverged = -1;\n";
	print DRIVER "\t\t\tCf = nan(size(C0));\n";
	print DRIVER "\t\t\tC = C00*1e-6;\n";
	print DRIVER "\t\t\tT = T00;\n";
	print DRIVER "\t\t\tflux = nan;\n";
	print DRIVER "\t\t\tout = nan;\n";
	print DRIVER "\t\t\toutflux_matrix = nan;\n";
	print DRIVER "\t\t\tsources = nan;\n";
	print DRIVER "\t\t\tform = nan;\n";
	print DRIVER "\t\t\tJ_out = nan;\n";
	print DRIVER "\t\t\treturn;\n";
	print DRIVER "\t\tend\n";
	print DRIVER "\t\tif i>imax\n";
	print DRIVER "% Not converged but exiting anyway\n";
	if($l_old_output){
		print DRIVER "\t\t\tif nargout < 11\n";
	}else{
		print DRIVER "\t\t\tif nargout < 3\n";
	}
	print DRIVER "\t\t\t\tdisp('Not converging, try a larger Tmax');\n";
	print DRIVER "\t\t\tend\n";
	print DRIVER "\t\t\tbreak;\n";
	print DRIVER "\t\telseif i>1\n";
	print DRIVER "% Adding more time in order to reach convergence\n";
    print DRIVER "\t\t\tTfin = T00(end)+Tadd;\n";
	print DRIVER "\t\tend\n";	
	print DRIVER "\t\t[T,C] = ode15s(@(t,c) dcdt($input_for_eq), [T00(end) Tfin], C0, options);\n";
	print DRIVER "\t\tif i == 1 && min(min(C)) < 0;\n";
	print DRIVER "\t\t\t[T5,C5] = ode15s(@(t,c) dcdt($input_for_eq), [T00(end) T00(end)+(Tfin(1)-T00(end))*1e-6], C0, options);\n";
	print DRIVER "\t\t\t[T,C] = ode15s(@(t,c) dcdt($input_for_eq), [T5(end) Tfin], C5(end,:), options);\n";
	print DRIVER "\t\t\tif numel(Tmax)==1\n";
	print DRIVER "\t\t\t\tT = [T5; T(2:end)];\n";	
	print DRIVER "\t\t\t\tC = [C5; C(2:end,:)];\n";
	print DRIVER "\t\t\tend\n";
	print DRIVER "\t\tend\n";
	print DRIVER "\t\tCf = C(end,:);\n";
	print DRIVER "\t\tif max(abs(Cf(1:$max_cluster_number)-C0(1:$max_cluster_number))./max(Cf(1:$max_cluster_number),1e-6))<1e-6\n\t\t\tconverged = 1;\n";
	print DRIVER "\t\telseif length(T)>2 && max(abs(Cf(1:$max_cluster_number)-C(floor(length(T)*.75),1:$max_cluster_number))./max(Cf(1:$max_cluster_number),1e-6))<1e-6\n\t\t\tconverged = 1;\n\t\tend\n";
	print DRIVER "\t\tif any(isfitted)\n";
	print DRIVER "\t\t\tfor j = 1:$max_cluster\n";
	print DRIVER "\t\t\t\tif isfitted(j)\n";
	print DRIVER "\t\t\t\t\tif abs(C_fit{j,1}-Cf(j)-sum(Cf(C_fit{j,2})))/C_fit{j,1} > 1e-6\n";
	print DRIVER "\t\t\t\t\t\tCf([j C_fit{j,2}]) = Cf([j C_fit{j,2}])*C_fit{j,1}/sum(Cf([j C_fit{j,2}]));\n";
	print DRIVER "\t\t\t\t\t\tconverged = 0;\n";
	print DRIVER "\t\t\t\t\tend\n";
	print DRIVER "\t\t\t\tend\n";
	print DRIVER "\t\t\tend\n";
	print DRIVER "\t\tend\n";
	if($include_generic_neg && $include_generic_pos){
		if($generic_positive_ions_corr ne ''){
			print DRIVER "\t\tif abs(Cf($n_generic_pos)+sum(Cf(positives))-Cf($n_generic_neg)-sum(Cf(negatives)))/max(Cf($n_generic_pos),1e-6) > 1e-6\n";
			print DRIVER "\t\t\tif Cf($n_generic_neg)+sum(Cf(negatives))-sum(Cf(positives))>0\n";
			print DRIVER "\t\t\t\tCf($n_generic_pos) = Cf($n_generic_neg)+sum(Cf(negatives))-sum(Cf(positives));\n";
			print DRIVER "\t\t\telse\n";
			print DRIVER "\t\t\t\tCf($n_generic_neg) = sum(Cf(positives))-sum(Cf(negatives));\n";
			print DRIVER "\t\t\t\tCf($n_generic_pos) = 0;\n";
			print DRIVER "\t\t\tend\n";
			print DRIVER "\t\t\tconverged = 0;\n";
			print DRIVER "\t\tend\n\n";
		}elsif($generic_negative_ions_corr ne ''){
			print DRIVER "\t\tif abs(Cf($n_generic_pos)+sum(Cf(positives))-Cf($n_generic_neg)-sum(Cf(negatives)))/max(Cf($n_generic_neg),1e-6) > 1e-6\n";
			print DRIVER "\t\t\tif Cf($n_generic_pos)+sum(Cf(positives))-sum(Cf(negatives))>0\n";
			print DRIVER "\t\t\t\tCf($n_generic_neg) = Cf($n_generic_pos)+sum(Cf(positives))-sum(Cf(negatives));\n";
			print DRIVER "\t\t\telse\n";
			print DRIVER "\t\t\t\tCf($n_generic_pos) = sum(Cf(negatives))-sum(Cf(positives));\n";
			print DRIVER "\t\t\t\tCf($n_generic_neg) = 0;\n";
			print DRIVER "\t\t\tend\n";
			print DRIVER "\t\t\tconverged = 0;\n";
			print DRIVER "\t\tend\n\n";
		}elsif($generic_positive_ions ne ''){
			print DRIVER "\t\tif abs(Cf($n_generic_pos)-Cf($n_generic_neg)-sum(Cf(negatives)))/max(Cf($n_generic_pos),1e-6) > 1e-6\n";
			print DRIVER "\t\t\tCf($n_generic_pos) = Cf($n_generic_neg)+sum(Cf(negatives));\n";
			print DRIVER "\t\t\tconverged = 0;\n";
			print DRIVER "\t\tend\n\n";
		}elsif($generic_negative_ions ne ''){
			print DRIVER "\t\tif abs(Cf($n_generic_pos)+sum(Cf(positives))-Cf($n_generic_neg))/max(Cf($n_generic_neg),1e-6) > 1e-6\n";
			print DRIVER "\t\t\tCf($n_generic_neg) = Cf($n_generic_pos)+sum(Cf(positives));\n";
			print DRIVER "\t\t\tconverged = 0;\n";
			print DRIVER "\t\tend\n\n";
		}
	}
	print DRIVER "\t\tC0 = Cf;\n";
	print DRIVER "\t\tC00 = [C00; C(2:end,:)];\n";
	print DRIVER "\t\tT00 = [T00; T(2:end)];\n";
	print DRIVER "\tend\n";
	print DRIVER "\tif ~C0_in && numel(Tmax)==1\n";
	print DRIVER "\t\tC = C00(2:end,:);\n";
	print DRIVER "\t\tT = T00(2:end);\n";
	print DRIVER "\telse\n";
	print DRIVER "\t\tC = C00;\n";
	print DRIVER "\t\tT = T00;\n";
	print DRIVER "\tend\n";
	print DRIVER "\tfor i=n_fun\n";
	print DRIVER "\t\tC(:,i) = eval(C_fun_str_3{i});\n";
	print DRIVER "\tend\n\n";

	print DRIVER "\% Printing table of neutral concentrations\n\tif(c_neutr)\n";
	$ncols = 1;
	for($i=1;$i<$n_neutral_mols;$i++){
		$ncols = $ncols*($n_max_type[$neutral_mols[$i]]+1);
	}
	print DRIVER "\t\tt0 = zeros($n_max_type[0],$ncols);\n";
	print DRIVER "\t\tfor i = 1:length(neutrals)\n\t\t\tt0(table_0{i}(1),table_0{i}(2)) = neutrals(i);\n\t\tend\n";
	print DRIVER "\t\tcfid = fopen(fname3, 'w');\n";
	print DRIVER "\t\tfprintf(cfid,'\%\% ";
	$line = '';
	for($icol=1;$icol<=$ncols;$icol++){
		$kmol = $icol-1;
		$iline = '';
		for($i=1;$i<$n_neutral_mols;$i++){
			$imol = $neutral_mols[$n_neutral_mols-$i];
			$nmol = 1;
			for($j=$i+1;$j<$n_neutral_mols;$j++){
				$jmol = $neutral_mols[$n_neutral_mols-$j];
				$nmol = $nmol*($n_max_type[$jmol]+1);
			}
			$iline = int($kmol/$nmol) . "$molecule_name[$imol]$iline";
			$kmol = $kmol-int($kmol/$nmol)*$nmol;
		}
		$line = "$line, '$iline'";
		print DRIVER "\%-9s ";
	}
	print DRIVER "\\n'$line);\n";
	print DRIVER "\t\tfor i=1:size(t0,1)\n\t\t\tfor j=1:size(t0,2)\n";
	print DRIVER "\t\t\t\tif t0(i,j)<1\n\t\t\t\t\tfprintf(cfid,'	  nan ');\n";
	print DRIVER "\t\t\t\telse\n\t\t\t\t\tfprintf(cfid,'\%.3e ',Cf(t0(i,j))*1e-6);\n";
	print DRIVER "\t\t\t\tend\n\t\t\tend\n\t\t\tfprintf(cfid,'\\n');\n\t\tend\n\t\tfclose(cfid);\n\tend\n\n";
	
	if($include_negative_ion_terms){
		print DRIVER "\% Printing table of negative concentrations\n\tif(c_neg)\n";
		$ncols = 1;
		for($i=1;$i<$n_neutral_mols;$i++){
			$ncols = $ncols*($n_max_type_neg[$neutral_mols[$i]]+1);
		}
		print DRIVER "\t\ttn = zeros($n_max_type_neg[0],$ncols);\n";
		print DRIVER "\t\tfor i = 1:length(negatives)\n\t\t\ttn(table_neg{i}(1),table_neg{i}(2)) = negatives(i);\n\t\tend\n";
		print DRIVER "\t\tcfid = fopen(fname4, 'w');\n";
		print DRIVER "\t\tfprintf(cfid,'\%\% ";
		$line = '';
		for($icol=1;$icol<=$ncols;$icol++){
			$kmol = $icol-1;
			$iline = '';
			for($i=1;$i<$n_neutral_mols;$i++){
				$imol = $neutral_mols[$n_neutral_mols-$i];
				$nmol = 1;
				for($j=$i+1;$j<$n_neutral_mols;$j++){
					$jmol = $neutral_mols[$n_neutral_mols-$j];
					$nmol = $nmol*($n_max_type_neg[$jmol]+1);
				}
				$iline = int($kmol/$nmol) . "$molecule_name[$imol]$iline";
				$kmol = $kmol-int($kmol/$nmol)*$nmol;
			}
			$line = "$line, '$iline'";
			print DRIVER "\%-9s ";
		}
		print DRIVER "\\n'$line);\n";
		print DRIVER "\t\tfor i=1:size(tn,1)\n\t\t\tfor j=1:size(tn,2)\n";
		print DRIVER "\t\t\t\tif tn(i,j)<1\n\t\t\t\t\tfprintf(cfid,'	  nan ');\n";
		print DRIVER "\t\t\t\telse\n\t\t\t\t\tfprintf(cfid,'\%.3e ',Cf(tn(i,j))*1e-6);\n";
		print DRIVER "\t\t\t\tend\n\t\t\tend\n\t\t\tfprintf(cfid,'\\n');\n\t\tend\n\t\tfclose(cfid);\n\tend\n\n";
	}
	
	if($include_positive_ion_terms){
		print DRIVER "\% Printing table of positive concentrations\n\tif(c_pos)\n";
		$ncols = 1;
		for($i=1;$i<$n_neutral_mols;$i++){
			$ncols = $ncols*($n_max_type_pos[$neutral_mols[$i]]+1);
		}
		print DRIVER "\t\ttp = zeros($n_max_type_pos[0],$ncols);\n";
		print DRIVER "\t\tfor i = 1:length(positives)\n\t\t\ttp(table_pos{i}(1),table_pos{i}(2)) = positives(i);\n\t\tend\n";
		print DRIVER "\t\tcfid = fopen(fname5, 'w');\n";
		print DRIVER "\t\tfprintf(cfid,'\%\% ";
		$line = '';
		for($icol=1;$icol<=$ncols;$icol++){
			$kmol = $icol-1;
			$iline = '';
			for($i=1;$i<$n_neutral_mols;$i++){
				$imol = $neutral_mols[$n_neutral_mols-$i];
				$nmol = 1;
				for($j=$i+1;$j<$n_neutral_mols;$j++){
					$jmol = $neutral_mols[$n_neutral_mols-$j];
					$nmol = $nmol*($n_max_type_pos[$jmol]+1);
				}
				$iline = int($kmol/$nmol) . "$molecule_name[$imol]$iline";
				$kmol = $kmol-int($kmol/$nmol)*$nmol;
			}
			$line = "$line, '$iline'";
			print DRIVER "\%-9s ";
		}
		print DRIVER "\\n'$line);\n";
		print DRIVER "\t\tfor i=1:size(tp,1)\n\t\t\tfor j=1:size(tp,2)\n";
		print DRIVER "\t\t\t\tif tp(i,j)<1\n\t\t\t\t\tfprintf(cfid,'	  nan ');\n";
		print DRIVER "\t\t\t\telse\n\t\t\t\t\tfprintf(cfid,'\%.3e ',Cf(tp(i,j))*1e-6);\n";
		print DRIVER "\t\t\t\tend\n\t\t\tend\n\t\t\tfprintf(cfid,'\\n');\n\t\tend\n\t\tfclose(cfid);\n\tend\n\n";
	}
	
	print DRIVER "\tif ~no_j\n";
	print DRIVER "\t\tJ_out = nan(1,length(T));\n";
	print DRIVER "\t\tfor i=1:length(T)\n";
	if(!$variable_temp){
		print DRIVER "\t\t\tJ_out(i) = $j_function_name(C(i,:)*1e-6,K);\n";
	}else{
		print DRIVER "\t\t\tJ_out(i) = $j_function_name(C(i,:)*1e-6,~,K);\n";
	}
	print DRIVER "\t\tend\n";
	if($l_old_output){
		print DRIVER "\telse\n";
		print DRIVER "\t\tJ_out = [];\n";
	}
	print DRIVER "\tend\n";
	print DRIVER "\tif ~no_fluxes\n";
	$input_for_flux_temp=$input_for_flux;
	$input_for_flux_temp =~ s/^c,/C(end,:),/;
	print DRIVER "\t\t[flux, ~, sources, ~, flux_out, ~] = $flux_function_name($input_for_flux_temp);\n";
	print DRIVER "\t\t\% Back to cubic centimeters\n";
	print DRIVER "\t\tflux = flux*1e-6;\n\t\tsources = sources*1e-6;\n\t\tflux_out = flux_out*1e-6;\n";
	if($l_old_output){
		print DRIVER "\telse\n";
		print DRIVER "\t\tflux = [];\n\t\tout = [];\n\t\toutflux_matrix = [];\n\t\tsources = [];\n\t\tform = [];\n";
	}
	print DRIVER "\tend\n\n";
	print DRIVER "\tC = C*1e-6;\n\tCf = Cf*1e-6;\n";
	
	
	print DRIVER "\t% printing out source file for constant monomer concentrations\n";
	print DRIVER "\tif consts_out && converged\n";
	print DRIVER "\t\tcfid = fopen(fname2c, 'w');\n";
	print DRIVER "\t\tfor i=monomers\n";
	print DRIVER "\t\t\tfprintf(cfid,'constant \%s \%f\\n',clust{i}, Cf(i));\n";
	print DRIVER "\t\tend\n";
	print DRIVER "\t\tfclose(cfid);\n";
	print DRIVER "\tend\n\n";
	print DRIVER "\tif ~no_fluxes\n";
	print DRIVER "\t\t% printing out source file for monomer sources\n";
	print DRIVER "\t\tif sources_out && converged\n";
	print DRIVER "\t\t\tsfid = fopen(fname2, 'w');\n";
	print DRIVER "\t\t\tfor i=monomers\n";
	print DRIVER "\t\t\t\tfprintf(sfid,'source \%s \%f\\n',clust{i}, sources(i));\n";
	print DRIVER "\t\t\tend\n";
	print DRIVER "\t\t\tfclose(sfid);\n";
	print DRIVER "\t\tend\n\n";
	
	$temp_val = max($n_max_type[0],$n_max_type_neg[0],$n_max_type_pos[0]);
	if(@cs_cluster){
		$temp_val = max($temp_val,$cs_cluster[0]);
	}
	if(@cs_cluster_neg){
		$temp_val = max($temp_val,$cs_cluster_neg[0]);
	}
	if(@cs_cluster_pos){
		$temp_val = max($temp_val,$cs_cluster_pos[0]);
	}
	$temp_val = $temp_val*2+1;
	print DRIVER "\t\tout = zeros($max_cluster_number,$temp_val,4);\n";
	print DRIVER "\t\toutflux_matrix = zeros($max_cluster_number,$max_cluster_number);\n";
	if($use_coag_to_outgoing && ($include_generic_neg || $include_generic_pos)){
		print DRIVER "\t\tfor i=[1:" . ($max_cluster-1) . "," . ($max_cluster+1) . ":$max_cluster_number]\n";
	}elsif($use_coag_to_outgoing){
		print DRIVER "\t\tfor i=1:" . ($max_cluster-1) . "\n";
	}else{
		print DRIVER "\t\tfor i=1:$max_cluster_number\n";
	}
	if($use_coag_to_outgoing && ($include_generic_neg || $include_generic_pos)){
		print DRIVER "\t\t\tfor j=[1:" . ($max_cluster-1) . "," . ($max_cluster+1) . ":$max_cluster_number]\n";
	}elsif($use_coag_to_outgoing){
		print DRIVER "\t\t\tfor j=1:" . ($max_cluster-1) . "\n";
	}else{
		print DRIVER "\t\t\tfor j=1:$max_cluster_number\n";
	}
	print DRIVER "\t\t\t\tk = " . $molecule_name[0] . "_in_clust(i)+" . $molecule_name[0] . "_in_clust(j)+1;\n";
	print DRIVER "\t\t\t\tout(i,k,:) = out(i,k,:)+flux_out(i,j,:);\n";
	print DRIVER "\t\t\t\toutflux_matrix(i,j) = sum(flux_out(i,j,:));\n";
	print DRIVER "\t\t\tend\n";
	print DRIVER "\t\tend\n";

	$temp_val = max($n_max_type[0],$n_max_type_neg[0],$n_max_type_pos[0]);
	if(@cs_cluster){
		$temp_val = max($temp_val,$cs_cluster[0]);
	}
	if(@cs_cluster_neg){
		$temp_val = max($temp_val,$cs_cluster_neg[0]);
	}
	if(@cs_cluster_pos){
		$temp_val = max($temp_val,$cs_cluster_pos[0]);
	}
	print DRIVER "\t\tform=zeros(1,$temp_val);\n";
	if($use_coag_to_outgoing && ($include_generic_neg || $include_generic_pos)){
		print DRIVER "\t\tfor i=[1:" . ($max_cluster-1) . "," . ($max_cluster+1) . ":$max_cluster_number]\n";
	}elsif($use_coag_to_outgoing){
		print DRIVER "\t\tfor i=1:" . ($max_cluster-1) . "\n";
	}else{
		print DRIVER "\t\tfor i=1:$max_cluster_number\n";
	}
	if($use_coag_to_outgoing && ($include_generic_neg || $include_generic_pos)){
		print DRIVER "\t\t\tfor j=[1:" . ($max_cluster-1) . "," . ($max_cluster+1) . ":$max_cluster_number]\n";
	}elsif($use_coag_to_outgoing){
		print DRIVER "\t\t\tfor j=1:" . ($max_cluster-1) . "\n";
	}else{
		print DRIVER "\t\t\tfor j=1:$max_cluster_number\n";
	}
	print DRIVER "\t\t\t\tif " . $molecule_name[0] . "_in_clust(i)<" . $molecule_name[0] . "_in_clust(j)\n";
	print DRIVER "\t\t\t\t\tform(" . $molecule_name[0] . "_in_clust(j))=form(" . $molecule_name[0] . "_in_clust(j))+flux(i,j)-flux(j,i);\n";
	print DRIVER "\t\t\t\tend\n";
	print DRIVER "\t\t\tend\n";
	print DRIVER "\t\tend\n";
	print DRIVER "\t\tform = form/2;\n\n";

	print DRIVER "\t\% Printing the flux matrix\n\t\tif(fl_out)\n";
	print DRIVER "\t\t\tfid=fopen(fname6,'w');\n\t\t\tfprintf(fid,'\%\% ');\n\t\t\tfor i=1:length(clust_flux)\n\t\t\t\tfprintf(fid,'\%s ',clust_flux{i});\n\t\t\tend\n";
	print DRIVER "\t\t\tfprintf(fid,'\\n');\n";
	if($keep_boundary_clusters){
		# Print also the "true" cluster names to distinguish the boundary clusters when the fluxes are tracked
		print DRIVER "\t\t\tfprintf(fid,'\%\% ');\n\t\t\tfor i=1:length(clust)\n\t\t\t\tfprintf(fid,'\%s ',clust{i});\n\t\t\tend\n";
		print DRIVER "\t\t\tfprintf(fid,'\\n');\n";
	}
	print DRIVER "\t\t\tfor i=1:size(flux,1)\n\t\t\t\tfor j=1:size(flux,2)\n\t\t\t\t\tfprintf(fid,'\%.6e ',flux(i,j));\n\t\t\t\tend\n\t\t\t\tfprintf(fid,'\\n');\n\t\t\tend\n\t\t\tfclose(fid);\n\t\tend\n\n";

	print DRIVER "\t\% Printing the outflux matrix\n\t\tif(outmat)\n";
	print DRIVER "\t\t\tfid=fopen(fname7,'w');\n\t\t\tfprintf(fid,'\%\% ');\n\t\t\tfor i=1:length(clust)\n\t\t\t\tfprintf(fid,'\%s ',clust{i});\n\t\t\tend\n\t\t\tfprintf(fid,'\\n');\n";
	print DRIVER "\t\t\tfor i=1:size(outflux_matrix,1)\n\t\t\t\tfor j=1:size(outflux_matrix,2)\n\t\t\t\t\tfprintf(fid,'\%.6e ',outflux_matrix(i,j));\n\t\t\t\tend\n\t\t\t\tfprintf(fid,'\\n');\n\t\t\tend\n\t\t\tfclose(fid);\n\t\tend\n\n";
	print DRIVER "\tend\n";
	

	if(! $l_save_outgoing){
		print DRIVER "\tC = C(:,1:$size_C_out);\n";
		print DRIVER "\tCf = Cf(1:$size_C_out);\n";
	}
	if(!$l_old_output){
		print DRIVER "\tif ~no_j && nargout>=6\n";
		print DRIVER "\t\tvarargout{1} = J_out;\n";
		print DRIVER "\tend\n";
		print DRIVER "\tif ~no_fluxes\n";
		print DRIVER "\t\tif nargout>=7\n";
		print DRIVER "\t\t\tvarargout{2} = flux;\n";
		print DRIVER "\t\t\tif nargout>=8\n";
		print DRIVER "\t\t\t\tvarargout{3} = out;\n";
		print DRIVER "\t\t\t\tif nargout>=9\n";
		print DRIVER "\t\t\t\t\tvarargout{4} = outflux_matrix;\n";
		print DRIVER "\t\t\t\t\tif nargout>=10\n";
		print DRIVER "\t\t\t\t\t\tvarargout{5} = sources;\n";
		print DRIVER "\t\t\t\t\t\tif nargout==11\n";
		print DRIVER "\t\t\t\t\t\t\tvarargout{6} = form;\n";
		print DRIVER "\t\t\t\t\t\tend\n";
		print DRIVER "\t\t\t\t\tend\n";
		print DRIVER "\t\t\t\tend\n";
		print DRIVER "\t\t\tend\n";
		print DRIVER "\t\tend\n";
		print DRIVER "\tend\n";
	}
	print DRIVER "end";


# Done with the driver routine
}

close(DRIVER);

################################################################
####                 END MAIN PROGRAM                       ####
################################################################


################################################################################
################################################################################
##   15. Subroutines                                                          ##
################################################################################
################################################################################

#######################################################################################################
# This subroutine combines two labels into one; if one of the clusters is positive and one is negative,
# it assumes that they recombine into a neutral molecule. If they are ions of the same sign, they can't
# be combined.
# The subroutine returns 3 values: 1) the label of the combined cluster, 2) a logical called valid_coll,
# which is 1 if the collision can happen, and 0 if it can't, which is the case when the colliding clusters
# are ions of the same sign, or we are disabling certain collisions, and 3) a logical called valid_evap,
# which is 1 if the cluster can evaporate back to the 2 original clusters, and 0 if it can't, which is
# true in the case where ions of opposite sign recombine or if cluster evaporations are forbidden.
#######################################################################################################

# combine_labels: ($value,$valid_coll,$valid_evap,$valid) = combine_labels($cluster1,$cluster2)
sub combine_labels(){
	my($cluster1,$cluster2,@nmols,$value,$valid_coll,$valid_evap);
	my($lcontinue,$lrec,$neg_clus,$pos_clus,$neutral_clus,@temp_array_neg,@temp_array_pos, @temp_array_neutral,@nmols1,@nmols2,$imol,$lout,$valid,$lfound,$i);
	$cluster1=$_[0];
	$cluster2=$_[1];

#	die "\nCluster numbers not defined! $cluster1 $cluster2\n\n" if ($cluster1 !~ /^\d+$/ || $cluster2 !~ /^\d+$/);

	$value='Undefined';
	$valid_coll=0;
	$valid_evap=0;
	$valid=0;

	$lcontinue=1;
	
	
	
	# First see if this is a collision between two charger ions -> that doesn't lead to a cluster
	if(($cluster1 > $max_cluster) &&  ($cluster2 > $max_cluster)){
		# the collision can still happen if the chargers have different polarities
		if($cluster1 != $cluster2){
			$valid_coll=1;
		}
		$lcontinue=0;
	
	# Then see if we are disabling (some) non-monomer collisions
	}elsif($disable_nonmonomers || $#nonmonomer_collisions_only>=0){
		# If so, check if this is a recombination in case non-monomer collisions are still allowed for those
		$lrec = 0;
		if($enable_rec_nonmonomers){
			if(($lnegative[$cluster1] && $lpositive[$cluster2]) || ($lpositive[$cluster1] && $lnegative[$cluster2])){
				$lrec = 1;
			}
		}
		# Now, check if this is a non-monomer process if this is is a collision type where they are disabled
		if(!$lrec){
			if(($cluster1 <= $max_cluster) &&  ($cluster2 <= $max_cluster)){
				$lmonomer1=$lmonomer[$cluster1];
				$lmonomer2=$lmonomer[$cluster2];
				if(! $lmonomer1 && ! $lmonomer2){
					$lcontinue=0;
					for($i=0;$i<=$#nonmonomer_collisions_only;$i++){
						if($i==$cluster1 || $i==$cluster2){
							$lcontinue=1;
							last;
						}
					}
				}
			}
		}
	}

	# Now some less trivial cases:
	if($lcontinue){

		# first the case where one of the clusters is negative and one is positive
		if(($lnegative[$cluster1] && $lpositive[$cluster2]) || ($lpositive[$cluster1] && $lnegative[$cluster2])){
			$valid_coll=1;
			# this can't evaporate back to the original pieces any more!
			$valid_evap=0;

			# which is which?
			if($lnegative[$cluster1]){
				$neg_clus=$cluster1;
				$pos_clus=$cluster2;
			}else{
				$neg_clus=$cluster2;
				$pos_clus=$cluster1;
			}

			# find the ions and change them into neutral molecules
			
			# first the negative ion
			@temp_array_neg=@{ $clust_comp[$neg_clus] };
			$wanted_mol=0;
			$wanted_mol++ until($wanted_mol > $nmol_types || ($temp_array_neg[$wanted_mol] > 0 && $charge[$wanted_mol] < 0));
			if($wanted_mol > $nmol_types){
				die "\nCan't find the charged molecule in the negative cluster $label[$neg_clus]!\n\n";
			}
			# is this ion described by a missing proton?
			if($missing_proton>0 && $temp_array_neg[$missing_proton] == 1){
				$temp_array_neg[$missing_proton] = 0;
			}else{
				# otherwise find the corresponding neutral molecule
				$corr_mol=$mol_num_from_name{$corr_neutral_mol[$wanted_mol]};
				# change this into a neutral cluster
				$temp_array_neg[$wanted_mol]-=1;
				$temp_array_neg[$corr_mol]+=1;
			}

			# then the positive ion
			@temp_array_pos=@{ $clust_comp[$pos_clus] };
			$temp_array_pos[$proton]-=1;

			# then combine the labels

			for($imol=0;$imol<$nmol_types;$imol++){
				$nmols[$imol]=$temp_array_neg[$imol]+$temp_array_pos[$imol];
				# If we get a negative number for some molecule type, something is wrong!
				if($nmols[$imol]<0){
					$value='Undefined';
					$valid_coll=0;
					$valid_evap=0;
					$lcontinue = 0;
					last;
				}
			}

		# if they both are ions of the same sign, they are not allowed to collide
		}elsif(($lnegative[$cluster1] && $lnegative[$cluster2]) || ($lpositive[$cluster1] && $lpositive[$cluster2])){

			$value='Undefined';
			$valid_coll=0;
			$valid_evap=0;
			$lcontinue = 0;

		# if the collision involves either two neutral clusters or a neutral cluster and an ion, 
		# the labels are combined by simply counting the sum of all the molecules
		}else{
			
			$valid_coll=1;
			$valid_evap=1;

			@nmols1=@{ $clust_comp[$cluster1] };
			@nmols2=@{ $clust_comp[$cluster2] };
			

			for($imol=0;$imol<$nmol_types;$imol++){
			
				$nmols[$imol]=$nmols1[$imol]+$nmols2[$imol];
				# If we try to charge something that can't be charged, we might get negative molecule numbers
				if($nmols[$imol]<0){
					$value='Undefined';
					$valid_coll=0;
					$valid_evap=0;
					$lcontinue = 0;
				}
			}

			# If there's a proton or missing proton, there needs to be a molecule that can accept or donate one
			if($proton>0 && $nmols[$proton]>0){
				$lfound = 0;
				for($imol=0;$imol<$nmol_types;$imol++){
					if($nmols[$imol] > 0){
						if($charge[$imol] == 0 && $corr_positive_mol[$imol] ne ''){
							$lfound=1;
						}
					}
				}
				if(!$lfound){
					$value='Undefined';
					$valid_coll=0;
					$valid_evap=0;
					$lcontinue = 0;
				}
			}elsif($missing_proton>0 && $nmols[$missing_proton]>0){
				$lfound = 0;
				for($imol=0;$imol<$nmol_types;$imol++){
					if($nmols[$imol] > 0){
						if($charge[$imol] == 0 && $corr_negative_mol[$imol] ne ''){
							$lfound=1;
						}
					}
				}
				if(!$lfound){
					$value='Undefined';
					$valid_coll=0;
					$valid_evap=0;
					$lcontinue = 0;
				}
			}
		}
	}
		
	# If the collision can happen, make a label for the product.
	if($lcontinue){
		$value=&create_label(\@nmols);

		# now, let's see if we are disabling collisions resulting in clusters outside of the system
		$valid=&check_validity($value);

		if(!$valid){
			$valid_evap=0;
			if($disable_excluded){
				$value='Undefined';
				$valid_coll=0;
				if($l_print_boundary){
					print "Not using boundary collision $label[$cluster1] + $label[$cluster2]\n";
				}
			}elsif($disable_excluded_ionization && (!$lneutral_cluster[$cluster1] || !$lneutral_cluster[$cluster2])){
				$value='Undefined';
				$valid_coll=0;
				if($l_print_boundary){
					print "Not using boundary collision $label[$cluster1] + $label[$cluster2]\n";
				}
			}
		}
	}
	
	if($valid_evap){
		# If one of the collision partners is a charger ion, the product can't evaporate back
		if($cluster1 > $max_cluster ||  $cluster2 > $max_cluster){
			$valid_evap = 0;
		}
	
		# Check if evaporations are disabled completely ...
		if($disable_evap){
			$valid_evap=0;
		# ... or just for this specific cluster ...
		}elsif($disable_evap_only[&get_cluster_number($value)]){
			$valid_evap=0;
		# ... or for fragmentation into two clusters ...
		}elsif($disable_nonmonomer_evaps){
			$lmonomer1=$lmonomer[$cluster1];
			$lmonomer2=$lmonomer[$cluster2];
			if(! $lmonomer1 && ! $lmonomer2){
				$valid_evap = 0;
			}
		}
	}


#	print"Combined $label[$cluster1] and $label[$cluster2] into $value, valid_evap is $valid_evap, valid_coll is $valid_coll\n";
	#print"Prevented collision of $label[$cluster1] and $label[$cluster2]\n" if(!$valid_coll);
	#print"Prevented evaporation of $value into $label[$cluster1] and $label[$cluster2]\n" if(!$valid_evap && $valid_coll);
	
	return ($value,$valid_coll,$valid_evap,$valid);
}


######################################################################################################################
## This subroutine checks if a cluster resulting from a collision or an ionization can get out of the system 
## (given how we have defined to let stuff get out), or if we will bring it back to the system boundary.
## It takes in the label of the cluster that is to be checked and returns the label of the cluster 
## after doing the check; if the cluster was allowed to get out, the label is the same as the input label,
## and if it is brougth back to the system, the label is the label for the cluster on the boundary.
## This also returns a logical $lout, which is 1 if the cluster got out, and 0 if it was brought back,
## and a hash array of how many molecules of each type were removed when bringing the cluster back to the boundary.
######################################################################################################################

# check_boundary: ($new_label,$lout,%nmonomers) = check_boundary($label0)
sub check_boundary(){
	my($label0,$valid,@temp_array,$n_first,@n_max,$lout,$new_label,$valid_new,%nmonomers);
	my($imol,$jmol,$element_name,@rules,$n_rules,$irule,$i);
	my($n_acid,$n_acid_other,$n_base,$n_base_other,$n_neutral,$l_remove_acids,$acid_min,$base_min,$lneg,$lpos);
	my($l_tried_acid_loop,$l_tried_base_loop,$diff,$n_neutral_old,$charge_clus);
	my($nmols_removable_max,$imol_max,$nmols_removable,$nmols_nonremovable);

	$label0=$_[0];
	$new_label=$label0;
	$lout=0;

	for($imol=0;$imol<$nmol_types;$imol++){
		$nmonomers{$molecule_name[$imol]}=0;
	}

	$valid=&check_validity($label0);
	if($valid){
		die "Cluster $label0 is inside the system, why is it being brought to the boundary?\n";
	}
	
	# see if we are allowing this to get out
	($ref,$comp_ok)=&determine_cluster_composition($label0);
	@temp_array=@$ref;
	# these are the criteria for neutrals getting out of the simulation
	@rules = @nucl_rules;
	$n_rules = $n_nucl_rules;
	@n_max=@n_max_type;
	$charge_clus = 0;
	#  but if this is a charged cluster, the criteria are different 
	for($imol=0;$imol<$nmol_types;$imol++){
		if($charge[$imol] == -1 && $temp_array[$imol] > 0){
			@rules = @nucl_rules_neg;
			$n_rules = $n_nucl_rules_neg;
			@n_max=@n_max_type_neg;
			$charge_clus -= $temp_array[$imol];
			last;   
		}elsif($charge[$imol] == 1 && $temp_array[$imol] > 0){
			@rules = @nucl_rules_pos;
			$n_rules = $n_nucl_rules_pos;
			@n_max=@n_max_type_pos;
			$charge_clus += $temp_array[$imol];
			last;
		}
	}

	# now we go through the rules to see if any of them are satisfied
	for($irule=0;$irule<$n_rules;$irule++){
		if($charge_clus<0){
			$lout = 2;
#			print "charge:$charge_clus,\t@n_max_type_neg\n;"
		}elsif($charge_clus>0){
			$lout = 3;
#			print "charge:$charge_clus,\t@n_max_type_pos\n;"
		}else{
			$lout = 1;
#			print "charge:$charge_clus,\t@n_max_type\n;"
		}
		for($imol=0;$imol<$nmol_types;$imol++){
			if($temp_array[$imol]<$rules[$irule][$imol]){
				$lout = 0;
			}
		}
		if($lout>0){
			last;
		}
	}
	
	# If the cluster didn't get out, we need to bring it back to the system
	if($lout==0){
		# First check if there are more molecules of some type than in ANY allowed clusters
		
		for($imol=0;$imol<$nmol_types;$imol++){
			if($charge_clus==0 && $temp_array[$imol] > $n_max_type[$imol]){
				$nmonomers{$molecule_name[$imol]} += $temp_array[$imol]-$n_max_type[$imol];
				$temp_array[$imol] = $n_max_type[$imol];
			}elsif($charge_clus<0 && $temp_array[$imol] > $n_max_type_neg[$imol]){
				$nmonomers{$molecule_name[$imol]} += $temp_array[$imol]-$n_max_type_neg[$imol];
				$temp_array[$imol] = $n_max_type_neg[$imol];
			}elsif($charge_clus>0 && $temp_array[$imol] > $n_max_type_pos[$imol]){
				$nmonomers{$molecule_name[$imol]} += $temp_array[$imol]-$n_max_type_pos[$imol];
				$temp_array[$imol] = $n_max_type_pos[$imol];
			}
		}
		$new_label = &create_label(\@temp_array);
		$valid=&check_validity($new_label);
	}
	
	# If this didn't help, try to bring it to the boundary one molecule at a time
	
	# First try using the acid and base strengths if they are given
	if($lout==0 && !$valid && $l_strength_neutral == 1){
		# First find out the number of acids and bases in the cluster
		$n_acid = 0;
		$n_base = 0;
		$n_acid_other = 0;
		$n_base_other = 0;
		$n_neutral = 0;
		for($imol=0;$imol<$nmol_types;$imol++){
			if($temp_array[$imol] > 0){
				# Number of neutral molecules in total
				if($charge[$imol]==0){
					$n_neutral += $temp_array[$imol];
				}
				if($acid_strength[$imol]>0){
					# Number of acids that can ...
					if($l_can_be_lost_bound[$imol] || $imol == $proton){
						$n_acid += $temp_array[$imol];
					}else{
					# ... and can't be lost
						$n_acid_other += $temp_array[$imol];
					}
					if($imol == $proton){
						# the proton forms an acid together with a molecule that would otherwise be a base
						$n_base_other -= 1;
						$n_neutral -= 1;
					}
				}
				if($base_strength[$imol]>0){
					# Number of bases that can ...
					if($l_can_be_lost_bound[$imol] || $imol == $missing_proton){
						$n_base += $temp_array[$imol];
					}else{
					# ... and can't be lost
						$n_base_other += $temp_array[$imol];
					}
					if($imol == $missing_proton){
						# the "missing proton" forms a base together with a molecule that would otherwise be an acid
						$n_acid_other -= 1;
						$n_neutral -= 1;
					}
				}
			}
		}
#		print "a:$n_acid, ao:$n_acid_other, b:$n_base, bo:$n_base_other\n";
		if($n_acid_other<0){
			$n_acid += $n_acid_other;
			$n_acid_other = 0;
		}
		if($n_base_other<0){
			$n_base += $n_base_other;
			$n_base_other = 0;
		}
#		print "a:$n_acid, ao:$n_acid_other, b:$n_base, bo:$n_base_other\n";
		if($n_acid+$n_acid_other > $n_base+$n_base_other){
			$l_remove_acids = 1;
		}else{
			$l_remove_acids = 0;
		}
		$acid_min = 0;
		$base_min = 0;
#		print "\n\n$label0\n";
		$n_neutral_old = $n_neutral;
		$l_tried_acid_loop = 0;
		$l_tried_base_loop = 0;
		$diff = 0;
		while(!$valid && $n_neutral>$diff){
#			print "$n_acid+$n_acid_other+$diff > $n_base+$n_base_other: ".($n_acid+$n_acid_other+$diff > $n_base+$n_base_other)."\n";
			# If we have more acids than bases, start removing acids ...
			if($l_remove_acids){
				$l_tried_acid_loop = 1;
				# ... if we have any neutral ones left
				if($acid_min>$acid_max){
					last;
				}
				ACID_MOL_LOOP: for($imol=0;$imol<$nmol_types;$imol++){
#					print "acid loop:$molecule_name[$imol]\n";
					if(!$l_can_be_lost_bound[$imol]){
						next;
					}
					# Start from the weakest acid in the cluster, but only if it is neutral
					if(($temp_array[$imol]>0) && ($acid_strength[$imol]==$acid_min) && $charge[$imol]==0){
						# Remove molecules one by one ...
						while($temp_array[$imol]>0){
#							print "removing $molecule_name[$imol]\n";
							$temp_array[$imol] -= 1;
							$n_neutral -= 1;
							$nmonomers{$molecule_name[$imol]} += 1;
							if($acid_strength[$imol]>0){
								$n_acid -= 1;
							}
							if($base_strength[$imol]>0){
								$n_base -= 1;
							}
							if($n_acid+$n_acid_other <= $n_base+$n_base_other+$diff){
								$l_remove_acids = 0;
							}
							$new_label = &create_label(\@temp_array);
#							print "->$new_label\n";
							$valid=&check_validity($new_label);
							# ... until we have a valid cluster or there are no longer more acids than bases left
							if($valid || !$l_remove_acids){
								last ACID_MOL_LOOP;
							}
						}
					}
				}
				# Move to the next stronger acid, if we are continuing to remove acids
				if($l_remove_acids){
					$acid_min += 1;
				}else{
					$base_min = 0;
				}
			# If we have more or as many bases than acids, start removing bases ...
			}else{
				$l_tried_base_loop = 1;
				# ... if we have any neutral ones left
				if($base_min>$base_max){
					last;
				}
				BASE_MOL_LOOP: for($imol=0;$imol<$nmol_types;$imol++){
#					print "base loop: $molecule_name[$imol]\n";
					if(!$l_can_be_lost_bound[$imol]){
						next;
					}
					# Start from the weakest base in the cluster, but only if it is neutral
					if(($temp_array[$imol]>0) && ($base_strength[$imol]==$base_min) && $charge[$imol]==0){
						# Remove molecules one by one ...
						while($temp_array[$imol]>0){
#							print "removing $molecule_name[$imol]\n";
							$temp_array[$imol] -= 1;
							$n_neutral -= 1;
							$nmonomers{$molecule_name[$imol]} += 1;
							if($base_strength[$imol]>0){
								$n_base -= 1;
							}
							if($acid_strength[$imol]>0){
								$n_acid -= 1;
							}
							if($n_acid+$n_acid_other+$diff > $n_base+$n_base_other){
								$l_remove_acids = 1;
							}
							$new_label = &create_label(\@temp_array);
#							print "->$new_label\n";
							$valid=&check_validity($new_label);
							# ... until we have a valid cluster or there are more acids than bases left
							if($valid || $l_remove_acids){
								last BASE_MOL_LOOP;
							}
						}
					}
				}
				if(!$l_remove_acids){
					$base_min += 1;
				}else{
					$acid_min = 0;
				}
			}
#			if($n_neutral<$n_neutral_old){
#				$n_neutral_old = $n_neutral;
#			}elsif(($l_remove_acids && $l_tried_acid_loop && $acid_min>$acid_max) || (!$l_remove_acids && $l_tried_base_loop && $base_min>$base_max)){
#				$diff += 1;
#				if($l_remove_acids){
#					$l_remove_acids = 0;
#				}else{
#					$l_remove_acids = 1;
#				}
#				$l_tried_acid_loop = 0;
#				$l_tried_base_loop = 0;
#			}
#			print "$valid, $n_neutral, $diff\n";
		}
	}
	
	# If no strengths are given, or the previous attempt didn't work, remove neutral molecules
	# according to the molecular content, starting from the compound with the highest molecule number
	if ($lout==0 && !$valid){
		$nmols_removable_max = 0;
		$imol_max = -1;
		while(!$valid){
			if($nmols_removable_max>0){
				#print "$new_label, removing $molecule_name[$imol_max]\n";
				$nmols_removable_max -= 1;
				$temp_array[$imol_max] -= 1;
				$nmonomers{$molecule_name[$imol_max]} += 1;
				$new_label = &create_label(\@temp_array);
				#print "->$new_label\n";
				$valid=&check_validity($new_label);
			}
			if(!$valid){
				# Find the most numerous removable neutral molecule
				for($imol=0;$imol<$nmol_types;$imol++){
					if(!$l_can_be_lost_bound[$imol] || $charge[$imol]!=0){
						next;
					}
					if($temp_array[$imol]>0){
						$nmols_removable = $temp_array[$imol];
						# If the cluster is charged and this molecule forms an ion with a proton or a missing proton,
						# subtract it from the number of removable neutral molecules of type $imol
						if($charge_clus<0 && $missing_proton >=0 && &compare_clusters($corr_negative_mol[$imol],"1$molecule_name[$imol]1$label_missing_proton")){
							if ($temp_array[$missing_proton]>0){
								$nmols_nonremovable = $temp_array[$missing_proton];
								$nmols_removable -= $nmols_nonremovable;
							}
							# Subtract other possible chargeable molecule types
							for($jmol=0;$jmol<$nmol_types;$jmol++){
								if($jmol==$imol || $charge[$imol]!=0){
									next;
								}
								if($temp_array[$jmol]>0 && &compare_clusters($corr_negative_mol[$jmol],"1$molecule_name[$jmol]1$label_missing_proton")){
									$nmols_removable += min($temp_array[$jmol],$nmols_nonremovable);
									$nmols_nonremovable -= min($temp_array[$jmol],$nmols_nonremovable);
								}
								if($nmols_nonremovable==0){
									last;
								}
							}
						}elsif($charge_clus>0 && $proton >=0 && &compare_clusters($corr_positive_mol[$imol],"1$molecule_name[$imol]1$label_proton")){
							if ($temp_array[$proton]>0){
								$nmols_nonremovable = $temp_array[$proton];
								$nmols_removable -= $nmols_nonremovable;
							}
							# Subtract other possible chargeable molecule types
							for($jmol=0;$jmol<$nmol_types;$jmol++){
								if($jmol==$imol || $charge[$imol]!=0){
									next;
								}
								if($temp_array[$jmol]>0 && &compare_clusters($corr_positive_mol[$jmol],"1$molecule_name[$jmol]1$label_proton")){
									$nmols_removable += min($temp_array[$jmol],$nmols_nonremovable);
									$nmols_nonremovable -= min($temp_array[$jmol],$nmols_nonremovable);
								}
								if($nmols_nonremovable==0){
									last;
								}
							}
						}
						if($nmols_removable > $nmols_removable_max){
							$nmols_removable_max = $nmols_removable;
							$imol_max = $imol;
						}
					}
				}
				if($nmols_removable_max==0){
					last;
				}
			}
		}
	}
	
	# If the cluster is still not back in the system, give up
	if ($lout==0 && !$valid){
		print "\nCan't bring back to boundary: $label[$iclus] + $label[$jclus] -> $label0 -> ?\n";
		die;
		$new_label = '';
	}
	return ($new_label,$lout,%nmonomers);


}


##############################################################################
## This subroutine finds the number of a cluster if it is in the system
## even if the molecules are in a different order than in the name label
##############################################################################

# get_cluster_number: $num = get_cluster_number($label1)
sub get_cluster_number(){
	my($label1,$label_new,$comp1,$comp_ok,$num);

	$label1=$_[0];
	if($use_coag_to_outgoing && $label1 eq 'cs_cluster'){
		$label1 = $cs_cluster_label;
	}elsif($use_coag_to_outgoing && $label1 eq 'cs_cluster_neg'){
		$label1 = $cs_cluster_neg_label;
	}elsif($use_coag_to_outgoing && $label1 eq 'cs_cluster_pos'){
		$label1 = $cs_cluster_pos_label;
	}

	$num='';
	# First check if we have this label in %num_from_label
	if(exists $num_from_label{"$label1"}){
		$num = $num_from_label{"$label1"};
	}else{
	# otherwise we can try to reorder the molecules in the label
		($comp1,$comp_ok)=&determine_cluster_composition($label1);
		if($comp_ok eq ''){
			$label_new=&create_label($comp1);
			if(exists $num_from_label{"$label_new"}){
				$num_from_label{"$label1"}=$num_from_label{"$label_new"};
				$num=$num_from_label{"$label1"};
			}
		}
	}
	
	return $num;
}


##############################################################################
## This subroutine finds the number of the corresponding dry cluster if this
## is a hydrate of a cluster in the system, the order of molecules is free
##############################################################################

# get_corresponding_dry: ($iclus,$nwaters) = get_corresponding_dry($label_in)
sub get_corresponding_dry(){
	my($iclus,$label_in,$imol,@temp_composition,$iwater,$nwaters,$temp_label);
	
	$label_in=$_[0];

	$nwaters=0;
	# Checking if we have the cluster in the system
	$iclus=&get_cluster_number($label_in);
	if($iclus eq '' && $lhydrate_average){
	# See if this is a hydrate of some cluster that we have
		@temp_composition=&find_composition($label_in);
		$iwater=-1;
		for($imol=0;$imol<=$#temp_composition;$imol++){
			if($temp_composition[$imol] =~ /^(\d+)$name_water$/){
				$iwater=$imol;
				$nwaters=$1;
			}
		}
		if($iwater>=0){
			# Remove the water and see if we have the dry cluster
			splice(@temp_composition,$iwater,1);
			$temp_label='';
			for($imol=0;$imol<=$#temp_composition;$imol++){
				$temp_label.=$temp_composition[$imol];
			}
			$iclus=&get_cluster_number($temp_label);
		}
	}
	
	return ($iclus,$nwaters)
}

##############################################################################
## This subroutine finds out if two clusters have the same composition,
## even if the molecules are in a different order in the name labels
##############################################################################

# compare_clusters: $lsame = compare_clusters($label1,$label2)

sub compare_clusters(){
	my($label1,$label2,$label_new,$comp1,$comp2,@composition1,@composition2,$lsame,$imol,$comp_ok);

	$label1=$_[0];
	if($use_coag_to_outgoing && $label1 eq 'cs_cluster'){
		$label1 = $cs_cluster_label;
	}elsif($use_coag_to_outgoing && $label1 eq 'cs_cluster_neg'){
		$label1 = $cs_cluster_neg_label;
	}elsif($use_coag_to_outgoing && $label1 eq 'cs_cluster_pos'){
		$label1 = $cs_cluster_pos_label;
	}
	$label2=$_[1];
	if($use_coag_to_outgoing && $label2 eq 'cs_cluster'){
		$label2 = $cs_cluster_label;
	}elsif($use_coag_to_outgoing && $label2 eq 'cs_cluster_neg'){
		$label2 = $cs_cluster_neg_label;
	}elsif($use_coag_to_outgoing && $label2 eq 'cs_cluster_pos'){
		$label2 = $cs_cluster_pos_label;
	}
	
	# First reordering the molecules in the labels if necessary
	if(!(exists $num_from_label{"$label1"})){
		($comp1,$comp_ok)=&determine_cluster_composition($label1);
		if($comp_ok eq ''){
			$label_new=&create_label($comp1);
			if(exists $num_from_label{"$label_new"}){
				$num_from_label{"$label1"}=$num_from_label{"$label_new"};
			}
		}
	}
	if(!(exists $num_from_label{"$label2"})){
		($comp2,$comp_ok)=&determine_cluster_composition($label2);
		if($comp_ok eq ''){
			$label_new=&create_label($comp2);
			if(exists $num_from_label{"$label_new"}){
				$num_from_label{"$label2"}=$num_from_label{"$label_new"};
			}
		}
	}
	
	# Now the labels point to the same cluster number if they have the same composition
	if((exists $num_from_label{"$label1"}) && (exists $num_from_label{"$label2"})){
		$lsame = ($num_from_label{"$label1"} ==  $num_from_label{"$label2"});
	}else{
		# This is just in case something funny happens - not sure if this is ever 
		@composition1=&find_composition($label1);
		@composition2=&find_composition($label2);
		$lsame=0;

		if($#composition1 == $#composition2){
			$lsame=1;
			for($imol=0;$imol<=$#composition1;$imol++){
				# see if all the elements from the 1st array belong to the 2nd
				$lsame=0 unless($composition1[$imol] ~~ @composition2);
			}
		}else{
			# these don't even have the same amount of molecule types
			$lsame=0; 
		}
			
	}
	return($lsame);
}

##############################################################################
## This subroutine finds out the composition of the cluster in the form e.g. 2A,1B,1D
## (i.e. both the amounts and names of the molecules)
##############################################################################

# find_composition: @composition = find_composition($ilabel)
sub find_composition(){
	my($ilabel,@composition,@nmols,@mol_name,$imol);

	$ilabel=$_[0];
	if($use_coag_to_outgoing && $ilabel eq 'cs_cluster'){
		$ilabel = $cs_cluster_label;
	}elsif($use_coag_to_outgoing && $ilabel eq 'cs_cluster_neg'){
		$ilabel = $cs_cluster_neg_label;
	}elsif($use_coag_to_outgoing && $ilabel eq 'cs_cluster_pos'){
		$ilabel = $cs_cluster_pos_label;
	}
	
	@nmols=split(/[a-zA-Z]+/,$ilabel);
	@mol_name=split(/-?\d+/,$ilabel);
	# do the shift because otherwise the first element will be empty
	shift(@mol_name);
	die "\nWhy doesn't the amount of molecule names match to the amount of molecule numbers in $label?\n\n" if ($#nmols != $#mol_name);
	for($imol=0;$imol<=$#nmols;$imol++){
		$composition[$imol]=$nmols[$imol] . $mol_name[$imol];
	}

	return(@composition);
}


###########################################################
## This subroutine checks to see if a given cluster is included
## in the clusters whose concentrations we are keeping track of;
## this returns 1 (true) if the cluster is valid, 0 (false) if it is not
###########################################################

# check_validity: $valid = check_validity($ilabel)
sub check_validity(){
	my ($ilabel,$iclus,$valid);

	$ilabel=$_[0];
	
	$valid=0;
	$iclus=&get_cluster_number($ilabel);
	if(($iclus ne '') && $iclus!=$n_flux_out && $iclus!=$n_flux_out_neg && $iclus!=$n_flux_out_pos){
		# if this is a valid cluster, it has a number
		# but if this is the cluster representing the coagulation sink, it's not really a valid cluster
		$valid=1;
	}

	return $valid;

}



####################################################################
## This subroutine takes in a label and sees if this is a monomer;
## 1 (true) is yes, 0 (false) is no
####################################################################

# check_monomer: $lmonomer = check_monomer($array_ref)
sub check_monomer(){
	my($total,$array_ref);

	$array_ref=$_[0];

	$total=&calculate_molecules($array_ref);

	if($total != 1){
		return 0;
	}else{
		return 1;
	}
}


####################################################################
## This subroutine takes in a label and calculates the total number
## of molecules in the cluster
####################################################################

# calculate_molecules: $total = calculate_molecules(\@temp_array)
sub calculate_molecules(){
	my(@temp_array,$total,$imol);

	@temp_array=@{ $_[0] };
	$total=0;
	for($imol=0;$imol<$nmol_types;$imol++){
		# Don't take the proton or the missing proton
		if($imol != $proton && $imol != $missing_proton){
			$total+=$temp_array[$imol];
		}
	}

	return $total;

}

####################################################################
## This calculates the mass of the cluster
####################################################################

# calculate_mass: $total_mass = calculate_mass(\@temp_array)
sub calculate_mass(){
	my(@temp_array,$total_mass,$imol);

	@temp_array=@{ $_[0] };
	$total_mass=0;
	for($imol=0;$imol<$nmol_types;$imol++){
		if($temp_array[$imol] != 0){
			$total_mass+=$temp_array[$imol]*$mass[$imol];
		}
	}

	return $total_mass;

}

####################################################################
## This calculates the radius from the total volume of the cluster
####################################################################

# calculate_radius: $radius = calculate_radius($array_ref)
sub calculate_radius(){
	my($total_volume,$radius,$array_ref);

	$array_ref=$_[0];

	($total_volume,$radius)=&calculate_volume_and_radius($array_ref);

	return $radius;

}

####################################################################
## This calculates the volume and radius of the cluster, assuming
## that it's just the sum of all the molecular volumes (which is
## assumed to be the mass divided by the density)
####################################################################

# calculate_volume_and_radius: ($total_volume, $radius) = calculate_volume_and_radius(\@temp_array)
sub calculate_volume_and_radius(){
	my(@temp_array,$total_volume,$volume,$imol);

	@temp_array=@{ $_[0] };
	$total_volume=0.0;
	for($imol=0;$imol<$nmol_types;$imol++){
		# Don't take the proton or the missing proton
		if($imol != $proton && $imol != $missing_proton){
			if($temp_array[$imol] != 0){
				$volume=$mass[$imol]/$density[$imol];
				$total_volume+=$temp_array[$imol]*$volume;
			}
		}
	}

	$radius=(3*$total_volume/(4*$pi))**(1.0/3.0);

	return ($total_volume, $radius);

}


####################################################################
## This subroutine takes the number of molecules in the cluster and
## creates a label for it
####################################################################

# create_label: $value = create_label(\@molecules_in_clust)
sub create_label(){
	my($value,$imol);
	my @molecules_in_clust=@{ $_[0] };

	# Create the label
	$value='';
	for($imol=0;$imol<$nmol_types;$imol++){
		if($molecules_in_clust[$imol] > 0){
			$value = $value . $molecules_in_clust[$imol] . $molecule_name[$imol];
		}
	}

	return $value;
}

####################################################################
## This subroutine takes in a label and figures out how many of each
## molecule type are present
####################################################################

# determine_cluster_composition: (\@molecules,$message) = determine_cluster_composition($ilabel)
sub determine_cluster_composition(){
	my($ilabel,@molecules,$imol,@nmols,@mol_name,$message);

	$ilabel=$_[0];
	if($use_coag_to_outgoing && $ilabel eq 'cs_cluster'){
		$ilabel = $cs_cluster_label;
	}elsif($use_coag_to_outgoing && $ilabel eq 'cs_cluster_neg'){
		$ilabel = $cs_cluster_neg_label;
	}elsif($use_coag_to_outgoing && $ilabel eq 'cs_cluster_pos'){
		$ilabel = $cs_cluster_pos_label;
	}
	
	# denote the charger ions as proton and removed proton, so they can easily be added and substracted from clusters
	if($ilabel eq $label_generic_neg){
		$ilabel = "1$label_missing_proton";
	}elsif($ilabel eq $label_generic_pos){
		$ilabel = "1$label_proton";
	}
	
	
	@nmols=split(/[a-zA-Z]+/,$ilabel);
	@mol_name=split(/-?\d+/,$ilabel);
	# shift the indices because otherwise the first element will be empty
	shift(@mol_name);
	# Going through the allowed molecule types to make a zero array
	for($imol=0;$imol<$nmol_types;$imol++){
		$molecules[$imol] = 0;
	}
	$message='';
	
	# Going through the molecule types in the cluster name, and if they are valid ones filling the array accordingly
	for($imol=0;$imol<=$#nmols;$imol++){
		if(exists $mol_num_from_name{$mol_name[$imol]}){
			$molecules[$mol_num_from_name{$mol_name[$imol]}]=$nmols[$imol];
		}else{
			$message="+ other molecule types";
		}
	}
	return (\@molecules,$message);
}


####################################################################
## This subroutine takes in a scalar and determines if it is a 
## positive number, possibly padded with spaces
####################################################################

# is_pos_number: $l_num = is_pos_number($a)
sub is_pos_number(){
	my($a,$l_num);
	
	$a=$_[0];
	$l_num = 0;
	
	if($a =~ m/^\s*\d+\.?\d*\s*$/){
		# 1 or 1. or 1.2
		$l_num = 1;
	}elsif($a =~ m/^\s*\.\d+\s*$/){
		# .1
		$l_num = 1;
	}elsif($a =~ m/^\s*\d+\.?\d*e[+-]?\d+\s*$/i){
		# 1e2, 1.e2, 1.2e3, 1.2e+3, 1.2e-3, 1.2E-3
		$l_num = 1;
	}elsif($a =~ m/^\s*\.\d+e[+-]?\d+\s*$/i){
		# .1e2, .1e+2, .1E-2
		$l_num = 1;
	}
	
	return $l_num;
}


####################################################################
####################################################################
## 
## List of subroutines
## 
####################################################################
####################################################################



# combine_labels: ($value,$valid_coll,$valid_evap,$valid) = combine_labels($cluster1,$cluster2)
# check_boundary: ($new_label,$lout,%nmonomers) = check_boundary($label0)

# get_cluster_number: $num = get_cluster_number($label1)
# get_corresponding_dry: ($iclus,$nwaters) = get_corresponding_dry($label_in)
# compare_clusters: $lsame = compare_clusters($label1,$label2)
# find_composition: @composition = find_composition($ilabel)

# check_validity: $valid = check_validity($ilabel)
# check_monomer: $lmonomer = check_monomer($array_ref)
# calculate_molecules: $total = calculate_molecules(\@temp_array)

# calculate_mass: $total_mass = calculate_mass(\@temp_array)
# calculate_radius: $radius = calculate_radius($array_ref)
# calculate_volume_and_radius: ($total_volume, $radius) = calculate_volume_and_radius(\@temp_array)

# create_label: $value = create_label(\@molecules_in_clust)
# determine_cluster_composition: (\@molecules,$message) = determine_cluster_composition($ilabel)

# is_pos_number: $l_num = is_pos_number($a)

