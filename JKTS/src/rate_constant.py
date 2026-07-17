from output import logger
import os
import runtime
from results import record_rate, format_rate
from metadata import load_metadata
from classes import Molecule


def eckart(SP_TS, SP_reactant, SP_product, imag, T=[298.15]):
    from numpy import pi, sqrt, arange, exp, max, cosh

    def Gcalc(Emin, GSize, V, A, B, L, h, kB, m, T, v1):
        E = arange(Emin, 2*max(V), GSize)
        K = [0 for x in range(len(E))]
        EK = [0 for x in range(len(E))]
        EKa = [0 for x in range(len(E))]
        for i in range(len(E)):
            C = (h**2)/(8*m*(L**2))
            a = 0.5*sqrt(E[i]/C)
            b = 0.5*sqrt((E[i]-A)/C)
            d = 0.5*sqrt((B-C)/C)
            K[i] = 1-((cosh(2*pi*(a-b))+cosh(2*pi*d))/(cosh(2*pi*(a+b))+cosh(2*pi*d)))
            EK[i] = E[i]
            EKa[i] = E[i]*627.509

        G = [0 for x in range(len(T))]
        for q in range(len(T)):
            GK = [0 for x in range(len(K))]
            for j in range(len(K)):
                GK[j] = K[j]*exp(-EK[j]/(kB*T[q]))
            GI = [0 for x in range(len(K))]
            for l in range(len(EK)-1):
                GI[l] = (0.5*(GK[l]+GK[l+1])*abs(EK[l]-EK[l+1]))

            GI = sum(GI)
            GI = GI*(exp(v1/(kB*T[q]))/(kB*T[q]))
            G[q] = GI+exp(v1/(kB*T[q]))*exp(-EK[len(EK)-1]/(kB*T[q]))

        return G, EKa, K, GK

    try:
        c = 2.99792458e+8
        kB = 3.1668152e-6
        h = 2*pi
        Na = 6.0221409e+23

        E1 = SP_TS - SP_reactant
        E2 = SP_TS - SP_product
        mu = 1
        v1 = ((E1*4184)/Na)/4.3597447222071e-18
        v2 = ((E2*4184)/Na)/4.3597447222071e-18
        wau = (imag*100)*c*2.418884326509e-17
        m = mu*1822.888479

        F = -4*(pi**2)*(wau**2)*m
        F2 = -4*(pi**2)*(wau**2)*1
        A = v1-v2
        B = (sqrt(v2)+sqrt(v1))**2
        L = -pi*(A-B)*(B+A)/(sqrt(-2*F*B)*B)

        x = arange(-3, 3, 0.01)
        x = x/(sqrt(mu))

        y = [0 for i in range(len(x))]
        V = [0 for i in range(len(x))]
        xa = [0 for i in range(len(x))]
        Va = [0 for i in range(len(x))]

        for i in range(len(x)):
            y[i] = -exp((2*pi*x[i])/L)
            V[i] = ((-(y[i]*A)/(1-y[i])) - ((y[i]*B)/((1-y[i])**2)))
            xa[i] = 0.529177*x[i]*sqrt(mu)
            Va[i] = V[i]*627.509

        VB = [0, 0]
        VB[0] = V[0]
        VB[1] = V[len(x)-1]
        Emin = max(VB)
        Gdiff = 1
        GSize = max(V)/50
        [Gold, EKa, K, GK] = Gcalc(Emin, GSize, V, A, B, L, h, kB, m, T, v1)
        GSize = GSize/10
        runs = 0
        while Gdiff >= 0.001:
            [Gnew, EKa, K, GK] = Gcalc(Emin, GSize, V, A, B, L, h, kB, m, T, v1)
            GSize = GSize/10
            Gdiffcalc = [0 for x in range(len(T))]
            for j in range(len(T)):
                Gdiffcalc[j] = abs(Gnew[j]-Gold[j])/Gold[j]
            Gdiff = max(Gdiffcalc)
            Gold = Gnew
            runs = runs+1

        [G, EKa, K, GK] = Gcalc(Emin, GSize, V, A, B, L, h, kB, m, T, v1)

        kappa = G[0]
        return kappa
    except Exception as e:
        logger.warning(f"Error calculating the Eckart tunneling ({e}). Returning tunneling coefficient 1")
        return 1


def rate_constant(TS_conformers, reactant_conformers, product_conformers, T=298.15, symmetry=1):
    from numpy import exp, sum
    from results import RateResult
    k_b = 1.380649e-23
    h = 6.62607015e-34
    HtoJ = 43.597447222e-19
    Htokcalmol = 627.509
    P_ref = 101325  # Pa (1 atm); must match P in Molecule.partition_function() and G16's log Q
    p_ref = P_ref / (k_b * T) / 1e6  # standard-state number density, molecules cm^-3
    kappa = 1

    if reactant_conformers and TS_conformers:
        # Identify the abstracting radical (OH or Cl) from the reactant list
        _radical_names = ('OH', 'OH_DLPNO', 'Cl', 'Cl_DLPNO', 'NO3', 'NO3_DLPNO')
        radical = next((mol for mol in reactant_conformers if mol.name in _radical_names), None)
        reactant_molecules = [mol for mol in reactant_conformers if mol.name not in _radical_names]

        # Drop molecules whose energies never parsed (would crash the min() below)
        dropped = [m.name for m in TS_conformers + reactant_molecules if m.zero_point_corrected is None]
        if radical is not None and radical.zero_point_corrected is None:
            dropped.append(radical.name)
        if dropped:
            logger.warning(f"Skipping molecules with missing ZPE/energy in rate step: {dropped}")
        TS_conformers = [m for m in TS_conformers if m.zero_point_corrected is not None]
        reactant_molecules = [m for m in reactant_molecules if m.zero_point_corrected is not None]
        if not TS_conformers or not reactant_molecules or (radical is not None and radical.zero_point_corrected is None):
            logger.error("Cannot compute rate constant - required TS/reactant/radical energies are missing.")
            return RateResult(sigma=symmetry, T=T)

        # Partition functions are stored at the 298.15 K parse-time default; recompute from
        # stored spectroscopic data if the rate is requested at a different temperature.
        if abs(T - 298.15) > 1e-9:
            to_update = list(TS_conformers) + list(reactant_molecules)
            if radical is not None:
                to_update.append(radical)
            for mol in to_update:
                if mol.vibrational_frequencies and getattr(mol, 'rot_temps', None):
                    try:
                        mol.partition_function(T)
                    except Exception as e:
                        logger.warning(f"Could not recompute Q for {mol.name} at {T} K ({e}); using 298.15 K value")

        lowest_reactant = min(reactant_molecules, key=lambda molecule: molecule.zero_point_corrected)
        lowest_TS = min(TS_conformers, key=lambda molecule: molecule.zero_point_corrected)

        lowest_ZP_TS_J = lowest_TS.zero_point_corrected * HtoJ
        lowest_ZP_reactant_J = lowest_reactant.zero_point_corrected * HtoJ
        lowest_ZP_TS_kcalmol = lowest_TS.zero_point_corrected * Htokcalmol

        if radical:
            sum_reactant_ZP_J = lowest_ZP_reactant_J + (radical.zero_point_corrected * HtoJ)
            lowest_ZP_reactant_kcalmol = (lowest_reactant.zero_point_corrected + radical.zero_point_corrected) * Htokcalmol
            Q_reactant = sum([exp(-((mol.zero_point_corrected - lowest_reactant.zero_point_corrected) * HtoJ) / (k_b * T)) * mol.Q for mol in reactant_molecules]) * radical.Q
        else:
            sum_reactant_ZP_J = lowest_ZP_reactant_J
            lowest_ZP_reactant_kcalmol = lowest_reactant.zero_point_corrected * Htokcalmol
            Q_reactant = sum([exp(-((mol.zero_point_corrected - lowest_reactant.zero_point_corrected) * HtoJ) / (k_b * T)) * mol.Q for mol in reactant_molecules])

        Q_TS = sum([exp(-(mol.zero_point_corrected - lowest_TS.zero_point_corrected) * HtoJ / (k_b * T)) * mol.Q for mol in TS_conformers])

        imag = None
        if product_conformers:
            # Identify the small product molecule (H2O or HCl)
            _product_names = ('H2O', 'H2O_DLPNO', 'HCl', 'HCl_DLPNO', 'HNO3', 'HNO3_DLPNO')
            product_small = next((mol for mol in product_conformers if mol.name in _product_names), None)
            product_molecules = [mol for mol in product_conformers if mol.name not in _product_names]

            # Incomplete products -> fall back to no tunneling rather than crash
            product_molecules = [mol for mol in product_molecules if mol.electronic_energy is not None]
            if product_small is not None and product_small.electronic_energy is None:
                product_small = None

            if product_molecules and product_small and radical:
                lowest_product = min(product_molecules, key=lambda molecule: molecule.electronic_energy)
                lowest_EE_product_kcalmol = (lowest_product.electronic_energy + product_small.electronic_energy) * Htokcalmol
                lowest_EE_TS_kcalmol = lowest_TS.electronic_energy * Htokcalmol
                lowest_EE_reactant_kcalmol = (lowest_reactant.electronic_energy + radical.electronic_energy) * Htokcalmol

                imag = abs(lowest_TS.vibrational_frequencies[0])
                kappa = eckart(lowest_EE_TS_kcalmol, lowest_EE_reactant_kcalmol, lowest_EE_product_kcalmol, imag, T=[T])
                k = kappa * (k_b*T)/(h*p_ref) * (Q_TS/Q_reactant) * exp(-(lowest_ZP_TS_J - sum_reactant_ZP_J) / (k_b * T))
            else:
                logger.warning("No small product molecule (H2O/HCl/HNO3) found. No tunneling correction will be calculated")
                k = kappa * (k_b*T)/(h*p_ref) * (Q_TS/Q_reactant) * exp(-(lowest_ZP_TS_J - sum_reactant_ZP_J) / (k_b * T))
        else:
            k = kappa * (k_b*T)/(h*p_ref) * (Q_TS/Q_reactant) * exp(-(lowest_ZP_TS_J - sum_reactant_ZP_J) / (k_b * T))

        k *= symmetry  # σᵢ reaction-path degeneracy (k only, not κ)

        Ea = lowest_ZP_TS_kcalmol - lowest_ZP_reactant_kcalmol
        return RateResult(k=k, kappa=kappa, Ea=Ea, Q_TS=Q_TS, Q_reactant=Q_reactant, sigma=symmetry,
                          n_ts=len(TS_conformers), n_reactant=len(reactant_molecules), imag=imag, T=T)

    return RateResult(sigma=symmetry, T=T)


def assemble_and_record_rate(TS_molecules):
    # Locate the finished reactant/product pickles for this TS channel directory,
    channel_name = os.path.basename(runtime.start_dir)
    symmetry = read_reaction_path_degeneracy(runtime.start_dir)
    reactant_pkl_name = channel_name.split("_")[0]
    for rates_dir in (os.path.dirname(runtime.start_dir), runtime.start_dir):
        reactant_pkl_path = os.path.join(rates_dir, f'reactants/Final_reactants_{reactant_pkl_name}.pkl')
        product_pkl_path = os.path.join(rates_dir, f'products/Final_products_{channel_name}.pkl')
        logger.info(f"Looking for reactant and product pickles: {reactant_pkl_path}, {product_pkl_path}")
        if os.path.exists(reactant_pkl_path):
            final_reactants = Molecule.load_molecules_from_pickle(reactant_pkl_path)
            final_products = Molecule.load_molecules_from_pickle(product_pkl_path) if os.path.exists(product_pkl_path) else []
            result = rate_constant(TS_molecules, final_reactants, final_products, T=runtime.args.T, symmetry=symmetry)
            record_rate(rates_dir, channel_name, result, method=runtime.args.method)
            logger.results(format_rate(result))
            return True
    logger.error(f"Could not find Final_reactants_{reactant_pkl_name}.pkl in ../reactants/ or ./reactants/")
    return False


def read_reaction_path_degeneracy(directory):
    # σᵢ recorded per TS directory in .metadata by mkdir(); defaults to 1
    try:
        return int(load_metadata(directory).get('reaction_path_degeneracy', 1))
    except (TypeError, ValueError):
        return 1