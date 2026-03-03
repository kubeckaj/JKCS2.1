import re

def log2vib(molecule):
    with open(molecule.log_file_path, 'r') as file:
        content = file.read()
        if molecule.program.lower() == "g16":
            vibrations = re.findall(r"Frequencies --\s+(-?\d+\.\d+)", content)
        elif molecule.program.lower() == "orca":
            vibrations = []
            vib = re.search(r'[-+]?\d*\.\d+\s*cm\*\*-1', content)
            if vib:
                vibration = float(vib.group().split()[0])
                vibrations.append(vibration)
        else:
            return 'No vibrations found'
    return vibrations


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
        print("Error in calculating the eckart tunneling. Returning tunneling coefficient 1")
        return 1


def rate_constant(TS_conformers, reactant_conformers, product_conformers, T=298.15):
    from numpy import exp, sum, float64, NaN
    k_b = 1.380649e-23
    h = 6.62607015e-34
    HtoJ = 43.597447222e-19
    Htokcalmol = 627.509
    Na = 6.022e23
    liters_to_cm3 = 1000
    mol_per_liter = 1/(0.0821*T)
    p_ref = (mol_per_liter*Na)/liters_to_cm3
    kappa = 1

    if reactant_conformers and TS_conformers:
        reactant_molecules = [mol for mol in reactant_conformers if 'OH' not in mol.name]
        OH = next((mol for mol in reactant_conformers if 'OH' in mol.name))

        lowest_reactant = min(reactant_molecules, key=lambda molecule: molecule.zero_point_corrected)
        lowest_TS = min(TS_conformers, key=lambda molecule: molecule.zero_point_corrected)

        lowest_ZP_TS_J = lowest_TS.zero_point_corrected * HtoJ
        lowest_ZP_reactant_J = lowest_reactant.zero_point_corrected * HtoJ
        lowest_ZP_TS_kcalmol = lowest_TS.zero_point_corrected * Htokcalmol
        lowest_ZP_reactant_kcalmol = (lowest_reactant.zero_point_corrected + OH.zero_point_corrected) * Htokcalmol

        if OH:
            sum_reactant_ZP_J = lowest_ZP_reactant_J + (OH.zero_point_corrected * HtoJ)
            Q_reactant = sum([exp(-(lowest_ZP_reactant_J - (mol.zero_point_corrected * HtoJ)) / (k_b * T)) * mol.Q for mol in reactant_molecules]) * OH.Q
        else:
            sum_reactant_ZP_J = lowest_ZP_reactant_J
            Q_reactant = sum([exp(-(lowest_ZP_reactant_J - (mol.zero_point_corrected * HtoJ)) / (k_b * T)) * mol.Q for mol in reactant_molecules])

        Q_TS = sum([exp(-(mol.zero_point_corrected - lowest_TS.zero_point_corrected) * HtoJ / (k_b * T)) * mol.Q for mol in TS_conformers])
        Q_reactant = sum([exp(-((mol.zero_point_corrected - lowest_reactant.zero_point_corrected) * HtoJ) / (k_b * T)) * mol.Q for mol in reactant_molecules]) * OH.Q

        if product_conformers:
            product_molecules = [mol for mol in product_conformers if mol.name not in ('H2O', 'H2O_DLPNO')]
            H2O = next((mol for mol in product_conformers if mol.name in ("H2O", "H2O_DLPNO")), None)

            if product_molecules and H2O:
                lowest_product = min(product_molecules, key=lambda molecule: molecule.electronic_energy)
                lowest_EE_product_kcalmol = (lowest_product.electronic_energy + H2O.electronic_energy) * Htokcalmol
                lowest_EE_TS_kcalmol = lowest_TS.electronic_energy * Htokcalmol
                lowest_EE_reactant_kcalmol = (lowest_reactant.electronic_energy + OH.electronic_energy) * Htokcalmol

                imag = abs(lowest_TS.vibrational_frequencies[0])
                kappa = eckart(lowest_EE_TS_kcalmol, lowest_EE_reactant_kcalmol, lowest_EE_product_kcalmol, imag)
                k = kappa * (k_b*T)/(h*p_ref) * (Q_TS/Q_reactant) * exp(-(lowest_ZP_TS_J - sum_reactant_ZP_J) / (k_b * T))
            else:
                print("No H2O in product molecules. No tunneling correction will be calculated")
                k = kappa * (k_b*T)/(h*p_ref) * (Q_TS/Q_reactant) * exp(-(lowest_ZP_TS_J - sum_reactant_ZP_J) / (k_b * T))
        else:
            k = kappa * (k_b*T)/(h*p_ref) * (Q_TS/Q_reactant) * exp(-(lowest_ZP_TS_J - sum_reactant_ZP_J) / (k_b * T))

        print(f"Ea: {lowest_ZP_TS_kcalmol - lowest_ZP_reactant_kcalmol:.4f} kcal/mol  Q_TS: {Q_TS}  k: {k}")
        return k, kappa

    from numpy import NaN
    return NaN, NaN
