def infer_bonds(atomic_numbers, positions, scale_factor=1.2):
  """Infer bonds based on interatomic distances and covalent radii."""
  from numpy import linalg
  from ase.data import covalent_radii
  bonds = []
  num_atoms = len(atomic_numbers)
  # Loop over all atom pairs
  for i in range(num_atoms):
      for j in range(i + 1, num_atoms):
          dist = linalg.norm(positions[i] - positions[j])  # Interatomic distance
          r_i = covalent_radii[atomic_numbers[i]]  # Covalent radius of atom i
          r_j = covalent_radii[atomic_numbers[j]]  # Covalent radius of atom j
          # If dist. is less than the sum of the covalent radii (with a scale factor)
          if dist < scale_factor * (r_i + r_j):
              bonds.append((i, j))  # Bond exists between atoms i and j
  return bonds

def read_xyz(file_i_XYZ):
  from ase.io import read
  try:
    out = read(file_i_XYZ)
  except:
    out = float("nan")
  return out 

def identify_1(out):
  from numpy import isnan
  if type(out) == type(float(0)):
    if isnan(out):
      return out
  
  try:
    from rdkit import Chem
    rdkit_mol = Chem.RWMol()
    atomic_numbers = out.get_atomic_numbers()
    positions = out.get_positions()[atomic_numbers != 1]
    atomic_numbers = atomic_numbers[atomic_numbers != 1]
    atom_map = {}  # To store ASE atom index to RDKit atom index mapping
    for i, atomic_num in enumerate(atomic_numbers):
      idx = rdkit_mol.AddAtom(Chem.Atom(int(atomic_num)))
      atom_map[i] = idx
    # Add bonds to the RDKit molecule
    bonds = infer_bonds(atomic_numbers, positions)
    for i, j in bonds:
      rdkit_mol.AddBond(atom_map[i], atom_map[j], Chem.BondType.SINGLE)
    rdkit_mol = rdkit_mol.GetMol()
    return Chem.MolToSmiles(rdkit_mol)
  except:
    return float("nan")

def identify_2(out):
  from numpy import isnan
  if type(out) == type(float(0)):
    if isnan(out):
      return out

  try:
    from rdkit import Chem
    rdkit_mol = Chem.RWMol()
    atomic_numbers = out.get_atomic_numbers()
    positions = out.get_positions()
    atom_map = {}  # To store ASE atom index to RDKit atom index mapping
    for i, atomic_num in enumerate(atomic_numbers):
      idx = rdkit_mol.AddAtom(Chem.Atom(int(atomic_num)))
      atom_map[i] = idx
    # Add bonds to the RDKit molecule
    bonds = infer_bonds(atomic_numbers, positions)
    for i, j in bonds:
      rdkit_mol.AddBond(atom_map[i], atom_map[j], Chem.BondType.SINGLE)
    rdkit_mol = rdkit_mol.GetMol()
   
    #TODO: WE NEED TO BE ABLE LOCATE RADICALS FIRST!! 
    #for a in rdkit_mol.GetAtoms():
    #  if a.GetNumImplicitHs():
    #    a.SetNumRadicalElectrons(a.GetNumImplicitHs())
    #    a.SetNoImplicit(True)
    #    a.UpdatePropertyCache()
    #rdkit_mol = Chem.RemoveAllHs(mol)
    return Chem.MolToSmiles(rdkit_mol,kekuleSmiles=True)
  except:
    return float("nan")

def identify_3(out):
  from numpy import isnan
  if isnan(out):
    return out
  
  #TODO somebody could try to do the 2nd level identification 
  # Step 5: Finalize the RDKit molecule
  #from rdkit.Chem import rdDetermineBonds
  #rdDetermineBonds.DetermineBonds(rdkit_mol,charge=0)
  #rdDetermineBonds.DetermineBonds(rdkit_mol,charge=0)
  #rdkit_mol = rdkit_mol.GetMol()  # Finalize the molecule
  #print(Chem.MolToSmiles(rdkit_mol))
  #rdkit_mol = Chem.RemoveHs(rdkit_mol)
  #print(Chem.MolToSmiles(rdkit_mol))
  #print(Chem.MolToSmiles(rdkit_mol,allHsExplicit=True))
  #print(Chem.MolToMolBlock(rdkit_mol)) 
  #print(out)
  #print(Chem.MolToSmiles(rdkit_mol,isomericSmiles=False))
  #print(Chem.MolToSmiles(rdkit_mol,kekuleSmiles=True))
  #params = AllChem.ETKDGv3()
  ##params.randomSeed = "0xf00d" # optional random seed for reproducibility
  ##AllChem.EmbedMolecule(rdkit_mol,useBasicKnowledge=False)
  ##AllChem.EmbedMolecule(rdkit_mol)
  #print(Chem.MolToMolBlock(rdkit_mol)) 
  #print(Chem.MolToSmiles(rdkit_mol))
  #print(out)
  #print(Chem.MolToSmiles(rdkit_mol,isomericSmiles=False))a
  return float("nan")
