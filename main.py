import networkx as nx

def read_mol_file(filename):
    """
    Reads an OEChem mol file in V2000 format and returns the atom and bond
    information as separate lists.

    Args:
        filename (str): the name of the mol file

    Returns:
        atoms (list): a list of dictionaries representing the atoms in the molecule,
            with keys for element symbol, charge, and coordinates
        bonds (list): a list of tuples representing the bonds in the molecule,
            where each tuple contains the indices of the two atoms involved in the bond
            and the bond order
    """
    with open(filename, 'r') as f:
        lines = f.readlines()
    num_atoms = int(lines[3].split()[0])
    num_bonds = int(lines[3].split()[1])
    atoms = []
    bonds = []
    for i in range(num_atoms):
        atom_line = lines[i+4].split()
        coords = tuple(map(float, [atom_line[0], atom_line[1], atom_line[2]]))
        atom = {'element': atom_line[3], 'charge': int(atom_line[4]), 'coords': coords}
        atoms.append(atom)
    for i in range(num_bonds):
        bond_line = lines[i+num_atoms+4].split()
        atom1 = int(bond_line[0])
        atom2 = int(bond_line[1])
        bond_order = int(bond_line[2])
        bonds.append((atom1, atom2, bond_order))
    return atoms, bonds
  
  def generate_molecular_graph(atoms, bonds):
    """
    Generates a molecular graph from a list of atoms and bonds.

    Args:
        atoms (list): a list of dictionaries representing the atoms in the molecule,
            with keys for element symbol, charge, and coordinates
        bonds (list): a list of tuples representing the bonds in the molecule,
            where each tuple contains the indices of the two atoms involved in the bond
            and the bond order

    Returns:
        G (networkx.Graph): a NetworkX graph object representing the molecule,
            where each node represents an atom and each edge represents a bond,
            and each node and edge has attributes for element symbol, charge,
            and bond order, respectively
    """
    G = nx.Graph()
    for i, atom in enumerate(atoms):
        # add atom node to graph
        node_label = str(i + 1)
        G.add_node(node_label, element=atom['element'])
    for bond in bonds:
        # add bond edge to graph
        atom1 = str(bond[0])
        atom2 = str(bond[1])
        bond_type = bond[2]
        G.add_edge(atom1, atom2, order=bond_type)
    return G
  
  def canonical_Morgan(G):
    # initialize invariants and labels
    inv = {i: 1 for i in G.nodes()}
    lab = {i: 0 for i in G.nodes()}
    # compute invariants
    inv = compute_invariant(G, inv)

    # Continue
    
    return inv, lab
  
  def compute_invariant(G, inv):
    inv_length = len(set(inv.values()))
    # compute invariant of each atom
    for x in G.nodes():
        pass
