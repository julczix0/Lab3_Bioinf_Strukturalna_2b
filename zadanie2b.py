import os
from Bio import PDB
from Bio.PDB import PDBParser
from Bio.PDB import Polypeptide
import numpy as np
import matplotlib.pyplot as plt

def calculate_angles(structure):
    val_phi = []
    val_psi = []
    secondary_structure = []

    for model in structure:
        for chain in model:
            polypeptide = Polypeptide.Polypeptide(chain)  # obiekt klasy Polypeptide, analiza sekwencji aminokwasów
            angles = polypeptide.get_phi_psi_list()  # lista katów phi i psi
            dssp = PDB.DSSP(model, fr)  # obiekt klasy DSSP, analiza struktury wtórnej

            residues = list(chain.get_residues())  # lista reszt aminokwasowych
            for phi, psi in angles:  # iteracja po parach kątów
                if phi is not None and psi is not None:  # warunek jeśli kąty są zdefiniowane
                    residue_index = polypeptide.get_phi_psi_list().index((phi, psi))  # indeks reszty aminokwasowej w polypeptide
                    residue_id = (chain.id, residues[residue_index].get_id()[1])  # identyfikator reszty aminokwasowej, identyfikator łańcucha i indeks reszty
                    val_phi.append(np.degrees(phi))
                    val_psi.append(np.degrees(psi))  # dodanie katów w stopniach do list
                    secondary_structure.append(dssp[residue_id][2])  # dodanie informacji o strukturze wtórnej do listy (trzeci element na liscie)

    return val_phi, val_psi, secondary_structure


def plot_ramachandran(val_phi, val_psi, secondary_structure):
    plt.figure(figsize=(8, 6))  # stworzenie plot 8x6

    # kolory dla typów struktury wtórnej
    color_map = {'H': 'red',
                 'B': 'orange',
                 'E': 'yellow',
                 'G': 'green',
                 'I': 'blue',
                 'T': 'purple',
                 'S': 'blue',
                 '-': 'grey'}

    # mapowanie kodów struktur wtórnych do ich nazw (legenda)
    full_names = {
        'H': 'Alpha helix',
        'B': 'Beta bridge',
        'E': 'Strand',
        'G': 'Helix-3',
        'I': 'Helix-5',
        'T': 'Turn',
        'S': 'Bend',
        '-': 'Unknown'
    }

    handles = []

    for label, color in color_map.items():
        handles.append(plt.Line2D([0], [0],
                                  marker='o', color='w', markerfacecolor=color, markersize=10,
                                  # punkty na legendzie
                                  label=f'{full_names[label]} ({label})'))  # etykieta dla legendy

    for phi, psi, ss in zip(val_phi, val_psi, secondary_structure):  # lączenie w pary elementów
        plt.scatter(phi, psi, s=10, color=color_map[ss])  # rysowanie punktów

    plt.title('Ramachandran Plot')
    plt.xlabel('Phi Angle (degrees)')
    plt.ylabel('Psi Angle (degrees)')
    plt.xlim(-180, 180)
    plt.ylim(-180, 180)
    plt.axhline(0, color='black', linewidth=0.5)
    plt.axvline(0, color='black', linewidth=0.5)
    plt.grid(color='gray', linestyle='--', linewidth=0.5, alpha=0.5)

    plt.legend(handles=handles, title='Secondary Structure', loc='lower right')  # pozycja legendy

    plt.savefig('ramachandran_plot.pdf')
    plt.show()


fr = input("PDB FILE NAME: ")
structure = PDBParser(PERMISSIVE=False, QUIET=False).get_structure(fr[-8:-4], os.path.basename(fr))

phi_values, psi_values, secondary_structure = calculate_angles(structure)

plot_ramachandran(phi_values, psi_values, secondary_structure)
