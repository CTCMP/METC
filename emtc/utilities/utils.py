# -*- coding: utf-8 -*-
from math import pi, sqrt

QE_atom_mass_factor = 1822.8884862 / 2  # amu
kb_J = 1.3806504e-23  # [J/K]
PlanckConstant = 4.13566733e-15  # [eV s]
Hbar = PlanckConstant / (2 * pi)  # [eV s]
Avogadro = 6.02214179e23
SpeedOfLight = 299792458  # [m/s]
AMU = 1.6605402e-27  # [kg]
Newton = 1.0  # [kg m / s^2]
Joule = 1.0  # [kg m^2 / s^2]
EV = 1.60217733e-19  # [J]
Angstrom = 1.0e-10  # [m]
THz = 1.0e12  # [/s]
Mu0 = 4.0e-7 * pi  # [Hartree/m]
Epsilon0 = 1.0 / Mu0 / SpeedOfLight**2  # [C^2 / N m^2]
Me = 9.10938215e-31


Bohr = 4e10 * pi * Epsilon0 * Hbar**2 / Me  # Bohr radius [A] 0.5291772
AngToBohr = 1.8897259886
Hartree = Me * EV / 16 / pi**2 / Epsilon0**2 / Hbar**2  # Hartree [eV] 27.211398
Rydberg = Hartree / 2  # Rydberg [eV] 13.6056991

THzToEv = PlanckConstant * 1e12  # [eV]
Kb = kb_J / EV  # [eV/K] 8.6173383e-05
THzToCm = 1.0e12 / (SpeedOfLight * 100)  # [cm^-1] 33.356410
CmToEv = THzToEv / THzToCm  # [eV] 1.2398419e-4
VaspToEv = sqrt(EV / AMU) / Angstrom / (2 * pi) * PlanckConstant  # [eV] 6.46541380e-2
VaspToTHz = sqrt(EV / AMU) / Angstrom / (2 * pi) / 1e12  # [THz] 15.633302
VaspToCm = VaspToTHz * THzToCm  # [cm^-1] 521.47083
EvTokJmol = EV / 1000 * Avogadro  # [kJ/mol] 96.4853910
Wien2kToTHz = (
    sqrt(Rydberg / 1000 * EV / AMU) / (Bohr * 1e-10) / (2 * pi) / 1e12
)  # [THz] 3.44595837
AbinitToTHz = sqrt(EV / (AMU * Bohr)) / Angstrom / (2 * pi) / 1e12  # [THz] 21.49068
PwscfToTHz = (
    sqrt(Rydberg * EV / AMU) / (Bohr * 1e-10) / (2 * pi) / 1e12
)  # [THz] 108.97077
ElkToTHz = (
    sqrt(Hartree * EV / AMU) / (Bohr * 1e-10) / (2 * pi) / 1e12
)  # [THz] 154.10794
SiestaToTHz = sqrt(EV / (AMU * Bohr)) / Angstrom / (2 * pi) / 1e12  # [THz] 21.49068
CP2KToTHz = (
    sqrt(Hartree * EV / (AMU * Bohr)) / Angstrom / (2 * pi) / 1e12
)  # CP2K uses a.u. for forces but Angstrom for distances
CrystalToTHz = VaspToTHz
CastepToTHz = VaspToTHz
DftbpToTHz = (
    sqrt(Hartree * EV / AMU) / (Bohr * 1e-10) / (2 * pi) / 1e12
)  # [THz] 154.10794344
dftbpToBohr = 0.188972598857892e01
TurbomoleToTHz = ElkToTHz  # Turbomole uses atomic units (Hartree/Bohr)
EVAngstromToGPa = EV * 1e21
FleurToTHz = ElkToTHz  # Fleur uses atomic units (Hartree/Bohr)

atom_data = [
    [0, "X", "X", None],  # 0
    [1, "H", "Hydrogen", 1.00794],  # 1
    [2, "He", "Helium", 4.002602],  # 2
    [3, "Li", "Lithium", 6.941],  # 3
    [4, "Be", "Beryllium", 9.012182],  # 4
    [5, "B", "Boron", 10.811],  # 5
    [6, "C", "Carbon", 12.0107],  # 6
    [7, "N", "Nitrogen", 14.0067],  # 7
    [8, "O", "Oxygen", 15.9994],  # 8
    [9, "F", "Fluorine", 18.9984032],  # 9
    [10, "Ne", "Neon", 20.1797],  # 10
    [11, "Na", "Sodium", 22.98976928],  # 11
    [12, "Mg", "Magnesium", 24.3050],  # 12
    [13, "Al", "Aluminium", 26.9815386],  # 13
    [14, "Si", "Silicon", 28.0855],  # 14
    [15, "P", "Phosphorus", 30.973762],  # 15
    [16, "S", "Sulfur", 32.065],  # 16
    [17, "Cl", "Chlorine", 35.453],  # 17
    [18, "Ar", "Argon", 39.948],  # 18
    [19, "K", "Potassium", 39.0983],  # 19
    [20, "Ca", "Calcium", 40.078],  # 20
    [21, "Sc", "Scandium", 44.955912],  # 21
    [22, "Ti", "Titanium", 47.867],  # 22
    [23, "V", "Vanadium", 50.9415],  # 23
    [24, "Cr", "Chromium", 51.9961],  # 24
    [25, "Mn", "Manganese", 54.938045],  # 25
    [26, "Fe", "Iron", 55.845],  # 26
    [27, "Co", "Cobalt", 58.933195],  # 27
    [28, "Ni", "Nickel", 58.6934],  # 28
    [29, "Cu", "Copper", 63.546],  # 29
    [30, "Zn", "Zinc", 65.38],  # 30
    [31, "Ga", "Gallium", 69.723],  # 31
    [32, "Ge", "Germanium", 72.64],  # 32
    [33, "As", "Arsenic", 74.92160],  # 33
    [34, "Se", "Selenium", 78.96],  # 34
    [35, "Br", "Bromine", 79.904],  # 35
    [36, "Kr", "Krypton", 83.798],  # 36
    [37, "Rb", "Rubidium", 85.4678],  # 37
    [38, "Sr", "Strontium", 87.62],  # 38
    [39, "Y", "Yttrium", 88.90585],  # 39
    [40, "Zr", "Zirconium", 91.224],  # 40
    [41, "Nb", "Niobium", 92.90638],  # 41
    [42, "Mo", "Molybdenum", 95.96],  # 42
    [43, "Tc", "Technetium", 98],  # 43 (mass is from wikipedia)
    [44, "Ru", "Ruthenium", 101.07],  # 44
    [45, "Rh", "Rhodium", 102.90550],  # 45
    [46, "Pd", "Palladium", 106.42],  # 46
    [47, "Ag", "Silver", 107.8682],  # 47
    [48, "Cd", "Cadmium", 112.411],  # 48
    [49, "In", "Indium", 114.818],  # 49
    [50, "Sn", "Tin", 118.710],  # 50
    [51, "Sb", "Antimony", 121.760],  # 51
    [52, "Te", "Tellurium", 127.60],  # 52
    [53, "I", "Iodine", 126.90447],  # 53
    [54, "Xe", "Xenon", 131.293],  # 54
    [55, "Cs", "Caesium", 132.9054519],  # 55
    [56, "Ba", "Barium", 137.327],  # 56
    [57, "La", "Lanthanum", 138.90547],  # 57
    [58, "Ce", "Cerium", 140.116],  # 58
    [59, "Pr", "Praseodymium", 140.90765],  # 59
    [60, "Nd", "Neodymium", 144.242],  # 60
    [61, "Pm", "Promethium", 145],  # 61 (mass is from wikipedia)
    [62, "Sm", "Samarium", 150.36],  # 62
    [63, "Eu", "Europium", 151.964],  # 63
    [64, "Gd", "Gadolinium", 157.25],  # 64
    [65, "Tb", "Terbium", 158.92535],  # 65
    [66, "Dy", "Dysprosium", 162.500],  # 66
    [67, "Ho", "Holmium", 164.93032],  # 67
    [68, "Er", "Erbium", 167.259],  # 68
    [69, "Tm", "Thulium", 168.93421],  # 69
    [70, "Yb", "Ytterbium", 173.054],  # 70
    [71, "Lu", "Lutetium", 174.9668],  # 71
    [72, "Hf", "Hafnium", 178.49],  # 72
    [73, "Ta", "Tantalum", 180.94788],  # 73
    [74, "W", "Tungsten", 183.84],  # 74
    [75, "Re", "Rhenium", 186.207],  # 75
    [76, "Os", "Osmium", 190.23],  # 76
    [77, "Ir", "Iridium", 192.217],  # 77
    [78, "Pt", "Platinum", 195.084],  # 78
    [79, "Au", "Gold", 196.966569],  # 79
    [80, "Hg", "Mercury", 200.59],  # 80
    [81, "Tl", "Thallium", 204.3833],  # 81
    [82, "Pb", "Lead", 207.2],  # 82
    [83, "Bi", "Bismuth", 208.98040],  # 83
    [84, "Po", "Polonium", None],  # 84
    [85, "At", "Astatine", None],  # 85
    [86, "Rn", "Radon", None],  # 86
    [87, "Fr", "Francium", None],  # 87
    [88, "Ra", "Radium", None],  # 88
    [89, "Ac", "Actinium", 227],  # 89 (mass is from wikipedia)
    [90, "Th", "Thorium", 232.03806],  # 90
    [91, "Pa", "Protactinium", 231.03588],  # 91
    [92, "U", "Uranium", 238.02891],  # 92
    [93, "Np", "Neptunium", 237],  # 93 (mass is from wikipedia)
    [94, "Pu", "Plutonium", None],  # 94
    [95, "Am", "Americium", None],  # 95
    [96, "Cm", "Curium", None],  # 96
    [97, "Bk", "Berkelium", None],  # 97
    [98, "Cf", "Californium", None],  # 98
    [99, "Es", "Einsteinium", None],  # 99
    [100, "Fm", "Fermium", None],  # 100
    [101, "Md", "Mendelevium", None],  # 101
    [102, "No", "Nobelium", None],  # 102
    [103, "Lr", "Lawrencium", None],  # 103
    [104, "Rf", "Rutherfordium", None],  # 104
    [105, "Db", "Dubnium", None],  # 105
    [106, "Sg", "Seaborgium", None],  # 106
    [107, "Bh", "Bohrium", None],  # 107
    [108, "Hs", "Hassium", None],  # 108
    [109, "Mt", "Meitnerium", None],  # 109
    [110, "Ds", "Darmstadtium", None],  # 110
    [111, "Rg", "Roentgenium", None],  # 111
    [112, "Cn", "Copernicium", None],  # 112
    [113, "Uut", "Ununtrium", None],  # 113
    [114, "Uuq", "Ununquadium", None],  # 114
    [115, "Uup", "Ununpentium", None],  # 115
    [116, "Uuh", "Ununhexium", None],  # 116
    [117, "Uus", "Ununseptium", None],  # 117
    [118, "Uuo", "Ununoctium", None],  # 118
]

symbol_data = {
    "H": 1,
    "He": 2,
    "Li": 3,
    "Be": 4,
    "B": 5,
    "C": 6,
    "N": 7,
    "O": 8,
    "F": 9,
    "Ne": 10,
    "Na": 11,
    "Mg": 12,
    "Al": 13,
    "Si": 14,
    "P": 15,
    "S": 16,
    "Cl": 17,
    "Ar": 18,
    "K": 19,
    "Ca": 20,
    "Sc": 21,
    "Ti": 22,
    "V": 23,
    "Cr": 24,
    "Mn": 25,
    "Fe": 26,
    "Co": 27,
    "Ni": 28,
    "Cu": 29,
    "Zn": 30,
    "Ga": 31,
    "Ge": 32,
    "As": 33,
    "Se": 34,
    "Br": 35,
    "Kr": 36,
    "Rb": 37,
    "Sr": 38,
    "Y": 39,
    "Zr": 40,
    "Nb": 41,
    "Mo": 42,
    "Tc": 43,
    "Ru": 44,
    "Rh": 45,
    "Pd": 46,
    "Ag": 47,
    "Cd": 48,
    "In": 49,
    "Sn": 50,
    "Sb": 51,
    "Te": 52,
    "I": 53,
    "Xe": 54,
    "Cs": 55,
    "Ba": 56,
    "La": 57,
    "Ce": 58,
    "Pr": 59,
    "Nd": 60,
    "Pm": 61,
    "Sm": 62,
    "Eu": 63,
    "Gd": 64,
    "Tb": 65,
    "Dy": 66,
    "Ho": 67,
    "Er": 68,
    "Tm": 69,
    "Yb": 70,
    "Lu": 71,
    "Hf": 72,
    "Ta": 73,
    "W": 74,
    "Re": 75,
    "Os": 76,
    "Ir": 77,
    "Pt": 78,
    "Au": 79,
    "Hg": 80,
    "Tl": 81,
    "Pb": 82,
    "Bi": 83,
    "Po": 84,
    "At": 85,
    "Rn": 86,
    "Fr": 87,
    "Ra": 88,
    "Ac": 89,
    "Th": 90,
    "Pa": 91,
    "U": 92,
    "Np": 93,
    "Pu": 94,
    "Am": 95,
    "Cm": 96,
    "Bk": 97,
    "Cf": 98,
    "Es": 99,
    "Fm": 100,
    "Md": 101,
    "No": 102,
    "Lr": 103,
    "Rf": 104,
    "Db": 105,
    "Sg": 106,
    "Bh": 107,
    "Hs": 108,
    "Mt": 109,
    "Ds": 110,
    "Rg": 111,
    "Cn": 112,
    "Uut": 113,
    "Uuq": 114,
    "Uup": 115,
    "Uuh": 116,
    "Uus": 117,
    "Uuo": 118,
}

atom_mass = {
  "H"   : 1.00794,
  "He"  : 4.0026,
  "Li"  : 6.941,
  "Be"  : 9.01218,
  "B"   : 10.811,
  "C"   : 12.01070,
  "N"   : 14.0067,
  "O"   : 15.9994,
  "F"   : 18.9984,
  "Ne"  : 20.1797, 
  "Na"  : 22.9898, 
  "Mg"  : 24.305,
  "Al"  : 26.9815, 
  "Si"  : 28.0850, 
  "P"   : 30.9738,
  "S"   : 32.066, 
  "Cl"  : 35.4527, 
  "Ar"  : 39.948, 
  "K"   : 39.0983,
  "Ca"  : 40.078,
  "Sc"  : 44.955912,
  "Ti"  : 47.88,
  "V"   : 50.9415,
  "Cr"  : 51.996,
  "Mn"  : 54.938,
  "Fe"  : 55.847,
  "Co"  : 58.9332,
  "Ni"  : 58.6934, 
  "Cu"  : 63.546,
  "Zn"  : 65.38, 
  "Ga"  : 69.723, 
  "Ge"  : 72.61,
  "As"  : 74.9216,
  "Se"  : 78.96, 
  "Br"  : 79.904,
  "Kr"  : 83.8, 
  "Rb"  : 85.4678,
  "Sr"  : 87.62, 
  "Y"   : 88.9059,
  "Zr"  : 91.224, 
  "Nb"  : 92.9064,
  "Mo"  : 95.94,
  "Tc"  : 98, 
  "Ru"  : 101.07,
  "Rh"  : 102.906, 
  "Pd"  : 106.42,
  "Ag"  : 107.868, 
  "Cd"  : 112.41,
  "In"  : 114.82,
  "Sn"  : 118.71, 
  "Sb"  : 121.757, 
  "Te"  : 127.6, 
  "I"   : 126.904, 
  "Xe"  : 131.29, 
  "Cs"  : 132.905, 
  "Ba"  : 137.33, 
  "La"  : 138.905, 
  "Ce"  : 140.12, 
  "Pr"  : 140.908, 
  "Nd"  : 144.24, 
  "Pm"  : 145, 
  "Sm"  : 150.36, 
  "Eu"  : 151.965, 
  "Gd"  : 157.25, 
  "Tb"  : 158.925, 
  "Dy"  : 162.5, 
  "Ho"  : 164.93, 
  "Er"  : 167.26, 
  "Tm"  : 168.934, 
  "Yb"  : 173.04, 
  "Lu"  : 174.967, 
  "Hf"  : 178.49, 
  "Ta"  : 180.948, 
  "W"   : 183.85, 
  "Re"  : 186.207, 
  "Os"  : 190.2, 
  "Ir"  : 192.22, 
  "Pt"  : 195.08, 
  "Au"  : 196.966, 
  "Hg"  : 200.59, 
  "Tl"  : 204.383, 
  "Pb"  : 207.2, 
  "Bi"  : 208.98, 
  "Po"  : 209, 
  "At"  : 210, 
  "Rn"  : 222, 
  "Fr"  : 223, 
  "Ra"  : 226.025, 
  "Ac"  : 227, 
  "Th"  : 232.038, 
  "Pa"  : 231.036, 
  "U"   : 238.029, 
  "Np"  : 237.048, 
  "Pu"  : 244, 
  "Am"  : 243, 
  "Cm"  : 247, 
  "Bk"  : 247, 
  "Cf"  : 251, 
  "Es"  : 252, 
  "Fm"  : 257, 
  "Md"  : 258, 
  "No"  : 259, 
  "Lr"  : 262, 
  "Rf"  : 261, 
  "Db"  : 262, 
  "Sg"  : 263,
  "Bh"  : 262, 
  "Hs"  : 265, 
  "Mt"  : 266, 
  "Uun" : 269, 
  "Uuu" : 272,
  "Uub" : 277,
}

def distance_Q1Q2(q1, q2):
  dis = 0
  for i in range(len(q1)):
    dis = dis + q2[i]**2 - q1[i]**2
  return dis