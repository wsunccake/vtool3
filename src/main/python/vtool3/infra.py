#!/usr/bin/env python3

import math

__author__ = ""
__date__ = "2016/09/05"
__version__ = "$Revision: 0.3$"


class Vector(object):
    def __init__(self, x, y, z):
        self.__x = x
        self.__y = y
        self.__z = z

    def get_basis(self):
        return self.__x, self.__y, self.__z

    def set_basis(self, x, y, z):
        self.__x = x
        self.__y = y
        self.__z = z

    basis = property(get_basis, 'basis property')

    def normalized(self):
        length = self.length()
        if length == 0.0:
            print("Don't normalized")
            return Vector(self.__x, self.__y, self.__z)
        else:
            return Vector(self.__x / length, self.__y / length, self.__z / length)

    def length(self):
        length = self.__x * self.__x + self.__y * self.__y + self.__z * self.__z
        return math.sqrt(length)

    def angle(self, vector):
        angle = math.acos(self.dot(vector) / self.length() / vector.length())
        return math.degrees(angle)

    def dot(self, vector):
        return self.__x * vector.__x + self.__y * vector.__y + self.__z * vector.__z

    def cross(self, vector):
        x = self.__y * vector.__z - self.__z * vector.__y
        y = self.__z * vector.__x - self.__x * vector.__z
        z = self.__x * vector.__y - self.__y * vector.__x
        return Vector(x, y, z)

    def __str__(self):
        return "%f, %f, %f" %(self.__x, self.__y, self.__z)

    def __add__(self, vector):
        return Vector(self.__x + vector.__x, self.__y + vector.__y, self.__z + vector.__z)

    def __sub__(self, vector):
        return Vector(self.__x - vector.__x, self.__y - vector.__y, self.__z - vector.__z)

    def __mul__(self, n):
        return Vector(n * self.__x, n * self.__y, n * self.__z)

    def __eq__(self, vector):
        x, y, z = vector.basis
        return True if (self.__x == x and self.__y == y and self.__z == z) else False

    def rotate(self, axis_vector, angle):
        """ rotating vector around arbitrary axis
                 cosQ+Nx^2(1-cosQ)        NxNy(1-cosQ) - Nz*sinQ   NxNz(1-cosQ) + Ny*sinQ     x1   x2   x3
            R =  NyNz(1-cosQ) + Nz*sinQ   cosQ + Ny^2(1-cosQ)      NyNz(1-cosQ) - Nz*sinQ  =  y1   y2   y3
                 NzNx(1-cosQ) - Ny*sinQ   NzNy(1-cosQ) + NxsinQ    cosQ + Nz^2(1-csoQ)a       z1   z2   z3
            ref: http://en.wikipedia.org/wiki/Rotation_matrix#Conversion_from_and_to_axis-angle
        """
        angle = math.radians(angle)
        x1 = math.cos(angle) + axis_vector.__x * axis_vector.__x * (1 - math.cos(angle))
        x2 = axis_vector.__x * axis_vector.__y * (1 - math.cos(angle)) - axis_vector.__z * math.sin(angle)
        x3 = axis_vector.__x * axis_vector.__z * (1 - math.cos(angle)) + axis_vector.__y * math.sin(angle)

        y1 = axis_vector.__y * axis_vector.__z * (1 - math.cos(angle)) + axis_vector.__z * math.sin(angle)
        y2 = math.cos(angle) + axis_vector.__y * axis_vector.__y * (1 - math.cos(angle))
        y3 = axis_vector.__y * axis_vector.__z * (1 - math.cos(angle)) - axis_vector.__z * math.sin(angle)

        z1 = axis_vector.__z * axis_vector.__x * (1 - math.cos(angle)) - axis_vector.__y * math.sin(angle)
        z2 = axis_vector.__z * axis_vector.__y * (1 - math.cos(angle)) + axis_vector.__x * math.sin(angle)
        z3 = math.cos(angle) + axis_vector.__z * axis_vector.__z * (1 - math.cos(angle))

        x = x1 * self.__x + x2 * self.__y + x3 * self.__z
        y = y1 * self.__x + y2 * self.__y + y3 * self.__z
        z = z1 * self.__x + z2 * self.__y + z3 * self.__z
        return Vector(x, y, z)


class Lattice(object):
    def __init__(self, v1, v2, v3, constant=1.0):
        """ set default argment
            v1:       lattice vector1  {Vector}
            v2:       lattice vector2  {Vector}
            v3:       lattice vector3  {Vector}
            constant: lattice constant {Number}
        """
        self.__vectors = [v1, v2, v3]
        self.__constant = constant

    @property
    def constant(self):
        return self.__constant

    @constant.setter
    def constant(self, const):
        self.__constant = const

    def get_vectors(self):
        return self.__vectors

    def set_vectors(self, v1, v2, v3):
        self.__vectors = [v1, v2, v3]

    vectors = property(get_vectors, 'vectors property')


class Element(object):
    def __init__(self, symbol='X', name='Dummy', number=0, mass=0.0):
        """ set default
            symbol:   chemical symbol {String}  [X]
            name:     name of element {String}  [Dummy]
            number:   atomic number   {Int}     [0]
            mass:     atomic mass     {Float}   [0.0]
        """
        self.__symbol = symbol
        self.__name = name
        self.__atomic_number = number
        self.__atomic_mass = mass

    @property
    def symbol(self):
        return self.__symbol

    @symbol.setter
    def symbol(self, s='X'):
        self.__symbol = s

    @property
    def name(self):
        return self.__name

    @name.setter
    def name(self, n='Dummy'):
        self.__name = n

    @property
    def atomic_number(self):
        return self.__atomic_number

    @atomic_number.setter
    def atomic_number(self, an=0):
        self.__atomic_number = an

    @property
    def atomic_mass(self):
        return self.__atomic_mass

    @atomic_mass.setter
    def atomic_mass(self, am=0.0):
        self.__atomic_mass = am

    def copy_element(self, e):
        self.symbol = e.symbol
        self.name = e.name
        self.atomic_number = e.atomic_number
        self.atomic_mass = e.atomic_mass

    def __str__(self):
        return '%s, %s, %d, %d' %(self.__symbol, self.__name, self.__atomic_number, self.__atomic_mass)


class Atom(Element):
    """ Atom basic info """
    def __init__(self, element_symbol='X',
                 x_coordinate=0.0, y_coordinate=0.0, z_coordinate=0.0,
                 x_dynamic='T', y_dynamic='T', z_dynamic='T',
                 x_displace=0.0, y_displace=0.0, z_displace=0.0):
        """ set default argments
            elelment:    atomic element
            x_coordinate: x-axix coordinate
            y_coordinate: y-axix coordinate
            z_coordinate: z-axix coordinate
            x_dynamic:    x-axix (T)ranslate/(F)reeze
            y_dynamic:    y-axix (T)ranslate/(F)reeze
            z_dynamic:    z-axix (T)ranslate/(F)reeze
            x_displace:   x-axix displacement
            y_displace:   y-axix displacement
            z_Displace:   z-axix displacement
        """

        self.symbol = element_symbol
        check_element_by_periodic_table(self)
        self.__x_coordinate = x_coordinate
        self.__y_coordinate = y_coordinate
        self.__z_coordinate = z_coordinate
        self.__x_dynamic = x_dynamic
        self.__y_dynamic = y_dynamic
        self.__z_dynamic = z_dynamic
        self.__x_displace = x_displace
        self.__y_displace = y_displace
        self.__z_displace = z_displace

    def get_coordinate(self):
        return self.__x_coordinate, self.__y_coordinate, self.__z_coordinate

    def set_coordinate(self, x, y, z):
        self.__x_coordinate = x
        self.__y_coordinate = y
        self.__z_coordinate = z

    coordinate = property(get_coordinate, 'coordinate property')

    def get_dynamic(self):
        return self.__x_dynamic, self.__y_dynamic, self.__z_dynamic

    def set_dynamic(self, x_dynamic, y_dynamic, z_dynamic):
        self.__x_dynamic = x_dynamic
        self.__y_dynamic = y_dynamic
        self.__z_dynamic = z_dynamic

    dynamic = property(get_dynamic, 'dynamic property')

    def get_displace(self):
        return self.__x_displace, self.__y_displace, self.__z_displace

    def set_displace(self, x_displace, y_displace, z_displace):
        self.__x_displace = x_displace
        self.__y_displace = y_displace
        self.__z_displace = z_displace

    displace = property(get_displace, 'displace property')

    def __repr__(self):  
        return repr((self.symbol, self.__x_coordinate, self.__y_coordinate, self.__z_coordinate))

    def __eq__(self, atom):
        if (self.symbol == atom.symbol):
            if (self.coordinate == atom.coordinate):
                return True
        return False


    def show_atom(self):
        print("%3s, %+14.10f, %+14.10f, %+14.10f, %4s, %4s, %4s" %(self.symbol, self.__x_coordinate, self.__y_coordinate, self.__z_coordinate, self.__x_dynamic, self.__y_dynamic, self.__z_dynamic))

    def copy_atom(self, atom):
        self.copy_element(atom)
        self.set_coordinate(atom.__x_coordinate, atom.__y_coordinate, atom.__z_coordinate)
        self.set_dynamic(atom.__x_dynamic, atom.__y_dynamic, atom.__z_dynamic)
        self.set_displace(atom.__x_displace, atom.__y_displace, atom.__z_displace)
        pass

    def add_coordinate(self, x, y, z):
        new_x = self.__x_coordinate + x
        new_y = self.__y_coordinate + y
        new_z = self.__z_coordinate + z
        return Atom(self.symbol, new_x, new_y, new_z, self.__x_dynamic, self.__y_dynamic, self.__z_dynamic)

    def sub_coordinate(self, x, y, z):
        new_x = self.__x_coordinate - x
        new_y = self.__y_coordinate - y
        new_z = self.__z_coordinate - z
        return Atom(self.symbol, new_x, new_y, new_z, self.__x_dynamic, self.__y_dynamic, self.__z_dynamic)

    def mul_coordinate(self, f):
        new_x = self.__x_coordinate * f
        new_y = self.__y_coordinate * f
        new_z = self.__z_coordinate * f
        return Atom(self.symbol, new_x, new_y, new_z, self.__x_dynamic, self.__y_dynamic, self.__z_dynamic)

    def div_coordinate(self, f):
        new_x = self.__x_coordinate / f
        new_y = self.__y_coordinate / f
        new_z = self.__z_coordinate / f
        return Atom(self.symbol, new_x, new_y, new_z, self.__x_dynamic, self.__y_dynamic, self.__z_dynamic)


""" default periodic table element:
    global constant: PERIODIC_TABLE_ElEMENTS
"""
PERIODIC_TABLE_ElEMENTS = [Element(),
                           Element('H',   'Hydrogen',      1,   1.00794),
                           Element('He',  'Helium',        2,   4.002602),
                           Element('Li',  'Lithium',       3,   6.941),
                           Element('Be',  'Beryllium',     4,   9.0121831),
                           Element('B',   'Boron',         5,   10.81),
                           Element('C',   'Carbon',        6,   12.011),
                           Element('N',   'Nitrogen',      7,   14.007),
                           Element('O',   'Oxygen',        8,   15.999),
                           Element('F',   'Fluorine',      9,   18.998403163),
                           Element('Ne',  'Neon',          10,  20.1797),
                           Element('Na',  'Sodium',        11,  22.98976928),
                           Element('Mg',  'Magnesium',     12,  24.305),
                           Element('Al',  'Aluminium',     13,  26.9815385),
                           Element('Si',  'Silicon',       14,  28.085),
                           Element('P',   'Phosphorus',    15,  30.973761998),
                           Element('S',   'Sulfur',        16,  32.066),
                           Element('Cl'   'Chlorine',      17,  35.45),
                           Element('Ar',  'Argon',         18,  39.948),
                           Element('K',   'Potassium',     19,  39.0983),
                           Element('Ca',  'Calcium',       20,  40.078),
                           Element('Sc',  'Scandium',      21,  44.955908),
                           Element('Ti',  'Titanium',      22,  47.867),
                           Element('V',   'Vanadium',      23,  50.9415),
                           Element('Cr',  'Chromium',      24,  51.9961),
                           Element('Mn',  'Manganese',     25,  54.938044),
                           Element('Fe',  'Iron',          26,  55.845),
                           Element('Co',  'Cobalt',        27,  58.933194),
                           Element('Ni',  'Nickel',        28,  58.6934),
                           Element('Cu',  'Copper',        29,  63.546),
                           Element('Zn',  'Zinc',          30,  65.38),
                           Element('Ga',  'Gallium',       31,  69.723),
                           Element('Ge',  'Germanium',     32,  72.630),
                           Element('As',  'Arsenic',       33,  74.921595),
                           Element('Se',  'Selenium',      34,  78.971),
                           Element('Br',  'Bromine',       35,  79.904),
                           Element('Kr',  'Krypton',       36,  83.798),
                           Element('Rb',  'Rubidium',      37,  85.4678),
                           Element('Sr',  'Strontium',     38,  87.62),
                           Element('Y',   'Yttrium',       39,  88.90584),
                           Element('Zr',  'Zirconium',     40,  91.224),
                           Element('Nb',  'Niobium',       41,  92.90637),
                           Element('Mo',  'Molybdenum',    42,  95.95),
                           Element('Tc',  'Technetium',    43,  98),
                           Element('Ru',  'Ruthenium',     44,  101.07),
                           Element('Rh',  'Rhodium',       45,  102.90550),
                           Element('Pd',  'Palladium',     46,  106.42),
                           Element('Ag',  'Silver',        47,  107.8682),
                           Element('Cd',  'Cadmium',       48,  112.414),
                           Element('In',  'Indium',        49,  114.818),
                           Element('Sn',  'Tin',           50,  118.710),
                           Element('Sb',  'Antimony',      51,  121.760),
                           Element('Te',  'Tellurium',     52,  127.60),
                           Element('I',   'Iodine',        53,  126.90447),
                           Element('Xe',  'Xenon',         54,  131.293),
                           Element('Cs',  'Caesium',       55,  132.90545196),
                           Element('Ba',  'Barium',        56,  137.327),
                           Element('La',  'Lanthanum',     57,  138.90547),
                           Element('Ce',  'Cerium',        58,  140.116),
                           Element('Pr',  'Praseodymium',  59,  140.90766),
                           Element('Nd',  'Neodymium',     60,  144.242),
                           Element('Pm',  'Promethium',    61,  145),
                           Element('Sm',  'Samarium',      62,  150.36),
                           Element('Eu',  'Europium',      63,  151.964),
                           Element('Gd',  'Gadolinium',    64,  157.25),
                           Element('Tb',  'Terbium',       65,  158.92535),
                           Element('Dy',  'Dysprosium',    66,  162.500),
                           Element('Ho',  'Holmium',       67,  164.93033),
                           Element('Er',  'Erbium',        68,  167.259),
                           Element('Tm',  'Thulium',       69,  168.93422),
                           Element('Yb',  'Ytterbium',     70,  173.054),
                           Element('Lu',  'Lutetium',      71,  174.9668),
                           Element('Hf',  'Hafnium',       72,  178.49),
                           Element('Ta',  'Tantalum',      73,  180.94788),
                           Element('W',   'Tungsten',      74,  183.84),
                           Element('Re',  'Rhenium',       75,  186.207),
                           Element('Os',  'Osmium',        76,  190.23),
                           Element('Ir',  'Iridium',       77,  192.217),
                           Element('Pt',  'Platinum',      78,  195.084),
                           Element('Au',  'Gold',          79,  196.966569),
                           Element('Hg',  'Mercury',       80,  200.592),
                           Element('Tl',  'Thallium',      81,  204.38),
                           Element('Pb',  'Lead',          82,  207.2),
                           Element('Bi',  'Bismuth',       83,  208.98040),
                           Element('Po',  'Polonium',      84,  209),
                           Element('At',  'Astatine',      85,  210),
                           Element('Rn',  'Radon',         86,  222),
                           Element('Fr',  'Francium',      87,  223),
                           Element('Ra',  'Radium',        88,  226),
                           Element('Ac',  'Actinium',      89,  227),
                           Element('Th',  'Thorium',       90,  232.0377),
                           Element('Pa',  'Protactinium',  91,  231.03588),
                           Element('U',   'Uranium',       92,  238.02891),
                           Element('Np',  'Neptunium',     93,  237),
                           Element('Pu',  'Plutonium',     94,  244),
                           Element('Am',  'Americium',     95,  243),
                           Element('Cm',  'Curium',        96,  247),
                           Element('Bk',  'Berkelium',     97,  247),
                           Element('Cf',  'Californium',   98,  251),
                           Element('Es',  'Einsteinium',   99,  252),
                           Element('Fm',  'Fermium',       100, 257),
                           Element('Md',  'Mendelevium',   101, 258),
                           Element('No',  'Nobelium',      102, 259),
                           Element('Lr',  'Lawrencium',    103, 262),
                           Element('Rf',  'Rutherfordium', 104, 267),
                           Element('Db',  'Dubnium',       105, 268),
                           Element('Sg',  'Seaborgium',    106, 269),
                           Element('Bh',  'Bohrium',       107, 270),
                           Element('Hs',  'Hassium',       108, 269),
                           Element('Mt',  'Meitnerium',    109, 278),
                           Element('Ds',  'Darmstadtium',  110, 281),
                           Element('Rg',  'Roentgenium',   111, 281),
                           Element('Cn',  'Copernicium',   112, 285),
                           Element('Uut', 'Ununtrium',     113, 286),
                           Element('Fl',  'Flerovium',     114, 289),
                           Element('Uup', 'Ununpentium',   115, 288),
                           Element('Lv',  'Livermorium',   116, 293),
                           Element('Uus', 'Ununseptium',   117, 294),
                           Element('Uuo', 'Ununoctium',    118, 294),
                           Element('Tv',  'Translation vectors', -2, 0), # Defined in Gaussian
                          ]


def check_element_by_periodic_table(element, method='symbol'):
    """ setup element all properties in peridic table
    """
    if method == 'symbol':
        for e in PERIODIC_TABLE_ElEMENTS:
            if e.symbol == element.symbol:
                element.copy_element(e)
                break
            element.name = 'dummy'
            element.atomic_number = 0
            element.atomic_mass = 0.0
    elif method == 'name':
        for e in PERIODIC_TABLE_ElEMENTS:
            if e.name == element.name:
                element.copy_element(e)
                break
            element.symbol = 'X'
            element.atomic_number = 0
            element.atomic_mass = 0.0
    elif method == 'number':
        for e in PERIODIC_TABLE_ElEMENTS:
            if e.atomic_number == element.atomic_number:
                element.copy_element(e)
                break
            element.symbol = 'X'
            element.name = 'dummy'
            element.atomic_mass = 0.0


if __name__ == "__main__":
    pass
