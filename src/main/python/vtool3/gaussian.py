#!/usr/bin/env python3

from vtool3.infra import *

__author__ = ""
__date__ = "2016/09/05"
__version__ = "$Revision: 0.3$"


class GJF(object):
    """ create/read Gaussian input file """
    def __init__(self, filename=None, comment=None):
        """ 
        """
        self.__specs = []
#        self._option_ = ''
        self.__option = '# opt freq hf/3-21g'
        self.__comment = 'This line is comment'
        self.__charge = 0
        self.__spin = 1
        self.__atoms = []
        self._elements_ = []
        self._element_lists_ = []
        self._numbers_ = []
        self.__lattice = None
        if not(filename is None):
            self.read_gjf(filename)

    def get_lattice(self):
        return self.__lattice

    def set_lattice(self, vectors, lattice_constant=1.0):
        """ set lattice
            vectors: {Vector array}
        """
        v1 = vectors[0]
        v2 = vectors[1]
        v3 = vectors[2]
        self.__lattice = Lattice(v1, v2, v3, lattice_constant)

    lattice = property(get_lattice, 'lattice property')

    def get_option(self):
        return self.__option

    def set_option(self, opt='# opt freq hf/3-21g'):
        self.__option = opt

    option = property(get_option, set_option, 'option property')

    def get_comment(self):
        return self.__comment

    def set_comment(self, comment='This line is comment'):
        self.__comment = comment

    comment = property(get_comment, set_comment, 'comment property')

    def get_charge(self):
        return self.__charge

    def set_charge(self, charge=0):
        self.__charge = charge

    charge = property(get_charge, set_charge, 'charge property')

    def get_spin(self):
        return self.__spin

    def set_spin(self, spin=1):
        self.__spin = spin

    spin = property(get_spin, set_spin, 'spin property')

    def get_atoms(self):
        return self.__atoms

    def add_atom(self, atom):
        self.__atoms.append(atom)

    def sort_atoms(self):
        self.__atoms = sorted(self.__atoms, key=lambda atom: (atom.symbol, atom.coordinate[2]))

    def read_gjf(self, filename):
        """ read Gaussian input file """
        f = open(filename, "r")
        i = 0
        opt_flag = -5
        for l in f.readlines():
            if l[0] == "%":
                self.__specs.append(l.rstrip() )
            elif l[0] == "#":
                self.option = l.rstrip()
                opt_flag = i
            elif l.strip() == "":
                 pass
            elif opt_flag + 2 == i:
                 self.comment = l.rstrip()
            elif len(l.split()) == 4:
                 atom = Atom(l.split()[0], float(l.split()[1]), float(l.split()[2]), float(l.split()[3]))
                 check_element_by_periodic_table(atom, 'symbol')
                 self.add_atom(atom)
            i += 1
        self.sort_atoms()
        self._check_lattice()

    def write_gjf(self, filename=None):
        """ write Gaussian input file """
        output = ''
        for l in self.__specs:
            output += l + "\n"
        output += self.__option + "\n\n"
        output += self.__comment + "\n\n"
        output += '%i %i\n' %(self.__charge, self.__spin)
        for a in self.__atoms:
            output += "%-2s" % a.symbol + "       %13.8f    %13.8f    %13.8f\n" % a.coordinate

        if isinstance(self.__lattice, Lattice):
            for v in self.__lattice.vectors:
                v1, v2, v3 = v.basis
                output += "Tv       %13.8f    %13.8f    %13.8f\n" %(v1, v2, v3)

        # setup output
        if filename is None:
            print(output)
        else:
            f = open(filename, "w")
            f.write(output)
            f.close()

    def _check_lattice(self):
        i = 0
        lattice_indexes = []
        vectors = []
        for atom in self.__atoms:
            if atom.symbol == "Tv":
                tc1, tc2, tc3 = atom.coordinate
                tmp_vec = Vector(tc1, tc2, tc3)
                vectors.append(tmp_vec)
                lattice_indexes.append(i)
            i += 1
        if len(lattice_indexes) == 3:
            lattice_indexes.reverse()
            for l in lattice_indexes:
                self.__atoms.pop(l)
            self.set_lattice(vectors)


class Log(object):
    """ create/read Gaussian output log """
    def __init__(self, filename=None):
        pass

    def set_atom(self):
        pass

    def set_frequency(self):
        self.__frequencies = []

    def read_frequency(self):
        pass

    def write_frequency(self):
        pass

    def write_log(self):
        pass

if __name__ == "__main__":
#    g = GJF('g.gjf')
#    g.writeGJF()
    pass
