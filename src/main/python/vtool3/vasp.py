#!/usr/bin/env python3

import copy
import re

from vtool3.infra import *

__author__ = ""
__date__ = "2016/10/22"
__version__ = "$Revision: 0.3$"


class POSCAR(object):
    """ create/read VASP POSCAR """
    def __init__(self, filename=None, comment='Comment line', lattice=None,
                 select='Selective', coordinate='Cartesian'):

        self.__comment = comment
        self.__lattice = lattice
        self.__selective_mode = select
        self.__coordinate_type = coordinate
        self.__atoms = []

        if not(filename is None):
            self.read_poscar(filename)

    def get_comment(self):
        return self.__comment

    def set_comment(self, comment):
        self.__comment = comment

    comment = property(get_comment, set_comment, 'comment property')

    def get_lattice(self):
        return self.__lattice

    def set_lattice(self, vectors, lattice_constant=1.0):
        """ set lattice
            vectors:         lattice vectors {vector array}
            lattice_constant:                 {number}
        """
        v1 = vectors[0]
        v2 = vectors[1]
        v3 = vectors[2]
        self.__lattice = Lattice(v1, v2, v3, lattice_constant)

    lattice = property(get_lattice, 'lattice property')

    def get_selective_mode(self):
        return self.__selective_mode

    def set_selective_mode(self, mode):
        self.__selective_mode = mode

    selective_mode = property(get_selective_mode, set_selective_mode, 'selective mode property')

    def get_coordinate_type(self):
        return self.__coordinate_type

    def set_coordinate_type(self, coordinate_type):
        self.__coordinate_type = coordinate_type

    coordinate_type = property(get_coordinate_type, set_coordinate_type, 'coordinate type property')

    def get_atoms(self):
        return self.__atoms

    def add_atom(self, atom):
        a = Atom()
        a.copy_atom(atom)
        self.__atoms.append(a)

    def del_atom(self):
        pass

    def set_atom_element(self, index, element_symbol):
        self.__atoms[index].symbol = element_symbol

    def set_atom_coordinate(self, index, x_coordinate, y_coordinate, z_coordinate):
        self.__atoms[index].set_coordinate(x_coordinate, y_coordinate, z_coordinate)

    def set_atom_dynamic(self, index, x_dynamic, y_dynamic, z_dynamic):
        self.__atoms[index].set_dynamic(x_dynamic, y_dynamic, z_dynamic)

    def set_atom(self, index, atom):
        self.__atoms[index] = atom

    def list_atom(self):
        for a in self.__atoms:
            a.show_atom()
    
    def _check__elements(self):
        tmp_element_type = self.__atoms[0].symbol
        tmp_element_number = 0
        elements = []
        for a in self.__atoms:
            e = a.symbol
            if tmp_element_type == e:
                tmp_element_number += 1
            else:
                elements.append({'element_type': tmp_element_type, 'element_number': tmp_element_number})
                tmp_element_type = e
                tmp_element_number = 1
        elements.append({'element_type': tmp_element_type, 'element_number': tmp_element_number})
        return elements

    def set_elements_type(self, elements):
        es = self._check__elements()
        n = 0
        for i in range(len(es)):
            for j in range(es[i]['element_number']):
                self.__atoms[n].symbol = elements[i]
                n += 1

    def read_poscar(self, filename):
        f = open(filename)
        self.__comment = f.readline().rstrip()
        # setup lattice constant
        lattice_constant = float(f.readline().split()[0])
        vectors = []

        # setup lattice vector
        for i in range(3):
            l = f.readline().split()
            tmp_vec = Vector(float(l[0]), float(l[1]), float(l[2]) )
            vectors.append(tmp_vec)
        self.__lattice = Lattice(vectors[0], vectors[1], vectors[2], lattice_constant)

        # setup number of element type
        tmp_atom_types = []
        tmp_atom_numbers = []
        l = f.readline().split()

        # for VASP5 POSCAR
        if l[0].isalpha():
            tmp_atom_types = l
            l = f.readline().split()

        total_atom_number = 0
        for n in l:
            tmp_atom_numbers.append(int(n) )

        # setup selective mode and coordinate type
        l = f.readline().rstrip()
        if l.split()[0][0].upper() == "S":
            # setup coordinate type
            self.__selective_mode = l
            l = f.readline().rstrip()
            # cartesian coordinates
            if l.split()[0][0].upper() == 'C':
                self.__coordinate_type = l
            # cartesian coordinates
            elif l.split()[0][0].upper() == 'K':
                self.__coordinate_type = l
            # direct/fractional coordinates)
            elif l.split()[0][0].upper() == 'D':
                self.__coordinate_type = l
            else:
                pass

        # setup atom coordinate
        for i in range(len(tmp_atom_numbers) ):
            if len(tmp_atom_types) is not 0:
                symbol = tmp_atom_types[i]
            else:
                symbol = str(i)
                pass
            for j in range(tmp_atom_numbers[i]):
                l = f.readline().split()
                if len(l) == 3:
                    a = Atom(symbol ,float(l[0]) ,float(l[1]), float(l[2]) )
                    check_element_by_periodic_table(a, 'symbol')
                elif len(l) == 6:
                    a = Atom(symbol ,float(l[0]) ,float(l[1]), float(l[2]), l[3], l[4], l[5])
                    check_element_by_periodic_table(a, 'symbol')
                else:
                    print('error format')
                self.add_atom(a)

        f.close()

    def write_poscar(self, filename=None):
        comment = ''
        number_of_atom = ''
        atom_of_type = ''

        for es in self._check__elements():
            atom_of_type += ' ' + str(es['element_type'])
            number_of_atom += ' ' + str(es['element_number'])

        l = self.__lattice.vectors
        v11, v12, v13 = l[0].basis
        v21, v22, v23 = l[1].basis
        v31, v32, v33 = l[2].basis
        format1 = '''%s
%.10f
%+14.10f %+14.10f %+14.10f
%+14.10f %+14.10f %+14.10f
%+14.10f %+14.10f %+14.10f
'''
        output1 = format1 % (self.comment + comment, self.__lattice.constant,
                             v11, v12, v13,
                             v21, v22, v23,
                             v31, v32, v33)


        if atom_of_type.split()[0].istitle():
            output2 = atom_of_type + "\n"
            output2 += number_of_atom + "\n"
        else:
            output2 = number_of_atom + "\n"
#        output2 = number_of_atom + "\n"
#        output += number_of_atom + "\n"

        output3 = ""
        if self.selective_mode is None:
            output3 = "%s\n" % (self.coordinate_type)
            for a in self.__atoms:
                output3 += "%+14.10f %+14.10f %+14.10f\n" % a.coordinate
        else:
            output3 = "%s\n%s\n" % (self.selective_mode, self.coordinate_type)
            for a in self.__atoms:
                output3 += "%+14.10f %+14.10f %+14.10f" % a.coordinate + "  %s %s %s\n" % a.dynamic
        
        # setup output
        if filename is None:
            print(output1 + output2 + output3)
        else:
            f = open(filename, 'w')
            f.write(output1 + output2 + output3)
            f.close()
#        self.listAtom()

    def line_scan(self, distance, nstep, ref_indexes, mot_indexes, grp_indexes):
        """ line scan
            distance:              {float}
            nstep:                 {int}
            reference atom index:  {int array}
            motion atom index:     {int array}
            group atom index:      {int array}
        """
        ref_atm = self.__atoms[ref_indexes[0]]
        mot_atm = self.__atoms[mot_indexes[0]]
        x1, y1, z1 = ref_atm.coordinate
        x2, y2, z2 = mot_atm.coordinate
        vec = Vector(x2 - x1, y2 - y1, z2 - z1)
        vec = vec.normalized()
    
        poscars = []
        tmp_indexes = mot_indexes + grp_indexes
    
        for i in range(nstep + 1):
            v = vec * i * distance
            tmp_poscar = POSCAR(comment = self.comment, select = self.selective_mode, coordinate = self.coordinate_type)
            tmp_poscar.__lattice = copy.deepcopy(self.lattice)
            tmp_poscar.__atoms = copy.deepcopy(self.__atoms)
            x, y, z = v.basis

            for j in tmp_indexes:
                tmp_x, tmp_y, tmp_z = self.__atoms[j].coordinate
                tmp_poscar.set_atom_coordinate(j, tmp_x + x, tmp_y + y, tmp_z + z)
            poscars.append(tmp_poscar)
        return poscars

    def angle_scan(self, angle, nstep, ref_indexes, mot_indexes, grp_indexes):
        """ angle scan
            angle:                 {float}
            nstep:                 {int}
            reference atom index:  {int array}
            motion atom index:     {int array}
            group atom index:      {int array}
        """
    # ref atm, fix/basic atm, mot atm
        ref_atm = self.__atoms[ref_indexes[0]]
        bas_atm = self.__atoms[ref_indexes[1]]
        mot_atm = self.__atoms[mot_indexes[0]]
    
        x1, y1, z1 = ref_atm.coordinate
        x2, y2, z2 = bas_atm.coordinate
        x3, y3, z3 = mot_atm.coordinate
        vec1 = Vector(x1 - x2, y1 - y2, z1 - z2)
        vec2 = Vector(x3 - x2, y3 - y2, z3 - z2)
        normal_vector = vec1.cross(vec2)
    
        poscars = []
        tmp_indexes = mot_indexes + grp_indexes
    
        for i in range(nstep):
            a = i * angle
            tmp_poscar = POSCAR(comment = self.comment, select = self.selective_mode, coordinate = self.coordinate_type)
            tmp_poscar.__lattice = copy.deepcopy(self.lattice)
            tmp_poscar.__atoms = copy.deepcopy(self.__atoms)
    
            for j in tmp_indexes:
                tmp_x, tmp_y, tmp_z = self.__atoms[j].coordinate
                tmp_v = Vector(tmp_x - x2, tmp_y - y2, tmp_z - z2)
                tmp_v = tmp_v.rotate(normal_vector, a)
                tmp_x, tmp_y, tmp_z = tmp_v.basis
                tmp_poscar.set_atom_coordinate(j, tmp_x + x2, tmp_y + y2, tmp_z + z2)
            poscars.append(tmp_poscar)
        return poscars

    def dihedral_scan(self, angle, nstep, ref_indexes, mot_indexes, grp_indexes):
        """ dihedral scan
            angle:                 {float}
            nstep:                 {int}
            reference atom index:  {int array}
            motion atom index:     {int array}
            group atom index:      {int array}
        """
    # ref atm, fix/basic atm, mot atm
        ref1_atm = self.__atoms[ref_indexes[0]]
        ref2_atm = self.__atoms[ref_indexes[1]]
        ref3_atm = self.__atoms[ref_indexes[2]]
        mot_atm  = self.__atoms[mot_indexes[0]]
    
        x1, y1, z1 = ref1_atm.coordinate
        x2, y2, z2 = ref2_atm.coordinate
        x3, y3, z3 = ref3_atm.coordinate
        x4, y4, z4 = mot_atm.coordinate
    
        vec1 = Vector(x1 - x2, y1 - y2, z1 - z2)
        vec2 = Vector(x3 - x2, y3 - y2, z3 - z2)
        vec3 = Vector(x2 - x3, y2 - y3, z2 - z3)
        vec4 = Vector(x4 - x3, y4 - y3, z4 - z3)
    
        normal_vec1 = vec1.cross(vec2)
        normal_vec2 = vec3.cross(vec4)
        normal_vec3 = normal_vec1.cross(normal_vec2)
    
        poscars = []
        tmp_indexes = mot_indexes + grp_indexes
    
        for i in range(nstep):
            a = i * angle
            tmp_poscar = POSCAR(comment = self.comment, select = self.selective_mode, coordinate = self.coordinate_type)
            tmp_poscar.__lattice = copy.deepcopy(self.lattice)
            tmp_poscar.__atoms = copy.deepcopy(self.__atoms)
    
            for j in tmp_indexes:
                tmp_x, tmp_y, tmp_z = self.__atoms[j].coordinate
                tmp_v = Vector(tmp_x - x3, tmp_y - y3, tmp_z - z3)
                tmp_v = tmp_v.rotate(normal_vec3, a)
                tmp_x, tmp_y, tmp_z = tmp_v.basis
                tmp_poscar.set_atom_coordinate(j, tmp_x + x3, tmp_y + y3, tmp_z + z3)
            poscars.append(tmp_poscar)
        return poscars

    def add_poscar_coordinate(self, atoms):
        pass

    def direct_to_cartesian(self):
        """ direct coordinate convert to cartesian coordinate
            ref: http://en.wikipedia.org/wiki/Fractional_coordinates#Conversion_to_cartesian_coordinates
        """
        if self.coordinate_type[0].upper() == 'D':
            self.coordinate_type = 'Cartesian'
            l = self.lattice.vectors
            lc = self.lattice.constant
            v1 = l[0]
            v2 = l[1]
            v3 = l[2]
            alpha = v2.angle(v3)
            beta = v1.angle(v3)
            gamma = v1.angle(v2)
            v1_len = v1.length()
            v2_len = v2.length()
            v3_len = v3.length()
            cos_alpha = math.cos(math.radians(alpha))
            cos_beta = math.cos(math.radians(beta))
            cos_gamma = math.cos(math.radians(gamma))
            sin_alpha = math.sin(math.radians(alpha))
            sin_beta = math.sin(math.radians(beta))
            sin_gamma = math.sin(math.radians(gamma))
            volume = math.sqrt(1 - cos_alpha * cos_alpha - cos_beta * cos_beta - cos_gamma * cos_gamma + 2 * cos_alpha * cos_beta * cos_gamma)

            for a in self.__atoms:
                x, y, z = a.coordinate
                tmp_x = v1_len * x + v2_len * cos_gamma * y + v3_len * cos_beta * z
                tmp_y = v2_len * sin_gamma * y + v3_len * (cos_alpha - cos_beta * cos_gamma) / sin_gamma * z
                tmp_z = v3_len * volume/sin_gamma*z
                a.set_coordinate(tmp_x, tmp_y, tmp_z)

    def cartesian_to_direct(self):
        """ cartesian coordinate convert to direct coordinate
            ref: http://en.wikipedia.org/wiki/Fractional_coordinates#Conversion_from_cartesian_coordinates
        """
        if self.get_coordinate_type()[0].upper() == 'C':
            self.set_coordinate_type('Direct')
            l = self.get_lattice().vectors
            lc = self.get_lattice().constant
            v1 = l[0]
            v2 = l[1]
            v3 = l[2]
            alpha = v2.angle(v3)
            beta = v1.angle(v3)
            gamma = v1.angle(v2)
            v1_len = v1.length()
            v2_len = v2.length()
            v3_len = v3.length()
            cos_alpha = math.cos(math.radians(alpha))
            cos_beta = math.cos(math.radians(beta))
            cos_gamma = math.cos(math.radians(gamma))
            sin_alpha = math.sin(math.radians(alpha))
            sin_beta = math.sin(math.radians(beta))
            sin_gamma = math.sin(math.radians(gamma))
            volume = math.sqrt(1 - cos_alpha * cos_alpha - cos_beta * cos_beta - cos_gamma * cos_gamma + 2 * cos_alpha * cos_beta * cos_gamma)

            for a in self.__atoms:
                x, y, z = a.coordinate
                tmp_x = 1 / v1_len * x - cos_gamma / v1_len / sin_gamma * y + (cos_alpha * cos_gamma - cos_beta) / v1_len / volume / sin_gamma * z
                tmp_y = 1 / v2_len / sin_gamma * y + (cos_beta * cos_gamma - cos_alpha) / v2_len / volume / sin_gamma * z
                tmp_z = sin_gamma / v3_len/volume * z
                a.set_coordinate(tmp_x, tmp_y, tmp_z)
            

class OUTCAR(object):
    def __init__(self, filename='OUTCAR'):
        self.__elements = []
        self.__lattices = []
        self.__dynamic_matrixes = []
        self.read_outcar(filename)

    def read_outcar(self, filename):
        f = open(filename)
        l = ' '
        # global reg pattern
        re_space = re.compile('^\s+?$')
        re_division_line = re.compile('--------------------------------------------------------------------------------------------------------')

        # INCAR section pattern
        re_incar = re.compile('^\s+?INCAR:')
        re_potcar = re.compile('^\s+?POTCAR:\s+?(\w+)\s+?(\w+)\s+?(\w+)')
        re_title = re.compile('^\s+?TITEL\s+?=\s+?(\w+)\s+?(\w+)\s+?(\w+)')

        # finite differences section pattern
        re_finite_differences = re.compile('^\s+?finite differences')
        re_direct_lattice_vectors_reciprocal_lattice_vertors = re.compile('^\s+?direct\s+?lattice\s+?vectors\s+?reciprocal\s+?lattice\s+?vectors')

        # ios position section pattern
        re_ion_position = re.compile('^\s+?ion  position')
        re_lattice_vectors = re.compile('^\s+?Lattice vectors:')
        re_vector_form = re.compile('\s+?(\w+) = \(\s+?(-?(\d+)\.(\d+)),\s+?(-?(\d+)\.(\d+)),\s+?(-?(\d+)\.(\d+))\)')

        # Dimension of arrays section pattern
        re_ions_per_type = re.compile('^\s+?ions per type =')

        re_position = re.compile('^\s+?position of ions in cartesian coordinates')
        re_dyn_mat = re.compile('Eigenvectors and eigenvalues of the dynamical matrix')
        re_finite = re.compile('Finite differences POTIM')

        total_atom_number = 0
        tmp_elements = []
        while l:
            l = f.readline()

            # Get potcar
            if re_title.search(l):
                r = re_title.match(l)
                e1, e2, e3 = r.groups()
                element = {'potential': e1, 'element': e2, 'date': e3}
                self.__elements.append(element)

            #  finite differences section
            if re_finite_differences.search(l):
                l = f.readline()
                tmp_flag = False
                while not (tmp_flag and re_division_line.search(l) ):
                    if re_division_line.search(l):
                        tmp_flag = True
                    if re_direct_lattice_vectors_reciprocal_lattice_vertors.search(l):
                        l = f.readline()
                        tmp_array = l.split()
                        tmp_vector1 = Vector(float(tmp_array[0]), float(tmp_array[1]), float(tmp_array[2]))

                        l = f.readline()
                        tmp_array = l.split()
                        tmp_vector2 = Vector(float(tmp_array[0]), float(tmp_array[1]), float(tmp_array[2]))

                        l = f.readline()
                        tmp_array = l.split()
                        tmp_vector3 = Vector(float(tmp_array[0]), float(tmp_array[1]), float(tmp_array[2]) )

                        self.__lattices.append(Lattice(tmp_vector1, tmp_vector2, tmp_vector3) )
                    l = f.readline()

            # Get atom type number
            if re_ions_per_type.search(l):
                 tmp_array = l.split()
                 for i in range(len(tmp_array) - 4):
                     self.__elements[i]['number'] = int(tmp_array[i + 4])
                     for j in range(int(tmp_array[i+4])):
                         tmp_elements.append(self.__elements[i]['element'])

            # Get dynamical matrix
            if re_dyn_mat.search(l):
                l = f.readline()

                while not (re_dyn_mat.search(l) or re_finite.search(l) ):
                    tmp_array = l.split()
                    # Get image freq
                    if len(tmp_array) == 10:
                        i = 0
                        freq = {"THz": float(tmp_array[2]) * -1.0,
                                "2PiTHz": float(tmp_array[4]) * 1.0,
                                "cm-1": float(tmp_array[6]) * -1.0,
                                "meV": float(tmp_array[8]) * -1.0}
                        atoms = []
                        l = f.readline()
                        tmp_array = l.split()

                        while len(tmp_array) == 6:
                            if tmp_array[0] == 'X':
                                l = f.readline()
                                tmp_array = l.split()
                            else:
                                tmp_atom = Atom(element_symbol= tmp_elements[i],
                                                x_coordinate = float(tmp_array[0]), y_coordinate = float(tmp_array[1]), z_coordinate = float(tmp_array[2]),
                                                x_displace = float(tmp_array[3]), y_displace = float(tmp_array[4]), z_displace = float(tmp_array[3]) )
                                l = f.readline()
                                tmp_array = l.split()
                                check_element_by_periodic_table(tmp_atom)
                                atoms.append(tmp_atom)
                                i += 1
                        self.__dynamic_matrixes.append({"freq": freq, "atoms": atoms})

                    # Get real freq
                    elif len(tmp_array) == 11:
                        i = 0
                        freq = {"THz": float(tmp_array[3]),
                                "2PiTHz": float(tmp_array[5]),
                                "cm-1": float(tmp_array[7]),
                                "meV": float(tmp_array[9])}
                        atoms = []
                        l = f.readline()
                        tmp_array = l.split()

                        while len(tmp_array) == 6:
                            if tmp_array[0] == 'X':
                                l = f.readline()
                                tmp_array = l.split()
                            else:
                                tmp_atom = Atom(element_symbol = tmp_elements[i],
                                                x_coordinate = float(tmp_array[0]), y_coordinate = float(tmp_array[1]), z_coordinate = float(tmp_array[2]),
                                                x_displace = float(tmp_array[3]), y_displace = float(tmp_array[4]), z_displace = float(tmp_array[3]) )
                                l = f.readline()
                                tmp_array = l.split()
                                check_element_by_periodic_table(tmp_atom)
                                atoms.append(tmp_atom)
                                i += 1
                        self.__dynamic_matrixes.append({"freq": freq, "atoms": atoms})
 
                    l = f.readline()
                
        f.close()
        self.__dynamic_matrixes.reverse()

    def write_log(self, filename=None):
        out0orientation = """GradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGrad
 Number of steps in this run=   2 maximum allowed number of steps=   2.
 GradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGrad
                          Standard orientation:
 ---------------------------------------------------------------------
 Center     Atomic     Atomic              Coordinates (Angstroms)
 Number     Number      Type              X           Y           Z
 ---------------------------------------------------------------------\n"""
        out4orientation = """
  ***** Axes restored to original set *****
  -------------------------------------------------------------------
  Center     Atomic                   Forces (Hartrees/Bohr)
  Number     Number              X              Y              Z
  -------------------------------------------------------------------
          25          25
          26          25
          27          25
          28          25
          29          25
          30          25
          31          25
          32          25
          33          25
          34          25
          35          25
          36          25
          37          25
          38           8
          39           8
          40           8
          41           8
          42           8
          43           8
          44           8
          45           8
          46           8
          47           8
          48           8
          49           8
          50           8
          51           8
          52           8
          53           8
          54           8
          55           8
          56           8
          57           8
          58           8
          59           8
          60           8
          61           8
          62           8
          63           8
          64           8
          65           8
          66           8
          67           8
          68           8
          69           8
          70           8
          71           8
          72           8
          73           8
          74           8
          75           8
          76           8
          77           8
          78           8
          79           8
          80           8
  -------------------------------------------------------------------

 GradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGrad
 Step number   1 out of a maximum of   2
 GradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGrad"""
        out0coordinate = '     %3d       %3d        %3d       %10.6f   %10.6f %10.6f\n'
        out0dash = '---------------------------------------------------------------------\n'
        out0frequency = """ Harmonic frequencies (cm**-1), IR intensities (KM/Mole), Raman scattering
 activities (A**4/AMU), depolarization ratios for plane and unpolarized
 incident light, reduced masses (AMU), force constants (mDyne/A),
 and normal coordinates:\n"""
        out1number = '                   %3d\n'
        out2number = '                   %3d                    %3d\n'
        out3number = '                   %3d                    %3d                         %3d\n'
        out1frequency = ' Frequencies --   %10.4f\n'
        out2frequency = ' Frequencies --   %10.4f             %10.4f\n'
        out3frequency = ' Frequencies --   %10.4f             %10.4f             %10.4f\n'
        out1title = ' Atom AN      X      Y      Z\n'
        out2title = ' Atom AN      X      Y      Z        X      Y      Z\n'
        out3title = ' Atom AN      X      Y      Z        X      Y      Z        X      Y      Z\n'
        out1col = ' %3d %3d   %6.2f %6.2f %6.2f\n'
        out2col = ' %3d %3d   %6.2f %6.2f %6.2f   %6.2f %6.2f %6.2f\n'
        out3col = ' %3d %3d   %6.2f %6.2f %6.2f   %6.2f %6.2f %6.2f   %6.2f %6.2f %6.2f\n'
        sentances = [{'orientation': out0orientation, 'coordinates': out0coordinate, 'dash': out0dash, 'frequency': out0frequency},
                     {'number': out1number, 'frequency': out1frequency, 'title': out1title, 'col': out1col},
                     {'number': out2number, 'frequency': out2frequency, 'title': out2title, 'col': out2col},
                     {'number': out3number, 'frequency': out3frequency, 'title': out3title, 'col': out3col},
                     {'orientation': out4orientation},]

        number_of_freq = len(self.__dynamic_matrixes)
        quotient = number_of_freq // 3
        remainder = number_of_freq % 3
        number_of_atoms = len(self.__dynamic_matrixes[0]['atoms'])
        out = sentances[0]['orientation']

        for i in range(number_of_atoms):
            x, y, z = self.__dynamic_matrixes[0]['atoms'][i].coordinate
            atomic_number = self.__dynamic_matrixes[0]['atoms'][i].atomic_number
            out += sentances[0]['coordinates'] %(i + 1, atomic_number, 0, x, y, z)

        for i in range(3):
            x, y, z = self.__lattices[0].vectors[i].basis
            out += sentances[0]['coordinates'] %(number_of_atoms + i + 1, -2, 0, x, y, z)

        out += sentances[0]['dash']
 
        out += sentances[0]['frequency']
        for i in range(1, quotient + 1):
            out += sentances[3]['number'] %(3 * i - 2, 3 * i - 1, 3 * i)
            out += sentances[3]['frequency'] %(self.__dynamic_matrixes[3*i-3]['freq']['cm-1'], self.__dynamic_matrixes[3 * i - 2]['freq']['cm-1'], self.__dynamic_matrixes[3 * i - 1]['freq']['cm-1'])
            out += sentances[3]['title']
            for j in range(number_of_atoms):
                a11, a12, a13 = self.__dynamic_matrixes[3*i-3]['atoms'][j].displace
                a21, a22, a23 = self.__dynamic_matrixes[3*i-2]['atoms'][j].displace
                a31, a32, a33 = self.__dynamic_matrixes[3*i-1]['atoms'][j].displace
                out += sentances[3]['col'] %(j+1, 1, a11, a12, a13, a21, a22, a23, a31, a32, a33)

        if remainder == 1:
            out += sentances[1]['number'] %(number_of_freq)
            out += sentances[1]['frequency'] %(self.__dynamic_matrixes[number_of_freq - 1]['freq']['cm-1'])
            out += sentances[1]['title']
            for j in range(number_of_atoms):
                a11, a12, a13 = self.__dynamic_matrixes[number_of_freq-1]['atoms'][j].displace
                out += sentances[1]['col'] %(j+1, 1, a11, a12, a13)
        elif remainder == 2:
            out += sentances[2]['number'] %(number_of_freq-1, number_of_freq)
            out += sentances[2]['frequency'] %(self.__dynamic_matrixes[number_of_freq - 2]['freq']['cm-1'], self.__dynamic_matrixes[number_of_freq - 2]['freq']['cm-1'])
            out = out + sentances[2]['title']
            for j in range(number_of_atoms):
                a11, a12, a13 = self.__dynamic_matrixes[number_of_freq - 2]['atoms'][j].displace
                a21, a22, a23 = self.__dynamic_matrixes[number_of_freq - 1]['atoms'][j].displace
                out += sentances[2]['col'] %(j+1, 1, a11, a12, a13, a21, a22, a23)

        out += sentances[4]['orientation']

        # setup output
        if filename is None:
            print(out)
        else:
            f = open(filename, 'w')
            f.write(out)
            f.close()


if __name__ == "__main__":
    p = POSCAR('c1')
    p.cartesian_to_direct()
    p.write_poscar('d1')

    p = POSCAR('d1')
    p.direct_to_cartesian()
    p.write_poscar('c2')
