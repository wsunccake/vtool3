import unittest
from io import StringIO
from unittest.mock import patch

from vtool3.vasp import *

class VaspTest(unittest.TestCase):
    def test_poscar(self):

        # p = POSCAR('POSCAR')
        # p.direct_to_cartesian()
        # p.write_poscar()

        # p.cartesian_to_direct()
        # p.write_poscar()
        p1 = POSCAR()

        p1.set_comment('This is comment line.')
        self.assertEquals('This is comment line.', p1.comment)

        p1.set_lattice([1, 2, 3], 10)
        self.assertEquals([1, 2, 3], p1.lattice.vectors)
        self.assertEquals(10, p1.lattice.constant)

        p1.set_selective_mode('Selective')
        self.assertEquals('Selective', p1.selective_mode)

        p1.set_coordinate_type('Direct')
        self.assertEquals('Direct', p1.coordinate_type)

        h1 = Atom('H', 10, 20, 30, 'T', 'F', 'T', 100, 200, 300)
        h2 = Atom('H', 1, 2, 3, 'F', 'F', 'F', 10, 20, 30)
        f = Atom('F', 100, 200, 300, 'F', 'T', 'F', 1000, 2000, 3000)
        p1.add_atom(h1)
        p1.add_atom(h2)
        p1.add_atom(f)
        self.assertEquals("[('H', 10, 20, 30), ('H', 1, 2, 3), ('F', 100, 200, 300)]", str(p1.get_atoms()))

        p1.set_atom_element(1, 'h')
        self.assertEquals("[('H', 10, 20, 30), ('h', 1, 2, 3), ('F', 100, 200, 300)]", str(p1.get_atoms()))

        p1.set_atom_coordinate(1, 10.1, 20.2, 30.3)
        self.assertEquals("[('H', 10, 20, 30), ('h', 10.1, 20.2, 30.3), ('F', 100, 200, 300)]", str(p1.get_atoms()))

        p1.set_atom_dynamic(1, 'T', 'T', 'T')
        self.assertEquals(('T', 'T', 'T'), p1.get_atoms()[1].dynamic)

        n = Atom('N', 1.002, 2.001, 3.002, 'F', 'T', 'F', 11.1, 21.1, 31.1)
        p1.set_atom(1, n)
        self.assertEquals("[('H', 10, 20, 30), ('N', 1.002, 2.001, 3.002), ('F', 100, 200, 300)]", str(p1.get_atoms()))


        out = '''  H, +10.0000000000, +20.0000000000, +30.0000000000,    T,    F,    T
  N,  +1.0020000000,  +2.0010000000,  +3.0020000000,    F,    T,    F
  F, +100.0000000000, +200.0000000000, +300.0000000000,    F,    T,    F
'''
        with patch('sys.stdout', new=StringIO()) as fake_out:
            p1.list_atom()
            self.assertEqual(out, fake_out.getvalue())

