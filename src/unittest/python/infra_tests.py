import unittest

from vtool3.infra import *


class InfraTest(unittest.TestCase):
    def test_vector(self):
        v1 = Vector(1, 0, 0)
        self.assertEqual((1, 0, 0), v1.basis)
        self.assertEqual(1, v1.length())
        self.assertEqual(Vector(1, 0, 0), v1.normalized())

        v2 = Vector(0, 1, 0)
        self.assertEqual(90.0, v1.angle(v2))
        self.assertEqual(0, v1.dot(v2))
        v3 = v1.cross(v2)
        self.assertEqual(Vector(0, 0, 1), v3)

        v2.set_basis(0, 0, 1)
        self.assertEqual(v3, v2)

    def test_lattice_init(self):
        l = Lattice(1, 2, 3)
        self.assertEqual([1, 2, 3], l.vectors)
        l.constant = 2
        self.assertEqual(2, l.constant)

    def test_element(self):
        x = Element()
        self.assertEqual('X', x.symbol)
        self.assertEqual('Dummy', x.name)
        self.assertEqual(0, x.atomic_number)
        self.assertEqual(0.0, x.atomic_mass)

    def test_check_element_by_periodic_table(self):
        h = Element('H')
        check_element_by_periodic_table(h)
        self.assertEqual('Hydrogen', h.name)

    def test_atom(self):
        a1 = Atom('H', 1, 2, 3)
        self.assertEqual('H', a1.symbol)
        self.assertEquals((1, 2, 3), a1.coordinate)
        self.assertEqual(('T', 'T', 'T'), a1.dynamic)
        self.assertEqual((0, 0, 0), a1.displace)

        self.assertEquals((1, 2, 3), a1.coordinate)
        self.assertEqual(Atom('H', 2, 2, 3), a1.add_coordinate(1, 0, 0))
        self.assertEqual(Atom('H', 0, 0, 0), a1.sub_coordinate(1, 2, 3))
        self.assertEqual(Atom('H', 0, 0, 0), a1.mul_coordinate(0))
        self.assertEqual(Atom('H', 1, 2, 3), a1.div_coordinate(1))

        print(a1)
        a1.show_atom()

        a2 = Atom()
        a2.copy_atom(a1)


