import unittest

from vtool3.infra import *


class InfraTest(unittest.TestCase):
    def test_vector(self):
        v1 = Vector(1, 0, 0)
        self.assertEqual(1, v1.length())

        v2 = Vector(0, 1, 0)
        self.assertEqual(90.0, v1.angle(v2))
        self.assertEqual(0, v1.dot(v2))

    def test_lattice_init(self):
        l = Lattice(1, 2, 3)
        self.assertEqual([1, 2, 3], l.vectors)

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
        a = Atom('H', 1, 2, 3)
        self.assertEqual('H', a.symbol)
        self.assertEquals((1, 2, 3), a.coordinate)
