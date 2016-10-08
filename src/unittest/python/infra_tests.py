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

        v4 = Vector(1, 2, 3)
        self.assertEquals('1.000000, 2.000000, 3.000000', v4.__str__())
        self.assertEquals(Vector(4, 4, 4), v4.__add__(Vector(3, 2, 1)))
        self.assertEquals(Vector(0, 0, 0), v4.__sub__(Vector(1, 2, 3)))
        self.assertEquals(Vector(5, 10, 15), v4.__mul__(5))
        self.assertEquals(True, v4.__eq__(Vector(1, 2, 3)))

        v5 = Vector(10, 0, 0)
        self.assertEquals(v5, v5.rotate(Vector(1, 0, 0), 360.0))

    def test_lattice_init(self):
        l = Lattice(1, 2, 3)
        self.assertEqual([1, 2, 3], l.vectors)
        l.constant = 2
        self.assertEqual(2, l.constant)

        l1 = Lattice(4, 5, 6, 10)
        self.assertEquals([4, 5, 6], l1.vectors)
        self.assertEquals(10, l1.constant)
        l1.set_vectors(2, 2, 2)
        self.assertEquals([2, 2, 2], l1.vectors)
        l1.constant = 5
        self.assertEquals(5, l1.constant)

    def test_element(self):
        x = Element()
        self.assertEqual('X', x.symbol)
        self.assertEqual('Dummy', x.name)
        self.assertEqual(0, x.atomic_number)
        self.assertEqual(0.0, x.atomic_mass)

        x.symbol = 'Se'
        x.name = 'Selenium'
        x.atomic_number = 34
        x.atomic_mass = 78.971
        self.assertEqual('Se', x.symbol)
        self.assertEqual('Selenium', x.name)
        self.assertEqual(34, x.atomic_number)
        self.assertEqual(78.971, x.atomic_mass)
        self.assertEquals('Se, Selenium, 34, 78', x.__str__())

        y = Element('H', 'Hydrogen', 1, 1)
        y.copy_element(x)
        self.assertEquals('Se, Selenium, 34, 78', y.__str__())

    def test_check_element_by_periodic_table(self):
        h = Element('H')
        check_element_by_periodic_table(h)
        self.assertEqual('Hydrogen', h.name)
        self.assertEquals(1, h.atomic_number)
        self.assertEquals(1.00794, h.atomic_mass)

        he = Element(name='Helium')
        check_element_by_periodic_table(he, 'name')
        self.assertEquals('He', he.symbol)
        self.assertEquals(2, he.atomic_number)
        self.assertEquals(4.002602, he.atomic_mass)

        li = Element(number=3)
        check_element_by_periodic_table(li, 'number')
        self.assertEquals('Li', li.symbol)
        self.assertEquals('Lithium', li.name)
        self.assertEquals(6.941, li.atomic_mass)

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

        default = Atom()
        self.assertEquals('X', default.symbol)
        self.assertEquals((0, 0, 0), default.coordinate)
        self.assertEquals(('T', 'T', 'T'), default.dynamic)
        self.assertEquals((0, 0, 0), default.displace)

        default.symbol = 'C'
        default.set_coordinate(3.52, 4.67, 2.33)
        default.set_dynamic('F', 'T', 'F')
        default.set_displace(10.5223, 20.9889, 56.2837)
        self.assertEquals('C', default.symbol)
        self.assertEquals((3.52, 4.67, 2.33), default.coordinate)
        self.assertEquals(('F', 'T', 'F'), default.dynamic)
        self.assertEquals((10.5223, 20.9889, 56.2837), default.displace)
        self.assertEquals("('C', 3.52, 4.67, 2.33)", default.__repr__())

        copy_from_default = Atom()
        copy_from_default.copy_atom(default)
        self.assertEquals(True, copy_from_default.__eq__(default))
        self.assertEquals(Atom('C', 9.0, 9.0, 9.0), copy_from_default.add_coordinate(5.48, 4.33, 6.67))
        self.assertEquals(Atom('C', 3.0, 4.0, 2.0), copy_from_default.sub_coordinate(0.52, 0.67, 0.33))
        self.assertEquals(Atom('C', 10.56, 14.01, 6.99), copy_from_default.mul_coordinate(3))
        self.assertEquals(Atom('C', 1.76, 2.335, 1.165), copy_from_default.div_coordinate(2))

