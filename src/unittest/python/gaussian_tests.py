import unittest
from io import StringIO
from unittest.mock import patch

from vtool3.gaussian import *


class GaussianTest(unittest.TestCase):
    def test_gjf(self):
        pass
        atom = Atom('H', 0, 0, 0)
        gjf = GJF()
        gjf.add_atom(atom)
        out = '''# opt freq hf/3-21g

This line is comment

0 1
H           0.00000000       0.00000000       0.00000000

'''
        with patch('sys.stdout', new=StringIO()) as fake_out:
            gjf.write_gjf()
            self.assertEqual(out, fake_out.getvalue())

    def test_gjf_read_file(self):
        gjf_content = '''# opt freq hf/3-21g

This line is comment

0 1
H           0.00000000       0.00000000       0.00000000

'''
        gjf_file = 'tmp.gjf'
        with open(gjf_file, 'w') as f:
            f.write(gjf_content)

        gjf = GJF()
        gjf.read_gjf(gjf_file)


