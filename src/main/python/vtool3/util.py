#!/usr/bin/env python3

__author__ = ""
__date__ = "2016/09/05"
__version__ = "$Revision: 0.3$"


def gjf_to_poscar(gjf, poscar):
    """ Gaussain GJF convert VASP POSCAR
        gjf:    {GJF}
        poscar: {POSCAR}
    """
    poscar.set_lattice(gjf.lattice.vectors)
    for a in gjf.get_atoms():
        poscar.add_atom(a)


def poscar_to_gjf(poscar, gjf, elements=None):
    """ VASP POSCAR convert Gaussain GJF
        poscar:   {POSCAR}
        gjf:      {GJF}
        elements: {string array}
    """
    poscar.direct_to_cartesian()
    if not(elements is None):
        poscar.set_elements_type(elements)
    for a in poscar.get_atoms():
        gjf.add_atom(a)
    gjf.set_lattice(poscar.lattice.vectors)
#    for v in poscar.getLattice().getVectors():
#        tmp1, tmp2, tmp3 = v.getBasis()
#        a = Atom('Tv', tmp1, tmp2, tmp3)
#        gjf.addAtom(a)


if __name__ == "__main__":
    import sys
    import os

    pass
