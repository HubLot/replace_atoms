#!/usr/bin/env python

from __future__ import print_function, with_statement
import argparse
import random
import sys

__author__ = "Jonathan Barnoud"

# File format description. The key is the name of the field, the value is a
# tuple from which the first element is the first (included) and last
# (excluded) indices of the field in the line, and the second element the type
# of the field content.
GRO_FIELDS = {
    "resid": ((0, 5), int),
    "resname": ((5, 10), str),
    "atom_name": ((10, 15), str),
    "atomid": ((15, 20), int),
    "x": ((20, 28), float),
    "y": ((28, 36), float),
    "z": ((36, 44), float),
}


class FormatError(Exception):
    """
    Exception raised when the file format is wrong.
    """
    pass


def stop_at_empty_line(iterator):
    """
    Yield all item of an iterator but stop when the item is an empty line.

    An empty line is a string which is empty when stripped.
    """
    for line in iterator:
        if line.strip() == "":
            return
        yield line


def read_gro(lines):
    """
    Read the atoms, the header, and the box description from a gro file.

    Atoms are represented as dictionaries.

    :Parameters:
        - lines: an iterator over atom lines from the gro file. The two header
                 lines and the bottom line describing the box have to be
                 included.

    :Returns:
        - title: the title of the system as written on line 1 of the file
        - atoms: a list of atom, each atom is stored as a dictionary
        - box: the box description as written on the last line

    :Raise:
        - FormatError: raised if the file format does not fit.
    """
    # "lines" might be a list and not a proper iterator
    lines = iter(lines)
    # The two first lines are a header
    title = next(lines)
    next(lines)  # This is the number of atoms, we do not care
    # Loop over the lines but act on the previous one. We are reading atoms and
    # we do not want to consider the last line (the box description) as an
    # atom.
    atoms = []
    prev_line = next(lines)
    for line in stop_at_empty_line(lines):
        try:
            atoms.append(dict(((key, convert(prev_line[begin:end].strip()))
                               for key, ((begin, end), convert)
                               in GRO_FIELDS.items())))
            prev_line = line
        except ValueError:
            raise FormatError
    box = prev_line
    return (title, atoms, box)


def write_gro(title, atoms, box):
    """
    Yield lines of a GRO file from a title, a list of atoms and a box

    :Parameters:
        - title: the title of the system
        - atoms: a list of atom, each atom is stored as a dictionary
        - box: the box description as written on the last line
    """
    yield title
    yield '{0}\n'.format(len(atoms))
    for atom in atoms:
        yield ('{resid:>5}{resname:<5}{atom_name:>5}{atomid:>5}'
                '{x:8.3f}{y:8.3f}{z:8.3f}\n').format(**atom)
    yield box


def select(atoms, key, value):
    """
    Select the atoms for wich the given key match the given value.

    Return the selection as a generator on the atom index in the provided atom
    list.

    :Parameters:
        - atoms: a list of atom dictionaries
        - key: the key on which to filter
        - value: the value we want to keep

    :Return:
        - a generator of atom indices relative to the provided atom list
    """
    return (i for i, atom in enumerate(atoms) if atom[key] == value)


def alter_atoms(atoms, indices, resname, atom_name):
    for i in indices:
        atoms[i]['resname'] = resname
        atoms[i]['atom_name'] = atom_name
    

def reorder(length, insertions, index):
    for i in xrange(0, index):
        if not i in insertions:
            yield i
    for i in insertions:
        yield i
    for i in xrange(index, length):
        if not i in insertions:
            yield i


def iter_order(atoms, indices):
    for i in indices:
        yield atoms[i]


def renumber(atoms):
    resid = 0
    prev_resid = 0
    for atomid, atom in enumerate(atoms, start=1):
        if atom['resid'] != prev_resid:
            resid += 1
            prev_resid = atom['resid']
        atom['resid'] = resid%100000
        atom['atomid'] = atomid%100000
        yield atom


def replace_atom(atoms, key, value, resname, atom_name, nsub):
    candidates = list(select(atoms, key, value))
    selection = random.sample(candidates, nsub)
    alter_atoms(atoms, selection, resname, atom_name)
    new_order = reorder(len(atoms), selection, min(candidates))
    return renumber(iter_order(atoms, new_order))

def get_args(argv):
    usage = "%(prog)s [options] < input > output"
    parser = argparse.ArgumentParser(description=__doc__, usage=usage)
    parser.add_argument("--natoms", "-n", type=int, required=True,
                        help="Number of atom to modify")
    parser.add_argument("--oresname", "-o", type=str, required=True,
                        help="Original residue name of the atoms to replace")
    parser.add_argument("--nresname", "-r", type=str, required=True,
                        help="New residue name of the replaced atoms")
    parser.add_argument("--natname", "-a", type=str, required=True,
                        help="New atom name of the replaced atoms")
    args = parser.parse_args(argv)
    return args


def cli_main():
    args = get_args(sys.argv[1:])
    title, atoms, box = read_gro(sys.stdin)
    modified_atoms = list(replace_atom(atoms, 'resname', args.oresname,
                                       args.nresname, args.natname,
                                       args.natoms))
    for line in write_gro(title, modified_atoms, box):
        print(line, end='')


if __name__ == "__main__":
    cli_main()
