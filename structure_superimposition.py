#!/usr/bin/env python
#./structure_superimposition.py 3ZCF.pdb 3O20.pdb A A 1,2,3 1,2,3

import sys

def get_ca_atoms(pdbfile, chain, rlist, atom = 'CA'):
	list_coor = []
	file_pdb = open(pdbfile)
	for line in file_pdb:
		if line[:4] != 'ATOM' : continue # Si se cumple continua, no hace falta else
		if line[21] != chain : continue # Posicion 22 es caracter - chain identifier
		if line[22:26].strip() not in rlist: continue # Posicion 23-26 es Integer - Residue sequence number y strip para quitar espacios
		if line[12:16].strip() != atom: continue # Posicion 13-16 es Atom - Atom name
		x = float(line[30:38]) # Posicion 31-38 es Real(8.3) - Orthogonal coordinates for X in Angstroms
		y = float(line[38:46]) # Posicion 39-46  es Real(8.3) - Orthogonal coordinates for Y in Angstroms
		z = float(line[46:54]) # Posicion 47-54  es Real(8.3) - Orthogonal coordinates for Z in Angstroms
		list_coor.append([x,y,z])
	return list_coor

if __name__ == '__main__':
	pdbfile1 = sys.argv[1]
	pdbfile2 = sys.argv[2]
	
	chain1 = sys.argv[3]
	chain2 = sys.argv[4]
	
	list1 = sys.argv[5].split(',')
	list2 = sys.argv[6].split(',')
	
	list_c1 = get_ca_atoms(pdbfile1, chain1, list1)
	list_c2 = get_ca_atoms(pdbfile2, chain2, list2)
	
	print 'Coordination 1: ', list_c1
	print 'Coordination 2: ', list_c2
