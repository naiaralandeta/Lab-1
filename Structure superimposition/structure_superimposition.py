#!/usr/bin/env python
#./structure_superimposition.py 3ZCF.pdb 3O20.pdb A A 1,2,3 1,2,3 // si a la ultima lista se le anyade 4, da error de tamanyo
# URL: http://biopython.org/DIST/docs/api/Bio.SVDSuperimposer.SVDSuperimposer-class.html

import sys
from Bio.SVDSuperimposer import SVDSuperimposer
import numpy as np

def get_ca_atoms(pdbfile, chain, rlist, atom = 'CA'): # Funcion para obtener todos los carbonos alfa
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
	
def get_rmsd(coord1, coord2): # Funcion para hacer la superimposicion de las estructuras
	if len(coord1) != len(coord2):
		print >> sys.stderr, 'ERROR: The set of coordinates have different size.' # Estandar error
		sys.exit(1)
	svd = SVDSuperimposer()
	svd.set(np.array(coord1), np.array(coord2))
	svd.run() # run the superimposition
	rmsd = svd.get_rms() # Da la RMSD de la superimposicion
	rotation, translation = svd.get_rotran() # Matrix rotation and translation
	print 'ROTATION:', rotation
	print 'TRANSLATION:', translation
	return rmsd

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

	rmsd = get_rmsd(list_c1, list_c2)
	
	print 'RMSD:', rmsd
