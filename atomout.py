# python utility with input and output functions

import numpy as np
from numpy import cos, sin, arccos, arctan2
from numpy import linalg as la
from StringIO import StringIO
import sys

def parse_poscar(FILE):

	lines = np.genfromtxt(FILE, delimiter="\n", dtype=str)
	
	latconst = np.genfromtxt(StringIO(lines[1]))
	
	a1 = np.genfromtxt(StringIO(lines[2]))
	a2 = np.genfromtxt(StringIO(lines[3]))
	a3 = np.genfromtxt(StringIO(lines[4]))
	a = np.vstack((a1, a2, a3))

	nats = np.genfromtxt(StringIO(lines[5])).astype(int)
	form = lines[6]
	b = np.genfromtxt(StringIO(lines[7]))
	for i in range(8, 8 + nats - 1):
		badd = np.genfromtxt(StringIO(lines[i]))
		b = np.vstack((b,badd))	
	if nats == 1:
		b = np.array([b])

	return latconst, nats, form, a, b

def output_xyz(atoms, line):

	print np.size(atoms, axis=0), "atoms"
	print line

	for at in atoms:
		print '%1d\t %8f\t %8f\t %8f' % (0, at[0], at[1], at[2])
	
	return

def output_xyz_alloy(atoms, line, dec):

	print np.size(atoms, axis=0), "atoms"
	print line

	counter = 0
	for at in atoms:
		print '%1d\t %8f\t %8f\t %8f' % (dec[counter], at[0], at[1], at[2])
		counter += 1
	return

def output_xyz_color(atoms, mp, rel, boxbound):

	print np.size(atoms, axis=0), "atoms"
	print ""

	for atom in atoms:
		if atom[2] < mp + rel and atom[2] >= mp - rel:		#make sure that ALL gamma surfaces are able to relax
			typ = 1
		elif atom[2] >= boxbound - rel:
			typ = 1
		elif atom[2] < rel:
			typ = 1
		else:
			typ = 0
		print '%1s\t %12f\t %12f\t %12f' % (typ, atom[0], atom[1], atom[2])
	return

def output_car(a, lat, atoms, line):

	print line
	print lat
	for aa in a:
		print '%12f\t %12f\t %12f\t' % (aa[0], aa[1], aa[2])
	print atoms.shape[0]

	boxbound = a[2][2]

	#print "Selective dynamics"
	print "Cartesian"

	for atom in atoms:
		print '%12f\t %12f\t %12f\t' % (atom[0], atom[1], atom[2])
	
	return

def output_car_alloy(a, lat, atoms, line, dec, form):

	print line
	print lat
	for aa in a:
		print '%12f\t %12f\t %12f\t' % (aa[0], aa[1], aa[2])

	n0 = 0; n1 = 0
	for i in range(0, atoms.shape[0]):
		if dec[i] == 0:
			n0 += 1
		else:
			n1 += 1

	print n0, n1

	boxbound = a[2][2]

	#print "Selective dynamics"
	print form 

	# make sure atoms are printed in order of decoration
	counter = 0
	for atom in atoms:
		if dec[counter] == 0:
			print '%12f\t %12f\t %12f\t' % (atom[0], atom[1], atom[2])

	counter = 0
	for atom in atoms:
		if dec[counter] == 1:
			print '%12f\t %12f\t %12f\t' % (atom[0], atom[1], atom[2])
	
	return
	
def output_car_direct(a, lat, atoms, line):

	print line
	print lat
	for aa in a:
		print '%12f\t %12f\t %12f\t' % (aa[0], aa[1], aa[2])
	print atoms.shape[0]

	boxbound = a[2][2]

	#print "Selective dynamics"
	print "Direct"

	for atom in atoms:
		print '%12f\t %12f\t %12f\t' % (atom[0], atom[1], atom[2])
	
	return

def output_car_sd(a, atoms, line, flags):

	print line
	print LAT
	for aa in a:
		print '%12f\t %12f\t %12f\t' % (aa[0], aa[1], aa[2])
	print atoms.shape[0]

	boxbound = a[2][2]

	print "Selective dynamics"
	print "Cartesian"

	counter = 0
	for atom in atoms:
		flag = flags[counter]
		print '%12f\t %12f\t %12f\t %1s\t %1s\t %1s' % (atom[0], atom[1], atom[2], flag[0], flag[1], flag[2])
		counter += 1

	return

def output_lmp(box, atoms, nats, line):

	print line
	print '%4d %5s' % (nats, "atoms")
	print "1 atom types"


	if ( box[0][0] > 0. ):
		print '%8f %8f %3s %3s' % (0, box[0][0], "xlo", "xhi")
	else:
		print '%8f %8f %3s %3s' % (box[0][0], 0, "xlo", "xhi")
	if ( box[1][1] > 0. ):
		print '%8f %8f %3s %3s' % (0, box[1][1], "ylo", "yhi")
	else:
		print '%8f %8f %3s %3s' % (box[1][1], 0, "ylo", "yhi")
	if ( box[2][2] > 0. ):
		print '%8f %8f %3s %3s' % (0, box[2][2], "zlo", "zhi")
	else:
		print '%8f %8f %3s %3s' % (box[2][2], 0, "zlo", "zhi")

	print "Atoms"
	print ""
	counter=1
	for atom in atoms:
		print '%4d\t %1d\t %8f\t %8f\t %8f' % (counter, 1, atom[0], atom[1], atom[2])
		counter += 1
	return 

def output_lmp_alloy(box, atoms, nats, line, dec):

	ntypes = np.unique(dec).shape[0]

	print line
	print '%4d %5s' % (nats, "atoms")
	print '%2d atom types' % ntypes
	# need to make sure the lattice vectors form lower triangular matrix (row-wise ordering)

	# convert atoms into direct coordinates

	CtoD = np.linalg.inv(np.transpose(np.array([box[0,:], box[1,:], box[2,:]])))
	
	for i in range(len(atoms)):
		atoms[i,:] = np.around(np.dot(CtoD, atoms[i,:]), decimals=16)
		
	# first rotate a into x-y plane
	if box[0,1] != 0.0:
		tet = np.arctan(-box[0,2]/box[0,1])
	elif box[0,1] == 0.0 and box[0,2] != 0.0:
		sign = -box[0,2]/np.fabs(box[0,2])
		tet = sign*np.pi/2
	else:
		tet = 0.0

	if np.round(tet, decimals=8) != 0.:
		rot = np.array([[1.,	0.,	0.],
				[0.,	np.cos(tet), -np.sin(tet)],
				[0.,	np.sin(tet),  np.cos(tet)]])
		for i in range(3):
			box[i,:]   = np.around(np.dot(rot, box[i,:]), decimals=16)
	#	for i in range(len(atoms)):
	#		atoms[i,:] = np.around(np.dot(rot, atoms[i,:]), decimals=8)
	
	# next rotate a into x
	if box[0,0] != 0.0:
		tet = np.arctan(-box[0,1]/box[0,0])
	elif box[0,0] == 0.0 and box[0,1] != 0.0:
		sign = -box[0,1]/np.fabs(box[0,1])
		tet = sign*np.pi/2
	else:
		tet = 0.0

	if np.round(tet, decimals=8) != 0.:
		rot = np.array([[np.cos(tet),	-np.sin(tet),	0.],
				[np.sin(tet),	 np.cos(tet),	0.],
				[0.,		 0.,		1.]])
		for i in range(3):
			box[i,:]   = np.around(np.dot(rot, box[i,:]), decimals=16)
	#	for i in range(len(atoms)):
	#		atoms[i,:] = np.around(np.dot(rot, atoms[i,:]), decimals=8)
#	print box
	# lastly rotate around the new x axis such that b lies in the x-y plane
	if box[1,1] != 0.0:
		tet = np.arctan(-box[1,2]/box[1,1])
	elif box[1,1] == 0.0 and box[1,2] != 0.0:
		sign = -box[1,2]/np.fabs(box[1,2])
		tet = sign*np.pi/2
	else:
		tet = 0.0

	if np.round(tet, decimals=8) != 0.:
		rot = np.array([[1.,	0.,	0.],
				[0.,	np.cos(tet), -np.sin(tet)],
				[0.,	np.sin(tet),  np.cos(tet)]])
		for i in range(3):
			box[i,:]   = np.around(np.dot(rot, box[i,:]), decimals=8)
	#	for i in range(len(atoms)):
	#		atoms[i,:] = np.around(np.dot(rot, atoms[i,:]), decimals=8)

	# now do any necessary reflections
	if box[0,0] < 0.:

		for i in range(3):
			box[i,:] = np.dot(np.array([[-1., 0., 0.],
				[0., 1., 0.],
				[0., 0., 1.]]), box[i,:])

	if box[1,1] < 0.:

		for i in range(3):
			box[i,:] = np.dot(np.array([[1., 0., 0.],
				[0., -1., 0.],
				[0., 0., 1.]]), box[i,:])

	if box[2,2] < 0.:

		for i in range(3):
			box[i,:] = np.dot(np.array([[1., 0., 0.],
				[0., 1., 0.],
				[0., 0., -1.]]), box[i,:])

	# back to cartesian
	for i in range(len(atoms)):
		atoms[i,:] = atoms[i,0]*box[0,:] + atoms[i,1]*box[1,:] + atoms[i,2]*box[2,:]

	# now check if the box is triclinic
	av = np.array(box[0,:]); bv = np.array(box[1,:]); cv = np.array(box[2,:])
	a = np.linalg.norm(av); b = np.linalg.norm(bv); c = np.linalg.norm(cv)

	alpha = np.arccos(np.dot(bv, cv)/b/c) 
	beta  = np.arccos(np.dot(av, cv)/a/c) 
	gamma = np.arccos(np.dot(av, bv)/a/b) 

	if np.round(gamma, decimals=6) != 0.:
		xy = np.linalg.norm(box[1,:])*np.cos(gamma)
		ly = np.sqrt(b * b - xy * xy)
	else:
		xy = 0.
		ly = b
	if np.round(beta,  decimals=6) != 0.:
		xz = c * np.cos(beta)
	else:
		xz = 0.
	if np.round(alpha, decimals=6) != 0.:
		yz = (b * c * np.cos(alpha) - xy * xz)/ly
	if xz != 0. or yz != 0.:
		lz = np.sqrt(c * c - xz * xz - yz * yz)
	else:
		lz = c
	
	print '%8f %8f %3s %3s' % (min(a,0.) , max(a,0.),  "xlo", "xhi")
	print '%8f %8f %3s %3s' % (min(ly,0.), max(ly,0.), "ylo", "yhi")
	print '%8f %8f %3s %3s' % (min(lz,0.), max(lz,0.), "zlo", "zhi")

	if xy != 0. or xz != 0. or yz != 0:
		print '%8f %8f %8f %2s %2s %2s' % (xy, xz, yz, "xy", "xz", "yz")
		#print '%8f %3s' % (xy, "xy")
		#print '%8f %3s' % (xz, "xz")
		#print '%8f %3s' % (yz, "yz")

	print "Atoms"
	print ""
	counter=0
	for atom in atoms:
		typ = dec[counter] + 1	
		counter += 1
		print '%4d\t %1d\t %8f\t %8f\t %8f' % (counter, typ, atom[0], atom[1], atom[2])
	return
 
def output_pdb(a, atoms, nats):

	alpha = 180. * np.arccos( np.dot(a[1], a[2]) / ( np.linalg.norm(a[1]) * np.linalg.norm(a[2]) ) ) / np.pi
	beta  = 180. * np.arccos( np.dot(a[0], a[2]) / ( np.linalg.norm(a[0]) * np.linalg.norm(a[2]) ) ) / np.pi
	gamma = 180. * np.arccos( np.dot(a[0], a[1]) / ( np.linalg.norm(a[0]) * np.linalg.norm(a[1]) ) ) / np.pi

	print "%-6s%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f%11s%4d" % ("CRYST1", a[0][0], a[1][1], a[2][2], alpha, beta, gamma, "P 1", 1)
	count = 0
	symbol = "Ti"
	for atom in atoms:
		count += 1
        	print "ATOM  %5d %2s   MOL     1  %8.3f%8.3f%8.3f  1.00  0.00" % (count, symbol, atom[0], atom[1], atom[2])
	
	print "END"

def output_pdb_alloy(a, atoms, nats, dec):

	alpha = 180. * np.arccos( np.dot(a[1], a[2]) / ( np.linalg.norm(a[1]) * np.linalg.norm(a[2]) ) ) / np.pi
	beta  = 180. * np.arccos( np.dot(a[0], a[2]) / ( np.linalg.norm(a[0]) * np.linalg.norm(a[2]) ) ) / np.pi
	gamma = 180. * np.arccos( np.dot(a[0], a[1]) / ( np.linalg.norm(a[0]) * np.linalg.norm(a[1]) ) ) / np.pi

	chems = np.array(["Ti", "Nb"])

	print "%-6s%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f%11s%4d" % ("CRYST1", a[0][0], a[1][1], a[2][2], alpha, beta, gamma, "P 1", 1)
	count = 0
	for atom in atoms:
		symbol = chems[dec[count]]
		count += 1
        	print "ATOM  %5d %2s   MOL     1  %8.3f%8.3f%8.3f  1.00  0.00" % (count, symbol, atom[0], atom[1]*a[1][1]/a[0][0], atom[2]*a[2][2]/a[0][0])
	
	print "END"

def output_conf(box, atoms, nats):
	print "#N", nats, "1"
	print "## generated with gammaTilt.py"
	print '%2s %8f %8f %8f' % ("#X", box[0][0], box[0][1], box[0][2])
	print '%2s %8f %8f %8f' % ("#Y", box[1][0], box[1][1], box[1][2])
	print '%2s %8f %8f %8f' % ("#Z", box[2][0], box[2][1], box[2][2])
	print '%2s %8f'	% ("#E", 0.0)
	print '%2s %8f %8f %8f %8f %8f %8f' % ("#S", 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
	print '#F'

	
	for atom in atoms:
		print '%1d %8f %8f %8f %8f %8f %8f' % (0, atom[0], atom[1], atom[2], 0.0, 0.0, 0.0)

	return

def output_conf_alloy(box, atoms, nats, dec, line):

	ntypes = np.unique(dec).shape[0]

	print "#N", nats, ntypes
	print "##", line
	print '%2s %8f %8f %8f' % ("#X", box[0][0], box[0][1], box[0][2])
	print '%2s %8f %8f %8f' % ("#Y", box[1][0], box[1][1], box[1][2])
	print '%2s %8f %8f %8f' % ("#Z", box[2][0], box[2][1], box[2][2])
	print '%2s %8f'	% ("#E", 0.0)
	print '%2s %8f %8f %8f %8f %8f %8f' % ("#S", 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
	print '#F'

	boa = box[1][1]/box[0][0]; coa = box[2][2]/box[0][0]

	counter = 0	
	for atom in atoms:
		typ = dec[counter]

		print '%1d %8f %8f %8f %8f %8f %8f' % (typ, atom[0], atom[1]*boa, atom[2]*coa, 0.0, 0.0, 0.0)
		counter += 1

	return
def output_conf_alloy_cart(box, atoms, nats, dec, line):

	ntypes = np.unique(dec).shape[0]

	print "#N", nats, ntypes
	print "##", line
	print '%2s %8f %8f %8f' % ("#X", box[0][0], box[0][1], box[0][2])
	print '%2s %8f %8f %8f' % ("#Y", box[1][0], box[1][1], box[1][2])
	print '%2s %8f %8f %8f' % ("#Z", box[2][0], box[2][1], box[2][2])
	print '%2s %8f'	% ("#E", 0.0)
	print '%2s %8f %8f %8f %8f %8f %8f' % ("#S", 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
	print '#F'

	boa = box[1][1]/box[0][0]; coa = box[2][2]/box[0][0]

	counter = 0	
	for atom in atoms:
		typ = dec[counter]

		print '%1d %8f %8f %8f %8f %8f %8f' % (typ, atom[0], atom[1], atom[2], 0.0, 0.0, 0.0)
		counter += 1

	return
def output_tri(box, atoms, nats, XY, XZ, YZ, line):

	print line
	print '%4d %5s' % (nats, "atoms")
	print "1 atom types"

	if ( box[0][0] > 0. ):
		print '%8f %8f %3s %3s' % (0, box[0][0], "xlo", "xhi")
	else:
		print '%8f %8f %3s %3s' % (box[0][0], 0, "xlo", "xhi")
	if ( box[1][1] > 0. ):
		print '%8f %8f %3s %3s' % (0, box[1][1], "ylo", "yhi")
	else:
		print '%8f %8f %3s %3s' % (box[1][1], 0, "ylo", "yhi")
	if ( box[2][2] > 0. ):
		print '%8f %8f %3s %3s' % (0, box[2][2], "zlo", "zhi")
	else:
		print '%8f %8f %3s %3s' % (box[2][2], 0, "zlo", "zhi")

	lx = np.fabs(box[0][0])
	ly = np.fabs(box[1][1])
	lz = np.fabs(box[2][2])

	while XY > lx/2 or XY < -lx/2:
		if XY > lx/2:
			XY = XY - lx
		elif XY < -lx/2:
			XY = XY + lx

	while XZ > lx/2 or XZ < -lx/2:
		if XZ > lx/2:
			XZ = XZ - lx
		elif XZ < -lx/2:
			XZ = XZ + lx

	while YZ > ly/2 or YZ < -ly/2:
		if YZ > ly/2:
			YZ = YZ - ly
		elif YZ < -ly/2:
			YZ = YZ + ly

	print '%8f %8f %8f %2s %2s %2s' % (XY, XZ, YZ, "xy", "xz", "yz")

	print "Atoms"
	print ""
	counter=1
	for atom in atoms:
		print '%4d\t %1d\t %8f\t %8f\t %8f' % (counter, 1, atom[0], atom[1], atom[2])
		counter += 1
	return 

def output_atom(box, atoms, nats):
	print "ITEM: TIMESTEP"
	print "0"
	print "ITEM: NUMBER OF ATOMS"
	print nats
	print "ITEM: BOX BOUNDS pp pp pp"
	print '%8f %8f' % (0, box[0][0])
	print '%8f %8f' % (0, box[1][1])
	print '%8f %8f' % (0, box[2][2])
	print "ITEM: ATOMS id element x y z"
	counter = 1
	for atom in atoms:
		print '%4d %2s %8f %8f %8f' % (counter, 'C', atom[0], atom[1], atom[2])
		counter += 1	

