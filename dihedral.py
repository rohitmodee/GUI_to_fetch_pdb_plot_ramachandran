import numpy as np
import pdb
def dihedral(p):
	b = np.array([p[0]-p[1], p[1]-p[2], p[1]-p[3]])
	q1 = np.cross(b[0], b[1])
	q2 = np.cross(b[1], b[2])
	q1q2 = np.dot(q1, q2)
	mod_q1 = np.sqrt(sum([i**2 for i in q1]))
	mod_q2 = np.sqrt(sum([i**2 for i in q2]))
	theta = np.arccos(q1q2/(mod_q1*mod_q2))
	
	q3 = np.cross(q1, q2)
	sign = np.dot(q3, b[1])
	if sign < 0:
		theta = theta * -1
	else:
		theta = theta * 1
	return np.degrees(theta)
	
p0 = np.array([24.969, 13.428, 30.692]) # N
p1 = np.array([24.044, 12.661, 29.808]) # CA
p2 = np.array([22.785, 13.482, 29.543]) # C
p3 = np.array([21.951, 13.670, 30.431]) # O
p4 = np.array([23.672, 11.328, 30.466]) # CB
p5 = np.array([22.881, 10.326, 29.620]) # CG
p6 = np.array([23.691,  9.935, 28.389]) # CD1
p7 = np.array([22.557,  9.096, 30.459]) # CD2

dihed=dihedral(np.array([p0,p1,p2,p3]))
print(dihed)
