#!/home/rohit/anaconda3/bin/python

import tkinter as tk
import tkinter.filedialog as filedialog
import numpy as np
import urllib
import urllib.request
import matplotlib.pyplot as plt
import pdb

root = tk.Tk()

def fetch_pdb():
	pdb_id = str(e1.get()).split(".")
	pdb_url = "http://files.rcsb.org/download/" + pdb_id[0] + ".pdb"
	pdb_file_name=pdb_url.split("/")[-1]
	a, b = urllib.request.urlretrieve (pdb_url, pdb_file_name)
	try:
		with open(pdb_file_name,'r') as f:
			global data
			data = [i.strip().split() for i in f.readlines()]
			fetchsubmit = tk.Button(root, text="Submit", bg="green", command=Submit).grid(row=2, column=3, pady=2)
	except:
		print("No file exists")

def browsefunc():
	# pdb.set_trace()
	filename = filedialog.askopenfilename(title = "Select PDB file",filetypes = (("pdb files","*.pdb"),("all files","*.*")))
	try:
		with open(filename,'r') as f:
			global data
			data = [i.strip().split() for i in f.readlines()]
			browsesubmit = tk.Button(root, text="Submit", bg="green", command=Submit).grid(row=4, column=2, padx=10)
	except:
		print("No file exists")

def Submit():
	all_crd=[]
	des_crd=[]
	for line_items in data:
		if line_items[0] == "ATOM":
			all_crd.append([int(line_items[1]), line_items[2], line_items[3], int(line_items[5]), float(line_items[6]), float(line_items[7]), float(line_items[8])])
			#data1=np.array(crd, dtype=np.float)
			#x.append(float(line_items[6]))		
			#y.append(float(line_items[7]))		
			#z.append(float(line_items[8]))
		if line_items[0] == "TER":
			break

	for line_items in all_crd:
		if(line_items[1] == "N" or line_items[1] == "CA" or line_items[1] == "C"):
			des_crd.append(line_items)


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
	
	name_dict={}
	for i in des_crd:
		key=str(i[3])+i[1]
		name_dict[key]=i[4:7]
	data_di=[]
	for i in range(len(des_crd)-3):
		p0 = np.array(des_crd[i][4:7]) # N
		p1 = np.array(des_crd[i+1][4:7]) # CA
		p2 = np.array(des_crd[i+2][4:7]) # C
		p3 = np.array(des_crd[i+3][4:7]) # O
		# print(des_crd[i][0], des_crd[i+1][0], des_crd[i+2][0], des_crd[i+3][0])
		# pdb.set_trace()
		dihed=dihedral(np.array([p0,p1,p2,p3]))
		data_di.append(dihed*-1)
			
	print(data_di)
	psi, omega, phi=[],[],[]
	for i in range(0, len(data_di), 3):
		psi.append(data_di[i])
		omega.append(data_di[i+1])
		phi.append(data_di[i+2])
	#print(psi,'\n\n')
	#print(omega,'\n\n')
	#print(phi,'\n\n')
	
	fig = plt.figure()
	ax= fig.gca()
	#ax.set_axis_bgcolor("darkslategray")
	#ax.set_clip_on(False)
	ax.set_xticks(np.arange(-180.1,180.1,30))
	ax.set_yticks(np.arange(-180.1,180.1,30))
	plt.scatter(phi,psi,color='r')
	plt.xlim(-180,180)
	plt.ylim(-180,180)
	plt.grid()
	
	xline = np.array(list(range(-180,180)))
	yline = np.array([0]*360)
	plt.plot(xline,yline,color='b',linewidth=1.5)
	plt.plot(yline,xline,color='b',linewidth=1.5)
	plt.xlabel('Phi')
	plt.ylabel('Psi')
	plt.title('Ramachandran Plot')
	plt.show()
	
pathlabel = tk.Label(root, text="\n\nPlease submit PDB file to get Ramachandran Plot\n", font = "Times 20 italic").grid(row=1, columnspan=4)

tk.Label(root, text="PDB ID").grid(row=2)
e1 = tk.Entry(root)
e1.grid(row=2, column=1)

fetchsubmit = tk.Button(root, text="Fetch PDB file", command=fetch_pdb).grid(row=2, column=2, pady=2)
tk.Label(root, text="OR", font = "Times 20 italic").grid(row=3, columnspan=4, pady=2)

browsebutton = tk.Button(root, text="Browse", command=browsefunc).grid(row=4, column=1, pady=2)
button = tk.Button(root, text="QUIT", command=quit).grid(row=4, column=3, pady=2)

root.mainloop()

