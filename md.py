#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  boltzmann_factor.py
#  
#  Copyright 2019 Rafa <rafa@Frost>
#  
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#  
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#  
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.
#  
#  



# este programa gera a distribuição de normal de velocidade de um gas para uma unica componete

import numpy as np

import math
M  = 28.0
kb = 1.38E-23
T  = 273.15
NA = 6.022E23
m  = (M/1000)/ NA
dt = 1E-12
number_of_steps = 1000
r  = 2.0


#parameters
natoms =  5 
sigma  = math.sqrt((kb * T)/m)
mu     = 0 
box    = [(-50,50), (-50, 50), (-50, 50)]  


pos = np.random.rand(natoms, 3)

ux  = (np.random.normal(mu, sigma, natoms))*1E10
uy  = (np.random.normal(mu, sigma, natoms))*1E10
uz  = (np.random.normal(mu, sigma, natoms))*1E10

fx = np.zeros(natoms)
fy = np.zeros(natoms)
fz = np.zeros(natoms)



def export_cell (box):
	""" Function doc 
	ATOM  1      H   LIG  1        -50.00  -50.00  -50.000   1.00  300.00          H 0000  
	ATOM  2      H   LIG  1        -50.00   50.00  -50.000   1.00  300.00          H 0000  
	ATOM  3      H   LIG  1        -50.00  -50.00   50.000   1.00  300.00          H 0000  
	ATOM  4      H   LIG  1         50.00  -50.00  -50.000   1.00  300.00          H 0000  
	ATOM  5      H   LIG  1         50.00   50.00  -50.000   1.00  300.00          H 0000  
	ATOM  6      H   LIG  1         50.00  -50.00   50.000   1.00  300.00          H 0000  
	ATOM  7      H   LIG  1        -50.00   50.00   50.000   1.00  300.00          H 0000  
	ATOM  8      H   LIG  1         50.00   50.00   50.000   1.00  300.00          H 0000  
	p0_0_0 = "ATOM  1      H   LIG  1         %s  %s  %s   1.00  300.00          H 0000\n" %(str(box[0][0]), str(box[1][0]), str(box[2][0]))
	p0_1_0 = "ATOM  2      H   LIG  1         %s  %s  %s   1.00  300.00          H 0000\n" %(str(box[0][0]), str(box[1][1]), str(box[2][0]))
	p0_0_1 = "ATOM  3      H   LIG  1         %s  %s  %s   1.00  300.00          H 0000\n" %(str(box[0][0]), str(box[1][0]), str(box[2][1]))
	p1_0_1 = "ATOM  4      H   LIG  1         %s  %s  %s   1.00  300.00          H 0000\n" %(str(box[0][1]), str(box[1][0]), str(box[2][0]))
	p1_1_0 = "ATOM  5      H   LIG  1         %s  %s  %s   1.00  300.00          H 0000\n" %(str(box[0][1]), str(box[1][1]), str(box[2][0]))
	p1_0_1 = "ATOM  6      H   LIG  1         %s  %s  %s   1.00  300.00          H 0000\n" %(str(box[0][1]), str(box[1][0]), str(box[2][1]))
	p0_1_1 = "ATOM  7      H   LIG  1         %s  %s  %s   1.00  300.00          H 0000\n" %(str(box[0][0]), str(box[1][1]), str(box[2][1]))
	p1_1_1 = "ATOM  8      H   LIG  1         %s  %s  %s   1.00  300.00          H 0000\n" %(str(box[0][1]), str(box[1][1]), str(box[2][1]))
	"""

	cellout = open('cell_xyz.pdb', 'w')
	string  = "ATOM  1      H   LIG  1         %s  %s  %s   1.00  300.00          H 0000\n" %(str(box[0][0]), str(box[1][0]), str(box[2][0]))
	string += "ATOM  2      H   LIG  1         %s  %s  %s   1.00  300.00          H 0000\n" %(str(box[0][0]), str(box[1][1]), str(box[2][0]))
	string += "ATOM  3      H   LIG  1         %s  %s  %s   1.00  300.00          H 0000\n" %(str(box[0][0]), str(box[1][0]), str(box[2][1]))
	string += "ATOM  4      H   LIG  1         %s  %s  %s   1.00  300.00          H 0000\n" %(str(box[0][1]), str(box[1][0]), str(box[2][0]))
	string += "ATOM  5      H   LIG  1         %s  %s  %s   1.00  300.00          H 0000\n" %(str(box[0][1]), str(box[1][1]), str(box[2][0]))
	string += "ATOM  6      H   LIG  1         %s  %s  %s   1.00  300.00          H 0000\n" %(str(box[0][1]), str(box[1][0]), str(box[2][1]))
	string += "ATOM  7      H   LIG  1         %s  %s  %s   1.00  300.00          H 0000\n" %(str(box[0][0]), str(box[1][1]), str(box[2][1]))
	string += "ATOM  8      H   LIG  1         %s  %s  %s   1.00  300.00          H 0000\n" %(str(box[0][1]), str(box[1][1]), str(box[2][1]))
	string += "CONECT    1    2    3    4 \n"
	string += "CONECT    2    7    5      \n"
	string += "CONECT    3    7    6      \n"
	string += "CONECT    4    5    6      \n" 
	string += "CONECT    5    8           \n"
	string += "CONECT    6    8           \n"
	string += "CONECT    7    8           \n"
	cellout.write(string)







def write_xyz_file (pos):
	""" Function doc """
	
	trajout = open('traj.xyz', 'a')
	lines   = [] 
	lines.append("" + str(natoms)+"\n")
	lines.append("\n")

	for i in range(natoms):
		lines.append("Ar  "+str(pos[i][0])+ "   "+ str(pos[i][1])+"   "+str(pos[i][2])+"\n")
	
	trajout.writelines(lines)
	
	

def computeForce (pos, mass):
	""" Function doc """
	pass


def integrate (pos, ux, uy, uz, fx, fy, fz, mass, dt):
	""" Function doc """
	#usaremos a formula de euler - tem coisa melhor
	for i in range(natoms):
		print pos[i][0],  pos[i][1], pos[i][2]
		pos[i][0]  += ux[i] * dt
		ux[i]      += fx[i] * dt / m

		pos[i][1]  += uy[i] * dt
		ux[i]      += fy[i] * dt / m
	
		pos[i][2]  += uz[i] * dt
		ux[i]      += fz[i] * dt / m

def sphereHitCheck (pos, ux, uy, uz, fx, fy, fz):
	""" Function doc """
	for i in range(natoms):
		for j in range (i+1, natoms):
			#check X
			
			
			d_x = pos[i][0] - pos[j][0]
			d_y = pos[i][1] - pos[j][1]
			d_z = pos[i][2] - pos[j][2]

			dist = (d_x**2 + d_y**2 + d_z**2) 
			
			
			

			#if diff_x <= r:
			#	ux[i] = -ux[i]
			#	ux[j] = -ux[j]
			#
			#if diff_y <= r:
			#	uy[i] = -uy[i]
			#	uy[j] = -uy[j]
			#
			#if diff_z <= r:
			#	uz[i] = -uz[i]
			#	uz[j] = -uz[j]
			
			#print pos[i][0],  pos[i][1], pos[i][2]
			#pos[i][0]  += ux[i] * dt
			#ux[i]      += fx[i] * dt / m
			#
			#pos[i][1]  += uy[i] * dt
			#ux[i]      += fy[i] * dt / m
			#
			#pos[i][2]  += uz[i] * dt
			#ux[i]      += fz[i] * dt / m	


def wallHitCheck (pos, ux, uy, uz, fx, fy, fz):
	""" Function doc """
	
	for i in range(natoms):
		#print pos[i][0],  pos[i][1], pos[i][2]
		
		pos[i][0]  += ux[i] * dt
		#'''
		#walls in X
		if pos[i][0] >= box[0][1]:
			diff      =  pos[i][0] - box[0][1]
			pos[i][0] = pos[i][0]-diff
			ux[i]     =  -ux[i]
		else:
			pass
			
		if pos[i][0] <= box[0][0]:
			diff      =  pos[i][0] - box[0][0]
			pos[i][0] = pos[i][0]-diff
			ux[i]     =  -ux[i]
		else:
			pass
		#'''



		#'''
		#walls in Y
		if pos[i][1] >= box[1][1]:
			diff      =  pos[i][1] - box[1][1]
			pos[i][1] = pos[i][1]-diff
			uy[i]     =  -uy[i]
		else:
			pass
			
		if pos[i][1] <= box[1][0]:
			diff      =  pos[i][1] - box[1][0]
			pos[i][1] = pos[i][1]-diff
			uy[i]     =  -uy[i]
		else:
			pass

		#'''
		
		
		
		#'''
		#walls in Z
		if pos[i][2] >= box[2][1]:
			diff      =  pos[i][2] - box[2][1]
			pos[i][2] = pos[i][2]-diff
			uz[i]     =  -uz[i]
		else:
			pass
			
		if pos[i][2] <= box[2][0]:
			diff      =  pos[i][2] - box[2][0]
			pos[i][2] = pos[i][2]-diff
			uz[i]     =  -uz[i]
		else:
			pass
		#'''
	pass
	


def run (number_of_steps):
	""" Function doc """
	export_cell (box)
	trajout = open('traj.xyz', 'w')
	
	for i in range (number_of_steps):
		
	
		forces =  computeForce(pos, m)
		
		integrate (pos, ux, uy, uz, fx, fy, fz, m, dt)

		sphereHitCheck(pos, ux, uy, uz, fx, fy, fz)

		wallHitCheck(pos, ux, uy, uz, fx, fy, fz)
		
		write_xyz_file(pos)

run(number_of_steps)












	
