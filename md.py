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
number_of_steps = 100

#parameter
natoms =  100


sigma = math.sqrt((kb * T)/m)
mu    = 0 



pos = np.random.rand(natoms, 3)

ux  = (np.random.normal(mu, sigma, natoms))*1E10
uy  = (np.random.normal(mu, sigma, natoms))*1E10
uz  = (np.random.normal(mu, sigma, natoms))*1E10

fx = np.zeros(natoms)
fy = np.zeros(natoms)
fz = np.zeros(natoms)

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


def wallHitCheck (pos, ux, uy, uz, fx, fy, fz):
	""" Function doc """
	
	pass
	


def run (number_of_steps):
	""" Function doc """

	for i in range (number_of_steps):
		
	
		forces =  computeForce(pos, m)
		
		integrate (pos, ux, uy, uz, fx, fy, fz, m, dt)

		wallHitCheck(pos, ux, uy, uz, fx, fy, fz)
		
		write_xyz_file(pos)

run(number_of_steps)












	
