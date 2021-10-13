#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
torsion_analyzer - Program to calculate torsion angle for PDB, netcdf trajectory and xtc trajectory
"""

import sys, signal
signal.signal(signal.SIGINT, signal.SIG_DFL)

import argparse
import os
import parmed
import MDAnalysis
import netCDF4
import numpy as np
import tqdm
import itertools
import csv

from mods.func_prompt_io import check_exist, check_overwrite



# =============== class =============== #
class Structure:
	def __init__(self, topology_file, list_mask_axis, list_mask_dihedral):
		self._obj_topology = None
		self._list_atom_idx_axis = []
		self._list_atom_idx_dihedral = []
		self._is_midpoint = True
		self._apply_mass = True

		self._obj_topology = parmed.load_file(topology_file)
		self._list_atom_idx_axis = [list(parmed.amber.AmberMask(self._obj_topology, mask).Selected()) for mask in list_mask_axis]
		self._list_atom_idx_dihedral = [list(parmed.amber.AmberMask(self._obj_topology, mask).Selected()) for mask in list_mask_dihedral]
		if len(self._list_atom_idx_axis[0]) != 1:
			self._is_midpoint = False


	@property
	def topology(self):
		return self._obj_topology


	def set_apply_mass(self, apply_mass):
		"""
		Method to apply_mass flag (use in calculation of mid-point)

		Args:
			apply_mass (bool): apply mass to atom coordinates in calculation of mid-point

		Returns:
			self
		"""
		self._apply_mass = apply_mass
		return self


	def get_dihedrals(self, coordinates):
		"""
		Method to get dihedral angles

		Args:
			coordinates (ndarray): coordinates for structure

		Returns:
			list: [dihedral_for_mask1(float), ...]
		"""
		coords_axis = None
		if self._is_midpoint == False:
			mid_points = []
			for list_atom_idx in self._list_atom_idx_axis:
				coords = coordinates[list_atom_idx][:]
				n = len(list_atom_idx)
				if self._apply_mass:
					mass = np.array([self._obj_topology.atoms[atom_idx].mass for atom_idx in list_atom_idx])
					coords = coords * np.array([mass]).T
					n = np.sum(mass)
				mid_points.append(np.mean(coords, axis=0) / n)
			coords_axis = np.array(mid_points)
		else:
			coords_axis = coordinates[[atom_idx for list_atom_idx in self._list_atom_idx_axis for atom_idx in list_atom_idx]][:]
		obj_axis = Axis(coords_axis)
		dihedrals = [obj_axis.get_dihedral(*coordinates[list_atom_idx][:]) for list_atom_idx in self._list_atom_idx_dihedral]
		return dihedrals


class Axis:
	def __init__(self, coordinates=None):
		self._vector = None
		self._coords = []

		if coordinates is not None:
			self.determine_axis_vector(coordinates)

	@property
	def vector(self):
		return self._vector

	@property
	def coords(self):
		return self._coords

	@property
	def coord1(self):
		return self._coords[0]

	@property
	def coord2(self):
		return self._coords[1]


	def determine_axis_vector(self, coordinates):
		"""
		Method to determine axis vector and coordinates

		Args:
			coordinates (ndarray): [[x1, y1, z1], ...]

		Returns:
			tuple: ([x1(np.float), y1(np.float), z1(np.float)], [x2(np.float), y2(np.float), z2(np.float)])
			ndarray: [x(np.float), y(np.float), z(np.float)]
		"""
		results = []	# [[sum_L, coord_pair, axis_vector], ...]
		for coords in itertools.combinations(coordinates, 2):
			# 座標ペアで回す

			# 仮の中心軸ベクトル
			u = np.array(coords[1]) - np.array(coords[0])
			if np.linalg.norm(u) == 0.0:
				continue

			sum_L = 0
			for coord in coordinates:
				# 各点で回す

				# 点までのベクトル
				v = np.array(coord) - np.array(coords[0])

				# 中心軸から点までの最短距離
				L = np.linalg.norm(np.cross(u, v) / np.linalg.norm(u))
				sum_L += L
			results.append([sum_L, coords, u])
		results = sorted(results, key=lambda x:x[0])[0]
		self._coords = results[1]
		self._vector = results[2]
		return self


	def get_nvector(self, coordinate):
		"""
		Method to get normal vector for given coordinate

		Args:
			coordinate (ndarray or list): [x(np.float), y(np.float), z(np.float)]

		Returns:
			ndarray: [x(np.float), y(np.float), z(np.float)]
		"""
		v = np.array(coordinate) - self._coords[0]
		w = np.cross(self._vector, v)
		w_norm = np.linalg.norm(w)
		return w / w_norm


	def get_dihedral(self, coordinate1, coordinate2):
		"""
		Method to get dihedral angle for given two coordinates

		Args:
			coordinate1 (ndarray or list): [x1(np.float), y1(np.float), z1(np.float)]
			coordinate2 (ndarray or list): [x2(np.float), y2(np.float), z2(np.float)]

		Returns:
			float: dihedral angle
		"""
		vec1 = np.array(coordinate1)
		vec2 = np.array(coordinate2)
		n_vec1 = self.get_nvector(vec1)
		n_vec2 = self.get_nvector(vec2)

		ip = np.inner(n_vec1, n_vec2)
		norm = np.linalg.norm(n_vec1) * np.linalg.norm(n_vec2)
		angle = np.rad2deg(np.arccos(np.clip(ip / norm, -1.0, 1.0)))
		return angle



# =============== main =============== #
if __name__ == '__main__':
	parser = argparse.ArgumentParser(description="torsion_analyzer.py - Program to calculate torsion angle for PDB, netcdf trajectory and xtc trajectory", formatter_class=argparse.RawTextHelpFormatter)
	parser.add_argument("-p", dest="TOPOLOGY_FILE", metavar="TOPOLOGY_FILE", help="topology file for trajectory files (.pdb for .pdb, .prmtop for .nc, and .gro for .xtc)")
	parser.add_argument("-x", dest="COORDINATE_FILE", metavar="COORDINATE_FILE", required=True, help="coordinate files (.pdb, .nc and .xtc)")
	parser.add_argument("-ma", dest="MASK_AXIS", metavar="AMBERMASK_AXIS", required=True, nargs="+", help="Ambermask for axis (single mask of one arg: mid point atom; twi masks of one arg: pair atoms for mid point)")
	parser.add_argument("-mb", dest="MASK_TORSION", metavar="AMBERMASK_TORSION", required=True, nargs="+", help="torsion angle atoms")
	parser.add_argument("-o", dest="OUTPUT_FILE", metavar="OUTPUT_FILE", required=True, help="output file (.csv or .txt)")
	parser.add_argument("-b", dest="BEGIN", metavar="START_TIME", type=int, help="First frame index to read from trajectory (start from 0)")
	parser.add_argument("-e", dest="END", metavar="END_TIME", type=int, help="Last frame index to read from trajectory (start from 0)")
	parser.add_argument("--offset", dest="OFFSET", metavar="OFFSET", type=int, default=1, help="Output interval (Default: 1)")
	parser.add_argument("-O", dest="FLAG_OVERWRITE", action="store_true", default=False, help="overwrite forcibly")
	args = parser.parse_args()

	check_exist(args.TOPOLOGY_FILE, 2)
	check_exist(args.COORDINATE_FILE, 2)

	obj_trajectory = None
	obj_structure = Structure(args.TOPOLOGY_FILE, args.MASK_AXIS, args.MASK_TORSION)
	list_dihedrals = []
	if os.path.splitext(args.TOPOLOGY_FILE)[1].lower() == ".prmtop":
		if os.path.splitext(args.COORDINATE_FILE)[1].lower() != ".nc":
			sys.stderr.write("ERROR: .prmtop needs .nc file.\n")
			sys.exit(1)

		obj_trajectory = netCDF4.Dataset(args.COORDINATE_FILE)
		frame_begin = 0
		if args.BEGIN is not None:
			frame_begin = args.BEGIN
		frame_end = len(obj_trajectory.dimensions["frame"])
		if args.END is not None:
			frame_end = args.END
		frame_skip = args.OFFSET
		if frame_skip == 0:
			frame_skip = 1
		for frame in tqdm.tqdm(range(frame_begin, frame_end, frame_skip), ascii=True):
			list_dihedrals.append([frame] + obj_structure.get_dihedrals(obj_trajectory.variables["coordinates"][frame].data))

	elif os.path.splitext(args.TOPOLOGY_FILE)[1].lower() == ".gro":
		if os.path.splitext(args.COORDINATE_FILE)[1].lower() != ".xtc":
			sys.stderr.write("ERROR: .gro needs .xtc file.\n")
			sys.exit(1)
		obj_trajectory = MDAnalysis.Universe(args.TOPOLOGY_FILE, args.COORDINATE_FILE)
		frame = -1
		for obj_frame in tqdm.tqdm(obj_trajectory.trajectory, ascii=True):
			frame += 1
			if args.BEGIN is not None and args.BEGIN > frame:
				continue

			if args.END is not NOne and args.END < frame:
				break

			if args.OFFSET is not None:
				n_offset += 1
				if n_offset != args.OFFSET:
					continue
				n_offset = 0

			list_dihedrals.append([frame] + obj_structure.get_dihedrals(obj_frame.positions))

	elif os.path.splitext(args.TOPOLOGY_FILE)[1].lower() == ".pdb":
		if os.path.splitext(args.COORDINATE_FILE)[1].lower() != ".pdb":
			sys.stderr.write("ERROR: .pdb needs .pdb file.\n")
			sys.exit(1)

		list_dihedrals.append([1] + obj_structure.get_dihedrals(obj_structure.topology.coordinates))

	else:
		sys.stderr.write("ERROR: undefined file type for topology file\n")
		sys.exit(1)


	if args.FLAG_OVERWRITE == False:
		check_overwrite(args.OUTPUT_FILE)

	with open(args.OUTPUT_FILE, "w") as obj_output:
		writer = None
		if os.path.splitext(args.OUTPUT_FILE)[1].lower() == ".csv":
			writer = csv.writer(obj_output)
		else:
			writer = csv.writer(obj_output, delimiter="\t")

		writer.writerow(["Frame"] + args.MASK_TORSION)
		writer.writerows(list_dihedrals)
