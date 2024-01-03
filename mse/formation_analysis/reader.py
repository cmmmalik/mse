from ase import Atoms
from ase.io import read as aseread
from functools import partial 
import re

from preprocessing.scan import Searcher as _searcher , Read
from common_functions import sort_multiple
from preprocessing.atoms import  Aatoms
import os


def _extendatoms(ase_atoms):

	def decorator(self):
		
		atoms = ase_atoms(self)
		paths = self._paths

		for i,p in zip(atoms, paths):
			
			yield Aatoms(i, p)

	return decorator


def read(p):

	atoms = aseread(p)
	return Aatoms(atoms, p)


class Structures:

	"""Front-end class for structures of the same type"""

	def __init__(self, name, paths=None, atoms=None, iatoms=None):
		
		assert paths or atoms or iatoms
		self.name = name

		if paths:
			self.paths = paths

		if atoms:
			self.atoms = atoms

	def __repr__(self):

		display = "name={}".format(self.name)

		if hasattr(self, "_atoms"):

			display = display + ",{}".format(self.atoms)

		if hasattr(self, "_paths"):

			display = display + ",{}".format(self.paths)
		
		return "("+display+")"
	
	@_extendatoms
	def __getatoms__(self):
		""" Get atoms generator from paths """
		return Read.readatoms(self._paths)

	@property
	def atoms(self):

		return self._atoms

	@atoms.setter
	def atoms(self, value):

		if isinstance(value, list):
			self._atoms = value
		else:
			raise TypeError("{} does not have a {} type".format(value, list))

	@staticmethod
	def getpaths(filename, folder, *subfolder):

		scanner = _searcher(filename, folder)

		if subfolder:
			paths = [ k for j in subfolder for k in scanner.get_filespaths(tomatch=["%s" %j])]

	
		else:
			paths = scanner.get_filespaths()

		return paths

	@classmethod
	def getpaths_vfiles(cls, filenames, folder, *subfolder):

		paths = []	
		for f in filenames:
			paths += cls.getpaths(f, folder, *subfolder)
		
		return paths
	
	@staticmethod
	def select_paths_preference(filenames, paths):

		selpaths = []
		for p in paths:

			for f in filenames:

				if os.path.join(os.path.dirname(p), f) in selpaths:
					break	

				if f == os.path.basename(p):
					selpaths.append(p)
					break
						 
		return selpaths

	@classmethod
	def getpathsatoms(cls, filename, folder, *subfolder):  #change it to __call__ method
		"""Left for backward compatibility,
		returns dictionary 'paths' and 'atoms' """
	
		paths = cls.getpaths(filename, folder, *subfolder)
		
#		firstpa = scanner.get_filespaths(tomatch=["%s" %subfolder[0]])
#		secondpa = scanner.get_filespaths(tomatch=["%s" %subfolder[1]])
		
		# get atoms
		iatoms = cls._getatoms(paths)
	
		return  {"paths": paths, "atoms": iatoms}

	@classmethod
	def getstructures(cls, name, filename, folder, *subfolder):
		"""returns structures object"""
		if isinstance(filename, list):
			paths = cls.getpaths_vfiles(filename, folder, *subfolder) # will contain duplicates differ only by filename at the end
			paths = cls.select_paths_preference(filename, paths)
				
		else:

			paths  = cls.getpaths(filename, folder, *subfolder)
		return  Structures(name, paths, atoms=None)


	def get_formulas_path(self, tomatch):

		"""obtain formulas from the folder names, by matching the 'tomatch'"""
		return _searcher.searchstrings(self.paths, tomatch)
	
	def iget_energy_fu(self):

		for at in self._atoms:

			yield at.get_energy_per_fu()

	def get_energy_fu(self):

		return list(self.iget_energy_fu())

	@property
	def name(self):

		return self._name

	@name.setter
	def name(self, name):

		if isinstance(name, str):		
			self._name = name

		else:
			raise TypeError("'{}' type is not {}".format(name, str))

	@property
	def paths(self):

		return self._paths

	@paths.setter
	def paths(self, p):

		self._paths = p
		self._atoms = list(self.readatoms())

	def readatoms(self):

		return self.__getatoms__()


