import numpy as np
import fnmatch as fnm
import warnings
from colorama import Fore, Back, Style, init
import os
import re
from ase.io import read
from pymatgen.core.composition import Composition
from common_functions import get_Ael
from energy_calculations import create_uHF
import itertools as it
from utility_functions import ins_ch_number, replace_1_string
from sklearn import linear_model
init(autoreset=True)

__version__ = "1.0"

#ToDo: Implement Decoraters----------##

"""Contains classes and  functions to calculate reaction enthalpy(formation energy) of :
1.MXene
2.Competing reactions
3. orthomxene (MX2ene)

ToDo: as general as possible"""


class Reader:

	"""Object consisting of methods to read MXenes, orthomxene, competing phases as required """

	def __init__(self, maxphase, mxene, criteria="notstrict", verbose=True):

		self.maxphase = maxphase
		self.mxene = mxene	
		self.verbose = verbose
		self.criteria = criteria  # criterion with max stress 0.01 or 0.001 eV/ A^3

	def readMXenes(self, filename="rel-ions-1.txt", filename2="rel-ions.txt", folder="."):

		self.mxeneatomslist = []
		self.mxenepath = []

		for roots, dirs, files in os.walk(folder):


			if self.criteria == "strict":

				cond = fnm.fnmatch(roots, "*/%s/0,001relaxation" % self.mxene) or fnm.fnmatch(roots, "*%sF2-0,001relaxation/*" % self.mxene)
	
			else:

				cond = fnm.fnmatch(roots, "*/%sF2*" % self.mxene) or fnm.fnmatch(roots, "*%s"% self.mxene)

			if cond:
				if filename in files or filename2 in files:
	
					try:
						at = read("%s/%s" %(roots, filename), format="gpaw-out")
	
					except OSError:
	
						warnings.warn("OS Error faced\n Reading rel-ions.txt from folder %s" %roots)
						at = read("%s/%s" %(roots, filename2), format="gpaw-out")
	
					self.mxeneatomslist.append(at)
					self.mxenepath.append(roots)

		self.mxenepath = [i.replace("./", "") for i in self.mxenepath]

	# sorting
		self.mxenepath, self.mxeneatomslist = Reader._sort_multiple(self.mxenepath, self.mxeneatomslist)
	
		if self.verbose:

			print("MXene atoms list"+ Fore.RED + "%s" % self.mxeneatomslist)
			print("path" + Fore.CYAN + "%s" % self.mxenepath)

	def readcarbides(self):

		self.carbidesdict = {}
		mainpath = "/nfshome/malik/Documents/Thesis/"

		if self.criteria == "strict":
 
			warnings.warn("At the moment, only strict criteria is available for TiC")
			self.carbidesdict["TiC"] = read(mainpath + "TiC/0,001relaxation/500ev/6/static.txt")

		elif self.criteria == "notstrict":

			self.carbidesdict["Cr3C2"] = read(mainpath + "Cr3C2/relaxation/rel-all-1.txt")
			self.carbidesdict["TiC"] = read(mainpath + "TiC/500ev/6/static.txt")
		
		carbides_energies = self._get_energies_fu(verbose=self.verbose)
		self.en_carbides = carbides_energies
	
	def readorthomxene(self, filename="rel-ions-1.txt", filename2="rel-ions.txt", folder="."):

		"""Same as MX2ene"""

		atomslist = []
		path = []

		for roots, dirs, files in os.walk(folder): 
		
			# place conditions to skip directories here

			if fnm.fnmatch(roots, "*relaxation*"):

				continue

			if filename in files:

				try:
					at = read("%s/%s" %(roots, filename), format="gpaw-out")

				except OSError:

					warnings.warn("OS Error faced\n Reading rel-ions.txt from folder %s" %roots)
					at = read("%s/%s" % (roots, filename2))

				atomslist.append(at)
				path.append(roots)


		path = [i.replace("./", "") for i in path]

		self.orthomxeneatoms = atomslist
		self.orthomxene_paths = path

		if self.verbose:

			print("OrthoMXenes atoms list are" + Style.BRIGHT +Fore.YELLOW + "%s" % self.orthomxeneatoms)
			print("path" +Fore.CYAN+ "%s" %self.orthomxene_paths)

	def readMAXphase(self, retarn=False):

		if self.criteria == "notstrict":

			maxatom = read("/nfshome/malik/Documents/Thesis/Bulk-%s/relaxation/rel-all-1.txt" % self.maxphase)

		elif self.criteria == "strict":

			maxatom = read("/nfshome/malik/Documents/Thesis/Bulk-%s/0,001eVrelaxation/500ev/static.txt" % self.maxphase)
		else:

			warnings.warn("criteria %s is invalid, recheck again, may lead to inconsistencies" % self.criteria)

		if retarn:

			return maxatom

		self.maxatom = maxatom 
		self._get_MAXenergy() #	getting the energies also

	def _get_MAXenergy(self):
		# ToDo: change the Ael
		from orthomxene_functions import Comp
		maxcomp = Comp(self.maxphase)
		self.eMAXphase = {"%s" % self.maxphase : maxcomp.get_energy_fu(self.maxatom)}
		eldict = maxcomp.get_el_amt_dict()
		self.n_MAX = eldict["C"] # No. of C = 'n' # Will fail if "C" is not present
		self.Ael = get_Ael(maxcomp, verbose=self.verbose)

	def readallmaxphases(self, phases=["Ti2AlC", "Cr2AlC", "Cr2GaC", "Ti3AlC2"]):

		self.allmaxatoms = []	

		for i in phases:

			maxatom = self.readMAXphase(retarn=True)

			self.allmaxatoms.append(maxatom)

	def read_energies_mxene(self):

		energylist = []

		for i, at in enumerate(self.mxeneatomslist):

			comp = Composition(self.mxeneformulas_fold[i])
			num_perfu = comp.num_atoms # No. of atoms pe F.U
			en = at.get_potential_energy()/len(at) # eV/atom
			en = en * num_perfu # eV/F.U
			energylist.append(en)

			if self.verbose:
	
				print("No. of atoms per FU %s" % num_perfu)
				print("energy eV/cell %s" % at.get_potential_energy())
				print("energy eV/F.U %s" % en)

		self.mxeneenergies = energylist

	def _get_energies_fu(self, verbose=True):

		#ToDo: Need to be tested thoruglhy before using genrally

		atomsdict = self.carbidesdict

		encarbidesdict = {}	
	
		for carbide,at in atomsdict.items():
			
			comp = Composition(at.get_chemical_formula("metal", empirical=True))
			compkey = Composition(carbide)
	
			check_atomsperunit = compkey.num_atoms
			noatomsperunit = comp.num_atoms 
			assert noatomsperunit == check_atomsperunit 
	
			if verbose:
	
				print("No. of atoms in %s" % carbide + " per F.U :" + Fore.RED+ " %s" % noatomsperunit)
	
			en = at.get_potential_energy()/len(at) # ev/atom
			enperfu = en*noatomsperunit	# ev/fu
			
			encarbidesdict[carbide] = enperfu	
			
		return encarbidesdict 	
	
	def _sort_multiple(main, b):

		from common_functions import sort_multiple
		return sort_multiple(main, b)

	def get_mxenes_maxformulas(self):
		
		""" Get mxene formulas from the path i.e.(folder name)"""
		#ToDo : generalize this
		# what if folder name is incorrect

		mxeneformulas = []
		maxformulas = []

		for p in self.mxenepath:

			formula = re.findall("%sF?\d?" % self.mxene, p)
			mxeneformulas.append(formula)

			maxf = re.findall("%s" % self.maxphase, p)
			maxformulas.append(maxf)	

		self.mxeneformulas_fold = [i for j in mxeneformulas for i in j]
		self.mxeneF_fold = [re.findall("F\d?", i) for i in self.mxeneformulas_fold]

		self.maxformulas_fold = [i for j in maxformulas for i in j]
	
		if self.verbose:

			print("mxene formulas: %s" % self.mxeneformulas_fold)
			print("F %s" % self.mxeneF_fold)
			print()
			print("max formulas: %s" % self.maxformulas_fold)

	def orthomxene_formulas(self):

		maxformulas = []

		for p in self.orthomxene_paths:

			match = re.findall("%s" % self.maxphase, p)

			if not match:

				warnings.warn("MAX phase is not found in folder name of %s\nRecheck the folder name/input data" %p )

			maxformulas.append(match)

		self.ortho_maxformulas = [i for j in maxformulas for i in j]
	
	def create_legend(self):

		legend = []

		for p in self.mxenepath:

			match = re.findall("site-?[a-z]", p)

			if not match:

				match  = re.findall("%s" %self.mxene, p)

			legend.append(match)
	
		self.legend = [i for j in legend for i in j] # reduce one dimension-- nested list from re.findall

	def create_ortholegend(self):

		ortholegend = []

		for p in self.orthomxene_paths:

			match = re.findall("slab-\d", p)

			if not match:

				match = re.findall("relaxation", p)

			ortholegend.append(match)

		self.ortholegend = [i for j in ortholegend for i in j]

#ToDo: print information function

class Refenergies_handler:

	def __init__(self):

		self.ecomp = {}
		self.eEl = {}

	def own_dataenergies(self, criteria="notstrict" ):

		from data_functions import DFT_energies_001, DFT_energies
		if criteria is "notstrict":

			ecomp, eEl = DFT_energies()

		elif criteria is "strict":

			ecomp, eEl = DFT_energies_001()
		else:

			raise Exception("Unhandled criteria entered")

		# updating

		self.ecomp = {**self.ecomp, **ecomp} # this will overwrite the common key values in first
		self.eEl = {**self.eEl, **eEl}


	def etchantsenergies(self):

		from energy_calculations import get_eHF, get_eH2O, get_eOH
		self.ecomp["HF"] = get_eHF(self.eEl)
		self.ecomp["H2O"] = get_eH2O(self.eEl)
		self.ecomp["OH"] = get_eOH(self.eEl)

class Calculator:
	
	"""Class to get reaction energies
	Don't confuse it with ase calculator"""


	def __init__(self, Refenergies_handler, Reader):
		
		self.refenergies = Refenergies_handler
		self.reader = Reader # assuming that relevant attributes are already present

	def formationenergymxenes(self):

		deltaGlist = []
		deltauHFlist = []
		Limits = []
		Ael = None
		
		maxcomp = Composition(self.reader.maxphase)
		maxeldict = maxcomp.get_el_amt_dict()
		self.n_MAX = maxeldict["C"] # No of C = 'n' # Will fail if "C" is not present
		self.Ael = get_Ael(maxcomp, verbose=self.reader.verbose)
		warnings.warn("Assuming that MAX phase is same for all read MXenes")

		if self.reader.verbose:

			print("A element is : " + Fore.YELLOW + "%s" % self.Ael)
	
		e_H = self.refenergies.eEl["H"]
		eBulkAF3 = self.refenergies.ecomp["%sF3" %self.Ael]

		for i, en in enumerate(self.reader.mxeneenergies):
			
			e_Max = self.reader.eMAXphase[self.reader.maxformulas_fold[i]]

			mxenecomp = Composition(self.reader.mxeneformulas_fold[i])
			eldict = mxenecomp.get_el_amt_dict()

			try:
				noF = eldict["F"]

			except KeyError:
				
				noF = 0

			if self.reader.verbose:

				print("F " +Style.BRIGHT+ Fore.MAGENTA+  "%s" % noF)
				print("MXene " + Style.BRIGHT + Fore.YELLOW + "%s" % self.reader.mxeneformulas_fold[i] )

			uHF, limits = self._uHF_related(e_MXene=en, e_Max=e_Max, noF=noF)

			deltauHFlist.append(uHF-limits[0]) # delta uHF --need for plotting
			limits = [i-limits[0] for i in limits] # deltalimits --needed for plotting
			
			if noF == 0:

				DelG = Formulas_handler.mxene_related.get_deltaG(uHF, e_Mxene=en, e_Max=e_Max, Bulk_AF3=eBulkAF3, e_H=e_H)

			elif noF == 2:

				DelG = Formulas_handler.mxene_related.get_deltaG_Fsaturated(uHF, e_MxeneF2=en, e_Max=e_Max, Bulk_AF3=eBulkAF3, e_H=e_H)

			else:
				raise Exception("Unhandled/unexpected F atoms in MXene\nmake sure the 'reading' part is correct")

			deltaGlist.append(DelG)
			Limits.append(limits)
	
			self.mxenedeltaG = deltaGlist # nested list
			self.mxenedeltauHF = deltauHFlist # nested list
			self.mxenelimits = Limits # nested list
			
	def _uHF_related(self, e_MXene, e_Max, noF):

		"""Also includes the uHF at zero deltaG"""

		uHF, ulimit, m_minlimit = create_uHF(self.refenergies.ecomp, self.refenergies.eEl, extra=0.5)
		zerouHF = Formulas_handler.mxene_related.get_zerouHF(e_Mxene=e_MXene, e_Max=e_Max, Bulk_AF3=self.refenergies.ecomp["%sF3" %self.Ael], e_H=self.refenergies.eEl["H"] , noF=noF)
		uHF = np.append(uHF, zerouHF)
		uHF = np.sort(uHF)

		return uHF, [ulimit, m_minlimit]


	def mxene_minimum_functionalized(self):

		selection = []
		finaluHF = []
		finallimit = []
		finalenergymxene = [] # legend missing

		uniquemaxphase_fold = np.unique(self.reader.maxformulas_fold) # names taken from folder/path
		uniquemappingmax = {i : (np.asarray(self.reader.maxformulas_fold) == i).nonzero()[0] for i in uniquemaxphase_fold}
		formulas = []
		for i, name in enumerate(uniquemaxphase_fold):

			# to find minium only in functionalized mxenes

			ind = [i for i,j in enumerate(self.reader.mxeneF_fold) if j]
				
			indices = uniquemappingmax[name] #  indices from MAX phase

			indices = list(set(indices).intersection(ind)) # overlap
			zeroindices = (np.isclose([self.mxenedeltaG[j] for j in indices], 0)).nonzero()
			minindices = np.argmin(np.asarray([self.mxenedeltauHF[j] for j in indices])[zeroindices], axis=0)
			fmainindices = indices[minindices]
	
			selection.append(self.mxenedeltaG[fmainindices])
			finaluHF.append(self.mxenedeltauHF[fmainindices])
			finallimit.append(self.mxenelimits[fmainindices])
			finalenergymxene.append(self.reader.mxeneenergies[fmainindices])
			formulas.append(self.reader.mxeneformulas_fold[fmainindices])

		# updating
		self.min_deltaG_mxene = selection
		self.min_uHF_mxene = finaluHF
		self.min_limits_mxene = finallimit
		self.min_finalen_mxene = finalenergymxene
		self.min_formulas = formulas

	def pristine_mxene(self):

		indices = [i for i,j in enumerate(self.reader.mxeneF_fold) if not j]
		
		self.formula_prist = [self.reader.mxeneformulas_fold[i] for i in indices]
		self.deltaG_prist = [self.mxenedeltaG[i] for i in indices]
		self.deltauHF_prist = [self.mxenedeltauHF[i] for i in indices]
		self.limits_prist = [self.mxenelimits[i] for i in indices]

	def printout_energies(self, detailed=0, show_all=True):

		print(Style.BRIGHT + Fore.RED + "----------Pristine--------------")
		print("MXene : %s" %self.formula_prist)
		
		print("Energies")
		print(*self.deltaG_prist, sep=",")

		if detailed == 1:
	
			print("deltauHF")
			print(*self.deltauHF_prist, sep=",")
			print(*self.limits_prist, sep=",")

		print(Style.BRIGHT + Fore.RED + "-------------Functionalized-(min)------------")
		print("MXene: %s" %self.min_formulas)
		print("Energies:")
		print(*self.min_deltaG_mxene, sep=",")

		if detailed == 1:
		
			print("deltauHF")
			print(*self.min_uHF_mxene, sep=",")
			print(*self.min_limits_mxene, sep=",")

		if show_all:

			print(Style.BRIGHT + Fore.RED + "----all------")
			print("MXene formulas: %s" % self.reader.mxeneformulas_fold)
			print("Energies:")
			print(*self.mxenedeltaG, sep=",")
			
			if detailed == 1:

				print("deltauHF")
				print(*self.mxenedeltauHF, sep=",")
				print(*self.mxenelimits, sep=",")


class side_reactions:
	
	def __init__(self, Reader, Refenergies_handler ):
		
		self.reader = Reader
		self.refenergies = Refenergies_handler
		
		# update refenergies with carbides
		self.reader.readcarbides()

		for c, en in self.reader.en_carbides.items():

			if c in self.refenergies.ecomp:

				os.warnings.warn("Key value is already present in refenergies. It will be overwritten")

			self.refenergies.ecomp[c] = en

	def formationenergy(self, ret=False):

		n = self.reader.n_MAX
		e_Max = self.reader.eMAXphase[self.reader.maxphase]
	
		side_formulas = Formulas_handler.side_reactions(n=n, ecomp=self.refenergies.ecomp, eEl=self.refenergies.eEl, maxphase=self.reader.maxphase, emaxphase=e_Max, Ael=self.reader.Ael)

		DeltauHF = []
		DeltaG = []
		legends = []
		Llimits = []
	# reaction type a
		for r in ["a", "b"]:

			uHF, limits = self._uHF_related(r, side_formulas)
			
			if r == "a":

				deltag, leg = side_formulas._deltaG_a(uHF)

			elif r == "b":

				deltag, leg = side_formulas._deltaG_b(uHF)

			legends.append(leg)
			DeltaG.append(deltag)
			DeltauHF.append(uHF - limits[0])
			Llimits.append([i-limits[0] for i in limits])
		
		if ret:
			return legends, DeltaG, DeltauHF, Llimits

		self.legends = legends
		self.DeltaG = DeltaG
		self.DeltauHF = DeltauHF
		self.Llimits = Llimits


	def _uHF_related(self, reaction_type, formulas):

		uHF, ulimit, m_minlimit = create_uHF(self.refenergies.ecomp, self.refenergies.eEl, extra=0.5)

		if reaction_type == "a":

			zerouHF = formulas._zerouHF_a(Ael=self.reader.Ael)

		elif reaction_type == "b":

			zerouHF = formulas._zerouHF_b(Ael=self.reader.Ael)

		if zerouHF:
	
			if zerouHF > np.amin(uHF) and zerouHF < np.amax(uHF):
	
				uHF = np.append(uHF, zerouHF)
				uHF = np.sort(uHF)

		return uHF, [ulimit, m_minlimit] 

	def printout(self, detailed=0):

		print(Style.BRIGHT + Fore.RED + "-----------Side_reactions----------")

		for i,r in enumerate(["a", "b"]):

			print("Reaction-%s" % r)
			print("%s" % self.legends[i])

			print("Energies")
			print(self.DeltaG[i], sep=",")

			if detailed == 1:

				print("deltauHF")
				print(self.DeltauHF[i], sep=",")
				print(self.lLimits[i], sep=",")


class orthomxene_calculator:

	def __init__(self, Reader, Refenergies_handler):

		# input attributes
		self.refenergies = Refenergies_handler
		self.reader = Reader
	
		# calculated attributes	
		self.uHF = None
		self.limits = None
		self.DeltaG = None
		self.compositions = None

	def wrapper_calc(self): 

		from orthomxene_functions import calculator_function
		
		uHF, limits = self._uHF_related()
		DeltaG, compositionlist, uHFlist = calculator_function(self.reader.orthomxeneatoms, self.reader.maxatom, self.refenergies.ecomp, self.refenergies.eEl, uHF)
	
		self.deltauHF = uHFlist # nested lists
		self.limits = [i-limits[0] for i in limits]
		self.DeltaG = DeltaG
		self.compositions = compositionlist

	def _uHF_related(self):

		uHF, ulimit, minlimit = create_uHF(self.refenergies.ecomp, self.refenergies.eEl)
		return uHF, [ulimit, minlimit]

	def minimum_deltaG(self):

		selection = []
		finaluHF = []
		finallimit = []
		compositions = []
		
		uniquemaxphase_fold = np.unique(self.reader.ortho_maxformulas)
		uniquemappingmax = {i : (np.asarray(self.reader.ortho_maxformulas) == i).nonzero()[0] for i in uniquemaxphase_fold}

		for i, name in enumerate(uniquemaxphase_fold):

			# to find minimum in orthomxenes

			indices = uniquemappingmax[name]

			zeroindices = []

			for j in indices:

				zeroindices.append(np.isclose(self.DeltaG[j], 0).nonzero()[0])

			zeroindices = [i for j in zeroindices for i in j]

			zerouHF = []

			for i,j in enumerate(zeroindices): # Need testing

				zerouHF.append(self.deltauHF[i][j])
			
			minindices = np.argmin(zerouHF, axis=0)
			fmainindices = indices[minindices]
		
			selection.append(self.DeltaG[fmainindices])
			finaluHF.append(self.deltauHF[fmainindices])
			compositions.append(self.compositions[fmainindices])

		# creating attributes
		self.min_deltaG = selection
		self.min_deltauHF = finaluHF
		self.min_compositions = compositions

	def printout(self, detailed=0):

		print(Style.BRIGHT + Fore.RED + "--------Orthomxene(Functionalized)-------")
		print("MXene: %s" % self.min_compositions)

		print("Energies:")
		print(*self.min_deltaG, sep=",")
		
		if detailed == 1:

			print("deltauHF")
			print(*self.min_deltauHF, sep=",")
			print(*self.min_limit, sep=",")

class Formulas_handler:

	def __init__(self):
	
		pass

	class mxene_related:

		def get_zerouHF( e_Mxene, e_Max, Bulk_AF3, e_H, noF):

#			""" returns the x(uHF)-point at which deltaG is zero"""
			uHF_zero = (e_Mxene + Bulk_AF3  +  (3+noF)*e_H - e_Max)/(3.0+ noF)

			return uHF_zero

		def get_deltaG(uHF, e_Mxene, e_Max, Bulk_AF3, e_H):

#			""" Calculates Delta G of a Mn+1AXn type MXene. The AF compound must be of AF3"""

			 deltaG = e_Mxene + Bulk_AF3  +  3*e_H - e_Max - 3*uHF
			  
			 return deltaG

		def get_deltaG_Fsaturated(uHF, e_MxeneF2, e_Max, Bulk_AF3, e_H):

#			""" Calculates DeltaG for a saturated MXeneMn+1AXnF2 AF3 compounds forms during the reaction"""
		
			deltaG  = e_MxeneF2 + Bulk_AF3 + 5*e_H - e_Max - 5*uHF
		
			return deltaG

	class side_reactions:

		def __init__(self, n, ecomp, eEl, maxphase, emaxphase, Ael):

			self.n = n
			self.ecomp = ecomp
			self.eEl = eEl
			self.emaxphase = emaxphase
			self.maxphase = maxphase
			self.Ael = Ael
		

		def _deltaG_a(self, uHF):
			
			if self.maxphase == "Ti2AlC" or self.maxphase == "Ti3AlC2":

				deltaG = self.n*self.ecomp["TiC"] + self.ecomp["TiF3"] + self.ecomp["AlF3"] + 6*self.eEl["H"] - self.emaxphase - 6*uHF

				legend = ["TiC and TiF$_3$"]

			elif self.maxphase == "Cr2AlC" or self.maxphase == "Cr2GaC":

				deltaG = self.ecomp["Cr3C2"] + self.ecomp["CrF3"] + 2*self.ecomp["%sF3" %self.Ael] + 9*self.eEl["H"] - 2*self.emaxphase - 9*uHF
				legend = ["Cr$_3$C$_2$ and CrF$_3$" ]	

			return	deltaG, legend

		def _deltaG_b(self, uHF):

			if self.maxphase == "Ti2AlC" or self.maxphase == "Ti3AlC2":

				DeltaG = self.n*self.ecomp["TiC"] + self.eEl["Ti"] + self.ecomp["AlF3"] + 3*self.eEl["H"] - self.emaxphase - 3*uHF
				legend = ["TiC and Ti"]


			elif self.maxphase == "Cr2AlC" or self.maxphase == "Cr2GaC":

			
				DeltaG = self.ecomp["Cr3C2"] + self.eEl["Cr"] + 2*self.ecomp["%sF3" % self.Ael] + 6*self.eEl["H"] - 2*self.emaxphase - 6*uHF
				legend = ["Cr$_3$C$_2$ and Cr"]

			return DeltaG, legend

			
		def _zerouHF_a(self, Ael):
			
			if self.maxphase == "Ti2AlC" or self.maxphase == "Ti3AlC2":

				zerouHF = (self.n*self.ecomp["TiC"] + self.ecomp["TiF3"] + self.ecomp["AlF3"] + 6*self.eEl["H"] - self.emaxphase)/6.0

			elif self.maxphase == "Cr2AlC" or self.maxphase == "Cr2GaC":
				
		                zerouHF = (self.ecomp["Cr3C2"] + self.ecomp["CrF3"] + 2*self.ecomp["%sF3" %Ael] + 9*self.eEl["H"] - 2*self.emaxphase)/9.0

			return zerouHF


		def _zerouHF_b(self, Ael):

			if self.maxphase == "Ti2AlC" or self.maxphase == "Ti3AlC2":

				zerouHF = (self.n*self.ecomp["TiC"] + self.eEl["Ti"] + self.ecomp["AlF3"] + 3*self.eEl["H"] - self.emaxphase)/3.0

			elif self.maxphase == "Cr2AlC" or self.maxphase == "Cr2GaC":
				
				zerouHF = (self.ecomp["Cr3C2"] + self.eEl["Cr"] + 2*self.ecomp["%sF3" %Ael] + 6*self.eEl["H"] - 2*self.emaxphase)/6.0

			return zerouHF


class Plotter:
	

	def __init__(self):

		pass
	
	def oldplotter(x, energylist, legend, *xlimits, **kwargs):

		from orthomxene_functions import global_plot 

		"""Uses the older plotter function"""

		global_plot(x, energylist, legend, *xlimits, **kwargs)

class Analyzer:

	def find_intersect(x1, y1, x2, y2, verbose=True):

		"""Find the intersection of two lines(linear) by linear regression"""
	
		line1 = linear_model.LinearRegression().fit(x1.reshape(-1, 1), y1.reshape(-1,1))
		line2 = linear_model.LinearRegression().fit(x2.reshape(-1, 1), y2.reshape(-1, 1))
	
		determinent = line1.coef_ - line2.coef_
	
		if determinent != 0: # Non-linear
	
			x_intersect = (line2.intercept_ - line1.intercept_)/ determinent
			y_intersect = line1.predict(x_intersect)
	
		else:
	
			warnings.warn("The intersection point is not found. The two lines are parallel")
	
		return x_intersect[0, 0], y_intersect[0, 0]

if __name__ == "__main__":

#first part---related to directories/reading/inital setups-----

	verbose = False
	maxphase = "Cr2GaC"
	mxene = "Cr2C"
	criteria = "strict"
	save = False

# Reading part-----

	inp_reader = Reader(maxphase=maxphase, mxene=mxene, criteria=criteria, verbose=verbose,)
	
	# mxenes
	inp_reader.readMXenes(folder="%s" % maxphase) # read atoms
	
	inp_reader.get_mxenes_maxformulas() # get formulas from folder/paths
	inp_reader.read_energies_mxene() # get energies

	# read maxphase
	inp_reader.readMAXphase(False)

	# ref energies
	refeng = Refenergies_handler()
	refeng.own_dataenergies(criteria)
	refeng.etchantsenergies()

	# calculator
	tocalc = Calculator(refeng, inp_reader)	
	tocalc.formationenergymxenes()
	
	tocalc.mxene_minimum_functionalized()
	tocalc.pristine_mxene()

	# print
	tocalc.printout_energies()
#	x = tocalc.mxenedeltauHF
#	y = tocalc.mxenedeltaG
	
	# side_reactions
	tocalc_side = side_reactions(inp_reader, refeng )
	tocalc_side.formationenergy() # will exist as attributes---better use dictionary --too many attributes
	tocalc_side.printout()

	# orthomxene
	ortho_calc = orthomxene_calculator(inp_reader, refeng)
	ortho_calc.reader.readorthomxene(folder="%s/orthomxene" %maxphase)
	ortho_calc.reader.orthomxene_formulas()
	ortho_calc.wrapper_calc()
	ortho_calc.minimum_deltaG()
	ortho_calc.printout()

	ortho_calc.reader.create_legend()

	# combine_data for plotting ToDo: A very hard way!!! Change it

	y = [*tocalc.min_deltaG_mxene, *tocalc_side.DeltaG, *ortho_calc.min_deltaG]
	x = [*tocalc.min_uHF_mxene, *tocalc_side.DeltauHF, *ortho_calc.min_deltauHF]
	limits = [tocalc.min_limits_mxene, tocalc_side.Llimits, ortho_calc.limits] #ToDo: Use * to avoid on loop
	legends = [tocalc.min_formulas, *tocalc_side.legends, ortho_calc.min_compositions]

	limts = [i for i in it.chain(limits)]
	
	llimits = []

	for i in limits:
	
		for j in i:

			try:
				for k in j:

					llimits.append(k)

			except TypeError:

				llimits.append(k)

	llimits = np.unique(llimits)
	print(llimits)

# intersection

	x_int, y_int = Analyzer.find_intersect(x[0], y[0], x[1], y[1], verbose=verbose)
	extraverticals = [x_int]

# Plotter

	print(legends)
	legends = [ins_ch_number(replace_1_string(i)) if "$" not in i else i  for j in legends for i in j ]
	Plotter.oldplotter(x, y, legends, llimits, **{"filename":"Ti2AlC", "setcolors":True, "figsize":[5.0, 4.0], "show":True, "fontsize":14.5, "labelsize":10.7, "legendfontsize":12.2}, linewidth=2.3, extra_verticals=extraverticals)

	if save:

		from utility_functions import save_data
		outfilename = "test"
		alldata = save_data(x, y)
		alldata.create_header(tocalc.reader.legend)
		alldata.savedata(filename="%s" % outfilename)

