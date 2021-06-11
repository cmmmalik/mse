
class Reader:

	"""Object consisting of methods to read MXenes, orthomxene, competing phases as required """

	def __init__(self, maxphase, mxene, criteria="notstrict", verbose=True):

		self.maxphase = {"name": maxphase} 
		self.mxene = mxene	
		self.verbose = verbose
		self.criteria = criteria  # criterion with max stress 0.01 or 0.001 eV/ A^3

	def readMXenes(self, filename="rel-ions-1.txt", filename2="rel-ions.txt", folder="."):

		mxeneatomslist = []
		mxenepath = []

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
	
					mxeneatomslist.append(at)
					mxenepath.append(roots)

		mxenepath = [i.replace("./", "") for i in self.mxenepath]

	# sorting

		mxenepath, mxeneatomslist = Reader._sort_multiple(mxenepath, mxeneatomslist)
	
		if self.verbose:

			print("MXene atoms list"+ Fore.RED + "%s" % mxeneatomslist)
			print("path" + Fore.CYAN + "%s" % mxenepath)

		return mxenepath, mxeneatomslist

	def readcarbides(self, mainpath="/nfshome/malik/Documents/Thesis/"):

		carbidesdict = {}

		if self.criteria == "strict":

			warnings.warn("At the moment, only strict criteria is available for TiC")
			carbidesdict["TiC"] = read(mainpath + "TiC/0,001relaxation/500ev/6/static.txt")

		elif self.criteria == "notstrict":

			carbidesdict["Cr3C2"] = read(mainpath + "Cr3C2/relaxation/rel-all-1.txt")
			carbidesdict["TiC"] = read(mainpath + "TiC/500ev/6/static.txt")
		
		carbides_energies = self._get_energies_fu(verbose=self.verbose)

		return carbides_energies, carbidesdict
	
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

		if self.verbose:

			print("OrthoMXenes atoms list are" + Style.BRIGHT +Fore.YELLOW + "%s" % atomslist)
			print("path" +Fore.CYAN+ "%s" % path)

		return atomslist, path

	def readMAXphase(self):

		if self.criteria == "notstrict":

			maxatom = read("/nfshome/malik/Documents/Thesis/Bulk-%s/relaxation/rel-all-1.txt" % self.maxphase)

		elif self.criteria == "strict":

			maxatom = read("/nfshome/malik/Documents/Thesis/Bulk-%s/0,001eVrelaxation/500ev/static.txt" % self.maxphase)
		else:

			warnings.warn("criteria %s is invalid, recheck again, may lead to inconsistencies" % self.criteria)

	
		_get_MAXenergy() #	getting the energies also

	def _get_MAXenergy(self):

		# ToDo: change the Ael
		from orthomxene_functions import Comp
		
		maxcomp = Comp(self.maxphase)
		self.eMAXphase = {"%s" % self.maxphase : maxcomp.get_energy_fu(self.maxatom)}
		eldict = maxcomp.get_el_amt_dict()
		self.n_MAX = eldict["C"] # No. of C = 'n' # Will fail if "C" is not present
		self.Ael = get_Ael(maxcomp, verbose=self.verbose)


	def readallmaxphases(self, phases=["Ti2AlC", "Cr2AlC", "Cr2GaC", "Ti3AlC2"]):

		allmaxatoms = []	

		for i in phases:

			maxatom = self.readMAXphase(retarn=True)

			allmaxatoms.append(maxatom)

		return allmaxatoms

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

		return energylist

	def _get_energies_fu(self, atomsdict, verbose=True):

		#ToDo: Need to be tested thoruglhy before using genrally

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

