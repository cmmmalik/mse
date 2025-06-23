from collections import UserDict, abc, OrderedDict
from pymatgen.core.composition import Composition
# from preprocessing.util import remove_1, float_int
from pymatgen.analysis.reaction_calculator import Reaction
import re
import numpy as np
from ase.formula import Formula
from chempy import balance_stoichiometry
from colorama import Fore, Back, init, Style

# from energy_calculations import create_uHF, create_uHCl
# from data_functions import DFT_energies_001, DFT_energies
# import sympy

#from warnings import warn
import warnings
#from warnings import RuntimeWarning
import functools as func

init(autoreset=True)

print = func.partial(print, flush=True) # enabling flush on all the print statements

class Formation_enthalpy:

	def __init__(self, max_at, Refenergies_dict, select_ecomp=["H2O", "HF", "OH"], select_eEl=["H", "O", "F"], default_reactants=[]):

		try:
			self.ecomp = Refenergies_dict.get_reduced_dict("ecomp", *select_ecomp)
			self.eEl = Refenergies_dict.get_reduced_dict("eEl", *select_eEl)

		except AttributeError:

			warnings.warn("Attributes of {0} could not be passed to {1}".format(Refenergies_dict, Formation_enthalpy))

		self.max_at = max_at
		self.reactants = (max_at.formula_composition.reduced_formula, *default_reactants)


	@property
	def max_at(self):

		return self._max_at 

	@max_at.setter
	def max_at(self, at):

		self._max_at = at

	@property
	def reactants(self):
	
		return self._reactants

	@reactants.setter
	def reactants(self, value):

		self._reactants = value

	def __getattr__(self, attr):


		try:
			try:
				getattr(self, "max_at")

			except AttributeError:
				pass
			else:

				return	getattr(self, "max_at").attr

		except AttributeError:
			pass

		try:

			return getattr(self, "ecomp")[attr]

		except KeyError:
			pass

		try:

			return getattr(self, "eEl")[attr]
		except KeyError:

			raise AttributeError("'{}' object does not have the attribute '{}'".format(self.__class__.__name__, attr))

	def __repr__(self):

		toshow = [self.max_at.__repr__()]

		for i in ["ecomp", "eEl"]:

			d = getattr(self, i).__repr__()
			toshow.append(d)

		return "{0}({1})".format(self.__class__.__name__,"\n".join(toshow))

	def mxene_reaction_stable(self, mxen_at, el_adsorb="F", verbose=True):
		"""More stable method than the new one!"""

		from energy_calculations import create_uHF, create_uHCl

		sys_comps = self.get_selection_compounds()


		AelF = sys_comps["A%s" %el_adsorb]
		f_af = re.search("(?<=%s)[0-9]*.?[0-9]*" %el_adsorb, AelF).group()
		print(f_af)
		
		if not f_af or f_af != "3":

			raise ValueError("reactions involvin {} can't be handled at the moment; only AF3 types are possible".format(AelF))	
			
		e_mxene = mxen_at.get_energy_fu()
		e_Max = self.max_at.get_energy_fu()
		e_H = self.H
		e_AF3 = getattr(self, AelF)
		noF = mxen_at.formula_composition[el_adsorb]

	# uHF_related
		if el_adsorb == "F":

			uHF_x, ulimit, l_limit = create_uHF(self.ecomp, self.eEl) 

		elif el_adsorb == "Cl":
			
			uHF_x, ulimit, l_limit = create_uHCl(self.ecomp, self.eEl)

		else:

			raise ValueError("analysis for {} as adsorbant is not available".format(el_adsorb))

		zerouHF_x = (e_mxene + e_AF3 + (3+noF)*e_H - e_Max)/(3.0 + noF)
		uHF_x = np.sort(np.append(uHF_x, zerouHF_x))
	
		if verbose:
			print("eMAX: {}".format(e_Max))
			print("emxene: {}".format(e_mxene))
			print("No. F: {}".format(noF))
			print("AelF: {}".format(AelF))

			
		if noF == 0:

			deltaG = e_mxene + e_AF3  +  3*e_H - e_Max - 3*uHF_x

		elif noF == 2:

			deltaG = e_mxene + e_AF3 + 5*e_H - e_Max - 5*uHF_x
		else:

			raise Exception("%s are unexpected" %noF)

		return uHF_x, deltaG, ulimit, l_limit

	def orthomxene_stable(self, orthomxen_at, reaction_type, el_adsorb="F", out_specie="H"):
		"""Calculates Raction enthalpy of orthomxene(M2AlC2) formation from MAX phase
		orthomxen_at: 'M2AC2', atoms instance,
		Note: The M and A assumes to be forming MF3 and AF3 here.
		reaction_type: str, available possible reactions are following:
		a:
			2Ti2AlC + (3+x)HF ---> Ti2AlC2F_x + 2Ti + AlF3 + (3+x)H
			3/2Ti3AlC2 + (3/2+x)HF ---> Ti3AlC3F2 + 3/2Ti + 1/2AlF3 + (3/2 + x)H

		b:
			2Ti2AlC + (9+x)HF ---> Ti2AlC2F_x + 2TiF3 + AlF3 + (9+x)H
			3/2Ti3AlC2 + (6+x)HF ----> Ti3AlC3F2 + 3/2TiF3 + 1/2AlF3 + (6+x)H
		Can handle only two types of reactions and thus is not general.

		General expression in terms of 'n'ToDo: Need to implement this:
		a:
			[(n+1)/n]M(n+1)ACn + (3/n + x)HF ----> M(n+1)AC(n+1)Fx + [(n+1)/n]M + (1/n)AF3 + (3/n + x)H

		b:
			[(n+1)/n]M(n+1)ACn + (3(n+2)/n + x)HF ---> M(n+1)AC(n+1)Fx + [(n+1)/n]MF3 + (1/n)AF3 + (3(n+2)/n + x)H


		Returns uHF array, deltaG array, upper and lower limits of uHF"""

		from energy_calculations import create_uHF, create_uHCl
		from preprocessing.util import remove_1, float_int



		noF = orthomxen_at.formula_composition[el_adsorb] # no of adsorb atoms

		Max_mapp = self._max_at.get_mappings() # get M-A-X mapping of elements
		sys_comps = self.get_selection_compounds()

		#energies
		e_orthomxen = orthomxen_at.get_energy_fu()
		e_max = self.max_at.get_energy_fu()
		e_af = self.ecomp[sys_comps["A%s" % el_adsorb]]
		e_mf = self.ecomp[sys_comps["M%s" % el_adsorb]]
		e_m = self.eEl[Max_mapp["M"]]
		e_h = self.eEl[out_specie]
		
		print("Energies")
		print("ortho:{}, e_max:{}, e_af:{}, e_mf:{}, e_m:{}, e_h:{}".format(e_orthomxen, e_max, e_af, e_mf, e_m, e_h))
		
		f_mf = re.search('(?<=%s)[0-9]*.?[0-9]*' %el_adsorb,  sys_comps["M%s" %el_adsorb]).group()
		if not f_mf or f_mf != "3":
				raise ValueError('Reactions involving {} can\'t be handled at the momet; only MF3 types are possible'.format(sys_comps["M%s" %el_adsorb]))

		f_af = re.search('(?<=%s)[0-9]*.?[0-9]*' %el_adsorb, sys_comps["A%s" %el_adsorb]).group()
	
		if not f_af or f_af != "3":
				
				raise ValueError('Reactions involving {} can\'t be handled at the momet; only AF3 types are possible'.format(sys_comps["A%s" %el_adsorb]))


		formula = orthomxen_at.formula_composition.iupac_formula
		symbolic_formula = formula.replace(" ", "")
		print("Formula is: {}".format(symbolic_formula))

		symbolic_formula = symbolic_formula.replace(el_adsorb+"%s"% float_int(noF) , "")
		for i,j in orthomxen_at.get_mappings().items():
			symbolic_formula = symbolic_formula.replace(j, i)

		#removing adsorbspecie

		symbolic_formula = remove_1(symbolic_formula)
		print(symbolic_formula)
		#x-axis
		if el_adsorb == "F":

			uHF_x, ulimit, l_limit = create_uHF(self.ecomp, self.eEl)

		elif el_adsorb == "Cl":
			
			uHF_x, ulimit, l_limit = create_uHCl(self.ecomp, self.eEl)

		else:

			raise ValueError("Analysis for {} as adsorbant is not available".format(el_adsorb))
		
		if reaction_type == "a":

				mainlegend = "MAX2ene-a"
				if "M2AX2"  == symbolic_formula:

					zerouHF_x = (e_orthomxen + e_af + 2*e_m + (3+noF)*e_h - 2*e_max)/(3 + noF)
					uHF_x = np.sort(np.append(uHF_x, zerouHF_x))
					deltaG = e_orthomxen + e_af + 2*e_m + (3 + noF)*e_h - 2*e_max - (3 + noF)*uHF_x

					deltaG = deltaG/2

				elif "M3AX3" == symbolic_formula:
					zerouHF_x = (e_orthomxen + e_af/2 + 3*e_m/2 + (3/2 + noF)*e_h - 3*e_max/2)/(3/2 + noF)
					uHF_x = np.sort(np.append(uHF_x, zerouHF_x))
					

					deltaG = e_orthomxen + e_af/2 + 3*e_m/2 + (3/2 + noF)*e_h - 3*e_max/2 - (3/2 + noF)*uHF_x

					deltaG = deltaG*2/3 # dividing by MAX phase mole

				else:

					raise ValueError("deltaG energy calculator for {} compound is not available".format(symbolic_formula))

		elif reaction_type == "b":

				mainlegend = "MAX2ene-b"

				if "M2AX2" == symbolic_formula:

					zerouHF_x = (e_orthomxen + 2*e_mf + e_af + (9 + noF)*e_h - 2*e_max)/(9+noF)
					uHF_x = np.sort(np.append(uHF_x, zerouHF_x))

					deltaG = e_orthomxen + 2*e_mf + e_af + (9 + noF)*e_h - 2*e_max - (9 + noF)*uHF_x
					deltaG = deltaG/2

				elif "M3AX3" == symbolic_formula:

					zerouHF_x = (e_orthomxen + e_af/2 + 3*e_mf/2 + (6 + noF)*e_h - 3*e_max/2)/(6+noF)
					uHF_x = np.sort(np.append(uHF_x, zerouHF_x))

					deltaG = e_orthomxen + e_af/2 + 3*e_mf/2 + (6 + noF)*e_h - 3*e_max/2 - (6+noF)*uHF_x
					deltaG = deltaG*2/3

				else:
					raise ValueError("deltaG energy calculator for {} compound is not available".format(symbolic_formula))
		else:
			raise ValueError("Unexpected value of reaction type {}".format(reaction_type))

		print(Fore.RED+"I m scaling DeltaG to the moles of MAX")
		return uHF_x, deltaG, ulimit, l_limit, mainlegend


	def mxene_general(self, mxene_at, etchant="HF", el_adsorb="F", default_product="H"):
		"""Note: Still in testing phase
		Needs more testing and is likely not as stable as the old method"""

		warnings.warn('{} is still in testing phase.\n It requires intensive testing. Always cross-check the output'.format(self.mxene_general.__name__), cetagory=RuntimeWarning)

		reactants = set([*self.reactants, etchant])
		products = default_product.split(",")

		Max_mapp = self._max_at.get_mappings()
		system_compounds = self.get_selection_compounds()

		products = set([*products, system_compounds["A%s" %el_adsorb]])
		products = products | {mxene_at.formula_composition.reduced_formula}
		energies = {at.formula_composition.reduced_formula : at.get_energy_fu()} #

		return self.get_DeltaG_univerisal(reactants, products, u_xcomp=etchant, **energies)

	def orthomxene_general(self, orthomxen_at, reaction_type, etchant="HF", default_product="H", el_adsorb="F"):

		""" Warning: still in testing phase......
		Todo: Needs more testing; Use Chempy equation balancer
		Calculates Reaction enthalpy of orthomxene (M2AlC2) formation
		orthomxen_at: 'M2AlC2' atoms instance
		reaction_type: str, current possibilities are a, b,

		a:
			2Ti2AlC + HF ----> Ti2AlC2 + 2Ti + 2AlF3+  3/2H2

		b:
			2Ti2AlC + 12HF ----> Ti2AlC2 + 2TiF3 + 2AlF3 + 6H2

		"""

		warnings.warn('{} is still in testing phase.\n It requires intensive testing. Always cross-check the output'.format(self.orthomxene_general.__name__))

		if el_adsorb not in etchant:

			raise ValueError("Found incompatibility of {} and {}".format(etchant, el_adsorb))

		reactants = set([*self.reactants, etchant])
		products = default_product.split()

		Max_mapp = self._max_at.get_mappings()	

		sys_comps = self.get_selection_compounds()

		if reaction_type == "a":
			products = set([*products, sys_comps["A%s" %el_adsorb], Max_mapp["M"]])
		elif reaction_type == "b":
			products = set([*products, sys_comps["M%s" %el_adsorb], sys_comps["A%s" %el_adsorb]])	

		else:
			raise ValueError("The reaction type {} is not valid".format(reaction_type))

		legend = "MAX2ene-%s" % reaction_type
		
		#energies
		ortho_formula = orthomxen_at.formula_composition.reduced_formula
		products = products | {ortho_formula}
		energies = {ortho_formula : orthomxen_at.get_energy_fu()}

		return self.get_DeltaG_universal(reactants, products, u_xcomp=etchant, **energies), legend


	def side_reactions(self, reaction_type, out_specie="H", etchant="HF", anion="F", verbose=True ):
		"""Calculates reaction enthalpy of side phases against MXene formation from MAX phase under a given acidic 'AB' 
		type acidic etchant.
		reaction_type: str either 'a' or 'b', which correspond to following reactions:
		'a':
		
			Ti2AlC + 6HF ---> TiC + TiF3 + AlF3 + 6H
			Ti3AlC2 + 6HF ---> 2TiC + TiF3 + AlF3 + 6H
			2Cr2AlC + 9HF ---> Cr3C2 + CrF3 + 2AlF3 + 9H

		'b':
			Ti2AlC + 3HF ---> TiC + Ti + AlF3 + 3H
			Ti3AlC2 + 3HF ----> 2TiC + Ti + AlF3 + 3H
			2Cr2AlC + 6HF ----> Cr3C2 + Cr + 2AlF3 + 6H

		Internally uses equation balancer provided by chempy"""

		warnings.warn(u"\u001b[38;5;210m" + '{} is still in testing phase.\n It requires intensive testing. Always cross-check the output'.format(self.side_reactions.__name__), category=RuntimeWarning)

		if anion not in etchant:

			raise ValueError("Found incompatibility of {} and {}".format(etchant, anion))				

		reactants = set([*self.reactants, etchant])
		products = out_specie.split(",")

		Max_mapp = self._max_at.get_mappings()
		sys_comps = self.get_selection_compounds()

		if reaction_type == "a":
			products = set([*products, sys_comps["A%s" %anion], sys_comps["M%s" %anion], sys_comps["MC"]])
			mainlegend = [sys_comps[i] for i in ["MC" , "M%s" %anion]]
		elif reaction_type == "b":
			products = set([*products, sys_comps["A%s" %anion], Max_mapp["M"], sys_comps["MC"]])
			mainlegend = [sys_comps["MC"], Max_mapp["M"]]
		else:
			raise ValueError("The reaction type {} is not valid".format(reaction_type))

		mainlegend = " and ".join(mainlegend)
		return self.get_DeltaG_universal(reactants, products, u_xomp=etchant), mainlegend
	
	def _releva_energies(self, reactants, products):
	
		max_formula = self.max_at.formula_composition.reduced_formula
		e_max = self.max_at.get_energy_fu()
		energies = {max_formula : e_max}

		for i in reactants | products:

				

			if i in self.ecomp:

				energies[i] = self.ecomp[i]

			else:
				try:

					energies[i] = self.eEl[i]
				except KeyError:
					pass


		return energies


	def get_DeltaG_universal(self, reactants, products, u_xcomp="HF", verbose=True, **specificeng):

		"""Calculates deltaG of a reaction involving a variable chemical potential of a reactant i.e. HF. It balances the reaction using chempy class: 'balance_stoichiometry', 
		reactants: set of reactants
		products: set of products
		energies: dictionary of compound, energy as key, value respectively
		u_xcomp: (etchant), the compound whose chemical potential acts as a variable i.e. HF or HCl, or alternatively, independent variable x
		anion: element (F, Cl). The anion introducted by the etchant(u_xcomp) in the solution, also reacts with various metallic ions to form respective compound i.e AF3, MF3,"""

		warnings.warn('{} is still in testing phase\n. internally uses chempy package\n Cross check it please'.format(self.get_DeltaG_universal.__name__))

		from energy_calculations import create_uHF, create_uHCl

		energies = self._releva_energies({i for i in reactants if i != u_xcomp}, set(products))
		energies.update(specificeng)
		coeff_reac, coeff_pro = balance_stoichiometry(reactants, products) #Equation balancer
		coeff_reac = OrderedDict( (k, int(v)) if v.is_integer else (k, float(v)) for k,v in coeff_reac.items())
		coeff_pro = OrderedDict( (k, int(v)) if v.is_integer else (k, float(v)) for k,v in coeff_pro.items())
		
		if verbose:
			print("Reaction is:")	
			print(self._displayreaction(reactants, products, coeff_reac, coeff_pro))
			print("Energies :")
			print(energies)

		if u_xcomp == "HF":

			uHF_x, ulimit, llimit = create_uHF(self.ecomp, self.eEl)

		elif u_xcomp == "HCl":

			uHF_x, ulimit, llimit = create_uHCl(self.ecomp, self.eEl)

		else: 

			raise ValueError("u_xomp(etchant) {} is not available and invalid".format(u_xcomp))

		print(u_xcomp)
		zerouHF_x = sum([coeff_pro[p]*energies[p] for p in products]) - sum([coeff_reac[r]*energies[r] for r in reactants if u_xcomp != r])
		zerouHF_x = zerouHF_x/(coeff_reac[u_xcomp])
		uHF_x = np.sort(np.append(uHF_x, zerouHF_x))

		deltaG = sum([coeff_pro[p]*energies[p] for p in products]) - sum([coeff_reac[r]*energies[r] for r in reactants if u_xcomp != r]) - coeff_reac[u_xcomp]*uHF_x
		
		print(Fore.RED + "Scaling the DeltaG by No.of moles of MAX")
		print("Returing DeltaG eV/MAXmol")
		
		max_formula = self.max_at.formula_composition.reduced_formula
		deltaG = deltaG/coeff_reac[max_formula]
		
		return uHF_x, deltaG, ulimit, llimit

	@staticmethod
	def _displayreaction(reactants, products,coeff_reac, coeff_pro):

		reacstr = " + ".join(["{}{}".format(coeff_reac[r], r) for r in reactants])
		productstr = " + ".join(["{}{}".format(coeff_pro[p], p) for p in products])

		return reacstr + " ------> " + productstr

	def get_DeltaG_pymatgen(self, reactants, products, energies, etchant="HF", verbose=True):
		"""Note:Still in testing phase.......
		 DeltaG of the reaction given by reactants and products
		Todo: Needs more testing"""

		warnings.warn('{} is still in testing phase.\n It requires intensive testing. Always cross-check the output'.format(self.get_DeltaG.__name__), cetagory=RuntimeWarning)

		energies[self.max_at.reduced_formula] = self.max_at.get_energy_fu()
		energies.update(self.ecomp)
		energies.update(self.eEl)
		
		products = [Composition(i) for i in products]
		reactants = [Composition(i) for i in reactants]

		ceq = Reaction(reactants, products)

		print("Reaction:")
		print(ceq)
		print(energies)
		print("Coeffs are:")
		print([ceq.get_coeff(i) for i in reactants])
		print([ceq.get_coeff(i) for i in products])
		#related to uHF

		if etchant == "HF":

			uHF, ulimit, l_limit = 	create_uHF(self.ecomp, self.eEl)

			zerouHF_x = sum([ceq.get_coeff(i)*energies[i.reduced_formula] for i in products]) + \
						sum([ceq.get_coeff(j)*energies[j.reduced_formula] for j in reactants if j.reduced_formula != etchant])

			zerouHF_x = zerouHF_x/(-1*ceq.get_coeff(Composition(etchant)))

			uHF = np.sort(np.append(uHF, zerouHF_x))

		else:
			raise ValueError("Etchant is unknown")
		
		if verbose:

			print("Calculating deltaG of reaction:")
			print(ceq)

		deltaG = sum([c*energies[i.reduced_formula] for c,i in zip(ceq.coeffs, ceq.all_comp) if i.reduced_formula != etchant]) \
				 + ceq.get_coeff(Composition(etchant))*uHF
	
		return uHF, deltaG, ulimit, l_limit


	def get_selection_compounds(self):

		"""Returns dictionary of MC, MF, AF as keys and values containing compounds/compositions usable for a given reactions involving MAX phase"""

		MAX_map = self.max_at.get_mappings()
		els = sorted(MAX_map.values())

		if els == sorted(["Ti", "Al", "C"]):

			return {"MC": "TiC", "MF": "TiF3", "AF":"AlF3",}


		elif els == sorted(["Cr", "Al", "C"]) or els == sorted(["Cr", "Ga", "C"]):

			return {"MC": "Cr3C2", "AF": "%sF3" % MAX_map["A"], "MF": "CrF3"
				,"ACl": "%sCl3" % MAX_map["A"], "MCl": "CrCl3"}

		elif els == sorted(["V", "Al", "C"]):

			return {"MC": "V6C5", "AF": "AlF3", "MF": "VF3"}

		elif els == sorted(["Nb", "Al", "C"]):

			return {"MC": "Nb6C5", "AF": "AlF3", "MF": "NbF5"}

		elif els == sorted(["Zr", "Al", "C"]):

			return {"MC": "ZrC", "AF": "AlF3", "MF": "ZrF4"}

		else:

			raise ValueError("{} is not listed under given compounds".format(MAX_map))


class Formation_static:
	
	def __init__(self, reactants, products):

		self.reactants = reactants
		self.products = products

	@property
	def reactants(self):
		return self._reactants
	@property
	def products(self):
		return self._products

	@reactants.setter
	def reactants(self, value):
		self._reactants = value

	@products.setter
	def products(self, value):
		self._products = value

	def reaction_enthalpy(self, energies, coeffs):

		print(self.displayreaction(self._reactants, self._products, coeffs))
		deltaG = sum(coeffs[p]*energies[p] for p in self._products) - sum(coeffs[r]*energies[r] for r in self._reactants)
		return deltaG

	@staticmethod
	def generate_coeffs(names, coeffs):
		return {i : c for i, c in zip(names, coeffs)}

	@staticmethod	
	def displayreaction(reactants, products, coeffs):

		reacstr = " + ".join(["{}{}".format(coeffs[r], r) for r in reactants])
		productstr = " + ".join(["{}{}".format(coeffs[p], p) for p in products])

		return reacstr + "-------->" + productstr


class Refenergies_dict(UserDict):

	"""Dictioanary object for reference energies and etchant-specific methods"""

	def __init__(self, criteria="notstrict", *args, **kwargs):

		super().__init__(*args, **kwargs)
		self._extend_(criteria)	

	def __getattr__(self, key):

		if key in self:

			return self[key]

		else:

			for i in self.keys():

				try:
					return self[i][key]
				
				except KeyError:
					continue
				
			raise AttributeError("{0} does not have the attribute {1}".format(self, key) )

	def _extend_(self, criteria):

		from data_functions import DFT_energies_001, DFT_energies

		
		if criteria is "notstrict":
			ecomp, eEl = DFT_energies()

		elif criteria is "strict":
			ecomp, eEl = DFT_energies_001()

		else:

			raise Exception("Unhandled criteria entered")

		self["ecomp"] = ecomp
		self["eEl"] = eEl

		self.calculate_etchantenergies()		

	def calculate_etchantenergies(self):

		from energy_calculations import get_eHF, get_eH2O, get_eOH, get_eHCl

		self.ecomp["HF"] = get_eHF(self.eEl)
		self.ecomp["HCl"] = get_eHCl(self.eEl)
		self.ecomp["H2O"] = get_eH2O(self.eEl)
		self.ecomp["OH"] = get_eOH(self.eEl)


	def get_reduced_dict(self, attribute, *keys):

		try:

			return	{i: self[attribute][i] for i in keys}

		except KeyError:

			raise KeyError("{} does not have the keys {} ".format(attribute, keys))

class side_reactions:

	"""A class for finding out enthalpies of side ractions. Currently only two possible reaction types are available """

	def __init__(self, Formation_enthalpy):

		self.max_at = Formation_enthalpy.max_at
		self.ecomp = Formation_enthalpy.ecomp
		self.eEl = Formation_enthalpy.eEl
		self.get_selection_compounds = Formation_enthalpy.get_selection_compounds

	def deltaG(self, reaction_type, el_adsorb="F"):
		""" Reaction_type:
		a:
				
			Ti2AlC + 6HF ---> TiC + TiF3 + AlF3 + 6H
			Ti3AlC2 + 6HF ---> 2TiC + TiF3 + AlF3 + 6H

		b:

			Ti2AlC + 3HF ---> TiC + Ti + AlF3 + 3H"""
	
		from energy_calculations import create_uHF, create_uHCl


		max_formula = self.max_at.formula_composition.refined_iupac_formula
		print(max_formula)
		max_n = self.max_at.formula_composition["C"]
		
		Max_mapp = self.max_at.get_mappings()
		Ael = Max_mapp["A"]
		sys_comps = self.get_selection_compounds()
		
		#energies
		e_max = self.max_at.get_energy_fu()
		ecomp = self.ecomp
		eEl = self.eEl	

		#x-axis, uHF
		if el_adsorb == "F":
	
			uHF_x, ulimit, llimit = create_uHF(self.ecomp, self.eEl)

		elif el_adsorb == "Cl":
			
			uHF_x, ulimit, llimit = create_uHCl(self.ecomp, self.eEl)

		else:

			raise ValueError("analysis for {} as adsorbant is not available".format(el_adsorb))
	

		#using laymayn approach
		if reaction_type == "a":
			
			if max_formula == "Ti2AlC" or max_formula  == "Ti3AlC2":	
				zerouHF_x = (max_n*ecomp["TiC"] + ecomp["Ti%s3" %el_adsorb] + ecomp["Al%s3" %el_adsorb] + 6*eEl["H"] - e_max)/6.0
				uHF_x = np.sort(np.append(uHF_x, zerouHF_x))	

				deltaG = max_n*ecomp["TiC"] + ecomp["Ti%s3" %el_adsorb] + self.ecomp["Al%s3" %el_adsorb] + 6*self.eEl["H"] - e_max - 6*uHF_x
				legend = "TiC and Ti%s3" %el_adsorb
	
			elif max_formula == "Cr2AlC" or max_formula  == "Cr2GaC":

				zerouHF_x = (ecomp["Cr3C2"] + ecomp["Cr%s3" %el_adsorb] + 2*ecomp["%s%s3" %(Ael, el_adsorb)] + 9*eEl["H"] - 2*e_max)/9.0
				uHF_x = np.sort(np.append(uHF_x, zerouHF_x))

				deltaG = ecomp["Cr3C2"] + ecomp["Cr%s3" %el_adsorb] + 2*ecomp["%s%s3" %(Ael, el_adsorb)] + 9*eEl["H"] - 2*e_max - 9*uHF_x
				legend = "Cr3C2 and Cr%s3" %el_adsorb

				print(Fore.RED + "Aplying scaling to deltaG")
				deltaG = deltaG/2.0

			elif max_formula == "V2AlC":
				zerouHF_x = (ecomp["V6C5"] + 4*ecomp["V%s3" %el_adsorb] + 5*ecomp["Al%s3" %el_adsorb] + 27*eEl["H"] - 5*e_max)/27
				uHF_x = np.sort(np.append(uHF_x, zerouHF_x))
				
				deltaG =  ecomp["V6C5"] + 4*ecomp["V%s3" %el_adsorb] + 5*ecomp["Al%s3" %el_adsorb] + 27*eEl["H"] - 5*e_max - 27*uHF_x
				legend = "V6C5 and V%s3" %el_adsorb
				print(Fore.RED + "Im applying scaling to deltaG- to have same no of MAX moles")
				deltaG = deltaG/5.0

			else:
	
				raise ValueError("{} formation energy calculator is not available".format(max_formula))

		
		elif reaction_type == "b":

			if max_formula == "Ti2AlC" or max_formula == "Ti3AlC2":
				zerouHF_x = (max_n*ecomp["TiC"] + eEl["Ti"] + ecomp["Al%s3" %el_adsorb] + 3*eEl["H"] - e_max)/3.0
				uHF_x = np.sort(np.append(uHF_x, zerouHF_x))

				deltaG = max_n*ecomp["TiC"] + eEl["Ti"] + ecomp["Al%s3" %el_adsorb] + 3*eEl["H"] - e_max - 3*uHF_x
				legend = "TiC and Ti"
	
			elif max_formula == "Cr2AlC" or max_formula == "Cr2GaC":
				zerouHF_x = (max_n*ecomp["Cr3C2"] + eEl["Cr"] + 2*ecomp["%s%s3" %(Ael, el_adsorb)] + 6*eEl["H"] - 2*e_max)/6.0
				uHF_x = np.sort(np.append(uHF_x, zerouHF_x))
	
				deltaG = ecomp["Cr3C2"] + eEl["Cr"] + 2*ecomp["%s%s3" % (Ael,el_adsorb)] + 6*eEl["H"] - 2*e_max - 6*uHF_x
				legend = "Cr3C2 and Cr"
				print(Fore.RED + "Aplying scaling to deltaG")
				deltaG = deltaG/2.0

			elif max_formula == "V2AlC":

				zerouHF_x = (ecomp["V6C5"] + 4*eEl["V"] + 5*ecomp["Al%s3" %el_adsorb] + 15*eEl["H"] - 5*e_max)/15
				uHF_x = np.sort(np.append(uHF_x, zerouHF_x))
				
				deltaG = ecomp["V6C5"] + 4*eEl["V"] + 5*ecomp["Al%s3" %el_adsorb] + 15*eEl["H"] - 5*e_max - 15*uHF_x
				legend = "V6C5 and V"
				print(Fore.RED + "Aplying scaling to deltaG")

				deltaG = deltaG/5.0

		elif reaction_type == "c":
			
			if max_formula == "V2AlC":

				zerouHF_x = (ecomp["V2C"] + ecomp["Al%s3" %el_adsorb] + 3*eEl["H"] - e_max)/3.0
				uHF_x = np.sort(np.append(uHF_x, zerouHF_x))

				deltaG = ecomp["V2C"] + ecomp["Al%s3" %el_adsorb] + 3*eEl["H"] - e_max - 3*uHF_x
				legend = "V2C"


		return uHF_x, deltaG, ulimit, llimit, legend


#	def zerouHF_a(self, Ael):
#		
#		if self.maxphase == "Ti2AlC" or self.maxphase == "Ti3AlC2":
#
#			zerouHF = (self.n*self.ecomp["TiC"] + self.ecomp["TiF3"] + self.ecomp["AlF3"] + 6*self.eEl["H"] - self.emaxphase)/6.0
#
#		elif self.maxphase == "Cr2AlC" or self.maxphase == "Cr2GaC":
#			
#	                zerouHF = (self.ecomp["Cr3C2"] + self.ecomp["CrF3"] + 2*self.ecomp["%sF3" %Ael] + 9*self.eEl["H"] - 2*self.emaxphase)/9.0
#
#		return zerouHF
#
#
#	def zerouHF_b(self, Ael):
#
#		if self.maxphase == "Ti2AlC" or self.maxphase == "Ti3AlC2":
#
#			zerouHF = (self.n*self.ecomp["TiC"] + self.eEl["Ti"] + self.ecomp["AlF3"] + 3*self.eEl["H"] - self.emaxphase)/3.0
#
#		elif self.maxphase == "Cr2AlC" or self.maxphase == "Cr2GaC":
#			
#			zerouHF = (self.ecomp["Cr3C2"] + self.eEl["Cr"] + 2*self.ecomp["%sF3" %Ael] + 6*self.eEl["H"] - 2*self.emaxphase)/6.0
#
#		return zerouHF



class Formation_maxphases(Formation_enthalpy):

	def __init__(self, max_ats, Refenergies_dict, select_ecomp=["H2O", "HF", "OH"], select_eEl=["H", "O", "F"],  default_reactants=[]):

		try:
			super().__init__( max_ats, Refenergies_dict, select_ecomp, select_eEl, default_reactants)
		
		except AttributeError:
			self.reactants = default_reactants

		self.max_ats = max_ats
		del self._max_at

	@property
	def max_ats(self):
		return self._max_ats

	@max_ats.setter
	def max_ats(self, value):

		if isinstance(value, abc.Iterable):

			self._max_ats = value

		else:

			raise TypeError("'{}' is not an iterable".format(value))

	def imxenes_reactions_stable(self, mxen_ats, el_adsorb="F"):


		for i, mx_at in enumerate(mxen_ats):

			# maxphase check
			mname, mindex = self.selectmax_fromidentitystring(mx_at.path)
			maxat = self.max_ats[mindex]
	
			max_f = maxat.formula_composition.iupac_formula
			max_Ael = maxat.get_Ael() # Need to change it for solid solutions
			
			if re.sub(max_Ael + "\d?\.?\d?", "", max_f) != mx_at.formula_composition.get_withoutel_iupac_formula(el_adsorb):

				raise ValueError(" max phase '{}' is not the parent phase of the mxene '{}'. Probably the sequence of max and mxene lists/tuples are incompatible".format(maxat, mxen_ats[i]) )

			self.max_at = maxat
			print(Fore.RED + "MAx phase {}".format(self.max_at))
			print(Fore.MAGENTA + "Mxene phase {}".format(mx_at))

			yield super().mxene_reaction_stable(mx_at, el_adsorb), mname, mindex

				

	def selectmax_fromidentitystring(self, string):

		from preprocessing.util import remove_1

		
		max_fs= [i.formula_composition.iupac_formula for i in self.max_ats]
		max_fs = [remove_1(i) for i in max_fs]

		maxphase = [(name, index) for index, name in enumerate(max_fs) if re.match(name, string)]

		if len(maxphase) >= 2:

			raise ValueError("Search of a max phase in '{}' found more than one max phases '{}', bug or string is incorrect".format(string, maxphase ))
		elif not maxphase:
			raise ValueError("Could not match a max phases '{}' to a string '{}'".format(max_fs, string))

		return maxphase[0]
	

class Chemicalreaction(object):
		"""Class for dealing with and inspecting reactions."""

		def __init__(self, reactants, products):
			self.reactants = reactants
			self.products = products

		@property
		def reactants(self):
			return self._reactants
		@reactants.setter
		def reactants(self, value):
			self._reactants = value
			self._react_elements = self.get_el_lists(value)

		@property
		def products(self):
			return self._products
		@products.setter
		def products(self, value):
			self._products  = value
			self._product_elements = self.get_el_lists(value)

		@property
		def product_elements(self):
			return self._product_elements

		@property
		def react_elements(self):
			return	self._react_elements

		@staticmethod
		def get_elements(str):
			"""returns the elements present in the string.
			It is based on string parsing function i.e. regrex matching method.
			Intended to be used by higher order methods."""
			return re.findall("[A-Z][a-z]?]", str)

		def get_el_lists(self, lst):
			"""Higher level method that returns the elements present in the list of strings.
			Internally, calls the get_elements method"""
			el = map(self.get_elements, lst)
			return [j for i in el for j in i] # keeping it unsorted now

		def get_all_compositions(self):
			"""Returns all the compositions of both reactants and products"""
			return	[i for i in self._reactants + self._products]

		@staticmethod
		def count(str):
				"""Returns a dictionary with counts of each element in the given
				'str' composition. Internally, it is based on the ase.formula.Formula
				object"""

				f = Formula(str)
				return f.count()

		def inspect(self, coeffs):
			"""Returns True if the reaction is balanced otherwise False.
			 coeffs: A nested dict with 'reactants' and 'products' as keys"""

			allreactants = self.get_elements(coeffs["reactants"], "reactants")
			allproducts = self.get_elements(coeffs["products"], "products")

			return r_el_counts == p_el_counts

		def get_counts_reaction(self, coeffs, kind):
			"""Returns a dict of all elements constituting all the 'reactants' or 'products'
			controlled by the argument kind."""
			species = "".join(getattr(self, kind))
			r_counts = self.count(species)
			return {i: coeffs*value for i, value in r_counts.items()}
