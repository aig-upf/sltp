import sys
from myTask import MyTask
from timeit import default_timer as timer
import os

class MyCNFError(Exception):
  def __init__(self, value):
    self.value = value
  def __str__(self):
    return repr(self.value)

class CNF:
	type1 = 'Atom-Controller'
	type2 = 'Action-Controller'
	type3 = 'Triplet'
	type4 = 'Reachable-I'
	type5 = 'Reachable-G'
	type6 = 'Replacement-Goal'
	type7 = 'Controller-Controller'
	type8 = 'Replacement-Equality'
	type9 = 'Inequality-CSCS'
	type10 = 'Replacement-Goal'
	num_types = 18
	print_types = [11, 1, 2, 3, 7]

	def __init__(self, n_file, n_file_extra, fair, strong):
		self.disjunctions = [] # list of disjunctions
		self.maxKey = 1
		self.mapVariableNumber = {}
		self.mapNumberVariable = {}
		self.mapVariableType = {}
		self.clauseSizeCounter = {}
		self.name_file_formula = n_file
		self.name_file_formula_extra = n_file_extra
		self.file_formula = open(n_file, 'w')
		self.file_formula_extra = open(n_file_extra, 'w')
		self.file_formula.close()
		self.file_formula_extra.close()
		# open and close to delete the contents of the files
		self.number_clauses = 0
		self.fair   = fair
		self.strong = strong

	def reset(self):
		self.disjunctions = [] # list of disjunctions
		self.maxKey = 1
		self.mapVariableNumber = {}
		self.mapNumberVariable = {}
		self.mapVariableType = {}
		self.clauseSizeCounter = {}
		self.file_formula = open(self.name_file_formula, 'w')
		self.file_formula.write('p cnf 1 1\n')
		self.file_formula_extra = open(self.name_file_formula_extra, 'a')
		# File formula extra is not used, can be ignored
		self.number_clauses = 0

	###########################################
	################# GENERAL #################
	###########################################

	def generateAtomControllerState(self, atom, controllerState):
		var = atom + '(' + controllerState + ')'
		self.assignKey(var, 1)
		return var

	def generatePairActionControllerState(self, action, controllerState):
		var = '(' + controllerState + ',' + action + ')'
		self.assignKey(var, 2)
		return var

	def generateTripletCSACS(self, initialState, action, finalState):
		var = '(' + initialState + ',' + action + ',' + finalState + ')'
		self.assignKey(var, 3)
		return var

	def generateReachableI(self, controllerState):
		var = 'reachableI(' + controllerState + ')'
		self.assignKey(var, 4)
		return var

	def generateReachableI2(self, controllerState, j):
		var = 'reachableI(' + controllerState + ',' + str(j) + ')'
		self.assignKey(var, 4)
		return var

	def generateReachableG(self, controllerState, j):
		var = 'reachableG(' + controllerState + ',' + str(j) + ')'
		self.assignKey(var, 5)
		return var

	def generateReplacementGoalPropagation(self, controllerState1, controllerState2, i):
		var = 'YR1-' + controllerState1 + '-' + controllerState2 + '-' + str(i)
		self.assignKey(var, 6)
		return var

	def generateReplacementGoalPropagation3(self, controllerState1, controllerState2, i):
		var = 'YR1-FAIR-' + controllerState1 + '-' + controllerState2 + '-' + str(i)
		self.assignKey(var, 6)
		return var

	def generatePairCSCS(self, n1, n2):
		var = '(' + n1 + ',' + n2 + ')'
		self.assignKey(var, 7)
		return var

	def generateReplacementEquality(self, n1, n2, atom):
		var = 'YR2-' + n1 + '-' + n2 + '-' + atom
		self.assignKey(var, 8)
		return var

	def generateInequalityN(self, n1, n2): # n1 < n2
		var = n1 + '<' + n2
		self.assignKey(var, 9)
		return var	

	def generateReplacementIPropagation(self, controllerState1, controllerState2, i):
		var = 'YR3-' + controllerState1 + '-' + controllerState2 + '-' + str(i)
		self.assignKey(var, 10)
		return var

	def generateFirstG(self, n):
		var = 'FirstG(' + str(n) + ')'
		self.assignKey(var, 13)
		return var

	def generateAfterG(self, n):
		var = 'AfterG(' + str(n) + ')'
		self.assignKey(var, 14)
		return var

	def generateConn(self, CStates):
		var = 'conn('
		for n in CStates:
			var += str(n) + ','
		var += ')'
		self.assignKey(var, 15)
		return var

	def generatePairFairCS(self, n):
		var = 'F(' + n + ',fair)'
		self.assignKey(var, 16)
		return var

	def generatePairUnfairCS(self, n):
		var = 'U(' + n + ',unfair)'
		self.assignKey(var, 17)
		return var

	def generateLowerPredecessor(self, n1, n2):
		var = 'Lower(' + n1 + ', ' + n2 + ')'
		self.assignKey(var, 18)
		return var

	def generateInputSat(self, nameFile):
		self.file_formula.close()
		with open(nameFile, 'r') as formula:
			name_final = nameFile + 'header'
			with open(name_final, 'w') as final_formula:
				final_formula.write('p cnf %i %i\n' % (len(self.mapNumberVariable), self.number_clauses))
				first_line = True
				for line in formula:
					if first_line:
						first_line = False
						continue
					final_formula.write(line)
		# # name_final is the name of the file that contains the formula with header and everything
		return name_final

	def generateInputSatWSB(self, nameFile):
		self.file_formula.close()
		with open(nameFile, 'r') as formula:
			name_final = nameFile + 'header'
			with open(name_final, 'w') as final_formula:
				final_formula.write('p cnf %i %i\n' % (len(self.mapNumberVariable), self.number_clauses))
				first_line = True
				for line in formula:
					if first_line:
						first_line = False
						continue
					final_formula.write(line)
		# # name_final is the name of the file that contains the formula with header and everything
		cmd = './shatter.pl ' + name_final
		os.system(cmd)
		return name_final + '.Sym.cnf'

	def writeDisjunctions(self, file):
		for i in self.disjunctions:
			for j in i:
				if j[0] == '-':
					file.write('-' + str(self.mapVariableNumber[j[1:]]) + '\t')
				else:
					file.write(str(self.mapVariableNumber[j]) + '\t')
			file.write('0\n')

	def printDisjunctions(self):
		for i in self.disjunctions:
			for j in i:
				print j + ' ',
			print '\n',

	def printDisjunctionsNumbers(self, printExpanded):
		for i in self.disjunctions:
			if printExpanded:
				for j in i:
					print j + '\t',
				print ''
			for j in i:
				if j[0] == '-':
					print '-' + str(self.mapVariableNumber[j[1:]]), '\t',
				else:
					print self.mapVariableNumber[j], '\t',
			print ' 0'

	def printVariables(self):
		for i in self.mapVariableNumber:
			print i

	def parseOutput(self, nameFile, controllerStates, parser = None):
		fres = open(nameFile, 'r')
	 	res = fres.readlines()
	 	if 'UNSAT' in res[0]:
	 		return False
	 	print 'Solved with', len(controllerStates), 'states'
	 	return True

	def parseOutputLin(self, nameFile, controllerStates, parser = None):
		fres = open(nameFile, 'r')
	 	res = fres.readlines()
	 	for line in res:
	 		if 's UNSATISFIABLE' in line:
	 			return False
	 		if 's SATISFIABLE' in line:
	 			return True
	 	return False

	def parseOutput2(self, nameFile, controllerStates, parser):
		# print_types = [11, 1, 2, 3, 7] # 11: (cs,arg) 1: p(n) 2: (cs, act) 3: (cs, a, cs) 7: (n1, n2)
		sets = [set([]) for i in xrange(self.num_types)]
		fres = open(nameFile, 'r')
		res = fres.readlines()
		if 'UNSAT' in res[0]:
			return False
		res = res[1]
		res = res.split(' ')
		for r in res:
			if '\n' in res:
				continue
			var = int(r)
			if var > 0:
				varName = self.mapNumberVariable[var]
				t = self.mapVariableType[varName]
				sets[t - 1].add(varName)
		for i in xrange(len(sets)):
			if i + 1 in self.print_types:
				s = sets[i]
				print '=========================='
				if i == 0:
					# pair atom controller
					for n in controllerStates:
						print '----------'
						for j in s:
							ind = '(' + n + ')'
							if ind in j:
								print parser.get_var_string(j.split(ind)[0]), ind
				elif i == 1:
					# pair cs action
					for n in controllerStates:
						for j in s:
							if '(' + n + ',' in j:
								print j
				elif i == 2:
					# Triplet
					for n in controllerStates:
						for j in s:
							if '(' + n + ',' in j:
								print j
				elif i == 3:
					# reachable I
					for n in controllerStates:
						print '----------'
						for j in s:
							if n in j:
								print j
				elif i == 4:
					# reachable G
					for n in controllerStates:
						print '----------'
						for j in s:
							if n in j:
								print j
				elif i == 10:
					# cs arguments
					for n in controllerStates:
						for j in s:
							if '(' + n + ',' in j:
								print j
				else:
					for j in s:
						print j
		print 'Solved with', len(controllerStates), 'states'
		return True

	def getNumberVariables(self):
		return len(self.mapVariableNumber)

	def getNumberClauses(self):
		return self.number_clauses

	def printMapVarNumber(self):
		for i in self.mapVariableNumber:
			print i, '-->', self.mapVariableNumber[i]

	def alreadyUsed(self, var):
		if var in self.mapVariableNumber:
			return True
		return False

	def assignKey(self, var, typeVar = -1):
		if not self.alreadyUsed(var):
			self.mapVariableNumber[var] = self.maxKey
			self.mapNumberVariable[self.maxKey] = var
			self.mapVariableType[var] = typeVar
			self.maxKey += 1

	def get_num_cl_vars(self):
		return self.getNumberClauses(), self.getNumberVariables()

	def addClause(self, clause):
		self.number_clauses += 1
		for j in clause:
			negated = (j[0] == '-')
			if negated:
				var = j[1:]
			else:
				var = j
			# self.assignKey(var)
			if negated:
				self.file_formula.write('-' + str(self.mapVariableNumber[var]) + '\t')
			else:
				self.file_formula.write(str(self.mapVariableNumber[var]) + '\t')
		self.file_formula.write('0\n')

	def addClauseExtra(self, clause):
		for j in clause:
			self.file_formula_extra.write(j + '|')
		self.file_formula_extra.write('\n')
		# j1|j2|...|jn|\n

	def printClausesSizes(self, n):
		print self.clauseSizeCounter
		sum_all = 0
		sum_greater = 0
		for i in self.clauseSizeCounter:
			sum_all += self.clauseSizeCounter[i]
			if i >= n:
				sum_greater += self.clauseSizeCounter[i]
			else:
				print i, ':', self.clauseSizeCounter[i]
		print '>=', i, ':', sum_greater
		# print 'sum_all:', sum_all

	###########################################
	############## GENERATION #################
	###########################################

	def generate_clauses(self, task, initialCState, goalCState, controllerStates, k, parser = None, debug = False):
		self.generateInitial(task, initialCState, debug)
		self.generateGoal(task, goalCState, debug)
		self.generatePreconditions(task, controllerStates, debug)
		self.generatePossibleNonDet(task, controllerStates, debug)
		self.generateOneSuccessor(task, controllerStates, debug)
		self.generateTripletForcesBin(task, controllerStates, debug)
		self.generateAtLeastOneAction(task, controllerStates, debug)
		self.generateNegativeForwardPropagation(task, controllerStates, debug)
		self.generateGeneralizeConnection(task, controllerStates, debug)
		self.generateReachableIClauses(task, initialCState, controllerStates, k, debug)
		# Includes I-G propagation
		# 1: Original, no distances, unary
		# 2: With distances, like reachability G
		self.generateReachableGClauses(task, controllerStates, goalCState, k, debug)
		self.generateSymmetryBreakingOnlyOther(task, controllerStates, initialCState, goalCState, debug)
		# self.generateCardinalitySimplification(task, controllerStates, debug)
		self.generateMutexGroupsClauses(task, controllerStates, debug)
		# self.generateH1ContraintsClauses(task, controllerStates, debug)

	###########################################
	################ INITIAL ##################
	###########################################

	def generateInitial(self, task, initialCState, debug = False):
		# -p(n0) for all p not in initial state
		c1, v1 = self.get_num_cl_vars()
		start = timer()
		initial = task.get_initial()
		atoms = task.get_atoms()
		for a in atoms:
			if a not in initial:
				variable = self.generateAtomControllerState(a, initialCState)
				self.addClause(['-' + variable])

		c2, v2 = self.get_num_cl_vars()
		if debug:
			print 'Generation: Initial\t\t v : ', v2 - v1, '\t\t c : ', c2 - c1, '\t', timer() - start

	###########################################
	############## GOAL #######################
	###########################################

	def generateGoal(self, task, goalCState, debug = False):
		# p(ng) for all p in goal state
		c1, v1 = self.get_num_cl_vars()
		start = timer()
		goal = task.get_goal()
		atoms = task.get_atoms()
		for a in atoms:
			if a in goal:
				variable = self.generateAtomControllerState(a, goalCState)
				self.addClause([variable])

		c2, v2 = self.get_num_cl_vars()
		if debug:
			print 'Generation: Goal\t\t v : ', v2 - v1, '\t\t c : ', c2 - c1, '\t', timer() - start

	###########################################
	############## PRECONDITIONS ##############
	###########################################

	def generatePreconditions(self, task, controllerStates, debug = False):
		# (n, b) --> p(n) // if p \in prec(b)
		c1, v1 = self.get_num_cl_vars()
		start = timer()
		atoms = task.get_atoms()
		actions = task.get_actions()
		
		for n in controllerStates:
			for a in actions:
				pre = task.get_preconditions(a)
				var = self.generatePairActionControllerState(a, n)
				for p in pre:
					varPre = self.generateAtomControllerState(p, n)
					self.addClause(['-' + var, varPre])
		 			# print clause + [varPre] # DEBUG

		c2, v2 = self.get_num_cl_vars()
		if debug:
			print 'Generation: Precs\t\t v : ', v2 - v1, '\t\t c : ', c2 - c1, '\t', timer() - start


	###########################################
	############## NON-DET ####################
	###########################################

	def generatePossibleNonDet(self, task, controllerStates, debug = False):
		# 1. (n, act) --> (n, act') // act, act' are non det act of equal det action (action names)
		# 2. (n, act) --> -(n, act'') //
		# 3. (n, b) --> (n, b') // b, b' are siblings
		# 4. (n, b) --> -(n, b'') // b, b'' belong to same act
		c1, v1 = self.get_num_cl_vars()
		start = timer()
		acts = task.get_action_names()
		for n in controllerStates:
			for act in acts:
				var_pair = self.generatePairActionControllerState(act, n)
				other_acts = task.get_other_actions_name(act)
				for act2 in acts:
					if act2 == act:
						continue
					var_pair2 = self.generatePairActionControllerState(act2, n)
					if act2 in other_acts:
						self.addClause(['-' + var_pair, var_pair2]) # 1
					else:
						self.addClause(['-' + var_pair, '-' + var_pair2]) # 2

				for a1 in task.get_actions_with_name(act):
					var1 = self.generatePairActionControllerState(a1, n)
					for a2 in task.get_actions_with_name(act):
						if a2 == a1:
							continue
						if task.actions_are_compatible(a1, a2): # IMPORTANT!! ie. prec not mutex
							var2 = self.generatePairActionControllerState(a2, n)
							self.addClause(['-' + var1, '-' + var2]) # 4

			for a in task.get_actions():
				var1 = self.generatePairActionControllerState(a, n)
				other_actions = task.get_other_actions(a)
				for other in other_actions:
					if other == a:
						continue
					var2 = self.generatePairActionControllerState(other, n)
					self.addClause(['-' + var1, var2]) # 3

		c2, v2 = self.get_num_cl_vars()
		if debug:
			print 'Generation: NonDet\t\t v : ', v2 - v1, '\t\t c : ', c2 - c1, '\t', timer() - start

	###########################################
	############## ONE SUCC ###################
	###########################################

	def generateOneSuccessor(self, task, controllerStates, debug = False):
		# 1. (n, act) --> \OR_{b} (n, b) // for b with name act
		# 2. (n, act) --> \OR_{n'} (n, act, n')
		# 3. -(n, act, n') \lor -(n, act, n'')
		# 4. (n, b) --> (n, act)
		c1, v1 = self.get_num_cl_vars()
		start = timer()
		acts = task.get_action_names()

		for n1 in controllerStates:
			for a_name in acts:
				for n2 in controllerStates:
					for n3 in controllerStates:
						if n3 == n2:
							continue
						var1 = self.generateTripletCSACS(n1, a_name, n2)
						var2 = self.generateTripletCSACS(n1, a_name, n3)
						self.addClause(['-' + var1, '-' + var2]) # 3

				for a in task.get_actions_with_name(a_name):
					pair1 = self.generatePairActionControllerState(a, n1)
					pair2 = self.generatePairActionControllerState(a_name, n1)
					self.addClause(['-' + pair1, pair2]) # 4

				var1 = self.generatePairActionControllerState(a_name, n1)
				var_triplets = []
				for n2 in controllerStates:
					var_triplets.append(self.generateTripletCSACS(n1, a_name, n2))
				self.addClause(['-' + var1] + var_triplets) # 2

				var1 = self.generatePairActionControllerState(a_name, n1)
				var_bin = []
				for a in task.get_actions_with_name(a_name):
					var_bin.append(self.generatePairActionControllerState(a, n1))
				self.addClause(['-' + var1] + var_bin) # 1

		c2, v2 = self.get_num_cl_vars()
		if debug:
			print 'Generation: One succ\t\t v : ', v2 - v1, '\t\t c : ', c2 - c1, '\t', timer() - start

	###########################################
	############## TRIPLET-BIN ################
	###########################################

	def generateTripletForcesBin(self, task, controllerStates, debug = False):
		# (n, act, n') --> (n, act)
		c1, v1 = self.get_num_cl_vars()
		start = timer()
		actions = task.get_action_names()
		for n1 in controllerStates:
			for a in actions:
				for n2 in controllerStates:
					var1 = self.generateTripletCSACS(n1, a, n2)
					var2 = self.generatePairActionControllerState(a, n1)
					# print ['-' + var1, var2] # DEBUG
					self.addClause(['-' + var1, var2])

		c2, v2 = self.get_num_cl_vars()
		if debug:
			print 'Generation: Trip bin\t\t v : ', v2 - v1, '\t\t c : ', c2 - c1, '\t', timer() - start

	###########################################
	############## ONE ACTION #################
	###########################################

	def generateAtLeastOneAction(self, task, controllerStates, debug):
		# \OR_{act_name} (n, act) for all n except goal, ng
		c1, v1 = self.get_num_cl_vars()
		start = timer()
		actions = task.get_action_names()
		for n in controllerStates:
			if n == 'ng':
				continue
			disj = []
			for a in actions:
				pair = self.generatePairActionControllerState(a, n)
				disj.append(pair)
			self.addClause(disj)

		c2, v2 = self.get_num_cl_vars()
		if debug:
			print 'Generation: One act\t\t v : ', v2 - v1, '\t\t c : ', c2 - c1, '\t', timer() - start

	###########################################
	############## NEG-FORWARD PROP ###########
	###########################################

	def __sibling_action_adds_atom(self, task, action, atom):
		siblings = task.get_other_actions(action)
		for a in siblings:
			add_list = task.get_add_list(a)
			if atom in add_list:
				return True
		return False

	def generateNegativeForwardPropagation(self, task, controllerStates, debug):
		# 1. (n, act, n') \land (n, b) --> -p(n') // if p in delete list
		# 2. (n, n') \land -p(n) --> -p(n') \lor \OR_{b: p \in add(b)} (n, b)
		# 3. (n, act, n') \land (n, b) \land -p(n) --> -p(n') // if p \not \in add(b) 
		# and some sibling of b adds p
		c1, v1 = self.get_num_cl_vars()
		start = timer()
		atoms = task.get_atoms()
		actions = task.get_actions()

		#for a in actions:
		#	print(a, task.get_del_list(a))
		#	print(a, task.get_add_list(a))

		for n1 in controllerStates:
			for n2 in controllerStates:
				for p in atoms:
					var_atom_n1 = self.generateAtomControllerState(p, n1)
					var_pair_n1n2 = self.generatePairCSCS(n1, n2)
					var_atom_n2 = self.generateAtomControllerState(p, n2)
					disj_add_clause = ['-' + var_pair_n1n2, var_atom_n1, '-' + var_atom_n2]
					for a in task.get_relevant_actions(p):
						a_name = task.get_action_name(a)
						del_list = task.get_del_list(a)
						add_list = task.get_add_list(a)
						var_triplet = self.generateTripletCSACS(n1, a_name, n2)
						var_bin = self.generatePairActionControllerState(a, n1)
						if p in del_list:
							self.addClause(['-' + var_triplet, '-' + var_bin,'-' + var_atom_n2]) # 1
						if p in add_list:
							disj_add_clause.append(var_bin)
						if p not in add_list and self.__sibling_action_adds_atom(task, a, p):
							self.addClause(['-' + var_triplet, '-' + var_bin, var_atom_n1, '-' + var_atom_n2]) # 3
					self.addClause(disj_add_clause) # 2

		c2, v2 = self.get_num_cl_vars()
		if debug:
			print 'Generation: Neg Prop\t\t v : ', v2 - v1, '\t\t c : ', c2 - c1, '\t', timer() - start

	def generateNegativeForwardPropagation2(self, task, controllerStates, debug):
		# 1. (n, act, n') \land (n, b) --> -p(n') // del list
		# 2. (n, act, n') \land -p(n) --> -p(n') // if no action with name act adds p; else
		# 3. (n, n') \land -p(n) --> -p(n') \lor \OR_{b: p \in add(b)} (n, b)
		# 4. (n, act, n') \land (n, b) \land -p(n) --> -p(n') // if p \not \in add(b) and some sibling of b adds p
		c1, v1 = self.get_num_cl_vars()
		start = timer()
		atoms = task.get_atoms()
		actions = task.get_actions()
		action_names = task.get_action_names()

		for p in atoms:
			for n1 in controllerStates:
				var_atom_n1 = self.generateAtomControllerState(p, n1)
				for n2 in controllerStates:
					var_atom_n2 = self.generateAtomControllerState(p, n2)
					var_pair_cs = self.generatePairCSCS(n1, n2)
					clause = ['-' + var_pair_cs, var_atom_n1, '-' + var_atom_n2]
					for a in task.actions_del_p[p]:
						del_list = task.get_del_list(a)
						act_name = task.get_action_name(a)
						var_triplet = self.generateTripletCSACS(n1, act_name, n2)
						var_bin = self.generatePairActionControllerState(a, n1)
						clause = ['-' + var_triplet, '-' + var_bin, '-' + var_atom_n2]
						if p in del_list:
							self.addClause(clause) # 1
							# print clause # DEBUG

					for a_name in action_names:
						var_triplet = self.generateTripletCSACS(n1, a_name, n2)
						if p in task.atoms_aname_never_adds[a_name]:
							clause = ['-' + var_triplet, var_atom_n1, '-' + var_atom_n2]
							self.addClause(clause) # 2
						else:
							for a in task.get_actions_with_name(a_name):
								add_list = task.get_add_list(a)
								if p not in add_list and self.__sibling_action_adds_atom(task, a, p):
									var_bin = self.generatePairActionControllerState(a, n1)
									clause = ['-' + var_triplet, '-' + var_bin, var_atom_n1, '-' + var_atom_n2]
									self.addClause(clause) # 4

					clause = ['-' + var_pair_cs, var_atom_n1, '-' + var_atom_n2]
					for a in task.actions_add_p[p]:
						add_list = task.get_add_list(a)
						if p in add_list:
							clause.append(self.generatePairActionControllerState(a, n1))
					self.addClause(clause) # 3

		c2, v2 = self.get_num_cl_vars()
		if debug:
			print 'Generation: Neg Prop\t\t v : ', v2 - v1, '\t\t c : ', c2 - c1, '\t', timer() - start

	###########################################
	############## GEN-CONNS ##################
	###########################################

	def generateGeneralizeConnection(self, task, controllerStates, debug = False):
		# (n, n') <--> \OR_{act} (n, act, n')
		c1, v1 = self.get_num_cl_vars()
		start = timer()
		actions = task.get_action_names()
		for n1 in controllerStates:
			for n2 in controllerStates:
				varBin = self.generatePairCSCS(n1, n2)
				triplets = ['-' + varBin]
				for a in actions:
					triplet = self.generateTripletCSACS(n1, a, n2)
					self.addClause(['-' + triplet, varBin])
					# print ['-' + triplet, varBin] # DEBUG
					triplets.append(triplet)
				# print triplets # DEBUG
				self.addClause(triplets)

		c2, v2 = self.get_num_cl_vars()
		if debug:
			print 'Generation: Gen conn\t\t v : ', v2 - v1, '\t\t c : ', c2 - c1, '\t', timer() - start

	###########################################
	############## REACH-I ####################
	###########################################

	def generateReachableIClauses(self, task, initialCState, controllerStates, k, debug = False):
		self.generateReachableIinitial(initialCState, debug)
		# 1
		self.generatePropagationReachableI(task, controllerStates, debug)
		# 1: Uses generalized connections
		self.generatePropagationIG(task, controllerStates, k - 1, debug)
		# 1

	def generateReachableIinitial(self, initialCState, debug = False):
		c1, v1 = self.get_num_cl_vars()
		start = timer()
		self.addClause([self.generateReachableI(initialCState)])
		# print [self.generateReachableI(initialCState)]

		c2, v2 = self.get_num_cl_vars()
		if debug:
			print 'Generation: RI init\t\t v : ', v2 - v1, '\t\t c : ', c2 - c1, '\t', timer() - start

	def generatePropagationReachableI(self, task, controllerStates, debug = False):
		c1, v1 = self.get_num_cl_vars()
		start = timer()
		for n1 in controllerStates:
			for n2 in controllerStates:
				var1 = self.generateReachableI(n1)
				var2 = self.generateReachableI(n2)
				var3 = self.generatePairCSCS(n1, n2)
				self.addClause(['-' + var3, '-' + var1, var2])
				# print ['-' + var3, '-' + var1, var2]

		c2, v2 = self.get_num_cl_vars()
		if debug:
			print 'Generation: RI prop\t\t v : ', v2 - v1, '\t\t c : ', c2 - c1, '\t', timer() - start

	def generatePropagationIG(self, task, controllerStates, k, debug = False):
		c1, v1 = self.get_num_cl_vars()
		start = timer()
		for n in controllerStates:
			var1 = self.generateReachableI(n)
			var2 = self.generateReachableG(n, k)
			self.addClause(['-' + var1, var2])
			# print ['-' + var1, var2]

		c2, v2 = self.get_num_cl_vars()
		if debug:
			print 'Generation: IG prop\t\t v : ', v2 - v1, '\t\t c : ', c2 - c1, '\t', timer() - start

	def generateReachableIClauses2(self, task, initialCState, controllerStates, k, debug = False):
		self.generateReachableIinitial2(task, initialCState, controllerStates, k - 1, debug)
		self.generateCompletionReachabilityI(task, controllerStates, k - 1, debug)
		self.generatePropagationReachableI2(task, controllerStates, k - 1, debug)
		self.generatePropagationIG2(task, controllerStates, k - 1, debug)

	def generatePropagationIG2(self, task, controllerStates, k, debug = False):
		# ReachI(n,k) --> ReachG(n,k)
		c1, v1 = self.get_num_cl_vars()
                start = timer()
		for n in controllerStates:
			var1 = self.generateReachableI2(n, k)
			var2 = self.generateReachableG(n, k)
			self.addClause(['-' + var1, var2])
			# print ['-' + var1, var2]

		c2, v2 = self.get_num_cl_vars()
		if debug:
			print 'Generation: IG prop\t\t v : ', v2 - v1, '\t\t c : ', c2 - c1, '\t', timer() - start

	def generateReachableIinitial2(self, task, initialCState, controllerStates, numberControllerStates, debug = False):
		# ReachI(n0,0), ReachI(n0,1), ...
		# -ReachI(n,0) for all n != n0
		c1, v1 = self.get_num_cl_vars()
                start = timer()
		# set the initial values for initial state
		for i in xrange(numberControllerStates + 1):
			self.addClause([self.generateReachableI2(initialCState, i)])
			# print self.generateReachableI2(initialCState, i)
		# set the negation for non initial states
		for n in controllerStates:
			if n == initialCState:
				continue
			self.addClause(['-' + self.generateReachableI2(n, 0)])
			# print '-' + self.generateReachableI2(n, 0)

		c2, v2 = self.get_num_cl_vars()
		if debug:
			print 'Generation: RI init\t\t v : ', v2 - v1, '\t\t c : ', c2 - c1, '\t', timer() - start

	def generateCompletionReachabilityI(self, task, controllerStates, k, debug = False):
		# ReachI(n,j) --> ReachI(n,j+1)
		c1, v1 = self.get_num_cl_vars()
                start = timer()
		for n in controllerStates:
			for i in xrange(k):
				var1 = self.generateReachableI2(n, i)
				var2 = self.generateReachableI2(n, i + 1)
				self.addClause(['-' + var1, var2])
				# print ['-' + var1, var2]

		c2, v2 = self.get_num_cl_vars()
		if debug:
			print 'Generation: RI compl\t\t v : ', v2 - v1, '\t\t c : ', c2 - c1, '\t', timer() - start

	def generatePropagationReachableI2(self, task, controllerStates, k, debug = False):
		c1, v1 = self.get_num_cl_vars()
                start = timer()
		# Force the equivalences
		for i in xrange(k):
			for n in controllerStates:
				for np in controllerStates:
					if np == 'n0':
						continue
					var1 = self.generateReplacementIPropagation(n, np, i)
					var2 = self.generatePairCSCS(n, np)
					var3 = self.generateReachableI2(n, i)
					self.addClause(['-' + var1, var2])
					self.addClause(['-' + var1, var3])
					self.addClause([var1, '-' + var2, '-' + var3])
		# Right arrow
		for np in controllerStates:
			if np == 'n0':
				continue
			for i in xrange(k):
				var1 = self.generateReachableI2(np, i + 1)
				var2 = ['-' + var1]
				for n in controllerStates:
					var3 = self.generateReplacementIPropagation(n, np, i)
					var2.append(var3)
				self.addClause(var2)
		# Left arrow
		for np in controllerStates:
			if np == 'n0':
					continue
			for i in xrange(k):
				var1 = self.generateReachableI2(np, i + 1)
				for n in controllerStates:
					var2 = self.generateReplacementIPropagation(n, np, i)
					self.addClause([var1, '-' + var2])

		c2, v2 = self.get_num_cl_vars()
		if debug:
			print 'Generation: RI prop\t\t v : ', v2 - v1, '\t\t c : ', c2 - c1, '\t', timer() - start

	###########################################
	############## REACH-G ####################
	###########################################

	def generateReachableGClauses(self, task, controllerStates, goalCState, k, debug = False):
		self.generateReachableGInitial(task, goalCState, controllerStates, k - 1, debug)
		self.generateCompletionReachabilityG(task, controllerStates, k - 1, debug)
		if self.strong:
			self.generatePropagationReachableGStrong(task, controllerStates, k - 1, debug)
		else:
			if not self.fair:
				self.generatePropagationReachableGUnfair(task, controllerStates, k - 1, debug)
			else:
				self.generatePropagationReachableGCyclic(task, controllerStates, k - 1, debug)	

	def generateReachableGInitial(self, task, goalCState, controllerStates, numberControllerStates, debug = False):
		# ReachG(ng,0), ReachG(ng,1), ...
		# -ReachG(n, 0) for n != ng
		c1, v1 = self.get_num_cl_vars()
                start = timer()
		# set the initial values for goal state
		for i in xrange(numberControllerStates + 1):
			self.addClause([self.generateReachableG(goalCState, i)])
		# set the negation for non goal states
		for n in controllerStates:
			if n == goalCState:
				continue
			self.addClause(['-' + self.generateReachableG(n, 0)])

		c2, v2 = self.get_num_cl_vars()
		if debug:
			print 'Generation: RG init\t\t v : ', v2 - v1, '\t\t c : ', c2 - c1, '\t', timer() - start

	def generateCompletionReachabilityG(self, task, controllerStates, k, debug = False):
		# ReachG(n,j) --> ReachG(n, j+1)
		c1, v1 = self.get_num_cl_vars()
                start = timer()
		for n in controllerStates:
			for i in xrange(k):
				var1 = self.generateReachableG(n, i)
				var2 = self.generateReachableG(n, i + 1)
				self.addClause(['-' + var1, var2])

		c2, v2 = self.get_num_cl_vars()
		if debug:
			print 'Generation: RG compl\t\t v : ', v2 - v1, '\t\t c : ', c2 - c1, '\t', timer() - start

	def setFairUnfairActions(self, task, controllerStates):
		# 1: (n, unfair) <-> \OR_{b: unf} (n,b)
		# 2: (n, fair) <-> \OR_{b: unf} (n,b)
		# 3: -(n, fair) \lor -(n, unfair)
		actions = task.get_action_names()
		for n in controllerStates:
			varPairUnf = self.generatePairUnfairCS(n)
			disj = ['-' + varPairUnf]
			for a in actions:
				if '_unfair_' in a:
					varPair = self.generatePairActionControllerState(a, n)
					self.addClause(['-' + varPair, varPairUnf]) # 1
					# print(['-' + varPair, varPairUnf])
					disj.append(varPair)
			self.addClause(disj) # 1
			# print(disj)
		for n in controllerStates:
			varPairF = self.generatePairFairCS(n)
			disj = ['-' + varPairF]
			for a in actions:
				if '_unfair_' not in a:
					varPair = self.generatePairActionControllerState(a, n)
					self.addClause(['-' + varPair, varPairF]) # 2
					# print(['-' + varPair, varPairF])
					disj.append(varPair) 
			self.addClause(disj) # 2
			# print(disj)
		for n in controllerStates:
			varF = self.generatePairFairCS(n)
			varU = self.generatePairUnfairCS(n)
			self.addClause(['-' + varF, '-' + varU]) # 3
			# print(['-' + varF, '-' + varU])

	def generatePropagationReachableGUnfair(self, task, controllerStates, k, debug = False):
		# ReachG(n, j+1) <--> [1] 
		# [1] = [2] \lor [3]
		# [2] = (n, unfair) \land \AND_{n'} [(n,n') --> RG(n',j)]
		# [3] = \OR_{n'} [(n, fair) \land (n,n') \land RG(n',j)]
		c1, v1 = self.get_num_cl_vars()
		start = timer()
		# Set variables (n, fair) and (n, unfair)
		self.setFairUnfairActions(task, controllerStates)
		# Force the equivalences for [2]
		# [(n,n') --> RG(n',j)] <-> Repl(n,n',j)
		for i in xrange(k):
			for n1 in controllerStates:
				if n1 == 'ng':
					continue
				for n2 in controllerStates:
					varRepl = self.generateReplacementGoalPropagation(n1, n2, i)
					varPair = self.generatePairCSCS(n1, n2)
					varRG   = self.generateReachableG(n2, i)
					self.addClause(['-' + varRepl, '-' + varPair, varRG])
					self.addClause(['-' + varRG, varRepl])
					self.addClause([varPair, varRepl])
					# print(['-' + varRepl, '-' + varPair, varRG])
					# print(['-' + varRG, varRepl])
					# print([varPair, varRepl])
		# Force the equivalences for [3]
		# [(n, fair) \land (n,n') \land RG(n',j)] <-> Repl3(n,n',j)
		for i in xrange(k):
			for n1 in controllerStates:
				if n1 == 'ng':
					continue
				for n2 in controllerStates:
					varRepl = self.generateReplacementGoalPropagation3(n1, n2, i)
					varPair = self.generatePairCSCS(n1, n2)
					varFair = self.generatePairFairCS(n1)
					varRG   = self.generateReachableG(n2, i)
					self.addClause(['-' + varRepl, varPair])
					self.addClause(['-' + varRepl, varFair])
					self.addClause(['-' + varRepl, varRG])
					self.addClause([varRepl, '-' + varPair, '-' + varFair, '-' + varRG])
					# print(['-' + varRepl, varPair])
					# print(['-' + varRepl, varFair])
					# print(['-' + varRepl, varRG])
					# print([varRepl, '-' + varPair, '-' + varFair, '-' + varRG])
		# Right arrow
		for n1 in controllerStates:
			if n1 == 'ng':
				continue
			for i in xrange(k):
				varRG = self.generateReachableG(n1, i + 1)
				listStrong = [self.generatePairUnfairCS(n1)]
				listCyclic = []
				for n in controllerStates:
					listStrong.append(self.generateReplacementGoalPropagation(n1, n, i))
					listCyclic.append(self.generateReplacementGoalPropagation3(n1, n, i))
				clause = ['-' + varRG] + listCyclic
				for e in listStrong:
					self.addClause(clause + [e])
					# print(clause + [e])
		# Left arrow
		for n1 in controllerStates:
			if n1 == 'ng':
				continue
			for i in xrange(k):
				varRG = self.generateReachableG(n1, i + 1)
				listStrong = [self.generatePairUnfairCS(n1)]
				listCyclic = []
				for n in controllerStates:
					listStrong.append(self.generateReplacementGoalPropagation(n1, n, i))
					listCyclic.append(self.generateReplacementGoalPropagation3(n1, n, i))
				self.addClause([varRG] + ['-' + e for e in listStrong])
				# print([varRG] + ['-' + e for e in listStrong])
				for e in listCyclic:
					self.addClause([varRG, '-' + e])
					# print([varRG, '-' + e])

		c2, v2 = self.get_num_cl_vars()
		if debug:
			print 'Generation: RG prop\t\t v : ', v2 - v1, '\t\t c : ', c2 - c1, '\t', timer() - start

	def generatePropagationReachableGCyclic(self, task, controllerStates, k, debug = False):
		# ReachG(n, j+1) <--> \OR_{n'} [(n,n') \land ReachG(n', j)]
		c1, v1 = self.get_num_cl_vars()
		start = timer()
		# Force the equivalences
		for i in xrange(k):
			for n1 in controllerStates:
				if n1 == 'ng':
					continue
				for n2 in controllerStates:
					var1 = self.generateReplacementGoalPropagation(n1, n2, i)
					var2 = self.generatePairCSCS(n1, n2)
					var3 = self.generateReachableG(n2, i)
					self.addClause(['-' + var1, var2])
					self.addClause(['-' + var1, var3])
					self.addClause([var1, '-' + var2, '-' + var3])
		# Right arrow
		for n1 in controllerStates:
			if n1 == 'ng':
				continue
			for i in xrange(k):
				var1 = self.generateReachableG(n1, i + 1)
				var2 = ['-' + var1]
				for n2 in controllerStates:
					var3 = self.generateReplacementGoalPropagation(n1, n2, i)
					var2.append(var3)
				self.addClause(var2)
		# Left arrow
		for n1 in controllerStates:
			if n1 == 'ng':
				continue
			for i in xrange(k):
				var1 = self.generateReachableG(n1, i + 1)
				for n2 in controllerStates:
					var2 = self.generateReplacementGoalPropagation(n1, n2, i)
					self.addClause([var1, '-' + var2])

		c2, v2 = self.get_num_cl_vars()
		if debug:
			print 'Generation: RG prop\t\t v : ', v2 - v1, '\t\t c : ', c2 - c1, '\t', timer() - start

	def generatePropagationReachableGStrong(self, task, controllerStates, k, debug = False):
		# ReachG(n, j+1) <--> \AND_{n'} [(n,n') --> ReachG(n', j)]
		c1, v1 = self.get_num_cl_vars()
		start = timer()
		# Force the equivalences
		for i in xrange(k):
			for n1 in controllerStates:
				if n1 == 'ng':
					continue
				for n2 in controllerStates:
					varRepl = self.generateReplacementGoalPropagation(n1, n2, i)
					varPair = self.generatePairCSCS(n1, n2)
					varRG   = self.generateReachableG(n2, i)
					self.addClause(['-' + varRepl, '-' + varPair, varRG])
					self.addClause(['-' + varRG, varRepl])
					self.addClause([varPair, varRepl])
		# Right arrow
		for n1 in controllerStates:
			if n1 == 'ng':
				continue
			for i in xrange(k):
				varRG = self.generateReachableG(n1, i + 1)
				for n2 in controllerStates:
					varRepl = self.generateReplacementGoalPropagation(n1, n2, i)
					self.addClause(['-' + varRG, varRepl])
		# Left arrow
		for n1 in controllerStates:
			if n1 == 'ng':
				continue
			for i in xrange(k):
				disj = [self.generateReachableG(n1, i + 1)]
				for n2 in controllerStates:
					varRepl = self.generateReplacementGoalPropagation(n1, n2, i)
					disj.append('-' + varRepl)
				self.addClause(disj)

		c2, v2 = self.get_num_cl_vars()
		if debug:
			print 'Generation: RG prop\t\t v : ', v2 - v1, '\t\t c : ', c2 - c1, '\t', timer() - start

	###########################################
	############## SYMM-BREAKING ##############
	###########################################

	def generateSymmetryBreaking2(self, task, controllerStates, initialCState, goalCState, debug = False):
		# j > i: reachI(nj,k) --> reachI(ni,k) // CNF: -reachI(nj,k) OR reachI(ni,k)
		c1, v1 = self.get_num_cl_vars()
		start = timer()
		value = len(controllerStates) - 1
		for j, nj in enumerate(controllerStates[:-1]):
			for i, ni in enumerate(controllerStates[:-1]):
				# j has to be greater than i
				if j <= i:
					continue
				for k in xrange(value):
					# we have two controller states (ni, nj), such that j > i
					var1 = self.generateReachableI2(ni, k)
					var2 = self.generateReachableI2(nj, k)
					self.addClause(['-' + var2, var1])

		c2, v2 = self.get_num_cl_vars()
		if debug:
			print 'Generation: Sym Brk\t\t\t v : ', v2 - v1, '\t\t\t c : ', c2 - c1, '\t', timer() - start



	def generateSymmetryBreakingBis(self, task, controllerStates, initialCState, goalCState, debug = False):
		# [1] (ni, aj) --> -(nk, al) for all k<i, l>j
		# if aj is applied in ni, no action greater than aj can be applied in states nk before ni
		# [2] (ni, aj) --> \OR_{l <= j} (nk, al) for k<i
		c1, v1 = self.get_num_cl_vars()
		start = timer()
		
#		orderedActions = task.getOrderedActions()
#		if len(controllerStates) >= 5:
#			reducedCS = controllerStates[1:-1] # remove n0 and nG
#			for ni in reducedCS:
#				for nk in reducedCS:
#					if ni == nk:
#						break # the nk are the controller states such that k < i
#					for j in range(len(orderedActions)):
#						aj = orderedActions[j]
#						for al in orderedActions[j + 1:]:
#							pairni = self.generatePairActionControllerState(aj, ni)
#							pairnk = self.generatePairActionControllerState(al, nk)
#							self.addClause(['-' + pairni, '-' + pairnk]) [1]

		orderedActions = task.getOrderedActions()
		if len(controllerStates) >= 5:
			reducedCS = controllerStates[1:-1] # remove n0 and nG
			for ni in reducedCS:
				for nk in reducedCS:
					if ni == nk:
						break # the nk are the controller states such that k < i
					for j in range(len(orderedActions)):
						aj = orderedActions[j]
						pairni = self.generatePairActionControllerState(aj, ni)
						clause = ['-' + pairni]
						for al in orderedActions[:j + 1]:
							pairnk = self.generatePairActionControllerState(al, nk)
							clause.append(pairnk)
						self.addClause(clause)

			for j in xrange(len(reducedCS) - 1):
				nj   = reducedCS[j]
				njp1 = reducedCS[j + 1]
				for a in orderedActions:
					pairja   = self.generatePairActionControllerState(a, nj)
					pairjp1a = self.generatePairActionControllerState(a, njp1)
					for nl in controllerStates[:-1]:
						pairljp1 = self.generatePairCSCS(nl, njp1)
						clause = ['-' + pairja, '-' + pairjp1a, '-' + pairljp1]
						for ni in controllerStates:
							clause.append(self.generatePairCSCS(ni, nj))
							if ni == nl:
								break
						self.addClause(clause)


		c2, v2 = self.get_num_cl_vars()
		if debug:
			print 'Generation: Sym Brk\t\t v : ', v2 - v1, '\t\t c : ', c2 - c1, '\t', timer() - start		




	def generateSymmetryBreakingOnlyOther(self, task, controllerStates, initialCState, goalCState, debug = False):
		c1, v1 = self.get_num_cl_vars()
		start = timer()

		if len(controllerStates) >= 4:
			# (na, nb) --> \OR_{i<=a} (n_i, n_{b-1})
			for ia, na in enumerate(controllerStates[:-1]):
				for ib, nb in enumerate(controllerStates[:-1]):
					if ib in [0, 1]:
						continue
					nb_1 = controllerStates[ib - 1]
					var_pair_ab = self.generatePairCSCS(na, nb)
					disj = ['-' + var_pair_ab]
					for i in xrange(ia + 1):
						ni = controllerStates[i]
						disj.append(self.generatePairCSCS(ni, nb_1))
					self.addClause(disj)

					#for j in range(ia + 1, len(controllerStates) - 1):
					#	nj = controllerStates[j]
					#	pair_post = self.generatePairCSCS(nj, nb_1)
					#	self.addClause(['-' + var_pair_ab, '-' + var_afterG_nb_1, '-' + pair_post])



			# lower(ni, nj) <--> \AND_{k<j} -(nk, ni) \land (nj, ni)
			# (ni,nj) \land (ni, nk) \land lower(nj,ni) \land lower(nk,ni) \land (ni, a_name)
			# --> (ni, a_name_0, nj) \land (ni, a_name_1, nk) // nj != ni, nj != ng, nk != ni, nk != ng, j < k
			
			# equivalences for lower
			for ni in controllerStates:
				for nj in controllerStates[:-1]:
					lower_ni_nj = self.generateLowerPredecessor(ni, nj)
					pair_nj_ni  = self.generatePairCSCS(nj, ni)
					clause = [lower_ni_nj, '-' + pair_nj_ni]
					self.addClause(['-' + lower_ni_nj, pair_nj_ni]) # right arrow
					for nk in controllerStates:
						if nk == nj: # want k<j
							break
						pair_nk_ni = self.generatePairCSCS(nk, ni)
						clause.append(pair_nk_ni)
						self.addClause(['-' + lower_ni_nj, '-' + pair_nk_ni]) # right arrow
					self.addClause(clause) # left arrow

			# real clause
			if task.max_card == 2:
				for ni in controllerStates[1:-1]:
					for nk in controllerStates[1:-1]:
						lower_nk_ni = self.generateLowerPredecessor(nk, ni)
						if nk == ni:
							continue
						for index_j, nj in enumerate(controllerStates[:-1]):
							if nj == 'n0':
								continue
							lower_nj_ni = self.generateLowerPredecessor(nj, ni)
							if nj == ni:
								continue
							if nj == nk: # I want j<k
								break
							pair_ni_njp1 = self.generatePairCSCS(ni, controllerStates[index_j + 1])
							pair_ni_nj   = self.generatePairCSCS(ni, nj)
							pair_ni_nk   = self.generatePairCSCS(ni, nk)
							# self.addClause(['-' + lower_nj_ni, '-' + lower_nk_ni, pair_ni_njp1])
							# self.addClause(['-' + lower_nj_ni, '-' + lower_nk_ni, pair_ni_njp1])
							# if nk != controllerStates[index_j + 1]:
							#	self.addClause(['-' + lower_nj_ni, '-' + lower_nk_ni, '-' + pair_ni_nk])


			# real clause
			for ni in controllerStates[1:-1]:
				for nk in controllerStates[1:-1]:
					pair_ni_nk = self.generatePairCSCS(ni, nk)
					lower_nk_ni = self.generateLowerPredecessor(nk, ni)
					if nk == ni:
						continue
					for index_j, nj in enumerate(controllerStates[:-1]):
						if nj == 'n0':
							continue
						pair_ni_nj = self.generatePairCSCS(ni, nj)
						lower_nj_ni = self.generateLowerPredecessor(nj, ni)
						if nj == ni:
							continue
						if nj == nk: # I want j<k
							break
						for a_0 in task.get_action_names():
							if '_DETDUP_0' in a_0:
								other_actions = task.get_other_actions_name(a_0)
								if len(other_actions) == 1:
									a_1 = other_actions[0]
									pair_ni_a0 = self.generatePairActionControllerState(a_0, ni)
									triplet_ni_0_nj = self.generateTripletCSACS(ni, a_0, nj)
									triplet_ni_1_nk = self.generateTripletCSACS(ni, a_1, nk)
									triplet_ni_1_njp1 = self.generateTripletCSACS(ni, a_1, controllerStates[index_j + 1])
									clause = ['-' + pair_ni_nj, '-' + pair_ni_nk, '-' + lower_nj_ni, '-' + lower_nk_ni, '-' + pair_ni_a0]
									self.addClause(clause + [triplet_ni_0_nj])
									#self.addClause(clause + [triplet_ni_1_nk])
									self.addClause(clause + [triplet_ni_1_njp1])




		c2, v2 = self.get_num_cl_vars()
		if debug:
			print 'Generation: Sym Brk\t\t v : ', v2 - v1, '\t\t c : ', c2 - c1, '\t', timer() - start







	def generateSymmetryBreaking(self, task, controllerStates, initialCState, goalCState, debug = False):
		c1, v1 = self.get_num_cl_vars()
		start = timer()
		
		if len(controllerStates) >= 4:		
			self.generateSymmetryBreakingPathGoal(task, controllerStates, initialCState, goalCState, False)
			self.generateSymmetryBreakingFirstAfter(task, controllerStates, initialCState, goalCState, False)
			self.generateSymmetryBreakingOtherNodes(task, controllerStates, initialCState, goalCState, False)
			# 1 is the original, 2 is with RI
			#self.generateSymmetryBreakingOtherNodesOrdering(task, controllerStates, initialCState, goalCState, False)

		c2, v2 = self.get_num_cl_vars()
		if debug:
			print 'Generation: Sym Brk\t\t v : ', v2 - v1, '\t\t c : ', c2 - c1, '\t', timer() - start		

	def generateSymmetryBreakingPathGoal(self, task, controllerStates, initialCState, goalCState, debug = False):
		# \AND_{j <= i} -(nj, ng) --> (ni, n_{i+1})
		c1, v1 = self.get_num_cl_vars()
		start = timer()
		for i, n in enumerate(controllerStates[:-1]):
			# n is the state 'i'
			if n == goalCState:
				continue
			disj = []
			for i in xrange(i + 1):
				disj.append(self.generatePairCSCS(controllerStates[i], goalCState))
			disj.append(self.generatePairCSCS(n, controllerStates[i + 1]))
			self.addClause(disj)

		c2, v2 = self.get_num_cl_vars()
		if debug:
			print 'Generation: Sym Brk\t\t v : ', v2 - v1, '\t\t c : ', c2 - c1, '\t', timer() - start

	def generateSymmetryBreakingFirstAfter(self, task, controllerStates, initialCState, goalCState, debug = False):
		# FirstG(ni) <--> (ni, ng) \land \AND_{j<i} -(nj, ng)
		# AfterG(ni) <--> \OR_{j<i} FirstG(n_j)
		c1, v1 = self.get_num_cl_vars()
		start = timer()
		# First right arrow
		for i, n in enumerate(controllerStates):
			if n == goalCState:
				continue
			var_first = self.generateFirstG(n)
			self.addClause(['-' + var_first, self.generatePairCSCS(n, goalCState)])
			for j in xrange(i):
				nj = controllerStates[j]
				self.addClause(['-' + var_first, '-' + self.generatePairCSCS(nj, goalCState)])
		# First left arrow
		for i, n in enumerate(controllerStates):
			if n == goalCState:
				continue
			var_n_ng = self.generatePairCSCS(n, goalCState)
			clause = ['-' + var_n_ng]
			for j in xrange(i):
				nj = controllerStates[j]
				clause.append(self.generatePairCSCS(nj, goalCState))
			clause.append(self.generateFirstG(n))
			self.addClause(clause)
		# After right arrow
		for i, n in enumerate(controllerStates):
			if n == initialCState:
				continue
			var_after = self.generateAfterG(n)
			clause = ['-' + var_after]
			for j in xrange(i):
				nj = controllerStates[j]
				clause.append(self.generateFirstG(nj))
			self.addClause(clause)
		# After left arrow
		for i, n in enumerate(controllerStates):
			if n == initialCState:
				continue
			for j in xrange(i):
				nj = controllerStates[j]
				self.addClause(['-' + self.generateFirstG(nj), self.generateAfterG(n)])

		c2, v2 = self.get_num_cl_vars()
		if debug:
			print 'Generation: Sym Brk\t\t v : ', v2 - v1, '\t\t c : ', c2 - c1, '\t', timer() - start

	def generateSymmetryBreakingOtherNodes2(self, task, controllerStates, initialCState, goalCState, debug = False):
		# j > i: reachI(nj,k) --> reachI(ni,k) // CNF: -reachI(nj,k) OR reachI(ni,k)
		c1, v1 = self.get_num_cl_vars()
		start = timer()
		value = len(controllerStates) - 1
		for j, nj in enumerate(controllerStates[:-1]):
			for i, ni in enumerate(controllerStates[:-1]):
				# j has to be greater than i
				if j <= i:
					continue
				for k in xrange(value):
					# we have two controller states (ni, nj), such that j > i
					var1 = self.generateReachableI2(ni, k)
					var3 = self.generateAfterG(ni)
					var2 = self.generateReachableI2(nj, k)
					self.addClause(['-' + var2, '-' + var3, var1])

		c2, v2 = self.get_num_cl_vars()
		if debug:
			print 'Generation: Sym Brk\t\t\t v : ', v2 - v1, '\t\t\t c : ', c2 - c1, '\t', timer() - start

	def generateSymmetryBreakingOtherNodes(self, task, controllerStates, initialCState, goalCState, debug = False):
		# (na, nb) \land AfterG(n_{b-1}) --> \OR_{i<=a} (n_i, n_{b-1})
		c1, v1 = self.get_num_cl_vars()
		start = timer()
		if len(controllerStates) <= 4:
			return
		for ia, na in enumerate(controllerStates[:-1]):
			for ib, nb in enumerate(controllerStates[:-1]):
				if ib in [0, 1]:
					continue
				nb_1 = controllerStates[ib - 1]
				var_afterG_nb_1 = self.generateAfterG(nb_1)
				var_pair_ab = self.generatePairCSCS(na, nb)
				disj = ['-' + var_pair_ab, '-' + var_afterG_nb_1]
				for i in xrange(ia + 1):
					ni = controllerStates[i]
					disj.append(self.generatePairCSCS(ni, nb_1))
				self.addClause(disj)

				#for j in range(ia + 1, len(controllerStates) - 1):
				#	nj = controllerStates[j]
				#	pair_post = self.generatePairCSCS(nj, nb_1)
				#	self.addClause(['-' + var_pair_ab, '-' + var_afterG_nb_1, '-' + pair_post])

		c2, v2 = self.get_num_cl_vars()
		if debug:
			print 'Generation: Sym Brk\t\t v : ', v2 - v1, '\t\t c : ', c2 - c1, '\t', timer() - start


	def generateSymmetryBreakingOtherNodesOrdering(self, task, controllerStates, initialCState, goalCState, debug = False):
		# -lower(nk, ni) <--> \OR_{j<i} (nj, nk) // ni != n0
		# lower(nk,n0) <--> (n0, nk)
		# After(ni) \land (ni,nj) \land (ni, nk) \land lower(nj,ni) \land lower(nk,ni) \land (ni, a_name)
		# --> (ni, a_name_0, nj) \land (ni, a_name_1, nk) // nj != ni, nj != ng, nk != ni, nk != ng, j < k
		
		# equivalences for lower
		for ni in controllerStates:
			for nj in controllerStates[:-1]:
				lower_ni_nj = self.generateLowerPredecessor(ni, nj)
				pair_nj_ni  = self.generatePairCSCS(nj, ni)
				clause = [lower_ni_nj, '-' + pair_nj_ni]
				self.addClause(['-' + lower_ni_nj, pair_nj_ni]) # right arrow
				for nk in controllerStates:
					if nk == nj: # want k<j
						break
					pair_nk_ni = self.generatePairCSCS(nk, ni)
					clause.append(pair_nk_ni)
					self.addClause(['-' + lower_ni_nj, '-' + pair_nk_ni]) # right arrow
				self.addClause(clause) # left arrow

		# real clause
		for ni in controllerStates[:-1]:
			afterg_ni = self.generateAfterG(ni)
			for nk in controllerStates[:-1]:
				pair_ni_nk = self.generatePairCSCS(ni, nk)
				lower_nk_ni = self.generateLowerPredecessor(nk, ni)
				if nk == ni:
					continue
				for index_j, nj in enumerate(controllerStates[:-1]):
					pair_ni_nj = self.generatePairCSCS(ni, nj)
					lower_nj_ni = self.generateLowerPredecessor(nj, ni)
					if nj == ni:
						continue
					if nj == nk: # I want j<k
						break
					for a_0 in task.get_action_names():
						if '_DETDUP_0' in a_0:
							other_actions = task.get_other_actions_name(a_0)
							if len(other_actions) == 1:
								a_1 = other_actions[0]
								pair_ni_a0 = self.generatePairActionControllerState(a_0, ni)
								triplet_ni_0_nj = self.generateTripletCSACS(ni, a_0, nj)
								triplet_ni_1_nk = self.generateTripletCSACS(ni, a_1, nk)
								triplet_ni_1_njp1 = self.generateTripletCSACS(ni, a_1, controllerStates[index_j + 1])
								clause = ['-' + afterg_ni, '-' + pair_ni_nj, '-' + pair_ni_nk, '-' + lower_nj_ni, '-' + lower_nk_ni, '-' + pair_ni_a0]
								self.addClause(clause + [triplet_ni_0_nj])
								#self.addClause(clause + [triplet_ni_1_nk])
								self.addClause(clause + [triplet_ni_1_njp1])

	###########################################
	############## CARD SIMPL #################
	###########################################

	def generateCardinalitySimplification(self, task, controllerStates, debug = False):
		c1, v1 = self.get_num_cl_vars()
		start = timer()
		actions = task.get_actions()
		# If cardinality 1
		for a in task.get_action_names():
			card = task.get_action_cardinality(a)
			if card == 1:
				for n in controllerStates:
					for n2 in controllerStates:
						for n3 in controllerStates:
							if n2 == n3:
								continue
							var1 = self.generateTripletCSACS(n, a, n2)
							var2 = self.generatePairCSCS(n, n3)
							self.addClause(['-' + var1, '-' + var2])
			if card == 2:
				a2 = task.get_other_actions_name(a)[0]
				for n in controllerStates:
					for n1 in controllerStates:
						for n2 in controllerStates:
							for n3 in controllerStates:
								if n3 == n1 or n3 == n2:
									continue
								var1 = self.generateTripletCSACS(n, a, n1)
								var2 = self.generateTripletCSACS(n, a2, n2)
								var3 = self.generatePairCSCS(n, n3)
								self.addClause(['-' + var1, '-' + var2, '-' + var3])
			# Way too big, number of clauses grows to the power fo 4 in k
		# HI
		for n in controllerStates:
			for ni in controllerStates:
				for nj in controllerStates:
					if ni == nj:
						continue
					var1 = self.generatePairCSCS(n, ni)
					var2 = self.generatePairCSCS(n, nj)
					clause = ['-' + var1, '-' + var2]
					for a in task.get_action_names():
						card = task.get_action_cardinality(a)
						if card >= 2:
							clause.append(self.generatePairActionControllerState(a, n))
					self.addClause(clause)

		c2, v2 = self.get_num_cl_vars()
		if debug:
			print 'Generation: Card simpl\t\t v : ', v2 - v1, '\t\t c : ', c2 - c1, '\t', timer() - start

	###########################################
	############## MUTEX GROUPS ###############
	###########################################

	def generateMutexGroupsClauses(self, task, controllerStates, debug = False):
		c1, v1 = self.get_num_cl_vars()
		start = timer()
		mutex_groups = task.get_mutex_groups()
		for mg in mutex_groups:
			pairs = self.__get_all_pairs(mg)
			for (a1, a2) in pairs:
				for n in controllerStates:
					var1 = self.generateAtomControllerState(a1, n)
					var2 = self.generateAtomControllerState(a2, n)
					self.addClause(['-' + var1, '-' + var2])

		c2, v2 = self.get_num_cl_vars()
		if debug:
			print 'Generation: Mutex\t\t v : ', v2 - v1, '\t\t c : ', c2 - c1, '\t', timer() - start

	def __get_all_pairs(self, els):
		return [(e1, e2) for (i1, e1) in enumerate(els) for (i2, e2) in enumerate(els) if i2 > i1]

	###########################################
	############## HEUR CONSTR ################
	###########################################

	def generateH1ContraintsClauses(self, task, controllerStates, debug = False):
		c1, v1 = self.get_num_cl_vars()
                start = timer()
		#self._generateH1Constraints1(task, controllerStates)
		self._generateH1Constraints2(task, controllerStates)

		c2, v2 = self.get_num_cl_vars()
		if debug:
			print 'Generation: H1 Constr\t\t v : ', v2 - v1, '\t\t c : ', c2 - c1, '\t', timer() - start

	def _generateH1Constraints1(self, task, controllerStates):
		# Basically, ni can be reachI in at most i steps, no more; every atom p st h1(p) > i cannot be true in ni
		atoms = task.get_atoms()
		for a in atoms:
			h1_value = task.get_h1_value(a)
			for i, n in enumerate(controllerStates):
				if i < h1_value:
					var = self.generateAtomControllerState(a, n)
					self.addClause(['-' + var])

	def _generateH1Constraints2(self, task, controllerStates):
		# Basically, if p(n) is true in n, then -reachI(n, j) for all j < h1(p)
		atoms = task.get_atoms()
		for a in atoms:
			h1_value = task.get_h1_value(a)
			for j in xrange(len(controllerStates)):
				if h1_value > j:
					for n in controllerStates:
						var1 = self.generateReachableI2(n, j)
						var2 = self.generateAtomControllerState(a, n)
						self.addClause(['-' + var1, '-' + var2])