import sys

class MyTaskError(Exception):
  def __init__(self, value):
    self.value = value
  def __str__(self):
    return repr(self.value)

class MyTask():
	def __init__(self):
		self.atoms = None # Set
		self.initial = None # Set of atoms
		self.goal = None # Set of atoms
		self.actions = None # Dictionary, name --> [precs, adds, dels]
		self.action_nondet = {} # name --> list other actions
		self.action_cardinality = {} # name --> number
		self.mutex_groups = None # list of list; each sub-list has atoms belonging to the same mutex group
		self.h1_values = {} # atom p --> h1(p)
		self.compatible_actions = {} # action name --> set of compatible actions
		self.mutex_groups_set = [] # list of sets of mutex groups
		self.max_card = 0
		self.relevant_actions_atom = {} # dictionary: atom --> relevant actions
		self.relevant_actions_atom_aname = {} # dictionary: (atom, a_name) --> relevant actions
		self.max_num_args = None # int
		self.action_names = None # set
		self.other_actions_name = {} # Dictionary name --> list of other actions
		self.arguments = {} # i --> list; each list contains the possible arguments in position i (1,2,...)
		self.arguments_action_name = {} # action_name --> list of lists; each sublist contains the possible i-th arguments for the action
		self.action_name_to_actions = {} # action_name --> list of actions
		self.act_name_compatible_arg_i = {} # tuple (i, arg) --> list of act compatible with i-th arg
		self.possible_arguments_for_action = {} # act_name --> list of lists, each sublist is a combination of possible arguments
		self.invalid_next_arguments = {} # Dictionary; (pos, arg) --> not valid args in position pos + 1
		self.valid_next_arguments = {} # Dictionary; (pos, arg) --> possible arguments in position pos + 1
		self.valid_arguments = {} # dictionary tuple (act, arg1, ...) or (act), --> next arguments that are possible
		self.actions_del_p = {} # atom --> actions that delete atom
		self.actions_add_p = {} # atom --> actions that add atom
		self.atoms_aname_never_adds = {} # a_name --> set of atoms that are not added by any action with name a_name
		self.ordered_actions = [] # list of ordered actions

	def getOrderedActions(self):
		return self.ordered_actions

	def is_fair(self):
		for a in self.actions:
			print a
			if '_unfair_' in a:
				return False
		return True

	def get_initial(self):
		return self.initial

	def get_goal(self):
		return self.goal

	def get_atoms(self):
		return self.atoms

	def get_actions(self):
		return self.actions

	def get_action_names(self):
		return self.action_names

	def get_preconditions(self, action):
		return self.actions[action][0]

	def get_add_list(self, action):
		return self.actions[action][1]

	def get_del_list(self, action):
		return self.actions[action][2]

	def get_mutex_groups(self):
		return self.mutex_groups

	def get_h1_value(self, atom):
		return self.h1_values[atom]

	def get_action_cardinality(self, action):
		return self.action_cardinality[action]

	def get_all_possible_arguments(self):
		arguments = []
		for i in xrange(self.max_num_args):
			arguments.append(self.get_all_possible_arguments_i(i))
		return arguments

	def get_possible_arguments_for_action(self, name):
		return self.possible_arguments_for_action[name]

	def get_act_name_compatible_arg_i(self, i, arg):
		return self.act_name_compatible_arg_i[(i, arg)]

	def get_arguments_action_name(self, name):
		return self.arguments_action_name[name]

	def get_actions_with_name(self, name):
		return self.action_name_to_actions[name]

	def get_other_actions_name(self, action):
		return self.other_actions_name[action]

	def get_relevant_actions_2(self, atom, a_name):
		return self.relevant_actions_atom_aname[(atom, a_name)]

	def get_all_possible_arguments_i(self, i):
		if i >= 0 and i < self.max_num_args:
			return self.arguments[i]
		else:
			raise MyTaskError('Bad argument index...')

	def get_action_name(self, action):
		return action.split('(')[0]

	def get_action_args(self, action):
		num_args = self.get_number_args(action)
		max_args = self.max_num_args
		if num_args == 0:
			return ['---' for i in xrange(max_args)]
		if num_args == 1:
			args = [action.split('(')[1].split(')')[0]]
		else:
			args = action.split('(')[1].split(')')[0].split(',')
		return args + ['---' for i in xrange(max_args - num_args)]

	def get_number_args(self, action):
		if '()' in action:
			return 0
		return action.count(',') + 1

	def get_valid_next_arguments(self, i, arg):
		return self.valid_next_arguments[(i, arg)]

	def get_invalid_next_arguments(self, i, arg):
		return self.invalid_next_arguments[(i, arg)]

	def set_atoms(self, atoms, debug = False):
		"""Set of all atoms"""
		self.atoms = atoms
		print '# Atoms:', len(atoms)
		if debug:
			print 'set atoms:', self.atoms

	def set_initial(self, initial, debug = False):
		""""Set of initial"""
		self.initial = initial
		if debug:
			print 'set initial', self.initial

	def set_goal(self, goal, debug = False):
		"""Goal"""
		self.goal = goal
		if debug:
			print 'goal', self.goal

	def set_actions_atomic(self, actions, debug = False):
		"""Actions"""
		self.actions = actions
		print '# Actions:', len(actions)
		print '\tSetting other actions'
		self.generate_other_actions()
		print '\tSetting action card'
		self.generate_action_cardinality()
		if debug:
			for a in self.actions:
				print a, self.action_nondet[a]
				print actions[a]

	def set_mutex_groups(self, m_groups, debug = False):
		"""Mutex groups"""
		self.mutex_groups = m_groups
		if debug:
			for mg in self.mutex_groups:
				print mg

	def generate_action_cardinality(self):
		# self.action_nondet = {} # name --> list other actions
		# self.action_cardinality = {} # name --> number
		for a in self.actions:
			others = self.action_nondet[a]
			self.action_cardinality[a] = len(others) + 1
			name = self.get_action_name(a)
			self.action_cardinality[name] = len(others) + 1

	def generate_other_actions(self):
		counter = 0
		total = len(self.actions)
		for a in self.actions:
			if counter % 1000 == 0:
				print counter, '/', total
			counter += 1
			if '_DETDUP_' not in a:
				# No nondet effects
				self.action_nondet[a] = []
				continue
			other = []
			for a2 in self.actions:
				if a2 == a:
					continue
				if self.__check_if_belong_same_det(a, a2):
					other.append(a2)
			self.action_nondet[a] = other

	def __check_if_belong_same_det(self, a1, a2):
		# I know that a1 has '_DETDUP_' in name
		sep = '_DETDUP_'
		if sep not in a2:
			return False
		if a1.split(sep)[0] != a2.split(sep)[0]:
			return False
		if a1.split('(')[1] != a2.split('(')[1]:
			return False
		return True

	def get_other_actions(self, action):
		return self.action_nondet[action]

	def get_other_actions_ordered(self, a):
		group_actions = [a]
		for action in self.action_nondet[a]:
			group_actions.append(action)
		if len(group_actions) == 1:
			return group_actions
		else:
			ordered_actions = []
			for i in xrange(len(group_actions)):
				minimum = group_actions[0]
				for a in group_actions:
					minimum = self.__get_smaller_action(minimum, a)
				ordered_actions.append(minimum)
				group_actions.remove(minimum)
			return ordered_actions

	def __get_smaller_action(self, a1, a2):
		sep = '_DETDUP_'
		if sep not in a1 or sep not in a2:
			raise MyTaskError('Actions deterministic, no order!')
		i1 = self.__get_action_index(a1)
		i2 = self.__get_action_index(a2)
		if i1 <= i2:
			return a1
		return a2

	def __get_action_index(self, a):
		sep = '_DETDUP_'
		return int(a.split(sep)[1].split('(')[0])

	def set_h_values(self, debug):
		atoms_layers = self._compute_atoms_layers()
		# H1
		for a in self.atoms:
			self.h1_values[a] = self._get_minimum_layer1(a, atoms_layers)
		# H2
		# for a1 in self.atoms:
		# 	for a2 in self.atoms:
		# 		if a1 != a2:
		# 			pair = (a1, a2)
		# 			self.h2_values[pair] = self._get_minimum_layer2(a1, a2, atoms_layers)
		if debug:
			# for a1 in self.atoms:
			# 	for a2 in self.atoms:
			# 		pair = (a1, a2)
			# 		if a1 != a2:
			# 			print 'h2', pair, '=', self.h2_values[pair]
			for a in self.atoms:
				print 'h1(', a, ') =', self.h1_values[a]

	def set_maximum_cardinality(self, debug):
		max_card = 0
		for a in self.actions:
			if self.action_cardinality[a] > max_card:
				max_card = self.action_cardinality[a]
		self.max_card = max_card

	def _get_minimum_layer1(self, atom, layers):
		for i in xrange(len(layers)):
			layer = layers[i]
			if atom in layer:
				return i
		# if it is not reachable
		return -1

	def _get_minimum_layer2(self, a1, a2, layers):
		for i in xrange(len(layers)):
			layer = layers[i]
			if a1 in layer and a2 in layer:
				return i
		return -1

	def _compute_atoms_layers(self):
		added_atoms = set([])
		layer = set([])
		for a in self.initial:
			added_atoms.add(a)
			layer.add(a)
		layers = {}
		layers[0] = layer
		for i in xrange(1, len(self.atoms) + 1):
			prev_layer = layers[i - 1]
			curr_layer = set([])
			for a in prev_layer:
				curr_layer.add(a)
			for a in self.actions:
				if self._is_applicable(a, prev_layer):
					for eff in self.get_add_list(a):
						curr_layer.add(eff)
			layers[i] = curr_layer
		return layers

	def _is_applicable(self, action, atoms):
		for pre in self.get_preconditions(action):
			if pre not in atoms:
				return False
		return True

	def set_relevant_actions_split(self, debug = False):
		for p in self.atoms:
			self.actions_del_p[p] = self.compute_actions_del_atom(p)
		for p in self.atoms:
			self.actions_add_p[p] = self.compute_actions_add_atom(p)
		for a_name in self.action_names:
			self.atoms_aname_never_adds[a_name] = self.compute_atoms_aname_never_adds(a_name)

	def compute_actions_add_atom(self, atom):
		actions = set()
		for a in self.actions:
			add_list = self.get_add_list(a)
			if atom in add_list:
				actions.add(a)
		return actions

	def compute_actions_del_atom(self, atom):
		actions = set()
		for a in self.actions:
			del_list = self.get_del_list(a)
			if atom in del_list:
				actions.add(a)
		return actions

	def compute_atoms_aname_never_adds(self, a_name):
		atoms = set(self.atoms)
		for a in self.get_actions_with_name(a_name):
			add_list = set(self.get_add_list(a))
			atoms = atoms - add_list
		return atoms

	def set_relevant_actions(self, debug = False):
		# action a is relevant to atom p if:
		# - p \in del(b)
		# - p \in add(b)
		# - p \not \in add(b) and some sibling action adds p
		for p in self.atoms:
			rel_actions = self._compute_relevant_actions(p)
			self.relevant_actions_atom[p] = rel_actions

	def _compute_relevant_actions(self, atom):
		rel_actions = set([])
		for a in self.actions:
			add_list = self.get_add_list(a)
			del_list = self.get_del_list(a)
			if atom in del_list or atom in add_list:
				rel_actions.add(a)
			siblings = self.get_other_actions(a)
			for s in siblings:
				add_sibling = self.get_add_list(s)
				if atom in add_sibling and atom not in add_list:
					rel_actions.add(a)
		return rel_actions

	def get_relevant_actions(self, atom):
		return self.relevant_actions_atom[atom]

	def create_compatible_actions(self, debug = False):
		self.mutex_groups_set = [set(mg) for mg in self.mutex_groups]
		counter = 0
		total = len(self.actions)
		for a1 in self.actions:
			if counter % 1000 == 0:
				print counter, '/', total
			counter += 1
			compatible_acts = set([])
			other_actions = self.action_nondet[a1]
			#for a2 in self.actions:
			for a2 in self.get_actions_with_name(self.get_action_name(a1)):
				if a2 == a1 or a2 in other_actions:
					continue
				if self._actions_are_compatible(a1, a2):
					compatible_acts.add(a2)
			self.compatible_actions[a1] = compatible_acts

	def _actions_are_compatible(self, a1, a2):
		precs1 = self.get_preconditions(a1)
		precs2 = self.get_preconditions(a2)
		pairs = self._get_all_pairs(precs1, precs2)
		for p1, p2 in pairs:
			if p1 == p2:
				continue
			if self._atoms_belong_to_same_mutex(p1, p2):
				return False
		return True

	def _atoms_belong_to_same_mutex(self, p1, p2):
		for mg in self.mutex_groups_set:
			if p1 in mg and p2 in mg:
				return True
		return False

	def _get_all_pairs(self, coll1, coll2):
		pairs = []
		for i in coll1:
			for j in coll2:
				pairs.append([i, j])
		return pairs

	def actions_are_compatible(self, a1, a2):
		return a2 in self.compatible_actions[a1]

	def initialize_splitting(self, debug):
		# Everything to do with info regarding the splitting of actions
		self.max_num_args = self.compute_max_num_args()
		if debug:
			print 'PRINTING MAX AND DECOMPOSITION'
			print self.max_num_args
			print '================================'
			for a in self.actions:
				print a, '-------', self.decompose_action(a)
			print '================================'
		self.action_names = self.compute_action_names()
		if debug:
			print 'PRINTING ACTION NAMES'
			print self.action_names
			print '================================'
		self.other_actions_name = self.compute_other_actions_name()
		if debug:
			print 'PRINTING OTHER ACTIONS'
			for i in self.other_actions_name:
				print i, '---', self.get_other_actions_name(i)
			print '=============================='
		self.arguments = self.compute_all_possible_arguments()
		if debug:
			print 'PRINTING POSSIBLE ARGS IN POSITION i'
			for i in xrange(self.max_num_args):
				print i, '---', self.get_all_possible_arguments_i(i)
				print ' '
			print '=============================='
		self.action_name_to_actions = self.compute_actions_name_to_actions()
		if debug:
			print 'PRINTING ACTION NAME TO POSSIBLE ACTIONS'
			counter = 0
			for n in self.action_name_to_actions:
				print n, '---', self.get_actions_with_name(n)
				counter += len(self.get_actions_with_name(n))
				print counter, '/', len(self.actions)
				print ' '
			print '=============================='
		self.arguments_action_name = self.compute_possible_arguments_action_name()
		if debug:
			print 'PRINTING VALID ARGUMENTS FOR ACTION'
			for a_name in self.action_names:
				print a_name
				for valid_args in self.get_arguments_action_name(a_name):
					print valid_args
			print '=============================='
		self.act_name_compatible_arg_i = self.compute_compatible_act_names_with_arg_i()
		if debug:
			print 'PRINTING ACTIONS COMPAT WITH ARGUMENT IN POSITION i'
			for (i, arg) in self.act_name_compatible_arg_i:
				print (i, arg), '---', self.get_act_name_compatible_arg_i(i, arg)
			print '=============================='
		self.possible_arguments_for_action = self.compute_possible_arguments_actions()
		if debug:
			print 'PRINTING POSSIBLE SETS OF ARGUMENTS FOR ACTIONS'
			for a_name in self.possible_arguments_for_action:
				print a_name, '---', self.get_possible_arguments_for_action(a_name)
			print '=============================='
		self.relevant_actions_atom_aname = self.compute_relevant_actions_atom_aname()
		if debug:
			print 'PRINTING RELEVANT ACTIONS ANAME ATOM'
			for a in self.action_names:
				for p in self.atoms:
					print (p, a), '---', self.get_relevant_actions_2(p, a)
			print '=============================='
		self.valid_next_arguments = self.compute_valid_next_arguments() # self.valid__nextarguments = {} # Dictionary; (pos, arg) --> possible arguments in position pos + 1
		if debug:
			print 'PRINTING VALID ARGUMENTS'
			print self.max_num_args
			for i, arg in self.valid_next_arguments:
				print i, arg, '---', self.get_valid_next_arguments(i, arg)
			print '=============================='
		self.invalid_next_arguments = self.compute_invalid_next_arguments() # self.invalid_next_arguments = {} # Dictionary; (pos, arg) --> not valid args in position pos + 1
		if debug:
			print 'PRINTING INVALID ARGUMENTS'
			print self.max_num_args
			for i, arg in self.invalid_next_arguments:
				print i, arg, '---', self.get_invalid_next_arguments(i, arg)
			print '=============================='
		self.valid_arguments = self.compute_valid_arguments()
		if debug:
			print 'PRINTING VALID ARGUMENTS ALL'
			for key in self.valid_arguments:
				print key, '===', self.valid_arguments[key]
			print '=============================='
		self.ordered_actions = self.compute_ordered_actions() # list of ordered actions ##################################
		#if debug:
		#	print 'PRINTING VALID ARGUMENTS ALL'
		#	for key in self.valid_arguments:
		#		print key, '===', self.valid_arguments[key]
		#	print '==============================' #############################################
		if debug:
			print 'ACTION NAME CARDINALITY'
			for name in self.action_names:
				print name, '---', self.get_action_cardinality(name) ###
			print '=============================='
		# sys.exit() ############ To test splitting generation

	def compute_ordered_actions(self):
		actions = set()
		for a in self.actions:
			other = self.get_other_actions(a)
			already_present = False
			for o in other:
				if o in actions:
					already_present = True
					break
			if already_present:
				continue
			actions.add(a)
		return list(actions)

	def compute_valid_arguments(self):
		valid_arguments = {}
		for a in self.actions:
			split_action = self.decompose_action(a)
			for i in xrange(1, len(split_action)):
				# print split_action, len(split_action)
				key = tuple(split_action[:i])
				val = split_action[i]
				if key not in valid_arguments:
					valid_arguments[key] = set()
				valid_arguments[key].add(val)
		return valid_arguments

	def compute_invalid_next_arguments(self):
		if len(self.valid_next_arguments) == 0:
			self.compute_valid_next_arguments()
		invalid_args = {}
		for i, arg in self.valid_next_arguments:
			pair = (i, arg)
			all_args = self.get_all_possible_arguments_i(i + 1)
			invalid_args[pair] = set(all_args) - self.valid_next_arguments[pair]
		return invalid_args

	def compute_valid_next_arguments(self):
		valid_arguments = {}
		for a in self.actions:
			# a_name = self.get_action_name(a)
			a_args = self.get_action_args(a)
			for i, arg in enumerate(a_args[:-1]):
				pair = (i, arg)
				next_arg = a_args[i + 1]
				if pair not in valid_arguments:
					valid_arguments[pair] = set()
				valid_arguments[pair].add(next_arg)
		return valid_arguments

	def compute_relevant_actions_atom_aname(self):
		atom_aname_to_rel_actions = {}
		for a_name in self.action_names:
			actions = self.get_actions_with_name(a_name)
			for p in self.atoms:
				atom_aname_to_rel_actions[(p, a_name)] = []
				for a in actions:
					add_list = self.get_add_list(a)
					del_list = self.get_del_list(a)
					if p in add_list or p in del_list:
						atom_aname_to_rel_actions[(p, a_name)].append(a)
		return atom_aname_to_rel_actions

	def compute_possible_arguments_actions(self):
		act_to_possible_args = {}
		for a_name in self.action_names:
			possible_arguments = []
			for a in self.action_name_to_actions[a_name]:
				possible_arguments.append(self.get_action_args(a))
			act_to_possible_args[a_name] = possible_arguments
		return act_to_possible_args

	def compute_compatible_act_names_with_arg_i(self):
		arg_position_action = {}
		for i in xrange(self.max_num_args):
			args_i = self.get_all_possible_arguments_i(i)
			for arg in args_i:
				arg_position_action[(i, arg)] = self.compute_compatible_act_names_with_specific_arg(i, arg)
		return arg_position_action

	def compute_compatible_act_names_with_specific_arg(self, i, arg):
		compat_actions = set([])
		for a in self.actions:
			name = self.get_action_name(a)
			args = self.get_action_args(a)
			if args[i] == arg:
				compat_actions.add(name)
		return list(compat_actions)

	def compute_actions_name_to_actions(self):
		a_names_to_actions = {}
		for a_name in self.action_names:
			a_names_to_actions[a_name] = []
		for a in self.actions:
			name = self.get_action_name(a)
			a_names_to_actions[name].append(a)
		return a_names_to_actions

	def compute_possible_arguments_action_name(self):
		arguments = {}
		for a_name in self.action_names:
			arguments_a = []
			for i in xrange(self.max_num_args):
				args_a_i = set([])
				for a in self.action_name_to_actions[a_name]:
					args_a_i.add(self.get_action_args(a)[i])
				arguments_a.append(list(args_a_i))
			arguments[a_name] = arguments_a
		return arguments

	def compute_all_possible_arguments(self):
		args = {}
		for i in xrange(self.max_num_args):
			index = i
			args_i = set([])
			for a in self.actions:
				args_i.add(self.get_action_args(a)[index])
			args[index] = list(args_i)
		return args

	def compute_other_actions_name(self):
		other = {}
		for a in self.action_names:
			if '_DETDUP_' not in a:
				other[a] = []
				continue
			other_names = []
			for a2 in self.action_names:
				if a2 == a:
					continue
				if self.get_action_name(a).split('_DETDUP_')[0] == self.get_action_name(a2).split('_DETDUP_')[0]:
					other_names.append(a2)
			other[a] = other_names
		return other

	def compute_action_names(self):
		names = set([])
		for a in self.actions:
			name = self.get_action_name(a)
			names.add(name)
		return names

	def decompose_action(self, action):
		name = self.get_action_name(action)
		args = self.get_action_args(action)
		return [name] + args

	def compute_max_num_args(self):
		max_args = 0
		for a in self.actions:
			num_args = self.get_number_args(a)
			if num_args > max_args:
				max_args = num_args
		return max_args

	def print_task(self):
		print 'ATOMS ================================================'
		for a in self.atoms:
			print a
		print 'INITIAL =============================================='
		for a in self.initial:
			print a
		print 'GOAL ================================================='
		for a in self.goal:
			print a
		print 'ACTIONS =============================================='
		for a in self.actions:
			print a, self.get_other_actions(a)
			print a, self.actions[a]