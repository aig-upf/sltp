from CNF import CNF
import os
import sys
from parser import Parser
from myTask import MyTask
from timeit import default_timer as timer
import time

def generateControllerStates(i):
	controllerStates = ['n0']
	for j in xrange(i):
		controllerStates.append('n' + str(j + 1))
	controllerStates.append('ng')
	return controllerStates

start = timer()

p = Parser()
if len(sys.argv) == 4:
	p.set_domain(sys.argv[1])
	p.set_problem(sys.argv[2])
else:
	print 'Incorrect usage\npython main.py <domain> <instance> <auxiliary name for temp files>'
	sys.exit()

name_formula_file = 'formula-' + sys.argv[3] + '.txt'
name_formula_file_extra = 'formula-extra-' + sys.argv[3] + '.txt'
name_output_satsolver = 'outsat-' + sys.argv[3] + '.txt'
name_sas_file = 'outputtrans-' + sys.argv[3] + '.sas'
verbose = True
print("X.0")
p.generate_file(name_sas_file)
print("X.1")
p.generate_task(name_sas_file)
print("X.2")
my_task = p.translate_to_atomic()
p.print_task()
fair = my_task.is_fair()
strong = False
print 'Looking for strong plans: ', strong, '\nFair actions: ', fair
# exit()

cnf = CNF(name_formula_file, name_formula_file_extra, fair, strong)

solver_time = []
for i in xrange(100):
	controllerStates = generateControllerStates(i)

	print '================================================='
	print 'Trying with', len(controllerStates), 'states'
	print 'Number of atoms: ', len(my_task.get_atoms())
	print 'Number of actions: ', len(my_task.get_actions())

	cnf.reset()
	start_g = timer()
	cnf.generate_clauses(my_task, 'n0', 'ng', controllerStates, len(controllerStates), p, verbose)
	
	print 'Generation time = ', timer() - start_g
	print 'Done generation...'
	print '# Clauses =', cnf.getNumberClauses()
	print '# Variables =', cnf.getNumberVariables()

	print 'Creating formula...'
	name_final = cnf.generateInputSat(name_formula_file)
	# name_final = cnf.generateInputSatWSB(name_formula_file)
	
	print 'Done creating formula. Calling solver...'
	start_solv = timer()

	command = './minisat ' + ' ' + name_final + ' ' + name_output_satsolver
	if not verbose:
		command = command + ' | grep "noprint"'
	os.system(command)
	end_solv = timer()
	solver_time.append(end_solv - start_solv)
	print 'Done solver. Round time: ', end_solv - start_solv
	print 'Cumulated solver time: ', sum(solver_time)

	if cnf.parseOutput2(name_output_satsolver, controllerStates, p):
	# if cnf.parseOutputLin(name_output_satsolver, controllerStates, p):
		print '==========='
		print 'SAS!!!'
		break
	print 'UNSAT'

end = timer()
print 'Elapsed total time (s):', end - start
print 'Elapsed solver time (s):', sum(solver_time)
print 'Elapsed solver time (s):', solver_time
print 'Looking for strong plans: ', strong, '\nFair actions: ', fair
os.system('rm' + ' ' + name_formula_file + ' ' + name_output_satsolver + ' ' + name_sas_file + ' ' + name_formula_file_extra + ' ' + name_final)
