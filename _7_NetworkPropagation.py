# -*- coding: ms949 -*-

import datetime
import os
import _0_Preprocess
import copy
import random
import math
import PPI
import threading
import multiprocessing
import MyUtil_pypy
#import numpy
import MULTIPROCESS
from datetime import datetime
import time
import sys
import cPickle

import statistics

#_0_Preprocess.initialize()


RND_ID = _0_Preprocess.RND_ID
TEMP_PATH = _0_Preprocess.TEMP_PATH


print 'RND_ID = ', RND_ID

R_PARAMETER = 0.3
DISEASE_SET = []

Pnp_matrix = None
Pnp_dict = {}


USE_CACHE = True
fcache = None

#fcache = './7. Network Propagation/' + \
#    os.path.basename (_0_Preprocess.STRING_FILE) + \
#    '_' + str( _0_Preprocess.STRING_CUTOFF).strip() + '.cache.txt'




MP_PROCESSORS = 1
REMOVE_TEMP_FILES = True



P_VALUE = False
PRELOAD_NO_OF_RANDOM_NETWORKS = 100  # # of random networks to calculate p-values, default=100



HOW_MANY_TEST_MODIFERS = 100  # iteratio, default=100
MAX_THREAD_FOR_ITERATION = 1
NETWORK_TEMP_FILES = []
DEGREE = 1
# 0 - random
# 1 - random preserving interaction degrees
# 3 - interaction reliability
# =================================




















#class JOB(multiprocessing.Process):
class JOB(threading.Thread):

	values = None
	over = False
	result = None

	def __init__(self):
		self.values = None
		self.over = False
		self.result = None
		threading.Thread.__init__(self)
		#multiprocessing.Process.__init__(self)


	def set(self, values):
		self.values = values

	def run(self):

		self.over = False

		disease = self.values[0]
		mod_pos = self.values[1]
		total_genes = self.values[2]
		total_genes_without_modifiers = self.values[3]
		iteration_th = self.values[4]
		ppi_box = self.values[5]

		self.result = mainJob(disease, mod_pos, total_genes, total_genes_without_modifiers, iteration_th, ppi_box)

		self.over = True

	def getResult(self):
		return self.result




def __getModifiers():


	r = _0_Preprocess.getModifiers()
	return r


def __randomSampling (total_genes, pos_mod, basket, ppi):

	new_pos_mod = []


	global DEGREE  # default=3

	if DEGREE == 0:
		for g in random.sample( total_genes, len(pos_mod) ):
			new_pos_mod.append(g)

	elif DEGREE == 1:
		# preserving interaction degrees

		for g in pos_mod:

			deg = len(ppi.getPartners(g))
			deg_str1 = str(deg-1).strip()
			deg_str2 = str(deg).strip()
			deg_str3 = str(deg+1).strip()


			box = []
			if basket.has_key(deg_str1):
				box = box + basket[deg_str1]
			if basket.has_key(deg_str2):
				box = box + basket[deg_str2]        
			if basket.has_key(deg_str3):
				box = box + basket[deg_str3]        


			while(True):        
				v = random.sample( box, 1 ) [0]
				if not v in new_pos_mod:
					new_pos_mod.append(v)
					break



	elif DEGREE == 3:
		# preserving reliability sum
		new_pos_mod = __randomSamplingBasedOnReliability(total_genes, [], pos_mod, ppi)

		if len(new_pos_mod) != len(pos_mod):
			print 'Error: DEGREE=3'
			print 'pos=',len(pos_mod), ' new pos=', len(new_pos_mod)

	else:

		print '[NetProp] Degree Error'
		sys.exit(1)

	return new_pos_mod

'''

def __getRandomNetwork():

    # random networks

    global NETWORK_TEMP_FILES 

    [interaction_file, index1, index2, score_index, threshold] = NETWORK_TEMP_FILES

    ppi_filename_only = os.path.basename(interaction_file) # file names
    ppi_filename_only = os.path.splitext(ppi_filename_only)[0] + '_' + str(threshold).strip()

    files, folders = MyUtil.getFileListIn('./7. Network Propagation/random networks/')
    prebuilt = []
    for f in files:
        if f.find(ppi_filename_only)>=0:
            prebuilt.append('./7. Network Propagation/random networks/' + f)    

    fname = random.sample( prebuilt, 1 ) [0]



    ppi = PPI.PPIDB_STRING()
    ppi.load(fname, 0, 1, 2, -1)
    return ppi
'''

def __getRandomNetworks(number):
	# random networks

	global NETWORK_TEMP_FILES 



	#ppi_filename_only = os.path.basename(interaction_file) # file names
	#ppi_filename_only = os.path.splitext(ppi_filename_only)[0] + '_' + str(threshold).strip()
	path = './ppi/random/human/'  + str(_0_Preprocess.STRING_CUTOFF).strip() 
	files, folders = MyUtil_pypy.getFileListIn(path)
	prebuilt = []
	for f in files:
		#print f
		if f.find('.txt')>=0 and f.find('.cache')<0:
			prebuilt.append(path+   '/' + f)    

	fnames = random.sample( prebuilt, number ) 




	return fnames    



def loadInteraction(interaction_file, iteration, index1, index2, score_index, threshold, separator):



	print 'Loading interactions = ', interaction_file
	box = []

	ppi = PPI.PPIDB_STRING()
	ppi.load(interaction_file, index1, index2, score_index, threshold, separator, _0_Preprocess.STRING_ALIAS_FILE)
	box.append(ppi) # normal network

	#ppi_filename_only = os.path.basename(interaction_file) 
	#ppi_filename_only = os.path.splitext(ppi_filename_only)[0] + '_' + str(threshold).strip() 
	# add random networks?
	#path = './'


	return box



def run_mp(args):
	config_file, ppi_info, iteration, disease, interaction_file, index1, index2, score_index, threshold, separator, queue = args

	init(config_file)


	r = run2(ppi_info, iteration, disease, interaction_file, index1, index2, score_index, threshold, separator)
	queue.put(r)
	return r


def run(ppi_info, iteration, interaction_file, index1, index2, score_index, threshold, separator):


	global RND_ID, PRELOAD_NO_OF_RANDOM_NETWORKS, NETWORK_TEMP_FILES, DEGREE, DISEASE_SET, R_PARAMETER, Pnp_dict, Pnp_matrix
	global MAX_THREAD_FOR_ITERATION, MP_PROCESSORS
	global fcache


	#_0_Preprocess.init(config_file)

	Pnp_matrix = None
	Pnp_dict = {}


	#int_file = os.path.basename(interaction_file)
	#fcache = './7. Network Propagation/' + int_file + '_' + str(threshold).strip() + '.cache.txt'
	fcache = './7. Network Propagation/' + \
	        os.path.basename ( interaction_file ) + \
	        '_' + str( threshold ).strip() + '.cache.txt'




	mod_maps, mods = __getModifiers()


	summary_file = TEMP_PATH + '/NetProp' + RND_ID + '_' + os.path.basename(interaction_file) + '_summary.txt'
	f=open(summary_file, 'w')
	dd = mods.keys()
	dd.sort()


	NETWORK_TEMP_FILES = [interaction_file, index1, index2, score_index, threshold, separator]  

	info = '--------------------------------------------' + '\n' + \
	        'THREAD = ' + str( MP_PROCESSORS ) + '\n' + \
	        'Remove temporary files = ' + str( REMOVE_TEMP_FILES ) + '\n' + \
	        'Calculate p-values = ' + str( P_VALUE) + '\n' + \
	        '\t# of random networks = ' + str( PRELOAD_NO_OF_RANDOM_NETWORKS) + '\n' + \
	        'Iteration = ' + str( HOW_MANY_TEST_MODIFERS) + '\n' + \
	        'Degree (method for random networking) = ' + str(DEGREE) + '\n' + \
	        'r = ' + str(R_PARAMETER)
	print info





	ppi_box = loadInteraction(interaction_file, PRELOAD_NO_OF_RANDOM_NETWORKS, index1, index2, score_index, threshold, separator) # 100개 네트워크를 만든다.
	total_genes = __getTotalGeneList(ppi_box)





	print 'Initialization....',
	#global Pnp_matrix, Pnp_dict
	total_gene_names = ppi_box[0].getWholeGeneList()
	total_gene_names.sort()

	#Pnp_matrix, Pnp_dict = __buildPnpMatrix(total_genes, ppi_box[0], cache_file=fcache)

	# load cache
	#fcache = './7. Network Propagation/' + \
	#    os.path.basename (_0_Preprocess.STRING_FILE) + \
	#    '_' + str( _0_Preprocess.STRING_CUTOFF).strip() + '.cache.txt'
	#Pnp_matrix, Pnp_dict = __buildPnpMatrix(total_genes, ppi_box[0], cache_file=fcache)

	#Pnp_matrix, Pnp_dict = __buildPnpMatrix(total_gene_names, ppi_box[0])
	print 'done'



	# copy a cache file over a temp folder
	# copyCache()


	#USE_CACHE = False ###########################
	init = True


	for disease in dd:

		print '--------------------------------------'
		print disease


		ppi_box2 = []
		ppi_box2 = copy.deepcopy(ppi_box)


		auc_output =  TEMP_PATH + '/NetPropagation_' + RND_ID + '_' + disease + '_auc.txt'


		mod_pos = mods[disease] 

		new_mod_pos = __removeUnannotatedModifiers(total_genes, mod_pos)
		total_genes_without_modifiers = __getTotalGenesWithoutModifiers(total_genes, mod_pos)


		s = 'Interaction = ' + interaction_file + '\n' + \
		        'Disease = ' + disease + '\n' + \
		        'Total annotated proteins # = ' + str(len(total_genes)) + '\n' + \
		        'Positive modifiers # = ' + str( len(new_mod_pos) ) + '\n'


		f.write(info+'\n')
		f.write('-------------------------\n' +s+'\n')

		print s

		prank = []




		if MAX_THREAD_FOR_ITERATION > 1:


			v = [ ppi_info, disease, new_mod_pos, total_genes, total_genes_without_modifiers, 0 ,ppi_box2]

			box = []

			for i in range(iteration):
				th = MULTIPROCESS.JOB_THREAD()
				th.set_args(mainjob_4_mp, v)
				box.append(th)

			ret = MULTIPROCESS.runMultiprocesses(box, max_cpu=MP_PROCESSORS, sleep_interval_seconds=5)



			for o in ret:
				prank.append(o)

			print ''

		else:

			# no thread
			for i in range(iteration):


				print '\r\t', RND_ID, ' Disease = ', disease, ' # = ', i+1,

				pos_rank = mainJob( ppi_info, disease, new_mod_pos, total_genes, total_genes_without_modifiers, i ,ppi_box2)

				print '\t', pos_rank,

				prank.append(  pos_rank )


			print ''
		'''





        for i in range(iteration):


            print '\r\t', RND_ID, ' Disease = ', disease, ' # = ', i+1,

            pos_rank = mainJob( disease, new_mod_pos, total_genes, total_genes_without_modifiers, i ,ppi_box2)

            print '\t', pos_rank,

            prank.append(  pos_rank )


            print ''
        '''


		foname = TEMP_PATH + '/NetProp' +RND_ID+'_'+disease+'_ranks.txt'
		__saveRank(prank, foname)



		pos_auc = __calculateAUC(prank, auc_output)


		s = 'Avg rank = ' + str( statistics.average( prank  ) )  + ' +- ' + str( statistics.stdev(prank) ) + '\n' + \
		        'AUC + = ' + str( pos_auc ) + '\n' + \
		        'Output file = ' + auc_output + '\n' + \
		        'RND_ID = ' + RND_ID



		print s
		f.write(s+'\n')
	f.close()

	print 'Result summary = ', summary_file

	return summary_file








def run2(ppi_info, iteration, disease, interaction_file, index1, index2, score_index, threshold, separator):


	global RND_ID, PRELOAD_NO_OF_RANDOM_NETWORKS, NETWORK_TEMP_FILES, MP_PROCESSORS , DEGREE, DISEASE_SET, R_PARAMETER, Pnp_dict, Pnp_matrix
	global MAX_THREAD_FOR_ITERATION, MP_PROCESSORS
	global fcache


	#_0_Preprocess.init(config_file)

	Pnp_matrix = None
	Pnp_dict = {}


	#int_file = os.path.basename(interaction_file)
	#fcache = './7. Network Propagation/' + int_file + '_' + str(threshold).strip() + '.cache.txt'
	fcache = './7. Network Propagation/' + \
	        os.path.basename ( interaction_file ) + \
	        '_' + str( threshold ).strip() + '.cache.txt'




	mod_maps, mods = __getModifiers()


	summary_file = TEMP_PATH + '/NetProp' + RND_ID + '_' + os.path.basename(interaction_file) + '_summary.txt'
	f=open(summary_file, 'w')
	dd = mods.keys()
	dd.sort()


	NETWORK_TEMP_FILES = [interaction_file, index1, index2, score_index, threshold, separator]  

	info = '--------------------------------------------' + '\n' + \
	        'THREAD = ' + str( MP_PROCESSORS ) + '\n' + \
	        'Remove temporary files = ' + str( REMOVE_TEMP_FILES ) + '\n' + \
	        'Calculate p-values = ' + str( P_VALUE) + '\n' + \
	        '\t# of random networks = ' + str( PRELOAD_NO_OF_RANDOM_NETWORKS) + '\n' + \
	        'Iteration = ' + str( HOW_MANY_TEST_MODIFERS) + '\n' + \
	        'Degree (method for random networking) = ' + str(DEGREE) + '\n' + \
	        'r = ' + str(R_PARAMETER)
	print info





	ppi_box = loadInteraction(interaction_file, PRELOAD_NO_OF_RANDOM_NETWORKS, index1, index2, score_index, threshold, separator) # 100개 네트워크를 만든다.
	total_genes = __getTotalGeneList(ppi_box)





	print 'Initialization....',
	#global Pnp_matrix, Pnp_dict
	total_gene_names = ppi_box[0].getWholeGeneList()
	total_gene_names.sort()

	#Pnp_matrix, Pnp_dict = __buildPnpMatrix(total_genes, ppi_box[0], cache_file=fcache)

	# load cache
	#fcache = './7. Network Propagation/' + \
	#    os.path.basename (_0_Preprocess.STRING_FILE) + \
	#    '_' + str( _0_Preprocess.STRING_CUTOFF).strip() + '.cache.txt'
	#Pnp_matrix, Pnp_dict = __buildPnpMatrix(total_genes, ppi_box[0], cache_file=fcache)

	#Pnp_matrix, Pnp_dict = __buildPnpMatrix(total_gene_names, ppi_box[0])
	print 'done'



	# copy a cache file over a temp folder
	# copyCache()


	#USE_CACHE = False ###########################
	init = True




	print '--------------------------------------'
	print disease


	ppi_box2 = []
	ppi_box2 = copy.deepcopy(ppi_box)


	auc_output =  TEMP_PATH + '/NetPropagation_' + RND_ID + '_' + disease + '_auc.txt'


	mod_pos = mods[disease] 

	new_mod_pos = __removeUnannotatedModifiers(total_genes, mod_pos)
	total_genes_without_modifiers = __getTotalGenesWithoutModifiers(total_genes, mod_pos)


	s = 'Interaction = ' + interaction_file + '\n' + \
	        'Disease = ' + disease + '\n' + \
	        'Total annotated proteins # = ' + str(len(total_genes)) + '\n' + \
	        'Positive modifiers # = ' + str( len(new_mod_pos) ) + '\n'


	f.write(info+'\n')
	f.write('-------------------------\n' +s+'\n')

	print s

	prank = []




	if MAX_THREAD_FOR_ITERATION > 1:


		v = [ ppi_info, disease, new_mod_pos, total_genes, total_genes_without_modifiers, 0 ,ppi_box2]

		box = []

		for i in range(iteration):
			th = MULTIPROCESS.JOB_THREAD()
			th.set_args(mainjob_4_mp, v)
			box.append(th)

		ret = MULTIPROCESS.runMultiprocesses(box, max_cpu=MAX_THREAD_FOR_ITERATION, sleep_interval_seconds=5)



		for o in ret:
			prank.append(o)

		print ''

	else:

		# no thread
		for i in range(iteration):


			print '\r\t', RND_ID, ' Disease = ', disease, ' # = ', i+1,

			pos_rank = mainJob( ppi_info, disease, new_mod_pos, total_genes, total_genes_without_modifiers, i ,ppi_box2)

			print '\t', pos_rank,

			prank.append(  pos_rank )


		print ''


	foname = TEMP_PATH + '/NetProp' +RND_ID+'_'+disease+'_ranks.txt'
	__saveRank(prank, foname)



	pos_auc = __calculateAUC(prank, auc_output)


	s = 'Avg rank = ' + str( statistics.average( prank  ) )  + ' +- ' + str( statistics.stdev(prank) ) + '\n' + \
	        'AUC + = ' + str( pos_auc ) + '\n' + \
	        'Output file = ' + auc_output + '\n' + \
	        'RND_ID = ' + RND_ID



	print s
	f.write(s+'\n')
	f.close()

	print 'Result summary = ', summary_file

	return summary_file


def mainjob_4_mp(args):


	#global Pnp_dict, Pnp_matrix

	ppi_info, disease, new_mod_pos, total_genes, total_genes_without_modifiers, i ,ppi_box2, queue = args



	rank = mainJob( ppi_info, disease, new_mod_pos, total_genes, total_genes_without_modifiers, i ,ppi_box2)

	queue.put(rank)


def __saveRank(prank, output):



	f=open(output,'w')
	f.write('Pos rank\n')

	for index in range( len(prank)):
		s = str(prank[index]) 
		f.write(s+'\n')
	f.close()

def __int2strList(v):
	r = []
	for g in v:
		r.append( str(g).strip() )
	return r

def __countOverlaps(a, b):
	cnt = 0
	for v in a:
		if v in b:
			cnt = cnt + 1
	return cnt




def mainJob(ppi_info, disease, mod_pos, total_genes, total_genes_without_modifiers, iteration_th, ppi_box):





	global RND_ID, TEMP_PATH, P_VALUE, DEGREE

	ppi_box2 = []
	ppi_box2 = ppi_box

	total_gene_names = ppi_box[0].getWholeGeneList()
	total_gene_names.sort()

	new_mod_pos, test_pos, test_others = __sampleModifiers(total_genes_without_modifiers, mod_pos)




	if disease.upper().find('RANDOM')>0 :


		# build random networks
		'''
        a = PPI.PPIDB()
        a.interDB = copy.deepcopy( ppi_box[0].interDB )
        a.scoreDB = copy.deepcopy( ppi_box[0].scoreDB )

        #print 'randomize the network'
        a.randomizeNetwork()

        ppi_box2.append(a)


        for i in range( 1, len(ppi_box)):
            ppi_box2.append(ppi_box[i])
        '''


		#index =random.randint(1, len(ppi_box)-1)
		#ppi_box2.append(  ppi_box[index]    )


		#[ output, index1, index2, score_index, threshold   ] = ppi_box[index]        
		#p = PPI.PPIDB()
		#p.loadInteractionFile(output, index1, index2, score_index, -1) # no threshold
		#p = __getRandomNetwork()
		#ppi_box2.append(p)

		# ---------------------------------------------------------------------------------
		# select random modifiers
		basket = getBasket(total_genes, ppi_box[0])


		#new_mod_pos, new_mod_neg = __randomSampling (total_genes, len(new_mod_pos)+1, len(new_mod_neg)+1)    
		new_mod_pos = __randomSampling (total_genes, mod_pos, basket, ppi_box[0])

		total_genes_without_modifiers = __getTotalGenesWithoutModifiers(total_genes, new_mod_pos)
		while(True):
			x_mod_pos, X_test_pos, x_test_others = __sampleModifiers(total_genes_without_modifiers, new_mod_pos)
			if not (test_pos in x_mod_pos or test_pos in x_test_others):
				new_mod_pos = x_mod_pos
				test_others = x_test_others
				test_pos = X_test_pos

				break

	else:
		ppi_box2 = ppi_box


	train_output = TEMP_PATH + '/NetProp_' + RND_ID+'__trained_'+str(iteration_th).strip()+disease + '_' + str(iteration_th).strip()+'.txt'
	test_set_all = [test_pos ] +  test_others

	__train(ppi_info, total_genes, new_mod_pos, train_output, test_set_all, ppi_box2)


	train_output2 = train_output + '_gene.txt'
	r1 = __test( train_output2, test_pos, test_others  )


	global REMOVE_TEMP_FILES
	if REMOVE_TEMP_FILES:
		try:
			os.remove(train_output)        
		except:
			pass

		try:
			os.remove(train_output2)
		except:
			pass


	return r1


def getBasket(total_genes, ppi):

	box = {}

	for g in total_genes:
		deg = len( ppi.getPartners(g)  )
		deg_str = str(deg).strip()

		if not box.has_key(deg_str):
			box[deg_str] = []

		if not g in box[deg_str]:
			box[deg_str].append(g)

	return box

















def __removeUnannotatedModifiers(total_genes, mods):
	r = {}
	for k in mods:
		if k in total_genes:
			r[k] = None
	return r.keys()

def __getTotalGeneList(ppi_box):

	ppi = ppi_box[0]

	total_gene_list = []
	for g in ppi.getWholeGeneList():
		total_gene_list.append(  g )
	return total_gene_list

def __sampleModifiers(total_genes_without_modifiers, mod_pos):


	mod_pos_list = copy.deepcopy(mod_pos)


	pos_1 = random.sample( mod_pos_list, 1) [0]
	mod_pos_list.remove (pos_1)
	other_genes = random.sample( total_genes_without_modifiers, 99 )

	return mod_pos_list, pos_1, other_genes




def __getTotalGenesWithoutModifiers(total_genes, mod_pos ):
	r = []
	for t in total_genes:
		if not t in mod_pos and not t in r:
			r.append(t)
	return r



def __getAvgGeneNumberPerOneGO(go_box):

	r = []
	for goid in go_box:
		r.append(   len(go_box[goid]))

	avg = statistics.average(r)
	std = statistics.stdev(r)

	return avg, std

def __linkGOandGenes(go_box, genes):


	for g in genes.keys():
		for goid in genes[g]:
			if go_box.has_key(goid):
				if not g in go_box[goid]:
					go_box[goid].append(g)



	new_go_box = {}
	for goid in go_box.keys():
		if len(go_box[goid])>0:
			new_go_box [goid] = go_box[goid]

	return new_go_box



def __split(go_box):

	fly = {}
	yeast = {}
	worm = {}

	for gid in go_box.keys():
		for g in go_box[gid]:
			if g[0] == 'F':
				if not fly.has_key(gid):
					fly[gid] = []
				fly[gid].append(g)
			elif g[0] == 'S':
				if not yeast.has_key(gid):
					yeast[gid] = []
				yeast[gid].append(g)
			else:
				if not worm.has_key(gid):
					worm[gid] = []
				worm[gid].append(g)

	return fly, yeast, worm








def __initializeGO(go_class):
	r = {} # go_id as a key, values are gene id
	go_ids = go_class.go_terms.keys()
	for g in go_ids:
		r[g] = []
		#print g
	return r
















# =============================================


def __train(ppi_info, total_genes, mod_pos, output, test_set_all, ppi_box):




	global MP_PROCESSORS

	# score
	pvalues =None
	if MP_PROCESSORS>1:
		pvalues = __calculatePvalueMP(ppi_info, total_genes, test_set_all, mod_pos, output, ppi_box)
	else:
		pvalues = __calculatePvalue(ppi_info, total_genes, test_set_all, mod_pos, output, ppi_box)

	#__save(go, go_box, mod_pos, mod_neg, pvalues, output)

	__resave(total_genes, mod_pos, pvalues, output, test_set_all, '')




def __resave(total_genes, mod_pos, pvalues, output, test_set_all, md5):



	ofile = output + '_gene.txt'


	box = __rebox(test_set_all, mod_pos, pvalues, output )
	__save2(box, mod_pos, ofile, md5)

def __save2(box, mod_pos, output, md5):


	#print '[NetProp] Sorting results... '
	box3 = sorted( box.items(), key=lambda (k,v): v[0], reverse = True)

	#print '[NetProp] Saving final scores...'

	f=open(output, 'w')
	f.write('!'+md5+'\n')
	f.write('Gene\tscore\tDesc\n')


	for k, v in box3:
		score = v[0]
		desc = v[1]

		if k is None:
			print 'found None'
			continue

		#print k, score, desc

		s = __addTag(k, mod_pos)+k+'\t'+str(score)+'\t'+desc
		f.write(s+'\n')


	f.close()


def __addTag(gene, mod_pos):
	tag = ''
	if gene in mod_pos:
		tag = 'm_'

	return tag

def __rebox(total_genes, mod_pos, pvalues, output):

	#print '[NetProp] processing resulting scores...'

	# calculate score !
	# 

	# ==================================================================================================
	box = {}


	score = 0.0
	desc = []
	for g in total_genes:

		if g is None: continue  # ?????


		desc = [ _0_Preprocess.ID_BOX.getCommonNameOf(g)   ]

		p_pos = pvalues[g]

		desc.append( '[' + str(p_pos)  +  ' ] ')
		score = p_pos
		#score = score - 2.0 * math.log(p_pos) 


		box[g] = [ score, ','.join(desc)] 


	return box


def __save(go, go_box, mod_pos, mod_neg, pvalues, output):

	f=open(output,'w')

	for index in range(2):

		m_v = pvalues[index]


		s = 'GO_ID\tp-value\t#gene\t#modifier\tGO Desc\tList'
		f.write(s+'\n')

		# vx = sorted( go_box.items(), key=lambda (k,v): m_v[k], reverse=False)

		for gid in go_box.keys():
			genes = go_box[gid]

			go_t = go.getGOTermOf(gid)
			name = go_t.getGOName()

			lst = ','.join( __makeGeneList(genes, mod_pos, mod_neg) )

			tv = len(genes)

			mv = None

			if index == 0:
				mv = __countModifierNumber(mod_pos, genes)
			else:
				mv = __countModifierNumber(mod_neg, genes)

			pv = m_v[gid]

			s = gid + '\t' + str(pv) + '\t' + str(tv) + '\t' + str(mv) + '\t' + name + '\t' + lst
			f.write(s + '\n')

	f.close()


def __makeGeneList(gene_list, mod_pos, mod_neg):
	r = []
	for g in gene_list:
		if g in mod_pos:
			r.append( 'm_' + g )
		elif g in mod_neg:
			r.append( 'x_' + g )
		else:
			r.append ( g )
	return r        

def __getCommonValuesInTwoLists(l1, l2):

	r = []
	for l in l1:
		if l in l2:
			r.append(l)
	return r


def __calculate1(g1, g2, intersection, ppi ):

	r = 0.0
	for g in intersection:
		score1 = float(ppi.getScore(g1, g)) 
		score2 = float(ppi.getScore(g2, g)) 
		r = r + score1 * score2

	r=r*2.0
	return r

def __calculate2(g1, intersection, ppi):
	r = 0.0
	for g in intersection:

		score = float(ppi.getScore(g1, g)) # / 100.0
		r = r + score    
	return r
def __calculate3(g1, g2, intersection, ppi ):

	r = 0.0
	for g in intersection:
		score1 = float(ppi.getScore(g1, g)) # / 100.0
		score2 = float(1.0 - ppi.getScore(g2, g)) # / 100.0
		r = r + score1 * score2

	r=r*2.0
	return r    


def __diffList(l1, l2):
	# L1 - L2 
	r = []
	for g in l1:
		if not g in l2:
			r.append(g)
	return r

def __calculateLambda(u, v, u_partners, v_partners, intersection):

	n_avg = float((len(u_partners) - 1 + len(v_partners) - 1 ))/2.0  # - self
	v = n_avg - float( len(__diffList(u_partners, v_partners)) + len(intersection)  )
	return max(0.0, float(v))

def __calculateFunctionalSimilarity(u, v, confidence, ppi):


	u_partners = ppi.getPartners(u) #+ [ u ]
	v_partners = ppi.getPartners(v) #+ [ v ]


	intersection = __getCommonValuesInTwoLists(u_partners, v_partners)


	if len(intersection)>0:

		#print repr(u_partners)
		#print repr(v_partners)
		#print repr(intersection)

		v1 = __calculate1(u, v, intersection, ppi)
		v2 = __calculate2(u, intersection, ppi)
		v3 = __calculate3(u, v, intersection, ppi)

		v4 = __calculate2(v, intersection, ppi)
		v5 = __calculate3(v, u, intersection, ppi)

		lam1 = __calculateLambda(u, v, u_partners, v_partners, intersection)
		lam2 = __calculateLambda(v, u, v_partners, u_partners, intersection)


		#print v1, v2, v3, v4, v5, lam1, lam2
		score1 = v2+v3+v1+lam1 
		score2 = v4+v5+v1+lam2

		score = 0.0
		if score2>0 and score1>0:
			score = v1/score1*v1/score2

			# confidence
			score = score * confidence

			#print u, v, score

		'''
                  v1                          v1
        -----------------------  * --------------------------
         (v2 + v3) + v1 + lam1      (v4 + v5) + v1 + lam2

        '''

		return score
	else:
		return 0.0


def __getReliabilitySum(gene, ppi):

	v = 0.0
	for pg in ppi.getPartners(gene):
		v = v + ppi.getScore(gene, pg)
	return v



def __randomSamplingBasedOnReliability(total_genes, test_set_gene, mod_pos, ppi):

	# modifier
	#    1. test_set_gene --> X
	#    2. mod_pos -- > X

	new_mod = []

	resvoir = []
	for g in total_genes:
		if not g in test_set_gene and not g in mod_pos:
			resvoir.append(g)


	for m in mod_pos:

		m_reliability = __getReliabilitySum( m , ppi)

		other_min_reliability = 10000000.0 
		for m2 in mod_pos:

			if m2 == m: 
				continue

			m2_reliability = __getReliabilitySum(m2, ppi)
			if abs(m_reliability - m2_reliability)< abs(m_reliability - other_min_reliability):
				other_min_reliability = m2_reliability


		while(True):

			rg = random.sample(  resvoir , 1 ) [0] 

			if not rg in new_mod:

				rg_reliability = __getReliabilitySum(rg, ppi)
				if abs( rg_reliability - m_reliability  ) <= abs( rg_reliability - other_min_reliability  ):

					if not rg in new_mod:
						new_mod.append(rg)
						break


	return new_mod


def for_multi_processing(args):

	# MULTIPROCESS에서  사용할 함수 제일 끝에는 multiprocess.Queue() 변수를 넣어줘야함.
	fname, total_genes, mod_pos, cutoff, output, r , scores_queue = args

	global Pnp_dict, Pnp_matrix

	ppi = PPI.PPIDB_STRING()
	ppi.load(fname, 0, 1, 2, cutoff, ' ', None)

	fcache = fname + '.cache'
	Pnp_matrix, Pnp_dict = __buildPnpMatrix(total_genes, ppi, cache_file=fcache)




	sc = __calculatePvaluSUB(total_genes, mod_pos, output, ppi, r)
	scores_queue.put(sc)

	logmsg('Results added into the queue: ' + fcache)

	#return sc



def __calculatePvalueMP(total_genes, test_set_gene, mod_pos, output, ppi_box):




	global P_VALUE, PRELOAD_NO_OF_RANDOM_NETWORKS, Pnp_dict, Pnp_matrix, fcache

	global R_PARAMETER, MP_PROCESSORS
	r = R_PARAMETER


	logmsg('calculatePvalueMP ============')

	fcache = './7. Network Propagation/' + \
	        os.path.basename (_0_Preprocess.STRING_FILE) + \
	        '_' + str( _0_Preprocess.STRING_CUTOFF).strip() + '.cache'    
	Pnp_matrix, Pnp_dict = __buildPnpMatrix(total_genes, ppi_box[0], cache_file=fcache)


	cutoff = _0_Preprocess.STRING_CUTOFF




	logmsg('Calculate scores -------------')
	real_scores = __calculatePvaluSUB(total_genes, mod_pos, output, ppi_box[0], r)

	if P_VALUE == False:
		return real_scores

	else:

		logmsg('calculate random scores -----------')
		rand = {}

		r_nets = __getRandomNetworks(PRELOAD_NO_OF_RANDOM_NETWORKS)

		box = []

		for i in range(PRELOAD_NO_OF_RANDOM_NETWORKS):

			fname = r_nets[i]

			th = MULTIPROCESS.JOB_THREAD()
			#th.set_args(for_multi_processing, [fname, test_set_gene, mod_pos, cutoff, output, r, queues])
			th.set_args(for_multi_processing, [fname, test_set_gene, mod_pos, cutoff, output, r])
			box.append(th)


		# waiting....


		logmsg('Running multiprocesses = ' + str(len(box)))

		ret = MULTIPROCESS.runMultiprocesses(box, max_cpu=MP_PROCESSORS, sleep_interval_seconds=10)


		logmsg ("Retreiving results from the queue")
		for r in ret:

			scores = r

			for key in scores.keys():
				if not rand.has_key(key):
					rand[key] = []

				rand[key].append( scores[key]   )


		pvalues = {}
		for g in rand.keys():
			avg = statistics.average( rand[g])
			std = statistics.stdev( rand[g])
			p = 1.0
			try:
				if std != 0.0:
					# z score 이용
					#p = abs( float(( avg - real_scores[g] )/ float(std) ) )
					p = statistics.getNormalDistPvalue(avg, std, real_scores[g])
				else:
					p = 1.0
			except:
				p = 1.0

			pvalues[g]=p

		return pvalues # p-value based on z-score    



	'''
    elif P_VALUE == True:

        # r= 0.0, reference
        ppi = ppi_box[0].randomizeNetwork()
        zero_scores = __calculatePvaluSUB(test_set_gene, mod_pos, output, ppi, 0.0)

        r = {}

        for k in real_scores.keys():
            a = real_scores[k]
            b = zero_scores[k]

            if b == 0.0 or a==0.0:
                r[k] = 0.0
            else:
                v = math.log(  a/b    )
                r[k] = v
        return r
    '''


'''
def __getCacheFileName():


    fcache = None

    # load cache
    fname =  os.path.basename (_0_Preprocess.STRING_FILE) + \
        '_' + str( _0_Preprocess.STRING_CUTOFF).strip() + '.cache.txt'

    temp_file = _0_Preprocess.TEMP_PATH + '/' + fname

    if os.path.exists(temp_file):
        fcache = temp_file
    else:

'''


def __calculatePvalue(ppi_info, total_genes, test_set_gene, mod_pos, output, ppi_box):

	global P_VALUE, PRELOAD_NO_OF_RANDOM_NETWORKS, Pnp_dict, Pnp_matrix, fcache

	global R_PARAMETER
	r = R_PARAMETER


	interaction_file = ppi_info[0]
	index1 = ppi_info[1]
	index2 = ppi_info[2]
	score_index = ppi_info[3]
	threshold = ppi_info[4]
	separator = ppi_info[5]


	# load cache

	fcache = './7. Network Propagation/' + \
	        os.path.basename ( interaction_file ) + \
	        '_' + str( threshold ).strip() + '.cache.txt'

	Pnp_matrix, Pnp_dict = __buildPnpMatrix(total_genes, ppi_box[0], cache_file=fcache)



	#fcache = './7. Network Propagation/' + \
	#        os.path.basename (_0_Preprocess.STRING_FILE) + \
	#        '_' + str( _0_Preprocess.STRING_CUTOFF).strip() + '.cache.txt'

	#fcache = copyCache()
	#if fcache is None:
	#    fcache = os.path.basename (_0_Preprocess.STRING_FILE) + \
	#    '_' + str( _0_Preprocess.STRING_CUTOFF).strip() + '.cache.txt'


	print 'Cache file = ', fcache


	real_scores = __calculatePvaluSUB(test_set_gene, mod_pos, output, ppi_box[0], r)
	
	print 'zero=00000000'
	zero_scores = __calculatePvaluSUB(test_set_gene, mod_pos, output, ppi_box[0], 0.0)
	pvalues = {}
	
	for g in real_scores.keys():
		if zero_scores.has_key(g):
			pvalues[g] = real_scores[g] - zero_scores[g]
		else:
			pvalues[g] = real_scores[g]
	
	return pvalues



	'''
	if P_VALUE == False:
		
		
		
		#return real_scores


		# 만약 r=0일떄와 비교해서 score 차이를 계산한다면?
		# random network랑 비교하니 영 아니게 결과가 나온다.


		print 'zero=00000000'
		zero_scores = __calculatePvaluSUB(test_set_gene, mod_pos, output, ppi_box[0], 0.0)
		pvalues = {}

		for g in real_scores.keys():
			if zero_scores.has_key(g):
				pvalues[g] = real_scores[g] - zero_scores[g]
			else:
				pvalues[g] = real_scores[g]

		return pvalues




	else:

		rand = {}

		r_nets = __getRandomNetworks(PRELOAD_NO_OF_RANDOM_NETWORKS)

		for i in range(PRELOAD_NO_OF_RANDOM_NETWORKS):

			fname = r_nets[i]

			print '[RND NETWORK]', (i+1), '/', PRELOAD_NO_OF_RANDOM_NETWORKS, fname

			ppi = PPI.PPIDB_STRING()
			ppi.load(fname, 0, 1, 2, _0_Preprocess.STRING_CUTOFF, 
			         ' ', None)

			#tmp_ppi = ppi_box[0].randomizeNetwork()
			fcache = fname + '.cache'
			Pnp_matrix, Pnp_dict = __buildPnpMatrix(total_genes, ppi, cache_file = fcache)


			scores = __calculatePvaluSUB(test_set_gene, mod_pos, output, ppi, r)

			for key in scores.keys():
				if not rand.has_key(key):
					rand[key] = []

				rand[key].append( scores[key]   )


		pvalues = {}
		for g in rand.keys():
			avg = statistics.average( rand[g])
			std = statistics.stdev( rand[g])
			p = 1.0
			try:
				if std != 0.0:
					# z score 이용
					#p = abs( float(( avg - real_scores[g] )/ float(std) ) )
					p = statistics.getNormalDistPvalue(avg, std, real_scores[g])
				else:
					p = 1.0
			except:
				p = 1.0

			pvalues[g]=p

		return pvalues # p-value based on z-score    
	'''


	'''
    elif P_VALUE == True:

        # r= 0.0, reference
        ppi = ppi_box[0].randomizeNetwork()
        zero_scores = __calculatePvaluSUB(test_set_gene, mod_pos, output, ppi, 0.0)

        r = {}

        for k in real_scores.keys():
            a = real_scores[k]
            b = zero_scores[k]

            if b == 0.0 or a==0.0:
                r[k] = 0.0
            else:
                v = math.log(  a/b    )
                r[k] = v
        return r
    '''




def __getTrim(total_genes, mod_pos, ppi):

	total2 = []
	score = {}
	for g in total_genes:
		score[g] = 0.0 # init

	for g in mod_pos.keys():
		for g2 in ppi.getPartners(g):

			if not g2 in total2:
				total2.append(g2)

			for g3 in ppi.getPartners(g2):
				if not g3 in total2:
					total2.append(g3)


	return score, total2

def __getW(u, ppi):

	total = 0.0
	for k in ppi.getPartners(u):
		if k != u:
			total = total + ppi.getScore(k, u)
	return total  

def __sortStr(a, b):
	if a>b:
		return a+':'+b
	else:
		return b+':'+a


def __buildPnpMatrix(total_gene_names, ppi, cache_file = None):


	global USE_CACHE, MP_PROCESSORS

	dic = {}

	t_genes = ppi.getWholeGeneList()

	# if USE_CACHE is True, use a cache file. 
	if USE_CACHE and cache_file is not None:


		logmsg('cache_file = ' + cache_file)

		if os.path.exists(cache_file):

			if MP_PROCESSORS == 1:
				print 'Loading cached matrix: ', cache_file

			f=open(cache_file,'r')
			
			while(True):
				lines = f.readlines(5000)
				if not lines: break
				
				for s in lines:
					s=s.replace('\n','').strip()
					if len(s)>0:
						x = s.split('\t')
						key = x[0]
						v=eval(x[1])
						dic[key]=v
			f.close()

			#if MP_PROCESSORS == 1:
			print 'Cached data = ', len(dic)

			logmsg('cached_data = ' + str(len(dic)))

			return None, dic
		else:

			#if MP_PROCESSORS == 1:
			print 'Cache file missing: ', cache_file



	logmsg('Creating cache file...')
	# chaching
	i = 0.0
	q = float( len(t_genes) )

	for u in t_genes:

		i+=1.0

		#if MP_PROCESSORS == 1:
		print '\r[CACHE]', (i/q*100.0), '%',


		W_u = __getW(u, ppi)

		partners = ppi.getPartners(u)



		#for v in total_gene_names:
		for v in partners:


			key = __sortStr(u, v)

			if u == v:
				pass
			else:
				W_v = __getW(v, ppi)
				w_uv = ppi.getScore(u, v)
				s = w_uv/math.sqrt( W_u * W_v)

				dic[key]=s

	#if MP_PROCESSORS == 1:
	print ''
	print '# of cached scores = ', len(dic)


	logmsg('# of cached scores = ' + str(len(dic)))


	f=open(cache_file,'w')
	for k in dic.keys():
		s = k+'\t'+str(dic[k]).strip()
		f.write(s+'\n')
	f.close()


	#copyCache()

	logmsg('Cache saved = ' + cache_file)


	return None, dic

def runFinal(ppi_info, interaction_file, index1, index2, score_index, threshold, separator):



	global RND_ID, NETWORK_TEMP_FILES, MP_PROCESSORS, Pnp_dict, Pnp_matrix
	global fcache, MP_PROCESSORS, REMOVE_TEMP_FILES, P_VALUE, PRELOAD_NO_OF_RANDOM_NETWORKS
	global HOW_MANY_TEST_MODIFERS, DEGREE, TEMP_PATH



	Pnp_matrix = None
	Pnp_dict = {}


	#int_file = os.path.basename(interaction_file)
	#fcache = './7. Network Propagation/' + int_file + '_' + str(threshold).strip() + '.cache.txt'
	fcache = './7. Network Propagation/' + \
	        os.path.basename ( interaction_file ) + \
	        '_' + str( threshold ).strip() + '.cache.txt'




	mod_maps, mods = __getModifiers()


	dd = mods.keys()
	dd.sort()


	NETWORK_TEMP_FILES = [interaction_file, index1, index2, score_index, threshold]  


	info = '--------------------------------------------' + '\n' + \
	        'THREAD = ' + str( MP_PROCESSORS ) + '\n' + \
	        'Remove temporary files = ' + str( REMOVE_TEMP_FILES ) + '\n' + \
	        'Calculate p-values = ' + str( P_VALUE) + '\n' + \
	        '\t# of random networks = ' + str( PRELOAD_NO_OF_RANDOM_NETWORKS) + '\n' + \
	        'Iteration = ' + str( HOW_MANY_TEST_MODIFERS) + '\n' + \
	        'Method for randomization = ' + str(DEGREE) + '\n' + \
	        'Interaction file = ' + interaction_file + '\n' + \
	        '\tcutoff='+str(threshold) + '\n'
	print info


	# ------------------------------------------------------------------
	dd2 = []
	# if cache exists, skip it



	for disease in dd:
		if disease.find('RANDOM')>=0:
			dd2.append(disease)
			continue

		#output = _0_Preprocess.OUTPUT_FOLDER + '/NetPropagation_'+disease+'_'+  os.path.basename(interaction_file) + '_(' + str(threshold) + ')_' + RND_ID + '.txt'
		output = _0_Preprocess.TEMP_PATH + '/NetPropagation_' + disease + '_' + os.path.basename(interaction_file) + '_(' + str(threshold) + ')_' + RND_ID + '.txt'
		
		

		# md5 within the output file

		x = sorted( copy.deepcopy( mods[disease] ) )
		x.append( interaction_file )
		x.append( str(  threshold ))

		s = '\n'.join( x )
		md5 = _0_Preprocess.generateMD5(s)


		print 'Modifiers MD5=', md5
		if md5 == _0_Preprocess.getMD5(output+'_gene.txt') and _0_Preprocess.REFRESH_RESULTS == False:
			print '\tUse previous results: ', output+'_gene.txt'
			#print 'MD5 = ', md5[disease]

		else:
			dd2.append(disease)



	if len(dd2) == 0:
		print '======================================'
		print ' [NetProp]   All cached before '
		print '======================================'
		return

	# ------------------------------------------------------------------







	ppi_box = loadInteraction(interaction_file, PRELOAD_NO_OF_RANDOM_NETWORKS, index1, index2, score_index, threshold, separator) # 100개 네트워크를 만든다.
	# ppi_box[0] 

	total_genes = __getTotalGeneList(ppi_box)






	print 'init...',
	#global Pnp_matrix, Pnp_dict
	total_gene_names = ppi_box[0].getWholeGeneList()
	total_gene_names.sort()

	#Pnp_matrix, Pnp_dict = __buildPnpMatrix(total_genes, ppi_box[0], cache_file=fcache)

	#fcache = './7. Network Propagation/' + \
	#    os.path.basename (_0_Preprocess.STRING_FILE) + \
	#    '_' + str( _0_Preprocess.STRING_CUTOFF).strip() + '.cache.txt'
	#Pnp_matrix, Pnp_dict = __buildPnpMatrix(total_genes, ppi_box[0], cache_file=fcache)

	#Pnp_matrix, Pnp_dict = __buildPnpMatrix(total_gene_names, ppi_box[0])
	print 'done'

	init = True

	for disease in dd2:



		if disease.find('RANDOM')>=0:
			mods[disease] = random.sample(total_genes, len(mods[disease]))

		print '--------------------------------------'        
		print disease

		#output = _0_Preprocess.OUTPUT_FOLDER + '/NetPropagation_'+disease+'_'+  os.path.basename(interaction_file) + '_(' + str(threshold) + ').txt'
		output = _0_Preprocess.TEMP_PATH + '/NetPropagation_' + disease + '_' + os.path.basename(interaction_file) + '_(' + str(threshold) + ').txt'

		# md5 within the output file

		x = sorted( copy.deepcopy( mods[disease] ) )
		x.append( interaction_file )
		x.append( str(  threshold ))

		s = '\n'.join( x )
		md5 = _0_Preprocess.generateMD5(s)


		print 'Modifiers MD5=', md5
		if md5 == _0_Preprocess.getMD5(output+'_gene.txt') and _0_Preprocess.REFRESH_RESULTS == False:
			print '\tUse previous results: ', output+'_gene.txt'
			#print 'MD5 = ', md5[disease]

			continue


		#Pnp_matrix = copy.deepcopy(local_Pnp_matrix)
		#Pnp_dict = copy.deepcopy(local_Pnp_dict)


		#if init:
		#    USE_CACHE = False ##############################
		#    init = False
		#else:
		#    USE_CACHE = True




		ppi_box2 = []
		#ppi_box2 = copy.deepcopy( ppi_box )
		ppi_box2 = ppi_box


		mod_pos = mods[disease] 

		new_mod_pos = __removeUnannotatedModifiers(total_genes, mod_pos)

		total_genes_without_modifiers = __getTotalGenesWithoutModifiers(total_genes, mod_pos)


		s = 'Interaction = ' + interaction_file + '\n' + \
		        'Disease = ' + disease + '\n' + \
		        'Total annotated proteins # = ' + str(len(total_genes)) + '\n' + \
		        'Modifiers # = ' + str( len(new_mod_pos) ) + '\n' 

		print s





		test_set_all = total_genes

		__train2(ppi_info, total_genes, mod_pos, output, test_set_all, ppi_box, md5)

		print '[OUTPUT] = ', output



def saveVariableAsPickle(variable):

	global TEMP_PATH

	#print 'temp=', TEMP_PATH
	#print 'RND=',  MyUtil_pypy.getRandomString(10)

	fname = TEMP_PATH + '/pickle_' + MyUtil_pypy.getRandomString(10) + '.obj'

	#print fname

	with open(fname, "wb") as of:
		cPickle.dump(variable, of)

	return fname

def loadVariableFromPickle(fname):

	e = None

	with open(fname, 'rb') as of:
		e = cPickle.load(of)

	return e

def runFinal_mp(config_file, ppi_info, interaction_file, index1, index2, score_index, threshold, separator, max_thread= 5):




	global RND_ID, NETWORK_TEMP_FILES, MP_PROCESSORS, OPTIONS
	global fcache, REMOVE_TEMP_FILES, P_VALUE, PRELOAD_NO_OF_RANDOM_NETWORKS
	global HOW_MANY_TEST_MODIFERS, DEGREE

	#int_file = os.path.basename(interaction_file)
	#fcache = './7. Network Propagation/' + int_file + '_' + str(threshold).strip() + '.cache.txt'



	mod_maps, mods = __getModifiers()
	
	opt = OPTIONS


	dd = mods.keys()
	dd.sort()


	NETWORK_TEMP_FILES = [interaction_file, index1, index2, score_index, threshold]  


	info = '--------------------------------------------' + '\n' + \
	        'THREAD = ' + str( max_thread ) + '\n' + \
	        'Remove temporary files = ' + str( REMOVE_TEMP_FILES ) + '\n' + \
	        'Calculate p-values = ' + str( P_VALUE) + '\n' + \
	        '\t# of random networks = ' + str( PRELOAD_NO_OF_RANDOM_NETWORKS) + '\n' + \
	        'Iteration = ' + str( HOW_MANY_TEST_MODIFERS) + '\n' + \
	        'Method for randomization = ' + str(DEGREE) + '\n' + \
	        'Interaction file = ' + interaction_file + '\n' + \
	        '\tcutoff='+str(threshold) + '\n'
	print info





	dd2 = []

	for disease in dd:



		if disease.find('RANDOM')>=0:
			dd2.append(disease)
			continue

		#output = _0_Preprocess.OUTPUT_FOLDER + '/NetPropagation_'+disease+'_'+  os.path.basename(_0_Preprocess.STRING_FILE) + '_(' + str(_0_Preprocess.STRING_CUTOFF) + ')_' + RND_ID+'.txt'
		output = _0_Preprocess.TEMP_PATH + '/NetPropagation_' + disease + '_' + os.path.basename(_0_Preprocess.STRING_FILE) + '_(' + str(_0_Preprocess.STRING_CUTOFF) + ')_' + RND_ID + '.txt'
		
		

		# md5 within the output file

		x = sorted( copy.deepcopy( mods[disease] ) )
		x.append( str(_0_Preprocess.STRING_FILE))
		x.append( str(_0_Preprocess.STRING_CUTOFF))

		s = '\n'.join( x )
		md5 = _0_Preprocess.generateMD5(s)


		print 'Modifiers MD5=', md5
		if md5 == _0_Preprocess.getMD5(output+'_gene.txt') and _0_Preprocess.REFRESH_RESULTS == False:
			print '\tUse previous results: ', output+'_gene.txt'
			#print 'MD5 = ', md5[disease]
			continue
		else:
			dd2.append(disease)

	if len(dd2) == 0:
		print '======================================'
		print ' [NetProp]   All cached before '
		print '======================================'
		return











	ppi_box = loadInteraction(interaction_file, PRELOAD_NO_OF_RANDOM_NETWORKS, index1, index2, score_index, threshold, separator) # 100개 네트워크를 만든다.
	# ppi_box[0] 


	ppi_box_file = saveVariableAsPickle(ppi_box)


	total_genes = __getTotalGeneList(ppi_box)






	print 'init...',
	global Pnp_matrix, Pnp_dict
	total_gene_names = ppi_box[0].getWholeGeneList()
	total_gene_names.sort()

	#fcache = './7. Network Propagation/' + \
	#    os.path.basename (_0_Preprocess.STRING_FILE) + \
	#    '_' + str( _0_Preprocess.STRING_CUTOFF).strip() + '.cache.txt'
	#Pnp_matrix, Pnp_dict = __buildPnpMatrix(total_genes, ppi_box[0], cache_file=fcache)

	#Pnp_matrix, Pnp_dict = __buildPnpMatrix(total_gene_names, ppi_box[0])
	print 'done'

	init = True


	jobs = []
	outputs = []

	for disease in dd2:



		if disease.find('RANDOM')>=0:
			mods[disease] = random.sample(total_genes, len(mods[disease]))

		print '--------------------------------------'        
		print disease

		#output = _0_Preprocess.OUTPUT_FOLDER + '/NetPropagation_'+disease+'_'+  os.path.basename(_0_Preprocess.STRING_FILE) + '_(' + str(_0_Preprocess.STRING_CUTOFF) + ')_' + RND_ID + '.txt'
		output = _0_Preprocess.TEMP_PATH + '/NetPropagation_' + disease + '_' + os.path.basename(_0_Preprocess.STRING_FILE) + '_(' + str(_0_Preprocess.STRING_CUTOFF) + ')_' + RND_ID + '.txt'

		# md5 within the output file

		x = sorted( copy.deepcopy( mods[disease] ) )
		x.append( str(_0_Preprocess.STRING_FILE))
		x.append( str(_0_Preprocess.STRING_CUTOFF))

		s = '\n'.join( x )
		md5 = _0_Preprocess.generateMD5(s)


		print 'Modifiers MD5=', md5
		if md5 == _0_Preprocess.getMD5(output+'_gene.txt') and _0_Preprocess.REFRESH_RESULTS == False:
			print '\tUse previous results: ', output+'_gene.txt'
			#print 'MD5 = ', md5[disease]

			continue


		#Pnp_matrix = copy.deepcopy(local_Pnp_matrix)
		#Pnp_dict = copy.deepcopy(local_Pnp_dict)


		#if init:
		#    USE_CACHE = False ##############################
		#    init = False
		#else:
		#    USE_CACHE = True




		#ppi_box2 = []
		#ppi_box2 = copy.deepcopy( ppi_box )


		mod_pos = mods[disease] 

		new_mod_pos = __removeUnannotatedModifiers(total_genes, mod_pos)

		total_genes_without_modifiers = __getTotalGenesWithoutModifiers(total_genes, mod_pos)


		s = 'Interaction = ' + interaction_file + '\n' + \
		        'Disease = ' + disease + '\n' + \
		        'Total annotated proteins # = ' + str(len(total_genes)) + '\n' + \
		        'Modifiers # = ' + str( len(new_mod_pos) ) + '\n' 

		print s





		test_set_all = total_genes

		#__train2(total_genes, mod_pos, output, test_set_all, ppi_box, md5)
		th = MULTIPROCESS.JOB_THREAD()
		#th.set_args(__train2_mp, [config_file, ppi_info, total_genes, mod_pos, output, test_set_all, ppi_box, md5, opt])
		
		th.set_args(__train2_mp, [config_file, ppi_info, total_genes, mod_pos, output, test_set_all, ppi_box_file, md5, opt])

		jobs.append(th)


		outputs.append(output)
		#print '[OUTPUT] = ', output


	MULTIPROCESS.runMultiprocesses(jobs, max_cpu=max_thread, sleep_interval_seconds=10)


	try:
		os.remove(ppi_box_file)
	except:
		pass

	print '========================================='
	for o in outputs:
		print '[OUTPUT] = ', o
	print '========================================='


def __train2_mp(args):

	config_file, ppi_info, total_genes, mod_pos, output, test_set_all, ppi_box_file, md5, opt, queue = args

	_0_Preprocess.init(config_file)
	
	RND_ID = opt[_0_Preprocess.OPT_RND_ID]
	_0_Preprocess.RND_ID = RND_ID

	ppi_box = loadVariableFromPickle(ppi_box_file)

	# run train2
	# save this function stores all the results, it is not necessary to return any results
	__train2(ppi_info, total_genes, mod_pos, output, test_set_all, ppi_box, md5)
	queue.put(1)
	
	return 



def __train2(ppi_info, total_genes, mod_pos, output, test_set_all, ppi_box, md5):


	pvalues = __calculatePvalue(ppi_info, total_genes, test_set_all, mod_pos, output, ppi_box)


	'''
	# score
	pvalues =None
	if MAX_THREAD>1:
		pvalues = __calculatePvalueMP(ppi_info, total_genes, test_set_all, mod_pos, output, ppi_box)
	else:
		pvalues = __calculatePvalue(ppi_info, total_genes, test_set_all, mod_pos, output, ppi_box)
	'''

	#__save(go, go_box, mod_pos, mod_neg, pvalues, output)

	__resave(total_genes, mod_pos, pvalues, output, test_set_all, md5)



def __getSeedScore(total_gene_names, mod, ppi):


	r = {}
	arr = []

	for k in total_gene_names:
		r[k] = 0.0



	method = 2 # 1: neurogem confidence, 2: score=1

	if method == 1:
		total = 0.0
		for k in mod:
			total = total + mod[k]

		for k in mod:
			r[k] = mod[k]/total


		for k in total_gene_names:
			arr.append(  [r[k]] )
			#arr.append(  r[k] )

	elif method == 2:
		total = 0.0
		for k in mod:
			total = total + 1.0 # instead of confidence score

		for k in mod:
			r[k] = 1.0/total


		for k in total_gene_names:
			arr.append(  [r[k]] )
			#arr.append(  r[k] )            

	return arr, r


'''
def __calculatePvaluSUB(total_genes, mod_pos, output, ppi, r):    

    global Pnp_matrix

    #print 'p-value, sorting'
    total_gene_names = ppi.getWholeGeneList()
    total_gene_names.sort()

    #print 'init'
    tlist, t0_dict =  __getSeedScore(total_gene_names, mod_pos, ppi)
    tn = numpy.array( tlist )
    rho = numpy.copy(tn)
    #Pnp = __buildPnpMatrix(total_gene_names, ppi)

    #print 'rep'

    while(True):
        t0 = numpy.copy(tn)


        print tn.shape, 
        print Pnp_matrix.shape

        tn = (1.0-r)*Pnp_matrix*tn + r*rho

        t0 = t0 - tn
        #t0.tolist()
        total = 0.0

        for v in range(len(total_gene_names)):
            for k in range(len(total_gene_names)):
                total=total+t0[v,k]

        #print '__converge: ', i, ' total = ', total
        if total == 0.0:
            break

    p = {}
    for i in range(len(total_gene_names)):
        g = total_gene_names[i]
        if g in total_genes:
            p[ g  ] = tn[i]

    return p 

'''


def __calculatePvaluSUB(total_genes, mod_pos, output, ppi, r):

	global Pnp_matrix, Pnp_dict, MP_PROCESSORS


	# total_genes <- not used in this method

	logmsg('sub-initiated')
	cutoff = 1e-10
	max_iteration = 200


	#print 'p-value, sorting'
	total_gene_names = ppi.getWholeGeneList()
	#total_gene_names.sort()

	#print 'init'
	tlist, t0_dict =  __getSeedScore(total_gene_names, mod_pos, ppi)

	rho = copy.deepcopy(t0_dict)

	tn = t0_dict

	cnt = 0
	totals = []

	prev_tn = {}

	while(True):

		cnt = cnt + 1
		if max_iteration>200:
			break

		total = 0.0
		tn_1 = {}
		for k in total_gene_names:
			s = 0.0

			if r!=1.0:
				for pg in ppi.getPartners(k):

					if pg == k: continue

					kk = __sortStr(k, pg) 
					s = s + Pnp_dict[ kk ] * tn[ pg ]


					'''
                    kk = __sortStr(k, pg) 
                    if not dict_cache.has_key(kk):

                        W_u = __getW(k, ppi)
                        W_v = __getW(pg, ppi)
                        w_uv = ppi.getScore(k, pg)
                        ss = w_uv/math.sqrt( W_u * W_v) 

                        dict_cache[kk] = ss

                    s = s + dict_cache[ kk ] * tn[ pg ]
                    '''




				s = s * (1.0-r)

			s = s + r * rho[k]

			tn_1[k] = s


			if cnt == 1: 
				prev_tn [k] = s
				total = total + 1.0
			else:
				total = total + abs( prev_tn[k] - s  )
				prev_tn [k] = s

		'''        
        total = 0.0
        for v in tn_1.keys():
            total = total + tn_1[v]
        '''
		total_s = str( total ).strip()

		if MP_PROCESSORS == 1:
			print '\rDiff = ', total, 'count = ', cnt, '            ',

		tn = tn_1

		if total <= cutoff or cnt>=max_iteration or total_s in totals:
			break
		else:
			if len(totals) == 4:
				del totals[0]
			totals.append(total_s)

	p = {}
	#for i in range(len(total_gene_names)):
	#    g = total_gene_names[i]
	#    if g in total_genes:
	#        p[ g  ] = tn[i]
	for g in total_genes:
		if g is None:
			continue

		try:
			p[g] = tn[g]
		except:
			p[g] = 0.0

		'''
		if tn.has_key(g):
			p[g] = tn[g]
		else:
			p[g] = 0.0
		'''

	logmsg('sub- done')
	#Pnp_dict = dict_cache

	return p 




def __countModifierNumber(mod_list, gene_list):
	cnt = 0
	for g in gene_list:
		if g in mod_list:
			cnt = cnt + 1
	return cnt




def __test( train_output2, test_pos, test_others   ):


	f=open(train_output2, 'r')

	pos_rank = -1

	cnt = 0

	init = True
	for s in f.readlines():
		s = s.replace('\n', '')

		if len(s.strip()) == 0:
			continue
		else:
			if init == True:
				init = False
			else:
				x = s.split('\t')

				gid = x[0].replace('m_', '').replace('x_', '')

				if gid == test_pos:
					cnt += 1
					if pos_rank == -1:
						pos_rank = cnt
					else:
						print 'Error!!!! two or more positives'

				elif gid in test_others:
					cnt += 1





	f.close()


	if cnt != 100:
		print 'Error !!! Not 100 test set', cnt

	return pos_rank 

def __loadSavedFile( infile ):

	f=open( infile, 'r' )
	init = True

	r = {}


	for s in f.readlines():
		if init:
			init = False
		else:
			s = s.replace('\n','')
			x = s.split('\t')

			genes = x[5].split(',')
			score = x[1]
			goid = x[0]
			pvalue = x[2]

			for g in genes:
				if not r.has_key(g):
					r[g] = [ [], [] ]

				#r[g][0].append( eval(score))
				r[g][0].append( eval(pvalue))
				r[g][1].append( goid )




	f.close()

	return r


def __calculateAUC2(p_rank, n_rank, output):


	f = open(output+'.txt','w')
	f.write('False-positive rate(1-specificity)\tTrue-positive rate(sensitivity)\n')

	start_from = 0.0
	start_to = 1.0
	step = 0.01




	pos_false_positive_rate_x = [] # 1- specificity
	pos_true_positive_rate_y = [] # sensitivity


	threshold = start_from 
	while(threshold <= start_to  ):
		threshold = threshold + step


		specificity = 1.0 - threshold

		total_trial = float( len( p_rank  ))
		total_rank = 100.0

		success = 0.0

		for index in range(len(p_rank)):
			p_v = p_rank[index]
			n_v = n_rank[index]

			p_r = float(p_v) / total_rank  
			n_r = float(n_v) / total_rank  

			if p_r <= threshold < n_r:
				success = success + 1.0

		sensitivity = success / total_trial

		pos_false_positive_rate_x.append( 1.0 - specificity  )
		pos_true_positive_rate_y.append( sensitivity )

		f.write(  str(1.0-specificity) + '\t' + str(sensitivity) + '\n'  )


	auc = statistics.average( pos_true_positive_rate_y  )
	f.close()

	return auc

def __calculateAUC(p_rank, output):


	f = open(output+'_pos.txt','w')
	f.write('False-positive rate(1-specificity)\tTrue-positive rate(sensitivity)\n')

	start_from = 0.0
	start_to = 1.0
	step = 0.01




	pos_false_positive_rate_x = [] # 1- specificity
	pos_true_positive_rate_y = [] # sensitivity


	threshold = start_from 
	while(threshold <= start_to  ):
		threshold = threshold + step

		specificity = 1.0 - threshold

		total_trial = float( len( p_rank  ))
		total_rank = 100.0

		success = 0.0
		for v in p_rank:
			r = float(v) / total_rank  # rank ratio
			if r <= threshold:
				success = success + 1

		sensitivity = success / total_trial

		pos_false_positive_rate_x.append( 1.0 - specificity  )
		pos_true_positive_rate_y.append( sensitivity )

		f.write(  str(1.0-specificity) + '\t' + str(sensitivity) + '\n'  )


	auc_pos = statistics.average( pos_true_positive_rate_y  )
	f.close()




	return auc_pos








def copyCache(interaction_file, threshold):

	fname = os.path.basename (interaction_file) + \
	        '_' + str( threshold ).strip() + '.cache.txt'


	temp_file = _0_Preprocess.TEMP_PATH + '/' + fname
	cache_file = './7. Network Propagation/' + fname

	if not os.path.exists(temp_file) and os.path.exists(cache_file):
		import shutil
		shutil.copy(cache_file, temp_file)

		return temp_file


	elif os.path.exists(temp_file):
		return temp_file

	elif os.path.exists(cache_file):
		return cache_file

	else:
		return None





def logmsg(msg):

	'''
	f=open('r:/netpro.log', 'a')
	t = repr(time.ctime())
	f.write(t + '\t' + msg + '\n')
	f.close()
	'''

	return



def init(config_file, opt):

	global RND_ID, TEMP_PATH

	_0_Preprocess.init(config_file)



	RND_ID = opt[_0_Preprocess.OPT_RND_ID]
	TEMP_PATH = _0_Preprocess.TEMP_PATH


def main(opt):


	global MP_PROCESSORS, RND_ID, OPTIONS

	print '''
========================================================    
    [7] Network Propagation
========================================================          
'''
	start_time = datetime.now()	


	OPTIONS = opt

	config_file = opt[_0_Preprocess.OPT_CONFIG_FILE]
	MP_PROCESSORS = opt[_0_Preprocess.OPT_MULTIMP]
	rnd_or_not = opt[_0_Preprocess.OPT_RANDOM_OR_NOT]
	test_or_not = opt[_0_Preprocess.OPT_TEST]
	iteration = opt[_0_Preprocess.OPT_ITERATION]	


	# random network를 이용한 p-value를 사용하지 않음.
	# 단, r=0.0을 reference로 사용함.


	print ('Network Propagation CPU=' + str(MP_PROCESSORS))


	if config_file is None:
		print '!!!!!!!!! No config file'
		return

	init(config_file, opt)

	RND_ID = opt[_0_Preprocess.OPT_RND_ID]
	_0_Preprocess.RND_ID = RND_ID



	if test_or_not == False:

		for i in range(len(_0_Preprocess.PPI_FILES)):

			string_db_file = _0_Preprocess.PPI_FILES[i][0]
			index1 = _0_Preprocess.PPI_FILES[i][1]
			index2 = _0_Preprocess.PPI_FILES[i][2]
			score_index = _0_Preprocess.PPI_FILES[i][3]
			threshold = _0_Preprocess.PPI_FILES[i][4]
			separator = _0_Preprocess.PPI_FILES[i][5]

			print '-----------------------------------'
			print ' [PPI_FILE]', string_db_file
			print '-----------------------------------'

			if MP_PROCESSORS == 1:
				# single processor
				runFinal(_0_Preprocess.PPI_FILES[i], string_db_file, index1, index2, score_index, threshold, separator)
			else:
				# multiple processors
				runFinal_mp(config_file, _0_Preprocess.PPI_FILES[i], string_db_file, index1, index2, score_index, threshold, separator, max_thread=MP_PROCESSORS)




	else:


		r = []



		if MP_PROCESSORS == 1:

			for i in range(len(_0_Preprocess.PPI_FILES)):



				string_db_file = _0_Preprocess.PPI_FILES[i][0]
				index1 = _0_Preprocess.PPI_FILES[i][1]
				index2 = _0_Preprocess.PPI_FILES[i][2]
				score_index = _0_Preprocess.PPI_FILES[i][3]
				threshold = _0_Preprocess.PPI_FILES[i][4]
				separator = _0_Preprocess.PPI_FILES[i][5]

				print '-----------------------------------'
				print ' [PPI_FILE]', string_db_file
				print '-----------------------------------'



				#iteration = 100


				#iteration = HOW_MANY_TEST_MODIFERS


				x = run(config_file, _0_Preprocess.PPI_FILES[i], iteration, string_db_file, index1, index2, score_index, threshold, separator)
				r.append(x)
		else:
			jobs = []

			mod_maps, mods = __getModifiers()
			diseases = mods.keys()


			for i in range(len(_0_Preprocess.PPI_FILES)):

				string_db_file = _0_Preprocess.PPI_FILES[i][0]
				index1 = _0_Preprocess.PPI_FILES[i][1]
				index2 = _0_Preprocess.PPI_FILES[i][2]
				score_index = _0_Preprocess.PPI_FILES[i][3]
				threshold = _0_Preprocess.PPI_FILES[i][4]
				separator = _0_Preprocess.PPI_FILES[i][5]

				print '-----------------------------------'
				print ' [PPI_FILE]', string_db_file
				print '-----------------------------------'

				#iteration = HOW_MANY_TEST_MODIFERS



				for disease in diseases:

					th = MULTIPROCESS.JOB_THREAD()
					th.set_args(run_mp, [config_file, _0_Preprocess.PPI_FILES[i], iteration, disease, string_db_file, index1, index2, score_index, threshold, separator])
					jobs.append(th)
					#x = run(config_file, _0_Preprocess.PPI_FILES[i], iteration, string_db_file, index1, index2, score_index, threshold, separator)
					#r.append(x)

			ret = MULTIPROCESS.runMultiprocesses(jobs, max_cpu=MP_PROCESSORS)
			for th in jobs:
				r.append(  th.getResult()  )


		print '=================='
		for t in r:
			print '[SUMMARY_FILE]', t




	print '''
========================================================    
    [7] Network Propagation (End)
========================================================          
'''

	end_time = datetime.now()
	print ('Elapsed time: {}'.format(end_time - start_time))
	print time.ctime()



if __name__ == '__main__':


	opt = _0_Preprocess.process_args(sys.argv[1:])
	main(opt)



	#global MP_PROCESSORS

	#print '''
#========================================================    
	#[7] Network Propagation
#========================================================          
#'''
	#start_time = datetime.now()

	##MP_PROCESSORS = 5


	#if len(sys.argv) == 2:


		#config_file = sys.argv[1]
		#init(config_file)

		#for i in range(len(_0_Preprocess.PPI_FILES)):



			#string_db_file = _0_Preprocess.PPI_FILES[i][0]
			#index1 = _0_Preprocess.PPI_FILES[i][1]
			#index2 = _0_Preprocess.PPI_FILES[i][2]
			#score_index = _0_Preprocess.PPI_FILES[i][3]
			#threshold = _0_Preprocess.PPI_FILES[i][4]
			#separator = _0_Preprocess.PPI_FILES[i][5]

			#print '-----------------------------------'
			#print ' [PPI_FILE]', string_db_file
			#print '-----------------------------------'

			#if MP_PROCESSORS == 1:
				## single processor
				#runFinal(_0_Preprocess.PPI_FILES[i], string_db_file, index1, index2, score_index, threshold, separator)
			#else:
				## multiple processors
				#runFinal_mp(string_db_file, index1, index2, score_index, threshold, separator, max_thread=MP_PROCESSORS)




	#elif len(sys.argv) == 3:
		#if sys.argv[2].lower() == 'test':

			#config_file = sys.argv[1]
			#init(config_file)

			#r = []



			#if MP_PROCESSORS == 1:

				#for i in range(len(_0_Preprocess.PPI_FILES)):



					#string_db_file = _0_Preprocess.PPI_FILES[i][0]
					#index1 = _0_Preprocess.PPI_FILES[i][1]
					#index2 = _0_Preprocess.PPI_FILES[i][2]
					#score_index = _0_Preprocess.PPI_FILES[i][3]
					#threshold = _0_Preprocess.PPI_FILES[i][4]
					#separator = _0_Preprocess.PPI_FILES[i][5]

					#print '-----------------------------------'
					#print ' [PPI_FILE]', string_db_file
					#print '-----------------------------------'



					#iteration = 100


					#iteration = HOW_MANY_TEST_MODIFERS


					#x = run(config_file, _0_Preprocess.PPI_FILES[i], iteration, string_db_file, index1, index2, score_index, threshold, separator)
					#r.append(x)
			#else:
				#jobs = []

				#mod_maps, mods = __getModifiers()
				#diseases = mods.keys()


				#for i in range(len(_0_Preprocess.PPI_FILES)):

					#string_db_file = _0_Preprocess.PPI_FILES[i][0]
					#index1 = _0_Preprocess.PPI_FILES[i][1]
					#index2 = _0_Preprocess.PPI_FILES[i][2]
					#score_index = _0_Preprocess.PPI_FILES[i][3]
					#threshold = _0_Preprocess.PPI_FILES[i][4]
					#separator = _0_Preprocess.PPI_FILES[i][5]

					#print '-----------------------------------'
					#print ' [PPI_FILE]', string_db_file
					#print '-----------------------------------'

					#iteration = HOW_MANY_TEST_MODIFERS



					#for disease in diseases:

						#th = MULTIPROCESS.JOB_THREAD()
						#th.set_args(run_mp, [config_file, _0_Preprocess.PPI_FILES[i], iteration, disease, string_db_file, index1, index2, score_index, threshold, separator])
						#jobs.append(th)
						##x = run(config_file, _0_Preprocess.PPI_FILES[i], iteration, string_db_file, index1, index2, score_index, threshold, separator)
						##r.append(x)

				#ret = MULTIPROCESS.runMultiprocesses(jobs, max_cpu=MP_PROCESSORS)
				#for th in jobs:
					#r.append(  th.getResult()  )


			#print '=================='
			#for t in r:
				#print '[SUMMARY_FILE]', t




	#print '''
#========================================================    
	#[7] Network Propagation (End)
#========================================================          
#'''

	#end_time = datetime.now()
	#print ('Elapsed time: {}'.format(end_time - start_time))
	#print time.ctime()
