# -*- coding: ms949 -*-
import datetime
import os
import _0_Preprocess
import copy
import random
import threading
import MULTIPROCESS
import MyUtil
from datetime import datetime
import time
import sys
import pickle
from PPI import PPIDB_STRING
import statistics

RND_ID = _0_Preprocess.RND_ID
TEMP_PATH = _0_Preprocess.TEMP_PATH
OUTPUT_FOLDER = _0_Preprocess.OUTPUT_FOLDER

Pnp_matrix = None
Pnp_dict = {}
DIRECTED_PPI = False

# =================================
MP_PROCESSORS = 1
REMOVE_TEMP_FILES = True

P_VALUE = _0_Preprocess.NETWORK_EXPANSION_PVALUE
# PRELOAD_NO_OF_RANDOM_NETWORKS = 100
ITERATION = 100

HOW_MANY_TEST_MODIFERS = 100
NETWORK_TEMP_FILES = []

# 여기가 중요함.
DEGREE = 1 # 1로 하니까Akt가 자주 나옴...
#DEGREE = 1  #
# 0 - 그냥 진짜로 random
# 1 - degree가 +-1 범위내의 것만 골라서
# 3 - reliability 차이가 1% 이내인 것들만 사용
# =================================


JOBS_PER_CPU = 2
# __train_k에서 사용함.

PYPY = False

P_VALUE_CUTOFF = 1e-50

'''
def mean(data):
	"""Return the sample arithmetic mean of data."""
	n = len(data)
	if n < 1:
		raise ValueError('mean requires at least one data point')
	return sum(data)/float(n) # in Python 2 use sum(data)/float(n)

def _ss(data):
	"""Return sum of square deviations of sequence data."""
	c = mean(data)
	ss = sum((x-c)**2 for x in data)
	return ss

def pstdev(data):
	"""Calculates the population standard deviation."""
	n = len(data)
	if n < 2:
		raise ValueError('variance requires at least two data points')
	ss = _ss(data)
	pvar = ss/float(n) # the population variance
	return pvar**0.5


def getFileListIn(path):

	files = []
	folders = []

	temp = os.listdir(path)
	for a in temp:
		try:
			if os.path.isfile(path + '/' + a):
				files.append(a)
			elif os.path.isdir(path + '/' + a ):
				folders.append(a)

		except:
			print "Can't read ", a

	return files, folders



#---------------------------------------------
def executeDosCommand2(cmd, printout=True):
	f = subprocess.Popen(cmd, shell = True, stdout = subprocess.PIPE).stdout
	re = ''
	line = ''
	while True:
		last_line = line
		line = f.readline()
		if printout:
			print "[EXEC] " + line.replace('\n','')
		re += line
		if not line: break
	return re



def calculate_pvalues(values):

	new_scores = {}


	temp_file = _0_Preprocess.TEMP_PATH + '/' + getRandomString(10) + '.txt'
	temp_file2 = _0_Preprocess.TEMP_PATH + '/' + getRandomString(10) + '.txt'

	f=open(temp_file,'w')
	f.write('value\tmean\tstdev\n')
	for g in values.keys():
		val, avg, std = values[g]


		txt = g + '\t' + str(val) + '\t' + str(avg) + '\t' + str(std)
		f.write(txt+'\n')

	f.close()

	cmd = 'python _X_norm_test.py "' + temp_file + '" "' + temp_file2 + '"'
	executeDosCommand2(cmd)

	f=open(temp_file2,'r')
	for s in f.readlines():
		s = s.strip()
		if len(s) == 0: continue

		t = s.split('\t')
		gene = t[0].strip()
		score = eval(t[1])

		new_scores[gene] = score

	return new_scores

def getRandomString(length):
	str1 = 'abcdefghijklmnopqrstuvwxyz'
	str2 = str1.upper()
	numbers = '0123456789'

	box = str1+str2 + numbers

	r = ''

	for i in range(length):
		s = random.sample(box, 1) [0]
		r += s
	return r

'''


# class JOB(multiprocessing.Process):
class JOB(threading.Thread):
	values = None
	over = False
	result = None
	
	def __init__(self):
		self.values = None
		self.over = False
		self.result = None
		threading.Thread.__init__(self)
	
	# multiprocessing.Process.__init__(self)
	
	def set(self, values):
		self.values = values
	
	def run(self):
		self.over = False
		
		disease = self.values[0]
		mod_pos = self.values[1]
		mod_neg = self.values[2]
		total_genes = self.values[3]
		total_genes_without_modifiers = self.values[4]
		iteration_th = self.values[5]
		ppi_box = self.values[6]
		
		self.result = mainJob(disease, mod_pos, mod_neg, total_genes, total_genes_without_modifiers, iteration_th,
		                      ppi_box)
		
		self.over = True
	
	def getResult(self):
		return self.result


def chunks(seq, size):
	# list를 쪼갠다. 이때 n개가 들어가도록 쪼갠다. 리스트의 갯수가 아니라 item의 갯수를 정해준다.
	return [seq[i:i + size] for i in xrange(0, len(seq), size)]


def __randomSampling(total_genes, cross_set, mod_box, basket, ppi):
	new_pos_mod = []
	
	global DEGREE  # default
	
	if DEGREE == 0:
		
		# 완전 랜덤하게 고른다.
		new_pos_mod = random.sample(total_genes, len(mod_box))
	
	elif DEGREE == 1:
		# preserving interaction degrees
		basket = getBasket(total_genes, ppi)
		
		new_pos_mod = []
		
		for g in mod_box:
			
			deg = len(ppi.getPartners(g))
			deg_str1 = str(deg - 1).strip()
			deg_str2 = str(deg).strip()
			deg_str3 = str(deg + 1).strip()
			
			box = []
			if deg_str1 in basket:
				box = box + basket[deg_str1]
			if deg_str2 in basket:
				box = box + basket[deg_str2]
			if deg_str3 in basket:
				box = box + basket[deg_str3]
			
			random.shuffle(box)
			
			for v in box:
				if not v in new_pos_mod:
					new_pos_mod.append(v)
					break
		
		if len(mod_box) != len(new_pos_mod):
			print('Random sampling error: ', len(mod_box), len(new_pos_mod))
	
	
	
	
	
	elif DEGREE == 3:
		
		# reliability
		# test 및 cross val set에는 포함되지 않는 걸로 골라야하는데, 여기선 고려 안해도 되니 [] 빈 리스트로 넘김
		new_pos_mod = __randomSamplingBasedOnReliability(total_genes, [], mod_box, ppi, [])
		
		if len(mod_box) != len(new_pos_mod):
			print('Random sampling error: ', len(mod_box), len(new_pos_mod))
	
	return new_pos_mod


def __getOutputFilename(fname, threshold):
	while (True):
		rnd = str(random.randint(0, 1000000000)).strip()
		
		ppi_filename_only = os.path.basename(fname)
		ppi_filename_only = os.path.splitext(ppi_filename_only)[0]
		ppi_filename_only = './6. Functional Similarity/network/' + ppi_filename_only + '_' + str(threshold).strip()
		output = ppi_filename_only + '_' + rnd + '.net'
		if not os.path.exists(output):
			return output
	
	return None


def __getRandomNetwork(interaction_file, how_many):
	# random network .
	
	# interaction_file = _0_Preprocess.STRING_FILE
	index1 = 0
	index2 = 1
	score_index = 2
	threshold = -1  # 이미 cutoff로 잘라서 저장해 놓은 것이라서 상관없음
	
	organism = interaction_file[:  interaction_file.find('.')]
	path = './ppi/random/' + organism + '/-1/'
	
	print(path)
	
	files, folders = MyUtil.getFileListIn(path)
	prebuilt = []
	for f in files:
		if '.txt' in f and not 'cache' in f and not 'pickle' in f:
			prebuilt.append(path + f)
	# print path+f
	
	print('Random network # =', len(prebuilt))
	fnames = random.sample(prebuilt, how_many)
	
	# ppi = PPI.PPIDB_STRING()
	# ppi.load(fname, 0, 1, 2, -1)
	
	return fnames

'''
def loadInteraction_bak(interaction_file, iteration, index1, index2, index3, score_index, threshold):
	global P_VALUE
	
	# global NETWORK_TEMP_FILES
	
	print('Loading interactions = ', interaction_file)
	box = []
	cache_file_box = []
	
	ppi = PPIDB_STRING()
	ppi.load(interaction_file, index1, index2, score_index, threshold, ' ', _0_Preprocess.STRING_ALIAS_FILE)
	box.append(ppi)  # normal network
	
	cache_file = interaction_file + '_' + str(threshold) + '_cache.pickle'
	cache_file_box.append(cache_file)
	
	ppi_filename_only = os.path.basename(interaction_file)
	ppi_filename_only = os.path.splitext(ppi_filename_only)[0] + '_' + str(threshold).strip()
	# ppi_filename_only = './6. Functional Similarity/network/'+ppi_filename_only
	
	if P_VALUE:
		
		# random networks
		random_ppi_path = './ppi/random/fly/-1/'
		files, folders = MyUtil.getFileListIn(random_ppi_path)
		prebuilt = []
		for f in files:
			fname = random_ppi_path + '/' + f
			if f.find('.pickle') < 0:
				prebuilt.append(fname)
		
		if len(prebuilt) > iteration:
			prebuilt = random.sample(prebuilt, iteration)
		
		for i in range(iteration):
			
			print('\t', i + 1, '/', iteration, end='\r')
			p = PPI.PPIDB_STRING()
			
			if i < len(prebuilt):
				p.load(prebuilt[i], index1, index2, score_index, -1)
				box.append(p)
				
				cf = prebuilt[i] + '_' + str(threshold) + '_cache.pickle'
				cache_file_box.append(cf)
			
			# box.append(    [ prebuilt[i] , index1, index2, score_index, threshold   ] )
			
			else:
				
				# random network가 모자라네.
				# 더 만든다.
				output = __getOutputFilename(interaction_file, threshold)
				
				p.interDB = copy.deepcopy(ppi.interDB)
				p.scoreDB = copy.deepcopy(ppi.scoreDB)
				
				p.randomizeNetwork()
				p.saveInteractionIntoFile(output)
				
				box.append(p)
				cf = output + '_' + str(threshold) + '_cache.pickle'
				cache_file_box.append(cf)
	
	return box, cache_file_box
'''

def loadInteraction(interaction_file, index1, index2, index3, score_index, threshold):
	global P_VALUE
	
	# global NETWORK_TEMP_FILES
	
	print('Loading interactions = ', interaction_file)
	box = []
	cache_file_box = []
	
	ppi = PPIDB_STRING()
	ppi.load(interaction_file, index1, index2, score_index, threshold, ' ', _0_Preprocess.STRING_ALIAS_FILE)
	# box.append(ppi)  # normal network
	
	cache_file = interaction_file + '_' + str(threshold) + '_cache.pickle'
	x = {'ppi_file': interaction_file, 'index1': index1, 'index2': index2,
	     'score_index': score_index, 'threshold': threshold, 'alias_file': _0_Preprocess.STRING_ALIAS_FILE,
	     'cache_file': cache_file}
	box.append(x)
	
	ppi_filename_only = os.path.basename(interaction_file)
	ppi_filename_only = os.path.splitext(ppi_filename_only)[0] + '_' + str(threshold).strip()
	# ppi_filename_only = './6. Functional Similarity/network/'+ppi_filename_only
	
	if P_VALUE:
		
		# random networks
		random_ppi_path = './ppi/random/fly/-1/'
		files, folders = MyUtil.getFileListIn(random_ppi_path)
		
		for f in files:
			fname = random_ppi_path + '/' + f
			if f.find('.pickle') < 0:
				#print('Random ppi = ', f)
				cf = fname + '_' + str(threshold) + '_cache.pickle'
				x = {'ppi_file': fname, 'index1': index1, 'index2': index2, 'score_index': score_index,
				     'threshold': threshold, 'alias_file': _0_Preprocess.STRING_ALIAS_FILE,
				     'cache_file': cf}
				box.append(x)
		
		'''
		if len(box) < ITERATION * 1.2:
			
			print('making more random networks...')
			
			for i in range(int(ITERATION * 1.2 - len(box))):
				# random network가 모자라네.
				# 더 만든다.
				output = __getOutputFilename(interaction_file, threshold)
				
				p = PPIDB_STRING()
				p.interDB = copy.deepcopy(ppi.interDB)
				p.scoreDB = copy.deepcopy(ppi.scoreDB)
				
				p.randomizeNetwork()
				p.saveInteractionIntoFile(output)
				# box.append(p)
				
				cf = output + '_' + str(threshold) + '_cache.pickle'
				x = {'ppi_file': output, 'index1': index1, 'index2': index2, 'score_index': score_index,
				     'threshold': threshold, 'alias_file': _0_Preprocess.STRING_ALIAS_FILE,
				     'cache_file': cf}
				box.append(x)
		'''
		
	return box


def run(iteration, interaction_file, index1, index2, index3, score_index, threshold):
	global RND_ID, NETWORK_TEMP_FILES, DEGREE, DISEASE_SET, MP_PROCESSORS
	global P_VALUE, ITERATION
	
	print('ID=', RND_ID)
	
	mods = __getModifiers()
	
	summary_file = TEMP_PATH + '/NetExp' + RND_ID + '_' + '_summary.txt'
	f = open(summary_file, 'w')
	dd = list(mods)
	dd.sort()
	
	NETWORK_TEMP_FILES = [interaction_file, index1, index2, score_index, threshold]
	
	info = '--------------------------------------------' + '\n' + \
	       'THREAD = ' + str(MP_PROCESSORS) + '\n' + \
	       'Remove temporary files = ' + str(REMOVE_TEMP_FILES) + '\n' + \
	       'Calculate p-value = ' + str(P_VALUE) + '\n' + \
	       '\t# of random networks = ' + str(ITERATION) + '\n' + \
	       'Iteration = ' + str(HOW_MANY_TEST_MODIFERS) + '\n' + \
	       'Randomization method = ' + str(DEGREE) + '\n' + \
	       'r = ' + str(R_PARAMETER)
	print(info)
	
	ppi_box = loadInteraction(interaction_file, index1, index2, index3, score_index, threshold)  # 100개 네트워크를 만든다.
	total_genes = __getTotalGeneList(ppi_box)
	
	print('init...')
	global Pnp_matrix, Pnp_dict
	total_gene_names = ppi_box[0].getWholeGeneList()
	total_gene_names.sort()
	local_Pnp_matrix, local_Pnp_dict = __buildPnpMatrix(total_gene_names, ppi_box[0])
	print('done')
	
	ppi_box2 = copy.deepcopy(ppi_box)
	
	auc_output = TEMP_PATH + '/NetExp_' + RND_ID + '_' + '_'.join(dd) + '_auc.txt'
	
	mod_box = {}
	
	for key in dd:
		mod_pos = mods[key]
		new_mod_pos = __removeUnannotatedModifiers(total_genes, mod_pos)
		mod_box[key] = new_mod_pos
	
	total_genes_without_modifiers = __getTotalGenesWithoutModifiers(total_genes, mod_box)
	
	s = 'Interaction = ' + interaction_file + '\n' + \
	    'Disease = ' + repr(DISEASE_SET) + '\n' + \
	    'Total annotated gene # = ' + str(len(total_genes)) + '\n'
	
	f.write(info + '\n')
	f.write('-------------------------\n' + s + '\n')
	
	print(s)
	
	prank = []
	
	# 정해진 만큼 train과 test를 반복한다.
	
	random_or_not = 1  # 1-random. 0-not random  This must be RANDOM to calculate AUC
	
	for i in range(iteration):
		Pnp_dict = copy.deepcopy(local_Pnp_dict)
		Pnp_matrix = copy.deepcopy(local_Pnp_matrix)
		
		print('\r\t', RND_ID, ' Disease = ', repr(dd), ' # = ', i + 1)
		
		pos_rank = mainJob(dd, copy.deepcopy(mod_box), total_genes, total_genes_without_modifiers, i, ppi_box2,
		                   random_or_not)
		
		print('\t', pos_rank)
		
		prank.append(pos_rank)
	
	print('')
	
	foname = TEMP_PATH + '/' + RND_ID + '_' + '_'.join(DISEASE_SET) + '_ranks.txt'
	__saveRank(prank, foname)
	
	pos_auc = __calculateAUC(prank, auc_output)
	
	s = 'Avg rank = ' + str(statistics.average(prank)) + ' +- ' + str(statistics.stdev(prank)) + '\n' + \
	    'AUC + = ' + str(pos_auc) + '\n'
	
	print(s)
	f.write(s + '\n')
	f.close()


def __saveRank(prank, output):
	f = open(output, 'w')
	f.write('Pos rank\n')
	
	for index in range(len(prank)):
		s = str(prank[index])
		f.write(s + '\n')
	f.close()


def __int2strList(v):
	r = []
	for g in v:
		r.append(str(g).strip())
	return r


def __countOverlaps(a, b):
	cnt = 0
	for v in a:
		if v in b:
			cnt = cnt + 1
	return cnt


def __extractCommonModifiers(total_genes, mod_box):
	r = []
	
	print('Common modifiers = ')
	for g in total_genes:
		cnt = 0
		for key in mod_box:
			mp = mod_box[key]
			if not g in mp:
				cnt = 1
				break
		
		if cnt == 0:
			# print g
			r.append(g)
			
			for key in mod_box:
				mp = mod_box[key]
				
				# print 'deleting ', g, ', ', len(mp), ' -> ',
				
				del mp[g]
		# print mod_box[key][0].has_key(g), len(mp)
	return r


def mainJob(disease_set, mod_box, total_genes, total_genes_without_modifiers, iteration_th, ppi_box, random_or_not):
	global RND_ID, TEMP_PATH, P_VALUE, DEGREE
	
	ppi_box2 = []
	ppi_box2 = ppi_box
	
	cross_set = __extractCommonModifiers(total_genes, mod_box)
	
	test_pos, test_others = __sampleModifiers(total_genes_without_modifiers, cross_set, mod_box)
	
	new_mod_box = None
	
	if random_or_not == 1:
		
		# random network
		# degree
		basket = getBasket(total_genes, ppi_box[0])
		
		cross_set, new_mod_box = __randomSampling(total_genes, cross_set, mod_box, basket, ppi_box[0])
		
		test_pos, test_others = __sampleModifiers(total_genes, cross_set, new_mod_box)
	
	
	
	
	
	
	else:
		ppi_box2 = ppi_box
		new_mod_box = copy.deepcopy(mod_box)
	
	print('TEST gene = ', test_pos, ' out of ', len(cross_set))
	train_output = TEMP_PATH + '/NetExp' + RND_ID + '__trained_' + str(iteration_th).strip() + '_'.join(
		disease_set) + '_' + str(iteration_th).strip() + '.txt'
	
	test_set_all = [test_pos] + test_others
	
	__train(total_genes, new_mod_box, train_output, test_set_all, ppi_box2)
	
	train_output2 = train_output + '_gene.txt'
	r1 = __test(train_output2, test_pos, test_others)
	
	'''
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

    '''
	return r1


def getBasket(total_genes, ppi):
	box = {}
	
	for g in total_genes:
		# deg = len( ppi.getIncomingPartnersOf(g)  )
		deg = len(ppi.getPartners(g))
		deg_str = str(deg).strip()
		
		if not deg_str in box:
			box[deg_str] = []
		
		if not g in box[deg_str]:
			box[deg_str].append(g)
	
	return box


def __removeUnannotatedModifiers(total_genes, mods):
	r = {}
	for k in mods:
		if k in total_genes and not k in r:
			r[k] = mods[k]
	return r


def __getTotalGeneList(ppi):
	# ppi = ppi_box[0]
	
	total_gene_list = []
	for g in ppi.getWholeGeneList():
		total_gene_list.append(g)
	return total_gene_list


def __sampleModifiers(total_genes, cross_set, mod_box):
	print('common set # ', len(cross_set))
	
	other_genes = []
	test_pos = random.sample(cross_set, 1)[0]
	
	while (len(other_genes) < 99):
		o = random.sample(total_genes, 1)[0]
		has = False
		for d in mod_box:
			mp = mod_box[d]
			if o in mp:
				has = True
				break
		
		if has == False and o != test_pos:
			other_genes.append(o)
	
	return test_pos, other_genes


def __getTotalGenesWithoutModifiers(total_genes, mod_box):
	r = []
	for t in total_genes:
		
		pass_or_fail = False
		
		for k in mod_box:
			if t in mod_box[k]:
				pass_or_fail = True
				break
		
		if pass_or_fail == False:
			r.append(t)
	
	return r


# =============================================


def __train(total_genes, mod_box, output, test_set_all, ppi_box):
	#
	
	pvalues = __calculatePvalue(total_genes, test_set_all, mod_box, output, ppi_box)
	return __resave(total_genes, mod_box, pvalues, output, test_set_all, ppi_box)


def __resave(total_genes, mod_box, pvalues, output, test_set_all, ppi_box):
	global P_VALUE
	
	ofile = output + '_gene.txt'
	
	box = __rebox(test_set_all, mod_box, pvalues, output)
	ofile2 = __save2(box, mod_box, ofile, ppi_box)  # pvalue로 정렬, pvalue=False라면 __save3랑 결과가 동일할 것임.
	ofile3 = __save3(box, mod_box, ofile, ppi_box)  # real score (avg)로 정렬
	ofile4 = None
	
	if P_VALUE:
		ofile4 = __save4(box, mod_box, ofile, ppi_box)  # pvalue cutoff이하인 것은 동일하게 취급하고, 이 경우만 real score로 정렬함.
	
	# __save2 + __save3 공통 top100만 추려낸다.
	top = 100
	
	cfile = __saveCommonGenes(
		output + '_top100_both_in_zscore_and_realscore_' + MyUtil.getRandomString(15) + '.txt',
		ofile2, ofile3, top)
	print('Saved modifiers identified in top100 (both pvalue/realscores): ', cfile)
	
	# ofile = ofile2 이네...
	return [ofile2, ofile3, ofile4, cfile]


def __saveCommonGenes(output, ofile1, ofile2, top):
	title1, list1 = getModifiersRank2(ofile1, top)
	title2, list2 = getModifiersRank2(ofile2, top)
	
	f = open(output, 'w')
	
	f.write('rank_zscore\trank_real_score\t' + title1 + '\n')
	
	for g in list1:
		if g in list2:
			r1, text1 = list1[g]
			r2, text2 = list2[g]
			f.write(str(r1) + '\t' + str(r2) + '\t' + text1 + '\n')
	
	f.close()
	
	return output


def getModifiersRank2(fname, threshold_rank):
	f = open(fname, 'r')
	
	rank = 0.0
	
	r = {}
	
	init = True
	title = ''
	
	print('Loading a modifier file', fname, '....')
	
	for s in f.readlines():
		
		if init:
			
			init = False
			title = s.strip()
		
		else:
			s = s.replace('m_', '').replace('x_', '')
			
			x = s.split('\t')
			
			gene_id = x[0]
			# pvalue = eval(x[1])  # meaningless
			# real_score = eval(x[2]) # meaningless
			text = s.strip()
			
			rank = rank + 1.0
			r[gene_id] = [rank, text]
			
			if rank >= threshold_rank:
				break
	
	f.close()
	
	return title, r


def __save4(box, mod_box, output, ppi_box):
	global P_VALUE, P_VALUE_CUTOFF
	
	output = output + '____pvalue_cutoff_and_real_score_based_' + MyUtil.getRandomString(15) + '.txt'
	
	ppi = PPIDB_STRING()
	ppi = ppi_box[0]
	
	order = False
	
	# p-value를 기준으로 cutoff 이내인 것과, 이 외인 것으로 구분한다.
	
	# 단순히 score만 가지고 정렬한 후 rank를 저장한다.
	
	print(' ----------- Saving pvalue cutt + real score...: ', P_VALUE_CUTOFF)
	
	if P_VALUE:
		order = False
		#
		box3 = []
		for g in box:
			box3.append([g, box[g]])
		
		# box[g] = [p-value, score, desc]
		
		# sort by z-score
		for j in range(len(box3) - 1):
			for i in range(len(box3) - 1):
				if box3[i][1][0] < box3[i + 1][1][0]:  # if p-value is -1, it must have a low rank
					box3[i], box3[i + 1] = box3[i + 1], box3[i]
		
		'''
		# sort by real score if pvalue < cutoff
		for j in range(len(box3) - 1):
			for i in range(len(box3) - 1):
				
				# if box3[i][1][0] == box3[i + 1][1][0]:  # p-value
				if box3[i][1][0] <= P_VALUE_CUTOFF and box3[i + 1][1][0] <= P_VALUE_CUTOFF:  # p-value
					
					# below pvalue cutoff
					if box3[i][1][1] < box3[i + 1][1][1]:  # real score
						box3[i], box3[i + 1] = box3[i + 1], box3[i]
				
				else:
					break
		'''
		
		# elif box3[i][1][0] > box3[i + 1][1][0] or \
		#				box3[i][1][0] == -1:  # if p-value is -1, it must have a low rank
		#	box3[i], box3[i + 1] = box3[i + 1], box3[i]
	
	# else:
	#	# 사실 여긴 실행될 일이 없음.
	#	box3 = sorted(box.items(), key=lambda (k, v): v[0], reverse=order)
	
	maps, mods = _0_Preprocess.getModifiers()
	diseases = list(mods)
	diseases.sort()
	
	good_box = []
	bad_box = []
	
	f = open(output, 'w')
	f.write('Gene\tpvalue\treal_score\tIntx Degree\tName\tDesc\t' + '\t'.join(diseases) + '\n')
	
	all_avg_deg = []
	str_avg_deg = []
	str_proteins = 0
	non_str_proteins = 0
	
	for k, v in box3:
		
		b = []
		for d in diseases:
			if k in mods[d]:
				b.append('o')
			else:
				b.append('_')
		
		score = v[0]  # pvalue
		real_score = v[1]  # netexp score
		desc = v[-1]
		deg = 0
		
		try:
			deg = len(ppi.getPartners(k))
			
			str_proteins += 1
			str_avg_deg.append(deg)
		except:
			non_str_proteins += 1
		
		all_avg_deg.append(deg)
		
		s = k + '\t' + str(score) + \
		    '\t' + str(real_score) + \
		    '\t' + str(deg) + \
		    '\t' + _0_Preprocess.ID_BOX.getCommonNameOf(k) + \
		    '\t' + desc + \
		    '\t' + '\t'.join(b)
		
		# if box[k][0] <= P_VALUE_CUTOFF:
		#	good_box.append(s)
		# else:
		#	bad_box.append(s)
		
		f.write(s + '\n')
	
	# f.write('\n'.join(good_box) + '\n')
	# f.write('\n'.join(bad_box) + '\n')
	
	# f.write('---------- Additional information ----------------\n')
	# f.write(' STRING proteins = ' + str(str_proteins) + '\n')
	# f.write('        avg edges = ' + str(MyUtil_pypy.mean(str_avg_deg)) + '+-' + str(
	#	MyUtil_pypy.pstdev(str_avg_deg)) + '\n')
	# f.write(' Non-STRING proteins = ' + str(non_str_proteins) + '\n')
	
	f.close()
	
	return output


def __save2(box, mod_box, output, ppi_box):
	global P_VALUE
	
	ppi = PPIDB_STRING()
	ppi = ppi_box[0]
	
	order = True
	
	if P_VALUE:
		order = False
		#
		box3 = []
		for g in box:
			box3.append([g, box[g]])
		# box[g] = [p-value, score, desc]
		
		# sort
		for j in range(len(box3) - 1):
			for i in range(len(box3) - 1):
				
				if box3[i][1][0] == box3[i + 1][1][0]:  # p-value
					# same p-value
					if box3[i][1][1] < box3[i + 1][1][1]:  # real score
						box3[i], box3[i + 1] = box3[i + 1], box3[i]
				#elif box3[i][1][0] > box3[i + 1][1][0] or \
				#		box3[i][1][0] == -1:  # if p-value is -1, it must have a low rank
				elif box3[i][1][0] < box3[i + 1][1][0]:  # if p-value is -1, it must have a low rank
					box3[i], box3[i + 1] = box3[i + 1], box3[i]
	
	else:
		box3 = sorted(box.items(), key=lambda item: item[1][0], reverse=order)
	
	maps, mods = _0_Preprocess.getModifiers()
	diseases = list(mods)
	diseases.sort()
	
	f = open(output, 'w')
	f.write('Gene\tpvalue\treal_score\tIntx Degree\tName\tDesc\t' + '\t'.join(diseases) + '\n')
	
	all_avg_deg = []
	str_avg_deg = []
	str_proteins = 0
	non_str_proteins = 0
	
	for k, v in box3:
		
		b = []
		for d in diseases:
			if k in mods[d]:
				b.append('o')
			else:
				b.append('_')
		
		score = v[0]
		real_score = v[1]
		desc = v[-1]
		deg = 0
		try:
			deg = len(ppi.getPartners(k))
			
			str_proteins += 1
			str_avg_deg.append(deg)
		except:
			non_str_proteins += 1
		
		all_avg_deg.append(deg)
		
		s = k + '\t' + str(score) + \
		    '\t' + str(real_score) + \
		    '\t' + str(deg) + \
		    '\t' + _0_Preprocess.ID_BOX.getCommonNameOf(k) + \
		    '\t' + desc + \
		    '\t' + '\t'.join(b)
		f.write(s + '\n')
	
	# f.write('---------- Additional information ----------------\n')
	# f.write(' STRING proteins = ' + str(str_proteins) + '\n')
	# f.write('        avg edges = ' + str(MyUtil_pypy.mean(str_avg_deg)) + '+-' + str( MyUtil_pypy.pstdev(str_avg_deg)) + '\n')
	# f.write(' Non-STRING proteins = ' + str(non_str_proteins) + '\n')
	
	f.close()
	
	return output


def __save3(box, mod_box, output, ppi_box):
	global P_VALUE
	
	output = output + '_real_score_based_' + MyUtil.getRandomString(15) + '.txt'
	
	ppi = PPIDB_STRING()
	ppi = ppi_box[0]
	
	order = True
	
	# real score를 높은 순서대로 정렬한다.
	box3 = sorted(box.items(), key=lambda item: item[1][1], reverse=order)
	
	'''
	if P_VALUE:
		order = False
		#
		box3 = []
		for g in box.():
			box3.append([g, box[g]])
			# box[g] = [p-value, score, desc]

		# sort
		for j in range(len(box3)-1):
			for i in range(len(box3)-1):
				if box3[i][1][0] == box3[i+1][1][0]:  # p-value
					# same p-value
					if box3[i][1][1] < box3[i+1][1][1]: # real score
						box3[i], box3[i+1] = box3[i+1], box3[i]
				elif box3[i][1][0] > box3[i+1][1][0] or \
				     box3[i][1][0] == -1: # if p-value is -1, it must have a low rank
					box3[i], box3[i+1] = box3[i+1], box3[i]

	else:
		# real score를 높은 순서대로 정렬한다.
		box3 = sorted( box.items(), key=lambda (k,v): v[1], reverse = order)
	'''
	
	maps, mods = _0_Preprocess.getModifiers()
	diseases = list(mods)
	diseases.sort()
	
	f = open(output, 'w')
	f.write('Gene\tpvalue\treal_score\tIntx Degree\tName\tDesc\t' + '\t'.join(diseases) + '\n')
	
	all_avg_deg = []
	str_avg_deg = []
	str_proteins = 0
	non_str_proteins = 0
	
	for k, v in box3:
		
		b = []
		for d in diseases:
			if k in mods[d]:
				b.append('o')
			else:
				b.append('_')
		
		score = v[0]
		real_score = v[1]
		desc = v[-1]
		deg = 0
		try:
			deg = len(ppi.getPartners(k))
			
			str_proteins += 1
			str_avg_deg.append(deg)
		except:
			non_str_proteins += 1
		
		all_avg_deg.append(deg)
		
		s = k + '\t' + str(score) + \
		    '\t' + str(real_score) + \
		    '\t' + str(deg) + \
		    '\t' + _0_Preprocess.ID_BOX.getCommonNameOf(k) + \
		    '\t' + desc + \
		    '\t' + '\t'.join(b)
		f.write(s + '\n')
	
	# f.write('---------- Additional information ----------------\n')
	# f.write(' STRING proteins = ' + str(str_proteins) + '\n')
	# f.write('        avg edges = ' + str(MyUtil_pypy.mean(str_avg_deg)) + '+-' + str( MyUtil_pypy.pstdev(str_avg_deg)) + '\n')
	# f.write(' Non-STRING proteins = ' + str(non_str_proteins) + '\n')
	
	f.close()
	
	return output


def __addTag(gene, mod_pos):
	tag = ''
	if gene in mod_pos:
		tag = 'm_'
	
	return tag


def __rebox(total_genes, mod_box, pvalues, output):
	global P_VALUE
	
	# calculate score !
	#
	
	# ==================================================================================================
	box = {}
	
	desc = []
	for g in total_genes:
		
		if P_VALUE:
			desc = []
			
			if g in pvalues:
				p_pos, score, txt = pvalues[g]
				
				desc.append(_0_Preprocess.ID_BOX.getCommonNameOf(g))
				desc.append('[' + str(score) + ' ] ')
				desc.append(txt)
				
				box[g] = [p_pos, score, ','.join(desc)]
			else:
				box[g] = [0.0, 0.0, _0_Preprocess.ID_BOX.getCommonNameOf(g)]
		
		else:
			desc = []
			
			p_pos = pvalues[g]
			
			desc.append(_0_Preprocess.ID_BOX.getCommonNameOf(g))
			desc.append('[' + str(p_pos) + ' ] ')
			
			box[g] = [p_pos, ','.join(desc)]
	
	return box


def __save(go, go_box, mod_pos, mod_neg, pvalues, output):
	f = open(output, 'w')
	
	for index in range(2):
		
		m_v = pvalues[index]
		
		s = 'GO_ID\tp-value\t#gene\t#modifier\tGO Desc\tList'
		f.write(s + '\n')
		
		# vx = sorted( go_box.items(), key=lambda (k,v): m_v[k], reverse=False)
		
		for gid in go_box:
			genes = go_box[gid]
			
			# Desc
			go_t = go.getGOTermOf(gid)
			name = go_t.getGOName()
			
			lst = ','.join(__makeGeneList(genes, mod_pos, mod_neg))
			
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
			r.append('m_' + g)
		else:
			r.append(g)
	return r


def __getCommonValuesInTwoLists(l1, l2):
	r = []
	for l in l1:
		if l in l2:
			r.append(l)
	return r


def __calculate1(g1, g2, intersection, ppi):
	r = 0.0
	for g in intersection:
		score1 = float(ppi.getScore(g1, g))
		score2 = float(ppi.getScore(g2, g))
		r = r + score1 * score2
	
	r = r * 2.0
	return r


def __calculate2(g1, intersection, ppi):
	r = 0.0
	for g in intersection:
		score = float(ppi.getScore(g1, g))  # / 100.0
		r = r + score
	return r


def __calculate3(g1, g2, intersection, ppi):
	r = 0.0
	for g in intersection:
		score1 = float(ppi.getScore(g1, g))  # / 100.0
		score2 = float(1.0 - ppi.getScore(g2, g))  # / 100.0
		r = r + score1 * score2
	
	r = r * 2.0
	return r


def __diffList(l1, l2):
	# L1 - L2
	r = []
	for g in l1:
		if not g in l2:
			r.append(g)
	return r


def __calculateLambda(u, v, u_partners, v_partners, intersection):
	n_avg = float((len(u_partners) - 1 + len(v_partners) - 1)) / 2.0
	v = n_avg - float(len(__diffList(u_partners, v_partners)) + len(intersection))
	return max(0.0, float(v))


def __calculateFunctionalSimilarity(u, v, confidence, ppi):
	u_partners = ppi.getPartners(u)  # + [ u ]
	v_partners = ppi.getPartners(v)  # + [ v ]
	
	intersection = __getCommonValuesInTwoLists(u_partners, v_partners)
	
	if len(intersection) > 0:
		
		# print repr(u_partners)
		# print repr(v_partners)
		# print repr(intersection)
		
		v1 = __calculate1(u, v, intersection, ppi)
		v2 = __calculate2(u, intersection, ppi)
		v3 = __calculate3(u, v, intersection, ppi)
		
		v4 = __calculate2(v, intersection, ppi)
		v5 = __calculate3(v, u, intersection, ppi)
		
		lam1 = __calculateLambda(u, v, u_partners, v_partners, intersection)
		lam2 = __calculateLambda(v, u, v_partners, u_partners, intersection)
		
		# print v1, v2, v3, v4, v5, lam1, lam2
		score1 = v2 + v3 + v1 + lam1
		score2 = v4 + v5 + v1 + lam2
		
		score = 0.0
		if score2 > 0 and score1 > 0:
			score = v1 / score1 * v1 / score2
			
			# confidence .
			score = score * confidence
		
		# print u, v, score
		
		'''
                  v1                          v1
        -----------------------  * --------------------------
         (v2 + v3) + v1 + lam1      (v4 + v5) + v1 + lam2

        '''
		
		return score
	else:
		return 0.0


def __getReliabilitySum(gene, ppi):
	return Pnp_dict[gene]  # no direction


def __randomSamplingBasedOnReliability(total_genes, test_set_gene, mod_pos, ppi, cross_set):
	#    1. test_set_gene
	#    2. mod_pos
	
	new_mod = []
	
	resvoir = []
	# for g in total_genes:
	#	if not g in test_set_gene and not g in mod_pos.keys() and not g in cross_set:
	#		resvoir.append(g)
	
	resvoir = copy.deepcopy(total_genes)
	
	for m in mod_pos:
		
		m_reliability = __getReliabilitySum(m, ppi)
		
		'''
		other_min_reliability = 10000000.0



		# m protein과 modifier 리스트 내의 다른 것들과 비교했을때 reliability 차기 최소 값을 고른다.
		# 이걸 cutoff로 정한다.
		for m2 in mod_pos:

			if m2 == m:
				continue

			m2_reliability = __getReliabilitySum(m2, ppi)
			if abs(m_reliability - m2_reliability)< abs(m_reliability - other_min_reliability):
				other_min_reliability = m2_reliability

		# 전체 protein에서 reliability를 계산했을 때, (random - m) <= (random - 최소 reliability값) 을 만족하면 고른다.
		'''
		
		random.shuffle(resvoir)
		ratio = 0.01  # reliability 값 차이가 원래 reliability의 10% 이내면 받아들인다.
		
		for rg in resvoir:
			
			if not rg in new_mod:
				
				rg_reliability = __getReliabilitySum(rg, ppi)
				
				if abs(rg_reliability - m_reliability) <= ratio * m_reliability:
					new_mod.append(rg)
					break
	
	return new_mod


def __randomSamplingBasedOnReliability_old(total_genes, test_set_gene, mod_pos, ppi, cross_set):
	#    1. test_set_gene
	#    2. mod_pos
	
	new_mod = {}
	
	resvoir = []
	for g in total_genes:
		if not g in test_set_gene and not g in mod_pos and not g in cross_set:
			resvoir.append(g)
	
	# resvoir = total_genes
	
	for m in mod_pos:
		
		m_reliability = __getReliabilitySum(m, ppi)
		
		other_min_reliability = 10000000.0
		for m2 in mod_pos:
			
			if m2 == m:
				continue
			
			m2_reliability = __getReliabilitySum(m2, ppi)
			if abs(m_reliability - m2_reliability) < abs(m_reliability - other_min_reliability):
				other_min_reliability = m2_reliability
		
		while (True):
			
			rg = random.sample(resvoir, 1)[0]
			
			if not rg in new_mod:
				
				rg_reliability = __getReliabilitySum(rg, ppi)
				if abs(rg_reliability - m_reliability) <= abs(rg_reliability - other_min_reliability):
					new_mod[rg] = mod_pos[m]
					break
	
	return new_mod


def logmsg(msg):
	f = open('r:/netexp.log', 'a')
	t = repr(time.ctime())
	f.write(t + '\t' + msg + '\n')
	f.close()


def for_multi_processing(args):
	# MULTIPROCESS에서  사용할 함수 제일 끝에는 multiprocess.Queue() 변수를 넣어줘야함.
	total_genes, new_mod_pos, output, ppi, r, pnp_matrix, pnp_dict, score_queue = args
	
	global Pnp_dict, Pnp_matrix
	Pnp_dict = pnp_dict
	Pnp_matrix = pnp_matrix
	
	'''
    ppi = PPI.PPIDB_STRING()
    ppi.load(fname, 0, 1, 2, cutoff, ' ', None)

    fcache = fname + '.cache'
    Pnp_matrix, Pnp_dict = __buildPnpMatrix(total_genes, ppi, cache_file=fcache)
    '''
	
	# r <-- meaningless
	sc = __calculatePvalueSUB2(total_genes, new_mod_pos, output, ppi, r)
	score_queue.put(sc)
	
	logmsg('Results added into the queue... ' + str(len(sc)))


# return sc


def saveVariableAsPickle(variable):
	global TEMP_PATH
	
	# print 'temp=', TEMP_PATH
	# print 'RND=',  MyUtil_pypy.getRandomString(10)
	
	fname = TEMP_PATH + '/pickle_' + MyUtil.getRandomString(10) + '.obj'
	
	# print fname
	
	with open(fname, "wb") as of:
		pickle.dump(variable, of)
	
	return fname


def loadVariableFromPickle(fname):
	e = None
	
	with open(fname, 'rb') as of:
		e = pickle.load(of)
	
	return e


def __train_k(args):
	
	global TEMP_PATH, Pnp_dict, Pnp_matrix
	global JOBS_PER_CPU
	
	r, total_genes, test_set_gene, mod_box, ppi_file_box, temp_path, job_per_cpu, rnd_id, queue = args
	TEMP_PATH = temp_path
	
	output = None
	
	ret_box = {}
	
	for i in range(JOBS_PER_CPU):
		
		
		# random network 중에서 하나를 그냥 랜덤하게 뽑자....
		
		#ppi_info_dict = random.sample(ppi_file_box[1:], 1)[0]
		ppi_info_dict = ppi_file_box[0] # seed를 randomize하니 network는 그대로 둔다.
		
		cache_file = ppi_info_dict['cache_file']
		ppi_file = ppi_info_dict['ppi_file']
		index1 = ppi_info_dict['index1']
		index2 = ppi_info_dict['index2']
		score_index = ppi_info_dict['score_index']
		threshold = ppi_info_dict['threshold']
		alias_file = ppi_info_dict['alias_file']
		
		ppi = PPIDB_STRING()
		ppi.load(ppi_file, index1, index2, score_index, threshold, ' ', alias_file)
		
		cache_file = ppi_file + '_' + str(threshold) + '_cache.pickle'
		Pnp_matrix, Pnp_dict = __buildPnpMatrix(total_genes, ppi_info_dict)
		
		## Purely random selection
		new_mod_pos = {}
		
		for k in mod_box:
			# new_mod_pos[k] = random.sample(ttt, len(mod_box[k]))
			
			gene_list = __randomSampling(total_genes, None, mod_box[k], None, ppi)
			new_mod_pos[k] = {}
			for g in gene_list:
				new_mod_pos[k][g] = random.random()  # g는 원래 STRING ID임.
		
		pr = __calculatePvalueSUB2(test_set_gene, new_mod_pos, output, ppi_info_dict, r)
		
		for p in pr:
			if not p in ret_box:
				ret_box[p] = []
			ret_box[p].append(pr[p])
		
		#ret_box.append(ret)
	
	r_file = temp_path + '/netexp_' + MyUtil.getRandomString(15) + '_' + rnd_id + '.pickle'
	
	MyUtil.saveVariableAsPickle(ret_box, filename=r_file)
	
	queue.put(r_file)


def remove(fname):
	try:
		os.remove(fname)
	except:
		print('Failed to remove: ', fname)


def __calculatePvalue_mpi(total_genes, test_set_gene, mod_box, output, ppi_file_box):
	
	global P_VALUE, TEMP_PATH, Pnp_matrix, Pnp_dict
	global JOBS_PER_CPU, MP_PROCESSORS, PYPY, ITERATION
	
	time_out = _0_Preprocess.get_default_opt()[_0_Preprocess.OPT_TIMEOUT] * 2
	
	# MP_PROCESSORS = 10
	
	# ------------------------------------
	# PRELOAD_NO_OF_RANDOM_NETWORKS = iteration
	
	r = -1
	
	# my algorithm
	real_scores = __calculatePvalueSUB2(test_set_gene, mod_box, output, ppi_file_box[0], r)
	
	pvalues = {}
	
	if P_VALUE == False:
		print('No p-value calculation: use scores....')
		return real_scores
	
	
	else:
		
		rand = {}
		box = []
		
		# ppi_box_file = ''
		# ppi_box_file = saveVariableAsPickle( ppi_box[0] )
		
		# ppi_file = _0_Preprocess.STRING_FILE
		# ppi_alias_file = _0_Preprocess.STRING_ALIAS_FILE
		# ppi_cutoff = _0_Preprocess.STRING_CUTOFF
		
		# mod_box_file = saveVariableAsPickle( mod_box )
		# total_gene_file = saveVariableAsPickle( total_genes )
		# test_gene_file = saveVariableAsPickle( test_set_gene)
		
		# pnp_dict_file = saveVariableAsPickle(Pnp_dict)
		# pnp_matrix_file = saveVariableAsPickle(Pnp_matrix)
		
		th = []
		
		cpu = MP_PROCESSORS
		jobs_per_cpu = int(ITERATION / cpu)
		
		#for i in range( int( cpu * 1.1) ):
		for i in range(int(cpu)):
			# print ('[p-value] Iteration: ', (i + 1), '/', PRELOAD_NO_OF_RANDOM_NETWORKS)
			
			job = MULTIPROCESS.JOB_THREAD()
			
			args = [r, total_genes, test_set_gene, mod_box, ppi_file_box, TEMP_PATH, jobs_per_cpu, RND_ID]
			job.set_args(__train_k, args)
			job.set_timeout(time_out)
			
			th.append(job)
		
		#ret_files = MULTIPROCESS.runMultiprocesses(th, max_cpu=int(cpu * 1.1), cutoff=cpu, title="[P_VALUE] CrossNetExp...")
		ret_files = MULTIPROCESS.runMultiprocesses(th, max_cpu=cpu, cutoff=cpu,
		                                           title="[P_VALUE] CrossNetExp...")
		
		print('[NetExp] Merging results... \n\n\n')
		
		for rfile in ret_files:
			
			print('loading netexp result: ', rfile, time.ctime())
			ret = MyUtil.loadVariableFromPickle(rfile, remove_after_load=True)
			
			for gene in ret:
			
				if not gene in rand:
					rand[gene] = []
				
				rand[gene] += ret[gene]
	
		print('Calculating p-values...')
		
		temp_scores = {}
		scores = {}
		pvalues = {}
		
		for g in rand:
			avg = statistics.average(rand[g])
			std = statistics.stdev(rand[g])
			p = 1.0
			try:
				if std != 0.0:
					# z score 이용
					p = float(( real_scores[g] - avg )/ float(std) )
					#p = statistics.getNormalDistPvalue(avg, std, real_scores[g])
				else:
					p = -10000000
			except:
				p = -10000000
			
			scores[g] = p
			
			txt = 'Score:' + str(real_scores[g]) + ', avg=' + str(avg) + ' std=' + str(std)
			pvalues[g] = [scores[g], real_scores[g], txt]  # pvalue, real score, comment
		
		return pvalues  # p-value based on z-score


'''
def saveVariableAsPickle(variable):

	global TEMP_PATH

	print ('temp=', TEMP_PATH)
	print ('RND=',  MyUtil.getRandomString(10))

	fname = TEMP_PATH + '/pickle_' + MyUtil.getRandomString(10) + '.obj'

	print (fname)

	with open(fname, "wb") as of:
		pickle.dump(variable, of)

	return fname

def loadVariableFromPickle(fname):

	e = None

	with open(fname, 'rb') as of:
		e = pickle.load(of)

	return e
'''


def __cal_rnd_network(args):
	global Pnp_matrix, Pnp_dict
	
	interaction_file, alias_file, test_set_gene, mod_box, output, r, queue = args
	interaction_pickle = interaction_file + '.pickle'
	matrix_pickle = interaction_file + '.matrix.pickle'
	dict_pickle = interaction_file + '.dict.pickle'
	
	ppi_rnd = None
	
	'''
	if os.path.exists(interaction_pickle) and \
			os.path.exists(matrix_pickle) and \
			os.path.exists(dict_pickle):

		ppi_rnd = loadVariableFromPickle(interaction_pickle)
		Pnp_matrix = loadVariableFromPickle(matrix_pickle)
		Pnp_dict = loadVariableFromPickle(dict_pickle)

	else:

		ppi_rnd = PPI.PPIDB_STRING()

		# ppi.load(interaction_file, index1, index2, score_index, threshold, ' ', _0_Preprocess.STRING_ALIAS_FILE)
		ppi_rnd.load(interaction_file, 0, 1, 2, -1, ' ', alias_file)
		Pnp_matrix, Pnp_dict = __buildPnpMatrix(ppi_rnd.getWholeGeneList(), ppi_rnd, cache_file=interaction_file+'.cache.txt')


		saveVariableAsPickle(ppi_rnd, interaction_pickle)
		saveVariableAsPickle(Pnp_matrix,matrix_pickle)
		saveVariableAsPickle(Pnp_dict, dict_pickle)
	'''
	
	ppi_rnd = PPIDB_STRING()
	# ppi.load(interaction_file, index1, index2, score_index, threshold, ' ', _0_Preprocess.STRING_ALIAS_FILE)
	ppi_rnd.load(interaction_file, 0, 1, 2, -1, ' ', alias_file)
	Pnp_matrix, Pnp_dict = __buildPnpMatrix(ppi_rnd.getWholeGeneList(), ppi_rnd,
	                                        cache_file=interaction_file + '.cache.txt')
	
	ret = __calculatePvalueSUB2(test_set_gene, mod_box, output, ppi_rnd, r)
	
	fname = 'R:/TEMP/' + MyUtil.getRandomString(15) + '.obj'
	saveVariableAsPickle(ret, fname)
	
	queue.put(fname)


def __getTrim(total_genes, mod_pos, ppi):
	total2 = []
	score = {}
	for g in total_genes:
		score[g] = 0.0  #
	
	for g in mod_pos:
		for g2 in ppi.getPartners(g):
			
			if not g2 in total2:
				total2.append(g2)
			
			for g3 in ppi.getPartners(g2):
				if not g3 in total2:
					total2.append(g3)
	
	return score, total2


def __getOutgoingW(u, ppi):
	total = 0.0
	for k in ppi.getOutgoingPartnersOf(u):
		total = total + ppi.getScore(u, k)
	return total


def __getIncomingW(u, ppi):
	total = 0.0
	for k in ppi.getIncomingPartnersOf(u):
		total = total + ppi.getScore(k, u)
	return total


def __sortStr(a, b):
	if a > b:
		return a + ':' + b
	else:
		return b + ':' + a


def __getW(u, ppi):
	total = 0.0
	for k in ppi.getPartners(u):
		total = total + ppi.getScore(k, u)
	return total


def __buildPnpMatrix(total_gene_names, ppi_info_dict):
	
	dic = {}
	
	cache_file = ppi_info_dict['cache_file']
	ppi_file = ppi_info_dict['ppi_file']
	index1 = ppi_info_dict['index1']
	index2 = ppi_info_dict['index2']
	score_index = ppi_info_dict['score_index']
	threshold = ppi_info_dict['threshold']
	alias_file = ppi_info_dict['alias_file']
	
	if os.path.exists(cache_file):
		print('[12_NetExp] Loading pnp matrix from cache = ', cache_file)
		dic = MyUtil.loadVariableFromPickle(cache_file)
		return None, dic
	else:
		
		ppi = PPIDB_STRING()
		ppi.load(ppi_file, index1, index2, score_index, threshold, ' ', alias_file)
		
		c = 0.0
		q = float(len(total_gene_names))
		for u in total_gene_names:
			c += 1.0
			# print ('\r', (c/q*100.0), '%')
			
			dic[u] = __getW(u, ppi)
		
		MyUtil.saveVariableAsPickle(dic, filename=cache_file)
		print(['_12_NetExp] pnp matrix cache saved = ', cache_file])
		return None, dic


def runFinal(iteration, interaction_file, index1, index2, index3, score_index, threshold, modI=None, test_set_all=[]):
	global P_VALUE
	global RND_ID, ITERATION, NETWORK_TEMP_FILES, MP_PROCESSORS, DEGREE, DISEASE_SET, R_PARAMETER
	
	mod_maps, mods = _0_Preprocess.getModifiers()
	
	summary_file = TEMP_PATH + '/NetExp' + RND_ID + '_' + '_summary.txt'
	f = open(summary_file, 'w')
	dd = list(mods)
	dd.sort()
	
	NETWORK_TEMP_FILES = [interaction_file, index1, index2, score_index, threshold]
	
	info = '--------------------------------------------' + '\n' + \
	       'THREAD = ' + str(MP_PROCESSORS) + '\n' + \
	       'Remove temporary files = ' + str(REMOVE_TEMP_FILES) + '\n' + \
	       'Calculate p-values = ' + str(P_VALUE) + '\n' + \
	       '\t# of random networks = ' + str(ITERATION) + '\n' + \
	       'Iteration = ' + str(HOW_MANY_TEST_MODIFERS) + '\n' + \
	       'Randomization method = ' + str(DEGREE) + '\n'
	print(info)
	
	ppi_box = loadInteraction(interaction_file, ITERATION, index1, index2, index3, score_index,
	                          threshold)  # 100개 네트워크를 만든다.
	total_genes = __getTotalGeneList(ppi_box)
	
	cache_file = copyCache()
	
	print('init...')
	global Pnp_matrix, Pnp_dict
	total_gene_names = ppi_box[0].getWholeGeneList()
	total_gene_names.sort()
	Pnp_matrix, Pnp_dict = __buildPnpMatrix(total_gene_names, ppi_box[0], cache_file=cache_file)
	print('done')
	
	ppi_box2 = []
	ppi_box2 = ppi_box
	
	mod_box = {}
	
	for key in dd:
		mod_pos = mods[key]
		
		new_mod_pos = __removeUnannotatedModifiers(total_genes, mod_pos)
		
		mod_box[key] = new_mod_pos
	
	s = 'Interaction = ' + interaction_file + '\n' + \
	    'Disease = ' + repr(dd) + '\n' + \
	    'Total annotated gene # = ' + str(len(total_genes)) + '\n'
	
	f.write(info + '\n')
	f.write('-------------------------\n' + s + '\n')
	
	print(s)
	
	fname = '__Common_modifiers.txt'
	
	if _0_Preprocess.OUTPUT_FILENAME is not None:
		fname = _0_Preprocess.OUTPUT_FILENAME
	
	output = _0_Preprocess.OUTPUT_FOLDER + '/' + fname
	
	ofiles = __train2(total_genes, mod_box, output, total_genes, ppi_box)
	
	# pvalue= ofiles[1]
	# real scores = ofiles[2]
	# pvalue + real scores = ofiles[3]
	
	return ofiles


def copyCache():
	fname = '_12_' + os.path.basename(_0_Preprocess.STRING_FILE) + \
	        '_' + str(_0_Preprocess.STRING_CUTOFF).strip() + '.cache.txt'
	
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


def __train2(total_genes, mod_box, output, test_set_all, ppi_box):
	# mpi ????
	pvalues = __calculatePvalue_mpi(total_genes, test_set_all, mod_box, output, ppi_box)
	
	# pvalues = __calculatePvalue(total_genes, test_set_all, mod_box, output, ppi_box)
	
	ofiles = __resave(total_genes, mod_box, pvalues, output, total_genes, ppi_box)
	return ofiles


def __getSeedScore(total_gene_names, mod, ppi):
	r = {}
	for k in total_gene_names:
		r[k] = 0.0
	
	total = 0.0
	for k in mod:
		total = total + mod[k]
	
	for k in mod:
		r[k] = mod[k] / total
	
	arr = []
	for k in total_gene_names:
		arr.append([r[k]])
	# arr.append(  r[k] )
	
	return arr, r


def __calculateScore(p, k, j, new_box, r, i, ppi):
	global Pnp_matrix, Pnp_dict
	
	key = p + ':' + __sortStr(k, j)
	w = Pnp_dict[key]
	
	s = r * (new_box[k] * new_box[j]) * w
	
	return s


def __categorizePartners(mod_box, partners):
	r_box = {}
	dis = list(mod_box)
	dis.sort()
	
	for p in partners:
		
		for disease in dis:
			mp = mod_box[disease]
			
			if p in mp:
				if not disease in r_box:
					r_box[disease] = []
				r_box[disease].append(p)
	
	q_box = []
	cbu = 0
	cbs = ''
	
	for d in dis:
		if not d in r_box:
			return None
		else:
			
			cbu = cbu + len(r_box[d]) - 1
			if cbs == '':
				cbs = r_box[d][0]
			else:
				if cbs != r_box[d][0]:
					cbu = cbu + 1
	
	if cbu == 0:
		return None
	
	return r_box


def __calculatePvalueSUB2(total_genes, mod_box, output, ppi_info_dict, r):
	
	global Pnp_matrix, Pnp_dict
	
	cache_file = ppi_info_dict['cache_file']
	ppi_file = ppi_info_dict['ppi_file']
	index1 = ppi_info_dict['index1']
	index2 = ppi_info_dict['index2']
	score_index = ppi_info_dict['score_index']
	threshold = ppi_info_dict['threshold']
	alias_file = ppi_info_dict['alias_file']
	
	ppi = PPIDB_STRING()
	ppi.load(ppi_file, index1, index2, score_index, threshold, ' ', alias_file)
	
	total_gene_names = ppi.getWholeGeneList()
	Pnp_matrix, Pnp_dict = __buildPnpMatrix(total_gene_names, ppi_info_dict)
	
	# global R_THRESHOLD
	
	box = {}
	# threshold = R_THRESHOLD
	prev_mod_box = {}
	
	print('Seed genes normalization')  # rank ratio
	for disease in mod_box:
		total = 0.0
		mod_pos = mod_box[disease]
		
		prev_mod_box[disease] = {}
		
		print(disease, '-->', len(mod_pos))
		
		for k in mod_pos:
			total = total + mod_pos[k]
		
		for k in mod_pos:
			prev_mod_box[disease][k] = mod_pos[k] / total
	
	prev_cache = {}
	
	cnt = 0
	
	dis = list(mod_box)
	dis.sort()
	
	identified_common_modifiers = {}
	
	while (True):
		
		cnt = cnt + 1
		if cnt > 100:
			break
		
		new_mod_box = {}
		# new_cache = copy.deepcopy( prev_cache  )
		
		for g in ppi.getWholeGeneList():
			
			if g in identified_common_modifiers:
				# print g, 'is a generic modifier'
				continue
			
			partners = ppi.getPartners(g) + [g]
			
			''' common modifier '''
			cat_partners = __categorizePartners(prev_mod_box, partners)
			
			if cat_partners == None:
				# print 'skip ,, ', g
				continue
			
			sc = 1.0 / (Pnp_dict[g])
			
			# score = 1/(incombin wet sum) * sigma[ score(p)*reliability(from,to)/(outgoing w sum)  ]
			for disease_name in cat_partners:
				
				disease_modifiers = cat_partners[disease_name]
				
				qsco = 0.0
				
				for dm in disease_modifiers:
					
					q = None
					if g == dm:
						# self
						# weight = 1
						# ppi.getScore(dm,g)/Pnp_dcit[dm][0] = 1
						q = prev_mod_box[disease_name][dm]
					
					else:
						
						q = prev_mod_box[disease_name][dm] * ppi.getScore(dm, g) / Pnp_dict[dm]
					
					qsco = qsco + q
				
				sc = sc * qsco
				
				'''
                for di in dis:
                    if not new_mod_box.has_key(di):
                        new_mod_box[di]={}

                    new_mod_box[di][g] = sc
                '''
				
				# print 'found ', g, sc
				new_mod_box[g] = sc
		# new_cache[g]=sc
		
		for g in new_mod_box:
			
			identified_common_modifiers[g] = new_mod_box[g]
			
			for d in dis:
				prev_mod_box[d][g] = new_mod_box[g]
		
		print('step = ', cnt, ' identified at this step = ', len(new_mod_box), ' identified common modifiers so far = ',
		      len(identified_common_modifiers))
		if len(new_mod_box) == 0:
			print('no more modifiers found')
			break
	
	# new_mod_box = {}
	
	for g in total_genes:
		if not g in identified_common_modifiers:
			identified_common_modifiers[g] = 0.0
	
	return identified_common_modifiers


def __calculatePvaluSUB(total_genes, mod_box, output, ppi, r):
	prev_mod_box = {}
	
	# seed gene normalization
	for disease in mod_box:
		total = 0.0
		mod_pos, mod_neg = mod_box[disease]
		
		prev_mod_box[disease] = {}
		
		for k in mod_pos:
			total = total + mod_pos[k]
		
		for k in mod_pos:
			prev_mod_box[disease][k] = mod_pos[k] / total
	
	prev_cache = {}  # common modifier
	
	# print 'mod # = ', len(mod_pos)
	
	cnt = 0
	
	dis = list(mod_box)
	dis.sort()
	
	while (True):
		
		cnt = cnt + 1
		if cnt > 100:
			break
		
		new_mod_box = {}
		new_cache = copy.deepcopy(prev_cache)
		
		for g in ppi.getWholeGeneList():
			
			if g in prev_cache:
				continue
			
			partners = ppi.getPartners(g) + [g]
			ptx = [] + partners
			
			for p in partners:
				psh = ppi.getPartners(p)
				for p2 in psh:
					if not p2 in ptx:
						ptx.append(p2)
			
			partners = ptx
			
			cat_partners = __categorizePartners(prev_mod_box, partners)
			
			if cat_partners == None:
				# print 'skip ,, ', g
				continue
			
			sc = 1.0
			
			for dindex in range(len(cat_partners)):
				
				disease_modifiers = cat_partners[dindex]
				disease_name = dis[dindex]
				
				fscore = 1.0
				
				for dm in disease_modifiers:
					key = __sortStr(dm, g)
					fscore = fscore * (1.0 - Pnp_dict[key] * prev_mod_box[disease_name][dm])
				
				fscore = 1.0 - fscore
				
				sc = sc * fscore
			
			for di in dis:
				if not di in new_mod_box:
					new_mod_box[di] = {}
				
				new_mod_box[di][g] = sc
			
			new_cache[g] = sc
		
		for di in mod_box:
			mp, np = mod_box[di]
			
			if not di in new_mod_box:
				new_mod_box[di] = {}
			
			nn = new_mod_box[di]
			
			for gk in mp:
				if not gk in nn:
					nn[gk] = mp[gk]
		
		print('step = ', cnt, ' # = ', len(new_cache))
		
		if len(prev_cache) == len(new_cache):
			prev_mod_box = new_mod_box
			prev_cache = new_cache
			break
		
		prev_mod_box = new_mod_box
		prev_cache = new_cache
	
	ret = {}
	for g in total_genes:
		if not g in new_cache:
			ret[g] = 0.0
		else:
			ret[g] = new_cache[g]
	
	return ret


def __countModifierNumber(mod_list, gene_list):
	cnt = 0
	for g in gene_list:
		if g in mod_list:
			cnt = cnt + 1
	return cnt


def __test(train_output2, test_pos, test_others):
	f = open(train_output2, 'r')
	
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
				
				gid = x[0].replace('m_', '')
				
				if gid == test_pos:
					cnt = cnt + 1
					if pos_rank == -1:
						pos_rank = cnt
					else:
						print('Error!!!! two or more positives')
				
				elif gid in test_others:
					cnt = cnt + 1
	
	f.close()
	
	if cnt != 100:
		print('Error !!! Not 100 test set', cnt)
	
	return pos_rank


def __loadSavedFile(infile):
	f = open(infile, 'r')
	init = True
	
	r = {}
	
	for s in f.readlines():
		if init:
			init = False
		else:
			s = s.replace('\n', '')
			x = s.split('\t')
			
			genes = x[5].split(',')
			score = x[1]
			goid = x[0]
			pvalue = x[2]
			
			for g in genes:
				if not g in r:
					r[g] = [[], []]
				
				# r[g][0].append( eval(score))
				r[g][0].append(eval(pvalue))
				r[g][1].append(goid)
	
	f.close()
	
	return r


def __calculateAUC2(p_rank, n_rank, output):
	f = open(output + '.txt', 'w')
	f.write('False-positive rate(1-specificity)\tTrue-positive rate(sensitivity)\n')
	
	start_from = 0.0
	start_to = 1.0
	step = 0.01
	
	pos_false_positive_rate_x = []  # 1- specificity
	pos_true_positive_rate_y = []  # sensitivity
	
	threshold = start_from
	while (threshold <= start_to):
		threshold = threshold + step
		
		specificity = 1.0 - threshold
		
		total_trial = float(len(p_rank))
		total_rank = 100.0
		
		success = 0.0
		
		for index in range(len(p_rank)):
			p_v = p_rank[index]
			n_v = n_rank[index]
			
			p_r = float(p_v) / total_rank  # rank ratio
			n_r = float(n_v) / total_rank
			
			if p_r <= threshold < n_r:
				success = success + 1.0
		
		sensitivity = success / total_trial
		
		pos_false_positive_rate_x.append(1.0 - specificity)
		pos_true_positive_rate_y.append(sensitivity)
		
		f.write(str(1.0 - specificity) + '\t' + str(sensitivity) + '\n')
	
	auc = MyUtil.mean(pos_true_positive_rate_y)
	f.close()
	
	return auc


def __calculateAUC(p_rank, output):
	f = open(output + '_pos.txt', 'w')
	f.write('False-positive rate(1-specificity)\tTrue-positive rate(sensitivity)\n')
	
	start_from = 0.0
	start_to = 1.0
	step = 0.01
	
	pos_false_positive_rate_x = []  # 1- specificity
	pos_true_positive_rate_y = []  # sensitivity
	
	threshold = start_from
	while (threshold <= start_to):
		threshold = threshold + step
		
		specificity = 1.0 - threshold
		
		total_trial = float(len(p_rank))
		total_rank = 100.0
		
		success = 0.0
		for v in p_rank:
			r = float(v) / total_rank  # rank ratio
			if r <= threshold:
				success = success + 1
		
		sensitivity = success / total_trial
		
		pos_false_positive_rate_x.append(1.0 - specificity)
		pos_true_positive_rate_y.append(sensitivity)
		
		f.write(str(1.0 - specificity) + '\t' + str(sensitivity) + '\n')
	
	auc_pos = MyUtil.mean(pos_true_positive_rate_y)
	f.close()
	
	return auc_pos


def getModifiersRank(fname, threshold_rank, ignore_rank=False, ignore_first_line=True):
	
	
	print('getModifiersRank: ', fname)
	
	f = open(fname, 'r')
	
	rank = 0.0
	total = float(threshold_rank)
	
	r = {}
	
	if ignore_first_line:
		f.readline()
	
	#print('Loading a modifier file', fname, '....')
	
	for s in f.readlines():
		
		s = s.strip()
		
		s = s.replace('m_', '').replace('x_', '')
		
		x = s.split('\t')
		
		gene_id = x[0].strip()
		# pvalue = eval( x[1] ) # meaningless
		string_id = _0_Preprocess.ID_BOX.getSTRING_ID_of(gene_id)
		
		if string_id is not None:
			
			if ignore_rank:
				r[string_id] = 1.0
				rank = rank + 1.0
			
			else:
				r[string_id] = 1.0 - rank / total
				rank = rank + 1.0
			
			# print gene_id, '\t', r[gene_id]
			
			if rank >= total:
				break
		else:
			print('\t[_12_CrossNetworkExpansion.getModifiersRank] Unrecognized protein name: ', gene_id)

	f.close()
	
	print('getModifiersRank [done]: ', fname, len(r))
	
	
	return r


def getPredictedModifiers(diseases, threshold_rank):
	global RND_ID
	
	# 여기에 리스트 형태로 들어가 있는 질병에 대해서만 disease-specific modifiers를 예측한다.
	predict_disease_specific_modifiers = _0_Preprocess.PREDICT_DISEASE_MODIFIERS
	
	r = {}
	
	for d in diseases:
		
		if d in predict_disease_specific_modifiers:
			# disease-specific modifiers를 예측했으므로 그 결과를 갖다 쓴다.
			fname = OUTPUT_FOLDER + '/Integrated_' + d + '_' + RND_ID + '.txt'
			v = getModifiersRank(fname, threshold_rank, ignore_first_line=True)
			# r[d] = [ v, {} ]
			r[d] = v
			print(d, fname, len(v), 'genes loaded...')
		else:
			# query로 들어온 유전자를 그대로 사용한다.
			# 정상적인 경우엔 아래 코드를 실행할 일이 없음.
			#============================================
			TEST = False
			#==========================================
			
			threshold_rank = 100000 # 사용자가 입력한거라서 그냥 다 이용한다.
			
			
			if TEST:
				fname = _0_Preprocess.MODIFIER_FILES[d]
				v = getModifiersRank(fname, threshold_rank, ignore_first_line=False)
				# r[d] = [ v, {} ]
				r[d] = v
				print(d, fname, len(v), 'ranked genes loaded...')
			
			else:
				# ============================================================
				# Query 들어온 걸 순서 상관없이 사용하기 위해 아래 코드 사용해야함.
				fname = _0_Preprocess.MODIFIER_FILES[d]
				x, filtered_genes = __loadModifier(fname)
				x = {}
				for g in filtered_genes:
					x[g] = 1.0  # 이 경우 rank 가 의미없으므로 모두 seed 값의 max인 1을 부여함.
					if len(x) >= threshold_rank:
						# seed 숫자보다 많으면 그냥 끊는다.
						break
				r[d] = x
				print(d, fname, len(x), 'genes loaded... (no rank)')
	
	return r


def __loadModifier(fname):
	r = {}
	l = {}
	
	f = open(fname, 'r', encoding='utf-8')
	
	# if DISPLAY:
	print('Loading ' + fname)
	
	hash_md5 = ''
	
	for s in f.readlines():
		s = s.replace('\n', '')
		if len(s) == 0:
			continue
		
		if s[0] == '#':
			continue  # comment
		
		x = s.split('\t')
		gene_id = x[0].strip()  # gene id
		string_id = _0_Preprocess.ID_BOX.getSTRING_ID_of(gene_id)
		if string_id is None:
			print('Unrecognized ID = ' + gene_id)
			r[gene_id] = None
		else:
			r[gene_id] = string_id
	
	f.close()
	
	for user_id in list(r):
		
		string_id = r[user_id]
		if string_id is not None:
			l[string_id] = None
	
	print('Queried proteins, # = %d ' % (len(r)))
	print('\tLoaded proteins, # = %d ' % (len(l)))
	
	return r, list(l)


def FinalizeAll(opt):
	
	global OPTIONS, P_VALUE
	
	# P_VALUE = True # 임시로 삽입한 것임.
	
	# seed_number = _0_Preprocess.SEED_NUMBER
	seed_number = _0_Preprocess.getSEEDNumber()
	output_pickle = opt[_0_Preprocess.OPT_OUTPUT]
	
	print('OUTPUT_PICKLE = ', output_pickle)
	
	
	
	
	string_db_file = _0_Preprocess.STRING_FILE
	index1 = 0
	index2 = 1
	index3 = 3
	score_index = 2
	threshold = _0_Preprocess.STRING_CUTOFF
	alias_file = _0_Preprocess.STRING_ALIAS_FILE
	
	maps, mods = _0_Preprocess.getModifiers()
	diseases = list(mods)
	
	
	
	
	#print('9999999999999999999999')
	#print(_0_Preprocess.ID_BOX.getSTRING_ID_of('Akt1'))
	#print('9999999999999999999999')
	
	
	
	# seeds 가져올떄 disease-specific modifiers 예측해서 가져올지, user query를 그대로 쓸지 결정함.
	seeds = getPredictedModifiers(diseases, seed_number)  # dict, key=protein, value=score(rank order)
	# cache_file = string_db_file + '_' + str(threshold) + '_cache.pickle'
	
	NETWORK_TEMP_FILES = [string_db_file, index1, index2, score_index, threshold]
	
	# 원래 PPI, random PPI까지 모두 파일이름을 읽어온다. PPI를 만드는 건 아님.
	ppi_file_box = loadInteraction(string_db_file, index1, index2, index3, score_index, threshold)  # 100개 네트워크를 만든다.
	
	print('loading basic ppi')
	ppi = PPIDB_STRING()
	ppi.load(string_db_file, index1, index2, score_index, threshold, ' ', alias_file)
	
	total_genes = __getTotalGeneList(ppi)
	
	print('init...')
	global Pnp_matrix, Pnp_dict
	total_gene_names = ppi.getWholeGeneList()
	total_gene_names.sort()
	# local_Pnp_matrix, local_Pnp_dict = __buildPnpMatrix(total_gene_names, ppi_box[0])
	Pnp_matrix, Pnp_dict = __buildPnpMatrix(total_gene_names, ppi_file_box[0])
	# Pnp_matrix=None이 리턴되네...
	print('done')
	
	ppi_box2 = []
	# ppi_box2 = copy.deepcopy(ppi_box)
	
	mod_box = {}
	
	for key in diseases:
		# if key.find('RANDOM') > 0:
		#    continue
		
		mod_pos = seeds[key]  # [ v{}, {} ]  # v contains key=gene, value=score --> v{} only
		
		print('------------------------')
		print(key)
		
		# Pnp_matrix, Pnp_dict = copy.deepcopy(local_Pnp_dict), copy.deepcopy(local_Pnp_matrix)
		
		#print(mod_pos)
		new_mod_pos = __removeUnannotatedModifiers(total_genes, mod_pos)
		
		mod_box[key] = new_mod_pos
		
		print("Filtered seeds: ", key, len(new_mod_pos))
	
	total_genes_without_modifiers = __getTotalGenesWithoutModifiers(total_genes, mod_box)
	
	fname = '__Common_modifiers.txt'
	if _0_Preprocess.OUTPUT_FILENAME is not None:
		fname = _0_Preprocess.OUTPUT_FILENAME
	# if '[RND]' in fname:
	#	fname = fname.replace('[RND]', RND_ID)
	
	output = _0_Preprocess.OUTPUT_FOLDER + '/' + fname
	
	ofiles = __train2(total_genes, mod_box, output, total_genes, ppi_file_box)
	
	print('######################################################')
	print(' Result files are : ')
	for o in ofiles:
		print('==>', repr(o))
	print('######################################################')
	
	# 여러 종류의 결과를 넘긴다.
	
	# pickle로 저장한다.
	MyUtil.saveVariableAsPickle(ofiles, filename=output_pickle)
	
	return ofiles


def FinalizeAll_random():
	# seed_number = _0_Preprocess.SEED_NUMBER
	seed_number = _0_Preprocess.getSEEDNumber()
	string_db_file = _0_Preprocess.STRING_FILE
	index1 = 0
	index2 = 1
	index3 = 3
	score_index = 2
	threshold = _0_Preprocess.STRING_CUTOFF
	
	maps, mods = _0_Preprocess.getModifiers()
	diseases = list(mods)
	
	NETWORK_TEMP_FILES = [string_db_file, index1, index2, score_index, threshold]
	
	ppi_box = loadInteraction(string_db_file, PRELOAD_NO_OF_RANDOM_NETWORKS, index1, index2, index3, score_index,
	                          threshold)  # 100개 네트워크를 만든다.
	total_genes = __getTotalGeneList(ppi_box)
	
	# seeds = getPredictedModifiers(diseases, seed_number) # dict, key=protein, value=score(rank order)
	
	# randomization -------------------------------------------------------
	seeds = {}
	total = float(len(total_genes))
	
	for k in diseases:
		seeds[k] = {}
		
		temp = random.sample(total_genes, seed_number)
		
		print('[', k, '] Randomly selected = ', len(temp))
		
		for i in range(len(temp)):
			seeds[k][temp[i]] = 1.0 - float(i) / total
	
	print('init...')
	global Pnp_matrix, Pnp_dict
	total_gene_names = ppi_box[0].getWholeGeneList()
	total_gene_names.sort()
	# local_Pnp_matrix, local_Pnp_dict = __buildPnpMatrix(total_gene_names, ppi_box[0])
	Pnp_matrix, Pnp_dict = __buildPnpMatrix(total_gene_names, ppi_box[0])
	print('done')
	
	ppi_box2 = []
	# ppi_box2 = copy.deepcopy(ppi_box)
	
	mod_box = {}
	
	for key in diseases:
		# if key.find('RANDOM') > 0:
		#    continue
		
		mod_pos = seeds[key]  # [ v{}, {} ]  # v contains key=gene, value=score --> v{} only
		
		print('------------------------')
		print(key)
		
		# Pnp_matrix, Pnp_dict = copy.deepcopy(local_Pnp_dict), copy.deepcopy(local_Pnp_matrix)
		
		new_mod_pos = __removeUnannotatedModifiers(total_genes, mod_pos)
		
		mod_box[key] = new_mod_pos
	
	total_genes_without_modifiers = __getTotalGenesWithoutModifiers(total_genes, mod_box)
	
	fname = '__Common_modifiers_random.txt'
	if _0_Preprocess.OUTPUT_FILENAME is not None:
		fname = _0_Preprocess.OUTPUT_FILENAME + '_random.txt'
	# if '[RND]' in fname:
	#	fname = fname.replace('[RND]', RND_ID)
	
	output = _0_Preprocess.OUTPUT_FOLDER + '/' + fname
	
	ofiles = __train2(total_genes, mod_box, output, total_genes, ppi_box)
	
	# pvalue= ofiles[1]
	# real scores = ofiles[2]
	# pvalue + real scores = ofiles[3]
	
	print('[OUTPUT] = ', output)
	
	return output


def init(config_file, opt):
	global RND_ID, TEMP_PATH, OUTPUT_FOLDER, P_VALUE
	
	_0_Preprocess.init(config_file, rnd_id=RND_ID)
	
	# RND_ID = _0_Preprocess.RND_ID
	# RND_ID = opt[_0_Preprocess.OPT_RND_ID]
	# _0_Preprocess.RND_ID = RND_ID
	
	TEMP_PATH = _0_Preprocess.TEMP_PATH
	OUTPUT_FOLDER = _0_Preprocess.OUTPUT_FOLDER
	P_VALUE = _0_Preprocess.NETWORK_EXPANSION_PVALUE
	# P_VALUE = _0_Preprocess.P_VALUE
	
	print('***** Loaded P-VALUE= ', P_VALUE)


def main(opt):
	
	global PYPY, PRELOAD_NO_OF_RANDOM_NETWORKS, MP_PROCESSORS, P_VALUE_CUTOFF, OPTIONS
	global RND_ID, P_VALUE, ITERATION
	
	OPTIONS = opt
	
	config_file = opt[_0_Preprocess.OPT_CONFIG_FILE]
	MP_PROCESSORS = opt[_0_Preprocess.OPT_MULTIMP]
	rnd_or_not = opt[_0_Preprocess.OPT_RANDOM_OR_NOT]
	test_or_not = opt[_0_Preprocess.OPT_TEST]
	ITERATION = opt[_0_Preprocess.OPT_ITERATION]
	
	RND_ID = opt[_0_Preprocess.OPT_RND_ID]
	_0_Preprocess.RND_ID = RND_ID
	
	# 이게 iteration이다.
	# PRELOAD_NO_OF_RANDOM_NETWORKS = opt[_0_Preprocess.OPT_ITERATION]
	
	PYPY = opt[_0_Preprocess.OPT_PYPY]
	P_VALUE_CUTOFF = opt[_0_Preprocess.OPT_PVALUE_CUTOFF]
	
	# print ('---- [CROSS] ', P_VALUE_CUTOFF)
	
	init(config_file, opt)
	
	start_time = datetime.now()
	
	# python _1_GO.py config.txt test
	print('''
========================================================
    [12] NETWORK EXPANSION
========================================================
''')
	
	ofiles = None
	
	if config_file is not None:
		ofiles = FinalizeAll(opt)
	# output = ofiles[0]
	# p-valie = ofiles[1]
	# real score = ofiles[2]
	# pvalie+real score = ofiles[3]
	
	else:
		
		if test_or_not == True:
			run()  # iteration 의미 없음.
		elif rnd_or_not == True:
			FinalizeAll_random()
	
	print('''
========================================================    
    [12] NETOWKR EXPANSION (End)
========================================================          
''')
	
	end_time = datetime.now()
	print('Elapsed time: {}'.format(end_time - start_time))
	print(time.ctime())
	
	for o in ofiles:
		print('____  -> ' + repr(o))
	
	return ofiles


if __name__ == '__main__':
	multiple_try = True
	attempt = 3
	
	opt = _0_Preprocess.process_args(sys.argv[1:])
	ofiles = main(opt)
	
	'''
	else:

		ofiles_box = []
		for i in range(attempt):
			ofiles=main(opt)
			ofiles_box.append(ofiles)

			# pvalue+real
			ofiles_box.append(ofiles[3])
			# 0 - original
			# 1 - pvalue
			# 2 - real
			# 3 - pvalue+real


		# order statistics
	'''



