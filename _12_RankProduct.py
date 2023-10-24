# -*- coding: ms949 -*-
import time
import datetime
import os
import rankproduct
#import statistics
import _0_Preprocess
import copy
import random
import math
import PPI
import threading
import multiprocessing
import MULTIPROCESS
import _100_NetworkExpansionModel
import MyUtil_pypy
from datetime import datetime
import time
import subprocess
import sys

try:
	import statistics
except:
	pass


RND_ID = _0_Preprocess.RND_ID
TEMP_PATH = _0_Preprocess.TEMP_PATH
OUTPUT_FOLDER = _0_Preprocess.OUTPUT_FOLDER


Pnp_matrix = None
Pnp_dict = {}
DIRECTED_PPI = False




# =================================
MAX_THREAD = 5
REMOVE_TEMP_FILES = True

P_VALUE = _0_Preprocess.NETWORK_EXPANSION_PVALUE
PRELOAD_NO_OF_RANDOM_NETWORKS = 100

HOW_MANY_TEST_MODIFERS = 100
NETWORK_TEMP_FILES = []

# 여기가 중요함.
DEGREE = 0 #
# 0 - 그냥 진짜로 random
# 1 - degree가 +-1 범위내의 것만 골라서
# 3 - reliability 차이가 10% 이내인 것들만 사용
# =================================




PYPY = False







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
		mod_neg = self.values[2]
		total_genes = self.values[3]
		total_genes_without_modifiers = self.values[4]
		iteration_th = self.values[5]
		ppi_box = self.values[6]

		self.result = mainJob(disease, mod_pos, mod_neg, total_genes, total_genes_without_modifiers, iteration_th, ppi_box)

		self.over = True
	def getResult(self):
		return self.result


def __randomSampling (total_genes, cross_set, mod_box, basket, ppi):

	new_pos_mod = []


	global DEGREE # default 

	if DEGREE == 0:
		
		# 완전 랜덤하게 고른다.
		new_pos_mod =  random.sample( total_genes, len(mod_box) )

	elif DEGREE == 1:
		# preserving interaction degrees
		basket = getBasket(total_genes, ppi)

		new_pos_mod = []


		for g in mod_box:

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


			random.shuffle(box)

			for v in box:
				if not v in new_pos_mod:
					new_pos_mod.append(v)
					break

		if len(mod_box) != len(new_pos_mod):
			print ('Random sampling error: ', len(mod_box), len(new_pos_mod))





	elif DEGREE == 3:




		# reliability
		# test 및 cross val set에는 포함되지 않는 걸로 골라야하는데, 여기선 고려 안해도 되니 [] 빈 리스트로 넘김
		new_pos_mod = __randomSamplingBasedOnReliability(total_genes, [], mod_box, ppi, [])
		
		if len(mod_box) != len(new_pos_mod):
			print ('Random sampling error: ', len(mod_box), len(new_pos_mod))



	return new_pos_mod


def __getOutputFilename(fname, threshold):


	while(True):
		rnd = str(random.randint(0,1000000000)).strip()

		ppi_filename_only = os.path.basename(fname)
		ppi_filename_only = os.path.splitext(ppi_filename_only)[0]
		ppi_filename_only = './6. Functional Similarity/network/'+ppi_filename_only + '_' + str(threshold).strip()
		output = ppi_filename_only + '_' + rnd +'.net'
		if not os.path.exists(output):
			return output


	return None


def __getRandomNetwork():

	# random network .

	global NETWORK_TEMP_FILES 

	[interaction_file, index1, index2, score_index, threshold] = NETWORK_TEMP_FILES

	ppi_filename_only = os.path.basename(interaction_file) 
	ppi_filename_only = os.path.splitext(ppi_filename_only)[0] + '_' + str(threshold).strip() 

	files, folders = getFileListIn('./6. Functional Similarity/network/')
	prebuilt = []
	for f in files:
		if f.find(ppi_filename_only)>=0:
			prebuilt.append('./6. Functional Similarity/network/' + f)    

	fname = random.sample( prebuilt, 1 ) [0]



	ppi = PPI.PPIDB_STRING()
	ppi.load(fname, 0, 1, 2, -1)
	return ppi

def loadInteraction(interaction_file, iteration, index1, index2, index3, score_index, threshold):


	#global NETWORK_TEMP_FILES

	print ('Loading interactions = ', interaction_file)
	box = []

	ppi = PPI.PPIDB_STRING()
	ppi.load(interaction_file, index1, index2, score_index, threshold, ' ', _0_Preprocess.STRING_ALIAS_FILE)
	box.append(ppi) # normal network

	ppi_filename_only = os.path.basename(interaction_file) 
	ppi_filename_only = os.path.splitext(ppi_filename_only)[0] + '_' + str(threshold).strip() 
	#ppi_filename_only = './6. Functional Similarity/network/'+ppi_filename_only


	'''
    # random networks
    files, folders = MyUtil.getFileListIn('./6. Functional Similarity/network/')
    prebuilt = []
    for f in files:
        if f.find(ppi_filename_only)>=0:
            prebuilt.append('./6. Functional Similarity/network/' + f)

    if len(prebuilt)>iteration:
        prebuilt = random.sample( prebuilt, iteration )


    for i in range(iteration):

        print '\r\t', i+1, '/', iteration,
        p = PPI.PPIDB_STRING()

        if i < len(prebuilt):
            #p.load(prebuilt[i], index1, index2, score_index, -1)  
            #box.append(p)

            box.append(    [ prebuilt[i] , index1, index2, score_index, threshold   ] ) 

        else:
            output = __getOutputFilename(interaction_file, threshold)

            p.interDB = copy.deepcopy( ppi.interDB )
            p.scoreDB = copy.deepcopy( ppi.scoreDB )

            p.randomizeNetwork()
            p.saveInteractionIntoFile(output)



            #box.append(p)

            box.append( [ output, index1, index2, score_index, threshold   ] ) 
            del p

    '''

	return box


def run(iteration, interaction_file, index1, index2, index3, score_index, threshold):


	global RND_ID, PRELOAD_NO_OF_RANDOM_NETWORKS, NETWORK_TEMP_FILES, MAX_THREAD , DEGREE, DISEASE_SET
	global P_VALUE
	
	print ('ID=', RND_ID)

	mods = __getModifiers()

	summary_file = TEMP_PATH + '/NetExp' + RND_ID + '_' + '_summary.txt'
	f=open(summary_file, 'w')
	dd = mods.keys()
	dd.sort()


	NETWORK_TEMP_FILES = [interaction_file, index1, index2, score_index, threshold]  


	info = '--------------------------------------------' + '\n' + \
	        'THREAD = ' + str( MAX_THREAD ) + '\n' + \
	        'Remove temporary files = ' + str( REMOVE_TEMP_FILES ) + '\n' + \
	        'Calculate p-value = ' + str( P_VALUE) + '\n' + \
	        '\t# of random networks = ' + str( PRELOAD_NO_OF_RANDOM_NETWORKS) + '\n' + \
	        'Iteration = ' + str( HOW_MANY_TEST_MODIFERS) + '\n' + \
	        'Randomization method = ' + str(DEGREE) + '\n' + \
	        'r = ' + str(R_PARAMETER)
	print (info)





	ppi_box = loadInteraction(interaction_file, PRELOAD_NO_OF_RANDOM_NETWORKS, index1, index2, index3, score_index, threshold) # 100개 네트워크를 만든다.
	total_genes = __getTotalGeneList(ppi_box)




	print ('init...')
	global Pnp_matrix, Pnp_dict
	total_gene_names = ppi_box[0].getWholeGeneList()
	total_gene_names.sort()     
	local_Pnp_matrix, local_Pnp_dict = __buildPnpMatrix(total_gene_names, ppi_box[0])
	print ('done')





	ppi_box2 = copy.deepcopy(ppi_box)



	auc_output = TEMP_PATH + '/NetExp_' + RND_ID + '_' + '_'.join(dd) + '_auc.txt'


	mod_box = {}


	for key in dd:
		mod_pos = mods[key] 
		new_mod_pos = __removeUnannotatedModifiers(total_genes, mod_pos)
		mod_box[key] = new_mod_pos

	total_genes_without_modifiers = __getTotalGenesWithoutModifiers(total_genes, mod_box)


	s = 'Interaction = ' + interaction_file + '\n' + \
	        'Disease = ' + repr(DISEASE_SET)+ '\n' + \
	        'Total annotated gene # = ' + str(len(total_genes)) + '\n' 

	f.write(info+'\n')
	f.write('-------------------------\n' +s+'\n')

	print (s)

	prank = []

	# 정해진 만큼 train과 test를 반복한다.


	random_or_not = 1 # 1-random. 0-not random  This must be RANDOM to calculate AUC


	for i in range(iteration):

		Pnp_dict = copy.deepcopy(local_Pnp_dict)
		Pnp_matrix = copy.deepcopy(local_Pnp_matrix)


		print ('\r\t', RND_ID, ' Disease = ', repr(dd), ' # = ', i+1)


		pos_rank = mainJob( dd, copy.deepcopy(mod_box), total_genes, total_genes_without_modifiers, i ,ppi_box2, random_or_not)

		print ('\t', pos_rank)

		prank.append(  pos_rank )


	print ('')





	foname = TEMP_PATH + '/'+RND_ID+'_'+'_'.join(DISEASE_SET)+'_ranks.txt'
	__saveRank(prank, foname)



	pos_auc = __calculateAUC(prank,  auc_output)


	s = 'Avg rank = ' + str( MyUtil_pypy.mean( prank  ) )  + ' +- ' + str( MyUtil_pypy.pstdev(prank) ) + '\n' + \
	        'AUC + = ' + str( pos_auc ) + '\n' 

	print (s)
	f.write(s+'\n')
	f.close()


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


def __extractCommonModifiers(total_genes, mod_box):

	r = []

	print ('Common modifiers = ')
	for g in total_genes:
		cnt = 0    
		for key in mod_box.keys():
			mp = mod_box[key]
			if not mp.has_key(g):
				cnt = 1
				break

		if cnt == 0:
			#print g
			r.append(g)

			for key in mod_box.keys():
				mp = mod_box[key]

				#print 'deleting ', g, ', ', len(mp), ' -> ', 

				del mp[g]
				#print mod_box[key][0].has_key(g), len(mp)
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

		cross_set, new_mod_box = __randomSampling (total_genes, cross_set, mod_box, basket, ppi_box[0])

		test_pos, test_others = __sampleModifiers(total_genes, cross_set, new_mod_box)






	else:
		ppi_box2 = ppi_box
		new_mod_box = copy.deepcopy(mod_box)



	print ('TEST gene = ', test_pos, ' out of ', len(cross_set))
	train_output = TEMP_PATH + '/NetExp' + RND_ID+'__trained_'+str(iteration_th).strip()+ '_'.join(disease_set) + '_' + str(iteration_th).strip()+'.txt'


	test_set_all = [test_pos ] + test_others

	__train(total_genes, new_mod_box, train_output, test_set_all, ppi_box2)



	train_output2 = train_output + '_gene.txt'
	r1 = __test( train_output2, test_pos, test_others  )

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
		#deg = len( ppi.getIncomingPartnersOf(g)  )
		deg = len( ppi.getPartners(g)  )
		deg_str = str(deg).strip()

		if not box.has_key(deg_str):
			box[deg_str] = []

		if not g in box[deg_str]:
			box[deg_str].append(g)

	return box

















def __removeUnannotatedModifiers(total_genes, mods):
	r = {}
	for k in mods.keys():
		if k in total_genes and not r.has_key(k):
			r[k] = mods[k]
	return r

def __getTotalGeneList(ppi_box):

	ppi = ppi_box[0]

	total_gene_list = []
	for g in ppi.getWholeGeneList():
		total_gene_list.append(  g )
	return total_gene_list

def __sampleModifiers(total_genes, cross_set, mod_box):



	print ('common set # ', len(cross_set))

	other_genes = []
	test_pos = random.sample( cross_set, 1) [0]



	while ( len(other_genes) < 99 ):
		o = random.sample(total_genes, 1) [0]
		has = False
		for d in mod_box.keys():
			mp = mod_box[d]
			if mp.has_key(o):
				has = True
				break

		if has==False and o != test_pos:
			other_genes.append(o)





	return test_pos, other_genes


def __getTotalGenesWithoutModifiers(total_genes, mod_box ):
	r = []
	for t in total_genes:

		pass_or_fail = False

		for k in mod_box.keys():
			if mod_box[k].has_key(t):
				pass_or_fail = True
				break

		if pass_or_fail == False:
			r.append(t)

	return r


















# =============================================


def __train(total_genes, mod_box, output, test_set_all, ppi_box):





	# 

	pvalues = __calculatePvalue(total_genes, test_set_all, mod_box, output, ppi_box)
	__resave(total_genes, mod_box, pvalues, output, test_set_all)




def __resave(total_genes, mod_box,  pvalues, output, test_set_all):



	ofile = output + '_gene.txt'


	box = __rebox(test_set_all, mod_box, pvalues, output )
	__save2(box, mod_box, ofile)



def __save2(box, mod_box, output):

	global P_VALUE

	order = True

	if P_VALUE:
		order = False
		#
		box3 = []
		for g in box.keys():
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
		box3 = sorted( box.items(), key=lambda (k,v): v[0], reverse = order)



	maps, mods = _0_Preprocess.getModifiers()
	diseases = mods.keys()
	diseases.sort()

	f=open(output, 'w')
	f.write('Gene\tscore\tName\tDesc\t' + '\t'.join(diseases) + '\n')


	for k, v in box3:

		b = []
		for d in diseases:
			if k in mods[d]:
				b.append('o')
			else:
				b.append('_')


		score = v[0]
		desc = v[-1]
		s = k + '\t' + str(score) + '\t' + _0_Preprocess.ID_BOX.getCommonNameOf(k) + '\t' + desc + '\t' + '\t'.join(b)
		f.write(s+'\n')

	f.close()


def __addTag(gene, mod_pos):
	tag = ''
	if gene in mod_pos.has_key:
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

			p_pos, score, txt = pvalues[g]

			desc.append( _0_Preprocess.ID_BOX.getCommonNameOf(g))
			desc.append( '[' + str(score)  +  ' ] ')
			desc.append(txt)

			box[g] = [ p_pos, score, ','.join(desc)]
		else:
			desc = []

			p_pos = pvalues[g]

			desc.append( _0_Preprocess.ID_BOX.getCommonNameOf(g))
			desc.append( '[' + str(p_pos)  +  ' ] ')


			box[g] = [ p_pos, ','.join(desc)]


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

			# Desc
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

	n_avg = float((len(u_partners) - 1 + len(v_partners) - 1 ))/2.0  
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

			# confidence .
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

	return Pnp_dict[gene] # no direction



def __randomSamplingBasedOnReliability(total_genes, test_set_gene, mod_pos, ppi, cross_set):

	#    1. test_set_gene 
	#    2. mod_pos 

	new_mod = []

	resvoir = []
	#for g in total_genes:
	#	if not g in test_set_gene and not g in mod_pos.keys() and not g in cross_set:
	#		resvoir.append(g)

	resvoir = copy.deepcopy(total_genes)


	for m in mod_pos:

		m_reliability = __getReliabilitySum( m , ppi)

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
		
		
		random.shuffle(  resvoir )
		ratio = 0.1  # reliability 값 차이가 원래 reliability의 10% 이내면 받아들인다.
		
		for rg in resvoir :

			if not rg in new_mod:

				rg_reliability = __getReliabilitySum(rg, ppi)
				
				if abs( rg_reliability - m_reliability  ) <= ratio * m_reliability:

					new_mod.append(rg)
					break





	return new_mod


def __randomSamplingBasedOnReliability_old(total_genes, test_set_gene, mod_pos, ppi, cross_set):

	#    1. test_set_gene 
	#    2. mod_pos 

	new_mod = {}

	resvoir = []
	for g in total_genes:
		if not g in test_set_gene and not g in mod_pos.keys() and not g in cross_set:
			resvoir.append(g)

	#resvoir = total_genes


	for m in mod_pos:

		m_reliability = __getReliabilitySum( m , ppi)

		other_min_reliability = 10000000.0 
		for m2 in mod_pos.keys():

			if m2 == m: 
				continue

			m2_reliability = __getReliabilitySum(m2, ppi)
			if abs(m_reliability - m2_reliability)< abs(m_reliability - other_min_reliability):
				other_min_reliability = m2_reliability


		while(True):

			rg = random.sample(  resvoir , 1 ) [0] 

			if not rg in new_mod.keys():

				rg_reliability = __getReliabilitySum(rg, ppi)
				if abs( rg_reliability - m_reliability  ) <= abs( rg_reliability - other_min_reliability  ):

					new_mod[rg] = mod_pos[m]
					break


	return new_mod



def logmsg(msg):




	f=open('r:/netexp.log', 'a')
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

	#return sc


def __calculatePvalue(total_genes, test_set_gene, mod_box, output, ppi_box):

	global P_VALUE, PRELOAD_NO_OF_RANDOM_NETWORKS

	r = -1

	# my algorithm
	real_scores = __calculatePvalueSUB2(test_set_gene, mod_box, output, ppi_box[0], r)

	if P_VALUE == False:
		print ('No p-value calculation: use scores....')
		return real_scores

		'''
        elif P_VALUE == True:

            # r= 0.0, reference

            zero_scores = __calculatePvalueSUB2(test_set_gene, mod_box, output, ppi_box[0], 0.0)

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

	else:    

		rand = {}
		box = []

		

		for i in range(PRELOAD_NO_OF_RANDOM_NETWORKS):


			print ('[p-value] Iteration: ', (i+1), '/', PRELOAD_NO_OF_RANDOM_NETWORKS)

			## Purely random selection
			new_mod_pos = {}
			for k in mod_box.keys():
				#new_mod_pos[k] = random.sample(ttt, len(mod_box[k]))

				gene_list = __randomSampling (total_genes, None, mod_box[k], None, ppi_box[0])
				new_mod_pos[k] = {}
				for g in gene_list:
					new_mod_pos[k][g] = random.random()  # g는 원래 STRING ID임.


			ret = __calculatePvalueSUB2(test_set_gene, new_mod_pos, output, ppi_box[0], r)

			for key in ret.keys():
				if not rand.has_key(key):
					rand[key] = []

				rand[key].append( ret[key]   )

		

		print ('Calculating p-values...')
		
		temp_scores = {}
		scores = {}
		
		if PYPY:
			for g in rand.keys():
				avg = MyUtil_pypy.mean( rand[g])
				std = MyUtil_pypy.pstdev( rand[g])	
				
				temp_scores[g] = [ real_scores[g], avg, std   ]
			
			
			scores = MyUtil_pypy.calculate_pvalues(temp_scores, _0_Preprocess.TEMP_PATH)
		
		else:

			for g in rand.keys():
				avg = MyUtil_pypy.mean( rand[g])
				std = MyUtil_pypy.pstdev( rand[g])
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
					
				scores [g] = p
		
		pvalues = {}
		for g in rand.keys():
			avg = MyUtil_pypy.mean( rand[g])
			std = MyUtil_pypy.pstdev( rand[g])				
			txt = 'Score:'+str(real_scores[g]) + ', avg='+str(avg) + ' std=' + str(std)
			pvalues[g]=[ scores[g], real_scores[g], txt]

		return pvalues # p-value based on z-score


def __getTrim(total_genes, mod_pos, ppi):

	total2 = []
	score = {}
	for g in total_genes:
		score[g] = 0.0 # 

	for g in mod_pos.keys():
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
	if a>b:
		return a+':'+b
	else:
		return b+':'+a



def __getW(u, ppi):

	total = 0.0
	for k in ppi.getPartners(u):
		total = total + ppi.getScore(k, u)
	return total  



def __buildPnpMatrix(total_gene_names, ppi, cache_file = None):

	dic = {}


	method = 2


	if method == 1:


		c = 0.0
		q = len(total_gene_names) 

		for ig1 in range(len(total_gene_names)):

			g1 = total_gene_names[ig1]

			c+=1.0
			print ('\r', (c/q*100.0), '%')

			for ig2 in range(ig1, len(total_gene_names)):

				g2 = total_gene_names[ig2]

				key = __sortStr(g1, g2)

				if not dic.has_key(key):
					score = 1.0
					if g1!=g2:
						score = __calculateFunctionalSimilarity(g1, g2, 1.0, ppi)
					dic[key] = score

					ss = key + '\t' + str(score).strip()
					#f.write(ss+'\n')


	else:






		c = 0.0
		q = float( len(total_gene_names) )
		for u in total_gene_names:

			c+=1.0
			print ('\r', (c/q*100.0), '%')

			dic[u] = __getW(u, ppi)



			'''
            #ppi.getScore(k,p) * ppi.getScore(p, j) / ( __getW(k,ppi) * __getW(p,ppi) * __getW(j, ppi)  ) ** (1.0/3.0)

            partners = ppi.getPartners(u)


            for i in range( len(partners) - 1 ):
                for j in range( i+1, len(partners)):

                    pi = partners[i]
                    pj = partners[j]



                    #s = ppi.getScore(pi, u) * ppi.getScore(pj, u)

                    #s = s / ( __getW(pi,ppi) * __getW(pj, ppi) * __getW(u, ppi) )

                    #s = s / ( __getW(pi,ppi) * __getW(pj, ppi) * __getW(u, ppi) )  ** (1.0/3.0)

                    key = u + ':' + __sortStr(pi, pj)        
                    dic[key] = s
            '''

		if cache_file is not None:
			f=open(cache_file,'w')
			for k in dic.keys():
				s = k+'\t'+str(dic[k]).strip()
				f.write(s+'\n')
			f.close()

	print ('')


	#f.close()


	return None, dic

def runFinal(iteration, interaction_file, index1, index2, index3, score_index, threshold, modI = None, test_set_all = []):



	global P_VALUE
	global RND_ID, PRELOAD_NO_OF_RANDOM_NETWORKS, NETWORK_TEMP_FILES, MAX_THREAD , DEGREE, DISEASE_SET, R_PARAMETER


	mod_maps, mods = _0_Preprocess.getModifiers()

	summary_file = TEMP_PATH + '/NetExp' + RND_ID + '_' + '_summary.txt'
	f=open(summary_file, 'w')
	dd = mods.keys()
	dd.sort()


	NETWORK_TEMP_FILES = [interaction_file, index1, index2, score_index, threshold]  


	info = '--------------------------------------------' + '\n' + \
	        'THREAD = ' + str( MAX_THREAD ) + '\n' + \
	        'Remove temporary files = ' + str( REMOVE_TEMP_FILES ) + '\n' + \
	        'Calculate p-values = ' + str( P_VALUE) + '\n' + \
	        '\t# of random networks = ' + str( PRELOAD_NO_OF_RANDOM_NETWORKS) + '\n' + \
	        'Iteration = ' + str( HOW_MANY_TEST_MODIFERS) + '\n' + \
	        'Randomization method = ' + str(DEGREE) + '\n'
	print (info)





	ppi_box = loadInteraction(interaction_file, PRELOAD_NO_OF_RANDOM_NETWORKS, index1, index2, index3, score_index, threshold) # 100개 네트워크를 만든다.
	total_genes = __getTotalGeneList(ppi_box)


	cache_file = copyCache()


	print ('init...')
	global Pnp_matrix, Pnp_dict
	total_gene_names = ppi_box[0].getWholeGeneList()
	total_gene_names.sort()     
	Pnp_matrix, Pnp_dict = __buildPnpMatrix(total_gene_names, ppi_box[0], cache_file=cache_file)
	print ('done')


	ppi_box2 = []
	ppi_box2 = ppi_box


	mod_box = {}

	for key in dd:
		mod_pos = mods[key] 

		new_mod_pos = __removeUnannotatedModifiers(total_genes, mod_pos)

		mod_box[key] = new_mod_pos









	s = 'Interaction = ' + interaction_file + '\n' + \
	        'Disease = ' + repr(dd)+ '\n' + \
	        'Total annotated gene # = ' + str(len(total_genes)) + '\n' 

	f.write(info+'\n')
	f.write('-------------------------\n' +s+'\n')

	print (s)


	fname = '__Common_modifiers.txt'

	if _0_Preprocess.OUTPUT_FILENAME is not None:
		fname = _0_Preprocess.OUTPUT_FILENAME


	output = _0_Preprocess.OUTPUT_FOLDER + '/' + fname

	__train2(total_genes, mod_box, output, total_genes, ppi_box)



def copyCache():

	fname = '_12_' + os.path.basename (_0_Preprocess.STRING_FILE) + \
	        '_' + str( _0_Preprocess.STRING_CUTOFF).strip() + '.cache.txt'


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


def __train2(total_genes, pvalues, output, test_set_all, ppi_box):
	__resave(total_genes, None, pvalues, output, total_genes)



def __getSeedScore(total_gene_names, mod, ppi):

	r = {}
	for k in total_gene_names:
		r[k] = 0.0

	total = 0.0
	for k in mod.keys():
		total = total + mod[k]

	for k in mod.keys():
		r[k] = mod[k]/total

	arr = []
	for k in total_gene_names:
		arr.append(  [r[k]] )
		#arr.append(  r[k] )

	return arr, r



def __calculateScore(p, k, j, new_box, r, i, ppi):
	global Pnp_matrix, Pnp_dict

	key = p+':'+__sortStr(k, j)
	w = Pnp_dict[key]



	s = r * (new_box[k] * new_box[j]) * w



	return s

def __categorizePartners( mod_box, partners ):

	r_box = {}
	dis = mod_box.keys()
	dis.sort()

	for p in partners:


		for disease in dis:
			mp = mod_box[disease]

			if mp.has_key(p):
				if not r_box.has_key(disease):
					r_box[disease]=[]
				r_box[disease].append(p)

	q_box = []
	cbu = 0
	cbs = ''

	for d in dis:
		if not r_box.has_key(d):
			return None
		else:

			cbu = cbu + len(r_box[d] ) - 1
			if cbs == '':
				cbs = r_box[d][0]
			else:
				if cbs != r_box[d][0]:
					cbu = cbu + 1



	if cbu == 0:
		return None

	return r_box




def __calculatePvalueSUB2(total_genes, mod_box, output, ppi, r):    



	#global R_THRESHOLD

	box = {} 
	#threshold = R_THRESHOLD
	prev_mod_box = {}


	print ('Seed genes normalization') # rank ratio
	for disease in mod_box.keys():
		total = 0.0
		mod_pos = mod_box[disease]

		prev_mod_box[disease] ={}

		print (disease, '-->', len(mod_pos))

		for k in mod_pos:
			total = total + mod_pos[k]

		for k in mod_pos:
			prev_mod_box[disease][k]=mod_pos[k]/total



	prev_cache = {} 


	cnt = 0


	dis = mod_box.keys()
	dis.sort()

	identified_common_modifiers = {}

	while(True):

		cnt = cnt + 1
		if cnt>100:
			break

		new_mod_box = {}
		#new_cache = copy.deepcopy( prev_cache  )




		for g in ppi.getWholeGeneList():

			if  identified_common_modifiers.has_key(g):
				#print g, 'is a generic modifier'
				continue


			partners = ppi.getPartners(g) + [g] 

			''' common modifier '''
			cat_partners = __categorizePartners( prev_mod_box, partners )



			if cat_partners == None:
				#print 'skip ,, ', g
				continue



			sc = 1.0 / (Pnp_dict[g])

			# score = 1/(incombin wet sum) * sigma[ score(p)*reliability(from,to)/(outgoing w sum)  ]
			for disease_name in cat_partners.keys():

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


				#print 'found ', g, sc
				new_mod_box[g] = sc
				#new_cache[g]=sc









		for g in new_mod_box.keys():

			identified_common_modifiers[g] = new_mod_box[g]

			for d in dis:

				prev_mod_box[d][g] = new_mod_box[g]



		print ('step = ', cnt, ' identified at this step = ', len(new_mod_box), ' identified common modifiers so far = ', len(identified_common_modifiers))
		if len(new_mod_box) == 0:
			print ('no more modifiers found')
			break

		#new_mod_box = {}



	for g in total_genes:
		if not identified_common_modifiers.has_key(g):
			identified_common_modifiers[g] = 0.0

	return identified_common_modifiers


def __calculatePvaluSUB(total_genes, mod_box, output, ppi, r):    




	prev_mod_box = {}

	# seed gene normalization
	for disease in mod_box.keys():
		total = 0.0
		mod_pos, mod_neg = mod_box[disease]

		prev_mod_box[disease] ={}

		for k in mod_pos.keys():
			total = total + mod_pos[k]

		for k in mod_pos.keys():
			prev_mod_box[disease][k]=mod_pos[k]/total

	prev_cache = {} #  common modifier 

	#print 'mod # = ', len(mod_pos)

	cnt = 0


	dis = mod_box.keys()
	dis.sort()

	while(True):

		cnt = cnt + 1
		if cnt>100:
			break

		new_mod_box = {}
		new_cache = copy.deepcopy( prev_cache  )




		for g in ppi.getWholeGeneList():

			if  prev_cache.has_key(g):
				continue



			partners = ppi.getPartners(g) + [g]  
			ptx = [] + partners

			for p in partners:
				psh = ppi.getPartners(p)
				for p2 in psh:
					if not p2 in ptx:
						ptx.append(p2)

			partners = ptx


			cat_partners = __categorizePartners( prev_mod_box, partners )



			if cat_partners == None:
				#print 'skip ,, ', g
				continue




			sc = 1.0



			for dindex in range(len(cat_partners)):

				disease_modifiers = cat_partners[dindex]
				disease_name = dis[dindex]

				fscore = 1.0

				for dm in disease_modifiers:


					key = __sortStr(dm, g)
					fscore=fscore*(1.0 - Pnp_dict[key] * prev_mod_box[disease_name][dm])   

				fscore = 1.0 - fscore

				sc = sc  * fscore




			for di in dis:
				if not new_mod_box.has_key(di):
					new_mod_box[di]={}

				new_mod_box[di][g] = sc


			new_cache[g]=sc


		for di in mod_box.keys():
			mp, np = mod_box[di]


			if not new_mod_box.has_key(di):
				new_mod_box[di] = {}

			nn = new_mod_box[di]

			for gk in mp.keys():
				if not nn.has_key(gk):
					nn[gk] = mp[gk]


		print ('step = ', cnt, ' # = ', len(new_cache))


		if len(prev_cache) == len(new_cache):
			prev_mod_box = new_mod_box
			prev_cache = new_cache
			break

		prev_mod_box = new_mod_box
		prev_cache = new_cache



	ret = {}
	for g in total_genes:
		if not new_cache.has_key(g):
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



				gid = x[0].replace('m_', '')

				if gid == test_pos:
					cnt = cnt + 1
					if pos_rank == -1:
						pos_rank = cnt
					else:
						print ('Error!!!! two or more positives')

				elif gid in test_others:
					cnt = cnt + 1


	f.close()


	if cnt != 100:
		print ('Error !!! Not 100 test set', cnt)

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

			p_r = float(p_v) / total_rank  # rank ratio 
			n_r = float(n_v) / total_rank  

			if p_r <= threshold < n_r:
				success = success + 1.0

		sensitivity = success / total_trial

		pos_false_positive_rate_x.append( 1.0 - specificity  )
		pos_true_positive_rate_y.append( sensitivity )

		f.write(  str(1.0-specificity) + '\t' + str(sensitivity) + '\n'  )


	auc = MyUtil_pypy.mean( pos_true_positive_rate_y  )
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


	auc_pos = MyUtil_pypy.mean( pos_true_positive_rate_y  )
	f.close()




	return auc_pos

















def getModifiersRank(fname):

	f=open(fname,'r')

	rank = 0.0

	r = {}

	init = True

	print ('Loading a modifier file', fname, '....')

	for s in f.readlines():
		if init:
			init = False
		else:
			s = s.replace('m_', '').replace('x_', '')

			x = s.split('\t')

			gene_id = x[0]
			pvalue = eval( x[1] ) # meaningless
			string_id = _0_Preprocess.ID_BOX.getSTRING_ID_of(gene_id)

			if string_id is not None:
				rank = rank + 1.0
				r[gene_id] = rank


				#print gene_id, '\t', r[gene_id]

			else:
				print ('\tUnrecognized protein name: ', gene_id)

	f.close()

	return r


def getPredictedModifiers(diseases):


	fnames = {}
	for d in diseases:
		fnames[d] = OUTPUT_FOLDER + '/Integrated_'+d+'.txt'

	r = {}

	for d in fnames.keys():
		v = getModifiersRank(fnames[d] )
		#r[d] = [ v, {} ]
		r[d] = v

		print (d, fnames[d], len(v), 'genes loaded...')

	return r



def FinalizeAll():



	# raqnk product를 기반으로 한다.


	#seed_number = _0_Preprocess.SEED_NUMBER
	seed_number = _0_Preprocess.getSEEDNumber()
	string_db_file = _0_Preprocess.STRING_FILE
	index1 = 0
	index2 = 1
	index3 = 3
	score_index = 2
	threshold = _0_Preprocess.STRING_CUTOFF

	maps, mods = _0_Preprocess.getModifiers()
	diseases = mods.keys()

	seeds = getPredictedModifiers(diseases) # dict, key=protein, value=score(rank order)
	# seeds['disease'] --> dㅕ기에 들어가는 value는, STRING id가 key, value는 rank (정수)


	protein_ids = seeds[diseases[0] ].keys()
	rho = []

	for p in protein_ids:

		rp = 1

		for d in diseases:
			rp = rp * seeds[d][p]

		rho.append(rp)


	n = len(protein_ids)
	k = len(diseases)
	pvalues_0 = rankproduct.calculate_rankproduct_p_value(rho, n, k)

	pvalues = {}
	for p_id, p_val in zip(protein_ids, pvalues_0):
		pvalues[p_id] = [p_val, 0.0, 'None' ]


	ppi_box = loadInteraction(string_db_file, PRELOAD_NO_OF_RANDOM_NETWORKS, index1, index2, index3, score_index, threshold) # 100개 네트워크를 만든다.
	total_genes = __getTotalGeneList(ppi_box)

	fname = '__Common_modifiers.txt'
	if _0_Preprocess.OUTPUT_FILENAME is not None:
		fname = _0_Preprocess.OUTPUT_FILENAME


	output = _0_Preprocess.OUTPUT_FOLDER + '/' + fname

	__train2(total_genes, pvalues, output, total_genes, ppi_box)
	
	print ('######################################################')
	print (' Result file is ', output)
	print ('######################################################')








def FinalizeAll_random():


	#seed_number = _0_Preprocess.SEED_NUMBER
	seed_number = _0_Preprocess.getSEEDNumber()
	string_db_file = _0_Preprocess.STRING_FILE
	index1 = 0
	index2 = 1
	index3 = 3
	score_index = 2
	threshold = _0_Preprocess.STRING_CUTOFF

	maps, mods = _0_Preprocess.getModifiers()
	diseases = mods.keys()










	NETWORK_TEMP_FILES = [string_db_file, index1, index2, score_index, threshold]


	ppi_box = loadInteraction(string_db_file, PRELOAD_NO_OF_RANDOM_NETWORKS, index1, index2, index3, score_index, threshold) # 100개 네트워크를 만든다.
	total_genes = __getTotalGeneList(ppi_box)


	#seeds = getPredictedModifiers(diseases, seed_number) # dict, key=protein, value=score(rank order)

	# randomization -------------------------------------------------------
	seeds = {}
	total = float(len(total_genes))

	for k in diseases:
		seeds[k] = {}

		temp = random.sample(total_genes, seed_number)
		
		print ('[', k, '] Randomly selected = ', len(temp))
		
		for i in range(len(temp)):
			seeds[k][temp[i]] = 1.0 - float(i)/total






	print ('init...')
	global Pnp_matrix, Pnp_dict
	total_gene_names = ppi_box[0].getWholeGeneList()
	total_gene_names.sort()
	#local_Pnp_matrix, local_Pnp_dict = __buildPnpMatrix(total_gene_names, ppi_box[0])
	Pnp_matrix, Pnp_dict = __buildPnpMatrix(total_gene_names, ppi_box[0])
	print ('done')






	ppi_box2 = []
	#ppi_box2 = copy.deepcopy(ppi_box)


	mod_box = {}

	for key in diseases:

		#if key.find('RANDOM') > 0:
		#    continue

		mod_pos = seeds[key]   #    [ v{}, {} ]  # v contains key=gene, value=score --> v{} only

		print ('------------------------')
		print (key)

		#Pnp_matrix, Pnp_dict = copy.deepcopy(local_Pnp_dict), copy.deepcopy(local_Pnp_matrix)

		new_mod_pos = __removeUnannotatedModifiers(total_genes, mod_pos)


		mod_box[key] = new_mod_pos


	total_genes_without_modifiers = __getTotalGenesWithoutModifiers(total_genes, mod_box)


	fname = '__Common_modifiers_random.txt'
	if _0_Preprocess.OUTPUT_FILENAME is not None:
		fname = _0_Preprocess.OUTPUT_FILENAME + '_random.txt'
	#if '[RND]' in fname:
	#	fname = fname.replace('[RND]', RND_ID)

	output = _0_Preprocess.OUTPUT_FOLDER + '/' + fname


	__train2(total_genes, mod_box, output, total_genes, ppi_box)
	
	print ('[OUTPUT] = ', output)


def init(config_file):

	global RND_ID, TEMP_PATH, OUTPUT_FOLDER, P_VALUE

	_0_Preprocess.init(config_file)

	RND_ID = _0_Preprocess.RND_ID
	TEMP_PATH = _0_Preprocess.TEMP_PATH
	OUTPUT_FOLDER = _0_Preprocess.OUTPUT_FOLDER
	P_VALUE = _0_Preprocess.NETWORK_EXPANSION_PVALUE 
	#P_VALUE = _0_Preprocess.P_VALUE

	print ('***** Loaded P-VALUE= ', P_VALUE)


def main(opt):
	
	global PYPY, PRELOAD_NO_OF_RANDOM_NETWORKS, MP_PROCESSORS

	config_file = opt[_0_Preprocess.OPT_CONFIG_FILE]
	MP_PROCESSORS = opt[_0_Preprocess.OPT_MULTIMP]
	rnd_or_not = opt[_0_Preprocess.OPT_RANDOM_OR_NOT]
	test_or_not = opt[_0_Preprocess.OPT_TEST]
	PRELOAD_NO_OF_RANDOM_NETWORKS = opt[_0_Preprocess.OPT_ITERATION]
	PYPY = opt[_0_Preprocess.OPT_PYPY]


	start_time = datetime.now()


	# python _1_GO.py config.txt test
	print ('''
========================================================    
    [12] Rank Products
========================================================          
''')


	if config_file is not None:
		
		init(config_file)
		FinalizeAll()
	
	else:
		init(config_file)


		if test_or_not == True:
			run() # iteration 의미 없음.
		elif rnd_or_not == True:
			FinalizeAll_random()




	print ('''
========================================================    
    [12] Rank Products (End)
========================================================          
''')


	end_time = datetime.now()
	print ('Elapsed time: {}'.format(end_time - start_time))
	print (time.ctime())


if __name__ == '__main__':


	opt = _0_Preprocess.process_args(sys.argv[1:])
	main(opt)
	
	
