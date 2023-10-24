# -*- coding: utf-8 -*-

'''
Randomization을 해서 p-value 계산하는 코드임.
_5_Usearch.py는 그냥 score만 계산하는 코드임.
'''


import sys
import os
import statistics
import _0_Preprocess
import copy
import random
import math
import USEARCH
from datetime import datetime
import time
import MULTIPROCESS
import MyUtil
from itertools import product
from itertools import permutations
from itertools import combinations
import shutil

RND_ID = _0_Preprocess.RND_ID
TEMP_PATH = _0_Preprocess.TEMP_PATH
#TEMP_PATH = '~/netexp/TEMP'
OUTPUT_PATH = _0_Preprocess.OUTPUT_FOLDER
P_VALUE = True
ITERATION = 10000

USE_EVALUE_OR_BIT_SCORE = 2 # 1 - e-value, 2 - bit-score
MP_PROCESSORS = 1


def __getModifiers():
	r, l = _0_Preprocess.getModifiers()
	return r, l


def __randomSampling (total_genes, pos_mod_num):
	
	pos_mod = random.sample( total_genes, pos_mod_num )
	return pos_mod
	

def run(iteration):

	
	
	RND_ID = _0_Preprocess.RND_ID

	mod_maps, mods = __getModifiers()
	
	summary_file = TEMP_PATH + '/USEARCH_' + RND_ID + '_summary.txt'
	f=open(summary_file, 'w')

	dd = list(mods)
	
	dd.sort()

	for disease in dd:
		
		
		auc_output = './TEMP/USearch_'+RND_ID + '_' + disease + '_auc.txt'

		total_genes = _0_Preprocess.ID_BOX.getAllStringIDs()
		
		mod_pos = mods[disease] 
		new_mod_pos = __removeUnannotatedModifiers(total_genes, mod_pos)
		

		if len(new_mod_pos) == 0 :
			print (disease, ' has no recognized genes=',len(new_mod_pos))
			continue

		
		total_genes_without_modifiers = __getTotalGenesWithoutModifiers(total_genes, mod_pos)
		
		
		s = '==================' + '\n' + \
			'Disease = ' + disease + '\n' + \
			'Total annotated gene # = ' + str(len(total_genes)) + '\n' + \
			'Modifiers # = ' + str( len(new_mod_pos) ) 
		
		f.write(s+'\n')  
		
		
		prank = []
		
		for i in range(iteration):
			
			
			print ('\r\t' + RND_ID + ' Disease = ', disease, ' # = ', i+1)
			pos_rank = __mainJob( disease, new_mod_pos, total_genes, total_genes_without_modifiers, i )
			print ('\t', pos_rank)
			prank.append(  pos_rank )
		
		
		print ('')
		
		foname = TEMP_PATH + '/'+RND_ID+'_'+disease+'_ranks.txt'
		__saveRank(prank, foname)
		
		pos_auc = __calculateAUC(prank, auc_output)
		
		
		s = 'Average rank = ' + str( statistics.average( prank  ) )  + ' +- ' + str( statistics.stdev(prank) ) + '\n' + \
			'AUC = ' + str( pos_auc ) + '\n' 
		
		print (s)
		f.write(s+'\n')
	f.close()
	
	print ('Result summary = ', summary_file)
	

		

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

def __mainJob(disease, mod_pos, total_genes, total_genes_without_modifiers, iteration_th):
	
	global TEMP_PATH
	
	new_mod_pos, test_pos, test_others = __sampleModifiers(total_genes_without_modifiers, mod_pos)
	
	if disease.upper().find('RANDOM')>0:
		new_mod_pos = __randomSampling (total_genes, len(new_mod_pos)+1)
		total_genes_without_modifiers = __getTotalGenesWithoutModifiers(total_genes, new_mod_pos)
		while(True):
			x_mod_pos, X_test_pos, x_test_others = __sampleModifiers(total_genes_without_modifiers, new_mod_pos)
			if not (test_pos in x_mod_pos or test_pos in x_test_others):
				new_mod_pos = x_mod_pos
				test_others = x_test_others

				test_pos = X_test_pos

				break

	train_output = TEMP_PATH+ './USEARCH_trained_' + RND_ID+'_'+disease + '_' + str(iteration_th).strip()+'.txt'
	
	test_set_all = [test_pos ] + test_others
 
	#__train(total_genes, new_mod_pos, train_output)
	__train(test_set_all, new_mod_pos, train_output)
	
	train_output2 = train_output + '_gene.txt'
	r1 = __test( train_output2, test_pos, test_others  )
	
	
	try:
		os.remove(train_output)		
	except:
		pass
	
	try:
		os.remove(train_output2)
	except:
		pass
	
	
	return r1



def __removeUnannotatedModifiers(total_genes, mods):
	r = []
	for k in mods:
		if k in total_genes and not k in r:
			r.append(k)
	return r
	

	
def __sampleModifiers(total_genes_without_modifiers, mod_pos):
	
	# randomly select 1 modifier and 99 non-modifiers

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

def preprocess(pathVisio):

	fly_box = {} # key: TF, values: target genes
	fly_genes = {} # key: target gene, values: TFs
	
	
	
	for p in pathVisio.getAllPathwayNames():
		
		pathway_id = p
		genes = pathVisio.getGenesInvolvedIn(p)
		
		for g in genes:
			if not g in fly_genes:
				fly_genes[g] = []
			fly_genes[g].append(pathway_id)

			
		
		fly_box[pathway_id] = genes
	


	return fly_box, fly_genes

def __linkGOandGenes(go_box, genes):
	
	
	for g in genes:
		for goid in genes[g]:
			if goid in go_box:
				if not g in go_box[goid]:
					go_box[goid].append(g)
	
	
	new_go_box = {}
	for goid in go_box:
		if len(go_box[goid])>0:
			new_go_box [goid] = go_box[goid]
	
	return new_go_box
	
	
	
def __split(go_box):
	
	fly = {}
	yeast = {}
	worm = {}
	
	for gid in go_box:
		for g in go_box[gid]:
			if g[0] == 'F':
				if not gid in fly:
					fly[gid] = []
				fly[gid].append(g)
			elif g[0] == 'S':
				if not gid in yeast:
					yeast[gid] = []
				yeast[gid].append(g)
			else:
				if not gid in worm:
					worm[gid] = []
				worm[gid].append(g)
	
	return fly, yeast, worm


def __selectLeafsOnly( go_box, go_class ):
	
	new_go_box = {}
	for gid in go_box:
		go = go_class.getGOTermOf(gid)

		#if len( go.getChildRelationshipGOIDs() )== 0: # leaf node
		if len( go.getChildRelationshipGOIDsOfISA() )== 0 and len( go.getChildRelationshipGOIDsOfPARTOF() ) == 0: # leaf node
		
			# add a leaf node
			new_go_box[gid] = go_box[gid]

	# replace with new one
	print ('GO reduced from ', len(go_box) , ' to ', len(new_go_box))
	return new_go_box



def __loadGeneInfo(flybase, category):
	
	r = {}
	for g in flybase.getAllGeneIDs():
		data = flybase.getDataByUniqueID(g)
		
		gu = g 
		
		for gt in data.getGOTerms():
			
			x = gt.split(';') 
			if x[0].strip() in category:
				if not gu in r:
					r[gu] = []
				
				r[gu].append( x[1].strip() )
			
	return r






def __initializeGO(go_class):
	r = {} # go_id as a key, values are gene id
	go_ids = list(go_class.go_terms)
	for g in go_ids:
		r[g] = []
		#print g
	return r
















# =============================================


def __train(total_genes, mod_pos, output):

	# p value
	#pvalues = __calculatePvalue(total_genes, mod_pos, output)
	pvalues = __calculatePvalue2(total_genes, mod_pos, output)

	__resave(total_genes, mod_pos, pvalues, output, '')
	
	
	

def __resave(total_genes, mod_pos, pvalues, output, md5):
	
	
	
	ofile = output + '_gene.txt'

	box = __rebox(total_genes, mod_pos, pvalues, output)
	__save2(box, mod_pos, ofile, md5)

def __save2(box, mod_pos, output, md5, order = True):

	global USE_EVALUE_OR_BIT_SCORE, P_VALUE

	if P_VALUE:
		order = False

	print('[_5_Usearchy] Saving result file = ' + output)
	

	
	box3 = None





	if P_VALUE == False:
		if USE_EVALUE_OR_BIT_SCORE == 1: # evalue
			box3 = sorted( box.items(), key=lambda item: item[1][0], reverse = False)
		elif USE_EVALUE_OR_BIT_SCORE == 2: # bit-score
			box3 = sorted( box.items(), key=lambda item: item[1][0], reverse = True)
	else:
		box3 = sorted(box.items(), key=lambda item: item[1][0], reverse=order)
	
	
	f=open(output, 'w')
	
	f.write('Gene\tscore\tDesc\n')
	f.write('!'+md5+'\n')
	
	for k, v in box3:
		score = v[0]
		desc = v[1]
		
		s = __addTag(k, mod_pos)+k+'\t'+str(score)+'\t'+desc
		f.write(s+'\n')
		
	
	f.close()


def __addTag(gene, mod_pos):
	tag = ''
	if gene in mod_pos:
		tag = 'm_'
	return tag

def __rebox(total_genes, mod_pos, pvalues, output):

	
	# calculate score !
	box = {}
	for g in total_genes:
		
		p = 0.0
		desc = ''
		
		p_pos = pvalues[g]
		p = p_pos 
			
		desc = _0_Preprocess.ID_BOX.getCommonNameOf(g)
		box[g] = [ p, desc] 


	return box


def __makeGeneList(gene_list, mod_pos):
	r = []
	for g in gene_list:
		if g in mod_pos:
			r.append( 'm_' + g )
		else:
			r.append ( g )
	return r		


def __makeDBFile_with_sequence(db_file, total_genes, mod_pos):
	# making a db for Blast
	f = open(db_file, 'w')

	for m in mod_pos:

		seq = _0_Preprocess.ID_BOX.getSequenceOf(m)

		if seq is not None:
			s = '>' + m + '\n' + seq
			f.write(s + '\n')

	f.close()


def __makeDBFile(db_file, total_genes, mod_pos):

	# making a db for Blast 
	f=open(db_file,'w')

	for m in mod_pos:

		seq = _0_Preprocess.ID_BOX.getSequenceOf(m)
		
		if seq is not None:
			s = '>'+m+'\n'+seq
			f.write(s+'\n')


	f.close()

def __calculatePvalue(total_genes, mod_pos, output):

	evalues = {}
	blast = USEARCH.USEARCH()
	
	RND_ID = _0_Preprocess.RND_ID


	db_file = TEMP_PATH+'/' + RND_ID + '_usearch_db.fasta'
	__makeDBFile(db_file, total_genes, mod_pos)

	new_db_file = blast.makeDB(db_file)
   
	for g in total_genes:
		
		protein_seq = _0_Preprocess.ID_BOX.getSequenceOf(g)

		r = blast.runUsingSequence(protein_seq, new_db_file)


		if USE_EVALUE_OR_BIT_SCORE == 1: # use e-values
			e = 10000000.0

			for gene, bits, evalue in r:
				if evalue < e:
					e = evalue

			evalues[g] = e
		elif USE_EVALUE_OR_BIT_SCORE == 2: # bit-score
			# larger is better
			e = -1000000.0

			for gene, bits, evalue in r:
				if bits > e:
					e = bits

			evalues[g] = e




		#print g, e

	try:
		os.remove(db_file)
		os.remove(new_db_file)
	except:
		pass
		

	return evalues




def __calculatePvalue2(total_genes, mod_pos, output):


	rnd_str = MyUtil.getRandomString(20)
	evalues = {}
	#Blast.DEBUG = True

	blast = USEARCH.USEARCH()
	RND_ID = _0_Preprocess.RND_ID


	db_file = TEMP_PATH+'/usearch_db_' + RND_ID + '_' + rnd_str + '.fasta'
	#db_file = '/home/ssbio/netexp/TEMP/' + RND_ID + '_' + rnd_str + '_usearch_db.fasta'
	
	__makeDBFile(db_file, total_genes, mod_pos)

	new_db_file = blast.makeDB(db_file)
	sequences_as_dict = {}
	trash = []
	for g in total_genes:
		p = _0_Preprocess.ID_BOX.getSequenceOf(g)
		if p is None:
			print (g,'not in UniProt')
			trash.append(g)
			continue
		else:
			if len(p) == 0:
				print (g, 'len = 0')
				trash.append(g)
				continue
			else:
				sequences_as_dict[g] = p
		
	# run blast
	#print ('Running Usearch...')

	ret = blast.runUsingSequences(sequences_as_dict, new_db_file)
	
	#print(ret)
	
	
	for g in total_genes:
		
		if g in trash:
			continue
		
		# no results
		if not g in ret:
			if USE_EVALUE_OR_BIT_SCORE == 1: # use e-values
				evalues[g] = 10000000.0
			elif USE_EVALUE_OR_BIT_SCORE == 2: # bit-score
				evalues[g] = -10000000.0
			continue
			
		# result found
			
		r = ret[g]

		if USE_EVALUE_OR_BIT_SCORE == 1: # use e-values
			e = 10000000.0

			for gene, bits, evalue in r:
				if evalue < e:
					e = evalue

			evalues[g] = e
		elif USE_EVALUE_OR_BIT_SCORE == 2: # bit-score
			# larger is better
			e = -1000000.0

			for gene, bits, evalue in r:
				if bits > e:
					e = bits

			evalues[g] = e




		#print g, e

	#try:
	#	os.remove(db_file)
	#	os.remove(new_db_file)
	#except:
	#	pass
	remove(db_file)
	remove(new_db_file)



	return evalues


def __calculatePvalue3(cache_scores, total_genes, mod_pos, output):

	new_scores = {}

	for g1 in total_genes:

		score = 0
		if USE_EVALUE_OR_BIT_SCORE == 1:  # use e-values
			score = 10000000.0
		elif USE_EVALUE_OR_BIT_SCORE == 2:
			score = -110000000.0

		for g2 in mod_pos:

			try:
				key = __sort(g1, g2)
				e = cache_scores[key]

				if USE_EVALUE_OR_BIT_SCORE == 1:
					if e < score:
						score = e
				elif USE_EVALUE_OR_BIT_SCORE == 2:
					if e > score:
						score = e
			except:
				pass

		new_scores[g1] = score

	return new_scores


def remove(fname):

	try:
		os.remove(fname)
	except:
		pass


def __countModifierNumber(mod_list, gene_list):
	cnt = 0
	for g in gene_list:
		if g in mod_list:
			cnt = cnt + 1
	return cnt




def __test( train_output2, test_pos,  test_others   ):

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
					cnt = cnt + 1
					if pos_rank == -1:
						pos_rank = cnt
					else:
						print ('Error!!!! two or more positives')
				elif gid in test_others:
					cnt = cnt + 1





	f.close()


	if cnt != 100:
		print ('Error !!! Not in the 100 test set', cnt)

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
				if not g in r:
					r[g] = [ [], [] ]
				
				#r[g][0].append( eval(score))
				r[g][0].append( eval(pvalue))
				r[g][1].append( goid )
				
			
	
	
	f.close()

	return r

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
			r = float(v) / total_rank  
			if r <= threshold:
				success = success + 1

		sensitivity = success / total_trial

		pos_false_positive_rate_x.append( 1.0 - specificity  )
		pos_true_positive_rate_y.append( sensitivity )

		f.write(  str(1.0-specificity) + '\t' + str(sensitivity) + '\n'  )


	auc_pos = statistics.average( pos_true_positive_rate_y  )
	f.close()


	return auc_pos




def calculate_mpi(args):

	global TEMP_PATH, RND_ID

	cache_file, total_genes, mod_pos, output, per, queue = args

	cache_scores = MyUtil.loadVariableFromPickle(cache_file)

	box = {}

	for i in range(per):
		rnd_mods = random.sample(total_genes, len(mod_pos))
		#sc = __calculatePvalue2(total_genes, rnd_mods, output)
		sc = __calculatePvalue3(cache_scores, total_genes, rnd_mods, output)
		for g in sc:
			if not g in box:
				box[g] = []
			box[g].append(sc[g])

	# save
	r_file = TEMP_PATH + '/pickle_Usearch_RET_' + MyUtil.getRandomString(20) + '_' + RND_ID + '.obj'
	MyUtil.saveVariableAsPickle(box, filename=r_file)
	#print('[InterPro] ', time.ctime(), ' result file was saved into ', r_file)
	queue.put(r_file)






def __train2_pvalue(cache_scores, cache_file, total_genes, mod_pos, output, md5):

	global MP_PROCESSORS, ITERATION

	# p value
	#pvalues = __calculatePvalue(total_genes, mod_pos, output)
	#print('*****************************')
	print('[_5_USEARCH2.py] *** USearch basic calculation....')
	#print('*****************************')

	scores = __calculatePvalue3(cache_scores, total_genes, mod_pos, output)

	#__resave(total_genes, mod_pos, pvalues, output, md5)

	# randomization해서 p-value 구한다.
	box = {}

	time_out = _0_Preprocess.get_default_opt()[_0_Preprocess.OPT_TIMEOUT]

	if MP_PROCESSORS <= 1:

		#print('*****************************')
		print('[_5_USEARCH2.py] *** P-value calculation (single process)....')
		#print('*****************************')

		for i in range(ITERATION):

			# print('------> interprot iteration: ', i)

			rnd_mods = random.sample(total_genes, len(mod_pos))
			ps = __calculatePvalue3(cache_scores, total_genes, rnd_mods, output)
			for g in ps:
				if not g in box:
					box[g] = []
				box[g].append(ps[g])

	else:

		# multi-process
		#print('*****************************')
		print('[_5_USEARCH2.py] *** P-value calculation (multiprocess)....')
		#print('*****************************')

		thread_box = []
		#steps = int(MP_PROCESSORS * 1.1)  # cpu 10개. 혹시 몰라 10% 계산 추가함.
		steps = int(MP_PROCESSORS)  # cpu 10개. 혹시 몰라 10% 계산 추가함.
		#if steps == MP_PROCESSORS:    steps = MP_PROCESSORS + 1
		per = int(ITERATION / steps)  # 1/10 씩 나눠서 계산한다.
		for i in range(steps):
			args = [ cache_file, total_genes, mod_pos, output, per]
			th = MULTIPROCESS.JOB_THREAD()
			th.set_args(calculate_mpi, args)
			th.set_timeout(time_out)
			thread_box.append(th)

		ret_files = MULTIPROCESS.runMultiprocesses(thread_box, max_cpu=steps, sleep_interval_seconds=10,
		                                           title='[USearch calculation]', cutoff=MP_PROCESSORS)

		for r_file in ret_files:
			ps = MyUtil.loadVariableFromPickle(r_file, remove_after_load=True)
			for g in ps:
				if not g in box:
					box[g] = []
				box[g] += ps[g]

	new_scores = {}
	for g in total_genes:
		mean = statistics.average(box[g])
		std = statistics.stdev(box[g])
		if std == 0:
			new_scores[g] = 1
			continue

		p = statistics.getNormalDistPvalue(mean, std, scores[g])
		if p == -1:
			p = 1
		new_scores[g] = p

	# save results
	#__resave(interpro, interpro_box, interpro_box_rev, total_genes, mod_pos, pvalues, scores, output, test_set_all, md5)
	__resave(total_genes, mod_pos, new_scores, output, md5)

def __train2(cache_scores, cache_file, total_genes, mod_pos, output, md5):

	# p value
	#pvalues = __calculatePvalue(total_genes, mod_pos, output)
	pvalues = __calculatePvalue2(total_genes, mod_pos, output)

	__resave(total_genes, mod_pos, pvalues, output, md5)	


def __sort(a, b):
	if a > b:
		return a + ':' + b
	else:
		return b + ':' + a




def cache_processing_mpi(total_genes):

	global MP_PROCESSORS

	print('-------------------------')
	print(' USearch caching... ')
	print(' CPU = ', MP_PROCESSORS)
	print('-------------------------')


	# pairwise 조합을 모두 만든다.
	print('\tGenerating gene pairs....')
	combi = list( combinations(total_genes, 2) )

	print('\tDiviging gene pairs in groups....')
	# cpu 숫자만큼 list를 나눈다.

	per = int( len(combi) / (MP_PROCESSORS * 100 ) )
	#per = int(len(combi) / (MP_PROCESSORS * 10000))
	divided_genes = MyUtil.divideListByElementNumber(combi, per)


	# 테스트용
	# -- 실전에서는 여기 지워야 함.
	#divided_genes = divided_genes[ : 10]
	# ==========

	print('\tPreparing MPI runs....  jobs per CPU=', per)

	thread_box = []
	output_filenames = []

	total = len(divided_genes)
	cnt = 0


	string_alias_file = _0_Preprocess.STRING_ALIAS_FILE
	string_sequence_file = _0_Preprocess.STRING_SEQUENCE_FILE


	for gene_pairs in divided_genes:

		cnt += 1
		print('\tCreating threads = ', round(cnt / total * 100.0, 2), '%    ', end='\r')

		output = TEMP_PATH + '/usearch_output_temp_' + MyUtil.getRandomString(20) + '.pickle'


		gene_pairs_file = TEMP_PATH + '/usearch_gene_pairs_' + MyUtil.getRandomString(30)

		MyUtil.saveVariableAsPickle(gene_pairs, filename = gene_pairs_file)


		args = [gene_pairs_file, string_alias_file, string_sequence_file, output]
		th = MULTIPROCESS.JOB_THREAD_NO_RETURN()
		th.set_args( calculate_cache_score_mpi, args  )
		thread_box.append(th)
		output_filenames.append(output)


	MULTIPROCESS.runMultiprocesses_no_return(thread_box, max_cpu=MP_PROCESSORS, sleep_interval_seconds=10,
		                               title='[USearch cache score calculation]')



	box = {}
	for output in output_filenames:
		# score 값을 모두 box에 합침.
		var = MyUtil.loadVariableFromPickle(output, remove_after_load=True)
		box.update(var)


	return box



def cache_processing_mpi2(total_genes):

	global MP_PROCESSORS

	print('-------------------------')
	print(' USearch caching... ')
	print(' CPU = ', MP_PROCESSORS)
	print('-------------------------')

	thread_box = []
	output_filenames = []



	indexes = [i for i in range(len(total_genes)-1)]
	per = int( len(indexes) / (MP_PROCESSORS * 2) )
	#per = int(len(combi) / (MP_PROCESSORS * 10000))
	divided_indexes = MyUtil.divideListByElementNumber(indexes, per)
	#divided_genes = MyUtil.divideList(indexes,  )

	string_alias_file = _0_Preprocess.STRING_ALIAS_FILE
	string_sequence_file = _0_Preprocess.STRING_SEQUENCE_FILE

	total = len(divided_indexes)
	cnt = 0

	total_genes_file = TEMP_PATH + '/usearch_total_genes_' + MyUtil.getRandomString(30)
	MyUtil.saveVariableAsPickle(total_genes, filename=total_genes_file)

	for group in reversed(divided_indexes):

		cnt += 1
		print('\tCreating threads = ', round(cnt / total * 100.0, 2), '%    ', end='\r')

		output = TEMP_PATH + '/usearch_output_temp_' + MyUtil.getRandomString(20) + '.pickle'
		index_file = TEMP_PATH + '/usearch_indexes_' + MyUtil.getRandomString(30)

		MyUtil.saveVariableAsPickle(group, filename = index_file)


		args = [index_file, total_genes_file, string_alias_file, string_sequence_file, output]
		th = MULTIPROCESS.JOB_THREAD_NO_RETURN()
		th.set_args( calculate_cache_score_mpi2, args  )
		thread_box.append(th)
		output_filenames.append(output)


	MULTIPROCESS.runMultiprocesses_no_return(thread_box, max_cpu=MP_PROCESSORS, sleep_interval_seconds=10,
		                               title='[USearch cache score calculation]')


	remove(total_genes_file)

	box = {}
	for output in output_filenames:
		# score 값을 모두 box에 합침.
		var = MyUtil.loadVariableFromPickle(output, remove_after_load=True)
		box.update(var)


	return box



def cache_processing_mpi3(total_genes):

	global MP_PROCESSORS

	print('-------------------------')
	print(' USearch caching... ')
	print(' CPU = ', MP_PROCESSORS)
	print('-------------------------')

	thread_box = []
	output_filenames = []




	string_alias_file = _0_Preprocess.STRING_ALIAS_FILE
	string_sequence_file = _0_Preprocess.STRING_SEQUENCE_FILE


	print('Making a db of total genes')


	cache_db_file = '_cache_db_for_usearch_' + os.path.basename(string_sequence_file) + '.db'

	if os.path.exists('./'+cache_db_file):
		print('** use cached db file...')
		new_db_file = TEMP_PATH + '/' + cache_db_file
		shutil.copyfile('./'+cache_db_file, new_db_file)
	else:
		db_fasta_file = TEMP_PATH + '/usearch_db_' + MyUtil.getRandomString(30)
		__makeDBFile(db_fasta_file, None, total_genes)
		print('Making a db file for USEARCH')
		blast  = USEARCH.USEARCH()
		new_db_file = blast.makeDB(db_fasta_file)
		# 만약을 대비해 파일을 현재 폴더로 이동해 놓는다.
		shutil.copyfile(new_db_file, './' + cache_db_file)

		remove(db_fasta_file)


	total = len(total_genes)
	cnt = 0

	time_out = 60 # 15분 안에 종료되면 재실행한다.

	#groups = MyUtil.divideListByElementNumber(total_genes, 100)  # 100개씩


	for gene in total_genes:

		cnt += 1
		print('\tCreating threads = ', round(cnt / total * 100.0, 2), '%    ', end='\r')

		if _0_Preprocess.ID_BOX.getSequenceOf(gene) is not None:

			output = TEMP_PATH + '/usearch_output_temp_' + MyUtil.getRandomString(20) + '.pickle'

			args = [gene, _0_Preprocess.ID_BOX.getSequenceOf(gene), new_db_file, output]

			th = MULTIPROCESS.JOB_THREAD_NO_RETURN()
			th.set_args( calculate_cache_score_mpi3, args  )
			th.set_timeout(time_out)

			thread_box.append(th)
			output_filenames.append(output)


	MULTIPROCESS.runMultiprocesses_no_return(thread_box, max_cpu=MP_PROCESSORS, sleep_interval_seconds=10,
		                               title='[USearch cache score calculation]')


	remove(new_db_file)

	fx=open('usearch_failure.log','w')

	box = {}
	for i in range(len(output_filenames)):

		gene = total_genes[i]
		output = output_filenames[i]

		if os.path.exists(output):
			# score 값을 모두 box에 합침.
			var = MyUtil.loadVariableFromPickle(output, remove_after_load=True)
			box.update(var)
		else:
			print('[USearch caching error gene] = ', gene)
			fx.write(gene+'\n')
	fx.close()

	return box



def cache_processing_mpi4(total_genes):

	global MP_PROCESSORS

	print('-------------------------')
	print(' USearch caching... ')
	print(' CPU = ', MP_PROCESSORS)
	print('-------------------------')

	thread_box = []
	output_filenames = []




	string_alias_file = _0_Preprocess.STRING_ALIAS_FILE
	string_sequence_file = _0_Preprocess.STRING_SEQUENCE_FILE


	print('Making a db of total genes')


	cache_db_file = '_cache_db_for_usearch_' + os.path.basename(string_sequence_file) + '.db'

	if os.path.exists('./'+cache_db_file):
		print('** use cached db file...')
		new_db_file = TEMP_PATH + '/' + cache_db_file
		shutil.copyfile('./'+cache_db_file, new_db_file)
	else:
		db_fasta_file = TEMP_PATH + '/usearch_db_' + MyUtil.getRandomString(30)
		__makeDBFile(db_fasta_file, None, total_genes)
		print('Making a db file for USEARCH')
		blast  = USEARCH.USEARCH()
		new_db_file = blast.makeDB(db_fasta_file)
		# 만약을 대비해 파일을 현재 폴더로 이동해 놓는다.
		shutil.copyfile(new_db_file, './' + cache_db_file)

		remove(db_fasta_file)


	total = len(total_genes)
	cnt = 0

	time_out = 60 # 15분 안에 종료되면 재실행한다.

	groups = MyUtil.divideListByElementNumber(total_genes, 500)  # 100개씩

	for genes in groups:

		cnt += 1
		print('\tCreating threads = ', round(cnt / total * 100.0, 2), '%    ', end='\r')


		output = TEMP_PATH + '/usearch_output_temp_' + MyUtil.getRandomString(20) + '.pickle'
		genes_file = TEMP_PATH + '/pickle_InterPro_' + MyUtil.getRandomString(20) + '_' + RND_ID + '.obj'
		genes_file = MyUtil.saveVariableAsPickle(genes, filename=genes_file)

		args = [genes_file, string_alias_file,string_sequence_file, new_db_file, output]

		th = MULTIPROCESS.JOB_THREAD_NO_RETURN()
		th.set_args( calculate_cache_score_mpi4, args  )
		th.set_timeout(time_out)

		thread_box.append(th)
		output_filenames.append(output)


	MULTIPROCESS.runMultiprocesses_no_return(thread_box, max_cpu=MP_PROCESSORS, sleep_interval_seconds=10,
		                               title='[USearch cache score calculation]')


	remove(new_db_file)

	fx=open('usearch_failure.log','w')

	box = {}
	for i in range(len(output_filenames)):

		gene = total_genes[i]
		output = output_filenames[i]

		if os.path.exists(output):
			# score 값을 모두 box에 합침.
			var = MyUtil.loadVariableFromPickle(output, remove_after_load=True)
			box.update(var)
		else:
			print('[USearch caching error gene] = ', gene)
			fx.write(gene+'\n')
	fx.close()

	return box


def calculate_cache_score_mpi(args):

	import ID

	gene_pair_file,  string_alias_file, string_sequence_file , output = args

	gene_pairs = MyUtil.loadVariableFromPickle(gene_pair_file)
	#sequences = MyUtil.loadVariableFromPickle(sequence_file, remove_after_load=True)

	id_box = ID.ID_CLASS(string_alias_file, string_sequence_file)


	scores = {}
	for g1, g2 in gene_pairs:

		seq1 = id_box.getSequenceOf(g1)
		seq2 = id_box.getSequenceOf(g2)

		# Usearch로 bit score 계산한다.
		key = __sort(g1, g2)
		scores[key] = get_a_score_of_two_sequences(g1, g2, seq1, seq2)

	# 파일로 저장한다.
	MyUtil.saveVariableAsPickle(scores, filename=output)

	# 리턴은 따로 안 한다.
	remove(gene_pair_file)


def calculate_cache_score_mpi2(args):

	import ID
	global USE_EVALUE_OR_BIT_SCORE

	index_file, total_genes_file, string_alias_file, string_sequence_file, output = args

	total_genes = MyUtil.loadVariableFromPickle(total_genes_file)
	indexes = MyUtil.loadVariableFromPickle(index_file, remove_after_load=True)

	id_box = ID.ID_CLASS(string_alias_file, string_sequence_file)


	scores = {}

	for i in indexes:

		blast = USEARCH.USEARCH()
		rnd_str = MyUtil.getRandomString(30)

		gene = total_genes[i]

		print('making a db fasta')
		db_file = TEMP_PATH + '/usearch_db_' + RND_ID + '_' + rnd_str + '.fasta'
		for g in total_genes[i+1:]:
			f = open(db_file, 'w')
			s = '>' + g + '\n' + id_box.getSequenceOf(g)
			f.write(s + '\n')
			f.close()

		new_db_file = blast.makeDB(db_file)

		#print('usearching...')
		#sequence2 = _0_Preprocess.ID_BOX.getSequenceOf(gene2)
		seq_dict = { gene: id_box.getSequenceOf(gene) }
		ret = blast.runUsingSequences(seq_dict, new_db_file)

		#print('removing...')
		remove(db_file)
		remove(new_db_file)

		if gene in ret:

			for gene_t, bits, evalue in ret[gene]:

				score = 0

				if USE_EVALUE_OR_BIT_SCORE == 1:  # use e-values
					score = evalue
				elif USE_EVALUE_OR_BIT_SCORE == 2:  # bit-score
					score = bits

				key = __sort(gene, gene_t)
				scores[key] = score

	# 파일로 저장한다.
	MyUtil.saveVariableAsPickle(scores, filename=output)

	# 리턴은 따로 안 한다.




def calculate_cache_score_mpi3(args):

	global USE_EVALUE_OR_BIT_SCORE

	gene, sequence, new_db_file, output = args

	scores = {}



	blast = USEARCH.USEARCH()
	seq_dict = { gene: sequence }
	ret = blast.runUsingSequences(seq_dict, new_db_file)

	if gene in ret:

		for gene_t, bits, evalue in ret[gene]:

			score = 0

			if USE_EVALUE_OR_BIT_SCORE == 1:  # use e-values
				score = evalue
			elif USE_EVALUE_OR_BIT_SCORE == 2:  # bit-score
				score = bits


			key = __sort(gene, gene_t)
			scores[key] = score

	# 파일로 저장한다.
	MyUtil.saveVariableAsPickle(scores, filename=output)

	# 리턴은 따로 안 한다.




def calculate_cache_score_mpi4(args):

	import ID
	global USE_EVALUE_OR_BIT_SCORE

	genes_file, string_alias_file,string_sequence_file, new_db_file, output = args

	scores = {}

	genes = MyUtil.loadVariableFromPickle(genes_file, remove_after_load=True)
	id_box = ID.ID_CLASS( string_alias_file, string_sequence_file  )

	blast = USEARCH.USEARCH()
	seq_dict = {}
	for g in genes:
		seq_dict [g] = id_box.getSequenceOf(g)

	ret = blast.runUsingSequences(seq_dict, new_db_file)

	for gene in genes:

		if gene in ret:

			for gene_t, bits, evalue in ret[gene]:

				score = 0

				if USE_EVALUE_OR_BIT_SCORE == 1:  # use e-values
					score = evalue
				elif USE_EVALUE_OR_BIT_SCORE == 2:  # bit-score
					score = bits


				key = __sort(gene, gene_t)
				scores[key] = score

				#print(key, score)

	# 파일로 저장한다.
	MyUtil.saveVariableAsPickle(scores, filename=output)

	# 리턴은 따로 안 한다.


def cache_processing(total_genes):

	global MP_PROCESSORS

	print('-------------------------')
	print(' USearch caching... ')
	print('-------------------------')

	box = {}


	total = float(len(total_genes) ** 2 ) / 2.0
	show = 0
	cnt = 0


	# pairwise 조합을 모두 만든다.
	#combi = combinations(total_genes, 2)
	# cpu 숫자만큼 list를 나눈다.
	#divided_genes = MyUtil.divideList(combi, MP_PROCESSORS)



	for i in range(len(total_genes)-1):

		for j in range(i+1, len(total_genes)):

			cnt += 1
			show += 1

			#if show == 100:
			print('Usearch cache calculation = ', round(cnt / total * 100.0, 2), '%    ', end='\r')
			show = 0

			g1 = _0_Preprocess.ID_BOX.getSTRING_ID_of( total_genes [i] )
			g2 = _0_Preprocess.ID_BOX.getSTRING_ID_of (total_genes[j] )

			# Usearch로 bit score 계산한다.
			key = __sort(g1, g2)
			try:
				# 에러가 없으면 이미 계산한 것이므로 건너 뜀.
				x = box[key]
				continue
			except:
				pass

			box[key] = get_a_score_of_two_sequences(g1, g2)
			print("Score = ", box[key])

	return box

def get_a_score_of_two_sequences(gene1, gene2, sequence1, sequence2):

	global USE_EVALUE_OR_BIT_SCORE

	blast = USEARCH.USEARCH()
	rnd_str = MyUtil.getRandomString(30)

	#print('making a db fasta')
	db_file = TEMP_PATH + '/usearch_db_' + RND_ID + '_' + rnd_str + '.fasta'
	f = open(db_file, 'w')
	s = '>' + gene1 + '\n' + sequence1
	f.write(s + '\n')
	f.close()


	#__makeDBFile(db_file, None, [gene1] ) # gene1을 DB로 만든다.
	new_db_file = blast.makeDB(db_file)

	#print('usearching...')
	#sequence2 = _0_Preprocess.ID_BOX.getSequenceOf(gene2)
	seq_dict = { gene2: sequence2 }
	ret = blast.runUsingSequences(seq_dict, new_db_file)

	#print('removing...')
	remove(db_file)
	remove(new_db_file)

	score = 0
	if USE_EVALUE_OR_BIT_SCORE == 1:
		score = 10000000
	elif USE_EVALUE_OR_BIT_SCORE == 2:
		score = -10000000

	if not gene2 in ret:
		return score

	else:

		for gene, bits, evalue in ret[gene2]:
			if USE_EVALUE_OR_BIT_SCORE == 1:  # use e-values
				if evalue < score:
					score = evalue
			elif USE_EVALUE_OR_BIT_SCORE == 2:  # bit-score
				if bits > score:
					score = bits

		return score


def runFinal():
	
	global OUTPUT_PATH, RND_ID, TEMP_PATH, P_VALUE, ITERATION
	
	print ('Usearch...')
	rnd_str = MyUtil.getRandomString(30)

	mod_maps, mods = _0_Preprocess.getModifiers()
	# mod_maps={}, key=user id, value=string id
	# mods = [], list of converted string IDs

	# 여기에 리스트 형태로 들어가 있는 질병에 대해서만 disease-specific modifiers를 예측한다.
	predict_disease_specific_modifiers = _0_Preprocess.PREDICT_DISEASE_MODIFIERS
	if len(predict_disease_specific_modifiers) == 0:
		# 예측할게 하나도 없으면 그냥 건너뛴다.
		print("[_5_USearch.py] 예측할 disease-specific modifiers가 없음 -> 건너뜀")
		return


	total_genes = _0_Preprocess.ID_BOX.getAllStringIDs()

	# =============================================================
	# cache file
	ppi_file = os.path.basename(_0_Preprocess.PPI_FILES[0][0])
	cache_file = './_cache.usearch_'+ppi_file+'.pickle'
	cache_scores = {}
	if os.path.exists(cache_file):
		cache_scores = MyUtil.loadVariableFromPickle(cache_file)
		print('Loaded caches # = ', len(cache_scores), cache_file)

	else:
		# caching
		# microarray 처럼 gene1:gene2 가 key
		cache_scores = cache_processing_mpi4(total_genes)
		MyUtil.saveVariableAsPickle( cache_scores, cache_file)
	# =============================================================



	
	for disease in mods:

		if not disease in predict_disease_specific_modifiers:
			print('[_5_USearch.py] disease-specific modifiers 예측 안하고 건너뜀: ', disease)
			continue

		if disease.upper().find('RANDOM') > 0:
			mods[disease] = random.sample(total_genes, len(mods[disease]))

		#output = OUTPUT_PATH + '/USearch_' + disease+'.txt'
		output = TEMP_PATH + '/USearch_' + disease + '_' + RND_ID + '.txt'
		#output = '/home/ssbio/netexp/TEMP/USearch_' + disease + '_' + RND_ID + '.txt'
		
		
		print (disease, output)

		md5 = MyUtil.getRandomString(30)
		'''
		# md5 within the output file
		x = sorted( copy.deepcopy( mods[disease] ) )
		s = '\n'.join( x )
		md5 = _0_Preprocess.generateMD5(s)
		print ('Modifiers MD5=', md5)
		if md5 == _0_Preprocess.getMD5(output+'_gene.txt') and _0_Preprocess.REFRESH_RESULTS == False:
			print ('\tUse previous results: ', output+'_gene.txt')
			#print 'MD5 = ', md5[disease]
	
			continue
		'''
		
		mod_pos = mods[disease]
		
		new_mod_pos = __removeUnannotatedModifiers(total_genes, mod_pos)
		
		if P_VALUE:
			#__train2_pvalue(cache_scores, cache_file, total_genes, new_mod_pos, output, md5)
			__train2_pvalue(cache_scores, cache_file, total_genes, new_mod_pos, output, md5)
		else:
			__train2(cache_scores, cache_file, total_genes, new_mod_pos, output, md5)




def init(config_file):
	_0_Preprocess.init(config_file)


	global RND_ID, TEMP_PATH, OUTPUT_PATH

	RND_ID = _0_Preprocess.RND_ID
	TEMP_PATH = _0_Preprocess.TEMP_PATH
	OUTPUT_PATH = _0_Preprocess.OUTPUT_FOLDER


def runFinal_mp(args):
	
	global OUTPUT_PATH, RND_ID, TEMP_PATH
	rnd_str = MyUtil.getRandomString(20)
	
	config_file, disease, opt, queue = args
	init(config_file)
	
	RND_ID = opt[_0_Preprocess.OPT_RND_ID]


	print ('[USEARCH] Init ...', disease)



	# 여기에 리스트 형태로 들어가 있는 질병에 대해서만 disease-specific modifiers를 예측한다.
	predict_disease_specific_modifiers = _0_Preprocess.PREDICT_DISEASE_MODIFIERS
	#if len(predict_disease_specific_modifiers) == 0:
	#	# 예측할게 하나도 없으면 그냥 건너뛴다.
	#	print("[_5_USearch.py (mp)] 예측할 disease-specific modifiers가 없음 -> 건너뜀")
	#	return

	if not disease in predict_disease_specific_modifiers:
		print('[_5_USearch.py (mp)] disease-specific modifiers 예측 안하고 건너뜀: ', disease)
		return


	mod_maps, mods = _0_Preprocess.getModifiers()
	# mod_maps={}, key=user id, value=string id
	# mods = [], list of converted string IDs


	total_genes = _0_Preprocess.ID_BOX.getAllStringIDs()
		
	if disease.upper().find('RANDOM') > 0:
		mods[disease] = random.sample(total_genes, len(mods[disease]))

	#output = OUTPUT_PATH + '/Usearch_' + disease+'_' + RND_ID + '.txt'
	#output = '/home/ssbio/netexp/TEMP/USearch_' + disease + '_' + RND_ID + '.txt'
	output = TEMP_PATH + '/USearch_' + disease+'_' + RND_ID + '.txt'
	
	
	print (disease, output)
	# md5 within the output file
	x = sorted( copy.deepcopy( mods[disease] ) )
	s = '\n'.join( x )
	md5 = _0_Preprocess.generateMD5(s)
	print ('Modifiers MD5=', md5)
	
	
	if md5 == _0_Preprocess.getMD5(output+'_gene.txt') and _0_Preprocess.REFRESH_RESULTS == False:
		print ('\tUse previous results: ', output+'_gene.txt')
		#print 'MD5 = ', md5[disease]

	else:		
	
		mod_pos = mods[disease]
		new_mod_pos = __removeUnannotatedModifiers(total_genes, mod_pos)
		__train2(total_genes, new_mod_pos, output, md5)
		
	queue.put(output)
	
	return output
	
	

def main(opt):
	
	global MP_PROCESSORS, ITERATION, RND_ID, P_VALUE


	config_file = opt[_0_Preprocess.OPT_CONFIG_FILE]
	MP_PROCESSORS = opt[_0_Preprocess.OPT_MULTIMP]
	rnd_or_not = opt[_0_Preprocess.OPT_RANDOM_OR_NOT]
	test_or_not = opt[_0_Preprocess.OPT_TEST]
	ITERATION = opt[_0_Preprocess.OPT_ITERATION]

	# USearch의 경우 P_VALUE 계산이 너무 오래걸림.
	P_VALUE = _0_Preprocess.P_VALUE
	

	init(config_file)
	
	RND_ID = opt[_0_Preprocess.OPT_RND_ID]
	_0_Preprocess.RND_ID = RND_ID
		


	start_time = datetime.now()

	print ('''
========================================================	
	[6] USearch
========================================================		  
''')

	
	if test_or_not == False:
		# one per disease
		
		runFinal()
		
		'''
		if MP_PROCESSORS == 1:

			# 실제로는 이게 실행됨.
			runFinal()



		else:
			mod_maps, mods = __getModifiers()
			diseases = list(mods)
			
			jobs = []
			
			for disease in diseases:


				th = MULTIPROCESS.JOB_THREAD()
				th.set_args(runFinal_mp, [config_file, disease, opt])
				jobs.append(th)
				
				# 그냥 하나씩 돌리자.
				th.run()
				
			#MULTIPROCESS.runMultiprocesses(jobs, max_cpu=MP_PROCESSORS, sleep_interval_seconds=10, title = '[_5_Usearch]')
			#MULTIPROCESS.runMultiprocesses(jobs, max_cpu=1, sleep_interval_seconds=10, title = '[_5_Usearch]')
			
		'''
			
			
	else:
		# support multi-processsors later
		run(ITERATION)
	

	print ('''
========================================================	
	[6] USearch (End)
========================================================		  
''')

	end_time = datetime.now()
	print ('Elapsed time: {}'.format(end_time - start_time))
	print (time.ctime())
	
	
	

if __name__ == '__main__':


	opt = _0_Preprocess.process_args(sys.argv[1:])
	main(opt)
	


	#global MP_PROCESSORS

	#start_time = datetime.now()

	#print '''
#========================================================	
	#[5] Blast
#========================================================		  
#'''

	#import sys
	#if len(sys.argv) == 2:
		#config_file = sys.argv[1]
		#init(config_file)
		
		## one per disease
		
		#MP_PROCESSORS = 1
		
		
		#if MP_PROCESSORS == 1:
			#runFinal()
		#else:
			#mod_maps, mods = __getModifiers()
			#diseases = mods.keys()
			
			#jobs = []
			#import MULTIPROCESS
			#for disease in diseases:
				#th = MULTIPROCESS.JOB_THREAD()
				#th.set_args(runFinal_mp, [config_file, disease])
				#jobs.append(th)
				
			#MULTIPROCESS.runMultiprocesses(jobs, max_cpu=MP_PROCESSORS, sleep_interval_seconds=10)
			
			
			
	#elif len(sys.argv) == 3:
		#if sys.argv[2].lower() == 'test':
			#config_file = sys.argv[1]
			#init(config_file)
			#iteration = 100
			#run(iteration)
	

	#print '''
#========================================================	
	#[5] Blast (End)
#========================================================		  
#'''

	#end_time = datetime.now()
	#print ('Elapsed time: {}'.format(end_time - start_time))
	#print time.ctime()
