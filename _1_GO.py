# -*- coding: ms949 -*-
import random
import os
import statistics
import GO
import _0_Preprocess
import copy
import math
import numpy
import MyUtil
import MULTIPROCESS
from datetime import datetime
import time





# files from config.txt
GO_FILE = _0_Preprocess.GO_FILE
ANNOTATION_FILE = _0_Preprocess.GO_ANNOTATION_FILE
TEMP_PATH = _0_Preprocess.TEMP_PATH
OUTPUT_PATH = _0_Preprocess.OUTPUT_FOLDER
RND_ID = _0_Preprocess.RND_ID
GENE_LIST = _0_Preprocess.TOTAL_GENE_FILE_TO_CALCULATE
P_VALUE = True
ITERATION = 100
CPU = 1

USE_LEAFNODES_ONLY = False


def __getModifiers():
	r, l = _0_Preprocess.getModifiers()
	return r, l



def __randomSampling (total_genes, pos_mod_num):

	pos_mod = random.sample( total_genes, pos_mod_num )
	return pos_mod


def run(iteration,  category ):


	global ITERATION, RND_ID
	ITERATION = iteration


	'''
	category = ['P', 'F', 'C' ]  
	'''

	global GO_FILE, ANNOTATION_FILE

	go_file = GO_FILE
	annotation_file = ANNOTATION_FILE


	mod_maps, mods = __getModifiers()

	summary_file = TEMP_PATH + '/GO_' + repr(category) + _0_Preprocess.RND_ID + '_summary.txt'
	f=open(summary_file, 'w')
	dd = list(mods)
	dd.sort()
	
	
	go1, go_box1, go_box_rev1 = preprocess(go_file, annotation_file, category)
	
	

	for disease in dd:


		go = copy.deepcopy(go1)
		go_box = copy.deepcopy(go_box1)
		go_box_rev = copy.deepcopy(go_box_rev1)
		

		auc_output = TEMP_PATH + '/GO_'+ repr(category) + RND_ID + '_' + disease + '_auc.txt'


		# remove unannotated proteins
		total_genes = __getTotalGeneList(go_box)

		mod_pos = mods[disease] 


		new_mod_pos = __removeUnannotatedModifiers(total_genes, mod_pos)

		# get a list of proteins except modifiers
		total_genes_without_modifiers = __getTotalGenesWithoutModifiers(total_genes, mod_pos)

		s = '==================' + '\n' + \
		        'Disease = ' + disease + '\n' + \
		        'Total annotated gene # = ' + str(len(total_genes)) + '\n' + \
		        'Modifiers # after filtering  = ' + str( len(new_mod_pos) ) + '\n' + \
		        '# of GO terms = ' + str( len(go_box) )

		f.write(s+'\n')

		print (s)

		prank = []

		for i in range(iteration):

			print ('\r\t', repr(category),  RND_ID + ' Disease = ', disease, ' # = ', i+1)

			pos_rank = __mainJob( disease, go, go_box, go_box_rev, new_mod_pos, total_genes, total_genes_without_modifiers, i , category)

			print ('\t', pos_rank)

			prank.append(  pos_rank )

		print ('')

		foname = TEMP_PATH + '/GO_'+ repr(category)+RND_ID+'_'+disease+'_ranks.txt'
		__saveRank(prank, foname)



		pos_auc = __calculateAUC(prank, auc_output)


		s = 'Average rank = ' + str( numpy.mean( prank  ) )  + ' +- ' + str( numpy.array(prank).std() ) + '\n' + \
		        'AUC = ' + str( pos_auc ) + '\n' 

		print (s)
		f.write(s+'\n')
	f.close()
	
	print ('Result summary = ', summary_file)
	
	print ('done')



def run_mpi(iteration, category):

	# 아직 작업안함.


	global ITERATION, RND_ID
	ITERATION = iteration

	'''
	category = ['P', 'F', 'C' ]  
	'''

	global GO_FILE, ANNOTATION_FILE

	go_file = GO_FILE
	annotation_file = ANNOTATION_FILE

	mod_maps, mods = __getModifiers()

	summary_file = TEMP_PATH + '/GO_' + repr(category) + _0_Preprocess.RND_ID + '_summary.txt'
	f = open(summary_file, 'w')
	dd = list(mods)
	dd.sort()

	go1, go_box1, go_box_rev1 = preprocess(go_file, annotation_file, category)

	for disease in dd:

		go = copy.deepcopy(go1)
		go_box = copy.deepcopy(go_box1)
		go_box_rev = copy.deepcopy(go_box_rev1)

		auc_output = TEMP_PATH + '/GO_' + repr(category) + RND_ID + '_' + disease + '_auc.txt'

		# remove unannotated proteins
		total_genes = __getTotalGeneList(go_box)

		mod_pos = mods[disease]

		new_mod_pos = __removeUnannotatedModifiers(total_genes, mod_pos)

		# get a list of proteins except modifiers
		total_genes_without_modifiers = __getTotalGenesWithoutModifiers(total_genes, mod_pos)

		s = '==================' + '\n' + \
		    'Disease = ' + disease + '\n' + \
		    'Total annotated gene # = ' + str(len(total_genes)) + '\n' + \
		    'Modifiers # after filtering  = ' + str(len(new_mod_pos)) + '\n' + \
		    '# of GO terms = ' + str(len(go_box))

		f.write(s + '\n')

		print(s)

		prank = []

		for i in range(iteration):
			print('\r\t', repr(category), RND_ID + ' Disease = ', disease, ' # = ', i + 1)

			pos_rank = __mainJob(disease, go, go_box, go_box_rev, new_mod_pos, total_genes,
			                     total_genes_without_modifiers, i, category)

			print('\t', pos_rank)

			prank.append(pos_rank)

		print('')

		foname = TEMP_PATH + '/GO_' + repr(category) + RND_ID + '_' + disease + '_ranks.txt'
		__saveRank(prank, foname)

		pos_auc = __calculateAUC(prank, auc_output)

		s = 'Average rank = ' + str(numpy.mean(prank)) + ' +- ' + str(numpy.array(prank).std()) + '\n' + \
		    'AUC = ' + str(pos_auc) + '\n'

		print(s)
		f.write(s + '\n')
	f.close()

	print('Result summary = ', summary_file)

	print('done')


def __saveRank(prank, output):

	f=open(output,'w')
	f.write('Pos rank\n')

	for index in range( len(prank)):
		s = str(prank[index])
		f.write(s+'\n')
	f.close()




def __mainJob(disease, go, go_box, go_box_rev, mod_pos, total_genes, total_genes_without_modifiers, iteration_th, category):

	global TEMP_PATH

	new_mod_pos, test_pos, test_others = __sampleModifiers(total_genes_without_modifiers, mod_pos)

	# randomly generate a list of modifiers
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


	train_output = TEMP_PATH + '/GO__trained_' + repr(category)+ '_'+ disease + '_' + str(iteration_th).strip()+'.txt'


	test_set_all = [test_pos ] + test_others

	__train(go, go_box, go_box_rev, total_genes, new_mod_pos, train_output, test_set_all)




	train_output2 = train_output + '_gene.txt'
	r1 = __test( train_output2, test_pos,  test_others  )



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
	# remove proteins with no annotations
	r = []
	for k in mods:
		if k in total_genes and not k in r:
			r.append(k)
	return r

def __getTotalGeneList(go_box):

	total_genes = []

	for goid in go_box:
		for gid in go_box[goid]:
			if not gid in total_genes:
				total_genes.append(gid) 
	return total_genes


def __sampleModifiers(total_genes_without_modifiers, mod_pos):


	# randomly select 1 modifiers and 99 other genes
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



def preprocess(go_file, annotation_file, category):

	go = GO.GO(obo_file=go_file)
	genes = __loadAnnotationFile(annotation_file, category)
	go_box = __initializeGO(go)

	# select only leaf nodes
	if USE_LEAFNODES_ONLY:
		go_box = __selectLeavesOnly( go_box, go )

	# remove terms with no genes
	go_box, go_box_rev = __linkGOandGenes( go_box, genes )
	# key=GO ID, value=annotated genes


	avg, std = __getAvgGeneNumberPerOneGO(go_box)
	print ('avg gene # per one GO term = ',  avg, '+-', std)

	return go, go_box, go_box_rev



def __getAvgGeneNumberPerOneGO(go_box):

	r = []
	for goid in go_box:
		r.append(   len(go_box[goid]))

	avg = numpy.mean(r)
	std = numpy.array(r).std()

	return avg, std




def __linkGOandGenes(go_box, genes):

	# genes ==> GO terms
	for g in genes:
		for goid in genes[g]:
			if goid in go_box:
				if  (g in go_box[goid]) == False:
					go_box[goid].append(g)

	# remove gene IDs without GO terms

	new_go_box = {}
	for goid in go_box:
		if len(go_box[goid])>0:
			new_go_box [goid] = go_box[goid]
			
	go_box_rev = {}
	for goid in new_go_box:
		for g in new_go_box[goid]:
			if not g in go_box_rev:
				go_box_rev[g] = []
			go_box_rev[g].append(goid)

	return new_go_box, go_box_rev





def __selectLeavesOnly( go_box, go_class ):

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



def __loadAnnotationFile(annotation_file, category):

	r = {}

	id_box = _0_Preprocess.ID_BOX

	f = open(annotation_file, 'r')

	for s in f.readlines():
		s = s.strip()
		if len(s) == 0: continue
		if s[0] == '!' : continue

		t = s.split('\t')

		gid = t[1].strip()
		go_id = t[4].strip()
		cat = t[8].strip()

		if cat in category:
			string_id = id_box.getSTRING_ID_of(gid)
			if string_id is not None:
				# add a GO term
				if not string_id in r:
					r[string_id] = []
				if not go_id in r[string_id]:
					r[string_id].append(go_id)





	f.close()


	return r




def __initializeGO(go_class):
	r = {} # go_id as a key, values are gene id
	go_ids = list(go_class.go_terms)
	for g in go_ids:
		r[g] = []
		#print g
	return r
















# =============================================

def __calculateScores(go_box, go_box_rev, total_genes, pvalues):


	# calculate score !
	# 

	# ==================================================================================================
	box = {}


	'''
	for g in total_genes: 

		score = 0.0
		desc = []
		for goid in go_box.keys():
			if g in go_box[goid]:
				p_pos = pvalues[goid]
				desc.append( goid + '[' + str(p_pos)  +  ' ] (' + go.getGOTermOf(goid).getGOName() + ')')

				score = score - math.log10(p_pos) 

		box[g] = [ score, ','.join(desc)] # score, description
	'''

	method = 1
	for g in total_genes: 

		score = 0.0
		desc = []
		for goid in go_box_rev[g]:
			
			p_pos = pvalues[goid]

			if method == 1:
				score = score - math.log10(p_pos) 
			elif method == 2:			
				p = - math.log10(p_pos)
				if score < p:
					score = p
			elif method == 3:
				# auc 제일 안 좋음.
				score = score - math.log10(p_pos)

		if method == 3:
			# norm
			#
			score = score / float(len(go_box_rev[g]))
			
		box[g] = score  # score
		
	return box

def __train(go, go_box, go_box_rev, total_genes, mod_pos, output, test_set_all):

	global P_VALUE, ITERATION


	# calculate p-values
	pos_pvalues = __calculatePvalue(go_box, total_genes, mod_pos, output)
	scores = __calculateScores(go_box, go_box_rev, total_genes, pos_pvalues)

	if P_VALUE == True:
		box = {}
		for i in range(ITERATION):
			rnd_mods = random.sample(total_genes, len(mod_pos))
			ps = __calculatePvalue(go_box, total_genes, rnd_mods, output)
			sc = __calculateScores(go_box, go_box_rev, total_genes, ps)
			for g in sc:
				if not g in box:
					box[g] = []
				box[g].append(sc[g])
			
		new_scores = {}
		for g in total_genes:
			mean = numpy.mean(box[g])
			std = numpy.array(box[g]).std()
			if std == 0:
				new_scores[g] = 1
				continue
			
			p = statistics.getNormalDistPvalue(mean, std, scores[g])
			if p == -1:
				p = 1
			new_scores[g] = p
		
		scores = new_scores

	#__save(go, go_box, mod_pos, mod_neg, pvalues, output)

	__resave(go, go_box, go_box_rev, total_genes, mod_pos, pos_pvalues, scores, output, test_set_all, '')




def __resave(go, go_box, go_box_rev, total_genes, mod_pos, pvalues, scores, output, test_set_all, md5):


	global P_VALUE
	ofile = output + '_gene.txt'

	# calculate modifier scores
	box = __rebox(go, go_box, go_box_rev, total_genes, mod_pos, pvalues, scores, output )
	
	order = True
	if P_VALUE:
		order = False
		
	__save2(box, mod_pos, ofile, md5, order)


def __save2(box, mod_pos, output, md5, order = True):

	box3 = sorted( box.items(), key=lambda item: item[1][0], reverse = order)

	f=open(output, 'w')
	f.write('!' + md5 + '\n')
	f.write('Gene\tscore\tDesc\n')


	for k, v in box3:
		score = v[0]
		desc = v[1]

		s = __addTag(k, mod_pos)+k+'\t'+str(score)+'\t'+desc
		f.write(s+'\n')

	'''
	for g in _0_Preprocess.ID_BOX.getAllStringIDs():
		if not box.has_key(g):
			s = __addTag(g, mod_pos)+g+'\t'+str(-1)+'\t'+ 'Missing in the DB'
			f.write(s+'\n')
	'''
	
	f.close()


def __addTag(gene, mod_pos):
	tag = ''
	if gene in mod_pos:
		tag = 'm_'

	return tag



def __rebox(go, go_box, go_box_rev, total_genes, mod_pos, pvalues, scores, output):


	# calculate score !
	# 

	# ==================================================================================================
	box = {}

	for g in total_genes: 

		score = scores[g]
		desc = []
		
		for goid in go_box_rev[g]:
			
			p_pos = pvalues[goid]
			desc.append( goid + '[' + str(p_pos)  +  ' ] (' + go.getGOTermOf(goid).getGOName() + ')')
			
		box[g] = [ score, ','.join(desc)] # score, description
		
	return box


def __save(go, go_box, mod_pos, pvalues, output):

	f=open(output,'w')

	mod_pos_dict = {}
	for m in mod_pos:
		mod_pos_dict[m] = 0

	m_v = pvalues


	s = 'GO_ID\tp-value\t#gene\t#modifier\tGO Desc\tList'
	f.write(s+'\n')

	# sort by p-value
	vx = sorted( go_box.items(), key=lambda item: m_v[ item[0] ], reverse=False)


	for gid, genes in vx:
		genes = go_box[gid]

		# Desc
		go_t = go.getGOTermOf(gid)
		name = go_t.getGOName()

		# gene list
		lst = ','.join( __makeGeneList(genes, mod_pos) )

		tv = len(genes)

		mv = __countModifierNumber(mod_pos_dict, genes)

		pv = m_v[gid]

		s = gid + '\t' + str(pv) + '\t' + str(tv) + '\t' + str(mv) + '\t' + name + '\t' + lst
		f.write(s + '\n')

	f.close()


def __makeGeneList(gene_list, mod_pos):
	r = []
	for g in gene_list:
		if g in mod_pos:
			r.append( 'm_' + g )
		else:
			r.append ( g )
	return r        



def __calculatePvalue(go_box, total_genes, mod_pos, output):


	
	pr = float(len(mod_pos)) / float(len(total_genes))
	p = {}
	
	mod_pos_dict = {}
	for m in mod_pos:
		mod_pos_dict[m] = 0

	

	for goid in go_box:
		# count the number of genes annotated with goid
		mod_number = __countModifierNumber(mod_pos_dict, go_box[goid])
		total_number = len( go_box[goid] )


		# hypergeometric distribution

		sample_pop = len(mod_pos)
		sample_succ = mod_number

		pop_pop = len(total_genes)
		pop_succ = total_number

		pv = statistics.getHypergeomPvalue(sample_succ, sample_pop, pop_succ, pop_pop  )
		if pv == -1:
			pv = 1

		p [goid] = pv

	return p



def __countModifierNumber(mod_dict, gene_list):
	cnt = 0
	for g in gene_list:
		if g in mod_dict:
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



def __calculateAUC(p_rank, output):

	f = open(output+'_pos.txt','w')
	f.write('False-positive rate(1-specificity)\tTrue-positive rate(sensitivity)\n')

	start_from = 0.0
	start_to = 1.0
	step = 0.01

	pos_false_positive_rate_x = [] # 1- specificity
	pos_true_positive_rate_y = [] # sensitivity

	# positive 계산함.
	threshold = start_from 
	while(threshold <= start_to  ):
		threshold = threshold + step

		# threshold 값 자체가 specificity를 의미한다.
		specificity = 1.0 - threshold

		total_trial = float( len( p_rank  ))
		total_rank = 100.0

		success = 0.0
		for v in p_rank:
			r = float(v) / total_rank  # rank ratio 값이 나오고 바로 threshold랑 비교 가능하다.
			if r <= threshold:
				success = success + 1

		# 다 했다.
		sensitivity = success / total_trial

		pos_false_positive_rate_x.append( 1.0 - specificity  )
		pos_true_positive_rate_y.append( sensitivity )

		f.write(  str(1.0-specificity) + '\t' + str(sensitivity) + '\n'  )


	# auc 값은 true positive 값을 다 더하면 되는 구나.
	auc_pos = numpy.mean( pos_true_positive_rate_y  )
	f.close()


	return auc_pos


def calculate_mpi(args):

	global TEMP_PATH, RND_ID

	go, go_box, go_box_rev, total_genes, mod_pos, output, per, TEMP_PATH, RND_ID, queue = args

	box = {}

	for i in range(per):
		rnd_mods = random.sample(total_genes, len(mod_pos))
		ps = __calculatePvalue(go_box, total_genes, rnd_mods, output)
		sc = __calculateScores(go_box, go_box_rev, total_genes, ps)

		for g in sc:
			if not g in box:
				box[g] = []
			box[g].append(sc[g])

	# save
	temp_file = TEMP_PATH + '/pickle_GO_'+MyUtil.getRandomString(20)+'_'+RND_ID+'.obj'
	temp_file = MyUtil.saveVariableAsPickle(box, filename=temp_file)
	#print('[GO] ', time.ctime(), ' result file was saved into ', r_file)
	queue.put(temp_file)


def __train2(go, go_box, go_box_rev, total_genes, mod_pos, output, test_set_all, md5):

	global P_VALUE, ITERATION, CPU
	
	# p value

	pvalues = __calculatePvalue(go_box, total_genes, mod_pos, output)
	scores = __calculateScores(go_box, go_box_rev, total_genes, pvalues)



	# save p-values of GO-terms (accesory data)
	__save(go, go_box, mod_pos, pvalues, output)


	time_out = _0_Preprocess.get_default_opt( ) [ _0_Preprocess.OPT_TIMEOUT ]

	if P_VALUE == True:

		print ('[GO ITERATION] ', ITERATION)


		box = {}


		# 여길 수정하면 됨.



		if CPU <= 1:
			# single process
			for i in range(ITERATION):
				rnd_mods = random.sample(total_genes, len(mod_pos))
				ps = __calculatePvalue(go_box, total_genes, rnd_mods, output)
				sc = __calculateScores(go_box, go_box_rev, total_genes, ps)
				for g in sc:
					if not g in box:
						box[g] = []
					box[g].append(sc[g])
		else:
			# multi-process

			thread_box = []
			#steps = int(CPU * 1.1) # cpu 10개. 혹시 몰라 10% 계산 추가함.
			steps = int(CPU )  # cpu 10개. 혹시 몰라 10% 계산 추가함.
			#if steps == CPU:
			#	steps = CPU + 1

			per = int(ITERATION / steps) # 1/10 씩 나눠서 계산한다.
			for i in range(steps):
				args = [go, go_box, go_box_rev, total_genes, mod_pos, output, per, TEMP_PATH, RND_ID]
				th = MULTIPROCESS.JOB_THREAD()
				th.set_args(calculate_mpi, args)
				th.set_timeout(time_out) # 20분 안에 안 끝나면 재실행
				thread_box.append(th)

			ret_files = MULTIPROCESS.runMultiprocesses(thread_box, max_cpu = steps, sleep_interval_seconds=10,
			                                   title='[GO calculation]', cutoff=CPU )

			for r_file in ret_files:
				ps = MyUtil.loadVariableFromPickle(r_file, remove_after_load=True)
				for g in ps:
					if not g in box:
						box[g] = []
					box[g] += ps[g]


		new_scores = {}
		for g in total_genes:
			mean = numpy.mean(box[g])
			std = numpy.array(box[g]).std()
			if std == 0:
				new_scores[g] = 1
				continue
			
			p = statistics.getNormalDistPvalue(mean, std, scores[g])
			if p == -1:
				p = 1
			new_scores[g] = p
		
		scores = new_scores


	# save results
	__resave(go, go_box, go_box_rev, total_genes, mod_pos, pvalues, scores, output, test_set_all, md5)
	

def runFinal( category ):
	
	global ITERATION, RND_ID
	global GO_FILE, ANNOTATION_FILE

	# 여기에 리스트 형태로 들어가 있는 질병에 대해서만 disease-specific modifiers를 예측한다.
	predict_disease_specific_modifiers = _0_Preprocess.PREDICT_DISEASE_MODIFIERS
	if len(predict_disease_specific_modifiers) == 0:
		# 예측할게 하나도 없으면 그냥 건너뛴다.
		print("[_1_GO.py] 예측할 disease-specific modifiers가 없음 -> 건너뜀")
		return

	ofiles = []

	print ('GO finalizing...')

	#category = [ 'C', 'F', 'P' ]



	go_file = GO_FILE
	annotation_file = ANNOTATION_FILE


	mod_maps, mods = __getModifiers()
	# mod_maps={}, key=user id, value=string id
	# mods = [], list of converted string IDs

	go, go_box, go_box_rev = preprocess(go_file, annotation_file, category)
	total_genes = __getTotalGeneList(go_box)

	test_set_genes = None

	for disease in mods:

		if not disease in predict_disease_specific_modifiers:
			print('[_1_GO.py] disease-specific modifiers 예측 안하고 건너뜀: ', disease)
			continue



		# run with random genes
		if disease.find('RANDOM') > 0:
			mods[disease] = random.sample(total_genes, len(mods[disease]))


		
		#output = OUTPUT_PATH +'/GO_' + disease+'('+'-'.join(category)+')_' + str(ITERATION) + '_' + RND_ID + '.txt'
		output = TEMP_PATH + '/GO_' + disease + '(' + '-'.join(category) + ')_' + str(ITERATION) + '_' + RND_ID + '.txt'
		
		# md5 within the output file
		x = sorted( copy.deepcopy( mods[disease] ) )
		x.append(_0_Preprocess.GO_ANNOTATION_FILE)
		x.append(_0_Preprocess.GO_FILE)
		x.append( str(ITERATION) )

		'''
		s = '\n'.join( x )
		md5 = _0_Preprocess.generateMD5(s)
		print ('Modifiers MD5=', md5)
		if md5 == _0_Preprocess.getMD5(output+'_gene.txt') and _0_Preprocess.REFRESH_RESULTS == False:
			print ('\tUse previous results: ', output+'_gene.txt')
			#print 'MD5 = ', md5[disease]
			continue
		'''

		md5 = MyUtil.getRandomString(40)
		

		print (disease, output)

		mod_pos = mods[disease]

		new_mod_pos = __removeUnannotatedModifiers(total_genes, mod_pos)


		__train2(go, go_box, go_box_rev, total_genes, new_mod_pos, output, test_set_genes, md5)
		
		ofiles.append(output)

	print ('done')
	
	return ofiles

def remove(fname):
	try:
		os.remove(fname)
	except:
		pass

def init(config_file, iteration_for_pvalue):


	global GO_FILE, ANNOTATION_FILE, TEMP_PATH
	global OUTPUT_PATH, RND_ID, GENE_LIST, P_VALUE, ITERATION

	_0_Preprocess.init(config_file)

	GO_FILE = _0_Preprocess.GO_FILE
	ANNOTATION_FILE = _0_Preprocess.GO_ANNOTATION_FILE
	TEMP_PATH = _0_Preprocess.TEMP_PATH
	OUTPUT_PATH = _0_Preprocess.OUTPUT_FOLDER
	RND_ID = _0_Preprocess.RND_ID
	GENE_LIST = _0_Preprocess.TOTAL_GENE_FILE_TO_CALCULATE
	
	#
	P_VALUE = _0_Preprocess.P_VALUE

	ITERATION = iteration_for_pvalue



if __name__ == '__main__':



	start_time = datetime.now()
	

	# python _1_GO.py config.txt test
	print (''')
========================================================    
    [1] GO
========================================================          
''')

	import sys
	if len(sys.argv) == 3:

		config_file = sys.argv[1]
		iteration = eval(sys.argv[2])

		init(config_file, iteration)

		#runFinal(['C'])
		#runFinal(['P'])
		#runFinal(['F'])
		
		runFinal(['C', 'F', 'P'])
		
		

	elif len(sys.argv) == 4:

		if sys.argv[3].lower() == 'test':

			config_file = sys.argv[1]
			iter = eval(sys.argv[2])
			iteration = 100 # test #

			init(config_file, iter)

			#iteration = 100
			#run(iteration, ['C'])
			#run(iteration, ['P'])
			#run(iteration, ['F'])
			
			run(iteration, ['C', 'F', 'P'])
	


	print ('''
========================================================    
    [1] GO (End)
========================================================          
''')


	end_time = datetime.now()
	print ('Elapsed time: {}'.format(end_time - start_time))
	print (time.ctime())