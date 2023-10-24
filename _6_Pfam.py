# -*- coding: ms949 -*-
import os
import statistics
import _0_Preprocess
import copy
import random
import math
from datetime import datetime
import time




# files from config.txt
PFAM_FILE = _0_Preprocess.PFAM_FILE
TEMP_PATH = _0_Preprocess.TEMP_PATH
OUTPUT_PATH = _0_Preprocess.OUTPUT_FOLDER
RND_ID = _0_Preprocess.RND_ID

P_VALUE = True
ITERATION = 100

def __getModifiers():
	r, l = _0_Preprocess.getModifiers()
	return r, l


def __randomSampling (total_genes, pos_mod_num):

	pos_mod = random.sample( total_genes, pos_mod_num )
	return pos_mod


def run(iteration):


	

	tf_file = _0_Preprocess.PFAM_FILE
	
	RND_ID = _0_Preprocess.RND_ID

	mod_maps, mods = __getModifiers()

	summary_file = TEMP_PATH + '/PFAM_' + _0_Preprocess.RND_ID + '_summary.txt'
	f=open(summary_file, 'w')
	dd = mods.keys()
	dd.sort()
	
	
	tf_box1, tf_box_rev1 = preprocess(tf_file)
	
	

	for disease in dd:


		
		box = copy.deepcopy(tf_box1)
		box_rev = copy.deepcopy(tf_box_rev1)

		auc_output = TEMP_PATH + '/Pfam_'+_0_Preprocess.RND_ID + '_' + disease + '_auc.txt'



		total_genes = __getTotalGeneList(box)
		mod_pos = mods[disease] 
		new_mod_pos = __removeUnannotatedModifiers(total_genes, mod_pos)
		total_genes_without_modifiers = __getTotalGenesWithoutModifiers(total_genes, mod_pos)

		s = '==================' + '\n' + \
		        'Disease = ' + disease + '\n' + \
		        'Total annotated gene # = ' + str(len(total_genes)) + '\n' + \
		        'Modifiers # = ' + str( len(new_mod_pos) ) + '\n' + \
		        '# of INTERPRO terms = ' + str( len(box) )


		f.write(s+'\n')


		if len(new_mod_pos) == 0:
			print ('===> NO AVAILABLE PROTEINS: quit')
			continue

		print (s)

		prank = []

		for i in range(iteration):

			print ('\r\t' + RND_ID + ' Disease = ', disease, ' # = ', i+1,)

			pos_rank = __mainJob( disease, None, box, box_rev, new_mod_pos, total_genes, total_genes_without_modifiers, i )

			print ('\t', pos_rank,)

			prank.append(  pos_rank )

		print ('')

		foname = TEMP_PATH + '/TF_'+RND_ID+'_'+disease+'_ranks.txt'
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




def __mainJob(disease, none, box, box_rev, mod_pos, total_genes, total_genes_without_modifiers, iteration_th):

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



	train_output = TEMP_PATH + '/GO__trained_'+disease + '_' + str(iteration_th).strip()+'.txt'


	test_set_all = [test_pos ] + test_others

	__train(None, box, box_rev, total_genes, new_mod_pos, train_output, test_set_all)



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
	r = []
	for k in mods:
		if k in total_genes and not k in r:
			r.append(k)
	return r

def __getTotalGeneList(box):

	total_genes = []

	for goid in box.keys(): 
		for gid in box[goid]:
			if not gid in total_genes:
				total_genes.append(gid) 
	return total_genes


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

def preprocess(tf_file):



	tf_box = __loadAnnotationFile(tf_file)
	
	tf_box_rev = {}
	for tf in tf_box.keys():
		for g in tf_box[tf]:
			if not g in tf_box_rev:
				tf_box_rev[g] = []
			tf_box_rev[g].append(tf)
	
	

	avg, std = __getAvgGeneNumberPerOneGO(tf_box)
	print ('avg gene # per one TF = ',  avg, '+-', std)

	return tf_box, tf_box_rev

def __getAvgGeneNumberPerOneGO(go_box):

	r = []
	for goid in go_box:
		r.append(   len(go_box[goid]))

	avg = statistics.average(r)
	std = statistics.stdev(r)

	return avg, std

def __linkInterProandGenes(annotations):

	
	interpro_box = {}
	
	for g in annotations.keys():
		for inter in annotations[g]:
			if not inter in interpro_box:
				interpro_box[inter] = []
			if not g in interpro_box[inter]:
				interpro_box[inter].append(g)
				
	return interpro_box




def __loadAnnotationFile(annotation_file):

	r = {}

	id_box = _0_Preprocess.ID_BOX

	f = open(annotation_file, 'r')


	while(True):
		lines = f.readlines(10000)
	
		if not lines:
			break
		
		for s in lines:
			s = s.strip()
			if len(s) == 0: continue
			
	
			t = s.split('\t')
			
			gene_id = t[0].strip()
			pfam_id = t[1].strip()
	
			
			string_id = id_box.getSTRING_ID_of(gene_id)
			
			
			
			if string_id is not None:
				# add an TF
				if not pfam_id in r:
					r[pfam_id] = []
				if not pfam_id in r[pfam_id]:
					r[pfam_id].append(string_id)
	
	
	f.close()


	return r




def __initializeGO(go_class):
	r = {} # go_id as a key, values are gene id
	go_ids = go_class.go_terms.keys()
	for g in go_ids:
		r[g] = []
		#print g
	return r









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

	min_p_value = 1
	for goid in pvalues.keys():
		p = pvalues[goid]
		if p>0.0 and p<min_p_value:
			min_p_value = p

	method = 1
	for g in total_genes: 

		score = 0.0
		desc = []
		for goid in go_box_rev[g]:
			
			p_pos = pvalues[goid]
			
			# p-value=0.0 --> causes an error
			if p_pos == 0.0:
				p_pos = min_p_value

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








# =============================================


def __train(interpro, tf_box, tf_box_rev, total_genes, mod_pos, output, test_set_all):


	global P_VALUE

	# p value
	pos_pvalues = __calculatePvalue(tf_box, total_genes, mod_pos, output)
	scores = __calculateScores(tf_box, tf_box_rev, total_genes, pos_pvalues)	
	
	#__save(go, go_box, mod_pos, mod_neg, pvalues, output)


	if P_VALUE == True:
		xbox = {}
		for i in range(ITERATION):
			rnd_mods = random.sample(total_genes, len(mod_pos))
			ps = __calculatePvalue(tf_box, total_genes, rnd_mods, output)
			sc = __calculateScores(tf_box, tf_box_rev, total_genes, ps)
			for g in sc.keys():
				if not g in xbox:
					xbox[g] = []
				xbox[g].append(sc[g])
			
		new_scores = {}
		for g in total_genes:
			mean = statistics.average(xbox[g])
			std = statistics.stdev(xbox[g])
			if std == 0:
				new_scores[g] = 1
				continue
			
			p = statistics.getNormalDistPvalue(mean, std, scores[g])
			if p == -1:
				p = 1
			new_scores[g] = p
		
		scores = new_scores
		
	__resave(interpro, tf_box, tf_box_rev, total_genes, mod_pos, pos_pvalues, scores, output, test_set_all, '')




def __resave(desc, tf_box, tf_box_rev, total_genes, mod_pos, pvalues, scores, output, test_set_all, md5):


	ofile = output + '_gene.txt'

	# calculate modifier scores
	box = __rebox(desc, tf_box, tf_box_rev, total_genes, mod_pos, pvalues, scores, output )
	__save2(box, mod_pos, ofile, md5)

def __save2(box, mod_pos, output, md5, order = True):

	if P_VALUE:
		order = False
		
	#box3 = sorted( box.items(), key=lambda (k,v): v[0], reverse = order)
	box3 = sorted(box.items(), key=lambda item: item[1][0], reverse=order)

	f=open(output, 'w')
	f.write('!'+md5+'\n')
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

def __rebox(description, tf_box, tf_box_rev, total_genes, mod_pos, pvalues, scores, output):


	# ==================================================================================================
	r = {}

	for g in total_genes: 


		score = scores[g]
		desc = []
		for o_id in tf_box_rev[g]:
			p_pos = pvalues[o_id]
			name = '---'
			desc.append( o_id + '[' + str(p_pos)  +  ' ] (' + name + ')')

		r[g] = [ score, ','.join(desc)] 
		

	return r


def __save(description, box, mod_pos, pvalues, output):

	f=open(output,'w')


	m_v = pvalues


	s = 'TF_ID\tp-value\t#gene\t#modifier\tDesc\tList'
	f.write(s+'\n')

	# order by p value
	#vx = sorted( box.items(), key=lambda (k,v): m_v[k], reverse=False)
	vx = sorted(box.items(), key=lambda item: m_v[item[0]], reverse=False)

	for gid, genes in vx:
		genes = box[gid]

		
		
		name = '----'

		lst = ','.join( __makeGeneList(genes, mod_pos) )

		tv = len(genes)

		mv = __countModifierNumber(mod_pos, genes)

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

def __calculatePvalue(box, total_genes, mod_pos, output):


	
	pr = float(len(mod_pos)) / float(len(total_genes))
	p = {}

	

	for tf in box.keys():
		# count the number of genes annotated with goid
		mod_number = __countModifierNumber(mod_pos, box[tf])
		total_number = len( box[tf] )


		# hypergeometric distribution

		sample_pop = len(mod_pos)
		sample_succ = mod_number

		pop_pop = len(total_genes)
		pop_succ = total_number
		
		

		pv = statistics.getHypergeomPvalue(sample_succ, sample_pop, pop_succ, pop_pop  )
		if pv == -1:
			pv = 1

		#print pv, sample_pop, sample_succ, pop_pop, pop_succ
		p [tf] = pv

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

def __train2(description, tf_box, tf_box_rev, total_genes, mod_pos, output, test_set_all, md5):


	global P_VALUE

	# p value
	pos_pvalues = __calculatePvalue(tf_box, total_genes, mod_pos, output)
	scores = __calculateScores(tf_box, tf_box_rev, total_genes, pos_pvalues)	
	
	#__save(go, go_box, mod_pos, mod_neg, pvalues, output)


	if P_VALUE == True:

		print ('[PFAM ITERATION] ', ITERATION)

		xbox = {}
		for i in range(ITERATION):
			rnd_mods = random.sample(total_genes, len(mod_pos))
			ps = __calculatePvalue(tf_box, total_genes, rnd_mods, output)
			sc = __calculateScores(tf_box, tf_box_rev, total_genes, ps)
			for g in sc.keys():
				if not g in xbox:
					xbox[g] = []
				xbox[g].append(sc[g])
			
		new_scores = {}
		for g in total_genes:
			mean = statistics.average(xbox[g])
			std = statistics.stdev(xbox[g])
			if std == 0:
				new_scores[g] = 1
				continue
			
			p = statistics.getNormalDistPvalue(mean, std, scores[g])
			if p == -1:
				p = 1
			new_scores[g] = p
		
		scores = new_scores
		
	# save p-values of GO-terms (accesory data)
	__save(description, tf_box, mod_pos, pos_pvalues, output)

	# save results
	__resave(description, tf_box, tf_box_rev, total_genes, mod_pos, pos_pvalues, scores, output, test_set_all, md5)    


def runFinal():

	global RND_ID

	tf_file = _0_Preprocess.PFAM_FILE

	#RND_ID = _0_Preprocess.RND_ID

	# 여기에 리스트 형태로 들어가 있는 질병에 대해서만 disease-specific modifiers를 예측한다.
	predict_disease_specific_modifiers = _0_Preprocess.PREDICT_DISEASE_MODIFIERS
	if len(predict_disease_specific_modifiers) == 0:
		# 예측할게 하나도 없으면 그냥 건너뛴다.
		print("[_6_PFam] 예측할 disease-specific modifiers가 없음 -> 건너뜀")
		return




	mod_maps, mods = __getModifiers()
	# mod_maps={}, key=user id, value=string id
	# mods = [], list of converted string IDs

	
	box1, box_rev1 = preprocess(tf_file)
	


	test_set_genes = None

	for disease in mods.keys():
		
		if not disease in predict_disease_specific_modifiers:
			print('[_6_PFam.py] disease-specific modifiers 예측 안하고 건너뜀: ', disease)
			continue
		
		
		if disease.find('RANDOM') > 0:
			mods[disease] = random.sample(total_genes, len(mods[disease]))




		box = copy.deepcopy(box1)
		box_rev = copy.deepcopy(box_rev1)
		total_genes = __getTotalGeneList(box1)


		#output = OUTPUT_PATH + '/Pfam_' + disease+'_' + str(ITERATION) + '_' + RND_ID + '.txt'
		output = TEMP_PATH + '/Pfam_' + disease + '_' + str(ITERATION) + '_' + RND_ID + '.txt'
		
		
		# md5 within the output file
		x = sorted( copy.deepcopy( mods[disease] ) )
		x.append(_0_Preprocess.TF_FILE)
		x.append(str(ITERATION))
		s = '\n'.join( x )
		md5 = _0_Preprocess.generateMD5(s)
		print ('Modifiers MD5=', md5)
		if md5 == _0_Preprocess.getMD5(output+'_gene.txt')  and _0_Preprocess.REFRESH_RESULTS == False:
			print ('\tUse previous results: ', output+'_gene.txt')
			#print 'MD5 = ', md5[disease]
			
			continue



		print (disease, output)

		mod_pos = mods[disease]

		new_mod_pos = __removeUnannotatedModifiers(total_genes, mod_pos)


		__train2(None, box, box_rev, total_genes, new_mod_pos, output, test_set_genes, md5)


	print ('Done')



def init(config_file, iter_for_pvalue):

	_0_Preprocess.init(config_file)


	global TF_FILE, TEMP_PATH, OUTPUT_PATH, P_VALUE, ITERATION

	PFAM_FILE = _0_Preprocess.PFAM_FILE
	TEMP_PATH = _0_Preprocess.TEMP_PATH
	OUTPUT_PATH = _0_Preprocess.OUTPUT_FOLDER
	P_VALUE = _0_Preprocess.P_VALUE
	ITERATION = iter_for_pvalue



if __name__ == '__main__':


	import sys
	
	
	
	
	
	start_time = datetime.now()


	# python _1_GO.py config.txt test
	print ('''
========================================================    
    [6] PFam
========================================================          
''')
	
	
	if len(sys.argv) == 3:
		config_file = sys.argv[1]
		iter_for_pvalue = eval(sys.argv[2])
		init(config_file, iter_for_pvalue)
		runFinal()
	elif len(sys.argv) == 4:
		if sys.argv[3].lower() == 'test':
			config_file = sys.argv[1]
			iter_for_pvalue = eval(sys.argv[2])
			iteration = 100
			init(config_file, iter_for_pvalue)
			run(iteration)
	

	print ('''
========================================================    
    [6] PFam (End)
========================================================          
''')


	end_time = datetime.now()
	print ('Elapsed time: {}'.format(end_time - start_time))
	print (time.ctime())