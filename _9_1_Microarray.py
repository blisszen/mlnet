# -*- coding: utf-8 -*-
import MULTIPROCESS
import os
import _0_Preprocess
import copy
import random
import math
import sys
from datetime import datetime
import time
import subprocess
import MyUtil
import hashlib
# import dill
import pickle
import sqlite3
# import queue



USE_SQLITE = False
RAMDISK_PATH = './'
#RAMDISK_PATH = 'C:/Users/blisszen/Desktop'

# 'C:/Users/blisszen/Desktop'
# 'Z:/'

# import numpy
# import statistics

# from numba import jit


RND_ID = _0_Preprocess.RND_ID
TEMP_PATH = _0_Preprocess.TEMP_PATH
OUTPUT_PATH = _0_Preprocess.OUTPUT_FOLDER
MONGO_DB_NAME = _0_Preprocess.MONGO_DB
TEMP_PATH = _0_Preprocess.TEMP_PATH

CAL_PER_PROCESSOR = 1  # 10 for fly

mongo_db = None
mongo_col = 'cache'

P_VALUE = True
ITERATION = 1000

THREAD = 0  # 여러 질병이 들어오면 한번에 다 실행한다. CPU 성능이 되면 이걸 4로 올린다. 그러면 4가지 질병을 동시에 계산함.
CACHE_BOX = {}

USE_CACHE = False  # this requires a lot of memory !!! be cautious
MP_PROCESSORS = 2

PYPY = False

DISPLAY = _0_Preprocess.DISPLAY

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
'''


def chunks(seq, size):
	# list�� �ɰ���. �̶� n���� ������ �ɰ���. ����Ʈ�� ������ �ƴ϶� item�� ������ �����ش�.
	return [seq[i:i + size] for i in range(0, len(seq), size)]


def remove(fname):
	try:
		os.remove(fname)
	except Exception as ex:
		print('Failed to remove: ', fname, repr(ex))


# @jit(nopython=True)
def __randomSampling(total_genes, pos_mod_num):
	pos_mod = random.sample(total_genes, pos_mod_num)
	return pos_mod


def saveVariableAsPickle(variable, filename=None):
	global TEMP_PATH, RND_ID

	# print 'temp=', TEMP_PATH
	# print 'RND=',  MyUtil.getRandomString(10)
	# fname = None

	while (True):

		fname = TEMP_PATH + '/pickle_microarray_' + MyUtil.getRandomString(30) + '_' + RND_ID + '.obj'
		if filename is not None:
			fname = filename

		if os.path.exists(fname):
			print('[MicroArray.saveVariableAsPickle] duplicate file name: ' + fname)
			time.sleep(1)
		else:
			# print fname
			try:
				# with open(fname, "wb") as of:
				#	pickle.dump(variable, of)
				of = open(fname, 'wb')
				pickle.dump(variable, of)
				of.close()
				break

			except Exception as ex:
				print('[microarray] Error in saving result pickle. Going to try again: ', time.ctime())
				print(print(repr(ex)))
				time.sleep(1)

	return fname


def loadVariableFromPickle(fname, remove_after_load=False):
	e = None

	try:
		# with open(fname, 'rb') as of:
		#	e = pickle.load(of)
		of = open(fname, 'rb')
		e = pickle.load(of)
		of.close()
	except Exception as ex:
		print('Error in LoadVariableFromPickle: ', fname, repr(ex))

	if remove_after_load:
		try:
			os.remove(fname)
		except Exception as ex:
			print('[loadVariableFromPickle] Error to delete ', fname, repr(ex))

	return e


def __train_k(args):
	global TEMP_PATH, CAL_PER_PROCESSOR

	array_box_fname, total_genes_fname, mod_pos_fname, output, temp_path, queue = args

	TEMP_PATH = temp_path

	array_box = loadVariableFromPickle(array_box_fname)
	total_genes = loadVariableFromPickle(total_genes_fname)
	rnd_mods = loadVariableFromPickle(mod_pos_fname, remove_after_load=True)

	# output = loadVariableFromPickle(output_fname)

	ps_set = {}
	for i in range(CAL_PER_PROCESSOR):
		rm = random.sample(total_genes, len(rnd_mods))
		ps = __calculatePvalue2(array_box, total_genes, rm, output)

		for g in ps:
			if not g in ps_set:
				ps_set[g] = []
			ps_set[g].append(ps[g])

	r_file = saveVariableAsPickle(ps_set)

	# print('[microarray] ', time.ctime(), ' result file was saved into ', r_file)

	queue.put(r_file)


def __train_k2(args):
	
	global TEMP_PATH, USE_SQLITE, RND_ID

	cache_scores_fname, cache_db_file, total_genes_fname, rnd_mods, output, TEMP_PATH, RND_ID, per, queue = args

	cache_scores = None

	if USE_SQLITE == False:
		cache_scores = loadVariableFromPickle(cache_scores_fname)
		
	total_genes = loadVariableFromPickle(total_genes_fname)
	#rnd_mods = loadVariableFromPickle(mod_pos_fname)

	# output = loadVariableFromPickle(output_fname)

	ps_set = {}
	for i in range(per):

		rm = random.sample(total_genes, len(rnd_mods))
		ps = __calculatePvalue3(cache_scores, cache_db_file, total_genes, rm, output)

		for g in ps:
			if not g in ps_set:
				ps_set[g] = []
			ps_set[g].append(ps[g])

	r_file = saveVariableAsPickle(ps_set)

	# print('[microarray] ', time.ctime(), ' result file was saved into ', r_file)


	#remove(mod_pos_fname)

	queue.put(r_file)


def __train2_mpi(cache_scores, cache_file, cache_db, cache_total_genes, total_genes, mod_pos, output, test_set_all, md5):

	global P_VALUE, ITERATION, PYPY, MP_PROCESSORS, TEMP_PATH, CAL_PER_PROCESSOR, RND_ID


	time_out = _0_Preprocess.get_default_opt( ) [ _0_Preprocess.OPT_TIMEOUT ]

	# MP_PROCESSORS = 3

	# return as a dict

	# print 'Gene # without data = ', len(genes_without_data)

	# score
	print('Calculating first scores....')
	pvalues = __calculatePvalue3(cache_scores, cache_db,  total_genes, mod_pos, output)
	ofile = output + '_gene.txt'
	
	# 메모리 비우기
	#----------------------------------------
	cache_scores = None
	# ----------------------------------------
	

	if P_VALUE == True:
	
		print('Calculating reference scores...')

		box = {}

		# cache_file = saveVariableAsPickle(cache_scores)
		#total_genes_file = saveVariableAsPickle(total_genes)
		# output_file = saveVariableAsPickle(output)
		if USE_SQLITE == False:
			print('[Array score pickle] = ', cache_file)
			

		th = []

		print('[Microarray] cpu=', MP_PROCESSORS, 'Iteration = ', ITERATION, ' Rnd per cpu=', CAL_PER_PROCESSOR)

		#steps = int(MP_PROCESSORS * 1.1)  # 작업을 나눠서 처리한다.
		steps = int(MP_PROCESSORS)  # 작업을 나눠서 처리한다.
		#if steps == MP_PROCESSORS:
		#	steps = MP_PROCESSORS + 1

		per = int(ITERATION / steps)

		for i in range(steps):
			job = MULTIPROCESS.JOB_THREAD()
			rnd_mods = random.sample(total_genes, len(mod_pos))
			#rnd_mods_file = saveVariableAsPickle(rnd_mods)
			args = [cache_file, cache_db, cache_total_genes, rnd_mods, output, TEMP_PATH, RND_ID, per]
			job.set_args(__train_k2, args)
			job.set_timeout(time_out) # 30분 안에 종료 안되면 재실행.

			th.append(job)

		print('[Microarray] Running multi-processes')

		ret = MULTIPROCESS.runMultiprocesses(th, max_cpu=steps, sleep_interval_seconds=10,
		                                     title='[Microarray Calculation]', cutoff=MP_PROCESSORS)

		print('[Microarray] Merging resulting data...')

		# ===== 작업해야함.

		# pypy �޸� ���� �ذ� �ڵ�

		box_f = {}
		box_t = {}

		for ps_file in ret:

			ps = loadVariableFromPickle(ps_file, remove_after_load=True)

			for g in total_genes:

				try:
					x = box_t[g]
				except:
					box_t[g] = []

				for p in ps[g]:
					box_t[g].append(p)

		# calculate average and stdev
		for g in box_t:
			avg = average(box_t[g])
			std = stdev(box_t[g])

			box_f[g] = [avg, std]

		# print '[Microarray] calculating p-values...'
		# process를 여분으로 더 돌려서 작업하는데, 이게 끝나면 결과 파일이 생성됨. 어쨌든 지워야 하네.
		for j in th:
			r_file = j.getResult()
			if r_file is not None:
				try:
					if os.path.exists(r_file):
						os.remove(r_file)
				except Exception as ex:
					pass

		# try ----------------------------------
		# $remove(cach)
		#remove(total_genes_file)
		# remove(output_file)

		new_scores = {}
		temp_scores = {}

		print('[Microarray] calculating p-values...')

		if PYPY:

			for g in total_genes:

				# avg = MyUtil.mean(box[g])
				# std = MyUtil.pstdev(box[g])

				if g in box_f:
					avg, std = box_f[g]

					temp_scores[g] = [pvalues[g], avg, std]
				else:
					# ���� �����Ͱ� ������ pvalue=1 �������� �����.
					temp_scores[g] = [pvalues[g], pvalues[g], 0]

			pvalues = MyUtil.calculate_pvalues(temp_scores, _0_Preprocess.TEMP_PATH)

		else:

			import statistics

			for g in total_genes:
				if g in box_f:
					mean, std = box_f[g]
					if std == 0:
						new_scores[g] = 1.0
						print('std=0', g)
						continue

					p = statistics.getNormalDistPvalue(mean, std, pvalues[g])
					if p == -1:
						p = 1.0
					new_scores[g] = p
				else:
					print('No gene: ', g)
					new_scores[g] = 1.0

			pvalues = new_scores
	else:
		print('No p-value calculation')
		
	box2 = __rebox(total_genes, mod_pos, pvalues, output)
	__save2(box2, mod_pos, ofile, md5)


def __train2(array_box, total_genes, mod_pos, output, test_set_all, md5):
	global P_VALUE, ITERATION, PYPY

	# return as a dict

	# print 'Gene # without data = ', len(genes_without_data)

	# score
	print('Calculating....')
	pvalues = __calculatePvalue2(array_box, total_genes, mod_pos, output)
	ofile = output + '_gene.txt'

	if P_VALUE == True:
		box = {}

		for i in range(ITERATION):

			print('[Microarray Iteration]', i, '/', ITERATION)

			rnd_mods = random.sample(total_genes, len(mod_pos))
			ps = __calculatePvalue2(array_box, total_genes, rnd_mods, output)
			for g in ps:
				try:
					x = box[g]
				except:
					box[g] = []
				box[g].append(ps[g])

		print('[Microarray] calculating p-values...')

		new_scores = {}
		temp_scores = {}

		if PYPY:

			for g in total_genes:
				avg = average(box[g])
				std = stdev(box[g])

				temp_scores[g] = [pvalues[g], avg, std]

			pvalues = MyUtil.calculate_pvalues(temp_scores, _0_Preprocess.TEMP_PATH)

		else:
			import statistics
			for g in total_genes:
				mean = average(box[g])
				std = stdev(box[g])
				if std == 0:
					new_scores[g] = 1
					continue

				p = statistics.getNormalDistPvalue(mean, std, pvalues[g])
				if p == -1:
					p = 1
				new_scores[g] = p

			pvalues = new_scores

	box2 = __rebox(array_box, total_genes, mod_pos, pvalues, output)
	__save2(box2, mod_pos, ofile, md5)


def runFinal():
	# load cache
	global THREAD, P_VALUE, ITERATION, RND_ID, TEMP_PATH, USE_SQLITE
	global RAMDISK_PATH
	
	# 여기에 리스트 형태로 들어가 있는 질병에 대해서만 disease-specific modifiers를 예측한다.
	predict_disease_specific_modifiers = _0_Preprocess.PREDICT_DISEASE_MODIFIERS
	if len(predict_disease_specific_modifiers) == 0:
		# 예측할게 하나도 없으면 그냥 건너뛴다.
		print("[_9_1_Microarray_pypy.py] 예측할 disease-specific modifiers가 없음 -> 건너뜀")
		return

	maps, mods = _0_Preprocess.getModifiers()
	dd = list(mods)
	dd.sort()

	data_folder = _0_Preprocess.MICROARRAY_FOLDER

	array_files = []
	array_file_string = data_folder

	files, folders = MyUtil.getFileListIn(data_folder)
	for f in files:
		if f.find('norm.txt') < 0:
			#print('Skipped microarray file = ', f)
			continue

		array_files.append(f)
		array_file_string += ',' + f
		
	

	dd2 = []
	for disease in dd:

		if not disease in predict_disease_specific_modifiers:
			print('[_9_1_Microarray_pypy.py] disease-specific modifiers 예측 안하고 건너뜀: ', disease)
			continue

		if disease.find('RANDOM') > 0:
			dd2.append(disease)
			continue

		# train_output = OUTPUT_PATH + '/MicroArray_'+disease+ '_' + str(P_VALUE) + '_' + str(ITERATION) + '_' + RND_ID + '.txt'
		train_output = TEMP_PATH + '/MicroArray_' + disease + '_' + str(P_VALUE) + '_' + str(
			ITERATION) + '_' + RND_ID + '.txt'

		'''
		# md5 within the output file
		x = sorted( copy.deepcopy( mods[disease] ) )
		x.append(_0_Preprocess.MICROARRAY_FOLDER)
		x.append(str(ITERATION))

		for fi in array_files:
			x.append(fi)

		s = '\n'.join( x )
		md5 = _0_Preprocess.generateMD5(s)
		print ('Modifiers MD5=', md5)
		if md5 == _0_Preprocess.getMD5(train_output+'_gene.txt') and _0_Preprocess.REFRESH_RESULTS == False:
			print ('\tUse previous results: ', train_output+'_gene.txt')
			#print 'MD5 = ', md5[disease]

		else:
			dd2.append(disease)
		'''
		dd2.append(disease)

	if len(dd2) == 0:
		print('======================================')
		print(' [Microarray]   All cached before ')
		print('======================================')
		return






	#md5 = _0_Preprocess.generateMD5(array_file_string)
	md5 = os.path.basename(data_folder)
	
	
	cache_db = RAMDISK_PATH + '/_cache_' + md5 + '.microarray.db'
	# REST_APP.py에서 cache_db 파일을 RAMDISK_PATH로 옮겨서 사용함.
	cache_file = './_cache_' + md5 + '.microarray.pickle'

	
	#if not os.path.exists(cache_db):
	#	cache_db = './_cache_' + md5 + '.microarray.pickle'
		
	
	
	cache_total_genes = './_cache_' + md5 + '.microarray_genes.pickle'
	
	
	
	
	#cache_file_md5 = cache_file + '.md5.txt'
	print('cache file = ', cache_file)
	print('cache db = ', cache_db)
	print('cache genes = ', cache_total_genes)
	
	
	total_genes = []
	flag = False

	#if os.path.exists(cache_db) and os.path.exists(cache_file) and os.path.exists(cache_total_genes):  # and _0_Preprocess.REFRESH_RESULTS == False:
	if os.path.exists(cache_file) and os.path.exists(cache_total_genes):  # and _0_Preprocess.REFRESH_RESULTS == False:

		# md5 읽어서 비교한다.
		#f = open(cache_file_md5, 'r')
		#cache_md5 = f.readline().strip()
		#f.close()

		#if md5 == cache_md5:
		print('[Microarray] Loading cache scores = ' + cache_file)
		
		if USE_SQLITE == False:
			cache_scores = loadVariableFromPickle(cache_file)
		total_genes = loadVariableFromPickle(cache_total_genes)
		

	else:
		print('[Microarray] Cache not found: ' + cache_file)


		# multiprocess 없이 cache 만들기.
		cache_scores = preprocess2(data_folder)  # list
		
		# multiprocess 이용하려면 아래 함수로 교체하거나,
		# _9_X_make_microarray_cache.py를 사용하면 됨. <- 이걸 추천함.
		#cache_scores = preprocess3(data_folder)  # list
		
		
		# [  [file, value]   ]    <- valie is dict

		print('[_9 Microarray.py] saving cache...', cache_file)
		# save md5

		# saveVariableAsPickle(cache_scores, filename = cache_file)

		# if os.path.exists(cache_file) == False:
		f = open(cache_file, 'wb')
		pickle.dump(cache_scores, f)
		f.close()
	
		# ------------------------------
		
		box = {}
		for key in cache_scores:
			g1, g2 = key.split(':')
			box[g1] = 0
			box[g2] = 0
		
		total_genes = list(box)  # microarray 안에 포함된 모든 유전자 (STRING ID)
		f=open(cache_total_genes, 'wb')
		pickle.dump(total_genes, f)
		f.close()
		
		
		# -=----------------------------------------------------
		# db에 저장한다.
		if os.path.exists(cache_db):
			remove(cache_db)
		
		conn = sqlite3.connect(cache_db)
		cur = conn.cursor()
		
		table = 'CREATE TABLE expression_scores(id TEXT, score REAL)'
		cur.execute(table)
		
		i = 0
		step = 1000
		cnt = 0
		total = len(cache_scores)
		
		for key in cache_scores:
			
			# display progress
			
			i += 1
			cnt += 1
			if i == step:
				i = 0
				per = round(cnt / total * 100, 2)
				print('Progress : ', per, '%', end='\r')
			
			cur.execute('INSERT INTO expression_scores VALUES(?, ?)',
			            [key, cache_scores[key]])
		
		# indexing
		conn.commit()
		
		conn.execute('CREATE INDEX id_index ON expression_scores (id)')
		conn.commit()
		
		conn.close()
		# md5도 저장한다.
		#f = open(cache_file_md5, 'w')
		#f.write(md5)
		#f.close()

	
	
	if USE_SQLITE:
		cache_scores = None # 메모리 비운다.

	jobs = []

	for disease in dd2:

		if disease.find('RANDOM') > 0:
			mods[disease] = random.sample(total_genes, len(mods[disease]))

		print('-----------------------')
		print(disease)
		# train_output = OUTPUT_PATH + '/MicroArray_'+disease+ '_' + str(P_VALUE) + '_' + str(ITERATION) + '_' + RND_ID + '.txt'
		train_output = TEMP_PATH + '/MicroArray_' + disease + '_' + str(P_VALUE) + '_' + str(
			ITERATION) + '_' + RND_ID + '.txt'

		'''
		# md5 within the output file
		x = sorted( copy.deepcopy( mods[disease] ) )
		x.append(_0_Preprocess.MICROARRAY_FOLDER)
		x.append(str(ITERATION))


		for fi, val in array_box:
			x.append(fi)


		s = '\n'.join( x )
		md5 = _0_Preprocess.generateMD5(s)
		print ('Modifiers MD5=', md5)
		if md5 == _0_Preprocess.getMD5(train_output+'_gene.txt') and _0_Preprocess.REFRESH_RESULTS == False:
			print ('\tUse previous results: ', train_output+'_gene.txt')
			#print 'MD5 = ', md5[disease]
			continue

		'''
		md5 = MyUtil.getRandomString(30)

		mod_pos = mods[disease]

		new_mod_pos = __removeUnannotatedModifiers(total_genes, mod_pos)

		print('=======> [Microarray] Total gene # = ', len(total_genes), ' Modifiers = ', len(new_mod_pos))

		if len(new_mod_pos) == 0:
			print('[ERROR!!!! Microarray] ', disease, ' has no expression data=', len(new_mod_pos))
			continue

		total_genes_without_modifiers = __getTotalGenesWithoutModifiers(total_genes, mod_pos)

		s = '==================' + '\n' + \
		    'Disease = ' + disease + '\n' + \
		    'Total annotated gene # = ' + str(len(total_genes)) + '\n' + \
		    'Positive modifiers # = ' + str(len(new_mod_pos))
		# '# of array data = ' + str( len(cache_scores) )

		print(s)

		prank = []

		# mpi

		# print ("[Microarray MPI !!!!!!!!!!!!!!!!!!!!]")
		__train2_mpi(cache_scores, cache_file, cache_db, cache_total_genes, total_genes, new_mod_pos, train_output, None, md5)

		# non-mpi
		# __train2(array_box, total_genes, new_mod_pos, train_output, None, md5)

		'''
		if THREAD == 0:
			__train2(array_box, total_genes, new_mod_pos, train_output, None, md5)
		else:

			th = MULTIPROCESS.JOB_THREAD_ONLY()
			th.set_args( __train2_mp, [array_box, total_genes, new_mod_pos, train_output, None, md5])
			jobs.append(th)
		'''

		print('[_9_Microarray] OUTPUT =', train_output)

	if THREAD > 0:
		print('!!!!!!!!!! Running multiprocessors: ', THREAD)
		MULTIPROCESS.runMultiprocesses(jobs, max_cpu=THREAD, sleep_interval_seconds=10)


# save cache


def run_mp_func(args):
	config_file, disease, iteration, queue = args
	init(config_file)
	r = run_mp(disease, iteration)
	queue.put(r)
	return r


def run_mp(disease, iteration):
	maps, mods = _0_Preprocess.getModifiers()
	data_folder = _0_Preprocess.MICROARRAY_FOLDER

	array_box = preprocess2(data_folder)

	summary_file = TEMP_PATH + '/Microarray_' + RND_ID + '_summary.txt'
	f = open(summary_file, 'w')

	dd = list(mods)
	dd.sort()

	jobs = []
	texts = []

	auc_output = TEMP_PATH + '/Microarray_' + RND_ID + '_' + disease + '_auc.txt'

	total_genes = __getTotalGeneList(array_box)

	mod_pos = mods[disease]

	new_mod_pos = __removeUnannotatedModifiers(total_genes, mod_pos)

	if len(new_mod_pos) == 0:
		print(disease, ' has no expression data=', len(new_mod_pos))
		return 'No modifiers'

	total_genes_without_modifiers = __getTotalGenesWithoutModifiers(total_genes, mod_pos)

	s = '==================' + '\n' + \
	    'Disease = ' + disease + '\n' + \
	    'Total annotated gene # = ' + str(len(total_genes)) + '\n' + \
	    'Positive modifiers # = ' + str(len(new_mod_pos)) + '\n' + \
	    '# of array data = ' + str(len(array_box))

	print(s)

	f.write(s + '\n')

	prank = []

	for i in range(iteration):
		print('\r\t' + RND_ID + ' Disease = ', disease, ' # = ', i + 1)

		pos_rank = __mainJob(disease, array_box, new_mod_pos, total_genes, total_genes_without_modifiers, i)

		print('\t', pos_rank)

		prank.append(pos_rank)

	print('')

	foname = TEMP_PATH + '/Microarray_' + RND_ID + '_' + disease + '_ranks.txt'
	__saveRank(prank, foname)

	pos_auc = __calculateAUC(prank, auc_output)

	s2 = 'Avg rank = ' + str(average(prank)) + ' +- ' + str(stdev(prank)) + '\n' + \
	     'AUC + = ' + str(pos_auc) + '\n'

	print(s2)
	f.write(s2 + '\n')
	f.close()

	return s + '\n' + s2


def run(iteration):
	maps, mods = _0_Preprocess.getModifiers()
	data_folder = _0_Preprocess.MICROARRAY_FOLDER

	array_box = preprocess2(data_folder)

	summary_file = TEMP_PATH + '/Microarray_' + RND_ID + '_summary.txt'
	f = open(summary_file, 'w')

	dd = list(mods)
	dd.sort()
	for disease in dd:

		auc_output = TEMP_PATH + '/Microarray_' + RND_ID + '_' + disease + '_auc.txt'

		total_genes = __getTotalGeneList(array_box)

		mod_pos = mods[disease]

		new_mod_pos = __removeUnannotatedModifiers(total_genes, mod_pos)

		if len(new_mod_pos) == 0:
			print(disease, ' has no expression data=', len(new_mod_pos))
			continue

		total_genes_without_modifiers = __getTotalGenesWithoutModifiers(total_genes, mod_pos)

		s = '==================' + '\n' + \
		    'Disease = ' + disease + '\n' + \
		    'Total annotated gene # = ' + str(len(total_genes)) + '\n' + \
		    'Positive modifiers # = ' + str(len(new_mod_pos)) + '\n' + \
		    '# of array data = ' + str(len(array_box))

		print(s)

		f.write(s + '\n')

		prank = []

		for i in range(iteration):
			print('\r\t' + RND_ID + ' Disease = ', disease, ' # = ', i + 1)

			pos_rank = __mainJob(disease, array_box, new_mod_pos, total_genes, total_genes_without_modifiers, i)

			print('\t', pos_rank)

			prank.append(pos_rank)

		print('')

		foname = TEMP_PATH + '/Microarray_' + RND_ID + '_' + disease + '_ranks.txt'
		__saveRank(prank, foname)

		pos_auc = __calculateAUC(prank, auc_output)

		s = 'Avg rank = ' + str(average(prank)) + ' +- ' + str(stdev(prank)) + '\n' + \
		    'AUC + = ' + str(pos_auc) + '\n'

		print(s)
		f.write(s + '\n')
	f.close()

	print('Result summary = ', summary_file)


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

	auc = average(pos_true_positive_rate_y)
	f.close()

	return auc


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


def __mainJob(disease, array_box, mod_pos, total_genes, total_genes_without_modifiers, iteration_th):
	new_mod_pos, test_pos, test_others = __sampleModifiers(total_genes_without_modifiers, mod_pos)

	if disease.find('RANDOM') > 0:
		new_mod_pos = __randomSampling(total_genes, len(new_mod_pos) + 1)
		total_genes_without_modifiers = __getTotalGenesWithoutModifiers(total_genes, new_mod_pos)
		while (True):
			x_mod_pos, X_test_pos, x_test_others = __sampleModifiers(total_genes_without_modifiers, new_mod_pos)
			if not (test_pos in x_mod_pos or test_pos in x_test_others):
				new_mod_pos = x_mod_pos
				test_others = x_test_others
				test_pos = X_test_pos

				break

	train_output = TEMP_PATH + '/Microarray_' + RND_ID + '__trained_' + disease + '_' + str(
		iteration_th).strip() + '.txt'
	test_set_all = [test_pos] + test_others
	__train(array_box, total_genes, new_mod_pos, train_output, test_set_all)

	train_output2 = train_output + '_gene.txt'
	r1 = __test(train_output2, test_pos, test_others)

	remove(train_output)
	remove(train_output2)

	return r1


def __removeUnannotatedModifiers(total_genes, mods):
	r = []
	for k in mods:
		if k in total_genes and not k in r:
			r.append(k)
	return r


def __getTotalGeneList(array_box):
	total_genes = []

	for fi, pat in array_box:
		for string_id in pat:
			if not string_id in total_genes:
				total_genes.append(string_id)

	return total_genes


def __sampleModifiers(total_genes_without_modifiers, mod_pos):
	mod_pos_list = copy.deepcopy(mod_pos)
	pos_1 = random.sample(mod_pos_list, 1)[0]

	mod_pos_list.remove(pos_1)

	other_genes = random.sample(total_genes_without_modifiers, 99)

	return mod_pos_list, pos_1, other_genes


def __getTotalGenesWithoutModifiers(total_genes, mod_pos):
	r = []
	for t in total_genes:
		if not t in mod_pos:
			r.append(t)
	return r


def hasCache(path):
	m = hashlib.md5()

	stored_md5 = None

	files, folders = MyUtil.getFileListIn(path)
	files2 = []
	for f in files:
		if f.find('norm.txt') < 0:
			#print('Skipped microarray file = ', f)
			continue
		files2.append(f)

	for f in files2:
		m.update(f)

	hash_md5 = m.hexdigest()

	cache_file = path.replace(':', '_').replace('/', '_').replace(',', '_')
	if os.path.exists(cache_file + '.txt'):
		# get md5
		f = open(cache_file + '.txt', 'r')
		for s in f.readlines():
			s = s.strip()
			if len(s) == 0: continue
			if s[0] == '!':
				stored_md5 = s[1:].strip()
			break
		f.close()

	if hash_md5 == stored_md5:
		return True, cache_file, hash_md5
	else:
		return False, cache_file, hash_md5







def preprocess2(path):
	use_cache_for_norm_data = False

	array_box = []

	files, folders = MyUtil.getFileListIn(path)

	files2 = []
	for f in files:
		if f.find('norm.txt') < 0:
			#print('Skipped microarray file = ', f)
			continue
		files2.append(f)

	cnt = 0
	for fi in files2:
		cnt = cnt + 1

		print('loading ', fi, cnt, '/', len(files2))

		f = open(path + '/' + fi, 'r')

		box = {}

		init = True
		# IGNORE the first line in the file
		for s in f.readlines():

			if init:
				init = False
				continue

			s = s.replace('\n', '')
			if len(s) == 0:
				continue

			s = s.upper()

			if 'NA' in s or 'NULL' in s or 'NONE' in s:
				continue

			x = s.split('\t')

			if len(x) > 3:

				gene_id = x[1].strip()
				if gene_id.find('FBGN') >= 0:
					gene_id = gene_id.replace('FBGN', 'FBgn')

				string_id = _0_Preprocess.ID_BOX.getSTRING_ID_of(gene_id)

				# if string_id is None:
				#	print('String ID not found: ', gene_id)

				if string_id is not None:

					del x[0]  # probe id
					del x[0]  # gene id

					exp_pattern = []
					failed = False
					for v in x:
						try:
							exp_pattern.append(eval(v))
						except:
							failed = True
							break

					if not failed:
						# avg, std
						# avg = MyUtil.mean( exp_pattern )
						# std = MyUtil.pstdev( exp_pattern )

						avg = average(exp_pattern)
						std = stdev(exp_pattern)

						# print(gene_id, avg, std)

						if std != 0.0:
							c_exp_pattern = []
							for ev in exp_pattern:
								vx = (ev - avg) / std
								c_exp_pattern.append(vx)

							# box[string_id]=[ numpy.array(exp_pattern), numpy.array(c_exp_pattern)]
							if len(exp_pattern) > 1:
								box[string_id] = [exp_pattern, c_exp_pattern]
						# box[string_id]=[ numpy.array(exp_pattern), numpy.array(c_exp_pattern)]
						# print('Microarray ', gene_id, string_id, 'was added')

		f.close()

		print('[Microarray.Preprocess2.for] box len = ', len(box))
		array_box.append([fi, box])

	print('')
	print("Preparation of expression data is done")
	print("Preprocessing and caching expression data.")

	# ==========================================
	# 이제 cache score 를 구해서 모두 저장한다.

	cache_scores = {}
	total_genes = __getTotalGeneList(array_box)

	total = float(len(total_genes))
	show = 0
	cnt = 0

	print('Total number of genes = ', len(total_genes))

	for i in range(len(total_genes) - 1):

		cnt += 1
		show += 1

		if show == 100:
			print('Microarray cache calculation = ', round(cnt / total * 100.0, 2), '%    ', end='\r')
			show = 0

		for j in range(i + 1, len(total_genes)):

			g1 = total_genes[i]  # 유전자 1
			g2 = total_genes[j]  # 유전자 2

			# score 계산함.
			key = __sort(g1, g2)

			# 이미 존재하는 것이면 건너 뛴다.
			if key in cache_scores:
				continue

			score = 0

			# array data에 유전자 g1, g2가 있으면 계산한다.
			for fi, exp_set in array_box:

				try:
					x = exp_set[g1]
					y = exp_set[g2]
				except:
					# not available for g and m
					continue

				p = 0.0

				for a, b in zip(exp_set[g1][1], exp_set[g2][1]):
					p += a * b
				p = p / float(len(exp_set[g1][1]) - 1.0)

				score += abs(p)  #

			# 계산이 안되는것은 0으로 처리한다.
			if math.isnan(score) == True:
				score = 0.0

			cache_scores[key] = score

	# return array_box
	return cache_scores




'''
def preprocess2_sqlite(conn, path):
	
	
	import sqlite3
	# table 부터 만든다.
	

	
	use_cache_for_norm_data = False

	array_box = []

	files, folders = MyUtil.getFileListIn(path)

	files2 = []
	for f in files:
		if f.find('norm.txt') < 0:
			#print('Skipped microarray file = ', f)
			continue
		files2.append(f)

	cnt = 0
	for fi in files2:
		cnt = cnt + 1

		print('loading ', fi, cnt, '/', len(files2))

		f = open(path + '/' + fi, 'r')

		box = {}

		init = True
		# IGNORE the first line in the file
		for s in f.readlines():

			if init:
				init = False
				continue

			s = s.replace('\n', '')
			if len(s) == 0:
				continue

			s = s.upper()

			if 'NA' in s or 'NULL' in s or 'NONE' in s:
				continue

			x = s.split('\t')

			if len(x) > 3:

				gene_id = x[1].strip()
				if gene_id.find('FBGN') >= 0:
					gene_id = gene_id.replace('FBGN', 'FBgn')

				string_id = _0_Preprocess.ID_BOX.getSTRING_ID_of(gene_id)

				# if string_id is None:
				#	print('String ID not found: ', gene_id)

				if string_id is not None:

					del x[0]  # probe id
					del x[0]  # gene id

					exp_pattern = []
					failed = False
					for v in x:
						try:
							exp_pattern.append(eval(v))
						except:
							failed = True
							break

					if not failed:
						# avg, std
						# avg = MyUtil.mean( exp_pattern )
						# std = MyUtil.pstdev( exp_pattern )

						avg = average(exp_pattern)
						std = stdev(exp_pattern)

						# print(gene_id, avg, std)

						if std != 0.0:
							c_exp_pattern = []
							for ev in exp_pattern:
								vx = (ev - avg) / std
								c_exp_pattern.append(vx)

							# box[string_id]=[ numpy.array(exp_pattern), numpy.array(c_exp_pattern)]
							if len(exp_pattern) > 1:
								box[string_id] = [exp_pattern, c_exp_pattern]
						# box[string_id]=[ numpy.array(exp_pattern), numpy.array(c_exp_pattern)]
						# print('Microarray ', gene_id, string_id, 'was added')

		f.close()

		print('[Microarray.Preprocess2.for] box len = ', len(box))
		array_box.append([fi, box])

	print('')
	print("Preparation of expression data is done")
	print("Preprocessing and caching expression data.")

	# ==========================================
	# 이제 cache score 를 구해서 모두 저장한다.
	
	cur = conn.cursor()
	
	table = 'CREATE TABLE expression_scores(id TEXT, score REAL)'
	
	conn.execute(table)


	cache_scores = {}
	total_genes = __getTotalGeneList(array_box)

	total = float(len(total_genes))
	show = 0
	cnt = 0

	print('Total number of genes = ', len(total_genes))

	for i in range(len(total_genes) - 1):

		cnt += 1
		show += 1

		if show == 100:
			print('Microarray cache calculation = ', round(cnt / total * 100.0, 2), '%    ', end='\r')
			show = 0

		for j in range(i + 1, len(total_genes)):

			g1 = total_genes[i]  # 유전자 1
			g2 = total_genes[j]  # 유전자 2

			# score 계산함.
			key = __sort(g1, g2)

			# 이미 존재하는 것이면 건너 뛴다.
			if key in cache_scores:
				continue

			score = 0

			# array data에 유전자 g1, g2가 있으면 계산한다.
			for fi, exp_set in array_box:

				try:
					x = exp_set[g1]
					y = exp_set[g2]
				except:
					# not available for g and m
					continue

				p = 0.0

				for a, b in zip(exp_set[g1][1], exp_set[g2][1]):
					p += a * b
				p = p / float(len(exp_set[g1][1]) - 1.0)

				score += abs(p)  #

			# 계산이 안되는것은 0으로 처리한다.
			if math.isnan(score) == True:
				score = 0.0

			#cache_scores[key] = score
			
			cur.execute('INSERT INTO expression_scores VALUES(?, ?)',
			            [key, score])
			
	
	# index 생성
	conn.execte('CREATE INDEX id_index ON expression_scores (id)')
	
	# return nothing
	return
'''

def preprocess3(path):
	
	global TEMP_PATH
	
	
	use_cache_for_norm_data = False

	array_box = []

	files, folders = MyUtil.getFileListIn(path)

	files2 = []
	for f in files:
		if f.find('norm.txt') < 0:
			#print('Skipped microarray file = ', f)
			continue
		files2.append(f)

	cnt = 0
	for fi in files2:
		cnt = cnt + 1

		print('loading ', fi, cnt, '/', len(files2))

		f = open(path + '/' + fi, 'r')

		box = {}

		init = True
		# IGNORE the first line in the file
		for s in f.readlines():

			if init:
				init = False
				continue

			s = s.replace('\n', '')
			if len(s) == 0:
				continue

			s = s.upper()

			if 'NA' in s or 'NULL' in s or 'NONE' in s:
				continue

			x = s.split('\t')

			if len(x) > 3:

				gene_id = x[1].strip()
				if gene_id.find('FBGN') >= 0:
					gene_id = gene_id.replace('FBGN', 'FBgn')

				string_id = _0_Preprocess.ID_BOX.getSTRING_ID_of(gene_id)

				# if string_id is None:
				#	print('String ID not found: ', gene_id)

				if string_id is not None:

					del x[0]  # probe id
					del x[0]  # gene id

					exp_pattern = []
					failed = False
					for v in x:
						try:
							exp_pattern.append(eval(v))
						except:
							failed = True
							break

					if not failed:
						# avg, std
						# avg = MyUtil.mean( exp_pattern )
						# std = MyUtil.pstdev( exp_pattern )

						avg = average(exp_pattern)
						std = stdev(exp_pattern)

						# print(gene_id, avg, std)

						if std != 0.0:
							c_exp_pattern = []
							for ev in exp_pattern:
								vx = (ev - avg) / std
								c_exp_pattern.append(vx)

							# box[string_id]=[ numpy.array(exp_pattern), numpy.array(c_exp_pattern)]
							if len(exp_pattern) > 1:
								box[string_id] = [exp_pattern, c_exp_pattern]
						# box[string_id]=[ numpy.array(exp_pattern), numpy.array(c_exp_pattern)]
						# print('Microarray ', gene_id, string_id, 'was added')

		f.close()

		print('[Microarray.Preprocess3.for] box len = ', len(box))
		array_box.append([fi, box])

	print('')
	print("Preparation of expression data is done")
	print("Preprocessing and caching expression data.")

	# ==========================================
	# 이제 cache score 를 구해서 모두 저장한다.
	global MP_PROCESSORS
	cpu = MP_PROCESSORS

	total_genes = __getTotalGeneList(array_box)

	total_genes_file = TEMP_PATH + '/' + MyUtil.getRandomString(20) + '_total_genes.pickle'
	array_box_file = TEMP_PATH + '/' + MyUtil.getRandomString(20) + '_array_box.pickle'
	
	print('Save data into pickles....')
	MyUtil.saveVariableAsPickle(total_genes, filename = total_genes_file)
	MyUtil.saveVariableAsPickle(array_box, filename = array_box_file)
	
	print('Total gene # = ', len(total_genes))
	
	thread_jobs = []
	pickle_files = []
	
	array_box = None # 메모리를 비워준다.
	
	groups = MyUtil.divideList( list(range(len(total_genes)))  , cpu )
	
	
	for gene_indexes in groups:
		
		#if i==50:  break
		
		#print('Creating multiprocesses ', (i+1), '/', len(total_genes))
		
		pickle_file = TEMP_PATH + '/'+MyUtil.getRandomString(20) + '_ret.pickle'
		pickle_files.append(pickle_file)
		
		#th = MULTIPROCESS.JOB_THREAD_NO_RETURN()
		th = MULTIPROCESS.JOB_THREAD()
		th.set_args(  preprocess3_mpi, [  gene_indexes, total_genes_file, array_box_file, pickle_file   ]     )
		th.set_timeout( 9999 ) # 20분 안에 안끝나면 재시작
		thread_jobs.append(th)
		
	

	#MULTIPROCESS.runMultiprocesses_no_return(thread_jobs, max_cpu = 3, sleep_interval_seconds=10)
	r = MULTIPROCESS.runMultiprocesses(thread_jobs, max_cpu=cpu, sleep_interval_seconds=10)
	


	# 결과를 받아와서 하나로 합친다.
	cache_scores = {}

	for pick in pickle_files:
		ret = MyUtil.loadVariableFromPickle(pick, remove_after_load=True)
		cache_scores.update(ret)
		
	
	remove(total_genes_file)
	remove(array_box_file)

	# return array_box
	return cache_scores

def preprocess3_mpi(args):
	
	time.sleep( random.random() * 10 )
	
	gene_indexes, total_genes_file, array_box_file, pickle_file, fake_queue = args
	
	total_genes = MyUtil.loadVariableFromPickle(total_genes_file)
	array_box = MyUtil.loadVariableFromPickle(array_box_file)
	
	cache_scores = {}
	
	
	for i in gene_indexes:
	 
		for j in range(i + 1, len(total_genes)):
			
			g1 = total_genes[i]  # 유전자 1
			g2 = total_genes[j]  # 유전자 2
			
			# score 계산함.
			key = __sort(g1, g2)
			
			# 이미 존재하는 것이면 건너 뛴다.
			if key in cache_scores:
				continue
			
			score = 0
			
			# array data에 유전자 g1, g2가 있으면 계산한다.
			for fi, exp_set in array_box:
				
				try:
					x = exp_set[g1]
					y = exp_set[g2]
				except:
					# not available for g and m
					continue
				
				p = 0.0
				
				for a, b in zip(exp_set[g1][1], exp_set[g2][1]):
					p += a * b
				p = p / float(len(exp_set[g1][1]) - 1.0)
				
				score += abs(p)  #
			
			# 계산이 안되는것은 0으로 처리한다.
			if math.isnan(score) == True:
				score = 0.0
			
			cache_scores[key] = score
		
	array_box=None
		
	MyUtil.saveVariableAsPickle(cache_scores, filename = pickle_file)


# =============================================

# @jit(nopython=True)
def shuffle_array_box(array_box):
	print('------ Shuffling expression data ----------')
	# genes = array_box.keys()
	# random.shuffle(genes)
	ret = []

	for fi, expr in array_box:

		genes = list(expr)
		random.shuffle(genes)

		new_expr = {}
		i = 0
		for g in expr:
			new_expr[genes[i]] = expr[g]
			i += 1

		ret.append([fi, new_expr])

	print('------ Shuffling expression data DONE ----------')
	return ret



def __save2(box, mod_pos, output, md5):
	global P_VALUE
	
	order = True
	if P_VALUE:
		order = False

	box3 = sorted(box.items(), key=(lambda x: x[1]), reverse=order)

	f = open(output, 'w')
	f.write('!' + md5 + '\n')
	f.write('Gene\tscore\tDesc\n')

	for k, v in box3:
		score = v[0]
		desc = v[1]

		s = __addTag(k, mod_pos) + k + '\t' + str(score) + '\t' + desc
		f.write(s + '\n')

	f.close()

	print('[_9_Microarray] Result saved = ' + output)


def __addTag(gene, mod_pos):
	tag = ''
	if gene in mod_pos:
		tag = 'm_'

	return tag


def __rebox(total_genes, mod_pos, pvalues, output):
	# calculate score !
	box = {}
	for g in total_genes:

		if not g in pvalues:
			continue

		desc = _0_Preprocess.ID_BOX.getCommonNameOf(g)

		p_pos = pvalues[g]

		score = p_pos
		box[g] = [score, desc]

	return box


def __save(go_box, mod_pos, mod_neg, pvalues, output):
	f = open(output, 'w')

	for index in range(2):

		if pvalues[index] == None:
			continue

		m_v = pvalues[index]

		s = 'TF_ID\tp-value\t#target genes\t#modifier\tList'
		f.write(s + '\n')

		# p value
		# vx = sorted( go_box.items(), key=lambda (k,v): pvalues[index][k]), reverse=False)
		vx = sorted(go_box.items(), key=(lambda k: pvalues[index][k[0]]), reverse=False)

		for gid, genes in vx:

			lst = ','.join(__makeGeneList(genes, mod_pos, mod_neg))

			tv = len(genes)

			mv = None

			if index == 0:
				mv = __countModifierNumber(mod_pos, genes)
			else:
				mv = __countModifierNumber(mod_neg, genes)

			pv = m_v[gid]

			s = gid + '\t' + str(pv) + '\t' + str(tv) + '\t' + str(mv) + '\t' + lst
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


# @jit(nopython=True)
def __calculatePvalue(go_box, total_genes, mod_pos, output):
	r = {}

	for g in total_genes:

		p_total = -100000.0

		for m in mod_pos:

			pm = 1.0
			cnt = 0

			for exp_set in go_box:

				if g in exp_set and m in exp_set:
					g_pattern = exp_set[g]
					m_pattern = exp_set[m]

					# p = statistics.getPearsonCorrelation(g_pattern, m_pattern)
					p = getPearsonCorrelation(g_pattern, m_pattern)
					pm = p * pm

					cnt = cnt + 1

			if abs(pm) > abs(p_total) and cnt:
				p_total = pm

		r[g] = p_total

	return r


# @jit(nopython=True)
def _gene_with_no_data(array_box, total_genes):
	r = {}
	for g in total_genes:
		r[g] = None

	for fi, exp in array_box:
		for g in r:
			try:
				x = exp[g]
				del r[g]
			except:
				pass
	return r  # return as a dict


def __calculatePvalue3(cache_scores, cache_db_file, total_genes, mod_pos, output):
	
	global CACHE_BOX, USE_CACE, mongo_db, USE_SQLITE
	# USE_CACHE = False


	conn = None
	cur = None
	
	if USE_SQLITE:
		print('Opening SQLITE Database: ', cache_db_file)
		conn = sqlite3.connect(cache_db_file, isolation_level=None)
		#conn.execute('pragma journal_mode=MEMORY')
		conn.execute('pragma synchronous=OFF')
		conn.execute('pragma cache_size=10000')
		#cur.execute('pragma mmap_size=10000')
		
		cur = conn.cursor()




	r = {}
	i = 0.0
	to = float(len(total_genes))
	show = 0

	for g in total_genes:

		i += 1.0
		show += 1

		'''
		if show == 1000:
			print('p-value calculation = ', round(i / to * 100.0, 2), '%    ', end='\r')
			show = 0
		'''

		b = 0.0
		score = 0.0

		for m in mod_pos:

			key = __sort(g, m)

			try:
				
				if USE_SQLITE:
					# db에서 읽어온다.
				
					try:
						cur.execute("SELECT score FROM expression_scores WHERE expression_scores.id = '" + key + "' LIMIT 1")
						
						ret = cur.fetchone()
						
						if ret is not None:
							score = ret[0]
							#print('score = ', score)
						
					except Exception as ex:
						print ('[_9_1_Microarray] Error in SQLITE. ', repr(ex))
				
				else:
					score = cache_scores[key]
				
				
				
				b += abs(score)
			except:
				# data가 없다.
				continue

		if math.isnan(b) == False:
			r[g] = b
		else:
			r[g] = 0.0

	return r


def __calculatePvalue2(array_box, total_genes, mod_pos, output):
	global CACHE_BOX, USE_CACHE, mongo_db
	# USE_CACHE = False

	r = {}
	i = 0.0
	to = float(len(total_genes))
	show = 0

	for g in total_genes:

		i += 1.0
		show += 1

		# if P_VALUE == False:
		#	if show == 100:
		#		print '\r', (i/to*100.0), '%          ',
		#		show = 0

		if show == 1000:
			print('p-value calculation = ', round(i / to * 100.0, 2), '%    ', end='\r')
			show = 0

		b = 0.0

		p_total = -100000.0
		max_score = -1000000.0

		for m in mod_pos:

			pm = 1.0

			cnt = 0
			p = 0
			tmp = []
			# ind = 0

			cache_score = 0.0
			avg_correlation = 0.0

			'''
			key = __sort(g,m)
			if USE_CACHE and mongo_db is not None:

				ret = mongo_db.find_one( { 'key': key})
				if ret is not None:
					b += ret['score']
					continue
			'''

			key = __sort(g, m)

			if USE_CACHE:

				try:
					# has a cached score
					v = CACHE_BOX[key]
					if v is not None:
						b += abs(v)
					continue
				except:
					pass

			for fi, exp_set in array_box:

				# ind += 1

				# if exp_set.has_key(g) and exp_set.has_key(m):

				try:
					x = exp_set[g]
					y = exp_set[m]
				except:
					# not available for g and m
					continue

				# modified pearson
				# p = numpy.sum (  numpy.multiply(exp_set[g][1], exp_set[m][1] )   ) / float( len(exp_set[g][1]) - 1.0 )

				p = 0.0

				'''
				g_pattern = exp_set[g][1]
				m_pattern = exp_set[m][1]

				for ind in range(len(g_pattern)):
					p = p + g_pattern[ind]*m_pattern[ind]

				p =  p / float( len(g_pattern) - 1.0 )
				'''

				for a, b in zip(exp_set[g][1], exp_set[m][1]):
					p += a * b
				p = p / float(len(exp_set[g][1]) - 1.0)

				# @ numpy 사용 가능하면 아래 코드 쓰면 됨.
				# p = numpy.dot(exp_set[g][1], exp_set[m][1]) / float( len(exp_set[g][1]) - 1.0 )

				'''
				method=2  # 1: perason correlation coefficient, 2: modified pearson
				p = 0.0





				if method == 1:
					g_pattern = exp_set[g][1] # numpy array
					m_pattern = exp_set[m][1] # numpy array

					p = 0.0
					cnt = cnt + 1

					p = statistics.getPearsonCorrelation(g_pattern, m_pattern)

				elif method == 2:

					g_pattern = exp_set[g][1]
					m_pattern = exp_set[m][1]


					p = numpy.sum (  numpy.multiply(g_pattern, m_pattern)   ) / float( len(g_pattern) - 1.0 )


					#p = 0.0
					#cnt = cnt + 1

					#for ind in range(len(g_pattern)):
					#    p = p + g_pattern[ind]*m_pattern[ind]

					#p =  p / float( len(g_pattern) - 1.0 )


					#p = __cal(g_pattern, m_pattern)
				'''

				# passed = False

				b = b + abs(p)
				cache_score += abs(p)

				'''
				tt = type(p)
				if tt is int or tt is float:
					if math.isnan(p) == False:
						#passed = True
						# ===========================
						b = b + abs(p)
						cache_score += abs(p)
						# ===========================
				'''

			# else:
			#    passed = False
			#    #print 'nan -> ', p
			# else:
			#
			#    #print "Error: ", p, type(p)
			#    passed = False
			#
			# if passed:
			#    b = b + abs(p)
			#    #b += p

			#    #cache_score = cache_score + abs(p)
			#    #tmp.append(abs(p))
			#    #tmp.append(p)

			# average correlations

			'''
			if len(tmp)>0:
				avg_correlation = scipy.stats.mstats.gmean(tmp)

				#avg_correlation = numpy.mean(tmp)

				if max_score < avg_correlation:
					max_score = avg_correlation
			'''

			if USE_CACHE:
				try:
					v = CACHE_BOX[key]
				# print '@##@#$@#$@#$@#$@#$ ??????'
				except:
					# if math.isnan(cache_score) == False:
					CACHE_BOX[key] = cache_score
			# else:
			#    CACHE_BOX[key] = None

			'''
			if USE_CACHE and mongo_db is not None:
				doc = { 'key': key,
						'score': cache_score}
				mongo_db.insert(doc)
			'''

		# r[g] = p_total
		# r[g] = statistics.average(b)
		# r[g] = statistics.sum(b)

		# r[g] = max_score

		if math.isnan(b) == False:
			r[g] = b
		else:
			r[g] = 0.0

		'''
		m = 1

		if m == 1:

			if numpy.isnan(b) == False:
				r[g] = b
		'''

		'''
		elif m ==2:
			if len(tmp)>0:
				#vv = statistics.average(tmp)
				vv = statistics.geometric_mean(tmp)
				#if math.isnan(vv) == False:
				if numpy.isnan(vv) == False:
					r[g] = vv
				else:
					print 'nan error: ', vv, repr(tmp)
		'''

	# print ''
	return r


# @jit(nopython=True)
def __cal(g_pattern, m_pattern):
	p = 0.0

	for ind in range(len(g_pattern)):
		p = p + g_pattern[ind] * m_pattern[ind]

	# p = numpy.sum (  numpy.multiply(g_pattern, m_pattern)   )

	p = p / float(len(g_pattern) - 1.0)

	return p


def __sort(a, b):
	if a > b:
		return a + ':' + b
	else:
		return b + ':' + a


# @jit(nopython=True)
def __countModifierNumber(mod_list, gene_list):
	cnt = 0
	for g in gene_list:
		if g in mod_list:
			cnt = cnt + 1
	return cnt


'''
def saveCache():
	cache_file = './9. Microarray/' +

	global CACHE_BOX
	f=open('')
'''


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

				gid = x[0].replace('m_', '').replace('x_', '')

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


def __calculateAUC(p_rank, output):
	f = open(output + '_pos.txt', 'w')
	f.write('False-positive rate(1-specificity)\tTrue-positive rate(sensitivity)\n')

	start_from = 0.0
	start_to = 1.0
	step = 0.01

	pos_false_positive_rate_x = []  # 1- specificity
	pos_true_positive_rate_y = []  # sensitivity

	neg_false_positive_rate_x = []
	neg_true_positive_rate_y = []

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

	auc_pos = average(pos_true_positive_rate_y)
	f.close()

	return auc_pos


def init(config_file):
	global RND_ID, TEMP_PATH, OUTPUT_PATH, MONGO_DB_NAME, mongo_db, mongo_col, USE_CACHE, P_VALUE

	_0_Preprocess.init(config_file)

	RND_ID = _0_Preprocess.RND_ID
	TEMP_PATH = _0_Preprocess.TEMP_PATH
	OUTPUT_PATH = _0_Preprocess.OUTPUT_FOLDER
	MONGO_DB_NAME = _0_Preprocess.MONGO_DB

	P_VALUE = _0_Preprocess.P_VALUE

	'''
	if MONGO_DB_NAME is not None and USE_CACHE:
		mongo_db = MyMongoDB.NoSQL_MongoDB(host = 'localhost')
		mongo_db.loadDB(MONGO_DB_NAME)
		mongo_db.loadCollection(mongo_col)
	'''


def run_mp_main(config_file, iteration):
	maps, mods = _0_Preprocess.getModifiers()

	jobs = []

	for disease in mods:
		th = MULTIPROCESS.JOB_THREAD()
		th.set_args(run_mp_func, [config_file, disease, iteration])
		jobs.append(th)

	r = MULTIPROCESS.runMultiprocesses(jobs, max_cpu=5,
	                                   sleep_interval_seconds=10)

	s_file = _0_Preprocess.TEMP_PATH + '/' + RND_ID + '_Microarray_summary.txt'
	f = open(s_file, 'w')
	print('------------------------')
	for s in r:
		f.write(s + '\n')
		print(s)
	f.close()

	print('------------------------')
	print(' [SUMMARY_FILE]', s_file)



def runFinal_mp_main(config_file):
	global OUTPUT_PATH, P_VALUE, ITERATION, MP_PROCESSORS

	maps, mods = _0_Preprocess.getModifiers()
	dd = list(mods)
	dd.sort()
	data_folder = _0_Preprocess.MICROARRAY_FOLDER
	array_files = []

	files, folders = MyUtil.getFileListIn(data_folder)
	for f in files:
		if f.find('norm.txt') < 0:
			#print('Skipped microarray file = ', f)
			continue
		#print('Microarray File = ', f)
		array_files.append(f)

	dd2 = []

	jobs = []

	for disease in dd:

		if disease.find('RANDOM') > 0:
			dd2.append(disease)
			continue

		# train_output = OUTPUT_PATH + '/MicroArray_'+disease+ '_' + str(P_VALUE) + '_' + str(ITERATION) + '.txt'
		train_output = TEMP_PATH + '/MicroArray_' + disease + '_' + str(P_VALUE) + '_' + str(ITERATION) + '.txt'

		# md5 within the output file
		x = sorted(copy.deepcopy(mods[disease]))
		x.append(_0_Preprocess.MICROARRAY_FOLDER)
		for fi in array_files:
			x.append(fi)

		s = '\n'.join(x)
		md5 = _0_Preprocess.generateMD5(s)
		print('Modifiers MD5=', md5)
		if md5 == _0_Preprocess.getMD5(train_output + '_gene.txt') and _0_Preprocess.REFRESH_RESULTS == False:
			print('\tUse previous results: ', train_output + '_gene.txt')
		# print 'MD5 = ', md5[disease]

		else:
			dd2.append(disease)
	if len(dd2) == 0:
		print('======================================')
		print(' [Microarray]   All cached before ')
		print('======================================')
		return

	for disease in dd2:
		th = MULTIPROCESS.JOB_THREAD()
		th.set_args(runFinal_mp_func, [config_file, disease])
		jobs.append(th)

	r = MULTIPROCESS.runMultiprocesses(jobs, max_cpu=MP_PROCESSORS,
	                                   sleep_interval_seconds=10)


def getPearsonCorrelation(a_list, b_list):
	'''
	Pearson-correlation 계산한다. List 2개 주면 숫자 아닌 것은 알아서 빼고 계산한다.
	'''

	if len(a_list) != len(b_list):
		raise "Different length of two lists"

	# r, p = scipy.pearsonr(a_list, b_list)
	# return r

	# a2_list, b2_list = __removeAlpha(a_list, b_list) # 안에 숫자 아닌게 있으면 없애버린다.
	a2_list = a_list
	b2_list = b_list

	a_avg = average(a2_list)
	a_std = stdev(a2_list)

	b_avg = average(b2_list)
	b_std = stdev(b2_list)

	r = 0.0
	cnt = 0.0

	v1 = v2 = v3 = v4 = 0.0

	for index in range(len(a_list)):
		a = a2_list[index]
		b = b2_list[index]

		# if ((type(a) is int) or (type(a) is long) or (type(a) is float)) and ((type(b) is int) or (type(b) is long) or (type(b) is float)):
		cnt = cnt + 1.0

		v1 += (a - a_avg) * (b - b_avg)
		v3 += (a - a_avg) ** 2
		v4 += (b - b_avg) ** 2

	# s = ( float(a) - a_avg ) * ( float(b) - b_avg )
	# r = r + s

	# r = r / (cnt-1.0) / (a_std*b_std)
	r = v1 / (v3 * v4) ** 0.5

	return r


def average(a_list):
	# return numpy.mean(a_list)

	return float(sum(a_list)) / float(len(a_list))


def stdev(a_list):
	# a = numpy.array(a_list)
	# return a.std()	
	# return stdev_quick(a_list)
	return stdev_manual(a_list)


def stdev_manual(values):
	option = 1

	if len(values) < 2:
		return 0

	sd = 0.0
	sum = 0.0
	meanValue = average(values)

	for v in values:
		diff = v - meanValue
		sum += diff * diff
	vv = sum / (len(values) - option)
	if vv < 0:
		vv = 0
	sd = math.sqrt(vv)
	return sd


def stdev_quick(data):
	# http://en.wikipedia.org/wiki/Algorithms_for_calculating_variance
	# 위에 소개된 naive algorithm이다.
	# 실제 numpy의 stdev 함수와 비교했을 때 오차가 매우 적다.
	# online-algorithm 이 더 나을 것 같아 해봤더니, 그건 오차가 너무 크다 (아래의 stdev_quick)
	n = 0
	Sum = 0
	Sum_sqr = 0

	for x in data:
		n += 1
		Sum += x
		Sum_sqr += x * x

	# variance = (Sum_sqr - ((Sum*Sum)/n))/(n - 1)
	variance = (Sum_sqr - ((Sum * Sum) / n)) / (n - 1)
	if variance < 1:
		variance = 0
	return math.sqrt(variance)


def main(opt):
	# global mongo_db
	global MP_PROCESSORS, ITERATION, PYPY, P_VALUE, TEMP_PATH, RND_ID

	print('''
========================================================    
	[9] Microarray
========================================================          
''')

	start_time = datetime.now()

	config_file = opt[_0_Preprocess.OPT_CONFIG_FILE]
	MP_PROCESSORS = opt[_0_Preprocess.OPT_MULTIMP]
	rnd_or_not = opt[_0_Preprocess.OPT_RANDOM_OR_NOT]
	test_or_not = opt[_0_Preprocess.OPT_TEST]
	ITERATION = opt[_0_Preprocess.OPT_ITERATION]
	PYPY = opt[_0_Preprocess.OPT_PYPY]

	TEMP_PATH = _0_Preprocess.TEMP_PATH

	RND_ID = opt[_0_Preprocess.OPT_RND_ID]
	_0_Preprocess.RND_ID = RND_ID

	print('[CPU=', MP_PROCESSORS, ']')

	if config_file is None:
		print('!!!!!!!!! No config file')
		return

	init(config_file)

	if test_or_not == False:

		runFinal()

		'''
		# final
		if MP_PROCESSORS == 1:
			runFinal()
		else:
			runFinal_mp_main(config_file)
		'''

	else:
		# calculate auc
		# iteration = 100
		# run(iteration)
		run_mp_main(config_file, ITERATION)

	'''
	if mongo_db is not None:
		mongo_db.close()
	'''

	print('''
========================================================    
	[9] Micriarray (End)
========================================================          
''')

	end_time = datetime.now()
	print('Elapsed time: {}'.format(end_time - start_time))
	print(time.ctime())


if __name__ == '__main__':
	opt = _0_Preprocess.process_args(sys.argv[1:])
	main(opt)

# global mongo_db
# global MP_PROCESSORS

# print '''
# ========================================================    
# [9] Microarray
# ========================================================          
# '''

# start_time = datetime.now()


# if len(sys.argv) == 2:
# config_file = sys.argv[1]
# init(config_file)

# if MP_PROCESSORS == 1:
# runFinal()
# else:
# runFinal_mp_main(config_file)


# elif len(sys.argv) == 3:
# if sys.argv[2].lower() == 'test':

# config_file = sys.argv[1]
# init(config_file)

# iteration = 100
##run(iteration)
# run_mp_main(config_file, iteration)

# '''
# if mongo_db is not None:
# mongo_db.close()
# '''

# print '''
# ========================================================    
# [9] Micriarray (End)
# ========================================================          
# '''


# end_time = datetime.now()
# print ('Elapsed time: {}'.format(end_time - start_time))
# print time.ctime()
