import random
import os
import subprocess
import math
import time
import sys
#import statistics

def which_os():
	
	'''
	
	
	:return: one of [ 'linux', 'windows', 'OSX ]
	'''
	
	
	if 'linux' in sys.platform.lower():
		return 'linux'
	elif 'win' in sys.platform.lower():
		return 'windows'
	elif 'darwin' in sys.platform.lower():
		return 'osx'
	else:
		return None


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
			print("Can't read ", a)

	return files, folders



#---------------------------------------------
def executeDosCommand2(cmd, printout=True):

	os.system(cmd)
	return ''

	'''
	f = subprocess.Popen(cmd, shell = True, stdout = subprocess.PIPE).stdout
	re = ''
	line = ''
	while True:
		last_line = line
		line = f.readline()
		if printout:
			print("[EXEC] " + line.replace('\n',''))
		re += line
		if not line: break
	return re
	'''


def calculate_pvalues(values, temp_path):

	new_scores = {}

	
	temp_file = temp_path + '/' + getRandomString(10) + '.txt'
	temp_file2 = temp_path + '/' + getRandomString(10) + '.txt'

	f=open(temp_file,'w')
	f.write('value\tmean\tstdev\n')
	for g in list(values):
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
	
	
	
	
	'''
	for g in list(values.keys()):
		val, avg, std = values[g]
		if std != 0.0:
			p = statistics.getNormalDistPvalue(mean, std, value)
			if p == -1:
				p = 1.0
		else:
			p = 1.0
			
		new_scores[g] = p
	'''

	return new_scores

def getMyFileName(user_file):
	return os.path.basename(user_file)

def getMyPath(user_path):
	f = getMyFileName(user_path)
	return os.path.abspath(user_path).replace(f,'')

	#return os.path.dirname(user_path)

def calculate_hypergeom_pvalues(values, temp_path):

	new_scores = {}

	
	temp_file = temp_path + '/' + getRandomString(10) + '.txt'
	temp_file2 = temp_path + '/' + getRandomString(10) + '.txt'

	f=open(temp_file,'w')
	f.write('value\tmean\tstdev\n')
	for g in list(values):
		success_in_sample, sample_number, success_in_pop , pop_number = values[g]


		txt = g + '\t' + str(success_in_sample) + '\t' + str(sample_number) + '\t' + str(success_in_pop) + '\t' + str(pop_number)
		f.write(txt+'\n')

	f.close()

	cmd = 'python _X_hyper_test.py "' + temp_file + '" "' + temp_file2 + '"'
	executeDosCommand2(cmd)

	f=open(temp_file2,'r')
	for s in f.readlines():
		s = s.strip()
		if len(s) == 0: continue

		t = s.split('\t')
		gene = t[0].strip()
		score = eval(t[1])

		new_scores[gene] = score
	
	
	'''
	for g in list(values.keys()):
		val, avg, std = values[g]
		if std != 0.0:
			p = statistics.getHypergeomPvalue( success_in_sample, sample_number, success_in_pop, pop_number)
			if p == -1:
				p = 1.0
		else:
			p = 1.0
			
		new_scores[g] = p
	'''

	return new_scores

def getRandomString(length, string_set =None):
	str1 = 'abcdefghijklmnopqrstuvwxyz'
	str2 = str1.upper()
	numbers = '0123456789'
	
	box = str1 + str2 + numbers
	
	if string_set is not None:
		box = string_set
	

	r = ''

	for i in range(length):
		s = random.sample(box, 1) [0]
		r += s
	return r
