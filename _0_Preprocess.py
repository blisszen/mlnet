# -*- coding: utf-8 -*-
import os
import FASTA
import ID
import sys
import random
import hashlib
import MyUtil_pypy



TIMEOUT_FOR_EACH_OMICS = 180 # 각 omics 계산은 30분 안에 끝나야함.


DEBUG_MODE = False
#----------------------------------

MODULE_GO = True
MODULE_INTERPRO = True
MODULE_TF = True
MODULE_PATHWAY = True

MODULE_MICROARRAY =  True
MODULE_USEARCH = True

MODULE_PPI = False    # network propagation을 할지 안할지 정하는 것임. 너무 오래 걸려서 안되겠음.
MODULE_BLAST = False
MODULE_PFAM = False # InterPro를 사용하고 있으므로 이건 제외함.


















# ALLOW PYPY or NOT
MODULE_GO_PYPY = False
MODULE_INTERPRO_PYPY = False
MODULE_TF_PYPY = False
MODULE_PATHWAY_PYPY = False
MODULE_BLAST_PYPY = False
MODULE_PFAM_PYPY = False
MODULE_PPI_PYPY = True # network propagation
MODULE_MICROARRAY_PYPY =  True
MODULE_USEARCH_PYPY = False
MODULE_NETWORK_EXPANSION_PY = False

#----------------------------------


if DEBUG_MODE:
	print('_0_Preprocess.py: Debug mode.......')
	MODULE_GO = True
	MODULE_INTERPRO = False
	MODULE_TF = False
	MODULE_PATHWAY = True
	MODULE_BLAST = False
	MODULE_PFAM = False
	MODULE_PPI = False
	MODULE_MICROARRAY = False
	MODULE_USEARCH = False

TEMP_PATH = None
OUTPUT_FOLDER = None
OUTPUT_FILENAME = None

STRING_FILE = None
STRING_ALIAS_FILE = None #'./ppi/fly.protein.aliases.v10.txt'
STRING_SEQUENCE_FILE = None #'./ppi/fly.protein.sequences.v10.fa'
STRING_CUTOFF = -1

MIN_SEED_NUMBER = 100

PPI_FILES = []  # filename, indexes (protein1, protein2), cutoff


GO_FILE = None
GO_ANNOTATION_FILE = None

INTERPRO_ENTRY_FILE = None
INTERPRO_ANNOTATION_FILE = None
INTERPRO_TO_GO_FILE = None
PFAM_FILE = None

TF_FILE = None

PATHWAY_FILES = []

MODIFIER_FILES = {}
PREDICT_DISEASE_MODIFIERS = []

SEED_NUMBER = None

MICROARRAY_FOLDER = None
RNA_SEQ_FOLDER = None
MONGO_DB = None

RND_ID = str(int(random.random()*1000000)).strip()

DISPLAY = False

ID_BOX = None

P_VALUE = True
P_VALUE_CUTOFF = 1e-10

REFRESH_RESULTS = True

TOTAL_GENE_FILE_TO_CALCULATE = None

# config 파일 이름을 제공해줘야 함.

#--------------------------------------

NETWORK_EXPANSION_PVALUE = False

OPT_MULTIMP = '-cpu'
OPT_MULTIMP_FOR_MICROARRAY = '-cpu_microarray'
OPT_TEST = '-test'
OPT_CONFIG_FILE = '-config'
OPT_RANDOM_OR_NOT = '-random'
OPT_ITERATION = '-iteration'
OPT_PYPY = '-pypy'
OPT_MODIFIER_FILE = '-mod'
OPT_GROUP = '-group'
OPT_NO_EXPANSION = '-no_expansion'
OPT_PVALUE_CUTOFF = '-pvalue_cutoff'
OPT_RND_ID = '-rnd_id'
OPT_DEBUG = '-debug'
OPT_OUTPUT = '-o'
OPT_RELIABLE_STRING_ONLY = '-reliable_only'  # 현재는 모든 STRING ID를 다 보여주는데, 그러지 말자.
OPT_TIMEOUT = '-timeout'
OPT_DISEASE_ITERATION = '-disease_iteration' # single disease일 경우, 몇개의 그룹으로 나눠서 계산할지 옵션
# Conf가 0.8 이상이고, interaction이 orphant가 아닌 것만 결과값으로 보여주자.

#print 'RND_ID = ', RND_ID



def init(config_filename, rnd_id = None):

	global TEMP_PATH, GO_FILE, GO_ANNOTATION_FILE
	global INTERPRO_TO_GO_FILE, INTERPRO_ENTRY_FILE, INTERPRO_ANNOTATION_FILE
	global STRING_ALIAS_FILE, STRING_CUTOFF, STRING_FILE, STRING_SEQUENCE_FILE
	global TF_FILE, PATHWAY_FILE, OUTPUT_FILENAME, OUTPUT_FOLDER
	global SEED_NUMBER, MICROARRAY_FOLDER, MONGO_DB
	global ID_BOX, MODIFIER_FILES, DISPLAY, REFRESH_RESULTS, PPI_FILES
	global PFAM_FILE, RNA_SEQ_FOLDER, PATHWAY_FILES, TOTAL_GENE_FILE_TO_CALCULATE, P_VALUE, NETWORK_EXPANSION_PVALUE
	global P_VALUE_CUTOFF, RND_ID
	global PREDICT_DISEASE_MODIFIERS, OPT_DISEASE_ITERATION


	MODIFIER_FILES = {}


	if rnd_id is not None:
		RND_ID = rnd_id


	PPI_FILES = []
	PATHWAY_FILES = []


	if DISPLAY:
		print('-----------------------------------')
		print('Loading ' + config_filename)

	f=open(config_filename, 'r')
	for s in f.readlines():
		s = s.strip()
		if len(s) == 0: continue
		if s[0]=='#': continue

		t = s.split('=')
		if len(t) == 2:

			#if DISPLAY:
			#	print(s)

			tag = t[0].strip()
			val = t[1].strip()

			if tag == 'TEMP_FOLDER':
				TEMP_PATH = val
			elif tag == 'GO_FILE':
				GO_FILE = val
			elif tag == 'GO_ANNOTATION':
				GO_ANNOTATION_FILE = val
			elif tag == 'INTER_PRO_ENTRY_FILE':
				INTERPRO_ENTRY_FILE = val
			elif tag == 'INTERPRO_ANNOTATION_FILE':
				INTERPRO_ANNOTATION_FILE = val
			elif tag == 'STRING_PPI':
				STRING_FILE = val
			elif tag == 'STRING_ALIAS':
				STRING_ALIAS_FILE = val
			elif tag == 'STRING_SEQUENCE':
				STRING_SEQUENCE_FILE = val
			elif tag == 'INTERPRO_TO_GO_FILE':
				INTERPRO_TO_GO_FILE = val
			elif tag == 'TF_FILE':
				TF_FILE = val
			elif tag == 'RNA_SEQ_FOLDER':
				RNA_SEQ_FOLDER = val
			elif tag == 'P_VALUE':
				# use z-scores for every calculations
				P_VALUE = eval(val)

			#elif tag == 'P_VALUE_CUTOFF':
			#	P_VALUE_CUTOFF = eval(val)

			elif tag == 'PATHWAY_FILE':

				# add ppi files
				# format:  ppi_file: index1(p1): index2(p2): index3(score): cutoff
				t2 = val.split(':')

				PATHWAY_FILES.append(
					[  t2[0].strip(),   # pathway file
					   eval(t2[1]),     # path name index
					   eval(t2[2])     # proteins index (separated by comma)
					   ]

				)

				#print 'pathway=', t2[0]



			elif tag == 'OUTPUT_FOLDER':
				OUTPUT_FOLDER = val
			elif tag == 'STRING_CUTOFF':
				STRING_CUTOFF = eval(val)
			elif tag == 'SEED_NUMBER':
				SEED_NUMBER = eval(val)
			elif tag == 'MICROARRAY_FOLDER':
				MICROARRAY_FOLDER = val
			elif tag == 'DB':
				MONGO_DB = val
			elif tag == 'OUTPUT_FILENAME':
				OUTPUT_FILENAME = val
			elif tag == 'REFRESH_RESULTS':
				REFRESH_RESULTS = eval(val)
			elif tag == 'NETWORK_EXPANSION_PVALUE':
				NETWORK_EXPANSION_PVALUE = eval(val)
				
				'''
				print "-----------------------------------------------------------------------------------------"
				print "Network expansion pvalue = ", eval(val)
				print "-----------------------------------------------------------------------------------------"
				'''
				
			elif tag == 'PFAM_FILE':
				PFAM_FILE = val.strip()

			elif tag == 'PPI_FILE':
				# add ppi files
				# format:  ppi_file: index1(p1): index2(p2): index3(score): cutoff
				t2 = val.split(':')

				separator = '\t'
				sep = t2[5].strip().upper()
				if sep == 'TAB':
					separator = '\t'
				elif sep == 'SPACE':
					separator = ' '
				else:
					print('Error: ' +  sep)


				PPI_FILES.append(
					[  t2[0].strip(),   # ppi file name
					   eval(t2[1]),     # protein1
					   eval(t2[2]),     # protein2
					   eval(t2[3]),     # score
					   eval(t2[4]),      # cutoff
					   separator         # separator
					   ]

				)

			elif tag == 'GENE_LIST':
				TOTAL_GENE_FILE_TO_CALCULATE = val.strip()

			elif tag == 'PREDICT_DISEASE_SPECIFIC_MODIFIERS':
				PREDICT_DISEASE_MODIFIERS = val.strip().replace(' ', '').split(',')   # [AD,HD,SCA1] 이런 형태로 저장됨. 여기에 등록된 질병은 disease-specific modifiers 예측함.

			else:

				if tag[0] == '[':
					disease = tag[1:-1].strip()
					MODIFIER_FILES[disease] = val

				else:
					print(' --> Error: "' + s + '"', tag[0])
					sys.exit(1)

		else:
			print('Error: ' + s)
			sys.exit(1)


			
	# if run on linux, change TEMP_PATH to /media/RamDisk
	# this code should be removed later
	# --------------------------------------------------------------------
	
	if MyUtil_pypy.which_os() == 'linux':
		TEMP_PATH = './TEMP'
	
	# --------------------------------------------------------------------
	



	f.close()
	if DISPLAY:
		print('-----------------------------------')



	if OUTPUT_FILENAME is not None:
		if '[RND]' in OUTPUT_FILENAME:
			OUTPUT_FILENAME = OUTPUT_FILENAME.replace('[RND]', RND_ID)

	global ID_BOX
	print('ALIAS FILE = ', STRING_ALIAS_FILE)
	ID_BOX = ID.ID_CLASS(STRING_ALIAS_FILE, STRING_SEQUENCE_FILE)


def debug_mode_on():
	
	
	global MODULE_USEARCH, MODULE_PPI, MODULE_MICROARRAY, MODULE_BLAST
	
	MODULE_BLAST = False
	MODULE_MICROARRAY = False
	MODULE_PPI = False
	MODULE_USEARCH = False
	
	
	
def getModifiers():

	'''
	Load modifiers 
    return a, b
    # a: dict, key=질병, value가 dict인데, 안에는 key=사용자입력 ID, value는 STRING ID
    # b: dict, key=질병, value가 list인데, 변환된 STRING ID가 들어 있음.
	'''
	
	global MODIFIER_FILES

	m = hashlib.md5()


	r = {}
	l = {}
	md5 = {}

	if DISPLAY:
		print(' Loading modifiers =====================================')
	
	
	
	for d in MODIFIER_FILES:
	
		print(d, '............')
		r[d], l[d] = __loadModifier(MODIFIER_FILES[d])

		if DISPLAY:
			print('User queries = ' + str(  len(r[d]) ) + ' converted STRING proteins = ' + str ( len(l[d])) )

	return r, l



def getSEEDNumber():

	global SEED_NUMBER, MIN_SEED_NUMBER

   
	# SEED_NUMBER = -1이면, modifier 개수 중 max를 취함
	# -1.5 이면 max*2 를 취함.
	# -2 이면, min을 취함
	# -3이면, average
	# -4이면 모든 modifier의 합



	if SEED_NUMBER > 0:
		return SEED_NUMBER
	else:

		r, l = getModifiers()

		# max(no of modifiers)
		if SEED_NUMBER == -1 or SEED_NUMBER == -1.5:
			max_no = -1

			for d in list(l):
				if max_no < len(l[d]):
					max_no = len(l[d])

			if SEED_NUMBER == -1.5:
				max_no = max_no * 2


			max_no = max(max_no, MIN_SEED_NUMBER)



			return max_no




		elif SEED_NUMBER == -2:
			min_no = 100000000000000000

			for d in list(l):
				if min_no > len(l[d]):
					min_no = len(l[d])

			return min_no

		elif SEED_NUMBER == -3:

			# average

			b = []
			for d in list(l):
				b.append(len(l[d]))

			return mean(b)

		elif SEED_NUMBER == -4:
			total = -1

			for d in list(l):
				total +=  len(l[d])

			max_no = max(total, MIN_SEED_NUMBER)

			return max_no

def mean(data):
	"""Return the sample arithmetic mean of data."""
	n = len(data)
	if n < 1:
		raise ValueError('mean requires at least one data point')
	return sum(data)/float(n) # in Python 2 use sum(data)/float(n)

def __loadModifier(fname):

	# return a, b
	# a: dict이며, key는 사용자 입력 ID, value는 STRING ID
	# b: list이며, 변환된 STRING ID



	r = {}
	l = []
	m = hashlib.md5()



	f=open(fname,'r',encoding='utf-8')

	#if DISPLAY:
	print('Loading ' + fname)

	hash_md5 = ''

	for s in f.readlines():
		s = s.replace('\n','')
		if len(s) == 0:
			continue

		if s[0] == '#':
			continue # comment

		x = s.split('\t')
		gene_id = x[0].strip() # gene id
		string_id = ID_BOX.getSTRING_ID_of(gene_id)

		if string_id is None:
			print('Unrecognized ID = '+ gene_id)
			r[gene_id] = None
		else:
			r[gene_id] = string_id
			if not string_id in l:
				l.append(string_id)

			m.update(string_id.encode())


	hash_md5 = m.hexdigest()
	#hash_md5 = 'XXXXXXXXXXXX'

	f.close()


	'''
	for user_id in list(r):

		string_id = r[user_id]
		if string_id is not None:
			l[string_id] = None
	'''

	print('Queried proteins, # = %d ' % ( len(r)) )
	print('\tLoaded proteins, # = %d ' %  (len(l)))
	print('\tMD5 = ' + hash_md5)


	return r, l


def getMD5(fname):
	md5=''

	if not os.path.exists(fname):
		return md5

	f=open(fname,'r')
	for s in f.readlines():
		s = s.strip()
		if len(s) == 0: continue
		if s[0] == '!':
			md5 = s[1:].strip()
			break
	f.close()
	print('Saved MD5 = ' + md5)

	return md5

def generateMD5(string):
	md5 = hashlib.md5()
	md5.update(string.encode())
	return md5.hexdigest()

def regenerate_rnd_id():
	global RND_ID
	RND_ID = str(int(random.random() * 1000000)).strip()
	return RND_ID


def get_default_opt():
	
	global TIMEOUT_FOR_EACH_OMICS
	
	
	opt= { OPT_MULTIMP: 1,
               OPT_TEST: False,
               OPT_CONFIG_FILE: None,
               OPT_RANDOM_OR_NOT: False,
               OPT_ITERATION: 100,
               OPT_PYPY: False,  # 무조건 pypy 사용한다
               OPT_MODIFIER_FILE: None,
               OPT_GROUP: 2,
               OPT_NO_EXPANSION: False,
               OPT_PVALUE_CUTOFF: 1e-100,
               OPT_RND_ID: RND_ID,
               OPT_DEBUG: False,
	       OPT_MODIFIER_FILE: None,
	       OPT_OUTPUT: None,
	       OPT_RELIABLE_STRING_ONLY: False,
	       OPT_MULTIMP_FOR_MICROARRAY: 5,
	       OPT_DISEASE_ITERATION: 4, # single disease일 경우 몇 개의 랜덤 그룹으로 나눌지 정함.
	       OPT_TIMEOUT: TIMEOUT_FOR_EACH_OMICS  # 30분 안에 프로세스 끝나야하는데, 안 끝나면 재실행
               }	
	
	return opt

def process_args(args):

	global OPT_MULTIMP, OPT_TEST, OPT_CONFIG_FILE, OPT_RANDOM_OR_NOT, OPT_ITERATION, OPT_MODIFIER_FILE
	global OPT_GROUP, OPT_NO_EXPANSION, OPT_PVALUE_CUTOFF, RND_ID, OPT_DEBUG, DEBUG_MODE, OPT_OUTPUT
	global OPT_RELIABLE_STRING_ONLY, OPT_MULTIMP_FOR_MICROARRAY, OPT_DISEASE_ITERATION

	opt = get_default_opt()

	for a in args:
		v = a.split('=')

		print(v)

		if v[0] == OPT_MULTIMP:
			opt[OPT_MULTIMP] = eval(v[1])
		elif v[0] == OPT_TEST:
			opt[OPT_TEST] = True
		elif v[0] == OPT_RANDOM_OR_NOT:
			opt[OPT_RANDOM_OR_NOT] = True
		elif v[0] == OPT_CONFIG_FILE:
			opt[OPT_CONFIG_FILE] = v[1].strip()
		elif v[0] == OPT_ITERATION:
			opt[OPT_ITERATION] = eval(v[1])
		elif v[0] == OPT_PYPY:
			opt[OPT_PYPY] = True
		elif v[0] == OPT_MODIFIER_FILE:
			opt[OPT_MODIFIER_FILE] = v[1].strip()
		elif v[0] == OPT_GROUP:
			opt[OPT_GROUP] = eval(v[1])
		elif v[0] == OPT_NO_EXPANSION:
			opt[OPT_NO_EXPANSION] = True
		elif v[0] == OPT_RELIABLE_STRING_ONLY:
			opt[OPT_RELIABLE_STRING_ONLY] = True
		elif v[0] == OPT_PVALUE_CUTOFF:
			opt[OPT_PVALUE_CUTOFF] = eval(v[1])
			print ('#### pvalue cutoff = ', eval(v[1]))
		elif v[0] == OPT_RND_ID:
			opt[OPT_RND_ID] = v[1]
		elif v[0] == OPT_MULTIMP_FOR_MICROARRAY:
			opt[OPT_MULTIMP_FOR_MICROARRAY] = eval(v[1])
		elif v[0] == OPT_DEBUG:
			if eval(v[1]) == True:
				debug_mode_on()
				opt[OPT_DEBUG] = True
				DEBUG_MODE = True
		elif v[0] == OPT_OUTPUT:
			opt[OPT_OUTPUT] = v[1].strip()
		elif v[0] == OPT_DISEASE_ITERATION:
			opt[OPT_DISEASE_ITERATION] = eval(v[1])


	return opt