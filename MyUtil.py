# -*- coding: ms949 -*-
import time
import copy
import os
from urllib.request import urlopen
import subprocess
import itertools
import random
import math
import MULTIPROCESS
import multiprocessing
import smtplib
#from sklearn import metrics
import smtplib
from contextlib import closing
from email.mime.text import MIMEText
from email.mime.multipart import MIMEMultipart
from email.header import Header
import sys
import pickle


TEMP_PATH = './TEMP'

def saveVariableAsPickle(variable, filename = None):

	global TEMP_PATH

	# print 'temp=', TEMP_PATH
	# print 'RND=',  MyUtil_pypy.getRandomString(10)



	while (True):
		fname = TEMP_PATH + '/pickle_obj_' + getRandomString(20) + '.obj'
		if filename is not None:
			fname = filename

		if os.path.exists(fname):
			time.sleep(10)
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
				print('[pickle] Error in saving result pickle. Going to try again: ', time.ctime())
				print(print(repr(ex)))
				time.sleep(30)

	return fname


def loadVariableFromPickle(fname, remove_after_load = False):
	e = None

	try:
		# with open(fname, 'rb') as of:
		#	e = pickle.load(of)
		of = open(fname, 'rb')
		e = pickle.load(of)
		of.close()
	except Exception as ex:
		print('Error in LoadVariableFromPickle: ', fname, repr(ex))


	if remove_after_load :
		try:
			os.remove(fname)
		except Exception as ex:
			print('Error in deleting a pickle file: ', fname, repr(ex))

	return e


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
	


def calculateAUC2(class_0_or_1, prediction_0_to_1, pos_label = 1):
	'''
	input: predicted scores [0-1], their classes [0 and 1]
	output: AUC value	
	'''
	
	fpr, tpr, thresholds = metrics.roc_curve(class_0_or_1, prediction_0_to_1, pos_label = 1)
	auc = metrics.auc(fpr, tpr)
	return auc
	
	


def calculateAUC(specificities, sensitivities):
	
	'''
	input: specificity [list], sensitivities [list]
	return: auc
	'''
	
	
	fpr = []
	for s in specificities:
		fpr.append( 1.0 - s )
	tpr = sensitivities
	
	auc = metrics.auc(fpr, tpr)
	return auc

	
	

def getListWithRankAsDict(fname, index, delim, ignore_firstLine = False):
	
	box = []
	box_dict = {}
	
	init = True
	f=open(fname, 'r')
	

	rank = 0
	for s in f.readlines():
		
		
		s = s.strip()
		if len(s) == 0: continue
		
		if init and ignore_firstLine:
			init = False
			continue
		
		value = s.split(delim)[index].strip()
		rank += 1
		
		box.append([value, rank])
		box_dict[value] = rank
		
		#print value, rank
		

	f.close()
	
	
	
	return box_dict


# ------------- LIST 관련 함수 -------------------- #
def getList(fname, index, delim, ignore_firstLine = False):

	# 파일 속 문자열을 읽고
	# delim으로 쪼갠 다음
	# index 번째 내용을 저장한다.

	r = []
	

	init = True

	f = open(fname,'r')

	for s in f.readlines():

		# 주석처리된 줄은 무시한다.
		if len(s)>0:
			if s[0] == '#':
				continue


		if ignore_firstLine and init:
			init = False
			continue


		import os
		s = s.replace('\n','').rstrip( os.linesep )



		if len(s)>0:
			x = s.split(delim)
			r.append( x [ index] )




	f.close()

	return r



def getKeyValueAsDict(fname, key_index, value_index, delim, ignore_firstLine = False):

	# 파일 속 문자열을 읽고
	# delim으로 쪼갠 다음
	# index 번째 내용을 저장한다.

	r = {}


	init = True

	f = open(fname,'r')

	for s in f.readlines():

		# 주석처리된 줄은 무시한다.
		if len(s)>0:
			if s[0] == '#':
				continue

		if ignore_firstLine and init:
			init = False
			continue

		s = s.replace('\n','')



		if len(s)>0:
			x = s.split(delim)

			#print repr(x)

			key = x[key_index]
			value = x[value_index]

			r[key] = value

	f.close()

	return r


def removeRedundances(your_list):
	# 리스트 내부에서 중복된 내용은 제거한다.
	r = []
	for v in your_list:
		if not v in r:
			r.append(v)

	return r

def mergeList(list1, list2):
	# 2개의 리스트 내용을 하나로 합친다.
	# 이때 중복되는 내용은 제거한다.
	y_list = list1 + list2
	return removeRedundances(y_list)

def subtractList(list1, list2):
	# list1에서 list2와 일치하는 내용이 있으면 제거한다.
	r = [        ]
	for v in list1:
		if not v in list2:
			r.append(v)

	return r

#---------------파일 관련 함수
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
			print ("Can't read ", a)

	return files, folders




def isNumeric(value):

	try:
		x=eval(value)
		return True
	except:
		return False



def system(cmd, display=False):
	if display == False:
		#cmd += ' > /dev/null 2>&1'
		cmd += ' 1> /dev/null'
	
	os.system(cmd)
	

def executeDosCommand3(cmd, display = False):
	
	if display == False:
		#cmd += ' > /dev/null 2>&1'
		cmd += ' 1> /dev/null'
		
	os.system(cmd)
	return ""
	
def executeDosCommand4(cmd):
	
	#p = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True)
	p = subprocess.Popen(cmd, stdout=subprocess.PIPE)
	(output, err) = p.communicate()
	return output	
	
#---------------------------------------------
'''
def executeDosCommand2(cmd, printout=True):
	f = subprocess.Popen(cmd, shell = True, stdout = subprocess.PIPE).stdout
	re = ''
	line = ''
	while True:
		last_line = line
		line = f.readline()
		if printout:
			print ("[EXEC] " + line.replace('\n',''))
		re += line
		if not line: break
	return re
'''

def executeDosCommand2(cmd, printout=True):
	return executeDosCommand3(cmd)

def executeDosCommand(cmd):
	'''도스명령을 실행하고, 화면에 뿌려지는 내용을 넘겨준다'''
	dummy, stdout = os.popen4(cmd)
	return stdout.read()




# ===============================================
def fetchURL(url, max_trial=10):
	'''사용자가 입력한 URL을 날리고, 리턴되는 값을 돌려준다.'''


	count = 0
	
	while(count < max_trial):
		try:

			#req = urllib2.Request(url)
			#f = urllib2.urlopen(req)
			f = urlopen(url)
			result = f.read()

			return result
		except:
			print ('Fetching Error')
			time.sleep(5)
			count += 1
			
	if count == max_trial:
		print ('Failed to fetch. Discard this fetch = ' + url)



def combinationN(list_set, n):
	'''list_set에 있는 것 중에서 n개를 골라서 조합만든 후 넘긴다'''
	ret = []
	for s in itertools.combinations(list_set, n):
		ret.append(s)
	return ret

def combinationAll(ss):
	ret = []
	for s in itertools.chain(*map(lambda x: itertools.combinations(ss, x), range(1, len(ss)+1))):
		ret.append(s)
	return ret



#url='http://daum.net'
#print fetchURL(url)
#a = [ '1', '2', '3', '4']
#print combinationN(a, 2)
#print combinationAll(a)














#=============================================================


def divideListByElementNumber(user_list, element_number):
	# element_number 만큼 가지는 리스트로 쪼개기
	return [user_list[i:i+element_number] for i in range(0, len(user_list), element_number  )  ]

def divideList(user_list, group_number, shuffle = False):

	'''
	List를 쪼개준다.
	'''


	if group_number > len(user_list):
		group_number = len(user_list)

	# 사용자가 입력한 리스트를 number 개수로 쪼갠다.
	ret = []
	for i in range(group_number):
		ret.append( [] )

	indexes = list( range(len(user_list)) ) # index만 이용하자.

	if shuffle:
		random.shuffle(indexes) # 섞어서 쓴다.

	g = 0

	for i in range( len(indexes)):
		real_index = indexes[i]
		ret[ g ].append(  user_list [real_index ]  )

		g += 1
		if g == group_number:
			g = 0


	#for i in range(len(ret)):
	#	print i+1, len(ret[i]), '/', len(user_list)

	return ret



# -----------------------------------------------------------------
def sortDict( user_dict, descending_order = False ):
	'''dict 값을 이용해서 정렬한다. 값이 여러개면... 이거 못 쓴다.
	for key, value in .... 이런 식으로 키와 값을 다 받는다.
	'''
	return sorted(user_dict, key=user_dict.get, reverse=descending_order)

def sortDictgetKeys(user_dict,  descending_order = False ):
	'''Dict를 값에 따라 정렬한다.
	return: 키만 순서대로 넘겨준다.'''
	
	
	return sorted(user_dict, key=user_dict.get, reverse=descending_order)
	
	

def sortDictHavingValue(user_dict, descending_order=False):
	'''
	Dictionary안에 list로 값이 들어있는 경우, index에 해당하는 값으로 정렬한 후 키 값만 넘겨준다.
	만약 list가 아니라면 오류가 발생할 것임...
	return: 키 값을 순서대로 넘긴다.
	'''


	r = []
	a = sorted( user_dict.items(), key=lambda x:x[1], reverse=descending_order)
	for key, values in a:
		r.append(key)
	return r


def sortDictHavingList(user_dict, index_in_list, descending_order = False):
	'''
	Dictionary안에 list로 값이 들어있는 경우, index에 해당하는 값으로 정렬한 후 키 값만 넘겨준다.
	만약 list가 아니라면 오류가 발생할 것임...
	키 값만 순서대로 넘긴다.
	'''
	'''
	from operator import itemgetter
	r = []
	for k in sorted( user_dict.iteritems(), key=itemgetter(index_in_values+1), reverse=descending_order):
		r.append(k[0]) # key 값만 저장한다.

	return r
	'''

	r = []
	a = sorted( user_dict.items(), key=lambda x:x[1][index_in_list], reverse=descending_order)
	for key, values in a:
		r.append(key)
	return r



#-----------------------------------------------
def frange(start, end, step):
	'''
	python의 range는 정수만 가능하기에 실수도 가능하게 변경
	'''
	box = []

	h = start

	while(h < end):
		box.append(h)
		h += step

	return box



def getMyFileName(user_file):
	return os.path.basename(user_file)

def getMyPath(user_path):
	f = getMyFileName(user_path)
	return os.path.abspath(user_path).replace(f,'')

	#return os.path.dirname(user_path)



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


#print getRandomString(15)


def getDictSum(dct):

	# dict에 든 값의 합을 넘겨준다.
	total = 0.0
	for i in dct:
		total += dct[i]

	return total



def crossTwoFiles(fname_for_x, index_x_for_key, index_x_for_value, fname_for_y, index_y_for_key, index_y_for_value, output_filename, ignore_firstLine_for_x=False, ignore_firstLine_for_y=False):

	'''
	아래의 crossTwoDicts를 이용하는 건 동일한데,
	대신 편하게 쓰려고 아예 파일 이름을 정해주면 거기서부터 시작한다.
	'''

	dict_x = getKeyValueAsDict(fname_for_x, index_x_for_key, index_x_for_value, ignore_firstLine=ignore_firstLine_for_x)
	dict_y = getKeyValueAsDict(fname_for_y, index_y_for_key, index_y_for_value, ignore_firstLine=ignore_firstLine_for_y)
	crossTwoDicts(dict_x, dict_y, output_filename )



def crossTwoDicts(dict_for_x, dict_for_y, output_file):

	'''
	DICT의 KEY 값이 같은 것만 추려서 파일로 저장한다.
	dict_for_x의 값은 X로
	dict_for_y의 값은 Y로 저장
	'''

	f=open(output_file,'w')
	for k in dict_for_x:
		if k in dict_for_y:
			s = k + '\t' + repr(dict_for_x[k]).strip() + '\t' + repr(dict_for_y[k]).strip()
			f.write(s+'\n')

	f.close()

def getFileSize(fname):
	'''
	파일 용량을 넘겨준다. String으로 된 텍스트파일의 경우 char 1개씩 계산하면 동일해진다.
	'''

	return float( os.path.getsize(fname) )




def sendHtmlEmail(target_email, subject, msg,
                  smtp_server='ssbio.cau.ac.kr',
                  smtp_port=465,
                  user_id = 'web@ssbio.cau.ac.kr',
                  passwd = 'myserver@0987',
                  your_email = 'web@ssbio.cau.ac.kr'):







	'''
	Email 보내는 코드이다.
	현재 Gmail에 최적화되어있음
	'''

	box = MIMEMultipart('alternative')
	box['Subject'] = Header(subject, 'utf-8')
	box['From'] = your_email
	box['To'] = target_email

	part2 = MIMEText(msg, 'html', 'utf-8')
	box.attach(part2)


	#server = smtplib.SMTP_SSL(smtp_server, smtp_port)
	server = smtplib.SMTP_SSL(smtp_server + ':' + str(smtp_port).strip())
	#server.ehlo()

	#server.starttls()


	#server.ehlo()
	server.login(user_id, passwd)
	p = server.sendmail(your_email, target_email, box.as_string())

	server.quit()


def sendEmail(target_email, subject, msg,
              smtp_server='ssbio.cau.ac.kr',
              smtp_port=465,
              user_id = 'web@ssbio.cau.ac.kr',
              passwd = 'myserver@0987',
              your_email = 'web@ssbio.cau.ac.kr',
              content_type = 'plain'):


	#print smtp_server

	msgMIME = MIMEText(msg, content_type)
	msgMIME['To'] = target_email
	msgMIME['From'] = your_email
	msgMIME['Subject'] = Header(subject, 'utf-8')

	#server = smtplib.SMTP_SSL(smtp_server, smtp_port)
	server = smtplib.SMTP_SSL(smtp_server + ':' + str(smtp_port).strip())
	#server.ehlo()

	#server.starttls()


	#server.ehlo()
	server.login(user_id, passwd)
	p = server.sendmail(your_email, target_email, msgMIME.as_string())

	server.quit()



def sendEmailViaGmail(target_email, subject, msg,
                      smtp_server='smtp.gmail.com',
                      smtp_port=587,
                      user_id = 'blisszentmp@gmail.com',
                      passwd = 'freetempgmail',
                      your_email = 'blisszentmp@gmail.com'):



	msgMIME = MIMEText(msg, 'plain')
	msgMIME['To'] = target_email
	msgMIME['From'] = your_email
	msgMIME['Subject'] = Header(subject, 'utf-8')

	server = smtplib.SMTP(smtp_server, smtp_port)
	server.ehlo()

	server.starttls()

	server.ehlo()
	server.login(user_id, passwd)
	p = server.sendmail(your_email, target_email, msgMIME.as_string())

	server.quit()


def sendEmail2(target_email, subject, msg,
               smtp_server='smtp.gmail.com:587',
               user_id = 'blisszentmp@gmail.com',
               passwd = 'freetempgmail',
               your_email = 'blisszentmp@gmail.com'):

	'''
	Email 보내는 코드이다.
	현재 Gmail에 최적화되어있음
	'''


	header = 'From: %s\n' % your_email
	header += 'To: %s\n' % target_email
	header += 'Subject: %s\n' % subject
	content = header + msg

	server = smtplib.SMTP(smtp_server)
	server.starttls()
	server.login(user_id, passwd)
	p = server.sendmail(your_email, target_email, content)
	server.quit()
