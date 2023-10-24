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


# ------------- LIST ���� �Լ� -------------------- #
def getList(fname, index, delim, ignore_firstLine = False):

	# ���� �� ���ڿ��� �а�
	# delim���� �ɰ� ����
	# index ��° ������ �����Ѵ�.

	r = []
	

	init = True

	f = open(fname,'r')

	for s in f.readlines():

		# �ּ�ó���� ���� �����Ѵ�.
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

	# ���� �� ���ڿ��� �а�
	# delim���� �ɰ� ����
	# index ��° ������ �����Ѵ�.

	r = {}


	init = True

	f = open(fname,'r')

	for s in f.readlines():

		# �ּ�ó���� ���� �����Ѵ�.
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
	# ����Ʈ ���ο��� �ߺ��� ������ �����Ѵ�.
	r = []
	for v in your_list:
		if not v in r:
			r.append(v)

	return r

def mergeList(list1, list2):
	# 2���� ����Ʈ ������ �ϳ��� ��ģ��.
	# �̶� �ߺ��Ǵ� ������ �����Ѵ�.
	y_list = list1 + list2
	return removeRedundances(y_list)

def subtractList(list1, list2):
	# list1���� list2�� ��ġ�ϴ� ������ ������ �����Ѵ�.
	r = [        ]
	for v in list1:
		if not v in list2:
			r.append(v)

	return r

#---------------���� ���� �Լ�
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
	'''��������� �����ϰ�, ȭ�鿡 �ѷ����� ������ �Ѱ��ش�'''
	dummy, stdout = os.popen4(cmd)
	return stdout.read()




# ===============================================
def fetchURL(url, max_trial=10):
	'''����ڰ� �Է��� URL�� ������, ���ϵǴ� ���� �����ش�.'''


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
	'''list_set�� �ִ� �� �߿��� n���� ��� ���ո��� �� �ѱ��'''
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
	# element_number ��ŭ ������ ����Ʈ�� �ɰ���
	return [user_list[i:i+element_number] for i in range(0, len(user_list), element_number  )  ]

def divideList(user_list, group_number, shuffle = False):

	'''
	List�� �ɰ��ش�.
	'''


	if group_number > len(user_list):
		group_number = len(user_list)

	# ����ڰ� �Է��� ����Ʈ�� number ������ �ɰ���.
	ret = []
	for i in range(group_number):
		ret.append( [] )

	indexes = list( range(len(user_list)) ) # index�� �̿�����.

	if shuffle:
		random.shuffle(indexes) # ��� ����.

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
	'''dict ���� �̿��ؼ� �����Ѵ�. ���� ��������... �̰� �� ����.
	for key, value in .... �̷� ������ Ű�� ���� �� �޴´�.
	'''
	return sorted(user_dict, key=user_dict.get, reverse=descending_order)

def sortDictgetKeys(user_dict,  descending_order = False ):
	'''Dict�� ���� ���� �����Ѵ�.
	return: Ű�� ������� �Ѱ��ش�.'''
	
	
	return sorted(user_dict, key=user_dict.get, reverse=descending_order)
	
	

def sortDictHavingValue(user_dict, descending_order=False):
	'''
	Dictionary�ȿ� list�� ���� ����ִ� ���, index�� �ش��ϴ� ������ ������ �� Ű ���� �Ѱ��ش�.
	���� list�� �ƴ϶�� ������ �߻��� ����...
	return: Ű ���� ������� �ѱ��.
	'''


	r = []
	a = sorted( user_dict.items(), key=lambda x:x[1], reverse=descending_order)
	for key, values in a:
		r.append(key)
	return r


def sortDictHavingList(user_dict, index_in_list, descending_order = False):
	'''
	Dictionary�ȿ� list�� ���� ����ִ� ���, index�� �ش��ϴ� ������ ������ �� Ű ���� �Ѱ��ش�.
	���� list�� �ƴ϶�� ������ �߻��� ����...
	Ű ���� ������� �ѱ��.
	'''
	'''
	from operator import itemgetter
	r = []
	for k in sorted( user_dict.iteritems(), key=itemgetter(index_in_values+1), reverse=descending_order):
		r.append(k[0]) # key ���� �����Ѵ�.

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
	python�� range�� ������ �����ϱ⿡ �Ǽ��� �����ϰ� ����
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

	# dict�� �� ���� ���� �Ѱ��ش�.
	total = 0.0
	for i in dct:
		total += dct[i]

	return total



def crossTwoFiles(fname_for_x, index_x_for_key, index_x_for_value, fname_for_y, index_y_for_key, index_y_for_value, output_filename, ignore_firstLine_for_x=False, ignore_firstLine_for_y=False):

	'''
	�Ʒ��� crossTwoDicts�� �̿��ϴ� �� �����ѵ�,
	��� ���ϰ� ������ �ƿ� ���� �̸��� �����ָ� �ű⼭���� �����Ѵ�.
	'''

	dict_x = getKeyValueAsDict(fname_for_x, index_x_for_key, index_x_for_value, ignore_firstLine=ignore_firstLine_for_x)
	dict_y = getKeyValueAsDict(fname_for_y, index_y_for_key, index_y_for_value, ignore_firstLine=ignore_firstLine_for_y)
	crossTwoDicts(dict_x, dict_y, output_filename )



def crossTwoDicts(dict_for_x, dict_for_y, output_file):

	'''
	DICT�� KEY ���� ���� �͸� �߷��� ���Ϸ� �����Ѵ�.
	dict_for_x�� ���� X��
	dict_for_y�� ���� Y�� ����
	'''

	f=open(output_file,'w')
	for k in dict_for_x:
		if k in dict_for_y:
			s = k + '\t' + repr(dict_for_x[k]).strip() + '\t' + repr(dict_for_y[k]).strip()
			f.write(s+'\n')

	f.close()

def getFileSize(fname):
	'''
	���� �뷮�� �Ѱ��ش�. String���� �� �ؽ�Ʈ������ ��� char 1���� ����ϸ� ����������.
	'''

	return float( os.path.getsize(fname) )




def sendHtmlEmail(target_email, subject, msg,
                  smtp_server='ssbio.cau.ac.kr',
                  smtp_port=465,
                  user_id = 'web@ssbio.cau.ac.kr',
                  passwd = 'myserver@0987',
                  your_email = 'web@ssbio.cau.ac.kr'):







	'''
	Email ������ �ڵ��̴�.
	���� Gmail�� ����ȭ�Ǿ�����
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
	Email ������ �ڵ��̴�.
	���� Gmail�� ����ȭ�Ǿ�����
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
