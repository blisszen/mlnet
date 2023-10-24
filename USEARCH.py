# -*- coding: utf-8 -*-'''


# query sequence to blast (NCBI)

import os
import re
import random
import MyUtil



'''
http://www.drive5.com/usearch/manual/cmd_usearch_global.html
'''

D_MELANOGASTER = 'txid7227'
S_SEREVISIAE = 'txid4932'
ALL_ORGANISM = ''
DEBUG = False
DISPLAY = False




class USEARCH():
	'''
	local에 설치된 usearch를 실행한다.
	여기서 지원하는 것은 FASTA로 된 파일로 db를 만드는 것과
	설치된 DB를 대상으로 search를 실행하는 것
	그리고 결과를 추출하는 것이다.

	'''

	result = None
	bin_path = 'D:\\NetworkExpansionModel\\Network_Expansion_Model(human)\\usearch'
	db_path = './TEMP'
	ofile = './TEMP/usearch_output_'+ MyUtil.getRandomString(30) + '.txt'

	if MyUtil.which_os() == 'linux':
		bin_path = './usearch'
		db_path = './TEMP'
		ofile = './TEMP/usearch_output_' + MyUtil.getRandomString(30) + '.txt'

	def __init__(self, bin_path = None, db_path = None):
		self.result = None


		# 사용자가 bin 폴더랑 db 폴더를 정해줘야 한다. 안그러면 기본으로 세팅한다.
		if bin_path != None:
			self.bin_path = bin_path

		if db_path != None:
			self.db_path = db_path



	def parseOutput(self, output_file):
		'''
		# blast output 파일을 받아서 결과를 분석한다.

		#만약 이전 결과가 있으면 그걸 바로 읽어들일까?
		# 현재는 [gene, bits, e-value]를 리스트로 넘겨준다.
		'''


		# 여기 수정해야 함.....

		r = [  ]
		# 들어가는 값은 gene info, bits, e-value

		if not os.path.exists(output_file):
			return r


		f=open(output_file,'r')

		recording = False
		blank = 0

		for s in f.readlines():
			s = s.strip()
			if len(s) == 0: continue

			t = s.split('\t')

			gene_name = t[0]
			gene = t[1]
			bits = eval( t[-1] )
			evalue  = eval( t[-2] )

			r.append( [gene, bits, evalue] )

		f.close()
		
		return r


	def parseOutput2(self, output_file):
		'''
		# blast output 파일을 받아서 결과를 분석한다.

		#만약 이전 결과가 있으면 그걸 바로 읽어들일까?
		# 현재는 [gene, bits, e-value]를 리스트로 넘겨준다.
		'''
		#print ('[USearch.py] Parsing ' + output_file)

		# 여기 수정해야 함.....

		r = { }
		
		# 들어가는 값은 gene info, bits, e-value

		#if not os.path.exists(output_file):
		#	return r


		f=open(output_file,'r')

		recording = False
		blank = 0
		gene_name = None

		for s in f.readlines():

			s = s.strip()
			if len(s) == 0: continue

			t = s.split('\t')

			try:
				gene_name = t[0]
				gene = t[1]
				bits = eval( t[-1] )
				evalue  = eval( t[-2] )
	
				if not gene_name in r:
					r[gene_name] = []
				r[gene_name].append([gene, bits, evalue])
			except:
				if DISPLAY:
					print('Error in USEARCH.parseOutput2: Line=' + s)
				# 이건 USearch 프로그램 자체의 문제임. 가끔 제대로 값이 저장안되는 부분이 있음.
				# 이건 그냥 넘어가야겠다.
				#a=raw_input("Press Enter to quit")
				#sys.exit(1)
			
		f.close()

		return r



	def runUsingSequence(self, sequence, db_file, tool = 'usearch11.0.667_win32.exe', option = '-usearch_local'):
		'''
		sequence 받으면 그걸로 FASTA 만들고 blast돌린다.
		'''

		#
		TEMP_PATH = './TEMP/usearch_input_fasta_'+MyUtil.getRandomString(30)
		if MyUtil.which_os() == 'linux':
			#TEMP_PATH = './TEMP/'+MyUtil.getRandomString(30)
			#tool = 'usearch8.1.1861_i86linux32'
			tool = 'usearch11.0.667_i86linux32'


		ifile = TEMP_PATH + '_test.txt'
		f=open(ifile, 'w')
		f.write('>TEST\n')
		f.write(sequence+'\n')
		f.close()

		ret = self.runUsingFASTA(ifile, db_file, tool=tool, option=option)

		#
		self.remove(ifile)

		return ret



	def runUsingSequences(self, sequences_as_dict, db_file, tool = 'usearch11.0.667_win32.exe', option = '-usearch_local'):
		'''
		sequence 받으면 그걸로 FASTA 만들고 blast돌린다.
		이 경우 sequence는 { name: seq } 형태로 받아야 하며, n 개의 seq를 이용함.
		'''

		#
		TEMP_PATH = './TEMP/'+str( random.randint(0,100000) ).strip()
		if MyUtil.which_os() == 'linux':
			TEMP_PATH = './TEMP/'+str( random.randint(0,100000) ).strip()
			#tool = 'usearch8.1.1861_i86linux32'
			tool = 'usearch11.0.667_i86linux32'


		ifile = TEMP_PATH + '_test.txt'
		f=open(ifile, 'w')
		
		for name in sequences_as_dict:
			f.write('>' + name + '\n')
			f.write(sequences_as_dict[name]+'\n')
		f.close()

		ret = self.runUsingFASTA2(ifile, db_file, tool=tool, option=option)
		# 결과는 ret = { name: score}   score=[gene info, bits, e-value] gene info는 찾아낸 seq gene과 동일.

		self.remove(ifile)



		return ret

	def remove(self, fname):
		try:
			os.remove(fname)
		except:
			pass
	
	#def runUsingFASTA2(self, query_fasta_file, db_file, tool = 'usearch8.1.1861_win32.exe', option = '-usearch_local'):
	def runUsingFASTA2(self, query_fasta_file, db_file, tool='usearch11.0.667_win32.exe', option='-usearch_local'):

		# blast를 돌리고
		# 결과 파일 이름을 넘겨준다.


		ofile = './TEMP/usearch_output_' + MyUtil.getRandomString(30) + '.txt'

		#cmd = self.bin_path + '\\blastp -db '+db_file+' -query '+query_fasta_file+' -out ' + ofile
		cmd = '"' + self.bin_path + '/' + tool + '" ' + option + ' ' + query_fasta_file + ' -db "'+db_file+'" -id 0.3 -blast6out "' + ofile + '" -maxaccepts 16'
		
		#print ('[USearch.runUsingFASTA2] ' + cmd)
		
		#os.system(cmd)

		# 출력 내용을 안 보여줌.
		if DISPLAY == False:
			cmd += ' >/dev/null 2>&1'

		os.system(cmd  )

		#MyUtil.executeDosCommand2(cmd, printout=False)
		#MyUtil.executeDosCommand3(cmd)

		#os.popen4(cmd)

		ret_values = self.parseOutput2(ofile)
		
		self.remove(ofile)

		return ret_values
	


	def runUsingFASTA(self, query_fasta_file, db_file, tool = 'usearch11.0.667_win32.exe', option = '-usearch_local'):

		# blast를 돌리고
		# 결과 파일 이름을 넘겨준다.
		
		#option = '-usearch_local'
		
		
		ofile = './TEMP/usearch_output_' + MyUtil.getRandomString(30) + '.txt'

		#cmd = self.bin_path + '\\blastp -db '+db_file+' -query '+query_fasta_file+' -out ' + ofile
		cmd = '"' + self.bin_path + '/' + tool + '" ' + option + ' ' + query_fasta_file + ' -db "'+db_file+'" -id 0.3 -blast6out "' + ofile + '" -maxaccepts 16'
		# -evalue 1e-3
		#cmd = '"' + self.bin_path + '/' + tool + '" ' + option + ' ' + query_fasta_file + ' -db "'+db_file+'" -evalue 1e-3 -blast6out ' + ofile 

		
		#print ('[USearch.runUsingFASTA] ' + cmd)
		
		#os.system(cmd)

		
		
		#MyUtil.executeDosCommand2(cmd, printout=False)
		#MyUtil.executeDosCommand3(cmd)
		# 출력 내용을 안 보여줌.
		if DISPLAY == False:
			cmd += ' >/dev/null 2>&1'
		os.system(cmd )

		#os.popen4(cmd)

		ret_values = self.parseOutput(ofile)

		self.remove(ofile)

		return ret_values


	def makeDB(self, fasta_file):

		# fasta 파일을 db로 만든다.
		# 만들어진 파일을 리턴해준다.
		# 우선 fasta를 복사해서
		# db 폴더로 옮긴다.




		#db_file = self.db_path + '\\' + fasta_file[ : fasta_file.rfind('.')] + '.db'

		# 파일 이름만 추출한다.
		base = os.path.basename(fasta_file)


		db_file = self.db_path + '/' + base[ : base.rfind('.')] + '.db'

		#cmd = '"' + self.bin_path+'/usearch8.1.1861_win32.exe" -makeudb_usearch "' + fasta_file + '" -output "' + db_file + '"'
		cmd = '"' + self.bin_path + '/usearch11.0.667_win32.exe" -makeudb_usearch "' + fasta_file + '" -output "' + db_file + '"'

		if MyUtil.which_os()=='linux':
			#cmd = '"' + self.bin_path + '/usearch8.1.1861_i86linux32" -makeudb_usearch "' + fasta_file + '" -output "' + db_file + '"'
			cmd = '"' + self.bin_path + '/usearch11.0.667_i86linux32" -makeudb_usearch "' + fasta_file + '" -output "' + db_file + '"'
		
		#print ('[USearch.makeDB] ' + cmd)

		# 실행이 끝날때까지 기다려야 하는데. 나중에 테스트 한다.
		#os.system(cmd)
		#os.popen4(cmd)
		#a=raw_input('sdsss')
		
		#MyUtil.executeDosCommand2(cmd)
		#MyUtil.executeDosCommand3(cmd)
		#print ('-------------------------------')
		if DISPLAY == False:
			cmd +=  ' >/dev/null 2>&1'
		os.system(cmd )
		
		return db_file 



#test()

if __name__ == '__main__': 

	#global DEBUG
	#DEBUG = True
	
	seq = '''
	MDISKVDSTRALRLLDKRLGCADERQIIMAGIERAEFIFRTIFRGLACTVVLGIIYISASSEPTLMYPT
	WIPWNWKDSTSAYLATAMLHTTALMANATLVLNLSSYPGTYLILVSVHTKALALRVSKLGYGAPLPAVRMQAILVGYIHD
	HQIILRLFKSLERSLSMTCFLQFFSTACAQCTICYFLLFGNVGIMRFMNMLFLL'''
	
	seq1 = 'ATCGTGTACGACTgttatcgtacgtacgt'
	seq2 = 'GGCAGAAGAA AATGTTCAGC TTGAGCAACT TCTGGCGCTCTCTTCGCCGA TTTACCGCCT TGAGCTGGCG ATGGTCGGGC'
	
	fasta_db_file= './temp/h.fa'
	f=open(fasta_db_file, 'w')
	f.write('>AA\n'+seq+'\n')
	f.write('>BB\n' + seq1 + '\n')
	f.write('>CC\n' + seq2 + '\n')
	f.close()

	
	b = USEARCH()
	#db_file= 'D:\\Work\\PythonLib\\NCBI\\db\\ecoli'
	#fasta_db_file= 'r:\\db.fa'
	# db conversion
	db_file = b.makeDB(fasta_db_file)
	
	
	print ('---------')
	
	
	#d = b.runBLASTUsingFASTA('R\\test2.txt', db_file)


	# word size는 query seq보다 최소 반보다 작아야한다.
	
	#d = b.runBLASTUsingSequence( seq, db_file, tool= 'blastn', option=option)
	#print repr(d)

	#c = b.parserroOutput('c:/NCBI/blast-2.2.26+/test2out.txt')
	#print repr(c)
	print ('Searching...')
	
	
	c = b.runUsingSequence(seq, db_file)
	print (repr(c))



	r = {'a1': seq1, 'a2':seq2, 'a3': seq}
	c = b.runUsingSequences(r, db_file)
	print (repr(c))
	
	
	print ('DONE')