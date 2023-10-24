# -*- coding: utf-8 -*-


'''
Integrated rank ratios predicted using various information

'''

import sys
import statistics
import math
import _0_Preprocess
import _100_NetworkExpansionModel
import os


from datetime import datetime
import time

class IntegrateRanks():

	files = []
	loaded_ranks = {}
	ordered_ranks = {}
	order = []
	info = None
	last_step = False
	cutoff = None
	fisher = False


	def __init__(self):
		self.files = []
		self.loaded_ranks = {}
		self.ordered_ranks = {}
		self.order = []
		self.info = None
		self.last_step = False
		self.cutoff = None
		self.fisher = False
		self.cutoff = 1


	def __load(self, files_list, exclude_set):


		# 여기는 정보가 하나라도 누락되면 아예 없애버린다.
		'''
        ranks = self.__loadFiles(files_list)

        r_ranks = {}

        for gid in ranks[0].keys():

            r = [  ranks[0][gid]   ]

            for i in range( 1, len(ranks)):
                if ranks[i].has_key(gid):
                    r.append( ranks[i][gid] )
                else:
                    break

            if len(r) == len(ranks):
                r_ranks[gid] = r

        return r_ranks
        '''

		# 여기는 없으면 없는대로 아쉬운대로 rank를 만든다.
		#ranks, aucs = self.__loadFiles(files_list, exclude_set)
		ranks = self.__loadFiles(files_list, exclude_set)

		r_ranks = {}


		for gid in self.__getWholeGenes(ranks):

			r = [    ]

			for i in range( 0, len(ranks)):

				#auc = aucs[i]

				if gid in ranks[i]:
					#r.append( [ ranks[i][gid], auc ] )
					r.append( ranks[i][gid] )

			r_ranks[gid] = r

		return r_ranks




	def __getWeights3 (self, rank_ratios):
		r = {}        
		if self.fisher == False:
			# use weights

			for g in rank_ratios:

				score = 0.0
				for ratio, auc in rank_ratios[g]:
					s = (1.0-ratio)*auc
					if s>score:
						score = s

				r[g] = score 

		return r

	def __getExcludeSet(self, ranks):

		ex = []

		for gid in ranks:
			if len( ranks[gid] ) < self.cutoff:
				ex.append( gid )

		return ex


	def __getWholeGenes(self, ranks):

		r = {}
		for rank in ranks:
			for gid in rank:
				if not gid in r:
					r[gid]=None

		return list(r)

	def __getPvalues(self, rank_ratios):
		'''
		r = {}        
		if self.fisher == False:
		    # use order statistic
		    O = statistics.OrderStatistic()

		    for g in rank_ratios.keys():
		        r[g] = O.rankUsingOrderStatistics( rank_ratios[g] )
		else:
		    # use fisher's score
		    for g in rank_ratios.keys():
		        v = 0.0
		        for pv in rank_ratios[g]:
		            v = v - 2*math.log(pv)
		        r[g] = v

		return r
		'''



		r = {}        
		if self.fisher == False:
			# use order statistic
			O = statistics.OrderStatistic()

			for g in rank_ratios:
				q = []
				#for ratio, auc in rank_ratios[g]:
				for ratio in rank_ratios[g]:
					q.append(ratio)
					# rank만 저장해서 계산한다.

				r[g] = O.rankUsingOrderStatistics( q )  # rank_ratios[g][1]은 auc        else:

		return r    

	def __order(self, pvalues):

		'''
        Returns a tuple of ordered (genes, pvalue) from top to bottom  
        '''
		r = []
		box3 = sorted( pvalues.items(), key=lambda item: item[1], reverse = False) # smaller is better
		return box3

	def __getWeights (self, rank_ratios):
		r = {}        
		if self.fisher == False:
			# use weights


			for g in rank_ratios:


				score = 0.0
				for ratio, auc in rank_ratios[g]:
					score = score + (1.0-ratio)*auc

				r[g] = score 

		return r

	def __write(self, tup, output):


		f=open(output,'w')
		f.write('GeneID\tp-value\n')

		for gid, pvalue in tup:
			s = gid + '\t' + str(pvalue).strip() + '\t'

			'''
            gg = gid.replace('m_', '').replace('x_','')


            if info.has_key(gg):
                s = s + info[gg]
            else:
                s = s + 'not tested ever'
            '''

			f.write(s + '\n')

		f.close()


	def run(self, files_list, output):
		'''
		Returns a tuple of ordered (gene, pvalue) from top rank to bottom 
		'''

		# load calculated rank ratios
		self.loaded_ranks = self.__load(files_list, [])


		# feature 수가 적은 것은 제거한다.
		ex_set = self.__getExcludeSet(self.loaded_ranks)





		# feature 감안해서 새로 읽는다.
		self.loaded_ranks = self.__load(files_list, ex_set)



		# calculate order statistic
		self.ordered_ranks = self.__getPvalues( self.loaded_ranks  )
		#self.ordered_ranks = self.__getWeights( self.loaded_ranks  )

		# order
		self.order = self.__order(  self.ordered_ranks )


		# save into a file
		self.__write(self.order, output)



		return self.order


	def __loadFiles(self, files_list, exclude_set):


		global RANK_METHOD

		r = []
		a = []

		for f in files_list:

			#auc = self.__getAUC(f)
			#a.append( auc )

			r.append( self.__loadFile(f, exclude_set)  )


		return r #, a


	def loadRankRatio(self, fname):
		# Rank ratio를 구하고자 하는 사람들을 위해 이 함수를 외부로 노출시킴
		return self.__loadFile(fname, [])

	def __loadFile(self, fname, exclude_set):
		'''
		Reads xxx_gene.txt containing gene id and its rank
		Returns rank
		'''

		fname = fname.replace("'", '')

		r = {}

		f=open(fname,'r')

		init = True


		rank = 1.0
		counter = 0.0


		#print ('Loading ', fname )

		prev_score = -10000000



		for s in f.readlines():



			s = s.replace('\n','')
			if len(s) == 0:
				continue

			# 주석 제거함
			if s[0] == '#' or s[0] == '!': continue


			if init:
				init=False
				continue

			# remove a tag such as m_ or x_
			if self.last_step:
				s = s.replace('m_', '').replace( 'x_', '')

			# get score and rank
			x = s.split('\t')

			gid = x[0].strip()
			score = eval( x[1] )

			# exclude_set에 있으면 제거한다.
			if gid in exclude_set:
				continue

			if self.fisher:
				# get Fisher's score from p-values
				r[gid] = score
			else:


				if score != prev_score:
					prev_score = score
					rank = counter

					r[gid]=rank + 1.0
				else:
					r[gid]=rank + 1.0

				#print gid,rank

			counter = counter + 1.0
		f.close()


		# normalize
		if self.fisher == False:
			for g in r:
				r[g]=r[g]/counter


		return r


def rewriteForPlot(ofiles):

	# just for plot

	for f1 in ofiles:
		for f2 in ofiles:
			if f1 != f2:
				d1r = getRankedScores(f1)
				d2r = getRankedScores(f2)

				d1 = f1.replace('__','').replace('.txt', '')
				d2 = f2.replace('__','').replace('.txt', '')

				output = '___'+d1+'_'+d2+'.txt'
				f=open(output,'w')
				f.write(d1+'\t'+d2+'\n')
				for g in d1r:
					if g in d2r:
						s = str(d1r[g]).strip() + '\t' + str(d2r[g]).strip()
						f.write(s+'\n')
				f.close()

def getRankedScores(of):

	r = {}

	f=open(of,'r')


	cnt = 0.0

	for s in f.readlines():
		s = s.replace('\n','')
		x = s.split('\t')

		cnt = cnt + 1.0

		gid = x[0]
		rank = cnt
		r[gid]=rank  


	# normalize
	for g in r:
		r[g]=r[g]/cnt
	f.close()

	return r






def run(iteration, rnd_id):

	FILES = _100_NetworkExpansionModel.getFiles(iteration, rnd_id)


	# 여기에 리스트 형태로 들어가 있는 질병에 대해서만 disease-specific modifiers를 예측한다.

	predict_disease_specific_modifiers = _0_Preprocess.PREDICT_DISEASE_MODIFIERS

	ofiles = []

	x, query_modifiers = _0_Preprocess.getModifiers()

	for disease in FILES:

		print (disease, '...')

		output = _0_Preprocess.OUTPUT_FOLDER + '/Integrated_'+disease+'_' + rnd_id + '.txt'
		ofiles.append(output)


		if not disease in predict_disease_specific_modifiers:

			print('[_10_Integration.py] disease-specific modifiers 예측 안하고 건너뜀: ', disease)
			# 대신 결과 파일은 query 리스트를 그대로 넣는다.
			f=open(output, 'w')
			f.write('GeneID	p-value\n')
			# p-value는 아무 의미 없음.
			cnt = 0
			for gene in query_modifiers[disease]:
				cnt += 1
				f.write(gene+'\t'+str(cnt).strip()+'\n') # score가 다르도록 넣어줘야 rank 계산이 제대로 됨.
			f.close()

			print('[[[ Query를 그냥 복사한 integrated file = ', output)
			# add common names
			addGeneNames(output)



		else:

			# disease-specific modifiers를 예측한다.

			flist = []
			for f in FILES[disease]:
				flist.append(f)



			ir = IntegrateRanks()
			ir.cutoff = 1 # feature 가 1개만 있으면 된다.

			ir.run(flist, output)

			#ofiles.append(output) # rank ratio 내용이 여기에 파일로 저장된다.

			print ('[[[ Integrated file = ', output)


			# add common names
			addGeneNames(output)



	print ('done')
	
	
	# 기존의 결과 파일은 모두 삭제하자.
	for f in FILES:
		delete(f)
		
	
	
	return ofiles


def addGeneNames(output):

	foname = output + '_genes.txt'
	fi = open(output, 'r')
	fo = open(foname, 'w')

	init = True
	for s in fi.readlines():
		s = s.strip()
		if len(s) == 0: continue

		if init:
			init = False
			fo.write(s + '\n')
		else:
			t = s.split('\t')
			n = t[0].replace('m_', '')
			c = _0_Preprocess.ID_BOX.getCommonNameOf(n)
			fo.write(s + '\t' + c + '\n')

	fi.close()
	fo.close()




def mergeTwoResult(fnam1, fname2, output):

	ir = IntegrateRanks()

	f1 = ir.loadRankRatio(fnam1)
	f2 = ir.loadRankRatio(fname2)


	f=open(output,'w')

	for g in f1:
		if g in f2:
			f.write(str(f1[g]).strip() + '\t' + str(f2[g]).strip()+'\n' )

	f.close()

def delete(fname):
	try:
		os.remove(fname)
	except:
		pass

def main(config, iteration, rnd_id):

	print ('''
========================================================    
    [10] Integration
========================================================          
''')


	_0_Preprocess.init(config)
	_100_NetworkExpansionModel.init(config)

	x = run(iteration, rnd_id)
	
	
	print ('''
========================================================    
    [10] Integration (End)
========================================================          
''')
	
	return x



if __name__ == '__main__':



	start_time = datetime.now()


	# python _1_GO.py config.txt test
	print ('''
========================================================    
    [10] Integration
========================================================          
''')

	if len(sys.argv) == 3:
		config_file = sys.argv[1]
		iteration = eval(sys.argv[2])

		_0_Preprocess.init(config_file)
		_100_NetworkExpansionModel.init(config_file)

		run(iteration=iteration)

	'''
    #fnam1 = './0. Integration/__ALL_fisher.txt'
    #fname2 = './0. Integration/__ALL_rank.txt'
    fnam1 = './0. Integration/f.txt'
    fname2 = './0. Integration/n.txt'
    output = './0. Integration/__ALL_merged.txt'

    mergeTwoResult(fnam1, fname2, output)
    '''

	print ('''
========================================================    
    [10] Integration (End)
========================================================          
''')


	end_time = datetime.now()
	print ('Elapsed time: {}'.format(end_time - start_time))
	print (time.ctime())