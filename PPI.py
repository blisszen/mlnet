# -*- coding:cp949 -*-




# interaction network�� ������� �� protein ������ path�� ã���ش�.
import sys
import random
import copy
#import ShortestPath
import networkx as nx
import os
#import gzip
import shutil







class DirectedPPIDB():


	outgoing_interDB = {}
	incoming_interDB = {}

	#outgoing_scoreDB = {}
	#incoming_scoreDB= {}

	whole_gene_list = []

	score_DB = {}

	def __init__(self):
		self.outgoing_interDB = {}
		self.incoming_interDB = {}
		#self.outgoing_scoreDB = {}
		#self.incoming_scoreDB = {}
		self.score_DB = {}

		self.whole_gene_list = []


	'''
	def saveInteractionIntoFile(self, output_file):


		f = open(output_file,'w')

		for g in self.interDB.keys():

			s = g + '\t'

			for ig in self.interDB[g]:
				key = ig+':'+g
				ss = s + ig + '\t' + str( self.scoreDB[key]).strip()
				f.write(ss+'\n')

		f.close()
	'''


	def getWholeGeneList(self):
		return self.whole_gene_list

	'''
	def __makeIndex(self):


		r = {}
		for g in self.interDB.keys():
			r[g] = []
			for index in range( len( self.interDB[g])):
				r[g].append( index )

		return r

	'''


	'''
	def randomizeNetwork(self):
		# ���� ����� ���� �������� �� �׷���.
		# ���⼭�� �����ڸ� �ϳ��� �������� ����.
		# �׷��� ����� ���� �ξ� �پ���.
		genes = self.getWholeGeneList()

		for g in genes:

			iteration = len( self.interDB[g] )
			box = self.interDB[g]

			n1 = len(box) # network degree�� ���ϴ��� Ȯ���Ϸ��� �Ѵ�.

			for index in range(iteration):

				g1 = None
				g2 = None
				g3 = None
				g4 = None


				while(True):
					# partner�� ���ϸ� �ȴ�.

					g1 = g
					g2 = box[iteration - index - 1] # g ������ �����ϴ� �� �� �ϳ� �̴´�.
					g3, g4 = self.__pickOneInteraction()

					n2 = len(self.getPartners(g2))
					n3 = len(self.getPartners(g3))
					n4 = len(self.getPartners(g4))


					gA = self.getPartners(g1) + self.getPartners(g2)
					gB = self.getPartners(g3) + self.getPartners(g4)

					#gA = [g1, g2]
					#gB = [g3, g4]


					if (not g3 in gA) and (not g4 in gA) and (not g1 in gB) and (not g2 in gB):
						# ���� �ߺ��� ���� ��쿡�� �����Ѵ�.

						break


				# interaction�� ��ü�Ѵ�.
				self.__replaceInteraction(g1, g2, g3, g4)

				# randomize �� degree Ȯ���ϱ� ���� �����ϴ� ����
				_n1 = len( self.getPartners(g1))
				_n2 = len(self.getPartners(g2))
				_n3 = len(self.getPartners(g3))
				_n4 = len(self.getPartners(g4))

				if not ( n1 == _n1 and n2 == _n2 and n3 == _n3 and n4 == _n4 ):
					print 'Degree screw up! ', g1, g2, g3, g4
					sys.exit(1)


		# score�� �����ϰ� �ؾ��ϴµ�...
		self.score = {}
		for g in self.interDB.keys():
			for pg in self.getPartners(g):
				key1 = g+':'+pg
				key2 = pg + ':' + g
				ns = random.random()

				self.scoreDB[key1]=ns
				self.scoreDB[key2]=ns
	'''

	'''
	def __replaceInteraction(self, g1, g2, g3, g4):
		# g1�� g3�� ���� �ٲ۴�.


		# �켱 ���� interaction �����
		self.__delInteraction(g1, g2)
		self.__delInteraction(g3, g4)

		# ���ο� interaction �߰��Ѵ�.
		self.__addInteraction(g1, g3)
		self.__addInteraction(g2, g4)


	def __addInteraction(self, g1, g2):

		# intxn �߰��Ѵ�.


		self.interDB[g1].append(g2)
		self.interDB[g2].append(g1)

	def __delInteraction(self, g1, g2):
		# g1, g2���� interaction�� �����Ѵ�.
		self.__delGene2FromGene1(g1, g2)
		self.__delGene2FromGene1(g2, g1)

	def __delGene2FromGene1(self, g1, g2):
		# g1-g2 interaction���� g2 ����


		v = self.interDB[g1] # g1 �� ��Ʈ�ʸ� �޾ƿ���, �� �� g2 �����Ѵ�.
		for index in range( len(v)):
			if v[index] == g2:
				del v[index]
				break

	def __pickOneInteraction(self):


		gene1 = random.choice( self.interDB.keys() )
		gene2 = random.choice( self.interDB[gene1] )
		return gene1, gene2
	'''


	def getScore(self, gene_from, gene_to):

		gene1 = gene_from
		gene2 = gene_to

		if gene1==gene2:
			return 1.0
		else:
			key=gene1+':'+gene2
			return self.score_DB[ gene1+':'+gene2]

	'''
	def getIncomingScoreOf(self, gene_to, gene_from):
		# score�� �Ѱ��ش�.

		gene1 = gene_to.upper()
		gene2 = gene_from.upper()


		if gene1==gene2:
			return 1.0
		else:
			key=gene1+':'+gene2
			return self.incoming_scoreDB[ gene1+':'+gene2]


	def getOutgoingScoreOf(self, gene_from, gene_to):
		# score�� �Ѱ��ش�.

		gene1 = gene_from.upper()
		gene2 = gene_to.upper()

		if gene1==gene2:
			return 1.0
		else:
			key=gene1+':'+gene2
			return self.outgoing_scoreDB[ gene1+':'+gene2]

	'''

	def getIncomingPartnersOf(self, gene_id):

		#gene_id = gene_id.upper()

		if not gene_id in self.incoming_interDB:
			return []  # interaction�� ����. �׳� ���ٰ� �˷��ش�.
		else:
			return self.incoming_interDB[gene_id] # interacting partner�̴�.

	def getOutgoingPartnersOf(self, gene_id):

		#gene_id = gene_id.upper()

		if not gene_id in self.outgoing_interDB:
			return []  # interaction�� ����. �׳� ���ٰ� �˷��ش�.
		else:
			return self.outgoing_interDB[gene_id] # interacting partner�̴�.

	def __checkKeys(self, outgoing_r, incoming_r, gS, gT):

		if not gS in outgoing_r:
			outgoing_r[gS] = []
		if not gT in incoming_r:
			incoming_r[gT] = []

	'''
	def addInteractionFile( self, fname, index1, index2, index3, threshold  ):
		# �̹� interaction ������ �о��µ�,
		# �ٸ� ������ �о ���θ� ���������.


		# index3�� StringDB threshold�̴�.
		# -1�̸� �׳� �� �д´�.


		f=open(fname,'r')

		for s in f.readlines():
			s = s.replace('\n','').strip().upper()

			if len(s) > 0:
				x = s.split('\t')

				#print s
				g1 = x[index1]
				g2 = x[index2]
				s =  x[index3]

				if threshold != -1:
					# threshold�� �ִ�.
					if eval(s) < threshold:
						continue

				if g1 != g2: # self interaction�� ������� �ʴ´�.
					self.__checkKeys(self.interDB, g1, g2)

					if not g2 in self.interDB[g1]:
						self.interDB[g1].append(g2)

					if not g1 in self.interDB[g2]:
						self.interDB[g2].append(g1)



		f.close()
	'''

	def __addInteraction(self, outgoing_interDB, incoming_interDB, score_DB, gS, gT, score):

		self.__checkKeys(outgoing_interDB, incoming_interDB, gS, gT)

		if not gT in outgoing_interDB[gS]:
			outgoing_interDB[gS].append(gT)


		if not gS in incoming_interDB[gT]:
			incoming_interDB[gT].append(gS)


		score_DB[gS+':'+gT] = score



	def loadInteractionFile( self, fname, index1, index2, index3, index4, threshold  ):

		'''
		PPI ���� �д´�.

		input:  fname = �����̸�
				index1 = gene 1
				index2 = gene 2
				index3 = score
				index4 = direction ( > (g1->g2 �Ϲ�), < (g2<-g1 �Ϲ�) , <> (�����) )
				threshold = score ���� threshold ���ϸ� �����Ѵ�.
		'''



		# index3�� StringDB threshold�̴�.
		# -1�̸� �׳� �� �д´�.

		outgoing_r = {}
		incoming_r = {}
		#outgoing_c = {}
		#incoming_c = {}
		c = {}

		f=open(fname,'r', encoding='utf-8')

		# �Ѳ����� �� �������� �ʹ� ������.
		# ���پ� �е��� �����ؼ� �ӵ��� ������Ŵ.
		while(True):
			s = f.readline()
			if not s:
				break

			s = s.replace('\n','').strip()

			if len(s) > 0:
				x = s.split('\t')

				#print s

				#print repr(x)

				g1 = x[index1]
				g2 = x[index2]
				s = '1.0'  # conf ������ ������ ������ 1�̴�.
				dir = x[index4]

				if index3 != -1:
					s = x[index3]

				if threshold != -1:
					# threshold�� �ִ�.
					if eval(s) < threshold:
						continue

				if g1 != g2: # self-interaction�� ������� �ʴ´�.


					score = eval(s)

					# gS (source), gT (target) ����� �����Ѵ�.
					if dir == '>':
						# --> ����
						#self.__addInteraction(outgoing_r, incoming_r, outgoing_c, incoming_c, g1, g2, score)
						self.__addInteraction(outgoing_r, incoming_r, c, g1, g2, score)
					elif dir == '<':
						# <-- ����
						#self.__addInteraction(outgoing_r, incoming_r, outgoing_c, incoming_c, g2, g1, score)
						self.__addInteraction(outgoing_r, incoming_r, c, g2, g1, score)
					elif dir == '<>':
						# �����
						#self.__addInteraction(outgoing_r, incoming_r, outgoing_c, incoming_c, g1, g2, score)
						#self.__addInteraction(outgoing_r, incoming_r, outgoing_c, incoming_c, g2, g1, score)

						self.__addInteraction(outgoing_r, incoming_r, c, g1, g2, score)
						self.__addInteraction(outgoing_r, incoming_r, c, g2, g1, score)


					# ��ü gene list�� �ִ´�.
					if not g1 in self.whole_gene_list:
						self.whole_gene_list.append(g1)
					if not g2 in self.whole_gene_list:
						self.whole_gene_list.append(g2)

		f.close()

		# �̰� ����.
		self.incoming_interDB = incoming_r
		self.outgoing_interDB = outgoing_r
		#self.outgoing_scoreDB = outgoing_c
		#self.incoming_scoreDB = incoming_c
		self.score_DB = c



	def __trace(self, c_gene, final_gene, path_length, max_path_length, path):




		# ���� �����ߴ�.
		if c_gene == final_gene:
			return [ path + '\t' + c_gene ]

		# ���̸� �Ѿ� ��������...
		#print 'len = ', path_length
		if path_length > max_path_length:
			return None


		# ==============================================================
		# �������� path ������ �����Ѵ�.
		if len(path) != 0:
			new_path = path + '\t' + c_gene
		else:
			new_path = c_gene

		#print new_path

		new_path_length = path_length + 1
		# =============================================================


		partners = []
		partners = self.getOutgoingPartnersOf(c_gene)

		r = []


		for p in partners:

			if new_path.find(p)<0: # �̹� ������ �����ڴ� ������� �ʴ´�.

				pstrs = self.__trace( p, final_gene, new_path_length, max_path_length, new_path)
				if pstrs != None: # ���� ������ ��������װ�, �̰� string���� ���°� �ƴ϶� list�� ����´�.
					for ps in pstrs:
						r.append( ps )


		return r



	def findAllPaths(self, gene1, gene2, max_path_length):
		# �� �Լ��� ����ϱ� ���ؼ��� ���� loadInteractionFile(�����̸�, gene1 index, g2 index) �� ����������Ѵ�.
		paths = self.__trace(gene1, gene2, 1, max_path_length, '')
		return paths





























































class PPIDB_STRING(object):

	# networkx class �̿��Ѵ�.
	
	
	# STRING ���� ��, aliases, link �� ���� ������ �о PPI network�� ������ش�.
	G = None
	aliases = None

	def __init__(self):
		self.G = nx.Graph()
		self.aliases = {}


	def setEdgeAttribute(self, gene1, gene2, attr_key, attr_value):

		# edge attribute�� �����Ѵ�.
		# ���⼭�� ������ ���õǹǷ� gene1->gene2, gene2->gene1 ��� �����Ѵ�.
		self.G[gene1][gene2][attr_key] = attr_value
		self.G[gene2][gene1][attr_key] = attr_value

	def calculatePathwayDensity(self, path_list, attribute_key = 'confidence'):

		# �־��� pathway�� density�� ���Ѵ�.

		score = 0.0
		p_len = len(path_list)
		edges = 0.0

		for i in range( p_len - 1 ):
			for j in range( i+1, p_len):

				g1 = path_list[i]
				g2 = path_list[j]

				s = self.getScore(g1, g2, attribute_key)
				if s>0:
					edges += 1.0

				score = score + s # confidence�� �� ���Ѵ�.


		#den = score / float( p_len )
		den = score/ float( edges )
		return den


	def getWholeGeneList(self):
		# ��ü node �� ����Ʈ ���·� ��ȯ�Ѵ�.
		#return self.G.nodes()
		return list(self.G)


	def getStringIDOf(self, your_id):
		if your_id in self.aliases:
			return self.aliases[your_id]
		else:
			return None

	def getAliasesOf(self, your_id):
		#
		box = []
		
		if not your_id in self.aliases:
			return []
		
		string_id = self.aliases[your_id]
		
		for k in self.aliases:
			
			if self.aliases[k] == string_id:
				if not k in box:
					box.append(k)
				
		return box

	def getScore(self, gene1, gene2, attribute_key = None):
		# score�� �Ѱ��ش�.

		if attribute_key == None:
			attribute_key = 'confidence'

		try:
			gene11 = self.__convertID(gene1)
			gene22 = self.__convertID(gene2)
			
			
			return self.G[gene11][gene22][ attribute_key ]
		except:
			# gene1-gene2 ���谡 ����.
			return -1


	def getPartners(self, gene_id):
		
		
		try:
			gene = self.__convertID(gene_id)
			return list(self.G.neighbors(gene))
		except:
			return []


	
	def __convertID(self, your_id):
		if your_id in self.aliases:
			return self.aliases[your_id]
		else:
			return None
		
	
	
	def load( self, interaction_file, index1, index2, index3, threshold, delimiter, alias_file  ):
		
		if alias_file is not None:
			#print ('[PPI] Loading aliases...')
			self.aliases = self.__getAliases(alias_file)
			
		#print ('[PPI] Loading interactions...')
		self.__loadInteraction(interaction_file, index1, index2, index3, threshold, delimiter)
		#print '    ---> done'
		
	def __getAliases(self, alias_file):
		box = {}
		f=open(alias_file, 'r', encoding='utf-8')
		
		init = True
		
		while(True):
			line = f.readline(10000)
			if not line:
				break
			
			line = line.strip()
			if len(line) == 0: continue
			
			if init:
				init = False
				continue
			
			content = line.split('\t')
			string_id = content[0]
			alias_id = content[1]
			alias_db = content[2]
			
			box[alias_id]=string_id
			box[string_id]=string_id
			
			if alias_id in box:
				if box[alias_id] != string_id:
					print ('[PPI, INCONSISTENT ID]', alias_id, alias_db)
			
		
		f.close()
		
		return box


	def __loadInteraction( self, interaction_file, index1, index2, index3, threshold, delimiter, ignore_first_line = True ):



		IGNORE_UNKNOW_ID = True

		cnt = 0

		f=open(interaction_file,'r', encoding='utf-8')
		
		init = ignore_first_line

		# �Ѳ����� �� �������� �ʹ� ������.
		# ���پ� �е��� �����ؼ� �ӵ��� ������Ŵ.
		
		cnt = 0.0
		cnt2 = 0.0
		
		total = os.path.getsize(interaction_file)
		
		while(True):
			#lines = f.readlines(5000)
			#if not lines:
			#	break

			#for s in lines:
			
			s = f.readline()
			if not s:
				break
			else:
				cnt2 += 1.0
				
				
				
				
				s = s.replace('\n','').strip()
			
	
				if len(s) > 0:
					
					if init:
						init = False
						continue
					
					#x = s.split('\t')
					x = s.split(delimiter)
	
					m = max( [index1, index2, index3])
					if len(x) < m + 1:
						print ('[PPI, Error] ', s)
						continue
	
					g1 = x[index1].strip()
					g2 = x[index2].strip()
					
					v = '1.0'  # conf ������ ������ ������ 1�̴�.
					if index3 != -1:
						v =  x[index3].strip()
						
						# for string
						if eval(v) > 1.0:
							v = '0.' + v
						
						
	
					#print v, threshold
					if threshold != -1:
						#
						
						# threshold�� �ִ�.
						if eval(v) < threshold:
							continue
						
						
						
					
	
					if g1 != g2: # self-interaction�� ������� �ʴ´�.
	
						if IGNORE_UNKNOW_ID:
							try:
								x = self.aliases[g1]
								x = self.aliases[g2]
							except:
								# unknown ID
								#print '\t', g1, g2
								continue


						if not g1 in self.aliases :
							self.aliases[g1] = g1
						if not g2 in self.aliases:
							self.aliases[g2] = g2
						
						sg1 = self.aliases[g1]
						sg2 = self.aliases[g2]
						
						self.G.add_node( sg1 )
						self.G.add_node( sg2 )
						#self.G.add_edge( sg1, sg2, { 'confidence': eval(v), 'distance':1.0-eval(v)}   )
						self.G.add_edge( sg1, sg2 )
						self.G[sg1][sg2].update( { 'confidence': eval(v), 'distance':1.0-eval(v)}   )
						
						cnt += 1
	
		f.close()
		
		#print ('[PPI] Loaded ', cnt, 'interactions from ', interaction_file)



	def findAllPaths(self, gene1, gene2, max_path_length):
		# �� �Լ��� ����ϱ� ���ؼ��� ���� loadInteractionFile(�����̸�, gene1 index, g2 index) �� ����������Ѵ�.

		gene11 = self.__convertID(gene1)
		gene22 = self.__convertID(gene2)
		return nx.dijkstra_path(self.G, gene11, gene22, 'distance')


	def getSTRING_ID(self, your_id):
		return self.__convertID(your_id)

	def findAllShortestPaths(self, gene1, gene2, attribute=None):

		# �ִ� �Ÿ��� ������, �� node�� ���� list�� �����ϰ�
		# ������ None�� ���ϵ�.
		# �Ÿ� �Ի��� attribute�� �̿���. distance���� �⺻�̸� �̴� 1-confidence�̴�.
		# setEdgeAttribute() �Լ��� �̿��� attribute�� ������ �� ������,
		# ���� �Լ����� attribute = 'xxx' ��� �������ָ� �װ� �̿��ؼ� ã�� �� �ִ�.
		# ��� shortest path�� ����Ʈ�� ��Ƽ� �����Ѵ�.

		path = None

		#if attribute is None:
		#	attribute = 'distance'

		try:
			gene11 = self.__convertID(gene1)
			gene22 = self.__convertID(gene2)
			
			path = nx.all_shortest_paths(self.G, source=gene11, target=gene22, weight = attribute)


			#path = nx.dijkstra_path(self.G, gene1, gene2, attribute)
			#path = nx.all_shortest_paths(self.G, gene1, gene2, attribute)

			path = [p for p in path]


		except:
			#print 'Nodes not connected'
			path = None


		return path



	def findShortestPath(self, gene1, gene2, attribute=None):

		# �ִ� �Ÿ��� ������, �� node�� ���� list�� �����ϰ�
		# ������ None�� ���ϵ�.
		# �Ÿ� �Ի��� attribute�� �̿���. distance���� �⺻�̸� �̴� 1-confidence�̴�.
		# setEdgeAttribute() �Լ��� �̿��� attribute�� ������ �� ������,
		# ���� �Լ����� attribute = 'xxx' ��� �������ָ� �װ� �̿��ؼ� ã�� �� �ִ�.

		path = None

		#if attribute is None:
		#	attribute = 'distance'

		try:
			
			gene11 = self.__convertID(gene1)
			gene22 = self.__convertID(gene2)			
			
			#path = nx.astar_path(self.G, gene1, gene2, attribute)
			# a star algorithm�� ��Ȯ������ �ʴ�. �޸���ƽ�̶� �׷���?
			path = nx.dijkstra_path(self.G, gene11, gene22, attribute)
			#path = nx.shortest_path(self.G, gene1, gene2, attribute)
		except:
			#print 'Nodes not connected'
			path = None

		return path


	def getAverageShortestPathLength(self, fname):

		
		#return nx.average_shortest_path_length(self.G, 'confidence')
		return nx.average_shortest_path_length(self.G) # �׳� �ܼ� �Ÿ��� ����Ѵ�. confidence ����.


	def randomizeNetwork(self):
		
		'''
		return a randomized PPIDB_STRING() class
		'''

		id_map = {}
		degree = {}
		ppi = PPIDB_STRING()
		
		
		
		for g in self.getWholeGeneList():
			
			id_map[g] = None
			
			deg = len( self.getPartners(g))
			deg_str = str(deg).strip()
			if not deg_str in degree:
				degree[deg_str] = []
			degree[deg_str].append(g)
			
		# randomize
		for d in degree:
			app = copy.deepcopy( degree[d] )
			random.shuffle(app)
			
			for i in range(len(app)):
				before = degree[d][i]
				after = app[i]
				id_map[before] = after
				

		for g in self.getWholeGeneList():
			
			for p in self.getPartners(g):
							
				g1 = id_map[g]
				g2 = id_map[p]
				v = self.getScore(g, p)
				
				if g1 != g2:
					ppi.G.add_node( g1 )
					ppi.G.add_node( g2 )
					ppi.G.add_edge( g1, g2, { 'confidence': v, 'distance':1.0-v}   )
					
				
			
		
		return ppi
		
	def randomizeAndSaveNetwork(self,fname, gzip=False):
		
		id_map = {}
		degree = {}
		
		
		for g in self.getWholeGeneList():
			
			id_map[g] = None
			
			deg = len( self.getPartners(g))
			deg_str = str(deg).strip()
			if not deg_str in degree:
				degree[deg_str] = []
			degree[deg_str].append(g)
			
		# randomize
		for d in degree:
			app = copy.deepcopy( degree[d] )
			random.shuffle(app)
			
			for i in range(len(app)):
				before = degree[d][i]
				after = app[i]
				id_map[before] = after
				
		# save into a file
		f=open(fname, 'w', encoding='utf-8')
		f.write('gene1 gene2 conf\n')
		for g in self.getWholeGeneList():
			
			for p in self.getPartners(g):
				
				if id_map[g] != id_map[p]:
					txt = id_map[g] + ' ' + id_map[p] + ' ' + str(int( 1000 * self.getScore(g, p)))
					f.write(txt+'\n')
			
		
		f.close()

		# zip the file
		if gzip:
			with open(fname, 'rb') as f_in, gzip.open(fname+'.gz', 'wb') as f_out:
				shutil.copyfileobj(f_in, f_out)




		








def test4():


	string_db_file = 'C:/Users/dkna/Desktop/workspace/Gsponer/NetworkAnalysis/Predictor/ppi/string_interaction_db(yeast).txt'
	index1 = 0
	index2 = 1
	score_index = 2
	threshold = 0.4 #0.4 # -1


	'''
	string_db_file = 'C:/Users/dkna/Desktop/workspace/Gsponer/NetworkAnalysis/Predictor/ppi/biogrid_converted.txt'
	index1 = 0
	index2 = 2
	score_index = -1
	threshold = -1
	'''

	'''
	string_db_file = 'C:/Users/dkna/Desktop/workspace/Gsponer/NetworkAnalysis/Predictor/ppi/yeast_interaction(sgd_id).txt'
	index1 = 0
	index2 = 3
	score_index = 6
	threshold = -1
	'''


	ppi = PPIDB2()
	ppi.loadInteractionFile(string_db_file, index1, index2, score_index, threshold)


	#g1 = 'S000001661' # ste3
	#g2 = 'S000001126' # ste12


	#g1 = 'S000006369' # rho1
	#g2 = 'S000006010' # rlm1


	'''
	g1 = 'S000005042' # Ras2
	g2 = 'S000001126' # STE12

	g1 = 'STE3'
	g2 = 'STE12'
	'''
	g1 = 'S000001409' # sln1
	g2 = 'S000004103' # hog1


	#g1 = 'S000000920' # sho1
	#g2 = 'S000000287' # tec1

	path = ppi.findShortestPath(g1, g2)
	print (path)

	# pathway Ȯ��
	path2 = expandPathway(ppi, path)
	print ('\n'.join( path2 ))

def test3():


	string_db_file = 'C:/Users/dkna/Desktop/workspace/Gsponer/NetworkAnalysis/Predictor/ppi/string_interaction_db(fly).txt'
	index1 = 0
	index2 = 1
	score_index = 2
	threshold = 0.8 #0.4 # -1

	'''

	string_db_file = 'C:/Users/dkna/Desktop/workspace/Gsponer/NetworkAnalysis/Predictor/ppi/biogrid_converted.txt'
	index1 = 0
	index2 = 2
	score_index = -1
	threshold = -1
	'''

	ppi = PPIDB2()
	ppi.loadInteractionFile(string_db_file, index1, index2, score_index, threshold)

	# test pathway
	# tsl = FBgn0003867
	# tor = FBgn0021796
	# phl = FBgn0003079
	# http://www.kegg.jp/kegg/pathway/dme/dme04013.html

	#g1 = 'FBgn0004638'
	#g2 = 'FBgn0003079'

	g1 = 'FBgn0005648' # pAbp2
	g2 = 'FBgn0041188' # Atx2
	path = ppi.findShortestPath(g1, g2)
	print ('******', path)


	'''
	g1 = 'FBgn0001137' # Grk
	g2 = 'FBgn0003205' # Ras85D

	path = ppi.findShortestPath(g1, g2)
	print '\n'.join(path)
	 '''



def test2():

	ppi = PPIDB2()
	ppi.loadInteractionFile("ppi_test.txt", 0, 1, 2, -1)
	print (ppi.getWholeGeneList())
	print (ppi.getScore('A', 'B'))
	print (ppi.getScore('B', 'A'))

def test():

	dppi = DirectedPPIDB()
	dppi.loadInteractionFile("ppi_test.txt", 0, 1, 2, 3, -1)

	a = 'A'
	for g in dppi.getIncomingPartnersOf(a):
		print (a, '<-', g, 'with', dppi.getScoreOf(g, a))

	for g in dppi.getOutgoingPartnersOf(a):
		print (a, '->', g,'with', dppi.getScoreOf(a, g))

	a = 'C'
	for g in dppi.getIncomingPartnersOf(a):
		print (a, '<-', g, 'with', dppi.getScoreOf(g, a))


	for g in dppi.getOutgoingPartnersOf(a):
		print (a, '->', g,'with', dppi.getScoreOf(a, g))




if __name__ == '__main__':
	test3()
