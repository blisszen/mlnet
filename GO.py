# -*- coding: ms949 -*-
'''

GO ������ �о ������ �ִ´�.
�̸� �̿��� �ٸ� ��꿡 ����.

'''

import DBPATH

EXE_PATH = DBPATH.DB_PATH +'/GO'

GO_FILE = EXE_PATH + '/gene_ontology.obo'
GO_EXT_FILE = EXE_PATH + '/gene_ontology_ext.obo'

class GOTerm(object):

	go_id = ''

	RELATION_IS_A = 'is_a'
	RELATION_PART_OF = 'part_of'
	RELATION_REGULATES = 'regulates'
	RELATION_NEGATIVELY_REGULATES = 'negatively_regulates'
	RELATION_POSITIVELY_REGULATES = 'positively_regulates'

	CATEGORY_PROCESS = 'biological_process'
	CATEGORY_FUNCTION = 'molecular_function'
	CATEGORY_COMPONENT = 'cellular_component'

	category = ''

	is_obsolete = False


	name = ''
	definition = ''
	name_space = ''  # process/function/structure ū �����̴�.


	relation = { }  # ���⿣ relation[GO ID] = IS_A  �̷� �� �� ����ȴ�. Parent node ���̴�.
	child_relation = {}

	alt_ids = []

	# =============
	count = {}
	flyids = {}
	flyids_cnt = {}


	def __init__(self):
		self.go_id = ''
		self.name = ''
		self.name_space = ''
		self.definition = ''
		self.category = ''

		self.relation = {}
		self.child_relation = {}

		self.count = {}
		self.flyids = {}
		self.flyids_cnt = {}

		self.is_obsolete = False
		self.alt_ids = []



	def setCategory(self, category):
		self.category = category

	def getCategory(self):
		return self.category

	def addAlternativeID(self, aid):
		self.alt_ids.append(aid)
	def getAlternativeIDs(self):
		return self.alt_ids

	def setObsolete(self, obsolete):
		self.is_obsolete = obsolete
	def isObsolete(self):
		return self.is_obsolete

	def setGOID( self, go_id ):
		self.go_id = go_id
	def getGOID(self):
		return self.go_id

	def setGOName(self, name):
		self.name = name
	def getGOName(self): return self.name
	def setGODefinition(self, definition):
		self.definition = definition
	def getGODefinition(self):
		return self.definition

	def setGONameSpace(self, namespace):
		self.name_space = namespace
	def getGONameSpace(self):
		return self.name_space

	def addParentRelationship(self, go_id, relation):
		self.relation [ go_id ] = relation
	def getParentRelationshipSize(self):
		''' ����� relationship ��ü ������ �Ѱ��ش� '''
		return len( self.relation )


	def getParentRelationshipGOIDs(self):
		''' relationship �� ����� ��� GO id�� �Ѱ��ش�. List '''
		return list( self.relation )

	def getParentRelationshipGOIDsOfISA(self):
		# is_a ���迡 �ִ� parent id�� �޾ƿ´�.
		r=[]
		for k in self.relation:
			if self.relation[k] == self.RELATION_IS_A :
				r.append(k)
		return r


	def getParentRelationshipGOIDsOfREGULATE(self):
		r=[]
		for k in self.relation:
			if self.relation[k] == self.RELATION_REGULATES :
				r.append(k)
		return r

	def getParentRelationshipGOIDsOfNEGATIVELYREGULATE(self):
		r=[]
		for k in self.relation:
			if self.relation[k] == self.RELATION_NEGATIVELY_REGULATES :
				r.append(k)
		return r


	def getParentRelationshipGOIDsOfPOSITIVELYREGULATE(self):
		r=[]
		for k in self.relation:
			if self.relation[k] == self.RELATION_POSITIVELY_REGULATES :
				r.append(k)
		return r





	def getChildRelationshipGOIDsOfISA(self):
		# is_a ���迡 �ִ� parent id�� �޾ƿ´�.
		r=[]
		for k in self.child_relation:
			if self.child_relation[k] == self.RELATION_IS_A :
				r.append(k)
		return r


	def getChildRelationshipGOIDsOfPARTOF(self):
		# is_a ���迡 �ִ� parent id�� �޾ƿ´�.
		r=[]
		for k in self.child_relation:
			if self.child_relation[k] == self.RELATION_PART_OF :
				r.append(k)
		return r


	def getChildRelationshipGOIDsOfREGULATE(self):
		# is_a ���迡 �ִ� parent id�� �޾ƿ´�.
		r=[]
		for k in self.child_relation:
			if self.child_relation[k] == self.RELATION_REGULATES :
				r.append(k)
		return r


	def getChildRelationshipGOIDsOfNEGATIVELYREGULATE(self):
		# is_a ���迡 �ִ� parent id�� �޾ƿ´�.
		r=[]
		for k in self.child_relation:
			if self.child_relation[k] == self.RELATION_NEGATIVELY_REGULATES :
				r.append(k)
		return r

	def getChildRelationshipGOIDsOfPOSITIVELYREGULATE(self):
		# is_a ���迡 �ִ� parent id�� �޾ƿ´�.
		r=[]
		for k in self.child_relation:
			if self.child_relation[k] == self.RELATION_POSITIVELY_REGULATES :
				r.append(k)
		return r

	def addChildRelationship(self, go_id, relation):
		self.child_relation [ go_id ] = relation
	def getChildRelationshipSize(self):
		''' ����� relationship ��ü ������ �Ѱ��ش� '''
		return len( self.child_relation )
	def getChildRelationshipGOIDs(self):
		''' relationship �� ����� ��� GO id�� �Ѱ��ش�. List '''
		return list(self.child_relation)


	def getParentRelationshipOf(self, go_id):
		''' GO ID���� �־��ָ� �ش��ϴ� relationship ���� �����ش� '''
		return self.relation[go_id]

	def getChildRelationshipOf(self, go_id):
		''' GO ID���� �־��ָ� �ش��ϴ� relationship ���� �����ش� '''
		return self.child_relation[go_id]


	def toString(self):
		r = 'ID:\t' + self.getGOID() + '\t' + \
		  'Name:\t' + self.getGOName() + '\t' + \
		  'NameSpace\t'+ self.getGONameSpace() + '\t' + \
		  'Definition\t' + self.getGODefinition() + '\t' + \
		  'Relation\t'

		for ids in self.getParentRelationshipGOIDs():
			r = r + '\t' + ids + '\t' + self.getParentRelationshipOf(ids) + '\n'

		for k in self.count:
			r = r + k + '\t' + str( self.count[k] ) + '\n'

		r = r + '====================================='
		return r


	''' �Ʒ� ���� �м��� ���� Ư���� �߰��� �Լ��� '''
	def increasePoint(self, cnt_id, point, flyDB_id):
		if not cnt_id in self.count:
			self.count [cnt_id ] = 0

		if not cnt_id in self.flyids:
			self.flyids[cnt_id] = []
			self.flyids_cnt[cnt_id] = 0

		self.count [cnt_id] = self.count [cnt_id] + point

		if not flyDB_id in self.flyids[cnt_id]:
			self.flyids[cnt_id].append( flyDB_id )
			self.flyids_cnt[cnt_id] = self.flyids_cnt[cnt_id] + 1


	def getCount(self, cnt_id):

		if cnt_id in self.count :
			return self.count [cnt_id]
		else:
			return 0.0

	def getFlyDBIDSizeOf(self, cnt_id):
		if cnt_id in self.flyids_cnt:
			return self.flyids_cnt[cnt_id]
		else:
			return 0

	def getFlyDBIDs(self, cnt_id):
		if cnt_id in self.flyids:
			#print repr( self.flyids )

			return self.flyids[cnt_id]
		else:
			return []




	# =====================================

class GO(object):

	go_terms = {}

	def __init__(self, obo_file = GO_EXT_FILE):
		self.go_terms = {}
		self.loadFile(obo_file)
		self.__setCategory() # ����� GO ID�� �ִ�. �̰� category�� �ִ´�.

	def __setCategory(self):
		# molecular function: GO:0003674
		# biological process: GO:0008150
		# cellular component: GO:0005575

		for go_id in self.go_terms:
			d1 = self.getDistance(go_id, 'GO:0003674') # function
			d2 = self.getDistance(go_id, 'GO:0008150') # process
			d3 = self.getDistance(go_id, 'GO:0005575') # component

			#print go_id, d1, d2, d3

			if d1>=0:
				#self.go_terms[go_id].setCategory( self.go_terms[go_id].CATEGORY_FUNCTION  )
				self.go_terms[go_id].setCategory( 'F'  )
			elif d2>=0:
				#self.go_terms[go_id].setCategory( self.go_terms[go_id].CATEGORY_PROCESS  )
				self.go_terms[go_id].setCategory( 'P' )
			elif d3>=0:
				self.go_terms[go_id].setCategory( self.go_terms[go_id].CATEGORY_COMPONENT  )
				self.go_terms[go_id].setCategory( 'C'  )




	def size(self):
		''' ����� GO_TERM ������ �Ѱ��ش� '''
		return len(self.go_terms)

	def getGOTermOf(self, go_id):
		''' GO ID�� ������ �ش��ϴ� GOTERM class�� �Ѱ��ش� '''
		if go_id in self.go_terms:
			return self.go_terms[go_id]
		else:
			# ���� �� �� ���⼭ ã��.
			# print '[GO class] unknown GO ID = ', go_id

			# alternative ID�� ���� �ִ�.
			for k in self.go_terms:
				aids = self.go_terms[k].getAlternativeIDs()
				if go_id in aids:
					return self.go_terms[k]

			return None

	def loadFile(self, fname):
		''' GO ������ �д´�. ������ FlyBase OBO ���ϸ� �д´� '''


		print ('[GO class] Loading ', fname)

		f=open(fname,'r')

		flag = 0
		goterm = None

		for s in f.readlines():
			s=s.replace('\n','')

			#print s

			if len(s)>0:
				if s == '[Term]' or s == '[Typedef]':
					flag = 1
					if goterm != None:
						# ����Ȱ� ������ dict �� �ִ´�.
						#print goterm.getGOID()
						self.go_terms [ goterm.getGOID() ] = goterm

						#print goterm.getGOID(), '========='
						#a=raw_input("XXX")

						#if goterm.getGOID() == 'GO:0000221':
						#	print goterm.toString()

					goterm = GOTerm()

				else:
					if flag == 1: # �ϴ� ó�� ���ߵ��� ���� comment�� ������ �� ���� �͸� �Ű澴��.
						x = s.split(":")
						if x[0] == 'id': # GO_ID�� ���´�.
							y = s.split(' ')
							goterm.setGOID( y[1].strip() )
						elif x[0] == 'name':
							goterm.setGOName( x[1].strip() )
						elif x[0] == 'namespace':
							goterm.setGONameSpace( x[1] )
						elif x[0] == 'def':
							goterm.setGODefinition( x[1].strip().replace('"','') )
						elif x[0] == 'alt_id':
							y = s.split(' ')
							aid = y[1].strip() # alternative go id
							goterm.addAlternativeID(aid)

						elif x[0] == 'is_obsolete':
							if x[1].strip() == 'true':
								goterm.setObsolete(True)
						elif x[0] == 'is_a':
							# IS_A GO ID�� ���� �̴´�.
							y = s.split(' ')
							goterm.addParentRelationship( y[1], goterm.RELATION_IS_A )
						elif x[0] == 'relationship':
							y = s.split(' ')
							rel = y[1]
							gid = y[2]

							if rel == 'regulates' :
								goterm.addParentRelationship( gid, goterm.RELATION_REGULATES )
							elif rel == 'part_of':
								goterm.addParentRelationship( gid, goterm.RELATION_PART_OF )
							elif rel == 'negatively_regulates':
								goterm.addParentRelationship( gid, goterm.RELATION_NEGATIVELY_REGULATES)
							elif rel == 'positively_regulates':
								goterm.addParentRelationship( gid, goterm.RELATION_POSITIVELY_REGULATES)

							else:
								#print '�𸣴� RELATION ', rel
								pass





		f.close()

		# parent - child ���踦 ����� �����Ѵ�.
		self.assignRelation()

	def assignRelation(self):

		#print "���� ������"
		# parent-child ���谡 ����� parent�� �� �� �ִ�.
		# �̰� �����ؼ� child ���赵 �� �� �ְ� �����Ѵ�.
		
		# IS_A�� ���� �Ϸ��� �Ʒ� getParentRelationShipGOIDsOfISA()�� ����Ѵ�.
		# �׸��� �Ʒ� ������ ��� children/Ȥ�� ��� parent�� �ƴ϶� 1�ܰ��� ���踸 �����ϴ� ���̴�.
		# �� �� �����ؾ��Ѵ�.
		
		g = GOTerm()

		for g_id in self.go_terms:

			g = self.go_terms[g_id]

			v = g.getParentRelationshipGOIDs()
			#print repr(v)

			for p_gi in v:
				#pg = GOTerm()
				#print p_gi
				pg = self.getGOTermOf( p_gi )

				#print 'P: ', p_gi, '   .... ', g_id
				pg.addChildRelationship( g_id, g.getParentRelationshipOf(p_gi) )



	def printAll(self):
		for k in self.go_terms:
			g = GOTerm()
			g = self.go_terms[k]
			print (g.toString())

	'''  �м��� ���� �߰��� def '''
	def increaseCounter (self, go_id, cnt_id, flyDB_id):
		g = GOTerm()
		g = self.go_terms[go_id]

		g.increaseCount()


		parent_ids = g.getParentRelationshipGOIDs()

		for p in parent_ids:
			i = g.getParentRelationshipOf(p)
			if i == g.RELATION_IS_A: # or i == g.RELATION_PART_OF:
				self.increaseCounter( p , cnt_id )  # ���� GO term�� ī���� ������Ų��.

	def increasePoint (self, go_id, cnt_id, point, flyDBID):
		g = GOTerm()
		g = self.go_terms[go_id]

		g.increasePoint( cnt_id, point, flyDBID )


		parent_ids = g.getParentRelationshipGOIDs()

		# IS_A ���� ���� �͵� ������ ����.
		cnt = 0
		for p in parent_ids:
			i = g.getParentRelationshipOf(p)
			if i == g.RELATION_IS_A: # or i == g.RELATION_PART_OF:
				cnt = cnt + 1

		# ������ ������ �ø���.
		for p in parent_ids:
			i = g.getParentRelationshipOf(p)
			if i == g.RELATION_IS_A: # or i == g.RELATION_PART_OF:

				#self.increasePoint( p, cnt_id, point / float(cnt) )

				self.increasePoint( p, cnt_id, point, flyDBID  )

	def findRelation(self, go1_child, go2_parent):
		# go1�� go2 ������ ���踦 ã�´�.
		# ã�� ������ PARENT-CHILD ����, �׸��� distance�� �󸶳� �Ǵ���.
		# ����� ���� �������� �𸣹Ƿ� go1->go2, go2->go1 ��� ã�´�.
		storate = self.__findRelation(go1_child, go2_parent)
		# b = self.__findRelation(go2, go1) # �̰� ���� �ܰ迡�� �����Ѵ�.


		return storate


	def getDistance(self, go_child, go_parent):

		# CHild-Parent node���� �Ÿ��� �Ի��Ѵ�.
		# ���谡 ������, -1
		# �ڱ��ڽ��̸� 0
		# �ƴϸ� �Ÿ�
		# ���� �Լ����� ISA�� �����.

		box = self.findRelation(go_child, go_parent)
		if len(box) == 0:
			return -1 # parent-child���谡 �ƴϴ�.
		else:
			distance = 0
			for d, goid in box:
				if d>distance:
					distance = d
			return distance



	def __findRelation(self, go1, go2):

		# go2�� parent node��� ��������.
		storage = []

		if go1 == go2:
			# �Ȱ���.
			#print 'same'
			storage =  [ [ 0, go2 ] ]
			return storage # distance�� 0�̵ȴ�. �ٷ� �ڱ� �ڽ��̴ϱ�.

		self.__traceUp(0, go1, go2, storage)

		return storage  # [  [distance, go2 id] , [ ], [ ], .... ] �̷� ������ �ѱ��.

	def __traceUp(self, cnt, go1, go2, storage):

		# go2�� parent
		go1_term = self.getGOTermOf(go1)

		if go1_term is None:
			return


		#parent_ids = go1_term.getParentRelationshipGOIDs()
		# IS_A ������ parent�� ������ �´�.
		parent_ids = go1_term.getParentRelationshipGOIDsOfISA()

		cnt += 1

		for p in parent_ids:

			#print '�t��: ', cnt, p, go2

			if p == go2:
				# ���ϴ� ���� ã�Ҵٸ�,  node���� �Ÿ��� GO ������ �����Ѵ�.
				storage.append( [ cnt, go2 ] )
			else:
				# parent node ������ ���ϴ� GO_ID�� �ƴ϶�� �� �ܰ� ���� �ö󰣴�.

				self.__traceUp( cnt, p, go2, storage)






if __name__ == '__main__':
	print ('----')
	go = GO()
	#go.loadFile('./flybase/gene_ontology.obo')
	print (go.size())
	print ('===')
	#go.printAll()
