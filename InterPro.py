# -*- coding: utf-8 -*-

'''
InterPro 내용을 읽어온다.
현재는 ID, description, GOmapping 정도만 가지고 있다.
'''


class InterPro():
	
	__fname_entry = 'entry.list'
	__fname_gomap = 'interpro2go.sdx'
	
	__interpro_desc = {}
	__go_mapping = {}
	
	
	#def __init__(self, path = './DBfiles/InterPro/'):
	def __init__(self, entry_file = './2. InterPro/entry.list', interpro2go_file = './2. InterPro/interpro2go.sdx'):
		
		self.__go_mapping = {}
		self.__interpro_desc = {}
		
		self.__fname_entry = entry_file
		self.__fname_gomap = interpro2go_file
		
		self.load()
		
		
	def load(self):
		''' 
		INTERPRO 파일 중 entry.list 파일과 interpro2go 파일이 있는 경로를 정해주면 된다.
		'''
		
		self.__go_mapping = {}
		self.__interpro_desc = {}

		self.__loadEntryList( self.__fname_entry)
		self.__loadGOMap( self.__fname_gomap)

		print ('Loaded Interpro IDs = ', len(self.__interpro_desc))
		print ('Loaded GO mapping data = ', len(self.__go_mapping))

	def getAllInterProIDs(self):
		return list(self.__interpro_desc)

	def getDescriptionOf(self, interpro_id):
		if interpro_id in self.__interpro_desc:
			return self.__interpro_desc[interpro_id]
		else:
			return interpro_id

	def getGOIDofInterProID(self, interpro_id):
		'''
		IPR 1개에 여러개의 GO가 들어있는 것들도 있다.
		그래서 리턴 값은 리스트로 나간다.
		'''
		return self.__go_mapping[interpro_id]
	
	def getInterProIDofGOID(self, go_id):
		
		for iid in self.__go_mapping:
			
			if go_id in self.__go_mapping[iid]:
				return iid
		return None

	def __loadEntryList(self, fname):
		
		f = open(fname)
		
		for s in f.readlines():
			
			s = s.replace('\n', '')
			
			if len(s) == 0:
				continue
			
			
			if s[:3] == 'IPR':
				iid = s[:9].strip()
				desc = s[10:].strip()
				
				self.__interpro_desc [ iid ] = desc
		
		f.close()
		
	def __loadGOMap(self,fname):

		f = open(fname, 'r', encoding='utf-8')
		
		for s in f.readlines():
			
			s = s.replace('\n', '')
			
			if len(s) == 0:
				continue
			
			
			if s[0] != '!':
				iid = s[9:18].strip()
				goid = s[ s.rfind(';') + 1: ].strip()

				if not iid in self.__go_mapping:
					self.__go_mapping[iid] = []
				self.__go_mapping[iid].append( goid )
				
		
		f.close()
		
if __name__ == '__main__':
	ipr = InterPro()
	print (ipr.getDescriptionOf('IPR000067'))
	#print ipr.getInterProIDofGOID('GO:0009431')