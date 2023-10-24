
import FASTA

class ID_CLASS:
	
	
	alias_file = ''
	sequence_file = ''
	fasta = None
	aliases = {}
	common_name = {}  # key=string id, value=common name
	string_ids = []

	def __init__ (self, alias_file, sequence_file):
		self.alias_file = alias_file
		self.sequence_file = sequence_file
		
		self.aliases = {}
		self.common_name = {}
		self.string_ids = []
		
		#print ('Loading alias file: ', alias_file)
		self.aliases, self.common_name, self.string_ids = self.__getAliases(alias_file)
		#print('Loaded alias names = ', len(self.aliases))
		
		#print 'Loading sequence file: ', sequence_file
		if sequence_file is not None:
			self.fasta = FASTA.FASTA()
			self.fasta.load(sequence_file)
		
	def getSequenceOf(self, query_id):
		string_id = self.getSTRING_ID_of(query_id)
		if string_id is None:
			return None
		sequence = self.fasta.getSequenceOf(string_id)
		return sequence
		
	def getAllStringIDs(self):
		return self.string_ids
	
	def __getAliases(self, alias_file):
		box = {} # aliases
		pri_cbox = {} # common names
		sec_cbox = {} # secondary common names
		pri_dbox = {} # number of dbs using the alias
		sec_dbox = {}
		ter_cbox = {}
		ter_dbox = {}

		sbox = {} # string ids
		
		f=open(alias_file, 'r')
		
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
			
			sbox[string_id] = None
			
			box[alias_id]=string_id
			box[string_id]=string_id
			
			
			if alias_id in box:
				if box[alias_id] != string_id:
					print ('Inconsistent data ', alias_id, alias_db)
			
			
			# common names
			key = 'Ensembl_EntrezGene'
			key2 = 'BLAST_UniProt_GN'

			dbs = alias_db.split(' ')
			dbs_number = len(dbs)

			if key in dbs:

				if string_id in pri_cbox:
					# if there is already a common name,
					# replace it with one used in more databases

					if pri_dbox[string_id] < dbs_number:
						pri_cbox[string_id] = alias_id
						pri_dbox[string_id] = dbs_number
				else:
					pri_cbox[string_id] = alias_id
					pri_dbox[string_id] = dbs_number

			if key2 in dbs:

				if string_id in sec_cbox:

					if sec_dbox[string_id] < dbs_number:
						sec_cbox[string_id] = alias_id
						sec_dbox[string_id] = dbs_number

				else:
					sec_cbox[string_id] = alias_id
					sec_dbox[string_id] = dbs_number
			
			# tertiary
			if string_id in ter_cbox:
				if ter_dbox[string_id] < dbs_number:
					ter_cbox[string_id] = alias_id
					ter_dbox[string_id] = dbs_number
			else:
				ter_cbox[string_id] = alias_id
				ter_dbox[string_id] = dbs_number
	
		f.close()
		
		# if there is a string id without common name,
		# assign string id to it.

		for s_id in sbox:
			if not s_id in pri_cbox:
				if s_id in sec_cbox:
					pri_cbox[s_id] = sec_cbox[s_id]
				else:
					if s_id in ter_cbox:
						pri_cbox[s_id] = ter_cbox[s_id]
					else:
						pri_cbox[s_id] = s_id
		
		'''
		for s in cbox.keys():
			print s, cbox[s]
		'''
		
		
		return box, pri_cbox, list(sbox)
		
		
	def getSTRING_ID_of(self, your_id):
		if your_id in self.aliases:
			return self.aliases[your_id]
		else:
			return None		
	
	
	def getAliasesOf(self, your_id):
		# return all aliases of your_id including string ID
		box = []
		
		if not your_id in self.aliases:
			return []
		
		string_id = self.aliases[your_id]
		
		for k in self.aliases:
			
			if self.aliases[k] == string_id:
				if not k in box:
					box.append(k)
				
		return box
	
	
	def getCommonNameOf(self, your_id):
		
		
		
		if your_id in self.common_name:
			return self.common_name[your_id]
		else:
			string_id = self.getSTRING_ID_of(your_id)
			
			if string_id is not None:
				
				if string_id in self.common_name:
					return self.common_name[string_id]
				else:
					return None
			else:
			
				return None