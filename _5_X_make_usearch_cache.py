import os
import _5_Usearch2
import _0_Preprocess
import pickle
import sqlite3


'''
Microarray의 계산값을 미리 만들어 놓는 코드임.
아래에 STRING 파일과 microarray 데이터가 들어 있는 폴더를 기입하면 됨.
Multiprocess를 사용함.

_9_1_Microarray.py 내부에서 자동으로 생성하기는 하지만,
single process로 돌아가므로 여기서 미리 만드는게 좋음.
참고로, runFinal() 함수 내부 (Ln 607-612)에 preprocess2 (single process) 혹은 preprocess3 (multi)호출하는 부분이 있음.
'''


def make_for_human():
	
	import MyUtil
	
	_5_Usearch2.TEMP_PATH = './TEMP'
	_5_Usearch2.MP_PROCESSORS = 10 # 10 processes
	
	cache_scores = None
	
	ppi_file = os.path.basename ('./ppi/9606.protein.links.v11.0_2021.txt')
	
	_0_Preprocess.STRING_ALIAS_FILE = './ppi/9606.protein.aliases.v11.0_2021.txt'
	_0_Preprocess.STRING_SEQUENCE_FILE = './ppi/9606.protein.sequences.v11.0_2021.fa'
	_0_Preprocess.ID_BOX = _0_Preprocess.ID.ID_CLASS(
		_0_Preprocess.STRING_ALIAS_FILE,
		_0_Preprocess.STRING_SEQUENCE_FILE
	)
	
	cache_file = './_cache.usearch_' + ppi_file + '.pickle'
	cache_db = './_cache.usearch_' + ppi_file + '.db'
	total_genes = _0_Preprocess.ID_BOX.getAllStringIDs()
	
	print('cache file = ', cache_file)
	print('cache db = ', cache_db)
	print('# of genes = ', len(total_genes))
	
	if not os.path.exists(cache_file):
		print('[USearch] Cache not found: ' + cache_file)
		
		# multiprocess로 해도... pickle된 expression data를 읽어오느라 시간이 많이 걸림
		# 그래서 실제 연산 속도가 single process랑 큰 차이가 없음. ㅠㅠ
		# cache_scores = _9_1_Microarray.preprocess2(data_folder)  # list
		cache_scores = _5_Usearch2.cache_processing_mpi4(total_genes)
		# [  [file, value]   ]    <- valie is dict
		
		print('[_5_Usearch.py] saving cache...', cache_file)
		# save md5
		
		# saveVariableAsPickle(cache_scores, filename = cache_file)
		
		# if os.path.exists(cache_file) == False:
		f = open(cache_file, 'wb')
		pickle.dump(cache_scores, f)
		f.close()
	
	else:
		cache_scores = MyUtil.loadVariableFromPickle(cache_file)
		print('Loaded cache file = ', cache_file)
	
	if os.path.exists(cache_db):
		print('Remove previous cache db: ', cache_db)
		os.remove(cache_db)
	
	conn = sqlite3.connect(cache_db)
	cur = conn.cursor()
	
	table = 'CREATE TABLE usearch_scores(id TEXT, score REAL)'
	cur.execute(table)
	
	i = 0
	step = 100
	cnt = 0
	total = len(cache_scores)
	
	for key in cache_scores:
		
		# display progress
		
		i += 1
		cnt += 1
		if i == step:
			i = 0
			per = round(cnt / total * 100, 2)
			print('Progress : ', per, '%', end='\r')
		
		cur.execute('INSERT INTO usearch_scores VALUES(?, ?)',
		            [key, cache_scores[key]])
	
	# indexing
	conn.commit()
	
	conn.execute('CREATE INDEX id_index ON usearch_scores (id)')
	conn.commit()
	
	conn.close()
	
	print('done')



def make_for_fly():
	
	import MyUtil
	
	_5_Usearch2.TEMP_PATH = './TEMP'
	_5_Usearch2.MP_PROCESSORS = 10  # 10 processes
	
	cache_scores = None
	
	ppi_file = os.path.basename ('./ppi/7227.protein.links.v11.0_2021.txt')
	
	_0_Preprocess.STRING_ALIAS_FILE = './ppi/7227.protein.aliases.v11.0_2021.txt'
	_0_Preprocess.STRING_SEQUENCE_FILE = './ppi/7227.protein.sequences.v11.0_2021.fa'
	_0_Preprocess.ID_BOX = _0_Preprocess.ID.ID_CLASS(
		_0_Preprocess.STRING_ALIAS_FILE,
		_0_Preprocess.STRING_SEQUENCE_FILE
	)
	
	cache_file = './_cache.usearch_' + ppi_file + '.pickle'
	cache_db = './_cache.usearch_' + ppi_file + '.db'
	total_genes = _0_Preprocess.ID_BOX.getAllStringIDs()
	
	print('cache file = ', cache_file)
	print('cache db = ', cache_db)
	print('# of genes = ', len(total_genes))
	
	if not os.path.exists(cache_file):
		print('[USearch] Cache not found: ' + cache_file)
		
		# multiprocess로 해도... pickle된 expression data를 읽어오느라 시간이 많이 걸림
		# 그래서 실제 연산 속도가 single process랑 큰 차이가 없음. ㅠㅠ
		# cache_scores = _9_1_Microarray.preprocess2(data_folder)  # list
		cache_scores = _5_Usearch2.cache_processing_mpi4(total_genes)
		# [  [file, value]   ]    <- valie is dict
		
		print('[_5_Usearch.py] saving cache...', cache_file)
		# save md5
		
		# saveVariableAsPickle(cache_scores, filename = cache_file)
		
		# if os.path.exists(cache_file) == False:
		f = open(cache_file, 'wb')
		pickle.dump(cache_scores, f)
		f.close()
	
	else:
		cache_scores = MyUtil.loadVariableFromPickle(cache_file)
		print('Loaded cache file = ', cache_file)
	
	
	if os.path.exists(cache_db):
		print('Remove previous cache db: ', cache_db)
		os.remove(cache_db)
		
	
	
	conn = sqlite3.connect(cache_db)
	cur = conn.cursor()
	
	table = 'CREATE TABLE usearch_scores(id TEXT, score REAL)'
	cur.execute(table)
	
	i = 0
	step = 100
	cnt = 0
	total = len(cache_scores)
	
	for key in cache_scores:
		
		# display progress
		
		i += 1
		cnt += 1
		if i == step:
			i = 0
			per = round(cnt / total * 100, 2)
			print('Progress : ', per, '%', end='\r')
		
		cur.execute('INSERT INTO usearch_scores VALUES(?, ?)',
		            [key, cache_scores[key]])
	
	# indexing
	conn.commit()
	
	conn.execute('CREATE INDEX id_index ON usearch_scores (id)')
	conn.commit()
	
	conn.close()
	
	print('done')

if __name__ == '__main__':
	make_for_fly()
	make_for_human()
	