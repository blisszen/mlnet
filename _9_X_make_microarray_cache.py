import os
import _9_1_Microarray
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
	
	_9_1_Microarray.MP_PROCESSORS = 10  # process 10개를 동시에 돌리면서 계산한다.
	_9_1_Microarray.TEMP_PATH = './TEMP'
	
	_0_Preprocess.STRING_ALIAS_FILE = './ppi/9606.protein.aliases.v11.0_2021.txt'
	_0_Preprocess.STRING_SEQUENCE_FILE = './ppi/9606.protein.sequences.v11.0_2021.fa'
	_0_Preprocess.ID_BOX = _0_Preprocess.ID.ID_CLASS(
		_0_Preprocess.STRING_ALIAS_FILE,
		_0_Preprocess.STRING_SEQUENCE_FILE
	)
	
	#data_folder = './9. Microarray/human_2022_all_tissues'
	data_folder = './9. Microarray//human_2022_microarray_brain_only'
	# data_folder = './9. Microarray/temp'
	
	array_files = []
	array_file_string = data_folder
	
	files, folders = MyUtil.getFileListIn(data_folder)
	for f in files:
		if f.find('norm.txt') < 0:
			# print('Skipped microarray file = ', f)
			continue
		
		array_files.append(f)
		array_file_string += ',' + f
	
	#md5 = _0_Preprocess.generateMD5(array_file_string)
	md5 = os.path.basename(data_folder)
	cache_file = './_cache_' + md5 + '.microarray.pickle'
	cache_db = './_cache_' + md5 + '.microarray.db'
	cache_genes = './_cache_' + md5 + '.microarray_genes.pickle'
	
	print('cache file = ', cache_file)
	print('cache db = ', cache_db)
	print('cache genes = ', cache_genes)
	
	cache_scores = None
	
	if not os.path.exists(cache_file):
		print('[Microarray] Cache not found: ' + cache_file)
		
		# multiprocess로 해도... pickle된 expression data를 읽어오느라 시간이 많이 걸림
		# 그래서 실제 연산 속도가 single process랑 큰 차이가 없음. ㅠㅠ
		# cache_scores = _9_1_Microarray.preprocess2(data_folder)  # list
		cache_scores = _9_1_Microarray.preprocess3(data_folder)  # list
		# [  [file, value]   ]    <- valie is dict
		
		print('[_9 Microarray.py] saving cache...', cache_file)
		# save md5
		
		# saveVariableAsPickle(cache_scores, filename = cache_file)
		
		# if os.path.exists(cache_file) == False:
		f = open(cache_file, 'wb')
		pickle.dump(cache_scores, f)
		f.close()
	
	if not (os.path.exists(cache_db) and cache_scores is None):
		
		if cache_scores is not None:
			# 새로 계산된 것이므로
			if os.path.exists(cache_db):
				print('Remove previous cache db: ', cache_db)
				os.remove(cache_db)
		
		# cache db 파일도 새로 만든다.
		if not os.path.exists(cache_db) and cache_scores is None:
			# pickle에서 가져온다.
			import MyUtil
			cache_scores = MyUtil.loadVariableFromPickle(cache_file)
		
		conn = sqlite3.connect(cache_db)
		cur = conn.cursor()
		
		table = 'CREATE TABLE expression_scores(id TEXT, score REAL)'
		cur.execute(table)
		
		i = 0
		step = 1000
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
			
			cur.execute('INSERT INTO expression_scores VALUES(?, ?)',
			            [key, cache_scores[key]])
		
		# indexing
		conn.commit()
		
		conn.execute('CREATE INDEX id_index ON expression_scores (id)')
		conn.commit()
		
		conn.close()
		
		# gene 리스트도 저장한다.
		box = {}
		for key in cache_scores:
			g1, g2 = key.split(':')
			box[g1] = 0
			box[g2] = 0
		
		total_genes = list(box)
		MyUtil.saveVariableAsPickle(total_genes, filename=cache_genes)
	
	else:
		print('Skipping saving cache_db: ', cache_db)
	# 이제 저장한다.
	
	print('done')



def make_for_fly():
	
	import MyUtil
	
	_9_1_Microarray.MP_PROCESSORS = 10  # process 10개를 동시에 돌리면서 계산한다.
	
	_0_Preprocess.STRING_ALIAS_FILE = './ppi/7227.protein.aliases.v11.0_2021.txt'
	_0_Preprocess.STRING_SEQUENCE_FILE = './ppi/7227.protein.sequences.v11.0_2021.fa'
	_0_Preprocess.ID_BOX = _0_Preprocess.ID.ID_CLASS(
		_0_Preprocess.STRING_ALIAS_FILE,
		_0_Preprocess.STRING_SEQUENCE_FILE
	)
	
	data_folder = './9. Microarray/fly_2021_all_tissues'
	# data_folder = './9. Microarray/temp'
	
	array_files = []
	array_file_string = data_folder
	
	files, folders = MyUtil.getFileListIn(data_folder)
	for f in files:
		if f.find('norm.txt') < 0:
			# print('Skipped microarray file = ', f)
			continue
		
		array_files.append(f)
		array_file_string += ',' + f
	
	#md5 = _0_Preprocess.generateMD5(array_file_string)
	md5 = os.path.basename(data_folder)
	
	cache_file = './_cache_' + md5 + '.microarray.pickle'
	cache_db = './_cache_' + md5 + '.microarray.db'
	cache_genes = './_cache_' + md5 + '.microarray_genes.pickle'
	
	print('cache file = ', cache_file)
	print('cache db = ', cache_db)
	print('cache genes = ', cache_genes)
	
	cache_scores = None
	
	if not os.path.exists(cache_file):
		print('[Microarray] Cache not found: ' + cache_file)
		
		# multiprocess로 해도... pickle된 expression data를 읽어오느라 시간이 많이 걸림
		# 그래서 실제 연산 속도가 single process랑 큰 차이가 없음. ㅠㅠ
		# cache_scores = _9_1_Microarray.preprocess2(data_folder)  # list
		cache_scores = _9_1_Microarray.preprocess3(data_folder)  # list
		# [  [file, value]   ]    <- valie is dict
		
		print('[_9 Microarray.py] saving cache...', cache_file)
		# save md5
		
		# saveVariableAsPickle(cache_scores, filename = cache_file)
		
		# if os.path.exists(cache_file) == False:
		f = open(cache_file, 'wb')
		pickle.dump(cache_scores, f)
		f.close()
	
	if not (os.path.exists(cache_db) and cache_scores is None):
		
		if cache_scores is not None:
			# 새로 계산된 것이므로
			if os.path.exists(cache_db):
				print('Remove previous cache db: ', cache_db)
				os.remove(cache_db)
		
		# cache db 파일도 새로 만든다.
		if not os.path.exists(cache_db) and cache_scores is None:
			# pickle에서 가져온다.
			import MyUtil
			cache_scores = MyUtil.loadVariableFromPickle(cache_file)
		
		conn = sqlite3.connect(cache_db)
		cur = conn.cursor()


		table = 'CREATE TABLE expression_scores(id TEXT, score REAL)'
		cur.execute(table)


		i = 0
		step = 1000
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
			
			cur.execute('INSERT INTO expression_scores VALUES(?, ?)',
			            [key, cache_scores[key]])
		
		# indexing
		conn.commit()
		
		conn.execute('CREATE INDEX id_index ON expression_scores (id)')
		conn.commit()
		
		conn.close()
		
		
		# gene 리스트도 저장한다.
		box = {}
		for key in cache_scores:
			g1, g2 = key.split(':')
			box[g1] = 0
			box[g2] = 0
		
		total_genes = list(box)
		MyUtil.saveVariableAsPickle(total_genes, filename = cache_genes)
	
	else:
		print('Skipping saving cache_db: ', cache_db)
	# 이제 저장한다.
	
	print('done')

def hhh():
	

	_0_Preprocess.STRING_ALIAS_FILE = './ppi/7227.protein.aliases.v11.0_2021.txt'
	_0_Preprocess.STRING_SEQUENCE_FILE = './ppi/7227.protein.sequences.v11.0_2021.fa'
	_0_Preprocess.ID_BOX = _0_Preprocess.ID.ID_CLASS(
		_0_Preprocess.STRING_ALIAS_FILE,
		_0_Preprocess.STRING_SEQUENCE_FILE
	)
	
	x = _0_Preprocess.ID_BOX.getSTRING_ID_of('CHIP')
	print(x)
	
if __name__ == '__main__':
	
	
	make_for_human()
	make_for_fly()