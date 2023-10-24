# -*- coding: utf-8 -*-

import MyUtil
import os
import sys
#import twill
#import sympy
from scipy.stats import hypergeom
from scipy.stats import norm
from scipy.stats import binom
from scipy.stats.mstats import gmean
from scipy.stats import pearsonr
import scipy.stats
import copy
import math
import numpy
from math import log, exp
#from mpmath import loggamma

#from numba import autojit, jit
#from numba import float32, int32




CONTENT='''
4919	48
617	0
659	2
173	3
326	0
1914	13
443	5
421	4
784	7
1261	2
1317	6

'''

#@jit(float32(float32[:]))
def average(a_list):
	a = numpy.array( a_list )
	#return float(sum(a_list))/float(len(a_list))
	return a.mean()
	

#@jit(float32(float32[:]))
def stdev(a_list):
	a = numpy.array(a_list)
	return a.std()

#@jit(float32(float32[:]))
def stderr(a_list):
	sd = stdev(a_list)
	c = math.sqrt(len(a_list))
	return sd / c
	
def harmonic_mean(a_list):
	return scipy.stats.mstats.hmean(a_list)

def geometric_mean(a_list):
	return scipy.stats.mstats.gmean(a_list)
'''
def geometric_mean(a_list):

	# geometric mean을 구한다.

	total = 1.0
	length = float( len(a_list) )

	for v in a_list:
		total *= v

	gm = total ** (1.0/length)

	return gm
'''

'''
def stdev_naive(data):


	mean = 0.0
	total = 0.0
	cnt = 0.0

	for d in data:
		total += d
		cnt += 1.0

	mean = (total/cnt)

	var = 0.0
	total = 0.0


	for d in data:
		total += (d - mean) ** 2

	var = total / cnt

	return math.sqrt(var)
'''




def stdev_quick(data):

	#http://en.wikipedia.org/wiki/Algorithms_for_calculating_variance
	# 위에 소개된 naive algorithm이다.
	# 실제 numpy의 stdev 함수와 비교했을 때 오차가 매우 적다.
	# online-algorithm 이 더 나을 것 같아 해봤더니, 그건 오차가 너무 크다 (아래의 stdev_quick)
	n = 0
	Sum = 0
	Sum_sqr = 0

	for x in data:
		n += 1
		Sum += x
		Sum_sqr += x*x

	#variance = (Sum_sqr - ((Sum*Sum)/n))/(n - 1)
	variance = (Sum_sqr - ((Sum*Sum)/n))/(n - 1 )
	return math.sqrt(variance)

'''
def stdev_quick(data):

	#online_variance():

	#http://en.wikipedia.org/wiki/Algorithms_for_calculating_variance
	#Variance를 빠르게 계산하는 방법이다.
	#stdev보다는 약간 오차가 발생하지만, mean 구하고 stdev 구하는 2번의 과정이 필요없다.


	n = 0
	mean = 0
	M2 = 0

	for x in data:
		n += 1
		delta = x - mean
		mean += delta/n
		M2 += delta*(x - mean)

	variance = M2/(n - 1)

	return math.sqrt(variance)

'''

'''
def calculateAll(content):
	x = content.split('\n')

	total_box = []
	init = True

	r = ''

	for s in x:


		s=s.strip()
		if len(s)>0:
			y = s.split('\t')

			if init==True:
				for ss in y:
					total_box.append( eval(ss) )
				init = False
			else:
				background_succ = eval( y[0] )
				background_num = total_box[0]

				for index in range( 1, len(y) ):
					sample_success = eval(y[index])
					sample_num = total_box[index]

					p = hypergeometric_distribution(sample_success, sample_num, background_succ, background_num)
					p2 = str(p)

					if len(r) == 0:
						r = p2
					else:
						r = r + '\t' + p2
				r = r + '\n'

	return r
'''

def sum(a_list):
	r = 0.0
	for v in a_list:
		r+=v
	return r



def get_T_test_p_value(list_a, list_b, axis = 0, equal_var = False):
	t_value, p_value = scipy.stats.ttest_ind( list_a, list_b, equal_var)
	return p_value

def __removeAlpha(a_list, b_list):
	a=[]
	b=[]

	for i in range(len(a_list)):
		v1 = a_list[i]
		v2 = b_list[i]

		if ((type(v1) is int) or (type(v1) is long) or ( type(v1) is float)) and ((type(v2) is int) or (type(v2) is long) or ( type(v2) is float)):
			a.append(v1)
			b.append(v2)

	return a,b

def getPearsonCorrelation_old(a_list, b_list):

	#https://en.wikipedia.org/wiki/Pearson_product-moment_correlation_coefficient

	'''
	Pearson-correlation 계산한다. List 2개 주면 숫자 아닌 것은 알아서 빼고 계산한다.
	'''

	if len( a_list ) != len(b_list):
		raise "Different length of two lists"

	#a2_list, b2_list = __removeAlpha(a_list, b_list) # 안에 숫자 아닌게 있으면 없애버린다.
	a2_list = a_list
	b2_list = b_list


	r, p = pearsonr(a2_list, b2_list)
	return r

	'''
	a_avg = average(a2_list)
	a_std = stdev(a2_list)
	b_avg = average(b2_list)
	b_std = stdev(b2_list)

	r = 0.0
	cnt = 0.0
	for index in range( len(a_list )):
		a = a2_list[index]
		b = b2_list[index]

		#if ((type(a) is int) or (type(a) is long) or (type(a) is float)) and ((type(b) is int) or (type(b) is long) or (type(b) is float)):
		cnt = cnt + 1.0
		s = ( float(a) - a_avg ) * ( float(b) - b_avg )
		r = r + s

	r = r / (cnt-1.0) / (a_std*b_std)
	

	return r
	'''

	
	
	
def getPearsonCorrelation(a_list, b_list):

	'''
	Pearson-correlation 계산한다. List 2개 주면 숫자 아닌 것은 알아서 빼고 계산한다.
	'''

	if len( a_list ) != len(b_list):
		raise "Different length of two lists"

	#a2_list, b2_list = __removeAlpha(a_list, b_list) # 안에 숫자 아닌게 있으면 없애버린다.
	a2_list = a_list
	b2_list = b_list


	r, p = pearsonr(a2_list, b2_list)
	return r

	'''
	a_avg = average(a2_list)
	a_std = stdev(a2_list)
	b_avg = average(b2_list)
	b_std = stdev(b2_list)

	r = 0.0
	cnt = 0.0
	
	v1 = v2 = v3 = v4 = 0.0
	
	for index in range( len(a_list )):
		a = a2_list[index]
		b = b2_list[index]

		#if ((type(a) is int) or (type(a) is long) or (type(a) is float)) and ((type(b) is int) or (type(b) is long) or (type(b) is float)):
		cnt = cnt + 1.0
		
		v1 += (a - a_avg)*(b - b_avg)
		v3 += (a - a_avg) ** 2
		v4 += (b - b_avg) ** 2
		
		
		#s = ( float(a) - a_avg ) * ( float(b) - b_avg )
		#r = r + s

	#r = r / (cnt-1.0) / (a_std*b_std)
	r = v1/ (v3*v4)**0.5

	return r	
	'''
	
	
def getBinomialDistribution(trials, success, pr):
	#r = binom(n, pr)
	r = binom.pmf(success, trials, pr)
	return r



#@jit(float32(float32, float32, float32))
def getNormalDistPvalue(mean, std, value):
	"""

	:rtype: object
	"""
	p_value = -1

	try:
		z = abs( float(( mean - value )/float(std) ) )
		#print z

		p_value = norm.sf(z) * 2
	except:
		p_value = -1



	#if p_value == 0.0: # might be range error
	#    p_value = normal_distributionWeb(mean, std, value)
	#    # eval을 씌우니 매우 작은 수는 0이 되어 버린다.

	return p_value

#@jit(float32(int32,int32,float32))
def getBinomialDistPvalue(x, n, pr):
	# n total trials, x success,with probability of pr
	prb = binom.cdf(x-1,n,pr)
	h = 1.0 - prb
	return h

#@jit(float32(int32,int32,int32,int32))
def getHypergeomPMF(success_in_sample, sample_number, success_in_pop, pop_number):
	pmf = hypergeom.pmf(success_in_sample, pop_number, success_in_pop, sample_number)
	return pmf



def getHypergeomCDF(success_in_sample, sample_number, success_in_pop, pop_number):
	cdf = hypergeom.cdf(success_in_sample, pop_number, success_in_pop, sample_number)

	return cdf

'''
def getHypergeomPvalue3(success_in_sample, sample_number, success_in_pop, pop_number):
	cdf = 0.0

	for x in range( success_in_sample + 1, sample_number + 1 ):
		cdf += getHypergeomPMF( x, sample_number, success_in_pop, pop_number)

	return cdf
'''

#@jit(float32(int32,int32,int32,int32))
def getHypergeomPvalue(success_in_sample, sample_number, success_in_pop, pop_number):


	if success_in_sample == 0:
		#print '!!! Success in sample = 0, returning -1'
		return -1

	p = hypergeom.pmf( numpy.arange(success_in_sample, sample_number+1), pop_number, success_in_pop , sample_number).sum()
	return p












'''
def logchoose(ni, ki):
	try:
		lgn1 = loggamma(ni+1)
		lgk1 = loggamma(ki+1)
		lgnk1 = loggamma(ni-ki+1)
	except ValueError:
		#print ni,ki
		raise ValueError
	return lgn1 - (lgnk1 + lgk1)




def gauss_hypergeom(X, n, m, N):
	"""Returns the probability of drawing X successes of m marked items
	 in n draws from a bin of N total items."""

	assert N >= m, 'Number of items %i must be larger than the number of marked items %i' % (N, m)
	assert m >= X, 'Number of marked items %i must be larger than the number of sucesses %i' % (m, X)
	assert n >= X, 'Number of draws %i must be larger than the number of sucesses %i' % (n, X)
	assert N >= n, 'Number of draws %i must be smaller than the total number of items %i' % (n, N)


	r1 = logchoose(m, X)
	try:
		r2 = logchoose(N-m, n-X)
	except ValueError:
		return 0
	r3 = logchoose(N,n)

	return exp(r1 + r2 - r3)

def hypergeom_cdf(success_in_sample, sample_number, success_in_pop, pop_number):
	#(X, n, m, N):

	N = pop_number
	m = success_in_pop
	X = success_in_sample
	n = sample_number


	assert N >= m, 'Number of items %i must be larger than the number of marked items %i' % (N, m)
	assert m >= X, 'Number of marked items %i must be larger than the number of sucesses %i' % (m, X)
	assert n >= X, 'Number of draws %i must be larger than the number of sucesses %i' % (n, X)
	assert N >= n, 'Number of draws %i must be smaller than the total number of items %i' % (n, N)
	assert N-m >= n-X, 'There are more failures %i than unmarked items %i' % (N-m, n-X)

	s = 0
	for i in range(0, X+1):
		s += max(gauss_hypergeom(i, n, m, N), 0.0)
	return min(max(s,0.0), 1)

def hypergeom_sf(success_in_sample, sample_number, success_in_pop, pop_number):

	N = pop_number
	m = success_in_pop
	X = success_in_sample
	n = sample_number

	assert N >= m, 'Number of items %i must be larger than the number of marked items %i' % (N, m)
	assert m >= X, 'Number of marked items %i must be larger than the number of sucesses %i' % (m, X)
	assert n >= X, 'Number of draws %i must be larger than the number of sucesses %i' % (n, X)
	assert N >= n, 'Number of draws %i must be smaller than the total number of items %i' % (n, N)
	assert N-m >= n-X, 'There are more failures %i than unmarked items %i' % (N-m, n-X)

	s = 0
	for i in range(X, min(m,n)+1):
		s += max(gauss_hypergeom(i, n, m, N), 0.0)
	return min(max(s,0.0), 1)



'''


'''
def getHypergeomPvalue2(success_in_sample, sample_number, success_in_pop, pop_number):
	# calculate p value with a good precision
	n_len = 50 # precision, 10^-30

	s_all = sympy.Symbol('s_all')
	s_all = sympy.N(1.0, n_len)

	s_tmp = sympy.Symbol('s_tmp')
	s_tmp = sympy.N( 0.0, n_len)

	s_pmf = sympy.Symbol('s_pmf')

	for i in range( success_in_sample   ):


		pmf = sympy.S (hypergeom.pmf(i, pop_number, success_in_pop, sample_number))

		s_pmf = sympy.N ( pmf , n_len)
		s_tmp = sympy.N( s_tmp + s_pmf , n_len)



	p = sympy.N(   s_all  - s_tmp , n_len )
	if p<1.0e-11:
		# inaccurate, raise an error
		print "Maybe inaccurate result: ", p
		# in this case, use the web.
		p = hypergeometric_distributionWeb(success_in_sample, sample_number, success_in_pop, pop_number)
	return p
'''


'''
def hypergeometric_distributionWeb(success_in_sample, sample_number, success_in_pop, pop_number):
	# calculate p_value using a web-tool
	# hypergeometric distribution-based calculation

	try:
		#site = 'http://keisan.casio.com/has10/SpecExec.cgi?path=07000000%2eProbability%20Function%2f01021000%2eHypergeometric%20distribution%2f12211000%2eHypergeometric%20distribution%2fdefault%2exml&charset=utf-8'
		site = 'http://keisan.casio.com/exec/system/1180573201'
		go(site)

		#showforms()

		fv("2", "var_x", str(success_in_sample) )
		fv("2", "var_no", str(sample_number) )
		fv("2", "var_M", str(success_in_pop) )
		fv("2", "var_N", str(pop_number) )


		submit()


		#showforms()
		#from twill.namespaces import get_twill_glocals
		#global_dict, local_dict = get_twill_glocals()



		temp_file = './'+MyUtil.getRandomString(20) + '.html'
		save_html(temp_file)


		# read from html

		key_q1 = "ans2"
		key_q3 = '>'
		key_q2 = "</div>"
		p_value = -1

		f=open('./a.html','r')
		for s in f.readlines():
			#print s
			index = s.find(key_q1)
			if index>0:
				# found the first hit
				#print s

				index = index + len(key_q1)
				index2 = s.find(key_q2, index)

				p_value = s[index:index2]

				index3 = p_value.find(key_q3)
				p_value = p_value[index3+1:]
				p_value = p_value.replace('<b>','').replace('</b>','')

				break
		f.close()

		try:
			os.remove(temp_file)
		except:
			pass

		print success_in_sample, sample_number, success_in_pop, pop_number, '->', p_value
		return eval(p_value)

	except:
		print 'Web site error'
		print '\tsuccess in sample=', success_in_sample
		print '\tsample number=', sample_number
		print '\tsuccess_in_pop=', success_in_pop
		print '\tpop_number=', pop_number
		return -1

'''



'''
def normal_distributionWeb(mean, std, value):
	# calculate p_value using a web-tool
	# hypergeometric distribution-based calculation


	site = 'http://keisan.casio.com/has10/SpecExec.cgi?id=system/2006/1180573188'
	go(site)

	#showforms()

	fv("2", "var_x", str(value) )
	fv("2", "var_μ", str(mean) )
	fv("2", "var_σ", str(std) )

	submit()


	#showforms()
	#from twill.namespaces import get_twill_glocals
	#global_dict, local_dict = get_twill_glocals()



	save_html('./a.html')


	# read from html

	key_q1 = "ans2"
	key_q3 = '>'
	key_q2 = "</div>"
	p_value = -1

	f=open('./a.html','r')
	for s in f.readlines():
		#print s
		index = s.find(key_q1)
		if index>0:
			# found the first hit
			#print s

			index = index + len(key_q1)
			index2 = s.find(key_q2, index)

			p_value = s[index:index2]

			index3 = p_value.find(key_q3)
			p_value = p_value[index3+1:]

			break
	f.close()



	#print success_in_sample, sample_number, success_in_pop, pop_number, '->', p_value
	return p_value
'''

# ---------------------------------------------------------------

class OrderStatistic:

	def __init__(self):
		pass

	def __factorial(self, n):
		r = 1.0
		for i in range(n):
			r=r*(i+1.0)
		return float(r)

	def rankUsingOrderStatistics(self, rank_ratios):
		'''
		the smaller the value is, the higher its rank is. !!!!!!!!!!!
		this function returns a p-value
		input:
			rank_ratios = [ 0.1 ,0.2, ...]
		return:
			p-value (smaller is better)
		'''

		N = len(rank_ratios)
		rank_ratios.sort()

		v =  [ 1.0 ]

		for k in range(1, N+1):
			vk = 0.0
			for i in range(1, k+1):
				vi = (-1.0)**(i-1.0) * v[ k-i ] / ( self.__factorial(i)  ) * rank_ratios[ N-k ] ** (i)
				#print '\t',vi
				vk=vk+vi
			v.append(vk)

		vn = v[-1] * self.__factorial(N)
		#p = 1.0 - vn
		p = vn #
		return p # the last element is Vn








def get_pvalue_of_fisher_exact_test(table, option = 'two-sided'):

	'''
	option=two-sided, less, greater
	table은 2x2 여야하고


		A   B
	C   1   5
	D   4   8

	위의 table이라면 table = [ [1,5], [4,8] ] 이렇게 넣어야 한다.
	'''

	oddsration, pvalue = scipy.stats.fisher_exact(table, option)
	return pvalue


def get_spearman_rank_correlation(things_in_order):

	'''

	:param things_in_order: list 이며, list 안에는 rank 순서대로 protein ID 등 식별할 수 있는 ID가 들어있어야 한다.
			list에는 기본적으로 2개가 들어가는데, 만약 3개 이상이라면 각각의 조합을 계산해서 평균 coefficient 값을 돌려준다.
	:return: Spearman's rank correlation coefficient
	'''


	box = []
	genes = copy.deepcopy(things_in_order[0])
	genes.sort()

	# 중복된 데이터만 사용한다.
	new_genes = []
	for g in genes:
		included = True
		for p_list in things_in_order:
			if not g in p_list:
				included = False

		if included:
			new_genes.append(g)

	genes = new_genes
	genes.sort()

	print ('Entered data = ', len(genes))

	for p_list in things_in_order:
		ranks = {}
		for i in range(len(p_list)):
			ranks[p_list[i]] = i + 1

		o_list = []
		for g in genes:
			o_list.append(ranks[g])
		box.append(o_list)

	ret = []
	for i in range(len(box) - 1):
		for j in range(i + 1, len(box)):
			p1 = box[i]
			p2 = box[j]

			r, p = scipy.stats.spearmanr(p1, p2)

			ret.append(r)

	return average(ret)



if __name__ == '__main__':
	#print getNormalDistPvalue(14, 3, 500)
	#print getBinomialDistPvalue(2,3,0.01)
	#print getBinomialDistribution(5,4,0.05)

	'''
	a = [ 1, 2, 3, 4, 5, 6, 7]
	b = [ 14, 15 ,16, 17 ]
	print get_T_test_p_value(a, b)
	print 'done'
	'''

	P = 100
	p = 150
	S = 20
	s = 5

	print ('-->',getHypergeomPvalue( 5, 20, 100, 150 ))
	#print hypergeometric_distributionWeb(s, S, p, P)
	#print  getHypergeomPvalue3( s, S, p, P )

	

	'''
	a = []
	import random
	for i in range(100):
		a.append( random.randint(0, 10))
	print average(a), stdev(a), stdev_naive(a)
	print stdev_quick(a)
	'''

	a = 124
	b = 96
	#table = [ [a,b], [14788-a,14809-b] ]
	table = [ [a,14788-a], [b,14809-b] ]
	print (get_pvalue_of_fisher_exact_test(table))

