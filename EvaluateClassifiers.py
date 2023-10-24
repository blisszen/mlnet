# -*- coding: ms949 -*-

import math
import numpy

'''
Machine learning으로 만든 classifier의 효과를 측정한다.
이때 측정값을 뭘로 보여줄 것인지를 계산해주는 코드이다.


kappa statistics도 있긴한데, 이건 random 요소를 가미해야한다.
장점은 MCC와 같이 class 크기가 달라도 상관없다.
실제 계산하면 MCC와 값이 비슷하게 나온다.

kappa = (total_accuracy - random_accuracy) / ( 1- random_accuracy)
자세한 내용은 http://standardwisdom.com/softwarejournal/2011/12/confusion-matrix-another-single-value-metric-kappa-statistic/


'''


def get_F_score(TP, TN, FP, FN):
	F = -1
	try:
		precision = getPositivePredictiveValue(TP, TN, FP, FN)
		recall = getSensitivity(TP, TN, FP, FN)
		F = 2.0 * precision * recall / ( precision + recall )
	except:
		pass

	return F

def getMCC(TP, TN, FP, FN):
	''' Matthews correaltion coefficient '''

	base = math.sqrt ( (TP+FP) * (TP+FN) * (TN+FP) * (TN+FN) )
	if base == 0.0:
		base = 1.0

	cab = float (TP * TN - FP * FN)

	MCC = cab / base
	return MCC


def getSensitivity(TP, TN, FP, FN):
	try:
		s = float(TP) / float(TP+FN)
		return s
	except:
		return -1

def getSpecificity(TP, TN, FP, FN):
	try:
		s = float(TN) / float(FP+TN)
		return s
	except:
		return -1

def getPositivePredictiveValue(TP, TN, FP, FN):
	try:
		s = float(TP) /  float(TP+FP)
		return s
	except:
		return -1

def getNegativePredictiveValue(TP, TN, FP, FN):
	try:
		s= float(TN) / float(FN+TN)
		return s
	except:
		return -1


def calculateAUC_for_Ranks(predicted_rank_ratios, output = None):

	# rank ratio 값을 가지고 AUC 계산하는 코드임.
	# 주의: rank ratio는 positive data만 가지고 계산함.
	# rank ratio 값이 작을 수록 top rank 임

	if output is not None:
		f = open(output+'_pos.txt','w')
		f.write('False-positive rate(1-specificity)\tTrue-positive rate(sensitivity)\n')

	# rank ratio 범위는 0-1
	start_from = 0.0
	start_to = 1.0
	step = 0.01


	pos_false_positive_rate_x = [] # 1- specificity
	pos_true_positive_rate_y = [] # sensitivity


	threshold = start_from
	while(threshold <= start_to  ):
		threshold = threshold + step

		specificity = 1.0 - threshold

		total_trial = float( len( predicted_rank_ratios  ))

		success = 0.0
		for v in predicted_rank_ratios:
			if v <= threshold:
				success = success + 1

		sensitivity = success / total_trial

		pos_false_positive_rate_x.append( 1.0 - specificity  )
		pos_true_positive_rate_y.append( sensitivity )

		if output is not None:
			f.write(  str(1.0-specificity) + '\t' + str(sensitivity) + '\n'  )


	auc_pos = numpy.mean ( pos_true_positive_rate_y  )

	if output is not None:
		f.close()




	return auc_pos
