# -*- coding: ms949 -*-

'''
Rank Order�� ����ϴ� ��� �ڵ��̴�.

����,

1. loadRankRatio�� �̿��ؼ� rank order�� dict ���·� �޾ƿ´�.
2. �������� rank order�� ����Ʈ ���·� ���� getRankOrder �� parameter�� �Ѱ��ش�.
3. (key, p-value) ���·� ���� ����� �޴´�.

'''

from Stat import statistics
import MyUtil_pypy


def __LoadFileAndSort(fname, reverse = False, ignore_first_line=False):

    '''
    ������ �о �̸��� score�� ���´�.
    ��, score�� ������ �ȵǾ����� �� ������, �����Ѵ�.
    ������ (key, score) �� Ʃ�� ���·� �ѱ��.
    '''


    r = {}
    f = open(fname, 'r', encoding='utf-8')

    init = True

    for s in f.readlines():
        s = s.replace('\n', '')
        if len(s)>0:

            if init and ignore_first_line:
                init = False
                continue
            
            x = s.split('\t')

            name = x[0].strip()
            score = eval( x[1].strip() )

            r[name] = score

    f.close()


    rr = sorted( r.items(), key = lambda item: item[1], reverse = reverse)
    return rr



def loadRankRatio(fname, reverse, ignore_first_line=False):

    '''
    ������ �о rank ratio�� dict ���·� ��ȯ�Ѵ�.
    '''

    rank = 1.0
    counter = 1.0

    prev_score = -10000000

    r = {}

    for name, score in __LoadFileAndSort( fname, reverse, ignore_first_line):


        # score�� ������ �͵��� ���� rank order�� ����.
        if score != prev_score:
            prev_score = score
            rank = counter

            r[name]=rank
        else:
            r[name]=rank


        counter +=  1.0

    # normalize
    for g in r.keys():
        r[g]=r[g]/counter


    return r





def loadRank(fname, reverse, ignore_first_line=False):

    '''
    ������ �о rank�� dict ���·� ��ȯ�Ѵ�.
    Spearman rank correlation ���� Ȱ�� ����
    '''

    rank = 1.0
    counter = 1.0

    prev_score = -10000000

    r = {}
    temp_names = []
    temp_ranks = []
    

    for name, score in __LoadFileAndSort( fname, reverse, ignore_first_line):


        # score�� ������ �͵��� rank�� ����� ����.
        if score != prev_score:
            prev_score = score
            rank = counter

            r[name]=rank
            
            if len(temp_names) > 0:
                avg = MyUtil_pypy.mean(temp_ranks)
                for n in temp_names:
                    r[n] = avg
                    
                temp_names = []
                temp_ranks = []
                    
        else:
            
            temp_names.append(name)
            temp_ranks.append(counter)
            
            #r[name]=rank


        counter +=  1.0



    if len(temp_names) > 0:
        avg = MyUtil_pypy.mean(temp_ranks)
        for n in temp_names:
            r[n] = avg
            
        temp_names = []
        temp_ranks = []    


    return r


def __getWholeNames( ratio_dicts_as_list ):
    # Rank order �� �Ѱܹ��� �����Ϳ��� ��ü �̸��� �̾ƿ´�.
    r = []

    for ro in ratio_dicts_as_list:
        for name in ro.keys():
            if not name in r:
                r.append(name)
    return r

def getRankOrder( ratio_dicts_as_list ):
    '''loadRankRatio �Լ��� �޾ƿ� �����͸� ����Ʈ ���·� ����ָ�, �̰� ����ؼ�
    (name, p-value) ���·� �����͸� �Ѱ��ش�.
    '''

    ret = {}
    o = statistics.OrderStatistic()
    names = __getWholeNames( ratio_dicts_as_list )

    for name in names:
        q = []

        for ro in ratio_dicts_as_list:
            # �ش� rank ratio�� �����ϸ� q�� ��´�.
            if name in ro:
                q.append( ro[name] )


        ret[name] = o.rankUsingOrderStatistics( q ) # name�� �ش��ϴ� p-value�� ��´�.


    # ������ �� ����� tuple (name, p-value) ���·� �Ѱ��ش�.
    ret_sorted = sorted( ret.items(), key = lambda item: item[1], reverse = False )

    return ret_sorted



def writeRankOrder( rank_result, fname ):
    '''
    Rank ����� ���Ϸ� �����Ѵ�.
    '''

    f=open(fname,'w', encoding='utf-8')

    for name, pvalue in rank_result:
        s = name + '\t' + str(pvalue).strip()
        f.write(s+'\n')

    f.close()


if __name__ == '__main__':

    r1 = loadRankRatio('__test_rank1.txt', False) # ���� score�� ������ rank
    r2 = loadRankRatio('__test_rank2.txt', False)

    r_box = [ r1, r2 ]

    ret = getRankOrder( r_box )

    writeRankOrder( ret , '__test_order.txt')

    print 'done'


