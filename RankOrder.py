# -*- coding: ms949 -*-

'''
Rank Order를 사용하는 통계 코드이다.

사용법,

1. loadRankRatio를 이용해서 rank order를 dict 형태로 받아온다.
2. 여러개의 rank order를 리스트 형태로 만들어서 getRankOrder 의 parameter로 넘겨준다.
3. (key, p-value) 형태로 최종 결과를 받는다.

'''

from Stat import statistics
import MyUtil_pypy


def __LoadFileAndSort(fname, reverse = False, ignore_first_line=False):

    '''
    파일을 읽어서 이름과 score를 얻어온다.
    단, score가 정렬이 안되어있을 수 있으니, 정렬한다.
    리턴은 (key, score) 의 튜플 형태로 넘긴다.
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
    파일을 읽어서 rank ratio를 dict 형태로 반환한다.
    '''

    rank = 1.0
    counter = 1.0

    prev_score = -10000000

    r = {}

    for name, score in __LoadFileAndSort( fname, reverse, ignore_first_line):


        # score가 동일한 것들은 같은 rank order를 가짐.
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
    파일을 읽어서 rank를 dict 형태로 반환한다.
    Spearman rank correlation 계산시 활용 가능
    '''

    rank = 1.0
    counter = 1.0

    prev_score = -10000000

    r = {}
    temp_names = []
    temp_ranks = []
    

    for name, score in __LoadFileAndSort( fname, reverse, ignore_first_line):


        # score가 동일한 것들은 rank의 평균을 가짐.
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
    # Rank order 로 넘겨받은 데이터에서 전체 이름을 뽑아온다.
    r = []

    for ro in ratio_dicts_as_list:
        for name in ro.keys():
            if not name in r:
                r.append(name)
    return r

def getRankOrder( ratio_dicts_as_list ):
    '''loadRankRatio 함수로 받아온 데이터를 리스트 형태로 담아주면, 이걸 계산해서
    (name, p-value) 형태로 데이터를 넘겨준다.
    '''

    ret = {}
    o = statistics.OrderStatistic()
    names = __getWholeNames( ratio_dicts_as_list )

    for name in names:
        q = []

        for ro in ratio_dicts_as_list:
            # 해당 rank ratio가 존재하면 q에 담는다.
            if name in ro:
                q.append( ro[name] )


        ret[name] = o.rankUsingOrderStatistics( q ) # name에 해당하는 p-value를 얻는다.


    # 정렬한 후 결과를 tuple (name, p-value) 형태로 넘겨준다.
    ret_sorted = sorted( ret.items(), key = lambda item: item[1], reverse = False )

    return ret_sorted



def writeRankOrder( rank_result, fname ):
    '''
    Rank 결과를 파일로 저장한다.
    '''

    f=open(fname,'w', encoding='utf-8')

    for name, pvalue in rank_result:
        s = name + '\t' + str(pvalue).strip()
        f.write(s+'\n')

    f.close()


if __name__ == '__main__':

    r1 = loadRankRatio('__test_rank1.txt', False) # 작은 score를 상위로 rank
    r2 = loadRankRatio('__test_rank2.txt', False)

    r_box = [ r1, r2 ]

    ret = getRankOrder( r_box )

    writeRankOrder( ret , '__test_order.txt')

    print 'done'


