# -*- coding: utf-8 -*-

'''
����ڰ� �ʿ��� �Լ��� arg�� �Ѱ��ָ�
�˾Ƽ� ��Ƽ���μ����� �۾��� �ش�.


https://pymotw.com/2/multiprocessing/communication.html


����: Queue ũ�� ������ �־ ��� ũ�Ⱑ ū �͵��� Queue�� �ѱ��� ���Ѵ�.
�׷��� Queue�� �� thread (job)���� ���� �����, �̰� ���� �޾Ƽ� list�� ���� �� �Ѱ��ִ� ����� ����

multiprocess�� ������ �Լ��� ���� �� ������ multiprocess.Queue() �����ϰ�, �� ����� ���⿡ ���������.

user_arg�� �׻� ����Ʈ�� �޾ƾ� ��.

'''

import multiprocessing
import threading
import time
import sys




DEBUG = False


class JOB_THREAD_ONLY(threading.Thread):

	p = None
	over = 0
	args = None
	method = None
	result = None


	def __init__(self):

		threading.Thread.__init__(self)
		self.p = None
		self.over =0
		self.args = None
		self.method = None
		self.result = None


	def set_args(self, user_method, user_args):
		self.method = user_method
		self.args = user_args

	def is_over(self):
		return self.over

	def run(self):
		self.over = 1
		self.result = self.method ( self.args )
		self.over = 2

	def getResult(self):
		return self.result


# 이건 deadlock 피하기 위한 코드인데, 아직 테스트 못함.
class JOB_THREAD(threading.Thread):
    p = None
    over = 0
    args = None
    method = None
    result_queue = None
    quit_signal = False
    time_out = -1
    trial = 0
    ctx = None
    
    def __init__(self):
        
        threading.Thread.__init__(self)
        
        self.ctx = multiprocessing.get_context('spawn')
        
        self.p = None
        self.over = 0
        self.args = None
        self.method = None
        self.quit_signal = False
        self.result = None
        self.time_out = -1
        self.result_queue = None
    
    def set_timeout(self, minutes):
        self.time_out = minutes * 60  # 실제로는 초 단위로 계산함.
    
    def clone(self):
        t = JOB_THREAD()
        t.over = 0
        t.args = self.args
        t.method = self.method
        t.quit_signal = False
        t.time_out = self.time_out
        t.trial = self.trial + 1
        t.result = None
        
        return t
    
    def set_args(self, user_method, user_args):
        
        '''
        :param user_method: multiprocess�� ������ �Լ�
        :param user_args: �Լ����� �޾Ƶ��̴� args (list���� ��). ��, ���� ���������� queue ������ �߰��� �޾ƾ��ϰ�, ����� ���⿡ �־�� ��.
        :return:
        '''
        
        self.method = user_method
        if user_args is not None:
            self.result_queue = self.ctx.Queue()
            u_2 = user_args + [self.result_queue]
            # self.args = [user_args]
            self.args = [u_2]
    
    def getResult(self):
        return self.result
    
    def is_over(self):
        return self.over
    
    def run(self):
        self.over = 1
        
        # print self.args
        if self.args is not None:
            self.p = self.ctx.Process(target=self.method, args=self.args)
        else:
            self.p = self.ctx.Process(target=self.method)
        
        if DEBUG:
            self.p.run()
        else:
            
            self.p.start()
            
            counter = 0
            sleep_interval = 10
            
            while (self.p.is_alive()):
                time.sleep(sleep_interval)
                
                if self.quit_signal:
                    # print('Terminating an unfinished job')
                    self.p.terminate()
                    print('Terminated an unfinished job')
                    break
                
                counter += 10  # 지나간 시간 측정함.
                if counter >= self.time_out and self.time_out != -1:
                    self.p.terminate()
                    print('[MULTIPROCESS] Timeout! Process terminated !')
                    
                    self.over = 3  # failed process
                    return  # 그냥 종료해버림.
        
        self.over = 2
        
        # result 저장
        
        if self.result_queue is not None:
            if not self.result_queue.empty():
                self.result = self.result_queue.get()
            else:
                self.result = None
        else:
            self.result = None
    
    def terminate(self):
        # if self.over == 1: # running process
        #    self.p.terminate()
        self.quit_signal = True


'''
class JOB_THREAD_BAK(threading.Thread):

    p = None
    over = 0
    args = None
    method = None
    result_queue = multiprocessing.Queue()
    quit_signal = False
    time_out = -1
    trial = 0

    def __init__(self):

        threading.Thread.__init__(self)
        self.p = None
        self.over = 0
        self.args = None
        self.method = None
        self.quit_signal = False
        self.result = None
        self.time_out = -1

    def set_timeout(self, minutes):
        self.time_out = minutes * 60 # 실제로는 초 단위로 계산함.


    def clone(self):
        t = JOB_THREAD()
        t.over = 0
        t.args = self.args
        t.method = self.method
        t.quit_signal = False
        t.time_out = self.time_out
        t.trial = self.trial + 1
        t.result = None

        return t

    def set_args(self, user_method, user_args):

        self.method = user_method
        if user_args is not None:
            u_2 = user_args + [ self.result_queue ]
            #self.args = [user_args]
            self.args = [ u_2 ]

    def getResult(self):
        return self.result

    def is_over(self):
        return self.over



    def run(self):
        self.over = 1

        # print self.args
        if self.args is not None:
            self.p = multiprocessing.Process(target=self.method, args=self.args)
        else:
            self.p = multiprocessing.Process(target=self.method)

        if DEBUG:
            self.p.run()
        else:

            self.p.start()

            counter = 0
            sleep_interval = 10
            
            while(self.p.is_alive()):
                time.sleep(sleep_interval)

                if self.quit_signal:
                    #print('Terminating an unfinished job')
                    self.p.terminate()
                    print('Terminated an unfinished job')
                    break
            
                counter += 10 # 지나간 시간 측정함.
                if counter >= self.time_out and self.time_out != -1:
                    self.p.terminate()
                    print('[MULTIPROCESS] Timeout! Process terminated !')

                    self.over = 3 # failed process
                    return # 그냥 종료해버림.



        self.over = 2

        # result 저장

        if not self.result_queue.empty():
            self.result = self.result_queue.get()
        else:
            self.result = None

    def terminate(self):
        #if self.over == 1: # running process
        #    self.p.terminate()
        self.quit_signal = True
'''

'''
class JOB_THREAD_NO_RETURN_BAK(threading.Thread):
    p = None
    over = 0
    args = None
    method = None
    quit_signal = False
    time_out = -1
    trial = 0

    def __init__(self):

        threading.Thread.__init__(self)
        self.p = None
        self.over = 0
        self.args = None
        self.method = None
        self.quit_signal = False
        self.result = None
        self.time_out = -1

    def clone(self):
        t = JOB_THREAD_NO_RETURN()
        t.over = 0
        t.args = self.args
        t.method = self.method
        t.quit_signal = False
        t.time_out = self.time_out
        t.trial = self.trial + 1
        return t

    def set_args(self, user_method, user_args):
        self.method = user_method
        if user_args is not None:
            #u_2 = user_args + [self.result_queue]
            self.args = [user_args]
            #self.args = [u_2]

    def set_timeout(self, minutes):
        self.time_out = minutes * 60 # 실제로는 초 단위로 계산함.


    def is_over(self):
        return self.over

    def run(self):
        self.over = 1

        self.trial += 1 # time-out 후 재실할때 마다 1씩 증가한다.

        # print self.args
        if self.args is not None:
            self.p = multiprocessing.Process(target=self.method, args=self.args)
        else:
            self.p = multiprocessing.Process(target=self.method)

        if DEBUG:
            self.p.run()
        else:
            self.p.start()

            counter = 0
            sleep_interval = 10
            while (self.p.is_alive()):
                time.sleep(sleep_interval)

                if self.quit_signal:
                    # print('Terminating an unfinished job')
                    self.p.terminate()
                    print('Terminated an unfinished job')
                    break

                counter += 10 # 지나간 시간 측정함.
                if counter >= self.time_out and self.time_out != -1:
                    self.p.terminate()
                    print('[MULTIPROCESS] Timeout! Process terminated !')

                    self.over = 3 # failed process
                    return # 그냥 종료해버림.


            # self.p.join()

        self.over = 2


    def terminate(self):
        # if self.over == 1: # running process
        #    self.p.terminate()
        self.quit_signal = True
'''

# deadlock 피하기 위한 코드임.
class JOB_THREAD_NO_RETURN(threading.Thread):
    p = None
    over = 0
    args = None
    method = None
    quit_signal = False
    time_out = -1
    trial = 0

    def __init__(self):

        threading.Thread.__init__(self)
        self.p = None
        self.over = 0
        self.args = None
        self.method = None
        self.quit_signal = False
        self.result = None
        self.time_out = -1

    def clone(self):
        t = JOB_THREAD_NO_RETURN()
        t.over = 0
        t.args = self.args
        t.method = self.method
        t.quit_signal = False
        t.time_out = self.time_out
        t.trial = self.trial + 1
        return t

    def set_args(self, user_method, user_args):
        self.method = user_method
        if user_args is not None:
            #u_2 = user_args + [self.result_queue]
            self.args = [user_args]
            #self.args = [u_2]

    def set_timeout(self, minutes):
        self.time_out = minutes * 60 # 실제로는 초 단위로 계산함.


    def is_over(self):
        return self.over

    def run(self):
        self.over = 1

        self.trial += 1 # time-out 후 재실할때 마다 1씩 증가한다.
        
        ctx = multiprocessing.get_context('spawn')

        # print self.args
        if self.args is not None:
            self.p = ctx.Process(target=self.method, args=self.args)
        else:
            self.p = ctx.Process(target=self.method)

        if DEBUG:
            self.p.run()
        else:
            self.p.start()

            counter = 0
            sleep_interval = 10
            while (self.p.is_alive()):
                time.sleep(sleep_interval)

                if self.quit_signal:
                    # print('Terminating an unfinished job')
                    self.p.terminate()
                    print('Terminated an unfinished job')
                    break

                counter += 10 # 지나간 시간 측정함.
                if counter >= self.time_out and self.time_out != -1:
                    self.p.terminate()
                    print('[MULTIPROCESS] Timeout! Process terminated !')

                    self.over = 3 # failed process
                    return # 그냥 종료해버림.


            # self.p.join()

        self.over = 2


    def terminate(self):
        # if self.over == 1: # running process
        #    self.p.terminate()
        self.quit_signal = True


def runMultiprocesses(jobs_list, max_cpu=1, sleep_interval_seconds=1, title = '', cutoff=9999999):
    '''
    ���� JOB�̶�� Ŭ������ �̿��� �ʿ��� ��� �۾��� �߰��� ����,
    ���� Ŭ������ ����Ʈ�� ��Ƽ� �ѱ�� �ȴ�.
    ������ ��: ������ ��ȯ�� ����� �� �ǹǷ� ���� ���� ���� ���� ����.
    �׷��� ���� ���� �� Queue�� �޾ƾ� �Ѵ�.


    Job_thread Ŭ������ �����, �װ� multi-thread�� ���ư���.
    �׳� job class�� ����� multi-process�� ���ư���.


    runMultiprocesses() �Լ��� ����� return��. list��.

	cutoff: job_list�� n���� �������?, cutoff ��ŭ�� ����� �Ǹ� �������� ����� �ǵ縻�� �׳� �����Ѵ�.

    '''

    MAX_TRIAL = 3 # 오류가 발생하면 3번까지는 실행함.

    results = []

    while (True):

        cnt = 0
        done = 0

        # ���� ���ư��� �ִ� process ���� ����.
        for m in range(len(jobs_list)):
            if jobs_list[m].is_over() == 1:  # Running processes
                cnt += 1
            elif jobs_list[m].is_over() == 2:  # Completed processes
                done += 1
            elif jobs_list[m].is_over() == 3:  # time-out processes
                # 이건 재실행해야함.

                if jobs_list[m].trial >= MAX_TRIAL:
                    # 여러번 실행했으나 실패함.
                    # 그냥 포기하자.
                    print('[MULTIPROCESS] && Despite of max trials, failed!')
                    done += 1
                else:
                    print('[MULTIPROCESS] Re-run timed out process again!')
                    j = jobs_list[m].clone()
                    jobs_list[m] = j
                    jobs_list[m].start()
                    cnt += 1

        if done >= len(jobs_list) or done >= cutoff:
            # send quit signal to all threads...
            #print('Killing processes....')
            for j in jobs_list:
                j.terminate()

            break



        new_run = max_cpu - cnt  # ���� �����ؾ��� process �� �̴�.
        running = cnt
        for m in range(len(jobs_list)):
            if jobs_list[m].is_over() == 0 and new_run > 0:  # ���� ���� �� �ߴ�. �׷� �޷�����.
                jobs_list[m].start()
                #print('Start ', m,'th process', new_run)

                new_run -= 1  # 1�� �� �ش�.
                running += 1

                #if new_run == 0:
                #    break

            # print '\nrunning ', m, '/', new_run, time.ctime(),


        #print '\r' + title + ', DONE=', done, '/', len(jobs_list), 'Running=', running, '::', time.ctime(),
        print (title + ', DONE=', done, '/', len(jobs_list), ', Running=', running, ', Cutoff=', cutoff, '::', time.ctime() , end='\r')

        time.sleep(sleep_interval_seconds)


    #if len(results) != len(jobs_list):
    #    print '[Multiprocess] ERROR results # != job #'
    #    sys.exit(1)

    # 작업 순서대로 결과를 담아서 넘긴다.
    results = []
    for j in jobs_list:
        if j.is_over() == 2:  # �۾� ����

            rt = j.getResult()
            if rt is not None:
                results.append(rt)

            if len(results) >= len(jobs_list) or len(results) >= cutoff:
                break



    return results



def runMultiprocesses_no_return(jobs_list, max_cpu=1, sleep_interval_seconds=1, title = '', cutoff=9999999, display=True):


    # 이 함수는 JOB_THREAD_NO_RETURN 하고 같이 사용해야 함.
    # JOB_THREAD 와 함께 사용시 오류 발생 가능.
    MAX_TRIAL = 4 # time-out 시간을 증가시키면서 확인.

    while (True):

        cnt = 0
        done = 0

        for m in range(len(jobs_list)):
            if jobs_list[m].is_over() == 1:  # Running processes
                cnt += 1
            elif jobs_list[m].is_over() == 2:  # Completed processes
                done += 1
            elif jobs_list[m].is_over() == 3: # time-out processes
                # 이건 재실행해야함.

                if jobs_list[m].trial >= MAX_TRIAL:
                    # 여러번 실행했으나 실패함.
                    # 그냥 포기하자.
                    print('[MULTIPROCESS] && Despite of max trials, failed!')
                    done += 1
                else:
                    print('[MULTIPROCESS] Re-run timed out process again!')

                    j = jobs_list[m].clone()
                    jobs_list[m] = j
                    jobs_list[m].start()
                    cnt += 1


        if done >= len(jobs_list) or done >= cutoff:
            # send quit signal to all threads...
            #print('Killing processes....')
            for j in jobs_list:
                j.terminate()

            break



        new_run = max_cpu - cnt
        running = cnt
        for m in range(len(jobs_list)):
            if jobs_list[m].is_over() == 0 and new_run > 0:
                jobs_list[m].start()

                new_run -= 1  # 1�� �� �ش�.
                running += 1

                #if new_run == 0:
                #    break

            # print '\nrunning ', m, '/', new_run, time.ctime(),


        #print '\r' + title + ', DONE=', done, '/', len(jobs_list), 'Running=', running, '::', time.ctime(),
        if display:
            print ('[**]', title + ', DONE=', done, '/', len(jobs_list), ', Running=', running, ', Cutoff=', cutoff, '::', time.ctime() , end='\r')

        time.sleep(sleep_interval_seconds)


    # 리턴하지 않음.

def _test(args):
    i, j, q = args  # q = multiprocess.Queue()
    for c in range(i, j):
        time.sleep(1)
        print (c)
    q.put(j)


if __name__ == '__main__':

    '''
    job_box = []
    for i in range(10):
        j = JOB()

        j.set_args( _test, (i, i+10) )
        job_box.append(j)


    runMultiprocesses(job_box, max_cpu=2, sleep_interval_seconds=1)

    print 'Done'
    '''


    job_box = []
    for i in range(10):
        j = JOB_THREAD()

        j.set_args(_test, [i, i + 10])
        job_box.append(j)

    ret = runMultiprocesses(job_box, max_cpu=2, sleep_interval_seconds=1)

    print ('-=------------------------------')
    for r in ret:
        print (r)
