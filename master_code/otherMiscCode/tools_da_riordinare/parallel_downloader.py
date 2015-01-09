#!/usr/bin/env python

import sys, os
from subprocess import Popen, PIPE
from multiprocessing import Process, Queue

# Set the max number of processes
n_procs = 3

print "Defining workers..."

def downloader(input, output):
    while input.qsize() != 0:
        item = input.get()
        if item[0] == 60:
            step = 2.5
        else:
            step = 10
        cmd = 'wget  --user=bziosi --password=t4n4GJB  --cookies=on '\
            +'--keep-session-cookies --save-cookies=cookie.txt '\
            +'--load-cookies=cookie.txt -O ./downloads/'+str(item[0])\
            +'.csv -o ./downloads/'+str(item[0])+'.log '\
            +'"http://www.g-vo.org/MyMillennium3?action=doQuery&SQL=select '\
            +'x, y, z from mmsnapshots..millimilsnapshots where snapnum=63 '\
            'and x between '+`item[0]`+' and '+`item[0]+step`\
            +' order by x"'
        try:
#            p = Popen(cmd, shell=True, stdin=PIPE, stdout=PIPE, 
#                      stderr=PIPE, close_fds=True).wait()
            pid = os.getpid()
            pid_cmd = 'echo "'+str(item[0])+'" >> '+str(pid)+'.log'
            os.system(pid_cmd)
            os.system(cmd)
#            print "popen exit code ", p
        except:
            print "Popen not done, exit..."
            sys.exit()

def fill_queue(task_queue):
    for i in range(0, 70, 10):
        task_queue.put([i])
    return task_queue

def status(proc):
    if proc.is_alive==True:
        return 'alive'
    elif proc.is_alive==False:
        return 'dead'
    else:
        return proc.is_alive()

print "Define queues..."

input_queue = Queue()
output_queue = Queue()

try:
    input_queue = fill_queue(input_queue)
except:
    print "Queue not filled, exit..."
    sys.exit()

procs = []

try:
    for i in range(n_procs):
        procs.append(Process(target=downloader, args=(input_queue, output_queue)))
except:
    print "Creating processes not complete, exit..."
    sys.exit()

try:
    for i in procs:
        i.start()
except:
    print "Start processes not complete, exit..."
    sys.exit()

for i in procs:
    print "Process ", i," @ " , i.pid, " is ", status(i)

print "Done."
