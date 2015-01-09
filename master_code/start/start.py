#!/bin/env python

import sys, os
from subprocess import Popen, PIPE
from multiprocessing import Process, Queue

"""Start the serial.py jobs. Start this with (python start.py &> start.log)&
"""

n_procs = 10
file_1 = None
file_2 =  "mill2_fof_sorted.h5" #"sub_haloes_global_coords_sorted_indexed.h5" #"all_haloes_gif2.h5"
m_factor = 1    # How many random more than data
start_slice = 0
end_slice = 999    # Number of data slices

def starter(input, output):
  """Start the code processes one after the other"""
  while input.qsize() != 0:
    item = input.get()
    file_1 = "mill2sort-"+str(item[0])+"-extracted.h5"       #"gif2_"+str(item[0])+".h5"
    cmd = "nice -n 19 /opt/epd-7.0-2-rh5-x86_64/bin/python -u ../serial.py --file_1 "+file_1+" --file_2 "+file_2+\
            " -l 400 -n 0 --m_factor "+str(m_factor)+" --slicing --log ../logs/"+file_1+"-"+file_2

      try:
        pid = os.getpid()
        pid_cmd = 'echo "'+str(item[0])+'" >> '+str(pid)+'.log'
        os.system(pid_cmd)
        p = Popen(cmd, shell=True, stdin=PIPE, stdout=PIPE, 
                      stderr=PIPE, close_fds=True).wait()
#            os.system(cmd)
      except:
        print "Popen/os.system not done, exit..."
        sys.exit()

def fill_queue(task_queue, start_slice, end_slice):
  """Fill the queue"""
  for i in range(start_slice, end_slice+1, 1):
    task_queue.put([i])
  return task_queue

def status(proc):
  """Check for processes status"""
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
  input_queue = fill_queue(input_queue, start_slice, end_slice)
except:
  print "Queue not filled, exit..."
  sys.exit()

procs = []    # processes container

# Start the starter processes
try:
  for i in range(n_procs):
    print "Process loop ", i
    procs.append(Process(target=starter, args=(input_queue, output_queue)))
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

# Wait processes to finish
while len(procs) != 0:
  for i in procs:
    if status(i) == False:
      procs.pop(procs.index(i))
  time.sleep(10) # loose 10 seconds
    
print "Done."
