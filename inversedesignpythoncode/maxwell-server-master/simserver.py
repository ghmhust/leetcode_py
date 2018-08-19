""" Simulation server for Maxwell.

    Executes uploaded jobs.

    Consists of an infinite loop which does the following:
    1.  Find the oldest job.
    2.  Run the solver on it.
    3.  Repeat.
"""

import os, time
import maxwell_config
import subprocess, shlex

def find_oldest_job():
    req = maxwell_config.list_requests() # Get the requests.
    if not req:
        return None

    req_with_time = {}
    for r in req:
        req_with_time[r] = os.stat(maxwell_config.path + r).st_ctime

    oldest_req = min(req_with_time) # Run this job.
    os.remove(maxwell_config.path + oldest_req)
    return oldest_req.rstrip('.request')

if __name__ == '__main__':
    path_to_solver_dir = os.path.abspath(__file__).replace( \
                            __file__.split('/')[-1], 'maxwell-solver') + '/'
    while True:
        job = find_oldest_job()

        if not job:
            time.sleep(1)
            continue

        else:
            print "Solving %s..." % job
            # os.remove(maxwell_config.path + job)
            out_file = open(maxwell_config.path + job + '.log', 'w')
            return_code = subprocess.call(shlex.split( \
                                           "mpirun -n 3 python " + \
                                           path_to_solver_dir + "fdfd.py " + \
                                          maxwell_config.path + job), \
                                         stdout=out_file, stderr=subprocess.STDOUT)

            out_file.close()

            # Used to let user know that files can be downloaded.
            f = open(maxwell_config.path + job + '.finished', 'w')
            f.write('Finished')
            f.close()
               


