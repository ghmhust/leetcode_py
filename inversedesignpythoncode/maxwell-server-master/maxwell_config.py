""" Configuration file for Maxwell.

    Holds constants and such...
"""

import os 

# my_dir = tempfile.mkdtemp() # Temporary directory that we will use.
path = "/home/tmp/"


if not os.path.exists(path):
    os.makedirs(path)

def list_requests():
    return [f for f in os.listdir(path) \
                if f[-len('.request'):] == '.request']

