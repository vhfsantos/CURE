from utilities import *
from functions_metrics import *
import os
from time import asctime
import sys


dataset_path = sys.argv[1]
try:
	output_path = sys.argv[2]
except:
	output_path = os.path.dirname(dataset_path)

print ("\n")
print ("analysing %s", dataset_path)

if not os.path.exists(output_path):
    os.makedirs(output_path)

os.chdir(output_path)

name = os.path.basename(dataset_path).rstrip(".nex")

print(asctime())
process_dataset_metrics(dataset_path, ['entropy'], minimum_window_size = 50, outfilename = '%s.csv' % (name))

print ("Done. Output is here %s", output_path)
