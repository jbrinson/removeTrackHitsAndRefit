#!/usr/bin/env python

# Purpose: Goes through log files from running DemoTrackAnalyzer and creates a script to remove root files associated with jobs that have failed

# Usage: python makeRemoveFailedJobsScript.py
#        source removeFailedAnalyzerJobs.src


from sys import argv
import os

#os.system('grep -r "Begin Fatal Exception" analyzerOutput_2/*.err > failedAnalyzerJobs.log')
os.system('grep -r "Begin Fatal Exception" analyzerOutput/*.err > failedAnalyzerJobs.log')

f = open('failedAnalyzerJobs.log')
f2 = open('removeFailedAnalyzerJobs.src', 'w')
for line in f:
    newline = 'rm analyzerOutput/hist' + line[21:-81] + '.root'
#    newline = 'rm analyzerOutput_2/hist' + line[23:-81] + '.root'
    f2.write(newline+'\n')
f.close()
f2.close()


## os.system('grep -r "Begin Fatal Exception" analyzerOutput_2/*.err > failedAnalyzerJobs_2.log')

## f3 = open('failedAnalyzerJobs_2.log')
## f4 = open('removeFailedAnalyzerJobs_2.src', 'w')
## for line in f3:
##     newline = 'rm analyzerOutput_2/hist' + line[21:-81] + '.root'
##     f4.write(newline+'\n')
## f3.close()
## f4.close()
    
    
