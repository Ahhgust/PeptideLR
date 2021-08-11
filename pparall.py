#!/usr/bin/python3
# Written by August Woerner

# MIT License

# Copyright (c) [2019] [August E. Woerner]

# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

from joblib import Parallel, delayed
import os


def inParallel(jobs, nProcs=2, prefer="threading", check=False, expectedReturn=0):
  """
  Trivial wrapper of joblib
  This takes an iterable of unix commands (jobs)
  and it runs each unix command in parallel (using os.system)
  Jobs are spread across nProcs
  and the exit codes are returned (check==False)
  if check is true the function returns True/False 
  if all processes exited normally (exit code == expectedReturn, True)
  """
  out = Parallel(n_jobs=nProcs, backend=prefer, verbose=0)( delayed( os.system)(i) for i in jobs)
  
  if check:
    return any(out) != expectedReturn
  
  return out
  


