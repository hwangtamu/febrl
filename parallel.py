# =============================================================================
# AUSTRALIAN NATIONAL UNIVERSITY OPEN SOURCE LICENSE (ANUOS LICENSE)
# VERSION 1.2
# 
# The contents of this file are subject to the ANUOS License Version 1.2
# (the "License"); you may not use this file except in compliance with
# the License. You may obtain a copy of the License at:
# 
#   http://datamining.anu.edu.au/linkage.html
# 
# Software distributed under the License is distributed on an "AS IS"
# basis, WITHOUT WARRANTY OF ANY KIND, either express or implied. See
# the License for the specific language governing rights and limitations
# under the License.
# 
# The Original Software is: "parallel.py"
# 
# The Initial Developers of the Original Software are:
#   Dr Ole Nielsen (Mathematical Sciences Institute, Australian National 
#                   University)
#   Dr Peter Christen (Department of Computer Science, Australian National
#                      University)
# 
# Copyright (C) 2002 - 2005 the Australian National University and
# others. All Rights Reserved.
# 
# Contributors:
# 
# Alternatively, the contents of this file may be used under the terms
# of the GNU General Public License Version 2 or later (the "GPL"), in
# which case the provisions of the GPL are applicable instead of those
# above. The GPL is available at the following URL: http://www.gnu.org/
# If you wish to allow use of your version of this file only under the
# terms of the GPL, and not to allow others to use your version of this
# file under the terms of the ANUOS License, indicate your decision by
# deleting the provisions above and replace them with the notice and
# other provisions required by the GPL. If you do not delete the
# provisions above, a recipient may use your version of this file under
# the terms of any one of the ANUOS License or the GPL.
# =============================================================================
#
# Freely extensible biomedical record linkage (Febrl) - Version 0.3
#
# See: http://datamining.anu.edu.au/linkage.html
#
# =============================================================================

"""Module parallel.py - Module for parallel imports and definitions.

   Written by Ole M. Nielsen, January 2003, ANU/MSI
   Extended by Peter Christen, February 2003

   Use PyPar for parallelism if it is installed, otherwise define a
   rudimentary interface for sequential execution.
"""

# =============================================================================

# Set the mode for printing ('host' or 'all')
#
printmode = 'host'  # 'all'

# Set mode for saving (writing) files (data sets) ('host' or 'all')
#
writemode = 'host'  # 'all'

# =============================================================================
# Imports go here

import logging

# -----------------------------------------------------------------------------
# Conditional import of PyPar
#
try:
  import pypar

# =============================================================================
# PyPar is not installed, so define sequential interface for parallel functions
#
except:
  logging.warn('Could not import module "PyPar", defining sequential '+ \
               'interface')
  def size(): return 1
  def rank(): return 0

  def Get_processor_name():  # - - - - - - - - - - - - - - - - - - - - - - - -
    import os
    try:
      hostname = os.environ['HOST']
    except:
      try:
        hostname = os.environ['HOSTNAME']
      except:
        hostname = 'Unknown'

    return hostname

  def Abort():  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    import sys
    sys.exit()

  def Finalize():  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    pass

  def Barrier():  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    pass

  def Wtime():  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    import time
    return time.time()

# -----------------------------------------------------------------------------
# PyPar is installed, so import it
#
else:
  from pypar import *

# -----------------------------------------------------------------------------
# Now define processor information for print statements (prompt)

if (size() > 1):
  prompt = 'P%d/%d: ' % (rank(),size())
else:
  prompt = ''  # Empty prompt for sequential runs

# =============================================================================
