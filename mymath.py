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
# The Original Software is: "mymath.py"
# 
# The Initial Developer of the Original Software is:
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

"""Module mymath.py - Various mathematical routines.

   See doc strings of individual functions for detailed documentation.
"""

# =============================================================================
# Imports go here

import logging
import math

# =============================================================================

def mean(x):
  """Compute the mean (average)  of a list of numbers.
  """

  if (len(x) == 1):  # Only one element in list
    return float(x[0])

  elif (len(x) == 0):  # Empty list
    logging.info('Empty list given: %s' % (str(x)))
    return None

  else:  # Calculate average
    sum = 0.0
    for i in x:
      sum += i

    res = sum / float(len(x))

    return res

# =============================================================================

def stddev(x):
  """Compute the standard deviation of a list of numbers.
  """

  if (len(x) == 1):  # Only one element in list
    return 0.0

  elif (len(x) == 0):  # Empty list
    logging.info('Empty list given: %s' % (str(x)))
    return None

  else:
    sum = 0.0
    for i in x:
      sum += i

    avrg = sum / float(len(x))

    sum = 0.0
    for i in x:
      sum = sum + (i - avrg) * (i - avrg)

    res = math.sqrt(sum / float(len(x)))

    return res

# =============================================================================

def log2(x):
  """Compute binary logarithm (log2) for a floating-point number.

  USAGE:
    y = log2(x)

  ARGUMENT:
    x  An positive integer or floating-point number

  DESCRIPTION:
    This routine computes and returns the binary logarithm of a positive
    number.
  """

  return math.log(x) / 0.69314718055994529  # = math.log(2.0)

# =============================================================================

def perm_tag_sequence(in_tag_seq):
  """Create all permuations of a tag sequence.

  USAGE:
    seq_list = perm_tag_sequence(in_tag_seq)

  ARGUMENT:
    in_tag_seq  Input sequence (list) with tags

  DESCRIPTION:
    This routine computes all permutations of the given input sequence. More
    than one permutation is created if at least one element in the input
    sequence contains more than one tag.

    Returns a list containing tag sequences (lists).
  """

  if (not isinstance(in_tag_seq, list)):
    logging.exception('Input tag sequence is not a list: %s' % \
                      (str(in_tag_seq)))
    raise Exception

  list_len = len(in_tag_seq)
  out_tag_seq = [[]]  # List of output tag sequences, start with one empty list

  for elem in in_tag_seq:
    if ('/' in elem):  # Element contains more than one tag, covert into a list
      elem = elem.split('/')

    tmp_tag_seq = []

    if (isinstance(elem,str)):  # Append a simple string
      for t in out_tag_seq:
        tmp_tag_seq.append(t + [elem])  # Append string to all tag sequences

    else:  # Process a list (that contains more than one tags)
      for tag in elem:  # Add each tag in the list to the temporary tag list
        for t in out_tag_seq:
          tmp_tag_seq.append(t+[tag])  # Append string to all tag sequences

    out_tag_seq = tmp_tag_seq

  # A log message for high volume log output (level 3) - - - - - - - - - - - -
  #
  logging.debug('Input tag sequence: %s' % (str(in_tag_seq)))
  logging.debug('Output permutations:')
  for p in out_tag_seq:
    logging.debug('    %s' % (str(p)))

  return out_tag_seq

# =============================================================================
