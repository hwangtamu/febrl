# =============================================================================
# output.py - Classes for extracting output data sets for matches, non-matches,
#             and clerical reviews.
#
# Freely extensible biomedical record linkage (Febrl) Version 0.2
# See http://datamining.anu.edu.au/projects/linkage.html
#
# =============================================================================
# AUSTRALIAN NATIONAL UNIVERSITY OPEN SOURCE LICENSE (ANUOS LICENSE)
# VERSION 1.0
#
# The contents of this file are subject to the ANUOS License Version 1.0 (the
# "License"); you may not use this file except in compliance with the License.
# Software distributed under the License is distributed on an "AS IS" basis,
# WITHOUT WARRANTY OF ANY KIND, either express or implied. See the License for
# the specific language governing rights and limitations under the License.
# The Original Software is "output.py".
# The Initial Developers of the Original Software are Dr Peter Christen
# (Department of Computer Science, Australian National University), Dr Tim
# Churches (Centre for Epidemiology and Research, New South Wales Department
# of Health) and Drs Markus Hegland, Stephen Roberts and Ole Nielsen
# (Mathematical Sciences Insitute, Australian National University). Copyright
# (C) 2002 the Australian National University and others. All Rights Reserved.
# Contributors:
#
# =============================================================================

"""Module output.py - Classes for preparing output data sets.

   This module contains functions for printing of record pairs, histograms,
   saving of results files, compiling output data sets, etc.

   Still under development, the current version only allows printing of record
   pairs and histograms, and saving of simple text files containing record pair
   numbers and corresponding weights.
"""

# =============================================================================
# Imports go here

import os

# =============================================================================

def print_histogram(results_dict):
  """Print a histogram for the data stored in the results dictionary.

     The routine sums up the number of record pairs with matching weights for
     each integer value of matching weights.
  """

  histo = {}  # Start with an empty histogram dictionary

  min_histo_x_val = 9999  # Minimum X axis value in the histogram
  max_histo_x_val = -999  # Maximum X axis value in the histogram
  max_histo_y_val = -999  # Maximal Y axis value in the histogram

  # Loop over all record dictionaries - - - - - - - - - - - - - - - - - - - - -
  #
  for (rec_num, rec_dict) in results_dict.items():

    weights = rec_dict.values()  # Only the weight values are needed

    for w in weights:
      w_int = int(round(w))  # Round to closest integer

      w_count = histo.get(w_int,0)
      histo[w_int] = w_count+1  # Increase count by one

      if (w_count > max_histo_y_val):
        max_histo_y_val = w_count

  # Sort histogram according to X axis values
  #
  x_vals = histo.keys()
  y_vals = histo.values()

  histo_list = map(None, x_vals, y_vals)
  histo_list.sort()

  scale_factor_y = 70.0 / max_histo_y_val

  # Print the histogram using ASCII characters rotated 90 degrees - - - - - - -
  # so the X axis is going downwards
  #
  print '1:'
  print '1:Weight histogram:'
  print '1:-----------------'
  print '1:'

  for (x,y) in histo_list:
    print '1:'+str(x).rjust(5)+' *'+'*'*int(y*scale_factor_y)+'  %i' % (y)

  print '1:'

# =============================================================================

def print_record_pairs(dataset_a, dataset_b, results_dict, assigned_dict,
                       threshold):
  """Print the record pairs as stored in the given results dictionary.

     Takes as input references to two data sets (the ones that were used for
     the classifcation task), a results dictionary, a threshold and a
     dictionary that contains the assigned record pairs.

     Both the threshold and the assigned record pairs dictionary can be set to
     None.

     If a threshold is given, only record pairs with weights equal to or above
     this threshold are printed.

     If assigned record pairs are given then these record pairs are printed
     with a comment "[assigned]".
  """

  key_print_length = 12  # Maximal number of characters for keys
                         # (should be an even number)
  rec_print_length = (70 - key_print_length) / 2

  print '1:'
  print '1:Resulting record pairs:'
  print '1:-----------------------'
  print '1:  Threshold:  %f' % (threshold)
  print '1:  Data set A: %s' % (dataset_a.name)
  print '1:  Data set B: %s' % (dataset_b.name)

  if (threshold == None):  # Set threshold to a very negative number
    threshold == -9999999.9

  rec_numbers = results_dict.keys()
  rec_numbers.sort()

  for rec_num in rec_numbers:  # Loop over all records in results - - - - - - -

    rec_dict = results_dict[rec_num]
    rec_nums = rec_dict.keys()
    weights  = rec_dict.values()

    rec_list = map(None, weights, rec_nums)  # Sort according to weights
    rec_list.sort()
    rec_list.reverse()  # High weights first

    rec_a = dataset_a.read_record(rec_num)  # Get the first record

    rec_a_id = '[RecID A: %s/%s]' % \
               (str(rec_a['_rec_num_']), str(rec_a['_dataset_name_']))

    for (weight, rec) in rec_list:

      # Filter out record pairs with weight less than the threshold
      #
      if (weight >= threshold):

        if (assigned_dict != None) and (assigned_dict.has_key((rec_num, rec))):
          assigned = '[assigned]'
        else:
          assigned = ''

        rec_b = dataset_b.read_record(rec)  # get the second record

        # Build the union of all keys in both record dictionaries
        #
        tmp_dict = rec_a.copy()
        tmp_dict.update(rec_b)
        key_union = tmp_dict.keys()
        key_union.sort()

        print '1:'
        print '1:'+'-'*77
        print '1:Weight: %f %s' % (weight, assigned)

        rec_b_id = '[RecID B: %s/%s]' % \
                   (str(rec_b['_rec_num_']), str(rec_b['_dataset_name_']))

        print '1:Fields '+(key_print_length-7)*' '+' | '+ \
              rec_a_id.ljust(rec_print_length)[:rec_print_length]+' | '+ \
              rec_b_id.ljust(rec_print_length)[:rec_print_length]

        for key in key_union:
          key_str = str(key)
          if (key_str[0] != '_'):
            str_a = str(rec_a.get(key,''))
            str_b = str(rec_b.get(key,''))
            print '1:'+str(key).ljust(key_print_length)[:key_print_length] + \
                  ' | '+str_a.ljust(rec_print_length)[:rec_print_length] + \
                  ' | '+str_b.ljust(rec_print_length)[:rec_print_length]

  print '1:'

# =============================================================================

def save(dataset_a_name, dataset_b_name, file_name, results_dict,
         assigned_dict, threshold):
  """Save the results in a text file as record number pairs with corresponding
     weights.

     Takes as input the names of the two datasets, a file name, a results
     dictionary, a threshold and a dictionary that contains the assigned
     record pairs.

     Both the threshold and the assigned record pairs dictionary can be set to
     None.

     If a threshold is given, only record pairs with weights equal to or above
     this threshold are saved.

     If assigned record pairs are given then these record pairs are saved with
     a comment "[assigned]".
  """

  if (not isinstance(file_name, str)) or (file_name == ''):
    print 'error:Illegal file name: "%s"' % (str(file_name))
    raise Exception

  # Try to open the file in write mode
  #
  try:
    fp = open(file_name,'w')
  except:
    print 'error:Can not open file: "%s" for writing' % (str(file_name))
    raise IOError

  fp.write('Resulting record pairs:'+os.linesep)
  fp.write('-----------------------'+os.linesep)
  fp.write('  Threshold:  %f' % (threshold)+os.linesep)
  fp.write('  Data set A: %s' % (dataset_a_name)+os.linesep)
  fp.write('  Data set B: %s' % (dataset_b_name)+os.linesep)
  fp.write(os.linesep)

  if (threshold == None):  # Set threshold to a very negative number
    threshold == -9999999.9

  num_rec_pairs_saved = 0   # Count the number of record pairs saved

  rec_numbers = results_dict.keys()
  rec_numbers.sort()

  # For nice output, get the length (in digits) of the largest number
  #
  max_len = len(str(rec_numbers[-1]))

  for rec_num in rec_numbers:  # Loop over all records in results - - - - - - -

    # Save all print lines for this record into a list
    #
    rec_line_print = [str(rec_num)+os.linesep]

    rec_dict = results_dict[rec_num]
    rec_nums = rec_dict.keys()
    weights  = rec_dict.values()

    rec_list = map(None, weights, rec_nums)  # Sort according to weights
    rec_list.sort()
    rec_list.reverse()  # High weights first

    for (weight, rec) in rec_list:  # Filter out records

      # Filter out record pairs with weight less than the threshold
      #
      if (weight >= threshold):

        if (assigned_dict != None) and (assigned_dict.has_key((rec_num, rec))):
          assigned = '  [assigned]'
        else:
          assigned = ''

        rec_line_print.append('  '+str(rec).rjust(max_len)+':  '+ \
                                str(weight)+assigned+os.linesep)
        num_rec_pairs_saved += 1

    if (len(rec_line_print) > 1):  # There were record pairs for this record
      for l in rec_line_print:
        fp.write(l)

  fp.write(os.linesep)
  fp.write('  Number of record pairs saved: '+str(num_rec_pairs_saved)+ \
           os.linesep)
  fp.write(os.linesep)
  fp.write(os.linesep)
  fp.close()

# =============================================================================
