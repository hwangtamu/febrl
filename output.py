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
# The Original Software is: "output.py"
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

"""Module output.py - Classes for extracting output data sets for matches,
                      non-matches, and clerical reviews.

   This module contains functions for printing and savinf to text files of
   record pairs (detailed and weights), and histograms,  compiling output data
   sets, etc.

   Still under development, the current version allows printing and saving of
   record pairs and histograms, and saving of simple text files containing
   record pair numbers and corresponding weights, as well as comma separated
   file (CSV) output.
"""

# =============================================================================
# Imports go here

import logging
import os

# =============================================================================

def histogram(results_dict, file_name=None):
  """Print or save a histogram for the data stored in the results dictionary.

     The routine sums up the number of record pairs with matching weights for
     each integer value of matching weights.

     If a file name is given, the output will be written into this text file.
  """

  histo = {}  # Start with an empty histogram dictionary

  min_histo_x_val = 9999  # Minimum X axis value in the histogram
  max_histo_x_val = -999  # Maximum X axis value in the histogram
  max_histo_y_val =    1  # Maximal Y axis value in the histogram
                          # Bug-fix by Marion Sturtevant (thanks!)

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

  # If a file name is given open it for writing - - - - - - - - - - - - - - - -
  #
  if (file_name != None):
    try:
      f = open(file_name, 'w')
    except:
      logging.exception('Can not open file "%s" for writing' % \
                        (str(file_name)))
      raise IOError

    f.write('Weight histogram:' + os.linesep)
    f.write('-----------------' + os.linesep)
    f.write(os.linesep)

  else:  # Print a header for the histogram - - - - - - - - - - - - - - - - - -
    print
    print 'Weight histogram:'
    print '-----------------'
    print

  # Print or save the histogram (using ASCII characters) rotated 90 degrees - -
  # so the X axis is going downwards
  #
  for (x,y) in histo_list:
    hist_str = str(x).rjust(5)+' *'+'*'*int(y*scale_factor_y)+'  %i' % (y)

    if (file_name != None):
      f.write(hist_str + os.linesep)
    else:
      print hist_str

  if (file_name != None):
    f.close()
    print 'Histogram written to file "%s"' % (file_name)

  else:
    print

# =============================================================================

def rec_pair_details(dataset_a, dataset_b, results_dict, assigned_dict,
                     threshold, file_name=None):
  """Print or save the record pairs in details (all stored fields) as stored in
     the given results dictionary.

     Takes as input references to two data sets (the ones that were used for
     the classifcation task), a results dictionary, a threshold and a
     dictionary that contains the assigned record pairs.

     Both the threshold and the assigned record pairs dictionary can be set to
     None.

     If a threshold is given, only record pairs with weights equal to or above
     this threshold are printed.

     If assigned record pairs are given then these record pairs are printed
     with a comment "[assigned]".

     If a file name is given, the output will be written into this text file.
  """

  key_print_length = 12  # Maximal number of characters for keys
                         # (should be an even number)
  rec_print_length = (70 - key_print_length) / 2

  # If a file name is given open it for writing - - - - - - - - - - - - - - - -
  #
  if (file_name != None):
    try:
      f = open(file_name, 'w')
    except:
      logging.exception('Can not open file "%s" for writing' % \
                        (str(file_name)))
      raise IOError

    f.write('Resulting record pairs:' + os.linesep)
    f.write('-----------------------' + os.linesep)
    f.write('  Output threshold:  %f' % (threshold) + os.linesep)
    f.write('  Data set A:        %s' % (dataset_a.name) + os.linesep)
    f.write('  Data set B:        %s' % (dataset_b.name) + os.linesep)

  else:  # Print a header for the histogram - - - - - - - - - - - - - - - - - -
    print
    print 'Resulting record pairs:'
    print '-----------------------'
    print '  Output threshold:  %f' % (threshold)
    print '  Data set A:        %s' % (dataset_a.name)
    print '  Data set B:        %s' % (dataset_b.name)

  if (threshold == None):  # Set threshold to a very negative number
    threshold == -9999999.9

  # Create a dictionary with record pairs and sort them according to weight - -
  #
  rec_pair_dict = {}

  for rec_num in results_dict:
    rec_dict = results_dict[rec_num]

    for rec_num2 in rec_dict:
      weight = rec_dict[rec_num2]
      rec_pair = (rec_num, rec_num2)

      if (rec_pair_dict.has_key(rec_pair)):
        logging.warn('Record pair %s more than once in results dictionary' % \
              (str(rec_pair)))
      else:
        rec_pair_dict[rec_pair] = weight

  # Now sort record pairs according to their weight
  #
  rec_pair_numbers = rec_pair_dict.keys()
  rec_pair_weights = rec_pair_dict.values()

  rec_pair_list = map(None, rec_pair_weights, rec_pair_numbers)
  rec_pair_list.sort()
  rec_pair_list.reverse()  # Large weights first

  # Loop over all record pairs  - - - - - - - - - - - - - - - - - - - - - - - -
  #
  for (rec_pair_weight, rec_pair_numbers) in rec_pair_list:

    # Filter out record pairs with weight less than the threshold
    #
    if (rec_pair_weight >= threshold):

      rec_num_a = rec_pair_numbers[0]
      rec_num_b = rec_pair_numbers[1]

      rec_a = dataset_a.read_record(rec_num_a)  # Get the first record
      rec_b = dataset_b.read_record(rec_num_b)  # Get the second record

      rec_a_id = '[RecID A: %s/%s]' % \
                 (str(rec_a['_rec_num_']), str(rec_a['_dataset_name_']))
      rec_b_id = '[RecID B: %s/%s]' % \
                 (str(rec_b['_rec_num_']), str(rec_b['_dataset_name_']))

      # Check if this record pair is assigned
      #
      if (assigned_dict != None) and \
         (assigned_dict.has_key(rec_pair_numbers)):
        assigned = '[assigned]'
      else:
        assigned = ''

      # Build the union of all keys in both record dictionaries
      #
      tmp_dict = rec_a.copy()
      tmp_dict.update(rec_b)
      key_union = tmp_dict.keys()
      key_union.sort()

      if (file_name != None):
        f.write(os.linesep)
        f.write('-'*77 + os.linesep)
        f.write('Weight: %f %s' % (rec_pair_weight, assigned) + os.linesep)
        f.write('Fields '+(key_print_length-7)*' '+' | '+ \
                rec_a_id.ljust(rec_print_length)[:rec_print_length]+' | '+ \
                rec_b_id.ljust(rec_print_length)[:rec_print_length] + \
                os.linesep)
      else:
        print
        print '-'*77
        print 'Weight: %f %s' % (rec_pair_weight, assigned)

        print 'Fields '+(key_print_length-7)*' '+' | '+ \
              rec_a_id.ljust(rec_print_length)[:rec_print_length]+' | '+ \
              rec_b_id.ljust(rec_print_length)[:rec_print_length]

      for key in key_union:
        key_str = str(key)
        if (key_str[0] != '_'):
          str_a = str(rec_a.get(key,''))
          str_b = str(rec_b.get(key,''))

          if (file_name != None):
            f.write(str(key).ljust(key_print_length)[:key_print_length] + \
                  ' | '+str_a.ljust(rec_print_length)[:rec_print_length] + \
                  ' | '+str_b.ljust(rec_print_length)[:rec_print_length] + \
                  os.linesep)
          else:
            print str(key).ljust(key_print_length)[:key_print_length] + \
                  ' | '+str_a.ljust(rec_print_length)[:rec_print_length] + \
                  ' | '+str_b.ljust(rec_print_length)[:rec_print_length]

  if (file_name != None):
    f.close()
    print 'Detailed record pairs written to file "%s"' % (file_name)

  else:
    print

# =============================================================================

def rec_pair_weights(dataset_a_name, dataset_b_name, results_dict,
         assigned_dict, threshold, file_name=None):
  """Print or save the results as record number pairs with corresponding
     weights.

     Takes as input the names of the two datasets, a results dictionary, a
     dictionary that contains the assigned record pairs, a threshold and a
     file name (set to None as default).

     Both the threshold and the assigned record pairs dictionary can be set to
     None.

     If a threshold is given, only record pairs with weights equal to or above
     this threshold are saved.

     If assigned record pairs are given then these record pairs are saved with
     a comment "[assigned]" (or a "assigned" in a CSV file).

     If a file name is given, the output will be written into this text file.
     If the file extension of the given file name is '.csv', then the format of
     the output file is comma separated values, otherwise it's a basic text
     file.
  """

  # If a file name is given open it for writing - - - - - - - - - - - - - - - -
  #
  if (file_name != None):
    try:
      f = open(file_name, 'w')
    except:
      logging.exception('Can not open file "%s" for writing' % (file_name))
      raise IOError

    file_ext = file_name.split('.')[-1]  # Get the file extension

    # Check the file type (CSV or not), set a flag and write a header
    #
    if (file_ext.lower().strip() == 'csv'):
      do_csv_file = True  # Flag for CSV file output

      f.write('Rec_ID_A, Rec_ID_B, Weight, Assigned' + os.linesep)

    else:
      do_csv_file = False

      f.write('Resulting record pairs:' + os.linesep)
      f.write('-----------------------' + os.linesep)
      f.write('  Output threshold:  %f' % (threshold) + os.linesep)
      f.write('  Data set A:        %s' % (dataset_a_name) + os.linesep)
      f.write('  Data set B:        %s' % (dataset_b_name) + os.linesep)
      f.write(os.linesep)

  else:  # Print a header - - - - - - - - - - - - - - - - - - - - - - - - - - -

    do_csv_file = False

    print
    print 'Resulting record pairs:'
    print '-----------------------'
    print '  Output threshold:  %f' % (threshold)
    print '  Data set A:        %s' % (dataset_a_name)
    print '  Data set B:        %s' % (dataset_b_name)
    print

  if (threshold == None):  # Set threshold to a very negative number
    threshold == -9999999.9

  num_rec_pairs_saved = 0  # Count the number of record pairs saved
  num_rec_pairs_assig = 0  # Number of assigned record pairs

  rec_numbers = results_dict.keys()
  rec_numbers.sort()

  # For nice output, get the length (in digits) of the largest number
  #
  max_len = len(str(rec_numbers[-1])) + 1

  for rec_num in rec_numbers:  # Loop over all records in results - - - - - - -

    if (do_csv_file == True):
      csv_line = str(rec_num) + ','
      csv_line_list = []

    else:
      # Save all print lines for this record into a list
      #
      rec_line_print = [str(rec_num)]

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
          assigned = 'assigned'
          num_rec_pairs_assig += 1
        else:
          assigned = ''

        if (do_csv_file == True):
          csv_line_list.append(csv_line+'%s,%f,%s' % (rec, weight, assigned))

        else:
          rec_line_print.append(max_len*' '+' '+str(rec).rjust(max_len)+ \
                                ':  %f  %s' % (weight, assigned))
        num_rec_pairs_saved += 1

    if (do_csv_file == False):
      rec_line_print.append('')  # Append an empty line

      if (len(rec_line_print) > 2):  # There were record pairs for this record
        for l in rec_line_print:

          if (file_name != None):
            f.write(l + os.linesep)
          else:
            print l

    elif ((do_csv_file == True) and (len(csv_line_list) > 0)):

      for l in csv_line_list:  # There were record pairs for this record
        f.write(l + os.linesep)

  # Write final statistics (not for CSV files) - - - - - - - - - - - - - - - -
  #
  if (do_csv_file == False):
    if (file_name != None):
      f.write(os.linesep)
      f.write('  Number of record pairs saved:    %i' % \
              (num_rec_pairs_saved) + os.linesep)
      f.write('  Number of record pairs assigned: %i' % \
              (num_rec_pairs_assig) + os.linesep)
      f.close()
      print 'Record pair with weights written to file "%s"' % (file_name)

    else:
      print
      print '  Number of record pairs printed:  %i' % (num_rec_pairs_saved)
      print '  Number of record pairs assigned: %i' % (num_rec_pairs_assig)
      print

  if (file_name != None):
    f.close()

# =============================================================================

def time_string(seconds):
  """Helper function which returns a time in hours, minutes or seconds
     according to the value of the argument 'seconds':

     - in milliseconds if less than one second
     - in seconds if the value is less than 60 seconds
     - in minutes and seconds if the value is less than one hour
     - in hours and minutes otherwiese

     Returned is a string.
  """

  if (seconds < 0.01):  # Less than 10 milli seconds
    stringtime = '%.2f milli sec' % (seconds*1000)
  elif (seconds < 1.0):
    stringtime = '%i milli sec' % (int(seconds*1000))
  elif (seconds < 10):
    stringtime = '%.2f sec' % (seconds)
  elif (seconds < 60):
    stringtime = '%i sec' % (int(seconds))
  elif (seconds < 3600):
    min = int(seconds / 60)
    sec = seconds - min*60
    stringtime = '%i min and %i sec' % (min, sec)
  else:
    hrs = int(seconds / 3600)
    min = int((seconds - hrs *3600) / 60)
    sec = seconds - hrs *3600 - min*60
    stringtime = '%i hrs, %i min and %i sec' % (hrs, min, sec)

  return stringtime

# =============================================================================
