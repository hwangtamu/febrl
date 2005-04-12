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
# The Original Software is: "qgramindex.py"
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

"""Module qgramindex.py - Class for positional q-gram inverted indexing.

   This module contains a class for a positional q-gram inverted index, and
   includes methods for building and querying such a index.
"""

# =============================================================================

# Imports go here

import logging
import cPickle as pickle
import sets
import time

import stringcmp

# =============================================================================

class PosQGramIndex:
  """Main class for the positional q-gram inverted index. Implements methods
     for building and querying such an index.

     An positional q-gram inverted index can contain q-grams of different
     length as specified when constructed by the argument 'q_gram_len_list',
     which must be a list with tuples containing values for q and corresponding
     minimum and maximum lengths (which can overlap). For example

       q_gram_len_list = [(1,(1,5)),(2,(5,20)),(3,(18,999))]

     would result in 1-grams being inserted into the inverted index for strings
     with length 1 to 5 (inclusive), 2-grams (bigrams) for strings of length 5
     to 20, and 3-grams (trigrams) for strings between 18 and 999 characters
     long.

     If a string contains more than one word, it is split and each word is
     inserted into the index separately.
  """

  # ---------------------------------------------------------------------------

  def __init__(self, **kwargs):
    """Constructor.

       Arguments needed:
       - name              A short name for the q-gram index
       - description       A longer description (can be left empty)
       - field_name        The name of a the field in the (geocode) reference
                           data set for which this index is built
       - q_gram_len_list   A list of 2-element tuples containing q and
                           corresponding minimum and maximum string lengths (as
                           shown in the example above)
       - max_edit_dist     The maximal edit distance tolerated when searching
                           for matches (this is per word if a value contains
                           more than one word)
       - pickle_file_name  If set to a string the index will be saved as pickle
                           file with the given name
       - load_pickle_file  A flag (True/False), if set to True and a pickle
                           with the given name is available then it is loaded,
                           otherwise the index is newly created
    """

    self.name =             ''
    self.description =      ''

    self.field_name = None  # The name of a the field in the geocode
                            # reference data set for which this index is built

    self.q_gram_len_list = None  # A list of 2-element tuples containing q and
                                 # corresponding minimum and maximum string
                                 # lengths

    self.max_edit_dist = None  # The maximal edit distance tolerated when
                               # searching for matches (this is per word if a
                               # value contains more than one word)

    self.inv_index =  None  # The inverted positional q-gram index
    self.value_list = None  # The list of values that were used for the index

    self.start_char = chr(1)  # Character to be added to the beginning of
                              # values before calculating positional bigrams
    self.end_char =   chr(2)  # Character to be added to the end of values
                              # before calculating positional bigrams
    self.pickle_file_name = None
    self.load_pickle_file = False

    # Process all keyword arguments
    #
    for (keyword, value) in kwargs.items():
      if (keyword == 'name'):
        self.name = value
      elif (keyword == 'description'):
        self.description = value

      elif (keyword == 'field_name'):
        if (not isinstance(value, str)):
          logging.exception('Argument "field_name" is not a string: %s' % \
                            str(typeof(value)))
          raise Exception
        self.field_name = value

      elif (keyword in ['q_gram_len_list','q_gram_list']):
        if (not isinstance(value, list)):
          logging.exception('Argument "q_gram_len_list" is not a list')
          raise Exception
        for q_elem in value:
          if (not isinstance(q_elem, tuple)) or (len(q_elem) != 2):
            logging.exception('An element of argument "q_gram_len_list" is' + \
                              'not a tuple with 2-elements: %s' % (str(value)))
            raise Exception
        self.q_gram_len_list = value

      elif (keyword in ['max_edit_dist','max_editdist','edit_dist']):
        if (not isinstance(value, int)) or (value <= 0):
          logging.exception('Argument "max_edit_dist" is not a positive ' + \
                            'integer: %s' % (str(value)))
          raise Exception
        self.max_edit_dist = value

      elif (keyword == 'pickle_file_name'):
        if (value != None):
          if (not isinstance(value, str)):
            logging.exception('Argument "pickle_file_name" is not a string' + \
                              ': %s' % (str(value)))
            raise Exception
          self.pickle_file_name = value

      elif (keyword == 'load_pickle_file'):
        if (value != None):
          if (value not in [True, False]):
            logging.exception('Argument "load_pickle_file" is not a True' + \
                              'or False: %s' % (str(value)))
            raise Exception
          self.load_pickle_file = value

      else:
        logging.exception('Illegal constructor argument keyword: "%s"' % \
                          (keyword))
        raise Exception

    # Check if the needed attributes are set  - - - - - - - - - - - - - - - - -
    #
    if (self.name == None):
      logging.exception('Positional q-gram index "name" not defined')
      raise Exception

    if (self.field_name == None):
      logging.exception('Positional q-gram index "field_name" not defined')
      raise Exception

    if (self.q_gram_len_list == None):
      logging.exception('Positional q-gram index "q_gram_len_list" not ' + \
                        'defined')
      raise Exception

    if (self.max_edit_dist == None):
      logging.exception('Maximal tolerated edit distance "max_edit_dist" ' + \
                        'not defined')
      raise Exception

    # Initialise the main data structures of the index - - - - - - - - - - - -
    #
    self.inv_index =  {}  # The inverted positional q-gram index
    self.value_list = []  # The list of values that were used for the index

    # A log message - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    #
    logging.info('Initialised positional q-gram index')
    logging.info('  Name:                            %s' % (str(self.name)))
    logging.info('  Description:                     %s' % \
                 (str(self.description)))
    logging.info('  Data set field name:             %s' % \
                 (str(self.field_name)))
    logging.info('  Maximal tolerated edit distance: %s' % \
                 (self.max_edit_dist))
    logging.info('  Build inverted q-gram index with:')
    for (q,(min_len, max_len)) in self.q_gram_len_list:
      logging.info('    %d-grams for strings with length between %d to %d' % \
                   (q, min_len, max_len))
    logging.info('  Pickle file name:                %s' % \
                 (str(self.pickle_file_name)))
    logging.info('  Load pickle file:                %s' % \
                 (str(self.load_pickle_file)))

  # ---------------------------------------------------------------------------

  def build_index(self, value_list):
    """Method to build an inverted index for the given values and q-grams.

       This method takes the input value list and creates an inverted index (a
       dictionary) with positional q-grams according to the given q-gram
       length list.

       Duplicate values in the input value list are removed.

       It stores the inverted index and a list of the cleaned input values. The
       numbers in the inverted index correspond to the index in this value
       list.

       Input values containing more than one word as split into words and these
       are processed separately.

       Creates an inverted index with q-grams as keys and values being
       dictionaries with q-gram positions as keys and sets of numbers as
       values. For example:

       self.inv_index = {'pe':{0:set(14,654,78),1:set(43,789,90)}}

       For q-grams with q > 1, (q-1) special characters are added at the start
       (ascii 1) and end (ascii 2) of string values, for example
         q = 2 -> 'peter' -> '1peter2'    (with '1' and '2' designating the
         q = 3 -> 'peter' -> '11peter22'  (correspnding ASCII characters)

       If the attribute 'self.load_pickle_file' is True and a pickle file name
       is provided in 'self.pickle_file_name' then instead of building a new
       index it is loaded from file.

       If 'self.load_pickle_file' is set to False and a pickle file name is
       given the created index will be written into the file.
    """

    # Initialise the main data structures of the index - - - - - - - - - - - -
    #
    self.inv_index =  {}  # The inverted positional q-gram index
    self.value_list = []  # The list of values that were used for the index

    # Load the index from a pickle file - - - - - - - - - - - - - - - - - - - -
    #
    if ((self.load_pickle_file == True) and (self.pickle_file_name != None)):
      try:
        pickle_file = open(self.pickle_file_name)
      except IOError:
        logging.warn('Index pickle file "%s" does not exist' % \
                     (self.pickle_file_name))
        create_index = True
        pickle_file = None
      except:
        raise Exception

      if (pickle_file != None):
        pickle_data = pickle.load(pickle_file)
        self.inv_index = pickle_data[0]
        self.value_list = pickle_data[1]
        create_index = False  # Index has been loaded, no need to create it

        logging.info('  Loaded q-gram index "%s" with %d entries' % \
                     (self.pickle_file_name, len(self.inv_index)))

    else:
      create_index = True

    # Index not loaded, so we have to create it - - - - - - - - - - - - - - - -
    #
    if (create_index == True):

      logging.info('  Build inverted q-gram index for field "%s"' % \
                 (self.field_name))

      # Collect various statistics
      #
      max_q_gram_count = 0  # Maximum number of occurences of a q-gram
      longest_value = ''    # The longest input value

      val_ind = 0  # Index number for cleaned values

      for val in value_list:  # Main loop over input values - - - - - - - - - -

        val = val.strip()

        while ('  ' in val):
          val = val.replace('  ', ' ')  # Make sure there are no double spaces
        while ('__' in val):
          val = val.replace('__', '_')  # Also remove double underlines

        if val in self.value_list:  #  Duplicate value, don't append
          logging.info('Duplicate cleaned value in input list: %s' % (val))

        else:  # Insert into value list and process for q-gram index

          self.value_list.append(val)  # Append to list of cleaned values

          val_len = len(val)

          if (val_len > len(longest_value)):
            longest_value = val

          if (' ' in val):  # More than one word in the input value
            val_list = val.split(' ')
          elif ('_' in val):  # More than one word in the input value
            val_list = val.split('_')
          else:
            val_list = [val]  # One word only

          for v in val_list:  # Process all words separately - - - - - - - - -

            v_len = len(v)

            # Loop over the q-grams in the q-gram length list, process if
            # lengths is in the corresponding interval
            #
            for (q, (min_len, max_len)) in self.q_gram_len_list:
              if ((min_len <= v_len) and (v_len <= max_len)):

                if (q > 1):  # Add leading and trailing special characters
                  q_v = (q-1)*self.start_char + v + (q-1)*self.end_char
                else:
                  q_v = v

                # Get positional q-gram list
                #
                pos_q_gram_list = self.str_to_pos_q_grams(q_v, q)

                # Insert positional q-grams into inverted index
                #
                for (q_gram, pos) in pos_q_gram_list:
                  ind_dict = self.inv_index.get(q_gram, {})
                    # Start with a new empty dictionary for a new q-gram
                  pos_set = ind_dict.get(pos, sets.Set())
                    # Start with a new empty set for a new position
                  pos_set.add(val_ind)
                  ind_dict[pos] = pos_set
                  self.inv_index[q_gram] = ind_dict
                  max_q_gram_count = max(max_q_gram_count, len(pos_set))

          val_ind += 1
          if ((val_ind % 10000) == 0):
            logging.info('  Processed %d values' % (val_ind))

      ## Test: Print complete inverted index
      ##
      #for bigr in self.inv_index:
      #  print 'Index for bigram "%s":' % (bigr)
      #  bigr_dict = self.inv_index[bigr]
      #  for pos in bigr_dict:
      #    print '  Position %d: %d' % (pos, len(bigr_dict[pos]))

      logging.info('    Number of values processed:               %d' % \
                   (val_ind))
      logging.info('    Number of different q-grams:              %d' % \
                   (len(self.inv_index)))
      logging.info('    Maximum number of occurences of a q-gram: %d' % \
                   (max_q_gram_count))
      logging.info('    Longest value (with length %d):           "%s"' % \
                   (len(longest_value), longest_value))

      # Write the created index into a pickle file - - - - - - - - - - - - - -
      #
      if ((self.load_pickle_file == False) and \
          (self.pickle_file_name != None)):
        try:
          pickle_file = open(self.pickle_file_name, 'w')
        except IOError:
          logging.warn('Can not write index pickle file "%s"' % \
                       (self.pickle_file_name))
          pickle_file = None
        except:
          raise Exception

        if (pickle_file != None):
          pickle_data = [self.inv_index, self.value_list]

          pickle.dump(pickle_data, pickle_file, pickle.HIGHEST_PROTOCOL)
          pickle_file.close()

  # ---------------------------------------------------------------------------

  def get_matches(self, in_val):
    """Method to find all matches for the given value in the inverted index up
       to the defined maximal edit distance.

       Returns a list with one or more values or an empty list if no match was
       found.
    """

    start_time = time.time()

    found_match_set_list = [] # List with sets for found value number matches

    while ('  ' in in_val):
      in_val = in_val.replace('  ', ' ') # Make sure there are no double spaces
    while ('__' in in_val):
      in_val = in_val.replace('__', '_')  # Also remove double underlines

    if (' ' in in_val):  # More than one word in the input value
      in_val_list = in_val.split(' ')
    elif ('_' in in_val):  # More than one word in the input value
      in_val_list = in_val.split('_')
    else:
      in_val_list = [in_val]  # One word only

    for in_v in in_val_list:  # Process all words separately - - - - - - - - -

      this_set = sets.Set()

      in_v_len = len(in_v)

      for (q, (min_len, max_len)) in self.q_gram_len_list:
        if ((min_len <= in_v_len) and (in_v_len <= max_len)):

          if (q > 1):  # Add leading and trailing whitespaces
            q_in_v = (q-1)*self.start_char + in_v + (q-1)*self.end_char
          else:
            q_in_v = in_v

          # Convert the input into positional q-grams
          #
          in_pos_q_gram_list = self.str_to_pos_q_grams(q_in_v, q)

          logging.debug('  Input positional %d-gram list: %s' % \
                        (q, str(in_pos_q_gram_list)))

          # Calculate minimum number of q-grams that must be in common ('count
          # filtering' according to Gravano et.al.)
          #
          #min_common_q_gram = 2*(q-1)+in_str_len - 1 - (max_edit_dist-1)*q
          #
          min_common_q_gram = max(0, (in_v_len - 1 - (self.max_edit_dist-1)*q))
          logging.debug('    Minimum number of %d-grams that must be in ' % \
                        (q) + 'common: %d' % (min_common_q_gram))

          num_in_q_gram = len(in_pos_q_gram_list)

          check_num_q_gram = num_in_q_gram

          logging.debug('    Check for %d-gram lists with length %d down to' \
                        % (q, check_num_q_gram) + ' %d' % (min_common_q_gram))

          for check_num_q_gram in range(min_common_q_gram, num_in_q_gram+1):

            # Get q-gram sub-lists of length 'check_num_q_gram'
            #
            q_gram_sub_lists = self.get_sub_lists(in_pos_q_gram_list[:],
                                                  check_num_q_gram)

            logging.debug('    Checked %d sub-lists of length %d' % \
                          (len(q_gram_sub_lists), len(q_gram_sub_lists[0])))

            for q_gram_list in q_gram_sub_lists:

              val_set = self.pos_q_gram_set_intersect(q_gram_list)

              if (val_set != None):
                logging.debug('      Found %d match(es) with score %.3f for ' \
                              % (len(val_set), \
                              len(q_gram_list)/float(num_in_q_gram)) + \
                              '%d-gram sub-list: %s' % (q, q_gram_list))

                this_set = this_set.union(val_set)

      found_match_set_list.append(this_set)

    # If there were more than one word in input, find the set intersection
    #
    if (len(found_match_set_list) > 1):
      for i in range(1,len(found_match_set_list)):
        found_match_set_list[0] = \
          found_match_set_list[0].intersection(found_match_set_list[i])

    # Calculate minimum and maximum possible lengths of matching strings
    # ('length filtering' according to Gravano et.al.)
    #
    in_val_len = len(in_val)
    min_match_len = in_val_len - self.max_edit_dist
    max_match_len = in_val_len + self.max_edit_dist
    logging.debug('  Minimum and maximum matching string lengths: %d / %d' \
                  % (min_match_len, max_match_len))

    # Retrieve all original values and do length filtering - - - - - - - - - -
    #
    match_set = sets.Set()

    for val_ind in found_match_set_list[0]:
      found_val = self.value_list[val_ind]
      found_val_len = len(found_val)

      # Length filtering by checking minimum and maximum value length
      #
      if ((found_val_len >= min_match_len) and \
          (found_val_len <= max_match_len)):
        match_set.add(found_val)
        logging.debug('      Found value "%s"' % (found_val))
      else:
        logging.debug('      Lengths filter value "%s"' % (found_val))

    logging.info('  Time used: %.3f milli-seconds' % \
                 (1000.0*(time.time()- start_time)))

    if (len(match_set) > 0):
      return list(match_set)
    else:
      return []

  # ---------------------------------------------------------------------------

  def select_best_match(self, org_val, approx_match_list, str_cmp_method):
    """Method to find the best match out of a list of approximate matches using
       an approximate string comparison method.

       Uses one of the methods available in the module stringcmp.py

       Returns a list of one or more best matches with the highest approximate
       matching score. Each match is a 2-element tuple containing the matched
       value and it's approximate matching score.

       All elements in the returned list will have the same (highest) matching
       score.
    """

    best_match_score = -1.0

    for match_val in approx_match_list:

      match_score = stringcmp.do_stringcmp(str_cmp_method, org_val, match_val)

      if (match_score > best_match_score):
        best_match_score = match_score
        best_match_list = [(match_val, match_score)]
      elif (match_score == best_match_score):
        best_match_list.append((match_val,match_score))

    logging.debug('  Found %d (of %d input) matches with best score: %s' % \
                  (len(best_match_list), len(approx_match_list), \
                   str(best_match_list)))

    return best_match_list

  # ---------------------------------------------------------------------------

  def select_best_matches(self, org_val, approx_match_list, str_cmp_method,
                          min_score):
    """Method to find the best matches with an approximate match score of at
       least 'min_score' out of a list of approximate matches using an
       approximate string comparison method.

       Uses one of the methods available in the module stringcmp.py

       Returns a list of one or more best matches with match scores equal to or
       larger than the given minimal score. Each match is a 2-element tuple
       containing the matched value and it's approximate matching score.

       The elements in the returned list are sorted according to their
       approximate matching scores, with highest scores first. If there are
       more than n best matches (all with the highest match score) then all of
       them will be returned (thus more than n matches will be returned).
    """

    best_match_list = []

    for match_val in approx_match_list:

      match_score = stringcmp.do_stringcmp(str_cmp_method, org_val, match_val)

      if (match_score >= min_score):
        best_match_list.append((match_score, match_val))

    best_match_list.sort()
    best_match_list.reverse()

    return_list = []

    for (match_score, match_val) in best_match_list:
      return_list.append((match_val, match_score))

    logging.debug('  Found %d (of %d input) best matches: %s' % \
                  (len(return_list), len(approx_match_list), str(return_list)))

    return return_list

  # ---------------------------------------------------------------------------

  def str_to_pos_q_grams(self, in_str, q):
    """Method that returns a list of positional q-grams for the given input
       string (positions starting with 0).

       q is the desired length of q-grams.

       If q is larger than the length of the input string an empty list is
       returned.
    """

    q_gram_list = []

    for i in range(0,len(in_str)-q+1):
      q_gram_list.append((in_str[i:i+q],i))

    return q_gram_list

  # ---------------------------------------------------------------------------

  def pos_q_grams_to_str(self, q_gram_list):
    """Method that returns a string compiled from the given positional q-gram
       list.

       If an empty q-gram list is given an empty string is returned.
    """

    out_str = ''

    if (q_gram_list != []):

      for i in range(0,len(q_gram_list)-1):
        out_str += q_gram_list[i][0][0]

      out_str += q_gram_list[-1][0]

    return out_str

  # ---------------------------------------------------------------------------

  def get_sub_lists(self, in_list, length):
    """Method to compute all sub-lists of the given length for the given input
       list.

       Returns a list of sub-lists, or an empty list if no sub-list of lengths
       'length' can be computed.
    """

    this_len = len(in_list)  # Get current length

    if (this_len == length):  # Desired length is input length
      return [in_list]

    elif (this_len-1 == length):  # Desired sub-lists are one element shorter

      res_list = []  # Create the desired sub-lists

      for i in range(this_len):
        res_list.append(in_list[:i]+in_list[i+1:])  # Without element i

      return res_list

    else:  # Current length is longer

      res_list = []

      for i in range(this_len):
        sub_list = in_list[:i]+in_list[i+1:]  # Remove element i from list

        this_res_list = self.get_sub_lists(sub_list, length)

        for this_sub_list in this_res_list:
          if (this_sub_list not in res_list):
            res_list = res_list + [this_sub_list]

      return res_list

  # ---------------------------------------------------------------------------

  def pos_q_gram_set_intersect(self, q_gram_list):
    """Method which finds the set intersection for the given positional q-gram
       list in the class' inverted index of values.

       Returns a set of values or None.

       Uses position filtering as described in:
         Approximate String Joins in a Database (almost) for free,
         L. Gravano, P.G. Ipeirotis, H.V. Jagadish, N. Koudas, S. Muthukrishnan
         and D. Srivastava, Proceedings of the 27th VLDB Conference, Rome, 2001
    """

    if (q_gram_list == []):  # Empty q-gram list given as input
      return None

    intersect_set = None  # Start with no set

    q_gram_list.reverse()  # Works faster for longer strings (q-grams are less
                           # common at later positions, therefore shorter sets)

    num_q_grams = len(q_gram_list)

    for (q_gram, pos) in q_gram_list:  # Loop over all positional q-grams

      # Get the dictionary for this q-gram or None if it doesn't exist in index
      #
      q_gram_dict = self.inv_index.get(q_gram, None)

      if (q_gram_dict == None):
        intersect_set = None
        break  # Doesn't exist, exit loop

      else:
        # Calculate the range of positions to consider (independent of q)
        #
        min_pos = max(pos - self.max_edit_dist, 0)
        max_pos = pos + self.max_edit_dist

        q_gram_set = sets.Set()  # Make union of record sets over positions

        for i in range(min_pos, max_pos+1):  # Include maximum position

          if (i in q_gram_dict):
            q_gram_set = q_gram_set.union(q_gram_dict[i])

        if (q_gram_set == sets.Set()):  # No value number found for this q-gram
          intersect_set = None
          break  # Empty intersection, exit loop

        # Now do intersection with final set
        #
        if (intersect_set == None):  # First time, don't intersect
          intersect_set = q_gram_set.copy()  # Make a copy

        else:  # Do intersection
          intersect_set = q_gram_set.intersection(intersect_set)

        if (len(intersect_set) == 0):
          intersect_set = None
          break  # Empty intersection, exit loop

    return intersect_set

# -----------------------------------------------------------------------------

