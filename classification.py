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
# The Original Software is: "classification.py"
# 
# The Initial Developers of the Original Software are:
#   Dr Tim Churches (Centre for Epidemiology and Research, New South Wales
#                    Department of Health)
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

"""Module classification.py - Record linkage classifiers.

   See doc strings of individual functions for detailed documentation.
"""

# =============================================================================
# Imports go here

import logging

# =============================================================================

class Classifier:
  """Base class.

     Classifiers classify weight vectors as produced by record comparators. The
     format of weight vectors is:

     [dataset_a, rec_num_a, dataset_b, rec_num_b, w_1, w_2, ... w_x]

     A comparison vector has been calculated by comparing two records. Thus the
     first four components are the identifiers of these two records. The first
     two components are the name of the data set of first record and the record
     number of the first record, while the third and fourth components are the
     data set name and record number of the second record.

     The remainder of a weight vector are the weights of the various fields
     comparisons, stored as floating point numbers.

     A classifier then uses these weights to calculate a final weight (e.g.
     the Fellegi and Sunter classifier simply adds all theweights together).

     The results of a classifier are stored in one dictionary, with record
     numbers from the first data set being the keys, and the values being the
     record numbers in the other data set and the corresponding weights. So a
     classfier results data structure can be viewed as a row oriented sparse
     matrix.
  """

  # ---------------------------------------------------------------------------

  def __init__(self, base_kwargs):
    """Constructor
    """

    # General attributes for all data sets
    #
    self.name =        ''
    self.description = ''
    self.dataset_a =   None
    self.dataset_b =   None

    self.results = {}  # Result data structure

    # Process base keyword arguments (all data set specific keywords were
    # processed in the derived class constructor)
    #
    for (keyword, value) in base_kwargs.items():
      if (keyword == 'name'):
        self.name = value
      elif (keyword == 'description'):
        self.description = value

      elif (keyword == 'dataset_a'):
        self.dataset_a = value
      elif (keyword == 'dataset_b'):
        self.dataset_b = value

      else:
        logging.exception('Illegal constructor argument keyword: %s' \
                          (str(keyword)))
        raise Exception

    # Check if the needed attributes are set  - - - - - - - - - - - - - - - - -
    #
    if (self.dataset_a == None) or (self.dataset_b == None):
      logging.exception('One or both data sets are not defined')
      raise Exception

  # ---------------------------------------------------------------------------

  def classify(self, vector):
    """Classify one weight vector.
       See implementations in derived classes for details.
    """

    logging.exception('Override abstract method in derived class')
    raise Exception

  # ---------------------------------------------------------------------------

  def classify_block(self, vector_list):
    """Classify a list of weight vectors.
    """

    if (not isinstance(vector_list, list)):
      logging.exception('Weight vector list is not a list: %s' % (str(vector)))
      raise Exception

    for vector in vector_list:
      self.classify(vector)

  # ---------------------------------------------------------------------------

  def merge(self, other_classifier_results):
    """Merge the results of another classifier into the classfier.
    """

    logging.info('  Merging classifiers result dictionaries')

    # Update dictionary - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    #
    for rec_num_a in other_classifier_results:

      rec_data = other_classifier_results[rec_num_a]

      dict_a = self.results.get(rec_num_a, {})

      # Now check if a record pair is already in the results dictionary
      #
      for (rec_num_b, weight) in rec_data.items():
        if (dict_a.has_key(rec_num_b)):
          if (dict_a[rec_num_b] == weight):  # Same weight
            logging.debug('    Record pair (%s,%s) in both result ' +\
                          'dictionaries' % (str(rec_num_a), str(rec_num_b)) + \
                          ' with same weight')
          else:
            logging.warn('Record pair (%s,%s) in both result dictionaries' % \
                  (str(rec_num_a), str(rec_num_b)) + ' with different ' + \
                  'weights: %f / %f' % (dict_a[rec_num_b], weight) + \
                  ' (keep old weight value)')
        else:
          dict_a[rec_num_b] = weight  # Insert into this results dictionary

      self.results[rec_num_a] = dict_a  # And update it

# =============================================================================

class FellegiSunterClassifier(Classifier):
  """The classical Fellegi and Sunter approach of summing all weights in a
     weight vector into the final weight and then classifying them using two
     thresholds.

     When record pairs are classified, only comparisons which have a final
     weight higher than the lower threshold are stored in the result data
     structures, i.e. only possible links and links are stored, non-links are
     discarded (but their number is counted).

     The methods 'classify_block' and 'merge' are done by the base class.
  """

  # ---------------------------------------------------------------------------

  def __init__(self, **kwargs):
    """Constructor.
    """

    # Initialise attributes
    #
    self.lower_threshold = None
    self.upper_threshold = None
    self.non_link_count =  0     # Number of comparisons resulting in non-links

    # Process all keyword arguments
    #
    base_kwargs = {}  # Dictionary, will contain unprocessed arguments for base
                      # class constructor

    for (keyword, value) in kwargs.items():
      if (keyword in ['upper','upper_threshold']):
        if (not (isinstance(value, int) or isinstance(value, float))):
          logging.exception('Argument "upper_threshold" is not a number')
          raise Exception
        self.upper_threshold = value

      elif (keyword in ['lower','lower_threshold']):
        if (not (isinstance(value, int) or isinstance(value, float))):
          logging.exception('Argument "lower_threshold" is not a number')
          raise Exception
        self.lower_threshold = value

      else:
        base_kwargs[keyword] = value

    Classifier.__init__(self, base_kwargs)  # Process base arguments

    # Check if thresholds are defined - - - - - - - - - - - - - - - - - - - - -
    #
    if (self.lower_threshold == None) or (self.upper_threshold == None):
      logging.exception('Lower and/or upper threshold not defined')
      raise Exception

    # Check if lower threshold is smaller than upper threshold  - - - - - - - -
    #
    if (self.lower_threshold >= self.upper_threshold):
      logging.exception('Lower threshold is equal to or larger than upper ' + \
            'threshold: Lower=%f, upper: %f' % \
            (self.lower_threshold, self.upper_threshold))
      raise Exception

    # A log message - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    #
    logging.info('Initialised Felligi and Sunter classifier')
    logging.info('  Name:            %s' % (str(self.name)))
    logging.info('  Data set A:      %s' % (str(self.dataset_a.name)))
    logging.info('  Data set B:      %s' % (str(self.dataset_b.name)))
    logging.info('  Lower threshold: %f' % (self.lower_threshold))
    logging.info('  Upper threshold: %f' % (self.upper_threshold))

  # ---------------------------------------------------------------------------

  def classify(self, vector):
    """Classify one weight vector.
    """

    if (not isinstance(vector, list)):
      logging.exception('Weight vector is not a list: %s' % (str(vector)))
      raise Exception

    if (vector[0] != self.dataset_a.name):
      logging.exception('Wrong data set A name in weight vector (should '+ \
            'be: %s): %s' % (str(self.dataset_a.name), str(vector[0])))

    if (vector[2] != self.dataset_b.name):
      logging.exception('Wrong data set B name in weight vector (should '+ \
            'be: %s): %s' % (str(self.dataset_b.name), str(vector[0])))

    if (not isinstance(vector[1], int)) or (vector[1] < 0):
      logging.exception('Record identifier A is not a valid number: %s' % \
            (str(vector[1])))
      raise Exception

    if (not isinstance(vector[3], int)) or (vector[3] < 0):
      logging.exception('Record identifier B is not a valid number: %s' % \
            (str(vector[3])))
      raise Exception

    rec_num_a = vector[1]
    rec_num_b = vector[3]

    # Do the Fellegi and Sunter summing - - - - - - - - - - - - - - - - - - - -
    #
    final = 0.0
    for w in vector[4:]:
      final += w

    # Now insert into the two result dictionaries if the sum is - - - - - - - -
    # higher than the lower threshold
    #
    if (final >= self.lower_threshold):

      dict_a = self.results.get(rec_num_a, {})
      dict_a[rec_num_b] = final
      self.results[rec_num_a] = dict_a

      if (final >= self.upper_threshold):
        link_status = 'link'
      else:
        link_status = 'possible link'

    else:  # A non-link, just count it
      self.non_link_count += 1
      link_status = 'non-link'

    logging.debug('    Weight vector:  %s' % (str(vector)))
    logging.debug('      Final weight: %f (%s)' % (final, link_status))

# =============================================================================

class FlexibleClassifier(Classifier):
  """This classifier allows a flexible definition of the final weight
     calculation by defining tuples containing a function and elements of the
     weight vector upon which the function is applied. The final weight is then
     calculated using another function.

     The following functions are possible:
     - min   Take the minimum value in the selected weight vector elements
     - max   Take the maximum value
     - add   Add the values in the selected weight vector elements
     - mult  Multiply the values
     - avrg  Take the average of the selected weight vector elements

     The argument 'calculate' needs to be set to a list made of tuples
     (function, weight vector elements) as given in the following example.
     The 'weight vector elements' must be a list of the integer numbers that
     correspond to the selected weight elements.

     The argument 'final_funct' must be set to one of the aboce given functions
     and it set the function that will be used to calculate the final weight
     using all the intermediate weights calculated with the tuples defined with
     the 'calculate' argument.

     my_classifier = FlexibleClassifier(name = 'My example classifier',
                                   dataset_a = my_first_dataset
                                   dataset_b = my_second_dataset,
                                   calculate = [('max',[0,1,2]),
                                                ('add',[5,6,7,8]),
                                                ('avrg',[10,11,12])],
                                 final_funct = 'sum',
                             lower_threshold = 10.0,
                             upper_threshold = 20.0)

     Similar to the classical Fellegi and Sunter approach two thresholds are
     used to classify the final weights.

     When record pairs are classified, only comparisons which have a final
     weight higher than the lower threshold are stored in the result data
     structures, i.e. only possible links and links are stored, non-links are
     discarded (but their number is counted).

     The methods 'classify_block' and 'merge' are done by the base class.
  """

  # ---------------------------------------------------------------------------

  def __init__(self, **kwargs):
    """Constructor.
    """

    # Initialise attributes
    #
    self.lower_threshold = None
    self.upper_threshold = None
    self.calculate =       None  # Definitions of how to calculate weights
    self.final_funct = None      # Definition of the final weights calculation
    self.non_link_count =  0     # Number of comparisons resulting in non-links

    # Process all keyword arguments
    #
    base_kwargs = {}  # Dictionary, will contain unprocessed arguments for base
                      # class constructor

    for (keyword, value) in kwargs.items():
      if (keyword in ['upper','upper_threshold']):
        if (not (isinstance(value, int) or isinstance(value, float))):
          logging.exception('Argument "upper_threshold" is not a number')
          raise Exception
        self.upper_threshold = value

      elif (keyword in ['lower','lower_threshold']):
        if (not (isinstance(value, int) or isinstance(value, float))):
          logging.exception('Argument "lower_threshold" is not a number')
          raise Exception
        self.lower_threshold = value

      elif (keyword in ['calc', 'calculate']):
        if (not isinstance(value, list)):
          logging.exception('Argument "calculate" is not a list')
          raise Exception
        self.calculate = value

      elif (keyword == 'final_funct'):
        if (value not in ['min','max','add','mult','avrg']):
          logging.exception('Illegal final function definition: %s' % \
                            (str(value)))
          raise Exception
        self.final_funct = value

      else:
        base_kwargs[keyword] = value

    Classifier.__init__(self, base_kwargs)  # Process base arguments

    # Check if thresholds are defined - - - - - - - - - - - - - - - - - - - - -
    #
    if (self.lower_threshold == None) or (self.upper_threshold == None):
      logging.exception('Lower and/or upper threshold not defined')
      raise Exception

    # Check if lower threshold is smaller than upper threshold  - - - - - - - -
    #
    if (self.lower_threshold >= self.upper_threshold):
      logging.exception('Lower threshold is equal to or larger than upper ' + \
            'threshold: Lower=%f, upper: %f' % \
            (self.lower_threshold, self.upper_threshold))
      raise Exception

    # Check if the calculate definition is correct  - - - - - - - - - - - - - -
    #
    if (self.calculate == None):
      logging.exception('No "calculate" functions defined')
      raise Exception

    for calc_tup in self.calculate:
      if (not (isinstance(calc_tup, tuple) or isinstance(calc_tup, list))) or \
         (len(calc_tup) != 2):
        logging.exception('Definition in "calculate" argument is not a ' + \
              'valid tuple: %s' % (str(calc_tup)))
        raise Exception
      if (calc_tup[0] not in ['min','max','avrg','add','mult']):
        logging.exception('Illegal function definition in tuple "%s"' % \
              (str(calc_tup)))
        raise Exception
      if (not isinstance(calc_tup[1], list)):
        logging.exception('Weight vector elements definition is not a list' + \
                          ': %s' % (str(calc_tup[1])))
        raise Exception

    # Check if the final function is defined  - - - - - - - - - - - - - - - - -
    #
    if (self.final_funct == None):
      logging.exception('Final function is not defined')
      raise Exception

    # A log message - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    #
    logging.info('Initialised flexible classifier')
    logging.info('  Name:                    %s' % (str(self.name)))
    logging.info('  Data set A:              %s' % (str(self.dataset_a.name)))
    logging.info('  Data set B:              %s' % (str(self.dataset_b.name)))
    logging.info('  Lower threshold:         %f' % (self.lower_threshold))
    logging.info('  Upper threshold:         %f' % (self.upper_threshold))
    logging.info('  Calculate function list: %s' % (str(self.calculate)))
    logging.info('  Final function:          %s' % (str(self.final_funct)))

  # ---------------------------------------------------------------------------

  def classify(self, vector):
    """Classify one weight vector.
    """

    if (not isinstance(vector, list)):
      logging.exception('Weight vector is not a list: %s' % (str(vector)))
      raise Exception

    if (vector[0] != self.dataset_a.name):
      logging.exception('Wrong data set A name in weight vector (should '+ \
            'be: %s): %s' % (str(self.dataset_a.name), str(vector[0])))

    if (vector[2] != self.dataset_b.name):
      logging.exception('Wrong data set B name in weight vector (should '+ \
            'be: %s): %s' % (str(self.dataset_b.name), str(vector[0])))

    if (not isinstance(vector[1], int)) or (vector[1] < 0):
      logging.exception('Record identifier A is not a valid number: %s' % \
            (str(vector[1])))
      raise Exception

    if (not isinstance(vector[3], int)) or (vector[3] < 0):
      logging.exception('Record identifier B is not a valid number: %s' % \
            (str(vector[3])))
      raise Exception

    rec_num_a = vector[1]
    rec_num_b = vector[3]

    # Do the flexible weight calculation  - - - - - - - - - - - - - - - - - - -
    #
    calc_res = []  # List of immediate results

    for (calc_funct, calc_elem) in self.calculate:
      calc_val = 0.0
      calc_vector = []
      for i in calc_elem:
        calc_vector.append(vector[4+i])

      if (calc_funct == 'min'):
        calc_res.append(min(calc_vector))
      elif (calc_funct == 'max'):
        calc_res.append(max(calc_vector))
      elif (calc_funct in ['add', 'avrg']):
        tmp_res = 0.0
        for calc_elem in calc_vector:
          tmp_res += calc_elem
        if (calc_funct == 'add'):
          calc_res.append(tmp_res)
        else:  # Calculate average
          calc_res.append((tmp_res / len(calc_vector)))
      elif (calc_funct == 'mult'):
        tmp_res = 1.0
        for calc_elem in calc_vector:
          tmp_res *= calc_elem
        calc_res.append(tmp_res)

    # Calculate final weight  - - - - - - - - - - - - - - - - - - - - - - - - -
    #
    if (self.final_funct == 'min'):
      final = min(calc_res)
    elif (self.final_funct == 'max'):
      final = max(calc_res)
    elif (self.final_funct in ['add', 'avrg']):
      final = 0.0
      for calc_elem in calc_res:
        final += calc_elem
      if (self.final_funct == 'avrg'):  # Calculate average
        final /= len(calc_res)
    elif (self.final_funct == 'mult'):
      final = 1.0
      for calc_elem in calc_res:
        final *= calc_elem   # Bug fixed 29/10/2003
                             # Thanks to (Thoetsak T.) thoetsak@ar.co.th

    # Now insert into the two result dictionaries if the sum is - - - - - - - -
    # higher than the lower threshold
    #
    if (final >= self.lower_threshold):

      dict_a = self.results.get(rec_num_a, {})
      dict_a[rec_num_b] = final
      self.results[rec_num_a] = dict_a

      if (final >= self.upper_threshold):
        link_status = 'link'
      else:
        link_status = 'possible link'

    else:  # A non-link, just count it
      self.non_link_count += 1
      link_status = 'non-link'

    logging.debug('    Weight vector %s' % (str(vector)))
    logging.debug('      Final weight: %f (%s)' % (final, link_status))

# =============================================================================
