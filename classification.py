# =============================================================================
# classification.py - Record linkage classifiers
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
# The Original Software is "classification.py".
# The Initial Developers of the Original Software are Dr Peter Christen
# (Department of Computer Science, Australian National University), Dr Tim
# Churches (Centre for Epidemiology and Research, New South Wales Department
# of Health) and Drs Markus Hegland, Stephen Roberts and Ole Nielsen
# (Mathematical Sciences Insitute, Australian National University). Copyright
# (C) 2002 the Australian National University and others. All Rights Reserved.
# Contributors:
#
# =============================================================================

"""Module classification.py - Record linkage classifiers.

   See doc strings of individual functions for detailed documentation.
"""

# =============================================================================
# Imports go here

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
        print 'error:Illegal constructor argument keyword: %s' (str(keyword))
        raise Exception

    # Check if the needed attributes are set  - - - - - - - - - - - - - - - - -
    #
    if (self.dataset_a == None) or (self.dataset_b == None):
      print 'error:One or both data sets are not defined'
      raise Exception

  # ---------------------------------------------------------------------------

  def classify(self, vector):
    """Classify one weight vector.
       See implementations in derived classes for details.
    """

    print 'error:Override abstract method in derived class'
    raise Exception

  # ---------------------------------------------------------------------------

  def classify_block(self, vector_list):
    """Classify a list of vectors.
       See implementations in derived classes for details.
    """

    print 'error:Override abstract method in derived class'
    raise Exception

  # ---------------------------------------------------------------------------

  def merge(self, other_classifier):
    """Merge another classifier (of the same sub-class) into the classfier.
       See implementations in derived classes for details.
    """

    print 'error:Override abstract method in derived class'
    raise Exception

# =============================================================================

class FellegiSunterClassifier(Classifier):
  """The classical Fellegi and Sunter approach of summing weight vectors into
     one number and then classifying them using two thresholds.

     When record pairs are classified, only comparisons which have a final
     weight higher than the lower threshold are stored in the result data
     structures.
  """

  # ---------------------------------------------------------------------------

  def __init__(self, **kwargs):
    """Constructor.
    """

    # Initialise attributes
    #
    self.lower_threshold = None
    self.upper_threshold = None

    # Process all keyword arguments
    #
    base_kwargs = {}  # Dictionary, will contain unprocessed arguments for base
                      # class constructor

    for (keyword, value) in kwargs.items():
      if (keyword in ['upper','upper_threshold']):
        if (not (isinstance(value, int) or isinstance(value, float))):
          print 'error:Argument "upper_threshold" is not a number'
          raise Exception
        self.upper_threshold = value

      elif (keyword in ['lower','lower_threshold']):
        if (not (isinstance(value, int) or isinstance(value, float))):
          print 'error:Argument "lower_threshold" is not a number'
          raise Exception
        self.lower_threshold = value

      else:
        base_kwargs[keyword] = value

    Classifier.__init__(self, base_kwargs)  # Process base arguments

    # Check if thresholds are defined - - - - - - - - - - - - - - - - - - - - -
    #
    if (self.lower_threshold == None) or (self.upper_threshold == None):
      print 'error:Lower and/or upper threshold not defined'
      raise Exception

    # Check if lower threshold is smaller than upper threshold  - - - - - - - -
    #
    if (self.lower_threshold >= self.upper_threshold):
      print 'error:Lower threshold is equal to or larger than upper '+ \
            'threshold: Lower=%f, upper: %f' % \
            (self.lower_threshold, self.upper_threshold)
      raise Exception

  # ---------------------------------------------------------------------------

  def classify(self, vector):
    """Classify one weight vector.
    """

    if (not isinstance(vector, list)):
      print 'error:Weight vector is not a list: %s' % (str(vector))
      raise Exception

    if (vector[0] != self.dataset_a.name):
      print 'error:Wrong data set A name in weight vector '+ \
            '(should be: %s): %s' % (str(self.dataset_a.name), str(vector[0]))

    if (vector[2] != self.dataset_b.name):
      print 'error:Wrong data set B name in weight vector '+ \
            '(should be: %s): %s' % (str(self.dataset_b.name), str(vector[0]))

    if (not isinstance(vector[1], int)) or (vector[1] < 0):
      print 'error:Record identifier A is not a valid number: %s' % \
            (str(vector[1]))
      raise Exception

    if (not isinstance(vector[3], int)) or (vector[3] < 0):
      print 'error:Record identifier B is not a valid number: %s' % \
            (str(vector[3]))
      raise Exception

    rec_num_a = vector[1]
    rec_num_b = vector[3]

    # Do the Fellegi and Sunter summing - - - - - - - - - - - - - - - - - - - -
    #
    sum = 0.0
    for w in vector[4:]:
      sum += w

    # Now insert into the two result dictionaries if the sum is - - - - - - - -
    # higher than the lower threshold
    #
    if (sum >= self.lower_threshold):

      dict_a = self.results.get(rec_num_a, {})
      dict_a[rec_num_b] = sum
      self.results[rec_num_a] = dict_a

    print '3:    Weight vector %s' % (str(vector))
    print '3:      Sum: %f' % (sum)

  # ---------------------------------------------------------------------------

  def classify_block(self, vector_list):
    """Classify a list of weight vectors.
    """

    if (not isinstance(vector_list, list)):
      print 'error:Weight vector list is not a list: %s' % (str(vector))
      raise Exception

    for vector in vector_list:
      self.classify(vector)

  # ---------------------------------------------------------------------------

  def merge(self, other_classifier):
    """Merge another classifier (of the same sub-class) into the classfier.
       See implementations in derived classes for details.
    """

    if type(self) != type(other_classifier):
      print 'Error:Different classifiers, merging not possible'
      raise Exception

    # Get references to the results from the other classifier - - - - - - - - -
    #
    other_results = other_classifier.results

    # Update dictionary - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    #
    for (rec_num_a, rec_data) in other_results.items():

      dict_a = self.results.get(rec_num_a, {})
      dict_a.update(rec_data)
      self.results[rec_num_a] = dict_a

# =============================================================================
