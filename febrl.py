# =============================================================================
# febrl.py - Main module and classes for febrl projects.
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
# The Original Software is "febrl.py".
# The Initial Developers of the Original Software are Dr Peter Christen
# (Department of Computer Science, Australian National University), Dr Tim
# Churches (Centre for Epidemiology and Research, New South Wales Department
# of Health) and Drs Markus Hegland, Stephen Roberts and Ole Nielsen
# (Mathematical Sciences Insitute, Australian National University). Copyright
# (C) 2002 the Australian National University and others. All Rights Reserved.
# Contributors:
#
# =============================================================================

"""Module febrl.py - Main module and classes for febrl projects.

   TODO
   - Allow saving and loading not only of pickled project files, but also XML,
     text, compressed and uncompressed, etc.
     Then, do not restrict project files to have a '.fbl' extension, but list
     all files found in a project directory (incl. their tyes)
     Changes in Project.save() and Febrl.__init__() needed.
"""

# =============================================================================
# The following flags can be set to True or False
# They are used for testing the parallel functionalities of Febrl and should
# be set to False for normal use.

DO_PARALLEL_TEST =         False  # Perform several parallel tests of
                                  # intermediate results
SAVE_PARALLEL_TEST_FILES = False  # Write intermediate results to file for
                                  # inspection

# =============================================================================

import cPickle
import os
import sys
import time
import traceback
import types
import copy

import parallel
import indexing
import output
import lap

# =============================================================================

class Febrl:
  """Class Febrl - Main class for Febrl projects.
  """

  def __init__(self, **kwargs):
    """Constructor - Set attributes and load list of available projects.
    """

    self.version_major =      '0.2'
    self.version_minor =      ''
    self.version =            self.version_major+'.'+self.version_minor
    self.license =            'ANUOS Version 1.0'
    self.copyright =          '(C) 2002 the Australian National University' + \
                              ' and others'
    self.initial_developers = 'Dr Peter Christen (Department of Computer ' + \
                              'Science, Australian National University), ' + \
                              'Dr Tim Churches (Centre for Epidemiology '+ \
                              'and Research, New South Wales Department ' + \
                              'of Health) and Drs Markus Hegland, Stephen ' + \
                              'Roberts and Ole Nielsen (Mathematical ' + \
                              'Sciences Insitute, Australian National ' + \
                              'University)'
    self.contributors =       ''

    self.description = None
    self.febrl_path =  os.curdir  # Default set to current directory
    self.project_names = []

    # Process all keyword arguments
    #
    for (keyword, value) in kwargs.items():
      if (keyword == 'febrl_path'):
        self.febrl_path = value
      elif (keyword == 'description'):
        self.description = value

      else:
        print 'error:Illegal constructor argument keyword: "%s"' % \
              (str(keyword))
        raise Exception

    # Check if Febrl projects are available in the 'project_path' directory by
    # scanning the directory for files with '.fbr' file extension
    #
    file_list = os.listdir(self.febrl_path)
    for fn in file_list:
      file_name = fn.strip().lower()
      if (file_name[-4:] == '.fbr'):
        self.project_names.append(file_name)

  # ---------------------------------------------------------------------------

  def __str__(self):
    """Create a string representation of the Febrl object.
    """

    linesep = os.linesep

    rep = 'Febrl (Freely extensible biomedical record linkage)' + linesep
    rep += '---------------------------------------------------' + linesep
    rep += linesep
    rep += '  Version: ' + self.version + linesep
    rep += '  License: ' + self.license + linesep
    rep += '  Copyright: ' + self.copyright + linesep + linesep

    rep += '  Initial developers: ' + self.initial_developers + linesep
    rep += '  Contributors: ' + self.contributors + linesep + linesep

    rep += 'Description: ' + str(self.description) + linesep
    rep += 'Febrl path: ' + str(self.febrl_path) + linesep
    rep += 'Avaliable projects:' + linesep
    if (len(self.project_names) == 0):
      rep += '  ' + str(None) + linesep
    else:
      i = 0
      for pn in self.project_names:
        rep += str(i).rjust(3) + ': ' + pn + linesep
        i += 1
      rep += linesep

    return rep                      

  # ---------------------------------------------------------------------------

  def load_project(self, project, project_path=None):
    """Load a project from file.

       A project can either be the file name as a string or the project
       number (as stored in the list of project names).
    """

    # Check project path  - - - - - - - - - - - - - - - - - - - - - - - - - - -
    #
    if (project_path == None):  # No project path given, use default
      file_name = self.febrl_path
    else:
      file_name = project_path

    if (file_name[-1] != os.sep):
      file_name += os.sep  # Make sure path ends with a directory separator

    # Check input argument type - - - - - - - - - - - - - - - - - - - - - - - -
    #
    if (type(project) == types.IntType):
      try:
        file_name += self.project_names[project]  # Get project file name
      except:
        print 'error:Illegal project number: %s' % (str(project))
        raise Exception

    elif (type(project) == types.StringType):
      file_name += project
    else:
      print 'error:Illegal type for "project" argument, must be either of ' + \
            'type string or integer'
      raise Exception

    # Open project file and load it - - - - - - - - - - - - - - - - - - - - - -
    #
    f = open(file_name,'r')
    loaded_project = cPickle.loads(f.read())
    f.close()
    loaded_project.febrl = self

    return loaded_project

  # ---------------------------------------------------------------------------

  def new_project(self, **kwargs):
    """Create a new project object and populate it.
    """

    new_project = Project(self, **kwargs)

    return new_project

# =============================================================================

class Project:
  """Class for record linkage projects.
  """

  # ---------------------------------------------------------------------------

  def __init__(self, febrl, **kwargs):
    """Constructor - Create a new project object.
    """

    self.febrl =            febrl
    self.name =             ''
    self.description =      ''
    self.file_name =        None
    self.project_path =     febrl.febrl_path  # Inherit path
    self.block_size =       10000  # File blocking size (in number of records)

    for (keyword, value) in kwargs.items():
      if (keyword == 'name'):
        self.name = value
      elif (keyword == 'description'):
        self.description = value

      elif (keyword == 'file_name'):
        self.file_name = value
      elif (keyword == 'project_path'):
        self.project_path = value

      elif (keyword == 'block_size'):
        if (not isinstance(value, int)) and (value > 0):
          print 'error:Argument "block_size" is not a positive integer'
          raise Esception
        self.block_size = value

      else:
        print 'error:Illegal constructor argument keyword: "%s"' % \
              (str(keyword))
        raise Exception

  # ---------------------------------------------------------------------------

  def __str__(self):
    """Create a string representation for this project.
    """

    linesep = os.linesep

    rep = linesep + 'Febrl project: "'+self.name + '"' + linesep
    rep += '  Description:  ' + self.description + linesep
    rep += '  Filename:     ' + self.file_name + linesep
    rep += '  Project path: ' + self.project_path + linesep

    return rep

  # ---------------------------------------------------------------------------

  def save(self, path=None):
    """Save the project into a file (currently a pickled file).
    """

    if (path is None):
      path = self.project_path  # Take the project's path

    if (path[-1] != os.sep):
      path += os.sep  # Make sure path ends with a directory separator

    # Unset the current febrl object
    #
    save_febrl = copy.copy(self.febrl)  # Make a deep copy first
    self.febrl = None

    file_name = path + self.file_name
    f = open(file_name, 'w+')
    f.write(cPickle.dumps(self, 1))
    f.close()

    # Restore febrl object
    #
    self.febrl = save_febrl

  # ---------------------------------------------------------------------------

  def link(self, **kwargs):
    """Link the given two data set using the defined record standardisers,
       record comparators, blocking indexes and classifiers.

       Records are loaded block wise from the input data sets, then
       standardised (if the record standardisers are defined, otherwise an
       input data set is directly taken for the linkage process), linked and
       the results are printed and/or saved into the result file.

       If the argument 'first_record' is not given, it is assumed to be the
       first record in the data set (i.e. record number 0).
       Similarly, if the argument 'number_records' is not given, it is assumed
       to be all records in the input data set.

       Currently, the output can be a printed list of record pairs (if the
       argument 'output_print' is set to 'True' and/or a text file with the
       linked record numbers, if the argument 'output_file' is set to a file
       name). The output can be filtered by setting the 'output_threshold'
       (meaning all record pairs with a weight less then this threshold are not
       printed or saved).
       If future versions, it will be possible to compile an output data set.

       A histogram can be printed by setting the argument 'output_histogram' to
       'True'.

       It is also possible to apply a one-to-one assignment procedure by
       setting the argument 'output_assignment' to 'one2one'.
    """

    self.input_dataset_a =    None   # A reference to the first input data set
    self.tmp_dataset_a =      None   # A reference to the first temporary
                                     # (random access) data set
    self.input_dataset_b =    None   # A reference ot the second input data set
    self.tmp_dataset_b =      None   # A reference to the second temporary
                                     # (random access) data set
#    self.output_dataset =     None   # A reference to the output data set

    self.rec_standardiser_a = None   # Reference to a record standardiser for
                                     # the first data set (A)
    self.rec_standardiser_b = None   # Reference to a record standardiser for
                                     # the second data set (B)
    self.blocking_index_a =   None   # Reference to a blocking index for data
                                     # the first set (A)
    self.blocking_index_b =   None   # Reference to a blocking index for data
                                     # the second set (B)
    self.rec_comparator =     None   # Reference to a record comparator
    self.classifier =         None   # Reference to a weight vector classifier

    self.first_record_a =     None   # Number of the first record to process
                                     # in the first data set (A)
    self.number_records_a =   None   # Number of records to process in the
                                     # first data set (A)
    self.first_record_b =     None   # Number of the first record to process
                                     # in the second data set (B)
    self.number_records_b =   None   # Number of records to process in the
                                     # second data set (B)

    self.output_print =     False    # Flag, set to True or False (default) if
                                     # record pairs should be printed
    self.output_histogram = False    # Flag, set to True or False (default) if 
                                     # a histogram of weights should be printed
    self.output_file =      None     # Set to a file name if results are to be
                                     # saved
    self.output_threshold = None     # Set to a weight threshold (only record
                                     # pairs with weights equal to or above
                                     # will be saved and or printed)
    self.output_assignment = None    # Set to 'one2one' if one-to-one
                                     # assignment should be forced
                                     # (default: None)

    for (keyword, value) in kwargs.items():
      if (keyword == 'input_dataset_a'):
        self.input_dataset_a = value
      elif (keyword == 'input_dataset_b'):
        self.input_dataset_b = value
      elif (keyword == 'tmp_dataset_a'):
        self.tmp_dataset_a = value
      elif (keyword == 'tmp_dataset_b'):
        self.tmp_dataset_b = value
#      elif (keyword == 'output_dataset'):
#        self.output_dataset = value

      elif (keyword == 'rec_standardiser_a'):
        self.rec_standardiser_a = value
      elif (keyword == 'rec_standardiser_b'):
        self.rec_standardiser_b = value
      elif (keyword == 'rec_comparator'):
        self.rec_comparator = value
      elif (keyword == 'blocking_index_a'):
        self.blocking_index_a = value
      elif (keyword == 'blocking_index_b'):
        self.blocking_index_b = value
      elif (keyword == 'classifier'):
        self.classifier = value

      elif (keyword == 'first_record_a'):
        if (not isinstance(value, int)) or (value < 0):
          print 'error:Argument "first_record_a" is not an integer number'
          raise Exception
        self.first_record_a = value
      elif (keyword == 'first_record_b'):
        if (not isinstance(value, int)) or (value < 0):
          print 'error:Argument "first_record_b" is not an integer number'
          raise Exception
        self.first_record_b = value
      elif (keyword == 'number_records_a'):
        if (not isinstance(value, int)) or (value <= 0):
          print 'error:Argument "number_records_a" is not a positive '+ \
                'integer number'
          raise Exception
        self.number_records_a = value
      elif (keyword == 'number_records_b'):
        if (not isinstance(value, int)) or (value <= 0):
          print 'error:Argument "number_records_b" is not a positive '+ \
                'integer number'
          raise Exception
        self.number_records_b = value

      elif (keyword == 'output_print'):
        if (value not in [True, False]):
          print 'error:Argument "output_print" must be "True" or "False"'
          raise Exception
        self.output_print = value
      elif (keyword == 'output_histogram'):
        if (value not in [True, False]):
          print 'error:Argument "output_histogram" must be "True" or "False"'
          raise Exception
        self.output_histogram = value
      elif (keyword == 'output_file'):
        if (not isinstance(value, str)) or (value == ''):
          print 'error:Argument "output_file" is not a valid string: %s' % \
                (str(value))
          raise Exception
        self.output_file = value
      elif (keyword == 'output_threshold'):
        if (not (isinstance(value, int) or isinstance(value, float))):
          print 'error:Argument "output_threshold" is not a number: %s' % \
                (str(value))
        self.output_threshold = value
      elif (keyword == 'output_assignment'):
        if (value not in ['one2one', None]):
          print 'error:Illegal value for argument "output_assignment": %s' % \
                (str(value))
          raise Exception
        else:
          self.output_assignment = value

      else:
        print 'error:Illegal constructor argument keyword: "%s"' % \
              (str(keyword))
        raise Exception

    # Check if the needed attributes are set  - - - - - - - - - - - - - - - - -
    #
    if (self.input_dataset_a == None):
      print 'error:Input data set A is not defined'
      raise Exception
    if (self.input_dataset_b == None):
      print 'error:Input data set B is not defined'
      raise Exception

    if (self.tmp_dataset_a == None):
      print 'error:Temporary data set A is not defined'
      raise Exception
    if (self.tmp_dataset_b == None):
      print 'error:Temporary data set B is not defined'
      raise Exception

#    if (self.output_dataset == None):
#      print 'error:Output data set is not defined'
#      raise Exception

    # Make sure either output_print is True or the output_file is defined
    #
    if (self.output_print == False) and (self.output_file == None):
      print 'error:No ouput of results (record pairs/results file) is defined.'
      raise Exception
    #
    # Code above to be removed once output data set functionality implemented

    if (self.first_record_a == None):
      self.first_record_a = 0  # Take default first record in data set
    if (self.first_record_b == None):
      self.first_record_b = 0  # Take default first record in data set

    if (self.number_records_a == None):
      self.number_records_a = input_dataset_a.num_records  # Take all records
    if (self.number_records_b == None):
      self.number_records_b = input_dataset_b.num_records  # Take all records

    if (self.rec_comparator == None):
      print 'error:No record comparator defined'
      raise Exception

    if (self.blocking_index_a == None):
      print 'error:No blocking index for data set A defined'
      raise Exception
    if (self.blocking_index_b == None):
      print 'error:No blocking index for data set B defined'
      raise Exception

    if (self.classifier == None):
      print 'error:No classifier defined'
      raise Exception

    if (self.rec_standardiser_a != None):
      if (self.rec_standardiser_a.input_dataset != self.input_dataset_a):
        print 'error:Illegal input data set definition in record '+ \
              'standardiser A: %s (should be: %s)' % \
              (str(self.rec_standardiser_a.input_dataset.name), \
               str(self.input_dataset_a.name))
        raise Exception
      if (self.rec_standardiser_a.output_dataset != self.tmp_dataset_a):
        print 'error:Illegal output data set definition in record '+ \
              'standardiser A: %s (should be: %s)' % \
              (str(self.rec_standardiser_a.output_dataset.name), \
               str(self.tmp_dataset_a.name))
        raise Exception

    else:  # No standardiser for data set A defined, so field names in input
           # and temporary data sets must be the same
      input_field_name_list = self.input_dataset_a.fields.keys()
      input_field_name_list.sort()
      tmp_field_name_list = self.tmp_dataset_a.fields.keys()
      tmp_field_name_list.sort()

      if (input_field_name_list != tmp_field_name_list):
        print 'error:Field names differ in input and temporary data sets ' + \
              '(with no record standardiser for data set A defined)'

    if (self.rec_standardiser_b != None):
      if (self.rec_standardiser_b.input_dataset != self.input_dataset_b):
        print 'error:Illegal input data set definition in record '+ \
              'standardiser B: %s (should be: %s)' % \
              (str(self.rec_standardiser_b.input_dataset.name), \
               str(self.input_dataset_b.name))
        raise Exception
      if (self.rec_standardiser_b.output_dataset != self.tmp_dataset_b):
        print 'error:Illegal output data set definition in record '+ \
              'standardiser B: %s (should be: %s)' % \
              (str(self.rec_standardiser_b.output_dataset.name), \
               str(self.tmp_dataset_b.name))
        raise Exception

    else:  # No standardiser for data set B defined, so field names in input
           # and temporary data sets must be the same
      input_field_name_list = self.input_dataset_b.fields.keys()
      input_field_name_list.sort()
      tmp_field_name_list = self.tmp_dataset_b.fields.keys()
      tmp_field_name_list.sort()

      if (input_field_name_list != tmp_field_name_list):
        print 'error:Field names differ in input and temporary data sets ' + \
              '(with no record standardiser for data set B defined)'

    if (self.blocking_index_a.dataset != self.tmp_dataset_a):
      print 'error:Illegal data set definition in blocking index A'
      raise Exception
    if (self.blocking_index_b.dataset != self.tmp_dataset_b):
      print 'error:Illegal data set definition in blocking index B'
      raise Exception

    if (self.rec_comparator.dataset_a != self.tmp_dataset_a) or \
       (self.rec_comparator.dataset_b != self.tmp_dataset_b):
      print 'error:Illegal data set definition in record comparator'
      raise Exception

    if (self.classifier.dataset_a != self.tmp_dataset_a) or \
       (self.classifier.dataset_b != self.tmp_dataset_b):
      print 'error:Illegal data set definition in classifier'
      raise Exception

    total_time = time.time()  # Get current time

    print '1:'
    print '1:*** Link data set: %s with data set: %s ***' % \
          (self.input_dataset_a.name, self.input_dataset_b.name)
    print '1:'

    print '1:'
    print '1:  Step 1: Load and standardise records, build blocking indexes'

    step_1_time = time.time()  # Get current time
    step_1_comm_time = 0.0  # Time for communication in step 1

    print '1:    Step 1a: Process data set A (%s)' % \
          (self.input_dataset_a.name)
    print '1:'

    # Do cleaning and standardisation for data set A  - - - - - - - - - - - - -
    #
    input_rec_counter = self.first_record_a  # Current record pointer

    block_cnt = 0  # A round robin block counter, used for parallelism

    # Load records in a blocked fashion - - - - - - - - - - - - - - - - - - - -

    while (input_rec_counter < (self.first_record_a + self.number_records_a)):

      block_size = min(self.block_size,
               ((self.first_record_a + self.number_records_a) - \
                input_rec_counter))

      # Distribute blocks equally to all processors
      #
      if ((block_cnt % parallel.size()) == parallel.rank()):

        # Load original records from input data set
        #
        in_recs = self.input_dataset_a.read_records(input_rec_counter,
                                                  block_size)
        print '1:    Loaded records %i to %i' % \
              (input_rec_counter, input_rec_counter+block_size)

        # Standardise them if a standardiser is defined
        #
        if (self.rec_standardiser_a != None):
          clean_recs = self.rec_standardiser_a.standardise_block(in_recs)
          print '1:    Standardised records %i to %i' % \
                (input_rec_counter, input_rec_counter+block_size)
        else:
          clean_recs = in_recs  # Take the original records directly

        # Store records in temporary data set
        #
        self.tmp_dataset_a.write_records(clean_recs)

        # Insert records into the blocking index
        #
        self.blocking_index_a.build(clean_recs)

        # If Febrl is run in parallel, send cleaned records to process 0
        #
        if (parallel.rank() > 0):
          tmp_time = time.time()
          parallel.send(clean_recs, 0)
          step_1_comm_time += (time.time() - tmp_time)

      # If Febrl is run in parallel, process 0 receives cleaned records
      #
      if (parallel.rank() == 0) and (block_cnt % parallel.size() != 0):

        p = (block_cnt % parallel.size())  # Process number to receive from
        tmp_time = time.time()
        tmp_recs = parallel.receive(p)
        step_1_comm_time += (time.time() - tmp_time)

        # Store the received records into the temporary data set
        #
        self.tmp_dataset_a.write_records(tmp_recs)

      input_rec_counter += block_size  # Increment current record pointer
      block_cnt += 1

    # If Febrl is run in parallel, collect blocking index in process 0  - - - -
    #
    if (parallel.rank() == 0):
      for p in range(1, parallel.size()):
        tmp_time = time.time()
        tmp_index = parallel.receive(p)
        step_1_comm_time += (time.time() - tmp_time)
        self.blocking_index_a.merge(tmp_index)
        print '1:  Received index from process %i' % (p)

    else:
      tmp_time = time.time()
      parallel.send(self.blocking_index_a, 0) # Send local index to process 0
      step_1_comm_time += (time.time() - tmp_time)
      print '1:  Sent index to process 0'

    # Compact the blocking index on process 0 - - - - - - - - - - - - - - - - -
    #
    if (parallel.rank() == 0):
      self.blocking_index_a.compact()

    # If run in parallel, broadcast the index from process 0  - - - - - - - - -
    #
    if (parallel.size() > 1):
      if (parallel.rank() == 0):
        for p in range(1, parallel.size()):
          tmp_time = time.time()
          parallel.send(self.blocking_index_a, p)
          step_1_comm_time += (time.time() - tmp_time)
          print '1:    Sent compacted index to process %i' % (p)

      else:
        tmp_time = time.time()
        self.blocking_index_a = parallel.receive(0)
        step_1_comm_time += (time.time() - tmp_time)
        print '1:    Received compacted index from process 0'

    # If run in parallel temporary data set needs to be sent to all processes -
    #
    if (parallel.size() > 1):
      if (parallel.rank() == 0):
        for p in range(1, parallel.size()):
          tmp_time = time.time()
          parallel.send(self.tmp_dataset_a, p)
          step_1_comm_time += (time.time() - tmp_time)
          print '1:    Sent temporary data set %i' % (p)

      else:
        tmp_time = time.time()
        self.tmp_dataset_a = parallel.receive(0)
        step_1_comm_time += (time.time() - tmp_time)
        print '1:    Received temporary data set from process 0'

    #################### START PARALLEL TEST CODE #############################
    # Save temporary data sets and indexes to files (on all processes)
    #
    if (SAVE_PARALLEL_TEST_FILES == True):
      f = open('tmp_data_set_a-link-'+str(parallel.rank())+'-'+ \
               str(parallel.size()),'w')
      tmp_list = self.tmp_dataset_a.dict.keys()
      tmp_list.sort()

      for r in tmp_list:
        rec = self.tmp_dataset_a.dict[r]
        rec_items = rec.items()
        rec_items.sort()
        rec = str(r)+': '+str(rec_items)
        f.write(rec+os.linesep)
      f.close()

      f = open('indexes_a-link-'+str(parallel.rank())+'-'+ \
               str(parallel.size()),'w')

      for i in range (self.blocking_index_a.num_indexes):
        tmp_index = self.blocking_index_a.index[i].keys()
        tmp_index.sort()

        for bi in tmp_index:
          ind = self.blocking_index_a.index[i][bi]
          ind_list = ind.items()
          ind_list.sort()

          ii = str(i)+'_'+str(bi)+': '+str(ind_list)
          f.write(ii+os.linesep)
      f.close()

    #################### END PARALLEL TEST CODE ###############################

    if (parallel.size() > 1):
      parallel.Barrier()  # Make sure all processes are here

    # Do cleaning and standardisation for data set B  - - - - - - - - - - - - -
    #
    print '1:'
    print '1:    Step 1b: Process data set B (%s)' % \
          (self.input_dataset_b.name)
    print '1:'

    input_rec_counter = self.first_record_b  # Current record pointer

    block_cnt = 0  # A round robin block counter, used for parallelism

    # Load records in a blocked fashion - - - - - - - - - - - - - - - - - - - -

    while (input_rec_counter < (self.first_record_b + self.number_records_b)):

      block_size = min(self.block_size,
               ((self.first_record_b + self.number_records_b) - \
                input_rec_counter))

      # Distribute blocks equally to all processors
      #
      if ((block_cnt % parallel.size()) == parallel.rank()):

        # Load original records from input data set
        #
        in_recs = self.input_dataset_b.read_records(input_rec_counter,
                                                  block_size)
        print '1:    Loaded records %i to %i' % \
              (input_rec_counter, input_rec_counter+block_size)

        # Standardise them if a standardiser is defined
        #
        if (self.rec_standardiser_b != None):
          clean_recs = self.rec_standardiser_b.standardise_block(in_recs)
          print '1:    Standardised records %i to %i' % \
                (input_rec_counter, input_rec_counter+block_size)
        else:
          clean_recs = in_recs  # Take the original records directly

        # Store records in temporary data set
        #
        self.tmp_dataset_b.write_records(clean_recs)

        # Insert records into the blocking index
        #
        self.blocking_index_b.build(clean_recs)

        # If Febrl is run in parallel, send cleaned records to process 0
        #
        if (parallel.rank() > 0):
          tmp_time = time.time()
          parallel.send(clean_recs, 0)
          step_1_comm_time += (time.time() - tmp_time)

      # If Febrl is run in parallel, process 0 receives cleaned records
      #
      if (parallel.rank() == 0) and (block_cnt % parallel.size() != 0):

        p = (block_cnt % parallel.size())  # Process number to receive from
        tmp_time = time.time()
        tmp_recs = parallel.receive(p)
        step_1_comm_time += (time.time() - tmp_time)

        # Store the received records into the temporary data set
        #
        self.tmp_dataset_b.write_records(tmp_recs)

      input_rec_counter += block_size  # Increment current record pointer

      block_cnt += 1

    # If Febrl is run in parallel, collect blocking index in process 0  - - - -
    #
    if (parallel.rank() == 0):
      for p in range(1, parallel.size()):
        tmp_time = time.time()
        tmp_index = parallel.receive(p)
        step_1_comm_time += (time.time() - tmp_time)
        self.blocking_index_b.merge(tmp_index)
        print '1:  Received index from process %i' % (p)

    else:
      tmp_time = time.time()
      parallel.send(self.blocking_index_b, 0) # Send local index to process 0
      step_1_comm_time += (time.time() - tmp_time)
      print '1:  Sent index to process 0'

    # Compact the blocking index on process 0 - - - - - - - - - - - - - - - - -
    #
    if (parallel.rank() == 0):
      self.blocking_index_b.compact()

    # If run in parallel, broadcast the index from process 0  - - - - - - - - -
    #
    if (parallel.size() > 1):
      if (parallel.rank() == 0):
        for p in range(1, parallel.size()):
          tmp_time = time.time()
          parallel.send(self.blocking_index_b, p)
          step_1_comm_time += (time.time() - tmp_time)
          print '1:    Sent compacted index to process %i' % (p)

      else:
        tmp_time = time.time()
        self.blocking_index_b = parallel.receive(0)
        step_1_comm_time += (time.time() - tmp_time)
        print '1:    Received compacted index from process 0'

    # If run in parallel temporary data set needs to be sent to all processes -
    #
    if (parallel.size() > 1):
      if (parallel.rank() == 0):
        for p in range(1, parallel.size()):
          tmp_time = time.time()
          parallel.send(self.tmp_dataset_b, p)
          step_1_comm_time += (time.time() - tmp_time)
          print '1:    Sent temporary data set %i' % (p)

      else:
        tmp_time = time.time()
        self.tmp_dataset_b = parallel.receive(0)
        step_1_comm_time += (time.time() - tmp_time)
        print '1:    Received temporary data set from process 0'

    step_1_time = time.time() - step_1_time  # Calculate time for step 1

    #################### START PARALLEL TEST CODE #############################
    # Save temporary data sets and indexes to files (on all processes)
    #
    if (SAVE_PARALLEL_TEST_FILES == True):
      f = open('tmp_data_set_b-link-'+str(parallel.rank())+'-'+ \
               str(parallel.size()),'w')
      tmp_list = self.tmp_dataset_b.dict.keys()
      tmp_list.sort()

      for r in tmp_list:
        rec = self.tmp_dataset_b.dict[r]
        rec_items = rec.items()
        rec_items.sort()
        rec = str(r)+': '+str(rec_items)
        f.write(rec+os.linesep)
      f.close()

      f = open('indexes_b-link-'+str(parallel.rank())+'-'+ \
               str(parallel.size()),'w')

      for i in range (self.blocking_index_b.num_indexes):
        tmp_index = self.blocking_index_b.index[i].keys()
        tmp_index.sort()

        for bi in tmp_index:
          ind = self.blocking_index_b.index[i][bi]
          ind_list = ind.items()
          ind_list.sort()

          ii = str(i)+'_'+str(bi)+': '+str(ind_list)
          f.write(ii+os.linesep)
      f.close()

    #################### END PARALLEL TEST CODE ###############################

    # End of step 1 - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    #
    if (parallel.size() > 1):
      parallel.Barrier()  # Make sure all processes are here

    print '1:'
    print '1:  Step 2: Perform linkage within blocks'

    step_2_time = time.time()  # Get current time
    step_2_comm_time = 0.0  # Time for communication in step 2

    # Get the record pairs which have to be compared  - - - - - - - - - - - - -
    #
    tmp_time = time.time()
    rec_pair_list = indexing.linkage_rec_pairs(self.blocking_index_a,
                                               self.blocking_index_b)

    print '1:    Built record pair list (time used %.2f sec)' % \
          (time.time()-tmp_time)

    rec_pair_list.sort()

#    if (parallel.size() > 1) or (DO_PARALLEL_TEST == True) or \
#       (SAVE_PARALLEL_TEST_FILES == True):
#      rec_pair_list.sort()  # Needed for parallel runs only because of round
#                            # robin fashion of work distribution (and list can
#                            # have different sequence on different processes

    #################### START PARALLEL TEST CODE #############################
    # Check if the rec_pair_list are the same on all processes
    #
    if (parallel.size() > 1) and (DO_PARALLEL_TEST == True):
      if (parallel.rank() == 0):

        for p in range(1, parallel.size()):
          tmp_list = parallel.receive(p)
          if (len(tmp_list) != len(rec_pair_list)):
            print 'warning:Record pair lists have differnt length on ' + \
                  'process 0 (%i) and %i (%i)' % \
                  (len(rec_pair_list), p, len(tmp_list))

          for i in range(len(rec_pair_list)):
            if (rec_pair_list[i] != tmp_list[i]):  # Element differs
              print 'warning:Record pair lists differ in position %i' % (i) + \
                    ' on process 0 (%s) and %i (%s)' % \
                    (p, rec_pair_list[i], tmp_list[i])
      else:
        parallel.send(rec_pair_list,0)

      parallel.Barrier()

    # Save record pair lists to files (all processes)
    #
    if (SAVE_PARALLEL_TEST_FILES == True):
      f = open('rec_pair_list-link-'+str(parallel.rank())+'-'+ \
               str(parallel.size()),'w')
      for r in rec_pair_list:
        f.write(r+os.linesep)
      f.close()

    #################### END PARALLEL TEST CODE ###############################

    weight_vector_dict = {}  # Dictionary with the comparison vectors

    compare_time =  0.0
    classify_time = 0.0

    num_rec_pairs = len(rec_pair_list)
    rec_pair_cnt = 0  # Loop counter

    # Compare records, and distribute comparisons equally to all processes  - -
    #
    for rec_pair in rec_pair_list:

      if ((rec_pair_cnt % parallel.size()) == parallel.rank()):

        [rec_num_a,rec_num_b] = rec_pair.split('_')
        print '2:      Compare records %s with %s' % (rec_num_a, rec_num_b)

        # Read the records from the data set
        #
        rec_a = self.tmp_dataset_a.read_record(int(rec_num_a))
        rec_b = self.tmp_dataset_b.read_record(int(rec_num_b))

        # Compare the two records
        #
        tmp_time = time.time()
        w_vector = self.rec_comparator.compare(rec_a, rec_b)
        compare_time += (time.time() - tmp_time)

        # Save the weight vector in a dictionary (used later for output)
        #
        if (weight_vector_dict.has_key(rec_pair)):
          print 'warning:This should never happen: Record pair %s ' % \
                (rec_pair) + ' already in weight vector dictionary'
        weight_vector_dict[rec_pair] = w_vector

        tmp_time = time.time()
        self.classifier.classify(w_vector)  # Classify the weight vector
        classify_time += (time.time() - tmp_time)

      rec_pair_cnt += 1

      # Progress report every 10 per cent
      #
      if (rec_pair_cnt > 0) and ((rec_pair_cnt % int(num_rec_pairs/10)) == 0):
        print '1:    %i/%i record pairs compared' % \
              (rec_pair_cnt, num_rec_pairs)
        print '1:      Average comparison time:     %.6f' % \
              (compare_time/rec_pair_cnt)
        print '1:      Average classification time: %.6f' % \
              (classify_time/rec_pair_cnt)

    # Now gather classifiers on process 0 and merge - - - - - - - - - - - - - -
    #
    if (parallel.size() > 1):
      if (parallel.rank() == 0):
        for p in range(1, parallel.size()):
          tmp_time = time.time()
          tmp_classifier = parallel.receive(p)
          step_2_comm_time += (time.time() - tmp_time)
          self.classifier.merge(tmp_classifier)
          print '1:  Received classifier from process %i and merged it' % (p)

      else:
        tmp_time = time.time()
        parallel.send(self.classifier, 0) # Send local classifier to process 0
        step_2_comm_time += (time.time() - tmp_time)
        print '1:  Sent classifier to process 0'

    # Also gather weight vectors on process 0 and merge - - - - - - - - - - - -
    #
    if (parallel.size() > 1):
      if (parallel.rank() == 0):
        for p in range(1, parallel.size()):
          tmp_time = time.time()
          tmp_weight_vector_dict = parallel.receive(p)
          step_2_comm_time += (time.time() - tmp_time)
          weight_vector_dict.update(tmp_weight_vector_dict)
          print '1:  Received weight vectors from process %i' % (p)

      else:
        tmp_time = time.time()
        parallel.send(weight_vector_dict, 0) # Send local weight vectors
        step_2_comm_time += (time.time() - tmp_time)
        print '1:  Sent weight vectors to process 0'

    #################### START PARALLEL TEST CODE #############################
    # Save classifiers and weight vectors to files (only process 0)
    #
    if (SAVE_PARALLEL_TEST_FILES == True) and (parallel.rank() == 0):
      tmp_list = self.classifier.results.keys()
      tmp_list.sort()
      f = open('classifier_results_dict-link-'+str(parallel.rank())+'-'+ \
               str(parallel.size()),'w')
      for c in tmp_list:
        res = self.classifier.results[c].items()
        res.sort()
        ce = str(c)+': '+str(res)
        f.write(ce+os.linesep)
      f.close()

      tmp_list = weight_vector_dict.keys()
      tmp_list.sort()
      f = open('weight_vector_dict-link-'+str(parallel.rank())+'-'+ \
               str(parallel.size()),'w')
      for v in tmp_list:
        wv = v+': '+str(weight_vector_dict[v])
        f.write(wv+os.linesep)
      f.close()

    #################### END PARALLEL TEST CODE ###############################

    print '1:Linkage done, totally %i comparisons' % (num_rec_pairs)

    step_2_time = time.time() - step_2_time  # Calculate time for step 2

    # Output the results  - - - - - - - - - - - - - - - - - - - - - - - - - - -

    if (parallel.rank() == 0):  # Only processor 0 prints results

      # Get the results dictionary with all the record pairs and their weights
      #
      results_dict = self.classifier.results

      if (self.output_histogram == True):  # Print a weights histogram

        output.print_histogram(results_dict)

      if (self.output_assignment != None):  # An output assignment is defined

        if (self.output_assignment == 'one2one'):

          # Do a one-to-one assignment on the classifier results dict
          #
          o2o_results_dict = lap.do_lap('lapmod', results_dict, \
                                        'linkage', self.output_threshold)
      else:  # No one-to-one assignment, set o2o result to None
        o2o_results_dict = None

      if (self.output_print == True):  # Print resulting record pairs

        output.print_record_pairs(self.tmp_dataset_a, self.tmp_dataset_b, \
                                  results_dict, o2o_results_dict, \
                                  self.output_threshold)

      if (self.output_file != None):  # Save results into a file

        output.save(self.tmp_dataset_a.name, self.tmp_dataset_b.name, \
                    self.output_file, results_dict, o2o_results_dict, \
                    self.output_threshold)
 
    print '1:Done.'

    total_time = time.time() - total_time  # Calculate total time

    parallel.Barrier()  # Wait here for all processes - - - - - - - - - - - - -

    print '1:Total time needed for linkage of %i' % (self.number_records_a) + \
          ' records with %i records: %.3f' % \
          (self.number_records_b, total_time) 

    print '1:  Time for step 1: %.3f' % (step_1_time)
    print '1:  Time for step 2: %.3f' % (step_2_time)
    print '1:  Time for communication in step 1: %.3f' % (step_1_comm_time)
    print '1:  Time for communication in step 2: %.3f' % (step_2_comm_time)

    parallel.Finalize()

  # ---------------------------------------------------------------------------

  def deduplicate(self, **kwargs):
    """Deduplicate the given data set using the defined record standardiser,
       record comparators, blocking indexes and classifiers.

       Records are loaded block wise from the input data set, then standardised
       (if the record standardiser is defined, otherwise the input data set is
       directly deduplicated), linked and the results are printed and/or saved
       into the result file.

       If the argument 'first_record' is not given, it is assumed to be the
       first record in the data set (i.e. record number 0).
       Similarly, if the argument 'number_records' is not given, it is assumed
       to be all records in the input data set.

       Currently, the output can be a printed list of record pairs (if the
       argument 'output_print' is set to 'True' and/or a text file with the
       linked record numbers, if the argument 'output_file' is set to a file
       name). The output can be filtered by setting the 'output_threshold'
       (meaning all record pairs with a weight less then this threshold are not
       printed or saved).
       If future versions, it will be possible to compile an output data set.

       A histogram can be printed by setting the argument 'output_histogram' to
       'True'.

       It is also possible to apply a one-to-one assignment procedure by
       setting the argument 'output_assignment' to 'one2one'.
    """

    self.input_dataset =    None   # A reference ot the (raw) input data set
    self.tmp_dataset =      None   # A reference to a temporary (random access)
                                   # data set
#    self.output_dataset =   None   # A reference to the output data set

    self.rec_standardiser = None   # Reference to a record standardiser
    self.rec_comparator =   None   # Reference to a record comparator
    self.blocking_index =   None   # Reference to a blocking index
    self.classifier =       None   # Reference to a weight vector classifier

    self.first_record =     None   # Number of the first record to process
    self.number_records =   None   # Number of records to process

    self.output_print =     False  # Flag, set to True or False (default) if
                                   # record pairs should be printed
    self.output_histogram = False  # Flag, set to True or False (default) if a
                                   # histogram of weights should be printed
    self.output_file =      None   # Set to a file name if results are to be
                                   # saved
    self.output_threshold = None   # Set to a weight threshold (only record
                                   # pairs with weights equal to or above will
                                   # be saved and or printed)
    self.output_assignment = None  # Set to 'one2one' if one-to-one assignment
                                   # should be forced (default: None)

    for (keyword, value) in kwargs.items():
      if (keyword == 'input_dataset'):
        self.input_dataset = value
      elif (keyword == 'tmp_dataset'):
        self.tmp_dataset = value
#      elif (keyword == 'output_dataset'):
#        self.output_dataset = value

      elif (keyword == 'rec_standardiser'):
        self.rec_standardiser = value
      elif (keyword == 'rec_comparator'):
        self.rec_comparator = value
      elif (keyword == 'blocking_index'):
        self.blocking_index = value
      elif (keyword == 'classifier'):
        self.classifier = value

      elif (keyword == 'first_record'):
        if (not isinstance(value, int)) or (value < 0):
          print 'error:Argument "first_record" is not an integer number'
          raise Exception
        self.first_record = value
      elif (keyword == 'number_records'):
        if (not isinstance(value, int)) or (value <= 0):
          print 'error:Argument "number_records" is not a positive integer '+ \
                'number'
          raise Exception
        self.number_records = value

      elif (keyword == 'output_print'):
        if (value not in [True, False]):
          print 'error:Argument "output_print" must be "True" or "False"'
          raise Exception
        self.output_print = value
      elif (keyword == 'output_histogram'):
        if (value not in [True, False]):
          print 'error:Argument "output_histogram" must be "True" or "False"'
          raise Exception
        self.output_histogram = value
      elif (keyword == 'output_file'):
        if (not isinstance(value, str)) or (value == ''):
          print 'error:Argument "output_file" is not a valid string: %s' % \
                (str(value))
          raise Exception
        self.output_file = value
      elif (keyword == 'output_threshold'):
        if (not (isinstance(value, int) or isinstance(value, float))):
          print 'error:Argument "output_threshold" is not a number: %s' % \
                (str(value))
        self.output_threshold = value
      elif (keyword == 'output_assignment'):
        if (value not in ['one2one', None]):
          print 'error:Illegal value for argument "output_assignment": %s' % \
                (str(value))
          raise Exception
        else:
          self.output_assignment = value

      else:
        print 'error:Illegal constructor argument keyword: "%s"' % \
              (str(keyword))
        raise Exception

    # Check if the needed attributes are set  - - - - - - - - - - - - - - - - -
    #
    if (self.input_dataset == None):
      print 'error:Input data set is not defined'
      raise Exception

    if (self.tmp_dataset == None):
      print 'error:Temporary data set is not defined'
      raise Exception

#    if (self.output_dataset == None):
#      print 'error:Output data set is not defined'
#      raise Exception

    # Make sure either output_print is True or the output_file is defined
    #
    if (self.output_print == False) and (self.output_file == None):
      print 'error:No ouput of results (record pairs/results file) is defined.'
      raise Exception
    #
    # Code above to be removed once output data set functionality implemented

    if (self.first_record == None):
      self.first_record = 0  # Take default first record in data set

    if (self.number_records == None):
      self.number_records = input_dataset.num_records  # Process all records

    if (self.rec_comparator == None):
      print 'error:No record comparator defined'
      raise Exception

    if (self.blocking_index == None):
      print 'error:No blocking index defined'
      raise Exception

    if (self.classifier == None):
      print 'error:No classifier defined'
      raise Exception

    if (self.rec_standardiser != None):
      if (self.rec_standardiser.input_dataset != self.input_dataset):
        print 'error:Illegal input data set definition in record '+ \
              'standardiser: %s (should be: %s)' % \
              (str(self.rec_standardiser.input_dataset.name), \
               str(self.input_dataset.name))
        raise Exception
      if (self.rec_standardiser.output_dataset != self.tmp_dataset):
        print 'error:Illegal output data set definition in record '+ \
              'standardiser: %s (should be: %s)' % \
              (str(self.rec_standardiser.output_dataset.name), \
               str(self.tmp_dataset.name))
        raise Exception

    else:  # No standardiser for data set defined, so field names in input
           # and temporary data sets must be the same
      input_field_name_list = self.input_dataset.fields.keys()
      input_field_name_list.sort()
      tmp_field_name_list = self.tmp_dataset.fields.keys()
      tmp_field_name_list.sort()

      if (input_field_name_list != tmp_field_name_list):
        print 'error:Field names differ in input and temporary data sets ' + \
              '(with no record standardiser defined)'

    if (self.blocking_index.dataset != self.tmp_dataset):
      print 'error:Illegal data set definition in blocking index'
      raise Exception

    if (self.rec_comparator.dataset_a != self.tmp_dataset) or \
       (self.rec_comparator.dataset_b != self.tmp_dataset):
      print 'error:Illegal data set definition in record comparator'
      raise Exception

    if (self.classifier.dataset_a != self.tmp_dataset) or \
       (self.classifier.dataset_b != self.tmp_dataset):
      print 'error:Illegal data set definition in classifier'
      raise Exception

    total_time = time.time()  # Get current time

    print '1:'
    print '1:*** Deduplicate data set: %s ***' % (self.input_dataset.name)
    print '1:'

    print '1:'
    print '1:  Step 1: Load and standardise records, build blocking indexes'

    step_1_time = time.time()  # Get current time
    step_1_comm_time = 0.0  # Time for communication in step 1

    input_rec_counter = self.first_record  # Current record pointer

    block_cnt = 0  # A round robin block counter, used for parallelism

    # Load records in a blocked fashion - - - - - - - - - - - - - - - - - - - -

    while (input_rec_counter < (self.first_record + self.number_records)):

      block_size = min(self.block_size,
               ((self.first_record + self.number_records) - input_rec_counter))

      # Distribute blocks equally to all processors
      #
      if ((block_cnt % parallel.size()) == parallel.rank()):

        # Load original records from input data set
        #
        in_recs = self.input_dataset.read_records(input_rec_counter,
                                                  block_size)
        print '1:    Loaded records %i to %i' % \
              (input_rec_counter, input_rec_counter+block_size)

        # Standardise them if a standardiser is defined
        #
        if (self.rec_standardiser != None):
          clean_recs = self.rec_standardiser.standardise_block(in_recs)
          print '1:    Standardised records %i to %i' % \
                (input_rec_counter, input_rec_counter+block_size)
        else:
          clean_recs = in_recs  # Take the original records directly

        # Store records in temporary data set
        #
        self.tmp_dataset.write_records(clean_recs)

        # Insert records into the blocking index
        #
        self.blocking_index.build(clean_recs)

        # If Febrl is run in parallel, send cleaned records to process 0
        #
        if (parallel.rank() > 0):
          tmp_time = time.time()
          parallel.send(clean_recs, 0)
          step_1_comm_time += (time.time() - tmp_time)

      # If Febrl is run in parallel, process 0 receives cleaned records
      #
      if (parallel.rank() == 0) and (block_cnt % parallel.size() != 0):

        p = (block_cnt % parallel.size())  # Process number to receive from
        tmp_time = time.time()
        tmp_recs = parallel.receive(p)
        step_1_comm_time += (time.time() - tmp_time)

        # Store the received records into the temporary data set
        #
        self.tmp_dataset.write_records(tmp_recs)

      input_rec_counter += block_size  # Increment current record pointer
      block_cnt += 1

    # If Febrl is run in parallel, collect blocking index in process 0  - - - -
    #
    if (parallel.rank() == 0):
      for p in range(1, parallel.size()):
        tmp_time = time.time()
        tmp_index = parallel.receive(p)
        step_1_comm_time += (time.time() - tmp_time)
        self.blocking_index.merge(tmp_index)
        print '1:  Received index from process %i' % (p)

    else:
      tmp_time = time.time()
      parallel.send(self.blocking_index, 0) # Send local index to process 0
      step_1_comm_time += (time.time() - tmp_time)
      print '1:  Sent index to process 0'

    # Compact the blocking index on process 0 - - - - - - - - - - - - - - - - -
    #
    if (parallel.rank() == 0):
      self.blocking_index.compact()

    # If run in parallel, broadcast the index from process 0  - - - - - - - - -
    #
    if (parallel.size() > 1):
      if (parallel.rank() == 0):
        for p in range(1, parallel.size()):
          tmp_time = time.time()
          parallel.send(self.blocking_index, p)
          step_1_comm_time += (time.time() - tmp_time)
          print '1:    Sent compacted index to process %i' % (p)

      else:
        tmp_time = time.time()
        self.blocking_index = parallel.receive(0)
        step_1_comm_time += (time.time() - tmp_time)
        print '1:    Received compacted index from process 0'

    # If run in parallel temporary data set needs to be sent to all processes -
    #
    if (parallel.size() > 1):
      if (parallel.rank() == 0):
        for p in range(1, parallel.size()):
          tmp_time = time.time()
          parallel.send(self.tmp_dataset, p)
          step_1_comm_time += (time.time() - tmp_time)
          print '1:    Sent temporary data set %i' % (p)

      else:
        tmp_time = time.time()
        self.tmp_dataset = parallel.receive(0)
        step_1_comm_time += (time.time() - tmp_time)
        print '1:    Received temporary data set from process 0'

    step_1_time = time.time() - step_1_time  # Calculate time for step 1

    #################### START PARALLEL TEST CODE #############################
    # Save temporary data sets and indexes to files (on all processes)
    #
    if (SAVE_PARALLEL_TEST_FILES == True):
      f = open('tmp_data_set-dedup-'+str(parallel.rank())+'-'+ \
               str(parallel.size()),'w')
      tmp_list = self.tmp_dataset.dict.keys()
      tmp_list.sort()

      for r in tmp_list:
        rec = self.tmp_dataset.dict[r]
        rec_items = rec.items()
        rec_items.sort()
        rec = str(r)+': '+str(rec_items)
        f.write(rec+os.linesep)
      f.close()

      f = open('indexes-dedup-'+str(parallel.rank())+'-'+ \
               str(parallel.size()),'w')

      for i in range (self.blocking_index.num_indexes):
        tmp_index = self.blocking_index.index[i].keys()
        tmp_index.sort()

        for bi in tmp_index:
          ind = self.blocking_index.index[i][bi]
          ind_list = ind.items()
          ind_list.sort()

          ii = str(i)+'_'+str(bi)+': '+str(ind_list)
          f.write(ii+os.linesep)
      f.close()

    #################### END PARALLEL TEST CODE ###############################

    # End of step 1 - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    #
    if (parallel.size() > 1):
      parallel.Barrier()  # Make sure all processes are here

    print '1:'
    print '1:  Step 2: Perform deduplication within blocks'

    step_2_time = time.time()  # Get current time
    step_2_comm_time = 0.0  # Time for communication in step 2

    # Get the record pairs which have to be compared  - - - - - - - - - - - - -
    #
    tmp_time = time.time()
    rec_pair_list = indexing.deduplication_rec_pairs(self.blocking_index)

    print '1:    Built record pair list (time used %.2f sec)' % \
          (time.time()-tmp_time)

    rec_pair_list.sort()

#    if (parallel.size() > 1) or (DO_PARALLEL_TEST == True) or \
#       (SAVE_PARALLEL_TEST_FILES == True):
#      rec_pair_list.sort()  # Needed for parallel runs only because of round
#                            # robin fashion of work distribution (and list can
#                            # have different sequence on different processes

    #################### START PARALLEL TEST CODE #############################
    # Check if the rec_pair_list are the same on all processes
    #
    if (parallel.size() > 1) and (DO_PARALLEL_TEST == True):
      if (parallel.rank() == 0):

        for p in range(1, parallel.size()):
          tmp_list = parallel.receive(p)
          if (len(tmp_list) != len(rec_pair_list)):
            print 'warning:Record pair lists have differnt length on ' + \
                  'process 0 (%i) and %i (%i)' % \
                  (len(rec_pair_list), p, len(tmp_list))

          for i in range(len(rec_pair_list)):
            if (rec_pair_list[i] != tmp_list[i]):  # Element differs
              print 'warning:Record pair lists differ in position %i' % (i) + \
                    ' on process 0 (%s) and %i (%s)' % \
                    (p, rec_pair_list[i], tmp_list[i])
      else:
        parallel.send(rec_pair_list,0)

      parallel.Barrier()

    # Save record pair lists to files (all processes)
    #
    if (SAVE_PARALLEL_TEST_FILES == True):
      f = open('rec_pair_list-dedup-'+str(parallel.rank())+'-'+ \
               str(parallel.size()),'w')
      for r in rec_pair_list:
        f.write(r+os.linesep)
      f.close()

    #################### END PARALLEL TEST CODE ###############################

    weight_vector_dict = {}  # Dictionary with the comparison vectors

    compare_time =  0.0
    classify_time = 0.0

    num_rec_pairs = len(rec_pair_list)
    rec_pair_cnt = 0  # Loop counter

    # Compare records, and distribute comparisons equally to all processes  - -
    #
    for rec_pair in rec_pair_list:

      if ((rec_pair_cnt % parallel.size()) == parallel.rank()):

        [rec_num_a,rec_num_b] = rec_pair.split('_')
        print '2:      Compare records %s with %s' % (rec_num_a, rec_num_b)

        # Read the records from the data set
        #
        rec_a = self.tmp_dataset.read_record(int(rec_num_a))
        rec_b = self.tmp_dataset.read_record(int(rec_num_b))

        #print '1: -------------------------------'  ########
        #print '1: record_a: %s' % (str(rec_a))      ########
        #print '1: record_b: %s' % (str(rec_b))      ########

        # Compare the two records
        #
        tmp_time = time.time()
        w_vector = self.rec_comparator.compare(rec_a, rec_b)
        compare_time += (time.time() - tmp_time)

        #print '1: weight vector: %s' % (str(w_vector))   ######
        #print '1: -------------------------------'       ######

        # Save the weight vector in a dictionary (used later for output)
        #
        if (weight_vector_dict.has_key(rec_pair)):
          print 'warning:This should never happen: Record pair %s ' % \
                (rec_pair) + ' already in weight vector dictionary'
        weight_vector_dict[rec_pair] = w_vector

        tmp_time = time.time()
        self.classifier.classify(w_vector)  # Classify the weight vector
        classify_time += (time.time() - tmp_time)

      rec_pair_cnt += 1

      # Progress report every 10 per cent
      #
      if (rec_pair_cnt > 0) and ((rec_pair_cnt % int(num_rec_pairs/10)) == 0):
        print '1:    %i/%i record pairs compared' % \
              (rec_pair_cnt, num_rec_pairs)
        print '1:      Average comparison time:     %.6f' % \
              (compare_time/rec_pair_cnt)
        print '1:      Average classification time: %.6f' % \
              (classify_time/rec_pair_cnt)

    # Now gather classifiers on process 0 and merge - - - - - - - - - - - - - -
    #
    if (parallel.size() > 1):
      if (parallel.rank() == 0):
        for p in range(1, parallel.size()):
          tmp_time = time.time()
          tmp_classifier = parallel.receive(p)
          step_2_comm_time += (time.time() - tmp_time)
          self.classifier.merge(tmp_classifier)
          print '1:  Received classifier from process %i and merged it' % (p)

      else:
        tmp_time = time.time()
        parallel.send(self.classifier, 0) # Send local classifier to process 0
        step_2_comm_time += (time.time() - tmp_time)
        print '1:  Sent classifier to process 0'

    # Also gather weight vectors on process 0 and merge - - - - - - - - - - - -
    #
    if (parallel.size() > 1):
      if (parallel.rank() == 0):
        for p in range(1, parallel.size()):
          tmp_time = time.time()
          tmp_weight_vector_dict = parallel.receive(p)
          step_2_comm_time += (time.time() - tmp_time)
          weight_vector_dict.update(tmp_weight_vector_dict)
          print '1:  Received weight vectors from process %i' % (p)

      else:
        tmp_time = time.time()
        parallel.send(weight_vector_dict, 0) # Send local weight vectors
        step_2_comm_time += (time.time() - tmp_time)
        print '1:  Sent weight vectors to process 0'

    #################### START PARALLEL TEST CODE #############################
    # Save classifiers and weight vectors to files (only process 0)
    #
    if (SAVE_PARALLEL_TEST_FILES == True) and (parallel.rank() == 0):
      tmp_list = self.classifier.results.keys()
      tmp_list.sort()
      f = open('classifier_results_dict-dedup-'+str(parallel.rank())+'-'+ \
               str(parallel.size()),'w')
      for c in tmp_list:
        res = self.classifier.results[c].items()
        res.sort()
        ce = str(c)+': '+str(res)
        f.write(ce+os.linesep)
      f.close()

      tmp_list = weight_vector_dict.keys()
      tmp_list.sort()
      f = open('weight_vector_dict-dedup-'+str(parallel.rank())+'-'+ \
               str(parallel.size()),'w')
      for v in tmp_list:
        wv = v+': '+str(weight_vector_dict[v])
        f.write(wv+os.linesep)
      f.close()

    #################### END PARALLEL TEST CODE ###############################

    print '1:Deduplication done, totally %i comparisons' % (num_rec_pairs)

    step_2_time = time.time() - step_2_time  # Calculate time for step 2

    # Output the results  - - - - - - - - - - - - - - - - - - - - - - - - - - -

    if (parallel.rank() == 0):  # Only processor 0 prints results

      # Get the results dictionary with all the record pairs and their weights
      #
      results_dict = self.classifier.results

      if (self.output_histogram == True):  # Print a weights histogram

        output.print_histogram(results_dict)

      if (self.output_assignment != None):  # An output assignment is defined

        if (self.output_assignment == 'one2one'):

          # Do a one-to-one assignment on the classifier results dict
          #
          o2o_results_dict = lap.do_lap('lapmod', results_dict, \
                                        'deduplication', self.output_threshold)
      else:  # No one-to-one assignment, set o2o result to None
        o2o_results_dict = None

      if (self.output_print == True):  # Print resulting record pairs

        output.print_record_pairs(self.tmp_dataset, self.tmp_dataset, \
                                  results_dict, o2o_results_dict, \
                                  self.output_threshold)

      if (self.output_file != None):  # Save results into a file

        output.save(self.tmp_dataset.name, self.tmp_dataset.name, \
                    self.output_file, results_dict, o2o_results_dict, \
                    self.output_threshold)

    print '1:Done.'

    total_time = time.time() - total_time  # Calculate total time

    parallel.Barrier()  # Wait here for all processes - - - - - - - - - - - - -

    print '1:Total time needed for deduplication of %i records: %.3f' % \
          (self.number_records, total_time)

    print '1:  Time for step 1: %.3f' % (step_1_time)
    print '1:  Time for step 2: %.3f' % (step_2_time)
    print '1:  Time for communication in step 1: %.3f' % (step_1_comm_time)
    print '1:  Time for communication in step 2: %.3f' % (step_2_comm_time)

    parallel.Finalize()

# =============================================================================

class ProjectLog:
  """Class for Febrl project logs.

     Project logs capture print statements with a special form as well as
     Python exceptions and warnings, and writes them to a log file and prints
     them to the standard terminal output.

     A log file name can be given with the attribute 'file_name' when a project
     log is initialised. If no file name is given, no logging information is
     written to a file (but still printed to standard output).

     When a project log is specified, the user needs to set both a verbose and
     a log level in the range 0 to 3 using the attributes 'log_level' and
     'verbose_level'. A level of 0 means nothing no log messages are logged or
     printed to standard output, level one means only important messages are
     logged/printed, level 2 means medium output and level 3 finally means a
     high volume output.

     Error messages are looged and printed in any case.

     Logging and printing of warning massages can be suppressed by setting the
     flag 'no_warnings' to 'True' when the project log file is initialised.
     The default is 'False' meaning that warning messages are printed.

     With the flag 'clear_log' the user can clear the contents of a log file
     when it is initialised by setting the value of this flag to 'True'. The
     default is 'False', in which case logging messages are appended to the
     given log file.

     Within Febrl modules, log messages are normal Python statements that must
     start with either of the following substrings followed by the message:

       '1:'        (a high priority log message)
       '2:'        (a medium priority log message)
       '3:'        (a low priority log message)
       'warning:'  (a warning message)
       'error:'    (an error message)

     All other print statements are handled like normal prints and sent to
     standard output 'sys.stdout'.

     The 'parallel_print' defines how printing and logging is done when Febrl
     is run in parallel. Possible are either 'host' (the default) in which case
     only the host node (where Febrl was started on) is printing and logging
     messages, or 'all' in which case all processes are printing and logging.
     Note that currently both error and warning messages are printed by all
     processes.
     The 'parallel_print' is used to set 'printmode' in the parallel.py module.

     The 'message' is the string that is printed to standard output and written
     to the log file if 'level' is equal to smaller than the 'verbose_level' or
     the 'log_level' respectively.
  """

  # ---------------------------------------------------------------------------

  def __init__(self, **kwargs):
    """Constructor - Create a new project log object.
    """

    self.febrl =          None    # A reference to the febrl object.
    self.project =        None    # A reference to the project object.
    self.log_level =      0       # Level of looging to file.
    self.verbose_level =  0       # Level of verbose printing.
    self.no_warnings =    False   # Supress warning messages
    self.clear_log =      False   # A flag (True/False), if set to True all
                                  # content of the log will be cleared first.
                                  # Default value is False, i.e. messages will
                                  # be appended to the existing log file.
    self.file_name =      None    # The log file name.
    self.log_file =       None    # The log file pointer
    self.parallel_print = 'host'  # Either set to 'host' (default) or 'all'

    # Process all keyword arguments
    #
    for (keyword, value) in kwargs.items():
      if (keyword == 'project'):
        self.project = value
        if (self.project.febrl != None):
          self.febrl = self.project.febrl  # Reference to the febrl object
        else:
          print 'error:Febrl object not defined'
          raise Exception

      elif (keyword == 'file_name'):
        if (not isinstance(value,str)):
          print 'error:Argument "file_name" is not a string'
          raise Exception
        self.file_name = value

      elif (keyword == 'clear_log'):
        if (value not in [True, False]):
          print 'error:Argument "clear_log" must be "True" or "False"'
          raise Exception
        self.clear_log = value

      elif (keyword in ['no_warnings','no_warn']):
        if (value not in [True, False]):
          print 'error:Argument "no_warnings" must be "True" or "False"'
          raise Exception
        self.no_warnings = value

      elif (keyword == 'parallel_print'):
        if (value not in ['host','all']):
          print 'error:Argument "parallel_print" must be set to "host"'+ \
                ' or "all"'
          raise Exception
        self.parallel_print = value

      elif (keyword == 'log_level'):
        if (not isinstance(value,int)) or (value < 0) or (value > 3):
          print 'error:Argument "log_level" must zero or a positive integer'
          raise Exception
        self.log_level = value

      elif (keyword == 'verbose_level'):
        if (not isinstance(value,int)) or (value < 0) or (value > 3):
          print 'error:Argument "verbose_level" must zero or a positive '+ \
                'integer'
          raise Exception
        self.verbose_level = value

      else:
        print 'error:Illegal constructor argument keyword: "%s"' % \
              (str(keyword))
        raise Exception

    # Check if the needed attributes are set  - - - - - - - - - - - - - - - - -
    #
    if (self.project == None):
      print 'error:Project not defined'
      raise Exception

    # Set the parallel printing mode  - - - - - - - - - - - - - - - - - - - - -
    #
    parallel.printmode = self.parallel_print

    # Check if file can be opened - - - - - - - - - - - - - - - - - - - - - - -
    #
    if (self.file_name != None):  # File logging activated

      if (self.clear_log == True):  # Open log file in 'write' mode
        try:
          self.log_file = open(self.file_name,'w')
        except:
          print 'error:Can not open log file "%s" for writing' % \
                (self.file_name)
          raise IOError

      else:  # Open log file in 'append' mode (don't delete old content)
        try:
          self.log_file = open(self.file_name,'a')
        except:
          print 'error:Can not open log file "%s" for appending' % \
                (self.file_name)
          raise IOError

      # Write a header to the log file - - - - - - - - - - - - - - - - - - - -
      #
      self.append('#'*75)
      self.append('# Febrl project log file')
      self.append('#')
      self.append('# Project name: '+self.project.name)
      self.append('# Project log file name: '+self.file_name)
      self.append('#')
      self.append('# Time and date: '+time.strftime('%d %b %Y %H:%M:%S', \
                                      time.localtime(time.time())))
      self.append('#')
      self.append('# Febrl version:   '+self.febrl.version)
      self.append('# Febrl license:   '+self.febrl.license)
      self.append('# Febrl copyright: '+self.febrl.copyright)
      self.append('#')
      self.append('# Febrl description: '+self.febrl.description)
      self.append('# Febrl path:        '+self.febrl.febrl_path)
      self.append('#'*75)
      self.append('')
      self.flush()

      # Redirect the system stdout so prints can be logged - - - - - - - - - -
      #
      sys.stdout = LogPrinter(self)

      # Set the system exception hook so exceptions and warnings are logged - -
      #
      sys.excepthook = self.except_hook

  # ---------------------------------------------------------------------------

  def except_hook(self, exc_type, value, trace_back):
    """A routine to catch exceptions and warnings and process them.
    """

    bug_report1 = 'Please submit an error report by sending an e-mail'+ \
                  ' to the Febrl authors'
    bug_report2 = 'and attach this error message.'

    time_stamp = time.strftime('%d %b %Y %H:%M:%S',time.localtime(time.time()))

    # Get complete trace stack
    #
    trace_list = traceback.extract_tb(trace_back)
    trace_stack_size = len(trace_list)

    self.append(parallel.prompt+'#'*75)
    self.append(parallel.prompt+'### Exception: '+str(exc_type))
    self.append(parallel.prompt+'###   Time:       '+time_stamp)
    self.append(parallel.prompt+'###   Message:    '+str(value))
    self.append(parallel.prompt+'###   Trace stack:')
    for lev in range(trace_stack_size):
      spc = '  '*lev
      self.append(parallel.prompt+'###     '+spc+'-'*(67-lev*2))
      self.append(parallel.prompt+'###     '+spc+'Module:   '+ \
                  str(trace_list[lev][0]))
      self.append(parallel.prompt+'###     '+spc+'Function: '+ \
                  str(trace_list[lev][2]))
      self.append(parallel.prompt+'###     '+spc+'Line:     '+ \
                  str(trace_list[lev][1]))
      self.append(parallel.prompt+'###     '+spc+'Text:     '+ \
                  str(trace_list[lev][3]))
    self.append(parallel.prompt+'#'*75)
    self.append(parallel.prompt+'### '+bug_report1)
    self.append(parallel.prompt+'### '+bug_report2)
    self.append(parallel.prompt+'#'*75)
    self.flush()

    sys.__stdout__.write(parallel.prompt+'#'*75+os.linesep)
    sys.__stdout__.write(parallel.prompt+'### Exception: '+str(exc_type)+ \
                         os.linesep)
    sys.__stdout__.write(parallel.prompt+'###   Time:       '+time_stamp+ \
                         os.linesep)
    sys.__stdout__.write(parallel.prompt+'###   Message:    '+str(value)+ \
                         os.linesep)
    sys.__stdout__.write(parallel.prompt+'###   Trace stack:'+os.linesep)
    for lev in range(trace_stack_size):
      spc = '  '*lev
      sys.__stdout__.write(parallel.prompt+'###     '+spc+'-'*(67-lev*2)+ \
                           os.linesep)
      sys.__stdout__.write(parallel.prompt+'###     '+spc+'Module:   '+ \
                           str(trace_list[lev][0])+os.linesep)
      sys.__stdout__.write(parallel.prompt+'###     '+spc+'Function: '+ \
                           str(trace_list[lev][2])+os.linesep)
      sys.__stdout__.write(parallel.prompt+'###     '+spc+'Line:     '+ \
                           str(trace_list[lev][1])+os.linesep)
      sys.__stdout__.write(parallel.prompt+'###     '+spc+'Text:     '+ \
                           str(trace_list[lev][3])+os.linesep)
    sys.__stdout__.write(parallel.prompt+'#'*75+os.linesep)
    sys.__stdout__.write(parallel.prompt+'### '+bug_report1+os.linesep)
    sys.__stdout__.write(parallel.prompt+'### '+bug_report2+os.linesep)
    sys.__stdout__.write(parallel.prompt+'#'*75+os.linesep)

  # ---------------------------------------------------------------------------

  def close(self):
    """Close the log file if it was opened.
    """

    if (self.log_file != None):
      self.log_file.close()
      self.log_file = None

  # ---------------------------------------------------------------------------

  def append(self, message):
    """Write the given message to the log file if logging is activated (i.e. if
       the log file is opened).

       The message can be a string or a list of strings, in which case each
       list element will be written as one line. A line separator is appended
       to the end of each line.
    """

    if (self.log_file != None):
      if (isinstance(message, str)):
        if (len(message) > 0) and (message[-1] == os.linesep):
          self.log_file.write(message)
        else:
          self.log_file.write(message+os.linesep)

      elif (isinstance(message, list)):
        for m in message:
          if (isinstance(m, str)):
            if (len(m) > 0) and (message[-1] == os.linesep):
              self.log_file.write(m)
            else:
              self.log_file.write(m+os.linesep)

          else:
            print 'error:Element in log file message list is not a ' + \
                  'string: %s' % (str(m))
            raise Exception
      else:
        print 'error:Log file message is not a string or a list: %s' % \
              (str(message))
        raise Exception

  # ---------------------------------------------------------------------------

  def flush(self):
    """Flush the log file if it is opened to make sure it is written out.
    """

    if (self.log_file != None):
      self.log_file.flush()

  # ---------------------------------------------------------------------------

  def print_log(self, print_type):
    """Print the contents of the log file in various types:
       Text, HTML, LaTex, XML
    """

    pass

    # Check if file is open

    # re-open in read mode

    # print all lines in format

    # close in read mode

# =============================================================================

class LogPrinter:
  """Class that replaces sys.stdout.

     Each normal print statement is analysed in the 'write' method, if it
     starts with a 'Error:', 'Warning:' or '1:', '2:', '3:', .. '9:' it will
     also be passed to the project loger.
     Otherwise, it will simple be given to the original standard out.

     For parallel runs, error and warning message will be printed (and logged)
     from all processes, but normal verbose message are printed acording to the
     value of 'printmode' ('all' or 'host') as defined in parallel.py
  """

  def __init__(self, project_log):
    self.project_log = project_log

  def write(self, msg): # - - - - - - - - - - - - - - - - - - - - - - - - - -

    #if (len(msg) > 0) and (msg[-1] != os.linesep):  # Append a line separator
    #  msg += os.linesep

    call_function =   sys._getframe(1).f_code.co_name
    call_linenumber = str(sys._getframe(1).f_lineno)
    call_module =     sys._getframe(1).f_code.co_filename
    time_stamp = time.strftime('%d %b %Y %H:%M:%S',time.localtime(time.time()))

    # Get message type and split message at line separators
    #
    if (msg[:6].lower() == 'error:'):  # An error message
      msg_type = 'error'
      msg = msg[6:]

    elif (msg[:8].lower() == 'warning:'):  # A warning message
      msg_type = 'warn'
      msg = msg[8:]

    elif (msg[0].isdigit()) and (msg[1] == ':'):  # A verbose level message
      msg_type = 'level'
      msg_level = int(msg[0])
      msg = msg[2:]

    else:  # Other message
      msg_type = 'other'

    msg_list = msg.split(os.linesep)

    if (msg_list == ['','']):  # An empty message, don't print or log it
      return

    if (msg_type == 'error'):   # - - - - - - - - - - - - - - - - - - - - - - -

      bug_report1 = 'Please submit an error report by sending an e-mail'+ \
                    ' to the Febrl authors'
      bug_report2 = 'and attach this error message.'

      self.project_log.append(parallel.prompt+'#'*75)
      self.project_log.append(parallel.prompt+'### Error')
      self.project_log.append(parallel.prompt+'###   Module:  '+ \
                              call_module+', function: '+ call_function+ \
                              ', line number: '+call_linenumber)
      self.project_log.append(parallel.prompt+'###   Time:    '+time_stamp)
      self.project_log.append(parallel.prompt+'###   Message: '+msg_list[0])
      for m in msg_list[1:]:
        self.project_log.append(parallel.prompt+'###            '+m)
      self.project_log.append(parallel.prompt+'#'*75)
      self.project_log.append(parallel.prompt+'### '+bug_report1)
      self.project_log.append(parallel.prompt+'### '+bug_report2)
      self.project_log.append(parallel.prompt+'#'*75)
      self.project_log.flush()

      sys.__stdout__.write(parallel.prompt+'#'*75+os.linesep)
      sys.__stdout__.write(parallel.prompt+'### Error'+os.linesep)
      sys.__stdout__.write(parallel.prompt+'###   Module:  '+call_module+ \
                           ', function: '+call_function+', line number: '+ \
                           call_linenumber+os.linesep)
      sys.__stdout__.write(parallel.prompt+'###   Time:    '+time_stamp+ \
                           os.linesep)
      sys.__stdout__.write(parallel.prompt+'###   Message: '+msg_list[0]+ \
                           os.linesep)
      for m in msg_list[1:]:
        self.project_log.append(parallel.prompt+'###            '+m+os.linesep)
      sys.__stdout__.write(parallel.prompt+'#'*75+os.linesep)
      sys.__stdout__.write(parallel.prompt+'### '+bug_report1+os.linesep)
      sys.__stdout__.write(parallel.prompt+'### '+bug_report2+os.linesep)
      sys.__stdout__.write(parallel.prompt+'#'*75+os.linesep)

    elif (msg_type == 'warn'):
      if (self.project_log.no_warnings == False):
        self.project_log.append(parallel.prompt)
        self.project_log.append(parallel.prompt+'### Warning')
        self.project_log.append(parallel.prompt+'###   Module:  '+ \
                                call_module+', function: '+ \
                                call_function+', line number: '+ \
                                call_linenumber)
        self.project_log.append(parallel.prompt+'###   Time:    '+time_stamp)
        self.project_log.append(parallel.prompt+'###   Message: '+msg_list[0])
        for m in msg_list[1:]:
          self.project_log.append(parallel.prompt+'###            '+m)
        self.project_log.flush()

        sys.__stdout__.write(parallel.prompt+os.linesep)
        sys.__stdout__.write(parallel.prompt+'### Warning'+os.linesep)
        sys.__stdout__.write(parallel.prompt+'###   Module:  '+call_module+ \
                             ', function: '+call_function+', line number: '+ \
                             call_linenumber+os.linesep)
        sys.__stdout__.write(parallel.prompt+'###   Time:    '+time_stamp+ \
                             os.linesep)
        sys.__stdout__.write(parallel.prompt+'###   Message: '+msg_list[0]+ \
                             os.linesep)
        for m in msg_list[1:]:
          self.project_log.append(parallel.prompt+'###            '+m+ \
                                  os.linesep)
      else:
        pass  # Don't print if 'no warn' flag is set

    elif (msg_type == 'level') and (msg_level in [1,2,3]):  # A 'level' message

      # Check if Febrl runs in parallel and adjust printing accordingly
      #
      if ((parallel.printmode == 'host') and (parallel.rank() == 0)) or \
          (parallel.printmode == 'all'):

        if (msg_level <= self.project_log.log_level):
          for m in msg_list:
            self.project_log.append(parallel.prompt+m)
          self.project_log.flush()

        if (msg_level <= self.project_log.verbose_level):
          for m in msg_list:
            sys.__stdout__.write(parallel.prompt+m+os.linesep)

    else:  # Print 'normal' print commands
      if (msg_list != [os.linesep]) and (msg_list != ['']): 

        if ((parallel.printmode == 'host') and (parallel.rank() == 0)) or \
            (parallel.printmode == 'all'):
          for m in msg_list:
            sys.__stdout__.write(parallel.prompt+m+os.linesep)

# =============================================================================
