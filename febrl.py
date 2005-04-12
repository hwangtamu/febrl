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
# The Original Software is: "febrl.py"
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

"""Module febrl.py - Main module and classes for Febrl projects.

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

import copy
import cPickle
import logging
import os
import sys
import time
import traceback
import types

import parallel
import indexing
import output
import lap

# =============================================================================

dirsep = os.sep  # Shorthand to the directory/folder separator character
                 # (mainly '/', ':', or '\')

# =============================================================================

class Febrl:
  """Class Febrl - Main class for Febrl projects.
  """

  def __init__(self, **kwargs):
    """Constructor - Set attributes and load list of available projects.
    """

    linesep = os.linesep

    self.version_major =      '0.3'
    self.version_minor =      ''
    self.version =            self.version_major+'.'+self.version_minor
    self.license =            'ANUOS Version 1.2'
    self.copyright =          '(C) 2002 - 2005 the Australian National ' + \
                              'University and others.'
    self.initial_developers = linesep + '    Dr Peter Christen (Department' + \
                              ' of Computer Science, Australian National' + \
                              linesep + \
                              '                       University) and' + \
                              linesep + \
                              '    Dr Tim Churches   (Centre for ' + \
                              'Epidemiology and Research, New South Wales ' + \
                              linesep + '                       Department' + \
                              ' of Health)'
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
        logging.exception('Illegal constructor argument keyword: "%s"' % \
                          (str(keyword)))
        raise Exception

    # Check if Febrl projects are available in the 'project_path' directory by
    # scanning the directory for files with '.fbr' file extension
    #
    file_list = os.listdir(self.febrl_path)
    for fn in file_list:
      file_name = fn.strip().lower()
      if (file_name[-4:] == '.fbr'):
        self.project_names.append(file_name)

    logging.info(self.__str__())  # Print header with license and copyright

  # ---------------------------------------------------------------------------

  def __str__(self):
    """Create a string representation of the Febrl object.
    """

    linesep = os.linesep

    rep = linesep  # Start with an empty line
    rep += '============================================================='
    rep += '====================' + linesep
    rep += linesep
    rep += 'Febrl (Freely extensible biomedical record linkage)' + linesep
    rep += '---------------------------------------------------' + linesep
    rep += linesep
    rep += '  Version:   ' + self.version + linesep
    rep += '  License:   ' + self.license + linesep
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
    rep += linesep
    rep += '============================================================='
    rep += '====================' + linesep

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
        logging.exception('Illegal project number: %s' % (str(project)))
        raise Exception

    elif (type(project) == types.StringType):
      file_name += project
    else:
      logging.exception('Illegal type for "project" argument, must be ' + \
                        'either of type string or integer')
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

  # ---------------------------------------------------------------------------

  def finalise(self):
    """Finalise a Febrl project.
    """

    logging.info('Febrl stopped.')

    parallel.Barrier()  # Make sure all processes are here

    parallel.Finalize()  # Finalise parallel environment

    sys.exit()

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
    self.parallel_write =   'host' # Set either to 'host' (default) or 'all'

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
        check_argument_is_integer(keyword, value)
        if (value == 0):
          logging.exception('Argument "block_size" is not a positive integer')
          raise Exception
        self.block_size = value

      elif (keyword == 'parallel_write'):
        if (value not in ['host','all']):
          logging.exception('Argument "parallel_write" must be set to ' + \
                            '"host" or "all"')
          raise Exception
        self.parallel_write = value

      else:
        logging.exception('Illegal constructor argument keyword: "%s"' % \
                          (str(keyword)))
        raise Exception

    # Set the parallel writing/saving mode  - - - - - - - - - - - - - - - - - -
    #
    parallel.writemode = self.parallel_write

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

  def standardise(self, **kwargs):
    """Clean and standardise the given data set using the defined record
       standardiser.

       Records are loaded block wise from the input data set, then standardised
       and written into the output data set.

       If the argument 'first_record' is not given, it will automatically be
       set to the first record in the data set (i.e. record number 0).
       Similarly, if the argument 'number_records' is not given, it will be set
       to the total number of records in the input data set.

       The output data set can be any data set type except a memory based data
       set (as all standardised records would lost once the program finishes).
       This output data set has to be initialised in 'write', 'append' or
       'readwrite' access mode.
    """

    self.input_dataset =    None   # A reference ot the (raw) input data set
                                   # data set
    self.output_dataset =   None   # A reference to the output data set

    self.rec_standardiser = None   # Reference to a record standardiser
    self.first_record =     None   # Number of the first record to process
    self.number_records =   None   # Number of records to process

    for (keyword, value) in kwargs.items():
      if (keyword == 'input_dataset'):
        self.input_dataset = value
      elif (keyword == 'output_dataset'):
        self.output_dataset = value

      elif (keyword == 'rec_standardiser'):
        self.rec_standardiser = value

      elif (keyword == 'first_record'):
        if (not isinstance(value, int)) or (value < 0):
          logging.exception('Argument "first_record" is not a valid ' + \
                            'integer number')
          raise Exception
        self.first_record = value
      elif (keyword == 'number_records'):
        if (not isinstance(value, int)) or (value <= 0):
          logging.exception('Argument "number_records" is not a positive' + \
                            ' integer number')
          raise Exception
        self.number_records = value

      else:
        logging.exception('Illegal constructor argument keyword: "%s"' % \
                          (str(keyword)))
        raise Exception

    # Do some checks on the input arguments - - - - - - - - - - - - - - - - - -
    #
    if (self.input_dataset == None):
      logging.exception('Input data set is not defined')
      raise Exception

    if (self.output_dataset == None):
      logging.exception('Output data set is not defined')
      raise Exception
    elif (self.output_dataset.dataset_type == 'MEMORY'):
      logging.exception('Output data set can not be a memory based data set')
      raise Exception
    if (self.output_dataset.access_mode not in ['write','append','readwrite']):
      logging.exception('Output dataset must be initialised in one of the' + \
                        ' access modes: "write", "append", or "readwrite"')
      raise Exception

    if (self.rec_standardiser == None):
      logging.exception('Record standardiser is not defined')
      raise Exception

    if (self.first_record == None):
      self.first_record = 0  # Take default first record in data set

    if (self.number_records == None):  # Process all records
      self.number_records = self.input_dataset.num_records

    logging.info('***** Standardise data set: "%s" (type "%s")' % \
                 (self.input_dataset.name, self.input_dataset.dataset_type))
    logging.info('*****        into data set: "%s" (type "%s")' % \
                 (self.output_dataset.name, self.output_dataset.dataset_type))

    # Call the main standardisation routine - - - - - - - - - - - - - - - - - -
    # (no indexing will be done)
    #
    [stand_time, comm_time] = do_load_standard_indexing(self.input_dataset,
                                                        self.output_dataset,
                                                        self.rec_standardiser,
                                                        None,
                                                        self.first_record,
                                                        self.number_records,
                                                        self.block_size)

  # ---------------------------------------------------------------------------

  def geocode(self, **kwargs):
    """Geocode a data set using referenced geocoded data. If a record
       standardiser is given, the input data set if first cleaned and
       standardised, otherwise the input data set is geocoded directly without
       cleaning and standardisation.

       Records are loaded block wise from the input data set, and then geocoded
       and written to the output data set.

       ***** this below will change *****
       For each record in the input data set, none, one or many matched
       geocoded reference records are written to the output data set. If more
       than one matching record are found, they are sorted with highest match
       weight first.
       *****

       If the argument 'first_record' is not given, it will automatically be
       set to the first record in the data set (i.e. record number 0).
       Similarly, if the argument 'number_records' is not given, it will be set
       to the total number of records in the input data set.

       The output data set can be any data set type except a memory based data
       set (as all standardised records would lost once the program finishes).
       This output data set has to be initialised in 'write', 'append' or
       'readwrite' access mode.
    """

    self.input_dataset =    None   # A reference ot the (raw) input data set
                                   # data set
    self.output_dataset =   None   # A reference to the output data set

    self.rec_standardiser = None   # Reference to a record standardiser
    self.rec_geocoder =     None   # Refernce to a record geocoder
    self.first_record =     None   # Number of the first record to process
    self.number_records =   None   # Number of records to process

    for (keyword, value) in kwargs.items():
      if (keyword == 'input_dataset'):
        self.input_dataset = value
      elif (keyword == 'output_dataset'):
        self.output_dataset = value

      elif (keyword == 'rec_standardiser'):
        self.rec_standardiser = value

      elif (keyword == 'rec_geocoder'):
        self.rec_geocoder = value

      elif (keyword == 'first_record'):
        if (not isinstance(value, int)) or (value < 0):
          logging.exception('Argument "first_record" is not a valid ' + \
                            'integer number')
          raise Exception
        self.first_record = value
      elif (keyword == 'number_records'):
        if (not isinstance(value, int)) or (value <= 0):
          logging.exception('Argument "number_records" is not a positive' + \
                            ' integer number')
          raise Exception
        self.number_records = value

      else:
        logging.exception('Illegal constructor argument keyword: "%s"' % \
                          (str(keyword)))
        raise Exception

    # Do some checks on the input arguments - - - - - - - - - - - - - - - - - -
    #
    if (self.input_dataset == None):
      logging.exception('Input data set is not defined')
      raise Exception

    if (self.output_dataset == None):
      logging.exception('Output data set is not defined')
      raise Exception
    elif (self.output_dataset.dataset_type == 'MEMORY'):
      logging.exception('Output data set can not be a memory based data set')
      raise Exception
    if (self.output_dataset.access_mode not in ['write','append','readwrite']):
      logging.exception('Output dataset must be initialised in one of the' + \
                        ' access modes: "write", "append", or "readwrite"')
      raise Exception

    if (self.rec_geocoder == None):
      logging.exception('Record geocoder is not defined')
      raise Exception

    if (self.first_record == None):
      self.first_record = 0  # Take default first record in data set

    if (self.number_records == None):  # Process all records
      self.number_records = self.input_dataset.num_records

    if (self.rec_standardiser == None):
      logging.info('***** Geocode data set: "%s" (type "%s")' % \
                   (self.input_dataset.name, self.input_dataset.dataset_type))
      logging.info('*****    into data set: "%s" (type "%s")' % \
                   (self.output_dataset.name,self.output_dataset.dataset_type))
    else:
      logging.info('***** Standardise and geocode data set: "%s" (type "%s")' \
                   % (self.input_dataset.name,self.input_dataset.dataset_type))
      logging.info('*****                    into data set: "%s" (type "%s")' \
                   % (self.output_dataset.name, \
                   self.output_dataset.dataset_type))

    # Call the main standardisation and geocoding routine - - - - - - - - - - -
    # (no indexing will be done)
    #
    [stand_time, comm_time] = do_load_standard_geocode(self.input_dataset,
                                                       self.output_dataset,
                                                       self.rec_standardiser,
                                                       self.rec_geocoder,
                                                       self.first_record,
                                                       self.number_records,
                                                       self.block_size)

  # ---------------------------------------------------------------------------

  def geocode_server_init(self, **kwargs):
    """Initialise a geocode server engine.
    """

    self.input_dataset =          None  # A reference ot the (raw) input data
                                        # set data set
    self.output_dataset =         None  # A reference to the output data set

    self.rec_standardiser =       None  # Reference to a record standardiser
    self.rec_geocoder =           None  # Reference to a record geocoder
    self.input_address_log_file = None  # A file for logging input addresses

    for (keyword, value) in kwargs.items():
      if (keyword == 'input_dataset'):
        self.input_dataset = value
      elif (keyword == 'output_dataset'):
        self.output_dataset = value

      elif (keyword == 'rec_standardiser'):
        self.rec_standardiser = value

      elif (keyword == 'rec_geocoder'):
        self.rec_geocoder = value

      elif (keyword == 'input_address_log_file'):
        self.input_address_log_file = value

      else:
        logging.exception('Illegal constructor argument keyword: "%s"' % \
                          (str(keyword)))
        raise Exception

    # Do some checks on the input arguments - - - - - - - - - - - - - - - - - -
    #
    if (self.input_dataset.dataset_type != 'MEMORY'):
      logging.exception('Input data set is not defined as memory data set')
      raise Exception

    if (self.output_dataset.dataset_type != 'MEMORY'):
      logging.exception('Output data set is not defined as memory data set')
      raise Exception

    if (self.input_dataset.access_mode != 'readwrite'):
      logging.exception('Input data set must be initialised in access' + \
                        ' mode: "readwrite"')
      raise Exception

    if (self.output_dataset.access_mode != 'readwrite'):
      logging.exception('output data set must be initialised in access' + \
                        ' mode: "readwrite"')
      raise Exception

    if (self.rec_geocoder == None):
      logging.exception('Record geocoder is not defined')
      raise Exception

    if (self.rec_standardiser == None):
      logging.exception('Record standardiser is not defined')
      raise Exception

    if (self.rec_standardiser == None):
      logging.info('***** Geocode data set: "%s" (type "%s")' % \
                   (self.input_dataset.name, self.input_dataset.dataset_type))
      logging.info('*****    into data set: "%s" (type "%s")' % \
                   (self.output_dataset.name, \
                   self.output_dataset.dataset_type))
    else:
      logging.info('***** Standardise and geocode data set: "%s" (type "%s")' \
                   % (self.input_dataset.name, \
                   self.input_dataset.dataset_type))
      logging.info('*****                    into data set: "%s" (type "%s")' \
                   % (self.output_dataset.name, \
                   self.output_dataset.dataset_type))

    if (self.input_address_log_file != None):
      logging.info('  Save received input addresses into log file "%s"' % \
                   (self.input_address_log_file))

    # -------------------------------------------------------------------------
    # The request handler function which does the standardisation
    # and geocoding.
    # The input arguments are:
    # - message_type ('text' or 'file')
    # - input_data   (non-empty string containing one or more addresses,
    #                 separated by '\n')
    # - neighbour_level (integer 0, 1, or 2)
    # - result_report   (integer 0, 1, or 2)

    def do_geocoding(message_type, input_data, neighbour_level, \
                     result_report_level):

      start_time = time.time()  # Get current time

      if (input_data == ''):
        return 'Error:Received empty input data'

      # Set the maximal neighbour level
      #
      self.rec_geocoder.maximal_neighbour_level = neighbour_level

      # Save input data into input data set
      #
      input_data_list = input_data.split('\n')
      org_input_records = {}

      i = 0
      for input_address in input_data_list:
        input_address = input_address.strip()  # Remove whitespaces, new-lines
        if (input_address != ''):
          org_input_records[i] = input_address
          input_record = {'_rec_num_':i, 'address':input_address.lower()}
          self.input_dataset.write_record(input_record)
          i += 1

      num_input_addresses = i

      logging.info('Received %d input address(es)' % (num_input_addresses))

      if (self.input_address_log_file != None):
        addr_log_file = open(self.input_address_log_file, 'a')

        header_str = 'Received %d addresses on %s' % \
                     (num_input_addresses, time.ctime(time.time()))
        addr_log_file.write(header_str+os.linesep)

        for i in range(num_input_addresses):
          addr_dict = self.input_dataset.read_record(i)
          addr_str = '  ' + addr_dict['address']
          addr_log_file.write(addr_str+os.linesep)

        addr_log_file.write(os.linesep)
        addr_log_file.close()

      if (message_type == 'text'):
        result_str = ''  # String to be returned

      elif (message_type == 'file'):

        # Create a CSV file header with the field names
        #
        field_dict =    self.output_dataset.fields
        field_names =   field_dict.keys()
        field_columns = field_dict.values()

        # Create sorted list of fields (according to column numbers)
        #
        field_list = map(None, field_columns, field_names)
        field_list.sort()

        csv_file_header = ''
        for (p,fn) in field_list:
          csv_file_header += fn + ','

        result_str = csv_file_header[:-1]+'\n'

      else:
        logging.exception('Illegale message type received: "%s"' % \
                          (message_type))
        raise Exception

      # Loop over all input addresses - - - - - - - - - - - - - - - - - - - - -
      #
      for i in range(num_input_addresses):

        # Standardise the input records
        #
        input_record = self.input_dataset.read_record(i)
        clean_record = self.rec_standardiser.standardise(input_record)

        logging.debug('Standardised input record %3d: "%s"' % \
                      (i, str(input_record)))
        logging.debug('                         into: "%s"' % \
                      (str(clean_record)))

        # Geocode the standardised record
        #
        geocode_result = self.rec_geocoder.match_record(clean_record)
        logging.debug('Geocoded standardised record: "%s"' % \
                      (str(geocode_result)))

        match_values = str(geocode_result[0][2])
        match_status = geocode_result[0][3]
        result_lati =  geocode_result[0][4]
        result_long =  geocode_result[0][5]

        clean_record['match_score'] =         str(geocode_result[0][0])
        clean_record['gnaf_pid'] =            str(geocode_result[0][1])
        clean_record['match_values'] =       str(match_values.replace(',',';'))
        clean_record['match_status'] =        str(match_status)
        clean_record['latitude'] =            str(result_lati)
        clean_record['longitude'] =           str(result_long)
        clean_record['collection_district'] = str(geocode_result[0][6])
        clean_record['max_avrg_distance'] =   str(geocode_result[0][7])
        clean_record['neighbour_level'] =     str(geocode_result[0][8])

        if (message_type == 'text'):  # Create output for text message  - - -

          # Compose result string (according to result_report_level)
          # The location is returned first (for all result_report_levels)
          #
          result_str = result_str + 'Input record %d:\n' % (i)
          result_str = result_str + '  Original input record:  "%s"\n' % \
                       org_input_records[i]
          result_str = result_str + '  Longitude and latitude: %s, %s\n' % \
                       (str(result_long), str(result_lati))

          # Add standardised record (level 1 and 2)
          #
          if (result_report_level > 0):
            del clean_record['_dataset_name_']  # Remove 'hidden record fields'
            del clean_record['_rec_num_']

            # Format the record dictionary in a nice way
            #
            max_key_lengh = 0
            for k in clean_record.keys():
              if (len(k) > max_key_lengh):
                max_key_lengh = len(k)

            result_str = result_str + '  Standardised record:\n'
            for (k, v) in clean_record.items():
              result_str = result_str + '    ' + ' '*(max_key_lengh-len(k)) + \
                           '"%s": "%s"\n' %(k,v)

          # For level 2, also add meta information
          #
          if  (result_report_level == 2):
            result_str = result_str + '  Matching information: ' + match_status

          result_str = result_str + '\n-------------------------------------\n'

        else:  # Create a CSV file  - - - - - - - - - - - - - - - - - - - - - -

          csv_rec = ''
          for (p,fn) in field_list:
            csv_rec += str(clean_record.get(fn,''))+','

          result_str += csv_rec[:-1]+'\n'

      return result_str

    # -------------------------------------------------------------------------

    # Finally link the request handler routine to the geocode server object
    #
    self.handler = do_geocoding

    return self

  # ---------------------------------------------------------------------------

  def deduplicate(self, **kwargs):
    """Deduplicate the given data set using the defined record standardiser,
       record comparators, blocking indexes and classifiers.

       Records are loaded block wise from the input data set, then standardised
       (if the record standardiser is defined, otherwise the input data set is
       directly deduplicated), linked and the results are printed and/or saved
       into the result file(s).

       If the argument 'first_record' is not given, it will automatically be
       set to the first record in the data set (i.e. record number 0).
       Similarly, if the argument 'number_records' is not given, it will be set
       to the total number of records in the input data set.

       The temporary data set must be a random acces data set implementation,
       i.e. either a Shelve or a Memory data set. For large data set it is
       recommended to use a Shelve data set. This temporary data set has to be
       initialised in access mode 'readwrite'.

       Currently, the output can be a printed or saved list of record pairs in
       both a detailed and condensed form (if the arguments
       'output_rec_pair_details' and 'output_rec_pair_weights' are set to
       'True' or to a file name (a string). The output can be filtered by
       setting the 'output_threshold' (meaning all record pairs with a weight
       less then this threshold are not printed or saved).

       If a file name is given for the 'output_rec_pair_weights' argument, and
       it's file extension is '.csv' then the output will be written as a
       comma separated file (otherwise as a normal text file).

       In future versions, it will be possible to compile an output data set.

       A histogram can be saved or printed by setting the argument
       'output_histogram' to 'True' or to a file name.

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

    self.output_histogram = False         # Set to True, a file name or False
                                          # (default) if a histogram of weights
                                          # should be printed or saved
    self.output_rec_pair_details = False  # Set to True, a file name or False
                                          # (default) if record pairs should
                                          # be printed or saved in details
    self.output_rec_pair_weights = False  # Set to True, a file name or False
                                          # (default) if record pairs should
                                          # be printed or saved with weights
    self.output_threshold = None          # Set to a weight threshold (only
                                          # record pairs with weights equal to
                                          # or above will be saved and or
                                          # printed)
    self.output_assignment = None         # Set to 'one2one' if one-to-one
                                          # assignment should be forced
                                          # (default: None)

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
          logging.exception('Argument "first_record" is not a valid ' + \
                            'integer number')
          raise Exception
        self.first_record = value
      elif (keyword == 'number_records'):
        if (not isinstance(value, int)) or (value <= 0):
          logging.exception('Argument "number_records" is not a positive' + \
                            ' integer number')
          raise Exception
        self.number_records = value

      elif (keyword == 'output_rec_pair_details'):
        if (not isinstance(value, str)) and (value not in [True, False]):
          logging.exception('Argument "output_rec_pair_details" must be ' + \
                            'a file name or "True" or "False"')
          raise Exception
        self.output_rec_pair_details = value
      elif (keyword == 'output_rec_pair_weights'):
        if (not isinstance(value, str)) and (value not in [True, False]):
          logging.exception('Argument "output_rec_pair_weights" must be ' + \
                            'a file name or "True" or "False"')
          raise Exception
        self.output_rec_pair_weights = value
      elif (keyword == 'output_histogram'):
        if (not isinstance(value, str)) and (value not in [True, False]):
          logging.exception('Argument "output_histogram" must be ' + \
                            'a file name or "True" or "False"')
          raise Exception
        self.output_histogram = value
      elif (keyword == 'output_threshold'):
        if (not (isinstance(value, int) or isinstance(value, float))):
          logging.exception('Argument "output_threshold" is not a number: %s' \
                            % (str(value)))
        self.output_threshold = value
      elif (keyword == 'output_assignment'):
        if (value not in ['one2one', None]):
          logging.exception('Illegal value for argument "output_assignment"' \
                            + ': %s' % (str(value)))
          raise Exception
        else:
          self.output_assignment = value

      else:
        logging.exception('Illegal constructor argument keyword: "%s"' % \
                          (str(keyword)))
        raise Exception

    # Do some checks on the input arguments - - - - - - - - - - - - - - - - - -
    #
    if (self.input_dataset == None):
      logging.exception('Input data set is not defined')
      raise Exception

    if (self.tmp_dataset == None):
      logging.exception('Temporary data set is not defined')
      raise Exception
    elif (self.tmp_dataset.dataset_type not in ['SHELVE', 'MEMORY']):
      logging.exception('Temporary data set must be a random access data ' + \
                        'set (either Shelve or Memory)')
      raise Exception
    if (self.tmp_dataset.access_mode not in ['write','append','readwrite']):
      logging.exception('Temporary data set must be initialised in one of ' + \
                        'the access  modes: "write", "append", or "readwrite"')
      raise Exception

    # Make sure at least one output is defined
    #
    if (self.output_rec_pair_weights == False) and \
       (self.output_rec_pair_details == False) and \
       (self.output_histogram == False):
      logging.exception('No ouput of results is defined')
      raise Exception
    #
    # Code above to be removed once output data set functionality implemented

    if (self.first_record == None):
      self.first_record = 0  # Take default first record in data set

    if (self.number_records == None):  # Process all records
      self.number_records = self.input_dataset.num_records

    if (self.rec_comparator == None):
      logging.exception('No record comparator defined')
      raise Exception
    if (self.rec_comparator.dataset_a != self.tmp_dataset) or \
       (self.rec_comparator.dataset_b != self.tmp_dataset):
      logging.exception('Illegal data set definition in record comparator')
      raise Exception

    if (self.blocking_index == None):
      logging.exception('No blocking index defined')
      raise Exception

    if (self.classifier == None):
      logging.exception('No classifier defined')
      raise Exception
    if (self.classifier.dataset_a != self.tmp_dataset) or \
       (self.classifier.dataset_b != self.tmp_dataset):
      logging.exception('Illegal data set definition in classifier')
      raise Exception

    total_time = time.time()  # Get current time

    logging.info('***** Deduplicate data set: "%s" (type "%s")' % \
                 (self.input_dataset.name, self.input_dataset.dataset_type))
    logging.info('*****   Temporary data set: "%s" (type "%s")' % \
                 (self.tmp_dataset.name, self.tmp_dataset.dataset_type))
    logging.info('Step 1: Loading, standardisation and indexing')
    logging.info('-------')

    step_1_time = time.time()  # Get current time

    # Call the main standardisation routine - - - - - - - - - - - - - - - - - -
    #
    [p, step_1_comm_time] = do_load_standard_indexing(self.input_dataset,
                                                      self.tmp_dataset,
                                                      self.rec_standardiser,
                                                      self.blocking_index,
                                                      self.first_record,
                                                      self.number_records,
                                                      self.block_size)

    # If Febrl is run in parallel, collect blocking index in process 0  - - - -
    #
    if (parallel.rank() == 0):
      for p in range(1, parallel.size()):
        tmp_time = time.time()
        tmp_indexes = parallel.receive(p)
        step_1_comm_time += (time.time() - tmp_time)
        logging.info('    Received index from process %i' % (p))

        self.blocking_index.merge(tmp_indexes)

    else:  # Send index to process 0
      tmp_time = time.time()
      parallel.send(self.blocking_index.index, 0) # Send indexes to process 0
      step_1_comm_time += (time.time() - tmp_time)
      logging.info('    Sent index to process 0')

    # If run in parallel, broadcast the blocking index from process 0 - - - - -
    #
    if (parallel.size() > 1):
      if (parallel.rank() == 0):
        for p in range(1, parallel.size()):
          tmp_time = time.time()
          parallel.send(self.blocking_index.index, p)
          step_1_comm_time += (time.time() - tmp_time)
          logging.info('    Sent index to process %i' % (p))

      else:
        tmp_time = time.time()
        tmp_indexes = parallel.receive(0)
        step_1_comm_time += (time.time() - tmp_time)
        logging.info('    Received index from process 0')

        self.blocking_index.merge(tmp_indexes)

    # Compact the blocking index  - - - - - - - - - - - - - - - - - - - - - - -
    #
    self.blocking_index.compact()

    step_1_time = time.time() - step_1_time  # Calculate time for step 1
    step_1_time_string = output.time_string(step_1_time)

    logging.info('Step 1 finished in %s' % (step_1_time_string))

    # End of step 1 - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    #
    parallel.Barrier()  # Make sure all processes are here

    # Now re-initialise the temporary data set in read access mode only - - - -
    #
    self.tmp_dataset.re_initialise('read')

    #################### START PARALLEL TEST CODE #############################
    # Save temporary data sets and indexes to files (on all processes)
    #
    if (SAVE_PARALLEL_TEST_FILES == True):
      #f = open('tmp_data_set-dedup-'+str(parallel.rank())+'-'+ \
      #         str(parallel.size()),'w')
      #tmp_list = self.tmp_dataset.dict.keys()
      #tmp_list.sort()

      #for r in tmp_list:
      #  rec = self.tmp_dataset.dict[r]
      #  rec_items = rec.items()
      #  rec_items.sort()
      #  rec = str(r)+': '+str(rec_items)
      #  f.write(rec+os.linesep)
      #f.close()

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

    logging.info('Step 2: Perform deduplication within blocks')
    logging.info('-------')

    step_2_time = time.time()  # Get current time
    step_2_comm_time = 0.0

    # Get the record pairs which have to be compared  - - - - - - - - - - - - -
    #
    tmp_time = time.time()
    [rec_pair_dict, rec_pair_cnt] = \
           indexing.deduplication_rec_pairs(self.blocking_index)
    rec_pair_time = time.time() - tmp_time
    rec_pair_time_string = output.time_string(rec_pair_time)

    logging.info('  Built record pair dictionary with %i entries in %s' % \
                 (rec_pair_cnt, rec_pair_time_string))

    # And do the comparisons of record pairs into classifer - - - - - - - - - -
    #
    [p] = do_comparison(self.tmp_dataset, self.tmp_dataset,
                        self.rec_comparator, self.classifier,
                        rec_pair_dict, rec_pair_cnt, self.block_size)

    # Now gather classifier results on process 0 and merge  - - - - - - - - - -
    #
    if (parallel.size() > 1):
      if (parallel.rank() == 0):
        for p in range(1, parallel.size()):
          tmp_time = time.time()
          tmp_classifier_results = parallel.receive(p)
          step_2_comm_time += (time.time() - tmp_time)
          logging.info('    Received classifier from process %i and merge it' \
                       % (p))
          self.classifier.merge(tmp_classifier_results)

      else:
        tmp_time = time.time()
        parallel.send(self.classifier.results, 0) # Send classifier results to
                                                  # process 0
        step_2_comm_time += (time.time() - tmp_time)
        logging.info('    Sent classifier to process 0')

    # If run in parallel, broadcast the classifier results from process 0 - - -
    #
    if (parallel.size() > 1):
      if (parallel.rank() == 0):
        for p in range(1, parallel.size()):
          tmp_time = time.time()
          parallel.send(self.classifier.results, p)
          step_2_comm_time += (time.time() - tmp_time)
          logging.info('    Sent classifier to process %i' % (p))

      else:
        tmp_time = time.time()
        self.classifier.results = parallel.receive(0)
        step_2_comm_time += (time.time() - tmp_time)
        logging.info('    Received classifier from process 0')

    #################### START PARALLEL TEST CODE #############################
    # Save classifiers and weight vectors to files (only process 0)
    #
    if (SAVE_PARALLEL_TEST_FILES == True):
      tmp_list = rec_pair_dict.keys()
      tmp_list.sort()
      f = open('rec-pair-dict-dedup-'+str(parallel.rank())+'-'+ \
               str(parallel.size()),'w')
      for rp in tmp_list:
        rec_list = rec_pair_dict[rp].items()
        rec_list.sort()
        r = str(rp)+': '+str(rec_list)
        f.write(r+os.linesep)
      f.close()

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

    #################### END PARALLEL TEST CODE ###############################

    step_2_time = time.time() - step_2_time  # Calculate time for step 2
    step_2_time_string = output.time_string(step_2_time)

    logging.info('Step 2 (deduplication) finished in %s' % \
                 (step_2_time_string))
    logging.info('  Totally %i record pair comparisons' % (rec_pair_cnt))

    # Output the results  - - - - - - - - - - - - - - - - - - - - - - - - - - -

    logging.info('Step 3: Output and assignment procedures')
    logging.info('-------')

    step_3_time = time.time()  # Get current time

    # Get the results dictionary with all the record pairs and their weights
    #
    results_dict = self.classifier.results

    if (results_dict == {}):
      logging.warn('Results dictionary empty')

    else:  # There are results

      # Do assignment restrictions if they are defined  - - - - - - - - - - - -
      #
      if (self.output_assignment != None):  # An output assignment is defined

        if (self.output_assignment == 'one2one'):

          # Do a one-to-one assignment on the classifier results dict
          #
          o2o_results_dict = lap.do_lap('auction', results_dict, \
                                        'deduplication', self.output_threshold)
      else:  # No one-to-one assignment, set one2one result to None
        o2o_results_dict = None

      #################### START PARALLEL TEST CODE ###########################

      if (SAVE_PARALLEL_TEST_FILES == True):
        tmp_list = o2o_results_dict.items()
        tmp_list.sort()
        f = open('one2one-dedup-'+str(parallel.rank())+'-'+ \
                 str(parallel.size()),'w')
        for c in tmp_list:
          f.write(str(c)+os.linesep)
        f.close()

      #################### END PARALLEL TEST CODE #############################

      if (parallel.rank() == 0):  # Only processor 0 prints results

        # Print or save weights histogram - - - - - - - - - - - - - - - - - - -
        #
        if (self.output_histogram == True):
          output.histogram(results_dict)
        elif (self.output_histogram != False):
          output.histogram(results_dict, self.output_histogram)

        # Print or save detailed record pairs - - - - - - - - - - - - - - - - -
        #
        if (self.output_rec_pair_details == True):
          output.rec_pair_details(self.tmp_dataset, self.tmp_dataset,
                                  results_dict, o2o_results_dict,
                                  self.output_threshold)
        elif (self.output_rec_pair_details != False):
          output.rec_pair_details(self.tmp_dataset, self.tmp_dataset,
                                  results_dict, o2o_results_dict,
                                  self.output_threshold,
                                  self.output_rec_pair_details)

        # Print or save record pairs with weights - - - - - - - - - - - - - - -
        #
        if (self.output_rec_pair_weights == True):
          output.rec_pair_weights(self.tmp_dataset.name,
                                  self.tmp_dataset.name,
                                  results_dict, o2o_results_dict,
                                  self.output_threshold)
        elif (self.output_rec_pair_weights != False):
          output.rec_pair_weights(self.tmp_dataset.name,
                                  self.tmp_dataset.name,
                                  results_dict, o2o_results_dict,
                                  self.output_threshold,
                                  self.output_rec_pair_weights)

    step_3_time = time.time() - step_3_time  # Calculate time for step 3
    step_3_time_string = output.time_string(step_3_time)

    logging.info('Step 3 (output and assignments) finished in %s' % \
                 (step_3_time_string))

    parallel.Barrier()  # Wait here for all processes - - - - - - - - - - - - -

    total_time = time.time() - total_time  # Calculate total time

    total_time_string =       output.time_string(total_time)
    step_1_comm_time_string = output.time_string(step_1_comm_time)
    step_2_comm_time_string = output.time_string(step_2_comm_time)

    logging.info('Total time needed for deduplication of %i records: %s' % \
                 (self.number_records, total_time_string))
    logging.info('  Time for step 1 (standardisation):       %s' % \
                 (step_1_time_string))
    logging.info('  Time for step 2 (deduplication):         %s' % \
                 (step_2_time_string))
    logging.info('  Time for step 3 (assignment and output): %s' % \
                 (step_3_time_string))
    logging.info('  Time for communication in step 1: %s' % \
                 (step_1_comm_time_string))
    logging.info('  Time for communication in step 2: %s' % \
                 (step_2_comm_time_string))

  # ---------------------------------------------------------------------------

  def link(self, **kwargs):
    """Link the given two data set using the defined record standardisers,
       record comparators, blocking indexes and classifiers.

       Records are loaded block wise from the input data sets, then
       standardised (if the record standardisers are defined, otherwise an
       input data set is directly taken for the linkage process), linked and
       the results are printed and/or saved into the result file(s).

       If the arguments 'first_record_a' and/or 'first_record_b' are not given,
       they will automatically be set to the first record in the data sets
       (i.e. record number 0). Similarly, if the arguments 'number_records_a'
       and/or 'number_records_b' are not given, they will be set to the total
       number of records in the input data sets.

       The temporary data sets must be random acces data set implementations,
       i.e. either Shelve or a Memory data sets. For large data set it is
       recommended to use Shelve data sets. This temporary data sets have to be
       initialised in access mode 'readwrite'.

       Currently, the output can be a printed or saved list of record pairs in
       both a detailed and condensed form (if the arguments
       'output_rec_pair_details' and 'output_rec_pair_weights' are set to
       'True' or to a file name (a string). The output can be filtered by
       setting the 'output_threshold' (meaning all record pairs with a weight
       less then this threshold are not printed or saved).

       If a file name is given for the 'output_rec_pair_weights' argument, and
       it's file extension is '.csv' then the output will be written as a
       comma separated file (otherwise as a normal text file).

       In future versions, it will be possible to compile an output data set.

       A histogram can be saved or printed by setting the argument
       'output_histogram' to 'True' or to a file name.

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

    self.output_histogram = False         # Set to True, a file name or False
                                          # (default) if a histogram of weights
                                          # should be printed or saved
    self.output_rec_pair_details = False  # Set to True, a file name or False
                                          # (default) if record pairs should
                                          # be printed or saved in details
    self.output_rec_pair_weights = False  # Set to True, a file name or False
                                          # (default) if record pairs should
                                          # be printed or saved with weights
    self.output_threshold = None          # Set to a weight threshold (only
                                          # record pairs with weights equal to
                                          # or above will be saved and or
                                          # printed)
    self.output_assignment = None         # Set to 'one2one' if one-to-one
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
          logging.exception('Argument "first_record_a" is not a valid ' + \
                            'integer number')
          raise Exception
        self.first_record_a = value
      elif (keyword == 'first_record_b'):
        if (not isinstance(value, int)) or (value < 0):
          logging.exception('Argument "first_record_b" is not a valid ' + \
                            'integer number')
          raise Exception
        self.first_record_b = value
      elif (keyword == 'number_records_a'):
        if (not isinstance(value, int)) or (value <= 0):
          logging.exception('Argument "number_records_a" is not a positive '+ \
                            'integer number')
          raise Exception
        self.number_records_a = value
      elif (keyword == 'number_records_b'):
        if (not isinstance(value, int)) or (value <= 0):
          logging.exception('Argument "number_records_b" is not a positive '+ \
                            'integer number')
          raise Exception
        self.number_records_b = value

      elif (keyword == 'output_rec_pair_details'):
        if (not isinstance(value, str)) and (value not in [True, False]):
          logging.exception('Argument "output_rec_pair_details" must be ' + \
                            'a file name or "True" or "False"')
          raise Exception
        self.output_rec_pair_details = value
      elif (keyword == 'output_rec_pair_weights'):
        if (not isinstance(value, str)) and (value not in [True, False]):
          logging.exception('Argument "output_rec_pair_weights" must be ' + \
                            'a file name or "True" or "False"')
          raise Exception
        self.output_rec_pair_weights = value
      elif (keyword == 'output_histogram'):
        if (not isinstance(value, str)) and (value not in [True, False]):
          logging.exception('Argument "output_histogram" must be ' + \
                            'a file name or "True" or "False"')
          raise Exception
        self.output_histogram = value
      elif (keyword == 'output_threshold'):
        if (not (isinstance(value, int) or isinstance(value, float))):
          logging.exception('Argument "output_threshold" is not a number: %s' \
                            % (str(value)))
        self.output_threshold = value
      elif (keyword == 'output_assignment'):
        if (value not in ['one2one', None]):
          logging.exception('Illegal value for argument "output_assignment"' \
                            + ': %s' % (str(value)))
          raise Exception
        else:
          self.output_assignment = value

      else:
        logging.exception('Illegal constructor argument keyword: "%s"' % \
                          (str(keyword)))
        raise Exception

    # Check if the needed attributes are set  - - - - - - - - - - - - - - - - -
    #
    if (self.input_dataset_a == None):
      logging.exception('Input data set A is not defined')
      raise Exception
    if (self.input_dataset_b == None):
      logging.exception('Input data set B is not defined')
      raise Exception

    if (self.tmp_dataset_a == None):
      logging.exception('Temporary data set A is not defined')
      raise Exception
    elif (self.tmp_dataset_a.dataset_type not in ['SHELVE', 'MEMORY']):
      logging.exception('Temporary data set A must be a random access data' + \
                        'set (either Shelve or Memory)')
      raise Exception
    if (self.tmp_dataset_a.access_mode not in ['write','append','readwrite']):
      logging.exception('Temporary data set A must be initialised in one ' + \
                        'of the access  modes: "write", "append", or ' + \
                        '"readwrite"')
      raise Exception

    if (self.tmp_dataset_b == None):
      logging.exception('Temporary data set B is not defined')
      raise Exception
    elif (self.tmp_dataset_b.dataset_type not in ['SHELVE', 'MEMORY']):
      logging.exception('Temporary data set B must be a random access ' + \
                        'data set (either Shelve or Memory)')
      raise Exception
    if (self.tmp_dataset_b.access_mode not in ['write','append','readwrite']):
      logging.exception('Temporary data set B must be initialised in one ' + \
                        'of the access  modes: "write", "append", or ' + \
                        '"readwrite"')
      raise Exception

    # Check if there are file names for the temporary data sets and - - - - - -
    # if they differ
    #
    tmp_file_name_a = getattr(self.tmp_dataset_a, 'file_name', None)
    tmp_file_name_b = getattr(self.tmp_dataset_b, 'file_name', None)
    if (tmp_file_name_a != None) and (tmp_file_name_b != None):
      if (tmp_file_name_a == tmp_file_name_b):
        logging.exception('The same file names for both temporary data sets')
        raise Exception

    # Make sure at least one output is defined
    #
    if (self.output_rec_pair_weights == False) and \
       (self.output_rec_pair_details == False) and \
       (self.output_histogram == False):
      logging.exception('No ouput of results is defined')
      raise Exception
    #
    # Code above to be removed once output data set functionality implemented

    if (self.first_record_a == None):
      self.first_record_a = 0  # Take default first record in data set
    if (self.first_record_b == None):
      self.first_record_b = 0  # Take default first record in data set

    if (self.number_records_a == None):  # Take all records
      self.number_records_a = self.input_dataset_a.num_records
    if (self.number_records_b == None):  # Take all records
      self.number_records_b = self.input_dataset_b.num_records

    if (self.rec_comparator == None):
      logging.exception('No record comparator defined')
      raise Exception
    if (self.rec_comparator.dataset_a != self.tmp_dataset_a) or \
       (self.rec_comparator.dataset_b != self.tmp_dataset_b):
      logging.exception('Illegal data set definition in record comparator')
      raise Exception

    if (self.blocking_index_a == None):
      logging.exception('No blocking index for data set A defined')
      raise Exception
    if (self.blocking_index_b == None):
      logging.exception('No blocking index for data set B defined')
      raise Exception

    if (self.classifier == None):
      logging.exception('No classifier defined')
      raise Exception
    if (self.classifier.dataset_a != self.tmp_dataset_a) or \
       (self.classifier.dataset_b != self.tmp_dataset_b):
      logging.exception('Illegal data set definition in classifier')
      raise Exception

    total_time = time.time()  # Get current time

    logging.info('***** Link data set: %s (type "%s") with data set: %s ' % \
                 (self.input_dataset_a.name, \
                 self.input_dataset_a.dataset_type, \
                 self.input_dataset_b.name) + '(type "%s")' % \
                 (self.input_dataset_b.dataset_type))
    logging.info('*****   Temporary data set A: "%s" (type "%s")' % \
                 (self.tmp_dataset_a.name, self.tmp_dataset_a.dataset_type))
    logging.info('*****   Temporary data set B: "%s" (type "%s")' % \
                 (self.tmp_dataset_b.name, self.tmp_dataset_b.dataset_type))
    logging.info('Step 1: Loading, standardisation and indexing')
    logging.info('-------')

    step_1_time = time.time()  # Get current time

    # Call the main standardisation routine for data set A  - - - - - - - - - -
    #
    [p, step_1_comm_time] = do_load_standard_indexing(self.input_dataset_a,
                                                      self.tmp_dataset_a,
                                                      self.rec_standardiser_a,
                                                      self.blocking_index_a,
                                                      self.first_record_a,
                                                      self.number_records_a,
                                                      self.block_size)

    # If Febrl is run in parallel, collect blocking index in process 0  - - - -
    #
    if (parallel.rank() == 0):
      for p in range(1, parallel.size()):
        tmp_time = time.time()
        tmp_indexes = parallel.receive(p)
        step_1_comm_time += (time.time() - tmp_time)
        logging.info('    Received index A from process %i' % (p))

        self.blocking_index_a.merge(tmp_indexes)

    else:  # Send index to process 0
      tmp_time = time.time()
      parallel.send(self.blocking_index_a.index, 0) # Send indexes to process 0
      step_1_comm_time += (time.time() - tmp_time)
      logging.info('    Sent index A to process 0')

    # If run in parallel, broadcast the blocking index from process 0 - - - - -
    #
    if (parallel.size() > 1):
      if (parallel.rank() == 0):
        for p in range(1, parallel.size()):
          tmp_time = time.time()
          parallel.send(self.blocking_index_a.index, p)
          step_1_comm_time += (time.time() - tmp_time)
          logging.info('    Sent index A to process %i' % (p))

      else:
        tmp_time = time.time()
        tmp_indexes = parallel.receive(0)
        step_1_comm_time += (time.time() - tmp_time)
        logging.info('    Received index A from process 0')

        self.blocking_index_a.merge(tmp_indexes)

    # Compact the blocking index  - - - - - - - - - - - - - - - - - - - - - - -
    #
    self.blocking_index_a.compact()

    # Call the main standardisation routine for data set B  - - - - - - - - - -
    #
    [p, tmp_time] = do_load_standard_indexing(self.input_dataset_b,
                                              self.tmp_dataset_b,
                                              self.rec_standardiser_b,
                                              self.blocking_index_b,
                                              self.first_record_b,
                                              self.number_records_b,
                                              self.block_size)
    step_1_comm_time += tmp_time

    # If Febrl is run in parallel, collect blocking index in process 0  - - - -
    #
    if (parallel.rank() == 0):
      for p in range(1, parallel.size()):
        tmp_time = time.time()
        tmp_indexes = parallel.receive(p)
        step_1_comm_time += (time.time() - tmp_time)
        logging.info('    Received index B from process %i' % (p))

        self.blocking_index_b.merge(tmp_indexes)

    else:  # Send index to process 0
      tmp_time = time.time()
      parallel.send(self.blocking_index_b.index, 0) # Send indexes to process 0
      step_1_comm_time += (time.time() - tmp_time)
      logging.info('    Sent index B to process 0')

    # If run in parallel, broadcast the blocking index from process 0 - - - - -
    #
    if (parallel.size() > 1):
      if (parallel.rank() == 0):
        for p in range(1, parallel.size()):
          tmp_time = time.time()
          parallel.send(self.blocking_index_b.index, p)
          step_1_comm_time += (time.time() - tmp_time)
          logging.info('    Sent index B to process %i' % (p))

      else:
        tmp_time = time.time()
        tmp_indexes = parallel.receive(0)
        step_1_comm_time += (time.time() - tmp_time)
        logging.info('    Received index B from process 0')

        self.blocking_index_b.merge(tmp_indexes)

    # Compact the blocking index  - - - - - - - - - - - - - - - - - - - - - - -
    #
    self.blocking_index_b.compact()

    step_1_time = time.time() - step_1_time  # Calculate time for step 1
    step_1_time_string = output.time_string(step_1_time)

    logging.info('Step 1 finished in %s' % (step_1_time_string))

    # End of step 1 - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    #
    parallel.Barrier()  # Make sure all processes are here

    # Now re-initialise the temporary data sets in read access mode only  - - -
    #
    self.tmp_dataset_a.re_initialise('read')
    self.tmp_dataset_b.re_initialise('read')

    #################### START PARALLEL TEST CODE #############################
    # Save temporary data sets and indexes to files (on all processes)
    #
    if (SAVE_PARALLEL_TEST_FILES == True):
      #f = open('tmp_data_set_a-link-'+str(parallel.rank())+'-'+ \
      #         str(parallel.size()),'w')
      #tmp_list = self.tmp_dataset_a.dict.keys()
      #tmp_list.sort()

      #for r in tmp_list:
      #  rec = self.tmp_dataset_a.dict[r]
      #  rec_items = rec.items()
      #  rec_items.sort()
      #  rec = str(r)+': '+str(rec_items)
      #  f.write(rec+os.linesep)
      #f.close()

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

    logging.info('Step 2: Perform linkage within blocks')
    logging.info('-------')

    step_2_time = time.time()  # Get current time
    step_2_comm_time = 0.0

    # Get the record pairs which have to be compared  - - - - - - - - - - - - -
    #
    tmp_time = time.time()
    [rec_pair_dict, rec_pair_cnt] = \
           indexing.linkage_rec_pairs(self.blocking_index_a,
                                      self.blocking_index_b)
    rec_pair_time = time.time() - tmp_time
    rec_pair_time_string = output.time_string(rec_pair_time)

    logging.info('  Built record pair dictionary with %i entries in %s' % \
                 (rec_pair_cnt, rec_pair_time_string))

    # And do the comparisons of record pairs into classifer - - - - - - - - - -
    #
    [p] = do_comparison(self.tmp_dataset_a, self.tmp_dataset_b,
                        self.rec_comparator, self.classifier,
                        rec_pair_dict, rec_pair_cnt, self.block_size)

    # Now gather classifier results on process 0 and merge  - - - - - - - - - -
    #
    if (parallel.size() > 1):
      if (parallel.rank() == 0):
        for p in range(1, parallel.size()):
          tmp_time = time.time()
          tmp_classifier_results = parallel.receive(p)
          step_2_comm_time += (time.time() - tmp_time)
          logging.info('    Received classifier from process %i and merge it' \
                       % (p))

          self.classifier.merge(tmp_classifier_results)

      else:
        tmp_time = time.time()
        parallel.send(self.classifier.results, 0) # Send classifier results to
                                                  # process 0
        step_2_comm_time += (time.time() - tmp_time)
        logging.info('    Sent classifier to process 0')

    # If run in parallel, broadcast the classifier results from process 0 - - -
    #
    if (parallel.size() > 1):
      if (parallel.rank() == 0):
        for p in range(1, parallel.size()):
          tmp_time = time.time()
          parallel.send(self.classifier.results, p)
          step_2_comm_time += (time.time() - tmp_time)
          logging.info('    Sent classifier to process %i' % (p))

      else:
        tmp_time = time.time()
        self.classifier.results = parallel.receive(0)
        step_2_comm_time += (time.time() - tmp_time)
        logging.info('    Received classifier from process 0')

    #################### START PARALLEL TEST CODE #############################
    # Save classifiers and weight vectors to files (only process 0)
    #
    if (SAVE_PARALLEL_TEST_FILES == True):
      tmp_list = rec_pair_dict.keys()
      tmp_list.sort()
      f = open('rec-pair-dict-link-'+str(parallel.rank())+'-'+ \
               str(parallel.size()),'w')
      for rp in tmp_list:
        rec_list = rec_pair_dict[rp].items()
        rec_list.sort()
        r = str(rp)+': '+str(rec_list)
        f.write(r+os.linesep)
      f.close()

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

    #################### END PARALLEL TEST CODE ###############################

    step_2_time = time.time() - step_2_time  # Calculate time for step 2
    step_2_time_string = output.time_string(step_2_time)

    logging.info('Step 2 (linkage) finished in %s' % (step_2_time_string))
    logging.info('  Totally %i record pair comparisons' % (rec_pair_cnt))

    # Output the results  - - - - - - - - - - - - - - - - - - - - - - - - - - -

    logging.info('Step 3: Output and assignment procedures')
    logging.info('-------')

    step_3_time = time.time()  # Get current time

    # Get the results dictionary with all the record pairs and their weights
    #
    results_dict = self.classifier.results

    if (results_dict == {}):
      logging.warn('Results dictionary empty')

    else:  # There are results

      # Do assignment restrictions if they are defined  - - - - - - - - - - - -
      #
      if (self.output_assignment != None):  # An output assignment is defined

        if (self.output_assignment == 'one2one'):

          # Do a one-to-one assignment on the classifier results dict
          #
          o2o_results_dict = lap.do_lap('auction', results_dict, \
                                        'linkage', self.output_threshold)
      else:  # No one-to-one assignment, set one2one result to None
        o2o_results_dict = None

      #################### START PARALLEL TEST CODE ###########################

      if (SAVE_PARALLEL_TEST_FILES == True):
        tmp_list = o2o_results_dict.items()
        tmp_list.sort()
        f = open('one2one-link-'+str(parallel.rank())+'-'+ \
                 str(parallel.size()),'w')
        for c in tmp_list:
          f.write(str(c)+os.linesep)
        f.close()

      #################### END PARALLEL TEST CODE #############################

      if (parallel.rank() == 0):  # Only processor 0 prints results

        # Print or save weights histogram - - - - - - - - - - - - - - - - - - -
        #
        if (self.output_histogram == True):
          output.histogram(results_dict)
        elif (self.output_histogram != False):
          output.histogram(results_dict, self.output_histogram)

        # Print or save detailed record pairs - - - - - - - - - - - - - - - - -
        #
        if (self.output_rec_pair_details == True):
          output.rec_pair_details(self.tmp_dataset_a, self.tmp_dataset_b, \
                                  results_dict, o2o_results_dict, \
                                  self.output_threshold)
        elif (self.output_rec_pair_details != False):
          output.rec_pair_details(self.tmp_dataset_a, self.tmp_dataset_b, \
                                  results_dict, o2o_results_dict, \
                                  self.output_threshold, \
                                  self.output_rec_pair_details)

        # Print or save record pairs with weights - - - - - - - - - - - - - - -
        #
        if (self.output_rec_pair_weights == True):
          output.rec_pair_weights(self.tmp_dataset_a.name, \
                                  self.tmp_dataset_b.name, \
                                  results_dict, o2o_results_dict, \
                                  self.output_threshold)
        elif (self.output_rec_pair_weights != False):
          output.rec_pair_weights(self.tmp_dataset_a.name, \
                                  self.tmp_dataset_b.name, \
                                  results_dict, o2o_results_dict, \
                                  self.output_threshold, \
                                  self.output_rec_pair_weights)

    step_3_time = time.time() - step_3_time  # Calculate time for step 3
    step_3_time_string = output.time_string(step_3_time)

    logging.info('Step 3 (output and assignments) finished in %s' % \
                 (step_3_time_string))

    parallel.Barrier()  # Wait here for all processes - - - - - - - - - - - - -

    total_time = time.time() - total_time  # Calculate total time

    total_time_string =       output.time_string(total_time)
    step_1_comm_time_string = output.time_string(step_1_comm_time)
    step_2_comm_time_string = output.time_string(step_2_comm_time)

    logging.info('Total time needed for linkage of %i records with %i ' % \
                 (self.number_records_a, self.number_records_b,) + \
                 'records: %s' % (total_time_string))
    logging.info('  Time for step 1 (standardisation):       %s' % \
                 (step_1_time_string))
    logging.info('  Time for step 2 (linkage):               %s' % \
                 (step_2_time_string))
    logging.info('  Time for step 3 (assignment and output): %s' % \
                 (step_3_time_string))
    logging.info('  Time for communication in step 1: %s' % \
                 (step_1_comm_time_string))
    logging.info('  Time for communication in step 2: %s' % \
                 (step_2_comm_time_string))

# =============================================================================
# The following are the main routines for standardisation, geocoding and record
# linkage as used withing the project methods 'standardise', 'geocode',
# 'deduplicate' and 'link'

def do_load_standard_indexing(input_dataset, output_dataset,
                              record_standardiser, blocking_index,
                              first_record, number_records,
                              febrl_block_size):

  """The main routine that does the loading of records from the input data set
     into the output data set, if a record standardiser is given these records
     are cleaned and standardised, and if a blocking index is given such a
     index will be built and returned as well.

     If no standardisation is needed the argument 'record_standardiser' has to
     be set to None, and if no indexing is needed the 'blocking_index' argument
     has to be set to None.

     The output dataset must be initialised in access mode "write", "append" or
     "readwrite".

     If run in parallel, all processes will open the input data set and read
     records from it (and clean and standardised them), and the processed
     records will then be sent to process 0 (the host process) and saved into
     the output data set.
  """

  # Check if the given record numbers are valid - - - - - - - - - - - - - - - -
  #
  last_record = first_record + number_records

  if (first_record < 0) or (last_record > input_dataset.num_records):
    logging.exception('Record range too large: (%i,%i)' % \
                      (first_record, last_record))
    raise Exception

  # Check if the data sets are set correctly within the record standardiser - -
  #
  if (record_standardiser != None):
    if (record_standardiser.input_dataset != input_dataset):
      logging.exception('Illegal input data set definition in record '+ \
                        'standardiser: %s (should be: %s)' % \
                        (str(record_standardiser.input_dataset.name), \
                        str(input_dataset.name)))
      raise Exception
    if (record_standardiser.output_dataset != output_dataset):
      logging.exception('Illegal output data set definition in record '+ \
                        'standardiser: %s (should be: %s)' % \
                        (str(record_standardiser.output_dataset.name), \
                        str(output_dataset.name)))
      raise Exception

  else:
    # If no record standardiser for the data set is defined, the field names in
    # the input and the output data sets must be the same
    #
    input_field_name_list = input_dataset.fields.keys()
    input_field_name_list.sort()
    output_field_name_list = output_dataset.fields.keys()
    output_field_name_list.sort()

    if (input_field_name_list != output_field_name_list):
      logging.exception('Field names differ in input and output data sets ' + \
                        '(with no record standardiser defined)')
      raise Exception

  # Check if the data set defined in the blocking index (if defined) is correct
  #
  if (blocking_index != None):
    if (blocking_index.dataset != output_dataset):
      logging.exception('Illegal data set definition in blocking index')
      raise Exception

  if (record_standardiser != None) and (blocking_index != None):
    do_string = 'Load, standardise and index records'
  elif (record_standardiser != None):
    do_string = 'Load and standardise records'
  elif (blocking_index != None):
    do_string = 'Load and index records'
  else:
    do_string = 'Load records'

  logging.info('  %s, then write them into output data set' % (do_string))
  logging.info('    First record: %i' % (first_record))
  logging.info('    Last record:  %i' % (last_record-1))

  start_time = time.time()  # Get current time
  comm_time  = 0.0          # Communication time

  input_rec_counter = first_record  # Current record pointer

  block_cnt = 0  # A round robin block counter, used for parallelism

  # Load records in a blocked fashion - - - - - - - - - - - - - - - - - - - - -

  while (input_rec_counter < last_record):

    block_size = min(febrl_block_size, (last_record - input_rec_counter))

    # Distribute blocks equally to all processors
    #
    if ((block_cnt % parallel.size()) == parallel.rank()):

      # Load original records from input data set
      #
      input_recs = input_dataset.read_records(input_rec_counter, block_size)
      logging.info('    Loaded records %i to %i' % \
                   (input_rec_counter, input_rec_counter+block_size))

      # Standardise them if a standardiser is defined - - - - - - - - - - - - -
      #
      if (record_standardiser != None):
        clean_recs = record_standardiser.standardise_block(input_recs)
        logging.info('      Standardised records %i to %i' % \
                     (input_rec_counter, input_rec_counter+block_size))
      else:
        clean_recs = input_recs  # Take the original records directly

      # Insert records into the blocking index if blocking index is defined - -
      #
      if (blocking_index != None):
        blocking_index.build(clean_recs)

      # If Febrl is run in parallel, send cleaned records to process 0
      #
      if (parallel.rank() > 0):
        tmp_time = time.time()
        parallel.send(clean_recs, 0)
        comm_time += (time.time() - tmp_time)

      else:  # Process 0, store standardised records

        # Store records in output data set
        #
        output_dataset.write_records(clean_recs)

    # If Febrl is run in parallel, process 0 receives cleaned records
    #
    if (parallel.rank() == 0) and (block_cnt % parallel.size() != 0):

      p = (block_cnt % parallel.size())  # Process number to receive from
      tmp_time = time.time()
      tmp_recs = parallel.receive(p)
      comm_time += (time.time() - tmp_time)

      # Store records in output data set
      #
      output_dataset.write_records(tmp_recs)

    input_rec_counter += block_size  # Increment current record pointer
    block_cnt += 1

    # Now determine timing and print progress report  - - - - - - - - - - - - -
    #
    if ((block_cnt % parallel.size()) == 0):
      used_time = time.time() - start_time
      recs_done = input_rec_counter - first_record
      perc_done = 100.0 * recs_done / number_records
      rec_time  = used_time / recs_done
      todo_time = (number_records - recs_done) * rec_time

      used_time_string = output.time_string(used_time)
      todo_time_string = output.time_string(todo_time)
      rec_time_string  = output.time_string(rec_time)

      logging.info('      Processed %.1f%% of records in %s (%s per record)' \
                   % (perc_done, used_time_string, rec_time_string))
      logging.info('        Estimated %s until finished' % (todo_time_string))

  # End of standardisation  - - - - - - - - - - - - - - - - - - - - - - - - - -
  #
  parallel.Barrier()  # Make sure all processes are here

  total_time = time.time() - start_time  # Calculate total time
  total_time_string = output.time_string(total_time)
  logging.info('  Total time needed for standardisation of %i records: %s' % \
               (number_records, total_time_string))

  if (parallel.size() > 1):
    comm_time_string = output.time_string(comm_time)
    logging.info('    Time for communication: %s' % (comm_time_string))

  return [total_time, comm_time]

# =============================================================================

def do_load_standard_geocode(input_dataset, output_dataset,
                             record_standardiser, record_geocoder,
                             first_record, number_records,
                             febrl_block_size):

  """The main routine that does the loading of records from the input data set
     into the output data set, if a record standardiser is given these records
     are cleaned and standardised, and then geocoded using the given geocoder.

     If no standardisation is needed the argument 'record_standardiser' has to
     be set to None.

     The output dataset must be initialised in access mode "write", "append" or
     "readwrite".

     If run in parallel, all processes will open the input data set and read
     records from it (and clean and standardised them if defined so), and then
     geocoded, before the processed records will be sent to process 0 (the host
     process) and saved into the output data set.
  """

  # Check if the given record numbers are valid - - - - - - - - - - - - - - - -
  #
  last_record = first_record + number_records

  if (first_record < 0) or (last_record > input_dataset.num_records):
    logging.exception('Record range too large: (%i,%i)' % \
                      (first_record, last_record))
    raise Exception

  # Check if the data sets are set correctly within the record standardiser - -
  #
  if (record_standardiser != None):
    if (record_standardiser.input_dataset != input_dataset):
      logging.exception('Illegal input data set definition in record '+ \
                        'standardiser: %s (should be: %s)' % \
                        (str(record_standardiser.input_dataset.name), \
                        str(input_dataset.name)))
      raise Exception
    if (record_standardiser.output_dataset != output_dataset):
      logging.exception('Illegal output data set definition in record '+ \
                        'standardiser: %s (should be: %s)' % \
                        (str(record_standardiser.output_dataset.name), \
                        str(output_dataset.name)))
      raise Exception

  else:
    # If no record standardiser for the data set is defined, the field names in
    # the input and the output data sets must be the same
    #
    input_field_name_list = input_dataset.fields.keys()
    input_field_name_list.sort()
    output_field_name_list = output_dataset.fields.keys()
    output_field_name_list.sort()

    if (input_field_name_list != output_field_name_list):
      logging.exception('Field names differ in input and output data sets ' + \
                        '(with no record standardiser defined)')
      raise Exception

  # Check if the geocoder is defined  - - - - - - - - - - - - - - - - - - - - -
  #
  if (record_geocoder == None):
    logging.exception('No record geocoder defined')
    raise Exception

  if (record_standardiser != None):
    do_string = 'Load, standardise and geocode records'
  else:
    do_string = 'Load and geocode records'

  logging.info('  %s, then write them into output data set' % (do_string))
  logging.info('    First record: %i' % (first_record))
  logging.info('    Last record:  %i' % (last_record-1))

  start_time = time.time()  # Get current time
  comm_time  = 0.0          # Communication time

  input_rec_counter = first_record  # Current record pointer

  block_cnt = 0  # A round robin block counter, used for parallelism

  # Load records in a blocked fashion - - - - - - - - - - - - - - - - - - - - -

  while (input_rec_counter < last_record):

    block_size = min(febrl_block_size, (last_record - input_rec_counter))

    # Distribute blocks equally to all processors
    #
    if ((block_cnt % parallel.size()) == parallel.rank()):

      # Load original records from input data set
      #
      input_recs = input_dataset.read_records(input_rec_counter, block_size)
      logging.info('    Loaded records %i to %i' % \
                   (input_rec_counter, input_rec_counter+block_size))

      # Standardise them if a standardiser is defined - - - - - - - - - - - - -
      #
      if (record_standardiser != None):
        clean_recs = record_standardiser.standardise_block(input_recs)
        logging.info('      Standardised records %i to %i' % \
                     (input_rec_counter, input_rec_counter+block_size))
      else:
        clean_recs = input_recs  # Take the original records directly

      # Geocode records - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      geocoded_recs = record_geocoder.match_records(clean_recs)

      # If Febrl is run in parallel, send cleaned records to process 0
      #
      if (parallel.rank() > 0):
        tmp_time = time.time()
        parallel.send(clean_recs, 0)
        comm_time += (time.time() - tmp_time)

      else:  # Process 0, store standardised records

        # Store records in output data set
        #
        output_dataset.write_records(geocoded_recs)

    # If Febrl is run in parallel, process 0 receives cleaned records
    #
    if (parallel.rank() == 0) and (block_cnt % parallel.size() != 0):

      p = (block_cnt % parallel.size())  # Process number to receive from
      tmp_time = time.time()
      tmp_recs = parallel.receive(p)
      comm_time += (time.time() - tmp_time)

      # Store records in output data set
      #
      output_dataset.write_records(tmp_recs)

    input_rec_counter += block_size  # Increment current record pointer
    block_cnt += 1

    # Now determine timing and print progress report  - - - - - - - - - - - - -
    #
    if ((block_cnt % parallel.size()) == 0):
      used_time = time.time() - start_time
      recs_done = input_rec_counter - first_record
      perc_done = 100.0 * recs_done / number_records
      rec_time  = used_time / recs_done
      todo_time = (number_records - recs_done) * rec_time

      used_time_string = output.time_string(used_time)
      todo_time_string = output.time_string(todo_time)
      rec_time_string  = output.time_string(rec_time)

      logging.info('      Processed %.1f%% of records in %s (%s per record)' \
                   % (perc_done, used_time_string, rec_time_string))
      logging.info('        Estimated %s until finished' % (todo_time_string))

  # End of standardisation  - - - - - - - - - - - - - - - - - - - - - - - - - -
  #
  parallel.Barrier()  # Make sure all processes are here

  total_time = time.time() - start_time  # Calculate total time
  total_time_string = output.time_string(total_time)
  logging.info('  Total time needed for standardisation of %i records: %s' % \
               (number_records, total_time_string))

  if (parallel.size() > 1):
    comm_time_string = output.time_string(comm_time)
    logging.info('    Time for communication: %s' % (comm_time_string))

  return [total_time, comm_time]

# =============================================================================

def do_comparison(dataset_a, dataset_b, record_comparator, classifier,
                  record_pair_dict, num_rec_pairs, febrl_block_size):
  """The main routine that does the comparison of record pairs given in the
     record pair dict and the two data sets using the given record comparator.
     The resulting weight vectors are then inserted into the given classifier.
  """

  start_time =    time.time()  # Get current time
  compare_time =  0.0          # Comparison time
  classify_time = 0.0          # Classification time

  rec_pair_cnt = 0  # Loop counter

  # Compare records, and distribute comparisons equally to all processes  - - -
  #
  for rec_num_a in record_pair_dict:

    rec_num_a_dict = record_pair_dict[rec_num_a]

    rec_a = dataset_a.read_record(int(rec_num_a))  # Read the first record

    for rec_num_b in rec_num_a_dict:

      logging.debug('          Compare records %i with %i' % \
                    (rec_num_a, rec_num_b))

      # Read the records from the data set
      #
      rec_b = dataset_b.read_record(rec_num_b)

      # Compare the two records
      #
      tmp_time = time.time()
      weight_vector = record_comparator.compare(rec_a, rec_b)
      compare_time += (time.time() - tmp_time)

      # And insert the weight vector into the classifier
      #
      tmp_time = time.time()
      classifier.classify(weight_vector)  # Classify the weight vector
      classify_time += (time.time() - tmp_time)

      rec_pair_cnt += 1

      # Now determine timing and print progress report  - - - - - - - - - - - -
      #
      if ((rec_pair_cnt % febrl_block_size) == 0):
        used_time =       time.time() - start_time
        if (num_rec_pairs > 0):
          perc_done =       100.0 * rec_pair_cnt / num_rec_pairs
        else:
          perc_done = 100.0  # No record pairs to process, so we're done
        rec_pair_time =   used_time / rec_pair_cnt
        todo_time =       (num_rec_pairs - rec_pair_cnt) * rec_pair_time
        avrg_comp_time =  (compare_time / rec_pair_cnt)
        avrg_class_time = (classify_time / rec_pair_cnt)

        used_time_string =       output.time_string(used_time)
        todo_time_string =       output.time_string(todo_time)
        rec_pair_time_string =   output.time_string(rec_pair_time)
        avrg_comp_time_string =  output.time_string(avrg_comp_time)
        avrg_class_time_string = output.time_string(avrg_class_time)

        logging.info('      Processed %.1f%% (%i/%i) of record pairs in %s' % \
                     (perc_done, rec_pair_cnt, num_rec_pairs, \
                     used_time_string) + ' (%s per record pair)' % \
                     (rec_pair_time_string))
        logging.info('        Average comparison time:     %s' % \
                     (avrg_comp_time_string))
        logging.info('        Average classification time: %s' % \
                     (avrg_class_time_string))
        logging.info('        Estimated %s until finished' % \
                     (todo_time_string))

  # Print final time for record pair comparison
  #
  total_time =    time.time() - start_time  # Calculate total time

  if (rec_pair_cnt > 0):  # Bug-fix (div by 0) by Marion Sturtevant (thanks!)
    rec_pair_time = (total_time / rec_pair_cnt)
  else:
    rec_pair_time = 0

  total_time_string =    output.time_string(total_time)
  rec_pair_time_string = output.time_string(rec_pair_time)

  logging.info('  Total time needed for comparison and classification of ' + \
               '%i record pairs: %s' % (rec_pair_cnt, total_time_string))
  logging.info('    (%s per record pair)' % (rec_pair_time_string))

  return [total_time]

# =============================================================================
# Initialise the logging system

def init_febrl_logger(log_file_name = None, file_level = None,
                      console_level = 'WARNING', clear_log = False,
                      parallel_output = 'host'):
  """Function to set up the Python logging system for both file logging and
     console (terminal) output.

     ARGUMENTS:
       log_file_name    Name of a log file, default is set to None, meaning no
                        log information is written to file.
       file_level       Set to one of the logging levels (see below) or None
                        (which is default).
       console_level    Set to one of the logging levels (see below), default
                        is 'WARNING'.
       clear_log        A flag, set to True if you want to clear (delete) the
                        old content of the log file. Set to False to keep the
                        old content.
       parallel_output  Defines how printing is done when Febrl is run in
                        parallel, can either be set to 'host' (only host
                        process is logging - except errors and exceptions), or
                        'all' (in which case all processes are logging). If
                        run in parallel each process writes into it's own log
                        file.

     DESCRIPTION:
       The following logging levels are possible:

       NOTSET, DEBUG, INFO, WARNING, ERROR, CRITICAL

       If set to a certain level only logging messages with a higher level will
       be logged.

       A current limitation is that the file log level can only be equal to or
       higher than the console level. If the console level is higher then the
       file level is set to this level as well by the Python logging system.
  """

  # Check and set the parallel printing mode  - - - - - - - - - - - - - - - - -
  #
  if (parallel_output not in ['host','all']):
    logging.exception('Illegal value for argument "parallel_output": %s' % \
                      (str(arallel_print)) + ', must beither "host" or "all"')
    raise Exception

  parallel.printmode = parallel_output

  # Define the formats for file and console logging, and date and time - - - -
  #
  file_format = parallel.prompt+'%(levelname)s [%(asctime)s; %(filename)s,' + \
                ' line %(lineno)s] %(message)s'

  #console_format = parallel.prompt+'%(levelname)s [%(filename)s, ' + \
  #                 'line %(lineno)s] %(message)s'
  console_format = parallel.prompt+'%(levelname)s [%(filename)s] %(message)s'

  datetime_format = '%d/%m/%y %H:%m'

  # A dictionary with mappings from string log levels to logging module levels
  # (including various aliases)
  #
  log_level_dict = {'NOTSET':logging.NOTSET, 'DEBUG':logging.DEBUG,
                    'INFO':logging.INFO, 'WARNING':logging.WARNING,
                    'ERROR':logging.ERROR, 'CRITICAL':logging.CRITICAL,
                    'WARN':logging.WARNING, 'ERR':logging.ERROR,
                    'FATAL':logging.CRITICAL,'ALL':logging.NOTSET}

  # Check and set console log level - - - - - - - - - - - - - - - - - - - - - -
  #
  if ((parallel_output == 'host') and (parallel.rank() > 0)):
    console_log_level = logging.ERROR  # Only log errors for parallel runs
  else:
    if (console_level in log_level_dict):
      console_log_level = log_level_dict[console_level]
    else:
      logging.warn('Illegal console logging level given: "%s"' % \
                   (console_level))
      console_log_level = logging.WARNING  # Set to default warning

  # Check and set file log level if given - - - - - - - - - - - - - - - - - - -
  #
  if (file_level != None):
    if ((parallel_output == 'host') and (parallel.rank() > 0)):
      file_log_level = logging.ERROR  # Only log errors for parallel runs
    else:
      if (file_level in log_level_dict):
        file_log_level = log_level_dict[file_level]
      else:
        logging.warn('Illegal file logging level given: "%s"' % (file_level))
        file_log_level = logging.WARNING  # Set to default warning

  # Intialise the root logger, which is also the console logger - - - - - - - -
  #
  febrl_logger = logging.getLogger()  # New logger at root level
  febrl_logger.setLevel(console_log_level)

  # Define a formater for it and set format in the logger
  #
  console_log_formatter = logging.Formatter(console_format, datetime_format)
  febrl_logger.handlers[0].setFormatter(console_log_formatter)

  print parallel.prompt + 'Initialised console logger with level "%s"' % \
        (logging.getLevelName(console_log_level))

  # Add file logger if defined - - - - - - - - - - - - - - - - - - - - - - - -
  #
  if (log_file_name != None):
    if (clear_log == True):
      file_mode = 'w'
    else:
      file_mode = 'a'

    # Modify file name, check for parallel runs and add '.log' if not given
    #
    if (parallel.size() > 1):
      if (log_file_name[-4:] in ['.log', '.LOG']):
        log_file_name = log_file_name[:-4]  # Remove log extension

      log_file_name += '_P%d' % (parallel.rank())  # Add process qualifier

    if (log_file_name[-4:] not in ['.log', '.LOG']):  # Add .log' extension
      log_file_name += '.log'

    file_log_handler =   logging.FileHandler(log_file_name, file_mode)
    file_log_formatter = logging.Formatter(file_format, datetime_format)
    file_log_handler.setFormatter(file_log_formatter)
    file_log_handler.setLevel(file_log_level)
    febrl_logger.addHandler(file_log_handler)
    print parallel.prompt + \
          'Initialised file logger with level "%s" into file "%s"' % \
          (logging.getLevelName(file_log_level), log_file_name)

# =============================================================================
# Various short functions

# -----------------------------------------------------------------------------

def check_argument_is_string(keyword, value):
  """Check if the type of the given value is a string, otherwise raise an
     exception.
  """

  if (not isinstance(value, str)):
    logging.exception('Value of argument "%s" is not a string: %s (%s)' % \
                      (keyword, str(value), type(value)))
    raise Exception

# -----------------------------------------------------------------------------

def check_argument_is_integer(keyword, value):
  """Check if the type of the given value is an integer, otherwise raise an
     exception.
  """

  if (not isinstance(value, int)):
    logging.exception('Value of argument "%s" is not an integer: %s (%s)' % \
                      (keyword, str(value), type(value)))
    raise Exception

# -----------------------------------------------------------------------------

def check_argument_is_float(keyword, value):
  """Check if the type of the given value is a floating point number,
     otherwise raise an exception.
  """

  if (not isinstance(value, float)):
    logging.exception('Value of argument "%s" is not a floating point ' % \
                      (keyword) + 'number: %s (%s)' % (str(value),type(value)))
    raise Exception

# -----------------------------------------------------------------------------

def check_argument_is_dictionary(keyword, value):
  """Check if the type of the given value is a dictionary, otherwise raise an
     exception.
  """

  if (not isinstance(value, dict)):
    logging.exception('Value of argument "%s" is not a dictionary: %s' % \
                      (keyword, type(value)))
    raise Exception

# -----------------------------------------------------------------------------

def check_argument_is_list(keyword, value):
  """Check if the type of the given value is a list, otherwise raise an
     exception.
  """

  if (not isinstance(value, list)):
    logging.exception('Value of argument "%s" is not a list: %s' % \
                      (keyword, type(value)))
    raise Exception

# -----------------------------------------------------------------------------

def check_argument_is_flag(keyword, value):
  """Check if the given value is either True or False, otherwise raise an
     exception.
  """

  if (value not in [True, False]):
    logging.exception('Value of argument "%s" is not True or False: %s' % \
                      (keyword, str(value)))
    raise Exception

# =============================================================================
