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
# The Original Software is: "geocoding.py"
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

"""Module geocoding.py - Classes for geocoding records using reference data.

   This module provides a class for geocoding of records using reference data
   based on GNAF data files that have previously been processed into inverted
   indices and stored in either binary pickled Python dictionary files, or
   Python shelve files.

   Preprocessing of the GNAF data files into cleaned and standardised inverted
   indices is done with the Febrl module 'process-gnaf.py' (which is available
   in the 'geocode/' directory)

   TODO:
   - improve processing of street numbers made of more than just a number
     (e.g. "1a")
   - same for unit numbers
   - if street mtaches and street+st. quali matches -> only use street+quali
"""

# =============================================================================
# Imports go here

import output
import febrl

import cPickle as pickle
import copy
import os
import logging
import sets
import shelve
import string
import time

# =============================================================================

class Geocoder:
  """Main class for the geocoding process. Implements routines for loading
     geocoding indices, and to do geocoding of records using exact as well as
     approximate matching, both block wise and for single records.
  """

  # ---------------------------------------------------------------------------

  def __init__(self, **kwargs):
    """Constructor.

       Initialise the geocoder and load inverted index and auxilliary files.
    """

    self.name =                     ''
    self.description =              ''

    self.geocode_file_directory = None    # Directory with GNAF index files
    self.pickle_files_extension = '.pik'  # File extension for pickle files
    self.shelve_files_extension = '.slv'  # File extension for shelve files

    self.address_site_geocode_file =     None  # Names of the three geocode
    self.street_locality_geocode_file =  None  # files
    self.locality_geocode_file =         None

    self.address_site_geocode_index =    None  # The three dictionaries with
    self.street_locality_geocode_index = None  # geocode information
    self.locality_geocode_index =        None

    self.collection_district_file =      None  # An file with CD information
    self.collection_district_index =     None

    self.input_data_set = None  # A reference to the data set to be geocoded

    # A dictionary which is used to assign fields from the (cleaned and
    # standardised) input data set to the corresponding processing fields
    # These field correspond to the Febrl address standardisation output fields
    #
    self.input_fields = {'building_name':None,
                  'location_description':None,
                         'locality_name':None,
                    'locality_qualifier':None,
                              'postcode':None,
                                 'state':None,
                          'wayfare_name':None,
                     'wayfare_qualifier':None,
                          'wayfare_type':None,
                        'wayfare_number':None,
                           'flat_number':None,
                        'flat_qualifier':None,
                             'flat_type':None,
                          'level_number':None,
                            'level_type':None,
                            'lot_number':None,
                  'lot_number_qualifier':None}

    self.match_threshold =         None
    self.best_match_only =         None
    self.maximal_neighbour_level = None
    self.max_average_address =     None  # For averaging addresses, in metres
    self.missing_value_weight =    0.0

    self.postcode_neighbours_1 = None  # Postcode level 1 region look-up table
    self.postcode_neighbours_2 = None  # Postcode level 2 region look-up table
    self.suburb_neighbours_1 =   None  # Suburb level 1 region look-up table
    self.suburb_neighbours_2 =   None  # Suburb level 2 region look-up table
    self.neighbour_score_adjust = [0.0, 0.0, 0.0]  # Score adjustments

    # Dictionary of all necessary index files (based on GNAF data files) used
    # (as created by the GNAF preprocessing module 'process-gnaf.py').
    #
    self.index_files = {'building_name':None,
                          'flat_number':None,
                   'flat_number_prefix':None,
                   'flat_number_suffix':None,
                            'flat_type':None,
                         'level_number':None,
                           'level_type':None,
                        'locality_name':None,
                       'location_descr':None,
                           'lot_number':None,
                    'lot_number_prefix':None,
                    'lot_number_suffix':None,
                         'number_first':None,
                  'number_first_prefix':None,
                  'number_first_suffix':None,
                          'number_last':None,
                   'number_last_prefix':None,
                   'number_last_suffix':None,
                             'postcode':None,
                         'state_abbrev':None,
                          'street_name':None,
                        'street_suffix':None,
                          'street_type':None}

    self.index_files_data = None  # Dictionary with GNAF inverted indices

    self.approx_index = None  # Approximate indices for selected fields

    # Some counters to collect statistics on geocode matching accuracy
    #
    self.match_stats = {'num_records':0,
                           'no match':0,
                      'exact address':0,
                    'average address':0,
                       'many address':0,
                       'exact street':0,
                     'average street':0,
                        'many street':0,
                     'exact locality':0,
                   'average locality':0,
                      'many locality':0,
             'no collection district':0}

    # Process all keyword arguments
    #
    for (keyword, value) in kwargs.items():
      if (keyword == 'name'):
        self.name = value
      elif (keyword == 'description'):
        self.description = value

      elif (keyword == 'geocode_file_directory'):
        febrl.check_argument_is_string(keyword, value)
        if (value[-1] != os.sep):  # Make sure there is a directory separator
          value = value+os.sep
        self.geocode_file_directory = value

      elif (keyword in ['pickle_files_ext', 'pickle_files_extension',
                        'pickle_file_ext', 'pickle_file_extension']):
        febrl.check_argument_is_string(keyword, value)
        self.pickle_file_extension = value

      elif (keyword in ['shelve_files_ext', 'shelve_files_extension',
                        'shelve_file_ext', 'shelve_file_extension']):
        febrl.check_argument_is_string(keyword, value)
        self.shelve_file_extension = value

      elif (keyword == 'address_site_geocode_file'):
        febrl.check_argument_is_string(keyword, value)
        self.address_site_geocode_file = value

      elif (keyword == 'street_locality_geocode_file'):
        febrl.check_argument_is_string(keyword, value)
        self.street_locality_geocode_file = value

      elif (keyword == 'locality_geocode_file'):
        febrl.check_argument_is_string(keyword, value)
        self.locality_geocode_file = value

      elif (keyword in ['collection_district_file', 'coll_dist_file']):
        febrl.check_argument_is_string(keyword, value)
        self.collection_district_file = value

      elif (keyword == 'input_data_set'):
        self.input_data_set = value

      elif (keyword == 'input_fields'):
        febrl.check_argument_is_dictionary(keyword, value)
        for (key, val) in value.items():  # Process the given dictionary
          if (key in self.input_fields):  # A known key
            self.input_fields[key] = val
          else:
            logging.exception('Illegal key word given in "input_fields" ' + \
                              'dictionary: "%s"' % (str(key)))
            raise Exception

      elif (keyword in ['match_threshold', 'threshold', 'match_threshold']):
        febrl.check_argument_is_float(keyword, value)
        self.match_threshold = value

      elif (keyword in ['best_match_only', 'best_match', 'best_only']):
        febrl.check_argument_is_flag(keyword, value)
        self.best_match_only = value

      elif (keyword in ['maximal_neighbour_level', 'max_neighbour_level']):
        febrl.check_argument_is_integer(keyword, value)
        if (value not in [0,1,2]):
          logging.exception('Argument "maximal_neighbour_level" is not one' + \
                            ' of 0,1,2')
          raise Exception
        self.maximal_neighbour_level = value

      elif (keyword in ['max_average_address', 'max_avrg_addr', \
                        'maximum_average_address']):
        febrl.check_argument_is_integer(keyword, value)
        if (value < 0):
          logging.exception('Argument "max_average_address" must be 0 or ' + \
                            'positive')
          raise Exception

        # Now convert metres into degrees
        # See: http://members.optusnet.com.au/fmet/main/degree.html
        #      http://www.gis.unbc.ca/courses/geog205/lectures/thegraticule/
        #
        # One degree of latitude is around 111 km, but one degree of longitude
        # differs from around 111 km at the equator to 0 km at the poles.
        #
        self.max_average_address = float(value) * (1.0 / 111000.0)

      elif (keyword in ['missing_value_weight', 'missing_weight', \
                        'miss_weight', 'miss_value_weight']):
        febrl.check_argument_is_float(keyword, value)
        self.missing_value_weight = value

      elif (keyword in ['postcode_neighbours_1', 'pc_neighbours_1']):
        self.postcode_neighbours_1 = value
      elif (keyword in ['postcode_neighbours_2', 'pc_neighbours_2']):
        self.postcode_neighbours_2 = value
      elif (keyword in ['suburb_neighbours_1', 'sub_neighbours_1']):
        self.suburb_neighbours_1 = value
      elif (keyword in ['suburb_neighbours_2', 'sub_neighbours_2']):
        self.suburb_neighbours_2 = value

      elif (keyword == 'index_files'):
        febrl.check_argument_is_dictionary(keyword, value)
        for (key, val) in value.items():  # Process the given dictionary
          if (key in self.index_files):  # A known key
            self.index_files[key] = val
          else:
            logging.exception('Illegal key word given in "index_files" ' + \
                  'dictionary: "%s"' % (str(key)))
            raise Exception

      elif (keyword == 'approx_index'):
        febrl.check_argument_is_dictionary(keyword, value)
        for (app_index_key, app_index) in value.items():
          if (app_index_key != app_index.field_name):  # Illegal field name
            logging.exception('Illegal field name for approximate index "%s": \
                              %s' % (app_index.name, app_index_key) + \
                              ' (must be the same as the "field_name" entry ' \
                              + ' of the index)')
            raise Exception
        self.approx_index = value

      else:
        logging.exception('Illegal constructor argument keyword: "%s"' % \
                          (keyword))
        raise Exception

    # Check if the needed attributes are set  - - - - - - - - - - - - - - - - -
    #
    if (self.name == None):
      logging.exception('Geocoder "name" not defined')
      raise Exception

    fields_defined = False  # Check if at least some fields are defined
    for (key,val) in self.input_fields.items():
      if (val != None):
        fields_defined = True
      #else:
      #  del self.input_fields[key]  # Delete the not used keys

    if (fields_defined == False):
      logging.exception('Geocoder "fields" not defined (at least one ' + \
                        'field must be different from "None")')
      raise Exception

    if (self.geocode_file_directory == None):
      logging.exception('Geocoder "geocode_file_directory" not defined')
      raise Exception
    if (self.pickle_files_extension == None):
      logging.exception('Geocoder "pickle_files_extension" not defined')
      raise Exception
    if (self.shelve_files_extension == None):
      logging.exception('Geocoder "shelve_files_extension" not defined')
      raise Exception

    if (self.address_site_geocode_file == None):
      logging.exception('Geocoder "address_site_geocode_file" not defined')
      raise Exception
    if (self.street_locality_geocode_file == None):
      logging.exception('Geocoder "street_locality_geocode_file" not defined')
      raise Exception
    if (self.locality_geocode_file == None):
      logging.exception('Geocoder "locality_geocode_file" not defined')
      raise Exception

    if (self.input_data_set == None):
      logging.exception('Geocoder "input_data_set" not defined')
      raise Exception

    # Check if the input fields definitions are of correct format (two element
    # tuples with a field name and a matching weight).
    #
    field_weight_sum = 0.0  # Get sum of field weights for normalisation later

    for (item_name, item_vals) in self.input_fields.items():
      if (item_vals != None):  # Check only if the field is defined
        if (not isinstance(item_vals, tuple)) and (len(item_vals) != 2):
          logging.exception('Input field entry "%s" is not of correct ' % \
                            (str(item_name)) + 'format: %s ' % \
                            (str(item_vals))+ '(must be a two-element tuple)')
          raise Exception
        if (item_vals[0] not in self.input_data_set.fields):
          logging.exception('Input data set "%s" does not have a field "%s"' \
                            % (str(self.input_data_set.name), \
                            str(item_vals[0])))
          raise Exception
        if ((not isinstance(item_vals[1], int)) and \
            (not isinstance(item_vals[1], float)) and (item_vals[1] <= 0.0)):
          logging.exception('Weight given for input field "%s" is not a ' % \
                (str(item_name)) + 'positive number: %s' % (str(item_vals[1])))
          raise Exception
        field_weight_sum += item_vals[1]

    for item_name in self.input_fields:  # Normalise field weights
      item_vals = self.input_fields[item_name]
      if (item_vals != None):
        item_vals = (item_vals[0], item_vals[1] / field_weight_sum)
        self.input_fields[item_name] = item_vals

    if (self.match_threshold == None):
      logging.exception('Geocoder "match_threshold" not defined')
      raise Exception

    elif ((self.match_threshold < 0.0) or (self.match_threshold > 1.0)):
      logging.exception('Illegal value of "match_threshold" (not in [0,1])')
      raise Exception

    if (self.best_match_only == None):
      logging.exception('Geocoder "best_match_only" not defined')
      raise Exception

    if (self.maximal_neighbour_level == None):
      logging.exception('Geocoder "maximal_neighbour_level" not defined')
      raise Exception

    if (self.max_average_address == None):
      logging.exception('Geocoder "max_average_address" not defined')
      raise Exception

    # Check references to neighbourhood look-up tables only if level > 0
    #
    if (self.maximal_neighbour_level > 0):
      if (self.postcode_neighbours_1 == None):
        logging.exception('Geocoder level 1 postcode look-up table not ' + \
                          'defined')
        raise Exception
      if (self.suburb_neighbours_1 == None):
        logging.exception('Geocoder level 1 suburb look-up table not defined')
        raise Exception
    if (self.maximal_neighbour_level > 1):  # Level 2 only if maximum larger 1
      if (self.postcode_neighbours_2 == None):
        logging.exception('Geocoder level 2 postcode look-up table not ' + \
                          'defined')
        raise Exception
      if (self.suburb_neighbours_2 == None):
        logging.exception('Geocoder level 2 suburb look-up table not defined')
        raise Exception

    # Check if the index files definitions are of correct format (two element
    # tuples with a file name and a 'pickle' or 'shelve' type).
    #
    for (item_name, item_vals) in self.index_files.items():
      if (item_vals != None):  # Check only if the field is defined
        if (not isinstance(item_vals, tuple)):
          logging.exception('Index file entry "%s" is not of correct ' % \
                            (str(item_name)) + 'format: %s ' % \
                            (str(item_vals))+ '(must be a two-element tuple)')
          raise Exception
        if (not isinstance(item_vals[0], str)):
          logging.exception('Index file entry "%s" is not of correct ' % \
                            (str(item_name)) + 'format: %s ' % \
                            (str(item_vals))+ '(first entry must be a ' + \
                            'file name)')
          raise Exception
        if (item_vals[1] not in ['pickle', 'shelve']):
          logging.exception('Index file entry "%s" is not of correct ' % \
                            (str(item_name)) + 'format: %s ' % \
                            (str(item_vals))+ '(second entry must be ' + \
                            'either "pickle" or "shelve")')
          raise Exception

    # Calculate a match score adjustment for neighbouring levels - - - - - - -
    # (interpolate beween 1.0 and match_threshold)
    #
    match_threshold_interval = (1.0 - self.match_threshold) / 3.0
    self.neighbour_score_adjust = [1.0, 1.0 - match_threshold_interval,
                                   1.0 - 2*match_threshold_interval]

    # A log message - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    #
    logging.info('Initialised geocoder')
    logging.info('  Name:                         %s' % (str(self.name)))
    logging.info('  Geocode file directory:       %s' % \
                 (self.geocode_file_directory))
    logging.info('  Pickle files extension:       %s' % \
                 (self.pickle_files_extension))
    logging.info('  Shelve files extension:       %s' % \
                 (self.shelve_files_extension))
    logging.info('  Address site geocode file:    %s' % \
                 (self.address_site_geocode_file))
    logging.info('  Street-locality geocode file: %s' % \
                 (self.street_locality_geocode_file))
    logging.info('  Locality geocode file:        %s' % \
                 (self.locality_geocode_file))

    if (self.collection_district_file != None):
      logging.info('  Collection district file:     %s' % \
                   (self.collection_district_file))
    else:
      logging.info('  No collection district file defined')

    logging.info('  Input data set:               %s' % \
                 (self.input_data_set.name))
    logging.info('  Input fields (and their matching weights):')
    for (field_name, field_vals) in self.input_fields.items():
      if (field_vals != None):
        logging.info('    %20s: %20s / %.3f' % \
                     (field_name, field_vals[0], field_vals[1]))
      else:
        logging.info('    %20s: %20s' % (field_name, 'Not defined'))

    logging.info('  Missing value weight:             %.3f' % \
                 (self.missing_value_weight))
    logging.info('  Maximal neighbour level:          %d' % \
                 (self.maximal_neighbour_level))
    logging.info('  Neighbour score adjustments:      %s' % \
                 (str(self.neighbour_score_adjust)))
    logging.info('  Match threshold:                  %.3f' % \
                 (self.match_threshold))
    logging.info('  Best match only:                  %s' % \
                 (str(self.best_match_only)))
    logging.info('  Maximum average address distance: %.8f' % \
                 (self.max_average_address))

    if (self.postcode_neighbours_1 != None):
      logging.info('  Postcode neighbours level 1 look-up table: %s' % \
                   (str(self.postcode_neighbours_1.name)))
    if (self.postcode_neighbours_2 != None):
      logging.info('  Postcode neighbours level 2 look-up table: %s' % \
                   (str(self.postcode_neighbours_2.name)))
    if (self.suburb_neighbours_1 != None):
      logging.info('  Suburb neighbours level 1 look-up table:   %s' % \
                   (str(self.suburb_neighbours_1.name)))
    if (self.suburb_neighbours_2 != None):
      logging.info('  Suburb neighbours level 2 look-up table:   %s' % \
                   (str(self.suburb_neighbours_2.name)))

    logging.info('  Index files (and their file type):')
    for (field_name, field_vals) in self.index_files.items():
      if (field_vals != None):
        logging.info('    %20s: %20s / %s' % \
                     (field_name, field_vals[0], field_vals[1]))
      else:
        logging.info('    %20s: Not defined' % (field_name))

    if (self.approx_index != None):  # Approximate indices are defined
      field_str = ''
      for field_name in self.approx_index.keys():
        field_str += field_name+', '
      field_str = field_str[:-2]
      logging.info('  Approximate indices defined for fields: %s' % \
                   (field_str))

    # Open the geocode reference pickles or shelves - - - - - - - - - - - - - -
    #
    logging.info('Load address site geocoding index')
    self.address_site_geocode_index = \
       self.load_shelve_pickle_file(self.address_site_geocode_file)
    self.street_locality_geocode_index = \
       self.load_shelve_pickle_file(self.street_locality_geocode_file)
    self.locality_geocode_index = \
       self.load_shelve_pickle_file(self.locality_geocode_file)

    # If defined open the collection district pickle or shelve - - - - - - - -
    #
    if (self.collection_district_file != None):
      logging.info('Load collection district index')
      self.collection_district_index = \
        self.load_shelve_pickle_file(self.collection_district_file)
    else:
      self.collection_district_index = {}  # For fast look-up

    # Load the index files - - - - - - - - - - - - - - - - - - - - - - - - - -
    #
    logging.info('Load geocoding index files')

    # A dictionary which will contain the inverted indices (a dictionary with
    # sets for each entry)
    #
    self.index_files_data = {}

    for (item_name, item_vals) in self.index_files.items():
      if (item_vals != None):
        logging.info('  Open index file for "%s":' % (item_name))

        if (item_vals[1] == 'pickle'):  # Pickle file type
          index_file_name = self.geocode_file_directory + item_vals[0] + \
                            self.pickle_files_extension
          logging.info('    Load index from pickle file "%s"' % \
                       (index_file_name))

          index_pickle_file = open(index_file_name)
          index_pickle_dict = pickle.load(index_pickle_file)
          index_pickle_file.close()

          self.index_files_data[item_name] = index_pickle_dict

        elif (item_vals[1] == 'shelve'):
          index_file_name = self.geocode_file_directory + item_vals[0] + \
                            self.shelve_files_extension
          logging.info('    Load index from shelve file "%s"' % \
                       (index_file_name))

          index_shelve_dict = shelve.open(index_file_name, flag='r')

          self.index_files_data[item_name] = index_shelve_dict

        #logging.info('    Index dictionary contains %d entries' % \
        #             (len(self.index_files_data[item_name].keys())))

    # Build the approximate indices - - - - - - - - - - - - - - - - - - - - - -
    #
    if (self.approx_index != None):

      for field_name in self.approx_index:
        field_data = self.index_files_data[field_name]
        self.approx_index[field_name].build_index(field_data)

    logging.info('Geocoder "%s" initialised on %s' % \
                 (self.name, time.ctime(time.time())))

  # ---------------------------------------------------------------------------

  def match_records(self, input_record_list):
    """Wrapper method to match a list of records. Simply calls 'match_record'
       for each input record given in the input record list.

       The following fields/attributes are added to each of the input records:
       - latitude
       - longitude
       - match_status
       - match_weight (numerical weight - the larger the better the match)
       - gnaf_pid (one or more G-NAF PIDs separated by whitespaces, for street
                   level matches combined STREET_PID/LOCALITY_PID are used,
                   separated with a '/')
       - collection_district

       If no match has been found all these fields will be empty except the
       match_status field which will be set to: 'No match'.

       Returns the list of output records (basically the input records with
       the above given fields attached).
    """

    start_time = time.time()

    if (not isinstance(input_record_list, list)):
      logging.exception('Input record list is not a list, but: %s' % \
                        (str(type(input_record_list))))
      raise Exception

    num_input_records = len(input_record_list)
    output_record_list = []

    for i in range(num_input_records):

      self.match_stats['num_records'] = self.match_stats['num_records'] + 1

      match_result_list = self.match_record(input_record_list[i])

      out_rec = copy.copy(input_record_list[i])

      if (len(match_result_list) == 1):
        out_rec['match_score'] =         str(match_result_list[0][0])
        out_rec['gnaf_pid'] =            str(match_result_list[0][1])
        out_rec['match_status'] =        str(match_result_list[0][3])
        out_rec['latitude'] =            str(match_result_list[0][4])
        out_rec['longitude'] =           str(match_result_list[0][5])
        out_rec['collection_district'] = str(match_result_list[0][6])
        out_rec['max_avrg_distance'] =   str(match_result_list[0][7])
        out_rec['neighbour_level'] =     str(match_result_list[0][8])

        if (match_result_list[0][6] == ''):
          self.match_stats['no collection district'] = \
            self.match_stats['no collection district'] + 1

#      else:  # How to handle this...? **************************************
#        pass

      output_record_list.append(out_rec)  # Add to output record list

    used_time = time.time() - start_time

    logging.warn('============================================')
    logging.warn('Statistics for geocode matching %d records:' % \
                 (self.match_stats['num_records']))
    for (match_type, stats_count) in self.match_stats.items():
      if (match_type != 'num_records'):
        logging.warn('   %20s: %d' % (match_type, stats_count))

    logging.warn('Time used: %s (for last %d records)' % \
                 (output.time_string(used_time), num_input_records))
    logging.warn('============================================')

    return output_record_list

  # ---------------------------------------------------------------------------

  def match_record(self, input_record):
    """Match input record with geocode reference data.

       This method matches the set input record (in self.in_rec) with G-NAF
       geocoded reference data. It sets the corresponding fields in
       self.out_rec (possibly with corrected field values), including the
       match status (describing what has been matched and how), the found
       longitude and latitude, the G-NAF PID (if a single match has been
       found, or a list of PIDs separated by whitespaces), the collection
       district, the maximum distance when average matches were found, and the
       neighbouring level used for matching, and

       If no match could be found this method returns False, otherwise True.

       If more than one match has been found, then for address level matches
       an average location is calculated and returned, while for multiple
       street and locality level matches no location is returned.

       The different components (postcode, locality, wayfare, wayfare name,
       unit, and property) of the input record are processed independently.
    """

    self.in_rec =  input_record  # The standardised input record

    if (not isinstance(self.in_rec, dict)) or (self.in_rec == {}):
      logging.warn('Input record is not a dictionary or is empty: "%s"' % \
                   (str(self.in_rec)))
      return False

    logging.info('----------------------------------------------------')

    logging.info('Geocode standardised input record:')
    for (field_name, field_val) in self.in_rec.items():
      logging.info('        %s: "%s"' % (field_name, field_val))

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Step 1: Find matches and their intersections

    # Get address (wayfare number) and street (wayfare) level matches - - - - -
    #
    addr_match_list = self.get_address_matches('exact')  # Currently only exact

    street_match_list = self.get_street_matches('exact')
    street_match_type = 'exact'
    if (street_match_list == []):
      street_match_type = 'approx'
      street_match_list = self.get_street_matches('approx')

    # Determine the neighbouring levels to be searched - - - - - - - - - - - -
    #
    if ((addr_match_list == []) and (street_match_list == [])):
      nb_level_list = [0]  # No neighbour search if no address and no street
    else:
      nb_level_list = range(self.maximal_neighbour_level+1)

    first_loc_match_list = None
    match_flag_list = []  # Will become a list of tuples with three flags

    # Loop over neighbouring levels - - - - - - - - - - - - - - - - - - - - - -
    #
    for nb_level in nb_level_list:

      # Get the locality level matches
      #
      loc_match_list = self.get_locality_matches(nb_level)

      # Save the first non-empty locality match list
      #
      if ((first_loc_match_list == None) and (loc_match_list != [])):
        first_loc_match_list = loc_match_list

      # Check for matches between locality and street level - - - - - - - - - -
      #
      if ((loc_match_list != []) and (street_match_list != [])):

        red_street_match_list = self.reduce_locality_matches(loc_match_list,
                                                             street_match_list)
        if (red_street_match_list != []):
          street_loc_match_flag = True

        elif (red_street_match_list == []) and (street_match_type == 'exact'):

          # So far only exact street matches, try approximate matching
          #
          street_approx_match_list = self.get_street_matches('approx')
          red_street_approx_match_list = \
            self.reduce_locality_matches(loc_match_list,
                                         street_approx_match_list)
          if (red_street_approx_match_list != []):
            street_loc_match_flag = True
            red_street_match_list = red_street_approx_match_list
            street_match_list = street_approx_match_list
            street_match_type = 'approx'
          else:
            red_street_match_list = street_match_list  # Keep unreduced matches
            street_loc_match_flag = False
        else:
          red_street_match_list = street_match_list  # Keep unreduced matches
          street_loc_match_flag = False
      else:
        street_loc_match_flag = False
        red_street_match_list = street_match_list  # Keep unreduced matches

      logging.info('street-loc matches: %s' % (str(street_loc_match_flag)))

      # Check for matches between locality and address level - - - - - - - - -
      #
      if ((loc_match_list != []) and (addr_match_list != [])):

        red_addr_match_list = self.reduce_locality_matches(loc_match_list,
                                                           addr_match_list)
        if (red_addr_match_list != []):
          addr_loc_match_flag = True
        else:
          red_addr_match_list = addr_match_list  # Keep unreduced matches
          addr_loc_match_flag = False
      else:
        addr_loc_match_flag = False
        red_addr_match_list = addr_match_list  # Keep unreduced matches

      logging.info('address-loc matches: %s' % (str(addr_loc_match_flag)))

      # Check for matches between street and address level - - - - - - - - - -
      #
      if ((red_addr_match_list != []) and (red_street_match_list != [])):

        re2_addr_match_list = self.reduce_street_matches(red_street_match_list,
                                                         red_addr_match_list)
        if (re2_addr_match_list != []):
          addr_street_match_flag = True
          red_addr_match_list = re2_addr_match_list
        else:
          addr_street_match_flag = False

      else:
        addr_street_match_flag = False

      logging.info('address-street matches: %s' %(str(addr_street_match_flag)))

      match_flag_list.append((street_loc_match_flag,addr_loc_match_flag,
                              addr_street_match_flag))

      logging.info(self.match_list_summary_string(addr_match_list, 'address'))
      logging.info(self.match_list_summary_string(red_addr_match_list,
                                                  'reduced address'))
      logging.info(self.match_list_summary_string(street_match_list, 'street'))
      logging.info(self.match_list_summary_string(red_street_match_list,
                                                  'reduced street'))
      logging.info(self.match_list_summary_string(loc_match_list, 'locality'))

      # Check if we're at neighbour level 1 or 2 and no changes in matches - -
      #
      if ((nb_level == 1) and (len(nb_level_list) == 2) and \
          (sum(match_flag_list[0]) >= sum(match_flag_list[1]))):

        # No improvement in matches
        #
        loc_match_list = first_loc_match_list  # Return nb_level 0 match lists
        if (loc_match_list == None):
          loc_match_list = []
        logging.warn('** use first locality match list')

      elif (nb_level == 2) and ((match_flag_list[2] != (True,True,True)) or \
                                (street_match_type == 'approx')):

        # Only return level 2 matches if all three flags are True and exact
        # street matches were used
        #
        red_street_match_list = street_match_list  # Original street match list
        loc_match_list = first_loc_match_list  # Return nb_level 0 match lists
        if (loc_match_list == None):
          loc_match_list = []
        logging.warn('** use first locality match list')

      # Check if intersection of matches has been successful - - - - - - - -
      #
      if (addr_street_match_flag == True) and (street_loc_match_flag == True) \
         and (addr_loc_match_flag == True):
        break  # Matches between all three levels: Exit loop

      if ((addr_match_list == []) and  (street_loc_match_flag == True)):
        break  # No address and matches between street and locality: Exit loop

      if ((loc_match_list == []) and  (addr_street_match_flag == True)):
        break  # No locality and matches between address and street: Exit loop

      if ((street_match_list == []) and  (addr_loc_match_flag == True)):
        break  # No street and matches between address and locality: Exit loop

    addr_match_list =   red_addr_match_list
    street_match_list = red_street_match_list

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Step 2: Try to refine the found address matches with unit and property

    logging.info('  Try to refine matches') # ********************

    if ((addr_match_list != []) and \
        (addr_street_match_flag == True) or (addr_loc_match_flag == True)):
      addr_match_list = self.get_address_refinement('exact',addr_match_list)

      logging.info(self.match_list_summary_string(addr_match_list,
                                                  'refined address'))

#    logging.info('address match list: %s' % (str(addr_match_list)))
#    logging.info('street match list: %s' % (str(street_match_list)))
#    logging.info('locality match list: %s' % (str(loc_match_list)))

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Step 3: Combine the matches and retrieve coordinates

    (final_match_dict, final_match_level) = \
      self.combine_matches(addr_match_list, street_match_list, loc_match_list)

    if (final_match_level == 'address'):
#      final_match_coordinates = self.get_address_coordinates(final_match_dict)
      final_match_coordinates = self.get_coordinates(final_match_dict,
                                                     'address')

      if (final_match_coordinates == []):  # No coordinates could be found

        # Redo matches with street and locality level matches only
        #
        (final_match_dict, final_match_level) = \
          self.combine_matches([], street_match_list, loc_match_list)

    if (final_match_level == 'street'):
#      final_match_coordinates = self.get_street_coordinates(final_match_dict)
      final_match_coordinates = self.get_coordinates(final_match_dict,
                                                     'street')

      if (final_match_coordinates == []):  # No coordinates could be found

        # Redo matches withlocality level matches only
        #
        (final_match_dict, final_match_level) = \
          self.combine_matches([], [], loc_match_list)
        logging.warn('*0*: %s / %s' % \
                     (final_match_level, str(final_match_dict)))

    if (final_match_level == 'locality'):

      if (nb_level > 0):  # Try the first non-empty locality match list first

        (first_final_match_dict, final_match_level) = \
          self.combine_matches([], [], first_loc_match_list)
        logging.warn('*1*: %s / %s' % \
                     (final_match_level, str(first_final_match_dict)))
#        final_match_coordinates = \
#          self.get_locality_coordinates(first_final_match_dict)
        final_match_coordinates = self.get_coordinates(first_final_match_dict,
                                                       'locality')

        if (final_match_coordinates == []):  # No coordinates could be found
          logging.warn('*2*: %s / %s' % \
                       (final_match_level, str(final_match_dict)))

#          final_match_coordinates = \
#            self.get_locality_coordinates(final_match_dict)
          final_match_coordinates = self.get_coordinates(final_match_dict,
                                                         'locality')
      else:
        final_match_coordinates = self.get_coordinates(final_match_dict,
                                                       'locality')

    if (final_match_level == 'no match') or (final_match_coordinates == []):
      logging.warn('  No final matches for input "%s"' % \
                   (self.in_rec['input_record']))
      final_match_coordinates = [(0.0, '', {}, 'no match', '', '', '', '0.0')]

    # Add the neighbouring level to all matches - - - - - - - - - - - - - - - -
    #
    i = 0
    for i in range(len(final_match_coordinates)):
      match_tuple = final_match_coordinates[i]
      match_tuple = tuple(list(match_tuple)+[str(nb_level)])
      final_match_coordinates[i] = match_tuple

    # Use the best match for the statistics - - - - - - - - - - - - - - - - - -
    #
    match_status = final_match_coordinates[0][3]
    self.match_stats[match_status] = self.match_stats[match_status] + 1

    logging.info('Number of matches found: %d' % \
                 (len(final_match_coordinates)))
##
    if (self.best_match_only == True):
      logging.info('    (only best match is returned)')
      final_match = [final_match_coordinates[0]]  # First list element only
    else:
      final_match = final_match_coordinates  # Return all matches (sorted)

    for match in final_match:
      logging.info('  Match score:               %f' % (match[0]))
      logging.info('  Match PIDs:                %s' % (match[1]))
      logging.info('  Match values:              %s' % (str(match[2])))
      logging.info('  Match status:              %s' % (match[3]))
      logging.info('  Match latitude:            %s' % (match[4]))
      logging.info('  Match longitude:           %s' % (match[5]))
      logging.info('  Match coll. district:      %s' % (match[6]))
      logging.info('  Match max. avrg. distance: %s' % (match[7]))
      logging.info('  Match neighbour level:     %s' % (match[8]))
      logging.info('')

    return final_match

  # ---------------------------------------------------------------------------

  def get_coordinates(self, match_dict, match_level):
    """This method finds coordinates for the given level (either 'address',
       'street', or 'locality') for the PIDs in the given amtch dictionary.

       It returns a list of tuples containing eight elements:

         (match_score, match_pids, match_descr, match_status, latitude,
          longitude, coll_district)

       with:
       - match_score       The sum of all contributing match scores.
       - match_pids        A string with one or more matched PIDs (either
                           ADDRESS_SITE_PIDs, pairs of
                           STREET_PID/LOCALITY_PIDs, or LOCALITY_PIDs) for
                           which coordinates were found (separated by a space
                           if more than one).
       - match_descr       A string of a dictionary with up to three keys
                           'address', 'street' and 'locality' and the
                           corresponding values that were matched.
       - match_status      A string which gives the match status, either
                           'exact', 'many', 'average' or 'no' (the match level
                           is part of the match status string, e.g. 'exact
                           address').
       - latitude          (as a string)
       - longitude         (as a string)
       - coll_district     (as a string)
       - max_avrg_distance  A numerical value (in metres) if the returned
                            matches were averaged, zero otherwise.
    """

    match_result_list = []  # Final match results

    match_score_list = match_dict.keys()
    match_score_list.sort()
    match_score_list.reverse()

    for match_score in match_score_list:

      match_list = match_dict[match_score] # Match(es) with this score

      score_result_list = []  # Results for this score

      for (pid_set, match_descr) in match_list:

        # Extract geocoding information for all matching records - - - - - - -
        #
        latitude_set =  sets.Set()
        longitude_set = sets.Set()
        match_pid_set = sets.Set()  # PIDs for which coordinates were found

        latitude_sum =  0.0
        longitude_sum = 0.0
        num_coordinates = 0

        # Process all PIDs for this match - - - - - - - - - - - - - - - - - - -
        #
        for pid in pid_set:

          # Get the coordinates from the corresponding geocode index - - - - -
          #
          if (match_level == 'address'):
            location = self.address_site_geocode_index.get(pid, None)

          elif (match_level == 'street'):
            (loc_pid, street_pid) = pid  # Street uses two PIDs
            loc_pid_dict = self.street_locality_geocode_index.get(loc_pid, {})
            location = loc_pid_dict.get(street_pid, None)

          elif (match_level == 'locality'):
            location = self.locality_geocode_index.get(pid, None)

          else:
            logging.exception('Illegal match level: %s' % (str(match_level)))
            raise Exception

          if (location != None):
            tmp_latitude, tmp_longitude = location.split(',')

            latitude_set.add(tmp_latitude)
            longitude_set.add(tmp_longitude)
            latitude_sum    += float(tmp_latitude)
            longitude_sum   += float(tmp_longitude)
            num_coordinates += 1

            if (match_level == 'street'):  # Street uses two PIDs
              match_pid_set.add(str(street_pid+'/'+loc_pid))
            else:
              match_pid_set.add(pid)

        # Get the collection district for these matches - - - - - - - - - - - -
        #
        coll_dist = self.get_collection_district(match_pid_set)

        # Check if there are coordinates - - - - - - - - - - - - - - - - - - -
        #
        if (num_coordinates == 0):
          logging.warn('No coordinates for %s level PIDs: %s' % \
                       (match_level, str(pid_set)))

        # One location only: Exact match  - - - - - - - - - - - - - - - - - - -
        #
        elif ((len(latitude_set) == 1) and (len(longitude_set) == 1)):

          score_result_list.append((match_score, match_pid_set,
                                    str(match_descr), 'exact '+match_level,
                                    latitude_set.pop(), longitude_set.pop(),
                                    coll_dist, '0.0'))

        else:  # Many different coordinates: Average or many match  - - - - - -

          # Find minimum and maximum latitudes and longitudes
          #
          lati_list = list(latitude_set)
          lati_list.sort()
          long_list = list(longitude_set)
          long_list.sort()

          maximum_dist = max((float(long_list[-1])-float(long_list[0])),
                             (float(lati_list[-1])-float(lati_list[0])))
          logging.info('(1)  Maximum distance: %f (%.1d metres)' % \
                       (maximum_dist, (maximum_dist*111000.0)))

          # Calculate average if smaller than allowed maximum - - - - - - - - -
          #
          if (maximum_dist <= self.max_average_address):

            latitude =  str(latitude_sum / float(num_coordinates))
            longitude = str(longitude_sum / float(num_coordinates))

            score_result_list.append((match_score, match_pid_set,
                                      str(match_descr), 'average '+match_level,
                                      latitude, longitude, coll_dist,
                                      str(maximum_dist*111000.0)))

          # If distances too large we have a 'many' match - - - - - - - - - - -
          #
          else:  # Directly append to match result lsit

            match_result_list += [(match_score, ' '.join(match_pid_set),
                                   str(match_descr), 'many '+match_level, '',
                                   '', coll_dist, str(maximum_dist*111000.0))]

      logging.info('** score result list: %s' % (str(score_result_list)))

      # Post-process the match list for this score - - - - - - - - - - - - - -
      #
      if (len(score_result_list) == 1):  # One match only
        match_tuple_list = list(score_result_list[0])
        match_tuple_list[1] = ' '.join(match_tuple_list[1])  # Set to string
        match_result_list += [tuple(match_tuple_list)]

      elif (len(score_result_list) > 1):  # Several matches

        # Check if the coordinates are the same, if so we have an exact match
        #
        match_pid_set =   sets.Set()
        match_descr_set = sets.Set()
        latitude_set =    sets.Set()
        longitude_set =   sets.Set()
        coll_dist_set =   sets.Set()

        latitude_sum =  0.0
        longitude_sum = 0.0
        num_coordinates = 0

        for match in score_result_list:
          match_pid_set = match_pid_set.union(match[1])
          match_descr_set.add(str(match[2]))
          latitude_set.add(match[4])
          longitude_set.add(match[5])
          coll_dist_set.add(match[6])

          latitude_sum    += float(match[4])
          longitude_sum   += float(match[5])
          num_coordinates += 1

        # One coordinate pair only, exact match
        #
        if ((len(latitude_set) == 1) and (len(longitude_set) == 1)):
          match_result_list += [(match_score, ' '.join(match_pid_set),
                                ' '.join(match_descr_set),
                                'exact '+match_level, latitude_set.pop(),
                                longitude_set.pop(), ' '.join(coll_dist_set),
                                '0.0')]

        else:  # Different coordinates, average or many match

          # Find minimum and maximum latitudes and longitudes
          #
          lati_list = list(latitude_set)
          lati_list.sort()
          long_list = list(longitude_set)
          long_list.sort()

          maximum_dist = max((float(long_list[-1])-float(long_list[0])),
                             (float(lati_list[-1])-float(lati_list[0])))
          logging.info('(2)  Maximum distance: %f (%.1d metres)' % \
                       (maximum_dist, (maximum_dist *111000.0)))

          # Calculate average if smaller than allowed maximum - - - - - - - - -
          #
          if (maximum_dist <= self.max_average_address):

            latitude =  str(latitude_sum / float(num_coordinates))
            longitude = str(longitude_sum / float(num_coordinates))

            match_result_list += [(match_score, ' '.join(match_pid_set),
                                ' '.join(match_descr_set),
                                'average '+match_level, latitude, longitude,
                                ' '.join(coll_dist_set),
                                str(maximum_dist*111000.0))]

          # If distances too large we have a 'many' match - - - - - - - - - - -
          #
          else:  # Directly append to match result list

            match_result_list += [(match_score, ' '.join(match_pid_set),
                                   str(match_descr), 'many '+match_level, '',
                                   '', coll_dist, str(maximum_dist*111000.0))]

      logging.info('match result list: %s' % (str(match_result_list))) # *****

      # Check if there were matches and if only best matches are needed - - - -
      #
      if ((match_result_list != []) and (self.best_match_only == True)):
        break  # Exit loop and return

    return match_result_list

  # ---------------------------------------------------------------------------

  def get_collection_district(self, pid_set):
    """Find the collection districts for the given PID set, and if they are all
       in the same collection district then return it (as a string), otherwise
       return an empty string.
    """

    # Make sure that the pid_set is a set (and not just one PID)
    #
    if (isinstance(pid_set, str)):
      pid_set = sets.Set([pid_set])

    coll_dist_set = sets.Set()
    for pid in pid_set:
      coll_dist_set.add(self.collection_district_index.get(pid, ''))
    coll_dist_set.discard('')  # Remove empty CD values

    if (len(coll_dist_set) == 1):  # One CD only
      coll_dist = str(coll_dist_set.pop())
    else:
      coll_dist = ''  # More than on CD

    return coll_dist

  # ---------------------------------------------------------------------------

  def reduce_locality_matches(self, locality_match_list, other_match_list):
    """A method that reduces the dictionary of LOCALITY_PIDs in the given
       other (either address or street) match list by only keeping the
       LOCALITY_PIDs which are in the given locality match list.

       The method returns the reduced other match list.
    """

    if (other_match_list == []):
      return other_match_list

    # Build union of all LOCALITY_PIDs in the locality match list - - - - - - -
    #
    loc_pid_union = sets.Set()
    for (loc_val, loc_match_score, loc_pid_set) in locality_match_list:
      loc_pid_union = loc_pid_union.union(loc_pid_set)

    logging.debug('  Union of all LOCALITY_PIDs has length %d' % \
                  (len(loc_pid_union)))

    # Only keep LOCALITY_PIDS in other matches which are in the union set - - -
    #
    red_other_match_list = []

    for (other_val, other_match_score, other_loc_pid_dict) in other_match_list:

      red_other_loc_dict = {}

      for loc_pid in loc_pid_union:

        if (loc_pid in other_loc_pid_dict):
          red_other_loc_dict[loc_pid] = other_loc_pid_dict[loc_pid]

      if (red_other_loc_dict != {}):  # Only keep if not empty
        red_other_match_list.append((other_val, other_match_score,
                                     red_other_loc_dict))

    logging.debug('  Reduced other match list from %d to %d entries' % \
                  (len(other_match_list), len(red_other_match_list)))

    return red_other_match_list

  # ---------------------------------------------------------------------------

  def reduce_street_matches(self, street_match_list, address_match_list):
    """A method that reduces the dictionaries of STREET_PIDs in the given
       address match list by only keeping the STREET_PIDs which are in the
       given street match list.

       The method returns the reduced address match list.
    """

    if (address_match_list == []):
      return address_match_list

    # Build dictionary of all LOCALITY_PIDs with union of all their STREET_PIDs
    # from the street match list.
    #
    street_loc_pid_dict = {}

    for (street_val, street_match_score, street_dict) in street_match_list:
      for loc_pid in street_dict:
        street_pid_set = street_loc_pid_dict.get(loc_pid, sets.Set())
        street_pid_set = street_pid_set.union(street_dict[loc_pid])

        street_loc_pid_dict[loc_pid] = street_pid_set

    logging.debug('  Number of different LOCALITY_PIDs from street match ' + \
                  'set: %d' % (len(street_loc_pid_dict)))

    # Only keep STREET_PIDs in address matches which are in union dictionary -
    #
    red_address_match_list = []

    for (addr_val, addr_match_score, addr_dict) in address_match_list:

      red_addr_dict = {}

      for loc_pid in street_loc_pid_dict:
        if (loc_pid in addr_dict):

          addr_street_dict = addr_dict[loc_pid]
          street_pid_set = street_loc_pid_dict[loc_pid]
          new_addr_street_dict = {}

          for street_pid in street_pid_set:
            if (street_pid in addr_street_dict):
              new_addr_street_dict[street_pid] = addr_street_dict[street_pid]

          if (new_addr_street_dict != {}):  # Only keep if not empty
            red_addr_dict[loc_pid] = new_addr_street_dict


      if (red_addr_dict != {}):  # Only keep if not empty
        red_address_match_list.append((addr_val, addr_match_score,
                                       red_addr_dict))

    logging.debug('  Reduced address match list from %d to %d entries' % \
                  (len(address_match_list), len(red_address_match_list)))

    return red_address_match_list

  # ---------------------------------------------------------------------------

  def combine_matches(self, addr_match_list, street_match_list,
                      loc_match_list):
    """Method which combines the three given match lists and returns a
       dictionary of matches, with keys being the match scores (the sum of all
       contributing match scores) and values being a list of tuples, each
       containing two elements:

         (match_pid_set, match_descr)

       with
       - match_pid_set  A set of one or more PIDs of the match (either
                        ADDRESS_SITE_PIDs, STREET_PIDs, or LOCALITY_PIDs,
                        according to the match level)
       - match_descr    A dictionary with up to three keys 'address', 'street'
                        and 'locality' an the corresponding values that were
                        matched.

       Also returned is a string 'match_level' (which will be either 'address',
       'street', 'locality', or 'no match') which indicates the level of the
       achieved matches.

       The list of matches returned is sorted (in descending order) according
       to the match scores.

       If no matches could be found an empty list will be returned.
    """

    final_match_dict = {}

    # Start with the address level - - - - - - - - - - - - - - - - - - - - - -
    #
    if (addr_match_list != []):

      # Case 1: There are also street and locality matches
      #
      if ((street_match_list != []) and (loc_match_list != [])):

        for (loc_val, loc_score, loc_pid_set) in loc_match_list:

          for loc_pid in loc_pid_set:

            for (street_val, street_score, street_loc_pid_dict) in \
              street_match_list:

              # Is LOCALITY_PID common in street and locality level
              #
              if (loc_pid in street_loc_pid_dict):

                street_street_pid_set = street_loc_pid_dict[loc_pid]

                for street_pid in street_street_pid_set:

                  for (addr_val, addr_score, addr_loc_pid_dict) in \
                    addr_match_list:

                    # Is LOCALITY_PID common in address and locality level
                    #
                    if (loc_pid in addr_loc_pid_dict):

                      addr_street_pid_dict = addr_loc_pid_dict[loc_pid]

                      # Is STREET_PID common in address and street level
                      #
                      if (street_pid in addr_street_pid_dict):

                        match_score =   addr_score+street_score+loc_score
                        match_pid_set = addr_street_pid_dict[street_pid]
                        match_descr = {'address':addr_val, 'street':street_val,
                                       'locality':loc_val}
                        match_tuple = (match_pid_set, match_descr)

                        match_list = final_match_dict.get(match_score, [])
                        if match_tuple in match_list:
                          logging.warn('** Tuple %s already in match list' % \
                                       (str(match_tuple)))
                        else:
                          match_list.append(match_tuple)
                          final_match_dict[match_score] = match_list

      # Case 2: There are also locality but no street matches
      #
      elif (loc_match_list != []):

        print loc_match_list # *********************

        for (loc_val, loc_score, loc_pid_set) in loc_match_list:

          for loc_pid in loc_pid_set:

            for (addr_val, addr_score, addr_loc_pid_dict) in addr_match_list:

              # Check if LOCALITY_PID is common on address and locality level
              #
              if (loc_pid in addr_loc_pid_dict):

                street_pid_dict = addr_loc_pid_dict[loc_pid]

                # Only take matches which identify a specific address site
                #
                if (len(street_pid_dict) == 1):
                  match_score =   addr_score+loc_score
                  match_pid_set = street_pid_dict.values()[0]
                  match_descr = {'address':addr_val, 'locality':loc_val}
                  match_tuple = (match_pid_set, match_descr)

                  match_list = final_match_dict.get(match_score, [])
                  if match_tuple in match_list:
                    logging.warn('** Tuple %s already in match list' % \
                                 (str(match_tuple)))
                  else:
                    match_list.append(match_tuple)
                    final_match_dict[match_score] = match_list

      # Case 3: There are also street matches but no locality matches
      #
      elif (street_match_list != []):

        for (street_val, street_score, street_loc_pid_dict) in \
          street_match_list:

          for loc_pid in street_loc_pid_dict:

            for (addr_val, addr_score, addr_loc_pid_dict) in addr_match_list:

              # Is LOCALITY_PID common in address and street level
              #
              if (loc_pid in addr_loc_pid_dict):
                addr_street_pid_dict = addr_loc_pid_dict[loc_pid]
                street_street_loc_set = street_loc_pid_dict[loc_pid]

                for street_pid in street_street_loc_set:

                  # Is STREET_PID common in address and street level
                  #
                  if (street_pid in addr_street_pid_dict):

                    match_score =   addr_score+street_score
                    match_pid_set = addr_street_pid_dict[street_pid]
                    match_descr = {'address':addr_val, 'street':street_val}
                    match_tuple = (match_pid_set, match_descr)

                    match_list = final_match_dict.get(match_score, [])
                    if match_tuple in match_list:
                      logging.warn('** Tuple %s already in match list' % \
                                   (str(match_tuple)))
                    else:
                      match_list.append(match_tuple)
                      final_match_dict[match_score] = match_list

      else:  # Case 4: There are no locality and no street matches

        for (addr_val, addr_score, addr_loc_pid_dict) in addr_match_list:

          if (len(addr_loc_pid_dict) == 1):  # Appears in one locality only

            loc_pid, street_pid_dict = addr_loc_pid_dict.items()[0]

            if (len(street_pid_dict) == 1):  # Appears in one street only

              street_pid, addr_site_pid_set = street_pid_dict.items()[0]

              # Only take matches which identify a specific address site
              #
              if (len(addr_site_pid_set) == 1):  # One address site only

                logging.warn('** Unique address site (no street/loc): %s' % \
                             (str(addr_match_list)))

                match_score =   addr_score
                match_pid_set = addr_site_pid_set.pop()
                match_descr = {'address':addr_val}
                match_tuple = (match_pid_set, match_descr)

                match_list = final_match_dict.get(match_score, [])
                if match_tuple in match_list:
                  logging.warn('** Tuple %s already in match list' % \
                               (str(match_tuple)))
                else:
                  match_list.append(match_tuple)
                  final_match_dict[match_score] = match_list

      if (final_match_dict != {}):
        match_level = 'address'

    # If no address matches try street matches next - - - - - - - - - - - - - -
    # (if no final address level matches were found,
    # and there are street matches)
    #
    if ((final_match_dict == {}) and (street_match_list != [])):

      # Case 1: There are also locality matches
      #
      if (loc_match_list != []):

        for (street_val, street_score, loc_pid_dict) in street_match_list:

          for loc_pid in loc_pid_dict:  # LOCALITY_PIDs from street

            for (loc_val, loc_score, loc_pid_set) in loc_match_list:

              if (loc_pid in loc_pid_set):  # LOCALITY_PID also on loc. level

                match_score =   street_score+loc_score
                match_pid_set = sets.Set()
                for street_pid in loc_pid_dict[loc_pid]:
                  match_pid_set.add((loc_pid, street_pid))
                match_descr = {'street':street_val, 'locality':loc_val}
                match_tuple = (match_pid_set, match_descr)

                match_list = final_match_dict.get(match_score, [])
                if match_tuple in match_list:
                  logging.warn('** Tuple %s already in match list' % \
                               (str(match_tuple)))
                else:
                  match_list.append(match_tuple)
                  final_match_dict[match_score] = match_list

      else:  # Case 2: No locality matches

        for (street_val, street_score, loc_pid_dict) in street_match_list:

          for (loc_pid, street_pid_set) in loc_pid_dict.items():

            # Only take matches which identify a specific street
            #
            if (len(street_pid_set) == 1):
              match_score =     street_score
              match_pid_tuple = sets.Set([(loc_pid, street_pid_set.pop())])
              match_descr =     {'street':street_val}
              match_tuple = (match_pid_tuple, match_descr)

              match_list = final_match_dict.get(match_score, [])
              if match_tuple in match_list:
                logging.warn('** Tuple %s already in match list' % \
                             (str(match_tuple)))
              else:
                match_list.append(match_tuple)
                final_match_dict[match_score] = match_list

      if (final_match_dict != {}):
        match_level = 'street'

    # If no street matches try locality matches - - - - - - - - - - - - - - - -
    # (if no final address or street level matches
    # were found, and there are locality matches)
    #
    if ((final_match_dict == {}) and (loc_match_list != [])):

      for (loc_val, loc_score, loc_pid_set) in loc_match_list:

        # Simply append the tuples
        #
        match_score =   loc_score
        match_pid_set = loc_pid_set
        match_descr = {'locality':loc_val}
        match_tuple = (match_pid_set, match_descr)

        match_list = final_match_dict.get(match_score, [])
        if match_tuple in match_list:
          logging.warn('** Tuple %s already in match list' % \
                       (str(match_tuple)))
        else:
          match_list.append(match_tuple)
          final_match_dict[match_score] = match_list

      if (final_match_dict != {}):
        match_level = 'locality'

    if (final_match_dict == {}):
      match_level = 'no match'

    comb_match_str = ''
    match_score_list = final_match_dict.keys()
    match_score_list.sort()
    match_score_list.reverse()
    for ms in match_score_list:
      comb_match_str += os.linesep + '    %15s: %s' % \
                        (ms, final_match_dict[ms])
    logging.info('** Combined matches: %s' % (comb_match_str))

    logging.info('Number of combined matches: %d, match level: %s' % \
                 (len(final_match_dict), match_level))

    return (final_match_dict, match_level)

  # ---------------------------------------------------------------------------
  # Process the locality level
  # ---------------------------------------------------------------------------

  def get_locality_matches(self, nb_level = 0):
    """Method to match the locality level part in the input record according to
       the given neighbour level.

       The neighbour level can be set to 0, 1 or 2.

       Returns a list (sorted according to descending match score) of tuples
       made of a found locality value together with its match score and its
       set of LOCALITY_PIDs:

         [(locality_val1, match_score1, (set of LOCALITY_PIDs)),
          (locality_val2, match_score2, (set of LOCALITY_PIDs)),
          (locality_val3, match_score3, (set of LOCALITY_PIDs))]

       Returns an empty list if no matches can be found.
    """

    if (nb_level not in [0,1,2]):
      logging.exception('Illegal value for locality neighbour region level' + \
                        ' (must be one of 0, 1 or 2): %s' % (str(nb_level)))
      raise Exception

    # First find postcode and locality matches separately - - - - - - - - - - -
    #
    self.get_postcode_pids('exact', nb_level)

    self.get_locality_pids('exact', nb_level)
    if (self.loc_match_list == []):
      self.get_locality_pids('approx', nb_level)
      logging.warn('** searched locality approx')

    # Now build a list of localities, postcodes, and their combinations (pairs)
    #
    loc_comb_match_list = self.loc_match_list + self.pc_match_list

    if ((self.loc_match_list != []) and (self.pc_match_list != [])):

      for (loc_val, loc_match_score, loc_loc_pid_set) in self.loc_match_list:

        for (pc_val, pc_match_score, pc_loc_pid_set) in self.pc_match_list:

          if (loc_loc_pid_set.issubset(pc_loc_pid_set)):

            # All LOCALITY_PIDs for this locality are also in this postcode set
            #
            loc_comb_match_list.append((loc_val+' '+pc_val,
                                        loc_match_score+pc_match_score,
                                        loc_loc_pid_set))

    # Calculate minimal threshold needed for these matches - - - - - - - - - -
    #
    min_threshold = self.match_threshold * \
                    min(self.input_fields['postcode'][1],
                        self.input_fields['locality_name'][1],
                        self.input_fields['locality_qualifier'][1])

    locality_match_list = self.filter_sort_match_list(loc_comb_match_list,
                                                      min_threshold)

    logging.info('  Length of locality level match list: %d (at neighbour' \
                 % (len(locality_match_list))+' level %d)' % (nb_level))

    return locality_match_list

  # ---------------------------------------------------------------------------

  def get_postcode_pids(self, match_mode, nb_level = 0):
    """A method that finds matches for the postcode value in the input record.

       The attribute 'self.pc_match_list' will be set to a list of tuples made
       of a found postcode together with its match score and its corresponding
       set of LOCALITY_PIDs.

       The search is performed according to the given neighbouring level.

       The attribute 'self.pc_match_list' will be set to an empty list if no
       postcode match can be found.
    """

    # Check if a postcode value is available in this record - - - - - - - - - -
    #
    if ('postcode' not in self.in_rec):
      self.pc_match_list = []
      return

    pc_val = self.in_rec['postcode']  # Get postcode value from input reccord

    # Get all postcode values according to given neighbour region level - - - -
    #
    if (nb_level == 0):
      pc_value_list = [pc_val]  # Original postcode value only
    elif (nb_level == 1):
      pc_value_list = self.postcode_neighbours_1.get(pc_val, [])
    else:
      pc_value_list = self.postcode_neighbours_2.get(pc_val, [])

    # Try to match each postcode in the postcode value list - - - - - - - - - -
    #
    pc_match_list = []

    for pc in pc_value_list:
      this_match_list = self.match_value_index(pc, 'postcode', match_mode)
      if (this_match_list != []):
        pc_match_list += this_match_list

    # Multiply all matching scores with postcode base weight - - - - - - - - -
    # (adjust for neighbouring level)
    #
    pc_base_weight = self.input_fields['postcode'][1] * \
                     self.neighbour_score_adjust[nb_level]

    self.pc_match_list = self.adjust_match_list_scores(pc_match_list, \
                                                       pc_base_weight)

  # ---------------------------------------------------------------------------

  def get_locality_pids(self, match_mode, nb_level = 0):
    """A method that finds matches for the locality values in the input record.

       The attribute 'self.loc_match_list' will be set to a list of tuples
       made of a found locality together with its match score and its
       corresponding set of LOCALITY_PIDs.

       Locality values include a locality name, possible locality qualifier and
       maybe state (currently not considered).

       The search is performed according to the given neighbouring level.
       - If a neighbouring level of 0 is given, then the possible variations of
         locality name, locality name plus locality qualifier, locality
         qualifer plus locality name, and wayfare qualifier plus locality name
         are used for matching.
       - If a neighbouring level of 1 or 2 is given then the list of
         neighbouring locality name values are used for matching.

       The attribute 'self.loc_match_list' will be set to an empty list if no
       locality match can be found.
    """

    # Check if locality name or qualifier values are available in this record -
    #
    if (('locality_name' not in self.in_rec) and \
        ('locality_qualifier' not in self.in_rec)):
      self.loc_match_list = []
      return

    # Get the locality values from the input record - - - - - - - - - - - - - -
    #
    loc_name_value =  self.in_rec.get('locality_name', '')
    loc_quali_value = self.in_rec.get('locality_qualifier', '')

    # Get the corresponding base match weights
    #
    loc_name_base_weight =  self.input_fields['locality_name'][1]
    loc_quali_base_weight = self.input_fields['locality_qualifier'][1]

    # Replace whitespaces with underscores to make it conform to look-up tables
    #
    loc_name_value =  loc_name_value.replace(' ', '_')
    loc_quali_value = loc_quali_value.replace(' ', '_')

    loc_value_list =  []  # List with possible locality values
    loc_weight_list = []  # Their corresponding base weights

    # Compile a set of possible locality values according to the given - - - -
    # neighbouring level
    #
    if (nb_level == 0):

      if (loc_name_value != ''):  # Start with locality name only
        loc_value_list.append(loc_name_value)  # Locality name only
        loc_weight_list.append(loc_name_base_weight)

      if (loc_quali_value != ''):

        if (loc_name_value == ''):  # Add locality qualifier only
          loc_value_list.append(loc_quali_value)
          loc_weight_list.append(loc_quali_base_weight)

        else:  # Add combinations with locality qualifier

          loc_value_list.append(loc_name_value+'_'+loc_quali_value)
          loc_value_list.append(loc_quali_value+'_'+loc_name_value)
          loc_weight_list.append(loc_name_base_weight+loc_quali_base_weight)
          loc_weight_list.append(loc_name_base_weight+loc_quali_base_weight)

      # Check if there is also a wayfare qualifier in the input record
      #
      if ('wayfare_qualifier' in self.in_rec):
        wf_quali_value = self.in_rec['wayfare_qualifier']
        wf_quali_value = wf_quali_value.replace(' ', '_')
        wf_quali_base_weight = self.input_fields['wayfare_qualifier'][1]

        loc_value_list.append(wf_quali_value+'_'+loc_name_value)
        loc_weight_list.append(loc_name_base_weight+wf_quali_base_weight)

    elif (nb_level == 1):
      loc_value_list =  self.suburb_neighbours_1.get(loc_name_value, [])
      loc_weight_list = [loc_name_base_weight*self.neighbour_score_adjust[1]] \
                        * len(loc_value_list)
    else:
      loc_value_list =  self.suburb_neighbours_2.get(loc_name_value, [])
      loc_weight_list = [loc_name_base_weight*self.neighbour_score_adjust[2]] \
                        * len(loc_value_list)

    # Try to match each locality in the locality value list - - - - - - - - - -
    #
    loc_match_list = []

    for i in range(len(loc_value_list)):
      this_loc = loc_value_list[i]
      this_base_weight = loc_weight_list[i]
      this_match_list = self.match_value_index(this_loc, 'locality_name', \
                                               match_mode)
      if (this_match_list != []):

        # Multiply all matching scores with base weight
        #
        loc_match_list += self.adjust_match_list_scores(this_match_list, \
                                                        this_base_weight)

    self.loc_match_list = loc_match_list

  # ---------------------------------------------------------------------------
  # Process the street level
  # ---------------------------------------------------------------------------

  def get_street_matches(self, match_mode):
    """Method to match the street level part in the input record. The street
       level includes a street name, type and qualifier (or suffix) (but not
       the street number, which is in the address level part).

       The match mode can be set to: 'exact', 'approx', or 'exact/approx'.

       Returns a list (sorted according to descending match score) of tuples
       made of a found street value together with its match score and its
       dictionary of LOCALITY_PIDs and their STREET_PID sets.

         [(street_val1, match_score1, {LOCALITY_PIDa: (set of STREET_PIDs),
                                       LOCALITY_PIDb: (set of STREET_PIDs),
                                       LOCALITY_PIDc: (set of STREET_PIDs)}),
          (street_val2, match_score2, {LOCALITY_PIDa: (set of STREET_PIDs),
                                       LOCALITY_PIDc: (set of STREET_PIDs),
                                       LOCALITY_PIDd: (set of STREET_PIDs)})]

       Returns an empty list if no matches can be found.
    """

    if (match_mode not in ['exact','approx','exact/approx','approx/exact']):
      logging.exception('Illegal argument for "match_mode" given: %s' %
                        (str(match_mode)))
      raise Exception

    self.get_street_pids(match_mode)

    # Calculate minimal threshold needed for these matches - - - - - - - - - -
    #
    min_threshold = self.match_threshold * \
                    min(self.input_fields['wayfare_name'][1],
                        self.input_fields['wayfare_type'][1],
                        self.input_fields['wayfare_qualifier'][1])

    street_match_list = self.filter_sort_match_list(self.street_match_list,
                                                    min_threshold)

    #logging.info('** street match list: %s' % \
    #             (str(street_match_list))) # *********************************

    logging.info('  Length of street level match list: %d' % \
                 (len(street_match_list)))

    return street_match_list

  # ---------------------------------------------------------------------------

  def get_street_pids(self, match_mode):
    """A method that finds matches for the street values in the input record.

       The attribute 'self.street_match_list' will be set to a list of tuples
       made of a found street together with its match score and its
       corresponding dictionary of LOCALITY_PIDs and their STREET_PID sets.

       Street values include a street name, type and qualifier.

       The attribute 'self.street_match_list' will be set to an empty list if
       no street match can be found.
    """

    # Check if there is any street information in the record - - - - - - - - -
    # (only a wayfare qualifier is not very useful)
    #
    if (('wayfare_name' not in self.in_rec) and \
        ('wayfare_type' not in self.in_rec)):
      self.street_match_list = []
      return

    # Get the street values from the input record - - - - - - - - - - - - - - -
    #
    st_name_value =  self.in_rec.get('wayfare_name', '')
    st_type_value =  self.in_rec.get('wayfare_type', '')
    st_quali_value = self.in_rec.get('wayfare_qualifier', '')

    # Get the corresponding base match weights
    #
    st_name_base_weight =  self.input_fields['wayfare_name'][1]
    st_type_base_weight =  self.input_fields['wayfare_type'][1]
    st_quali_base_weight = self.input_fields['wayfare_qualifier'][1]

    # Replace whitespaces with underscores to make it conform to look-up tables
    #
    st_name_value =  st_name_value.replace(' ', '_')
    st_quali_value = st_quali_value.replace(' ', '_')

    # Get the matches for the three street components - - - - - - - - - - - - -
    #
    st_name_match_list = self.match_value_index(st_name_value, 'street_name',
                                                match_mode)
    st_type_match_list = self.match_value_index(st_type_value, 'street_type',
                                                'exact')
    st_quali_match_list = self.match_value_index(st_quali_value,
                                                 'street_suffix', 'exact')

    # Multiply all matching scores with their base weights - - - - - - - - - -
    #
    st_name_match_list = self.adjust_match_list_scores(st_name_match_list,
                                                       st_name_base_weight)
    st_type_match_list = self.adjust_match_list_scores(st_type_match_list,
                                                       st_type_base_weight)
    st_quali_match_list = self.adjust_match_list_scores(st_quali_match_list,
                                                        st_quali_base_weight)

    # Find all non-empty intersection of the three street components - - - - -
    #
    street_match_list = []

    # Loop over all possible combinations (hopefully not many)
    #
    for st_name_match in st_name_match_list:

      street_match_list.append(st_name_match)  # Street name match only

      for st_type_match in st_type_match_list:

        isect_dict = self.isect_street_match_dicts(st_name_match[2],
                                                   st_type_match[2])
        if (isect_dict != None):

          # Street name and type combination
          #
          street_name_type_val =    st_name_match[0]+' '+st_type_match[0]
          street_name_type_weight = st_name_match[1]+st_type_match[1]
          street_match_list.append((street_name_type_val,
                                   street_name_type_weight, isect_dict))

          for st_quali in st_quali_match_list:

            isect_dict2 = self.isect_street_match_dicts(isec_dict,
                                                        st_quali_match[2])
            if (isect_dict2 != None):

               # Street name, type and qualifier combination
               #
               street_match_list.append((street_name_type_val+' '+ \
                                         st_quali_match[0],
                                         street_name_type_weight+ \
                                         st_quali_match[1], isect_dict2))

    # If no street name match has been found append street type - - - - - - - -
    # (only append if not too many matches for a street type)
    #
    if (st_name_match_list == []):
      for st_type_match in st_type_match_list:

        if (len(st_type_match[2]) < 100):
          street_match_list.append(st_type_match)  # Street type match only

    self.street_match_list = street_match_list

  # ---------------------------------------------------------------------------
  # Process the address level
  # ---------------------------------------------------------------------------

  def get_address_matches(self, match_mode):
    """Method to match the address level part in the input record. The address
       level includes a street number (incl. possibly pre- and suffix), and
       unit and property values (if available in the input record).

       The match mode can be set to: 'exact', 'approx', or 'exact/approx'.

       Returns a list (sorted according to descending match score) of tuples
       made of a found address value together with its match score and its
       dictionary of LOCALITY_PIDs and their STREET_PID dictionaries with
       ADDRESS_SITE_PID sets.

         [(addr_val1, match_score1,
              {LOCALITY_PIDa: {STREET_PID1: (set of ADDRESS_SITE_PIDs),
                               STREET_PID2: (set of ADDRESS_SITE_PIDs),
                               STREET_PID3: (set of ADDRESS_SITE_PIDs)},
               LOCALITY_PIDb: {STREET_PID5: (set of ADDRESS_SITE_PIDs),
                               STREET_PID6: (set of ADDRESS_SITE_PIDs),
                               STREET_PID7: (set of ADDRESS_SITE_PIDs)}}),
          (addr_val2, match_score2,
              {LOCALITY_PIDc: {STREET_PID1: (set of ADDRESS_SITE_PIDs),
                               STREET_PID2: (set of ADDRESS_SITE_PIDs),
                               STREET_PID3: (set of ADDRESS_SITE_PIDs)},
               LOCALITY_PIDd: {STREET_PID5: (set of ADDRESS_SITE_PIDs),
                               STREET_PID6: (set of ADDRESS_SITE_PIDs),
                               STREET_PID7: (set of ADDRESS_SITE_PIDs)}})]

       Returns an empty list if no matches can be found.
    """

    if (match_mode not in ['exact','approx','exact/approx','approx/exact']):
      logging.exception('Illegal argument for "match_mode" given: %s' %
                        (str(match_mode)))
      raise Exception

    self.get_street_num_pids(match_mode)

    # Calculate minimal threshold needed for these matches - - - - - - - - - -
    #
    min_threshold = self.match_threshold*self.input_fields['wayfare_number'][1]

    st_num_match_list = self.filter_sort_match_list(self.street_num_match_list,
                                                    min_threshold)

    #logging.info('** street number match list: %s' % \
    #             (str(st_num_match_list))) # *********************************

    logging.info('  Length of address level match list: %d' % \
                 (len(st_num_match_list)))

    return st_num_match_list

  # ---------------------------------------------------------------------------

  def get_address_refinement(self, match_mode, address_match_list):
    """Method which refines the address level match with unit and property
       values if available in the input record.

       It combines the given address match list with the match lists
       retrieved using unit and property information (if available), and adds
       the combinations to the give naddress match list.

       The match mode can be set to: 'exact', 'approx', or 'exact/approx'.

       *** only reduce the address_pid_sets of entries which have street and
       locality pids in common, otherwise don't modify them ***

       ***** What about 'flat' and 'lot' fields? *****
    """

    if (match_mode not in ['exact','approx','exact/approx','approx/exact']):
      logging.exception('Illegal argument for "match_mode" given: %s' %
                        (str(match_mode)))
      raise Exception

    # Get unit and property match sets
    #
    self.get_unit_pids(match_mode)
    self.get_property_pids(match_mode)

    # Check if unit and property matches are available, if not return - - - - -
    #
    if (self.unit_match_list == []) and (self.property_match_list == []):
      return address_match_list

    logging.warn(self.match_list_summary_string(self.unit_match_list, 'unit'))
    logging.warn(self.match_list_summary_string(self.property_match_list,
                 'property'))

    # Check if both unit and property matches are available - - - - - - - - - -
    #
    if (self.unit_match_list != []) and (self.property_match_list != []):

      addr_refine_match_list = []

      # Loop over all possible combinations (hopefully not many)
      #
      for unit_match in self.unit_match_list:

        addr_refine_match_list.append(unit_match)  # Unit match only

        for prop_match in self.property_match_list:

          isect_dict = self.isect_address_match_dicts(unit_match[2],
                                                      prop_match[2])
          if (isect_dict != None):

            # Unit and property combination
            #
            addr_refine_match_list.append(('"'+prop_match[0]+'" '+ \
                                           unit_match[0],
                                           unit_match[1]+prop_match[1],
                                           isect_dict))

      # Also add all property matches separately
      #
      for prop_match in self.property_match_list:
        addr_refine_match_list.append(prop_match)

    elif (self.unit_match_list != []):  # Only unit matches exist
      addr_refine_match_list = self.unit_match_list

    elif (self.property_match_list != []):  # Only property matches exist
      addr_refine_match_list = self.property_match_list

    # If address matches are already available find intersections of all pairs
    #
    if (address_match_list != []):

      new_addr_match_list = []

      for addr_match in address_match_list:

        for unit_prop_match in addr_refine_match_list:

          isect_dict = self.isect_address_match_dicts(addr_match[2],
                                                      unit_prop_match[2])
          if (isect_dict != None):

            # Addres, unit and property combination
            #
            new_addr_match_list.append((unit_prop_match[0]+' / '+addr_match[0],
                                        addr_match[1]+unit_prop_match[1],
                                        isect_dict))

      new_addr_match_list = address_match_list + new_addr_match_list

    else:  # No address matches available - - - - - - - - - - - - - - - - - - -

      new_addr_match_list = addr_refine_match_list

    # Calculate minimal threshold needed for these matches - - - - - - - - - -
    #
    min_threshold = self.match_threshold * \
                    min(self.input_fields['wayfare_number'][1],
                        self.input_fields['flat_number'][1],
                        self.input_fields['flat_type'][1],
                        self.input_fields['building_name'][1])

    new_address_match_list = self.filter_sort_match_list(new_addr_match_list,
                                                         min_threshold)

    #logging.info('** refined address match list: %s' % \
    #             (str(new_address_match_list))) # ****************************

    logging.debug('  Length of refined address match list: %d' % \
                  (len(new_address_match_list)))

    return new_address_match_list

  # ---------------------------------------------------------------------------

  def get_street_num_pids(self, match_mode):
    """A method that finds matches for the street number values in the input
       record.

       The attribute 'self.street_num_match_list' will be set to a list with
       tuples made of a street number value, its match score and its
       corresponding dictionary of LOCALITY_PIDs and their dictionaries with
       STREET_PID as keys and sets of ADDRESS_SITE_PIDs as values (see method
       get_address_matches() above for more details).

       The attribute 'self.street_num_match_list' will be set to an empty list
       if no strret number match can be found.
    """

    # Check if a street (wayfare) number value is available in this record - -
    #
    if ('wayfare_number' not in self.in_rec):
      self.street_num_match_list = []
      return

    # Get the street number values from the input record - - - - - - - - - - -
    #
    street_num_val = self.in_rec.get('wayfare_number', '')

    # Get the corresponding base match weight
    #
    street_num_base_weight = self.input_fields['wayfare_number'][1]

    # Check for most commen case with one number only - - - - - - - - - - - - -
    #
    if (street_num_val.isdigit() == True):

      street_num_match_list = self.match_value_index(street_num_val,
                                                     'number_first',
                                                     'exact')
      self.street_num_match_list = \
        self.adjust_match_list_scores(street_num_match_list,
                                      street_num_base_weight)
      return

    # Otherwise split the number into its components and match them - - - - - -
    #
    street_num_list = self.number_components(street_num_val)

    # If the list length is one element only, then it must be letters
    # (as digits have already been handled)
    #
    if (len(street_num_list) == 1):
      street_num_match_list = self.match_value_index(street_num_list[0],
                                                     'number_first_prefix',
                                                     'exact')
      self.street_num_match_list = \
        self.adjust_match_list_scores(street_num_match_list,
                                      street_num_base_weight)
      return

    # *************************************************************************
    # Problem: We now have more than one number component which we can match
    # with up to 6 indices:
    #
    index_name_list = ['number_first_prefix', 'number_first', \
                       'number_first_suffix', 'number_last_prefix', \
                       'number_last', 'number_last_suffix']

    logging.warn('  More than one street number component: %s' % \
                 (street_num_list))

    street_num_match_list = []

    candidate_index = 0  # A pointer into the candicate index name list
    done_num_first = False

    for num_comp in street_num_list:

      if (num_comp.isdigit() == True):  # A digit, thus number first or last
        if (candidate_index <= 1):
          candidate_index = 1  # Move to number first index
        else:
          candidate_index = 4  # Move to number last index

      elif (num_comp.isalpha() == True):  # A letter, thus pre- or suffix
        if (candidate_index in [1,4]):
          candidate_index += 1  # Move tonext index

      else:  # This should not happen
        logging.warn('    Strange number component: "%s"' % (num_comp))

      index_name = index_name_list[candidate_index]

      logging.warn('    Match value "%s" with index %s' % \
                   (num_comp, index_name))

      tmp_list = self.match_value_index(num_comp, index_name, 'exact')

      if (tmp_list != []):  # Only use non-empty match lists
        street_num_match_list += tmp_list

      candidate_index += 1

    if (len(street_num_match_list) == 1):
      self.street_num_match_list = \
        self.adjust_match_list_scores(street_num_match_list,
                                      street_num_base_weight)
      return

    # Now try to intersect as many as possible - - - - - - - - - - - - - - - -

    tmp_street_num_match_list = []
    for i in range(len(street_num_match_list)-1):

      street_match_tuple1 = street_num_match_list[i]
      street_match_tuple2 = street_num_match_list[i+1]

      isect_dict = self.isect_address_match_dicts(street_match_tuple1[2],
                                                  street_match_tuple2[2])
      if (isect_dict != None):
        street_num_val = street_num_match_list[i][0]+' '+ \
                         street_num_match_list[i+1][0]
        street_num_score = street_num_match_list[i][1]+ \
                         street_num_match_list[i+1][1]
        tmp_street_num_match_list.append((street_num_val, street_num_score,
                                          isect_dict))

    street_num_match_list = street_num_match_list + tmp_street_num_match_list

    self.street_num_match_list = \
      self.adjust_match_list_scores(street_num_match_list,
                                    street_num_base_weight)

    logging.warn(self.match_list_summary_string(self.street_num_match_list,
                                                'Final addess'))

  # ---------------------------------------------------------------------------

  def get_unit_pids(self, match_mode):
    """A method that finds matches for the unit values in the input record.

       The 'self.unit_match_list' will be set to a list with tuples made of a
       unit value, its match score and its corresponding dictionary of
       LOCALITY_PIDs and their dictionaries with STREET_PID as keys and sets of
       ADDRESS_SITE_PIDs as values (see method get_address_matches() above for
       more details).

       The attribute 'self.unit_match_list' will be set to an empty list if no
       unit match can be found.
    """

    # Check if a unit value is available in this record - - - - - - - - - - - -
    #
    if ('unit_number' not in self.in_rec) and ('unit_type' not in self.in_rec):
      self.unit_match_list = []
      return

    # Get the unit type and number values from the input record - - - - - - - -
    #
    unit_num_val =  self.in_rec.get('unit_number', '')
    unit_type_val = self.in_rec.get('unit_type', '')

    # Get the corresponding base match weight
    #
    unit_num_base_weight =  self.input_fields['flat_number'][1]
    unit_type_base_weight = self.input_fields['flat_type'][1]

    # Get the matches for the unit type - - - - - - - - - - - - - - - - - - - -
    #
    unit_type_match_list = self.match_value_index(unit_type_val, 'flat_type',
                                                  match_mode)
    # Adjust the match scores
    #
    unit_type_match_list = \
      self.adjust_match_list_scores(unit_type_match_list,
                                    unit_type_base_weight)

    # If there is no unit number return unit type match list - - - - - - - - -
    #
    if (unit_num_val == ''):
      self.unit_match_list = unit_type_match_list
      return

    # Check for most commen case with one number only - - - - - - - - - - - - -
    #
    if (unit_num_val.isdigit() == True):
      unit_num_match_list = self.match_value_index(unit_num_val, 'flat_number',
                                                   'exact')

    else:  # Split unit number into its components and match them - - - - - - -

      unit_num_list = self.number_components(unit_num_val)

      # If the list length is one element only, then it must be letters
      # (as digits have already been handled)
      #
      if (len(unit_num_list) == 1):
        unit_num_match_list = self.match_value_index(unit_num_list[0],
                                                     'flat_number_prefix',
                                                     'exact')
      else:

        # Unit number is assumed to be either 'xx42', '42yy' or 'xx42yy'

        # If first list element are digits then match with unit number - - - -
        #
        if (unit_num_list[0].isdigit() == True):
          unit_num_prefix_match_list = []
          unit_num_match_list = self.match_value_index(unit_num_list[0],
                                                       'flat_number',
                                                       'exact')
          unit_num_suffix_match_list = \
            self.match_value_index(unit_num_list[1], 'flat_number_suffix',
                                   'exact')

          if (len(unit_num_list) > 2):  # ******* How to handle this? ********
            logging.warn('** Strange unit number: "%s"' % (unit_num_val))

        else:  # Match with number prefix
          unit_num_prefix_match_list = \
            self.match_value_index(unit_num_list[0], 'flat_number_prefix',
                                   'exact')
          unit_num_match_list = self.match_value_index(unit_num_list[1],
                                                       'flat_number',
                                                       'exact')
          if (len(unit_num_list) > 2):
            unit_num_suffix_match_list = \
              self.match_value_index(unit_num_list[2], 'flat_number_suffix',
                                     'exact')
          else:
            unit_num_suffix_match_list = []

          if (len(unit_num_list) > 3):  # ******* How to handle this? ********
            logging.warn('** Strange unit number: "%s"' % (unit_num_val))

        # Now try to intersect as many as possible - - - - - - - - - - - - - -
        #
        if (unit_num_prefix_match_list != []):  # Intersection with prefix

          unit_num_match_list += unit_num_prefix_match_list # Add prefix itself

          unit_num_tuple =        unit_num_match_list[0]
          unit_num_prefix_tuple = unit_num_prefix_match_list[0]

          isect_dict = self.isect_address_match_dicts(unit_num_tuple[2],
                                                      unit_num_prefix_tuple[2])
          if (isect_dict != None):
            unit_num_val =   unit_num_prefix_tuple[0]+' '+unit_num_tuple[0]
            unit_num_score = unit_num_prefix_tuple[1]+unit_num_tuple[1]
            unit_num_match_list.append((unit_num_val, unit_num_score,
                                        isect_dict))

            # Try intersection with prefix and suffix
            #
            if (unit_num_suffix_match_list != []):
              unit_num_suffix_tuple = unit_num_suffix_match_list[0]

              isect_dict2 = self.isect_address_match_dicts(isect_dict,
                                                      unit_num_suffix_tuple[2])
              if (isect_dict2 != None):
                unit_num_val = unit_num_prefix_tuple[0]+' '+ \
                               unit_num_tuple[0]+' '+unit_num_suffix_tuple[0]
                unit_num_score = unit_num_prefix_tuple[1]+unit_num_tuple[1]+ \
                                 unit_num_suffix_tuple[1]
                unit_num_match_list.append((unit_num_val, unit_num_score,
                                            isect_dict2))

        if (unit_num_suffix_match_list != []):  # Intersection with suffix

          unit_num_match_list += unit_num_suffix_match_list # Add suffix itself

          unit_num_tuple =        unit_num_match_list[0]
          unit_num_suffix_tuple = unit_num_suffix_match_list[0]

          isect_dict = self.isect_address_match_dicts(unit_num_tuple[2],
                                                      unit_num_suffix_tuple[2])
          if (isect_dict != None):
            unit_num_val =   unit_num_tuple[0]+' '+unit_num_suffix_tuple[0]
            unit_num_score = unit_num_tuple[1]+unit_num_suffix_tuple[1]
            unit_num_match_list.append((unit_num_val, unit_num_score,
                                        isect_dict))

    # Adjust the match scores for unit number
    #
    unit_num_match_list = \
      self.adjust_match_list_scores(unit_num_match_list, unit_num_base_weight)

    # Last try to intersect with unit type - - - - - - - - - - - - - - - - - -
    #
    if (unit_type_match_list != []):

      isect_match_list = []

      for unit_num_match_tuple in unit_num_match_list:
        for unit_type_match_tuple in unit_type_match_list:

          isect_dict = self.isect_address_match_dicts(unit_num_match_tuple[2],
                                                      unit_type_match_tuple[2])
          if (isect_dict != None):
            unit_val = unit_num_match_tuple[0]+' '+unit_type_match_tuple[0]
            unit_score = unit_num_match_tuple[1]+unit_type_match_tuple[1]
            isect_match_list.append((unit_val, unit_score, isect_dict))

      if (isect_match_list != []):  # Non-empty intersections found
        self.unit_match_list = isect_match_list+unit_num_match_list

      else:
        self.unit_match_list = unit_num_match_list  # Unit numbers only

    else:
      self.unit_match_list = unit_num_match_list

  # ---------------------------------------------------------------------------

  def get_property_pids(self, match_mode):
    """A method that finds matches for the property value in the input record.

       The attribute 'self.property_match_list' will be set to a list with
       tuples made of a property value, its match score and its corresponding
       dictionary of LOCALITY_PIDs and their dictionaries with STREET_PID as
       keys and sets of ADDRESS_SITE_PIDs as values (see method
       get_address_matches() above for more details).

       The attribute 'self.property_match_list' will be set to an empty list if
       no property match can be found.
    """

    property_val = self.in_rec.get('property_name', '')  # Get property value

    property_match_list = self.match_value_index(property_val, 'building_name',
                                                 match_mode)
    # Get the corresponding base match weight
    #
    property_base_weight = self.input_fields['building_name'][1]

    self.property_match_list = \
      self.adjust_match_list_scores(property_match_list, property_base_weight)

  # ---------------------------------------------------------------------------
  # Various methods used by above methods for geocode matching
  # ---------------------------------------------------------------------------

  def isect_address_match_dicts(self, dict1, dict2):
    """Method that intersects the two given address level dictionaries which
       are assumed to contain dictionaires with sets (see method
       get_address_matches() above for more details about these dictionaries).

       Each dictionary is assumed to contain LOCALITY_PIDs as keys and
       dictionaries with STREET_PIDs as keys and sets of ADDRESS_SITE_PIDS as
       values.

       If the intersection results in an empty dictionary, None is returned.
    """

    intersect_dict = {}

    if (len(dict1) < len(dict2)):  # Select reference to shorter dictionary
      dict_a = dict1
      dict_b = dict2
    else:
      dict_a = dict2
      dict_b = dict1

    for loc_pid in dict_a:
      if loc_pid in dict_b:

        street_dict1 = dict_a[loc_pid]
        street_dict2 = dict_b[loc_pid]

        if (len(street_dict1) < len(street_dict2)):
          street_dict_a = street_dict1
          street_dict_b = street_dict2
        else:
          street_dict_a = street_dict2
          street_dict_b = street_dict1

        street_intersect_dict = {}

        for street_pid in street_dict_a:
          if street_pid in street_dict_b:
            intersect_set = \
              street_dict_a[street_pid].intersection(street_dict_b[street_pid])

            if (intersect_set != sets.Set()):  # Intersection is non-empty
              street_intersect_dict[street_pid] = intersect_set

        if (street_intersect_dict != {}):
          intersect_dict[loc_pid] = street_intersect_dict

    if (intersect_dict == {}):
      return None
    else:
      return intersect_dict

  # ---------------------------------------------------------------------------

  def isect_street_match_dicts(self, dict1, dict2):
    """Method that intersects the two given street level dictionaries which are
       assumed to contain sets (see method get_street_matches() above for more
       details about these dictionaries).

       Each dictionary is assumed to contain LOCALITY_PIDs as keys and sets of
       STREET_PIDs as values.

       If the intersection results in an empty dictionary, None is returned.
    """

    intersect_dict = {}

    if (len(dict1) < len(dict2)):  # Select reference to shorter dictionary
      dict_a = dict1
      dict_b = dict2
    else:
      dict_a = dict2
      dict_b = dict1

    for loc_pid in dict_a:
      if loc_pid in dict_b:
        intersect_set = dict_a[loc_pid].intersection(dict_b[loc_pid])

        if (intersect_set != sets.Set()):  # Intersection is non-empty
          intersect_dict[loc_pid] = intersect_set

    if (intersect_dict == {}):
      return None
    else:
      return intersect_dict

  # ---------------------------------------------------------------------------

  def match_value_index(self, in_val, index_name, match_mode):
    """Method to match the given input value in the input record (if available
       in the input record) with the inverted index of the given name using
       either exact or approximate matching mode, or both.

       The match mode can be set to: 'exact', 'approx', or 'exact/approx'.

       Returns a list of one or more matches, each being a tuple of the match
       value, its match score (1.0 for exact matches and less than 1.0 for
       approximate matches) and the corresponding match set or match
       dictionary.

       If the field is not available in the input record or no matches can be
       found the method returns an empty list.
    """

    if (in_val == ''):  # Empty input value
      return []

    in_val_dict = self.index_files_data[index_name]  # Dictionary for the index

    match_res_list = []  # The list of match tuples to be returned

    # Perform exact matching - - - - - - - - - - - - - - - - - - - - - - - - -
    #
    if (match_mode in ['exact','exact/approx','approx/exact']):

      if (in_val in in_val_dict):
        match_res_list.append((in_val, 1.0, in_val_dict[in_val]))
        logging.debug('  Found exact match for value "%s" in index "%s"' % \
                      (in_val, index_name))

    # Perform approximate matching - - - - - - - - - - - - - - - - - - - - - -
    #
    if ((match_mode in ['approx','exact/approx','approx/exact']) and \
        (index_name in self.approx_index)):

      approx_index = self.approx_index[index_name]
      approx_match_list = approx_index.get_matches(in_val)

      # Find best matches with a match score larger than the threshold using
      # the Winkler approximate comparison method
      #
      best_match_list = approx_index.select_best_matches(in_val,
                                                         approx_match_list,
                                                         'winkler',
                                                         self.match_threshold)

      # Append the approximate matches to the match results
      #
      for (match_val, match_score) in best_match_list:
        match_res_list.append((match_val, match_score, in_val_dict[match_val]))

      logging.debug('Found %d approximate matches for input value "%s"' % \
                    (len(best_match_list), in_val) + ' in index "%s": %s' % \
                    (index_name, str(best_match_list)))

    # Give warning if set to approximate matching but no index available  - - -
    #
    elif ((match_mode in ['approx','exact/approx','approx/exact']) and \
          (index_name not in self.approx_index)):
      logging.warn('Match mode includes "approximate" but no approximate ' + \
                   'index available for index "%s"' % (index_name))

    return match_res_list

  # ---------------------------------------------------------------------------

  def adjust_match_list_scores(self, in_match_list, base_match_weight):
    """Method which multiplies all match scores in the given list of matches
       with the given base match weight.

       This is needed to adjust the single match scores to the overall match
       weights.

       Returns a new list with the same tuples and adjusted weights.
    """

    out_match_list = []

    for (match_value, match_score, match_data) in in_match_list:
      new_match_score = match_score * base_match_weight
      out_match_list.append((match_value, new_match_score, match_data))

    return out_match_list

  # ---------------------------------------------------------------------------

  def filter_sort_match_list(self, in_match_list, min_threshold):
    """Method which removes all elements which have a match score less than the
       given minimal threshold 'min_threshold', and then sorts the remaining
       elements according to their match score in decreasing order (high match
       scores first).

       It is assumed that each element in the input match list contains a tuple
       (match_value, match_score, match_data), with the match data being either
       a list or dictionary containing the matched PIDs.

       Returns a reduced and sorted list.
    """

    tmp_match_list = []

    for (match_value, match_score, match_data) in in_match_list:
      if (match_score >= min_threshold):
        tmp_match_list.append((match_score, match_value, match_data))

    tmp_match_list.sort()
    tmp_match_list.reverse()

    out_match_list = []

    for (match_score, match_value, match_data) in tmp_match_list:
      out_match_list.append((match_value, match_score, match_data))

    return out_match_list

  # ---------------------------------------------------------------------------

  def number_components(self, in_str):
    """Method to process a number string (e.g. street or unit) into it's
       components, i.e. substrings which are digits or letters only. The

       Returns a list with the components, or an empty list if the input string
       is empty.

       For example, input "a1c" will be returned as ['a','1','c']
    """

    in_str = in_str.strip()  # Remove sourrounding whitespaces

    if (in_str == ''):
      return []  # Empty input

    # Check if the input string is digits or letters only
    #
    if (in_str.isdigit() == True) or (in_str.isalpha() == True):
      return [in_str]

    # First split a whitespace characters
    #
    in_list =   in_str.split(' ')
    res_list = []  # Result list to be returned

    for i in range(len(in_list)):  # Process each element

      # Check if the list component is only digits or only letters
      #
      if ((in_list[i].isdigit() == True) or (in_list[i].isalpha() == True)):
        res_list.append(in_list[i])  # Use input component directly

      else:  # List component is not digits nor letters only
        this_str = in_list[i]
        tmp_list = []
        tmp_str =  this_str[0]
        for j in range(1,len(this_str)):

          # Check if this character and next character are of different type
          #
          if (this_str[j-1].isdigit() != this_str[j].isdigit()) or \
             (this_str[j-1].isalpha() != this_str[j].isalpha()):
            tmp_list.append(tmp_str)
            tmp_str = ''

          tmp_str += this_str[j]  # Start a new temporary string

        tmp_list.append(tmp_str)  # Append last component

        res_list += tmp_list

    return res_list

  # ---------------------------------------------------------------------------

  def load_shelve_pickle_file(self, file_name):
    """Opens and loads a shelve or pickle file (according to the given file
       extension and returns the dictionary or shelve.
       - It is assumed the given file name has one of the defined pickle or
         shelve file extensions.
       - It is also assumed the file is in the geocode file directory.
    """

    pik_elen = len(self.pickle_files_extension)
    slv_elen = len(self.shelve_files_extension)

    # Construct the file name and check the file extension - - - - - - - - - -
    #
    full_file_name = self.geocode_file_directory + file_name

    logging.info('  Open file "%s"' % (full_file_name))

    # Check the file extension open it
    #
    if (full_file_name[-pik_elen:] == self.pickle_files_extension):

      # Open a pickle file
      #
      pickle_file = open(full_file_name)
      load_dict = pickle.load(pickle_file)
      pickle_file.close()

    elif (full_file_name[-slv_elen:] == self.shelve_files_extension):

      # Open a shelve file
      #
      load_dict = shelve.open(full_file_name, flag='r')

    else:  # Illegal file extension
      logging.exception('File name "%s" has an illegal file extension' % \
                        (self.geocode_ref_file))
      raise Exception

    return load_dict

  # ---------------------------------------------------------------------------

  def match_list_summary_string(self, in_match_list, match_level):
    """Returns a string representation of a summary for the given match list
       (i.e. match values and scores as well as length of their match sets or
       dictionaries (but not the values itself).

       The argument "match_level" can be either 'address', 'street' or
       'locality'.
    """

    sum_string = '  '+match_level.capitalize()+' level match list:'+os.linesep
    sum_string += '    Match value          Match score   Number of matches'+ \
                  os.linesep

    if (in_match_list == []):
      sum_string += '      *** No matches found *** '

    else:
      for match in in_match_list:
        sum_string += '    %20s %1.9f   %d' % (match[0].ljust(18), match[1],
                                                len(match[2]))+os.linesep

    return sum_string[:-1]  # Return without last line separator

  # ---------------------------------------------------------------------------
