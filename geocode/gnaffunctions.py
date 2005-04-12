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
# The Original Software is: "gnaffunctions.py"
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

"""Module gnaffunctions.py - Module with various functions needed to process
                             the GNAF data files.

   This module contains a class 'ProcessGNAF' with methods used to process the
   original G-NAF (Australian Geocoded National Address File) files, assuming
   they are available as comma separated values (CSV) files.

   While the methods in this module do the actual processing, the module
   'process-gnaf.py' is used configure the G-NAF processing (similar to the
   project.py modules used within Febrl).

   See the documentation within the methods for more details.
"""

# =============================================================================
# Imports go here

import logging

import cPickle as pickle
import csv
import shelve
import os
import sets
import time

from address import *          # Address tagging functionality
from standardisation import *  # Standardisation routines
import output

# =============================================================================

class ProcessGNAF:
  """Class ProcessGNAF - Methods and attributes used for processing a suite
                         of GNAF CSV files into inverted indices to be used
                         for the Febrl geocoding engine.
  """

  def __init__(self, **kwargs):
    """Constructor - Set various attributes.
    """

    self.gnaf_directory =   ''
    self.output_directory = ''
    self.gnaf_file_ext =    ''
    self.pickle_file_ext =  ''
    self.shelve_file_ext =  ''

    self.stats_max_num_val = 0

    self.address_correction_list = None
    self.address_lookup_table =    None

    self.aust_post_postcode_dict = None
    self.aust_post_suburb_dict =   None
    self.aust_post_pc_sub_set =    None

    for (keyword, value) in kwargs.items():
      if (keyword in ['gnaf_dir', 'gnaf_directory']):
        self.gnaf_directory = value
        if (self.gnaf_directory[-1] != os.sep):
          self.gnaf_directory = self.gnaf_directory + os.sep

      elif  (keyword in ['output_dir', 'output_directory', 'out_dir']):
        self.output_directory = value
        if (self.output_directory[-1] != os.sep):
          self.output_directory = self.output_directory + os.sep

      elif  (keyword in ['gnaf_file_ext','gnaf_file_extension']):
        self.gnaf_file_ext = value
      elif  (keyword in ['pickle_file_ext','pickle_file_extension']):
        self.pickle_file_ext = value
      elif  (keyword in ['shelve_file_ext','shelve_file_extension']):
        self.shelve_file_ext = value
      elif  (keyword in ['stats_max_num_val','stats_max_num_values']):
        self.stats_max_num_val = value

      elif  (keyword in ['address_correction_list','address_corr_list', \
                         'addr_correction_list','addr_corr_list']):
        self.address_correction_list = value
      elif  (keyword in ['address_lookup_table','addr_lookup_table']):
        self.address_lookup_table = value

      else:
        logging.exception('Illegal constructor argument keyword: "%s"' % \
                          (str(keyword)))
        raise Exception

  # ---------------------------------------------------------------------------

  def open_csv(self, gnaf_file_info):
    """Open a CSV file and return the CSV reader if successful.
       - 'gnaf_file_info' contains a file name and the names of all attributes.
       - All quotes are removed.
       - It is assumed the first line in the file contains the header line with
         the attribute (column) names.
       - The attribute names in the header line are checked against the names
         given in 'gnaf_file_info' (same number, same attribute names).
    """

    # Open the CSV file
    #
    csv_file_name = self.gnaf_directory+gnaf_file_info[0]+self.gnaf_file_ext
    try:
      csv_file = open(csv_file_name, 'r')
    except IOError:
      logging.exception('Cannot open CSV file "%s"' % (csv_file_name))
      raise IOError
    except:
      raise Exception

    logging.debug('Opened CSV file "%s" for reading' % (csv_file_name))

    # Define a CSV reader
    #
    csv_reader = csv.reader(csv_file)

    # It is assumed the first line in the file is the header - - - - - - - - -
    #
    csv_header_line = csv_reader.next()
    logging.debug('Read CVS file header: "%s"' % (csv_header_line))

    # Check if attribute names in the header line are the same as expected
    #
    if (len(csv_header_line) != len(gnaf_file_info[1])):
      logging.exception('Different number of attributes in file than ' + \
                        'expected.  File: %d, expected: %s' % \
                        (len(csv_header_line), len(gnaf_file_info[1])))
      raise Exception

    for i in range(len(csv_header_line)):
      if (csv_header_line[i] != gnaf_file_info[1][i]):
        logging.exception('Different attribute %d in file than expected.' % \
                          (i) + '  File: %s, expected: %s' % \
                          (csv_header_line[i], gnaf_file_info[1][i]))
        raise Exception

    # Return the CSV reader
    #
    return [csv_reader, csv_file]

  # ---------------------------------------------------------------------------

  def get_csv_statistics(self, csv_row_list, gnaf_file_info):
    """Analyse a list of CSV rows (lists) and collect statistics.
       - Count the number of rows with missing values for each column.
       - Make a dictionary of values for each column.
    """

    num_attr = len(gnaf_file_info[1])

    miss_values = []  # A list of missing value counters
    attr_values = []  # A list of attribute value dictionaries

    for i in range(num_attr):
      miss_values.append(0)
      attr_values.append({})

    # Main loop over rows (lists) in the list - - - - - - - - - - - - - - - - -
    #
    row_count = 0
    for r in csv_row_list:

      if (len(r) != num_attr):
        logging.warn('Wrong number of columns in row %d: %d (should be %d)' \
                     % (row_count, len(r), num_attr))

      col_count = 0

      for col_val in r:  # Loop over all columns in the row

        if (col_val == ''):  # One more missing values
          miss_values[col_count] = miss_values[col_count] + 1

        else:
          col_dict = attr_values[col_count]
          if (col_val in col_dict):
            val_count = col_dict[col_val] + 1  # Value already found in column
          else:
            val_count = 1  # A new value in this column

          col_dict[col_val] = val_count  # Insert into column dictionary
          attr_values[col_count] = col_dict

        col_count += 1

      row_count += 1

    # Print collected statistics - - - - - - - - - - - - - - - - - - - - - - -
    #
    logging.info('CVS row list statistics:')
    logging.info('  Number of rows: %d' % (len(csv_row_list)))
    for i in range(num_attr):
      logging.info('  Column name: %s' % (gnaf_file_info[1][i]))
      logging.info('    Number of rows with missing values: %d' % \
                   (miss_values[i]))
      if (len(attr_values[i]) > self.stats_max_num_val):
        values_comment = ' (showing only %d most frequent values)' % \
                         (self.stats_max_num_val)
      else:
        values_comment = ''
      logging.info('    Number of different values: %d' % \
                   (len(attr_values[i])) + values_comment)
      if (len(attr_values[i]) <= self.stats_max_num_val):
        for (k,v) in attr_values[i].items():
          logging.info('      "%s" : %d' % (k, v))
      else:  # More than stats_max_num_val values, sort for most frequent ones
        attr_values_keys = attr_values[i].keys()
        attr_values_vals = attr_values[i].values()

        values_list = map(None, attr_values_vals, attr_values_keys)
        values_list.sort()
        values_list.reverse()
        for (v,k) in values_list[:self.stats_max_num_val]:
          logging.info('      "%s" : %d' % (k, v))

  # ---------------------------------------------------------------------------

  def attr_list_to_dict(self, gnaf_file_info):
    """Convert an attribute name list (from 'gnaf_file_info') into a
       dictionary.
       - Assuming the attributes correspond to columns in a CSV file.
       - Attributes in the list are numbered starting from zero.
       - The dictionary returned contains the attribute names a keys, with the
         corresponding values being the column numbers.
    """

    attr_dict = {}

    col_num = 0
    for attr in gnaf_file_info[1]:
      attr_dict[attr] = col_num
      col_num += 1

    logging.debug('Input attribute name list:   %s' % (str(gnaf_file_info[1])))
    logging.debug('  Attribute name dictionary: %s' % (str(attr_dict)))

    return attr_dict

  # ---------------------------------------------------------------------------

  def csv_list_to_dict(self, csv_row_list, attrib_column):
    """Convert a list of CSV rows (lists) into a dictionary
       - Create a dictionary for the given attribute column (with the attribute
         column values as keys).
       - It is assumed that the values in the given attribute are unique
         (otherwise an error is triggered).
    """

    # Check if there are enough columns
    #
    if (len(csv_row_list[0]) <= attrib_column):
      logging.exception('Not enough columns in CVS row list ' + \
                        '(%d needed, only %d' % (attrib_column+1, \
                        len(csv_row_list[0])) + ' available)')
      raise Exception

    csv_dict = {}

    for r in csv_row_list:
      key = r[attrib_column]

      if (key in csv_dict):
        logging.exception('CSV dictionary already has key "%s"' % (key))
        raise Exception

      csv_dict[key] = r

    return csv_dict

  # ---------------------------------------------------------------------------

  def index_locality(self, g_locality_file_info,  g_locality_alias_file_info,
                     g_locality_geocode_file_info):
    """Clean, standardise and index the two GNAF files containing locality
       information as well as the locality geocode file.
       - Create an inverted index for the locality geocode, with LOCALITY_PIDs
         as keys and coordinates as values.
       - Use Febrl's clean_component() and tag_address_component() to clean and
         standardise the locality name and postcode.
       - Use lookup tables based on Australia Post lookup file for checking and
         possible imputing postcodes.
       - Create three separate inverted indices for locality names, postcodes
         and state abbreviations.
       - Don't insert entries which have no geocode into the inverted indices.
       - Returns inverted indices (dictionaries) for:
           locality_name
           postcode
           state_abbrev
           locality_geocode
         Values in the sets for these inverted indices are LOCALITY_PIDs.
    """

    logging.info('Build a dictionary for "G_LOCALITY_GEOCODE" and clean, ' + \
                 'standardise and index locality files "G_LOCALITY" and ' + \
                 '"G_LOCALITY_ALIAS"')

    locality_index = {'locality_geocode':{},
                      'locality_name':   {},
                      'postcode':        {},
                      'state_abbrev':    {}}

    # Counters for number of unknown / empty locality name and postcodes
    #
    unknown_geocode_count =     0  # Number of records with a LOCALITY_PID for
                                   # which no geocode is aavailable
    empty_postcode_count =      0
    unknown_postcode_count =    0
    empty_suburb_count =        0
    unknown_suburb_count =      0
    unknown_pc_sub_comb_count = 0

    min_locality_pid = '999999999'  # Find the smallest LOCALITY_PID
    max_locality_pid = '0'          # Find the largest LOCALITY_PID

    # Get direct references to Australia Post dictionaries and sets
    #
    aust_post_pc_dict =    self.aust_post_postcode_dict
    aust_post_sub_dict =   self.aust_post_suburb_dict
    aust_post_pc_sub_set = self.aust_post_pc_sub_set

    # Process the G_LOCALITY_GEOCODE file - - - - - - - - - - - - - - - - - - -

    # Get a dictionary of the column names
    #
    attr_dict = self.attr_list_to_dict(g_locality_geocode_file_info)

    # Get the column numbers
    #
    locality_pid_col = attr_dict['LOCALITY_PID']
    latitude_col =     attr_dict['LATITUDE']
    longitude_col =    attr_dict['LONGITUDE']

    # Open the CSV file
    #
    [csv_reader, csv_file] = self.open_csv(g_locality_geocode_file_info)

    logging.info('  Process "G_LOCALITY_GEOCODE" file')

    row_count =  0
    start_time = time.time()
    num_missing = 0  # Counter for number records with missing coordinates

    # Now read the file content and process it
    #
    for csv_line in csv_reader:

      locality_pid = csv_line[locality_pid_col].strip()
      latitude =     csv_line[latitude_col].strip()
      longitude =    csv_line[longitude_col].strip()

      locality_pid = self.check_pid(locality_pid, 'LOCALITY_PID', row_count)

      if (locality_pid != ''):
        min_locality_pid = min(min_locality_pid, locality_pid)
        max_locality_pid = max(max_locality_pid, locality_pid)

        if (latitude == '') or (longitude == ''):  # Coordinates missing
          logging.warn('Missing coordinates in row %d' % (row_count))
          num_missing += 1
        else:
          coordinates = latitude+','+longitude

          # Now insert into locality geocode inverted index
          #
          if (locality_pid in locality_index['locality_geocode']):
            logging.exception('Row %d: LOCALITY_PID "%s" already in locality' \
                              % (row_count, \
                              locality_geocode_index[locality_pid]) + \
                              ' geocode index')
            raise Exception
          else:
            locality_index['locality_geocode'][locality_pid] = coordinates

      row_count += 1
      if ((row_count % 10000) == 0):  # Process meter
        used_time = (time.time() - start_time)
        logging.info('    Processed %d rows in %s (time per row: %s)' % \
                     (row_count, output.time_string(used_time), \
                     output.time_string(used_time/row_count)))

    csv_file.close()
    logging.info('    Read %d data lines' % (row_count))

    logging.info('    Locality geocode index contains %d entries' % \
                 (len(locality_index['locality_geocode'])))
    logging.info('    Number of records with missing coordinates: %d' % \
                 (num_missing))
    logging.info('    Smallest LOCALITY_PID: %s' % (min_locality_pid))
    logging.info('    Largest LOCALITY_PID:  %s' % (max_locality_pid))

    # Process the G_LOCALITY and G_LOCALITY_ALIAS files - - - - - - - - - - - -

    for fi in [g_locality_file_info, g_locality_alias_file_info]:

      # Get a dictionary of the column names
      #
      attr_dict = self.attr_list_to_dict(fi)

      # Get the column numbers
      #
      locality_name_col = attr_dict['LOCALITY_NAME']
      state_abbrev_col =  attr_dict['STATE_ABBREVIATION']
      postcode_col =      attr_dict['POSTCODE']
      locality_pid_col =  attr_dict['LOCALITY_PID']

      # Open the CSV file
      #
      [csv_reader, csv_file] = self.open_csv(fi)

      logging.info('  Process "%s" file' % (fi[0]))  # The file name

      row_count =  0
      start_time = time.time()
      num_unknown = 0  # Counter for number of unknown geocodes

      # Now read the file content and process it
      #
      for csv_line in csv_reader:

        locality_pid =  csv_line[locality_pid_col].strip()
        locality_name = csv_line[locality_name_col].strip()
        state_abbrev =  csv_line[state_abbrev_col].strip()
        postcode =      csv_line[postcode_col].strip()

        locality_pid = self.check_pid(locality_pid, 'LOCALITY_PID', row_count)

        if (locality_pid != ''):

          # Check if LOCALITY_PID has a known geocode
          #
          if (locality_pid not in locality_index['locality_geocode']):
            logging.warn('LOCALITY_PID "%s" in row %d does not have a ' % \
                         (locality_pid, row_count) + 'known geocode')
            num_unknown += 1

          else:  # Process the record

            clean_locality_name = self.clean_value(locality_name, locality_pid)
            logging.debug('Cleaned locality name "%s" into "%s"' % \
                          (locality_name, clean_locality_name))

            state_abbrev =  state_abbrev.strip().lower()

            # Check locality name and try to impute or correct
            #
            if ((clean_locality_name not in aust_post_sub_dict) or \
                (clean_locality_name == '') and \
               (postcode in aust_post_pc_dict)):

              suburb_list = aust_post_pc_dict.get(postcode, [])

              if (len(suburb_list) == 1):  # One suburb name for this postcode
                suburb_name = suburb_list[0]

                if (clean_locality_name == ''):
                  logging.info('      Impute locality name "%s" in row %d ' % \
                               (clean_locality_name, row_count) + \
                               'from Australia Post look-up table')
                else:
                  logging.info('      Replace locality name "%s" with suburb' \
                               % (clean_locality_name) + \
                               ' name "%s" in row %d' % \
                               (suburb_name, row_count))

                clean_locality_name = suburb_name  ###### Make sure this is OK

            if (clean_locality_name == ''):  # Still an empty locality name
              logging.debug('      Locality name is empty in row %d' % \
                            (row_count))
              empty_suburb_count += 1

            elif (clean_locality_name not in aust_post_sub_dict):
              logging.debug('      Locality name "%s" in row %d is not a ' % \
                            (clean_locality_name, row_count) + \
                            'known Australia Post suburb')
              unknown_suburb_count += 1

            # Check postcode and try to impute or correct
            #
            if ((postcode not in aust_post_pc_dict) or \
                (postcode == '') and
               (clean_locality_name in aust_post_sub_dict)):

              postcode_list = aust_post_sub_dict.get(clean_locality_name, [])

              if (len(postcode_list) == 1):  # One postcode for this locality
                postcode_value = postcode_list[0]

                if (postcode == ''):
                  logging.info('      Impute postcode "%s" in row %d from ' % \
                               (postcode_value, row_count) + \
                               'Australia Post look-up table')
                else:
                  logging.info('      Replace postcode "%s" with postcode ' % \
                               (postcode) + '"%s" ' % (postcode_value) + \
                               'in row %d' % (row_count))

              postcode = postcode_value  ######### Make sure this is OK #####

            if (postcode == ''):  # Still an empty postcode
              logging.debug('      Postcode is empty in row %d' % (row_count))
              empty_postcode_count += 1

            elif (postcode not in aust_post_pc_dict):
              logging.debug('      Postcode "%s" in row %d is not a known ' % \
                            (postcode, row_count) +'Australia Post postcode')
              unknown_postcode_count += 1

            # Check postcode / locality name combination
            #
            postcode_suburb_comb = postcode+'/'+clean_locality_name

            if (postcode_suburb_comb not in aust_post_pc_sub_set):
              logging.debug('      Postcode / locality name combination "%s"' \
                            % (postcode_suburb_comb) + \
                            'in row %d is not a known ' % (row_count) + \
                            'Australia Post combination')
              unknown_pc_sub_comb_count += 1

            # Now insert into inverted indices - - - - - - - - - - - - - - - -
            #
            self.insert_locality_value(locality_index['locality_name'],
                                       clean_locality_name, locality_pid)
            self.insert_locality_value(locality_index['postcode'], postcode,
                                       locality_pid)
            self.insert_locality_value(locality_index['state_abbrev'],
                                       state_abbrev, locality_pid)

        row_count += 1
        if ((row_count % 1000) == 0):  # Process meter
          used_time = (time.time() - start_time)
          logging.info('    Processed %d rows in %s (time per row: %s)' % \
                       (row_count, output.time_string(used_time), \
                       output.time_string(used_time/row_count)))

      csv_file.close()
      logging.info('    Read %d data lines' % (row_count))

      logging.info('    Locality name index contains      %d entries' % \
                   (len(locality_index['locality_name'])))
      logging.info('    Postcode index contains           %d entries' % \
                   (len(locality_index['postcode'])))
      logging.info('    State abbreviation index contains %d entries' % \
                   (len(locality_index['state_abbrev'])))
      logging.info('    Number of unknown / empty postcodes values: %d / %d' \
                   % (unknown_postcode_count, empty_postcode_count))
      logging.info('    Number of unknown / empty suburb values:    %d / %d' \
                   % (unknown_suburb_count, empty_suburb_count))
      logging.info('    Number of unknown postcode / locality name ' + \
                   'combinations: %d' % (unknown_pc_sub_comb_count))
      logging.info('    Number of records with unknown geocodes: %d' % \
                   (num_unknown))

    return (locality_index['locality_name'], locality_index['postcode'],
            locality_index['state_abbrev'], locality_index['locality_geocode'])

  # ---------------------------------------------------------------------------

  def index_street(self, g_street_file_info,
                   g_street_locality_alias_file_info,
                   g_street_locality_geocode_file_info,
                   g_address_detail_file_info):
    """Clean, standardise and index the two GNAF files containing street
       information as well as the street-locality geocode file. Use the
       address detail file to build a dictionary of the possible
       street-locality pairs first.
       - Create an inverted index for the street-locality geocode, with
         LOCALITY_PID as keys and dictionaries of STREET_PID and coordinates.
       - Use Febrl's clean_component() and tag_address_component() to clean and
         standardise the street names and types.
       - Street suffix' only values are 'n','s','e', and 'w'.
       - Create three separate inverted indices for street names, types and
         suffixes.
       - Each inverted index is a dictionary with values as keys, and a
         dictionary with LOCALITY_PIDs as keys and sets of STREET_PIDs as
         values.
       - Don't insert entries which have no geocode into the inverted indices.
       - Returns inverted indices (dictionaries) for:
           street_name
           street_type
           street_suffix
           street_locality_geocode
    """

    logging.info('Build dictionaries for "G_STREET_LOCALITY_GEOCODE" and ' + \
                 'clean, standardise and index street files "G_STREET" and ' \
                 + '"G_STREET_LOCALITY_ALIAS"')

    street_index = {'street_locality_geocode': {},
                    'street_name':             {},
                    'street_type':             {},
                    'street_suffix':           {}}

    # A dictionary with all STREET_PID/LOCALITY_PID pairs, even if they don't
    # have coordinates (so we get at least the available pairs of PIDs)
    #
    all_street_loc_gecode_dict = {}

    min_locality_pid = '999999999'  # Find the smallest LOCALITY_PID
    max_locality_pid = '0'          # Find the largest LOCALITY_PID
    min_street_pid =   '999999999'  # Find the smallest STREET_PID
    max_street_pid =   '0'          # Find the largest STREET_PID

    # Process the G_STREET_LOCALITY_GEOCODE file - - - - - - - - - - - - - - -

    # Get a dictionary of the column names
    #
    attr_dict = self.attr_list_to_dict(g_street_locality_geocode_file_info)

    # Get the column numbers
    #
    street_pid_col =   attr_dict['STREET_PID']
    locality_pid_col = attr_dict['LOCALITY_PID']
    latitude_col =     attr_dict['LATITUDE']
    longitude_col =    attr_dict['LONGITUDE']

    # Open the CSV file
    #
    [csv_reader, csv_file] = self.open_csv(g_street_locality_geocode_file_info)

    logging.info('  Process "G_STREET_LOCALITY_GEOCODE" file')

    row_count =  0
    start_time = time.time()
    num_missing = 0  # Counter for number records with missing coordinates

    # Now read the file content and process it
    #
    for csv_line in csv_reader:

      street_pid =   csv_line[street_pid_col].strip()
      locality_pid = csv_line[locality_pid_col].strip()
      latitude =     csv_line[latitude_col].strip()
      longitude =    csv_line[longitude_col].strip()

      locality_pid = self.check_pid(locality_pid, 'LOCALITY_PID', row_count)
      street_pid =   self.check_pid(street_pid, 'STREET_PID', row_count)

      if (locality_pid == '') or (street_pid == ''):
        logging.warn('LOCALITY_PID or STREET_PID is empty in row %d (will' % \
                     (row_count) + ' not be inserted into locality geocode' + \
                     ' index)')

      else:
        min_locality_pid = min(min_locality_pid, locality_pid)
        max_locality_pid = max(max_locality_pid, locality_pid)
        min_street_pid =   min(min_street_pid,   street_pid)
        max_street_pid =   max(max_street_pid,   street_pid)

        # Insert into dictionary with all pairs of STREET_PID/LOCALITY_PID - -
        #
        st_pid_set = all_street_loc_gecode_dict.get(locality_pid, sets.Set())
        if (street_pid in st_pid_set):
          logging.exception('Row %d: STREET_PID "%s" already in set of ' \
                            (row_count, street_pid) + \
                            'streets for LOCALITY_PID "%s"' % \
                            (locality_pid))
          raise Exception
        else:
          st_pid_set.add(street_pid)
          all_street_loc_gecode_dict[locality_pid] = st_pid_set

        # Only insert into geocode index if coordinates available - - - - - - -
        #
        if (latitude == '') or (longitude == ''):
          logging.warn('Missing coordinates in row %d' % (row_count))
          num_missing += 1
        else:
          coordinates = latitude+','+longitude

          # Now insert into street-locality geocode inverted index
          #
          loc_dict = street_index['street_locality_geocode'].get(locality_pid,
                                                                 {})
          if (street_pid in loc_dict):
            logging.exception('Row %d: STREET_PID "%s" already in locality ' \
                              (row_count, street_pid) + \
                              'dictionary for LOCALITY_PID "%s"' % \
                              (locality_pid))
            raise Exception
          else:
            loc_dict[street_pid] = coordinates
            street_index['street_locality_geocode'][locality_pid] = loc_dict

      row_count += 1
      if ((row_count % 10000) == 0):  # Process meter
        used_time = (time.time() - start_time)
        logging.info('    Processed %d rows in %s (time per row: %s)' % \
                     (row_count, output.time_string(used_time), \
                     output.time_string(used_time/row_count)))

    csv_file.close()
    logging.info('    Read %d data lines' % (row_count))

    logging.info('    Street-locality geocode index contains %d entries' % \
                 (len(street_index['street_locality_geocode'])))
    logging.info('    Number of records with missing coordinates: %d' % \
                 (num_missing))
    logging.info('    Smallest LOCALITY_PID: %s' % (min_locality_pid))
    logging.info('    Largest LOCALITY_PID:  %s' % (max_locality_pid))
    logging.info('    Smallest STREET_PID:   %s' % (min_street_pid))
    logging.info('    Largest STREET_PID:    %s' % (max_street_pid))

    # Process the G_ADDRESS_DETAIL file - - - - - - - - - - - - - - - - - - - -
    #
    street_locality_address_index = {}  # Index with STREET_PIDs as keys and
                                        # sets of LOCALITY_PIDs as values

    # Get a dictionary of the column names
    #
    attr_dict = self.attr_list_to_dict(g_address_detail_file_info)

    # Get the column numbers
    #
    locality_pid_col =     attr_dict['LOCALITY_PID']
    street_pid_col =       attr_dict['STREET_PID']
    gnaf_pid_col =         attr_dict['GNAF_PID']
    address_site_pid_col = attr_dict['ADDRESS_SITE_PID']

    # Open the CSV file
    #
    [csv_reader, csv_file] = self.open_csv(g_address_detail_file_info)

    logging.info('  Process "G_ADDRESS_DETAIL" file')

    row_count =  0
    start_time = time.time()
    num_missing = 0  # Counter for number records with missing PIDs

    # Now read the file content and process it
    #
    for csv_line in csv_reader:

      locality_pid =     csv_line[locality_pid_col].strip()
      street_pid =       csv_line[street_pid_col].strip()
      gnaf_pid =         csv_line[gnaf_pid_col].strip()
      address_site_pid = csv_line[address_site_pid_col].strip()

      locality_pid =     self.check_pid(locality_pid, 'LOCALITY_PID',
                                        row_count)
      street_pid =       self.check_pid(street_pid, 'STREET_PID', row_count)
      address_site_pid = self.check_pid(address_site_pid, 'ADDRES_SITE_PID',
                                        row_count)

      # Only process street and locality PIDs if all PIDs are available
      #
      if ((gnaf_pid != '') and (address_site_pid != '') and \
          (locality_pid != '') and (street_pid != '')):

        # Check if street-locality PID pair is in corresponding dictionary - -
        #
        if (locality_pid in all_street_loc_gecode_dict):
          if (street_pid in all_street_loc_gecode_dict[locality_pid]):

            # Insert into street-locality address index
            #
            loc_set = street_locality_address_index.get(street_pid, sets.Set())
            loc_set.add(locality_pid)
            street_locality_address_index[street_pid] = loc_set

      else:
        num_missing += 1

      row_count += 1
      if ((row_count % 10000) == 0):  # Process meter
        used_time = (time.time() - start_time)
        logging.info('    Processed %d rows in %s (time per row: %s)' % \
                     (row_count, output.time_string(used_time), \
                     output.time_string(used_time/row_count)))

    csv_file.close()

    logging.info('    Read %d data lines' % (row_count))
    logging.info('  Number of STREET_PIDs with LOCALITY_PIDs: %d' % \
                 (len(street_locality_address_index)))
    logging.info('    Number of records with missing PIDs: %d' % (num_missing))

    # Process the G_STREET and G_STREET_LOCALITY_ALIAS files - - - - - - - - -

    for fi in [g_street_file_info, g_street_locality_alias_file_info]:

      if (fi == g_street_locality_alias_file_info):
        do_locality = True  # A flag for processing the locality PIDs as well
      else:
        do_locality = False

      # Get a dictionary of the column names
      #
      attr_dict = self.attr_list_to_dict(fi)

      # Get the column numbers
      #
      street_name_col =   attr_dict['STREET_NAME']
      street_type_col =   attr_dict['STREET_TYPE']
      street_suffix_col = attr_dict['STREET_SUFFIX']
      street_pid_col =    attr_dict['STREET_PID']
      if (do_locality == True):
        locality_pid_col = attr_dict['LOCALITY_PID']

      # Open the CSV file
      #
      [csv_reader, csv_file] = self.open_csv(fi)

      logging.info('  Process "%s" file' % (fi[0]))  # The file name

      row_count =  0
      start_time = time.time()
      num_unknown = 0  # Counter for number of unknown geocodes

      # Now read the file content and process it
      #
      for csv_line in csv_reader:

        street_pid =    csv_line[street_pid_col].strip()
        street_name =   csv_line[street_name_col].strip()
        street_type =   csv_line[street_type_col].strip()
        street_suffix = csv_line[street_suffix_col].strip()

        street_pid = self.check_pid(street_pid, 'STREET_PID', row_count)

        if (do_locality == True):
          locality_pid = csv_line[locality_pid_col].strip()
          locality_pid = self.check_pid(locality_pid, 'LOCALITY_PID',
                                        row_count)

        else:  # Get the set of LOCALITY_PIDs for this STREET_PID

          locality_pid = street_locality_address_index.get(street_pid, '')

        # Don't insert record if STREET_PID or LOCALITY_PID is empty - - - - -
        #
        if ((street_pid == '') or (locality_pid == '')):
          logging.warn('STREET_PID and/or LOCALITY_PID are empty in row %d' % \
                       (row_count) + \
                       ' (record will not be inserted into inverted index)')

        else:  # Process the record

          clean_street_name = self.clean_value(street_name, street_pid)
          logging.debug('Cleaned street name "%s" into "%s"' % (street_name, \
                        clean_street_name))

          clean_street_type = self.clean_value(street_type, street_pid)
          logging.debug('Cleaned street type "%s" into "%s"' % (street_type, \
                        clean_street_type))

          street_suffix = street_suffix.strip().lower()

          # Now insert into inverted indices - - - - - - - - - - - - - - - - -
          #
          self.insert_street_value(street_index['street_name'],
                                   clean_street_name, street_pid,
                                   locality_pid, do_locality)
          self.insert_street_value(street_index['street_type'],
                                   clean_street_type, street_pid,
                                   locality_pid, do_locality)
          self.insert_street_value(street_index['street_suffix'],
                                   street_suffix, street_pid,
                                   locality_pid, do_locality)

        row_count += 1
        if ((row_count % 1000) == 0):  # Process meter
          used_time = (time.time() - start_time)
          logging.info('    Processed %d rows in %s (time per row: %s)' % \
                       (row_count, output.time_string(used_time), \
                       output.time_string(used_time/row_count)))

      csv_file.close()
      logging.info('    Read %d data lines' % (row_count))

      logging.info('    Street name index contains   %d entries' % \
                   (len(street_index['street_name'])))
      logging.info('    Street type index contains   %d entries' % \
                   (len(street_index['street_type'])))
      logging.info('    Street suffix index contains %d entries' % \
                   (len(street_index['street_suffix'])))
      logging.info('    Number of records with unknown geocodes: %d' % \
                   (num_unknown))

    return (street_index['street_name'], street_index['street_type'],
            street_index['street_suffix'],
            street_index['street_locality_geocode'])

  # ---------------------------------------------------------------------------

  def index_address(self, g_address_detail_file_info,
                    g_address_site_geocode_file_info,
                    g_address_alias_file_info):
    """Clean, standardise and index the GNAF file containing address details
       information as well as the address site geocode file.
       - Create an inverted index for the address site geocode, with
         ADDRESS_SITE_PIDs as keys and coordinates as values.
       - Create a dictionary for the address alias file. This will be used to
         add address alias PIDs into the inverted indices.
       - Use Febrl's clean_component() and tag_address_component() to clean and
         standardise various attributes.
       - Each inverted index is a dictionary with values as keys, and a
         dictionary with LOCALITY_PIDs as keys and a dictionary with
         STREET_PIDs as keys and sets of ADDRESS_SITE_PIDs as values.
       - Don't insert entries which have no geocode into the inverted indices.
       - Create inverted indices for:
           flat_number_prefix
           flat_number
           flat_number_suffix
           flat_type
           level_number
           level_type
           building_name
           location_descr
           number_first_prefix
           number_first
           number_first_suffix
           number_last_prefix
           number_last
           number_last_suffix
           lot_number_prefix
           lot_number
           lot_number_suffix
           address_site_geocode
    """

    logging.info('Build dictionaries for "G_ADDRESS_SITE_GEOCODE" and ' + \
                 '"G_ADDRESS_ALIAS", and clean, standardise and index ' + \
                 'address file "G_ADDRESS_DETAIL"')

    # Two dictionaries with the principal and alias GNAF-PIDs
    #
    address_principal_alias_index = {}  # Alias PIDs for principal PID
    address_alias_principal_index = {}  # Principal PIDs for alias PID

    # Dictionaries for the geocode and inverted indices
    #
    address_index = {'address_site_geocode': {},
                     'flat_number_prefix':   {},
                     'flat_number':          {},
                     'flat_number_suffix':   {},
                     'flat_type':            {},
                     'level_number':         {},
                     'level_type':           {},
                     'building_name':        {},
                     'location_descr':       {},
                     'number_first_prefix':  {},
                     'number_first':         {},
                     'number_first_suffix':  {},
                     'number_last_prefix':   {},
                     'number_last':          {},
                     'number_last_suffix':   {},
                     'lot_number_prefix':    {},
                     'lot_number':           {},
                     'lot_number_suffix':    {}}

    min_address_site_pid = '999999999'  # Find the smallest ADDRESS_SITE_PID
    max_address_site_pid = '0'          # Find the largest ADDRESS_SITE_PID

    # Process the G_ADDRESS_ALIAS file - - - - - - - - - - - - - - - - - - - -

    # Get a dictionary of the column names
    #
    attr_dict = self.attr_list_to_dict(g_address_alias_file_info)

    # Get the column numbers
    #
    principal_pid_col = attr_dict['PRINCIPAL_PID']
    alias_pid_col =     attr_dict['ALIAS_PID']
    alias_type_col =    attr_dict['ALIAS_TYPE']  # Currently not really needed

    # Open the CSV file
    #
    [csv_reader, csv_file] = self.open_csv(g_address_alias_file_info)

    logging.info('  Process "G_ADDRESS_ALIAS" file')

    row_count =  0
    start_time = time.time()

    warn_principal_alias = 0
    warn_alias_principal = 0

    # Now read the file content and process it
    #
    for csv_line in csv_reader:

      principal_pid = csv_line[principal_pid_col].strip()
      alias_pid =     csv_line[alias_pid_col].strip()
      alias_type =    csv_line[alias_type_col].strip()

      principal_pid = self.check_pid(principal_pid, 'principal GNAF_PID',
                                     row_count)
      alias_pid = self.check_pid(alias_pid, 'alias GNAF_PID', row_count)


      if (principal_pid == '') or (alias_pid == ''):
        logging.warn('Principal or alias GNAF_PID is empty in row %d' % \
                     (row_count) + \
                     ' (will not be inserted into address alias index)')

      else:

        # Insert alias PID into set for principal PID dictionary
        #
        alias_pid_set = address_principal_alias_index.get(principal_pid,
                                                          sets.Set())
        if (alias_pid in alias_pid_set):
          logging.warn('Alias PID "%s" already in set for principal PID ' % \
                       (alias_pid) + '"%s" in row %d' % \
                       (principal_pid, row_count))
          warn_principal_alias += 1
        else:
          alias_pid_set.add(alias_pid)
          address_principal_alias_index[principal_pid] = alias_pid_set

        # Insert principal PID into set for alias PID dictionary
        #
        principal_pid_set = address_alias_principal_index.get(alias_pid,
                                                              sets.Set())
        if (principal_pid in principal_pid_set):
          logging.warn('Principal PID "%s" already in set for alias PID ' % \
                       (principal_pid) + '"%s" in row %d' % \
                       (alias_pid, row_count))
          warn_alias_principal += 1
        else:
          principal_pid_set.add(principal_pid)
          address_alias_principal_index[alias_pid] = principal_pid_set

      row_count += 1
      if ((row_count % 100000) == 0):  # Process meter
        used_time = (time.time() - start_time)
        logging.info('    Processed %d rows in %s (time per row: %s)' % \
                     (row_count, output.time_string(used_time), \
                     output.time_string(used_time/row_count)))

    csv_file.close()
    logging.info('  Read %d data lines' % (row_count))

    logging.info('Number of warnings %d / %d' % \
                 (warn_principal_alias, warn_alias_principal))

    logging.info('  Address alias PIDs for principal PID index contains ' + \
                 '%d entries' % (len(address_principal_alias_index)))
    logging.info('  Address principal PIDs for alias PID index contains ' + \
                 '%d entries' % (len(address_alias_principal_index)))

    # Process the G_ADDRESS_DETAIL file - - - - - - - - - - - - - - - - - - - -
    # (first time only to get a dictionary of all GNAF_PIDs (as key) with their
    # corresponding ADDRESS_SITE_PID (as value))
    #
    gnaf_address_site_index = {}

    # Get a dictionary of the column names
    #
    attr_dict = self.attr_list_to_dict(g_address_detail_file_info)

    # Get the column numbers
    #
    gnaf_pid_col =            attr_dict['GNAF_PID']
    address_site_pid_col =    attr_dict['ADDRESS_SITE_PID']

    # Open the CSV file
    #
    [csv_reader, csv_file] = self.open_csv(g_address_detail_file_info)

    logging.info('  Process "G_ADDRESS_DETAIL" file')

    row_count =  0
    start_time = time.time()
    num_unknown = 0  # Counter for number of unknown address site geocodes

    # Now read the file content and process it
    #
    for csv_line in csv_reader:

      gnaf_pid =            csv_line[gnaf_pid_col].strip()
      address_site_pid =    csv_line[address_site_pid_col].strip()

      gnaf_pid = self.check_pid(gnaf_pid, 'GNAF_PID', row_count)
      address_site_pid = self.check_pid(address_site_pid, 'ADDRESS_SITE_PID',
                                        row_count)

      if (gnaf_pid != '') and (address_site_pid != ''):

        if (gnaf_pid in gnaf_address_site_index):
          logging.warn('GNAF_PID "%s" already in index' % (gnaf_pid))

        address_site_set = gnaf_address_site_index.get(gnaf_pid, sets.Set())

        if (address_site_pid in address_site_set):
          logging.warn('ADRESS_SITE_PID "%s" already in address site set ' % \
                       (address_site_pid) + 'of GNAF_PID "%s" in row %d' \
                       % (gnaf_pid, row_count))
        else:
          address_site_set.add(address_site_pid)
          gnaf_address_site_index[gnaf_pid] = address_site_set

      row_count += 1
      if ((row_count % 100000) == 0):  # Process meter
        used_time = (time.time() - start_time)
        logging.info('    Processed %d rows in %s (time per row: %s)' % \
                     (row_count, output.time_string(used_time), \
                     output.time_string(used_time/row_count)))

    csv_file.close()
    logging.info('    Read %d data lines' % (row_count))

    logging.info('    G-NAF / Address site index contains %d entries' % \
                 (len(gnaf_address_site_index)))

    # Process the G_ADDRESS_SITE_GEOCODE file - - - - - - - - - - - - - - - - -

    # Get a dictionary of the column names
    #
    attr_dict = self.attr_list_to_dict(g_address_site_geocode_file_info)

    # Get the column numbers
    #
    address_site_pid_col = attr_dict['ADDRESS_SITE_PID']
    latitude_col =         attr_dict['LATITUDE']
    longitude_col =        attr_dict['LONGITUDE']

    # Open the CSV file
    #
    [csv_reader, csv_file] = self.open_csv(g_address_site_geocode_file_info)

    logging.info('  Process "G_ADDRESS_SITE_GEOCODE" file')

    row_count =  0
    start_time = time.time()
    num_missing = 0  # Counter for number records with missing coordinates

    # Now read the file content and process it
    #
    for csv_line in csv_reader:

      address_site_pid = csv_line[address_site_pid_col].strip()
      latitude =         csv_line[latitude_col].strip()
      longitude =        csv_line[longitude_col].strip()

      address_site_pid = self.check_pid(address_site_pid, 'ADDRESS_SITE_PID',
                                        row_count)

      if (address_site_pid != ''):
        min_address_site_pid = min(min_address_site_pid, address_site_pid)
        max_address_site_pid = max(max_address_site_pid, address_site_pid)

        if (latitude == '') or (longitude == ''):  # Coordinates missing
          logging.warn('Missing coordinates in row %d' % (row_count))
          num_missing += 1
        else:
          coordinates = latitude+','+longitude

          # Now insert into dictionary
          #
          if (address_site_pid in address_index['address_site_geocode']):
            logging.exception('Row %d: ADDRESS_SITE_PID "%s" already in ' % \
                              (row_count, \
                   address_index['address_site_geocode'][address_site_pid]) + \
                              'address site geocode index')
            raise Exception
          else:
            address_index['address_site_geocode'][address_site_pid] = \
              coordinates

      row_count += 1
      if ((row_count % 10000) == 0):  # Process meter
        used_time = (time.time() - start_time)
        logging.info('    Processed %d rows in %s (time per row: %s)' % \
                     (row_count, output.time_string(used_time), \
                     output.time_string(used_time/row_count)))

    csv_file.close()
    logging.info('    Read %d data lines' % (row_count))

    logging.info('    Address site geocode index contains %d entries' % \
                 (len(address_index['address_site_geocode'])))
    logging.info('    Number of records with missing coordinates: %d' % \
                 (num_missing))
    logging.info('    Smallest ADDRESS_SITE_PID: %s' % (min_address_site_pid))
    logging.info('    Largest  ADDRESS_SITE_PID: %s' % (max_address_site_pid))

    # Process the G_ADDRESS_DETAIL file - - - - - - - - - - - - - - - - - - - -
    # (second time process all fields and insert them into inverted indices,
    # also check for alias PIDs)

    # Get a dictionary of the column names
    #
    attr_dict = self.attr_list_to_dict(g_address_detail_file_info)

    # Get the column numbers
    #
    flat_number_prefix_col =  attr_dict['FLAT_NUMBER_PREFIX']
    flat_number_col =         attr_dict['FLAT_NUMBER']
    flat_number_suffix_col =  attr_dict['FLAT_NUMBER_SUFFIX']
    flat_type_col =           attr_dict['FLAT_TYPE']
    level_number_col =        attr_dict['LEVEL_NUMBER_PREFIX']
    level_type_col =          attr_dict['LEVEL_TYPE']
    building_name_col =       attr_dict['BUILDING_NAME']
    location_descr_col =      attr_dict['LOCATION_DESCRIPTION']
    number_first_prefix_col = attr_dict['NUMBER_FIRST_PREFIX']
    number_first_col =        attr_dict['NUMBER_FIRST']
    number_first_suffix_col = attr_dict['NUMBER_FIRST_SUFFIX']
    number_last_prefix_col =  attr_dict['NUMBER_LAST_PREFIX']
    number_last_col =         attr_dict['NUMBER_LAST']
    number_last_suffix_col =  attr_dict['NUMBER_LAST_SUFFIX']
    lot_number_prefix_col =   attr_dict['LOT_NUMBER_PREFIX']
    lot_number_col =          attr_dict['LOT_NUMBER']
    lot_number_suffix_col =   attr_dict['LOT_NUMBER_SUFFIX']
    gnaf_pid_col =            attr_dict['GNAF_PID']
    locality_pid_col =        attr_dict['LOCALITY_PID']
    street_pid_col =          attr_dict['STREET_PID']
    address_site_pid_col =    attr_dict['ADDRESS_SITE_PID']

    # Open the CSV file
    #
    [csv_reader, csv_file] = self.open_csv(g_address_detail_file_info)

    logging.info('  Process "G_ADDRESS_DETAIL" file')

    row_count =  0
    start_time = time.time()
    num_unknown = 0  # Counter for number of unknown address site geocodes
    num_missing = 0  # Counter for number of records with missing PIDs
    num_alias =   0  # Counter for number of records with corresponding aliases

    # Now read the file content and process it
    #
    for csv_line in csv_reader:

      gnaf_pid =            csv_line[gnaf_pid_col].strip()
      locality_pid =        csv_line[locality_pid_col].strip()
      street_pid =          csv_line[street_pid_col].strip()
      address_site_pid =    csv_line[address_site_pid_col].strip()
      flat_number_prefix =  csv_line[flat_number_prefix_col].strip()
      flat_number =         csv_line[flat_number_col].strip()
      flat_number_suffix =  csv_line[flat_number_suffix_col].strip()
      flat_type =           csv_line[flat_type_col].strip()
      level_number =        csv_line[level_number_col].strip()
      level_type =          csv_line[level_type_col].strip()
      building_name =       csv_line[building_name_col].strip()
      location_descr =      csv_line[location_descr_col].strip()
      number_first_prefix = csv_line[number_first_prefix_col].strip()
      number_first =        csv_line[number_first_col].strip()
      number_first_suffix = csv_line[number_first_suffix_col].strip()
      number_last_prefix =  csv_line[number_last_prefix_col].strip()
      number_last =         csv_line[number_last_col].strip()
      number_last_suffix =  csv_line[number_last_suffix_col].strip()
      lot_number_prefix =   csv_line[lot_number_prefix_col].strip()
      lot_number =          csv_line[lot_number_col].strip()
      lot_number_suffix =   csv_line[lot_number_suffix_col].strip()

      locality_pid =     self.check_pid(locality_pid, 'LOCALITY_PID',
                                        row_count)
      street_pid =       self.check_pid(street_pid, 'STREET_PID', row_count)
      gnaf_pid =         self.check_pid(gnaf_pid, 'GNAF_PID', row_count)
      address_site_pid = self.check_pid(address_site_pid, 'ADDRESS_SITE_PID',
                                        row_count)

      if (address_site_pid != ''):

        # Check if ADDRESS_SITE_PID has a known geocode
        #
        if (address_site_pid not in address_index['address_site_geocode']):
          logging.warn('ADDRESS_SITE_PID "%s" in row %d does not have a ' % \
                       (address_site_pid, row_count) + \
                       'known address site geocode')
          address_site_pid == '' # Set it to an empty value so it not used ####
          num_unknown += 1

      # Don't insert record if ADDRESS_SITE_PID and/or LOCALITY_PID and/or
      # STREET_PID are empty
      #
      if ((address_site_pid == '') or (locality_pid == '') or \
          (street_pid == '')):
        num_missing += 1
        logging.warn('ADDRESS_SITE_PID and/or LOCALITY_PID and/or STREET_PID' \
                     + ' are empty in row %d' % (row_count) + \
                     ' (record will not be inserted into inverted index)')

      else:  # Process the record

        #clean_flat_number_prefix = flat_number_prefix.lower()
        clean_flat_number_prefix = self.clean_value(flat_number_prefix,
                                                    gnaf_pid)
        #clean_flat_number = flat_number.lower()
        clean_flat_number = self.clean_value(flat_number, gnaf_pid)
        #clean_flat_number_suffix = flat_number_suffix.lower()
        clean_flat_number_suffix = self.clean_value(flat_number_suffix,
                                                    gnaf_pid)
        clean_flat_type =     self.clean_value(flat_type, gnaf_pid)

        #level_number = level_number.lower()
        clean_level_number = self.clean_value(level_number, gnaf_pid)
        clean_level_type =   self.clean_value(level_type, gnaf_pid)

        #number_first_prefix = number_first_prefix.lower()
        clean_number_first_prefix = self.clean_value(number_first_prefix,
                                                     gnaf_pid)
        #number_first = number_first.lower()
        clean_number_first = self.clean_value(number_first, gnaf_pid)
        #number_first_suffix = number_first_suffix.lower()
        clean_number_first_suffix = self.clean_value(number_first_suffix,
                                                     gnaf_pid)

        #number_last_prefix = number_last_prefix.lower()
        clean_number_last_prefix = self.clean_value(number_last_prefix,
                                                    gnaf_pid)
        #number_last = number_last.lower()
        clean_number_last = self.clean_value(number_last, gnaf_pid)
        #number_last_suffix = number_last_suffix.lower()
        clean_number_last_suffix = self.clean_value(number_last_suffix,
                                                    gnaf_pid)

        #clean_lot_number_prefix = lot_number_prefix.lower()
        clean_lot_number_prefix = self.clean_value(lot_number_prefix,
                                                   gnaf_pid)
        #clean_lot_number = lot_number.lower()
        clean_lot_number = self.clean_value(lot_number, gnaf_pid)
        #clean_lot_number_suffix = lot_number_suffix.lower()
        clean_lot_number_suffix = self.clean_value(lot_number_suffix,
                                                   gnaf_pid)

        clean_building_name = self.clean_value(building_name, gnaf_pid)
        clean_location_descr = self.clean_value(location_descr, gnaf_pid)

        # Check if there are alias GNAF_PIDs
        #
        if (gnaf_pid in address_principal_alias_index):
          num_alias += 1
          alias_gnaf_pid_set = address_principal_alias_index[gnaf_pid]
          #logging.info('== %s' % (str(alias_gnaf_pid_set))) #################
          address_site_pid_set = sets.Set()
          for pid in alias_gnaf_pid_set:
            #logging.info('  * %s' % (str(pid))) ###################3
            if (pid in gnaf_address_site_index):
              address_site_pid_set = \
                address_site_pid_set.union(gnaf_address_site_index[pid])

          #logging.info('** gnaf_pid/address_site_pid: %s / %s' % \
          #             (gnaf_pid, address_site_pid)) ###########
          #logging.info('** gnaf_pid_list:             %s' % \
          #             (str(alias_gnaf_pid_list))) ###############
          if (address_site_pid in address_site_pid_set):
            logging.warn('ADDRESS_SITE_PID already in alias list in row %d' \
                         % (row_count))
          else:
            address_site_pid_set.add(address_site_pid)

        else:
          address_site_pid_set = sets.Set([address_site_pid])

        #logging.info('** address_site_pid_list:     %s' % \
        #             (str(address_site_pid_list))) ###############

        # Now insert into inverted indices - - - - - - - - - - - - - - - - - -
        #
        self.insert_address_value(address_index['flat_number_prefix'],
                                  clean_flat_number_prefix, street_pid,
                                  locality_pid, address_site_pid_set)

        self.insert_address_value(address_index['flat_number'],
                                  clean_flat_number, street_pid, locality_pid,
                                  address_site_pid_set)
        self.insert_address_value(address_index['flat_number_suffix'],
                                  clean_flat_number_suffix, street_pid,
                                  locality_pid, address_site_pid_set)
        self.insert_address_value(address_index['flat_type'], clean_flat_type,
                                  street_pid, locality_pid,
                                  address_site_pid_set)

        self.insert_address_value(address_index['level_number'],
                                  clean_level_number, street_pid, locality_pid,
                                  address_site_pid_set)
        self.insert_address_value(address_index['level_type'],
                                  clean_level_type, street_pid, locality_pid,
                                  address_site_pid_set)

        self.insert_address_value(address_index['number_first_prefix'],
                                  clean_number_first_prefix, street_pid,
                                  locality_pid, address_site_pid_set)
        self.insert_address_value(address_index['number_first'],
                                  clean_number_first, street_pid, locality_pid,
                                  address_site_pid_set)
        self.insert_address_value(address_index['number_first_suffix'],
                                  clean_number_first_suffix, street_pid,
                                  locality_pid, address_site_pid_set)

        self.insert_address_value(address_index['number_last_prefix'],
                                  clean_number_last_prefix, street_pid,
                                  locality_pid, address_site_pid_set)
        self.insert_address_value(address_index['number_last'],
                                  clean_number_last, street_pid, locality_pid,
                                  address_site_pid_set)
        self.insert_address_value(address_index['number_last_suffix'],
                                  clean_number_last_suffix, street_pid,
                                  locality_pid, address_site_pid_set)

        self.insert_address_value(address_index['lot_number_prefix'],
                                  clean_lot_number_prefix, street_pid,
                                  locality_pid, address_site_pid_set)
        self.insert_address_value(address_index['lot_number'],
                                  clean_lot_number, street_pid, locality_pid,
                                  address_site_pid_set)
        self.insert_address_value(address_index['lot_number_suffix'],
                                  clean_lot_number_suffix, street_pid,
                                  locality_pid, address_site_pid_set)

        self.insert_address_value(address_index['building_name'],
                                  clean_building_name, street_pid,
                                  locality_pid, address_site_pid_set)
        self.insert_address_value(address_index['location_descr'],
                                  clean_location_descr, street_pid,
                                  locality_pid, address_site_pid_set)

      row_count += 1
      if ((row_count % 10000) == 0):  # Process meter
        used_time = (time.time() - start_time)
        logging.info('    Processed %d rows in %s (time per row: %s)' % \
                     (row_count, output.time_string(used_time), \
                     output.time_string(used_time/row_count)))

    csv_file.close()
    logging.info('    Read %d data lines' % (row_count))

    logging.info('    Address site geocode index contains   %d entries' % \
                 (len(address_index['address_site_geocode'])))
    logging.info('    Flat number prefix index contains     %d entries' % \
                 (len(address_index['flat_number_prefix'])))
    logging.info('    Flat number index contains            %d entries' % \
                 (len(address_index['flat_number'])))
    logging.info('    Flat number suffix index contains     %d entries' % \
                 (len(address_index['flat_number_suffix'])))
    logging.info('    Flat type index contains              %d entries' % \
                 (len(address_index['flat_type'])))
    logging.info('    Level number index contains           %d entries' % \
                 (len(address_index['level_number'])))
    logging.info('    Level type index contains             %d entries' % \
                 (len(address_index['level_type'])))
    logging.info('    Building name index contains          %d entries' % \
                 (len(address_index['building_name'])))
    logging.info('    Location description index contains   %d entries' % \
                 (len(address_index['location_descr'])))
    logging.info('    Number first prefix index contains    %d entries' % \
                 (len(address_index['number_first_prefix'])))
    logging.info('    Number first index contains           %d entries' % \
                 (len(address_index['number_first'])))
    logging.info('    Number first suffix index contains    %d entries' % \
                 (len(address_index['number_first_suffix'])))
    logging.info('    Number last prefix index contains     %d entries' % \
                 (len(address_index['number_last_prefix'])))
    logging.info('    Number last index contains            %d entries' % \
                 (len(address_index['number_last'])))
    logging.info('    Number last suffix index contains     %d entries' % \
                 (len(address_index['number_last_suffix'])))
    logging.info('    Lot number prefix index contains      %d entries' % \
                 (len(address_index['lot_number_prefix'])))
    logging.info('    Lot number index contains             %d entries' % \
                 (len(address_index['lot_number'])))
    logging.info('    Lot number suffix index contains      %d entries' % \
                 (len(address_index['lot_number_suffix'])))
    logging.info('    Number of records with unknown address site ' + \
                 'geocodes: %d' % (num_unknown))
    logging.info('    Number of records with missing PIDs:    %d ' % \
                 (num_missing))
    logging.info('    Number of records with alias addresses: %d' % \
                 (num_alias))

    return (address_index['flat_number_prefix'], address_index['flat_number'],
            address_index['flat_number_suffix'], address_index['flat_type'],
            address_index['level_number'], address_index['level_type'],
            address_index['building_name'], address_index['location_descr'],
            address_index['number_first_prefix'],
            address_index['number_first'],
            address_index['number_first_suffix'],
            address_index['number_last_prefix'], address_index['number_last'],
            address_index['number_last_suffix'],
            address_index['lot_number_prefix'], address_index['lot_number'],
            address_index['lot_number_suffix'],
            address_index['address_site_geocode'])

  # ---------------------------------------------------------------------------

  def index_collect_dist(self, address_cd_file_info,
                         street_locality_cd_file_info, locality_cd_file_info):
    """Index the the NSW Census Collection District boundaries files.
       - Create one inverted index with PIDs as keys and CD PIDs as values.
       - The key PIDs are one of:
         - ADDRESS_SITE_PID
         - STREET_PID/LOCALITY_PID
         - LOCALITY_PID
       - Don't insert entries which have no valid collection district (i.e.
         which have value 0)
    """

    logging.info('Build an index for collection district files "ADDRESS_CD",' \
                 + ' "STREET_LOCALITY_CD" and "LOCALITY_CD"')

    cd_index = {}

    # Counters for number of missing cd PIDs in the three files
    #
    missing_cd_count = [0,0,0]

    min_cd_pid = '999999999'  # Find the smallest CD PID
    max_cd_pid = '0'          # Find the largest CD PID

    # Process the three files - - - - - - - - - - - - - - - - - - - - - - - - -

    fi_num = 0

    for fi in [address_cd_file_info, street_locality_cd_file_info, \
               locality_cd_file_info]:

      # Get a dictionary of the column names
      #
      attr_dict = self.attr_list_to_dict(fi)

      # Get the column numbers
      #
      address_site_pid_col = attr_dict.get('ADDRESS_SITE_PID', None)
      street_pid_col =       attr_dict.get('STREET_PID', None)
      locality_pid_col =     attr_dict.get('LOCALITY_PID', None)
      cd_col =               attr_dict['CD']

      # Open the CSV file
      #
      [csv_reader, csv_file] = self.open_csv(fi)

      logging.info('  Process "%s" file' % (fi[0]))  # The file name

      row_count =  0
      start_time = time.time()

      # Now read the file content and process it
      #
      for csv_line in csv_reader:

        if (address_site_pid_col != None):
          address_site_pid = csv_line[address_site_pid_col].strip()

          if (address_site_pid == ''):
            logging.warn('ADDRESS_SITE_PID is empty in row %d' % (row_count) \
                         + ' (record will not be inserted into inverted index)')
          else:
            # Make sure there is no '.00' at the end of the ADDRESS_SITE_PID
            #
            if (address_site_pid[-3] == '.'):
              address_site_pid = address_site_pid[:-3]  # Remove '.00'
        else:
          address_site_pid = ''

        if (street_pid_col != None):
          street_pid = csv_line[street_pid_col].strip()

          if (street_pid == ''):
            logging.warn('STREET_PID is empty in row %d ' % (row_count) \
                         + '(record will not be inserted into inverted index)')
          else:
            # Make sure there is no '.00' at the end of the STREET_PID
            #
            if (street_pid[-3] == '.'):
              street_pid = street_pid[:-3]  # Remove '.00'
        else:
          street_pid = ''

        if (locality_pid_col != None):
          locality_pid = csv_line[locality_pid_col].strip()

          if (locality_pid == ''):
            logging.warn('LOCALITY_PID is empty in row %d ' % (row_count) \
                         + '(record will not be inserted into inverted index)')
          else:
            # Make sure there is no '.00' at the end of the LOCALITY_PID
            #
            if (locality_pid[-3] == '.'):
              locality_pid = locality_pid[:-3]  # Remove '.00'
        else:
          locality_pid = ''

        cd_val = csv_line[cd_col].strip()

        if (cd_val == ''):
          missing_cd_count[fi_num] = missing_cd_count[fi_num] + 1

          logging.warn('Missing CD value in row %d' % (row_count))

        elif (cd_val == '0'):
          missing_cd_count[fi_num] = missing_cd_count[fi_num] + 1

          logging.warn('CD value "0" in row %d' % (row_count))

        else:  # A valid CD, insert into index

          min_cd_pid = min(min_cd_pid, cd_val)
          max_cd_pid = max(max_cd_pid, cd_val)

          # Compile the PID key
          #
          if (address_site_pid != ''):
            pid = address_site_pid
          elif (street_pid != '') and (locality_pid != ''):
            pid = street_pid+'/'+locality_pid
          elif (locality_pid != ''):
            pid = locality_pid
          else:
            logging.warn('No valid PIDs in row %d' % (row_count))

          if (pid != ''):

            if (pid in cd_index):
              logging.warn('PID "%s" already in CD index with CD value: ' % \
                           (pid) + '"%s" (new value: "%s")' % \
                           (cd_index[pid], cd_val))
            else:
              cd_index[pid] = cd_val

        row_count += 1
        if ((row_count % 10000) == 0):  # Process meter
          used_time = (time.time() - start_time)
          logging.info('    Processed %d rows in %s (time per row: %s)' % \
                       (row_count, output.time_string(used_time), \
                       output.time_string(used_time/row_count)))

      fi_num += 1

    csv_file.close()
    logging.info('    Read %d data lines' % (row_count))

    logging.info('    CD index contains %d entries' % (len(cd_index)))
    logging.info('    Smallest CD value: %s' % (min_cd_pid))
    logging.info('    Largest CD value:  %s' % (max_cd_pid))

    logging.info('    Number of records with missing CD value in files')
    logging.info('      %d / %d / %d' % (missing_cd_count[0], \
                 missing_cd_count[1], missing_cd_count[2]))

    return cd_index

  # ---------------------------------------------------------------------------

  def dict_to_shelve(self, file_name, index_dict):
    """Copy a dictionary into a shelve (file) (into the output directory) for
       persistant storage.  Add a shelve file extension to the given file name.
    """

    # Open the shelve file for writing
    #
    shelve_file_name = self.output_directory+file_name+self.shelve_file_ext

    try:
      shelve_file = shelve.open(shelve_file_name, protocol=0) # -1
    except IOError:
      logging.exception('Cannot open shelve file "%s" for writing' % \
                        (shelve_file_name))
      raise IOError
    except:
      raise Exception

    logging.info('Write to shelve file "%s"' % (shelve_file_name))

    s = time.time()

    shelve_file.clear()  # Make sure the shelve is empty

    for key in index_dict:
      vals = index_dict[key]
      shelve_file[key] = vals  # Copy into shelve

    logging.info('  Time used: %f' % (time.time() - s))

    logging.info('  Saved shelved dictionary with %d entries' % \
                 (len(shelve_file)))

    shelve_file.close()

  # ---------------------------------------------------------------------------

  def dict_to_text(self, file_name, index_dict):
    """Write a dictionary into a text file (into the output directory) for
       persistant storage. Add a text file extension to the given file name.
    """

    # Open the text file for writing
    #
    text_file_name = self.output_directory+file_name+'.txt'

    try:
      text_file = open(text_file_name, 'w')
    except IOError:
      logging.exception('Cannot open text file "%s" for writing' % \
                        (text_file_name))
      raise IOError
    except:
      raise Exception

    logging.info('Write to text file "%s"' % (text_file_name))

    s = time.time()

    for key in index_dict:
      vals = index_dict[key]
      text_file.write('%s: %s' % (str(key),str(vals))+os.linesep)

    logging.info('  Time used: %f' % (time.time() - s))

    logging.info('  Saved dictionary into text file with %d entries' % \
                 (len(index_dict)))

    text_file.close()

  # ---------------------------------------------------------------------------

  def text_to_dict(self, file_name):
    """Read a text file containing a dictionary from the output directory and
       return it. Add a text file extension to the given file name.
    """

    Set = sets.Set  # Shorthand, so 'eval()' doesn't crash

    # Open the text file
    #
    text_file_name = self.output_directory+file_name+'.txt'

    try:
      text_file = open(text_file_name, 'r')
    except IOError:
      logging.exception('Cannot open text file "%s" for reading' % \
                        (text_file_name))
      raise IOError
    except:
      raise Exception

    s = time.time()

    logging.info('Read from text file "%s"' % (text_file_name))

    text_dict = {}
    line_count = 0

    for line in text_file:
      col_index = line.find(':')  # Find first :
      if (col_index < 0):
        logging.exception('Illegal format in line %d: %s' % (line_count, line))
        raise Exception

      line_key = line[:col_index].strip()
      line_val = line[(col_index+1):].strip()

      eval_val = eval(line_val)

      #logging.info('             %s' % (str(line_key)))
      #logging.info(' line %d, content type: %s' % \
      #             (line_count, str(type(eval_val))))
      #logging.info('             %s' % (str(eval_val)))

      text_dict[line_key] = eval_val

      line_count += 1

    text_file.close()

    logging.info('  Time used: %f' % (time.time() - s))

    logging.info('  Loaded dictionary from text file with %d entries' % \
                 (len(text_dict)))

    return text_dict

  # ---------------------------------------------------------------------------

  def dict_to_pickle(self, file_name, index_dict):
    """Save a dictionary into a pickled file (into the output directory) for
       persistent storage. Add a pickle file extension to the given file name.
    """

    # Open the pickle file for writing
    #
    pickle_file_name = self.output_directory+file_name+self.pickle_file_ext

    try:
      pickle_file = open(pickle_file_name, 'w')
    except IOError:
      logging.exception('Cannot open pickle file "%s" for writing' % \
                        (pickle_file_name))
      raise IOError
    except:
      raise Exception

    logging.info('Write to pickle file "%s"' % (pickle_file_name))

    s = time.time()

    pickle.dump(index_dict, pickle_file, pickle.HIGHEST_PROTOCOL)
    pickle_file.close()

    logging.info('  Time used: %f' % (time.time() - s))

    logging.info('  Saved pickled dictionary with %d entries' % \
                 (len(index_dict)))

  # ---------------------------------------------------------------------------

  def pickle_to_dict(self, file_name):
    """Read a pickled dictionary from the output directory and return it. Add a
       pickle file extension to the given file name.
  """

    # Open the pickle file
    #
    pickle_file_name = self.output_directory+file_name+self.pickle_file_ext

    try:
      pickle_file = open(pickle_file_name, 'r')
    except IOError:
      logging.exception('Cannot read from pickle file "%s"' % \
                        (pickle_file_name))
      raise IOError
    except:
      raise Exception

    s = time.time()

    logging.info('Read from pickle file "%s"' % (pickle_file_name))

    pickle_dict = pickle.load(pickle_file)
    pickle_file.close()

    logging.info('  Time used: %f' % (time.time() - s))

    logging.info('  Loaded pickled dictionary with %d entries' % \
                 (len(pickle_dict)))

    return pickle_dict

# -----------------------------------------------------------------------------

  def load_austpost_lookup(self, austpost_lookup_file_name, state_filter=None):
    """Load and process an Australia Post look-up table with postcodes and
       suburb names and build dictionaries (for postcodes and suburbs) and a
       set with all postcode/suburb combinations, which will be used to
       validate (and impute if possible) postcode and locality (suburb)
       attributes in GNAF locality files.
       An optional 'state_filter' value can be given to filter out all records
       with a state value different from the one given.

       The Australia Post lookup table file is available from:

         http://www1.auspost.com.au/postcodes/index.asp?sub=2

       The resulting two dictionaries and one set will be inserted into the
       main ProcessGNAF object and made available as:

         self.aust_post_postcode_dict
         self.aust_post_suburb_dict
         self.aust_post_pc_sub_set
    """

    if (state_filter != None):
      if (state_filter not in \
          ['aat','act','nsw','nt','qld','sa','tas','vic','wa']):
        logging.exception('Illegal "state_filter" value given: %s' % \
                          (state_filter))
        raise Exception

    try:
      austpost_file = open(austpost_lookup_file_name, 'r')
    except IOError:
      logging.exception('Cannot open Australia Post look-up file "%s"' \
                        % (austpost_lookup_file_name))
      raise IOError
    except:
      raise Exception

    logging.info('Opened Australia Post look-up file "%s" for reading' % \
                 (austpost_lookup_file_name))

    postcode_dict = {}  # A dictionary with postcodes as keys
    suburb_dict =   {}  # A dictionary with suburb names as keys

    postcode_suburb_set = sets.Set()  # A set of postcode/suburb combinations

    # Define a CSV reader
    #
    csv_reader = csv.reader(austpost_file)

    # It is assumed the first line in the file is the header (don't process it)
    #
    csv_header_line = csv_reader.next()
    logging.debug('Read file header: "%s"' % (csv_header_line))

    row_count =  0
    start_time = time.time()

    # Now read the file content and process it
    #
    for csv_line in csv_reader:

      # Extract the three necessary fields
      #
      postcode = csv_line[0]
      suburb =   csv_line[1].lower()
      state =    csv_line[2].lower()
      logging.debug('  Extracted: Postcode: %s, suburb: %s, state: %s' % \
                    (postcode, suburb, state))

      postcode = postcode.strip()
      suburb =   suburb.strip()
      state =    state.strip()

      # Process if either no state filter or approrpiate state in record
      #
      if (state_filter == None) or (state_filter == state):

        # Remove some special words from the suburb end
        #
        if suburb[-3:] in [' dc', ' bc',' mc', ' pc']:
          logging.info('    Extension "%s" removed from suburb "%s"' % \
                       (suburb[-3:],suburb[:-3]) + ' in row %d' % (row_count))
          suburb = suburb[:-3]

        # Remove all hyphens and remove all single quotes in suburb names
        #
        if ('-' in suburb):
          logging.info('    Replace "-" with " " in suburb "%s"' % (suburb) + \
                       ' in row %d' % (row_count))
          suburb = suburb.replace('-', ' ')
        if ("'" in suburb):
          logging.debug("1:    Remove ' from suburb " + '"%s"' % (suburb) + \
                        ' in row %d' % (row_count))
          suburb = suburb.replace("'", "")

        # Correct 'st into 'saint' and 'mt' into 'mount' in suburbs
        #
        suburb = ' '+suburb+' '
        if (' mt ' in suburb):
          logging.info('    Correct "mt" into "mount" in suburb "%s"' % \
                       (suburb) + ' in row %d' % (row_count))
          suburb = suburb.replace(' mt ', ' mount ')
        if (' st ' in suburb):
          logging.info('    Correct "st" into "saint" in suburb "%s"' % \
                       (suburb) + ' in row %d' % (row_count))
          suburb = suburb.replace(' st ', ' saint ')

        suburb = suburb.replace('   ', ' ')  # Replace triple spaces with one
        suburb = suburb.replace('  ', ' ')  # Replace double spaces with one
        suburb = suburb.replace('  ', ' ')  # Replace double spaces with one
        suburb = suburb.strip()  # Make sure no spaces at beginning or end

        # Make suburb values the same as within Febrl (after standardisation)
        #
        suburb = suburb.replace(' ', '_')  # Replace all spaces between words

        if (len(postcode) != 4):
          logging.info('    Australia Post look-up file contains a postcode ' \
                       + 'which does not have 4 digits: "%s"' % (postcode) + \
                       ' in row %d' % (row_count))
          raise Exception

        if (state not in ['aat','act','nsw','nt','qld','sa','tas','vic','wa']):
          logging.exception('Australia Post look-up file contains an ' + \
                            'unknown state value: "%s"' % (state))
          raise Exception

        # Insert into dictionaries if the values are not empty
        #
        if (suburb != '') and (postcode != ''):

          suburb_list = postcode_dict.get(postcode,[])
          if suburb in suburb_list:
            logging.info('    Suburb "%s" already in list of suburbs ' % \
                         (suburb) + 'for postcode "%s"' % (postcode) + \
                         ' in row %d' % (row_count))
          else:
            suburb_list.append(suburb)
            suburb_list.sort()
            postcode_dict.update({postcode:suburb_list})

          postcode_list = suburb_dict.get(suburb,[])
          if postcode in postcode_list:
            logging.info('    Postcode "%s" already in list of postcodes ' % \
                         (postcode) + 'for suburb "%s"' % (suburb) + \
                         ' in row %d' % (row_count))
          else:
            postcode_list.append(postcode)
            postcode_list.sort()
            suburb_dict.update({suburb:postcode_list})

          # Now create the postcode/suburb combination and insert into set
          #
          postcode_suburb_comb = postcode+'/'+suburb

          if (postcode_suburb_comb in postcode_suburb_set):
            logging.info('    Postcode/suburb combination "%s" is more than ' \
                         % (postcode_suburb_comb) + \
                         'once in Australia Post look-up table file in ' + \
                         'row %d' % (row_count))
          else:
            postcode_suburb_set.add(postcode_suburb_comb)

      # Increase line count and print status message
      #
      row_count += 1
      if ((row_count % 1000) == 0):  # Process meter
        used_time = (time.time() - start_time)
        logging.info('    Processed %d rows in %s (time per row: %s)' % \
                     (row_count, output.time_string(used_time), \
                     output.time_string(used_time/row_count)))

    austpost_file.close()
    logging.info('  Read %d data lines' % (row_count))

    logging.info('  Australia Post postcode dictionary contains %d unique ' % \
                 (len(postcode_dict)) + 'postcodes')
    logging.info('  Australia Post suburb dictionary contains %d unique ' % \
                 (len(suburb_dict)) + 'suburbs')
    logging.info('  Australia Post postcode/suburb combination set contains ' \
                 + '%d entries' % (len(postcode_suburb_set)))

    self.aust_post_postcode_dict = postcode_dict
    self.aust_post_suburb_dict =   suburb_dict
    self.aust_post_pc_sub_set =    postcode_suburb_set

  # ---------------------------------------------------------------------------

  def clean_value(self, in_val, in_pid):
    """Small routine which cleans the input value using Febrl's
       'clean_component' and 'tag_address_coponent' routines.

       Returns a cleaned string.
    """

    in_val = in_val.strip()
    in_val = in_val.lower()

    clean_val = clean_component(in_val, self.address_correction_list, in_pid)

    if (clean_val != ''):
      if ('|' in clean_val):
        clean_val = clean_val.replace(' | ','_')
        clean_val = clean_val.replace(' |','')
        clean_val = clean_val.replace('| ','')
        clean_val = clean_val.replace('_|_','_')
        clean_val = clean_val.replace('_|','')
        clean_val = clean_val.replace('|_','')
      clean_val_lists = tag_address_component(clean_val, \
                                             self.address_lookup_table, in_pid)
      clean_val = '_'.join(clean_val_lists[0])
      clean_val = clean_val.replace(' ', '_')

    return clean_val

  # ---------------------------------------------------------------------------

  def check_pid(self, in_pid, pid_name, row_count):
    """Check if the given PID is not empty and doesn't end with ".00"

       Return a 'cleaned' PID or give a warning and return an empty string.
    """

    if (in_pid == ''):
      logging.warn('%s is empty in row %d' % (pid_name, row_count))

    else:
      # Make sure there is no '.XX' at the end of the PID
      #
      if (in_pid[-3] == '.'):
        in_pid = in_pid[:-3]  # Remove '.XX' (X = any character)

      if (in_pid == ''):
        logging.warn('%s is empty in row %d' % (pid_name, row_count))

    return in_pid

  # ---------------------------------------------------------------------------

  def check_pid_uniqueness(self, gnaf_file_info_list):
    """Load the given GNAF files, and enter the PIDs (in the given primary key
       attribute) into one dctionary - report dublicate PIDs.
    """

    pid_dict = {}

    logging.info('Check uniqueness of GNAF PIDs')

    for (gnaf_file_info, pid_attr) in gnaf_file_info_list:

      logging.info('  Process "%s" file' % (gnaf_file_info[0]))

      attr_dict = self.attr_list_to_dict(gnaf_file_info)
      if (pid_attr not in attr_dict):
        logging.exception('PID attribute "%s" not an attribute in the ' % \
                          (pid_attr) + 'file "%s"' % (gnaf_file_info[0]))
        raise Exception

      logging.info('    Check PID attribute "%s"' % (pid_attr))

      pid_attr_col = attr_dict[pid_attr]  # Get the column number

      # Open the CSV file
      #
      [csv_reader, csv_file] = self.open_csv(gnaf_file_info)

      row_count =  0
      start_time = time.time()

      # Now read the file content and process it
      #
      for csv_line in csv_reader:

        pid_val = csv_line[pid_attr_col]

        if (pid_val in pid_dict):  # This PID already occured somewhere else
          pid_info = pid_dict[pid_val].split(',')
          logging.warn('PID "%s" already appeared in attribute "%s" file ' \
                       % (pid_val,pid_info[0]) + '"%s" in row %d' % \
                       (pid_info[1], row_count))

        else:  # A new PID
          pid_info = pid_attr+','+gnaf_file_info[0]

          pid_dict[pid_val] = pid_info  # Insert into dictionary

        row_count += 1
        if ((row_count % 10000) == 0):  # Process meter
          used_time = (time.time() - start_time)
          logging.info('    Processed %d rows in %s (time per row: %s)' % \
                       (row_count, output.time_string(used_time), \
                       output.time_string(used_time/row_count)))

      csv_file.close()

    logging.info('Number of entries in the PID dictionary: %d' % \
                 (len(pid_dict)))

    # Following lines are commented because it takes too long to get examples
    #
    # logging.info('  Some example entires:'))
    # a = pid_dict.items()
    # logging.info('    %s' % (str(a[:10])))

  # ---------------------------------------------------------------------------

  def insert_locality_value(self, inv_index_dict, in_val, in_loc_pid):
    """Small routine to insert the given locality value and PID into the given
       inverted index.

       The value is only inserted if it is not an empty string.
    """

    if (in_val != ''):
      val_set = inv_index_dict.get(in_val, sets.Set())
      val_set.add(in_loc_pid)
      inv_index_dict[in_val] = val_set

  # ---------------------------------------------------------------------------

  def insert_street_value(self, inv_index_dict, in_val, in_street_pid,
                          in_loc_pids, do_locality):
    """Small routine to insert the given street value and PID into the given
       inverted index.

       'in_loc_pid' can either be a single PID or a set of PIDs.

       The value is only inserted if it is not an empty string.
    """

    # One LOCALITY_PID only, make it a list
    #
    if (isinstance(in_loc_pids, str) == True):
      in_loc_pids = [in_loc_pids]

    if (in_val != ''):
      val_dict = inv_index_dict.get(in_val, {})

      for in_loc_pid in in_loc_pids:
        val_set = val_dict.get(in_loc_pid, sets.Set())
        val_set.add(in_street_pid)
        val_dict[in_loc_pid] = val_set

      inv_index_dict[in_val] = val_dict

  # ---------------------------------------------------------------------------

  def insert_address_value(self, inv_index_dict, in_val, in_street_pid,
                           in_loc_pid, in_addr_pid_set):
    """Small routine to insert the given value and PIDs into the given inverted
       index.

       The value is only inserted if it is not an empty string.
    """

    if (in_val != ''):
      loc_pid_dict = inv_index_dict.get(in_val, {})
      street_pid_dict = loc_pid_dict.get(in_loc_pid, {})

      addr_pid_set = street_pid_dict.get(in_street_pid, sets.Set())

      addr_pid_set = addr_pid_set.union(in_addr_pid_set)

      street_pid_dict[in_street_pid] = addr_pid_set
      loc_pid_dict[in_loc_pid] = street_pid_dict

      inv_index_dict[in_val] = loc_pid_dict

  # ---------------------------------------------------------------------------

  def gnaf_reverse_dict(self, g_locality_file_info, g_locality_alias_file_info,
                        g_street_file_info, g_street_locality_alias_file_info,
                        g_address_detail_file_info):
    """Routine to create a dictionary with PIDs as keys and the original values
       in the files as values - used for checking the original G-NAF values
       when a PID is given (reverse look-up).

       The input files are:
         G_LOCALITY
         G_LOCALITY_ALIAS
         G_STREET
         G_STREET_LOCALITY_ALIAS
         G_ADDRESS_DETAIL

       The routine creates one large shelve containing all the useful entires
       from the above listed G-NAF files.

       Note that specials entries in the shelve contain the attribute names of
       the processed G-NAF files. These keys are the above file names.

       Some example entries look as follows:

       502298712:{'G_STREET':'502298712,"ALLAMANDA","PL",,20040221,"}
       500207172:{'G_LOCALITY':'500211195,"SMITHS CREEK",4481,,"NSW","","S"',
                  'G_LOCALITY_ALIAS':'500211195,"SYN","WOMBAT CREEK","NSW",'}

       The routine returns one dictionary.
    """

    file_info_list = [g_locality_file_info, g_locality_alias_file_info,
                      g_street_file_info, g_street_locality_alias_file_info,
                      g_address_detail_file_info]

    # A list of lists with the PID attribute names (which will be used as keys)
    #
    file_pid_attr = [['LOCALITY_PID'], ['LOCALITY_PID'],
                     ['STREET_PID'], ['STREET_PID', 'LOCALITY_PID'],
                     ['LOCALITY_PID', 'STREET_PID', 'ADDRESS_SITE_PID']]

    reverse_dict = {}

    for i in range(5):  # Loop over all files - - - - - - - - - - - - - - - - -

      # First insert attribute name list
      #
      reverse_dict[file_info_list[i][0]] = file_info_list[i][1]

      logging.info('Process G-NAF file: %s' % (file_info_list[i][0]))

      # Get a dictionary of the column names
      #
      attr_dict = self.attr_list_to_dict(file_info_list[i])

      logging.info('  PID names: %s' % (str(file_pid_attr[i])))

      # Get the column numbers
      #
      pid_col = []  # make a list of PID column numbers
      for pid_name in file_pid_attr[i]:
        pid_col.append(attr_dict[pid_name])

      logging.info('  PID columns: %s' % (str(pid_col)))

      # Open the CSV file
      #
      [csv_reader, csv_file] = self.open_csv(file_info_list[i])

      row_count =   0
      start_time =  time.time()

      # Now read the file content and process it
      #
      for csv_line in csv_reader:

        #logging.info('%s' % (','.join(csv_line))) ##################

        for col in pid_col:  # Loop over all PIDs in this file

          # Only insert into dictionary if PID is not empty
          #
          pid = csv_line[col]

          # Make sure there is no '.00' at the end of the PID
          #
          if (pid[-3:] == '.00'):
            pid = pid[:-3].strip()  # Remove '.XX'

          #logging.info('  %s' % (pid))  #############################

          if (pid != ''):
            key = ','.join(csv_line)
            val = file_info_list[i][0]  # The file name

            pid_dict = reverse_dict.get(pid,{})

            if (key in pid_dict):
              logging.warn('Reverse dictionary already contains an entry' + \
                           ' {%s:%s} for PID %s in row %d: \n %s' % \
                           (key,val,pid, row_count, str(pid_dict)))
            else:
              pid_dict[key] = val
            reverse_dict[pid] = pid_dict

        row_count += 1
        if ((row_count % 10000) == 0):  # Process meter
          used_time = (time.time() - start_time)
          logging.info('    Processed %d rows in %s (time per row: %s)' % \
                       (row_count, output.time_string(used_time), \
                       output.time_string(used_time/row_count)))

      csv_file.close()
      logging.info('  Read %d data lines' % (row_count))

    return reverse_dict

# -----------------------------------------------------------------------------

  def save_gnaf_address_csv(self, g_locality_file_info, g_street_file_info,
                            g_address_detail_file_info,
                            g_address_site_geocode_file_info, csv_file_name):
    """Routine to compile a file of all addresses in G-NAF including latitude
       and longitude and save them into a comma separated values (CSV) text
       file.
    """

    # Define a list of the output fields (attributes) in the CSV file - - - - -
    # (these must be attribute names available in the G-NAF files used)
    #
    csv_attribute_list = ["GNAF_PID",
                          "LOCALITY_PID",
                          "STREET_PID",
                          "ADDRESS_SITE_PID",
                          "FLAT_NUMBER",
                          "FLAT_TYPE",
                          "LEVEL_TYPE",
                          "LEVEL_NUMBER",
                          "BUILDING_NAME",
                          "NUMBER_FIRST",
                          "NUMBER_LAST",
                          "LOT_NUMBER",
                          "RURAL_ADDRESS",
                          "STREET_NAME",
                          "STREET_TYPE",
                          "STREET_SUFFIX",
                          "LOCALITY_NAME",
                          "STATE_ABBREVIATION",
                          "POSTCODE",
                          "LATITUDE",
                          "LONGITUDE"]

    # Define a mapping from original G-NAF file attributes into CSV attributes
    # (set to -1 if attribute is not used)
    #
    address_detail_mapping =  [0,1,2,-1,-1,4,-1,5,6,-1,7,-1,8,-1,-1,-1,9,-1, \
                               -1,10,-1,-1,11,-1,-1,12,-1,3]
    address_site_pid_index =  27 # Index of ADDRESS_SITE_PID in address detail
    street_pid_index =        2  # Index of STREET_PID in address detail file
    locality_pid_index =      1  # Index of LOCALITY PID in address detail

    street_mapping =          [-1,13,14,15,-1,-1]
    locality_mapping =        [-1,16,-1,-1,17,18,-1]
    address_geocode_mapping = [19,20]

    csv_empty_row_list = []
    for i in range(len(csv_attribute_list)):
      csv_empty_row_list.append('')

    # Process the G_LOCALITY file an make a dictionary with LOCALITY_PID keys -
    #
    locality_pid_dict = {}

    # Get a dictionary of the column names
    #
    attr_dict = self.attr_list_to_dict(g_locality_file_info)
    locality_pid_col = attr_dict['LOCALITY_PID']

    # Open the CSV file
    #
    [csv_reader, csv_file] = self.open_csv(g_locality_file_info)

    logging.info('  Process "G_LOCALITY" file')

    row_count =  0
    start_time = time.time()
    num_missing = 0  # Counter for number records with missing coordinates

    # Now read the file content and process it
    #
    for csv_line in csv_reader:

      locality_pid = csv_line[locality_pid_col].strip()

      if (locality_pid != ''):
        if (locality_pid[-3:] == '.00'):
          locality_pid = locality_pid[:-3]

      if (locality_pid != ''):
        if (locality_pid in locality_pid_dict):
          logging.exception('LOCALITY_PID %s already in dictionary' % \
                            (locality_pid))
          raise Exception
        else:
          locality_pid_dict[locality_pid] = csv_line

      else:
        num_missing += 1

      row_count += 1
      if ((row_count % 10000) == 0):  # Process meter
        used_time = (time.time() - start_time)
        logging.info('    Processed %d rows in %s (time per row: %s)' % \
                     (row_count, output.time_string(used_time), \
                     output.time_string(used_time/row_count)))

    csv_file.close()
    logging.info('    Read %d data lines' % (row_count))

    logging.info('    Locality PID index contains %d entries' % \
                 (len(locality_pid_dict)))
    logging.info('    Number of records with missing coordinates: %d' % \
                 (num_missing))

    # Process the G_STREET file an make a dictionary with STREET_PID keys - - -
    #
    street_pid_dict = {}

    # Get a dictionary of the column names
    #
    attr_dict = self.attr_list_to_dict(g_street_file_info)
    street_pid_col = attr_dict['STREET_PID']

    # Open the CSV file
    #
    [csv_reader, csv_file] = self.open_csv(g_street_file_info)

    logging.info('  Process "G_STREET" file')

    row_count =  0
    start_time = time.time()
    num_missing = 0  # Counter for number records with missing coordinates

    # Now read the file content and process it
    #
    for csv_line in csv_reader:

      street_pid = csv_line[street_pid_col].strip()

      if (street_pid != ''):
        if (street_pid[-3:] == '.00'):
          street_pid = street_pid[:-3]

      if (street_pid != ''):
        if (street_pid in street_pid_dict):
          logging.exception('STREET_PID %s already in dictionary' % \
                            (street_pid))
          raise Exception
        else:
          street_pid_dict[street_pid] = csv_line

      else:
        num_missing += 1

      row_count += 1
      if ((row_count % 10000) == 0):  # Process meter
        used_time = (time.time() - start_time)
        logging.info('    Processed %d rows in %s (time per row: %s)' % \
                     (row_count, output.time_string(used_time), \
                     output.time_string(used_time/row_count)))

    csv_file.close()
    logging.info('    Read %d data lines' % (row_count))

    logging.info('    Street PID index contains %d entries' % \
                 (len(street_pid_dict)))
    logging.info('    Number of records with missing coordinates: %d' % \
                 (num_missing))

    # Process the G_ADDRESS_SITE_GEOCODE file to get coordinates  - - - - - - -
    #
    address_site_pid_dict = {}

    # Get a dictionary of the column names
    #
    attr_dict = self.attr_list_to_dict(g_address_site_geocode_file_info)

    # Get the column numbers
    #
    address_site_pid_col = attr_dict['ADDRESS_SITE_PID']
    latitude_col =         attr_dict['LATITUDE']
    longitude_col =        attr_dict['LONGITUDE']

    # Open the CSV file
    #
    [csv_reader, csv_file] = self.open_csv(g_address_site_geocode_file_info)

    logging.info('  Process "G_ADDRESS_SITE_GEOCODE" file')

    row_count =  0
    start_time = time.time()
    num_missing = 0  # Counter for number records with missing coordinates

    # Now read the file content and process it
    #
    for csv_line in csv_reader:

      address_site_pid = csv_line[address_site_pid_col].strip()
      latitude =         csv_line[latitude_col].strip()
      longitude =        csv_line[longitude_col].strip()

      if (address_site_pid == ''):
        logging.warn('ADDRESS_SITE_PID is empty in row %d' % (row_count) + \
                     ' (will not be inserted into address site geocode index)')

      else:
        # Make sure there is no '.00' at the end of the ADDRESS_SITE_PID
        #
        if (address_site_pid[-3:] == '.00'):
          address_site_pid = address_site_pid[:-3]  # Remove '.00'

        if (latitude == '') or (longitude == ''):  # Coordinates missing
          logging.warn('Missing coordinates in row %d' % (row_count))
          num_missing += 1
        else:
          coordinates = [latitude,longitude]

          # Now insert into dictionary
          #
          if (address_site_pid in address_site_pid_dict):
            logging.exception('Row %d: ADDRESS_SITE_PID "%s" already in ' + \
                              'address ' % (row_count, \
                              address_site_pid_dict[address_site_pid]) + \
                              ' site geocode index')
            raise Exception
          else:
            address_site_pid_dict[address_site_pid] = coordinates

      row_count += 1
      if ((row_count % 10000) == 0):  # Process meter
        used_time = (time.time() - start_time)
        logging.info('    Processed %d rows in %s (time per row: %s)' % \
                     (row_count, output.time_string(used_time), \
                     output.time_string(used_time/row_count)))

    csv_file.close()
    logging.info('    Read %d data lines' % (row_count))

    logging.info('    Address site geocode index contains %d entries' % \
                 (len(address_site_pid_dict)))
    logging.info('    Number of records with missing coordinates: %d' % \
                 (num_missing))

    # Process the G_ADDRESS_DETAIL file an create CSV file  - - - - - - - - - -

    # Open CSV file for writing - - - - - - - - - - - - - - - - - - - - - - - -
    #
    csv_file_name = self.output_directory+csv_file_name

    try:
      csv_file = open(csv_file_name, 'w')
    except IOError:
      logging.exception('Cannot write to CSV file "%s"' % (csv_file_name))
      raise IOError
    except:
      raise Exception

    logging.info('Opened CSV file "%s" for writing' % (csv_file_name))

    # Define a CSV writer
    #
    csv_writer = csv.writer(csv_file)

    # Write header line with attribute names
    #
    csv_writer.writerow(csv_attribute_list)

    csv_row_count =  0
    gnaf_row_count = 0
    start_time = time.time()
    num_missing = 0  # Counter for number records with missing coordinates

    # Get a dictionary of the column names
    #
    attr_dict = self.attr_list_to_dict(g_address_detail_file_info)

    address_site_pid_col = attr_dict['ADDRESS_SITE_PID']
    street_pid_col =       attr_dict['STREET_PID']
    locality_pid_col =     attr_dict['LOCALITY_PID']

    # Open the CSV file
    #
    [csv_reader, csv_file] = self.open_csv(g_address_detail_file_info)

    logging.info('  Process "G_ADDRESS_DETAIL" file')

    # Now read the file content and process it
    #
    for csv_line in csv_reader:

      address_site_pid = csv_line[address_site_pid_col].strip()
      street_pid =       csv_line[street_pid_col].strip()
      locality_pid =     csv_line[locality_pid_col].strip()

      if (address_site_pid != '') and (address_site_pid[-3:] == '.00'):
        address_site_pid = address_site_pid[:-3]
      if (street_pid != '') and (street_pid[-3:] == '.00'):
        street_pid = street_pid[:-3]
      if (locality_pid != '') and (locality_pid[-3:] == '.00'):
        locality_pid = locality_pid[:-3]

      # Get detailed information from dictionaries - - - - - - - - - - - - - -
      #
      locality_csv_list = locality_pid_dict.get(locality_pid, [])
      street_csv_list =   street_pid_dict.get(street_pid, [])
      coordinates =       address_site_pid_dict.get(address_site_pid, [])

      # Only process record if both locality and street PIDs are valid  - - - -
      #
      if ((locality_csv_list != []) and (street_csv_list != []) and \
          (coordinates != [])):

        new_csv_list = csv_empty_row_list[:]

        i = 0
        for attr in address_detail_mapping:  # Copy G-NAF values into CSV list
          if (attr > -1):
            new_csv_list[attr] = csv_line[i]

            if (attr < 4):  # Fix '.00' in certain PIDs
              tmp_pid = new_csv_list[attr]
              if (tmp_pid[-3:] == '.00'):
                tmp_pid = tmp_pid[:-3]
                new_csv_list[attr] = tmp_pid
          i += 1

        i = 0
        for attr in street_mapping:  # Copy street detail values into CSV list
          if (attr > -1):
            new_csv_list[attr] = street_csv_list[i]
          i += 1

        i = 0
        for attr in locality_mapping:  # Copy locality detail values into CSV
          if (attr > -1):
            new_csv_list[attr] = locality_csv_list[i]
          i += 1

        # Finally copy coordinates into CSV list
        #
        new_csv_list[address_geocode_mapping[0]] = coordinates[0]
        new_csv_list[address_geocode_mapping[1]] = coordinates[1]

        # logging.info('CSV list from row %d: %s' % \
        #              (gnaf_row_count, str(new_csv_list)))

        csv_writer.writerow(new_csv_list)

      else:
        num_missing += 1

      gnaf_row_count += 1
      if ((gnaf_row_count % 10000) == 0):  # Process meter
        used_time = (time.time() - start_time)
        logging.info('    Processed %d rows in %s (time per row: %s)' % \
                     (gnaf_row_count, output.time_string(used_time), \
                     output.time_string(used_time/gnaf_row_count)))

    csv_file.close()

    logging.info('  Saved G-NAF address details into CSV file with %d ' % \
                 (csv_row_count) + \
                 'entires  (from %d G-NAF address detail records)' % \
                 (gnaf_row_count))
    logging.info('    Number of records with missing coordinates: %d' % \
                 (num_missing))

# -----------------------------------------------------------------------------

# Clean, standardise, index and save the address detail file - - - - - - - - -
#
# 14/07/2004
#
#    Read 4145365 data lines
#    Address site index contains           3336778 entries
#    Flat number prefix index contains     57 entries
#    Flat number index contains            3050 entries
#    Flat number suffix index contains     42 entries
#    Flat type index contains              25 entries
#    Level number index contains           0 entries
#    Level type index contains             7 entries
#    Building name index contains          32729 entries
#    Location description index contains   4359 entries
#    Number first prefix index contains    90 entries
#    Number first index contains           7684 entries
#    Number first suffix index contains    53 entries
#    Number last prefix index contains     13 entries
#    Number last index contains            1641 entries
#    Number last suffix index contains     17 entries
#    Lot number prefix index contains      26 entries
#    Lot number index contains             9377 entries
#    Lot number suffix index contains      18 entries
#    Number of records with unknown address site geocodes: 774800
#
# With address aliases included (15/07/2004)
#    Read 4145365 data lines
#    Address site index contains           3336778 entries
#    Flat number prefix index contains     57 entries
#    Flat number index contains            3050 entries
#    Flat number suffix index contains     42 entries
#    Flat type index contains              25 entries
#    Level number index contains           0 entries
#    Level type index contains             7 entries
#    Building name index contains          32729 entries
#    Location description index contains   4359 entries
#    Number first prefix index contains    90 entries
#    Number first index contains           7684 entries
#    Number first suffix index contains    53 entries
#    Number last prefix index contains     13 entries
#    Number last index contains            1641 entries
#    Number last suffix index contains     17 entries
#    Lot number prefix index contains      26 entries
#    Lot number index contains             9377 entries
#    Lot number suffix index contains      18 entries
#    Number of records with unknown address site geocodes: 774800
#    Number of records with alias addresses:               287123
