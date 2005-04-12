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
# The Original Software is: "process-gnaf.py"
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

"""Module process-gnaf.py - Module to clean and index GNAF data files.

   This module controls how the original G-NAF files (assumed to be CSV files,
   i.e. comma separated values text files) are pre-processed into binary
   inverted index files to be used by the Febrl geocoding engine.

   Note that this pre-processing needs a machine with a large amount of main
   memory (almost 4 Gigabytes or more) in order to work (otherwise swapping
   will occur and slow down the pre-processing task tremendously).

   The main parts that need to be define are:
   - standard Febrl configuration options for  standardisation
   - input directory, i.e. where the original G-NAF files are
   - output directory, i.e. where to save the binary inverted indexes
     (note that this directory must exist - it will not be created)
   - description of all the G-NAF files (their file names and lists of their
     attributes)

   The processing of the G-NAF files into inverted indexes can be controlled by
   setting the appropriate flags below.

   The directory separator 'dirsep' is a shorthand to os.sep as defined in
   febrl.py.
"""

# =============================================================================
# Imports go here

import sys  # So we can add the febrl path to the Python search path
sys.path.append('..')

from febrl import *            # Main Febrl classes

from gnaffunctions import *    # Methods and attributes for processing GNAF
from lookup import *           # Look-up table routines

from dataset import *  # Otherwise we get an exception with logging system ;-(

# =============================================================================
# Flags that control the G-NAF pre-processing, set to either True or False

check_pid_uniqueness = False

save_pickle_files = True  # Save inverted indexes into binary Python pickles
save_shelve_files = True  # Save inverted indexes into binary Python shelves
save_text_files   = True # Save inverted indexes into text files

process_coll_dist_files = False  # Process collection district files

process_locality_files = True   # Process the G-NAF locality related files
process_street_files =  False    # Process the G-NAF street related files
process_address_files = False    # Process the G-NAF address related files

create_reverse_lookup_shelve = False  # Create one large shelve to be used for
                                     # reverse look-ups (i.e. given one or more
                                     # PID find the correspnding G-NAF records)

create_gnaf_address_csv_file = False # Create one large CSV file with all G-NAF
                                     # addresses (values merged from several
                                     # files) and their locations
gnaf_address_csv_file_name = 'gnaf_address_geocodes.csv'  # Corresponding file

# =============================================================================
# Define a project logger

init_febrl_logger(log_file_name = 'febrl-process-gnaf.log',
                     file_level = 'WARN',
                  console_level = 'INFO',
                      clear_log = True,
                parallel_output = 'host')

# =============================================================================
# Set up Febrl and create a new project (or load a saved project)

gnaf_febrl = Febrl(description = 'Process GNAF data files',
                   febrl_path = '.')

gnaf_project = gnaf_febrl.new_project(name = 'Process GNAF',
                         description = 'Clean and create indexes for GNAF',
                           file_name = 'processgnaf.fbr')

# =============================================================================
# Define and load address correction lists and tagging look-up tables

addr_lookup_table = TagLookupTable(name = 'Address tagging lookup table',
                                default = '')
addr_lookup_table.load(['..'+dirsep+'data'+dirsep+'country.tbl',
                        '..'+dirsep+'data'+dirsep+'address_misc.tbl',
                        '..'+dirsep+'data'+dirsep+'address_qual.tbl',
                        '..'+dirsep+'data'+dirsep+'institution_type.tbl',
                        '..'+dirsep+'data'+dirsep+'locality_name_act.tbl',
                        '..'+dirsep+'data'+dirsep+'locality_name_nsw.tbl',
                        '..'+dirsep+'data'+dirsep+'post_address.tbl',
                        '..'+dirsep+'data'+dirsep+'postcode_act.tbl',
                        '..'+dirsep+'data'+dirsep+'postcode_nsw.tbl',
                        '..'+dirsep+'data'+dirsep+'saints.tbl',
                        '..'+dirsep+'data'+dirsep+'territory.tbl',
                        '..'+dirsep+'data'+dirsep+'unit_type.tbl',
                        '..'+dirsep+'data'+dirsep+'wayfare_type.tbl'])

addr_correction_list = CorrectionList(name = 'Address correction list')
addr_correction_list.load('..'+dirsep+'data'+dirsep+'address_corr.lst')

# =============================================================================
# Define locations, file extensions for GNAF and index (pickle / shelve) files

gnaf_input_directory =   '..'+dirsep+'..'+dirsep+'..'+dirsep+'data' + \
                         dirsep+'gnaf'+dirsep
gnaf_output_directory =  '..'+dirsep+'..'+dirsep+'..'+dirsep+'data' + \
                         dirsep+'gnaf'+dirsep+'shelve_pickles'+dirsep
gnaf_input_file_ext =    '.csv'
pickle_output_file_ext = '.pik'
shelve_output_file_ext = '.slv'

# =============================================================================
# Define the GNAF data files (file names and attributes names)

gnaf_files = {}

# Address related data files - - - - - - - - - - - - - - - - - - - - - - - - -

gnaf_files['g_address_detail'] = ['G_ADDRESS_DETAIL',
                                  ["GNAF_PID",
                                   "LOCALITY_PID",
                                   "STREET_PID",
                                   "PROPERTY_PID",
                                   "FLAT_NUMBER_PREFIX",
                                   "FLAT_NUMBER",
                                   "FLAT_NUMBER_SUFFIX",
                                   "FLAT_TYPE",
                                   "LEVEL_TYPE",
                                   "LEVEL_NUMBER_PREFIX",
                                   "LEVEL_NUMBER",
                                   "LEVEL_NUMBER_SUFFIX",
                                   "BUILDING_NAME",
                                   "LOCATION_DESCRIPTION",
                                   "PRIVATE_ROAD",
                                   "NUMBER_FIRST_PREFIX",
                                   "NUMBER_FIRST",
                                   "NUMBER_FIRST_SUFFIX",
                                   "NUMBER_LAST_PREFIX",
                                   "NUMBER_LAST",
                                   "NUMBER_LAST_SUFFIX",
                                   "LOT_NUMBER_PREFIX",
                                   "LOT_NUMBER",
                                   "LOT_NUMBER_SUFFIX",
                                   "LEGAL_PARCEL_ID",
                                   "RURAL_ADDRESS",
                                   "CONFIDENCE",
                                   "ADDRESS_SITE_PID"]]

gnaf_files['g_address_type'] = ['G_ADDRESS_TYPE',
                                ["ADDRESS_TYPE",
                                 "DESCRIPTION"]]

gnaf_files['g_address_alias'] = ['G_ADDRESS_ALIAS',
                                 ["PRINCIPAL_PID",
                                  "ALIAS_PID",
                                  "ALIAS_TYPE",
                                  "ALIAS_COMMENT",
                                  "DATE_CREATED",
                                  "DATE_RETIRED"]]

gnaf_files['g_address_alias_type'] = ['G_ADDRESS_ALIAS_TYPE',
                                      ["ALIAS_TYPE",
                                       "DESCRIPTION"]]

gnaf_files['g_address_site'] = ['G_ADDRESS_SITE',
                                ["ADDRESS_SITE_PID",
                                 "DATE_CREATED",
                                 "DATE_RETIRED",
                                 "ADDRESS_TYPE",
                                 "ADDRESS_SITE_NAME"]]

gnaf_files['g_flat_type'] = ['G_FLAT_TYPE',
                             ["FLAT_TYPE",
                              "DESCRIPTION"]]

gnaf_files['g_level_type'] = ['G_LEVEL_TYPE',
                              ["LEVEL_TYPE",
                               "DESCRIPTION"]]

# Geocode data files - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

gnaf_files['geocode_reliability'] = ['GEOCODE_RELIABILITY',
                                     ["RELIABILITY",
                                      "DESCRIPTION"]]

gnaf_files['g_geocode_type'] = ['G_GEOCODE_TYPE',
                                ["GEOCODE_TYPE",
                                 "DESCRIPTION"]]

gnaf_files['g_address_site_geocode'] = ['G_ADDRESS_SITE_GEOCODE',
                                        ["ADDRESS_SITE_PID",
                                         "RELIABILITY",
                                         "BOUNDARY_EXTENT",
                                         "PLANIMETRIC_ACCURACY",
                                         "DATE_GEOCODE_CREATED",
                                         "DATE_GEOCODE_RETIRED",
                                         "LATITUDE",
                                         "LONGITUDE",
                                         "ELEVATION",
                                         "GEOCODE_SITE_NAME",
                                         "GEOCODE_SITE_DESCRIPTION",
                                         "GEOCODE_TYPE"]]

gnaf_files['g_locality_geocode'] = ['G_LOCALITY_GEOCODE',
                                    ["LOCALITY_PID",
                                     "BOUNDARY_EXTENT",
                                     "PLANIMETRIC_ACCURACY",
                                     "DATE_GEOCODE_CREATED",
                                     "DATE_GEOCODE_RETIRED",
                                     "LATITUDE",
                                     "LONGITUDE",
                                     "ELEVATION",
                                     "GEOCODE_SITE_NAME",
                                     "GEOCODE_SITE_DESCRIPTION",
                                     "GID"]]

gnaf_files['g_street_locality_geocode'] = ['G_STREET_LOCALITY_GEOCODE',
                                           ["STREET_PID",
                                            "LOCALITY_PID",
                                            "STREET_CONFIDENCE",
                                            "BOUNDARY_EXTENT",
                                            "PLANIMETRIC_ACCURACY",
                                            "DATE_GEOCODE_CREATED",
                                            "DATE_GEOCODE_RETIRED",
                                            "LATITUDE",
                                            "LONGITUDE",
                                            "ELEVATION",
                                            "GEOCODE_SITE_NAME",
                                            "GEOCODE_SITE_DESCRIPTION",
                                            "GID"]]

# Street related data files - - - - - - - - - - - - - - - - - - - - - - - - - -

gnaf_files['g_street'] = ['G_STREET',
                          ["STREET_PID",
                           "STREET_NAME",
                           "STREET_TYPE",
                           "STREET_SUFFIX",
                           "DATE_CREATED",
                           "DATE_ARCHIVED"]]

gnaf_files['g_street_locality_alias'] = ['G_STREET_LOCALITY_ALIAS',
                                         ["STREET_PID",
                                          "LOCALITY_PID",
                                          "ALIAS_TYPE",
                                          "STREET_NAME",
                                          "STREET_TYPE",
                                          "STREET_SUFFIX"]]

gnaf_files['g_street_locality_alias_type'] = ['G_STREET_LOCALITY_ALIAS_TYPE',
                                              ["ALIAS_TYPE",
                                               "DESCRIPTION"]]

gnaf_files['g_street_suffix'] = ['G_STREET_SUFFIX',
                                 ["STREET_SUFFIX",
                                  "DESCRIPTION"]]

gnaf_files['g_street_type'] = ['G_STREET_TYPE',
                               ["STREET_TYPE",
                                "DESCRIPTION"]]

# Locality related data files - - - - - - - - - - - - - - - - - - - - - - - - -

gnaf_files['g_locality'] = ['G_LOCALITY',
                            ["LOCALITY_PID",
                             "LOCALITY_NAME",
                             "SUBURB_PID",
                             "GAZETTER_PID",
                             "STATE_ABBREVIATION",
                             "POSTCODE",
                             "LOCALITY_TYPE"]]

gnaf_files['g_locality_type'] = ['G_LOCALITY_TYPE',
                                 ["LOCALITY_TYPE",
                                  "DESCRIPTION"]]

gnaf_files['g_locality_alias'] = ['G_LOCALITY_ALIAS',
                                  ["LOCALITY_PID",
                                   "ALIAS_TYPE",
                                   "LOCALITY_NAME",
                                   "STATE_ABBREVIATION",
                                   "POSTCODE"]]

gnaf_files['g_locality_alias_type'] = ['G_LOCALITY_ALIAS_TYPE',
                                       ["ALIAS_TYPE",
                                        "DESCRIPTION"]]

gnaf_files['g_state'] = ['G_STATE',
                         ["STATE_ABBREVIATION",
                          "STATE_NAME"]]

# NSW Census Collection District boundaries data files - - - - - - - - - - - -

gnaf_files['address_cd'] = ['ADDRESS_CD',
                            ["ADDRESS_SITE_PID",
                             "CD"]]

gnaf_files['street_locality_cd'] = ['STREET_LOCALITY_CD',
                                    ["STREET_PID",
                                     "LOCALITY_PID",
                                     "CD"]]

gnaf_files['locality_cd'] = ['LOCALITY_CD',
                            ["LOCALITY_PID",
                             "CD"]]

# =============================================================================
# Define an Australia Post postcode/suburb look-up file and a possible state
# filter (or None if no state filter should be used).
#
austpost_lookup_file_name = '..'+dirsep+'..'+dirsep+'..'+dirsep+'data' + \
                            dirsep+'gnaf'+dirsep+'austpost'+dirsep + \
                            'pc-full_20050222.csv'
austpost_lookup_state_filter = 'nsw'

# =============================================================================
# =============================================================================
# Do not change anything below here
# =============================================================================
# =============================================================================


# =============================================================================
# Define a GNAF process and set its attributes

gnafproc = ProcessGNAF(gnaf_directory = gnaf_input_directory,
                     output_directory = gnaf_output_directory,
                        gnaf_file_ext = gnaf_input_file_ext,
                      pickle_file_ext = pickle_output_file_ext,
                      shelve_file_ext = shelve_output_file_ext,
                    stats_max_num_val = 20,
              address_correction_list = addr_correction_list,
                 address_lookup_table = addr_lookup_table)

# =============================================================================
# Check PIDs for uniqueness

if (check_pid_uniqueness == True):
  gnafproc.check_pid_uniqueness([(gnaf_files['g_address_detail'], 'GNAF_PID'),
                           (gnaf_files['g_address_site'], 'ADDRESS_SITE_PID'),
                           (gnaf_files['g_street'],'STREET_PID'),
                           (gnaf_files['g_locality'],'LOCALITY_PID')])

# =============================================================================
# Load and process an Australia Post postcode/suburb look-up file

gnafproc.load_austpost_lookup(austpost_lookup_file_name,
                              austpost_lookup_state_filter)

# =============================================================================
# Now process the various GNAF files into indexes and save them as pickles and
# shelve files

# Process and index Census Collection District boundaries files - - - - - - - -
#
if (process_coll_dist_files == True):

  cd_ind = gnafproc.index_collect_dist(gnaf_files['address_cd'],
                                       gnaf_files['street_locality_cd'],
                                       gnaf_files['locality_cd'])

  if (save_pickle_files == True):
    gnafproc.dict_to_pickle('collection_district', cd_ind)
  if (save_shelve_files == True):
    gnafproc.dict_to_shelve('collection_district', cd_ind)
  if (save_text_files == True):
    gnafproc.dict_to_text('collection_district',   cd_ind)

  del cd_ind

# Clean, standardise, index and save the locality files - - - - - - - - - - - -
#
if (process_locality_files == True):

  loc_ind = gnafproc.index_locality(gnaf_files['g_locality'],
                                    gnaf_files['g_locality_alias'],
                                    gnaf_files['g_locality_geocode'])

  if (save_pickle_files == True):
    gnafproc.dict_to_pickle('locality_name',    loc_ind[0])
    gnafproc.dict_to_pickle('postcode',         loc_ind[1])
    gnafproc.dict_to_pickle('state_abbrev',     loc_ind[2])
    gnafproc.dict_to_pickle('locality_geocode', loc_ind[3])

  if (save_shelve_files == True):
    gnafproc.dict_to_shelve('locality_name',    loc_ind[0])
    gnafproc.dict_to_shelve('postcode',         loc_ind[1])
    gnafproc.dict_to_shelve('state_abbrev',     loc_ind[2])
    gnafproc.dict_to_shelve('locality_geocode', loc_ind[3])

  if (save_text_files == True):
    gnafproc.dict_to_text('locality_name',      loc_ind[0])
    gnafproc.dict_to_text('postcode',           loc_ind[1])
    gnafproc.dict_to_text('state_abbrev',       loc_ind[2])
    gnafproc.dict_to_text('locality_geocode',   loc_ind[3])

  del loc_ind  # Clean up memory

# Clean, standardise, index and save the street files - - - - - - - - - - - - -
#
if (process_street_files == True):

  street_ind = gnafproc.index_street(gnaf_files['g_street'],
                                     gnaf_files['g_street_locality_alias'],
                                     gnaf_files['g_street_locality_geocode'],
                                     gnaf_files['g_address_detail'])

  if (save_pickle_files == True):
    gnafproc.dict_to_pickle('street_name',             street_ind[0])
    gnafproc.dict_to_pickle('street_type',             street_ind[1])
    gnafproc.dict_to_pickle('street_suffix',           street_ind[2])
    gnafproc.dict_to_pickle('street_locality_geocode', street_ind[3])

  if (save_shelve_files == True):
    gnafproc.dict_to_shelve('street_name',             street_ind[0])
    gnafproc.dict_to_shelve('street_type',             street_ind[1])
    gnafproc.dict_to_shelve('street_suffix',           street_ind[2])
    gnafproc.dict_to_shelve('street_locality_geocode', street_ind[3])

  if (save_text_files == True):
    gnafproc.dict_to_text('street_name',               street_ind[0])
    gnafproc.dict_to_text('street_type',               street_ind[1])
    gnafproc.dict_to_text('street_suffix',             street_ind[2])
    gnafproc.dict_to_text('street_locality_geocode',   street_ind[3])

  del street_ind  # Clean up memory

# Clean, standardise, index and save the address detail file - - - - - - - - -
#
if (process_address_files == True):

  addr_ind = gnafproc.index_address(gnaf_files['g_address_detail'],
                                    gnaf_files['g_address_site_geocode'],
                                    gnaf_files['g_address_alias'])

  if (save_pickle_files == True):
    gnafproc.dict_to_pickle('flat_number_prefix',   addr_ind[0])
    gnafproc.dict_to_pickle('flat_number',          addr_ind[1])
    gnafproc.dict_to_pickle('flat_number_suffix',   addr_ind[2])
    gnafproc.dict_to_pickle('flat_type',            addr_ind[3])
    gnafproc.dict_to_pickle('level_number',         addr_ind[4])
    gnafproc.dict_to_pickle('level_type',           addr_ind[5])
    gnafproc.dict_to_pickle('building_name',        addr_ind[6])
    gnafproc.dict_to_pickle('location_descr',       addr_ind[7])
    gnafproc.dict_to_pickle('number_first_prefix',  addr_ind[8])
    gnafproc.dict_to_pickle('number_first',         addr_ind[9])
    gnafproc.dict_to_pickle('number_first_suffix',  addr_ind[10])
    gnafproc.dict_to_pickle('number_last_prefix',   addr_ind[11])
    gnafproc.dict_to_pickle('number_last',          addr_ind[12])
    gnafproc.dict_to_pickle('number_last_suffix',   addr_ind[13])
    gnafproc.dict_to_pickle('lot_number_prefix',    addr_ind[14])
    gnafproc.dict_to_pickle('lot_number',           addr_ind[15])
    gnafproc.dict_to_pickle('lot_number_suffix',    addr_ind[16])
    gnafproc.dict_to_pickle('address_site_geocode', addr_ind[17])

  if (save_shelve_files == True):
    gnafproc.dict_to_shelve('flat_number_prefix',   addr_ind[0])
    gnafproc.dict_to_shelve('flat_number',          addr_ind[1])
    gnafproc.dict_to_shelve('flat_number_suffix',   addr_ind[2])
    gnafproc.dict_to_shelve('flat_type',            addr_ind[3])
    gnafproc.dict_to_shelve('level_number',         addr_ind[4])
    gnafproc.dict_to_shelve('level_type',           addr_ind[5])
    gnafproc.dict_to_shelve('building_name',        addr_ind[6])
    gnafproc.dict_to_shelve('location_descr',       addr_ind[7])
    gnafproc.dict_to_shelve('number_first_prefix',  addr_ind[8])
    gnafproc.dict_to_shelve('number_first',         addr_ind[9])
    gnafproc.dict_to_shelve('number_first_suffix',  addr_ind[10])
    gnafproc.dict_to_shelve('number_last_prefix',   addr_ind[11])
    gnafproc.dict_to_shelve('number_last',          addr_ind[12])
    gnafproc.dict_to_shelve('number_last_suffix',   addr_ind[13])
    gnafproc.dict_to_shelve('lot_number_prefix',    addr_ind[14])
    gnafproc.dict_to_shelve('lot_number',           addr_ind[15])
    gnafproc.dict_to_shelve('lot_number_suffix',    addr_ind[16])
    gnafproc.dict_to_shelve('address_site_geocode', addr_ind[17])

  if (save_text_files == True):
    gnafproc.dict_to_text('flat_number_prefix',   addr_ind[0])
    gnafproc.dict_to_text('flat_number',          addr_ind[1])
    gnafproc.dict_to_text('flat_number_suffix',   addr_ind[2])
    gnafproc.dict_to_text('flat_type',            addr_ind[3])
    gnafproc.dict_to_text('level_number',         addr_ind[4])
    gnafproc.dict_to_text('level_type',           addr_ind[5])
    gnafproc.dict_to_text('building_name',        addr_ind[6])
    gnafproc.dict_to_text('location_descr',       addr_ind[7])
    gnafproc.dict_to_text('number_first_prefix',  addr_ind[8])
    gnafproc.dict_to_text('number_first',         addr_ind[9])
    gnafproc.dict_to_text('number_first_suffix',  addr_ind[10])
    gnafproc.dict_to_text('number_last_prefix',   addr_ind[11])
    gnafproc.dict_to_text('number_last',          addr_ind[12])
    gnafproc.dict_to_text('number_last_suffix',   addr_ind[13])
    gnafproc.dict_to_text('lot_number_prefix',    addr_ind[14])
    gnafproc.dict_to_text('lot_number',           addr_ind[15])
    gnafproc.dict_to_text('lot_number_suffix',    addr_ind[16])
    gnafproc.dict_to_text('address_site_geocode', addr_ind[17])

  del addr_ind  # Clean up memory

# Create the reverse look-up dictionary - - - - - - - - - - - - - - - - - - - -
#
if (create_reverse_lookup_shelve == True):

  rev_dict = gnafproc.gnaf_reverse_dict(gnaf_files['g_locality'],
                                        gnaf_files['g_locality_alias'],
                                        gnaf_files['g_street'],
                                        gnaf_files['g_street_locality_alias'],
                                        gnaf_files['g_address_detail'])

  # Save into a shelve file to be used with the auxilliary program
  # 'gnaf_reverse.py' (for quick interactive reverse look-ups)
  #
  gnafproc.dict_to_shelve('gnaf_reverse_index', rev_dict)

  del rev_dict  # Clean up memory

# Safe G-NAF address details and coordinates into a CSV text file - - - - - - -
#
if (create_gnaf_address_csv_file == True):

  gnafproc.save_gnaf_address_csv(gnaf_files['g_locality'],
                                 gnaf_files['g_street'],
                                 gnaf_files['g_address_detail'],
                                 gnaf_files['g_address_site_geocode'],
                                 gnaf_address_csv_file_name)

# =============================================================================
