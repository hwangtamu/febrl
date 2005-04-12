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
# The Original Software is: "project-geocode.py"
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

"""Module project-geocode.py - Configuration for a geocoding project.

   Briefly, what needs to be defined for a geocoding project is:
   - A Febrl object, a project, plus a project logger
   - One input data set (initialised in read access mode)
   - One output data set (initialised in write or append access mode)
   - Look-up tables for suburb and postcode neighbouring regions
   - A Geocoder object with all the information needed to perform the geocoding
     process, including the geocode reference file(s), the inverted index
     files, input field definitions, references to various look-up tables, etc.
   - Lookup tables to be used
   - Standardisers for names, addresses and dates

   and then the 'geocode' method can be called.

   For more information see chapter

   "Configuration and Running Febrl using a Module derived from 'project.py'"

   in the Febrl manual.

   This project module will geocode the example data set 'dataset1.csv'
   given in the directory 'dsgen' using the example geocode reference data
   given in the dictionary 'geocode'.

   The directory separator 'dirsep' is a shorthand to os.sep as defined in
   febrl.py.
"""

# =============================================================================
# Imports go here

import sys
import time

from febrl import *            # Main Febrl classes
from geocoding import *        # Geocoding routines
from dataset import *          # Data set routines
from standardisation import *  # Standardisation routines
from comparison import *       # Comparison functions
from lookup import *           # Look-up table routines
from indexing import *         # Indexing and blocking routines
from simplehmm import *        # Hidden Markov model (HMM) routines
from classification import *   # Classifiers for weight vectors
from qgramindex import *       # Q-Gram index for approximate matching

# =============================================================================
# Define a project logger

init_febrl_logger(log_file_name = 'febrl-example-geocode.log',
                     file_level = 'WARN',
                  console_level = 'INFO',
                      clear_log = True,
                parallel_output = 'host')

# =============================================================================
# Set up Febrl and create a new project

myfebrl = Febrl(description = 'Example geocoding Febrl instance',
                 febrl_path = '.')

myproject = myfebrl.new_project(name = 'example-geocode',
                         description = 'Geocode example data set 1',
                           file_name = 'example-geocode.fbr',
                          block_size = 100,
                      parallel_write = 'host')

# =============================================================================
# Define original input data set

indata = DataSetCSV(name = 'testaddressesIn',
             description = 'Test addresses input data set',
             access_mode = 'read',
            header_lines = 0,
               file_name = 'geocode'+dirsep+'testaddresses-small.txt',
                  fields = {'address':0,
                           },
          fields_default = '',
            strip_fields = True,
          missing_values = ['','missing'])

# =============================================================================
# Define the output data set

outdata = DataSetCSV(name = 'testaddressesGeocoded',
              description = 'Geocoded test addresses',
              access_mode = 'write',
             write_header = True,
                file_name = 'geocode'+dirsep+'testaddresses-geocoded.csv',
         write_quote_char = '',
           missing_values = ['', 'missing', 'n/a'],
                   fields = {'wayfare_number':0,
                             'wayfare_name':1,
                             'wayfare_qualifier':2,
                             'wayfare_type':3,
                             'unit_number':4,
                             'unit_type':5,
                             'property_name':6,
                             'institution_name':7,
                             'institution_type':8,
                             'postaddress_number':9,
                             'postaddress_type':10,
                             'locality_name':11,
                             'locality_qualifier':12,
                             'postcode':13,
                             'territory':14,
# The next five output fields contain the geocoding latitude and longitude,
# match status, match weight and the G-NAF idendentifier(s)
                             'latitude':15,
                             'longitude':16,
                             'match_status':17,
                             'match_weight':18,
                             'gnaf_pid':19,
# The NSW collection district will be written to the output record
                             'collection_district':20,
                             'neighbour_level':21,
                             'max_avrg_distance':22,
# Finally we pass the original input address into the output
                             'input_record':23})

# =============================================================================
# Define and load the neighbouring regions look-up tables
#
pc_level_1_table = NeighbourLookupTable(name = 'PC-1', default = [])
pc_level_1_table.load(file_names = 'geocode'+dirsep+'pc-neighbours-1.txt')

pc_level_2_table = NeighbourLookupTable(name = 'PC-2', default = [])
pc_level_2_table.load(file_names = 'geocode'+dirsep+'pc-neighbours-2.txt')

sub_level_1_table = NeighbourLookupTable(name = 'Sub-1', default = [])
sub_level_1_table.load(file_names = 'geocode'+dirsep+'suburb-neighbours-1.txt')

sub_level_2_table = NeighbourLookupTable(name = 'Sub-2', default = [])
sub_level_2_table.load(file_names = 'geocode'+dirsep+'suburb-neighbours-2.txt')

# =============================================================================
# Define the approximate q-gram index for the geocoder (the field names must
# correspond to a key in the 'index_files' dictionary of the geocoder)
#
loc_name_qgram_index = PosQGramIndex(name = 'loc name q-gram index',
                              description = 'Q-Gram index for locality name',
                               field_name = 'locality_name',
                          q_gram_len_list = [(1,(1,4)),(2,(5,99))],
                            max_edit_dist = 2,
                         load_pickle_file = True,
                pickle_file_name = 'geocode'+dirsep+'loc_name_qgram_index.pik')

street_name_qgram_index = PosQGramIndex(name = 'street name q-gram index',
                              description = 'Q-Gram index for street name',
                               field_name = 'street_name',
                          q_gram_len_list = [(1,(1,4)),(2,(5,99))],
                            max_edit_dist = 2,
                         load_pickle_file = True,
             pickle_file_name = 'geocode'+dirsep+'street_name_qgram_index.pik')

# Put all approximate q-gram indices into a dictionary with field names as keys
# (the keys must correspond to the 'field_name' entries in the corresponding
# approximate index)
#
approx_indices = {'locality_name':loc_name_qgram_index,
                    'street_name':street_name_qgram_index}

# =============================================================================
# The main Geocoder object
#
example_geocoder = Geocoder(name = 'example1geocode',
                     description = 'Example geocoder',
          geocode_file_directory = 'gnaf'+dirsep+'shelve_pickles'+dirsep,
          pickle_files_extension = '.pik',
          shelve_files_extension = '.slv',
       address_site_geocode_file = 'address_site_geocode.slv',
    street_locality_geocode_file = 'street_locality_geocode.pik',
           locality_geocode_file = 'locality_geocode.pik',
        collection_district_file = 'collection_district.slv',
                  input_data_set = outdata,
         input_fields = {'building_name':('property_name',      1.0),
                  'location_description':('institution_name',   1.0),
                         'locality_name':('locality_name',      6.0),
                    'locality_qualifier':('locality_qualifier', 1.0),
                              'postcode':('postcode',           5.0),
                                 'state':('territory',          1.0),
                          'wayfare_name':('wayfare_name',       5.0),
                     'wayfare_qualifier':('wayfare_qualifier',  1.0),
                          'wayfare_type':('wayfare_type',       3.0),
                        'wayfare_number':('wayfare_number',     3.0),
                           'flat_number':('unit_number',        2.0),
                        'flat_qualifier':None,
                             'flat_type':('unit_type',          2.0),
                          'level_number':None,
                            'level_type':None,
                            'lot_number':('postaddress_number', 2.0),
                  'lot_number_qualifier':('postaddress_type',   1.0)},
                 match_threshold = 0.8,
                 best_match_only = True,
            missing_value_weight = 0.0,
         maximal_neighbour_level = 2,
             max_average_address = 50,
           postcode_neighbours_1 = pc_level_1_table,
           postcode_neighbours_2 = pc_level_2_table,
             suburb_neighbours_1 = sub_level_1_table,
             suburb_neighbours_2 = sub_level_2_table,
       index_files = {'building_name':('building_name',       'shelve'),
                        'flat_number':('flat_number',         'shelve'),
                 'flat_number_prefix':('flat_number_prefix',  'shelve'),
                 'flat_number_suffix':('flat_number_suffix',  'shelve'),
                          'flat_type':('flat_type',           'shelve'),
                       'level_number':('level_number',        'shelve'),
                         'level_type':('level_type',          'shelve'),
                      'locality_name':('locality_name',       'shelve'),
                     'location_descr':('location_descr',      'shelve'),
                         'lot_number':('lot_number',          'shelve'),
                  'lot_number_prefix':('lot_number_prefix',   'shelve'),
                  'lot_number_suffix':('lot_number_suffix',   'shelve'),
                       'number_first':('number_first',        'shelve'),
                'number_first_prefix':('number_first_prefix', 'shelve'),
                'number_first_suffix':('number_first_suffix', 'shelve'),
                        'number_last':('number_last',         'shelve'),
                 'number_last_prefix':('number_last_prefix',  'shelve'),
                 'number_last_suffix':('number_last_suffix',  'shelve'),
                           'postcode':('postcode',            'shelve'),
                       'state_abbrev':('state_abbrev',        'shelve'),
                        'street_name':('street_name',         'shelve'),
                      'street_suffix':('street_suffix',       'shelve'),
                        'street_type':('street_type',         'shelve')},
                    approx_index = approx_indices)

# =============================================================================
# Define and load lookup tables

address_lookup_table = TagLookupTable(name = 'Address lookup table',
                                   default = '')
address_lookup_table.load(['data'+dirsep+'country.tbl',
                           'data'+dirsep+'address_misc.tbl',
                           'data'+dirsep+'address_qual.tbl',
                           'data'+dirsep+'institution_type.tbl',
                           'data'+dirsep+'locality_name_act.tbl',
                           'data'+dirsep+'locality_name_nsw.tbl',
                           'data'+dirsep+'post_address.tbl',
                           'data'+dirsep+'postcode_act.tbl',
                           'data'+dirsep+'postcode_nsw.tbl',
                           'data'+dirsep+'saints.tbl',
                           'data'+dirsep+'territory.tbl',
                           'data'+dirsep+'unit_type.tbl',
                           'data'+dirsep+'wayfare_type.tbl'])

address_correction_list = CorrectionList(name = 'Address correction list')
address_correction_list.load('data'+dirsep+'address_corr.lst')

# =============================================================================
# Define and load hidden Markov models (HMMs)

address_states = ['wfnu','wfna1','wfna2','wfql','wfty','unnu','unty','prna1',
                  'prna2','inna1','inna2','inty','panu','paty','hyph','sla',
                  'coma','opbr','clbr','loc1','loc2','locql','pc','ter1',
                  'ter2','cntr1','cntr2','rubb']
address_tags = ['PC','N4','NU','AN','TR','CR','LN','ST','IN','IT','LQ','WT',
                'WN','UT','HY','SL','CO','VB','PA','UN','RU']

myaddress_hmm = hmm('Address HMM', address_states, address_tags)
myaddress_hmm.load_hmm('hmm'+dirsep+'geocode-nsw-address.hmm')

# myaddress_hmm.load_hmm('hmm'+dirsep+'address-absdiscount.hmm')
# myaddress_hmm.load_hmm('hmm'+dirsep+'address.hmm')
# myaddress_hmm.load_hmm('hmm'+dirsep+'address-laplace.hmm')

# =============================================================================
# Define a standardiser for addresses based on HMM

address_hmm_std = AddressHMMStandardiser(name = 'Address-HMM',
                                 input_fields = ['address'],
                                output_fields = ['wayfare_number',
                                                 'wayfare_name',
                                                 'wayfare_qualifier',
                                                 'wayfare_type',
                                                 'unit_number',
                                                 'unit_type',
                                                 'property_name',
                                                 'institution_name',
                                                 'institution_type',
                                                 'postaddress_number',
                                                 'postaddress_type',
                                                 'locality_name',
                                                 'locality_qualifier',
                                                 'postcode',
                                                 'territory',
                                                 None,
                                                 None],
                            address_corr_list = address_correction_list,
                            address_tag_table = address_lookup_table,
                                  address_hmm = myaddress_hmm)

# =============================================================================

pass_fields = PassFieldStandardiser(name = 'Pass fields',
                            input_fields = ['address'],
                           output_fields = ['input_record'])

# =============================================================================
# Define record standardiser(s) (one for each data set)

comp_stand = [address_hmm_std, pass_fields]

example_standardiser = RecordStandardiser(name = 'Example-std',
                                   description = 'Example standardiser',
                                 input_dataset = indata,
                                output_dataset = outdata,
                                      comp_std = comp_stand)

# =============================================================================
# Start the geocoding task
# - If 'first_record' is set to 'None' then it will be set to 0
# - If 'number_records' is set to 'None' then it will be set to the total
#   number of records in the input data set.

myproject.geocode(input_dataset = indata,
                 output_dataset = outdata,
               rec_standardiser = example_standardiser,
                   rec_geocoder = example_geocoder,
                   first_record = 0,
                 number_records = 1000)

# =============================================================================

myfebrl.finalise()

# =============================================================================
