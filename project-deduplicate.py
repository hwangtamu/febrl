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
# The Original Software is: "project-deduplicate.py"
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

"""Module project-deduplicate.py - Configuration for a deduplication project

   Briefly, what needs to be defined for a deduplication project is:
   - A Febrl object, a project, plus a project logger
   - One input data set
   - One corresponding temporary data set (with 'readwrite' access)
   - Lookup tables to be used
   - Standardisers for names, addresses and dates
   - Field comparator functions and a record comparator
   - A blocking index
   - A classifier

   and then the 'deduplicate' method can be called.

   For more information see chapter

   "Configuration and Running Febrl using a Module derived from 'project.py'"

   in the Febrl manual.

   This project module will standardised and then deduplicate the example data
   set 'dataset2.csv' given in the 'dsgen' directory.

   The directory separator 'dirsep' is a shorthand to os.sep as defined in
   febrl.py.
"""

# =============================================================================
# Imports go here

import sys
import time

from febrl import *            # Main Febrl classes
from dataset import *          # Data set routines
from standardisation import *  # Standardisation routines
from comparison import *       # Comparison functions
from lookup import *           # Look-up table routines
from indexing import *         # Indexing and blocking routines
from simplehmm import *        # Hidden Markov model (HMM) routines
from classification import *   # Classifiers for weight vectors

# =============================================================================
# Define a project logger

init_febrl_logger(log_file_name = 'febrl-example-dedup.log',
                     file_level = 'WARN',
                  console_level = 'INFO',
                      clear_log = True,
                parallel_output = 'host')

# =============================================================================
# Set up Febrl and create a new project (or load a saved project)

myfebrl = Febrl(description = 'Example deduplication Febrl instance',
                 febrl_path = '.')

myproject = myfebrl.new_project(name = 'example-dedup',
                         description = 'Deduplicate example data set 2',
                           file_name = 'example-deduplicate.fbr',
                          block_size = 100,
                      parallel_write = 'host')

# =============================================================================
# Define original input data set(s)
# Only one data set is needed for deduplication

indata = DataSetCSV(name = 'example2in',
             description = 'Example data set 2',
             access_mode = 'read',
            header_lines = 1,
               file_name = 'dsgen'+dirsep+'dataset2.csv',
                  fields = {'rec_id':0,
                            'given_name':1,
                            'surname':2,
                            'street_num':3,
                            'address_part_1':4,
                            'address_part_2':5,
                            'suburb':6,
                            'postcode':7,
                            'state':8,
                            'date_of_birth':9,
                            'soc_sec_id':10},
          fields_default = '',
            strip_fields = True,
          missing_values = ['','missing'])

# =============================================================================
# Define temporary data set(s) (one per input data set)
# Commented lines are only needed for the disk based Shelve data set

tmpdata = DataSetMemory(name = 'example2tmp',
#tmpdata = DataSetShelve(name = 'example2tmp',
#                   file_name = 'example2-shelve',
#                       clear = True,
                 description = 'Temporary example 2 data set',
                 access_mode = 'readwrite',
                      fields = {'title':1,
                                'gender_guess':2,
                                'given_name':3,
                                'alt_given_name':4,
                                'surname':5,
                                'alt_surname':6,
                                'wayfare_number':7,
                                'wayfare_name':8,
                                'wayfare_qualifier':9,
                                'wayfare_type':10,
                                'unit_number':11,
                                'unit_type':12,
                                'property_name':13,
                                'institution_name':14,
                                'institution_type':15,
                                'postaddress_number':16,
                                'postaddress_type':17,
                                'locality_name':18,
                                'locality_qualifier':19,
                                'postcode':20,
                                'territory':21,
                                'country':22,
                                'dob_day':23,
                                'dob_month':24,
                                'dob_year':25,
# The following are output fields that are passed without standardisation
                                'rec_id':0,
                                'soc_sec_id':26,
# The last output field contains the probability of the address HMM
                                'address_hmm_prob':27,
                               },
              missing_values = ['','missing'])

# =============================================================================
# Define and load lookup tables

name_lookup_table = TagLookupTable(name = 'Name lookup table',
                                default = '')
name_lookup_table.load(['data'+dirsep+'givenname_f.tbl',
                        'data'+dirsep+'givenname_m.tbl',
                        'data'+dirsep+'name_prefix.tbl',
                        'data'+dirsep+'name_misc.tbl',
                        'data'+dirsep+'saints.tbl',
                        'data'+dirsep+'surname.tbl',
                        'data'+dirsep+'title.tbl'])

name_correction_list = CorrectionList(name = 'Name correction list')
name_correction_list.load('data'+dirsep+'name_corr.lst')

surname_freq_table = FrequencyLookupTable(name = 'Surname frequency table',
                                       default = 1)
surname_freq_table.load('data'+dirsep+'surname_nsw_freq.csv')

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

pc_geocode_table = GeocodeLookupTable(name = 'Postcode locations',
                                   default = [])
pc_geocode_table.load('data'+dirsep+'postcode_centroids.csv')

# =============================================================================
# Define and load hidden Markov models (HMMs)

name_states = ['titl','baby','knwn','andor','gname1','gname2','ghyph',
               'gopbr','gclbr','agname1','agname2','coma','sname1','sname2',
               'shyph','sopbr','sclbr','asname1','asname2','pref1','pref2',
               'rubb']
name_tags = ['NU','AN','TI','PR','GF','GM','SN','ST','SP','HY','CO','NE','II',
             'BO','VB','UN','RU']

myname_hmm = hmm('Name HMM', name_states, name_tags)
myname_hmm.load_hmm('hmm'+dirsep+'name-absdiscount.hmm')
# myname_hmm.load_hmm('hmm'+dirsep+'name.hmm')
# myname_hmm.load_hmm('hmm'+dirsep+'name-laplace.hmm')

address_states = ['wfnu','wfna1','wfna2','wfql','wfty','unnu','unty','prna1',
                  'prna2','inna1','inna2','inty','panu','paty','hyph','sla',
                  'coma','opbr','clbr','loc1','loc2','locql','pc','ter1',
                  'ter2','cntr1','cntr2','rubb']
address_tags = ['PC','N4','NU','AN','TR','CR','LN','ST','IN','IT','LQ','WT',
                'WN','UT','HY','SL','CO','VB','PA','UN','RU']

myaddress_hmm = hmm('Address HMM', address_states, address_tags)
myaddress_hmm.load_hmm('hmm'+dirsep+'address-absdiscount.hmm')
# myaddress_hmm.load_hmm('hmm'+dirsep+'address.hmm')
# myaddress_hmm.load_hmm('hmm'+dirsep+'address-laplace.hmm')

# =============================================================================
# Define a list of date parsing format strings

date_parse_formats = ['%d %m %Y',   # 24 04 2002  or  24 4 2002
                      '%d %B %Y',   # 24 Apr 2002 or  24 April 2002
                      '%m %d %Y',   # 04 24 2002  or  4 24 2002
                      '%B %d %Y',   # Apr 24 2002 or  April 24 2002
                      '%Y %m %d',   # 2002 04 24  or  2002 4 24
                      '%Y %B %d',   # 2002 Apr 24 or  2002 April 24
                      '%Y%m%d',     # 20020424                   ISO standard
                      '%d%m%Y',     # 24042002
                      '%m%d%Y',     # 04242002
                      '%d %m %y',   # 24 04 02    or  24 4 02
                      '%d %B %y',   # 24 Apr 02   or  24 April 02
                      '%y %m %d',   # 02 04 24    or  02 4 24
                      '%y %B %d',   # 02 Apr 24   or  02 April 24
                      '%m %d %y',   # 04 24 02    or  4 24 02
                      '%B %d %y',   # Apr 24 02   or  April 24 02
                      '%y%m%d',     # 020424
                      '%d%m%y',     # 240402
                      '%m%d%y',     # 042402
                     ]

# =============================================================================
# Define standardisers for dates

dob_std = DateStandardiser(name = 'DOB-std',
                    description = 'Date of birth standardiser',
                   input_fields = 'date_of_birth',
                  output_fields = ['dob_day','dob_month', 'dob_year'],
                  parse_formats = date_parse_formats)

# =============================================================================
# Define a standardiser for names based on rules

name_rules_std = NameRulesStandardiser(name = 'Name-Rules',
                               input_fields = ['given_name','surname'],
                              output_fields = ['title',
                                               'gender_guess',
                                               'given_name',
                                               'alt_given_name',
                                               'surname',
                                               'alt_surname'],
                             name_corr_list = name_correction_list,
                             name_tag_table = name_lookup_table,
                                male_titles = ['mr'],
                              female_titles = ['ms'],
                            field_separator = ' ',
                           check_word_spill = True)

# =============================================================================
# Define a standardiser for names based on HMM

name_hmm_std = NameHMMStandardiser(name = 'Name-HMM',
                           input_fields = ['given_name','surname'],
                          output_fields = ['title',
                                           'gender_guess',
                                           'given_name',
                                           'alt_given_name',
                                           'surname',
                                           'alt_surname'],
                         name_corr_list = name_correction_list,
                         name_tag_table = name_lookup_table,
                            male_titles = ['mr'],
                          female_titles = ['ms'],
                               name_hmm = myname_hmm,
                        field_separator = ' ',
                       check_word_spill = True)

# =============================================================================
# Define a standardiser for addresses based on HMM

address_hmm_std = AddressHMMStandardiser(name = 'Address-HMM',
                                 input_fields = ['street_num','address_part_1',
                                                 'address_part_2','suburb',
                                                 'postcode', 'state'],
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
                                                 'country',
                                                 'address_hmm_prob'],
                            address_corr_list = address_correction_list,
                            address_tag_table = address_lookup_table,
                                  address_hmm = myaddress_hmm)

# =============================================================================
# Define a pass field standardiser for all fields that should be passed from
# the input to the output data set without any cleaning or standardisdation.

pass_fields = PassFieldStandardiser(name = 'Pass fields',
                            input_fields = ['rec_id', 'soc_sec_id'],
                           output_fields = ['rec_id', 'soc_sec_id'])

# =============================================================================
# Define record standardiser(s) (one for each data set)

comp_stand = [dob_std, name_rules_std, address_hmm_std, pass_fields]

# The HMM based name standardisation is not used in the above standardiser,
# uncomment the lines below (and comment the ones above) to use HMM
# standardisation for names.
#
#comp_stand = [dob_std, name_hmm_std, address_hmm_std, pass_fields]

example_standardiser = RecordStandardiser(name = 'Example-std',
                                   description = 'Example standardiser',
                                 input_dataset = indata,
                                output_dataset = tmpdata,
                                      comp_std = comp_stand)

# =============================================================================
# Define blocking index(es) (one per temporary data set)

myblock_def = [[('surname','dmetaphone', 4),('dob_year','direct')],
               [('given_name','truncate', 3), ('postcode','direct')],
               [('locality_name','nysiis'),('dob_month','direct')],
              ]

# Define one or more indexes (to be used in the classifier further below)

example_index = BlockingIndex(name = 'Index-blocking',
                           dataset = tmpdata,
                         index_def = myblock_def)

example_sorting_index = SortingIndex(name = 'Index-sorting',
                                  dataset = tmpdata,
                                index_def = myblock_def,
                              window_size = 3)

example_bigram_index = BigramIndex(name = 'Index-bigram',
                                dataset = tmpdata,
                              index_def = myblock_def,
                              threshold = 0.75)

# =============================================================================
# Define comparison functions for deduplication

given_name_nysiis = FieldComparatorEncodeString(name = 'Given name NYSIIS',
                                            fields_a = 'given_name',
                                            fields_b = 'given_name',
                                              m_prob = 0.95,
                                              u_prob = 0.001,
                                      missing_weight = 0.0,
                                       encode_method = 'nysiis',
                                             reverse = False)

surname_dmetaphone = FieldComparatorEncodeString(name = 'Surname D-Metaphone',
                                             fields_a = 'surname',
                                             fields_b = 'surname',
                                               m_prob = 0.95,
                                               u_prob = 0.001,
                                       missing_weight = 0.0,
                                        encode_method = 'dmetaphone',
                                              reverse = False,
                                      frequency_table = surname_freq_table,
                                freq_table_max_weight = 9.9,
                                freq_table_min_weight = -4.3)

wayfare_name_winkler = FieldComparatorApproxString(name = 'Wayfare name ' + \
                                                          'Winkler',
                                               fields_a = 'wayfare_name',
                                               fields_b = 'wayfare_name',
                                                 m_prob = 0.95,
                                                 u_prob = 0.001,
                                         missing_weight = 0.0,
                                         compare_method = 'winkler',
                                       min_approx_value = 0.7)

locality_name_key = FieldComparatorKeyDiff(name = 'Locality name key diff',
                                       fields_a = 'locality_name',
                                       fields_b = 'locality_name',
                                         m_prob = 0.95,
                                         u_prob = 0.001,
                                 missing_weight = 0.0,
                                   max_key_diff = 2)

postcode_distance = FieldComparatorDistance(name = 'Postcode distance',
                                        fields_a = 'postcode',
                                        fields_b = 'postcode',
                                          m_prob = 0.95,
                                          u_prob = 0.001,
                                  missing_weight = 0.0,
                                   geocode_table = pc_geocode_table,
                                    max_distance = 50.0)

age = FieldComparatorAge(name = 'Age',
                     fields_a = ['dob_day','dob_month', 'dob_year'],
                     fields_b = ['dob_day','dob_month', 'dob_year'],
            m_probability_day = 0.95,
            u_probability_day = 0.03333,
          m_probability_month = 0.95,
          u_probability_month = 0.083,
           m_probability_year = 0.95,
           u_probability_year = 0.01,
                max_perc_diff = 10.0,
                     fix_date = 'today')

field_comparisons = [given_name_nysiis, surname_dmetaphone, \
                     wayfare_name_winkler, locality_name_key, \
                     postcode_distance, age]

example_comparator = RecordComparator(tmpdata, tmpdata, field_comparisons)

# =============================================================================
# Define a classifier for classifying the matching vectors

example_fs_classifier = FellegiSunterClassifier(name = 'Fellegi and Sunter',
                                           dataset_a = tmpdata,
                                           dataset_b = tmpdata,
                                     lower_threshold = 0.0,
                                     upper_threshold = 30.0)

example_flex_classifier = FlexibleClassifier(name = 'Example flex classifier',
                                        dataset_a = tmpdata,
                                        dataset_b = tmpdata,
                                  lower_threshold = 0.0,
                                  upper_threshold = 10.0,
                                        calculate = [('avrg', [0,1]),
                                                     ('max',  [2,3,4]),
                                                     ('min',  [5])],
                                      final_funct = 'avrg')

# =============================================================================
# Start the deduplication task

myproject.deduplicate(input_dataset = indata,
                        tmp_dataset = tmpdata,
                   rec_standardiser = example_standardiser,
                     rec_comparator = example_comparator,
                     blocking_index = example_index,
                         classifier = example_flex_classifier,
                       first_record = 0,
                     number_records = 5000,
                   output_histogram = 'dedup-example-histogram.res',
            output_rec_pair_details = 'dedup-example-details.res',
            output_rec_pair_weights = 'dedup-example-weights.csv',
                   output_threshold = 10.0,
                  output_assignment = 'one2one')

# =============================================================================

myfebrl.finalise()

# =============================================================================
