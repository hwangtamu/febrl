# =============================================================================
# project-linkage.py - Configuration for a record linkage project.
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
# The Original Software is "project-linkage.py".
# The Initial Developers of the Original Software are Dr Peter Christen
# (Department of Computer Science, Australian National University), Dr Tim
# Churches (Centre for Epidemiology and Research, New South Wales Department
# of Health) and Drs Markus Hegland, Stephen Roberts and Ole Nielsen
# (Mathematical Sciences Insitute, Australian National University). Copyright
# (C) 2002 the Australian National University and others. All Rights Reserved.
# Contributors:
#
# =============================================================================

"""Module project-linkage.py - Configuration for a record linkage project

   Briefly, what needs to be defined for a linkage project is:
   - A Febrl object, a project, plus a project logger
   - Two input data set
   - Two corresponding temporary data set (with read and write access)
   - Lookup tables to be used
   - Standardisers for names, addresses and dates for both input data sets
     (these standardisers can differ for thetwo input data sets)
   - Field comparator functions and a record comparator
   - Two blocking index (they should be the same, i.e. lead to similar blocking
     variables)
   - A classifier

   and then the 'link' method can be called.

   For more information see chapter

   "Configuration and Running Febrl using a Module derived from 'project.py'"

   in the Febrl manual.
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
# Set up Febrl and create a new project (or load a saved project)

myfebrl = Febrl(description = 'Linkage Febrl instance',
                 febrl_path = '.')

myproject = myfebrl.new_project(name = 'MDC',
                         description = 'Link Midwifes Data Collection',
                           file_name = 'mdc-linkage.fbr',
                          block_size = 1000)

# =============================================================================
# Define a project logger

mylog = ProjectLog(file_name = 'mdc-linkage.log',
                     project = myproject,
                   log_level = 1,
               verbose_level = 1,
                   clear_log = True,
                     no_warn = False,
              parallel_print = 'host')

print myfebrl
print myproject

# =============================================================================
# Define original input data set(s)
#
# Two data sets are needed for linkage
#
indata_a = DataSetCSV(name = 'MDC-a',
               description = 'Midwifes Data Collection, 1990-2000, ' + \
                             'source: NSW Department of Health',
               access_mode = 'read',
              header_lines = 1,
              write_header = True,
                 file_name = '../../data/nswhealth_mdc/mdc-100000.csv',
                    fields = {'year':0,
                              'rseqnum':1,
                              'ohoscode':2,
                              'omrn':3,
                              'gname':4,
                              'sname':5,
                              'omdob':6,
                              'ocob':7,
                              'wfarenum':8,
                              'wayfare':9,
                              'locality':10,
                              'pcode':11,
                              'state':12,
                              'obmrn':13,
                              'bdob':14,
                              'plural':15,
                              'plurnum':16},
             fields_default = '',
               strip_fields = True,
             missing_values = ['','missing'])

indata_b = DataSetCSV(name = 'MDC-b',
               description = 'Midwifes Data Collection, 1990-2000, ' + \
                             'source: NSW Department of Health',
               access_mode = 'read',
              header_lines = 1,
              write_header = True,
                 file_name = '../../data/nswhealth_mdc/mdc-200000.csv',
                    fields = {'year':0,
                              'rseqnum':1,
                              'ohoscode':2,
                              'omrn':3,
                              'gname':4,
                              'sname':5,
                              'omdob':6,
                              'ocob':7,
                              'wfarenum':8,
                              'wayfare':9,
                              'locality':10,
                              'pcode':11,
                              'state':12,
                              'obmrn':13,
                              'bdob':14,
                              'plural':15,
                              'plurnum':16},
             fields_default = '',
               strip_fields = True,
             missing_values = ['','missing'])

# =============================================================================
# Define temporary data set(s) (one per input data set)

tmpdata_a = DataSetMemory(name = 'MDCtmp-a',
                   description = 'Temporary MDC data set',
                   access_mode = 'readwrite',
                        fields = {'title':0,
                                  'gender_guess':1,
                                  'given_name':2,
                                  'alt_given_name':3,
                                  'surname':4,
                                  'alt_surname':5,
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
                                  'address_hmm_prob':23,
                                  'baby_day':24,
                                  'baby_month':25,
                                  'baby_year':26,
                                  'mother_day':27,
                                  'mother_month':28,
                                  'mother_year':29},
                missing_values = ['','missing'])

tmpdata_b = DataSetMemory(name = 'MDCtmp-b',
                   description = 'Temporary MDC data set',
                   access_mode = 'readwrite',
                        fields = {'title':0,
                                  'gender_guess':1,
                                  'given_name':2,
                                  'alt_given_name':3,
                                  'surname':4,
                                  'alt_surname':5,
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
                                  'address_hmm_prob':23,
                                  'baby_day':24,
                                  'baby_month':25,
                                  'baby_year':26,
                                  'mother_day':27,
                                  'mother_month':28,
                                  'mother_year':29},
                missing_values = ['','missing'])

# =============================================================================
# Define and load lookup tables

name_lookup_table = TagLookupTable(name = 'Name lookup table',
                                default = '')

name_lookup_table.load(['./data/givenname_f.tbl',
                        './data/givenname_m.tbl',
                        './data/name_prefix.tbl',
                        './data/name_misc.tbl',
                        './data/saints.tbl',
                        './data/surname.tbl',
                        './data/title.tbl'])

name_correction_list = CorrectionList(name = 'Name correction list')

name_correction_list.load('./data/name_corr.lst')

address_lookup_table = TagLookupTable(name = 'Geoloc lookup table',
                                   default = '')

address_lookup_table.load(['./data/country.tbl',
                          './data/geoloc_misc.tbl',
                          './data/geoloc_qual.tbl',
                          './data/institution_type.tbl',
                          './data/post_address.tbl',
                          './data/saints.tbl',
                          './data/territory.tbl',
                          './data/unit_type.tbl',
                          './data/wayfare_type.tbl'])

address_correction_list = CorrectionList(name = 'Address correction list')

address_correction_list.load('./data/geoloc_corr.lst')

pc_geocode_table = GeocodeLookupTable(name = 'Postcode centroids',
                                   default = [])

pc_geocode_table.load('./data/postcode_centroids.csv')

# =============================================================================
# Define and load hidden Markov models (HMMs)

name_states = ['titl','baby','knwn','andor','gname1','gname2','ghyph',
               'gopbr','gclbr','agname1','agname2','coma','sname1','sname2',
               'shyph','sopbr','sclbr','asname1','asname2','pref1','pref2',
               'rubb']
name_tags = ['NU','AN','TI','PR','GF','GM','SN','ST','SP','HY','CO','NE','II',
             'BO','VB','UN','RU']

myname_hmm = hmm('Name HMM', name_states, name_tags)
myname_hmm.load_hmm('./hmm/name-absdiscount.hmm')
# myname_hmm.load_hmm('./hmm/name.hmm')
# myname_hmm.load_hmm('./hmm/name-laplace.hmm')

address_states = ['wfnu','wfna1','wfna2','wfql','wfty','unnu','unty','prna1',
                  'prna2','inna1','inna2','inty','panu','paty','hyph','sla',
                  'coma','opbr','clbr','loc1','loc2','locql','pc','ter1',
                  'ter2','cntr1','cntr2','rubb']
address_tags = ['PC','N4','NU','AN','TR','CR','LN','ST','IN','IT','LQ','WT',
                'WN','UT','HY','SL','CO','VB','PA','UN','RU']

myaddress_hmm = hmm('Address HMM', address_states, address_tags)
myaddress_hmm.load_hmm('./hmm/geoloc-absdiscount.hmm')
# myaddress_hmm.load_hmm('./hmm/geoloc.hmm')
# myaddress_hmm.load_hmm('./hmm/geoloc-laplace.hmm')

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
# Define standardisers for dates (one for each data set needed)

mdc_baby_dob_std_a = DateStandardiser(name = 'MDC-BDOB-std-a',
                               description = 'MDC baby DOB standardiser',
                              input_fields = 'bdob',
                             output_fields = ['baby_day','baby_month',
                                              'baby_year'],
                             parse_formats = date_parse_formats)

mdc_baby_dob_std_b = DateStandardiser(name = 'MDC-BDOB-std-b',
                               description = 'MDC baby DOB standardiser',
                              input_fields = 'bdob',
                             output_fields = ['baby_day','baby_month',
                                              'baby_year'],
                             parse_formats = date_parse_formats)

mdc_mother_dob_std_a = DateStandardiser(name = 'MDC-MDOB-std-a',
                                 description = 'MDC mother DOB standardiser',
                                input_fields = 'omdob',
                               output_fields = ['mother_day','mother_month',
                                                'mother_year'],
                               parse_formats = date_parse_formats)

mdc_mother_dob_std_b = DateStandardiser(name = 'MDC-MDOB-std-b',
                                 description = 'MDC mother DOB standardiser',
                                input_fields = 'omdob',
                               output_fields = ['mother_day','mother_month',
                                                'mother_year'],
                               parse_formats = date_parse_formats)

# =============================================================================
# Define a standardiser for names based on rules (one for each data set needed)

mdc_name_rules_std_a = NameRulesStandardiser(name = 'MDC-Name-Rules-a',
                                     input_fields = ['gname','sname'],
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

mdc_name_rules_std_b = NameRulesStandardiser(name = 'MDC-Name-Rules-b',
                                     input_fields = ['gname','sname'],
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
# Define a standardiser for names based on HMM (one for each data set needed)

mdc_name_hmm_std_a = NameHMMStandardiser(name = 'MDC-Name-HMM-a',
                                 input_fields = ['gname','sname'],
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

mdc_name_hmm_std_b = NameHMMStandardiser(name = 'MDC-Name-HMM-b',
                                 input_fields = ['gname','sname'],
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
# Define a standardiser for addresses based on HMM (one for each data set)

mdc_address_hmm_std_a = AddressHMMStandardiser(name = 'MDC-Address-HMM-a',
                                       input_fields = ['wfarenum','wayfare',
                                                       'locality','pcode',
                                                       'state'],
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
                                                       None],
                                  address_corr_list = address_correction_list,
                                  address_tag_table = address_lookup_table,
                                        address_hmm = myaddress_hmm)

mdc_address_hmm_std_b = AddressHMMStandardiser(name = 'MDC-Address-HMM-b',
                                       input_fields = ['wfarenum','wayfare',
                                                       'locality','pcode',
                                                       'state'],
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
                                                       None],
                                  address_corr_list = address_correction_list,
                                  address_tag_table = address_lookup_table,
                                        address_hmm = myaddress_hmm)

# =============================================================================
# Define record standardiser(s) (one for each data set)

mdc_comp_stand_a = [mdc_baby_dob_std_a, mdc_mother_dob_std_a,
                    mdc_name_rules_std_a, mdc_address_hmm_std_a]

mdc_comp_stand_b = [mdc_baby_dob_std_b, mdc_mother_dob_std_b,
                    mdc_name_rules_std_b, mdc_address_hmm_std_b]

# The HMM based name standardisation is not used in the above standardisers,
# uncomment the lines below (and comment the ones above) to use HMM
# standardisation for names.
#
# mdc_comp_stand_a = [mdc_baby_dob_std_a, mdc_mother_dob_std_a,
#                     mdc_name_hmm_std_a, mdc_address_hmm_std_a]
#
# mdc_comp_stand_b = [mdc_baby_dob_std_b, mdc_mother_dob_std_b,
#                     mdc_name_hmm_std_b, mdc_address_hmm_std_b]

mdc_standardiser_a = RecordStandardiser(name = 'MDC-std-a',
                                 description = 'NSW MDC standardiser',
                               input_dataset = indata_a,
                              output_dataset = tmpdata_a,
                                    comp_std = mdc_comp_stand_a)

mdc_standardiser_b = RecordStandardiser(name = 'MDC-std-b',
                                 description = 'NSW MDC standardiser',
                               input_dataset = indata_b,
                              output_dataset = tmpdata_b,
                                    comp_std = mdc_comp_stand_b)

# =============================================================================
# Define blocking index(es) (one per temporary data set)

myblock_def_a = [[('surname','dmetaphone', 4),('mother_year','direct')],
                 [('given_name','truncate', 3), ('postcode','direct')],
                 [('locality_name','nysiis'),('mother_month','direct')],
                ]
myblock_def_b = [[('surname','dmetaphone', 4),('mother_year','direct')],
                 [('given_name','truncate', 3), ('postcode','direct')],
                 [('locality_name','nysiis'),('mother_month','direct')],
                ]

# Define one or more indexes (to be used in the classfier furter below)

mdc_index_a = BlockingIndex(name = 'MDC-blocking-a',
                         dataset = tmpdata_a,
                       index_def = myblock_def_a)

mdc_index_b = BlockingIndex(name = 'MDC-blocking-b',
                         dataset = tmpdata_b,
                       index_def = myblock_def_b)

# =============================================================================
# Define comparison functions for linkage

given_name_nysiis = FieldComparatorEncodeString(fields_a = 'surname',
                                                fields_b = 'surname',
                                                   m_prob = 0.95,
                                                   u_prob = 0.001,
                                           missing_weight = 0.0,
                                            encode_method = 'nysiis',
                                                  reverse = False)

surname_dmetaphone = FieldComparatorEncodeString(fields_a = 'surname',
                                                 fields_b = 'surname',
                                                   m_prob = 0.95,
                                                   u_prob = 0.001,
                                           missing_weight = 0.0,
                                            encode_method = 'dmetaphone',
                                                  reverse = False)

wayfare_name_winkler = FieldComparatorApproxString(fields_a = 'wayfare_name',
                                                   fields_b = 'wayfare_name',
                                                     m_prob = 0.95,
                                                     u_prob = 0.001,
                                             missing_weight = 0.0,
                                             compare_method = 'winkler',
                                           min_approx_value = 0.7)

locality_name_key = FieldComparatorKeyDiff(fields_a = 'locality_name',
                                           fields_b = 'locality_name',
                                             m_prob = 0.95,
                                             u_prob = 0.001,
                                     missing_weight = 0.0,
                                       max_key_diff = 2)

postcode_distance = FieldComparatorDistance(fields_a = 'postcode',
                                            fields_b = 'postcode',
                                              m_prob = 0.95,
                                              u_prob = 0.001,
                                      missing_weight = 0.0,
                                       geocode_table = pc_geocode_table,
                                        max_distance = 50.0)

mother_age = FieldComparatorAge(fields_a = ['mother_day','mother_month',
                                            'mother_year'],
                                fields_b = ['mother_day','mother_month',
                                            'mother_year'],
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
                     postcode_distance, mother_age]

mdc_comparator = RecordComparator(tmpdata_a, tmpdata_b, field_comparisons)

# =============================================================================
# Define a classifier for classifying the matching vectors

mdc_f_s_classifier = FellegiSunterClassifier(name = 'MDC Fellegi and Sunter',
                                        dataset_a = tmpdata_a,
                                        dataset_b = tmpdata_b,
                                  lower_threshold = 0.0,
                                  upper_threshold = 30.0)

# =============================================================================
# Start the linkage task

myproject.link(input_dataset_a = indata_a,
               input_dataset_b = indata_b,
                 tmp_dataset_a = tmpdata_a,
                 tmp_dataset_b = tmpdata_b,
            rec_standardiser_a = mdc_standardiser_a,
            rec_standardiser_b = mdc_standardiser_b,
                rec_comparator = mdc_comparator,
              blocking_index_a = mdc_index_a,
              blocking_index_b = mdc_index_b,
                    classifier = mdc_f_s_classifier,
                first_record_a = 0,
              number_records_a = 40000,
                first_record_b = 0,
              number_records_b = 20000,
                  output_print = True,
              output_histogram = True,
                   output_file = 'mdc-linkage.res',
              output_threshold = 30.0,
             output_assignment = None)

# =============================================================================
