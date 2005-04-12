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
# The Original Software is: "project-standardise.py"
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

"""Module project-standardise.py - Configuration for a standardisation project.

   Briefly, what needs to be defined for a standardisation project is:
   - A Febrl object, a project, plus a project logger
   - One input data set (initialised in read access mode)
   - One output data set (initialised in write or append access mode)
   - Lookup tables to be used
   - Standardisers for names, addresses and dates

   and then the 'standardise' method can be called.

   For more information see chapter

   "Configuration and Running Febrl using a Module derived from 'project.py'"

   in the Febrl manual.

   This project module will standardised the example data set 'dataset1.csv'
   given in the directory 'dsgen' into a data set 'dataset1clean.csv'.

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

init_febrl_logger(log_file_name = 'febrl-example-standard.log',
                     file_level = 'WARN',
                  console_level = 'INFO',
                      clear_log = True,
                parallel_output = 'host')

# =============================================================================
# Set up Febrl and create a new project (or load a saved project)

myfebrl = Febrl(description = 'Example standardisation Febrl instance',
                 febrl_path = '.')

myproject = myfebrl.new_project(name = 'example-std',
                         description = 'Standardise example data set 1',
                           file_name = 'example-standardise.fbr',
                          block_size = 100,
                      parallel_write = 'host')

# =============================================================================
# Define original input data set
#

indata = DataSetCSV(name = 'example1in',
             description = 'Example data set 1',
             access_mode = 'read',
            header_lines = 1,
               file_name = 'dsgen'+dirsep+'dataset1.csv',
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
# Define the output data set

outdata = DataSetCSV(name = 'example1out',
              description = 'Standardised example data set 1',
              access_mode = 'write',
             write_header = True,
                file_name = 'dsgen'+dirsep+'dataset1clean.csv',
         write_quote_char = '',
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
                             'phone_country_code':26,
                             'phone_country_name':27,
                             'phone_area_code':28,
                             'phone_number':29,
                             'phone_extension':30,
# The following are output fields that are passed without standardisation
                             'rec_id':0,
                             'soc_sec_id':31,
# The last output field contains the probability of the address HMM
                             'address_hmm_prob':32,
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
# Define a standardiser for dates

dob_std = DateStandardiser(name = 'DOB-std',
                    description = 'Date of birth standardiser',
                   input_fields = 'date_of_birth',
                  output_fields = ['dob_day','dob_month', 'dob_year'],
                  parse_formats = date_parse_formats)

# =============================================================================
# Define a standardiser for telephone numbers

phone_std = PhoneNumStandardiser(name = 'Phone-Num-std',
                          description = 'Phone number standardiser',
                         input_fields = 'soc_sec_id',
                        output_fields = ['phone_country_code',
                                         'phone_country_name',
                                         'phone_area_code',
                                         'phone_number',
                                         'phone_extension'],
                      default_country = 'australia')

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
                            input_fields = ['rec_id',
                                            'soc_sec_id'],
                           output_fields = ['rec_id',
                                            'soc_sec_id'])

# =============================================================================
# Define record standardiser(s) (one for each data set)

comp_stand = [dob_std, phone_std, name_rules_std, address_hmm_std, pass_fields]

# The HMM based name standardisation is not used in the above standardiser,
# uncomment the lines below (and comment the ones above) to use HMM
# standardisation for names.
#
#comp_stand = [dob_std, name_hmm_std, address_hmm_std, pass_fields]

example_standardiser = RecordStandardiser(name = 'Example-std',
                                   description = 'Example standardiser',
                                 input_dataset = indata,
                                output_dataset = outdata,
                                      comp_std = comp_stand)

# =============================================================================
# Start the standardisation task
# - If 'first_record' is set to 'None' then it will be set to 0
# - If 'number_records' is set to 'None' then it will be set to the total
#   number of records in the input data set.

myproject.standardise(input_dataset = indata,
                     output_dataset = outdata,
                   rec_standardiser = example_standardiser,
                       first_record = 0,
                     number_records = 1000)

# =============================================================================

myfebrl.finalise()

# =============================================================================
