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
# The Original Software is: "reverse-gnaf.py"
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

"""Module reverse-gnaf.py - Module for reverse G-NAF look-ups.

   This module allows reverse look-ups by giving one or more G-NAF persistent
   identifiers (PIDs) as arguments, retrieving the corresponding original G-NAF
   records.

   The name of the reverse look-up shelve (as created with process-gnaf.py) has
   to be defined in the variable 'shelve_name'.

   USAGE:

     python reverse-gnaf.py PID_1 [PID_2 PID_3 ...]
"""

# =============================================================================

import shelve
import sys

# =============================================================================
# Define the name of the reverse G-NAF look-up shelve to be used.
#
shelve_name = '..'+dirsep+'..'+dirsep+'..'+dirsep+'data'+dirsep+'gnaf' + \
              dirsep+'shelve_pickles'+dirsep+'gnaf_reverse_index.slv'

# =============================================================================

pid_list = sys.argv[1:]  # Get all the PIDs from the command line

gnaf_shelve = shelve.open(shelve_name,'r')

for pid in pid_list:

  if (pid[:3] == '500'):
    pid_type = 'Locality'
    search_keys = ['G_LOCALITY', 'G_LOCALITY_ALIAS']
  elif (pid[:3] == '502'):
    pid_type = 'Street'
    search_keys = ['G_STREET', 'G_STREET_LOCALITY_ALIAS']
  elif (pid[:3] >= '699'):
    pid_type = 'Address site'
    search_keys = ['G_LOCALITY', 'G_LOCALITY_ALIAS', 'G_STREET', \
                   'G_STREET_LOCALITY_ALIAS', 'G_ADDRESS_DETAIL']

  print 'Input PID %s (%s):' % (pid, pid_type)

  result_dict = gnaf_shelve.get(pid, None)

  if (result_dict == None):
    print '  Not in G-NAF reverse shelve'

  else:

    for (record, file_name) in result_dict.items():

      if (file_name in search_keys):  # Limit the files to search and print

        print '  ***** From file %s *****' % (file_name)

        file_attributes = gnaf_shelve[file_name]

        record_list = record.split(',')

        for i in range(len(record_list)):
          if (record_list[i].strip() != ''):
            attr_name = file_attributes[i]
            attr_valu = record_list[i].strip().lower()
            if ('PID' in attr_name) and (attr_valu[-3:] == '.00'):
              attr_valu = attr_valu[:-3]
            print '   %20s : %s' % (attr_name, attr_valu)

      if (pid_type == 'Address site') and (file_name == 'G_ADDRESS_DETAIL'):
        loc_pid = record_list[file_attributes.index('LOCALITY_PID')]
        if (loc_pid[-3:] == '.00'):
          loc_pid = loc_pid[:-3]

        loc_dict = gnaf_shelve.get(loc_pid, None)

        if (loc_dict == None):
          print '    No LOCALITY_PID detail information found'
        else:
          for (loc_record, loc_file_name) in loc_dict.items():
            if (loc_file_name in ['G_LOCALITY', 'G_LOCALITY_ALIAS']):

              print '     ***** From file %s *****' % (loc_file_name)

              loc_file_attributes = gnaf_shelve[loc_file_name]
              loc_record_list = loc_record.split(',')

              for i in range(len(loc_record_list)):
                if (loc_record_list[i].strip() != ''):
                  print '      %20s : %s' % (loc_file_attributes[i], \
                        loc_record_list[i].strip().lower())

        street_pid = record_list[file_attributes.index('STREET_PID')]
        if (street_pid[-3:] == '.00'):
          street_pid = street_pid[:-3]

        street_dict = gnaf_shelve.get(street_pid, None)

        if (street_dict == None):
          print '    No STREET_PID detail information found'

        else:
          for (street_record, street_file_name) in street_dict.items():
            if (street_file_name in ['G_STREET', 'G_STREET_LOCALITY_ALIAS']):

              print '     ***** From file %s *****' % (street_file_name)

              street_file_attributes = gnaf_shelve[street_file_name]
              street_record_list = street_record.split(',')

              for i in range(len(street_record_list)):
                if (street_record_list[i].strip() != ''):
                  print '      %20s : %s' % (street_file_attributes[i], \
                        street_record_list[i].strip().lower())

  print

gnaf_shelve.close()

# =============================================================================
