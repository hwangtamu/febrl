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
# The Original Software is: "get-neighbour-regions.py"
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

"""Module get-neighbour-regions.py - Find direct and indirect (level 2)
                                     neighbouring regions extracted from a
                                     pair of MapInfo 'MID' and 'MIF' files.

   This module can be used to create neighbouring region look-up table files
   for postcodes and localities (suburbs, towns) using a pair of MapInfo files.

   It needs three (or four) command line arguments:

     1) MID file name
     2) MIF file name
     3) Output file name with direct neighbours

     4) As an option, a second output file name can be given, into which second
        level neighbouring regions will be saved.

   It is assumed that the first (quoted) element in each line (assuming comma
   separated values) in the MID file are the region names.
"""

# =============================================================================

# For localities (suburb and town names) some special processing is needed
#
correction_dict = {"d aguilar":"d'aguilar",
               "d estrees bay":"d'estrees bay",
               "o briens hill":"o'briens hill",
                   "o connell":"o'connell",
                    "o connor":"o'connor",
             "o halloran hill":"o'halloran hill",
                    "o malley":"o'malley",
                   "o reillys":"o'reillys",
            "o sullivan beach":"o'sullivan beach",
              "stun sail boom":"stun'sail boom",
           "brighton le sands":"brighton-le-sands",
            "callaghans creek":"callaghans creeks",
                 "north nowra":"nowra north",
                "west gosford":"gosford west",
             "west wollongong":"wollongong west",
             "mount st thomas":"mount saint thomas",
               "lake tabourie":"tabourie lake",
            "north wollongong":"wollongong north",
              "south tamworth":"tamworth south",
           "north wagga wagga":"wagga wagga north",
          "south murwillumbah":"murwillumbah south",
               "glendon brook":"glendonbrook",
                "tacoma south":"south tacoma",
               "north dorrigo":"dorrigo north",
            "mount kuring gai":"mount kuring-gai",
                    "st clair":"saint clair",
            "st georges basin":"saint georges basin",
                   "st george":"saint george",
                   "st peters":"saint peters",
              "st helens park":"saint helens park",
                  "st andrews":"saint andrews",
               "st johns park":"saint johns park",
                 "st leonards":"saint leonards",
              "north st marys":"north saint marys",
                    "st marys":"saint marys",
               "st ives chase":"saint ives chase",
                     "st ives":"saint ives",
           "st huberts island":"saint huberts island",
                   "st albans":"saint albans",
                  "st fillans":"saint fillans",
           "st leonards creek":"saint leonards creek",
                 "st patricks":"saint patricks"}

QUOTED = True  # Set to False if region names in MID file are not quoted.

# =============================================================================

import os
import sets
import sys
import time

# =============================================================================

if (len(sys.argv) not in [4,5]):
  print 'Three arguments needed with %s:' % (sys.argv[0])
  print '  - Name of the MID file'
  print '  - Name of the MIF file'
  print '  - Output file name'
  print 'Optional argument:'
  print '  - Second output file name (for second level neighbours)'
  sys.exit()

mid_file_name =    sys.argv[1]
mif_file_name =    sys.argv[2]
output_file_name = sys.argv[3]
if (len(sys.argv) == 5):
  second_output_file_name = sys.argv[4]
else:
  second_output_file_name = None

# -----------------------------------------------------------------------------
# Step 1: Read MID and MIF files, and build sets of all polygon points for each
#         region.

# Step 1a: Open MID file (with region names) and read it into a list  - - - - -
#
mid_file = open(mid_file_name)

mid_data = []
for line in mid_file:
  line = line.lower()
  elem_list = []
  line_list = line.split(',')

  for elem in line_list:
    if (QUOTED == True):
      elem = elem[1:-1]  # Remove quotes
    elem = elem.strip()  # Remove whitespaces

    if elem in correction_dict:
      print '** Correct spelling "%s" into "%s"' % \
            (elem, correction_dict[elem])
      elem = correction_dict[elem]

    if ('-' in elem):  # Replace hyphens with spaces
      elem = elem.replace('-', ' ')

    elem_list.append(elem)

  mid_data.append(elem_list)

mid_file.close()

num_regions = len(mid_data)  # Number of regions to process

print 'MID file "%s" contains %i regions' % (mid_file_name, num_regions)

# Step 1b: Open MIF file and process it region wise - - - - - - - - - - - - - -
#
regions = []

mif_file = open(mif_file_name)

for region_count in range(num_regions):

  region_data = sets.Set()  # Regions are a set of points

  # Read lines until 'Region' keyword appears
  #
  line = mif_file.readline()
  while (line[:6] != 'Region'):
    line = mif_file.readline()

  # Get number of polygons in this region (number after 'Region' keyword)
  #
  num_polygon = int(line.split(' ')[-1])

  print 'Region %i contains %i polygons' % (region_count, num_polygon)

  # Loop through polygons
  #
  for p in range(num_polygon):

    # First line in each polygon contains the number of points in the polygon
    #
    num_polygon_points = int(mif_file.readline())

    print '  Polygon %i-%i contains %i points' % \
          (region_count, p, num_polygon_points)

    for point in range(num_polygon_points):
      point_line = mif_file.readline().strip()  # Get longitute and latitude
      lon, lat = point_line.split(' ')
      lon = float(lon)
      lat = float(lat)

      # print '    (%f,%f)' % (lon, lat)

      # Check if point is already in set for this region
      #
      if (lon,lat) in region_data:
        print '    Point (%f,%f) already in region %i' % \
              (lon, lat, region_count)
      else:
        region_data.add((lon,lat))  # Add new point into set of this region

  regions.append([mid_data[region_count], region_data])

mif_file.close()

print
print

# -----------------------------------------------------------------------------
# Step 2: Loop through all regions and find non-empty intersections of points
#         (i.e. regions that share points), which are neighbouring regions.

neigbours = {}  # A dictionary with the region names as keys

empty_set = sets.Set([])  # Define an empty set

i = 0  # A region counter

start_time = time.time()

for this_region in regions:

  this_region_name = this_region[0][0]

  print 'Finding neighbours for region %s' % (this_region_name)

  region_neighbours = sets.Set()  # Neighbourhoods are sets of region names

  for that_region in regions:
    if (this_region != that_region):  # Don's check same region

      # Find intersection of points of regions
      #
      common_points = this_region[1] & that_region[1]

      if (common_points != empty_set):  # Not an empty set

        that_region_name = that_region[0][0]

        # print 'that region name:', that_region_name

        region_neighbours.add(that_region_name)

  # Make a sorted list of all neighbours for a region and save it
  #
  region_neighbours_list = list(region_neighbours)
  region_neighbours_list.sort()

  neigbours[this_region_name] = region_neighbours_list

  print '  Found %i neighouring regions' % (len(region_neighbours_list))

  i += 1
  if (i % (num_regions / 10) == 0):
    print 'Processed %d regions of %d regions in %.1f seconds' % \
          (i, num_regions, (time.time()-start_time))
    print

print
print

# Print (and save) the list of region neighbours  - - - - - - - - - - - - - - -
#
region_names = neigbours.keys()
region_names.sort()

out_file = open(output_file_name, 'w')

out_file.write('# Region neighbour list' + os.linesep)
out_file.write('#' + os.linesep)
out_file.write('# MID file: %s' % (mid_file_name) + os.linesep)
out_file.write('# MIF file: %s' % (mif_file_name) + os.linesep)
out_file.write('# Output file (this file): %s' % (output_file_name) + \
               os.linesep)
out_file.write('#' + os.linesep)

for region_name in region_names:
  region_neighbours = neigbours[region_name]

  out_file.write('%s: %s' % (region_name, str(region_neighbours)) + \
                 os.linesep)

out_file.close()

print
print

# -----------------------------------------------------------------------------
# Step 3: If optional second output file is given, calculate second level
#         neighbours by looping through region neighbour dictionary.

if (second_output_file_name != None):
  print 'Calculate second level neighbours'

  out_file = open(second_output_file_name, 'w')
  out_file.write('# Region second level neighbour list' + os.linesep)
  out_file.write('#' + os.linesep)
  out_file.write('# MID file: %s' % (mid_file_name) + os.linesep)
  out_file.write('# MIF file: %s' % (mif_file_name) + os.linesep)
  out_file.write('# Output file (this file): %s' % \
                 (second_output_file_name) + os.linesep)
  out_file.write('#' + os.linesep)

  region_names = neigbours.keys()
  region_names.sort()

  for region_name in region_names:

    region_neighbours = neigbours[region_name]  # Direct neighbours

    second_region_neighbours = sets.Set()  # New empty set for second level

    for second_region in region_neighbours:
      that_neighbours = neigbours[second_region]
      for that_region in that_neighbours:
        second_region_neighbours.add(that_region)

    # Remove the central region from the second level neighbours
    #
    second_region_neighbours.discard(region_name)

    # Merge first and second level neighbours
    #
    all_neighbours = second_region_neighbours
    for that_region in region_neighbours:
      all_neighbours.add(that_region)

    all_neighbours_list = list(all_neighbours)
    all_neighbours_list.sort()

    print 'Region name:', region_name
    print '  Direct neighbours:', region_neighbours
    print '  Second neighbours:', list(second_region_neighbours)
    print '  All neighbours:   ', all_neighbours_list

    out_file.write('%s: %s' % (region_name, str(all_neighbours_list)) + \
                   os.linesep)
  out_file.close()

# =============================================================================
