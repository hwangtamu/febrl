# =============================================================================
# lap.py - Routines for linear assignment procedures
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
# The Original Software is "lap.py".
# The Initial Developers of the Original Software are Dr Peter Christen
# (Department of Computer Science, Australian National University), Dr Tim
# Churches (Centre for Epidemiology and Research, New South Wales Department
# of Health) and Drs Markus Hegland, Stephen Roberts and Ole Nielsen
# (Mathematical Sciences Insitute, Australian National University). Copyright
# (C) 2002 the Australian National University and others. All Rights Reserved.
# Contributors:
#
# =============================================================================

"""Module lap.py - Routines for linear assignment procedures

   Currently, the following lap algorithms are implemented:

   - lapmod

     Implementation of the Volgenant LAPMOD algorithm as describe in

     "Linear and Semi-Assignment Problems: A Core Based Approach"
     A. Volgenant, Computers Operation Research, Vol. 23, No. 10, pp. 917-932,
     1996, Elsevier Science Ltd.

   - auction (TO BE DONE)

     Implementation of the asymmtric auction algortihm by Bertsekas as
     described in:

    "..."
"""

# =============================================================================
# Imports go here

# =============================================================================

def do_lap(lap_method, results_dict, process_type, threshold):
  """Linear sum assignment procedure.

     This routine calculates a linear assignments for one-to-one matching.

     The routine do_lap does all kinds of preprocessing, including the
     extraction of unique record pairs (which can be removed before the lap is
     applied) and the extraction of sub-set which can be solved independently.
     These sub-sets are then given to the chosen lap routine.

     The routine takes as input a results dictionary, as produced by a
     classifier (see classification.py), and returns a disctionary with the
     assigned record pair numbers as keys and the corresponding weight as
     values.

     Possible methods are 'lapmod' ('auction' will follow soon, but is not yet
     implemented).

     The process_type attribute must either be set to 'deduplication' or to
     'linkage' in order to be able to preprocess the classifier data prior to
     the lap procedure.
  """

  if (lap_method not in ['lapmod']):
    print 'error:Illegal method for lap method: %s' % (str(lap_method))
    raise Exception

  if (process_type not in ['deduplication', 'linkage']):
    print 'error:Illegal value for attribute "process_type": %s' %\
          (str(process_type))
    raise Exception

  if (results_dict == {}):
    print 'error:Empty results dictionary'
    raise Exception

  # Make sure the threshold is a number if it is defined
  #
  if (threshold != None):
    if (not (isinstance(threshold, int) or isinstance(threshold, float))):
      print 'error:Threshold is not a number: %s' % (str(threshold))
      raise Exception

  lap_results = {}  # Result dictionary with final record pair matches

  print '1:Start linear assignment procedure using method: %s' % (lap_method)

  print '1:  Original length of results dictionary: %i' % (len(results_dict))

  # Step 1: Filter out record pairs with weight lower than the threshold  - - -
  #         (if a threshold is defined)
  #
  if (threshold != None):

    print '1:  Remove record pairs with weight less than: %f' % (threshold)

    work_dict = {}  # Make an empty working copy

    for rec_num in results_dict:  # Loop over all record numbers (keys)
      rec_dict = results_dict[rec_num]  # Get corresponding record dictionary

      new_rec_dict = {}  # Start a new record dictionary

      for rec_num2 in rec_dict:  # Loop over all records in this dictionary
        if (rec_dict[rec_num2] >= threshold):
          new_rec_dict[rec_num2] = rec_dict[rec_num2]  # Copy to new dictionary
        if (new_rec_dict != {}):  # Only insert non empty dictionaries
          work_dict[rec_num] = new_rec_dict

    print '1:  Length of results dictionary after filtering: %i' % \
          (len(work_dict))

  else:  # No threshold is set, use the original results dictionary
    work_dict = results_dict

  results_len = len(work_dict)  # Save results length (after filtering)

  # Step 2: Remove all matches (record pairs) which are unique  - - - - - - - -
  #         (i.e. which don't have matches with other records)
  #
  row_num_dict = {}  # Count occurences of record numbers in rows
  col_num_dict = {}  # Count occurences of record numbers in columns

  for (row_num, rec_dict) in work_dict.items():

    row_num_dict[row_num] = row_num_dict.get(row_num,0)+1

    for col_num in rec_dict:

      # Increase the count for a column number
      #
      col_num_dict[col_num] = col_num_dict.get(col_num,0)+1

    if (process_type == 'deduplication'):

      # For deduplication, insert symmetric record number as well
      #
      row_num_dict[col_num] = row_num_dict.get(col_num,0)+1
      col_num_dict[row_num] = col_num_dict.get(row_num,0)+1

  for rec_num in work_dict.keys():  # Loop over all record numbers (keys)
    rec_dict = work_dict[rec_num]  # Get corresponding record dictionary

    if (len(rec_dict) == 1):  # Only one record pair for this record
      rec_num2 = rec_dict.keys()[0]
      weight   = rec_dict.values()[0]

      if (row_num_dict[rec_num] == 1) and (col_num_dict[rec_num2] == 1):

        lap_results[(rec_num,rec_num2)] = weight  # Insert into results
        del work_dict[rec_num]  # And delete the record in the results

  print '1:  Found and extracted %i unique record ' % (len(lap_results)) + \
        'pairs in results dictionary'

  for rec_pair in lap_results:
    print '1:    %s: %f' % (str(rec_pair), lap_results[rec_pair])
  print '1:'

  # Make a test: Check if no record number appears doubled ##################
  #
  if (process_type == 'deduplication'):
    test_dict = {}
    for (rec_a,rec_b) in lap_results.keys():
      if (test_dict.has_key(rec_a)):
        print '1:########### record '+str(rec_a)+' already in test dictionary'
      else:
        test_dict[rec_a] = 1
      if (test_dict.has_key(rec_b)):
        print '1:########### record '+str(rec_b)+' already in test dictionary'
      else:
        test_dict[rec_b] = 1
  else:  # Linkage process
    test_dict_a = {}
    test_dict_b = {}
    for (rec_a,rec_b) in lap_results.keys():
      if (test_dict_a.has_key(rec_a)):
        print '1:######## record '+str(rec_a)+' already in test dictionary A'
      else:
        test_dict_a[rec_a] = 1
      if (test_dict_b.has_key(rec_b)):
        print '1:######## record '+str(rec_b)+' already in test dictionary B'
      else:
        test_dict_b[rec_b] = 1
  ####### end test ##########################################################

  print '1:  Remaining number of records in results dictionary: %i' % \
        (len(work_dict)) + ' (down from: %i)' % (results_len)

  # Step 3: Find 'connected sub-graphs' in the results dictionary - - - - - - -
  #         (using depth-first search)
  #
  visited =  {}  # Dictionary with will contain all visited rows at the end
  sub_sets = {}  # Dictionary which will contain the sub-sets extracted 

  print '1:  Find connected sub-graphs in results dictionary'

  max_sub_set_length = -1

  for rec_num in work_dict.keys():  # Loop over all rows

    if (not visited.has_key(rec_num)):  # This row has not been visited yet

      visited[rec_num] = rec_num  # Mark visited as 'seeding' row

      print '2:    Create sub-set with seeding record %i' % (rec_num)

      process_queue = [rec_num]  # Start a new queue of rows to process

      row_sub_set = {rec_num:1}  # Row numbers connected to this row

      while (process_queue != []): # Process rows until all connected rows done
        print '3:      Process queue: %s' % (str(process_queue))

        next_row = process_queue.pop(0)  # Get first row to process
        row_col_numbers = work_dict[next_row].keys()  # Get columns in this row

        # For deduplication, also insert row number into column numbers
        #
        if (process_type == 'deduplication'):
          row_col_numbers.append(next_row)   ################ NEW

        print '2:      Row %i with column numbers: %s' % \
              (next_row, str(row_col_numbers))

        # Get the row numbers from all column numbers
        #
        for col_num in row_col_numbers:

          row_num_list = []  # Create a list of rows having this column number

          if (process_type == 'deduplication') and (col_num in work_dict):
            row_num_list.append(col_num) ############ NEW

         # # For deduplication, also insert column number into sub-set
         # #
         # if (process_type == 'deduplication') and (col_num in work_dict):
         #   row_sub_set[col_num] = 1
         #   if (col_num not in visited):
         #     visited[col_num] = rec_num  # And mark is as visited
         #   if (col_num not in process_queue):
         #     process_queue.append(col_num)

          for rec_num2 in work_dict:  # Loop over all rows
            if (col_num in work_dict[rec_num2]):
              row_num_list.append(rec_num2)

          print '2:        Column: %i with row numbers: %s' % \
                (col_num, str(row_num_list))

          for row_num in row_num_list:
            row_sub_set[row_num] = 1
            if (not visited.has_key(row_num)):  # Check if it's a new row
              process_queue.append(row_num)
              print '3:        Appended row number %i to process queue' % \
                    (row_num)

              visited[row_num] = rec_num  # Mark row as visited by seeding row
              print '2:        Row %i connected to row %i' % (row_num, rec_num)

      sub_sets[rec_num] = row_sub_set.keys()  # Only store keys

      if (len(row_sub_set) > max_sub_set_length):
        max_sub_set_length = len(row_sub_set)

      print '2:      Sub-set contains records: %s' % (str(row_sub_set.keys()))

  print '1:  Number of sub-sets extracted: %i' % (len(sub_sets))
  print '1:    Longest sub-set contains %i rows' % (max_sub_set_length)

  # Test if all sub-sets are mutually exclusive ############################
  #
  for (seed_row, row_list) in sub_sets.items():

    for rec_num in row_list:

      for (seed_row2, row_list2) in sub_sets.items():

        if (seed_row != seed_row2):  # Don't test itself
          if (rec_num in row_list2):
            print '1:  ######  error: '+str(seed_row2)+' / '+str(row_list2)
  ########## End test ##############3

  # Now call the actual linear assignment method for each of the sub-sets - - -
  #
  for (seed_row, row_list) in sub_sets.items():

    print '1:    Sub-set seed row:  %s' % (str(seed_row))
    print '2:    Sub-set row list:  %s' % (str(row_list))

    # Construct dictionary for this sub-set and get minimal and maximal weights
    #
    min_weight = 9999.0
    max_weight = -999.0

    lap_dict = {}
    for row in row_list:  # Loop over rows
      lap_dict[row] = work_dict[row]
      row_dict_values = lap_dict[row].values()
      min_val = min(row_dict_values)
      max_val = max(row_dict_values)
      if (min_val < min_weight):
        min_weight = min_val
      if (max_val > max_weight):
        max_weight = max_val

    print '2:    Sub-set dictionary: %s' % (str(lap_dict))

    print '1:    Minimal and maximal weights in results: %f / %f' % \
          (min_weight, max_weight)

    if (lap_method == 'lapmod'):
      new_lap_results = lapmod(lap_dict, min_weight, max_weight, process_type)

    else:
       print 'error:LAP method %s not implemented' % (lap_method)
       raise Exception

    lap_results.update(new_lap_results)  # Insert into final results

  print '1:  Total number of assignments: %i' % (len(lap_results))
  print '1:  Number of rows in original results dictionary: %i' % \
        (len(results_dict))

  # Make a test: Check if no record number appears doubled
  #
  test_dict = {}
  for (rec_a,rec_b) in lap_results.keys():
    if (test_dict.has_key(rec_a)):
      print '1: ############ record '+str(rec_a)+' already in test dictionary'
    else:
      test_dict[rec_a] = 1
    if (test_dict.has_key(rec_b)):
      print '1: ############ record '+str(rec_b)+' already in test dictionary'
    else:
      test_dict[rec_b] = 1
  ########## end test ####################333

  return lap_results

# =============================================================================

def lapmod(result_dict, min_weight, max_weight, process_type):
  """Linear sum assignment procedure LAPMOD for rectangular problems.

     Re-implementation of a C++ and Pascal codes, taken from: 

     http://www.magiclogic.com/assignment.html

     and as described in:

     "Linear and Semi-Assignment Problems: A Core Based Approach"
     A. Volgenant, Computers Operation Research, Vol. 23, No. 10, pp. 917-932,
     1996, Elsevier Science Ltd.


     Takes as input a results dictionary, the minimum and maximum weights in
     this dictionary, and the process type (which must either be set to
     'deduplication' or 'linkage' in order to be able to process the results
     dictionary correctly).
  """

  max_cost = 1000000  # A value larger than max_weight

  lap_results = {}  # Result dictionary with final record pair matches

  # Make sure minimal weights don't become negative - - - - - - - - - - - - - -
  #
  max_weight += 1
  if (min_weight < 0):
    max_weight += -min_weight
  print '1:    Set maximal weight to: %f' % (max_weight)

  # For a deduplication process, insert symmetric elements  - - - - - - - - - -
  #
  if (process_type == 'deduplication'):
    for (rec_num_a, rec_dict) in result_dict.items():
      for (rec_num_b, weight) in rec_dict.items():  # All elements in a row
        rec_dict2 = result_dict.get(rec_num_b,{})

        if (rec_num_a not in rec_dict2):  # Only insert if not already there
          rec_dict2[rec_num_a] = weight  # Insert symmetric weight
          result_dict[rec_num_b] = rec_dict2
 
    print '1:  Results dictionaries made symmetric for deduplication process'
    print '1:    New results length: %i' % (len(result_dict))
    print '2:      %s' % (str(result_dict))

  else:  # Processing for linkage process - - - - - - - - - - - - - - - - - - -
    pass  # What to do here? #################

    # What happens if result_dict has more rows than columns ??? ####

  # Get the row and column numbers  - - - - - - - - - - - - - - - - - - - - - -
  #
  row_numbers = result_dict.keys()
  col_numbers = {}
  for row_dict in result_dict.values():
    col_numbers.update(row_dict)
  col_numbers = col_numbers.keys()

  # Create sparse cost matrix from results dictionary - - - - - - - - - - - - -
  #
  assign_cost_row = {}  # Row oriented sparse matrix

  for rec_num_row in result_dict:  # Loop over all keys in result dictionary

    rec_dict = result_dict[rec_num_row]  # All entries for this record
    row_dict = {}  # Dictionary with all non-zero elements in a row

    for (rec_num_col, w) in rec_dict.items():
      cost = max_weight - w  # Make smaller cost if original weight is large
      row_dict[rec_num_col] = cost

    # Insert a diagonal element if there is none
    #
    if (rec_num_row not in row_dict):  # No diagonal element in this row
      row_dict[rec_num_row] = 100.0*max_weight
      print '3:    Inserted diagonal element into row %i' % (rec_num_row)

      if (rec_num_row not in col_numbers):
        col_numbers.append(rec_num_row)  # Append new column number

    # Make sure the row has at least two elements - - - - - - - - - - - - - -
    # (insert an artificial element with a high value if necessary)
    #
    if (len(row_dict) == 1):
      i = rec_num_row+1  # Start with element after diagonal element
      max_col_number = max(col_numbers)

      while (i not in col_numbers):  # Find an available column number
        i += 1
        if (i > max_col_numbers):
          i = 0  # Start from first column

        if (not row_dict.has_key(i)):  # Should always be
          row_dict[i] = 50.0*max_weight
        else:
          print 'error:This should never happen'
          raise Exception

    assign_cost_row[rec_num_row] = row_dict  # Insert into row sparse matrix

  print '1:  Length of dictionary used for lapmod: %i' % (len(assign_cost_row))
  print '1:  Dictionary used for lapmod:'
  print '1:   %s' % (str(assign_cost_row))

  # Get sorted lists of all row and column numbers and dimensionality - - - - -
  #
  row_numbers.sort()
  col_numbers.sort()

  num_rows = len(row_numbers)
  num_cols = len(col_numbers)

  # Check if the number of rows equals to the number of columns
  #
  symmetric = (num_rows == num_cols)

  # Create several arrays needed  - - - - - - - - - - - - - - - - - - - - - - -
  #
  u = {}  # Dual variables, row reduction numbers, u[i] for row i
  v = {}  # Dual variables, column reduction numbers, v[j] for column j
  x = {}  # Column index assigned to row i, x[i] = 0 for unassigned row
  y = {}  # Row index assigned to column j, y[j] = 0 for unsassigne col

  lab =  {}
  todo = []  # List with row numbers to be processed (?)
  free_rows = []

  # Initialisation  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  #
  for rn in row_numbers:
    x[rn] = 0
    u[rn] = 0.0
    todo += [0,0]  # Must be longer than the number of rows
    free_rows += [0]

  for cn in col_numbers:
    y[cn] = 0
    v[cn] = max_cost

  # Find minimum costs over rows  - - - - - - - - - - - - - - - - - - - - - - -
  #
  for row_num in row_numbers:

    row_dict = assign_cost_row[row_num]  # Get all elements in this row

    for (col_num, cost) in row_dict.items():
      if (cost < v[col_num]):  # Find minimal cost row
         v[col_num] = cost
         y[col_num] = row_num

  # Column reduction (for symmetric problems only)  - - - - - - - - - - - - - -
  #
  if (symmetric):
    rev_col_numbers = col_numbers[:]  # Make a copy of the column numbers
    rev_col_numbers.reverse()  # Reverse order gives better results

    for col_num in rev_col_numbers:
      i = y[col_num]
      if (x[i] == 0):  # Row not assigned yet
        x[i] = col_num
      else:  # Row already assigned to another column
        y[col_num] = 0  # Mark column as not assigned
        x[i] = -abs(x[i])  # Mark row as 'double' assigned
  else:
    print '1:  Asymmetric problem, no column reduction is done'

  print '1:  Initialisation and column reduction done'
  print '2:    x= %s' % (str(x))
  print '2:    u= %s' % (str(u))
  print '2:    y= %s' % (str(y))
  print '2:    v= %s' % (str(v))

  # Reduction transfer (for symmetric problems only)  - - - - - - - - - - - - -
  #
  if (symmetric):
    l = 0  # Counter for free rows
    for i in row_numbers:
      if (x[i] < 0):
        x[i] = -x[i]
      elif (x[i] > 0):
        min_val = max_cost
        j1 = x[i]
    
        row_dict = assign_cost_row[i]  # Get all elements in this row

        for (col_num, cost) in row_dict.items():

          if (col_num != j1):
            if ((cost - v[col_num]) < min_val):
              min_val = cost - v[col_num]
        u[i] = min_val
        v[j1] = row_dict[j1] - min_val

      else:  # Free row
        free_rows[l] = i
        l += 1

  else:  # Asymmetric problem, add all rows into list of free rows
    print '1:  Asymmetric problem, all rows entered into free rows list'
    l = 0
    for rn in row_numbers:
      free_rows[l] = rn
      l +=1

  print '1:  Reduction transfer done'
  print '2:    x= %s' % (str(x))
  print '2:    u= %s' % (str(u))
  print '2:    y= %s' % (str(y))
  print '2:    v= %s' % (str(v))
  print '2:    l= %s' % (str(l))
  print '2:    free rows: %s' % (str(free_rows))

  # Augmenting row reduction  - - - - - - - - - - - - - - - - - - - - - - - - -
  #
  for arr in range(2):  # ARR performs best with two loops (see paper)

    h =  0
    l0 = l
    l =  -1

    while (h < l0):

      i = free_rows[h]  # Get row number of free row
      h += 1
      v0 = max_cost  # Minimum
      vj = max_cost  # Second minimum
      j0 = -1  # Column index of minimum (illegal initial value)

      row_dict = assign_cost_row[i]  # Get all elements in this row

      for (col_num, cost) in row_dict.items():  # Find min. and second min.

        dj = cost - v[col_num]  # Get element from row
        #print '1:  col_num:'+str(col_num)+', v[col_num] '+str(v[col_num])
        #print '1:  dj:'+str(dj)+', cost '+str(cost)
        if (dj < vj):  # It's smaller than the second minimum
          if (dj >= v0):  # It's larger than the minimum
            vj = dj
            j1 = col_num  # Update second minimum and its column index
          else:  # Found new minimum, move old minimum to second minimum
            vj = v0
            v0 = dj  # New minimum value
            j1 = j0
            j0 = col_num  # New index for minimum value

      #print '1:row '+str(i)+': '+str(row_dict)
      #print '1:min='+str(v0)+','+str(j0)
      #print '1:min2='+str(vj)+','+str(j1)

      i0 = y[j0]
      # u[i] = vj

      if (v0 < vj):  # Change the reduction of the minimum column to increase
                     # the minimum reduced cost in the row to the subminimum
        v[j0] = v[j0] - vj + v0
      elif (i0 > 0):  # Minimum and subminimum equal, and column i0 assigned
        j0 = j1  # Swap columns j0 and j1 as j1 may be unassigned
        i0 = y[j0]

      x[i] = j0  # (Re-)assign i to j0, possible de-assigning an i0
      y[j0] = i
      if (i0 > 0):  # Minimum column j2 assigned earlier
        if (v0 < vj):
          h -= 1  # Go back to previous free row
          free_rows[h] = i0
        else:
          l += 1
          free_rows[l] = i0  # Store i0 in free rows list for next phase

  print '1:  Augmenting row reduction done'
  print '2:    x= %s' % (str(x))
  print '2:    u= %s' % (str(u))
  print '2:    y= %s' % (str(y))
  print '2:    v= %s' % (str(v))
  print '2:    l= %s' % (str(l))
  print '2:    free rows: %s' % (str(free_rows))

  td1 = -1 # Initalise to do index counter 1
  goto_flag = ''  # Flag for goto statement (see original Pascal code)

  # Augment solution for each free row  - - - - - - - - - - - - - - - - - - - -
  #
  l0 = l
  for l in range(l0):
    d =  {}
    ok = {}

    min_val = max_cost
    i0 = free_rows[l]  # Number of the current free row

    row_dict = assign_cost_row[i0]  # Get all elements in this row

    for (col_num, cost) in row_dict.items():

      dj = cost - v[col_num]
      d[col_num] = dj
      lab[col_num] = i0

      if (dj <= min_val):
        if (dj < min_val):
          td1 = -1  # Initalise to do index counter
          min_val = dj
        td1 += 1
        todo[td1] = col_num  # Insert into to do list
        # print '1: td1:'+str(td1)+', todo:'+str(todo) ####

    # print '1:  td1:'+str(td1)+', todo:'+str(todo) ####

    for h in range(td1+1):  # Loop over rows to be done

      j = todo[h]  ####
      if (y[j] == 0):
        # print 'goto 2 statement :-('  ####
        goto_flag = '2'
        break
      ok[j] = True

    if (goto_flag != '2'):
      td2 = num_rows
      last = num_rows+1

    while (goto_flag == ''):
      # print '1: while 1: goto_flag: '+goto_flag
      if (goto_flag != ''):
        break  # Get out of here...

      # print '1:  td1:'+str(td1) ####
      j0 = todo[td1]
      td1 -= 1
      # print '1:  j0:'+str(j0)+', todo:'+str(todo)
      # print '1:  i:'+str(i)+', y:'+str(y)
      i = y[j0]
      todo[td2] = j0
      td2 -= 1

      row_dict = assign_cost_row[i]  # Get contents of this row

      h = row_dict[j0] - v[j0] - min_val

      for (col_num, cost) in row_dict.items():

        if (not ok.get(col_num, False)):
          vj = cost - v[col_num] - h

          if (vj < d.get(col_num, max_cost)):
            d[col_num] =   vj
            lab[col_num] = i

            if (vj == min_val):
              if (y[col_num] == 0):
                # print 'goto 1 statement :-(' ###
                goto_flag = '1'
                break
              td1 += 1
              todo[td1] = col_num
              ok[col_num] = True

      # print '1: after for: goto_flag: '+goto_flag
      if (goto_flag == '1'):
        break  # Get out of here...

      # print '1:td1: '+str(td1)+', td2: '+str(td2) ####
      if (td1 == -1):
        # print '1:I am here'
        #print '1:  '+str(d)
        #print '1:  '+str(ok)
        #print '1:  '+str(row_numbers)

        min_val = max_cost-1
        last = td2+1
        for j in row_numbers:  # Loop over all rows
          #print '1:*    '+str(j)
          #print '1:*    '+str(d.get(j, max_cost))
          #print '1:*    '+str(ok.get(j, False))
          #print '1:  td1='+str(td1)

          if (d.get(j, max_cost) <= min_val):
            if (not ok.get(j, False)):
              # print '1:And here'
              if (d.get(j, max_cost) < min_val):
                td1 = -1
                min_val = d.get(j, max_cost)
              td1 += 1
              todo[td1] = j

        # print '1:### '+str(td1)+', '+str(todo) #####

        for h in range(td1+1):
          j = todo[h]
          if (y[j] == 0):
            # print 'goto 1 statement :-('  #####
            goto_flag = '1'
            break
          ok[j] = True

    # End while 1 loop
    # print 'exit' #####

    if (goto_flag == '1'):
      for k in range(last, num_rows):
        j0 = todo[k]
        v[j0] = v[j0] + d[j0] - min_val
    if (goto_flag == '2'):
      i = lab[j]
      y[j] = i
      k = j
      j = x[i]
      x[i] = k
      while (i != i0):
        i = lab[j]
        y[j] = i
        k = j
        j = x[i]
        x[i] = k

#  h = 0.0
#  for i in row_numbers:
#    j = x[i]
#
#    row_dict = assign_cost_row[i]  # Get contents of this row
#   
#    u[i] = row_dict[j] - v[j]
#    h += row_dict[j]

  if (process_type == 'deduplication'):  # Some post-processing needed

    # A dictionary of potential assignments for this sub-set
    #
    assign_pairs = {}

    for (x_key, x_val) in x.items():  # Loop over all assignments in x
      if (x_key > x_val):  # Make sure x_key is smaller than x_val
        x_key, x_val = x_val, x_key

      # Get original weight from results dictionary
      #
      if ((x_key, x_val) not in assign_pairs) and \
         (result_dict.has_key(x_key)):
        rec_dict = result_dict[x_key]
        if (rec_dict.has_key(x_val)):  # Found a valid record pair
          weight = rec_dict[x_val]

          assign_pairs[(x_key, x_val)] = weight

    assign_weights = assign_pairs.values()  # Get the weights in a list
    assign_rec_pairs = assign_pairs.keys()  # And the record pairs

    assign_list = map(None,assign_weights,assign_rec_pairs)
    assign_list.sort()

    new_lap_results = {}
    rec_num_list = []  # A list with record numbers found in assigned pairs

    while (assign_list != []):  # Now get record pairs with largest weight
      check_pair = assign_list.pop()  # Get record pair with largest weight

      rec_a = check_pair[1][0]
      rec_b = check_pair[1][1]

      # Now make sure no record has been assigned previously
      #
      if (rec_a not in rec_num_list) and (rec_b not in rec_num_list):

        rec_num_list += [rec_a, rec_b]

        if (not lap_results.has_key((rec_a, rec_b))) or \
           (not new_lap_results.has_key((rec_a, rec_b))):
          new_lap_results[(rec_a, rec_b)] = check_pair[0]
        else:
          print 'warning:Record pair (%i,%i) already in LAP results' % \
                (rec_a,rec_b)

    print '1: new LAP results: %s' % (str(new_lap_results))
    # print '1: rec number list: '+str(rec_num_list)

  else:  # For linkage process

    print '1:    result_dict:' +str(result_dict)
    print '1:    x= %s' % (str(x))

    new_lap_results = {}
    for (x_key, x_val) in x.items():  # Loop over all assignments in x
      if ((x_key, x_val) in lap_results):
        print 'warning:Record pair (%i,%i) already in LAP results' % \
              (x_key, x_val)
      else:
        if (result_dict.has_key(x_key)):
          rec_dict = result_dict[x_key]
          if (rec_dict.has_key(x_val)):  # Found a valid record pair
            weight = rec_dict[x_val]

            new_lap_results[(x_key, x_val)] = weight  # Insert into LAP results
            print '1:      Inserted (%i,%i) with weight %f into results' % \
                  (x_key, x_val, weight)
          else:
            print 'warning:Key %i not in results dictionary for record %i' % \
                  (x_val, x_key)
        else:
          print 'warning:Key %i not in results dictionary' % (x_key)

  lap_results.update(new_lap_results)  # Insert into final results dictionary

  print '1:  Added %i pairs to assignments' % (len(new_lap_results))

  return lap_results

# =============================================================================
