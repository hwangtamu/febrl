# =============================================================================
# date.py - Routines for date conversions and parsing.
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
# The Original Software is "date.py".
# The Initial Developers of the Original Software are Dr Peter Christen
# (Department of Computer Science, Australian National University), Dr Tim
# Churches (Centre for Epidemiology and Research, New South Wales Department
# of Health) and Drs Markus Hegland, Stephen Roberts and Ole Nielsen
# (Mathematical Sciences Insitute, Australian National University). Copyright
# (C) 2002 the Australian National University and others. All Rights Reserved.
# Contributors:
#
# =============================================================================

"""Module date.py - Routines for date conversions and parsing.

   PUBLIC FUNCTIONS:
     epoch_to_date        Convert a Unix epoch day number into a date
     date_to_epoch        Convert a date into a Unix epoch day integer
     date_to_age          Convert a date into an age (relative to a fix date)
     str_to_date          A routine that converts a string into a date using
                          a format string
     get_today            Return today's date as a [day,month,year] tuple.

   Note that dates are returned as a tuple of strings, with day and month being
   of length 2 (i.e. '01' etc.), and year being of length 4 (e.g. '2003').

   See doc strings of individual functions for detailed documentation.

   TODO:
   - PC 30/11/2002  Fix date_to_epoch and epoch_to_date, which currently do not
                    allow processing of dates before 1901.
"""

# =============================================================================
# Imports go here

import string
import time

# =============================================================================

# A dictionary of month name abbreviations, used in date.str_to_date() routine

month_abbrev_dict = {'jan':'01', 'feb':'02', 'mar':'03', 'apr':'04', \
                     'may':'05', 'jun':'06', 'jul':'07', 'aug':'08', \
                     'sep':'09', 'oct':'10', 'nov':'11', 'dec':'12'}

# Define a character replace table for data strings - - - - - - - - - - - -
#
string_replace = ["'.,:-=_/\\", \
                  "         "]
# Characters in the first list are replaced by the corresponding character in
# the second list

replace_table = string.maketrans(string_replace[0], string_replace[1])

# =============================================================================

def epoch_to_date(daynum):
  """Convert a Unix epoch day number into a date [day, month, year], with
     day, month and year being strings of length 2, 2, and 4, respectively.

  USAGE:
    [year, month, day] = epoch_to_date(daynum)

  ARGUMENTS:
    daynum  A integer giving the Unix epoch day (0 = 1970-01-01)

  DESCRIPTION:
    Function for converting a number of days since Unix epoch time (integer
    value) into a date tuple [day, month, year].

  EXAMPLES:
    [day, month, year] = epoch_to_date(0)       # 1970-01-01
    [day, month, year] = epoch_to_date(11736)   # 2002-02-18
  """

  if (not (isinstance(daynum, int) or isinstance(daynum, long))):
    print 'error:Input value for "daynum" is not of integer type: %s' % \
          (str(daynum))
    raise Exception

  date_tuple  = time.gmtime(daynum*24*3600)

  # A log message for high volume log output (level 3)  - - - - - - - - - - - -
  #
  print '3:    Epoch: %i -> Date: %s' % (daynum, str(date_tuple))

  day =   string.zfill(str(date_tuple[2]),2)  # Add '0' if necessary
  month = string.zfill(str(date_tuple[1]),2)  # Add '0' if necessary
  year =  str(date_tuple[0])  # Is always four digits long

  return [day, month, year]

# =============================================================================

def date_to_epoch(day, month, year):
  """ Convert a date [day, month, year] into a Unix epoch day number.

  USAGE:
    daynum = date_to_epoch(year, month, day)

  ARGUMENTS:
    day    Day value (string or integer number)
    month  Month value (string or integer number)
    year   Year value (string or integer number)

  DESCRIPTION:
    Function for converting a date into a Unix epoch day number
    (integer value).

    Based on a Perl script... source unknown

  EXAMPLES:
    day1 = date_to_epoch(18,  2, 2002)  # 11736
    day2 = date_to_epoch('01', '01', '1970')  # 0
  """

  try:
    day_int = int(day)
  except:
    print 'error:"day" value is not an integer'
    raise Exception
  try:
    month_int = int(month)
  except:
    print 'error:"month" value is not an integer'
    raise Exception
  try:
    year_int = int(year)
  except:
    print 'error:"year" value is not an integer'
    raise Exception

  if (day_int <= 0) or (day_int > 31):
    print 'error:Input value for "day" is not a possible day number: %i' % \
          (day)
    raise Exception
  if (month_int <= 0) or (month_int > 12):
    print 'error:Input value for "month" is not a possible day number: %i' % \
          (month)
    raise Exception
  if (year_int <= 1000):
    print 'error:Input value for "year" is not a positive integer ' + \
          'number: %i' % (year)
    raise Exception

  # Note; mktime did not work on Tim's RedHat 8.0 Linux system running
  # Python 2.2.1, a date of 1/9/1968 was too early and mktime complained
  # about the values in the tuple.

  #epoch_time = time.mktime((year_int, month_int, day_int, 0,0,0,0,0,0)) - \
  #             time.timezone
  #epoch_date = int(epoch_time / (24 * 3600))

  # Do some adjustments for leap year etc.
  #
  if (month_int < 3):
    year_int = year_int - 1
  if (month_int > 2):
    month_int = month_int - 3
  else:
     month_int = month_int + 9

  c = year_int / 100.0
  ya = year_int - ( 100 * c )

  epoch_date = int(((146097*c)/4) + ((1461*ya)/4) + \
                  (((153*month_int)+2)/5) + day_int - 719469)

  if (epoch_date < 0):
    epoch_date -= 1  # Adjust for pre 1970 dates ##### Still wrong #####

  #print
  #print epoch_date
  #print epoch_date2

  # A log message for high volume log output (level 3)  - - - - - - - - - - - -
  #
  print '3:    Date: %s -> Epoch: %i' % \
        (str([day_int,month_int,year_int]), epoch_date)

  return epoch_date

# =============================================================================

def date_to_age(day, month, year, fix_date='today'):
  """Convert a date into an age (relative to a fix date)

  USAGE:
    age = date_to_age(day, month, year)
    age = date_to_age(day, month, year, fix_date)

  ARGUMENTS:
    day       Day value (integer number)
    month     Month value (integer number)
    year      Year value (integer number)
    fix_date  The date relative for which the age is computed. Can be a date
              tuple, the string 'today' (which is the default), or an integer
              (epoch day number)

  DESCRIPTION:
    Returns the age in years as a positive floating-point number.
    If the date is after the fix date a negative age is returned.
  """

  # Check if fix date is given, otherwise calculate it  - - - - - - - - - - - -
  #
  if (fix_date == 'today'):
    sys_time = time.localtime(time.time())  # Get current system date
    fix_day =   string.zfill(str(sys_time[2]),2)
    fix_month = string.zfill(str(sys_time[1]),2)
    fix_year =  str(sys_time[0])

  elif (isinstance(fix_date, list) and (len(fix_date) == 3)):
    fix_day =   string.zfill(str(fix_date[0]),2)
    fix_month = string.zfill(str(fix_date[1]),2)
    fix_year =  str(fix_date[2])

  elif (isinstance(fix_date, int)):
    fix_epoch = fix_date

  else:
    print 'error:"fix_date" is not in a valid form: %s' % (str(fix_date))
    raise Exception

  # Get epoch number for input date and fix date  - - - - - - - - - - - - - - -
  #
  date_epoch = date_to_epoch(day, month, year)

  if (not isinstance(fix_date, int)):
    fix_epoch  = date_to_epoch(fix_day, fix_month, fix_year)

  day_diff = fix_epoch - date_epoch  # Get day difference

  # Compute approximate age - - - - - - - - - - - - - - - - - - - - - - - - - -
  #
  age = float(day_diff) / 365.25  # Can be positive or negative

  # A log message for high volume log output (level 3)  - - - - - - - - - - - -
  #
  print '3:    Date: %s with fix date: %s -> Age: %.3f' % \
        (str([day,month,year]), str(fix_date), age)

  return age

# =============================================================================

def str_to_date(date_str, format_str, pivot_year):
  """A routine that converts a string into a date using a format string.

  USAGE:
    [day,month,year] = str_to_date(datestr, formatstr)

  ARGUMENTS:
    datestr     Input date as a string
    formatstr   A format string, must be made of three directives, which can
                either be written one after the other or being separated by
                a space between them (e.g. "%d%b%Y" or "%m %d %y")
                Possible format directives are (similar to Python 'strptime'
                format directives):
                  %b  Abbreviated month name (Jan, Feb, Mar, etc.)
                  %B  Full month name (January, February, etc.)
                  %d  Day of the month as a decimal number [01,31]
                  %m  Month as a decimal number [01,12]
                  %y  Year without century as a decimal number [00,99]
                  %Y  Year with century as a decimal number
    pivot_year  If a two-digit year is given, the pivot year is used to
                detemine if it is expanded into 19XX or 20XX. Two-digits years
                smaller than the pivot year will be expanded into 20XX, years
                larger and equal than the pivot year will be expanded into 19xx
                Example: pivot_year = 03:  68 -> 1968
                                           03 -> 1903
                                           02 -> 2002

  DESCRIPTION:
    This routine parses the input date string according to the given format
    and extracts a [day,month,year] triplet if possible. The output is a list
    of three strings, with both day and month having length 2, and year having
    length 4. Example: ['01','02','2003']

    If the routine can't parse the date string successfully an empty
    list is returned.

    Valid four digit year values must be in the interval [1850,2005].
  """

  # Apply character replace table - - - - - - - - - - - - - - - - - - - - - - -
  #
  date_str = date_str.translate(replace_table)

  # Remove leading and trailing whitespaces and make lower case - - - - - - - -
  #
  format_str = format_str.strip()
  date_str   = date_str.strip().lower()

  # Replace triple and double spaces with one space only  - - - - - - - - - - -
  #
  date_str = date_str.replace('   ',' ')
  date_str = date_str.replace('  ',' ')
  date_str = date_str.replace('  ',' ')

  # Now check the date and format strings - - - - - - - - - - - - - - - - - - -
  #
  if (date_str =='') or (format_str ==''):
    return []  # No date or no format string given

  if (date_str =='') or (format_str ==''):
    return []  # No date or no format string given

  elif (' ' in format_str) and (' ' in date_str):
    date_list   = date_str.split()
    format_list = format_str.split()

  elif (' ' not in format_str) and (' ' not in date_str):
    date_list   = []
    format_list = []

    work_format_str = format_str
    work_date_str   = date_str
    while (work_format_str != ''):
      if (work_format_str[:2] == '%Y'):  # Four digit year
        format_list.append(work_format_str[:2])
        work_format_str = work_format_str[2:]
        date_list.append(work_date_str[:4])
        work_date_str = work_date_str[4:]

      elif (work_format_str[:2] in ['%m','%d','%y']):  # 2-digit year/month/day
        format_list.append(work_format_str[:2])
        work_format_str = work_format_str[2:]
        date_list.append(work_date_str[:2])
        work_date_str = work_date_str[2:]

      else:
        print 'error:Illegal format string without spaces: "%s"' % (formatstr)
        raise Exception

  else:  # A space in either date or format string but not in both
    return []  # Illegal combination of format string and date string

  day   = -1  # Set initially to illegal values
  month = -1
  year  = -1

  if (len(format_list) != 3) or (len(date_list) != 3):
    return []  # No valid format or date

  while (format_list != []):  # process all elements in the format list - - - -

    directive = format_list[0]

    if (directive in ['%b','%B']):  # Month abbreviation or name
      if (len(date_list[0]) >= 3):
        # Get month number (between 1..12) or -1 if not a valid month
        month = month_abbrev_dict.get(date_list[0][:3],-1)

    elif (directive == '%d'):  # Day of the month number
      try:
        day = int(date_list[0])  # Convert into an integer number
      except:
        day = -1  # Day is no integer

    elif (directive == '%m'):  # Month number (between 1..12)
      try:
        month = int(date_list[0])  # Convert into an integer number
      except:
        month = -1  # Month is no integer

    elif (directive == '%y'):  # Two digit year without century
      if (len(date_list[0]) == 2):  # Parse only if length is 2 digits
        try:
          year = int(date_list[0])  # Convert into an integer number
          if (year < 0) or (year > 99):
            year = -1  # Illegal year value in two digit
        except:
          year = -1  # year is no integer

        # Convert year into four digit value according to value of pivot_year
        #
        if (year >= 0) and (year < pivot_year):
          year = year+2000
        elif (year >= pivot_year) and (year < 100):
          year = year+1900

    elif (directive == '%Y'):  # Four digit year with century
      if (len(date_list[0]) == 4):  # Parse only if length is 4 digits
        try:
          year = int(date_list[0])  # Convert into an integer number
          if (year < 1850) or (year > 2005):
            year = -1  # Illegal year value in four digit
        except:
          year = -1  # year is no integer

    else:
      print 'error:Illegal format directive: "%s"' % (directive)
      raise Exception

    date_list =   date_list[1:]    # Remove processed first element
    format_list = format_list[1:]

  if (day == -1) or (month == -1) or (year == -1):  # No date parsed
    return []

  # Now do some test on the range of the values - - - - - - - - - - - - - - - -
  #
  if ((year % 4) == 0):
    if ((year % 100) == 0):
      if ((year % 400) == 0):
        leap_year = True
      else:
        leap_year = False
    else:
      leap_year = True
  else:
    leap_year = False

  valid = True

  if (month == 2):
    if (leap_year == True) and (day > 29):
      valid = False  # Illegal day number in February in a leap year
    elif (leap_year == False) and (day > 28):
      valid = False  # Illegal day number in February in a normal year
  elif (month in [1,3,5,7,8,10,12]) and (day > 31):
      valid = False  # Illegal day number in 31-day months
  elif (month in [4,6,9,11]) and (day > 30):
      valid = False  # Illegal day number in 30-day months

  if (valid == False):
    return []

  else:
    day_str =   string.zfill(str(day),2)
    month_str = string.zfill(str(month),2)
    year_str =  str(year)

    return [day_str,month_str,year_str]

# =============================================================================

def get_today():
  """Return today's date as a [day,month,year] tuple, with three string (both
     day and month having length 2, and year having length 4).
  """

  sys_time = time.localtime(time.time())  # Get current system time and date
  today = [string.zfill(str(sys_time[2]),2), \
           string.zfill(str(sys_time[1]),2), \
           str(sys_time[0])]

  # A log message for high volume log output (level 3)  - - - - - - - - - - - -
  #
  print '3:    Today is %s' % (str(today))

  return today

# =============================================================================
