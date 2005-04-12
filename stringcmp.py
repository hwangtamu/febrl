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
# The Original Software is: "stringcmp.py"
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

"""Module stringcmp.py - Several approximate string comparison routines.

   Provides routines for several approximate string comparisons. All return
   a value between 1.0 (the strings are the same) and 0.0 (strings are totally
   different).

   ROUTINES
     jaro        Jaro
     winkler     Winkler (based on Jaro)
     bigram      Bigram based
     editdist    Edit-distance (or Levenshtein distance)
     bagdist     Bag distance (cheap distance based method)
     seqmatch    Uses Python's standard library 'difflib'
     compression Based on Zlib compression algorithm
     permwinkler Winkler combined with permutations of words, improves results
                 for swapped words
     sortwinkler Winkler with sorted words (if more than one),  improves
                 results for swapped words

   See doc strings of individual functions for detailed documentation.

   If called from command line, a test routine is run which prints example
   approximate string comparisons for various string pairs.
"""

# =============================================================================
# Imports go here

import difflib
import logging
import zlib

# =============================================================================
# Special character used in the Jaro and Winkler comparions functions.
# Thanks to Luca Montecchiani (luca.mon@aliceposta.it).
#
jaro_winkler_marker_char = chr(1)

# =============================================================================

def do_stringcmp(cmp_method, str1, str2):
  """A 'chooser' functions which performs the selected approximate string
     comparison method.
  """

  if (cmp_method == 'jaro'):
    match_score = jaro(str1, str2)
  elif (cmp_method == 'winkler'):
    match_score = winkler(str1, str2)
  elif (cmp_method == 'bigram'):
    match_score = bigram(str1, str2)
  elif (cmp_method == 'editdist'):
    match_score = editdist(str1, str2)
  elif (cmp_method == 'bagdist'):
    match_score = bagdist(str1, str2)
  elif (cmp_method == 'seqmatch'):
    match_score = seqmatch(str1, str2)
  elif (cmp_method == 'compression'):
    match_score = compression(str1, str2)
  elif (cmp_method == 'sortwinkler'):
    match_score = sortwinkler(str1, str2)
  elif (cmp_method == 'permwinkler'):
    match_score = permwinkler(str1, str2)
  else:
    logging.exception('Illegal approximate string comparison method: %s' \
                      % (str_cmp_method))
    raise Exception

  return match_score

# =============================================================================

# A function to create permutations of a list (from ASPN Python cookbook, see:
# http://aspn.activestate.com/ASPN/Cookbook/Python/Recipe/66463)

def getPermutations(a):
   if len(a)==1:
      yield a
   else:
      for i in range(len(a)):
         this = a[i]
         rest = a[:i] + a[i+1:]
         for p in getPermutations(rest):
            yield [this] + p

def permute(alist):
  reslist = []
  for l in getPermutations(alist):
    reslist.append(' '.join(l))

  return reslist

# =============================================================================

def jaro(str1, str2):
  """Return approximate string comparator measure (between 0.0 and 1.0)

  USAGE:
    score = jaro(str1, str2)

  ARGUMENTS:
    str1  The first string
    str2  The second string

  DESCRIPTION:
    As desribed in 'An Application of the Fellegi-Sunter Model of
    Record Linkage to the 1990 U.S. Decennial Census' by William E. Winkler
    and Yves Thibaudeau.
  """

  # Quick check if the strings are the same - - - - - - - - - - - - - - - - - -
  #
  if (str1 == str2):
    return 1.0

  len1 = len(str1)
  len2 = len(str2)
  halflen = max(len1, len2) / 2 + 1

  ass1 = ''  # Characters assigned in str1
  ass2 = ''  # Characters assigned in str2

  workstr1 = str1  # Copy of original string
  workstr2 = str2

  common1 = 0  # Number of common characters
  common2 = 0

  # Analyse the first string  - - - - - - - - - - - - - - - - - - - - - - - - -
  #
  for i in range(len1):
    start = max(0,i-halflen)
    end   = min(i+halflen+1,len2)
    index = workstr2.find(str1[i],start,end)
    if (index > -1):  # Found common character
      common1 += 1
      ass1 = ass1+str1[i]
      workstr2 = workstr2[:index]+jaro_winkler_marker_char+workstr2[index+1:]

  # Analyse the second string - - - - - - - - - - - - - - - - - - - - - - - - -
  #
  for i in range(len2):
    start = max(0,i-halflen)
    end   = min(i+halflen+1,len1)
    index = workstr1.find(str2[i],start,end)
    if (index > -1):  # Found common character
      common2 += 1
      ass2 = ass2 + str2[i]
      workstr1 = workstr1[:index]+jaro_winkler_marker_char+workstr1[index+1:]

  if (common1 != common2):
    logging.error('Jaro: Wrong common values for strings "%s" and "%s"' % \
          (str1, str2) + ', common1: %i, common2: %i' % (common1, common2) + \
          ', common should be the same.')
    common1 = float(common1+common2) / 2.0  ##### This is just a fix #####

  if (common1 == 0):
    return 0.0

  # Compute number of transpositions  - - - - - - - - - - - - - - - - - - - - -
  #
  transposition = 0
  for i in range(len(ass1)):
    if (ass1[i] != ass2[i]):
      transposition += 1
  transposition = transposition / 2.0

  common1 = float(common1)
  w = 1./3.*(common1 / float(len1) + common1 / float(len2) + \
           (common1-transposition) / common1)

  # A log message - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  #
  logging.debug('Jaro comparator string "%s" with "%s" value: %.3f' % \
                (str1, str2, w))
  return w

# =============================================================================

def winkler(str1, str2):
  """Return approximate string comparator measure (between 0.0 and 1.0)

  USAGE:
    score = winkler(str1, str2)

  ARGUMENTS:
    str1  The first string
    str2  The second string

  DESCRIPTION:
    As desribed in 'An Application of the Fellegi-Sunter Model of
    Record Linkage to the 1990 U.S. Decennial Census' by William E. Winkler
    and Yves Thibaudeau.

    Based on the 'jaro' string comparator, but modifies it according to wether
    the first few characters are the same or not.
  """

  # Quick check if the strings are the same - - - - - - - - - - - - - - - - - -
  #
  if (str1 == str2):
    return 1.0

  len1 = len(str1)
  len2 = len(str2)
  halflen = max(len1,len2) / 2 + 1

  ass1 = ''  # Characters assigned in str1
  ass2 = ''  # Characters assigned in str2
  workstr1 = str1
  workstr2 = str2

  common1 = 0  # Number of common characters
  common2 = 0

  # Analyse the first string  - - - - - - - - - - - - - - - - - - - - - - - - -
  #
  for i in range(len1):
    start = max(0,i-halflen)
    end   = min(i+halflen+1,len2)
    index = workstr2.find(str1[i],start,end)
    if (index > -1):  # Found common character
      common1 += 1
      ass1 = ass1 + str1[i]
      workstr2 = workstr2[:index]+jaro_winkler_marker_char+workstr2[index+1:]

  # Analyse the second string - - - - - - - - - - - - - - - - - - - - - - - - -
  #
  for i in range(len2):
    start = max(0,i-halflen)
    end   = min(i+halflen+1,len1)
    index = workstr1.find(str2[i],start,end)
    if (index > -1):  # Found common character
      common2 += 1
      ass2 = ass2 + str2[i]
      workstr1 = workstr1[:index]+jaro_winkler_marker_char+workstr1[index+1:]

  if (common1 != common2):
    logging.error('Winkler: Wrong common values for strings "%s" and "%s"' % \
          (str1, str2) + ', common1: %i, common2: %i' % (common1, common2) + \
          ', common should be the same.')
    common1 = float(common1+common2) / 2.0  ##### This is just a fix #####

  if (common1 == 0):
    return 0.0

  # Compute number of transpositions  - - - - - - - - - - - - - - - - - - - - -
  #
  transposition = 0
  for i in range(len(ass1)):
    if (ass1[i] != ass2[i]):
      transposition += 1
  transposition = transposition / 2.0

  # Now compute how many characters are common at beginning - - - - - - - - - -
  #
  minlen = min(len1,len2)
  for same in range(minlen+1):
    if (str1[:same] != str2[:same]):
      break
  same -= 1
  if (same > 4):
    same = 4

  common1 = float(common1)
  w = 1./3.*(common1 / float(len1) + common1 / float(len2) + \
           (common1-transposition) / common1)

  wn = w + same*0.1 * (1.0 - w)

  # A log message - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  #
  logging.debug('Winkler comparator string "%s" with "%s" value: %.3f' % \
                (str1, str2, w))
  return wn

# =============================================================================

def bigram(str1, str2):
  """Return approximate string comparator measure (between 0.0 and 1.0)
     using bigrams.

  USAGE:
    score = bigram(str1, str2)

  ARGUMENTS:
    str1  The first string
    str2  The second string

  DESCRIPTION:
    Bigrams are two-character sub-strings contained in a string. For example,
    'peter' contains the bigrams: pe,et,te,er.

    This routine counts the number of common bigrams and divides by the
    average number of bigrams. The resulting number is returned.
  """

  # Quick check if the strings are the same - - - - - - - - - - - - - - - - - -
  #
  if (str1 == str2):
    return 1.0

  bigr1 = []
  bigr2 = []

  # Make a list of bigrams for both strings - - - - - - - - - - - - - - - - - -
  #
  for i in range(1,len(str1)):
    bigr1.append(str1[i-1:i+1])

  for i in range(1,len(str2)):
    bigr2.append(str2[i-1:i+1])

  # Compute average number of bigrams - - - - - - - - - - - - - - - - - - - - -
  #
  average = (len(bigr1)+len(bigr2)) / 2.0
  if (average == 0.0):
    return 0.0

  # Get common bigrams  - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  #
  common = 0.0

  if (len(bigr1) < len(bigr2)):  # Count using the shorter bigram list
    short_bigr = bigr1
    long_bigr  = bigr2
  else:
    short_bigr = bigr2
    long_bigr  = bigr1

  for b in short_bigr:
    if (b in long_bigr):
      common += 1.0
      long_bigr[long_bigr.index(b)] = []  # Mark this bigram as counted

  w = common / average

  # A log message - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  #
  logging.debug('Bigram comparator string "%s" with "%s" value: %.3f' % \
                (str1, str2, w))
  return w

# =============================================================================

def editdist(str1, str2):
  """Return approximate string comparator measure (between 0.0 and 1.0)
     using the edit (or Levenshtein) distance.

  USAGE:
    score = editdist(str1, str2)

  ARGUMENTS:
    str1  The first string
    str2  The second string

  DESCRIPTION:
    The edit distance is the minimal number of insertions, deletions and
    substitutions needed to make two strings equal.

    For more information on the modified Soundex see:
    - http://www.nist.gov/dads/HTML/editdistance.html
  """

  # Quick check if the strings are the same - - - - - - - - - - - - - - - - - -
  #
  if (str1 == str2):
    return 1.0

  n = len(str1)
  m = len(str2)

  if (n == 0) or (m == 0):  # Check if strings are of length zero
    return 0.0

  d = range(n+1)  # Create matrix
  for i in range(n+1):
    d[i] = range(m+1)
    for j in range(m+1):
      d[i][j] = 0

  for i in range(n+1):  # Set initial values
    d[i][0] = i
  for j in range(m+1):
    d[0][j] = j

  for i in range(1,n+1):
    s = str1[i-1]

    for j in range(1,m+1):
      t = str2[j-1]
      if (s == t):
        cost = 0
      else:
        cost = 1

      d[i][j] = min(d[i-1][j]+1, d[i][j-1]+1,d[i-1][j-1]+cost)

  w = float(max(n,m) - d[n][m]) / float(max(n,m))

  # A log message - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  #
  logging.debug('Edit-distance comparator string "%s" with "%s" value: %.3f' \
                % (str1, str2, w))
  return w

# =============================================================================

def bagdist(str1, str2):
  """Return approximate string comparator measure (between 0.0 and 1.0)
     using the bag distance.

  USAGE:
    score = bagdist(str1, str2)

  ARGUMENTS:
    str1  The first string
    str2  The second string

  DESCRIPTION:
    Bag distance is a cheap method to calculate the distance between two
    strings. It is always smaller or equal to the edit distance, and therefore
    the similarity measure returned by the method is always larger than the
    edit distance similarity measure.

    For more details see for example:

      "String Matching with Metric Trees Using an Approximate Distance"
      Ilaria Bartolini, Paolo Ciaccia and Marco Patella,
      in Proceedings of the 9th International Symposium on String Processing
      and Information Retrieval, Lisbone, Purtugal, September 2002.
  """

  # Quick check if the strings are the same - - - - - - - - - - - - - - - - - -
  #
  if (str1 == str2):
    return 1.0

  n = len(str1)
  m = len(str2)

  if (n == 0) or (m == 0):  # Check if strings are of length zero
    return 0.0

  list1 = list(str1)
  list2 = list(str2)

  tmplist1 = list1[:]
  tmplist2 = list2[:]

  for ch in list2:
    if (ch in tmplist1):
      tmplist1.remove(ch)

  for ch in list1:
    if (ch in tmplist2):
      tmplist2.remove(ch)

  b = max(len(tmplist1),len(tmplist2))

  w = float(max(n,m) - b) / float(max(n,m))

  # A log message - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  #
  logging.debug('Bag-distance comparator string "%s" with "%s" value: %.3f' % \
                (str1, str2, w))
  return w

# =============================================================================

def seqmatch(str1, str2):
  """Return approximate string comparator measure (between 0.0 and 1.0)
     using the Python standard library 'difflib' sequence matcher.

     Because the matches are not commutative, the pair and the swapped pair are
     compared and the average is taken.

  USAGE:
    score = seqmatch(str1, str2)

  ARGUMENTS:
    str1  The first string
    str2  The second string

  DESCRIPTION:
    For more information on Python's 'difflib' library see:

      http://www.python.org/doc/current/lib/module-difflib.html
  """

  # Quick check if the strings are the same - - - - - - - - - - - - - - - - - -
  #
  if (str1 == str2):
    return 1.0

  seq_matcher_1 = difflib.SequenceMatcher(None, str1, str2)
  seq_matcher_2 = difflib.SequenceMatcher(None, str2, str1)

  w = (seq_matcher_1.ratio()+seq_matcher_2.ratio()) / 2.0  # Return average

  # A log message - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  #
  logging.debug('Seq-match comparator string "%s" with "%s" value: %.3f' % \
                (str1, str2, w))
  return w

# =============================================================================

def compression(str1, str2):
  """Return approximate string comparator measure (between 0.0 and 1.0)
     using the zlib compression library.

  USAGE:
    score = compression(str1, str2)

  ARGUMENTS:
    str1  The first string
    str2  The second string

  DESCRIPTION:
    For more information about using compression for similarity measures see:

    - Cilibrasi, R. and Vitanyi, P.: Clustering by compression. IEEE Trans.
      Infomat. Th. Submitted, 2004. See: http://arxiv.org/abs/cs.CV/0312044

    - Keogh, E., Lonardi, S. and Ratanamahatana, C.A.: Towards parameter-free
      data mining. Proceedings of the 2004 ACM SIGKDD international conference
      on Knowledge discovery and data mining, pp. 206-215, Seattle, 2004.
  """

  # Quick check if the strings are the same - - - - - - - - - - - - - - - - - -
  #
  if (str1 == str2):
    return 1.0

  c1 =  float(len(zlib.compress(str1)))
  c2 =  float(len(zlib.compress(str2)))
  c12 = 0.5 * (len(zlib.compress(str1+str2)) + len(zlib.compress(str2+str1)))

  if (c12 == 0.0):
    return 0.0  # Maximal distance

  w = 1.0 - (c12 - min(c1,c2)) / max(c1,c2)

  if (w < 0.0):
    print 'warning:Compression based comparison smaller than 0.0 with ' + \
          'strings "%s" and "%s": %.3f (cap to 1.0)' % (str1, str2, w)
    w = 0.0

  # A log message - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  #
  logging.debug('Compression comparator string "%s" with "%s" value: %.3f' % \
                (str1, str2, w))
  return w

# =============================================================================

def permwinkler(str1, str2):
  """Return approximate string comparator measure (between 0.0 and 1.0) using
     a combination of the Winkler string comparator on all permutations of
     words (ifd there are more than one in the input strings), which improves
     the results for swapped words.

  USAGE:
    score = permwinkler(str1, str2)

  ARGUMENTS:
    str1  The first string
    str2  The second string

  DESCRIPTION:
    If one or both of the input strings contain more than one words all
    possible permutations of are compared using the Winkler approximate string
    comparator, and the maximum value is returned.

    If both input strings contain one word only then the standard Winkler
    string comparator is used.
  """

  if (' ' not in str1) and (' ' not in str2):
    w = winkler(str1, str2)  # Standard Winkler

  else:  # At least one of the strings contains two words

    str_list1 = str1.split(' ')
    str_list2 = str2.split(' ')

    perm_list1 = permute(str_list1)
    perm_list2 = permute(str_list2)

    w =        -1.0  # Maximal similarity measure
    max_perm = None

    for perm1 in perm_list1:
      for perm2 in perm_list2:

        # Calculate standard winkler for this permutation
        #
        this_w = winkler(perm1, perm2)

        if (this_w > w):
          w        = this_w
          max_perm = [perm1, perm2]

    logging.debug('Permutation Winkler best permutation: %s' % (str(max_perm)))

  # A log message - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  #
  logging.debug('Permutation Winkler comparator string "%s" with "%s" value:' \
                % (str1, str2) + ' %.3f' % (w))
  return w

# =============================================================================

def sortwinkler(str1, str2):
  """Return approximate string comparator measure (between 0.0 and 1.0) using
     the Winkler string comparator on the word-sorted input strings (if there
     are more than one in the input strings), which improves the results for
     swapped words.

  USAGE:
    score = sortwinkler(str1, str2)

  ARGUMENTS:
    str1  The first string
    str2  The second string

  DESCRIPTION:
    If one or both of the input strings contain more than one words then the
    input string is word-sorted before the standard Winkler approximate string
    comparator is applied.

    If both input strings contain one word only then the standard Winkler
    string comparator is used.
  """

  if (' ' in str1):  # Sort string 1
    word_list = str1.split(' ')
    word_list.sort()
    str1 = ' '.join(word_list)

  if (' ' in str2):  # Sort string 2
    word_list = str2.split(' ')
    word_list.sort()
    str2 = ' '.join(word_list)

  w = winkler(str1, str2)  # Standard Winkler

  # A log message - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  #
  logging.debug('Sorted Winkler comparator string "%s" with "%s" value: %.3f' \
                % (str1, str2, w))
  return w

# =============================================================================
#
# Do some tests if called from command line
#
# Most test strings are taken from:
#   Approximate String Comparison and its Effect on an Advanced Record
#   Linkage System, Edward H. Porter and William W. Winkler, Bureau of
#   Census, 1997. Research report RR97/02.
#

if (__name__ == '__main__'):

  msg = []

  msg.append('Febrl module "stringcmp.py"')
  msg.append('---------------------------')
  msg.append('')

  strings = [['shackleford','dunningham','nichleson','jones','massey', \
              'abroms','hardin','itman','jeraldine','marhta','michelle', \
              'julies','tanya','dwayne','sean','jon','jon','brookhaven', \
              'brook hallow','decatur','fitzrureiter','higbee','higbee', \
              'lacura','iowa','1st','peter','abcde','yz','cunningham', \
              'campell','galloway','frederick','michele','jesse', \
              'jonathon','julies','yvette','dickson','dixon','peter', \
              'gondiwindi','delfini*','ein#','do','doe', \
              'louise marie', 'maria louisa', 'mighty joe', 'kim zhu', \
              'lim zhau kim'], \
             ['shackelford','cunnigham','nichulson','johnson','massie', \
              'abrams','martinez','smith','geraldine','martha','michael', \
              'julius','tonya','duane','susan','john','jan','brrokhaven', \
              'brook hllw','decatir','fitzenreiter','highee','higvee', \
              'locura','iona','ist','peter','fghij','abcdef', \
              'cunnigham','campbell','calloway','fredrick','michelle', \
              'jessie','jonathan','juluis','yevett','dixon','dickson', \
              'ole', 'gondiwindiro','delfini','eni','od','deo', \
              'marie lousie', 'louisa marie', 'joe mighty', 'zhou kim', \
              'kim lim zhao']]

  msg.append('     String 1      String 2   Jaro    Winkler'+ \
             ' Bigram  EditD   SeqMatch Compress BagD   PermWink  SortWink')

  for i in range(len(strings[0])):
    str1 = strings[0][i]
    str2 = strings[1][i]

    s = '%13s %13s   %.3f   %.3f' %(str1, str2, jaro(str1, str2), \
        winkler(str1, str2)) + \
        '   %.3f   %.3f   %.3f    %.3f    %.3f  %.3f     %.3f' % \
        (bigram(str1, str2), editdist(str1, str2), seqmatch(str1, str2), \
         compression(str1, str2),bagdist(str1,str2), permwinkler(str1,str2), \
         sortwinkler(str1,str2))
    msg.append(s)
  for m in msg:
    print m

# =============================================================================
