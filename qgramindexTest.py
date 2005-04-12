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
# The Original Software is: "qgramindexTest.py"
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

"""Module qgramindexTest.py - Test module for qgramindex.py.
"""

# -----------------------------------------------------------------------------

import sets
import unittest

import qgramindex

# -----------------------------------------------------------------------------

class TestCase(unittest.TestCase):

  # Initialise test case  - - - - - - - - - - - - - - - - - - - - - - - - - - -
  #
  def setUp(self):

    # Initialise the PosQGramIndex class
    #
    self.test_index = qgramindex.PosQGramIndex(name = 'test index',
                                        description = 'test Q-Gram index',
                                         field_name = 'test_field',
                                    q_gram_len_list = [(1,(1,4)),(2,(4,99))],
                                      max_edit_dist = 2)

    # List of strings and their positional q-gram, for testing the methods
    # 'str_to_pos_q_grams' and 'pos_q_grams_to_str'.
    #
    self.q_gram_names = [('peter',1,[('p',0),('e',1),('t',2),('e',3),('r',4)]),
                         ('peter',2,[('pe',0),('et',1),('te',2),('er',3)]),
                         ('peter',3,[('pet',0),('ete',1),('ter',2)]),
                         ('peter',4,[('pete',0),('eter',1)]),
                         ('peter',5,[('peter',0)]),
                         ('peter',6,[]),
                         ('peter',7,[]),
                         ('peter',99,[]),
                         ('tim',1,[('t',0),('i',1),('m',2)]),
                         ('tim',2,[('ti',0),('im',1)]),
                         ('tim',3,[('tim',0)]),
                         ('tim',4,[]),
                         ('x',1,[('x',0)]),
                         ('x',2,[]),
                         ('',1,[]),
                         ('',2,[]),
                         ('0123456789',1,[('0',0),('1',1),('2',2),('3',3),
                          ('4',4),('5',5),('6',6),('7',7),('8',8),('9',9)]),
                         ('0123456789',2,[('01',0),('12',1),('23',2),('34',3),
                          ('45',4),('56',5),('67',6),('78',7),('89',8)]),
                         ('0123456789',3,[('012',0),('123',1),('234',2),
                          ('345',3),('456',4),('567',5),('678',6),('789',7)]),
                         ('0123456789',4,[('0123',0),('1234',1),('2345',2),
                          ('3456',3),('4567',4),('5678',5),('6789',6)]),
                         ('0123456789',5,[('01234',0),('12345',1),('23456',2),
                          ('34567',3),('45678',4),('56789',5)]),
                         ('0123456789',6,[('012345',0),('123456',1),
                          ('234567',2),('345678',3),('456789',4)]),
                         ('0123456789',7,[('0123456',0),('1234567',1),
                          ('2345678',2),('3456789',3)]),
                         ('0123456789',8,[('01234567',0),('12345678',1),
                          ('23456789',2)]),
                         ('0123456789',9,[('012345678',0),('123456789',1)]),
                         ('0123456789',10,[('0123456789',0)]),
                         ('0123456789',11,[]),
                       ]
    # List of sub-lists for testing the method 'get_sub_lists'
    #
    self.sub_list_tuples = [([1],0,[[]]),
                            ([1],1,[[1]]),
                            ([1,2],0,[[]]),
                            ([1,2],1,[[2],[1]]),
                            ([1,2],2,[[1,2]]),
                            ([1,2,3],0,[[]]),
                            ([1,2,3],1,[[3],[2],[1]]),
                            ([1,2,3],2,[[2,3],[1,3],[1,2]]),
                            ([1,2,3],3,[[1,2,3]]),
                            ([1,2,3,4],0,[[]]),
                            ([1,2,3,4],4,[[1,2,3,4]]),
                            ([1,2,3,4],3,[[2,3,4],[1,3,4],[1,2,4],[1,2,3]]),
                            ([1,2,3,4],2,[[3,4],[2,4],[2,3],[1,4],[1,3],
                             [1,2]]),
                            ([1,2,3,4],1,[[4],[3],[2],[1]]),
                            ([1,2,3,4,5],0,[[]]),
                            ([1,2,3,4,5],5,[[1,2,3,4,5]]),
                            ([1,2,3,4,5],4,[[2,3,4,5],[1,3,4,5],[1,2,4,5],
                             [1,2,3,5],[1,2,3,4]]),
                            ([1,2,3,4,5],3,[[3,4,5],[2,4,5],[2,3,5],[2,3,4],
                             [1,4,5],[1,3,5],[1,3,4],[1,2,5],[1,2,4],[1,2,3]]),
                            ([1,2,3,4,5],2,[[4,5],[3,5],[3,4],[2,5],[2,4],
                             [2,3],[1,5],[1,4],[1,3],[1,2]]),
                            ([1,2,3,4,5],1,[[5],[4],[3],[2],[1]]),
                           ]

    # Values for building an inverted index
    #
    self.test_input_list = ['a','ab','abc','abcd','abcde','abcdef']

    # The expected inverted index created for the above input values
    #
    start_chr = self.test_index.start_char
    end_chr =   self.test_index.end_char

    self.test_inv_index = {'a':{0:sets.Set([0,1,2,3])},
                           'b':{1:sets.Set([1,2,3])},
                           'c':{2:sets.Set([2,3])},
                           'd':{3:sets.Set([3])},
                 start_chr+'a':{0:sets.Set([3,4,5])},
                          'ab':{1:sets.Set([3,4,5])},
                          'bc':{2:sets.Set([3,4,5])},
                          'cd':{3:sets.Set([3,4,5])},
                          'de':{4:sets.Set([4,5])},
                          'ef':{5:sets.Set([5])},
                   'd'+end_chr:{4:sets.Set([3])},
                   'e'+end_chr:{5:sets.Set([4])},
                   'f'+end_chr:{6:sets.Set([5])}}

    # The expected matches for the above 'test_input_list' (with edit_dist = 2)
    #
    self.test_matches_list = [['a','ab','abc'],
                              ['a','ab','abc','abcd'],
                              ['a','ab','abc','abcd'],
                              ['ab','abc','abcd','abcde','abcdef'],
                              ['abcd','abcde','abcdef'],
                              ['abcd','abcde','abcdef']]

    # Test input q-gram list for intersection with inverted index
    #
    self.test_qgrams = [[(0,'a'),(2,'c')]]

    # And the corresponding matches
    #
    self.test_matches = ['a','abc']

    # List of values for select best match methods
    #
    self.best_match_list = [('peter',['petra','pete','petar','pita'],1),
                            ('gail',['gayle','kyle','gaill','grail'],2)
                           ]

  # Clean up test case  - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  #
  def tearDown(self):
    pass  # Nothing to clean up

  # ---------------------------------------------------------------------------
  #
  # Start test cases

  def test_str_to_pos_q_grams(self):  # - - - - - - - - - - - - - - - - - - - -
    """Test 'str_to_pos_q_grams' method"""

    for (in_str, q, q_gram_list) in self.q_gram_names:

      test_list = self.test_index.str_to_pos_q_grams(in_str,q)

      assert (test_list == q_gram_list), \
             '"str_to_pos_q_grams" returns wrong list for test tuple %s: %s' \
             % (str((in_str, q, q_gram_list)),str(test_list))

  def test_pos_q_grams_to_str(self):  # - - - - - - - - - - - - - - - - - - - -
    """Test 'pos_q_grams_to_str' method"""

    for (in_str, q, q_gram_list) in self.q_gram_names:

      if (q_gram_list != []):
        test_str = self.test_index.pos_q_grams_to_str(q_gram_list)

        assert (test_str == in_str), \
               '"pos_q_grams_to_str" returns wrong string for test tuple %s: '\
               % (str((in_str, q, q_gram_list))) + '"%s"' % (str(test_str))

  def test_get_sub_lists(self):  # - - - - - - - - - - - - - - - - - - - - - -
    """Test 'get_sub_lists' method"""

    for (in_list, l, in_sub_lists) in self.sub_list_tuples:

      test_sub_lists = self.test_index.get_sub_lists(in_list,l)

      assert (test_sub_lists == in_sub_lists), \
             '"get_sub_lists" returns wrong sub-lists for test tuple %s: %s'\
              % (str((in_list, l, in_sub_lists)), str(test_sub_lists))

  def test_build_index(self):  # - - - - - - - - - - - - - - - - - - - - - - -
    """Test 'build_index' method"""

    self.test_index.build_index(self.test_input_list)

    assert (len(self.test_index.value_list) == len(self.test_input_list)), \
           'Wrong length of value list after "build_index": %d, should be %d' \
           % (len(self.test_index.value_list), len(self.test_input_set))

    assert (self.test_index.inv_index == self.test_inv_index), \
           'Wrong inverted index after "build_index": %s\nShould be: %s' % \
           (str(self.test_index.inv_index), str(self.test_inv_index))

  def test_get_matches(self): # - - - - - - - - - - - - - - - - - - - - - - - -
    """Test 'get_matches' method"""

    self.test_index.build_index(self.test_input_list)

    for i in range(len(self.test_input_list)):
      in_val = self.test_input_list[i]
      res_list = self.test_matches_list[i]  # Expected results

      test_list = self.test_index.get_matches(in_val)
      test_list.sort()

      assert (test_list == res_list), 'Wrong match list retrived from ' + \
             'index (in value "%s"): %s / %s' % \
             (in_val, str(res_list), str(test_list))

    # Some more tests
    #
    test_list = self.test_index.get_matches('b')
    assert (test_list == ['ab','abc']), 'Wrong match list retrived from ' + \
             'index (in value "b"): %s / ["ab","abc"]' % (str(test_list))

    test_list = self.test_index.get_matches('c')
    assert (test_list == ['abc']), 'Wrong match list retrived from ' + \
             'index (in value "c"): %s / ["abc"]' % (str(test_list))

    test_list = self.test_index.get_matches('d')
    assert (test_list == []), 'Wrong match list retrived from ' + \
             'index (in value "d"): %s / []' % (str(test_list))

  def test_select_best_match(self): # - - - - - - - - - - - - - - - - - - - - -
    """Test 'select_best_match' method"""

    for test_tuple in self.best_match_list:
      best_test = self.test_index.select_best_match(test_tuple[0],
                                                    test_tuple[1],
                                                    'winkler')

      assert (best_test[0][0] == test_tuple[1][test_tuple[2]]), \
             'Wrong best match returned: %s / %s' % (best_test[0][0],
             test_tuple[1][test_tuple[2]])

#  def pos_q_gram_set_intersect(self):  # - - - - - - - - - - - - - - - - - - -
#    """Test 'pos_q_gram_set_intersect' method"""
#
#    for (in_dict, max_ed, intersect_set) in self.intersect_lists:
#
#      test_set = self.test_index.pos_q_gram_set_intersect(in_dict, \
#                 self.inv_indices)

# -----------------------------------------------------------------------------
# Start tests when called from command line

if (__name__ == "__main__"):
  unittest.main()  # Run all test

  # The following code does the same as 'unittest.main()'
  #
  # mysuite = unittest.makeSuite(TestCase,'test')
  # testrunner = unittest.TextTestRunner(verbosity=1)
  # testrunner.run(mysuite)
