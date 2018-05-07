#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May  7 16:13:24 2018

@author: beccalove
"""

import unittest
import veffects

class ExonSequenceTestCase(unittest.TestCase):
    
    def setUp(self):
        
        self.exon = veffects.ExonSequence("2L", "foo", 3, "+", 109, 134)
        
    def test_length(self):
        
        self.assertEqual(self.exon.length, 26)
        
    def test_add_sequence(self):
        
        test_seq = "CGATATGAATATGACCAGATATGAGT"
        
        self.exon.add_sequence(test_seq)
        
        self.assertEqual(self.exon.sequence, test_seq)
        
    def test_sequence_is_right_length(self):
        
        test_seq = "A"

        self.assertRaises(ValueError, self.exon.add_sequence, test_seq)
    
    def test_change(self):
        
        pass
    
    def test_change_edge_case(self):
        
        pass

       