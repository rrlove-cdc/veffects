#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 16 15:48:51 2018

@author: beccalove
"""

import unittest
import requests
import veffects

class TestPOSTRequest(unittest.TestCase):
    
    def setUp(self):
        
        gene = "AGAP004687"
        
        self.basic_post = veffects.make_POST_request(gene)
        
    def test_basic(self):
        
        returned_id = "AGAP004687-RA"
        
        returned_molecule = 'dna'
        
        returned_seq = "ATGGTACACAACGGCATCGACTACGGTGATATGCAGCTGA" +\
        "TATGCGAAGCGTGTCACCTGATGCTGGCGCTGGGCATGACACGGAAGGAAATGGT" +\
        "ACAGGAGTTCGACGTGTGGAACAAGGGCGTATTGGATTCGTTCCTGATCGAGATC" +\
        "CCGCACGATTTTCTCAACCAACGCGACGTTGAGGGATAG" 

        self.assertEqual(len(self.basic_post), 1)
        
        self.assertEqual(self.basic_post[0]["id"], returned_id)
        
        self.assertEqual(self.basic_post[0]["molecule"], returned_molecule)
        
        self.assertEqual(self.basic_post[0]["seq"], returned_seq)

    def test_bad_gene_name(self):

        self.assertRaises(veffects.methods.RequestReturnError, 
                          veffects.make_POST_request, "foo")

    def test_bad_server(self):
        
        self.assertRaises(requests.exceptions.MissingSchema, 
                          veffects.make_POST_request, "bar", 
                          post_server = "foo")

class TestGETRequest(unittest.TestCase):
    
    def setUp(self):
        
        gene = "AGAP004687"
        
        self.basic_get = veffects.make_GET_request(gene)
    
    def test_basic(self):
        
        returned_start = 819113
        returned_end = 819301
        
        self.assertEqual(len(self.basic_get), 2)
        
        self.assertEqual(self.basic_get[0]["start"], returned_start)
        
        self.assertEqual(self.basic_get[1]["start"], returned_start)
        
        self.assertEqual(self.basic_get[0]["end"], returned_end)
        
        self.assertEqual(self.basic_get[1]["end"], returned_end)
        
    def test_bad_name(self):
        
        self.assertRaises(veffects.methods.RequestReturnError,
                          veffects.make_GET_request, "foo")
        
    def test_bad_server(self):
        
        self.assertRaises(requests.exceptions.MissingSchema,
                          veffects.make_GET_request, "foo",
                          get_server="bar")
        
class TestMakeTranscript(unittest.TestCase):
    
    def setUp(self):
        pass
    
class TestMakeExons(unittest.TestCase):
    
    def setUp(self):
        pass
    
class TestWorkflow(unittest.TestCase):
    
    def setUp(self):
        pass

if __name__ == '__main__':
	unittest.main()
