#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May  3 13:46:33 2018

@author: beccalove
"""

import unittest
import veffects

class TranscriptTestCase(unittest.TestCase):
    
    def setUp(self):
        
        test_seq = "ATGCTTCATCAGCAGCAGCCCGGGGCCTAG"
        self.transcript = veffects.Transcript("2L", test_seq)
        
        exon1 = veffects.ExonSequence("2L", "foo", 1, "+", 5, 18)#14
        exon2 = veffects.ExonSequence("2L", "foo", 2, "+", 50, 55)#6
        exon3 = veffects.ExonSequence("2L", "foo", 3, "+", 70, 74)#5
        exon4 = veffects.ExonSequence("2L", "foo", 4, "+", 103, 107)#5
        
        self.exons = [exon1, exon2, exon3, exon4]
        
        variant1 = veffects.VariantRecord("2L", 15, "A", "AA", 15)
        variant2 = veffects.VariantRecord("2L", 16, "GCA", "G", 18)
        variant3 = veffects.VariantRecord("2L", 52, "A", "AAAAA", 52)
        variant4 = veffects.VariantRecord("2L", 70, "CGG", "C", 72)
        
        self.variants = [variant1, variant2, variant3, variant4]
        
    def test_translate_seq(self):
        
        ##test the translated versions
        
        self.transcript.seq_changed = "ATGCGTCATCAGCAGCAGCAGTCCGGGGCCTAG"
        
        unchanged_translation = "MLHQQQPGA*"
        
        changed_translation = "MRHQQQQSGA*"
        
        self.transcript.translate_seqs()
        
        self.assertEqual(self.transcript.seq_translated, 
                         unchanged_translation)
        
        self.assertEqual(self.transcript.seq_changed_translated,
                         changed_translation)
        
    def test_populate_exon_seq_length_check(self):
        ##pass a set of exons whose collective length is less than the length
        ##of the sequence, and make sure an AssertionError is raised
        
        self.transcript.add_exon(self.exons[0])
        self.transcript.add_exon(self.exons[1])
        
        self.assertRaises(AssertionError, self.transcript.populate_exon_seq)
                
    def test_populate_exon_seq_correct_seqs(self):
        
        #ATGCTTCATCAGCAGCAGCCCGGGGCCTAG
        
        for exon in self.exons:
            self.transcript.add_exon(exon)
        
        self.transcript.populate_exon_seq()
        
        self.assertEqual(self.transcript.exons[0].sequence,
                         "ATGCTTCATCAGCA")
        
        self.assertEqual(self.transcript.exons[1].sequence,
                         "GCAGCC")
        
        self.assertEqual(self.transcript.exons[2].sequence,
                         "CGGGG")
                
        self.assertEqual(self.transcript.exons[3].sequence,
                         "CCTAG")
        
    def test_populate_exon_has_exons(self):
        
        self.assertRaises(AssertionError, self.transcript.populate_exon_seq)
        
    def test_catches_overlapping_variants(self):
        ##pass a set of variants that overlaps an exon boundary
        ##and make sure a VariantOverlapsExonBoundaryError is raised
        
        for exon in self.exons:
            
            self.transcript.add_exon(exon)

        variant5 = veffects.VariantRecord("2L", 4, "CAT", "C", 6)
        
        self.assertRaises(veffects.VariantCrossesExonBoundaryError,
                          self.transcript.check_for_overlapping_variants,
                          [variant5])
                
    def test_parse_variants_list(self):
        
        for exon in self.exons:
            self.transcript.add_exon(exon)

        self.transcript.parse_variants_list(self.variants)
        
        exon1_variants = [self.variants[0], self.variants[1]]
        exon2_variants = [self.variants[2]]
        exon3_variants = [self.variants[3]]
        
        self.assertEqual(self.transcript.exons[0].variants, exon1_variants)
        self.assertEqual(self.transcript.exons[1].variants, exon2_variants)
        self.assertEqual(self.transcript.exons[2].variants, exon3_variants)
        self.assertEqual(self.transcript.exons[3].variants, [])
    
    def test_assemble_changed_seq(self):
        
        for exon in self.exons:
            self.transcript.add_exon(exon)
                    
        self.transcript.exons[0].changed_sequence = "ATGCTTCATCAAGCACA"
        self.transcript.exons[1].changed_sequence = "GCAAAAAGCC"
        self.transcript.exons[2].changed_sequence = "CGG"
        self.transcript.exons[3].changed_sequence = "CCTAG"
        
        self.transcript.assemble_changed_seq()
        
        test_seq = "ATGCTTCATCAAGCACAGCAAAAAGCCCGGCCTAG"
        
        self.assertEqual(self.transcript.seq_changed, test_seq)
        
        