#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 16 15:48:51 2018

@author: beccalove
"""

import allel
import unittest
import requests
import veffects
from veffects import VariantRecord as vr

class TestPOSTRequest(unittest.TestCase):
    
    def setUp(self):
        
        gene = "AGAP004687-RA"
        
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

        self.assertRaises(veffects.RequestReturnError, 
                          veffects.make_POST_request, "foo-RZ")
        
    def test_gene_name_not_transcript_name(self):
        
        self.assertRaises(veffects.BadNameError, 
                          veffects.make_POST_request,
                          "AGAP004687")

    def test_bad_server(self):
        
        self.assertRaises(requests.exceptions.MissingSchema, 
                          veffects.make_POST_request, "foo-RZ", 
                          post_server = "bar")

class TestGETRequest(unittest.TestCase):
    
    def setUp(self):
        
        gene = "AGAP004687-RA"
        
        self.basic_get = veffects.make_GET_request(gene)
    
    def test_basic(self):
        
        returned_start = 819113
        returned_end = 819301
        
        self.assertEqual(len(self.basic_get), 2)
        
        self.assertEqual(self.basic_get[0]["start"], returned_start)
        
        self.assertEqual(self.basic_get[1]["start"], returned_start)
        
        self.assertEqual(self.basic_get[0]["end"], returned_end)
        
        self.assertEqual(self.basic_get[1]["end"], returned_end)
        
    def test_gene_name_not_transcript_name(self):
        
        self.assertRaises(veffects.BadNameError, veffects.make_GET_request, 
                          "AGAP004687")
        
    def test_bad_name(self):
        
        self.assertRaises(veffects.RequestReturnError,
                          veffects.make_GET_request, "foo-RZ")
        
    def test_bad_server(self):
        
        self.assertRaises(requests.exceptions.MissingSchema,
                          veffects.make_GET_request, "foo-RZ",
                          get_server="bar")
        
class TestMakeTranscript(unittest.TestCase):
    
    def setUp(self):
        
        self.feature_json = [{'desc': None,
                              'id': 'AGAP004687-RA',
                              'molecule': 'dna',
                              'query': 'AGAP004687',
                              'seq': 'ATGGTACACAACGGCATCGACTACGGTGATATGCA' +\
                              'GCTGATATGCGAAGCGTGTCACCTGATGCTGGCGCTGGGCAT' +\
                              'GACACGGAAGGAAATGGTACAGGAGTTCGACGTGTGGAACAA' +\
                              'GGGCGTATTGGATTCGTTCCTGATCGAGATCCCGCACGATTT' +\
                              'TCTCAACCAACGCGACGTTGAGGGATAG'}]
    
    def test_basic(self):
        
        self.transcript = veffects.make_transcript(self.feature_json)
        
        self.assertEqual(self.transcript.name, "AGAP004687-RA")
        
        test_seq = 'ATGGTACACAACGGCATCGACTACGGTGATATGCA' +\
        'GCTGATATGCGAAGCGTGTCACCTGATGCTGGCGCTGGGCAT' +\
        'GACACGGAAGGAAATGGTACAGGAGTTCGACGTGTGGAACAA' +\
        'GGGCGTATTGGATTCGTTCCTGATCGAGATCCCGCACGATTT' +\
        'TCTCAACCAACGCGACGTTGAGGGATAG'
    
        self.assertEqual(self.transcript.seq, test_seq)
        
class TestMakeExons(unittest.TestCase):
    
    def setUp(self):
                
        feature_table_dtype = [
                ('seqid','O'),
                ('source','O'),
                ('type','O'),
                ('start','<i8'),
                ('end','<i8'),
                ('score','<f8'),
                ('strand','O'),
                ('phase','<i8'),
                ('Parent','O'),
                ('Name','O'),
                ('ID','O'),
                ]
        
        feature_table_data = [
                ('2L', 'DB', 'mRNA', 819113, 819301, -1, '-', -1,
                 'AGAP004687', '.', 'AGAP004687-RA'),
                ('2L', 'DB', 'exon', 819113, 819301, -1, '-', -1,
                 'AGAP004687-RA', '.', 'E013854A'),
                ('2L', 'DB', 'CDS', 819113, 819301, -1, '-', 0,
                 'AGAP004687-RA', '.', '.')
                ]
            
        self.feature_table = \
            allel.FeatureTable(feature_table_data, dtype=feature_table_dtype)
        
    def test_basic(self):
        
        self.exons = veffects.make_exons(self.feature_table,
                                         "AGAP004687-RA")
        
        self.assertEqual(len(self.exons), 1)
        
        self.assertEqual(self.exons[0].chrom, "2L")
        
        self.assertEqual(self.exons[0].transcript, "AGAP004687-RA")
        
        self.assertEqual(self.exons[0].exon_number, 1)
        
        self.assertEqual(self.exons[0].strand, "-")
        
        self.assertEqual(self.exons[0].start, 819113)
        
        self.assertEqual(self.exons[0].end, 819301)
        
class TestValidateTranscriptName(unittest.TestCase):
    
    def test_basic(self):
        
        self.assertRaises(veffects.BadNameError,
                          veffects.validate_transcript_name,
                          "AGAP004687")
    
class TestWorkflowReverseStrand(unittest.TestCase):
        
    ##what about masked sequence?
    
    def setUp(self):
        
        self.gene = "AGAP004687-RA"
        
        self.variants_one = [vr("2L", 819183, "TTGTTCC", "T", 819189),
                             vr("2L", 819219, "C", "A", 819219),
                             vr("2L", 819235, "C", "A", 819235)]
        
        self.variants_two = [vr("2L", 819243, "A", "AG", 819243)]
        
        self.variants_three = [vr("2L", 819174, "A", "AA", 819174),
                               vr("2L", 819138, "AG", "A", 819139)]
        
        self.variants_four = [vr("2L", 819290, "G", "T", 819290),
                              vr("2L", 819233, "C", "G", 819233),
                              vr("2L", 819172, "C", "A", 819172),
                              vr("2L", 819195, "T", "C", 819195),
                              vr("2L", 819165, "A", "T", 819165),
                              vr("2L", 819201, "T", "G", 819201),
                              vr("2L", 819251, "C", "A", 819251),
                              vr("2L", 819134, "G", "A", 819134),
                              vr("2L", 819235, "C", "G", 819235)]
        
        self.variants_five = [vr("2L", 819216, "TTCCGTGTCATGC", "T", 819228),
                             vr("2L", 819185, "G", "GT", 819185),
                             vr("2L", 819182, "CT", "C", 819183),
                             vr("2L", 819160, "T", "G", 819160),
                             vr("2L", 819294, "T", "A", 819294),
                             vr("2L", 819225, "A", "G", 819225),
                             vr("2L", 819167, "C", "T", 819167),
                             vr("2L", 819203, "CT", "C", 819204),
                             vr("2L", 819196, "CGA", "C", 819198),
                             vr("2L", 819154, "T", "G", 819154),
                             vr("2L", 819279, "T", "G", 819279),
                             vr("2L", 819290, "G", "A", 819290),
                             vr("2L", 819145, "C", "G", 819145)]
        
        feature_table_dtype = [
                ('seqid','O'),
                ('source','O'),
                ('type','O'),
                ('start','<i8'),
                ('end','<i8'),
                ('score','<f8'),
                ('strand','O'),
                ('phase','<i8'),
                ('Parent','O'),
                ('Name','O'),
                ('ID','O'),
                ]
        
        feature_table_data = [
                ('2L', 'DB', 'mRNA', 819113, 819301, -1, '-', -1,
                 'AGAP004687', '.', 'AGAP004687-RA'),
                ('2L', 'DB', 'exon', 819113, 819301, -1, '-', -1,
                 'AGAP004687-RA', '.', 'E013854A'),
                ('2L', 'DB', 'CDS', 819113, 819301, -1, '-', 0,
                 'AGAP004687-RA', '.', '.')
                ]
            
        self.feature_table = \
            allel.FeatureTable(feature_table_data, dtype=feature_table_dtype)
         
    def test_variant_set_one(self):
        
        self.transcript = veffects.run_workflow(self.gene, self.feature_table,
                                                self.variants_one)
        
        output_one = "MVHNGIDYGDMQLICEACHLMLSLGMTLKEMVQEFDV*GVLDS" +\
        "FLIEIPHDFLNQRDVEG*"
        
        self.assertEqual(self.transcript.seq_changed_translated, output_one)
        
    def test_variant_set_two(self):
        
        self.transcript = veffects.run_workflow(self.gene, self.feature_table,
                                                self.variants_two)
        
        output_two = "MVHNGIDYGDMQLICEACHPDAGAGHDTEGNGTGVRRVEQGRIG" +\
        "FVPDRDPARFSQPTRR*GI"
                    
        self.assertEqual(self.transcript.seq_changed_translated, output_two)
        
    def test_variant_set_three(self):
        
        self.transcript = veffects.run_workflow(self.gene, self.feature_table,
                                                self.variants_three)
        
        output_three = "MVHNGIDYGDMQLICEACHLMLALGMTRKEMVQEFDVWNKGVF" +\
        "GFVPDRDPARFFNQRDVEG*"
        
        self.assertEqual(self.transcript.seq_changed_translated, output_three)
        
    def test_variant_set_four(self):
        
        self.transcript = veffects.run_workflow(self.gene, self.feature_table,
                                                self.variants_four)
        
        output_four = "MVHKGIDYGDMQLICEACHLMLPLGMTRKEMVQAFGVWNKGVLY" +\
        "SYLIEIPHDFLNQRDVEG*"
        
        self.assertEqual(self.transcript.seq_changed_translated, output_four)
    
    def test_variant_set_five(self):
        
        self.transcript = veffects.run_workflow(self.gene, self.feature_table,
                                                self.variants_five)
        
        output_five = "MVLNGIDSGDMQLICEACHLMLALEEMVRSDVWKQGVLDSFLL" +\
        "ELPHHFLNQRDVEG*"
        
        self.assertEqual(self.transcript.seq_changed_translated, output_five)
    
    
    def test_order_does_not_matter(self):

        v1 = vr("2L", 819219, "C", "A", 819219)
        v2 = vr("2L", 819235, "C", "A", 819235)
        
        transcript_1 = veffects.run_workflow(self.gene, self.feature_table, 
                                             [v1, v2])
        
        transcript_2 = veffects.run_workflow(self.gene, self.feature_table,
                                             [v2, v1])
        
        self.assertEqual(transcript_1.seq_changed_translated,
                         transcript_2.seq_changed_translated)
    
class TestWorkflowForwardStrand(unittest.TestCase):
        
    def setUp(self):
        
        self.gene = "AGAP013717-RA"
        
        self.variants_one = [vr("3R", 758683, "A", "G", 758683),
                             vr("3R", 758690, "CC", "C", 758691),
                             vr("3R", 758693, "A", "C", 758693),
                             vr("3R", 758697, "GAA", "G", 758699),
                             vr("3R", 758708, "C", "T", 758708),
                             vr("3R", 758718, "C", "T", 758718),
                             vr("3R", 758720, "A", "G", 758720),
                             vr("3R", 758727, "G", "A", 758727)]
        
        self.variants_two = [vr("3R", 758681, "G", "C", 758681),
                             vr("3R", 758697, "G", "A", 758697),
                             vr("3R", 758711, "A", "T", 758711),
                             vr("3R", 758713, "TACAGCAAA", "T", 758721),
                             vr("3R", 758722, "A", "C", 758722),
                             vr("3R", 758727, "G", "GTT", 758727)]
        
        feature_table_dtype = [
                ('seqid','O'),
                ('source','O'),
                ('type','O'),
                ('start','<i8'),
                ('end','<i8'),
                ('score','<f8'),
                ('strand','O'),
                ('phase','<i8'),
                ('Parent','O'),
                ('Name','O'),
                ('ID','O'),
                ]
        
        feature_table_data = [
                ('3R', 'DB', 'mRNA', 758261, 759686, -1, '+', -1,
                 'AGAP013717', '.', 'AGAP013717-RA'),
                 
                ('3R', 'DB', 'five_prime_UTR', 758261, 758676, -1, '+', -1,
                 'AGAP013717-RA', '.', '.'),
                 
                ('3R', 'DB', 'exon', 758261, 758733, -1, '+', -1,
                 'AGAP013717-RA', 'AGAP013717-RA-E1', '.'),
                 
                ('3R', 'DB', 'CDS', 758677, 758733, -1, '+', 0,
                 'AGAP013717-RA', '.', 'AGAP013717-PA'),
                 
                ('3R', 'DB', 'CDS', 759210, 759212, -1, '+', 0,
                 'AGAP013717-RA', '.', 'AGAP013717-PA'),
                 
                ('3R', 'DB', 'exon', 759210, 759686, -1, '+', -1,
                 'AGAP013717-RA', 'AGAP013717-RA-E2', '.'),
                 
                ('3R', 'DB', 'three_prime_UTR', 759213, 759686, -1, '+', -1,
                 'AGAP013717-RA', '.', '.')
                ]
            
        self.feature_table = \
            allel.FeatureTable(feature_table_data, dtype=feature_table_dtype)

    def test_variant_set_one(self):
        
        self.transcript = veffects.run_workflow(self.gene, self.feature_table,
                                                self.variants_one)
        
        output_one = "MGDAPPRQYIKYSRNTAE*"
        
        self.assertEqual(self.transcript.seq_changed_translated, output_one)

    def test_variant_set_two(self):
        
        self.transcript = veffects.run_workflow(self.gene, self.feature_table,
                                                self.variants_two)
        
        output_two = "MANAPNPKQYTMSHVAE*"
        
        self.assertEqual(self.transcript.seq_changed_translated, output_two)
        
if __name__ == '__main__':
	unittest.main()
