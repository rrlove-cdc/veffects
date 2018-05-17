#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 16 15:48:51 2018

@author: beccalove
"""

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

        self.assertRaises(veffects.methods.RequestReturnError, 
                          veffects.make_POST_request, "foo-RZ")
        
    def test_gene_name_not_transcript_name(self):
        
        self.assertRaises(ValueError, veffects.make_POST_request,
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
        
        self.assertRaises(ValueError, veffects.make_GET_request, 
                          "AGAP004687")
        
    def test_bad_name(self):
        
        self.assertRaises(veffects.methods.RequestReturnError,
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
        
        self.feature_coords_json = [{'Parent': 'AGAP004687-RA',
                          'assembly_name': 'AgamP4',
                          'constitutive': '1',
                          'end': 819301,
                          'ensembl_end_phase': '0',
                          'ensembl_phase': '0',
                          'exon_id': 'AGAP004687-RA-E1',
                          'feature_type': 'exon',
                          'id': 'AGAP004687-RA-E1',
                          'rank': 1,
                          'seq_region_name': '2L',
                          'source': 'VectorBase',
                          'start': 819113,
                          'strand': -1,
                          'version': 1},
                            {'Parent': 'AGAP004687-RA',
                             'assembly_name': 'AgamP4',
                             'end': 819301,
                             'feature_type': 'cds',
                             'id': 'AGAP004687-PA',
                             'phase': 0,
                             'protein_id': 'AGAP004687-PA',
                             'seq_region_name': '2L',
                             'source': 'VectorBase',
                             'start': 819113,
                             'strand': -1}]
        
    def test_basic(self):
        
        self.exons = veffects.make_exons(self.feature_coords_json,
                                         "AGAP004687-RA")
        
        self.assertEqual(len(self.exons), 1)
        
        self.assertEqual(self.exons[0].chrom, "2L")
        
        self.assertEqual(self.exons[0].transcript, "AGAP004687-RA")
        
        self.assertEqual(self.exons[0].exon_number, 1)
        
        self.assertEqual(self.exons[0].strand, -1)
        
        self.assertEqual(self.exons[0].start, 819113)
        
        self.assertEqual(self.exons[0].end, 819301)
        
class TestValidateTranscriptName(unittest.TestCase):
    
    def test_basic(self):
        
        self.assertRaises(ValueError,
                          veffects.validate_transcript_name,
                          "AGAP004687")
    
class TestWorkflow(unittest.TestCase):
        
    def setUp(self):
        
        self.gene = "AGAP004687-RA"
        
        self.variants_one = [vr("2L", 819190, "TGGAACA", "T"),
                             vr("2L", 819219, "G", "T"),
                             vr("2L", 819235, "G", "T")]
        
        self.variants_two = [vr("2L", 819244, "C", "CC")]
        
        self.variants_three = [vr("2L", 819175, "T", "TT"),
                               vr("2L", 819140, "TC", "T")]
        
        self.variants_four = [vr("2L", 819290, "C", "A"),
                              vr("2L", 819233, "G", "C"),
                              vr("2L", 819172, "G", "T"),
                              vr("2L", 819195, "A", "G"),
                              vr("2L", 819165, "T", "A"),
                              vr("2L", 819201, "A", "C"),
                              vr("2L", 819251, "G", "T"),
                              vr("2L", 819134, "C", "T"),
                              vr("2L", 819235, "G", "C")]
        
        self.variants_five = [vr("2L", 819229, "GGCATGACACGGA", "G"),
                             vr("2L", 819186, "A", "AA"),
                             vr("2L", 819184, "AA", "A"),
                             vr("2L", 819160, "A", "C"),
                             vr("2L", 819294, "A", "T"),
                             vr("2L", 819225, "T", "C"),
                             vr("2L", 819167, "G", "A"),
                             vr("2L", 819205, "CA", "C"),
                             vr("2L", 819199, "TTC", "T"),
                             vr("2L", 819154, "A", "C"),
                             vr("2L", 819279, "A", "C"),
                             vr("2L", 819290, "C", "T"),
                             vr("2L", 819145, "G", "C")]
        
    def test_variant_set_one(self):
        
        self.transcript = veffects.run_workflow(self.gene,
                                                self.variants_one)
        
        output_one = "MVHNGIDYGDMQLICEACHLMLSLGMTLKEMVQEFDV*GVLDS" +\
        "FLIEIPHDFLNQRDVEG*"
        
        self.assertEqual(self.transcript.seq_changed_translated, output_one)
        
    def test_variant_set_two(self):
        
        self.transcript = veffects.run_workflow(self.gene, 
                                                self.variants_two)
        
        output_two = "MVHNGIDYGDMQLICEACHPDAGAGHDTEGNGTGVRRVEQGRIG" +\
        "FVPDRDPARFSQPTRR*GI"
                    
        self.assertEqual(self.transcript.seq_changed_translated, output_two)
        
    def test_variant_set_three(self):
        
        self.transcript = veffects.run_workflow(self.gene,
                                                self.variants_three)
        
        output_three = "MVHNGIDYGDMQLICEACHLMLALGMTRKEMVQEFDVWNKGVF" +\
        "GFVPDRDPARFFNQRDVEG*"
        
        self.assertEqual(self.transcript.seq_changed_translated, output_three)
        
    def test_variant_set_four(self):
        
        self.transcript = veffects.run_workflow(self.gene, self.variants_four)
        
        output_four = "MVHKGIDYGDMQLICEACHLMLPLGMTRKEMVQAFGVWNKGVLY" +\
        "SYLIEIPHDFLNQRDVEG*"
        
        self.assertEqual(self.transcript.seq_changed_translated, output_four)
    
    def test_variant_set_five(self):
        
        self.transcript = veffects.run_workflow(self.gene, self.variants_five)
        
        output_five = "MVLNGIDSGDMQLICEACHLMLALEEMVRS*DVWKQGVLDSFLL" +\
        "ELPHHFLNQRDVEG*"
        
        self.assertEqual(self.transcript.seq_changed_translated, output_five)
    
    
    def test_positive_strand(self):
        
        raise ValueError("test not written")
        
        
        
if __name__ == '__main__':
	unittest.main()
