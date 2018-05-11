#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May  3 13:37:02 2018

@author: beccalove
"""

class ExonSequence:
    
    def __init__(self, chrom, transcript, exon_number, strand, start, end):
        
        self.chrom = chrom
        self.transcript = transcript
        self.exon_number = exon_number
        self.strand = strand
        self.start = start
        self.end = end
        self.length = abs(self.end - self.start) + 1
        self.variants = []
        self.sequence = ''
        self.changed_sequence = ''
        
    def add_sequence(self, sequence):
        
        if not len(sequence) == self.length:
            
           raise ValueError("Sequence and exon differ in length")
        
        self.sequence = sequence
        
    def change(self, variants):
        
        ##what to do about overlapping variants?
        ##if we are assuming these are haplotypes, there shouldn't be any
        ##make that assumption explicit?
        ##it seems like a case with overlapping variants would break the
        ##assumptions behind this code
        
        mutable_sequence_list = list(self.sequence)
        
        for variant in variants:
            
            if not variant.chrom == self.chrom:
                raise ValueError("Variant chrom and exon chrom don't match", 
                                 variant)
            
            if not self.start <= variant.pos <= self.end:
                raise ValueError("Variant is out of bounds of exon", variant)
            
            index = variant.pos - self.start
            
            if len(variant.alt) >= len(variant.ref):
                ##variant is insertion or SNP
                
                if not variant.ref == self.sequence[index]:
                    raise ValueError("Reference alleles don't match at" +\
                                     variant.pos)
                
                mutable_sequence_list[index] = variant.alt
                
            elif len(variant.ref) > len(variant.alt):##variant is a deletion
                
                test_string = self.sequence[(index):(index + len(variant.ref))]
                
                if not variant.ref == test_string:
                    raise ValueError("Reference alleles don't match at" +\
                                     variant.pos)
                
                mutable_sequence_list[index] = variant.alt
                
                for index in range((index+1),(index + len(variant.ref))):
                    mutable_sequence_list[index] = ''
                
        self.changed_sequence = ''.join(mutable_sequence_list)
