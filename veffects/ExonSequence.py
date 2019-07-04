#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May  3 13:37:02 2018

@author: beccalove
"""

from Bio import Seq

class OverlappingVariantError(Exception):
    pass

class ExonSequence:
    
    def __init__(self, chrom, transcript, exon_number, strand, start, end):
        
        self.chrom = chrom
        self.transcript = transcript
        self.exon_number = exon_number
        self.strand = str(strand)
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
        
        '''If this is on the reverse strand, take the reverse complement,
        because variants are reported in relation to the forward strand.
        However, we still want to retain the original sequence as well.
        '''
        
        self.seq_adjusted = sequence
            
        if self.strand == "-":
            
            reverse_complement = str(Seq.Seq(sequence).reverse_complement())
            
            self.seq_adjusted = reverse_complement
        
    def change(self):
        
        mutable_sequence_list = list(self.seq_adjusted)
                
        for variant in self.variants:
            
            # Does the chromosome match?
            if not variant.chrom == self.chrom:
                raise ValueError("Variant chrom and exon chrom don't match", 
                                 variant)

            index = variant.pos - self.start
            
            len_ref = len(variant.ref)
            
            len_alt = len(variant.alt)
            
            if not variant.ref == self.seq_adjusted[(index):(index + len_ref)]:
                    raise ValueError("Reference alleles don't match at " +\
                                     str(variant.pos))
                    
            # Check for spanning deletions
            
            affected = mutable_sequence_list[(index):(index + len_ref)]
            
            if '' in affected:
                
                raise OverlappingVariantError(variant.pos)
                        
            if len_ref == len_alt:
                
                for alt_index in range(len_alt):
                    
                    new_index = alt_index + index
                    
                    mutable_sequence_list[new_index] = variant.alt[alt_index]
                    
            else:
            
                mutable_sequence_list[index] = variant.alt
            
                if len_ref > 1:
                    
                    for new_index in range((index + 1), (index + len_ref)):
                        
                        mutable_sequence_list[new_index] = ''
        
        changed_sequence = ''.join(mutable_sequence_list)
        
        # Change the reverse strand back to the original orientation
        
        if self.strand == "-":
            
            changed_sequence =\
            str(Seq.Seq(changed_sequence).reverse_complement())
        
        self.changed_sequence = changed_sequence        