#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May  3 13:37:02 2018

@author: beccalove
"""

from Bio import Seq

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
        However, we still want to retain the original sequence as well.'''
        
        self.seq_adjusted = sequence
            
        if self.strand == "-1":
            
            reverse_complement = str(Seq.Seq(sequence).reverse_complement())
            
            self.seq_adjusted = reverse_complement
        
    def change(self):
        
        ##what to do about overlapping variants?
        ##if we are assuming these are haplotypes, there shouldn't be any
        ##make that assumption explicit?
        ##it seems like a case with overlapping variants would break the
        ##assumptions behind this code
        
        mutable_sequence_list = list(self.seq_adjusted)
                
        for variant in self.variants:
            
            if not variant.chrom == self.chrom:
                raise ValueError("Variant chrom and exon chrom don't match", 
                                 variant)
                
            if not self.start <= variant.pos <= self.end:
                raise ValueError("Variant is out of bounds of exon", variant)

            index = variant.pos - self.start
                
#            if self.strand == "-1":
                
#                index = self.end - variant.pos
            
            if len(variant.alt) >= len(variant.ref):
                ##variant is insertion or SNP
                
                if not variant.ref == self.seq_adjusted[index]:
                    raise ValueError("Reference alleles don't match at " +\
                                     str(variant.pos))

                if mutable_sequence_list[index] == '':
                    
                    print("Warning: position at {pos} is part of a \
                          spanning deletion. You may be analyzing data \
                          from multiple haplotypes or have an undetected \
                          error. Skipping this position.".format(pos =\
                          variant.pos))
                    
                    continue
                
                mutable_sequence_list[index] = variant.alt
                
            elif len(variant.ref) > len(variant.alt):##variant is a deletion
                
                test_string =\
                self.seq_adjusted[(index):(index + len(variant.ref))]
                
                if not variant.ref == test_string:
                    raise ValueError("Reference alleles don't match at " +\
                                     str(variant.pos))
                    
                to_be_deleted =\
                mutable_sequence_list[(index):(index + len(variant.ref))]
                    
                if '' in to_be_deleted:
                    
                    print("Warning: position at {pos} is part of a \
                          spanning deletion. You may be analyzing data \
                          from multiple haplotypes or have an undetected \
                          error. Skipping this position.".format(pos =\
                          variant.pos))
                    
                    continue
                
                mutable_sequence_list[index] = variant.alt
                
                for new_index in range((index+1),(index + len(variant.ref))):
                    mutable_sequence_list[new_index] = ''
        
        changed_sequence = ''.join(mutable_sequence_list)
        
        ##if this is on the reverse strand, we need to change back
        ##to the original orientation
        
        if self.strand == "-1":
            
            changed_sequence =\
            str(Seq.Seq(changed_sequence).reverse_complement())
        
        self.changed_sequence = changed_sequence        