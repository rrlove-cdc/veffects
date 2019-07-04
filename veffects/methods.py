#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May  3 13:37:17 2018

@author: beccalove
"""

from collections import namedtuple
import re
import requests

from .Transcript import Transcript
from .ExonSequence import ExonSequence

VariantRecord = namedtuple(
    'VariantRecord',['chrom', 'pos', 'ref', 'alt', 'end'])


class BadNameError(Exception):
    pass

class NumExonsAndCDSDifferError(Exception):
    pass

class RequestReturnError(Exception):
    pass

class TimeOutError(Exception):
    pass

def validate_transcript_name(name):
    
    if not re.search('-R[A-Z]', name):
        
        raise BadNameError("You may be using a gene name.\
                         Please use a transcript name instead.")

def make_POST_request(gene, 
                      post_server = "https://www.vectorbase.org/rest", 
                      post_ext = "/sequence/id/",
                     feature = "cds",
                     timeout = (0.1, None),
                     session = None):
    
    validate_transcript_name(gene)
    
    headers = {'Content-type' : 'application/json',
               'Accept' : 'application/json'}
    
    post_string = post_server + post_ext
    
    feature_seq_payload = '{"ids" : ["' + gene + '"],\
    "type" : "' + feature + '"}'
    
    if session:
        
        try:
        
            feature_seq = session.post(post_string,
                                       data = feature_seq_payload,
                                       headers = headers,
                                       timeout = timeout)
            
        except requests.exceptions.Timeout:
            
            raise TimeOutError()
        
    else:
        
        try:
    
            feature_seq = requests.post(post_string,
                                        data = feature_seq_payload,
                                        headers = headers,
                                        timeout = timeout)
            
        except requests.exceptions.Timeout:
            
            raise TimeOutError
    
    if not feature_seq.status_code == requests.codes.ok:
        
        raise RequestReturnError("POST request status: ", 
                                 feature_seq.status_code)
    
    return feature_seq.json()

def make_transcript(feature_seq_json):
    
    transcript = Transcript(name = feature_seq_json[0]["id"], 
                            seq = feature_seq_json[0]["seq"])
    
    return transcript

def make_GET_request(gene,
                    get_server = "https://www.vectorbase.org/rest",
                    feature_types = ["exon","cds"]):
    
    validate_transcript_name(gene)
    
    headers = {'Content-type' : 'application/json', \
               'Accept' : 'application/json'}
    
    get_ext = "/overlap/id/" + gene + "?"
    
    get_string =\
    get_server + get_ext +\
    "".join(["feature=" + feature +\
             ";" for feature in feature_types]).rstrip(";")
    
    feature_coords = requests.get(get_string, headers = headers)
    
    if not feature_coords.status_code == requests.codes.ok:
        
        raise RequestReturnError("GET request status: ",
                                 feature_coords.status_code)
        
    return feature_coords.json()

def make_exons(gff3, gene):
    
    chunk = gff3[(gff3["Parent"] == gene) | (gff3["ID"] == gene) |\
                 (gff3["Name"] == gene)]
    
    cds_list = []
    exon_object_list = []
    
    for item in chunk:
    
        if item["type"] == "CDS":

            cds_list.append(item)
                
    if cds_list[0]["strand"] == "-":
        
        cds_list.sort(key = lambda x: x["start"], reverse=True)
        
    else:
        
        cds_list.sort(key = lambda x: x["start"])
        
    for index, feature in enumerate(cds_list):
    
        exon = ExonSequence(chrom = feature["seqid"],
                        transcript = gene,
                       exon_number = index + 1,
                       strand = feature["strand"],
                       start = feature["start"],
                       end = feature["end"])
        
        exon_object_list.append(exon)
        
    return exon_object_list

def check_exon_order(exon_list):
    
    if exon_list[0].strand == "-":
        
        exon_list.reverse()
        
    for i, x in enumerate(exon_list[:-1]):
        
        y = exon_list[i+1]
        
        if not x.end <= y.start:
            
            raise ValueError("exons out of order")

def run_workflow(gene_name, gff3, variants, **kwargs):
    
    transcript = make_transcript(make_POST_request(gene = gene_name, **kwargs))
    
    for exon in make_exons(gff3, gene_name):
        transcript.add_exon(exon)
        
    transcript.populate_exon_seq()
    
    transcript.check_for_overlapping_variants(variants)
    
    transcript.parse_variants_list(variants)
    
    for exon in transcript.exons:
    
        exon.change()
        
    transcript.assemble_changed_seq()
    
    transcript.translate_seqs()
    
    return transcript