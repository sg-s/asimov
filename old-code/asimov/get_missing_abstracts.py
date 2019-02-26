#! /usr/bin/env python
# get_missing_abstracts.py
# reads a .bib file and scrapes Pubmed for abstracts
# usage: 
# bib2abs()
# 
# if you want to specify a .bib to work with, use 
# bib2abs(bibfile="/path/to/file.bib")

import os
import sys
import bibtexparser 
from Bio import Entrez
import configparser
import re
import pickle

def get_missing_abstracts():

    # load pickled data
    this_dir = os.path.dirname(os.path.realpath(__file__))
    if not os.path.isfile(os.path.join(this_dir,'bibtex.pickle')):
        raise AssertionError('Could not find any pickled data to load. Try running pickle_bib() first on your .bib file')

    with open(os.path.join(this_dir,'bibtex.pickle'),'rb') as f:
        all_abstracts, all_first_authors, all_titles, all_years = pickle.load(f)

	# read config file and get the email we need to use to connect to PubMed
    cf = configparser.RawConfigParser()   
    configFilePath = r'config.txt'
    cf.read(configFilePath)
    Entrez.email = cf.get('asimov-config', 'pubmed_email')
    pubmed_batch_size = int(cf.get('asimov-config','pubmed_batch_size'))

    # print out how many abstracts we have, and how many are missing. 
    n_missing_abs = len([i for i in all_abstracts if i == None]);
    print('Of the {} entries in this dataset, {} abstracts are missing.'.format(len(all_abstracts),n_missing_abs))

    # step 1: for papers without abstracts, search for them on pubmed and recover their all_pmids
    all_pmids = [None]*len(all_abstracts)

    for i in range(len(all_abstracts)):
        if all_abstracts[i] == None:
            # for the ones with no abstract, let's try scraping pubmed for the abstracts
            search_string_template = '((this_title [Title]) AND first_author [Author - First]) AND ("this_year"[Date - Publication] : "this_year"[Date - Publication])'

            # make the search string
            search_string = search_string_template
            search_string = search_string.replace('this_title',all_titles[i])
            search_string = search_string.replace('first_author',all_first_authors[i])
            search_string = search_string.replace('this_year',all_years[i])

            r = Entrez.read(Entrez.esearch(db="pubmed", term = search_string, rettype="xml"))
            if float(r['Count']) == 0:
                print("zero matches")
            elif float(r['Count']) == 1:
                print('found this paper on pubmed, adding PMID to database...')
                all_pmids[i] = r['IdList'][0]
                print(all_pmids[i])
            else : 
                print("more than 1 match, cannot uniquely identify this paper")


    print('done collecting PMIDs for articles with no abstracts. Now need to download abstracts from Pubmed...')

    idx = [i for i, x in enumerate(all_pmids) if x != None]
    ret_these = [x for i, x in enumerate(all_pmids) if x != None]

    # now retrieve abstracts from pubmed in batches using efetch 
    article_counter = 0
    pubmed_abstracts = [None]*len(ret_these)
    while article_counter < len(ret_these):
        pmid_batch = ret_these[article_counter : article_counter + pubmed_batch_size]
        a = Entrez.read(Entrez.efetch(db="pubmed", id = pmid_batch, rettype="xml"))
        articles = a['PubmedArticle']
        for i in range(len(articles)):
            # check if this article has an abstract
            if 'Abstract' in articles[i]['MedlineCitation']['Article'].keys():
                this_abstract = articles[i]['MedlineCitation']['Article']['Abstract']['AbstractText'][0]
                this_abstract = this_abstract[:]
            else:
                this_abstract = None
            pubmed_abstracts[article_counter] = this_abstract
            article_counter += 1

    print(pubmed_abstracts)

    # put these retrieved abstracts back in the global abstracts list 
    for i, pa in zip(idx,pubmed_abstracts):
        all_abstracts[i] = pa

    n_missing_abs = len([i for i in all_abstracts if i == None]);
    print('After searching Pubmed, of the {} entries in this dataset, {} abstracts are missing.'.format(len(all_abstracts),n_missing_abs))

    # finally, re-pickle all the data 
    print('saving the newly acquired abstracts....')
    with open(os.path.join(dir_path,'bibtex.pickle'), 'wb') as f:
        pickle.dump([all_abstracts, all_first_authors, all_titles, all_years], f)

if __name__ == "__main__":
    get_missing_abstracts()