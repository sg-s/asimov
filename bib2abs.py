#! /usr/bin/env python
# bib2abs.py
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

def bib2abs(**kwargs):
	if kwargs:
		bibfile = kwargs.get('bibfile')
	else:
		print("no optional argument specified, need to look for bib files")

		dir_path = os.path.dirname(os.path.realpath(__file__))
		bib_files = []
		for file in os.listdir(dir_path):
			if file.endswith(".bib"):
				bib_files.append(file)

		if len(bib_files) > 1:
			raise NameError("more than one bib file, cannot proceed");     
		else:
			bibfile = bib_files[0]

	print("Using file {}".format(bibfile))

	# read config file and get the email we need to use to connect to PubMed
	cf = configparser.RawConfigParser()   
	configFilePath = r'config.txt'
	cf.read(configFilePath)
	Entrez.email = cf.get('asimov-config', 'pubmed_email')
	pubmed_batch_size = int(cf.get('asimov-config','pubmed_batch_size'))


	# we make a set of lists for pubmed id, whether there is an abstract or not,
	# the abstract, the title, first author, year published, etc.

	all_years = []
	all_first_authors = []
	all_titles = []
	all_abstracts = []

	# first read the bibtext file and pull out abstracts when we can
	with open(bibfile, encoding='utf-8') as bibtex_file:
		bibtexfile = bibtexparser.load(bibtex_file)

	for i in range(len(bibtexfile.entries)):
		this_entry = bibtexfile.entries[i]

		this_title = None
		this_abstract = None
		this_first_author = None
		this_year = None

		# pull out entry title 
		if 'title' in this_entry.keys():
			this_title = this_entry['title']
			# remove things in curly braces from the title
			this_title = re.sub('{.*?}', '', this_title)

		# pull out year  
		if 'year' in this_entry.keys():
			this_year = this_entry['year']

		# pull out abstract
		if 'abstract' in this_entry.keys():
			this_abstract = this_entry['abstract'];
			# sometimes the word "abstract occurs in it"
			this_abstract = this_abstract.replace('Abstract','')
			# make sure it is long enough 
			if len(this_abstract) < 100:
				this_abstract = None

		# pull out authors
		if 'author' in this_entry.keys():
			all_authors = this_entry['author']
			if all_authors.find('and') != -1:
				# more than one author, let's ignore everyone but the first
				first_author = all_authors[0:all_authors.find('and')-1]
			else:
				# only one author
				first_author = all_authors

			# pubmed only works with the last name of the first author. providing the first name breaks the search because pubmed doesn't know first names
			first_author = first_author.replace(',','')
			first_author = first_author[0:first_author.find(' ') + 2]	

		# now append them all to the list
		all_abstracts.append(this_abstract)
		all_first_authors.append(first_author)
		all_titles.append(this_title)
		all_years.append(this_year)


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

if __name__ == "__main__":
	bib2abs()