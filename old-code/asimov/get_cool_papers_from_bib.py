#! /usr/bin/env python
# get_cool_papers_from_bib
# accepts a .bib file, and scrapes title, abstract information from it, and saves it to disk in a easy-to-read list format

import os
import sys
import bibtexparser 
import re
import pickle

def get_cool_papers_from_bib(**kwargs):
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
		elif len(bib_files) == 0:
			raise NameError("No .bib files found, cannot proceed");
		else:
			bibfile = bib_files[0]

	print("Using file {}".format(bibfile))

	# we make a set of lists for pubmed id, whether there is an abstract or not,
	# the abstract, the title, first author, year published, etc.

	all_years = []
	all_first_authors = []
	all_titles = []
	all_abstracts = []

	# first read the bibtext file and pull out abstracts when we can
	with open(os.path.join(dir_path,bibfile), encoding='utf-8') as bibtex_file:
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

			# pubmed only works with the last name of the first author.
			# providing the first name breaks the search because pubmed doesn't know first names
			first_author = first_author.replace(',', '')
			first_author = first_author[0:first_author.find(' ') + 2]

		# now append them all to the list
		all_abstracts.append(this_abstract)
		all_first_authors.append(first_author)
		all_titles.append(this_title)
		all_years.append(this_year)


	# now save them all to disk
	# Saving the objects:
	with open(os.path.join(dir_path, 'bibtex.pickle'), 'wb') as f:
		pickle.dump([all_abstracts, all_first_authors, all_titles, all_years], f)



if __name__ == "__main__":
	get_cool_papers_from_bib()
