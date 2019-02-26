#! /usr/bin/env python
# get_boring_papers_from_rss
# get a collection of abstracts from papers you find boring to use as negative training data

import os
import sys
import feedparser 
import html2text
import pickle

def get_boring_papers_from_rss():

    dir_path = os.path.dirname(os.path.realpath(__file__))
    this_file = os.path.join(dir_path,'boring-feeds.txt')

    # get all the feeds and the names of the feeds 
    with open(this_file) as f: 
        l = f.readlines()

    all_feed_urls = []
    all_feed_names = []
    for i in range(len(l)):
        this_line = l[i]
        all_feed_urls.append(this_line[this_line.find('|')+1:].strip())
        all_feed_names.append(this_line[:this_line.find('|')].strip())

    print('Successfully read list of RSS feeds...')
	# make a list to hold the abstracts (here, we lump tirles, 
    # abstracts, and anything else in together)

    all_abstracts = []
    for i in range(len(all_feed_urls)):
        print('Reading RSS feeds from {}'.format(all_feed_names[i]))
        d = feedparser.parse(all_feed_urls[i])
        entries = d['entries']
        print('read {} entries'.format(len(entries)))
        for j in range(len(entries)):
            this_txt = entries[j]['title'] + ' ' + entries[j]['summary']
            all_abstracts.append(this_txt)

	# save all thse things 
    with open(os.path.join(dir_path, 'boring_papers.pickle'), 'wb') as f:
        pickle.dump(all_abstracts, f)


if __name__ == "__main__":
	get_boring_papers_from_rss()
