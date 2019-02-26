# RSS module for asimov

import feedparser
import os
from asimov import asimov_root

# fetch all RSS feeds
def fetch():
    feed_list = getFeeds()

    for feed in feed_list:
        print("Getting feed:")
        print(feed)
        d = feedparser.parse(feed)
        items = d.entries
        if len(items) == 0:
            print("No items in this feed!")
            continue
        print("# items found:")
        print(len(items))

        for item in items:
            if len(item.summary) > 0:
                d = {'title' : item.title, 'link' : item.link, 'summary' : item.summary, 'updated' : item.updated}
                with open('test.json','a') as fp:
                 

# reads the cached entries from a particular feed
def readCache(feed):



# add RSS feed to list of feeds
def addFeed(feed_name):
    print("adding feed...")
    if not isinstance(feed_name,str):
        raise AssertionError("Expected a string, got something else")

    d = feedparser.parse(feed_name)
    if len(d.entries) == 0:
        raise AssertionError("Feed URL has no entries; and is probably wrong")

    # add to list of feeds if it doesn't already exist in the list of feeds
    feed_list = getFeeds()

    add_to_list = True
    if len(feed_list) > 0:
        for feed in feed_list:
            if feed == feed_name:
                add_to_list = False;
                print("Not adding feed -- alreacy exists in list")

    if add_to_list:
        with open(os.path.join(asimov_root,'feeds.txt'), "a") as myfile:
            myfile.write(feed_name)

def getFeeds():
    # returns a list of feeds 

    feed_list = []

    this_file = os.path.join(asimov_root,'feeds.txt')
    if os.path.isfile(this_file):
        # read this file
        print("feeds.txt already exists, need to read it")

        # get all the feeds and the names of the feeds 
        with open(this_file) as f: 
            l = f.readlines()

        for i in range(len(l)):
            this_line = l[i]
            feed_list.append(this_line)

    return feed_list



# show list of feeds
