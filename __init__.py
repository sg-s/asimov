import os

asimov_root = os.path.dirname(os.path.realpath(__file__))

list_of_feeds = []

# check that the config file exists 
if not os.path.isfile(os.path.join(asimov_root,'config.txt')):
    raise AssertionError('Could not find config file.')

if os.path.isfile(os.path.join(asimov_root,'feeds.txt')):
	print("List of feeds exists, need to load")