import os

this_dir = os.path.dirname(os.path.realpath(__file__))

# check that the config file exists 
if not os.path.isfile(os.path.join(os.path.dirname(this_dir),'config.txt')):
    raise AssertionError('Could not find config file.')

# from flask import Flask
# import os
# import time
# from app.pubmed_filter import fetch_and_filter
# from apscheduler.schedulers.background import BackgroundScheduler

# app = Flask(__name__)
# from app import views


# @app.before_first_request
# def schedule_fetch():
#     os.environ['TZ'] = 'US/Eastern'
#     time.tzset()
#     scheduler = BackgroundScheduler()
#     # scheduler.add_job(main, 'cron', day_of_week='*', hour='17', minute='21')
#     scheduler.add_job(fetch_and_filter, 'cron', hour='3')
#     scheduler.start()
#     print('Press Ctrl+{0} to exit'.format('Break' if os.name == 'nt' else 'C'))
