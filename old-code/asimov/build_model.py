#! /usr/bin/env python
# builds model, once you have constructed your training data

import os
import time
import pickle
from random import randint
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import LabelEncoder
from sklearn.linear_model import SGDClassifier
from sklearn.metrics import classification_report as clsr
from sklearn.feature_extraction.text import TfidfVectorizer
from sklearn.model_selection import train_test_split as tts



def timeit(func):
    """
    Simple timing decorator
    """
    def wrapper(*args, **kwargs):
        start = time.time()
        result = func(*args, **kwargs)
        delta = time.time() - start
        return result, delta
    return wrapper


def identity(arg):
    """
    Simple identity function works as a passthrough.
    """
    return arg


@timeit
def build_and_evaluate(X, y, classifier=SGDClassifier, outpath=None,
                       verbose=True):

    @timeit
    def build(classifier, X, y=None):
        """
        Inner build function that builds a single model.
        """

        model = Pipeline([
            ('vectorizer', TfidfVectorizer()),
            ('classifier', classifier()),
        ])

        model.fit(X, y)
        return model

    # Label encode the targets
    labels = LabelEncoder()
    y = labels.fit_transform(y)

    # Begin evaluation
    if verbose:
        print("Building for evaluation")
    X_train, X_test, y_train, y_test = tts(X, y, test_size=0.2)
    model, secs = build(classifier, X_train, y_train)

    if verbose:
        print("Evaluation model fit in {:0.3f} seconds".format(secs))
    if verbose:
        print("Classification Report:\n")

    y_pred = model.predict(X_test)
    print(clsr(y_test, y_pred, target_names=labels.classes_))

    if verbose:
        print("Building complete model and saving ...")
    model, secs = build(classifier, X, y)
    model.labels_ = labels

    if verbose:
        print("Complete model fit in {:0.3f} seconds".format(secs))

    if outpath:
        with open(outpath, 'wb') as f:
            pickle.dump(model, f)

        print("Model written out to {}".format(outpath))

    return model


def build_model():
    # load the data 
    this_dir = os.path.dirname(os.path.realpath(__file__))
    if not os.path.isfile(os.path.join(this_dir,'bibtex.pickle')):
        raise AssertionError('Could not find any pickled data to load. Try running pickle_bib() first on your .bib file')

    with open(os.path.join(this_dir,'bibtex.pickle'),'rb') as f:
        all_abstracts, all_first_authors, all_titles, all_years = pickle.load(f)

    # cut out all elements where we don't have the abstracts
    ok_values = [i for i,x in enumerate(all_abstracts) if x!=None]
    ok_abs = [all_abstracts[i] for i in ok_values]
    ok_titles = [all_titles[i] for i in ok_values]

    X = ["{} {}".format(a_, b_) for a_, b_ in zip(ok_titles, ok_abs)]
    y = ['pos']*len(X)

    print('done getting data from database')
    model, secs = build_and_evaluate(X, y, outpath='./model.pickle')


if __name__ == "__main__":
    build_model()