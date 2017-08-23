# README

Know everything, like [Asimov](https://en.wikipedia.org/wiki/Isaac_Asimov) did. 

`asimov` is a python module that uses machine learning to help you keep abrest of new scientific journal articles. 

# General principle

`asimov` uses a corpus of text, that you identify as being interesting or not, to train a model whose job it is to predict if a new piece of text is something that you're interested in. 

# Usage 

## 0. Configuration

Copy `default.txt` to `config.txt` and fill out the missing values. This information is needed for `asimov` to work. Then, import `asimov`:

```python
import asimov
```

## 1. Generate training data 

Training data consists of a list of abstracts (or any text that you want to train on), and a corresponding list of labels identifying whether each piece of text is to be trained for, or trained against. 

`asimov` can generate positive training data from a BibTex file containing papers that you find interesting. Once you have exported all the papers you are interested in a .bib file, use `asimov` to pull some useful text from the .bib file:

```python
asimov.get_cool_papers_from_bib()
```

You may not have abstracts for every paper in this .bib file. If this is the case, you can use PubMed to get these abtracts:

```python
asimov.get_missing_abstracts()
```
Hopefully, this dataset is mostly complete, and somewhat useable. 

### Generating negative training data 

`asimov` can also scan RSS feeds to determine a corpus of text to feed into the negative training data. To use this, you need to edit the file called `boring-feeds.txt` and add links to RSS feeds from sources that would never interest you. You might want to edit this file and remove entries that you might be interested in. Once this done, 

```python
asimov.get_boring_papers_from_rss()
```

# Installation 

Download this repository. 

# License 

[GPL v3](http://gplv3.fsf.org/)

