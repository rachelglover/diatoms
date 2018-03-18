from setuptools import setup, find_packages
from codecs import open
from os import path

here = path.abspath(path.dirname(__file__))

with open(path.join(here, 'README.rst'), encoding='utf-8') as f:
    long_description = f.read()

setup (
    name = 'diatompipeline',
    version = '1.0.0',
    description = 'Pipeline for diatom sample identification',
    author = 'Rachel Glover',
    author_email = 'rachel.glover@taxagenomics.com',
    classifiers = [
        'Development Status :: 4 - Beta',
        'Intended Audience :: Scientists',
        'Programming Lanugage :: Python :: 2',
        'Programming Lanugage :: Python :: 2.7',
    ],
    py_modules=["pipeline"],
    install_requires=[
        "pygal >= 2.4.0",
        "biopython >= 1.70",
        "cutadapt >= 1.16",
        "pandas",
    ],

)