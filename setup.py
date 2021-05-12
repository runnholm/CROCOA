import os
from setuptools import setup


# Utility function to read the README file.
# Used for the long_description.  It's nice, because now 1) we have a top level
# README file and 2) it's easier to type in the README file than to put a raw
# string in below ...
def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()


setup(
    name="CROCOA",
    version="0.0.1",
    author="Axel Runnholm",
    author_email="axel.runnholm@astro.su.se",
    description=("Tools for doing image alignment using 2D correlation"),
    license="GNU GPLv3",
    packages=['CROCOA'],
    long_description=read('README.md'),
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Topic :: Astronomy",
    ],
)
