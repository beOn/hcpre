#!/usr/bin/env python

import os
from setuptools import setup

def get_readme():
    md_path = os.path.join(os.path.dirname(__file__), "README.md")
    txt_path = os.path.join(os.path.dirname(__file__), "README.txt")
    if os.path.exists(txt_path):
        d = open(txt_path).read()
    elif os.path.exists(md_path):
        d = open(md_path).read()
    else:
        d = ""
    return d

setup(name='hcpre',
      version='0.5.2',
      author='Ben Acland',
      author_email='benacland@gmail.com',
      description='Generalized launcher for human connectome project BOLD preprocessing',
      license='BSD',
      keywords='connectome preprocessing fmri bold',
      url='https://github.com/beOn/hcpre',
      setup_requires=['nibabel','traits','networkx','pydicom'],
      install_requires=['pydicom','nibabel','networkx','traits','nipype','configobj'],
      packages=['hcpre'],
      scripts=['hcpre/hcpipe.py'],
      package_data={'hcpre':['duke_siemens/*']},
      long_description=get_readme(),
      classifiers=[
        'Development Status :: 3 - Alpha',
        'Topic :: Scientific/Engineering :: Information Analysis',
        'License :: OSI Approved :: BSD License',
        'Intended Audience :: Science/Research',
      ],
)
