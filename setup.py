#!/usr/bin/env python

from distutils.core import setup

setup(name='hc',
      version='0.1',
      description='Console-based RPN calculator',
      author='Vernon Mauery',
      author_email='vernon@mauery.com',
      url='',
      packages=['lhc'],
      package_dir={'lhc': 'lhc'},
      scripts=['hc'],
     )

