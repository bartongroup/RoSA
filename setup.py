#!/usr/bin/env python

from setuptools import setup, find_packages

setup(name='viewseq',
      version='0.1',
      description='RNA-Seq data summariser',
      keywords=['bioinformatics', 'biology', 'RNA-Seq'],
      url='',
      author='Kira Mourao',
      author_email='k.mourao@dundee.ac.uk',
      license='',
      packages=find_packages(),
      
      #package_data = {
        # Include files from websummary in the package data:
      #  '': ['websummary/*'],
      #},
      
      include_package_data=True,
      zip_safe=False,

      install_requires=[
          'pandas>=0.18.0',
          'numpy',
          'scipy',
          'pysam',
          'drmaa',
          'six',
          'pyyaml',
      ],
      
      entry_points={
        'console_scripts': [
            'viewseq = viewseq.summariser:main'
        ]
      },

      classifiers=[
            'Development Status :: 3 - Alpha',
            'Environment :: Console',
            'Environment :: Web Environment',
            'Intended Audience :: Science/Research',
            'Natural Language :: English',
            'Operating System :: MacOS :: MacOS X',
            'Operating System :: POSIX',
            'Operating System :: Unix',
            'Programming Language :: Python :: 2.7',
            'Programming Language :: JavaScript',
            'Topic :: Scientific/Engineering',
            'Topic :: Scientific/Engineering :: Bio-Informatics',
            'Topic :: Scientific/Engineering :: Visualization',
      ],
      )
