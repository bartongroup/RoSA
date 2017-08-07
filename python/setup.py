#!/usr/bin/env python

#  Copyright 2017 Kira Mourao, Nick Schurch
#
#  This file is part of RoSA.
#
#  RoSA is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  RoSA is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with RoSA.  If not, see <http://www.gnu.org/licenses/>.

from setuptools import setup, find_packages

setup(name='rosa',
      version='0.1',
      description='RNA-Seq antisense data analyser',
      keywords=['bioinformatics', 'biology', 'RNA-Seq'],
      url='',
      author='Kira Mourao',
      author_email='k.mourao@dundee.ac.uk',
      license='GPL3',
      packages=find_packages(),
      
      include_package_data=True,
      zip_safe=False,

      install_requires=[
          'pandas>=0.18.0,<=0.19.2',
          'numpy',
          'scipy',  ###>=0.16.1,<=0.17.1',
          'drmaa',
          'six',
      ],
      
      entry_points={
        'console_scripts': [
            'make_annotation = rosa.make_annotation:main',
            'count_spliced = rosa.antisense:main',
            'antisense_job = rosa.anti_reads_as_job:main'
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
