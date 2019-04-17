import os
from setuptools import setup

PACKAGE_DIR = os.path.abspath(os.path.dirname(__file__))
README_FP = os.path.join(PACKAGE_DIR, 'README.md')
with open(README_FP, encoding='utf-8') as f:
    long_description = f.read()

setup(name='unassigner',
      version='0.0.4',
      description='Bacterial species unassigner',
      long_description=long_description,
      long_description_content_type='text/markdown',
      author='Kyle Bittinger',
      author_email='kylebittinger@gmail.com',
      url='https://github.com/kylebittinger/unassigner',
      packages=['unassigner'],
      install_requires=[
          'biopython',
          'scipy',
      ],
      entry_points={
          'console_scripts': [
              'unassign=unassigner.command:main',
              'trimragged=unassigner.trim:main',
          ],
      },
      classifiers=[
          'Intended Audience :: Science/Research',
          'License :: OSI Approved :: GNU General Public License v2 or later (GPLv2+)',
          'Operating System :: POSIX :: Linux',
          'Operating System :: MacOS :: MacOS X',
          'Programming Language :: Python :: 3',
          'Topic :: Scientific/Engineering :: Bio-Informatics',
      ],
      license='GPLv2+',
)
