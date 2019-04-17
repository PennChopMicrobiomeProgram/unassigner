from setuptools import setup

setup(name='unassigner',
      version='0.0.3',
      description='Bacterial species unassigner',
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
