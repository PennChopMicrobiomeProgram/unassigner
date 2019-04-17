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
      entry_points = {
          'console_scripts': [
              'unassign=unassigner.command:main',
              'trimragged=unassigner.trim:main',
          ],
      }
)
