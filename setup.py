from setuptools import setup

setup(name='unassign',
      version='0.0.1',
      description='Bacterial species unassigner',
      author='Kyle Bittinger',
      author_email='kylebittinger@gmail.com',
      url='https://github.com/kylebittinger/unassigner',
      packages=['unassign'],
      entry_points = {
          'console_scripts': [
              'unassignseqs=unassign.command:main',
              'prepare_strain_data=unassign.prepare_strain_data:main',
              'trimragged=unassign.trim:main',
          ],
      }
)
