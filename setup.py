#!/usr/env python
from distutils.core import setup
setup(
	name='svangogh',
	version='1.0',
	description='Pixelate SV breakpoints',
	author='Danny Antaki',
	author_email='dantaki@ucsd.edu',
	url='https://github.com/dantaki/SVanGogh',
	packages=['svangogh'],
	scripts=['svangogh/svangogh'],
	license='LICENSE.txt',
	install_requires=[
		"pysam",
		"pybedtools",
		"scipy"
	]
	
     )
