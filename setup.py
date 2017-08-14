try:
	from setuptools import setup
except ImportError:
	from distutils.core import setup
setup(
    name='SVanGogh',
    version='0.1',
    author='Danny Antaki',
    author_email='dantaki@ucsd.edu',
    packages=['svangogh'],
    scripts=['svangogh/svangogh'],
    url='https://github.com/dantaki/SVanGogh',
    license='LICENSE',
    long_description=open('README').read(),
    install_requires=[
        "Pillow",
	"pysam",
        "pybedtools",
	"scipy",
    	"tqdm"
    ],
)
