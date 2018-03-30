# _*_ coding: utf-8 _*_

from setuptools import find_packages, setup
from codecs import open
from os import path

name = "dk_met_base"
author = __import__(name).__author__
version = __import__(name).__version__

here = path.abspath(path.dirname(__file__))

# Get the long description from the README file
with open(path.join(here, 'README.rst'), encoding='utf-8') as f:
    long_description = f.read()

setup(
    name=name,

    version=version,

    description='Collection of basic functions for array, mathematics, utility and so on.',
    long_description=long_description,

    # project home URL and download URL
    url='https://bitbucket.org/dave2017/'+name,
    download_url='https://bitbucket.org/dave2017/'+name+'/downloads/',

    # author
    author=author,
    author_email='daikan1998@163.com',

    # LICENSE
    license='MIT',

    classifiers=[
      'Development Status :: 3 - Alpha',
      'Intended Audience :: Developers',
      'Programming Language :: Python :: 3',
    ],

    keywords='meteorology data handling',

    packages=find_packages(exclude=['docs', 'tests', 'build', 'dist']),
    include_package_data=True,
    exclude_package_data={'': ['.gitignore']},

    install_requires=['numpy>=1.12.1',
                      'scipy>=0.19.0',
                      'pyproj>=1.9.5.1',
                      'python-dateutil',
                      'pandas>=0.20.0',
                      'numba>=0.36']
)

# development mode (DOS command):
#     setup fortran compiler environment.
#     python setup.py develop
#     python setup.py develop --uninstall

# build mode：
#     python setup.py build --build-base=D:/test/python/build

# distribution mode:
#     python setup.py sdist             # create source tar.gz file in /dist
#     python setup.py bdist_wheel       # create wheel binary in /dist
