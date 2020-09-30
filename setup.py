"""A setuptools based setup module.

See:
https://packaging.python.org/en/latest/distributing.html
https://github.com/pypa/sampleproject
"""

# Always prefer setuptools over distutils
from setuptools import setup, find_packages
# To use a consistent encoding
from codecs import open
from os import path
from xmipptomo import __version__

here = path.abspath(path.dirname(__file__))

# Get the long description from the README file
with open(path.join(here, 'README.rst'), encoding='utf-8') as f:
    long_description = f.read()

with open('requirements.txt') as f:
    requirements = f.read().splitlines()

# Arguments marked as "Required" below must be included for upload to PyPI.
# Fields marked as "Optional" may be commented out.

setup(
    name='scipion-em-xmipptomo',  # Required
    version=__version__,  # Required
    description='Scipion plugin to deal with xmipp tomography protocols.',  # Required

    long_description=long_description,  # Optional
    url='https://github.com/i2pc/scipion-em-xmipptomo',  # Optional
    author='I2PC',  # Optional
    author_email='scipion@cnb.csic.es',  # Optional
    keywords='scipion cryoem imageprocessing tomography scipion-3.0',  # Optional
    packages=find_packages(),
    install_requires=[requirements],
    entry_points={  # Optional
       'pyworkflow.plugin': 'xmipptomo = xmipptomo',
    },
    package_data={  # Optional
       'xmipptomo': ['xmipp_logo.png', 'protocols.conf'],
    }
)
