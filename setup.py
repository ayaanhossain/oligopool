from setuptools import setup
from os import path
from io import open

here = path.abspath(path.dirname(__file__))

# Get the long description from the README file
with open(path.join(here, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

# Arguments marked as "Required" below must be included for upload to PyPI.
# Fields marked as "Optional" may be commented out.

setup(
    name='oligopool',

    # Link: https://www.python.org/dev/peps/pep-0440/#version-scheme
    version='2024.10.24',

    description='Oligopool Calculator - Automated design and analysis of oligopool libraries',

    long_description=long_description,

    long_description_content_type='text/markdown',

    url='https://github.com/ayaanhossain/oligopool',

    author='Ayaan Hossain and Howard Salis',

    author_email='auh57@psu.edu, salis@psu.edu',  # Optional

    classifiers=[  # Optional
        # How mature is this project? Common values are
        #   3 - Alpha
        #   4 - Beta
        #   5 - Production/Stable
        'Development Status :: 4 - Beta',

        # Indicate who your project is intended for
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Topic :: Scientific/Engineering :: Chemistry',

        # Pick your license as you wish
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',

        # Specify the Python versions you support here. In particular, ensure
        # that you indicate whether you support Python 2, Python 3 or both.
        # These classifiers are *not* checked by 'pip install'. See instead
        # 'python_requires' below.
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
    ],
    keywords=' '.join([
        'synthetic',
        'computational',
        'biology',
        'nucleotide',
        'oligo',
        'pool',
        'calculator',
        'design',
        'analysis',
        'barcode',
        'primer',
        'spacer',
        'motif',
        'split',
        'pad',
        'assembly',
        'index',
        'pack',
        'scry',
        'classifier',
        'count',
        'acount',
        'xcount']),

    packages=['oligopool', 'oligopool.base'],

    package_dir={
        'oligopool': './oligopool'
    },

    python_requires='>=3.10, <4',

    install_requires=[
        'biopython>=1.84',
        'primer3-py>=2.0.3',
        'msgpack>=1.1.0',
        'pyfastx>=2.1.0',
        'edlib>=1.3.9.post1',
        'dinopy>=3.0.0',
        'bounter>=1.2.0',
        'parasail>=1.3.4',
        'nrpcalc>=1.7.0',
        'plyvel>=1.5.1',
        'numba>=0.60.0',
        'seaborn>=0.13.2',
    ],

    project_urls={  # Optional
        'Bug Reports': 'https://github.com/ayaanhossain/oligopool/issues',
        'Source'     : 'https://github.com/ayaanhossain/oligopool/tree/master/oligopool',
    },
)