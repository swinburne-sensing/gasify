#!/usr/bin/env python

import os
import re

from setuptools import setup, find_packages


_RE_URL_DEPENDENCY = re.compile(r'^[^:\s]+://[^#]+#egg=(.+)$')


# Read properties from __init__.py
with open(os.path.join(os.path.dirname(__file__), 'experimentdata', '__init__.py')) as file_init:
    content_init = file_init.read()

    version = re.search("__version__ = '([^']+)'", content_init).group(1)

    author = re.search("__author__ = '([^']+)'", content_init).group(1)

    maintainer = re.search("__maintainer__ = '([^']+)'", content_init).group(1)
    maintainer_email = re.search("__email__ = '([^']+)'", content_init).group(1)

# Read requirements from file
with open('requirements.txt', 'r') as file_init:
    requirements = [line.strip() for line in file_init if not line.startswith('#')]

# Fix dependencies from git
for n, requirement in enumerate(requirements):
    match_requirement = _RE_URL_DEPENDENCY.match(requirement)

    if match_requirement is not None:
        requirements[n] = match_requirement.group(1) + ' @ ' + requirement


setup(
    name='experimentdata',
    version=version,
    description='A module for handling data from experiments, primarily focused upon gas sensing results',
    long_description=open('README.md').read(),
    author=author,
    maintainer=maintainer,
    maintainer_email=maintainer_email,
    packages=find_packages(),
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Operating System :: OS Independent',
        'Programming Language :: Python',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
        'Topic :: Software Development :: Libraries',
        'Topic :: Software Development :: Libraries :: Python Modules',
    ],
    license='GPLv3',
    platforms='any',
    install_requires=requirements
)
