import os
from setuptools import find_packages, setup

with open(os.path.join(os.path.dirname(__file__), 'README.rst')) as readme:
    README = readme.read()

# allow setup.py to be run from any path
os.chdir(os.path.normpath(os.path.join(os.path.abspath(__file__), os.pardir)))

setup(
    name='mission-control',
    version='0.1',
    packages=find_packages(),
    include_package_data=True,
    license='Apache 2.0',
    description='A workflow toolkit',
    long_description=README,
    author='A. Dorsk',
    classifiers=[
        'Intended Audience :: Developers',
        'Operating System :: OS Independent',
        'Programming Language :: Python',
        'Programming Language :: Python :: 3',
    ],
    install_requires=[
        'dill',
        'jinja2',
        'marshmallow',
        'sqlalchemy',
    ],
)
