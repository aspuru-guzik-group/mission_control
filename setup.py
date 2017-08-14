import os
from setuptools import find_packages, setup

# allow setup.py to be run from any path
os.chdir(os.path.normpath(os.path.join(os.path.abspath(__file__), os.pardir)))

dependencies = [
    *[
        {'install_requires_value': simple_dependency}
        for simple_dependency in [
            'dill',
            'jinja2',
            'sqlalchemy',
        ]
    ],
    {
        'install_requires_value': 'jobman==0.0.1',
        'dependency_links_value': (
            'git+https://github.com/aspuru-guzik-group/jobman.git'
            '#egg=jobman-0.0.1'
        )
    },
]

setup(
    name='mission_control',
    version='0.0.1',
    packages=find_packages(),
    include_package_data=True,
    license='Apache 2.0',
    description='A workflow toolkit',
    author='A. Dorsk',
    classifiers=[
        'Intended Audience :: Developers',
        'Operating System :: OS Independent',
        'Programming Language :: Python',
        'Programming Language :: Python :: 3',
    ],
    install_requires=(
        [dependency['install_requires_value'] for dependency in dependencies
         if 'install_requires_value' in dependency]
    ),
    dependency_links=(
        [dependency['dependency_links_value'] for dependency in dependencies
         if 'dependency_links_value' in dependency]
    )
)
