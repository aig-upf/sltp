# coding=utf-8
from setuptools import setup, find_packages
from distutils.core import setup
from codecs import open
from os import path

here = path.abspath(path.dirname(__file__))
fsdir = path.abspath(path.join(here, '..', 'src'))
vendordir = path.abspath(path.join(here, '..', 'vendor'))
fs_builddir = path.abspath(path.join(here, '..', '.build'))


# Get the long description from the README file
with open(path.join(here, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()


def main():
    setup(
        name='sltp',
        version='0.1.0',
        description='sltp: The Sample+Learn+Transform+Plan Framework',
        long_description=long_description,
        url='https://github.com/aig-upf/features-generalized-planning',
        author='Blai Bonet and Guillem Francès',
        author_email='-',

        keywords='planning logic STRIPS',
        classifiers=[
            'Development Status :: 3 - Alpha',

            'Intended Audience :: Science/Research',
            'Topic :: Scientific/Engineering :: Artificial Intelligence',

            'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',

            'Programming Language :: Python :: 3',
            'Programming Language :: Python :: 3.4',
            'Programming Language :: Python :: 3.5',
            'Programming Language :: Python :: 3.6',
        ],

        packages=find_packages('src'),  # include all packages under src
        package_dir={'': 'src'},  # tell distutils packages are under src

        install_requires=[
            'setuptools',
            'tarski==0.1.0',
            'psutil',
            'bitarray',
            'numpy',
            'jsonpickle'
            # "pip @ git+ssh://git@github.com/aig-upf/tarski.git@cb92626051e79dca0439acd6b76a877adf30d497"
        ],

        # Git dependencies - egg names need to match whatever written in install_requires
        dependency_links=[
            'git+ssh://git@github.com/aig-upf/tarski.git@cb92626051e79dca0439acd6b76a877adf30d497#egg=tarski-0.1.0',
        ],


        extras_require={
            'dev': ['pytest', 'tox'],
            'test': ['pytest', 'tox'],
        },
    )


if __name__ == '__main__':
    main()
