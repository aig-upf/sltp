# coding=utf-8
from setuptools import setup, find_packages
from distutils.core import setup
from codecs import open
from os import path

here = path.abspath(path.dirname(__file__))


def get_description():
    # Get the long description from the README file
    with open(path.join(here, 'README.md'), encoding='utf-8') as f:
        return f.read()


def get_version():
    import sys
    sys.path.insert(0, path.join(here, "src", "sltp"))
    import version
    v = version.get_version()
    sys.path = sys.path[1:]
    return v


def main():
    setup(
        name='sltp',
        version=get_version(),
        description='The SLTP Generalized Planning Framework: Sample, Learn, Transform & Plan',
        long_description=get_description(),
        url='https://github.com/aig-upf/sltp',
        author='Blai Bonet and Guillem Francès',
        author_email='-',

        keywords='planning logic STRIPS generalized planning',
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
            'psutil',
            'bitarray',
            'numpy',
            "tarski @ git+ssh://git@github.com/aig-upf/tarski.git@c8c2e1b#egg=tarski-dev-0.1.0"
        ],

        extras_require={
            'dev': ['pytest', 'tox'],
            'test': ['pytest', 'tox'],
        },
    )


if __name__ == '__main__':
    main()
