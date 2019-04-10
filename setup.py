# coding=utf-8

import os
import pathlib

from setuptools import setup, find_packages, Extension
from setuptools.command.build_ext import build_ext as build_ext_orig

from distutils.core import setup


here = os.path.abspath(os.path.dirname(__file__))


def get_description():
    from codecs import open
    # Get the long description from the README file
    with open(os.path.join(here, 'README.md'), encoding='utf-8') as f:
        return f.read()


def get_version():
    import sys
    sys.path.insert(0, os.path.join(here, "src", "sltp"))
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
            "tarski @ git+ssh://git@github.com/aig-upf/tarski.git@bc3afb9#egg=tarski-dev-0.2.0"
        ],

#        ext_modules=[CMakeExtension('featuregen', os.path.join(here, "src", "features"))],
#        cmdclass={'build_ext': BuildExt, },

        extras_require={
            'dev': ['pytest', 'tox'],
            'test': ['pytest', 'tox'],
        },
    )


class BuildExt(build_ext_orig):
    """ A helper to build c++ code using CMake.
     @see https://stackoverflow.com/a/48015772 """

    def run(self):
        for ext in self.extensions:
            self.build_cmake(ext)
        super().run()

    def build_cmake(self, ext):
        cwd = pathlib.Path().absolute()

        # these dirs will be created in build_py, so if you don't have
        # any python sources to bundle, the dirs will be missing
        build_temp = pathlib.Path(self.build_temp)
        build_temp.mkdir(parents=True, exist_ok=True)
        extdir = pathlib.Path(self.get_ext_fullpath(ext.name))
        extdir.mkdir(parents=True, exist_ok=True)

        # example of cmake args
        config = 'Debug' if self.debug else 'Release'
        cmake_args = [
            # '-DCMAKE_RUNTIME_OUTPUT_DIRECTORY=' + str(extdir.parent.absolute()),
            '-DCMAKE_RUNTIME_OUTPUT_DIRECTORY={}'.format(ext.path),
            '-DCMAKE_BUILD_TYPE=' + config
        ]

        # example of build args
        build_args = [
            '--config', config,
            '--', '-j4'
        ]

        os.chdir(str(build_temp))
        self.spawn(['cmake', ext.path] + cmake_args)
        if not self.dry_run:
            self.spawn(['cmake', '--build', '.'] + build_args)
        os.chdir(str(cwd))


class CMakeExtension(Extension):
    def __init__(self, name, path):
        # don't invoke the original build_ext for this special extension
        super().__init__(name, sources=[])
        self.path = path


if __name__ == '__main__':
    main()
