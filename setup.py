import os
import pathlib
import sys

from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext as build_ext_orig

# Based on https://stackoverflow.com/a/48015772

class CMakeExtension(Extension):

    def __init__(self, name):
        # don't invoke the original build_ext for this special extension
        super().__init__(name, sources=[])


class build_ext(build_ext_orig):

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
            # Python_EXECUTABLE tells CMake's FindPython which Python executable to use.
            # Required for cibuildwheel to work.
            '-DPython_EXECUTABLE=' + sys.executable,
            '-DCMAKE_LIBRARY_OUTPUT_DIRECTORY=' + str(extdir.parent.absolute()),
            '-DCMAKE_BUILD_TYPE=' + config,
            '-DREQUIRE_PYTHON_WRAPPER=1',
            '-DENABLE_EXECUTABLES=0'
        ]

        # example of build args
        build_args = [
            '--', '-j2'
        ]

        os.chdir(str(build_temp))
        self.spawn(['cmake', str(cwd)] + cmake_args)
        if not self.dry_run:
            self.spawn(['cmake', '--build', '.'] + build_args)
        # Troubleshooting: if fail on line above then delete all possible
        # temporary CMake files including "CMakeCache.txt" in top level dir.
        os.chdir(str(cwd))


with open('readme', 'r') as fh:
    long_description = fh.read()

setup(
    name='varikn',
    setup_requires=['setuptools_scm>=3.4'],
    use_scm_version={
        'write_to': 'python-wrapper/varikn/_version.py'
    },
    author='Vesa Siivola',
    url='https://github.com/vsiivola/variKN',
    license='BSD 3-Clause License',
    description='VariKN language modeling toolkit',
    long_description=long_description,
    packages=['varikn'],
    package_dir={'': 'python-wrapper'},
    cmdclass={
        'build_ext': build_ext,
    },
    ext_modules=[CMakeExtension('varikn/wrapper')],
)
