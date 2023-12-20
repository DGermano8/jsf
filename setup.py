from setuptools import setup, find_packages, Extension
from Cython.Build import cythonize

extensions = [
    Extension(
        "jsf.hello",  # Name of the Python module
        ["cython/hello.pyx", "src/hello.c"]
    )
]

setup(
    name='jsf',
    version='0.1',
    packages=find_packages(),
    description='Implementation of the jump-switch-flow simulator',
    ext_modules=cythonize(extensions),
    author='Domenic Germano',
    author_email='germanod@student.unimelb.edu.au',
    url='https://github.com/DGermano8/jazzy-shrill-fart',
    zip_safe=False,
)
