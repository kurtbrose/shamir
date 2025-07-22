from setuptools import setup, find_packages

setup(
    name='shamir',
    version='17.12.0',
    url='https://github.com/kurtbrose/shamir',
    author='Kurt Rose',
    author_email='kurt@kurtrose.com',
    decription="fast, secure, pure python shamir's secret sharing",
    long_description = open('README.rst').read(),
    py_modules=['shamir'],
    extras_require={'numpy': ['numpy']},
)