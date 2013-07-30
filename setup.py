from distutils.core import setup

setup(
    name='taintedSwallow',
    version='0.0.1',
    author='Michael Imelfort',
    author_email='mike@mikeimelfort.com',
    packages=['taintedswallow', 'taintedswallow.test'],
    scripts=['bin/taintedSwallow'],
    url='http://pypi.python.org/pypi/taintedSwallow/',
    license='GPLv3',
    description='taintedSwallow',
    long_description=open('README.md').read(),
    install_requires=[],
)

