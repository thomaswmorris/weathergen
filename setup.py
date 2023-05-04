from os import path
import setuptools

here = path.abspath(path.dirname(__file__))

with open('README.rst', 'r') as fh:
    long_description = fh.read()

with open(path.join(here, "requirements.txt")) as requirements_file:
    # Parse requirements.txt, ignoring any commented-out lines.
    requirements = [line for line in requirements_file.read().splitlines() if not line.startswith("#")]

setuptools.setup(
    name='weathergen',
    version='2.3.0',
    description="Generates time-varying weather profiles using a synthesis of in-situ observations and satellite reanalysis estimates of meteorological parameters.",
    long_description=long_description,
    author="Thomas Morris",
    author_email='thomasmorris@princeton.edu',
    url='https://github.com/thomaswmorris/weathergen',
    python_requires='>=3.9.9',
    packages=setuptools.find_packages(exclude=['docs', 'tests']),
    include_package_data=True,
    package_data={
        'weathergen': []
    },
    install_requires=requirements,
    license="BSD (3-clause)",
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Natural Language :: English',
        'Programming Language :: Python :: 3',
    ],
)
