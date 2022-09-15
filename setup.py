
import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()


setuptools.setup(
    name='weathergen',
    version='2.0.1',
    description="Generates weather",
    long_description=long_description,
    author="Thomas Morris",
    author_email='thomasmorris@princeton.edu',
    url='https://github.com/thomaswmorris/weathergen',
    python_requires='>=3.6',
    packages=setuptools.find_packages(exclude=['docs', 'tests']),
    include_package_data=True,
    package_data={
        'weathergen': [
            # When adding files here, remember to update MANIFEST.in as well,
            # or else they will not be included in the distribution on PyPI!
            # 'path/to/data_file',
        ]
    },
    install_requires=['numpy', 'datetime'],
    license="BSD (3-clause)",
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Natural Language :: English',
        'Programming Language :: Python :: 3',
    ],
)