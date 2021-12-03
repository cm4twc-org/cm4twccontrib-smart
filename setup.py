from setuptools import setup, find_namespace_packages


pkg_name = 'SMART'

with open("README.rst", 'r') as fh:
    long_desc = fh.read()

with open("unifhycontrib/{}/version.py".format(pkg_name.lower()), 'r') as fv:
    exec(fv.read())


setup(
    name='unifhycontrib-{}'.format(pkg_name.lower()),

    version=__version__,

    description='unifhy components for the {} model'.format(pkg_name),
    long_description=long_desc,
    long_description_content_type="text/x-rst",

    author='Thibault Hallouin',

    project_urls={
        'Source Code':
            'https://github.com/unifhy-org/unifhycontrib-{}'.format(pkg_name.lower())
    },

    license='GPL-3',

    classifiers=[
        'Development Status :: 4 - Beta',
        'Natural Language :: English',
        'Topic :: Scientific/Engineering :: Hydrology',
    ],

    packages=find_namespace_packages(include=['unifhycontrib.*'],
                                     exclude=['tests']),

    namespace_packages=['unifhycontrib'],

    install_requires=[
        'unifhy'
    ],

    extras_require={
        'tests': [
            'cf-python'
        ]
    },

    zip_safe=False
)
