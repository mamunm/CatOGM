import setuptools

with open('readme.org') as f:
    readme = f.read()

with open('requirements.txt') as f:
    requirements = f.read()

setuptools.setup(
    name="CatOGM",
    version="0.0.1",
    url="https://github.com/mamunm/CatOGM",

    author="Osman Mamun",
    author_email="mamunm@stanford.edu",

    description="Fingerprint generator for atomic structures.",
    long_description=readme,

    license='None',

    packages=[
        'catogm',
        'catogm.fingerprint',
    ],

    install_requires=requirements,
    python_requires='<=3.5, <4',

    package_data={'catogm': ['data/*']},

    classifiers=[
        'Development Status :: 1 - Planning',
        'Intended Audience :: Developers',
        'Programming Language :: Python',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.6',
    ],
)
