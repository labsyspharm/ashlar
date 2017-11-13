import os
from setuptools import setup, find_packages
import versioneer

here = os.path.abspath(os.path.dirname(__file__))
with open(os.path.join(here, 'README.md')) as f:
    README = f.read()

requires = [
    'numpy>=1.13.0',
    'javabridge>=1.0.15',
    'python-bioformats==1.1.0',
    'matplotlib>=2.1.0',
    'ModestImage>=0.1',
    'networkx>=2.0',
    'pathlib2>=2.3.0',
    'pyfftw>=0.10.4',
    'scipy>=0.19.1',
    'scikit-image>=0.13.0',
    'scikit-learn>=0.19.1'
]


VERSION = versioneer.get_version()
DESCRIPTION = ('Alignment by Simultaneous Harmonization of Layer/Adjacency '
               'Registration')
AUTHOR = 'Jeremy Muhlich'
AUTHOR_EMAIL = 'jeremy_muhlich@hms.harvard.edu'
LICENSE = 'MIT License'
HOMEPAGE = 'https://github.com/sorgerlab/ashlar'

setup(
    name='ashlar',
    version=VERSION,
    description=DESCRIPTION,
    long_description=README,
    cmdclass=versioneer.get_cmdclass(),
    packages=find_packages(),
    include_package_data=True,
    install_requires=requires,
    entry_points={
        'console_scripts': [
            'mosaic=ashlar.scripts.mosaic:main',
            'preview_slide=ashlar.scripts.preview_slide:main',
            'make_alignment_movie=ashlar.scripts.make_alignment_movie:main'
        ]
    },
    classifiers=[
        'Development Status :: 5 - Production/Stable',
        'Intended Audience :: End Users/Desktop',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: %s' % LICENSE,
        'Natural Language :: English',
        'Operating System :: OS Independent',
        'Programming Language :: Python :: 2',
        'Topic :: Scientific/Engineering :: Visualization'
    ],
    author=AUTHOR,
    author_email=AUTHOR_EMAIL,
    license=LICENSE,
    url=HOMEPAGE,
    download_url='%s/archive/v%s.tar.gz' % (HOMEPAGE, VERSION),
    keywords=['scripts', 'microscopy', 'registration', 'stitching'],
    zip_safe=False,
)
