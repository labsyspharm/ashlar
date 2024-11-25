import os
from urllib.request import urlopen
import hashlib
from setuptools import setup, find_packages
from setuptools.command.develop import develop
from setuptools.command.sdist import sdist
from setuptools.command.build_py import build_py
import versioneer

requires = [
    'numpy>=1.18.1',
    'pyjnius>=1.2.1',
    'matplotlib>=3.1.2',
    'networkx>=2.4',
    'scipy>=1.4.1',
    'scikit-image>=0.19.2',
    'scikit-learn>=0.21.1',
    'tifffile>=2023.3.15',
    'zarr>=2.11.3',
    'blessed>=1.17',
]


VERSION = versioneer.get_version()
DESCRIPTION = ('Alignment by Simultaneous Harmonization of Layer/Adjacency '
               'Registration')
LONG_DESCRIPTION='''

ASHLAR: Alignment by Simultaneous Harmonization of Layer/Adjacency Registration

Ashlar implements efficient combined stitching and registration of multi-channel
image mosaics collected using the Tissue-CycIF microscopy protocol [1]_. Although
originally developed for CycIF, it may also be applicable to other tiled and/or
cyclic imaging approaches. The package offers both a command line script for the
most common use cases as well as an API for building more specialized tools.

.. [1] Tissue-CycIF is multi-round immunofluorescence microscopy on large fixed
   tissue samples. See https://doi.org/10.1101/151738 for details.

'''
AUTHOR = 'Jeremy Muhlich'
AUTHOR_EMAIL = 'jeremy_muhlich@hms.harvard.edu'
LICENSE = 'MIT License'
HOMEPAGE = 'https://github.com/sorgerlab/ashlar'

BIOFORMATS_JAR_URL = 'https://downloads.openmicroscopy.org/bio-formats/8.0.1/artifacts/bioformats_package.jar'
BIOFORMATS_JAR_SHA256 = '8c7557a9357a83bf40272292fbd676beb466a9a8bab34126e92a49d636c64bc2'

def download_bioformats():
    print("Ensuring latest bioformats is present:")
    dist_root = os.path.abspath(os.path.dirname(__file__))
    jar_dir = os.path.join(dist_root, 'ashlar', 'jars')
    lt_jar_path = os.path.join(jar_dir, 'bioformats_package.jar')
    if not os.access(jar_dir, os.F_OK):
        os.mkdir(jar_dir)
    try:
        with open(lt_jar_path, 'rb') as f:
            existing_sha256 = hashlib.sha256(f.read()).hexdigest()
            if existing_sha256 == BIOFORMATS_JAR_SHA256:
                print("    Up to date!")
                return
    except IOError:
        pass
    print("    Downloading BioFormats from %s ..." % BIOFORMATS_JAR_URL)
    # FIXME add progress bar
    content = urlopen(BIOFORMATS_JAR_URL).read()
    content_sha256 = hashlib.sha256(content).hexdigest()
    with open(lt_jar_path, 'wb') as f:
        f.write(content)
    if content_sha256 != BIOFORMATS_JAR_SHA256:
        raise RuntimeError("bioformats_package.jar hash mismatch")

# Define some distutils command subclasses for a few key commands to trigger
# downloading the BioFormats JAR before they run.

class PreDevelop(develop):
    def run(self):
        download_bioformats()
        develop.run(self)

class PreSdist(sdist):
    def run(self):
        download_bioformats()
        sdist.run(self)

class PreBuildPy(build_py):
    def run(self):
        download_bioformats()
        build_py.run(self)

cmdclass = {
    'develop': PreDevelop,
    'sdist': PreSdist,
    'build_py': PreBuildPy,
}

setup(
    name='ashlar',
    version=VERSION,
    description=DESCRIPTION,
    long_description=LONG_DESCRIPTION,
    long_description_content_type='text/x-rst',
    cmdclass=versioneer.get_cmdclass(cmdclass),
    packages=find_packages(),
    include_package_data=True,
    install_requires=requires,
    entry_points={
        'console_scripts': [
            'ashlar=ashlar.scripts.ashlar:main',
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
        'Programming Language :: Python :: 3',
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
