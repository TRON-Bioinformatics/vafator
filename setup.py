from setuptools import find_packages, setup
import vafator


VERSION = vafator.VERSION


# parses requirements from file
with open("requirements.txt") as f:
    required = f.read().splitlines()

with open("README.md", "r", encoding="utf-8") as f:
    long_description = f.read()

# Build the Python package
setup(
    name='vafator',
    version=VERSION,
    packages=find_packages(exclude=["legacy"]),
    entry_points={
        'console_scripts': [
            'vafator=vafator.command_line:annotator',
            'multiallelics-filter=vafator.command_line:multiallelics_filter',
            'vafator2decifer=vafator.command_line:vafator2decifer'
        ],
    },
    author="TRON - Translational Oncology at the University Medical Center of the Johannes Gutenberg University Mainz"
    "- Computational Medicine group",
    author_email='pablo.riesgoferreiro@tron-mainz.de',
    description='Annotate a VCF file with AF, AD and DP from tumor and normal BAMs',
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/tron-bioinformatics/vafator",
    requires=[],
    install_requires=required,
    classifiers=[
        'Development Status :: 4 - Beta',      # Chose either "3 - Alpha", "4 - Beta" or "5 - Production/Stable" as the current state of your package
        'Intended Audience :: Healthcare Industry',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
        'Programming Language :: Python :: 3 :: Only',
        "License :: OSI Approved :: MIT License",
        "Operating System :: Unix"
      ],
    python_requires='>=3.7',
    license='MIT'
)