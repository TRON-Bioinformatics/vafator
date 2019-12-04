from setuptools import find_packages, setup
import vafator


VERSION = vafator.VERSION


# Build the Python package
setup(
    name='vafator',
    version=VERSION,
    packages=find_packages(exclude=["legacy"]),
    entry_points = {
        'console_scripts': [
            'vafator=vafator.command_line:run'
        ],
    },
    author='Pablo Riesgo Ferreiro',
    author_email='pablo.riesgoferreiro@tron-mainz.de',
    description='Annotate a VCF with AF, AD and DP from tumor and normal BAMs',
    requires=[],
    install_requires=["pandas", "pysam", "cyvcf2"],
    classifiers=[
        'Development Status :: 3 - Alpha',      # Chose either "3 - Alpha", "4 - Beta" or "5 - Production/Stable" as the current state of your package
        'Intended Audience :: Healthcare Industry',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Programming Language :: Python :: 3.7'
      ]
)