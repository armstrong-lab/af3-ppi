#!/usr/bin/env python3
"""
Setup script for af3ppi - AlphaFold3 Protein-Protein Interaction CLI Tool
"""

from setuptools import setup, find_packages
import os

# Read the README file for long description
def read_readme():
    readme_path = os.path.join(os.path.dirname(__file__), 'README.md')
    if os.path.exists(readme_path):
        with open(readme_path, 'r', encoding='utf-8') as f:
            return f.read()
    return ''

# Read requirements
def read_requirements():
    req_path = os.path.join(os.path.dirname(__file__), 'requirements.txt')
    if os.path.exists(req_path):
        with open(req_path, 'r', encoding='utf-8') as f:
            return [line.strip() for line in f if line.strip() and not line.startswith('#')]
    return []

setup(
    name="af3ppi",
    version="0.1.0",
    author="Pablo Riera Freire",
    author_email="pablor_freire@dfci.harvard.edu",  # Add your email if desired
    description="CLI tool for generating AlphaFold3 server JSON inputs and analyzing outputs",
    long_description=read_readme(),
    long_description_content_type="text/markdown",
    url="https://github.com/armstrong-lab/af3-ppi.git",  # Replace with your actual GitHub URL
    packages=find_packages(),
    include_package_data=True,
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: OS Independent",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Programming Language :: Python :: 3.12",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    python_requires=">=3.9",
    install_requires=read_requirements(),
    entry_points={
        "console_scripts": [
            "af3ppi=af3ppi.__main__:main",
        ],
    },
    keywords="alphafold3 protein-protein-interaction bioinformatics",
    project_urls={
        "Bug Reports": "https://github.com/yourusername/af3ppi/issues",  # Update URL
        "Source": "https://github.com/yourusername/af3ppi",  # Update URL
    },
)