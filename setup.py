from setuptools import find_packages, setup
from pathlib import Path

setup(
    name="SCOT",
    version="2.0",
    description="Gromov-Wasserstein optimal transport for aligning single-cell multi-omics data",
    author="Pinar Demetci",
    packages=find_packages(),
    python_requires='>=3.6',
    install_requires=[
        l.strip() for l in Path('requirements.txt').read_text('utf-8').splitlines()
    ],
    license="MIT",
)
