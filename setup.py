from setuptools import setup, find_packages
import os

version = {}
with open(os.path.join("mags_visualization", "version.py"), encoding="utf-8") as f:
    exec(f.read(), version)

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

requirement_path = os.path.join(os.path.dirname(__file__), "requirements.txt")
install_requires = [
    "pandas",
    "numpy",
    "matplotlib",
    "seaborn",
    "plotly<7",
    "kaleido",
    "networkx",
    ]

if os.path.isfile(requirement_path):
    with open(requirement_path, encoding="utf-8") as f:
        install_requires = [line.strip() for line in f if line.strip() and not line.startswith("#")]

setup(
    name="MAGs-visualization",
    version=version["__version__"],
    description="MAGs visualization plots",
    long_description=long_description,
    long_description_content_type="text/markdown",
    author="Alexandra Hottmann",
    author_email="alexandra.hottm@gmx.de",
    license="GPL-3.0",
    url="https://github.com/usegalaxy-eu/MAGs-visualization",

    packages=find_packages(include=["mags_visualization", "mags_visualization.*"]),

    install_requires=install_requires,

    entry_points={
        "console_scripts": [
            "mags-visualization=mags_visualization.cli:cli",
        ]
    },

    include_package_data=True,
    python_requires=">=3.9",
)