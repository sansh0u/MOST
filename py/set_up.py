from setuptools import setup, find_packages

setup(
    name="toolkit",
    version="1.0.0",
    packages=find_packages(),
    install_requires=[
        "typer",
        "pyyaml",
        "biopython",
        "chromap",
    ],
    entry_points={
        "console_scripts": [
            "toolkit = toolkit.cli:app",  #修改toolkit名字
        ],
    },
)