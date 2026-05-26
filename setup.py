from setuptools import setup, find_packages

setup(
    name="toolkit",
    version="0.1.0",
    packages=find_packages(),
    install_requires=[
        "stpipeline"
    ],
    entry_points={
        "console_scripts": [
            "toolkit = toolkit.CLI:app",  #修改toolkit名字
        ],
    },
)