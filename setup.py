from setuptools import setup, find_packages
from most.version import __version__
setup(
    name="most", 
    version=__version__,
    packages=find_packages(),
    include_package_data=True,
    package_data={"most": ["barcode/*.txt","config/*.yaml",]},
    install_requires=[
        "stpipeline>=2.1.0", 
        "plotly==5.24.1",
        "kaleido==0.2.1"
    ],
    entry_points={
        "console_scripts": [
            "most = most.CLI:app",  #修改toolkit名字
        ],
    },
)