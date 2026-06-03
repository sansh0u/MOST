from setuptools import setup, find_packages

setup(
    name="most",
    version="0.1.1",
    packages=find_packages(),
    include_package_data=True,
    package_data={"most": ["barcode/*.txt","config/*.yaml",]},
    install_requires=[
        "stpipeline>=2.1.0" 
    ],
    entry_points={
        "console_scripts": [
            "most = most.CLI:app",  #修改toolkit名字
        ],
    },
)