from setuptools import setup, find_packages

setup(
    name="most", 
    version="0.2.1",
    packages=find_packages(),
    include_package_data=True,
    package_data={"most": ["barcode/*.txt","config/*.yaml","workflow/*.smk"]},
    install_requires=[],
    entry_points={
        "console_scripts": [
            "most = most.CLI:app",  
        ],
    },
)