from setuptools import setup, find_packages

setup(
    name="most",
    version="0.1.0",
    packages=find_packages(),
    package_dir= {},
    package_data={"": ["*.txt","*.yaml"]},
    use_scm_version=True,
    setup_requires=['setuptools_scm'],
    install_requires=[
        "stpipeline>=2.1.0" 
    ],
    entry_points={
        "console_scripts": [
            "most = most.CLI:app",  #修改toolkit名字
        ],
    },
)