from setuptools import setup, find_packages

setup(
    name="most",
    use_scm_version=True,
    setup_requires=["setuptools_scm"],
    packages=find_packages(),
    include_package_data=True,
    package_data={
        "most": ["barcode/*.txt", "config/*.yaml"]
    },
    install_requires=[],
    entry_points={
        "console_scripts": [
            "most = most.CLI:app",
        ],
    },
)