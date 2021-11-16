import setuptools


setuptools.setup(
    author="Abin Abraham",
    author_email="abin.abraham[AT]vanderbilt[dot]edu",
    name="gsel_vec",
    license="MIT",
    description="Detect signatures of evolution from GWAS data.",
    version="v0.0.1",
    long_description="long descrip",
    url="https://github.com/abraham-abin13/gsel_vec",
    scripts=["gsel_vec/run_evo_sig.py"],
    packages=setuptools.find_packages(),
    install_requires=[],
    # include_package_data=True,
    # package_data={'': ['data/*.csv']},
    python_requires=">=3.6",
    classifiers=[
        "Programming Language :: Python",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.6",
    ],
)
