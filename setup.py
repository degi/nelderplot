import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="nelderplot",
    version="0.0.8",
    author="Degi Harja Asmara",
    author_email="degiharja@gmail.com",
    description="A computer-aided tool for nelder plot design",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/degi/nelder-plot-designer",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: OS Independent",
    ],
    include_package_data=True,
    zip_safe=True,
    install_requires=[
          "matplotlib", "simplekml", "utm"
      ],
)
