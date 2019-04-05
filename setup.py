# Copyright (c) 2019 Leiden University Medical Center
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

from setuptools import setup

with open("README.md", "r") as readme_file:
    LONG_DESCRIPTION = readme_file.read()

setup(name="collect-columns",
      version="0.1.1",
      description="Retrieve a columns for each each in a set of tables, "
                  "placing them in a single output table.",
      long_description=LONG_DESCRIPTION,
      long_description_content_type='text/markdown',
      classifiers=[
        "Programming Language :: Python :: 3 :: Only",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.5",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "Development Status :: 3 - Alpha",
        "License :: OSI Approved :: MIT License",
      ],
      keywords="bioinformatics",
      url="https://github.com/biowdl/collect-columns",
      author="Leiden University Medical Center",
      author_email="sasc@lumc.nl",
      license="MIT",
      packages=["collect_columns"],
      package_dir={'': 'src'},
      install_requires=[
        "pandas>=0.23",
        "bcbio-gff",
        "biopython" #Required for bcbio-gff
      ],
      entry_points={
          "console_scripts":
              ["collect-columns=collect_columns.collect_columns:main"]
      })
