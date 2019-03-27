# Copyright (c) 2019 Leiden University Medical Center
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.


from pathlib import Path
import sys

import numpy as np
import pandas as pd
import pytest

from mergecounts import main, parse_args


datadir = Path(__file__).parent / Path("data")


def test_main_htseq(tmpdir):
    sample1 = str(datadir / Path("htseq") / Path("sample1.fragments_per_gene"))
    sample2 = str(datadir / Path("htseq") / Path("sample2.fragments_per_gene"))
    expected_result = pd.DataFrame(data={
        "feature": ["MSTRG.1", "MSTRG.2", "MSTRG.3", "MSTRG.4", "MSTRG.5",
                    "MSTRG.6", "__no_feature", "__ambiguous",
                    "__too_low_aQual", "__not_aligned",
                    "__alignment_not_unique"],
        "s1": [2371, 381, 741, 2361, 382, 706, 0, 2995, 0, 5, 131],
        "s2": [0, 1, 7, 2, 3, 7, 0, 295, 0, 51, 13]
    }, columns=["feature", "s1", "s2"])
    output_file = tmpdir.join("output.tsv")
    sys.argv = ["script", output_file.strpath, sample1, sample2, "-n", "s1",
                "s2"]
    main()
    result = pd.read_csv(output_file, sep="\t")
    assert result.equals(expected_result)


def test_main_stringtie(tmpdir):
    sample1 = str(datadir / Path("stringtie") / Path("sample1.abundance"))
    sample2 = str(datadir / Path("stringtie") / Path("sample2.abundance"))
    gtf = str(datadir / Path("merged.gtf"))
    expected_result = pd.DataFrame(data={
        "feature": ["MSTRG.1", "MSTRG.2", "MSTRG.3", "MSTRG.4", "MSTRG.5",
                    "MSTRG.6"],
        "ref_gene_id": ["g_1;g_7", "g_2", "g_3", "g_4", "g_5", "g_6"],
        "gene_name": ["gene_1;gene_7", "gene_2", "gene_3", "gene_4", "gene_5",
                      "gene_6"],
        "sample1.abundance": [185151.953125, 100160.070312, 91229.078125,
                              184648.109375, 104290.078125, 89926.898438],
        "sample2.abundance": [85151.953125, 160.070312, 1229.078125,
                              84648.109375, 4290.078125, 9926.898438],
    }, columns=["feature", "ref_gene_id", "gene_name", "sample1.abundance",
                "sample2.abundance"])
    output_file = tmpdir.join("output.tsv")
    sys.argv = ["script", output_file.strpath, sample1, sample2, "-c", "7",
                "-H", "-g", gtf, "-a", "ref_gene_id", "gene_name"]
    main()
    result = pd.read_csv(output_file, sep="\t")
    assert result.equals(expected_result)


def test_main_semicolon(tmpdir):
    sample1 = str(datadir / Path("semicolon") / Path("sample1.csv"))
    sample2 = str(datadir / Path("semicolon") / Path("sample2.csv"))
    gtf = str(datadir / Path("merged.gtf"))
    expected_result = pd.DataFrame(data={
        "feature": ["gene_1", "gene_2", "gene_3", "gene_4", "gene_5",
                    "gene_6"],
        "ref_gene_id": ["g_1", "g_2", "g_3", "g_4", "g_5", "g_6"],
        "gene_id": ["MSTRG.1", "MSTRG.2", "MSTRG.3", "MSTRG.4", "MSTRG.5",
                    "MSTRG.6"],
        "sample1.csv": [1, 2, 3, 4, 5, np.nan],
        "sample2.csv": [10, 20, 30, 40, np.nan, 60],
    }, columns=["feature", "ref_gene_id", "gene_id", "sample1.csv",
                "sample2.csv"])
    output_file = tmpdir.join("output.tsv")
    sys.argv = ["script", output_file.strpath, sample1, sample2, "-f", "1",
                "-c", "0", "-s", ";", "-H", "-g", gtf, "-a", "ref_gene_id",
                "gene_id", "-F", "gene_name"]
    main()
    result = pd.read_csv(output_file, sep=";")
    assert result.equals(expected_result)


def test_parse_args_a_but_no_g(capsys):
    with pytest.raises(SystemExit):
        sys.argv = ["script", "output", "input", "-a", "attribute"]
        parse_args()
    stderr = capsys.readouterr().err
    assert (
        "the following argument is required if -a is specified: -g" ==
        stderr[-58:-1])


def test_parse_args_defaults():
    sys.argv = ["script", "output", "input"]
    args = parse_args()
    assert args.table == [Path("input")]
    assert args.output == Path("output")
    assert args.feature_column == 0
    assert args.counts_column == 1
    assert args.sep == "\t"
    assert args.names is None
    assert args.header is False
    assert args.gtf is None
    assert args.additional_attributes is None
    assert args.feature_attribute == "gene_id"
