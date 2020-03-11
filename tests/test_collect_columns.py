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
from warnings import catch_warnings

from collect_columns.collect_columns import collect_columns


datadir = Path(__file__).parent / Path("data")


def test_collect_columns_htseq():
    tables = [datadir / Path("htseq") / Path("sample1.fragments_per_gene"),
              datadir / Path("htseq") / Path("sample2.fragments_per_gene")]
    expected_result = {
        "MSTRG.1": {"sample1": "2371", "sample2": "0"},
        "MSTRG.2": {"sample1": "381", "sample2": "1"},
        "MSTRG.3": {"sample1": "741", "sample2": "7"},
        "MSTRG.4": {"sample1": "2361", "sample2": "2"},
        "MSTRG.5": {"sample1": "382", "sample2": "3"},
        "MSTRG.6": {"sample1": "706", "sample2": "7"},
        "__alignment_not_unique": {"sample1": "131", "sample2": "13"},
        "__ambiguous": {"sample1": "2995", "sample2": "295"},
        "__no_feature": {"sample1": "0", "sample2": "0"},
        "__not_aligned": {"sample1": "5", "sample2": "51"},
        "__too_low_aQual": {"sample1": "0", "sample2": "0"}}
    result = collect_columns(tables, 0, 1, "\t", ["sample1", "sample2"], False)
    assert result == expected_result


def test_collect_columns_stringtie():
    tables = [datadir / Path("stringtie") / Path("sample1.abundance"),
              datadir / Path("stringtie") / Path("sample2.abundance")]
    expected_result = {
        "MSTRG.1": {"sample1": "185151.953125", "sample2": "85151.953125"},
        "MSTRG.2": {"sample1": "100160.070312", "sample2": "160.070312"},
        "MSTRG.3": {"sample1": "91229.078125", "sample2": "1229.078125"},
        "MSTRG.4": {"sample1": "184648.109375", "sample2": "84648.109375"},
        "MSTRG.5": {"sample1": "104290.078125", "sample2": "4290.078125"},
        "MSTRG.6": {"sample1": "89926.898438", "sample2": "9926.898438"}}
    with catch_warnings(record=True) as warnings:
        result = collect_columns(tables, 0, 7, "\t", ["sample1", "sample2"],
                                 True)
        assert "duplicate value for row MSTRG.6 in sample2, will overwrite " \
               "previous value" == str(warnings[0].message)
    assert result == expected_result


def test_collect_columns_semicolon():
    tables = [datadir / Path("semicolon") / Path("sample1.csv"),
              datadir / Path("semicolon") / Path("sample2.csv")]
    expected_result = {
        "gene_1": {"sample1": "1", "sample2": "10"},
        "gene_2": {"sample1": "2", "sample2": "20"},
        "gene_3": {"sample1": "3", "sample2": "30"},
        "gene_4": {"sample1": "4", "sample2": "40"},
        "gene_5": {"sample1": "5"},
        "gene_6": {"sample2": "60"}}
    result = collect_columns(tables, 1, 0, ";", ["sample1", "sample2"], True)
    assert result == expected_result
