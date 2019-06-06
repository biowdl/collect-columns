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

from collect_columns.collect_columns import add_additional_attributes


datadir = Path(__file__).parent / Path("data")


def test_add_additional_attributes():
    gtf = datadir / Path("merged.gtf")
    count_table = {
        "MSTRG.1": {"sample1": "1", "sample2": "10"},
        "MSTRG.2": {"sample1": "2", "sample2": "20"},
        "MSTRG.3": {"sample1": "3", "sample2": "30"},
        "MSTRG.4": {"sample1": "4", "sample2": "40"},
        "MSTRG.5": {"sample1": "5"},
        "MSTRG.6": {"sample2": "60"},
        "MSTRG.7": {"sample1": "7"}}

    expected_result = {
        "MSTRG.1": {"ref_gene_id": "g_1;g_7", "gene_name": "gene_1;gene_7",
                    "sample1": "1", "sample2": "10"},
        "MSTRG.2": {"ref_gene_id": "g_2", "gene_name": "gene_2",
                    "sample1": "2", "sample2": "20"},
        "MSTRG.3": {"ref_gene_id": "g_3", "gene_name": "gene_3",
                    "sample1": "3", "sample2": "30"},
        "MSTRG.4": {"ref_gene_id": "g_4", "gene_name": "gene_4",
                    "sample1": "4", "sample2": "40"},
        "MSTRG.5": {"ref_gene_id": "g_5", "gene_name": "gene_5",
                    "sample1": "5"},
        "MSTRG.6": {"ref_gene_id": "g_6", "gene_name": "gene_6",
                    "sample2": "60"},
        "MSTRG.7": {"sample1": "7"}}
    result = add_additional_attributes(count_table, gtf, "gene_id",
                                       ["ref_gene_id", "gene_name"])
    assert result == expected_result
