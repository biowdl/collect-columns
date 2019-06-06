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

import pytest

from collect_columns.collect_columns import main, parse_args


datadir = Path(__file__).parent / Path("data")


def test_main_htseq(tmpdir):
    sample1 = str(datadir / Path("htseq") / Path("sample1.fragments_per_gene"))
    sample2 = str(datadir / Path("htseq") / Path("sample2.fragments_per_gene"))
    expected_result = set([
        "feature\ts1\ts2\n",
        "MSTRG.1\t2371\t0\n",
        "MSTRG.2\t381\t1\n",
        "MSTRG.3\t741\t7\n",
        "MSTRG.4\t2361\t2\n",
        "MSTRG.5\t382\t3\n",
        "MSTRG.6\t706\t7\n",
        "__alignment_not_unique\t131\t13\n",
        "__ambiguous\t2995\t295\n",
        "__no_feature\t0\t0\n",
        "__not_aligned\t5\t51\n",
        "__too_low_aQual\t0\t0\n"])
    output_file = tmpdir.join("output.tsv")
    sys.argv = ["script", output_file.strpath, sample1, sample2, "-n", "s1",
                "s2"]
    main()
    result = set(output_file.open().readlines())
    assert result == expected_result


def test_main_stringtie(tmpdir):
    sample1 = str(datadir / Path("stringtie") / Path("sample1.abundance"))
    sample2 = str(datadir / Path("stringtie") / Path("sample2.abundance"))
    gtf = str(datadir / Path("merged.gtf"))
    expected_result = set([
        "feature\tref_gene_id\tgene_name\tsample1.abundance\t"
        "sample2.abundance\n",
        "MSTRG.1\tg_1;g_7\tgene_1;gene_7\t185151.953125\t85151.953125\n",
        "MSTRG.2\tg_2\tgene_2\t100160.070312\t160.070312\n",
        "MSTRG.3\tg_3\tgene_3\t91229.078125\t1229.078125\n",
        "MSTRG.4\tg_4\tgene_4\t184648.109375\t84648.109375\n",
        "MSTRG.5\tg_5\tgene_5\t104290.078125\t4290.078125\n",
        "MSTRG.6\tg_6\tgene_6\t89926.898438\t9926.898438\n"])

    output_file = tmpdir.join("output.tsv")
    sys.argv = ["script", output_file.strpath, sample1, sample2, "-c", "7",
                "-H", "-g", gtf, "-a", "ref_gene_id", "gene_name"]
    main()
    result = set(output_file.open().readlines())
    assert result == expected_result


def test_main_semicolon(tmpdir):
    sample1 = str(datadir / Path("semicolon") / Path("sample1.csv"))
    sample2 = str(datadir / Path("semicolon") / Path("sample2.csv"))
    gtf = str(datadir / Path("merged.gtf"))
    expected_result = set([
        "feature;ref_gene_id;transcript_id;sample1.csv;sample2.csv\n",
        "gene_1;g_1;t_1_2;1;10\n",
        "gene_2;g_2;t_2_1;2;20\n",
        "gene_3;g_3;t_3_1;3;30\n",
        "gene_4;g_4;\"t_4_2;t_4_1\";4;40\n",
        "gene_5;g_5;t_5_1;5;\n",
        "gene_6;g_6;t_6_1;;60\n"])

    output_file = tmpdir.join("output.tsv")
    sys.argv = ["script", output_file.strpath, sample1, sample2, "-f", "1",
                "-c", "0", "-s", ";", "-H", "-g", gtf, "-a", "ref_gene_id",
                "transcript_id", "-F", "gene_name"]
    main()
    result = set(output_file.open().readlines())
    assert result == expected_result


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
    assert args.value_column == 1
    assert args.sep == "\t"
    assert args.names is None
    assert args.header is False
    assert args.gtf is None
    assert args.additional_attributes is None
    assert args.feature_attribute == "gene_id"


def test_collect_columns_incorrect_number_of_names():
    sample1 = str(datadir / Path("semicolon") / Path("sample1.csv"))
    sample2 = str(datadir / Path("semicolon") / Path("sample2.csv"))
    sys.argv = ["script", "output", sample1, sample2, "-n", "sample1"]
    with pytest.raises(ValueError, match="The number of names did not match "
                                         "the number of inputs."):
        main()
