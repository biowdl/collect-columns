from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from mergecounts.mergecounts import mergecounts


datadir = Path(__file__).parent / Path("data")


def test_mergecounts_htseq():
    tables = [
        datadir / Path("htseq") / Path("sample1.fragments_per_gene").__str__(),
        datadir / Path("htseq") / Path("sample2.fragments_per_gene").__str__()
    ]
    expected_result = pd.DataFrame(data={
        "feature": ["MSTRG.1", "MSTRG.2", "MSTRG.3", "MSTRG.4", "MSTRG.5",
                    "MSTRG.6", "__no_feature", "__ambiguous",
                    "__too_low_aQual", "__not_aligned",
                    "__alignment_not_unique"],
        "sample1.fragments_per_gene": [2371, 381, 741, 2361, 382, 706, 0,
                                       2995, 0, 5, 131],
        "sample2.fragments_per_gene": [0, 1, 7, 2, 3, 7, 0, 295, 0, 51, 13]
    }).set_index("feature")
    result = mergecounts(tables, 0, 1, "\t", None, False)
    assert result.equals(expected_result)


def test_mergecounts_htseq_with_names():
    tables = [
        datadir / Path("htseq") / Path("sample1.fragments_per_gene").__str__(),
        datadir / Path("htseq") / Path("sample2.fragments_per_gene").__str__()
    ]
    expected_result = pd.DataFrame(data={
        "feature": ["MSTRG.1", "MSTRG.2", "MSTRG.3", "MSTRG.4", "MSTRG.5",
                    "MSTRG.6", "__no_feature", "__ambiguous",
                    "__too_low_aQual", "__not_aligned",
                    "__alignment_not_unique"],
        "sample1": [2371, 381, 741, 2361, 382, 706, 0, 2995, 0, 5, 131],
        "sample2": [0, 1, 7, 2, 3, 7, 0, 295, 0, 51, 13]
    }).set_index("feature")
    result = mergecounts(tables, 0, 1, "\t", ["sample1", "sample2"], False)
    assert result.equals(expected_result)


def test_mergecounts_stringtie():
    tables = [
        datadir / Path("stringtie") / Path("sample1.abundance").__str__(),
        datadir / Path("stringtie") / Path("sample2.abundance").__str__()
    ]
    expected_result = pd.DataFrame(data={
        "feature": ["MSTRG.1", "MSTRG.2", "MSTRG.3", "MSTRG.4", "MSTRG.5",
                    "MSTRG.6"],
        "sample1.abundance": [185151.953125, 100160.070312, 91229.078125,
                              184648.109375, 104290.078125, 89926.898438],
        "sample2.abundance": [85151.953125, 160.070312, 1229.078125,
                              84648.109375, 4290.078125, 9926.898438],
    }).set_index("feature")
    result = mergecounts(tables, 0, 7, "\t", None, True)
    assert result.equals(expected_result)


def test_mergecounts_stringtie_with_names():
    tables = [
        datadir / Path("stringtie") / Path("sample1.abundance").__str__(),
        datadir / Path("stringtie") / Path("sample2.abundance").__str__()
    ]
    expected_result = pd.DataFrame(data={
        "feature": ["MSTRG.1", "MSTRG.2", "MSTRG.3", "MSTRG.4", "MSTRG.5",
                    "MSTRG.6"],
        "sample1": [185151.953125, 100160.070312, 91229.078125, 184648.109375,
                    104290.078125, 89926.898438],
        "sample2": [85151.953125, 160.070312, 1229.078125, 84648.109375,
                    4290.078125, 9926.898438],
    }).set_index("feature")
    result = mergecounts(tables, 0, 7, "\t", ["sample1", "sample2"], True)
    assert result.equals(expected_result)


def test_mergecounts_semicolon():
    tables = [
        datadir / Path("semicolon") / Path("sample1.csv").__str__(),
        datadir / Path("semicolon") / Path("sample2.csv").__str__()
    ]
    expected_result = pd.DataFrame(data={
        "feature": ["gene_1", "gene_2", "gene_3", "gene_4", "gene_5",
                    "gene_6"],
        "sample1.csv": [1, 2, 3, 4, 5, np.nan],
        "sample2.csv": [10, 20, 30, 40, np.nan, 60],
    }).set_index("feature")
    result = mergecounts(tables, 1, 0, ";", None, True)
    assert result.equals(expected_result)
