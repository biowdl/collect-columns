from pathlib import Path

import pandas as pd

import pytest

from mergecounts.mergecounts import add_additional_attributes


datadir = Path(__file__).parent / Path("data")


def test_add_additional_attributes():
    gtf = datadir / Path("merged.gtf").__str__()
    count_table = pd.DataFrame(data={
        "feature": ["MSTRG.1", "MSTRG.2", "MSTRG.3", "MSTRG.4", "MSTRG.5",
                    "MSTRG.6", "__no_feature", "__ambiguous",
                    "__too_low_aQual", "__not_aligned",
                    "__alignment_not_unique"],
        "sample1.fragments_per_gene": [2371, 381, 741, 2361, 382, 706, 0,
                                       2995, 0, 5, 131]
    }).set_index("feature")
    expected_result = pd.DataFrame(data={
        "feature": ["MSTRG.1", "MSTRG.2", "MSTRG.3", "MSTRG.4", "MSTRG.5",
                    "MSTRG.6", "__no_feature", "__ambiguous",
                    "__too_low_aQual", "__not_aligned",
                    "__alignment_not_unique"],
        "ref_gene_id": ["g_1;g_7", "g_2", "g_3", "g_4", "g_5", "g_6",
                        None, None, None, None, None],
        "gene_name": ["gene_1;gene_7", "gene_2", "gene_3", "gene_4",
                      "gene_5", "gene_6", None, None, None, None, None],
        "sample1.fragments_per_gene": [2371, 381, 741, 2361, 382, 706, 0,
                                       2995, 0, 5, 131]
    }).set_index("feature")
    result = add_additional_attributes(count_table, gtf, "gene_id",
                              ["ref_gene_id", "gene_name"])
    assert result.equals(expected_result)
