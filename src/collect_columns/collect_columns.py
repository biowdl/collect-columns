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

import argparse
from pathlib import Path
from typing import List, Optional

from BCBio import GFF
import pandas as pd


def collect_columns(count_tables: List[Path], feature_column: int,
                    value_column: int, sep: str, names: Optional[List[str]],
                    tables_have_headers: bool) -> pd.DataFrame:
    """
    Retrieve a column from each in a set of tables and put them into a
    single table, mapping the rows based on other column.
    :param count_tables: A list of paths to the tables to be merged.
    :param feature_column: The position of the column with the (unique)
    feature/row ids.
    :param value_column: The position of the column with the values of
    interest.
    :param sep: The separator used in the tables.
    :param names: The names of the samples corresponding to the tables
    (in the same order as the tables). If None the basenames of tables
    will be used.
    :param tables_have_headers: Whether or not the tables have a header.
    """
    # If no sample names are given use table basenames instead.
    if names is None:
        names = [path.name for path in count_tables]
    # Convert tables_have_headers to conform with the pandas parameter.
    header = 0 if tables_have_headers else None
    # Load in the first table.
    merged_table = pd.read_csv(count_tables.pop(0), sep=sep, header=header,
                               index_col=feature_column,
                               usecols=[feature_column, value_column])
    merged_table.columns = [names.pop(0)]
    merged_table.index.names = ["feature"]
    # Merge the other tables into the first.
    for i, table in enumerate(count_tables):
        sample_table = pd.read_csv(table, sep=sep, header=header,
                                   index_col=feature_column,
                                   usecols=[feature_column, value_column])
        sample_table.columns = [names[i]]
        sample_table.index.names = ["feature"]
        merged_table = merged_table.merge(sample_table, left_index=True,
                                          right_index=True, how="outer")
    return merged_table


def add_additional_attributes(table: pd.DataFrame, gtf: Path,
                              feature_attribute: str,
                              additional_attributes: List[str]):
    """
    Retrieve additional attributes from the GTF/GFF and add them to the
    table.
    :param table: The pandas DataFrame to which the additional
    attributes will be added.
    :param gtf: The path to the the GTF/GFF file.
    :param feature_attribute: The attribute from the GTF/GFF used for
    matching the feature records with the rows in the table.
    :param additional_attributes: A list containing the keys of the
    attributes which will be added to the table.
    """
    # Create dictionary mapping attributes to features
    features_list = {feature: {attr: [] for attr in additional_attributes}
                     for feature in table.index}  # dict of dicts of list
    with gtf.open("r") as in_file:
        for rec in GFF.parse_simple(in_file):
            record_attributes = rec['quals']
            for attr in additional_attributes:
                for feature in record_attributes.get(feature_attribute, []):
                    if feature in table.index:
                        # Use lists to ensure the order stays the same.
                        # This way attributes which belong together
                        # will likely have to same position in their
                        # respective columns, assuming the available
                        # attributes for each record is consistent.
                        features_list[feature][attr].extend(
                            [a for a in record_attributes.get(attr, [])
                             if a not in features_list[feature][attr]])

    for attr in additional_attributes[::-1]:
        column = [";".join(features_list[feature][attr])
                  if len(features_list[feature][attr]) > 0 else None
                  for feature in table.index]
        table.insert(0, attr, column)
    return table


def main():
    args = parse_args()
    merged_table = collect_columns(args.table, args.feature_column,
                                   args.value_column, args.sep, args.names,
                                   args.header)
    if args.additional_attributes is not None:
        merged_table = add_additional_attributes(
            merged_table, args.gtf, args.feature_attribute,
            args.additional_attributes)
    merged_table.to_csv(args.output, sep=args.sep)


def parse_args():
    parser = argparse.ArgumentParser(
        description="Retrieves a column from a set of tables and puts "
                    "them into a single table.\n"
                    "Optionally, additional attributes may be retrieved from "
                    "a GTF or GFF file, which will be added as additional "
                    "column in the merged table as well.")
    # positional
    parser.add_argument("output", type=Path,
                        help="The path the output will be written to.")
    parser.add_argument("table", type=Path, nargs="+",
                        help="The tables to be merged.")
    # optional
    parser.add_argument("-f", "--feature-column", type=int, default=0,
                        metavar="I",
                        help="The position of the column with the "
                             "(unique) feature ids. Default to 0.")
    parser.add_argument("-c", "--value-column", type=int, default=1,
                        metavar="I",
                        help="The position of the column with the "
                             "values of interest. Defaults to 1.")
    parser.add_argument("-s", "--sep", "--separator", type=str, default="\t",
                        help="The separator used in the tables. This "
                             "will also be used in the output table. "
                             "Defaults to a tab.")
    parser.add_argument("-n", "--names", type=str, nargs="*", metavar="NAME",
                        help="The names of the samples corresponding "
                             "to the tables (in the same order as the "
                             "tables). These will be used as headers in the "
                             "merged table. If not specified the basenames "
                             "of tables will be used.")
    parser.add_argument("-H", "--header", action='store_true',
                        help="Whether or not the tables have a header. "
                             "Defaults to false.")
    parser.add_argument("-a", "--additional-attributes", type=str, nargs="+",
                        metavar="ATTR",
                        help="A list of attributes which will be added "
                             "to the merged table. These attributes "
                             "will be retrieved from the GTF or GFF file "
                             "specified with the -g option. Multiple values "
                             "will be separator by a ';'. Requires -g to be "
                             "specified.")
    parser.add_argument("-g", "--gtf", "--gff", type=Path, metavar="FILE",
                        help="The GTF or GFF file from which the "
                             "additional attributes (see -a) will be "
                             "retrieved. Ignored if -a is not specified. "
                             "Required if -a is specified.")
    parser.add_argument("-F", "--feature-attribute", type=str, metavar="ATTR",
                        default="gene_id",
                        help="The attribute from the GTF/GFF used for "
                             "matching the feature records with the rows in "
                             "the table. Ignored if -a is not "
                             "specified. Defaults to 'gene_id'.")
    args = parser.parse_args()
    if args.additional_attributes is not None and args.gtf is None:
        parser.error("the following argument is required if -a is "
                     "specified: -g")
    return args


if __name__ == "__main__":
    main()
