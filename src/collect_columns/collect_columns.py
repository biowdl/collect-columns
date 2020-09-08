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
import csv
from pathlib import Path
from typing import List
from warnings import warn

import gffutils


def collect_columns(count_tables: List[Path], feature_column: int,
                    value_column: int, sep: str, names: List[str],
                    tables_have_headers: bool,
                    sum_on_duplicate_id: bool) -> dict:
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
    :param sum_on_duplicate_id: Whether or not values should be added up
    if multiple rows exist with the same feature id.
    """
    # If no sample names are given use table basenames instead.
    merged_table = {}
    for i, table in enumerate(count_tables):
        reader = csv.reader(table.open(), delimiter=sep)
        column_name = names[i]
        if tables_have_headers is True:
            next(reader)
        for record in reader:
            feature = record[feature_column]
            value = (float(record[value_column]) if sum_on_duplicate_id
                     else record[value_column])
            try:
                if column_name in merged_table[feature].keys():
                    if sum_on_duplicate_id:
                        merged_table[feature][column_name] += value
                    else:
                        warn("duplicate value for row {} in {}, will "
                             "overwrite previous value".format(feature,
                                                               column_name))
                        merged_table[feature][column_name] = value
                else:
                    merged_table[feature][column_name] = value
            # If feature key does not have a dict yet, create the dict.
            except KeyError:
                merged_table[feature] = {column_name: value}

        if sum_on_duplicate_id:
            for feature_id in merged_table.keys():
                try:
                    merged_table[feature_id][column_name] = str(
                        merged_table[feature_id][column_name])
                except KeyError:
                    pass
    return merged_table


def add_additional_attributes(table: dict, gtf: Path,
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
    with gtf.open("r") as in_file:
        for line in in_file:
            record = gffutils.feature.feature_from_line(line)
            for attr in additional_attributes:
                for feature in record.attributes.get(feature_attribute, []):
                    if feature in table.keys():
                        # Use lists to ensure the order stays the same.
                        # This way attributes which belong together
                        # will likely have to same position in their
                        # respective columns, assuming the available
                        # attributes for each record is consistent.
                        try:
                            table[feature][attr] += (
                                [a for a in record.attributes.get(attr, [])
                                 if a not in table[feature][attr]])
                        except KeyError:
                            table[feature][attr] = (
                                [a for a in record.attributes.get(attr, [])])
    # Turn lists into strings
    for feature in table.keys():
        for attribute in additional_attributes:
            try:
                table[feature][attribute] = ";".join(table[feature][attribute])
            except KeyError:
                pass
    return table


def main():
    args = parse_args()
    if args.names is None:
        names = [path.name for path in args.table]
    elif len(args.names) == len(args.table):
        names = args.names
    else:
        raise ValueError(
            "The number of names did not match the number of inputs.")

    merged_table = collect_columns(args.table, args.feature_column,
                                   args.value_column, args.sep, names,
                                   args.header, args.sum_on_duplicate_id)

    if args.additional_attributes is not None:
        merged_table = add_additional_attributes(
            merged_table, args.gtf, args.feature_attribute,
            args.additional_attributes)
        names = args.additional_attributes + names

    with args.output.open("w", newline="") as output_file:
        writer = csv.writer(output_file, delimiter=args.sep)
        writer.writerow(["feature"] + names)
        for feature, values in merged_table.items():
            row = [feature] + [values.get(column_name)
                               for column_name in names]
            writer.writerow(row)


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
    parser.add_argument("-S", "--sum-on-duplicate-id", action="store_true",
                        help="Whether or not values should be added up if "
                             "multiple rows exist with the same feature id. "
                             "The values will become floats if this flag is "
                             "set. By default only the last value will be "
                             "taken and a warning will be give.")
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
