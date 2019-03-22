import argparse
from os.path import basename
import sys

from BCBio import GFF
import pandas as pd


def mergecounts(count_tables, feature_column, counts_column, sep, names,
                tables_have_headers):
    """
    Merge a set of count tables.
    :param count_tables: A list of paths (as strings) to the count
    tables to be merged.
    :param feature_column: The position of the column with the (unique)
    feature ids.
    :param counts_column: The position of the column with the values of
    interest.
    :param sep: The seperator used in the tables.
    :param names: The names of the samples corresponsing to the tables
    (in the same order as the tables). If None the basenames of tables
    will be used.
    :param tables_have_headers: Whether or not the tables have a header.
    """
    # If no sample names are given use table basenames instead.
    if names is None:
        names = [basename(path) for path in count_tables]
    # Convert tables_have_headers to conform with the pandas parameter.
    header = 0 if tables_have_headers else None
    # Load in the first table.
    merged_table = pd.read_csv(count_tables.pop(0), sep=sep, header=header,
                               index_col=feature_column,
                               usecols=[feature_column, counts_column])
    merged_table.columns = [names.pop(0)]
    merged_table.index.names = ["feature"]
    # Merge the other tables into the first.
    for i, table in enumerate(count_tables):
        sample_table = pd.read_csv(table, sep=sep, header=header,
                                   index_col=feature_column,
                                   usecols=[feature_column, counts_column])
        sample_table.columns = [names[i]]
        sample_table.index.names = ["feature"]
        merged_table = merged_table.merge(sample_table, left_index=True,
                                          right_index=True, how="outer")
    return merged_table


def add_additional_attributes(counts_table, gtf, feature_attribute,
                              additional_attributes):
    """
    Retrieve additional attributes from the GTF/GFF and add them to the
    counts table.
    :param counts_table: The pandas Dataframe to which the additional
    attributes will be added.
    :param gtf: The path (as a string) to the the GTF/GFF file.
    :param feature_attribute: The attribute from the GTF/GFF used for
    matching the feature records with the rows in the count table.
    :param additional_attributes: A list containing the keys of the
    attributes which will be added to the table.
    """
    # Add columns for the additional attributes to the counts table.
    for attribute in additional_attributes:
        counts_table.insert(0, attribute, None)
    # For each record in the GTF/GFF
    in_file = open(gtf)
    for i, rec in enumerate(GFF.parse_simple(in_file)):
        attributes = rec["quals"]
        # If the attribute signifying the feature/gene is present
        if feature_attribute in attributes.keys():
            gene = attributes[feature_attribute][0]
            # and the gene is present in the count table
            if gene in counts_table.index:
                # add the additional attributes
                for attribute in additional_attributes:
                    # if the attribute exists in the gtf record
                    if attribute in attributes.keys():
                        if counts_table[attribute][gene] is None:
                            counts_table.at[gene, attribute] = []
                        for value in attributes[attribute]:
                            # and the value isn't present already.
                            if value not in counts_table[attribute][gene]:
                                counts_table.at[gene, attribute].append(value)
    in_file.close()
    # Reformat the retreived atributes into strings
    for attribute in additional_attributes:
        counts_table[attribute] = counts_table[attribute].apply(
            lambda x: None if x is None else ";".join(x))
    return counts_table


def main(sample_tables, output_path, feature_column, counts_column, sep,
         names, tables_have_headers, gtf, additional_attributes,
         feature_attribute):
    merged_counts_table = mergecounts(sample_tables, feature_column,
                                      counts_column, sep=sep, names=names,
                                      tables_have_headers=tables_have_headers)
    if additional_attributes is not None:
        merged_counts_table = add_additional_attributes(merged_counts_table,
                                                        gtf, feature_attribute,
                                                        additional_attributes)
    merged_counts_table.to_csv(output_path, sep=sep)


def parse_args(argv):
    parser = argparse.ArgumentParser(
        description="Merges counts tables from multiple samples into a "
                    "single counts table. Optionally, additional "
                    "attributes may be retreived from a GTF or GFF "
                    "file, which will be added as additional column in "
                    "the merged count table as well.")
    # positional
    parser.add_argument("output", type=str,
                        help="The path the output will be written to.")
    parser.add_argument("table", type=str, nargs="+",
                        help="The count tables to be merged.")
    # optional
    parser.add_argument("-f", "--feature-column", type=int, default=0,
                        metavar="I",
                        help="The position of the column with the "
                        "(unique) feature ids. Default to 0.")
    parser.add_argument("-c", "--counts-column", type=int, default=1,
                        metavar="I",
                        help="The position of the column with the "
                        "values of interest. Defaults to 1.")
    parser.add_argument("-s", "--sep", "--seperator", type=str, default="\t",
                        help="The seperator used in the tables. This "
                        "will also be used in the output table. "
                        "Defaults to a tab.")
    parser.add_argument("-n", "--names", type=str, nargs="*", metavar="NAME",
                        help="The names of the samples corresponsing "
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
                        "to the merged count table. These attributes "
                        "will be retrieved from the GTF or GFF file "
                        "specified with the -g option. Multiple values "
                        "will be seperated by a ';'. Requires -g to be "
                        "specified.")
    parser.add_argument("-g", "--gtf", "--gff", type=str, metavar="FILE",
                        help="The GTF or GFF file from which the "
                        "additional attributes (see -a) will be "
                        "retrieved. Ignored if -a is not specified. "
                        "Required if -a is specified.")
    parser.add_argument("-F", "--feature-attribute", type=str, metavar="ATTR",
                        default="gene_id",
                        help="The attribute from the GTF/GFF used for "
                        "matching the feature records with the rows in "
                        "the count table. Ignored if -a is not "
                        "specified. Defaults to 'gene_id'.")
    args = parser.parse_args(argv)
    if args.additional_attributes is not None and args.gtf is None:
        parser.error("the following argument is required if -a is "
                     "specified: -g")
    return args


if __name__ == "__main__":
    args = parse_args(sys.argv[1:])
    main(args.table, args.output, args.feature_column, args.counts_column,
         args.sep, args.names, args.header, args.gtf,
         args.additional_attributes, args.feature_attribute)
