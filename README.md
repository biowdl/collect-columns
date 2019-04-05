# collect-columns

This tool retrieves a column from each in a set of tables and compiles into a single table.
Optionally, additional attributes from the associated GTF/GFF file may be added to the output
tables.

## Installation
* Clone the repository: `git clone https://github.com/biowdl/collect-columns.git`
* Enter the repository: `cd collect-columns`
* Install using pip: `pip install .`

## Usage
```
collect-columns output_path input_files...
```

It assumes that all input count tables are in the same format.
By default the format is assumed to be headerless and tab separated, with the
first column being the feature identifiers and the second the values of interest.
The output table will use the same separator as the input tables and contain
a header. The `feature` column will contain the feature identifiers, the value
columns will be named after the input files or according to the names given
through the `-n` option, which takes a list of names as argument.

In order to use a different input format the following options can be given:

| option | arguments | definition |
|:-:|:-:|:-|
| `-f` | a number | The index of the column containing the feature identifiers. |
| `-c` | a number | The index of the column containing the values/counts. |
| `-s` | a character | The separator.|
| `-H` | | Indicates that the table has a header. |

To add additional attributes from a GTF/GFF, the following options can be given:

| option | arguments | definition |
|:-:|:-:|:-|
| `-a` | a list of words | The attributes to be added to the output table. |
| `-g` | a path | The gtf file from which the attributes will be retrieved. |
| `-F` | a word | The attribute used to map rows in the input tables to gtf record. Defaults to `gene_id`. |

### Examples
#### HTSeq-count
Using the output from HTSeq-count as input the following command:
```
collect-columns all.tsv s1.tsv s2.tsv
```
will result in a table like:

| feature | s1.tsv | s2.tsv |
|:-------:|:------:|:------:|
| MSTRG.1 | 10     | 11     |
| MSTRG.2 | 60     | 12     |
| ...     | ...    | ...    |

#### Stringtie
Using stringtie abundance output as input, the following command:
```
collect-columns all.FPKM s1.abundance s2.abundance \
    -c 7 \
    -H \
    -a ref_gene_id gene_name \
    -g merged.gtf \
    -n sample1 sample2
```
will result in a table like:

| feature | ref_gene_id | gene_name | sample1       | sample2     |
|:-------:|:-----------:|:---------:|:-------------:|:-----------:|
| MSTRG.1 | g_1         | gene_1    | 185151.953125 | 151.964231  |
| MSTRG.2 | g_2         | gene_2    | 100160.070312 | 1160.030213 |
| ...     | ...         | ...       | ...           | ...         |
