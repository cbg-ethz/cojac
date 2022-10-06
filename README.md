
# COJAC - CoOccurrence adJusted Analysis and Calling

[![Bioconda package](https://img.shields.io/conda/dn/bioconda/cojac.svg?label=Bioconda)](https://bioconda.github.io/recipes/cojac/README.html)
[![Docker container](https://quay.io/repository/biocontainers/cojac/status)](https://quay.io/repository/biocontainers/cojac)

## Description

The _cojac_ package comprises a set of command-line tools to analyse co-occurrence of mutations on amplicons. It is useful, for example, for early detection of viral variants of concern (e.g. Alpha, Delta, Omicron) in environmental samples, and has been designed to scan for multiple SARS-CoV-2 variants in wastewater samples, as analyzed jointly by [ETH Zurich](https://bsse.ethz.ch/news-and-events/d-bsse-news/2021/01/sars-cov-2-variants-detected-in-wastewater-samples.html), [EPFL](https://actu.epfl.ch/news/covid-19-using-wastewater-to-track-the-pandemic/) and [Eawag](https://www.eawag.ch/en/department/sww/projects/sars-cov2-in-wastewater/).
Learn more about this project on [its Dashboard](https://bsse.ethz.ch/cbg/research/computational-virology/sarscov2-variants-wastewater-surveillance.html), amplicon cooccurrences measured with _cojac_ are visualized on the heatmaps available on per-station or per-variant subpages [displayed on CoV-Spectrum](https://cov-spectrum.ethz.ch/story/wastewater-in-switzerland).

The analysis requires the whole amplicon to be covered by sequencing read pairs. It currently works at the level of aligned reads, but [we plan](#upcoming-features) to be able to adjust confidence scores based on local (window) haplotypes (as generated, e.g., by [ShoRAH](https://github.com/cbg-ethz/shorah), [doi:10.1186/1471-2105-12-119](https://doi.org/10.1186/1471-2105-12-119)).

## Usage

Here are the available command-line tools:

| command                                    | purpose |
| :----------------------------------------- | :------ |
| [`cojac cooc-mutbamscan`](cooc-mutbamscan) | scan an alignment BAM/CRAM/SAM file for mutation co-occurrences and output a JSON or YAML file |
| [`cojac cooc-colourmut`](cooc-colourmut)   | display a JSON or YAML file as a coloured output on the terminal |
| [`cojac cooc-pubmut`](cooc-pubmut)         | render a JSON or YAML file to a table as in the publication |
| [`cojac cooc-tabmut`](cooc-tabmut)         | export a JSON or YAML file as a CSV/TSV table for downstream analysis (e.g.: RStudio) |
| [`cojac cooc-curate`](cooc-curate)         | an (experimental) tool to assist evaluating the quality of variant definitions by looking at mutations' or cooccurences' frequencies from [CoV-Spectrum](cov-spectrum.ethz.ch) |
| [`cojac phe2cojac`](phe2cojac)             | a tool to generate new variant definition YAMLs for cojac using YMLs available at [PHE Genomic's _Standardised Variant Definitions_](https://github.com/phe-genomics/variant_definitions) |

Use option `-h` / `--help` to see available command-line options:

```console
$ cojac cooc-mutbamscan --help
Usage: cojac cooc-mutbamscan [OPTIONS]

  Scan amplicon (covered by long read pairs) for mutation cooccurrence

Options:
  -a, --alignments BAM/CRAM       alignment files
  -n, --name NAME                 when using alignment files, name to use for
                                  the output
  -s, --samples TSV               V-pipe samples list tsv
  --batchname SEP                 concatenate samplename/batchname from
                                  samples tsv
  -p, --prefix PATH               V-pipe work directory prefix for where to
                                  look at align files when using TSV samples
                                  list
  -r, --reference REFID           reference to look for in alignment files
  -m, --vocdir DIR                directory containing the yamls defining the
                                  variant of concerns
  --rev, --with-revert / --no-rev, --without-revert
                                  also include reverts when compiling
                                  amplicons (requires VOC YAML files with
                                  revert category)
  -b, --bedfile BED               bedfile defining the amplicons, with format:
                                  ref\tstart\tstop\tamp_num\tpool\tstrand
  -#, --cooc COOC                 minimum number of cooccurences to search for
  -Q, --amplicons, --in-amp, --in-amplicons YAML
                                  use the supplied YAML file to query
                                  amplicons instead of building it from BED +
                                  voc's DIR
  -A, --out-amp, --out-amplicons YAML
                                  output amplicon query in a YAML file
  -j, --json JSON                 output results to as JSON file
  -y, --yaml YAML                 output results to as yaml file
  -t, --tsv TSV                   output results to as (raw) tsv file
  -d, --dump                      dump the python object to the terminal
  --help                          Show this message and exit.

  @listfile can be used to pass a long list of parameters (e.g.: a large
  number of BAMs) in a file instead of command line
```

```console
$ cojac cooc-colourmut --help
Usage: cojac cooc-colourmut [OPTIONS]

  Print coloured pretty table on terminal

Options:
  -a, --amplicons YAML  list of query amplicons, from mutbamscan  [required]
  -j, --json JSON       results generated by mutbamscan
  -y, --yaml YAML       results generated by mutbamscan
  --help                Show this message and exit.

  See cooc-pubmut for a CSV file that can be imported into an article
```

```console
$ cojac cooc-pubmut --help
Usage: cojac cooc-pubmut [OPTIONS]

  Make a pretty table

Options:
  -m, --vocdir DIR        directory containing the yamls defining the variant
                          of concerns
  -a, --amplicons YAML    list of query amplicons, from mutbamscan
  -j, --json JSON         results generated by mutbamscan
  -y, --yaml YAML         results generated by mutbamscan
  -o, --output CSV        name of (pretty) csv file to save the table into
  -e, --escape            use escape characters for newlines
  -e, --escape            use escape characters for newlines
  -x, --excel             use a semi-colon ';' instead of a comma ',' in the
                          comma-separated-files as required by Microsoft Excel
  -, --batchname BOOLEAN  split samplename/batchname (as in samples tsv)
  -, --batchname          split samplename/batchname (as in samples tsv)
  -q, --quiet             Run quietly: do not print the table
  --help                  Show this message and exit.

  You need to open the CSV in a spreadsheet that understands linebreaks
```

```console
$ cojac cooc-tabmut --help
Usage: cojac cooc-tabmut [OPTIONS]

  Make a table suitable for further processing: RStudio, etc

Options:
  -j, --json JSON         results generated by mutbamscan
  -y, --yaml YAML         results generated by mutbamscan
  -, --batchname BOOLEAN  split samplename/batchname (as in samples tsv)
  -, --batchname          split samplename/batchname (as in samples tsv)
  -o, --output CSV        name of (pretty) csv file to save the table into
  -l, --lines             Line-oriented table alternative
  -x, --excel             use a semi-colon ';' instead of a comma ',' in the
                          comma-separated-files as required by Microsoft Excel
  -m, --multiindex        Use multi-level indexing (amplicons and counts
                          categories)
  -q, --quiet             Run quietly: do not print the table
  --help                  Show this message and exit.
```

```console
$ cojac cooc-curate --help
Usage: cojac cooc-curate [OPTIONS] [VOC]...

  Helps determining specific mutations and cooccurences by querying CoV-
  Spectrum

Options:
  -a, --amplicons YAML        use the YAML file generated by mutbamscan to
                              query amplicons instead of mutations
  -m, --mutations             always do mutations (even if amplicons YAML
                              provided)
  -H, --high FLOAT            Fraction above which a mutation must be found
                              among seeked lineages
  -l, --low FLOAT             Fraction under which a mutation must be found
                              among other lineages
  --collapse / --no-collapse  combine counts of all sublineages together and
                              consider a signle value that corresponds to a
                              lineages family (e.g.: count all B.1.612.2*
                              together). This is especially useful for
                              assessing signature of old variants that have
                              branched out by now.
  --colour / --no-colour      use coloured output
  --help                      Show this message and exit.

  This tool queries LAPIS, see https://lapis.cov-spectrum.org/swagger/ and
  https://lapis.cov-spectrum.org/
```

```console
$ cojac phe2cojac  --help
usage: phe2cojac [-h] [-s SHRT] [-y [OUT_YAML]] IN_YAML

convert phe-genomics to cojac's dedicated variant YAML format

positional arguments:
  IN_YAML               phe-genomics variant YAML input file

options:
  -h, --help            show this help message and exit
  -s SHRT, --shortname SHRT
                        shortname to use (otherwise auto-build one based on phe-genomic's unique id)
  -y [OUT_YAML], --yaml [OUT_YAML]
                        write cojac variant to a YAML file instead of printing (if empty, build filename from shortname)
```

## Howto

### Input data requirements

Analysis needs to be performed on SARS-CoV-2 samples sequenced using a tiled multiplexed PCRs protocol for which you need a BED (Browser Extensible Data) file describing the amplified regions, and sequenced with read settings that covers the totality of an amplicon.

We provide BED files for the following examples:
 - [nCoV-2019.insert.V3.bed](nCoV-2019.insert.V3.bed) for [ARTIC V3](https://doi.org/10.17504/protocols.io.bibtkann)
 - [SARS-CoV-2.insert.V4.txt](SARS-CoV-2.insert.V4.txt) for [ARTIC V4](https://community.artic.network/t/sars-cov-2-version-4-scheme-release/312)

These protocols produces ~400bp long amplicons, and thus needs to be sequenced with, e.g., paired end sequencing with read length 250.

Select the desired bedfile using the `-b` / `--bedfile` option.

> **Note:**
> - this analysis method cannot work on read length much shorter than the amplicons (e.g.: it will not give reliable results for a read-length of 50).
> - to use different protocols (e.g. Nimagen), you need to provide a BED file describing the amplicons. Its columns "start" and "stop" are mandatory

Analysis will use variants description YAML that list mutation to be searched.

We provide several examples in the directory [`voc/`](voc/).

Select a directory containing a collection of virus definitions YAMLs using the `-m` / `--vocdir` option.

> **Note:**
> - you can create new YAML files if you need to look for new variants of concern.
> - e.g. it is possible to automatically generate YAMLs for cojac from  [PHE Genomic's _Standardised Variant Definitions_](https://github.com/phe-genomics/variant_definitions):
```bash
# fetch the repository of standardised variant definitions
git clone https://github.com/phe-genomics/variant_definitions.git
# generate a YAML for omicron subvariant BA.2 using the corresponding standardised variant definitions
phe2cojac --shortname 'om2' --yaml voc/omicron_ba2_mutations.yaml variant_definitions/variant_yaml/imagines-viewable.yml
# now have a look at the frequencies of mutations using CoV-Spectrum
cooc-curate voc/omicron_ba2_mutations.yaml
# adjust the content of the YAML files to your needs
```

### Collect the co-occurrence data

There are currently two modes to collect the data about co-occurring mutations in reads: analysing stand-alone BAM/CRAM/SAM alignment files, or analysing the output of a cohort analysed with [V-pipe](https://cbg-ethz.github.io/V-pipe/) ([doi:10.1093/bioinformatics/btab015](https://doi.org/10.1093/bioinformatics/btab015)).

#### Standalone files

Provide a list of BAM files using the `-a` / `--alignment` option. Run:

```bash
cojac cooc-mutbamscan -b nCoV-2019.insert.V3.bed -m voc/ -a sam1.bam sam2.bam -j cooc-test.json
```

> **Note:** you can also use the `-y` / `--yaml` option to write to a YAML file instead of a JSON.

#### Analyzing a cohort with V-pipe

You can learn how to analyse _fastq.gz_ files with V-pipe with this tutorial:

 - https://cbg-ethz.github.io/V-pipe/tutorial/sars-cov2/ ([video tutorial](https://youtu.be/pIby1UooK94))

Run:

```bash
cojac cooc-mutbamscan -b nCoV-2019.insert.V3.bed -m voc/ -t work/samples.tsv -p work/samples/ -j cooc-test.json
```

#### Number of cooccurences

By default `cooc-mutbamscan` will look for cooccurrences of at least 2 mutations on the same amplicon. You can change that number using option `-#`/`--cooc`:

 - you can increase it to e.g.: 3 if the variants you study requires more stringent identification
 - you can set it to 1, to also count isolated occurrences - in this case `cooc-mutbamscan` will also double as a generic (non coorcurrence-aware) variant caller, so you can get all counts with a single tool.

#### Store the amplicon query

Using the `-A` / `--out-amp` / `--out-amplicons` option, it is possible to store the exact request that was used to analyze samples.
You can then re-use the exact same request using the `-Q` / `--in-amp` / `--amplicons` option, or pass it to a visualisation tool.

```bash
# store the request in a YAML file
cojac cooc-mutbamscan -b nCoV-2019.insert.V3.bed -m voc/ -A amplicons.v3.yaml
# adjust the content of amplicons.v3.yaml

# now have a look at the frequencies of mutation cooccurences using CoV-Spectrum
cojac cooc-curate -a amplicons.v3.yaml voc/omicron_ba2_mutations.yaml voc/omicron_ba1_mutations.yaml voc/delta_mutations.yaml
# reuse the amplicon
cojac cooc-mutbamscan -Q amplicons.v3.yaml -a sam1.bam sam2.bam -j cooc-test.json
```

### Display data on terminal

The default `-d` / `--dump` option of `cooc-mutbamscan` is not a very user-friendly experience to display the data. You can instead pass a JSON or YAML file to the display script. Run:

```bash
cojac cooc-colourmut -a amplicons.v3.yaml -j cooc-test.json
```

![terminal screen shot](images/terminal.svg)

> **Notes:**
> - passing the `-a` / `--amplicons` parameter is currenlty mandatory, see [section _Store the amplicon query_ above](#store_the_amplicon_query)

### Render table for publication

And now, let’s go beyond our terminal and produce a table that can be included in a publication (see bibliography below for concrete example). Run:

```bash
cojac cooc-pubmut -m voc/ -a amplicons.v3.yaml -j cooc-test.json -o cooc-output.tsv
```

> **Note:**
> - if provided options `-m` / `--vocdir` and `-a` /  `--amplicons` can help generate human-friendly headers (_Amplicon 88, 26277-26635_) in the table instead of short names (`88_om`)
> - you can also output to comma-separated table (`-o cooc-output.csv`)
> - Microsoft Excel requires using option `-x`/`--excel` (using semi-colon instead of comma in comma-separated-value files). Some versions can also open TSV (but not the Office 365 web app).

You need to open the table with a spread-sheet that can understand line breaks, such as [LibreOffice Calc](https://www.libreoffice.org/discover/calc/), [Google Docs Spreadsheet](https://www.google.com/sheets/about/) or, using special options (see above), [Microsoft Excel](https://www.microsoft.com/en-us/microsoft-365/excel).

|          | 72_al               | 78_al            | 92_al              | 93_al               | 76_be             | 77_d614g               |
| :------- | ------------------: | ---------------: | -----------------: | ------------------: | ----------------: | ---------------------: |
| sam1.bam | 158 / 809<br>19.53% | 2 / 452<br>0.44% | 89 / 400<br>22.25% | 344 / 758<br>45.38% | 0 / 1090<br>0.00% | 371 / 371<br>100.00%   |
| sam2.bam | 0 / 1121<br> 0.00%  | 0 / 255<br>0.00% | 58 / 432<br>13.43% | 142 / 958<br>14.82% | 0 / 1005<br>0.00% | 1615 / 1615<br>100.00% |


It is also possible to use the software [_pandoc_](https://pandoc.org/) to further convert the CSV to other formats. Run:

```bash
cojac cooc-pubmut -j cooc-test.json -o cooc-output.csv
pandoc cooc-output.csv -o cooc-output.pdf
pandoc cooc-output.csv -o cooc-output.html
pandoc cooc-output.csv -o cooc-output.md
```

### Export table for downstream analysis

If you want to further analyse the data (e.g.: with RStudio), it's also possible to export the data into a more machine-readable CSV/TSV table. Run:

```bash
cojac cooc-tabmut -j cooc-test.json -o cooc-export.csv
```

You can try importing the resulting CSV in you favourite tool.

|          | A72_al.count | A72_al.mut_all | A72_al.mut_oneless | A72_al.frac | A72_al.cooc | A78_al.count | A78_al.mut_all | A78_al.mut_oneless | A78_al.frac | A78_al.cooc | ... |
| :------- | -----------: | -------------: | -----------------: | ----------: | ----------: | -----------: | -------------: | -----------------: | ----------: | ----------: | --- |
| sam1.bam |          809 |            158 |                234 |    0.195303 |           2 |          452 |              2 |                  7 |    0.004425 |           2 | ... |
| sam2.bam |         1121 |              0 |                  0 |    0.000000 |           2 |          255 |              0 |                 52 |    0.000000 |           2 | ... |

The columns are tagged as following:

 - **count**: total count of amplicons carrying the sites of interest
 - **mut_all**: amplicons carrying mutations on all site of interest (e.g.: variant mutations observed on all sites)
 - **mut_oneless**: amplicons where one mutation is missing (e.g.: only 2 out of 3 sites carried the variant mutation, 1 sites carries wild-type)
 - **frac**: fraction _(mut_all/count)_ or empty if no _counts_
 - **cooc**: number of considered site (e.g.: 2 sites of interests) or empty if no _counts_

If your tool supports multi-level indexing, use the `-m`/`--multiindex` option. The resulting table will be bilevel indexed:
the first level is the amplicon, the second is the category.

<table>
<thead>
<tr><th></th><th colspan="5">A72_al</th><th colspan="5">A78_al</th></tr>
<tr><th></th><th>count</th><th>mut_all</th><th>mut_oneless</th><th>frac</th><th>cooc</th><th>count</th><th>mut_all</th><th>mut_oneless</th><th>frac</th><th>cooc</th></tr>
</thead>
<tbody>
<tr><td>sam1.bam</td><td>809</td><td>158</td><td>234</td><td>0.195303</td><td>2</td><td>452</td><td>2</td><td>7</td><td>0.004425</td><td>2</td></tr>
<tr><td>sam2.bam</td><td>1121</td><td>0</td><td>0</td><td>0.0</td><td>2</td><td>255</td><td>0</td><td>52</td><td>0.0</td><td>2</td></tr>
</tbody>
</table>

Another different table orientation is provided by `-l`/`--lines`:

| sample   | amplicon | frac     | cooc | count | mut_all | mut_oneless | al | be | d614g |
| :------- | -------: | :------- | :--: | ----: | ------: | ----------: |:--:|:--:| :---: |
| sam1.bam |       72 | 0.195303 |    2 |   809 |     158 |         234 |  1 |
| sam1.bam |       78 | 0.004425 |    2 |   452 |       2 |           7 |  1 |
| sam1.bam |       92 | 0.222500 |    3 |   400 |      89 |           3 |  1 |
| sam1.bam |       93 | 0.453826 |    2 |   758 |     344 |         140 |  1 |
| sam1.bam |       76 | 0.000000 |    2 |  1090 |       0 |         377 |    |  1 |
| sam1.bam |       77 | 1.000000 |    1 |   371 |     371 |           0 |    |    |  1 |
| sam2.bam |       72 | 0.000000 |    2 |  1121 |       0 |           0 |  1 |
| sam2.bam |       78 | 0.000000 |    2 |   255 |       0 |          52 |  1 |
| sam2.bam |       92 | 0.134259 |    3 |   432 |      58 |           3 |  1 |
| sam2.bam |       93 | 0.148225 |    2 |   958 |     142 |          80 |  1 |
| sam2.bam |       76 | 0.000000 |    2 |  1005 |       0 |           0 |    |  1 |
| sam2.bam |       77 | 1.000000 |    1 |  1615 |    1615 |           0 |    |    |  1 |


## Installation

We recommend using [bioconda software repositories](https://bioconda.github.io/index.html) for easy installation.
You can find instruction to setup your bioconda environment at the following address:

 - https://bioconda.github.io/user/install.html

In those instructions, please follow carefully the [section _2. Set up channels_](https://bioconda.github.io/user/install.html#set-up-channels).

If you use [V-pipe’s `quick_install.sh`](https://github.com/cbg-ethz/V-pipe#using-quick-install-script), it will set up an environment that you can activate, e.g.:

```bash
bash quick_install.sh -b sars-cov2 -p testing -w work
. ./testing/miniconda3/bin/activate
```

### Prebuilt package

_cojac_ and its dependencies are all available in the bioconda repository. We strongly advise you to install [this pre-built package](https://bioconda.github.io/recipes/cojac/README.html) for a hassle-free experience.

You can install _cojac_ in its own environment and activate it:

```bash
conda create -n cojac cojac
conda activate cojac
# test it
cojac cooc-mutbamscan --help
```

And to update it to the latest version, run:

```bash
# activate the environment if not already active:
conda activate cojac
conda update cojac
```

Or you can add it to the current environment (e.g.: in environment _base_):

```bash
conda install cojac
```

### Dependencies

If you want to install the software yourself, you can see the list of dependencies in [`conda_cojac_env.yaml`](conda_cojac_env.yaml).

We recommend using conda to install them:

```bash
conda env create -f conda_cojac_env.yaml
conda activate cojac
```

Install _cojac_ using pip:
```bash
pip install .
# this will autodetect dependencies already installed by conda
```

cojac should now be accessible from your PATH

```bash
# activate the environment if not already active:
conda activate cojac
cojac cooc-mutbamscan --help
```

### Remove conda environment

You can remove the conda environment if you don't need it any more:

```bash
# exit the cojac environment first:
conda deactivate
conda env remove -n cojac
```

### Python Package Index

Alternatively, you can install `cojac` using pip:

```bash
pip install cojac
```

## Additional notebooks

The subdirectory [`notebooks/`](notebooks/) contains Jupyter and Rstudio notebooks used in the [publication](#citation).

## Upcoming features

- [x] ~~bioconda package~~
- [x] ~~further jupyter and rstudio code from the publication~~
- [x] ~~Move hard-coded amplicons to BED input file~~
- [x] ~~Move hard-coded mutations to YAML configuration~~
- [x] ~~Refactor code into proper Python package~~

Long term goal:

- [ ] Integration as part of V-pipe
- [ ] Integration with ShoRAH amplicons

## Contributions

#### Package developers:

- [David Dreifuss ![orcid]](https://orcid.org/0000-0002-5827-5387), [![github]](https://github.com/dr-david)
- [Kim Philipp Jablonski ![orcid]](https://orcid.org/0000-0002-4166-4343), [![github]](https://github.com/kpj)
- [Ivan Topolsky ![orcid]](https://orcid.org/0000-0002-7561-0810), [![github]](https://github.com/dryak)

#### Additional notebooks:

 - [Lara Fuhrmann ![orcid]](https://orcid.org/0000-0001-6405-0654), [![github]](https://github.com/LaraFuhrmann)
 - [Kim Philipp Jablonski ![orcid]](https://orcid.org/0000-0002-4166-4343), [![github]](https://github.com/kpj)
 - [Katharina Jahn ![orcid]](https://orcid.org/0000-0002-6983-4615), [![github]](https://github.com/jahnka)

#### Corresponding author:

 - [Niko Beerenwinkel ![orcid]](https://orcid.org/0000-0002-0573-6119)

[github]: images/mark-github.svg
[orcid]: images/ORCIDiD_iconvector.svg


## Citation

If you use this software in your research, please cite:

- Katharina Jahn, David Dreifuss, Ivan Topolsky, Anina Kull, Pravin Ganesanandamoorthy, Xavier Fernandez-Cassi, Carola Bänziger, Elyse Stachler, Lara Fuhrmann, Kim Philipp Jablonski, Chaoran Chen, Catharine Aquino, Tanja Stadler, Christoph Ort, Tamar Kohn, Timothy R. Julian, Niko Beerenwinkel

  "*Detection of SARS-CoV-2 variants in Switzerland by genomic analysis of wastewater samples*."

  medRxiv 2021.01.08.21249379; [doi:10.1101/2021.01.08.21249379](https://doi.org/10.1101/2021.01.08.21249379)

## Contacts

If you experience problems running the software:

- We encourage to use the [issue tracker on GitHub](https://github.com/cbg-ethz/cojac/issues)
- For further enquiries, you can also contact the [V-pipe Dev Team](https://cbg-ethz.github.io/V-pipe/contact/)
- You can contact the publication’s corresponding author
