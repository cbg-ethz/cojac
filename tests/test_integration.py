import subprocess


def test_workflow():
    # download data (TODO: use better sources)
    data_sources = [
        "https://polybox.ethz.ch/index.php/s/vncvlX4pAGo4h33/download?path=%2F&files=REF_aln.bam",
        "https://polybox.ethz.ch/index.php/s/vncvlX4pAGo4h33/download?path=%2F&files=sam1.bam",
        "https://polybox.ethz.ch/index.php/s/vncvlX4pAGo4h33/download?path=%2F&files=sam1.bam.bai",
        "https://polybox.ethz.ch/index.php/s/vncvlX4pAGo4h33/download?path=%2F&files=sam2.bam",
        "https://polybox.ethz.ch/index.php/s/vncvlX4pAGo4h33/download?path=%2F&files=sam2.bam.bai",
    ]

    for url in data_sources:
        name = url.split("=")[-1]
        print(name)
        subprocess.run(["wget", "--no-clobber", "-O", name, url])

    subprocess.check_call(
        [
            "wget",
            "--no-clobber",
            "https://raw.githubusercontent.com/ukhsa-collaboration/variant_definitions/main/variant_yaml/imagines-viewable.yml",
            "https://raw.githubusercontent.com/hodcroftlab/covariants/master/defining_mutations/21K.Omicron.tsv",
            "https://raw.githubusercontent.com/cbg-ethz/cowwid/master/voc/omicron_ba286_mutations_full.yaml",
        ]
    )

    subprocess.run(
        [
            "wget",
            "--no-clobber",
            "-O",
            "SARS-CoV-2.v532.primer.bed",
            "https://raw.githubusercontent.com/artic-network/primer-schemes/master/nCoV-2019/V5.3.2/SARS-CoV-2.primer.bed",
        ]
    )

    # dummy run
    subprocess.check_call(["cojac", "--version"])

    # retrieve variants
    # TODO: mock network access
    # alternative 1: CovSpectrum
    subprocess.check_call(
        [
            "cojac",
            "sig-generate",
            "--url",
            "https://lapis.cov-spectrum.org/open/v1",
            "--variant",
            "BA.1",
        ]
    )

    # alternative 2: covariants
    subprocess.check_call(
        [
            "cojac",
            "sig-generate",
            "--covariants",
            "21K.Omicron.tsv",
            "--variant",
            "BA.1",
        ]
    )

    # alternative 3: UKHSA
    subprocess.check_call(
        [
            "cojac",
            "phe2cojac",
            "--shortname",
            "om2",
            "--yaml",
            "omicron_ba2_mutations.yaml",
            "imagines-viewable.yml",
        ]
    )

    # check frequencies
    subprocess.check_call(["cojac", "cooc-curate", "omicron_ba2_mutations.yaml"])

    # amplicon queries
    subprocess.check_call(
        [
            "cojac",
            "cooc-mutbamscan",
            "-b",
            "nCoV-2019.insert.V3.bed",
            "-m",
            "voc/",
            "-A",
            "amplicons.v3.yaml",
        ]
    )
    subprocess.check_call(
        [
            "cojac",
            "cooc-curate",
            "-a",
            "amplicons.v3.yaml",
            "voc/omicron_ba2_mutations_full.yaml",
            "voc/omicron_ba1_mutations_full.yaml",
            "voc/delta_mutations_full.yaml",
        ]
    )
    subprocess.check_call(
        [
            "cojac",
            "cooc-mutbamscan",
            "-Q",
            "amplicons.v3.yaml",
            "-a",
            "sam1.bam",
            "-a",
            "sam2.bam",
            "-j",
            "cooc-test.json",
        ]
    )

    # display results
    subprocess.check_call(
        ["cojac", "cooc-colourmut", "-a", "amplicons.v3.yaml", "-j", "cooc-test.json"]
    )

    # create result tsv and csv
    subprocess.check_call(
        [
            "cojac",
            "cooc-pubmut",
            "-m",
            "voc/",
            "-a",
            "amplicons.v3.yaml",
            "-j",
            "cooc-test.json",
            "-o",
            "cooc-output.tsv",
        ]
    )
    subprocess.check_call(
        ["cojac", "cooc-tabmut", "-j", "cooc-test.json", "-o", "cooc-export.csv"]
    )

    # primers affected by mutations
    subprocess.check_call(
        [
            "cojac",
            "cooc-mutbamscan",
            "--voc",
            "omicron_ba286_mutations_full.yaml",
            "--bedfile",
            "SARS-CoV-2.v532.primer.bed",
            "--no-sort",
            "--cooc",
            "1",
            "--out-amp",
            "affected_primers.v532.yaml",
        ]
    )
