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
            "https://raw.githubusercontent.com/phe-genomics/variant_definitions/main/variant_yaml/imagines-viewable.yml",
        ]
    )

    # retrieve variants
    # TODO: mock network access
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
            "voc/omicron_ba2_mutations.yaml",
            "voc/omicron_ba1_mutations.yaml",
            "voc/delta_mutations.yaml",
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
