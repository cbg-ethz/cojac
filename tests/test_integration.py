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

    subprocess.run(
        [
            "wget",
            "--no-clobber",
            "https://raw.githubusercontent.com/phe-genomics/variant_definitions/main/variant_yaml/imagines-viewable.yml",
        ]
    )

    # retrieve variants
    # TODO: mock network access
    subprocess.run(
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
    subprocess.run(["cojac", "cooc-curate", "omicron_ba2_mutations.yaml"])
