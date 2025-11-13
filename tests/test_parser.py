from cojac.mut_parser import to_single


def test_to_single():
    voc = {
        193: "C>T",
        241: "T",
        515: "---",
        22577: "GC>TA",
        22881: "GG",
    }
    out_dels = [
        "C193T",
        "241T",
        "515-",
        "516-",
        "517-",
        "G22577T",
        "C22578A",
        "22881G",
        "22882G",
    ]
    out_nodels = ["C193T", "241T", "G22577T", "C22578A", "22881G", "22882G"]

    assert list(to_single(voc, nodels=False)) == out_dels
    assert list(to_single(voc, nodels=True)) == out_nodels
