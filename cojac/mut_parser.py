#!/usr/bin/env python3
import itertools
import re

rxmutdec = re.compile(
    r"^(?:(?:(?:(?P<ref>[ATCG]+)\>)?(?P<mut>[ATCG]+))|(?P<del>[\-]+)|(?:[\+](?P<ins>[ATGC]+)))$"
)


def mut_decode(mutstr):
    """function to decode the strings of a mutation
    mutation can specifie:
     - a single base (e.g.: 'C') or a run of bases ('CAT')
     - a substitution (e.g.: 'G>C'), can also be longer ('GTC>CAT')
     - a run of deletions ('---')
     - an insertions ('+TATA')

    (Basically, this function strips the reference part)

    Input:
            a string with any of the above
    Returns:
            a string with bases or deletions (either single or runs)
    """

    res = rxmutdec.match(mutstr)

    if res:
        match = res.groupdict()
        if match["mut"]:
            # print(f"{mutstr} : mutation {match['mut']}")
            return match["mut"]
        if match["del"]:
            # print(f"{mutstr} : deletion {match['del']}")
            return match["del"]
        if match["ins"]:
            print(f"insertions not supported (yet): {mutstr} : {match['ins']}")
            return None

    raise ValueError(f"cannot parse mutation f{mutstr}")


def filter_decode_vartiant(yam, categories=["mut", "extra", "shared", "subset"]):
    """function to filter mutations and decode the strings
    mutations can be classified as:
     - 'mut': defining mut
     - 'extra': officially not counted as defining but still specific
     - 'shared': present in this variant, but also occuring outside of it (e.g.: common ancestror)
     - 'subset': only present on some variants (i.e.: 'not fixed')
    mutation can be specifie:
     - a single base (e.g.: 'C') or a run of bases ('CAT')
     - a substitution (e.g.: 'G>C'), can also be longer ('GTC>CAT')
     - a run of deletions ('---')
     - an insertions ('+TATA')

    Input:
            yam: dictionnary loaded from a voc YAML
            categories: a list of categories to filter-in
    Returns:
            voc: dict() with
                    'name': short name of the variant
                    'mut': a an ordereddict() of 'position': 'bases or deletions'
    """

    return {
        "name": yam["variant"]["short"],
        # sort so that definitions with same mutations in different categories all yield identic result
        # (e.g.: YAML definitions with either "{mut:{123:A,567:C}}" or "{mut:{567:C},extra:{123:A}}"  or "{mut:{123:A},extra:{567:C}}"
        # will all always consistently yield a tmp_mut_dict with "{123:A,567:C}" )
        "mut": dict(
            sorted(
                (pos, mut_decode(mut))
                for c in categories
                if c in yam
                for (pos, mut) in yam[c].items()
                if mut_decode(mut)
            )
        ),
    }


def to_single(mutstrs, nodels=False):
    """function to convert to the single mutations format used by standard tools like covSpectrum, etc.
    YAML:
       1432: GTC>CAT
    output:
       G1432C, T1433A, C1434T
    """

    for pos, mstr in mutstrs.items():
        rxm = rxmutdec.match(mstr)
        ref = mut = None

        if rxm:
            match = rxm.groupdict()

            if match["ins"]:
                print(f"insertions not supported (yet): {mutstr} : {match['ins']}")
                continue

            if match["mut"]:
                ref = match.get("ref", "")
                mut = match["mut"]
            if match["del"]:
                if nodels:
                    continue
                ref = ""
                mut = match["del"]

            if mut is not None:
                for r, p, m in zip(
                    list(
                        [""] * len(mut)
                        if ref is None or len(mut) > len(ref)
                        else list(ref)
                    ),
                    itertools.count(int(pos)),
                    list(mut),
                ):
                    yield f"{r}{p}{m}"
                continue

        raise ValueError(f"cannot parse mutation {mstr}")
