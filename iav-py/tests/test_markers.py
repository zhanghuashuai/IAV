import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))

from markers import (  # noqa: E402
    is_applicable,
    marker_found,
    parse_deletions,
    parse_substitutions,
)


def test_parse_substitutions():
    subs = parse_substitutions("HA:Q226L,NA:H275Y")
    assert len(subs) == 2
    assert subs[0]["gene"] == "HA"
    assert subs[0]["pos"] == 226
    assert subs[1]["ref"] == "H"
    assert subs[1]["alt"] == "Y"


def test_substitution_marker_match():
    subs = parse_substitutions("NA:H275Y")
    marker = {
        "type": "substitution",
        "gene": "NA",
        "position": 275,
        "ref": "H",
        "alt": "Y",
    }
    found, observed = marker_found(marker, subs, [])
    assert found is True
    assert observed == "NA:H275Y"


def test_substitution_marker_mismatch():
    subs = parse_substitutions("NA:H275Y")
    marker = {
        "type": "substitution",
        "gene": "NA",
        "position": 275,
        "ref": "H",
        "alt": "N",
    }
    found, observed = marker_found(marker, subs, [])
    assert found is False
    assert observed == ""


def test_deletion_marker_match():
    dels = parse_deletions("NS1:80-84")
    marker = {
        "type": "deletion",
        "gene": "NS1",
        "start": 80,
        "end": 84,
    }
    found, observed = marker_found(marker, [], dels)
    assert found is True
    assert observed == "NS1:80-84"


def test_gene_prefix_mismatch():
    subs = parse_substitutions("M1:S31N")
    marker = {
        "type": "substitution",
        "gene": "M2",
        "position": 31,
        "ref": "S",
        "alt": "N",
    }
    found, _ = marker_found(marker, subs, [])
    assert found is False


def test_dataset_filter():
    marker = {"dataset_contains": "h1n1pdm"}
    assert is_applicable(marker, "flu_h1n1pdm_na") is True
    assert is_applicable(marker, "flu_h3n2_na") is False
    assert is_applicable(marker, "") is True
    assert is_applicable({}, "flu_h3n2_na") is True
