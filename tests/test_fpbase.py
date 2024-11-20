import pytest

import fpbase


def test_get_microscope() -> None:
    scope = fpbase.get_microscope("wKqWbgApvguSNDSRZNSfpN")
    repr(scope)
    assert scope.name == "Example Simple Widefield"


@pytest.mark.parametrize("name", ["EGFP", "Alexa Fluor 488"])
def test_get_fluor(name: str) -> None:
    fluor = fpbase.get_fluorophore(name)
    repr(fluor)
    assert fluor.name == name
    assert fluor.default_state
    assert fluor.default_state.excitation_spectrum is not None
    assert fluor.default_state.emission_spectrum is not None


@pytest.mark.parametrize("name", ["mEos3.2", "mScarlet-I"])
def test_get_protein(name: str) -> None:
    prot = fpbase.get_protein(name)
    repr(prot)
    assert prot.name == name
    assert prot.default_state.excitation_spectrum is not None
    assert prot.default_state.emission_spectrum is not None


def test_get_missing_protein() -> None:
    with pytest.raises(ValueError, match="Did you mean 'mscarlet'"):
        fpbase.get_protein("mScrlet")


@pytest.mark.parametrize("name", ["Chroma ET525/50m", "Semrock FF01-520/35"])
def test_get_filter(name: str) -> None:
    filt = fpbase.get_filter(name)
    repr(filt)
    assert filt.name == name


def test_lists() -> None:
    assert len(fpbase.list_microscopes()) > 0
    assert len(fpbase.list_fluorophores()) > 0
    assert len(fpbase.list_filters()) > 0
    assert len(fpbase.list_proteins()) > 0
