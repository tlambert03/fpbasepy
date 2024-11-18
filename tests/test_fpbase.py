import pytest

from fpbase import get_filter, get_fluorophore, get_microscope


def test_get_microscope() -> None:
    scope = get_microscope("wKqWbgApvguSNDSRZNSfpN")
    assert scope.name == "Example Simple Widefield"


@pytest.mark.parametrize("name", ["EGFP", "Alexa Fluor 488"])
def test_get_fluor(name: str) -> None:
    fluor = get_fluorophore(name)
    assert fluor.name == name
    assert fluor.default_state
    assert fluor.default_state.excitation_spectrum is not None
    assert fluor.default_state.emission_spectrum is not None


@pytest.mark.parametrize("name", ["Chroma ET525/50m", "Semrock FF01-520/35"])
def test_get_filter(name: str) -> None:
    filt = get_filter(name)
    assert filt.name == name
