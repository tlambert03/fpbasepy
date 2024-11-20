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


def test_get_camera() -> None:
    cam = fpbase.get_camera("Andor Zyla 5.5")
    repr(cam)
    assert cam.name == "Andor Zyla 5.5"


def test_get_light_source() -> None:
    light = fpbase.get_light_source("Lumencor Celesta UV")
    repr(light)
    assert light.name == "Lumencor Celesta UV"


def test_lists() -> None:
    assert len(fpbase.list_microscopes()) > 0
    assert len(fpbase.list_fluorophores()) > 0
    assert len(fpbase.list_filters()) > 0
    assert len(fpbase.list_cameras()) > 0
    assert len(fpbase.list_light_sources()) > 0
    assert len(fpbase.list_dyes()) > 0
    assert len(fpbase.list_proteins()) > 0


def test_generic_gql_query() -> None:
    data = fpbase.graphql_query("{proteins { name seq } }")
    EGFP = next(p for p in data["data"]["proteins"] if p["name"] == "EGFP")
    assert EGFP["seq"].startswith("MVSK")

    q = "query getProtein($id: String!){ protein(id: $id){ name } }"
    data = fpbase.graphql_query(q, {"id": "R9NL8"})
    assert data["data"]["protein"]["name"] == "EGFP"
