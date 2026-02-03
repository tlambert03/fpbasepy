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
    result = fpbase.list_microscopes()
    assert result
    assert fpbase.get_microscope(result[0])

    result = fpbase.list_fluorophores()
    assert result
    assert fpbase.get_fluorophore(result[0])

    result = fpbase.list_filters()
    assert result
    assert fpbase.get_filter(result[0])

    result = fpbase.list_cameras()
    assert result
    assert fpbase.get_camera(result[0])

    result = fpbase.list_light_sources()
    assert result
    assert fpbase.get_light_source(result[0])

    result = fpbase.list_dyes()
    assert result
    assert fpbase.get_fluorophore(result[0])

    result = fpbase.list_proteins()
    assert result
    assert fpbase.get_protein(result[0])


def test_generic_gql_query() -> None:
    data = fpbase.graphql_query("{proteins { name seq } }")
    EGFP = next(p for p in data["data"]["proteins"] if p["name"] == "EGFP")
    assert EGFP["seq"].startswith("MVSK")

    q = "query getProtein($id: String!){ protein(id: $id){ name } }"
    data = fpbase.graphql_query(q, {"id": "R9NL8"})
    assert data["data"]["protein"]["name"] == "EGFP"


@pytest.mark.parametrize("name", ["Clover1.5", "6C", "dClover2 A206K"])
def test_fluors_with_no_pdb(name: str) -> None:
    fpbase.get_fluorophore(name)


def test_get_multiple_proteins() -> None:
    proteins = fpbase.get_multiple_proteins(["EGFP", "mCherry", "mTurquoise2"])
    assert len(proteins) == 3
    assert proteins["EGFP"] is not None
    assert proteins["EGFP"].name == "EGFP"
    assert proteins["mCherry"] is not None
    assert proteins["mCherry"].name == "mCherry"
    assert proteins["mTurquoise2"] is not None
    assert proteins["mTurquoise2"].name == "mTurquoise2"
    # Check that spectral data is included
    assert proteins["EGFP"].default_state.spectra
    assert proteins["EGFP"].default_state.excitation_spectrum is not None
    assert proteins["EGFP"].default_state.emission_spectrum is not None


def test_get_multiple_proteins_with_invalid() -> None:
    """Test that invalid names return None in results."""
    proteins = fpbase.get_multiple_proteins(["EGFP", "NotAProtein"])
    assert len(proteins) == 2
    assert proteins["EGFP"] is not None
    assert proteins["EGFP"].name == "EGFP"
    assert proteins["NotAProtein"] is None


def test_get_multiple_proteins_empty() -> None:
    """Test empty input."""
    proteins = fpbase.get_multiple_proteins([])
    assert proteins == {}


def test_get_multiple_dyes() -> None:
    # Use actual dye names from the database
    dyes = fpbase.get_multiple_dyes(["DAPI", "Hoechst 33342"])
    assert len(dyes) == 2
    assert dyes["DAPI"] is not None
    assert dyes["DAPI"].name == "DAPI"
    assert dyes["Hoechst 33342"] is not None
    assert dyes["Hoechst 33342"].name == "Hoechst 33342"
    # Check that spectral data is included
    assert dyes["DAPI"].default_state.spectra


def test_get_multiple_dyes_with_invalid() -> None:
    """Test that invalid dye names return None in results."""
    dyes = fpbase.get_multiple_dyes(["DAPI", "NotADye"])
    assert len(dyes) == 2
    assert dyes["DAPI"] is not None
    assert dyes["NotADye"] is None


def test_get_multiple_microscopes() -> None:
    # Use actual microscope IDs
    microscopes = fpbase.get_multiple_microscopes(["wKqWbgApvguSNDSRZNSfpN"])
    assert len(microscopes) == 1
    assert microscopes["wKqWbgApvguSNDSRZNSfpN"] is not None
    assert microscopes["wKqWbgApvguSNDSRZNSfpN"].name == "Example Simple Widefield"
    # Check that optical configs are included
    assert microscopes["wKqWbgApvguSNDSRZNSfpN"].opticalConfigs


def test_get_multiple_microscopes_with_invalid() -> None:
    """Test that invalid microscope IDs return None in results."""
    microscopes = fpbase.get_multiple_microscopes(
        ["wKqWbgApvguSNDSRZNSfpN", "InvalidID"]
    )
    assert len(microscopes) == 2
    assert microscopes["wKqWbgApvguSNDSRZNSfpN"] is not None
    assert microscopes["InvalidID"] is None
