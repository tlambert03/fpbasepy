from __future__ import annotations

from typing import Any


def build_multiple_query(
    query_name: str,
    items: dict[str, dict[str, Any]],
    fields: str,
) -> str:
    """Build a GraphQL query to fetch multiple items using aliases.

    Parameters
    ----------
    query_name : str
        The GraphQL query name (e.g., "protein", "dye")
    items : dict[str, dict[str, Any]]
        Mapping of alias -> query arguments
        e.g., {"egfp": {"id": "R9NL8"}, "mcherry": {"id": "ZERB6"}}
    fields : str
        The fields to query for each item

    Returns
    -------
    str
        GraphQL query string with aliases
    """
    query_parts = []
    for alias, args in items.items():
        # Build argument string: id: "value", name: "value"
        args_str = ", ".join(
            f'{k}: "{v}"' if isinstance(v, str) else f"{k}: {v}"
            for k, v in args.items()
        )
        query_parts.append(f"{alias}: {query_name}({args_str}) {{ {fields} }}")

    return "{ " + " ".join(query_parts) + " }"


MICROSCOPE_QUERY = """
query getMicroscope($id: String!) {
    microscope(id: $id) {
        id
        name
        opticalConfigs {
            name
            filters {
                path
                reflects
                filter {
                  id
                  name
                  spectrum { id subtype data }
                }
            }
            camera { id name spectrum { id subtype data } }
            light { id name spectrum { id subtype data } }
            laser
        }
    }
}
"""

DYE_QUERY = """
query getDye($id: Int!) {
    dye(id: $id) {
        name
        id
        exMax
        emMax
        extCoeff
        qy
        spectra { id subtype data }
    }
}
"""

PROTEIN_QUERY = """
query getProtein($id: String!) {
    protein(id: $id) {
        name
        id
        seq
        genbank
        pdb
        uniprot
        mw
        agg
        switchType
        primaryReference { doi }
        references { doi }
        states {
            id
            name
            exMax
            emMax
            emhex
            exhex
            extCoeff
            qy
            lifetime
            spectra { id subtype data }
        }
        defaultState {
            id
            name
            exMax
            emMax
            emhex
            exhex
            extCoeff
            qy
            lifetime
            spectra { id subtype data }
         }
    }
}
"""

SPECTRUM_QUERY = """
query getSpectrum($id: Int!) {
    spectrum(id: $id) {
        id
        subtype
        data
        ownerFilter {
            id
            name
            manufacturer
            bandcenter
            bandwidth
            edge
            spectrum { id subtype data }
        }
        ownerCamera {
            id
            name
            manufacturer
            spectrum { id subtype data }
        }
        ownerLight {
            id
            name
            manufacturer
            spectrum { id subtype data }
        }
    }
}
"""
