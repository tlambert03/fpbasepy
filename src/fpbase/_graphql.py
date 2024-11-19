MICROSCOPE_QUERY = """
query getMicroscope($id: String!) {
    microscope(id: $id) {
        id
        name
        opticalConfigs {
            name
            filters {
                name
                path
                reflects
                spectrum { subtype data }
            }
            camera { name spectrum { subtype data } }
            light { name spectrum { subtype data } }
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
        spectra { subtype data }
    }
}
"""

PROTEIN_QUERY = """
query getProtein($id: String!) {
    protein(id: $id) {
        name
        id
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
            spectra { subtype data }
        }
        defaultState {
            id
        }
    }
}
"""

# slightly convuluted way to get to filter, since the API
# doesn't offer a top-level filter query
SPECTRUM_QUERY = """
query getSpectrum($id: Int!) {
    spectrum(id: $id) {
        subtype
        data
        ownerFilter {
            name
            manufacturer
            bandcenter
            bandwidth
            edge
            spectrum { subtype data }
        }
    }
}
"""
