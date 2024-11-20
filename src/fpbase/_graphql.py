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
