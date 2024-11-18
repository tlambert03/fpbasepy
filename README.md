# Python FPbase API

[![License](https://img.shields.io/pypi/l/fpbase.svg?color=green)](https://github.com/tlambert03/fpbasepy/raw/main/LICENSE)
[![PyPI](https://img.shields.io/pypi/v/fpbase.svg?color=green)](https://pypi.org/project/fpbase)
[![Python Version](https://img.shields.io/pypi/pyversions/fpbase.svg?color=green)](https://python.org)
[![CI](https://github.com/tlambert03/fpbasepy/actions/workflows/ci.yml/badge.svg)](https://github.com/tlambert03/fpbasepy/actions/workflows/ci.yml)
[![codecov](https://codecov.io/gh/tlambert03/fpbasepy/branch/main/graph/badge.svg)](https://codecov.io/gh/tlambert03/fpbasepy)

Python wrapper for FPBase.org GraphQL API.

See https://www.fpbase.org/graphql for full documentation on the graphql schema and an interactive playground.

This library provides simple Python access to commonly-accessed data.

## Installation

```
pip install fpbase
```

## Usage

```python
In [1]: from fpbase import get_fluorophore, get_microscope

In [2]: print(get_fluorophore("mCherry"))
Fluorophore(
    name='mCherry',
    id='ZERB6',
    states=[
        State(
            id=336,
            exMax=587.0,
            emMax=610.0,
            emhex='#f70000',
            exhex='#ff4600',
            extCoeff=72000.0,
            qy=0.22,
            spectra=[Spectrum(subtype='EX'), Spectrum(subtype='EM'), Spectrum(subtype='A_2P')],
            lifetime=1.4
        )
    ],
    defaultState=336
)

In [3]: print(get_fluorophore("DAPI"))
Fluorophore(
    name='DAPI',
    id='15',
    states=[
        State(
            id=15,
            exMax=359.0,
            emMax=461.0,
            emhex='',
            exhex='',
            extCoeff=None,
            qy=None,
            spectra=[Spectrum(subtype='AB'), Spectrum(subtype='EX'), Spectrum(subtype='EM')],
            lifetime=None
        )
    ],
    defaultState=None
)

# fetch info for https://www.fpbase.org/microscope/i6WL2W/
In [4]: print(get_microscope("i6WL2W"))
Microscope(
    id='i6WL2WdgcDMgJYtPrpZcaJ',
    name='Example Widefield (Sedat)',
    opticalConfigs=[
        OpticalConfig(
            name='Widefield Blue',
            filters=[
                FilterPlacement(name='Chroma ET395/25x', spectrum=Spectrum(subtype='BX'), path='EX', reflects=False),
                FilterPlacement(name='Chroma T425lpxr', spectrum=Spectrum(subtype='LP'), path='BS', reflects=False),
                FilterPlacement(name='Chroma ET460/50m', spectrum=Spectrum(subtype='BM'), path='EM', reflects=False)
            ],
            camera=SpectrumOwner(name='Andor Zyla 4.2 PLUS', spectrum=Spectrum(subtype='QE')),
            light=SpectrumOwner(name='SOLA 395', spectrum=Spectrum(subtype='PD')),
            laser=None
        ),
        OpticalConfig(
            name='Widefield Dual FRET',
            filters=[
                FilterPlacement(name='Lumencor 470/24x', spectrum=Spectrum(subtype='BX'), path='EX', reflects=False),
                FilterPlacement(name='Chroma 59022bs', spectrum=Spectrum(subtype='BS'), path='BS', reflects=False),
                FilterPlacement(name='Semrock FF02-641/75', spectrum=Spectrum(subtype='BP'), path='EM', reflects=False)
            ],
            camera=SpectrumOwner(name='Andor Zyla 4.2 PLUS', spectrum=Spectrum(subtype='QE')),
            light=SpectrumOwner(name='SOLA 395', spectrum=Spectrum(subtype='PD')),
            laser=None
        ),
        OpticalConfig(
            name='Widefield Dual Green',
            filters=[
                FilterPlacement(name='Lumencor 470/24x', spectrum=Spectrum(subtype='BX'), path='EX', reflects=False),
                FilterPlacement(name='Chroma 59022bs', spectrum=Spectrum(subtype='BS'), path='BS', reflects=False),
                FilterPlacement(name='Semrock FF03-525/50', spectrum=Spectrum(subtype='BP'), path='EM', reflects=False)
            ],
            camera=SpectrumOwner(name='Andor Zyla 4.2 PLUS', spectrum=Spectrum(subtype='QE')),
            light=SpectrumOwner(name='SOLA 395', spectrum=Spectrum(subtype='PD')),
            laser=None
        ),
        OpticalConfig(
            name='Widefield Dual Red',
            filters=[
                FilterPlacement(name='Lumencor 575/25x', spectrum=Spectrum(subtype='BX'), path='EX', reflects=False),
                FilterPlacement(name='Chroma 59022bs', spectrum=Spectrum(subtype='BS'), path='BS', reflects=False),
                FilterPlacement(name='Semrock FF02-641/75', spectrum=Spectrum(subtype='BP'), path='EM', reflects=False)
            ],
            camera=SpectrumOwner(name='Andor Zyla 4.2 PLUS', spectrum=Spectrum(subtype='QE')),
            light=SpectrumOwner(name='SOLA 395', spectrum=Spectrum(subtype='PD')),
            laser=None
        ),
        OpticalConfig(
            name='Widefield Far-Red',
            filters=[
                FilterPlacement(name='Chroma ET640/30x', spectrum=Spectrum(subtype='BX'), path='EX', reflects=False),
                FilterPlacement(name='Chroma T660lpxr', spectrum=Spectrum(subtype='LP'), path='BS', reflects=False),
                FilterPlacement(name='Semrock FF01-698/70', spectrum=Spectrum(subtype='BP'), path='EM', reflects=False)
            ],
            camera=SpectrumOwner(name='Andor Zyla 4.2 PLUS', spectrum=Spectrum(subtype='QE')),
            light=SpectrumOwner(name='SOLA 395', spectrum=Spectrum(subtype='PD')),
            laser=None
        ),
        OpticalConfig(
            name='Widefield Triple Cyan',
            filters=[
                FilterPlacement(name='Lumencor 440/20x', spectrum=Spectrum(subtype='BX'), path='EX', reflects=False),
                FilterPlacement(name='Chroma 69008bs', spectrum=Spectrum(subtype='BS'), path='BS', reflects=False),
                FilterPlacement(name='Chroma ET470/24m', spectrum=Spectrum(subtype='BM'), path='EM', reflects=False)
            ],
            camera=SpectrumOwner(name='Andor Zyla 4.2 PLUS', spectrum=Spectrum(subtype='QE')),
            light=SpectrumOwner(name='SOLA 395', spectrum=Spectrum(subtype='PD')),
            laser=None
        ),
        OpticalConfig(
            name='Widefield Triple FRET',
            filters=[
                FilterPlacement(name='Lumencor 440/20x', spectrum=Spectrum(subtype='BX'), path='EX', reflects=False),
                FilterPlacement(name='Chroma 69008bs', spectrum=Spectrum(subtype='BS'), path='BS', reflects=False),
                FilterPlacement(name='Chroma ET535/30m', spectrum=Spectrum(subtype='BM'), path='EM', reflects=False)
            ],
            camera=SpectrumOwner(name='Andor Zyla 4.2 PLUS', spectrum=Spectrum(subtype='QE')),
            light=SpectrumOwner(name='SOLA 395', spectrum=Spectrum(subtype='PD')),
            laser=None
        ),
        OpticalConfig(
            name='Widefield Triple Red',
            filters=[
                FilterPlacement(name='Lumencor 575/25x', spectrum=Spectrum(subtype='BX'), path='EX', reflects=False),
                FilterPlacement(name='Chroma 69008bs', spectrum=Spectrum(subtype='BS'), path='BS', reflects=False),
                FilterPlacement(name='Semrock FF02-641/75', spectrum=Spectrum(subtype='BP'), path='EM', reflects=False)
            ],
            camera=SpectrumOwner(name='Andor Zyla 4.2 PLUS', spectrum=Spectrum(subtype='QE')),
            light=SpectrumOwner(name='SOLA 395', spectrum=Spectrum(subtype='PD')),
            laser=None
        ),
        OpticalConfig(
            name='Widefield Triple Yellow',
            filters=[
                FilterPlacement(name='Chroma ET500/20x', spectrum=Spectrum(subtype='BX'), path='EX', reflects=False),
                FilterPlacement(name='Chroma 69008bs', spectrum=Spectrum(subtype='BS'), path='BS', reflects=False),
                FilterPlacement(name='Chroma ET535/30m', spectrum=Spectrum(subtype='BM'), path='EM', reflects=False)
            ],
            camera=SpectrumOwner(name='Andor Zyla 4.2 PLUS', spectrum=Spectrum(subtype='QE')),
            light=SpectrumOwner(name='SOLA 395', spectrum=Spectrum(subtype='PD')),
            laser=None
        )
    ]
)
```
