# Python FPbase API

[![License](https://img.shields.io/pypi/l/fpbase.svg?color=green)](https://github.com/tlambert03/fpbasepy/raw/main/LICENSE)
[![PyPI](https://img.shields.io/pypi/v/fpbase.svg?color=green)](https://pypi.org/project/fpbase)
[![Python Version](https://img.shields.io/pypi/pyversions/fpbase.svg?color=green)](https://python.org)
[![CI](https://github.com/tlambert03/fpbasepy/actions/workflows/ci.yml/badge.svg)](https://github.com/tlambert03/fpbasepy/actions/workflows/ci.yml)
[![codecov](https://codecov.io/gh/tlambert03/fpbasepy/branch/main/graph/badge.svg)](https://codecov.io/gh/tlambert03/fpbasepy)

Python wrapper for FPBase.org GraphQL API.

See <https://www.fpbase.org/graphql> for an interactive playground and full documentation on the graphql schema
(click the little book icon on the top left corner).

## Installation

```sh
pip install fpbase
```

## API

See all response model types in `fpbase.models`.

### Functions that return a single object

* `fpbase.get_fluorophore`  (can return either protein or dye)
* `fpbase.get_protein` (will only return proteins)
* `fpbase.get_microscope`
* `fpbase.get_filter`
* `fpbase.get_light_source`
* `fpbase.get_camera`

### Functions that return a list of available object keys

* `fpbase.list_fluorophores`  (includes both proteins and dyes)
* `fpbase.list_proteins`
* `fpbase.list_dyes`
* `fpbase.list_microscopes`
* `fpbase.list_filters`
* `fpbase.list_light_sources`
* `fpbase.list_cameras`

### Other

* `fpbase.graphql_query` : Send generic GraphQL query to FPbase (see <https://www.fpbase.org/graphql> for full documentation on the graphql schema and an interactive playground)

## Example Usage

```python
import fpbase

fpbase.get_fluorophore("mCherry")
# NOTE: you can also use `get_protein()` if you want to ensure that the
# fluorophore is a protein and the response is an `fpbase.models.Protein`
```

<details>

<summary>output</summary>

```python
Protein(
    name='mCherry',
    id='ZERB6',
    default_state=State(
        id=336,
        name='mCherry',
        exMax=587.0,
        emMax=610.0,
        emhex='#f70000',
        exhex='#ff4600',
        ext_coeff=72000.0,
        qy=0.22,
        spectra=[
            Spectrum(id=79, subtype='EX', owner_filter=None, owner_camera=None, owner_light=None),
            Spectrum(id=80, subtype='EM', owner_filter=None, owner_camera=None, owner_light=None),
            Spectrum(id=158, subtype='A_2P', owner_filter=None, owner_camera=None, owner_light=None)
        ],
        lifetime=1.4
    ),
    seq='MVSKGEEDNMAIIKEFMRFKVHMEGSVNGHEFEIEGEGEGRPYEGTQTAKLKVTKGGPLPFAWDILSPQFMYGSKAYVKHPADIPDYLKLSFPEGFKWERVMNFEDGGVVTVTQDSSLQDGEFIYKVKLRGTNFPSDGPVMQKKTMGWEASSERMYPEDGALKGEIKQRLKLKDGGHYDAEVKTTYKAKKPVQLPGAYNVNIKLDITSHNEDYTIVEQYERAEGRHSTGGMDELYK',
    pdb=['2H5Q'],
    genbank='AAV52164',
    uniprot='X5DSL3',
    agg=<Olig.MONOMER: 'M'>,
    switch_type=<SwitchType.BASIC: 'B'>,
    primary_reference=Reference(doi='10.1038/nbt1037', url='https://doi.org/10.1038/nbt1037'),
    references=[
        Reference(doi='10.1038/nbt1037', url='https://doi.org/10.1038/nbt1037'),
        Reference(doi='10.1021/bi060773l', url='https://doi.org/10.1021/bi060773l'),
        Reference(doi='10.1038/nmeth1062', url='https://doi.org/10.1038/nmeth1062'),
        Reference(doi='10.1038/nmeth.1596', url='https://doi.org/10.1038/nmeth.1596'),
        Reference(doi='10.1038/nmeth.1955', url='https://doi.org/10.1038/nmeth.1955'),
        Reference(doi='10.1091/mbc.e16-01-0063', url='https://doi.org/10.1091/mbc.e16-01-0063'),
        Reference(doi='10.1038/s41598-017-12212-x', url='https://doi.org/10.1038/s41598-017-12212-x'),
        Reference(doi='10.1038/s41598-018-28858-0', url='https://doi.org/10.1038/s41598-018-28858-0'),
        Reference(doi='10.1002/pld3.112', url='https://doi.org/10.1002/pld3.112'),
        Reference(doi='10.1371/journal.pone.0219886', url='https://doi.org/10.1371/journal.pone.0219886'),
        Reference(doi='10.1073/pnas.2017379117', url='https://doi.org/10.1073/pnas.2017379117'),
        Reference(doi='10.1038/s41467-023-40647-6', url='https://doi.org/10.1038/s41467-023-40647-6')
    ]
)
```

</details>

```python
fpbase.get_fluorophore("DAPI")
```

<details>

<summary>output</summary>

```python
Fluorophore(
    name='DAPI',
    id='15',
    default_state=State(
        id=15,
        name='DAPI',
        exMax=359.0,
        emMax=461.0,
        emhex='',
        exhex='',
        ext_coeff=None,
        qy=None,
        spectra=[
            Spectrum(id=7754, subtype='AB', owner_filter=None, owner_camera=None, owner_light=None),
            Spectrum(id=222, subtype='EX', owner_filter=None, owner_camera=None, owner_light=None),
            Spectrum(id=223, subtype='EM', owner_filter=None, owner_camera=None, owner_light=None)
        ],
        lifetime=None
    )
)
```

</details>

```python
# fetch info for <https://www.fpbase.org/microscope/i6WL2W/>
fpbase.get_microscope("i6WL2W")
```

<details>

<summary>output</summary>

```python
Microscope(
    id='i6WL2WdgcDMgJYtPrpZcaJ',
    name='Example Widefield (Sedat)',
    opticalConfigs=[
        OpticalConfig(
            name='Widefield Blue',
            filters=[
                FilterPlacement(
                    path='EX',
                    filter=Filter(
                        id=52,
                        name='Chroma ET395/25x',
                        spectrum=Spectrum(id=375, subtype='BX', owner_filter=None, owner_camera=None, owner_light=None),
                        manufacturer='',
                        bandcenter=None,
                        bandwidth=None,
                        edge=None
                    ),
                    reflects=False
                ),
                FilterPlacement(
                    path='BS',
                    filter=Filter(
                        id=10,
                        name='Chroma T425lpxr',
                        spectrum=Spectrum(id=333, subtype='LP', owner_filter=None, owner_camera=None, owner_light=None),
                        manufacturer='',
                        bandcenter=None,
                        bandwidth=None,
                        edge=None
                    ),
                    reflects=False
                ),
                FilterPlacement(
                    path='EM',
                    filter=Filter(
                        id=47,
                        name='Chroma ET460/50m',
                        spectrum=Spectrum(id=370, subtype='BM', owner_filter=None, owner_camera=None, owner_light=None),
                        manufacturer='',
                        bandcenter=None,
                        bandwidth=None,
                        edge=None
                    ),
                    reflects=False
                )
            ],
            camera=Camera(
                id=4,
                name='Andor Zyla 4.2 PLUS',
                spectrum=Spectrum(id=1328, subtype='QE', owner_filter=None, owner_camera=None, owner_light=None),
                manufacturer=''
            ),
            light=LightSource(
                id=9,
                name='SOLA 395',
                spectrum=Spectrum(id=394, subtype='PD', owner_filter=None, owner_camera=None, owner_light=None),
                manufacturer=''
            ),
            laser=None
        ),
        OpticalConfig(
            name='Widefield Dual FRET',
            filters=[
                FilterPlacement(
                    path='EX',
                    filter=Filter(
                        id=63,
                        name='Lumencor 470/24x',
                        spectrum=Spectrum(id=399, subtype='BX', owner_filter=None, owner_camera=None, owner_light=None),
                        manufacturer='',
                        bandcenter=None,
                        bandwidth=None,
                        edge=None
                    ),
                    reflects=False
                ),
                FilterPlacement(
                    path='BS',
                    filter=Filter(
                        id=62,
                        name='Chroma 59022bs',
                        spectrum=Spectrum(id=385, subtype='BS', owner_filter=None, owner_camera=None, owner_light=None),
                        manufacturer='',
                        bandcenter=None,
                        bandwidth=None,
                        edge=None
                    ),
                    reflects=False
                ),
                FilterPlacement(
                    path='EM',
                    filter=Filter(
                        id=689,
                        name='Semrock FF02-641/75',
                        spectrum=Spectrum(id=1025, subtype='BP', owner_filter=None, owner_camera=None, owner_light=None),
                        manufacturer='',
                        bandcenter=None,
                        bandwidth=None,
                        edge=None
                    ),
                    reflects=False
                )
            ],
            camera=Camera(
                id=4,
                name='Andor Zyla 4.2 PLUS',
                spectrum=Spectrum(id=1328, subtype='QE', owner_filter=None, owner_camera=None, owner_light=None),
                manufacturer=''
            ),
            light=LightSource(
                id=9,
                name='SOLA 395',
                spectrum=Spectrum(id=394, subtype='PD', owner_filter=None, owner_camera=None, owner_light=None),
                manufacturer=''
            ),
            laser=None
        ),
        OpticalConfig(
            name='Widefield Dual Green',
            filters=[
                FilterPlacement(
                    path='EX',
                    filter=Filter(
                        id=63,
                        name='Lumencor 470/24x',
                        spectrum=Spectrum(id=399, subtype='BX', owner_filter=None, owner_camera=None, owner_light=None),
                        manufacturer='',
                        bandcenter=None,
                        bandwidth=None,
                        edge=None
                    ),
                    reflects=False
                ),
                FilterPlacement(
                    path='BS',
                    filter=Filter(
                        id=62,
                        name='Chroma 59022bs',
                        spectrum=Spectrum(id=385, subtype='BS', owner_filter=None, owner_camera=None, owner_light=None),
                        manufacturer='',
                        bandcenter=None,
                        bandwidth=None,
                        edge=None
                    ),
                    reflects=False
                ),
                FilterPlacement(
                    path='EM',
                    filter=Filter(
                        id=569,
                        name='Semrock FF03-525/50',
                        spectrum=Spectrum(id=905, subtype='BP', owner_filter=None, owner_camera=None, owner_light=None),
                        manufacturer='',
                        bandcenter=None,
                        bandwidth=None,
                        edge=None
                    ),
                    reflects=False
                )
            ],
            camera=Camera(
                id=4,
                name='Andor Zyla 4.2 PLUS',
                spectrum=Spectrum(id=1328, subtype='QE', owner_filter=None, owner_camera=None, owner_light=None),
                manufacturer=''
            ),
            light=LightSource(
                id=9,
                name='SOLA 395',
                spectrum=Spectrum(id=394, subtype='PD', owner_filter=None, owner_camera=None, owner_light=None),
                manufacturer=''
            ),
            laser=None
        ),
        OpticalConfig(
            name='Widefield Dual Red',
            filters=[
                FilterPlacement(
                    path='EX',
                    filter=Filter(
                        id=67,
                        name='Lumencor 575/25x',
                        spectrum=Spectrum(id=403, subtype='BX', owner_filter=None, owner_camera=None, owner_light=None),
                        manufacturer='',
                        bandcenter=None,
                        bandwidth=None,
                        edge=None
                    ),
                    reflects=False
                ),
                FilterPlacement(
                    path='BS',
                    filter=Filter(
                        id=62,
                        name='Chroma 59022bs',
                        spectrum=Spectrum(id=385, subtype='BS', owner_filter=None, owner_camera=None, owner_light=None),
                        manufacturer='',
                        bandcenter=None,
                        bandwidth=None,
                        edge=None
                    ),
                    reflects=False
                ),
                FilterPlacement(
                    path='EM',
                    filter=Filter(
                        id=689,
                        name='Semrock FF02-641/75',
                        spectrum=Spectrum(id=1025, subtype='BP', owner_filter=None, owner_camera=None, owner_light=None),
                        manufacturer='',
                        bandcenter=None,
                        bandwidth=None,
                        edge=None
                    ),
                    reflects=False
                )
            ],
            camera=Camera(
                id=4,
                name='Andor Zyla 4.2 PLUS',
                spectrum=Spectrum(id=1328, subtype='QE', owner_filter=None, owner_camera=None, owner_light=None),
                manufacturer=''
            ),
            light=LightSource(
                id=9,
                name='SOLA 395',
                spectrum=Spectrum(id=394, subtype='PD', owner_filter=None, owner_camera=None, owner_light=None),
                manufacturer=''
            ),
            laser=None
        ),
        OpticalConfig(
            name='Widefield Far-Red',
            filters=[
                FilterPlacement(
                    path='EX',
                    filter=Filter(
                        id=445,
                        name='Chroma ET640/30x',
                        spectrum=Spectrum(id=781, subtype='BX', owner_filter=None, owner_camera=None, owner_light=None),
                        manufacturer='',
                        bandcenter=None,
                        bandwidth=None,
                        edge=None
                    ),
                    reflects=False
                ),
                FilterPlacement(
                    path='BS',
                    filter=Filter(
                        id=6,
                        name='Chroma T660lpxr',
                        spectrum=Spectrum(id=329, subtype='LP', owner_filter=None, owner_camera=None, owner_light=None),
                        manufacturer='',
                        bandcenter=None,
                        bandwidth=None,
                        edge=None
                    ),
                    reflects=False
                ),
                FilterPlacement(
                    path='EM',
                    filter=Filter(
                        id=719,
                        name='Semrock FF01-698/70',
                        spectrum=Spectrum(id=1055, subtype='BP', owner_filter=None, owner_camera=None, owner_light=None),
                        manufacturer='',
                        bandcenter=None,
                        bandwidth=None,
                        edge=None
                    ),
                    reflects=False
                )
            ],
            camera=Camera(
                id=4,
                name='Andor Zyla 4.2 PLUS',
                spectrum=Spectrum(id=1328, subtype='QE', owner_filter=None, owner_camera=None, owner_light=None),
                manufacturer=''
            ),
            light=LightSource(
                id=9,
                name='SOLA 395',
                spectrum=Spectrum(id=394, subtype='PD', owner_filter=None, owner_camera=None, owner_light=None),
                manufacturer=''
            ),
            laser=None
        ),
        OpticalConfig(
            name='Widefield Triple Cyan',
            filters=[
                FilterPlacement(
                    path='EX',
                    filter=Filter(
                        id=79,
                        name='Lumencor 440/20x',
                        spectrum=Spectrum(id=415, subtype='BX', owner_filter=None, owner_camera=None, owner_light=None),
                        manufacturer='',
                        bandcenter=None,
                        bandwidth=None,
                        edge=None
                    ),
                    reflects=False
                ),
                FilterPlacement(
                    path='BS',
                    filter=Filter(
                        id=60,
                        name='Chroma 69008bs',
                        spectrum=Spectrum(id=383, subtype='BS', owner_filter=None, owner_camera=None, owner_light=None),
                        manufacturer='',
                        bandcenter=None,
                        bandwidth=None,
                        edge=None
                    ),
                    reflects=False
                ),
                FilterPlacement(
                    path='EM',
                    filter=Filter(
                        id=46,
                        name='Chroma ET470/24m',
                        spectrum=Spectrum(id=369, subtype='BM', owner_filter=None, owner_camera=None, owner_light=None),
                        manufacturer='',
                        bandcenter=None,
                        bandwidth=None,
                        edge=None
                    ),
                    reflects=False
                )
            ],
            camera=Camera(
                id=4,
                name='Andor Zyla 4.2 PLUS',
                spectrum=Spectrum(id=1328, subtype='QE', owner_filter=None, owner_camera=None, owner_light=None),
                manufacturer=''
            ),
            light=LightSource(
                id=9,
                name='SOLA 395',
                spectrum=Spectrum(id=394, subtype='PD', owner_filter=None, owner_camera=None, owner_light=None),
                manufacturer=''
            ),
            laser=None
        ),
        OpticalConfig(
            name='Widefield Triple FRET',
            filters=[
                FilterPlacement(
                    path='EX',
                    filter=Filter(
                        id=79,
                        name='Lumencor 440/20x',
                        spectrum=Spectrum(id=415, subtype='BX', owner_filter=None, owner_camera=None, owner_light=None),
                        manufacturer='',
                        bandcenter=None,
                        bandwidth=None,
                        edge=None
                    ),
                    reflects=False
                ),
                FilterPlacement(
                    path='BS',
                    filter=Filter(
                        id=60,
                        name='Chroma 69008bs',
                        spectrum=Spectrum(id=383, subtype='BS', owner_filter=None, owner_camera=None, owner_light=None),
                        manufacturer='',
                        bandcenter=None,
                        bandwidth=None,
                        edge=None
                    ),
                    reflects=False
                ),
                FilterPlacement(
                    path='EM',
                    filter=Filter(
                        id=36,
                        name='Chroma ET535/30m',
                        spectrum=Spectrum(id=359, subtype='BM', owner_filter=None, owner_camera=None, owner_light=None),
                        manufacturer='',
                        bandcenter=None,
                        bandwidth=None,
                        edge=None
                    ),
                    reflects=False
                )
            ],
            camera=Camera(
                id=4,
                name='Andor Zyla 4.2 PLUS',
                spectrum=Spectrum(id=1328, subtype='QE', owner_filter=None, owner_camera=None, owner_light=None),
                manufacturer=''
            ),
            light=LightSource(
                id=9,
                name='SOLA 395',
                spectrum=Spectrum(id=394, subtype='PD', owner_filter=None, owner_camera=None, owner_light=None),
                manufacturer=''
            ),
            laser=None
        ),
        OpticalConfig(
            name='Widefield Triple Red',
            filters=[
                FilterPlacement(
                    path='EX',
                    filter=Filter(
                        id=67,
                        name='Lumencor 575/25x',
                        spectrum=Spectrum(id=403, subtype='BX', owner_filter=None, owner_camera=None, owner_light=None),
                        manufacturer='',
                        bandcenter=None,
                        bandwidth=None,
                        edge=None
                    ),
                    reflects=False
                ),
                FilterPlacement(
                    path='BS',
                    filter=Filter(
                        id=60,
                        name='Chroma 69008bs',
                        spectrum=Spectrum(id=383, subtype='BS', owner_filter=None, owner_camera=None, owner_light=None),
                        manufacturer='',
                        bandcenter=None,
                        bandwidth=None,
                        edge=None
                    ),
                    reflects=False
                ),
                FilterPlacement(
                    path='EM',
                    filter=Filter(
                        id=689,
                        name='Semrock FF02-641/75',
                        spectrum=Spectrum(id=1025, subtype='BP', owner_filter=None, owner_camera=None, owner_light=None),
                        manufacturer='',
                        bandcenter=None,
                        bandwidth=None,
                        edge=None
                    ),
                    reflects=False
                )
            ],
            camera=Camera(
                id=4,
                name='Andor Zyla 4.2 PLUS',
                spectrum=Spectrum(id=1328, subtype='QE', owner_filter=None, owner_camera=None, owner_light=None),
                manufacturer=''
            ),
            light=LightSource(
                id=9,
                name='SOLA 395',
                spectrum=Spectrum(id=394, subtype='PD', owner_filter=None, owner_camera=None, owner_light=None),
                manufacturer=''
            ),
            laser=None
        ),
        OpticalConfig(
            name='Widefield Triple Yellow',
            filters=[
                FilterPlacement(
                    path='EX',
                    filter=Filter(
                        id=41,
                        name='Chroma ET500/20x',
                        spectrum=Spectrum(id=364, subtype='BX', owner_filter=None, owner_camera=None, owner_light=None),
                        manufacturer='',
                        bandcenter=None,
                        bandwidth=None,
                        edge=None
                    ),
                    reflects=False
                ),
                FilterPlacement(
                    path='BS',
                    filter=Filter(
                        id=60,
                        name='Chroma 69008bs',
                        spectrum=Spectrum(id=383, subtype='BS', owner_filter=None, owner_camera=None, owner_light=None),
                        manufacturer='',
                        bandcenter=None,
                        bandwidth=None,
                        edge=None
                    ),
                    reflects=False
                ),
                FilterPlacement(
                    path='EM',
                    filter=Filter(
                        id=36,
                        name='Chroma ET535/30m',
                        spectrum=Spectrum(id=359, subtype='BM', owner_filter=None, owner_camera=None, owner_light=None),
                        manufacturer='',
                        bandcenter=None,
                        bandwidth=None,
                        edge=None
                    ),
                    reflects=False
                )
            ],
            camera=Camera(
                id=4,
                name='Andor Zyla 4.2 PLUS',
                spectrum=Spectrum(id=1328, subtype='QE', owner_filter=None, owner_camera=None, owner_light=None),
                manufacturer=''
            ),
            light=LightSource(
                id=9,
                name='SOLA 395',
                spectrum=Spectrum(id=394, subtype='PD', owner_filter=None, owner_camera=None, owner_light=None),
                manufacturer=''
            ),
            laser=None
        )
    ]
)
```

</details>
