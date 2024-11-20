"""Main fetching logic."""

from enum import Enum
from typing import TYPE_CHECKING, Any, Optional

from pydantic import BaseModel, Field, computed_field, model_validator

if TYPE_CHECKING:
    from collections.abc import Iterable

__all__ = [
    "Filter",
    "FilterPlacement",
    "Fluorophore",
    "Microscope",
    "OpticalConfig",
    "Spectrum",
    "SpectrumOwner",
    "SpectrumType",
    "State",
]


class SpectrumType(str, Enum):
    """Spectrum types."""

    A_2P = "A_2P"
    BM = "BM"
    BP = "BP"
    BS = "BS"
    BX = "BX"
    EM = "EM"
    EX = "EX"
    LP = "LP"
    PD = "PD"
    QE = "QE"
    AB = "AB"

    def __str__(self) -> str:  # pragma: no cover
        """Return the string representation of the enum."""
        return self.value

    def __repr__(self) -> str:
        """Return the repr of the enum."""
        return repr(self.value)


class FilterPath(str, Enum):
    """Placement of a filter in an optical config."""

    EX = "EX"
    EM = "EM"
    BS = "BS"

    def __str__(self) -> str:  # pragma: no cover
        """Return the string representation of the enum."""
        return self.value

    def __repr__(self) -> str:
        """Return the repr of the enum."""
        return repr(self.value)


class Olig(str, Enum):
    MONOMER = "M"
    DIMER = "D"
    TANDEM_DIMER = "TD"
    WEAK_DIMER = "WD"
    TETRAMER = "T"

    def __str__(self) -> str:  # pragma: no cover
        """Return the string representation of the enum."""
        return self.value


class SwitchType(str, Enum):
    BASIC = "B"
    PHOTOACTIVATABLE = "PA"
    PHOTOSWITCHABLE = "PS"
    PHOTOCONVERTIBLE = "PC"
    MULTIPHOTOCHROMIC = "MP"
    TIMER = "T"
    OTHER = "O"

    def __str__(self) -> str:  # pragma: no cover
        """Return the string representation of the enum."""
        return self.value


class Spectrum(BaseModel):
    """Spectrum with data."""

    id: int
    subtype: SpectrumType
    data: list[tuple[float, float]] = Field(..., repr=False)

    owner_filter: Optional["Filter"] = Field(None, alias="ownerFilter")
    owner_camera: Optional["Camera"] = Field(None, alias="ownerCamera")
    owner_light: Optional["LightSource"] = Field(None, alias="ownerLight")


class SpectrumOwner(BaseModel):
    """Something that can own a spectrum."""

    id: int
    name: str
    spectrum: Spectrum


class Filter(SpectrumOwner):
    """A filter with its properties."""

    manufacturer: str = ""
    bandcenter: Optional[float] = None
    bandwidth: Optional[float] = None
    edge: Optional[float] = None


class Camera(SpectrumOwner):
    manufacturer: str = ""


class LightSource(SpectrumOwner):
    manufacturer: str = ""


class State(BaseModel):
    """Fluorophore state."""

    id: int
    name: str
    exMax: Optional[float] = None  # nanometers
    emMax: Optional[float] = None  # nanometers
    emhex: str = ""
    exhex: str = ""
    ext_coeff: Optional[float] = Field(None, alias="extCoeff")  # M^-1 cm^-1
    qy: Optional[float] = None
    spectra: list[Spectrum] = Field(default_factory=list)
    lifetime: Optional[float] = None  # ns

    @property
    def excitation_spectrum(self) -> Optional[Spectrum]:
        """Return the excitation spectrum, absorption spectrum, or None."""
        spect = next((s for s in self.spectra if s.subtype == "EX"), None)
        if not spect:  # pragma: no cover
            spect = next((s for s in self.spectra if s.subtype == "AB"), None)
        return spect

    @property
    def emission_spectrum(self) -> Optional[Spectrum]:
        """Return the emission spectrum or None."""
        return next((s for s in self.spectra if s.subtype == "EM"), None)


class Fluorophore(BaseModel):
    """A fluorophore with its states."""

    name: str
    id: str
    default_state: Optional[State] = Field(None, alias="defaultState")
    states: list[State] = Field(default_factory=list)

    @model_validator(mode="before")
    @classmethod
    def _v_model(cls, v: Any) -> Any:
        if isinstance(v, dict):
            out = dict(v)
            if "states" not in v and "exMax" in v:
                # this is a single-state fluorophore. probably a Dye.
                # this is a bit of a hack around the fpbase API
                state = State(**v)
                out["states"] = [state]
                out["defaultState"] = state
            return out
        return v  # pragma: no cover

    def __repr_args__(self) -> "Iterable[tuple[str | None, Any]]":
        """Return the repr args, excluding the default state if it's the only one."""
        for key, val in super().__repr_args__():
            if key == "states" and len(val) == 1 and val[0] == self.default_state:
                continue
            yield key, val


class Reference(BaseModel):
    doi: str

    @computed_field
    def url(self) -> str:
        """Return the DOI URL."""
        return f"https://doi.org/{self.doi}"


class Protein(Fluorophore):
    seq: Optional[str] = None
    pdb: list[str] = Field(default_factory=list)
    genbank: Optional[str] = None
    uniprot: Optional[str] = None
    agg: Optional[Olig] = None
    switch_type: Optional[SwitchType] = Field(None, alias="switchType")
    primary_reference: Optional[Reference] = Field(None, alias="primaryReference")
    references: list[Reference] = Field(default_factory=list)
    states: list[State] = Field(default_factory=list)
    # default_state: Optional[State] = Field(None, alias="defaultState")


class FilterPlacement(BaseModel):
    """A filter placed in a microscope."""

    path: FilterPath
    filter: Filter
    reflects: bool = False


class OpticalConfig(BaseModel):
    """A collection of filters and light sources."""

    name: str
    filters: list[FilterPlacement]
    camera: Optional["Camera"]
    light: Optional["LightSource"]
    laser: Optional[int]


class Microscope(BaseModel):
    """A microscope with its optical configurations."""

    id: str
    name: str
    opticalConfigs: list[OpticalConfig]


class _MicroscopePayload(BaseModel):
    microscope: Microscope


class MicroscopeResponse(BaseModel):
    """Response for a microscope query."""

    data: _MicroscopePayload


class _ProteinPayload(BaseModel):
    protein: Protein


class ProteinResponse(BaseModel):
    """Response for a protein query."""

    data: _ProteinPayload


class _DyePayload(BaseModel):
    dye: Fluorophore


class DyeResponse(BaseModel):
    """Response for a dye query."""

    data: _DyePayload


class _SpectrumPayload(BaseModel):
    spectrum: Spectrum


class SpectrumResponse(BaseModel):
    """Response for a filter spectrum query."""

    data: _SpectrumPayload
