"""Main fetching logic."""

from enum import Enum
from typing import Any, Literal, Optional

from pydantic import BaseModel, Field, field_validator, model_validator

__all__ = [
    "Filter",
    "FilterPlacement",
    "FilterSpectrum",
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

    def __str__(self) -> str:
        """Return the string representation of the enum."""
        return self.value

    def __repr__(self) -> str:
        """Return the repr of the enum."""
        return repr(self.value)


class Spectrum(BaseModel):
    """Spectrum with data."""

    subtype: SpectrumType
    data: list[tuple[float, float]] = Field(..., repr=False)


class SpectrumOwner(BaseModel):
    """Something that can own a spectrum."""

    name: str
    spectrum: Spectrum


class Filter(SpectrumOwner):
    """A filter with its properties."""

    manufacturer: str
    bandcenter: Optional[float]
    bandwidth: Optional[float]
    edge: Optional[float]


class FilterSpectrum(Spectrum):
    """Spectrum owned by a filter."""

    ownerFilter: Filter


class State(BaseModel):
    """Fluorophore state."""

    id: int
    exMax: float  # nanometers
    emMax: float  # nanometers
    emhex: str = ""
    exhex: str = ""
    extCoeff: Optional[float] = None  # M^-1 cm^-1
    qy: Optional[float] = None
    spectra: list[Spectrum]
    lifetime: Optional[float] = None  # ns

    @property
    def excitation_spectrum(self) -> Optional[Spectrum]:
        """Return the excitation spectrum, absorption spectrum, or None."""
        spect = next((s for s in self.spectra if s.subtype == "EX"), None)
        if not spect:
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
    states: list[State] = Field(default_factory=list)
    defaultState: Optional[int] = None

    @model_validator(mode="before")
    @classmethod
    def _v_model(cls, v: Any) -> Any:
        if isinstance(v, dict):
            out = dict(v)
            if "states" not in v and "exMax" in v:
                out["states"] = [State(**v)]
            return out
        return v

    @field_validator("defaultState", mode="before")
    @classmethod
    def _v_default_state(cls, v: Any) -> int:
        if isinstance(v, dict) and "id" in v:
            return int(v["id"])
        return int(v)

    @property
    def default_state(self) -> Optional[State]:
        """Return the default state or the first state."""
        for state in self.states:
            if state.id == self.defaultState:
                return state
        return next(iter(self.states), None)


class FilterPlacement(SpectrumOwner):
    """A filter placed in a microscope."""

    path: Literal["EX", "EM", "BS"]
    reflects: bool = False


class OpticalConfig(BaseModel):
    """A collection of filters and light sources."""

    name: str
    filters: list[FilterPlacement]
    camera: Optional[SpectrumOwner]
    light: Optional[SpectrumOwner]
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
    protein: Fluorophore


class ProteinResponse(BaseModel):
    """Response for a protein query."""

    data: _ProteinPayload


class _DyePayload(BaseModel):
    dye: Fluorophore


class DyeResponse(BaseModel):
    """Response for a dye query."""

    data: _DyePayload


class _FilterSpectrumPayload(BaseModel):
    spectrum: FilterSpectrum


class FilterSpectrumResponse(BaseModel):
    """Response for a filter spectrum query."""

    data: _FilterSpectrumPayload


# WIP
# def generate_graphql_query(model: type[BaseModel], model_name: str = "") -> str:
#     def get_fields(model: type[BaseModel]) -> str:
#         fields = []
#         for name, field in model.model_fields.items():
#             annotation = field.annotation

#             if isinstance(annotation, type) and issubclass(annotation, BaseModel):
#                 sub_fields = get_fields(annotation)
#                 fields.append(f"{name} {{ {sub_fields} }}")
#             elif (
#                 get_origin(annotation) in (list, tuple)
#                 and isinstance(type_ := get_args(annotation)[0], type)
#                 and issubclass(type_, BaseModel)
#             ):
#                 sub_fields = get_fields(type_)
#                 fields.append(f"{name} {{ {sub_fields} }}")
#             else:
#                 fields.append(name)
#         return "\n".join(fields)

#     fields_str = get_fields(model)
#     model_name = model_name or model.__name__
#     return f"""
#     query get{model_name}($id: String!) {{
#         {model_name.lower()}(id: $id) {{
#             {fields_str}
#         }}
#     }}
#     """
