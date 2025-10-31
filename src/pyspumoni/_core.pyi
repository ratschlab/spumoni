from __future__ import annotations
import typing

__all__: list[str] = ['Alignment', 'Index']

class Alignment:
    @typing.overload
    def __init__(self) -> None:
        ...
    @typing.overload
    def __init__(self, header: str, fwd: bool, start: typing.SupportsInt, pres_frac: typing.SupportsFloat, qry_len: typing.SupportsInt) -> None:
        ...
    @property
    def ctg(self) -> str:
        ...
    @property
    def pres_frac(self) -> float:
        ...
    @property
    def r_en(self) -> int:
        ...
    @property
    def r_st(self) -> int:
        ...
    @property
    def strand(self) -> int:
        ...
class Index:
    def __init__(self, *args, **kwargs) -> None:
        ...
    def query(self, arg0: str) -> Alignment:
        ...
    def query_batch(self, arg0: list) -> list[Alignment]:
        ...