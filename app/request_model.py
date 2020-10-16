from pydantic import BaseModel
from typing import Optional, List


# Input model
class EntryId(BaseModel):
    name: str
    ownerUUID: str
    type: Optional[str] = None
    URI: str


class Features(BaseModel):
    name: str
    units: str
    conditions: dict
    category: str
    key: str
    uri: str


class DataEntry(BaseModel):
    entryId: Optional[EntryId] = None
    values: dict


class DataSet(BaseModel):
    dataEntry: Optional[List[DataEntry]] = None
    features: Optional[List[Features]] = None


class Parameters(BaseModel):
    categories: List[str]


class Request(BaseModel):
    dataset: Optional[DataSet] = None
    parameters: Optional[Parameters] = None
