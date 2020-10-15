from rdkit import Chem
from mordred import Calculator, descriptors
from fastapi import FastAPI
from pydantic import BaseModel
from fastapi.responses import JSONResponse
from typing import Optional, List
import uvicorn


# Calculation of mordred descriptors for a given smile
def calculate_single(smile: str):
    calc = Calculator(descriptors, ignore_3D=True)
    mol = Chem.MolFromSmiles(smile)
    return calc(mol)


# Converts calculations string to a python dictionary
def format_calculations(calculations: str):
    calculations = calculations.replace("'", "")
    dict1 = dict((x.strip(), y.strip())
                 for x, y in (element.split(':')
                              for element in calculations.split(', ')))
    return dict1


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


# Data Entry for output model
class DataEntryOut:

    def __init__(self, entryId, values):
        self.entryId = entryId
        self.values = values


class FeaturesOut:
    def __init__(self, name, conditions, category, uri):
        self.name = name
        self.conditions = conditions
        self.category = category
        self.uri = uri


app = FastAPI(
    debug=True,
    title="Mordred Descriptors API"
)


@app.post("/mordred")
async def apply_descriptor(request: Request):
    """
    Applies Mordred Descriptor on Dataset and creates a new Dataset
    """
    # Error check between the key at features and values
    for i in range(1, len(request.dataset.dataEntry) + 1):
        for key in request.dataset.dataEntry[i - 1].values:
            if key != request.dataset.features[0].key:
                return JSONResponse(status_code=422, content={422: "Mismatching Keys"})

    # Get all smiles from request
    smiles = []
    features_key = request.dataset.features[0].key
    for entry in request.dataset.dataEntry:
        smiles.append(entry.values[features_key])

    # Calculate all given smiles
    calculations = []
    for s in smiles:
        calculations.append(str(calculate_single(s)))

    # Strip the mordred generated string from unwanted brackets
    calculations = [item.replace('Result({', '') for item in calculations]
    calculations = [item.replace('})', '') for item in calculations]

    data_entry_out = [DataEntryOut(entryId={"name": s, "ownerUUID": None, "type": None, "URI": None},
                                   values=format_calculations(calculations[smiles.index(s)])) for s in smiles]
    # Get all descriptor names
    calc = Calculator(descriptors, ignore_3D=True)
    mol = Chem.MolFromSmiles('c1ccccc1')
    variable = calc(mol).asdict(False)
    descriptor_names_list = []
    for key, value in variable.items():
        descriptor_names_list.append(key)

    conditions_mordred = {
        "Implementation Vendor": "Molecular Descriptor Calculator",
        "Implementation Identifier": "1.2.1a1",
        "Implementation Title": "https://jcheminf.biomedcentral.com/articles/10.1186/s13321-018-0258-y",
        "Specification Reference": "http://mordred-descriptor.github.io/documentation/master/"
    }
    features_out = [
        FeaturesOut(name=d + " descriptor of feature with name smiles and URI " + request.dataset.features[0].uri,
                    conditions=conditions_mordred, category="CDK", uri=d) for d in descriptor_names_list]

    response_object = {
        "responseDataset": {
            "dataEntry": data_entry_out,
            # "features": [{}],
            "features": features_out,
            "descriptors": ["CDK"]
        }
    }

    return response_object


if __name__ == "__main__":
    uvicorn.run("main:app", host="0.0.0.0", port=8000, reload=True)
