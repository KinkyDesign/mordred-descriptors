from fastapi import FastAPI
from fastapi.responses import JSONResponse
import uvicorn
from calculations import calculate_single, format_calculations, get_descriptors
from request_model import *
from response_model import DataEntryOut, FeaturesOut


app = FastAPI(
    # debug=True,
    title="Mordred Descriptors API"
)


@app.post("/mordred/{smile}")
async def apply_descriptor_to_smile(smile: str):
    """
    Applies Mordred Descriptor on a given smile
    """
    # Convert Result to dict and return
    return calculate_single(smile).asdict(rawkey=False)


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

    descriptor_names_list = get_descriptors()

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
            "features": features_out,
            "descriptors": ["CDK"]
        }
    }

    return response_object


if __name__ == "__main__":
    uvicorn.run("main:app", host="0.0.0.0", port=8000, reload=True)
