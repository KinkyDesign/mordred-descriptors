from mordred import Calculator, descriptors
from rdkit import Chem


# Calculation of mordred descriptors for a given smile
def calculate_single(smile: str):
    calc = Calculator(descriptors, ignore_3D=True)
    mol = Chem.MolFromSmiles(smile)
    c = calc(mol).fill_missing(value="0")
    return c


# Converts calculations string to a python dictionary
def format_calculations(calculations: str):
    calculations = calculations.replace("'", "")
    dict1 = dict((x.strip(), y.strip())
                 for x, y in (element.split(':')
                              for element in calculations.split(', ')))
    return dict1


# Get all descriptor names
def get_descriptors():
    calc = Calculator(descriptors, ignore_3D=True)
    mol = Chem.MolFromSmiles('c1ccccc1')
    variable = calc(mol).asdict(False)
    descriptor_names_list = []
    for key, value in variable.items():
        descriptor_names_list.append(key)
    return descriptor_names_list
