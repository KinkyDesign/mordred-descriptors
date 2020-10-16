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
