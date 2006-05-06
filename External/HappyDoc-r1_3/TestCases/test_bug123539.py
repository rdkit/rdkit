"""member functions with equal names not recognized

If there are several classes in one file and two (or more) classes
contain a function with teh same name, only the last definition /
documentation ist recognized.
"""

class ClassOne:

    def __init__(self):
        "ClassOne __init__"
        pass

    def method_one(self, foo):
        "ClassOne method_one"
        pass

class ClassTwo:

    def __init__(self):
        "ClassTwo __init__"
        pass

    def method_one(self, bar):
        "ClassTwo method_one"
        pass
