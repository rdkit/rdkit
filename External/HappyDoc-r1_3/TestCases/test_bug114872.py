"""Test case for bug in exception detection."""

def parameterless_raise_in_function():
    "This function raises an exception that it catches."
    try:
        pass
    except:
        raise

class ParameterlessRaiseClass:
    "This class raises an exception that it catches."

    def __init__(self):
        "This method raises an exception that it catches."
        try:
            pass
        except:
            raise
        return

