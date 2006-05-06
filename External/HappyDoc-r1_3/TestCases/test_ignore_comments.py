"""This module is used to test the ignore comments flag.
"""

class WithComments:
    #
    # This class is documented only with comments.
    #
    # Any documentation which appears for this class with the
    # comment flag set to ignore comments indicates a bug.
    #

    def __init__(self):
        #
        # WithComments init method
        #
        
        pass


    

class WithoutComments:
    """This class is documented with __doc__ strings.

    The documentation for this class should always appear.
    """

    def __init__(self):
        "WithoutComments __init__ method."
        pass
    
        
