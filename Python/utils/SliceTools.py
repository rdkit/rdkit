#
#  Copyright (C) 2000  greg Landrum
#
""" tools for working with python slices

"""

def FindSliceBounds(slice,dim):
    """ makes the bounds of a slice explicit

      **Arguments**

       - slice: a Slice object

       - dim: the dimensions of the thing being sliced

      **Returns**

        a 3-tuple:

          1) the explicit beginning of the slice

          2) the explicit end of the slice

          3) the slice step

    """
    res = [slice.start,slice.stop,slice.step]
    if res[0] is None:
        res[0] = 0
    elif res[0] < 0:
        res[0] = res[0] + dim
    if res[1] is None:
        res[1] = dim
    elif res[1] < 0:
        res[1] = res[1] + dim
    if res[2] is None:
        res[2] = 1
    return tuple(res)
