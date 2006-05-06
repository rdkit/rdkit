"""Test case for function argument handling with lots of different types of arguments."""

def example_function_with_args(arg1, arg2,
                               arg3withDefault='hi there',
                               arg3aWithDefault="'hi again'",
                               arg3bWithDefault='"hi there again"',
                               arg4DefaultInt=101,
                               arg5DefaultTuple=(1,2),
                               arg6DefaultList=[3,4],
                               arg7DefaultNone=None,
                               arg8DefaultName=foo,
                               arg9DefaultInstance=DefaultClassInst(),
                               arg10DefaultInstanceWithParams= \
                               DefaultClassInstWithParams(1, 2,
                                                          ('tuple', 'param'),
                                                          ['list', 'param']
                                                          ),
                               negativeIntArg=-1,
                               floatArg=1.2,
                               negativeFloatArg=-3.4,
                               mathArg=1 + 2,
                               stringArgWithHTML='<h1>Hi, Dick & Jane!</h1>',
                               ):
    "This is an example function for testing purposes."
    if one:
        raise IOError('RAISE_class')
    else:
        raise 'RAISE_blah2'
    for i in range(1, 10):
        raise 'RAISE_loop'
    raise 'RAISE_main_level'
    return None
