"""Bug 126855: Expressions as argument defaults."""

def func(arg1=1+2,
         arg2=(3*4),
         arg3='string1' + 'string2',
         arg4=('string4' + 'string5'),
         arg5='str' * 5,
         ):
    """func

    Arguments

      arg1 -- Numerical expression, no parens.

      arg2 -- Numerical expression, with parens.

      arg3 -- String expression, no parens.

      arg4 -- String expression, with parens.

      arg5 -- String multiplication, no parens

    """
    pass
