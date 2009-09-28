#! /bin/sh
mkdir -p bindings/{libs,boost}/numeric
(cd bindings/boost/numeric;
  svn co http://svn.boost.org/svn/boost/sandbox/boost/numeric/bindings)
(cd bindings/libs/numeric;
  svn co http://svn.boost.org/svn/boost/sandbox/libs/numeric/bindings)

