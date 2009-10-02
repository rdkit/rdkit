##########################################################################
# Boost Utilities                                                        #
##########################################################################
# Copyright (C) 2007 Douglas Gregor <doug.gregor@gmail.com>              #
# Copyright (C) 2007 Troy Straszheim                                     #
#                                                                        #
# Distributed under the Boost Software License, Version 1.0.             #
# See accompanying file LICENSE_1_0.txt or copy at                       #
#   http://www.boost.org/LICENSE_1_0.txt                                 #
##########################################################################
# Macros in this module:                                                 #
#                                                                        #
#   list_contains: Determine whether a string value is in a list.        #
#                                                                        #
#   car: Return the first element in a list                              #
#                                                                        #
#   cdr: Return all but the first element in a list                      #
#                                                                        #
#   parse_arguments: Parse keyword arguments for use in other macros.    #
##########################################################################

# This utility macro determines whether a particular string value
# occurs within a list of strings:
#
#  list_contains(result string_to_find arg1 arg2 arg3 ... argn)
# 
# This macro sets the variable named by result equal to TRUE if
# string_to_find is found anywhere in the following arguments.
macro(list_contains var value)
  set(${var})
  foreach (value2 ${ARGN})
    if (${value} STREQUAL ${value2})
      set(${var} TRUE)
    endif (${value} STREQUAL ${value2})
  endforeach (value2)
endmacro(list_contains)

# This utility macro extracts the first argument from the list of
# arguments given, and places it into the variable named var.
#
#   car(var arg1 arg2 ...)
macro(car var)
  set(${var} ${ARGV1})
endmacro(car)

# This utility macro extracts all of the arguments given except the
# first, and places them into the variable named var.
#
#   car(var arg1 arg2 ...)
macro(cdr var junk)
  set(${var} ${ARGN})
endmacro(cdr)

# The PARSE_ARGUMENTS macro will take the arguments of another macro and
# define several variables. The first argument to PARSE_ARGUMENTS is a
# prefix to put on all variables it creates. The second argument is a
# list of names, and the third argument is a list of options. Both of
# these lists should be quoted. The rest of PARSE_ARGUMENTS are
# arguments from another macro to be parsed.
# 
#     PARSE_ARGUMENTS(prefix arg_names options arg1 arg2...) 
# 
# For each item in options, PARSE_ARGUMENTS will create a variable with
# that name, prefixed with prefix_. So, for example, if prefix is
# MY_MACRO and options is OPTION1;OPTION2, then PARSE_ARGUMENTS will
# create the variables MY_MACRO_OPTION1 and MY_MACRO_OPTION2. These
# variables will be set to true if the option exists in the command line
# or false otherwise.
# 
# For each item in arg_names, PARSE_ARGUMENTS will create a variable
# with that name, prefixed with prefix_. Each variable will be filled
# with the arguments that occur after the given arg_name is encountered
# up to the next arg_name or the end of the arguments. All options are
# removed from these lists. PARSE_ARGUMENTS also creates a
# prefix_DEFAULT_ARGS variable containing the list of all arguments up
# to the first arg_name encountered.
MACRO(PARSE_ARGUMENTS prefix arg_names option_names)
  SET(DEFAULT_ARGS)
  FOREACH(arg_name ${arg_names})
    SET(${prefix}_${arg_name})
  ENDFOREACH(arg_name)
  FOREACH(option ${option_names})
    SET(${prefix}_${option} FALSE)
  ENDFOREACH(option)

  SET(current_arg_name DEFAULT_ARGS)
  SET(current_arg_list)
  FOREACH(arg ${ARGN})
    LIST_CONTAINS(is_arg_name ${arg} ${arg_names})
    IF (is_arg_name)
      SET(${prefix}_${current_arg_name} ${current_arg_list})
      SET(current_arg_name ${arg})
      SET(current_arg_list)
    ELSE (is_arg_name)
      LIST_CONTAINS(is_option ${arg} ${option_names})
      IF (is_option)
      SET(${prefix}_${arg} TRUE)
      ELSE (is_option)
      SET(current_arg_list ${current_arg_list} ${arg})
      ENDIF (is_option)
    ENDIF (is_arg_name)
  ENDFOREACH(arg)
  SET(${prefix}_${current_arg_name} ${current_arg_list})
ENDMACRO(PARSE_ARGUMENTS)

# Perform a reverse topological sort on the given LIST. 
#   
#   topological_sort(my_list "MY_" "_EDGES")
#
# LIST is the name of a variable containing a list of elements to be
# sorted in reverse topological order. Each element in the list has a
# set of outgoing edges (for example, those other list elements that
# it depends on). In the resulting reverse topological ordering
# (written back into the variable named LIST), an element will come
# later in the list than any of the elements that can be reached by
# following its outgoing edges and the outgoing edges of any vertices
# they target, recursively. Thus, if the edges represent dependencies
# on build targets, for example, the reverse topological ordering is
# the order in which one would build those targets.
#
# For each element E in this list, the edges for E are contained in
# the variable named ${PREFIX}${E}${SUFFIX}, where E is the
# upper-cased version of the element in the list. If no such variable
# exists, then it is assumed that there are no edges. For example, if
# my_list contains a, b, and c, one could provide a dependency graph
# using the following variables:
#
#     MY_A_EDGES     b
#     MY_B_EDGES     
#     MY_C_EDGES     a b
#
#  With the involcation of topological_sort shown above and these
#  variables, the resulting reverse topological ordering will be b, a,
#  c.
function(topological_sort LIST PREFIX SUFFIX)
  # Clear the stack and output variable
  set(VERTICES "${${LIST}}")
  set(STACK)
  set(${LIST})

  # Loop over all of the vertices, starting the topological sort from
  # each one.
  foreach(VERTEX ${VERTICES})
    string(TOUPPER ${VERTEX} UPPER_VERTEX)

    # If we haven't already processed this vertex, start a depth-first
    # search from where.
    if (NOT FOUND_${UPPER_VERTEX})
      # Push this vertex onto the stack with all of its outgoing edges
      string(REPLACE ";" " " NEW_ELEMENT 
        "${VERTEX};${${PREFIX}${UPPER_VERTEX}${SUFFIX}}")
      list(APPEND STACK ${NEW_ELEMENT})

      # We've now seen this vertex
      set(FOUND_${UPPER_VERTEX} TRUE)

      # While the depth-first search stack is not empty
      list(LENGTH STACK STACK_LENGTH)
      while(STACK_LENGTH GREATER 0)
        # Remove the vertex and its remaining out-edges from the top
        # of the stack
        list(GET STACK -1 OUT_EDGES)
        list(REMOVE_AT STACK -1)

        # Get the source vertex and the list of out-edges
        separate_arguments(OUT_EDGES)
        list(GET OUT_EDGES 0 SOURCE)
        list(REMOVE_AT OUT_EDGES 0)

        # While there are still out-edges remaining
        list(LENGTH OUT_EDGES OUT_DEGREE)
        while (OUT_DEGREE GREATER 0)
          # Pull off the first outgoing edge
          list(GET OUT_EDGES 0 TARGET)
          list(REMOVE_AT OUT_EDGES 0)

          string(TOUPPER ${TARGET} UPPER_TARGET)
          if (NOT FOUND_${UPPER_TARGET})
            # We have not seen the target before, so we will traverse
            # its outgoing edges before coming back to our
            # source. This is the key to the depth-first traversal.

            # We've now seen this vertex
            set(FOUND_${UPPER_TARGET} TRUE)

            # Push the remaining edges for the current vertex onto the
            # stack
            string(REPLACE ";" " " NEW_ELEMENT 
              "${SOURCE};${OUT_EDGES}")
            list(APPEND STACK ${NEW_ELEMENT})

            # Setup the new source and outgoing edges
            set(SOURCE ${TARGET})
            string(TOUPPER ${SOURCE} UPPER_SOURCE)
            set(OUT_EDGES 
              ${${PREFIX}${UPPER_SOURCE}${SUFFIX}})
          endif(NOT FOUND_${UPPER_TARGET})

          list(LENGTH OUT_EDGES OUT_DEGREE)
        endwhile (OUT_DEGREE GREATER 0)

        # We have finished all of the outgoing edges for
        # SOURCE; add it to the resulting list.
        list(APPEND ${LIST} ${SOURCE})

        # Check the length of the stack
        list(LENGTH STACK STACK_LENGTH)
      endwhile(STACK_LENGTH GREATER 0)
    endif (NOT FOUND_${UPPER_VERTEX})
  endforeach(VERTEX)

  set(${LIST} ${${LIST}} PARENT_SCOPE)
endfunction(topological_sort)
