#if defined(__clang__)
/* Clang/LLVM. ---------------------------------------------- */
#pragma GCC diagnostic pop

#elif defined(__ICC) || defined(__INTEL_COMPILER)
/* Intel ICC/ICPC. ------------------------------------------ */

#elif (defined(__GNUC__) || defined(__GNUG__)) && \
    (__GNUC__ > 4 || (__GNUC__ == 4 && __GNUC_MINOR__ > 5))
/* GNU GCC/G++. these pragmas only work with >v4.1
 * --------------------------------------------- */
#pragma GCC diagnostic pop

#elif defined(__HP_cc) || defined(__HP_aCC)
/* Hewlett-Packard C/aC++. ---------------------------------- */

#elif defined(__IBMC__) || defined(__IBMCPP__)
/* IBM XL C/C++. -------------------------------------------- */

#elif defined(_MSC_VER)
/* Microsoft Visual Studio. --------------------------------- */
#pragma warning(pop)
#elif defined(__PGI)
/* Portland Group PGCC/PGCPP. ------------------------------- */

#elif defined(__SUNPRO_C) || defined(__SUNPRO_CC)
/* Oracle Solaris Studio. ----------------------------------- */

#endif
