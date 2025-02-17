// ---------------------------------------------------------------------
//
// Copyright (c) 2020 - 2020 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

#ifndef included_IBTK_config
#define included_IBTK_config

//
// Version information
//

// Major version of IBTK and IBAMR.
#define IBTK_VERSION_MAJOR 0

// Minor version of IBTK and IBAMR.
#define IBTK_VERSION_MINOR 12

// Subminor version of IBTK and IBAMR.
#define IBTK_VERSION_SUBMINOR 0

#define IBTK_VERSION_GTE(major, minor, subminor)                                          \
    ((IBTK_VERSION_MAJOR * 10000 + IBTK_VERSION_MINOR * 100 + IBTK_VERSION_SUBMINOR) >= \
     (major)*10000 + (minor)*100 + (subminor))

//
// Dependency information
//

// Whether or not _Pragma is available
#define IBTK_HAVE_PRAGMA_KEYWORD

// Whether or not libMesh is available
#define IBTK_HAVE_LIBMESH

// Whether or not Silo is available
#define IBTK_HAVE_SILO

//
// Utility macros
//

// Correctly mangle Fortran names
#define IBTK_FC_FUNC(name, NAME) name ## _

// Like IBTK_FC_FUNC, but for names containing underscores
#define IBTK_FC_FUNC_(name, NAME) name ## _

// Macro for disabling warnings
#ifdef IBTK_HAVE_PRAGMA_KEYWORD
// Prevent clang-format from doing strange things to this very long macro:
// clang-format off

// The first four warnings here should be left in that order: new warnings
// should be placed at the end.
#define IBTK_DISABLE_EXTRA_WARNINGS                             \
_Pragma("GCC diagnostic push")                                  \
_Pragma("GCC diagnostic ignored \"-Wunknown-pragmas\"")         \
_Pragma("GCC diagnostic ignored \"-Wpragmas\"")                 \
_Pragma("GCC diagnostic ignored \"-Wunknown-warning-option\"")  \
_Pragma("GCC diagnostic ignored \"-Wunknown-warning\"")         \
_Pragma("GCC diagnostic ignored \"-Wunused-variable\"")         \
_Pragma("GCC diagnostic ignored \"-Wignored-attributes\"")      \
_Pragma("GCC diagnostic ignored \"-Wdeprecated-declarations\"") \
_Pragma("GCC diagnostic ignored \"-Wmisleading-indentation\"")  \
_Pragma("GCC diagnostic ignored \"-Wint-in-bool-context\"")     \
_Pragma("GCC diagnostic ignored \"-Wmaybe-uninitialized\"")     \
_Pragma("GCC diagnostic ignored \"-Wunused-local-typedefs\"")   \
_Pragma("GCC diagnostic ignored \"-Wdeprecated-copy\"")         \
_Pragma("GCC diagnostic ignored \"-Wunused-parameter\"")        \
_Pragma("GCC diagnostic ignored \"-Wparentheses\"")             \
_Pragma("GCC diagnostic ignored \"-Wunneeded-internal-declaration\"")

#define IBTK_ENABLE_EXTRA_WARNINGS _Pragma("GCC diagnostic pop")

// clang-format on
#else

#define IBTK_DISABLE_EXTRA_WARNINGS
#define IBTK_ENABLE_EXTRA_WARNINGS

#endif // #ifdef IBTK_HAVE_PRAGMA_KEYWORD

#endif // included_IBTK_config
