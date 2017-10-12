/**
    This source file contains everything that is in freesasa.h
    interface and does not have a natural home in any of the other
    source files.
 */

#if HAVE_CONFIG_H
#include <config.h>
#endif
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include "freesasa_internal.h"

#ifdef PACKAGE_VERSION
const char *freesasa_version = PACKAGE_VERSION;
#else
const char *freesasa_version = "";
#endif

#ifdef PACKAGE_STRING
const char *freesasa_string = PACKAGE_STRING;
#else
const char *freesasa_string = "FreeSASA";
#endif

// Allows compilation with different defaults
// depending on USE_THREADS. but still exposing the value in a header
// that doesn't depend on USE_THREADS
#ifdef USE_THREADS
#define DEF_NUMBER_THREADS 2
#else
#define DEF_NUMBER_THREADS 1
#endif
const int FREESASA_DEF_NUMBER_THREADS = DEF_NUMBER_THREADS;

const freesasa_parameters freesasa_default_parameters = {
    .alg = FREESASA_DEF_ALGORITHM,
    .probe_radius = FREESASA_DEF_PROBE_RADIUS,
    .shrake_rupley_n_points = FREESASA_DEF_SR_N,
    .lee_richards_n_slices = FREESASA_DEF_LR_N,
    .n_threads = DEF_NUMBER_THREADS,
};

static freesasa_result *result_new(int n) {
  freesasa_result *result = malloc(sizeof(freesasa_result));

  if (result == NULL) {
    mem_fail();
    return NULL;
  }

  result->sasa = malloc(sizeof(double) * n);

  if (result->sasa == NULL) {
    mem_fail();
    freesasa_result_free(result);
    return NULL;
  }

  result->n_atoms = n;

  return result;
}

void freesasa_result_free(freesasa_result *r) {
  if (r) {
    free(r->sasa);
    free(r);
  }
}

freesasa_result *freesasa_calc(const coord_t *c, const double *radii,
                               const freesasa_parameters *parameters)

{
  assert(c);
  assert(radii);

  freesasa_result *result = result_new(freesasa_coord_n(c));
  int ret;

  if (result == NULL) {
    fail_msg("");
    return NULL;
  }

  if (parameters == NULL) parameters = &freesasa_default_parameters;

  switch (parameters->alg) {
    case FREESASA_SHRAKE_RUPLEY:
      ret = freesasa_shrake_rupley(result->sasa, c, radii, parameters);
      break;
    case FREESASA_LEE_RICHARDS:
      ret = freesasa_lee_richards(result->sasa, c, radii, parameters);
      break;
    default:
      assert(0);  // should never get here
      break;
  }
  if (ret == FREESASA_FAIL) {
    freesasa_result_free(result);
    return NULL;
  }

  result->total = 0;
  for (int i = 0; i < freesasa_coord_n(c); ++i) {
    result->total += result->sasa[i];
  }
  result->parameters = *parameters;

  return result;
}

freesasa_result *freesasa_calc_coord(const double *xyz, const double *radii,
                                     int n,
                                     const freesasa_parameters *parameters) {
  assert(xyz);
  assert(radii);
  assert(n > 0);

  coord_t *coord = NULL;
  freesasa_result *result = NULL;

  coord = freesasa_coord_new_linked(xyz, n);
  if (coord != NULL) result = freesasa_calc(coord, radii, parameters);
  if (result == NULL) fail_msg("");

  freesasa_coord_free(coord);

  return result;
}

freesasa_result *freesasa_calc_structure(
    const freesasa_structure *structure,
    const freesasa_parameters *parameters) {
  assert(structure);

  return freesasa_calc(freesasa_structure_xyz(structure),
                       freesasa_structure_radius(structure), parameters);
}

freesasa_node *freesasa_calc_tree(const freesasa_structure *structure,
                                  const freesasa_parameters *parameters,
                                  const char *name) {
  assert(structure);

  freesasa_node *tree = NULL;
  freesasa_result *result =
      freesasa_calc(freesasa_structure_xyz(structure),
                    freesasa_structure_radius(structure), parameters);

  if (result != NULL) {
    tree = freesasa_tree_init(result, structure, name);
  } else {
    fail_msg("");
  }

  if (tree == NULL) {
    fail_msg("");
  }

  freesasa_result_free(result);

  return tree;
}

static inline void count_err(int return_value, int *n_err) {
  if (return_value == FREESASA_FAIL) {
    (*n_err)++;
  }
}

int freesasa_tree_export(FILE *file, freesasa_node *root, int options) {
  assert(freesasa_node_type(root) == FREESASA_NODE_ROOT);
  int n_err = 0;
  if (options & FREESASA_LOG) {
    count_err(freesasa_write_log(file, root), &n_err);
  }
  if (options & FREESASA_RES) {
    count_err(freesasa_write_res(file, root), &n_err);
  }
  if (options & FREESASA_SEQ) {
    count_err(freesasa_write_seq(file, root), &n_err);
  }
  if (options & FREESASA_PDB) {
    count_err(freesasa_write_pdb(file, root), &n_err);
  }
  if (options & FREESASA_RSA) {
    count_err(freesasa_write_rsa(file, root, options), &n_err);
  }
  if (options & FREESASA_JSON) {
#if USE_JSON
    count_err(freesasa_write_json(file, root, options), &n_err);
#else
    return fail_msg("library was built without support for JSON output");
#endif
  }
  if (options & FREESASA_XML) {
#if USE_XML
    count_err(freesasa_write_xml(file, root, options), &n_err);
#else
    return fail_msg("library was built without support for XML output");
#endif
  }
  if (n_err > 0) {
    return fail_msg("there were errors when writing output");
  }
  return FREESASA_SUCCESS;
}

freesasa_result *freesasa_result_clone(const freesasa_result *result) {
  freesasa_result *clone = result_new(result->n_atoms);
  if (clone == NULL) {
    fail_msg("");
    return NULL;
  }

  clone->n_atoms = result->n_atoms;
  clone->total = result->total;
  clone->parameters = result->parameters;
  memcpy(clone->sasa, result->sasa, sizeof(double) * clone->n_atoms);

  return clone;
}

const char *freesasa_alg_name(freesasa_algorithm alg) {
  switch (alg) {
    case FREESASA_SHRAKE_RUPLEY:
      return "Shrake & Rupley";
    case FREESASA_LEE_RICHARDS:
      return "Lee & Richards";
  }
  assert(0 && "Illegal algorithm");
}
