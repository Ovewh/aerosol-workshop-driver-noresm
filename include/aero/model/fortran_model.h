#ifndef AERO_FORTRAN_MODEL_H
#define AERO_FORTRAN_MODEL_H

#include "model.h"

#ifdef __cplusplus
extern "C" {
#endif

/// Creates a new Fortran-backed C Model
/// The flag \c owns_model is used to indicate whether this wrapper owns
/// the underlying model
aero_model_t* aero_new_fortran_model(void *cpp_model, bool owns_model);

#ifdef __cplusplus
} // extern "C"
#endif

#endif
