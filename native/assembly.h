#pragma once
#include <stddef.h>
#ifdef _WIN32
# ifdef TM_BUILD
#  define TM_API __declspec(dllexport)
# else
#  define TM_API __declspec(dllimport)
# endif
#else
# define TM_API __attribute__((visibility("default")))
#endif
#ifdef __cplusplus
extern "C" {
#endif
TM_API int tm_abi_version(void);
/* Arrays: nodes, diagonal, load have n+1 entries; all others have n.
   Buffers must be valid and nonoverlapping. Outputs are usable only on success.
   SI inputs; temperatures are element means T-Tref. Signed alpha is allowed.
   Status: 0 success; 1 invalid count/null pointer; 2 invalid input;
   3 nonrepresentable arithmetic (including zero stiffness underflow).
   No allocation, exceptions or global mutable state cross this C ABI. */
TM_API int tm_assemble(size_t n, const double* nodes, const double* young,
                       const double* area, const double* alpha,
                       const double* temperature, double* diagonal,
                       double* off_diagonal, double* load);
#ifdef __cplusplus
}
#endif
