#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* .C calls */
extern void bilinearform_XAY(void *, void *, void *, void *, void *, void *, void *, void *);
extern void binit(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void cor_diag(void *, void *, void *, void *, void *, void *, void *);
extern void diag_quadraticform_XAX(void *, void *, void *, void *, void *, void *);
extern void diffpairs(void *, void *, void *, void *, void *, void *);
extern void distdiag(void *, void *, void *, void *);
extern void kb_sim_new(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void loccoords(void *, void *, void *, void *, void *, void *, void *);
extern void tgangle(void *, void *, void *, void *);
extern void veccorrval(void *, void *, void *, void *, void *, void *);

static const R_CMethodDef CEntries[] = {
    {"bilinearform_XAY",       (DL_FUNC) &bilinearform_XAY,        8},
    {"binit",                  (DL_FUNC) &binit,                  12},
    {"cor_diag",               (DL_FUNC) &cor_diag,                7},
    {"diag_quadraticform_XAX", (DL_FUNC) &diag_quadraticform_XAX,  6},
    {"diffpairs",              (DL_FUNC) &diffpairs,               6},
    {"distdiag",               (DL_FUNC) &distdiag,                4},
    {"kb_sim_new",             (DL_FUNC) &kb_sim_new,             21},
    {"loccoords",              (DL_FUNC) &loccoords,               7},
    {"tgangle",                (DL_FUNC) &tgangle,                 4},
    {"veccorrval",             (DL_FUNC) &veccorrval,              6},
    {NULL, NULL, 0}
};

void R_init_geoR(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
