#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: Check these declarations against the C/Fortran source code.  */

/* .C calls */
extern void aalen(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void addmult(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *) ;
extern void atriskindex(void *, void *, void *, void *, void *, void *, void *, void *);
extern void clusterindex(void *, void *, void *, void *, void *, void *, void *, void *);
extern void compSs(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void compSsforward(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void compSsrev(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void confBandBasePredict(void *, void *, void *, void *, void *, void *, void *);
extern void dynadd(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void Gtranssurv(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void itfit(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void localTimeReg(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void mgresid(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void nclusters(void *, void *, void *, void *, void *);
extern void OSbreslow(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *); 
extern void OSsemicox(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void OStimecox(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void * , void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void pava(void *, void *, void *);
extern void pes(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void posubdist2(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void resmean(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void robaalen(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void robaalenC(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void score(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void semiaalen(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void semibreslow(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, 
			void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, 
			void *, void *, void *, void *, void *, void *, void *, void *, void *, void *,
			void *,void *, void *, void *);
extern void semidynadd(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, 
		       void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, 
		       void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, 
		       void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, 
		       void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, 
		       void *, void *, void *, void *, void *);
extern void sindex(void *, void *, void *, void *, void *, void *);
extern void smooth2B(void *, void *, void *, void *, void *, void *, void *, void *);
extern void smoothB(void *, void *, void *, void *, void *, void *, void *, void *);
extern void transsurv(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void twostagereg(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);

static const R_CMethodDef CEntries[] = {
    {"aalen",               (DL_FUNC) &aalen,               11},
    {"addmult",             (DL_FUNC) &addmult,             29},
    {"atriskindex",         (DL_FUNC) &atriskindex,          8},
    {"clusterindex",        (DL_FUNC) &clusterindex,         8},
    {"compSs",              (DL_FUNC) &compSs,              18},
    {"compSsforward",       (DL_FUNC) &compSsforward,       18},
    {"compSsrev",           (DL_FUNC) &compSsrev,           18},
    {"confBandBasePredict", (DL_FUNC) &confBandBasePredict,  7},
    {"dynadd",              (DL_FUNC) &dynadd,              42},
    {"Gtranssurv",          (DL_FUNC) &Gtranssurv,          40},
    {"itfit",               (DL_FUNC) &itfit,               55},
    {"localTimeReg",        (DL_FUNC) &localTimeReg,        10},
    {"mgresid",             (DL_FUNC) &mgresid,             58},
    {"nclusters",           (DL_FUNC) &nclusters,            5},
    {"OSbreslow",           (DL_FUNC) &OSbreslow,           31},
    {"OSsemicox",           (DL_FUNC) &OSsemicox,           37},
    {"OStimecox",           (DL_FUNC) &OStimecox,           32},
    {"pava",                (DL_FUNC) &pava,                 3},
    {"pes",                 (DL_FUNC) &pes,                 29},
    {"posubdist2",          (DL_FUNC) &posubdist2,          48},
    {"resmean",             (DL_FUNC) &resmean,             50},
    {"robaalen",            (DL_FUNC) &robaalen,            37},
    {"robaalenC",           (DL_FUNC) &robaalenC,           32},
    {"score",               (DL_FUNC) &score,               65},
    {"semiaalen",           (DL_FUNC) &semiaalen,           52},
    {"semibreslow",         (DL_FUNC) &semibreslow,         34},
    {"semidynadd",          (DL_FUNC) &semidynadd,          55},
    {"sindex",              (DL_FUNC) &sindex,               6},
    {"smooth2B",            (DL_FUNC) &smooth2B,             8},
    {"smoothB",             (DL_FUNC) &smoothB,              8},
    {"transsurv",           (DL_FUNC) &transsurv,           41},
    {"twostagereg",         (DL_FUNC) &twostagereg,         45},
    {NULL, NULL, 0}
};

void R_init_timereg(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
