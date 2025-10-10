#ifndef _MOSEKTASK_P_H_
#define _MOSEKTASK_P_H_

#include <mosek.h>
#include <mutex>
#include <vector>
#include <memory>
#include <stdint.h>

#include "mosektask.h"

namespace mosek
{
  class Task
  {
  public:
    typedef std::function<void(const std::string &)>                                                      msghandler_t;
  private:
    bool         breakflag;
    MSKtask_t    task;

    msghandler_t      msguserfunc;
    datacbhandler_t   datacbuserfunc;
    cbhandler_t       cbuserfunc;

    static void MSKAPI msghandler(void * handle, const char * msg);
    static int  MSKAPI pgshandler(MSKtask_t         task,
                           void            * handle,
                           MSKcallbackcodee  caller,
                           const double    * dinf,
                           const int       * iinf,
                           const int64_t  * linf);

    //---------------------------------
    static MSKenv_t env;
    static int      env_refc;
    static bool     env_release;
    static std::mutex env_lock;

    //---------------------------------
    //static MSKenv_t getEnv();
    static void releaseEnv();
  public:
    static MSKenv_t getEnv(bool borrow);

    Task();
    Task(Task & other);
    ~Task();
    void dispose();
    MSKtask_t clone() const;
    MSKtask_t get_task() const;

    void setDataCallbackFunc(datacbhandler_t userfunc) { datacbuserfunc = userfunc; }
    void setCallbackFunc(cbhandler_t userfunc)         { cbuserfunc = userfunc; }
    void setStreamFunc  (msghandler_t userfunc)                           { msguserfunc = userfunc; }
    void removeStreamFunc()                     { msguserfunc = nullptr; }

    int getnumvar() const;
    int getnumcon() const;
    int getnumcone() const;
    int getnumbarvar() const;

    int getbarvardim(int idx) const;

    int appendvars(int num);
    int appendcons(int num);
    int appendbarvar(int dim);
    void appendcones(int numcone);
    void putcones(MSKconetypee ct, int * coneidxs, int conesize, int numcone, double * alpha, int * membs);

    // revert the task to this size
    void removevar(int idx);
    void revert(int numvar, int numcon, int numcone, int numbarvar);
    void fixvars(int numvar);

    void varbound     (int idx, MSKboundkeye bk, double lb, double ub);
    void varboundlist (size_t num, const int * sub,     const MSKboundkeye bk, const double * lb, const double * ub);
    void conbound     (int idx, MSKboundkeye bk, double lb, double ub);
    void conboundlist (size_t num, const int * sub,     const MSKboundkeye bk, const double * lb, const double * ub);

    void djcname   (int64_t idx, const std::string & name);
    void varname   (int idx, const std::string & name);
    void conname   (int idx, const std::string & name);
    void conename  (int idx, const std::string & name);
    void barvarname(int idx, const std::string & name);
    void objname   (         const std::string & name);
    void taskname  (         const std::string & name);

    void generatedjcnames(   int num, int64_t * sub,  const std::string & format, int ndims, int * dims, int64_t * sp, int numnamedaxis, const int * namedaxisidx, int64_t nameslength, const char ** names);
    void generateaccnames(   int num, int64_t * sub,  const std::string & format, int ndims, int * dims, int64_t * sp, int numnamedaxis, const int * namedaxisidx, int64_t nameslength, const char ** names);
    void generatevarnames(   int num, int     * subj, const std::string & format, int ndims, int * dims, int64_t * sp, int numnamedaxis, const int * namedaxisidx, int64_t nameslength, const char ** names);
    void generatebarvarnames(int num, int     * subj, const std::string & format, int ndims, int * dims, int64_t * sp, int numnamedaxis, const int * namedaxisidx, int64_t nameslength, const char ** names);
    void generateconnames(   int num, int     * subi, const std::string & format, int ndims, int * dims, int64_t * sp, int numnamedaxis, const int * namedaxisidx, int64_t nameslength, const char ** names);
    void generateconenames(  int num, int     * subk, const std::string & format, int ndims, int * dims, int64_t * sp, int numnamedaxis, const int * namedaxisidx, int64_t nameslength, const char ** names);

    void putarowlist( int num,
                      const int       * idxs,
                      const int64_t * ptrb,
                      const int       * subj,
                      const double    * cof);
    void putacollist( int num,
                      const int       * idxs,
                      const int64_t * ptrb,
                      const int       * subi,
                      const double    * cof);
    // putacolslice
    void putaijlist(const int    * subi,
                    const int    * subj,
                    const double * cof,
                    int64_t      num);



    void putbararowlist(int               num,
                        const int       * subi,
                        const int64_t * ptr,
                        const int       * subj,
                        const int64_t * matidx);

    void putbaraijlist(int64_t num,
                       const int       * subi,
                       const int       * subj,
                       const int64_t * alphaptrb,
                       const int64_t * alphaptre,
                       const int64_t * matidx,
                       const double    * weight);


    void putobjsense(MSKobjsensee sense);
    void putclist   (size_t num, const int * subj, const double * c);
    void putcfix    (double cfix);
    void putbarc (int64_t num,
                  const int *subj,
                  const int *subk,
                  const int *subl,
                  const double * val);

    //--------------------------
    void appendsymmatlist(int num,
                          const int    * dim,
                          const int64_t * nz,
                          const int    * subk,
                          const int    * subl,
                          const double * val,
                          int64_t    * res );

    //--------------------------

    void breakOptimize();
    MSKrescodee optimize();
    MSKrescodee optimizermt(const std::string & server, const std::string & port);
    void putoptserverhost(const std::string & addr);

    void solutionsummary(MSKstreamtypee strm);
    bool solutiondef(MSKsoltypee sol);
    void getbarxslice(MSKsoltypee sol, int first, int last, size_t barvarvecle, double* barxj);
    void getbarsslice(MSKsoltypee sol, int first, int last, size_t barvarvecle, double* barxj);
    void getxxslice  (MSKsoltypee sol, int first, int last, double *xx);
    void getslxslice (MSKsoltypee sol, int first, int last, double *slx);
    void getsuxslice (MSKsoltypee sol, int first, int last, double *sux);
    void getxcslice  (MSKsoltypee sol, int first, int last, double *xc);
    void getslcslice (MSKsoltypee sol, int first, int last, double *slc);
    void getsucslice (MSKsoltypee sol, int first, int last, double *suc);
    void evaluateaccs(MSKsoltypee sol, double *accx);
    MSKsolstae getsolsta(MSKsoltypee sol);
    MSKprostae getprosta(MSKsoltypee sol);
    void getaccdoty  (MSKsoltypee sol, MSKint64t idx, double *accy);
    int64_t getaccn  (MSKint64t idx);
    int64_t getaccntot();

    double getprimalobj(MSKsoltypee sol);
    double getdualobj(MSKsoltypee sol);



    //--------------------------

    int64_t appendquadraticconedomain(int conesize);
    int64_t appendrquadraticconedomain(int conesize);
    int64_t appendprimalexpconedomain();
    int64_t appenddualexpconedomain();
    int64_t appendprimalpowerconedomain(int conesize,size_t nleft,double * alpha);
    int64_t appenddualpowerconedomain(int conesize,size_t nleft,double * alpha);
    /* int64_t appendonenormconedomain(int conesize); */
    /* int64_t appendinfnormconedomain(int conesize); */
    int64_t appendprimalgeomeanconedomain(int conesize);
    int64_t appenddualgeomeanconedomain(int conesize);
    int64_t appendrminusdomain(int conesize);
    int64_t appendrplusdomain(int conesize);
    int64_t appendrdomain(int conesize);
    int64_t appendrzerodomain(int conesize);
    int64_t appendsvecpsdconedomain(int conesize);

    //--------------------------

    void putxxslice(MSKsoltypee which, int first, int last, double * xx);

    //--------------------------

    void putintlist(size_t num, const int * subj);
    void putcontlist(size_t num, const int * subj);

    //--------------------------

    void analyzeproblem(MSKstreamtypee stream);
    void write(const std::string & filename);
    void writedatahandle(MSKhwritefunc func, MSKuserhandle_t handle, MSKdataformate format, MSKcompresstypee compress);

    //--------------------------

    void putparam(const std::string & name, double value);
    void putparam(const std::string & name, int value);
    void putparam(const std::string & name, const std::string & value);

    double    getdinfitem (MSKdinfiteme  key);
    int       getiinfitem (MSKiinfiteme  key);
    int64_t getliinfitem(MSKliinfiteme key);

    static void env_syeig (int n, const double * a, double * w);
    //static void env_potrf (int n, double * a);
    static void env_syevd (int n, double * a, double * w);

    static void releaseGlobalEnv();

    //--------------------------
    int64_t getnumafe() const;
    int64_t getnumdjc() const;
    int64_t getnumacc() const;

    int64_t appendafes(int64_t num);
    void putafefrow(int64_t idx, int num, int * subj, double * cof);
    void putafefrowlist(int64_t num, int64_t * afeidxs, int32_t * nzrow, int64_t * ptrrow, int32_t * subj, double * cof);
    void putafefentrylist(int64_t num, int64_t * idx, int * subj, double * cof);
    void putafebarfentrylist(int64_t num,
                             const int64_t * afeidxlist,
                             const int32_t * barvaridxlist,
                             const int64_t * numtermslist,
                             const int64_t * ptrtermslist,
                             int64_t lenterms,
                             const int64_t * termidx,
                             const double * termweight);


    void putafeglist(int64_t num, const int64_t * idx, const double * g);
    void putacclist(int64_t numaccs, const int64_t * accidxs, const int64_t * domidxs, int64_t numafeidx, const int64_t * afeidxlist, const double * b);
    void appendaccs(int64_t numaccs, const int64_t * domidxs, int64_t numafeidx, const int64_t * afeidxlist, const double * b);
    void emptyafefrowlist(int64_t num, const int64_t * afeidx);
    void emptyafebarfrowlist(int64_t num, const int64_t * afeidx);

    int64_t appenddjcs(int64_t num);
    void putdjcslice(int64_t first, int64_t last,
                     int64_t numdom, int64_t * domidxlist,
                     int64_t numafe, int64_t * afeidxlist, double *b,
                     int64_t numtermsize, int64_t * termsizes,
                     int64_t * numterm);
 };

  void putlicensecode(int* code);
  void putlicensepath(const std::string & path);
  void putlicensewait(int wait);
  void optimizebatch(bool israce, double timelimit, int numthreads, int64_t numtasks, MSKtask_t *tasks, MSKrescodee *trm, MSKrescodee *res);
}
#endif
