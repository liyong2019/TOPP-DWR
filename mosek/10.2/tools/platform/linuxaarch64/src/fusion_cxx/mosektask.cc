#include "mosektask_p.h"
#include "mosektask.h"
#include "monty_iterator.h"

#include <iostream>
#define MSK(name,...) if (MSK_RES_OK != MSK_ ## name(task,__VA_ARGS__)) throw MosekException(task)

namespace mosek {
  MSKenv_t   Task::env;
  int        Task::env_refc;
  bool       Task::env_release;

  std::mutex Task::env_lock;

  MSKenv_t Task::getEnv(bool borrow)
  {
    // std::lock_guard<std::mutex> guard(env_lock);
    // if (env == NULL)
    // {
    //   MSKrescodee r = MSK_makeenv(&env, NULL);
    //   env_refc = 0;
    //   env_release = false;
    //   if (r != MSK_RES_OK || env == NULL)
    //   {
    //     std::stringstream ss;
    //     ss << r << ": Failed to create global environment";
    //     throw MosekException(ss.str());
    //   }
    // }
    // if (! borrow)
    //     ++env_refc;
    // return env;
      return NULL; // default to global env
  }
/*
  MSKenv_t Task::getEnv()
  {
    std::lock_guard<std::mutex> guard(env_lock);
    if (env == NULL)
    {
      MSKrescodee r = MSK_makeenv(&env, NULL);
      env_refc = 0;
      env_release = false;
      if (r != MSK_RES_OK || env == NULL)
      {
        std::stringstream ss;
        ss << r << ": Failed to create global environment";
        throw MosekException(ss.str());
      }
    }
    ++env_refc;
    return env;
  }
*/
  void Task::releaseEnv()
  {
    std::lock_guard<std::mutex> guard(env_lock);
    if (env_refc > 0)
    {
        --env_refc;
        if (env_refc == 0 && env_release)
            MSK_deleteenv(&env);
    }
  }

  void Task::releaseGlobalEnv()
  {
      std::lock_guard<std::mutex> guard(env_lock);
      env_release = true;
      if (env != NULL && env_refc == 0)
          MSK_deleteenv(&env);
  }
  void releaseGlobalEnv()
  {
      Task::releaseGlobalEnv();
  }


  void MSKAPI Task::msghandler(void * handle, const char * msg)
  {
    Task * self = (Task*)handle;
    if (self->msguserfunc)
      self->msguserfunc(msg);
  }

  int MSKAPI Task::pgshandler(MSKtask_t         task,
                        void            * handle,
                        MSKcallbackcodee  caller,
                        const double    * dinf,
                        const int       * iinf,
                        const int64_t * linf)
  {
    Task * self = (Task*)handle;

    if (self->datacbuserfunc)
        if (self->datacbuserfunc(caller,dinf,iinf,(int64_t*)linf))
            self->breakflag = 1;

    if (self->cbuserfunc)
        if (self->cbuserfunc(caller))
            self->breakflag = 1;

    return self->breakflag ? 1 : 0;
  }


  static std::string extract_task_error(MSKtask_t task)
  {
    MSKrescodee lastres;
    MSKint32t   len = 0;
    char        symname[MSK_MAX_STR_LEN];
    MSK_getlasterror(task, & lastres, 0, &len, NULL);
    MSK_getcodedesc(lastres, symname, NULL);
    std::vector<char> buf(len+1);
    MSK_getlasterror(task, & lastres, len, NULL, &buf[0]);
    std::string desc = std::string(symname) + "(" + std::to_string(lastres) + "): " + std::string(&buf[0]);
    return desc;
  }

  // mosektask.h
  MosekException::MosekException(MSKtask_t task) : std::runtime_error(extract_task_error(task).c_str())
  {
  }

  Task::Task() :
    breakflag(false),
    task(nullptr),
    msguserfunc(nullptr)
  {
    MSKrescodee r = MSK_maketask(getEnv(true),0,0,&task);
    if (r != MSK_RES_OK)
    {
      std::stringstream ss; ss << r << ": " << "Failed to create native task object";
      throw MosekException(ss.str());
    }
    MSK_linkfunctotaskstream(task,MSK_STREAM_LOG, this, msghandler);
    MSK_putcallbackfunc(task,pgshandler,this);
  }

  Task::Task(Task & other) :
    breakflag(false),
    task(nullptr),
    msguserfunc(nullptr)
  {
      MSKrescodee r = MSK_clonetask(other.task, &task);
      if (r != MSK_RES_OK)
      {
          std::stringstream ss; ss << r << ": " << "Failed to create native task object";
          throw MosekException(ss.str());
      }
      MSK_linkfunctotaskstream(task,MSK_STREAM_LOG, this, msghandler);
      MSK_putcallbackfunc(task,pgshandler,this);
  }

  Task::~Task()
  {
      dispose();
  }

  void Task::dispose()
  {
      if (task)
      {
          MSK_deletetask(&task);
          releaseEnv();
          task = NULL;
      }
  }


  MSKtask_t Task::clone() const
  {
    MSKrescodee r = MSK_RES_OK;
    MSKtask_t t2;
    r = MSK_clonetask(task,&t2);
    if (r != MSK_RES_OK)
    {
      std::stringstream ss; ss << r << ": " << "Failed to create native task object";
      throw MosekException(ss.str());
    }
    return t2;
  }

  MSKtask_t Task::get_task() const
  {
    return task;
  }






  int Task::getnumvar()    const { int num; MSK(getnumvar,&num);    return num; }
  int Task::getnumcon()    const { int num; MSK(getnumcon,&num);    return num; }
  int Task::getnumcone()   const { int num; MSK(getnumcone,&num);   return num; }
  int Task::getnumbarvar() const { int num; MSK(getnumbarvar,&num); return num; }

  int Task::getbarvardim(int j) const { int dim; MSK(getdimbarvarj, j, &dim); return dim; }

  int Task::appendvars  (int num) { int idx = getnumvar();    MSK(appendvars,num); return idx; }
  int Task::appendcons  (int num) { int idx = getnumcon();    MSK(appendcons,num); return idx; }
  int Task::appendbarvar(int dim) { int idx = getnumbarvar(); MSK(appendbarvars, 1, &dim); return idx; }

  int64_t Task::appendprimalpowerconedomain(int conesize,size_t nleft, double *alpha) { int64_t domidx; MSK(appendprimalpowerconedomain,conesize,nleft,alpha,&domidx); return domidx; }
  int64_t Task::appenddualpowerconedomain(int conesize,size_t nleft, double *alpha)   { int64_t domidx; MSK(appenddualpowerconedomain,conesize,nleft,alpha,&domidx); return domidx; }
  int64_t Task::appendquadraticconedomain(int conesize)     { int64_t domidx; MSK(appendquadraticconedomain,    conesize,&domidx); return domidx; }
  int64_t Task::appendrquadraticconedomain(int conesize)    { int64_t domidx; MSK(appendrquadraticconedomain,   conesize,&domidx); return domidx; }
  int64_t Task::appendprimalexpconedomain()                 { int64_t domidx; MSK(appendprimalexpconedomain,             &domidx); return domidx; }
  int64_t Task::appenddualexpconedomain()                   { int64_t domidx; MSK(appenddualexpconedomain,               &domidx); return domidx; }
  // int64_t Task::appendonenormconedomain(int conesize)       { int64_t domidx; MSK(appendonenormconedomain,      conesize,&domidx); return domidx; }
  // int64_t Task::appendinfnormconedomain(int conesize)       { int64_t domidx; MSK(appendinfnormconedomain,      conesize,&domidx); return domidx; }
  int64_t Task::appendprimalgeomeanconedomain(int conesize) { int64_t domidx; MSK(appendprimalgeomeanconedomain,conesize,&domidx); return domidx; }
  int64_t Task::appenddualgeomeanconedomain(int conesize)   { int64_t domidx; MSK(appenddualgeomeanconedomain,  conesize,&domidx); return domidx; }
  int64_t Task::appendrminusdomain(int conesize)            { int64_t domidx; MSK(appendrminusdomain,           conesize,&domidx); return domidx; }
  int64_t Task::appendrplusdomain(int conesize)             { int64_t domidx; MSK(appendrplusdomain,            conesize,&domidx); return domidx; }
  int64_t Task::appendrdomain(int conesize)                 { int64_t domidx; MSK(appendrdomain,                conesize,&domidx); return domidx; }
  int64_t Task::appendrzerodomain(int conesize)             { int64_t domidx; MSK(appendrzerodomain,            conesize,&domidx); return domidx; }
  int64_t Task::appendsvecpsdconedomain(int conesize)       { int64_t domidx; MSK(appendsvecpsdconedomain,      conesize,&domidx); return domidx; }


void Task::varbound(int idx, MSKboundkeye bk, double lb, double ub) { MSK(putvarbound, idx, bk, lb, ub); }

void Task::conbound(int idx, MSKboundkeye bk, double lb, double ub) { MSK(putconbound, idx, bk, lb, ub); }
  void Task::varboundlist
    ( size_t               num,
      const int          * sub,
      const MSKboundkeye   bk,
      const double       * lb,
      const double       * ub)
  {
      std::vector<MSKboundkeye> bks(num, bk);
      MSK(putvarboundlist,num,sub,bks.data(),lb,ub);
  }
  void Task::conboundlist
    ( size_t               num,
      const int          * sub,
      const MSKboundkeye   bk,
      const double       * lb,
      const double       * ub)
  {
      std::vector<MSKboundkeye> bks(num, bk);
      MSK(putconboundlist,num,sub,bks.data(),lb,ub);
  }

  void Task::djcname   (int64_t idx, const std::string & name) { MSK(putdjcname,    idx, name.c_str()); }
  void Task::varname   (int idx, const std::string & name) { MSK(putvarname,    idx, name.c_str()); }
  void Task::conname   (int idx, const std::string & name) { MSK(putconname,    idx, name.c_str()); }
  void Task::conename  (int idx, const std::string & name) { MSK(putconename,   idx, name.c_str()); }
  void Task::barvarname(int idx, const std::string & name) { MSK(putbarvarname, idx, name.c_str()); }
  void Task::objname   (         const std::string & name) { MSK(putobjname,         name.c_str()); }
  void Task::taskname  (         const std::string & name) { MSK(puttaskname,        name.c_str()); }

  void Task::generatedjcnames( int num, int64_t * sub,  const std::string & format, int ndims, int * dims, int64_t * sp, int numnamedaxis, const int * namedaxisidxs, int64_t nameslength, const char ** names) { MSK(generatedjcnames,    num,sub, format.c_str(),ndims,dims,sp,numnamedaxis,namedaxisidxs,nameslength,names); }
  void Task::generateaccnames( int num, int64_t * sub,  const std::string & format, int ndims, int * dims, int64_t * sp, int numnamedaxis, const int * namedaxisidxs, int64_t nameslength, const char ** names) { MSK(generateaccnames,    num,sub, format.c_str(),ndims,dims,sp,numnamedaxis,namedaxisidxs,nameslength,names); }
  void Task::generatevarnames( int num, int     * subj, const std::string & format, int ndims, int * dims, int64_t * sp, int numnamedaxis, const int * namedaxisidxs, int64_t nameslength, const char ** names) { MSK(generatevarnames,    num,subj,format.c_str(),ndims,dims,sp,numnamedaxis,namedaxisidxs,nameslength,names); }
  void Task::generatebarvarnames( int num, int  * subj, const std::string & format, int ndims, int * dims, int64_t * sp, int numnamedaxis, const int * namedaxisidxs, int64_t nameslength, const char ** names) { MSK(generatebarvarnames, num,subj,format.c_str(),ndims,dims,sp,numnamedaxis,namedaxisidxs,nameslength,names); }
  void Task::generateconnames( int num, int     * subi, const std::string & format, int ndims, int * dims, int64_t * sp, int numnamedaxis, const int * namedaxisidxs, int64_t nameslength, const char ** names) { MSK(generateconnames,    num,subi,format.c_str(),ndims,dims,sp,numnamedaxis,namedaxisidxs,nameslength,names); }

  void Task::putarowlist
  ( int num,
    const int *       idxs,
    const int64_t * ptrb,
    const int       * subj,
    const double    * cof)
  {
      MSK(putarowlist64,num,idxs, ptrb,ptrb+1,subj,cof);
  }

  void Task::putacollist
  ( int num,
    const int *       idxs,
    const int64_t * ptrb,
    const int       * subi,
    const double    * cof)
  {
      MSK(putacollist64,num,idxs, ptrb,ptrb+1,subi,cof);
  }


  void Task::putbararowlist(int num,
                            const int       * subi,
                            const int64_t * ptr,
                            const int       * subj,
                            const int64_t * matidx) {
      std::vector<int64_t> nummat(ptr[num]);
      std::vector<double> alpha(ptr[num]);
      for (ptrdiff_t i = 0; i < nummat.size(); ++i) {
          nummat[i] = 1;
          alpha[i] = 1.0;
      }

      MSK(putbararowlist,num,subi,ptr,ptr+1,subj,nummat.data(),matidx,alpha.data());
  }

  void Task::putbarc(int64_t num,
                     const int * subj,
                     const int * subk,
                     const int * subl,
                     const double * val) {
      MSK(putbarcblocktriplet,num,subj,subk,subl,val);
  }

  void Task::putbaraijlist(int64_t num,
                           const int       * subi,
                           const int       * subj,
                           const int64_t * alphaptrb,
                           const int64_t * alphaptre,
                           const int64_t * matidx,
                           const double    * weight)
  {
      if (num > 0)
          MSK(putbaraijlist,num,subi,subj,alphaptrb,alphaptre,matidx,weight);
  }

  void Task::putaijlist
  ( const int    * subi,
    const int    * subj,
    const double * cof,
    int64_t      num)
  {
      if (num > 0) // bug1989
          MSK(putaijlist64,num,subi,subj,cof);
  }

  void Task::putobjsense(MSKobjsensee sense) { MSK(putobjsense, sense); }
  void Task::putclist(size_t num, const int * subj,const double * c)
  {
    MSK(putclist, (int)num, &subj[0], &c[0]);
  }
  void Task::putcfix    (double cfix) { MSK(putcfix, cfix); }

  void Task::appendsymmatlist(int num,
                              const int    * dim,
                              const int64_t * nz,
                              const int    * subk,
                              const int    * subl,
                              const double * val,
                              int64_t    * res )
  {
    MSK(appendsparsesymmatlist,num,dim,nz,subk,subl,val,res);
  }

  void Task::breakOptimize()
  {
      breakflag = true;
  }

  MSKrescodee Task::optimize()
  {
    MSKrescodee r = MSK_RES_OK;
    breakflag = false;
    MSK(optimizetrm, &r);
    return r;
  }

  MSKrescodee Task::optimizermt(const std::string & server, const std::string & port)
  {
    MSKrescodee r = MSK_RES_OK;
    breakflag = false;
    MSK(optimizermt, server.c_str(), port.c_str(), &r);
    return r;
  }

  void Task::putoptserverhost(const std::string & addr)
  {
    MSK(putoptserverhost, addr.c_str());
  }

  void Task::solutionsummary(MSKstreamtypee strm)
  {
    if (msguserfunc)
      MSK(solutionsummary,strm);
  }

  bool Task::solutiondef(MSKsoltypee sol)
  {
    int soldef = 0;
    MSK(solutiondef,sol,&soldef);
    return soldef != 0;
  }

  void Task::getxxslice  (MSKsoltypee sol, int first, int last, double *xx)  { MSK(getxxslice,sol,first,last,xx); }
  void Task::getslxslice (MSKsoltypee sol, int first, int last, double *slx) { MSK(getslxslice,sol,first,last,slx); }
  void Task::getsuxslice (MSKsoltypee sol, int first, int last, double *sux) { MSK(getsuxslice,sol,first,last,sux); }
  void Task::getxcslice  (MSKsoltypee sol, int first, int last, double *xc)  { MSK(getxcslice,sol,first,last,xc);  }
  void Task::getslcslice (MSKsoltypee sol, int first, int last, double *slc) { MSK(getslcslice,sol,first,last,slc); }
  void Task::getsucslice (MSKsoltypee sol, int first, int last, double *suc) { MSK(getsucslice,sol,first,last,suc); }
  void Task::evaluateaccs(MSKsoltypee sol, double *accx ) { MSK(evaluateaccs,sol,accx); }
  void Task::getaccdoty  (MSKsoltypee sol, MSKint64t idx, double *accy) { MSK(getaccdoty,sol,idx,accy); }
  MSKsolstae Task::getsolsta(MSKsoltypee sol) { MSKsolstae sta; MSK(getsolsta,sol,&sta); return sta; }
  MSKprostae Task::getprosta(MSKsoltypee sol) { MSKprostae sta; MSK(getprosta,sol,&sta); return sta; }
  int64_t Task::getaccn  (MSKint64t idx) { int64_t n; MSK(getaccn,idx,&n); return n; }
  int64_t Task::getaccntot() { int64_t n; MSK(getaccntot,&n); return n; }


  void Task::getbarxslice(MSKsoltypee sol, int first, int last, size_t capacity, double* barxj) { MSK(getbarxslice,sol,first,last,capacity,barxj); }
  void Task::getbarsslice(MSKsoltypee sol, int first, int last, size_t capacity, double* barxj) { MSK(getbarsslice,sol,first,last,capacity,barxj); }

  double Task::getprimalobj(MSKsoltypee sol)
  {
    double obj;
    MSK(getprimalobj,sol,&obj);
    return obj;
  }

  double Task::getdualobj(MSKsoltypee sol)
  {
    double obj;
    MSK(getdualobj,sol,&obj);
    return obj;
  }

  void      Task::putparam(const std::string & name, double value)               { MSK(putnadouparam,name.c_str(),value); }
  void      Task::putparam(const std::string & name, int    value)               { MSK(putnaintparam,name.c_str(),value); }
  void      Task::putparam(const std::string & name, const  std::string & value) { MSK(putnastrparam,name.c_str(),value.c_str()); }
  double    Task::getdinfitem (MSKdinfiteme  key) { double     val; MSK(getdouinf,key,&val);  return val; }
  int       Task::getiinfitem (MSKiinfiteme  key) { int        val; MSK(getintinf,key,&val);  return val; }
  int64_t Task::getliinfitem(MSKliinfiteme key) { int64_t  val; MSK(getlintinf,key,&val); return val; }

  void Task::putintlist(size_t num, const int * subj) { for (int j = 0; j < num; ++j) MSK(putvartype, subj[j], MSK_VAR_TYPE_INT); }
  void Task::putcontlist(size_t num, const int * subj) { for (int j = 0; j < num; ++j) MSK(putvartype, subj[j], MSK_VAR_TYPE_CONT); }
  void Task::write(const std::string & filename) { MSK(writedata,filename.c_str()); }
  void Task::analyzeproblem(MSKstreamtypee stream) { MSK(analyzeproblem,stream); }

  void Task::writedatahandle(MSKhwritefunc func, MSKuserhandle_t handle, MSKdataformate format, MSKcompresstypee compress)
  {
      MSK(writedatahandle, func, handle, format, compress);
  }

  void Task::putxxslice(MSKsoltypee which, int first, int last, double * xx)
  {
    MSK(putxxslice,which,first,last,xx);
  }

  void Task::fixvars(int numvar)
  {
    int tasknumvar    = getnumvar();
    if (numvar < tasknumvar)
    {
        std::vector<MSKboundkeye> bk(tasknumvar-numvar);
        std::vector<double> bnd(tasknumvar-numvar);
        MSK(putvarboundslice, numvar, tasknumvar, &bk[0], &bnd[0],&bnd[0]);
    }
  }

  void Task::removevar(int idx)
  {
      MSK(removevars, 1, &idx);
  }

  void putlicensecode(int* code)
  {
      MSKrescodee r = MSK_putlicensecode(Task::getEnv(true), code);
      if (r != MSK_RES_OK) throw MosekException("Failed to put license code");
  }

  void putlicensepath(const std::string & path)
  {
      MSKrescodee r = MSK_putlicensepath(Task::getEnv(true), path.c_str());
      if (r != MSK_RES_OK) throw MosekException("Failed to put license path");
  }

  void putlicensewait(int wait)
  {
      MSK_putlicensewait(Task::getEnv(true), wait);
  }

  void optimizebatch(bool israce, double timelimit, int numthreads, int64_t numtasks, MSKtask_t *tasks, MSKrescodee *trm, MSKrescodee *res) {
      MSKrescodee r = MSK_optimizebatch(Task::getEnv(true), israce, timelimit, numthreads, numtasks, tasks, trm, res);
      if ( r != MSK_RES_OK )
      {
        char symname[MSK_MAX_STR_LEN];
        MSK_getcodedesc(r, symname, NULL);
        std::string desc = "Error in Model::solveBatch: " + std::string(symname) + " (" + std::to_string(r) + ")";
        throw MosekException(desc);
      }
  }

  void Task::revert(int numvar, int numcon, int numcone, int numbarvar)
  {
    int tasknumvar    = getnumvar();
    int tasknumcon    = getnumcon();
    int tasknumcone   = getnumcone();
    int tasknumbarvar = getnumbarvar();

    if (tasknumvar > numvar)
    {
      std::vector<int> subj(tasknumvar-numvar);
      for (int i = numvar; i < tasknumvar; ++i) subj[i-numvar] = i;
      MSK(removevars, tasknumvar-numvar, &subj[0]);
    }
    if (tasknumcone > numcone)
    {
      std::vector<int> subj(tasknumcone-numcone);
      for (int i = numcone; i < tasknumcone; ++i) subj[i-numcone] = i;
      MSK(removecones, tasknumcone-numcone, &subj[0]);
    }

    if (tasknumbarvar > numbarvar)
    {
      std::vector<int> subj(tasknumcon-numcon);
      for (int i = numbarvar; i < tasknumbarvar; ++i) subj[i-numbarvar] = i;
      MSK(removebarvars, tasknumbarvar-numbarvar, &subj[0]);
    }
    if (tasknumcon > numcon)
    {
      std::vector<int> subj(tasknumcon-numcon);
      for (int i = numcon; i < tasknumcon; ++i) subj[i-numcon] = i;
      MSK(removecons, tasknumcon-numcon, &subj[0]);
    }
  }

  int64_t Task::getnumacc() const { int64_t num; MSK(getnumacc,&num); return num; }
  int64_t Task::getnumafe() const { int64_t num; MSK(getnumafe,&num); return num; }
  int64_t Task::getnumdjc() const { int64_t num; MSK(getnumdjc,&num); return num; }
  int64_t Task::appendafes(int64_t num) { int64_t idx = getnumafe();  MSK(appendafes,num); return idx; }

  void Task::putafefrow(int64_t idx, int num, int * subj, double * cof) { MSK(putafefrow, idx, num, subj, cof); }

  void Task::emptyafefrowlist(int64_t num, const int64_t * afeidx) { MSK(emptyafefrowlist,num,afeidx); }
  void Task::emptyafebarfrowlist(int64_t num, const int64_t * afeidx) { MSK(emptyafebarfrowlist,num,afeidx); }
  void Task::putafeglist(int64_t num, const int64_t * idx, const double * g) { MSK(putafeglist,num,idx,g); }

  void Task::putacclist(int64_t numaccs, const int64_t * accidxs, const int64_t * domidxs, int64_t numafeidx, const int64_t * afeidxlist, const double * b) { MSK(putacclist,numaccs,accidxs,domidxs,numafeidx,afeidxlist,b); }
  void Task::appendaccs(int64_t numaccs, const int64_t * domidxs, int64_t numafeidx, const int64_t * afeidxlist, const double * b) { MSK(appendaccs,numaccs,domidxs,numafeidx,afeidxlist,b); }


  void Task::putafebarfentrylist(int64_t num,
                                 const int64_t * afeidxlist,
                                 const int32_t * barvaridxlist,
                                 const int64_t * numtermslist,
                                 const int64_t * ptrtermslist,
                                 int64_t lenterms,
                                 const int64_t * termidx,
                                 const double * termweight) {
      MSK(putafebarfentrylist,num, afeidxlist, barvaridxlist, numtermslist, ptrtermslist, lenterms, termidx, termweight);
  }



  void Task::putafefrowlist(int64_t num, int64_t * afeidxs, int32_t * nzrow, int64_t * ptrrow, int32_t * subj, double * cof) {
    int64_t lenidxval = 0; for (ptrdiff_t i = 0; i < num; ++i) if (lenidxval < ptrrow[i]+nzrow[i]) lenidxval = ptrrow[i]+nzrow[i];
    MSK(putafefrowlist,num,afeidxs,nzrow,ptrrow,lenidxval,subj,cof);
  }

  void Task::putafefentrylist(int64_t num, int64_t * idx, int32_t * subj, double * cof) {
      for (ptrdiff_t i = 0; i < num; ++i)
          MSK(putafefentry,idx[i],subj[i],cof[i]);
  }



  int64_t Task::appenddjcs(int64_t num) { int64_t idx = getnumdjc();  MSK(appenddjcs,num); return idx; }

  void Task::putdjcslice(int64_t first,
                         int64_t last,
                         int64_t numdom,
                         int64_t * domidxlist,
                         int64_t numafe,
                         int64_t * afeidxlist,
                         double *b,
                         int64_t numtermsize,
                         int64_t * termsizes,
                         int64_t * numterm)
  {
    MSK(putdjcslice, first, last, numdom, domidxlist, numafe, afeidxlist, b, numtermsize, termsizes, numterm);
  }

}




void mosek::LinAlg::axpy
( int n,
  double alpha,
  std::shared_ptr<monty::ndarray<double,1>> x,
  std::shared_ptr<monty::ndarray<double,1>> y)
{
    if (n > y->size() ||
        n > x->size())
        throw ArgumentError("Mismatching argument lengths");
    MSKrescodee r = MSK_axpy(Task::getEnv(true),n,alpha,x->raw(),y->raw());
    if (r != MSK_RES_OK)
        throw MosekException("Axpy() failed");
}

void mosek::LinAlg::dot
( int n,
  std::shared_ptr<monty::ndarray<double,1>> x,
  std::shared_ptr<monty::ndarray<double,1>> y,
  double & xty)
{
    if (n != y->size() ||
        n != x->size())
        throw ArgumentError("Mismatching argument lengths");
    MSKrescodee r = MSK_dot(Task::getEnv(true),n,x->raw(),y->raw(),&xty);
    if (r != MSK_RES_OK)
        throw MosekException("dot() failed");
}

void mosek::LinAlg::gemv
( bool     transa,
  int      m,
  int      n,
  double   alpha,
  std::shared_ptr<monty::ndarray<double,1>> a,
  std::shared_ptr<monty::ndarray<double,1>> x,
  double   beta,
  std::shared_ptr<monty::ndarray<double,1>> y)
{
    if ( m*n != a->size() ||
         m > (transa ? x->size() : y->size()) ||
         n > (transa ? y->size() : x->size()) )
        throw ArgumentError("Mismatching argument lengths");
    MSKrescodee r = MSK_gemv(Task::getEnv(true),
                             transa ? MSK_TRANSPOSE_YES : MSK_TRANSPOSE_NO,
                             m,
                             n,
                             alpha,
                             a->raw(),
                             x->raw(),
                             beta,
                             y->raw());
    if (r != MSK_RES_OK)
        throw MosekException("gemv() failed");

}

void mosek::LinAlg::gemm
( bool     transa,
  bool     transb,
  int      m,
  int      n,
  int      k,
  double   alpha,
  std::shared_ptr<monty::ndarray<double,1>> a,
  std::shared_ptr<monty::ndarray<double,1>> b,
  double   beta,
  std::shared_ptr<monty::ndarray<double,1>> c)
{
    if ( m*n != c->size() ||
         n*k != b->size() ||
         m*k != a->size() )
        throw ArgumentError("Mismatching argument lengths");

    MSKrescodee r = MSK_gemm(Task::getEnv(true),
                             transa ? MSK_TRANSPOSE_YES : MSK_TRANSPOSE_NO,
                             transb ? MSK_TRANSPOSE_YES : MSK_TRANSPOSE_NO,
                             m,
                             n,
                             k,
                             alpha,
                             a->raw(),
                             b->raw(),
                             beta,
                             c->raw());
    if (r != MSK_RES_OK)
        throw MosekException("gemm() failed");
}


void mosek::LinAlg::syrk
( MSKuploe uplo,
  bool           transa,
  int            n,
  int            k,
  double         alpha,
  std::shared_ptr<monty::ndarray<double,1>> a,
  double         beta,
  std::shared_ptr<monty::ndarray<double,1>> c)
{
    if ( n*k != a->size() ||
         n*n != c->size() )
        throw ArgumentError("Mismatching argument lengths");

    MSKrescodee r = MSK_syrk(Task::getEnv(true),
                             uplo,
                             transa ? MSK_TRANSPOSE_YES : MSK_TRANSPOSE_NO,
                             n,
                             k,
                             alpha,
                             a->raw(),
                             beta,
                             c->raw());
    if (r != MSK_RES_OK)
        throw MosekException("syrk() failed");
}

void mosek::LinAlg::syeig
( MSKuploe uplo,
  int      n,
  std::shared_ptr<monty::ndarray<double,1>> a,
  std::shared_ptr<monty::ndarray<double,1>> w)
{
    if ( n != a->size()/n ||
         n != w->size() )
        throw ArgumentError("Mismatching argument lengths");

    MSKrescodee r = MSK_syeig(Task::getEnv(true),
                              uplo,
                              n,
                              a->raw(),
                              w->raw());
    if (r != MSK_RES_OK)
        throw MosekException("syeig() failed");
}

void mosek::LinAlg::syevd
( MSKuploe uplo,
  int      n,
  std::shared_ptr<monty::ndarray<double,1>> a,
  std::shared_ptr<monty::ndarray<double,1>> w)
{
    if ( n != a->size()/n ||
         n != w->size() )
        throw ArgumentError("Mismatching argument lengths");

    MSKrescodee r = MSK_syevd(Task::getEnv(true),
                              uplo,
                              n,
                              a->raw(),
                              w->raw());
    if (r != MSK_RES_OK)
        throw MosekException("syevd() failed");
}

void mosek::LinAlg::potrf
( MSKuploe uplo,
  int      n,
  std::shared_ptr<monty::ndarray<double,1>> a)
{
    if ( n != a->size()/n )
        throw ArgumentError("Mismatching argument lengths");

    MSKrescodee r = MSK_potrf(Task::getEnv(true),
                              uplo,
                              n,
                              a->raw());
    if (r != MSK_RES_OK)
        throw MosekException("potrf() failed");
}
