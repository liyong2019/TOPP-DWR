#include "mosek.h"
#include "mosektask.h"
#include "mosektask_p.h"
#include "fusion_p.h"
#include <string>
#include <vector>
#include <map>
#include <memory>
#include <cassert>
#include <SolverInfo.h>

#include "mosektask_p.h"


namespace mosek
{
  namespace fusion
  {
    //-----------------------------
    BaseModel::BaseModel(p_BaseModel * _impl) : _impl(_impl)
    {
    }
    BaseModel::~BaseModel()
    {
        //std::cout << "BaseModel::~BaseModel()" << std::endl;
      delete _impl;
      _impl = nullptr;
    }

    void BaseModel::dispose()
    {
        //std::cout << "dispose(). ref count : " <<  this->refcount << std::endl;
        _impl->dispose();
    }

    //-----------------------------
    p_BaseModel::p_BaseModel(BaseModel * _pubthis) :
      synched(false),
      sol_itr(nullptr),
      sol_itg(nullptr),
      sol_bas(nullptr)
    {}

    void p_BaseModel::_initialize( const std::string & name,
                                   const std::string & licpath)
    {
        task = std::unique_ptr<Task>(new Task());
        taskname = name;
        task->taskname(name);
    }

    void p_BaseModel::_initialize( monty::rc_ptr<BaseModel> m)
    {
        task = std::unique_ptr<Task>(new Task(* (m->_impl->task)));
        task->putparam("MSK_IPAR_LOG_EXPAND",0);
        task->putparam("MSK_IPAR_REMOVE_UNUSED_SOLUTIONS",1);
        synched  = m->_impl->synched;
        taskname = m->_impl->taskname;
        sol_itr  = m->_impl->sol_itr.get() ? m->_impl->sol_itr->clone() : nullptr;
        sol_bas  = m->_impl->sol_bas.get() ? m->_impl->sol_bas->clone() : nullptr;
        sol_itg  = m->_impl->sol_itg.get() ? m->_impl->sol_itg->clone() : nullptr;
    }


    void p_BaseModel::report_solution_get_xx(array_t<double> v)  { task->getxxslice(  cursol,0,v->size(),v->raw()); }
    void p_BaseModel::report_solution_get_slx(array_t<double> v) { task->getslxslice( cursol,0,v->size(),v->raw()); }
    void p_BaseModel::report_solution_get_sux(array_t<double> v) { task->getsuxslice( cursol,0,v->size(),v->raw()); }
    void p_BaseModel::report_solution_get_xc(array_t<double> v)  { task->getxcslice(  cursol,0,v->size(),v->raw()); }
    void p_BaseModel::report_solution_get_slc(array_t<double> v) { task->getslcslice( cursol,0,v->size(),v->raw()); }
    void p_BaseModel::report_solution_get_suc(array_t<double> v) { task->getsucslice( cursol,0,v->size(),v->raw()); }
    void p_BaseModel::report_solution_get_barx(array_t<double> v){ task->getbarxslice(cursol,0,task->getnumbarvar(),v->size(),v->raw()); }
    void p_BaseModel::report_solution_get_bars(array_t<double> v){ task->getbarsslice(cursol,0,task->getnumbarvar(),v->size(),v->raw()); }
    void p_BaseModel::report_solution_get_accx(array_t<double> v){ task->evaluateaccs(cursol,v->raw()); }
    void p_BaseModel::report_solution_get_accy(array_t<double> v) {
      size_t ptr = 0;
      int64_t numacc = task->getnumacc();
      for (int64_t i = 0; i < numacc; ++i) {
          task->getaccdoty(cursol,i,v->raw()+ptr);
          ptr += task->getaccn(i);
      }
    }
    void p_BaseModel::report_solution_get_accptr(array_t<int32_t> v) {
      int64_t ptr = 0;
      int64_t numacc = task->getnumacc();
      for (int i = 0; i < numacc; ++i) {
          (*v)[i] = (int)ptr;
          ptr += task->getaccn(i);
      }
    }




    void p_BaseModel::task_setLogHandler (const msghandler_t & userfunc) { task->setStreamFunc(userfunc); }
    void p_BaseModel::task_setDataCallbackHandler (const datacbhandler_t & userfunc) { task->setDataCallbackFunc(userfunc); }
    void p_BaseModel::task_setCallbackHandler (const cbhandler_t & userfunc) { task->setCallbackFunc(userfunc); }


    int  p_BaseModel::task_append_barvar(int dim, int num)
    {
      int idx = task->getnumbarvar();
      for (int i = 0; i < num; ++i)
        task->appendbarvar(dim);
      return idx;
    }

    int  p_BaseModel::task_append_con(int num) { return task->appendcons(num); }

    int  p_BaseModel::task_append_var(int num)  { return task->appendvars(num); }
    int  p_BaseModel::task_barvardim(int index) { return task->getbarvardim(index);  }

    int  p_BaseModel::task_numbarvar()          { return task->getnumbarvar(); }
    int  p_BaseModel::task_numcon()             { return task->getnumcon(); }
    int  p_BaseModel::task_numvar()             { return task->getnumvar(); }
    int  p_BaseModel::task_numacc()             { return (int)task->getnumacc(); }
    int  p_BaseModel::task_numafe()             { return (int)task->getnumafe(); }
    int  p_BaseModel::task_numdjc()             { return (int)task->getnumdjc(); }

    void p_BaseModel::task_barvar_name(int idx, const std::string & name) { task->barvarname(idx,name); }
    void p_BaseModel::task_con_name   (int idx, const std::string & name) { task->conname   (idx,name); }
    void p_BaseModel::task_var_name   (int idx, const std::string & name) { task->varname   (idx,name); }
    void p_BaseModel::task_djc_name   (int64_t idx, const std::string & name) { task->djcname   (idx,name); }


    static void convert_names(const std::shared_ptr<monty::ndarray<std::shared_ptr<monty::ndarray<std::string,1>>>> names,
                         std::vector<int> & v_namedaxisidxs,
                         std::vector<const char *> & v_namelist) {
        v_namedaxisidxs.clear();
        v_namelist.clear();

        if (names.get() != NULL) {
            int axisi = 0;
            for (auto &nl: *names) {
                if (nl.get() != NULL) {
                    v_namedaxisidxs.push_back(axisi);
                    for (auto &s: *nl)
                        v_namelist.push_back(s.c_str());
                }
                ++axisi;
            }
        }
    }

    void p_BaseModel::task_format_djc_names(const std::shared_ptr<monty::ndarray<int64_t,1>> sub,
                                            const std::string & format,
                                            const std::shared_ptr<monty::ndarray<int,1>> dims,
                                            const std::shared_ptr<monty::ndarray<std::shared_ptr<monty::ndarray<std::string,1>>>> names)
    {
        std::vector<int64_t> v_sub(sub->flat_begin(),sub->flat_end());
        std::vector<int>     v_dims(dims->flat_begin(),dims->flat_end());
        std::vector<int>     v_namedaxisidxs;
        std::vector<const char *> v_namelist;
        convert_names(names,v_namedaxisidxs,v_namelist);
        if (v_namelist.size() > 0)
            task->generatedjcnames((int64_t) v_sub.size(), v_sub.data(), format, (int)v_dims.size(), v_dims.data(),NULL,(int32_t)v_namedaxisidxs.size(),v_namedaxisidxs.data(),(int64_t)v_namelist.size(),v_namelist.data());
        else
            task->generatedjcnames((int64_t) v_sub.size(), v_sub.data(), format, (int)v_dims.size(), v_dims.data(),NULL,0,NULL,0,NULL);
    }

    void p_BaseModel::task_format_acc_names(const std::shared_ptr<monty::ndarray<int64_t,1>> sub,
                                            const std::string & format,
                                            const std::shared_ptr<monty::ndarray<int,1>> dims,
                                            const std::shared_ptr<monty::ndarray<std::shared_ptr<monty::ndarray<std::string,1>>>> names )
    {
        std::vector<int64_t> v_sub(sub->flat_begin(),sub->flat_end());
        std::vector<int> v_dims(dims->flat_begin(),dims->flat_end());
        std::vector<int>     v_namedaxisidxs;
        std::vector<const char *> v_namelist;
        convert_names(names,v_namedaxisidxs,v_namelist);
        if (v_namelist.size() > 0)
            task->generateaccnames((int64_t) v_sub.size(), v_sub.data(), format, (int)v_dims.size(), v_dims.data(),NULL,(int32_t)v_namedaxisidxs.size(),v_namedaxisidxs.data(),(int64_t)v_namelist.size(),v_namelist.data());
        else
            task->generateaccnames((int64_t) v_sub.size(), v_sub.data(), format, (int)v_dims.size(), v_dims.data(),NULL,0,NULL,0,NULL);
    }

    void p_BaseModel::task_format_var_names(const std::shared_ptr<monty::ndarray<int,1>> subj,
                                            const std::string & format,
                                            const std::shared_ptr<monty::ndarray<int,1>> dims,
                                            const std::shared_ptr<monty::ndarray<int64_t,1>> sp,
                                            const std::shared_ptr<monty::ndarray<std::shared_ptr<monty::ndarray<std::string,1>>>> names )
    {
        std::vector<int> v_subj(subj->flat_begin(),subj->flat_end());
        std::vector<int> v_dims(dims->flat_begin(),dims->flat_end());
        std::vector<int64_t> v_sp;
        int64_t * p_sp = NULL;
        std::vector<int>     v_namedaxisidxs;
        std::vector<const char *> v_namelist;
        convert_names(names,v_namedaxisidxs,v_namelist);
        if (sp.get() != NULL) {
            v_sp = std::vector<int64_t> (sp->flat_begin(),sp->flat_end());
            p_sp = v_sp.data();
        }
        if (v_namelist.size() > 0)
            task->generatevarnames((int) v_subj.size(), v_subj.data(), format, (int)v_dims.size(), v_dims.data(),p_sp,(int32_t)v_namedaxisidxs.size(),v_namedaxisidxs.data(),(int64_t)v_namelist.size(),v_namelist.data());
        else
            task->generatevarnames((int) v_subj.size(), v_subj.data(), format, (int)v_dims.size(), v_dims.data(),p_sp,0,NULL,0,NULL);
    }

    void p_BaseModel::task_format_barvar_names(const std::shared_ptr<monty::ndarray<int,1>>     subj,
                                               const std::string                              & format,
                                               const std::shared_ptr<monty::ndarray<int,1>>     dims,
                                               const std::shared_ptr<monty::ndarray<std::shared_ptr<monty::ndarray<std::string,1>>>> names )
    {
        std::vector<int> v_subj(subj->flat_begin(),subj->flat_end());
        std::vector<int> v_dims(dims->flat_begin(),dims->flat_end());
        std::vector<int>     v_namedaxisidxs;
        std::vector<const char *> v_namelist;
        convert_names(names,v_namedaxisidxs,v_namelist);
        if (v_namelist.size() > 0)
            task->generatebarvarnames((int) v_subj.size(), v_subj.data(), format, (int)v_dims.size(), v_dims.data(),NULL,(int32_t)v_namedaxisidxs.size(),v_namedaxisidxs.data(),(int64_t)v_namelist.size(),v_namelist.data());
        else
            task->generatebarvarnames((int) v_subj.size(), v_subj.data(), format, (int)v_dims.size(), v_dims.data(),NULL,0,NULL,0,NULL);
    }

    void p_BaseModel::task_format_con_names(const std::shared_ptr<monty::ndarray<int,1>> subi,
                                            const std::string & format,
                                            const std::shared_ptr<monty::ndarray<int,1>> dims,
                                            const std::shared_ptr<monty::ndarray<int64_t,1>> sp,
                                            const std::shared_ptr<monty::ndarray<std::shared_ptr<monty::ndarray<std::string,1>>>> names )
    {
        std::vector<int> v_subi(subi->flat_begin(),subi->flat_end());
        std::vector<int> v_dims(dims->flat_begin(),dims->flat_end());
        std::vector<int64_t> v_sp;
        int64_t * p_sp = NULL;
        std::vector<int>     v_namedaxisidxs;
        std::vector<const char *> v_namelist;
        convert_names(names,v_namedaxisidxs,v_namelist);
        if (sp.get() != NULL) {
            v_sp = std::vector<int64_t> (sp->flat_begin(),sp->flat_end());
            p_sp = v_sp.data();
        }
        if (v_namelist.size() > 0)
            task->generateconnames((int) v_subi.size(), v_subi.data(), format, (int) v_dims.size(), v_dims.data(), p_sp,(int32_t)v_namedaxisidxs.size(),v_namedaxisidxs.data(),(int64_t)v_namelist.size(),v_namelist.data());
        else
            task->generateconnames((int) v_subi.size(), v_subi.data(), format, (int) v_dims.size(), v_dims.data(), p_sp,0,NULL,0,NULL);
    }

    void p_BaseModel::task_put_param(const std::string & name, double value) { task->putparam(name,value); }
    void p_BaseModel::task_put_param(const std::string & name, int    value) { task->putparam(name,value); }
    void p_BaseModel::task_put_param(const std::string & name, const std::string & value) { task->putparam(name,value); }

    double p_BaseModel::task_get_dinf (const std::string & name)
    {
      MSKdinfiteme key;
      if (SolverInfo::getdouinf(name,key))
        return task->getdinfitem(key);
      else
        throw NameError("Invalid double information item name");
    }

    int       p_BaseModel::task_get_iinf (const std::string & name)
    {
      MSKiinfiteme key;
      if (SolverInfo::getintinf(name,key))
        return task->getiinfitem(key);
      else
        throw NameError("Invalid integer information item name");
    }

    int64_t p_BaseModel::task_get_liinf(const std::string & name)
    {
      MSKliinfiteme key;
      if (SolverInfo::getlintinf(name,key))
        return task->getliinfitem(key);
      else
        throw NameError("Invalid long integer information item name");
    }


    void p_BaseModel::task_putaijlist
      ( const std::shared_ptr<monty::ndarray<int,1>>    & subi,
        const std::shared_ptr<monty::ndarray<int,1>>    & subj,
        const std::shared_ptr<monty::ndarray<double,1>> & cof,
        int64_t                        num )
    { task->putaijlist(subi->raw(),subj->raw(),cof->raw(),num); }

    void p_BaseModel::task_putarowlist
      ( const std::shared_ptr<monty::ndarray<int,1>>       & idxs,
        const std::shared_ptr<monty::ndarray<int64_t,1>> & ptrb,
        const std::shared_ptr<monty::ndarray<int,1>>       & subj,
        const std::shared_ptr<monty::ndarray<double,1>>    & cof)
    {
        assert(ptrb->size() == idxs->size()+1);
        for (auto i = ptrb->flat_begin(), e = ptrb->flat_end(); i != e; ++i)
        {
            assert(*i <= subj->size());
            assert(*i <= cof->size());
        }
        task->putarowlist((int)idxs->size(), idxs->raw(), ptrb->raw(),subj->raw(),cof->raw());
    }

    void p_BaseModel::task_cleararowlist ( const std::shared_ptr<monty::ndarray<int,1>> & idxs) {
        std::vector<int64_t> ptrb(idxs->size()+1); // zeros
        double cof;
        int subj;
        task->putarowlist((int)idxs->size(), idxs->raw(), ptrb.data(),&subj,&cof);
    }

    void p_BaseModel::task_clearacollist ( const std::shared_ptr<monty::ndarray<int,1>> & idxs) {
        std::vector<int64_t> ptrb(idxs->size()+1); // zeros
        double cof;
        int subj;
        task->putacollist((int)idxs->size(), idxs->raw(), ptrb.data(),&subj,&cof);
    }







    void p_BaseModel::task_putbararowlist
    ( const std::shared_ptr<monty::ndarray<int,1>>       subi,
      const std::shared_ptr<monty::ndarray<int64_t,1>> ptr,
      const std::shared_ptr<monty::ndarray<int,1>>       subj,
      const std::shared_ptr<monty::ndarray<int64_t,1>> matidx) {
        task->putbararowlist(subi->size(),
                             subi->raw(),
                             ptr->raw(),
                             subj->raw(),
                             matidx->raw());
    }

    void p_BaseModel::task_putbaraijlist
    ( const std::shared_ptr<monty::ndarray<int,1>>       subi,
      const std::shared_ptr<monty::ndarray<int,1>>       subj,
      const std::shared_ptr<monty::ndarray<int64_t,1>> matidx) {
        std::vector<int64_t> ptr(subi->size()+1);
        std::vector<double> weight(subi->size());
        for (ptrdiff_t i = 0, e = subi->size()+1; i < e; ++i) ptr[i] = i;
        for (ptrdiff_t i = 0, e = subi->size(); i < e; ++i) weight[i] = 1.0;
        task->putbaraijlist(subi->size(),
                            subi->raw(),
                            subj->raw(),
                            ptr.data(),
                            ptr.data()+1,
                            matidx->raw(),
                            weight.data());
    }

    void p_BaseModel::task_putbarc
    ( const std::shared_ptr<monty::ndarray<int,1>>    subj,
      const std::shared_ptr<monty::ndarray<int,1>>    subk,
      const std::shared_ptr<monty::ndarray<int,1>>    subl,
      const std::shared_ptr<monty::ndarray<double,1>> val) {
        task->putbarc(subj->size(),
                      subj->raw(),
                      subk->raw(),
                      subl->raw(),
                      val->raw());
    }

    std::shared_ptr<monty::ndarray<int64_t,1>> p_BaseModel::task_appendsymmatlist
      ( const std::shared_ptr<monty::ndarray<int,1>> & dim,
        const std::shared_ptr<monty::ndarray<int64_t,1>> & nz,
        const std::shared_ptr<monty::ndarray<int,1>> & subk,
        const std::shared_ptr<monty::ndarray<int,1>> & subl,
        const std::shared_ptr<monty::ndarray<double,1>> & val) {
        std::shared_ptr<monty::ndarray<int64_t,1>> r = std::shared_ptr<monty::ndarray<int64_t,1>>(new monty::ndarray<int64_t,1>(dim->size()));
        task->appendsymmatlist(dim->size(),
                               dim->raw(),
                               nz->raw(),
                               subk->raw(),
                               subl->raw(),
                               val->raw(),
                               r->raw());
        return r;
    }

    void p_BaseModel::task_con_putboundlist_fr(const std::shared_ptr<monty::ndarray<int,1>> idxs) { task->conboundlist(idxs->size(), idxs->raw(), MSK_BK_FR, NULL,NULL); }
    void p_BaseModel::task_con_putboundlist_lo(const std::shared_ptr<monty::ndarray<int,1>> idxs, const std::shared_ptr<monty::ndarray<double,1>> & rhs) { task->conboundlist(idxs->size(), idxs->raw(), MSK_BK_LO, rhs->raw(),rhs->raw()); }
    void p_BaseModel::task_con_putboundlist_up(const std::shared_ptr<monty::ndarray<int,1>> idxs, const std::shared_ptr<monty::ndarray<double,1>> & rhs) { task->conboundlist(idxs->size(), idxs->raw(), MSK_BK_UP, rhs->raw(),rhs->raw()); }
    void p_BaseModel::task_con_putboundlist_fx(const std::shared_ptr<monty::ndarray<int,1>> idxs, const std::shared_ptr<monty::ndarray<double,1>> & rhs) { task->conboundlist(idxs->size(), idxs->raw(), MSK_BK_FX, rhs->raw(),rhs->raw()); }
    void p_BaseModel::task_con_putboundlist_ra(const std::shared_ptr<monty::ndarray<int,1>> idxs,
                                               const std::shared_ptr<monty::ndarray<double,1>> & lb ,
                                               const std::shared_ptr<monty::ndarray<double,1>> & ub ) {
        task->conboundlist(idxs->size(), idxs->raw(), MSK_BK_RA,lb->raw(),ub->raw());
    }

    void p_BaseModel::task_var_putboundlist_fr(const std::shared_ptr<monty::ndarray<int,1>> idxs) {task->varboundlist(idxs->size(), idxs->raw(), MSK_BK_FR, NULL,NULL); }
    void p_BaseModel::task_var_putboundlist_lo(const std::shared_ptr<monty::ndarray<int,1>> idxs, const std::shared_ptr<monty::ndarray<double,1>> & rhs) {task->varboundlist(idxs->size(), idxs->raw(), MSK_BK_LO, rhs->raw(),rhs->raw()); }
    void p_BaseModel::task_var_putboundlist_up(const std::shared_ptr<monty::ndarray<int,1>> idxs, const std::shared_ptr<monty::ndarray<double,1>> & rhs) {task->varboundlist(idxs->size(), idxs->raw(), MSK_BK_UP, rhs->raw(),rhs->raw()); }
    void p_BaseModel::task_var_putboundlist_fx(const std::shared_ptr<monty::ndarray<int,1>> idxs, const std::shared_ptr<monty::ndarray<double,1>> & rhs) {task->varboundlist(idxs->size(), idxs->raw(), MSK_BK_FX, rhs->raw(),rhs->raw()); }
    void p_BaseModel::task_var_putboundlist_ra(const std::shared_ptr<monty::ndarray<int,1>> idxs,
                                               const std::shared_ptr<monty::ndarray<double,1>> & lb ,
                                               const std::shared_ptr<monty::ndarray<double,1>> & ub )
    {
        task->varboundlist(idxs->size(), idxs->raw(), MSK_BK_RA, lb->raw(),ub->raw());
    }


    void p_BaseModel::task_objectivename(const std::string & name) { task->objname(name); }

    void p_BaseModel::task_putobjective
      ( bool                                              maximize,
        const std::shared_ptr<monty::ndarray<int,1>>    & subj,
        const std::shared_ptr<monty::ndarray<double,1>> & cof,
        double                                            cfix)
    {
        std::vector<double> c(task->getnumvar());
        std::vector<int>    cidx(task->getnumvar());

        for (ptrdiff_t i = 0; i < cidx.size();  ++i) cidx[i] = i;
        for (ptrdiff_t i = 0; i < subj->size(); ++i) c[(*subj)[i]] += (*cof)[i];

        if (cidx.size() > 0)
          task->putclist(cidx.size(), & cidx[0], c.data());
        task->putcfix(cfix);
        task->putobjsense(maximize ? MSK_OBJECTIVE_SENSE_MAXIMIZE : MSK_OBJECTIVE_SENSE_MINIMIZE);
    }

    void p_BaseModel::task_putclist
      ( const std::shared_ptr<monty::ndarray<int,1>>    & subj,
        const std::shared_ptr<monty::ndarray<double,1>> & cof)
    {
        task->putclist(subj->size(), subj->raw(),cof->raw());
    }



    void p_BaseModel::task_putobjectivename(const std::string & name)
    {
      task->objname(name);
    }

    MSKtask_t p_BaseModel::task_get()
    {
      return task->get_task();
    }

    MSKtask_t p_BaseModel::__mosek_2fusion_2BaseModel__task_get() { return task_get(); }

    void p_BaseModel::dispose()
    {
        task->dispose();
    }

    //---------------------------
    void p_BaseModel::task_break_solve()
    {
        task->breakOptimize();
    }

    void p_BaseModel::task_putoptserver_host(const std::string & addr)
    {
      task->putoptserverhost(addr);
    }


    void p_BaseModel::report_task_solution(MSKsoltypee st, int numvar, int numcon, int numbarelm, int64_t numacc, int64_t numaccelm) {
        MSKsolstae ssta = task->getsolsta(st);
        MSKprostae psta = task->getprosta(st);
        SolutionStatus psolsta,dsolsta;
        bool
          hasprimal = false,
          hasdual   = false;

        switch (ssta) {
            case MSK_SOL_STA_OPTIMAL:
              psolsta = SolutionStatus::Optimal;
              dsolsta = SolutionStatus::Optimal;
              hasprimal = true;
              hasdual   = true;
              break;
            case MSK_SOL_STA_INTEGER_OPTIMAL:
              psolsta = SolutionStatus::Optimal;
              dsolsta = SolutionStatus::Undefined;
              hasprimal = true;
              break;
            case MSK_SOL_STA_PRIM_AND_DUAL_FEAS:
              psolsta = SolutionStatus::Feasible;
              dsolsta = SolutionStatus::Feasible;
              hasprimal = true;
              hasdual   = true;
              break;
            case MSK_SOL_STA_PRIM_FEAS:
              psolsta = SolutionStatus::Feasible;
              dsolsta = SolutionStatus::Undefined;
              hasprimal = true;
              break;
            case MSK_SOL_STA_DUAL_FEAS:
              psolsta = SolutionStatus::Undefined;
              dsolsta = SolutionStatus::Feasible;
              hasdual   = true;
              break;
            case MSK_SOL_STA_PRIM_INFEAS_CER:
              psolsta = SolutionStatus::Undefined;
              dsolsta = SolutionStatus::Certificate;
              hasdual   = true;
              break;
            case MSK_SOL_STA_DUAL_INFEAS_CER:
              psolsta = SolutionStatus::Certificate;
              dsolsta = SolutionStatus::Undefined;
              hasprimal = true;
              break;
            case MSK_SOL_STA_PRIM_ILLPOSED_CER:
              psolsta = SolutionStatus::Undefined;
              dsolsta = SolutionStatus::IllposedCert;
              hasdual   = true;
              break;
            case MSK_SOL_STA_DUAL_ILLPOSED_CER:
              psolsta = SolutionStatus::IllposedCert;
              dsolsta = SolutionStatus::Undefined;
              hasprimal = true;
              break;
            default:
              psolsta = SolutionStatus::Unknown;
              dsolsta = SolutionStatus::Unknown;
              hasprimal = true;
              hasdual   = true;
              break;
        }

        /* An unconditional override for integer problems */
        if (st == MSK_SOL_ITG)
        {
          dsolsta = SolutionStatus::Undefined;
          hasdual = false;
          if (psta == MSK_PRO_STA_PRIM_INFEAS)
          {
            psolsta = SolutionStatus::Undefined;
            hasprimal = false;
          }
        }

        double pobj     = hasprimal ? task->getprimalobj(st) : 0;
        double dobj     = hasdual ? task->getdualobj(st) : 0;

        ProblemStatus ps;
        switch (psta) {
            case MSK_PRO_STA_UNKNOWN:                    ps = ProblemStatus::Unknown; break;
            case MSK_PRO_STA_PRIM_AND_DUAL_FEAS:         ps = ProblemStatus::PrimalAndDualFeasible; break;
            case MSK_PRO_STA_PRIM_FEAS:                  ps = ProblemStatus::PrimalFeasible; break;
            case MSK_PRO_STA_DUAL_FEAS:                  ps = ProblemStatus::DualFeasible; break;
            case MSK_PRO_STA_PRIM_INFEAS:                ps = ProblemStatus::PrimalInfeasible; break;
            case MSK_PRO_STA_DUAL_INFEAS:                ps = ProblemStatus::DualInfeasible; break;
            case MSK_PRO_STA_PRIM_AND_DUAL_INFEAS:       ps = ProblemStatus::PrimalAndDualInfeasible; break;
            case MSK_PRO_STA_ILL_POSED:                  ps = ProblemStatus::IllPosed; break;
            case MSK_PRO_STA_PRIM_INFEAS_OR_UNBOUNDED:   ps = ProblemStatus::PrimalInfeasibleOrUnbounded; break;
            default:                                     ps = ProblemStatus::Unknown; break;
        }

        SolutionType soltp;
        switch (st) {
            case MSK_SOL_BAS: soltp = SolutionType::Basic;    break;
            case MSK_SOL_ITR: soltp = SolutionType::Interior; break;
            case MSK_SOL_ITG: soltp = SolutionType::Integer;  break;
            default:          soltp = SolutionType::Interior; break;
        }

        cursol = st;
        report_solution(soltp,ps,psolsta,dsolsta,pobj,dobj,numvar,numcon,numbarelm,(int)numacc,(int)numaccelm,hasprimal,hasdual);
    }

    
    std::shared_ptr<monty::ndarray<SolverStatus,1>> p_BaseModel::env_solve_batch(bool israce, 
                                                                                 double timelimit, 
                                                                                 int numthreads, 
                                                                                 std::shared_ptr<monty::ndarray<Model::t,1>> & models)
    {
        int n = models->size();
        std::vector<MSKrescodee> res(n);
        std::vector<MSKrescodee> trm(n);
        std::vector<MSKtask_t> tasks(n);
        std::vector<SolverStatus> wasOK(n);
        MSKrescodee r = MSK_RES_OK;
        bool ok = false;

        for(int i = 0; i < n; i++) 
        {
            res[i] = MSK_RES_ERR_UNKNOWN;
            trm[i] = MSK_RES_ERR_UNKNOWN;
            tasks[i] = (*models)[i]->getTask();
        }
        
        try 
        {
            optimizebatch(israce, timelimit, numthreads, n, tasks.data(), trm.data(), res.data());
        }
        catch (mosek::MosekException e)
        {
            throw OptimizeError(e.what());
        }

        for(int i = 0; i < n; i++)
            wasOK[i] = (res[i] != MSK_RES_OK) ? SolverStatus::Error : ((trm[i] == MSK_RES_TRM_LOST_RACE) ? SolverStatus::LostRace : SolverStatus::OK);

        return monty::new_array_ptr<SolverStatus>(wasOK);
    }

    void p_BaseModel::task_solve(bool remote, const std::string & server, const std::string & port)
    {
      synched = false;
      bool ok = false;

      monty::finally([&]()
        {
          if (! ok) // means exception before we reached end
          {
            clear_solutions();
          }
        });

      try
      {
        if (remote)
          task->optimizermt(server, port);
        else
          task->optimize();
        task->solutionsummary(MSK_STREAM_LOG);
        ok = true;
      }
      catch (mosek::MosekException e)
      {
        throw OptimizeError(e.what());
      }
    }

    void p_BaseModel::task_post_solve() 
    {
      int numcon    = task->getnumcon();
      int numvar    = task->getnumvar();
      int numcone   = task->getnumcone();
      int numbarvar = task->getnumbarvar();
      long numacc    = task->getnumacc();
      long numaccelm = task->getaccntot();
      size_t barvarveclen = 0;
      for (int j = 0; j < numbarvar; ++j) {
          int barjdim = task->getbarvardim(j);
          barvarveclen += barjdim*(1+barjdim)/2;
      }

      bool sol_bas_def = task->solutiondef(MSK_SOL_BAS);
      bool sol_itr_def = task->solutiondef(MSK_SOL_ITR);
      bool sol_itg_def = task->solutiondef(MSK_SOL_ITG);

      clear_solutions();
      if (sol_itr_def != 0)
          report_task_solution(MSK_SOL_ITR, numvar, numcon, (int)barvarveclen,numacc,numaccelm);
      if (sol_bas_def != 0)
          report_task_solution(MSK_SOL_BAS, numvar, numcon, (int)barvarveclen,numacc,numaccelm);
      if (sol_itg_def != 0)
          report_task_solution(MSK_SOL_ITG, numvar, numcon, (int)barvarveclen,numacc,numaccelm);
    }

    void p_BaseModel::task_var_putintlist(const std::shared_ptr<monty::ndarray<int,1>> & idxs) { task->putintlist(idxs->size(), idxs->raw()); }
    void p_BaseModel::task_var_putcontlist(const std::shared_ptr<monty::ndarray<int,1>> & idxs) { task->putcontlist(idxs->size(), idxs->raw()); }

    void p_BaseModel::task_analyze_problem(int detail)
    {
      task->putparam("MSK_IPAR_LOG_ANA_PRO",detail);
      task->analyzeproblem(MSK_STREAM_LOG);
    }

    void p_BaseModel::task_write(const std::string & filename)
    {
      task->putparam("MSK_IPAR_OPF_WRITE_SOLUTIONS",1);
      //task->putparam("MSK_IPAR_WRITE_IGNORE_INCOMPATIBLE_ITEMS",1);
      task->write(filename);
    }

    size_t MSKAPI datawritefunc (MSKuserhandle_t handle, const void * dest, const size_t count) {
        std::ostream * stream = (std::ostream*) handle;
        stream->write((char*) dest, count);
        return count;
    }

    void p_BaseModel::task_write_stream(const std::string & ext, std::ostream & stream)
    {
      task->putparam("MSK_IPAR_OPF_WRITE_SOLUTIONS",1);

      MSKdataformate format = MSK_DATA_FORMAT_MPS;
      MSKcompresstypee compress = MSK_COMPRESS_NONE;

      size_t p = ext.find_first_of('.');
      if (p == std::string::npos) {
          if      (ext.compare("lp") == 0)      format = MSK_DATA_FORMAT_LP;
          else if (ext.compare("ptf") == 0)     format = MSK_DATA_FORMAT_PTF;
          else if (ext.compare("cbf") == 0)     format = MSK_DATA_FORMAT_CB;
          else if (ext.compare("opf") == 0)     format = MSK_DATA_FORMAT_OP;
          else if (ext.compare("mps") == 0)     format = MSK_DATA_FORMAT_MPS;
          else if (ext.compare("jtask") == 0)   format = MSK_DATA_FORMAT_JSON_TASK;
          else if (ext.compare("task") == 0)    format = MSK_DATA_FORMAT_TASK;
      }
      else {
          std::string fmtext = ext.substr(0,p);
          std::string compressext = ext.substr(p+1,ext.size()-p-1);
          if      (fmtext.compare("lp") == 0)      format = MSK_DATA_FORMAT_LP;
          else if (fmtext.compare("ptf") == 0)     format = MSK_DATA_FORMAT_PTF;
          else if (fmtext.compare("cbf") == 0)     format = MSK_DATA_FORMAT_CB;
          else if (fmtext.compare("opf") == 0)     format = MSK_DATA_FORMAT_OP;
          else if (fmtext.compare("mps") == 0)     format = MSK_DATA_FORMAT_MPS;          
          else if (fmtext.compare("jtask") == 0)   format = MSK_DATA_FORMAT_JSON_TASK;
          else if (fmtext.compare("task") == 0)    format = MSK_DATA_FORMAT_TASK;

          if      (compressext.compare("gz") == 0) compress = MSK_COMPRESS_GZIP;
          else if (compressext.compare("zst") == 0) compress = MSK_COMPRESS_ZSTD;
      }

      task->writedatahandle(datawritefunc, &stream, format, compress);
    }

    void p_BaseModel::task_putxx_slice(SolutionType which, int first, int last, std::shared_ptr<monty::ndarray<double,1>> & xx)
    {
        switch (which)
        {
            case SolutionType::Default:
            case SolutionType::Interior: task->putxxslice(MSK_SOL_ITR, first,last,xx->raw()); break;
            case SolutionType::Integer:  task->putxxslice(MSK_SOL_ITG, first,last,xx->raw()); break;
            case SolutionType::Basic:    task->putxxslice(MSK_SOL_BAS, first,last,xx->raw()); break;
       }
    }

    void p_BaseModel::task_setnumvar(int num)
    {
        int numvar = task->getnumvar();
        if (numvar > num)
        {
            for (int i = numvar; i > num; --i)
                task->removevar(i-1);
        }
    }


    void p_BaseModel::task_cleanup (int inumvar, int inumcon, int inumcone, int inumbarvar)
    {
        task->fixvars(inumcone);
        task->revert(task->getnumvar(),inumcon,inumcone,inumbarvar);
    }




#if 0
    void p_BaseModel::env_syeig (int n, std::shared_ptr<monty::ndarray<double,1>> & a, std::shared_ptr<monty::ndarray<double,1>> & w)
    {
      if      (a->size() < n * (n+1) / 2)
        throw MosekException("Invalid length if a in call to syeig");
      else if (w->size() < n)
        throw MosekException("Invalid length if w in call to syeig");
      Task::env_syeig(n,a->raw(),w->raw());
    }

    void p_BaseModel::env_potrf (int n, std::shared_ptr<monty::ndarray<double,1>> & a)
    {
      if      (a->size() < n * (n+1) / 2)
        throw MosekException("Invalid length if a in call to potrf");
      Task::env_potrf(n,a->raw());
    }

    void p_BaseModel::env_syevd (int n, std::shared_ptr<monty::ndarray<double,1>> & a, std::shared_ptr<monty::ndarray<double,1>> & w)
    {
      if      (a->size() < n * (n+1) / 2)
        throw MosekException("Invalid length if a in call to syevd");
      else if (w->size() < n)
        throw MosekException("Invalid length if w in call to syevd");
      Task::env_syevd(n,a->raw(),w->raw());
    }
#endif


    void p_BaseModel::env_putlicensecode (std::shared_ptr<monty::ndarray<int,1>> code) { putlicensecode(code->raw()); }
    void p_BaseModel::env_putlicensepath (const std::string & path) { putlicensepath(path); }
    void p_BaseModel::env_putlicensewait (int wait) { putlicensewait(wait); }

    std::string p_BaseModel::env_getversion() {
       int major; int minor; int revision;
       char buf[20];
       MSK_getversion(&major, &minor, &revision);
       sprintf(buf, "%d.%d.%d", major, minor, revision);
       return std::string(buf);
    };


    int64_t p_BaseModel::task_append_afes(int64_t num)
    {
      return task->appendafes(num);
    }

    void p_BaseModel::task_putafeflist  (array_t<int64_t> idxs, array_t<int> ptr, array_t<int>subj, array_t<double>cof, array_t<double>g) {
        std::vector<int64_t> lptr(ptr->raw(),ptr->raw()+ptr->size()-1); //for (ptrdiff_t i = 0; i < lptr.size(); ++i) lptr[i] = (*ptr)[i];
        std::vector<int>     nzrow(lptr.size()); for (ptrdiff_t i = 0; i < nzrow.size(); ++i) nzrow[i] = (*ptr)[i+1] - (*ptr)[i];
        task->putafefrowlist(idxs->size(), idxs->raw(),nzrow.data(),lptr.data(), subj->raw(),cof->raw());
        task->putafeglist(g->size(),idxs->raw(),g->raw());
    }
    void p_BaseModel::task_putafebarfrowlist (array_t<int> idxs, array_t<int> ptr, array_t<int> barsubj, array_t<int64_t> symmatidx) {
        size_t num = barsubj->size();
        std::vector<int64_t> afeidxlist(idxs->raw(),idxs->raw()+idxs->size());
        for (ptrdiff_t i = 0,k = 0; i < idxs->size(); ++i)
            for (ptrdiff_t j = (*ptr)[i]; j < (*ptr)[i+1]; ++j, ++k)
                afeidxlist[k] = (*idxs)[i];


        std::vector<int64_t> numtermslist(num,1);
        std::vector<int64_t> ptrtermslist(num); for (ptrdiff_t i = 0; i < num; ++i) ptrtermslist[i] = i;
        std::vector<double>  termweights(num, 1.0);

        task->putafebarfentrylist(num,afeidxlist.data(),barsubj->raw(),numtermslist.data(),ptrtermslist.data(),num,symmatidx->raw(),termweights.data());

    }

    void p_BaseModel::task_putafefijlist(array_t<int32_t> & idxs, array_t<int32_t> & subj, array_t<double> & cof) {
        std::vector<int64_t> idxs64(idxs->raw(),idxs->raw()+idxs->size());
        task->putafefentrylist(idxs64.size(), idxs64.data(), subj->raw(),cof->raw());
    }

    void p_BaseModel::task_putafefglist (array_t<int64_t> idxs, array_t<double> g) {
        task->putafeglist(idxs->size(),idxs->raw(),g->raw());
    }
    void p_BaseModel::task_clearafelist (array_t<int64_t>idxs) {
        auto n = idxs->size();
        task->emptyafefrowlist(idxs->size(),idxs->raw());
        task->emptyafebarfrowlist(idxs->size(),idxs->raw());
        std::vector<double> g(n,0.0);
        task->putafeglist(idxs->size(),idxs->raw(),g.data());
    }
    void p_BaseModel::task_putacclist  (array_t<int64_t>idxs, array_t<int64_t>domidxs, array_t<int64_t> afeidxs, array_t<double>g) {
        task->putacclist(idxs->size(),idxs->raw(),domidxs->raw(),afeidxs->size(),afeidxs->raw(),g->raw());
    }
    void p_BaseModel::task_append_accs ( int64_t domidx, int numcone, array_t<int64_t> afeidxs,array_t<double> g) {
        std::vector<int64_t> domidxs(numcone,domidx);
        task->appendaccs(numcone,domidxs.data(),afeidxs->size(),afeidxs->raw(),g->raw());
    }




    int64_t p_BaseModel::task_append_domain_quad     (int conesize ){ return task->appendquadraticconedomain(conesize); }
    int64_t p_BaseModel::task_append_domain_rquad    (int conesize ){ return task->appendrquadraticconedomain(conesize); }
    int64_t p_BaseModel::task_append_domain_pexp     (){ return task->appendprimalexpconedomain(); }
    int64_t p_BaseModel::task_append_domain_dexp     (){ return task->appenddualexpconedomain(); }
    int64_t p_BaseModel::task_append_domain_ppow     (int conesize, array_t<double> alpha){ return task->appendprimalpowerconedomain(conesize,alpha->size(),alpha->raw()); }
    int64_t p_BaseModel::task_append_domain_dpow     (int conesize, array_t<double> alpha){ return task->appenddualpowerconedomain(conesize,alpha->size(),alpha->raw()); }
    // int64_t p_BaseModel::task_append_domain_onenorm  (int conesize ){ return task->appendonenormconedomain(conesize); }
    // int64_t p_BaseModel::task_append_domain_infnorm  (int conesize ){ return task->appendinfnormconedomain(conesize); }
    int64_t p_BaseModel::task_append_domain_pgeomean (int conesize ){ return task->appendprimalgeomeanconedomain(conesize); }
    int64_t p_BaseModel::task_append_domain_dgeomean (int conesize ){ return task->appenddualgeomeanconedomain(conesize); }
    int64_t p_BaseModel::task_append_domain_rpos     (int conesize ){ return task->appendrplusdomain(conesize); }
    int64_t p_BaseModel::task_append_domain_rneg     (int conesize ){ return task->appendrminusdomain(conesize); }
    int64_t p_BaseModel::task_append_domain_r        (int conesize ){ return task->appendrdomain(conesize); }
    int64_t p_BaseModel::task_append_domain_rzero    (int conesize ){ return task->appendrzerodomain(conesize); }
    int64_t p_BaseModel::task_append_domain_svec_psd (int conesize ){ return task->appendsvecpsdconedomain(conesize); }
    int64_t p_BaseModel::task_append_domain_empty    (){ return task->appendrdomain(0); }


    int64_t p_BaseModel::task_append_djc(int64_t num)
    {
      return task->appenddjcs(num);
    }

    void p_BaseModel::task_putdjcslice(int64_t first,
                                       int64_t last,
                                       const std::shared_ptr<monty::ndarray<int64_t,1>> numterm,
                                       const std::shared_ptr<monty::ndarray<int64_t,1>> termsizes,
                                       const std::shared_ptr<monty::ndarray<int64_t,1>> domidxlist,
                                       const std::shared_ptr<monty::ndarray<int64_t,1>> afeidxlist,
                                       const std::shared_ptr<monty::ndarray<double,1>> b)
    {
      task->putdjcslice(first, last,
                        domidxlist->size(), domidxlist->raw(),
                        afeidxlist->size(), afeidxlist->raw(), b->raw(),
                        termsizes->size(), termsizes->raw(),
                        numterm->raw());
    }

  }
}
