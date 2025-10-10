#ifndef __FUSION_P_H__
#define __FUSION_P_H__
#include "monty.h"
#include "mosektask_p.h"
#include "list"
#include "vector"
#include "unordered_map"
#include "fusion.h"
namespace mosek
{
  namespace fusion
  {
    struct p_Disjunction
    {
      Disjunction * _pubthis;
      static mosek::fusion::p_Disjunction* _get_impl(mosek::fusion::Disjunction * _inst){ assert(_inst); assert(_inst->_impl); return _inst->_impl; }
      static mosek::fusion::p_Disjunction * _get_impl(mosek::fusion::Disjunction::t _inst) { return _get_impl(_inst.get()); }
      p_Disjunction(Disjunction * _pubthis);
      virtual ~p_Disjunction() { /* std::cout << "~p_Disjunction" << std::endl;*/ };
      int64_t id{};

      virtual void destroy();

      static Disjunction::t _new_Disjunction(int64_t _7_id);
      void _initialize(int64_t _7_id);
    }; // struct Disjunction;

    struct p_DisjunctionTerms
    {
      DisjunctionTerms * _pubthis;
      static mosek::fusion::p_DisjunctionTerms* _get_impl(mosek::fusion::DisjunctionTerms * _inst){ assert(_inst); assert(_inst->_impl); return _inst->_impl; }
      static mosek::fusion::p_DisjunctionTerms * _get_impl(mosek::fusion::DisjunctionTerms::t _inst) { return _get_impl(_inst.get()); }
      p_DisjunctionTerms(DisjunctionTerms * _pubthis);
      virtual ~p_DisjunctionTerms() { /* std::cout << "~p_DisjunctionTerms" << std::endl;*/ };
      std::shared_ptr< monty::ndarray< monty::rc_ptr< ::mosek::fusion::Term >,1 > > terms{};

      virtual void destroy();

      static DisjunctionTerms::t _new_DisjunctionTerms(monty::rc_ptr< ::mosek::fusion::DisjunctionTerms > _8_terms1,monty::rc_ptr< ::mosek::fusion::ExprRangeDomain > _9_term);
      void _initialize(monty::rc_ptr< ::mosek::fusion::DisjunctionTerms > _8_terms1,monty::rc_ptr< ::mosek::fusion::ExprRangeDomain > _9_term);
      static DisjunctionTerms::t _new_DisjunctionTerms(monty::rc_ptr< ::mosek::fusion::DisjunctionTerms > _10_terms1,monty::rc_ptr< ::mosek::fusion::ExprLinearDomain > _11_term);
      void _initialize(monty::rc_ptr< ::mosek::fusion::DisjunctionTerms > _10_terms1,monty::rc_ptr< ::mosek::fusion::ExprLinearDomain > _11_term);
      static DisjunctionTerms::t _new_DisjunctionTerms(monty::rc_ptr< ::mosek::fusion::DisjunctionTerms > _12_terms1,monty::rc_ptr< ::mosek::fusion::DisjunctionTerms > _13_term2);
      void _initialize(monty::rc_ptr< ::mosek::fusion::DisjunctionTerms > _12_terms1,monty::rc_ptr< ::mosek::fusion::DisjunctionTerms > _13_term2);
      static DisjunctionTerms::t _new_DisjunctionTerms(monty::rc_ptr< ::mosek::fusion::DisjunctionTerms > _14_term1,std::shared_ptr< monty::ndarray< monty::rc_ptr< ::mosek::fusion::Term >,1 > > _15_term2);
      void _initialize(monty::rc_ptr< ::mosek::fusion::DisjunctionTerms > _14_term1,std::shared_ptr< monty::ndarray< monty::rc_ptr< ::mosek::fusion::Term >,1 > > _15_term2);
      static DisjunctionTerms::t _new_DisjunctionTerms(monty::rc_ptr< ::mosek::fusion::DisjunctionTerms > _19_term1,monty::rc_ptr< ::mosek::fusion::Term > _20_term2);
      void _initialize(monty::rc_ptr< ::mosek::fusion::DisjunctionTerms > _19_term1,monty::rc_ptr< ::mosek::fusion::Term > _20_term2);
      static DisjunctionTerms::t _new_DisjunctionTerms(std::shared_ptr< monty::ndarray< monty::rc_ptr< ::mosek::fusion::Term >,1 > > _22_terms);
      void _initialize(std::shared_ptr< monty::ndarray< monty::rc_ptr< ::mosek::fusion::Term >,1 > > _22_terms);
      static DisjunctionTerms::t _new_DisjunctionTerms(monty::rc_ptr< ::mosek::fusion::ExprRangeDomain > _24_term);
      void _initialize(monty::rc_ptr< ::mosek::fusion::ExprRangeDomain > _24_term);
      static DisjunctionTerms::t _new_DisjunctionTerms(monty::rc_ptr< ::mosek::fusion::ExprLinearDomain > _25_term);
      void _initialize(monty::rc_ptr< ::mosek::fusion::ExprLinearDomain > _25_term);
    }; // struct DisjunctionTerms;

    struct p_Term : public /*implements*/ virtual ::mosek::fusion::ExprDomain
    {
      Term * _pubthis;
      static mosek::fusion::p_Term* _get_impl(mosek::fusion::Term * _inst){ assert(_inst); assert(_inst->_impl); return _inst->_impl; }
      static mosek::fusion::p_Term * _get_impl(mosek::fusion::Term::t _inst) { return _get_impl(_inst.get()); }
      p_Term(Term * _pubthis);
      virtual ~p_Term() { /* std::cout << "~p_Term" << std::endl;*/ };
      std::shared_ptr< monty::ndarray< monty::rc_ptr< ::mosek::fusion::DJCDomain >,1 > > domains{};
      std::shared_ptr< monty::ndarray< monty::rc_ptr< ::mosek::fusion::Expression >,1 > > exprs{};

      virtual void destroy();

      static Term::t _new_Term(std::shared_ptr< monty::ndarray< monty::rc_ptr< ::mosek::fusion::Expression >,1 > > _26_elist,std::shared_ptr< monty::ndarray< monty::rc_ptr< ::mosek::fusion::DJCDomain >,1 > > _27_dlist);
      void _initialize(std::shared_ptr< monty::ndarray< monty::rc_ptr< ::mosek::fusion::Expression >,1 > > _26_elist,std::shared_ptr< monty::ndarray< monty::rc_ptr< ::mosek::fusion::DJCDomain >,1 > > _27_dlist);
      static Term::t _new_Term(monty::rc_ptr< ::mosek::fusion::Expression > _30_e,monty::rc_ptr< ::mosek::fusion::DJCDomain > _31_d);
      void _initialize(monty::rc_ptr< ::mosek::fusion::Expression > _30_e,monty::rc_ptr< ::mosek::fusion::DJCDomain > _31_d);
      static Term::t _new_Term(std::shared_ptr< monty::ndarray< monty::rc_ptr< ::mosek::fusion::SimpleTerm >,1 > > _32_t);
      void _initialize(std::shared_ptr< monty::ndarray< monty::rc_ptr< ::mosek::fusion::SimpleTerm >,1 > > _32_t);
      virtual int32_t numaccterms() ;
      virtual int32_t numaccrows() ;
      virtual int32_t num() ;
      virtual int32_t size() ;
      virtual monty::rc_ptr< ::mosek::fusion::Term > __mosek_2fusion_2Term__toDJCTerm() ;
      virtual monty::rc_ptr< ::mosek::fusion::Term > __mosek_2fusion_2ExprDomain__toDJCTerm() { return __mosek_2fusion_2Term__toDJCTerm(); }
    }; // struct Term;

    struct p_SimpleTerm : public ::mosek::fusion::p_Term
    {
      SimpleTerm * _pubthis;
      static mosek::fusion::p_SimpleTerm* _get_impl(mosek::fusion::SimpleTerm * _inst){ return static_cast< mosek::fusion::p_SimpleTerm* >(mosek::fusion::p_Term::_get_impl(_inst)); }
      static mosek::fusion::p_SimpleTerm * _get_impl(mosek::fusion::SimpleTerm::t _inst) { return _get_impl(_inst.get()); }
      p_SimpleTerm(SimpleTerm * _pubthis);
      virtual ~p_SimpleTerm() { /* std::cout << "~p_SimpleTerm" << std::endl;*/ };

      virtual void destroy();

      static SimpleTerm::t _new_SimpleTerm(monty::rc_ptr< ::mosek::fusion::Expression > _41_e,monty::rc_ptr< ::mosek::fusion::DJCDomain > _42_d);
      void _initialize(monty::rc_ptr< ::mosek::fusion::Expression > _41_e,monty::rc_ptr< ::mosek::fusion::DJCDomain > _42_d);
    }; // struct SimpleTerm;

    struct p_DJCDomain
    {
      DJCDomain * _pubthis;
      static mosek::fusion::p_DJCDomain* _get_impl(mosek::fusion::DJCDomain * _inst){ assert(_inst); assert(_inst->_impl); return _inst->_impl; }
      static mosek::fusion::p_DJCDomain * _get_impl(mosek::fusion::DJCDomain::t _inst) { return _get_impl(_inst.get()); }
      p_DJCDomain(DJCDomain * _pubthis);
      virtual ~p_DJCDomain() { /* std::cout << "~p_DJCDomain" << std::endl;*/ };
      mosek::fusion::DJCDomainType dom{};
      int32_t conedim{};
      std::shared_ptr< monty::ndarray< int32_t,1 > > shape{};
      std::shared_ptr< monty::ndarray< double,1 > > par{};
      std::shared_ptr< monty::ndarray< double,1 > > b{};

      virtual void destroy();

      static DJCDomain::t _new_DJCDomain(std::shared_ptr< monty::ndarray< double,1 > > _43_b_,std::shared_ptr< monty::ndarray< double,1 > > _44_par_,std::shared_ptr< monty::ndarray< int32_t,1 > > _45_shape_,mosek::fusion::DJCDomainType _46_dom_);
      void _initialize(std::shared_ptr< monty::ndarray< double,1 > > _43_b_,std::shared_ptr< monty::ndarray< double,1 > > _44_par_,std::shared_ptr< monty::ndarray< int32_t,1 > > _45_shape_,mosek::fusion::DJCDomainType _46_dom_);
      static DJCDomain::t _new_DJCDomain(std::shared_ptr< monty::ndarray< double,1 > > _47_b_,std::shared_ptr< monty::ndarray< double,1 > > _48_par_,std::shared_ptr< monty::ndarray< int32_t,1 > > _49_shape_,int32_t _50_conedim_,mosek::fusion::DJCDomainType _51_dom_);
      void _initialize(std::shared_ptr< monty::ndarray< double,1 > > _47_b_,std::shared_ptr< monty::ndarray< double,1 > > _48_par_,std::shared_ptr< monty::ndarray< int32_t,1 > > _49_shape_,int32_t _50_conedim_,mosek::fusion::DJCDomainType _51_dom_);
      virtual int32_t numaccterms() ;
      virtual int32_t numaccrows() ;
      virtual int32_t size() ;
    }; // struct DJCDomain;

    struct p_DJC
    {
      DJC * _pubthis;
      static mosek::fusion::p_DJC* _get_impl(mosek::fusion::DJC * _inst){ assert(_inst); assert(_inst->_impl); return _inst->_impl; }
      static mosek::fusion::p_DJC * _get_impl(mosek::fusion::DJC::t _inst) { return _get_impl(_inst.get()); }
      p_DJC(DJC * _pubthis);
      virtual ~p_DJC() { /* std::cout << "~p_DJC" << std::endl;*/ };

      virtual void destroy();

      static  monty::rc_ptr< ::mosek::fusion::Term > ANDFromTerms(std::shared_ptr< monty::ndarray< monty::rc_ptr< ::mosek::fusion::Term >,1 > > _54_tlist);
      static  monty::rc_ptr< ::mosek::fusion::Term > AND(std::shared_ptr< monty::ndarray< monty::rc_ptr< ::mosek::fusion::ExprDomain >,1 > > _61_elist);
      static  monty::rc_ptr< ::mosek::fusion::Term > AND(monty::rc_ptr< ::mosek::fusion::SimpleTerm > _63_s1,monty::rc_ptr< ::mosek::fusion::SimpleTerm > _64_s2,monty::rc_ptr< ::mosek::fusion::SimpleTerm > _65_s3);
      static  monty::rc_ptr< ::mosek::fusion::Term > AND(monty::rc_ptr< ::mosek::fusion::SimpleTerm > _66_s1,monty::rc_ptr< ::mosek::fusion::SimpleTerm > _67_s2);
      static  monty::rc_ptr< ::mosek::fusion::Term > AND(monty::rc_ptr< ::mosek::fusion::SimpleTerm > _68_s1);
      static  monty::rc_ptr< ::mosek::fusion::Term > AND(std::shared_ptr< monty::ndarray< monty::rc_ptr< ::mosek::fusion::SimpleTerm >,1 > > _69_slist);
      static  monty::rc_ptr< ::mosek::fusion::SimpleTerm > term(monty::rc_ptr< ::mosek::fusion::Expression > _70_expr,monty::rc_ptr< ::mosek::fusion::RangeDomain > _71_dom);
      static  monty::rc_ptr< ::mosek::fusion::SimpleTerm > term(monty::rc_ptr< ::mosek::fusion::Variable > _82_x,monty::rc_ptr< ::mosek::fusion::RangeDomain > _83_dom);
      static  monty::rc_ptr< ::mosek::fusion::SimpleTerm > term(monty::rc_ptr< ::mosek::fusion::ExprRangeDomain > _84_exprdom);
      static  monty::rc_ptr< ::mosek::fusion::SimpleTerm > term(monty::rc_ptr< ::mosek::fusion::Expression > _85_expr,monty::rc_ptr< ::mosek::fusion::LinearDomain > _86_dom);
      static  monty::rc_ptr< ::mosek::fusion::SimpleTerm > term(monty::rc_ptr< ::mosek::fusion::ExprLinearDomain > _94_exprdom);
      static  monty::rc_ptr< ::mosek::fusion::SimpleTerm > term(monty::rc_ptr< ::mosek::fusion::Variable > _95_x,monty::rc_ptr< ::mosek::fusion::LinearDomain > _96_dom);
    }; // struct DJC;

    // mosek.fusion.BaseModel from file 'src/fusion/cxx/BaseModel_p.h'
    // namespace mosek::fusion
    struct p_BaseModel
    {
      p_BaseModel(BaseModel * _pubthis);
    
      void _initialize( monty::rc_ptr<BaseModel> m);
      void _initialize( const std::string & name,
                        const std::string & licpath);
    
      virtual ~p_BaseModel() { /* std::cout << "~p_BaseModel()" << std::endl;*/  }
    
      static p_BaseModel * _get_impl(Model * _inst) { return _inst->_impl; }
    
      //----------------------
    
      bool synched;
      std::string taskname;
    
      monty::rc_ptr<SolutionStruct> sol_itr;
      monty::rc_ptr<SolutionStruct> sol_itg;
      monty::rc_ptr<SolutionStruct> sol_bas;
    
      MSKsoltypee cursol;
      //---------------------
    
      std::unique_ptr<Task> task;
    
      //---------------------
      template<class T>
      using array_t = std::shared_ptr<monty::ndarray<T,1>>;
    
      virtual void clear_solutions() = 0;
      virtual void report_solution(SolutionType soltype,
                                   ProblemStatus prosta,
                                   SolutionStatus psolsta,
                                   SolutionStatus dsolsta,
                                   double pobj,
                                   double dobj,
                                   int32_t numvar,
                                   int32_t numcon,
                                   int32_t numbarelm,
                                   int32_t numacc,
                                   int32_t numaccelm,
                                   bool hasprimal,
                                   bool hasdual) = 0;
    
      void report_solution_get_xx(array_t<double> v);
      void report_solution_get_slx(array_t<double> v);
      void report_solution_get_sux(array_t<double> v);
      void report_solution_get_xc(array_t<double> v);
      void report_solution_get_slc(array_t<double> v);
      void report_solution_get_suc(array_t<double> v);
      void report_solution_get_barx(array_t<double> v);
      void report_solution_get_bars(array_t<double> v);
      void report_solution_get_accx(array_t<double> v);
      void report_solution_get_accy(array_t<double> v);
      void report_solution_get_accptr(array_t<int32_t> v);
    
      //---------------------
      void task_setLogHandler (const msghandler_t & handler);
      void task_setDataCallbackHandler (const datacbhandler_t & handler);
      void task_setCallbackHandler (const cbhandler_t & handler);
    
      int task_append_barvar(int size, int num);
    
      void task_djc_name   (int64_t index, const std::string & name);
      void task_var_name   (int index, const std::string & name);
      void task_con_name   (int index, const std::string & name);
      void task_barvar_name(int index, const std::string & name);
      void task_objectivename(         const std::string & name);
    
      void task_format_djc_names
      ( const std::shared_ptr<monty::ndarray<int64_t,1>> sub,
        const std::string                              & format,
        const std::shared_ptr<monty::ndarray<int,1>>     dims,
        const std::shared_ptr<monty::ndarray<std::shared_ptr<monty::ndarray<std::string,1>>>> names);
      void task_format_acc_names
      ( const std::shared_ptr<monty::ndarray<int64_t,1>> sub,
        const std::string                              & format,
        const std::shared_ptr<monty::ndarray<int,1>>     dims,
        const std::shared_ptr<monty::ndarray<std::shared_ptr<monty::ndarray<std::string,1>>>> names);
      void task_format_var_names
      ( const std::shared_ptr<monty::ndarray<int,1>>     subj,
        const std::string                              & format,
        const std::shared_ptr<monty::ndarray<int,1>>     dims,
        const std::shared_ptr<monty::ndarray<int64_t,1>> sp,
        const std::shared_ptr<monty::ndarray<std::shared_ptr<monty::ndarray<std::string,1>>>> names);
      void task_format_con_names
      ( const std::shared_ptr<monty::ndarray<int,1>>     subi,
        const std::string                              & format,
        const std::shared_ptr<monty::ndarray<int,1>>     dims,
        const std::shared_ptr<monty::ndarray<int64_t,1>> sp,
        const std::shared_ptr<monty::ndarray<std::shared_ptr<monty::ndarray<std::string,1>>>> names);
      void task_format_barvar_names
      ( const std::shared_ptr<monty::ndarray<int,1>>     subj,
        const std::string                              & format,
        const std::shared_ptr<monty::ndarray<int,1>>     dims,
        const std::shared_ptr<monty::ndarray<std::shared_ptr<monty::ndarray<std::string,1>>>> names);
    
      void task_break_solve();
    
      //--------------------------
    
      int task_numvar();
      int task_numcon();
      int task_numbarvar();
      int task_numacc();
      int task_numdjc();
      int task_numafe();
    
      //--------------------------
    
      void task_put_param(const std::string & name, const std::string & value);
      void task_put_param(const std::string & name, int    value);
      void task_put_param(const std::string & name, double value);
    
      double    task_get_dinf (const std::string & name);
      int       task_get_iinf (const std::string & name);
      int64_t task_get_liinf(const std::string & name);
    
      //--------------------------
    
      void task_con_putboundlist_fr(const std::shared_ptr<monty::ndarray<int,1>> idxs);
      void task_con_putboundlist_lo(const std::shared_ptr<monty::ndarray<int,1>> idxs, const std::shared_ptr<monty::ndarray<double,1>> & rhs);
      void task_con_putboundlist_up(const std::shared_ptr<monty::ndarray<int,1>> idxs, const std::shared_ptr<monty::ndarray<double,1>> & rhs);
      void task_con_putboundlist_fx(const std::shared_ptr<monty::ndarray<int,1>> idxs, const std::shared_ptr<monty::ndarray<double,1>> & rhs);
      void task_con_putboundlist_ra(const std::shared_ptr<monty::ndarray<int,1>> idxs, const std::shared_ptr<monty::ndarray<double,1>> & lb ,
                                    const std::shared_ptr<monty::ndarray<double,1>> & ub );
    
      void task_var_putboundlist_fr(const std::shared_ptr<monty::ndarray<int,1>> idxs);
      void task_var_putboundlist_lo(const std::shared_ptr<monty::ndarray<int,1>> idxs, const std::shared_ptr<monty::ndarray<double,1>> & rhs);
      void task_var_putboundlist_up(const std::shared_ptr<monty::ndarray<int,1>> idxs, const std::shared_ptr<monty::ndarray<double,1>> & rhs);
      void task_var_putboundlist_fx(const std::shared_ptr<monty::ndarray<int,1>> idxs, const std::shared_ptr<monty::ndarray<double,1>> & rhs);
      void task_var_putboundlist_ra(const std::shared_ptr<monty::ndarray<int,1>> idxs, const std::shared_ptr<monty::ndarray<double,1>> & lb ,
                                    const std::shared_ptr<monty::ndarray<double,1>> & ub );
    
      void task_var_putintlist(const std::shared_ptr<monty::ndarray<int,1>> & idxs);
      void task_var_putcontlist(const std::shared_ptr<monty::ndarray<int,1>> & idxs);
    
      //--------------------------
    
      void task_putbararowlist(const std::shared_ptr<monty::ndarray<int,1>>       subi,
                               const std::shared_ptr<monty::ndarray<int64_t,1>> ptr,
                               const std::shared_ptr<monty::ndarray<int,1>>       subj,
                               const std::shared_ptr<monty::ndarray<int64_t,1>> matidx);
    
      void task_putbaraijlist(const std::shared_ptr<monty::ndarray<int,1>> subi,
                              const std::shared_ptr<monty::ndarray<int,1>> subj,
                              std::shared_ptr<monty::ndarray<int64_t,1>> matidx);
    
      void task_putbarc(const std::shared_ptr<monty::ndarray<int,1>> subj,
                        const std::shared_ptr<monty::ndarray<int,1>> subl,
                        const std::shared_ptr<monty::ndarray<int,1>> subk,
                        const std::shared_ptr<monty::ndarray<double,1>> val);
    
      std::shared_ptr<monty::ndarray<int64_t,1>> task_appendsymmatlist (const std::shared_ptr<monty::ndarray<int,1>>       & dim,
                                                                        const std::shared_ptr<monty::ndarray<int64_t,1>> & nz,
                                                                        const std::shared_ptr<monty::ndarray<int,1>>       & subk,
                                                                        const std::shared_ptr<monty::ndarray<int,1>>       & subl,
                                                                        const std::shared_ptr<monty::ndarray<double,1>>    & val);
      int  task_barvar_dim(int j);
      void task_putbaraij (int i, int j, int k);
      void task_putbaraij (int i, int j, const std::shared_ptr<monty::ndarray<int,1>> & k);
      void task_putbarcj  (int j, int k);
      void task_putbarcj  (int j,        const std::shared_ptr<monty::ndarray<int,1>> & k);
      int  task_barvardim (int index);
    
      int task_append_var(int num);
      int task_append_con(int num);
    
      void task_cleararowlist(const std::shared_ptr<monty::ndarray<int,1>> & idxs);
      void task_clearacollist(const std::shared_ptr<monty::ndarray<int,1>> & idxs);
    
      void task_putarowlist(
        const std::shared_ptr<monty::ndarray<int,1>>       & idxs,
        const std::shared_ptr<monty::ndarray<int64_t,1>> & ptrb,
        const std::shared_ptr<monty::ndarray<int,1>>       & subj,
        const std::shared_ptr<monty::ndarray<double,1>>    & cof);
      void task_putaijlist(
        const std::shared_ptr<monty::ndarray<int,1>>       & subi,
        const std::shared_ptr<monty::ndarray<int,1>>       & subj,
        const std::shared_ptr<monty::ndarray<double,1>>    & cof,
        int64_t                           num);
    
      void task_setnumvar(int num);
      void task_cleanup(int oldnum, int oldnumcon, int oldnumcone, int oldnumbarvar);
      void task_putoptserver_host(const std::string & addr);
      void report_task_solution(MSKsoltypee st, int numvar, int numcon, int numbarelm, int64_t numacc, int64_t numaccelm);
    
      void task_solve(bool remote, const std::string & server, const std::string & port);
      void task_post_solve();
      static std::shared_ptr<monty::ndarray<SolverStatus,1>>  env_solve_batch(bool israce, 
                                                                              double timelimit, 
                                                                              int numthreads, 
                                                                              std::shared_ptr<monty::ndarray<Model::t,1>> & models);
    
      void task_putobjective(
        bool                             maximize,
        const std::shared_ptr<monty::ndarray<int,1>>    & subj    ,
        const std::shared_ptr<monty::ndarray<double,1>> & cof     ,
        double                           cfix    );
    
      void task_putclist(
        const std::shared_ptr<monty::ndarray<int,1>>    & subj,
        const std::shared_ptr<monty::ndarray<double,1>> & cof);
    
    
      void task_putobjectivename(const std::string & name);
    
      void task_write(const std::string & filename);
      void task_write_stream(const std::string & ext, std::ostream & stream);
      void task_dump (const std::string & filename);
      void task_analyze_problem(int detail);
    
      MSKtask_t task_get();
      MSKtask_t __mosek_2fusion_2BaseModel__task_get();
    
      void dispose();
    
      void task_putxx_slice(SolutionType which, int first, int last, std::shared_ptr<monty::ndarray<double,1>> & xx);
    
      static void env_syeig (int n, std::shared_ptr<monty::ndarray<double,1>> & a, std::shared_ptr<monty::ndarray<double,1>> & w);
      static void env_potrf (int n, std::shared_ptr<monty::ndarray<double,1>> & a);
      static void env_syevd (int n, std::shared_ptr<monty::ndarray<double,1>> & a, std::shared_ptr<monty::ndarray<double,1>> & w);
    
      static void env_putlicensecode(std::shared_ptr<monty::ndarray<int,1>> code);
      static void env_putlicensepath(const std::string &licfile);
      static void env_putlicensewait(int wait);
    
      static std::string env_getversion();
    
      // void convertSolutionStatus(MSKsoltypee soltype, p_SolutionStruct * sol, MSKsolstae status, MSKprostae prosta);
    
      int64_t task_append_afes (int64_t n);
      void task_putafeflist  (array_t<int64_t> idxs, array_t<int> ptr, array_t<int>subj, array_t<double>cof, array_t<double>g);
      void task_putafebarfrowlist (array_t<int> idxs, array_t<int> ptr, array_t<int> barsubj, array_t<int64_t> symmatidx);
      void task_putafefijlist (array_t<int> &idxs, array_t<int> &subj, array_t<double> &cof);
      void task_putafefglist (array_t<int64_t> idxs, array_t<double> g);
      void task_clearafelist (array_t<int64_t>idxs);
      void task_putacclist  (array_t<int64_t>idxs, array_t<int64_t>domidxs, array_t<int64_t>afeidxs_t,array_t<double>g);
      void task_append_accs ( int64_t domidx, int numcone,array_t<int64_t> afeidxs,array_t<double> g);
    
      int64_t task_append_domain_quad     (int conesize);
      int64_t task_append_domain_rquad    (int conesize);
      int64_t task_append_domain_pexp     ();
      int64_t task_append_domain_ppow     (int conesize,array_t<double> alpha);
      int64_t task_append_domain_dexp     ();
      int64_t task_append_domain_dpow     (int conesize,array_t<double> alpha);
      /* int64_t task_append_domain_onenorm  (int conesize); */
      /* int64_t task_append_domain_infnorm  (int conesize); */
      int64_t task_append_domain_pgeomean (int conesize);
      int64_t task_append_domain_dgeomean (int conesize);
      int64_t task_append_domain_rpos     (int conesize);
      int64_t task_append_domain_rneg     (int conesize);
      int64_t task_append_domain_r        (int conesize);
      int64_t task_append_domain_rzero    (int conesize);
      int64_t task_append_domain_svec_psd (int conesize);
      int64_t task_append_domain_empty    ();
      int64_t task_append_djc             (int64_t n);
      void task_putdjcslice(int64_t first, int64_t last , array_t<int64_t> numterm_t, array_t<int64_t> termsizes, array_t<int64_t> domidxlist, array_t<int64_t> afeidxlist,  array_t<double> b);
    
    };
    // End of file 'src/fusion/cxx/BaseModel_p.h'
    struct p_Model : public ::mosek::fusion::p_BaseModel
    {
      Model * _pubthis;
      static mosek::fusion::p_Model* _get_impl(mosek::fusion::Model * _inst){ return static_cast< mosek::fusion::p_Model* >(mosek::fusion::p_BaseModel::_get_impl(_inst)); }
      static mosek::fusion::p_Model * _get_impl(mosek::fusion::Model::t _inst) { return _get_impl(_inst.get()); }
      p_Model(Model * _pubthis);
      virtual ~p_Model() { /* std::cout << "~p_Model" << std::endl;*/ };
      monty::rc_ptr< ::mosek::fusion::WorkStack > xs{};
      monty::rc_ptr< ::mosek::fusion::WorkStack > ws{};
      monty::rc_ptr< ::mosek::fusion::WorkStack > rs{};
      monty::rc_ptr< ::mosek::fusion::SolutionStruct > sol_itg{};
      monty::rc_ptr< ::mosek::fusion::SolutionStruct > sol_bas{};
      monty::rc_ptr< ::mosek::fusion::SolutionStruct > sol_itr{};
      monty::rc_ptr< ::mosek::fusion::Utils::StringIntMap > con_map{};
      std::shared_ptr< monty::ndarray< monty::rc_ptr< ::mosek::fusion::ModelConstraint >,1 > > acons{};
      std::shared_ptr< monty::ndarray< monty::rc_ptr< ::mosek::fusion::ModelConstraint >,1 > > cons{};
      int64_t task_numaferow{};
      std::shared_ptr< monty::ndarray< double,1 > > param_value{};
      int32_t param_num{};
      monty::rc_ptr< ::mosek::fusion::Utils::StringIntMap > par_map{};
      int32_t numparameter{};
      std::shared_ptr< monty::ndarray< monty::rc_ptr< ::mosek::fusion::Parameter >,1 > > parameters{};
      std::shared_ptr< monty::ndarray< bool,1 > > initsol_xx_flag{};
      std::shared_ptr< monty::ndarray< double,1 > > initsol_xx{};
      monty::rc_ptr< ::mosek::fusion::Utils::StringIntMap > var_map{};
      std::shared_ptr< monty::ndarray< monty::rc_ptr< ::mosek::fusion::ModelVariable >,1 > > barvars{};
      std::shared_ptr< monty::ndarray< monty::rc_ptr< ::mosek::fusion::ModelVariable >,1 > > vars{};
      int32_t bfixidx{};
      std::shared_ptr< monty::ndarray< int32_t,1 > > barvar_block_elm_j{};
      std::shared_ptr< monty::ndarray< int32_t,1 > > barvar_block_elm_i{};
      std::shared_ptr< monty::ndarray< int32_t,1 > > barvar_block_elm_barj{};
      std::shared_ptr< monty::ndarray< int32_t,1 > > barvar_block_elm_ptr{};
      std::shared_ptr< monty::ndarray< int32_t,1 > > barvar_block_dim{};
      std::shared_ptr< monty::ndarray< int32_t,1 > > barvar_block_ptr{};
      std::shared_ptr< monty::ndarray< int32_t,1 > > barvar_dim{};
      int32_t barvar_num{};
      std::shared_ptr< monty::ndarray< int32_t,1 > > var_elm_acc_ofs{};
      std::shared_ptr< monty::ndarray< int32_t,1 > > var_elm_acc_idx{};
      std::shared_ptr< monty::ndarray< int32_t,1 > > var_block_acc_id{};
      monty::rc_ptr< ::mosek::fusion::LinkedBlocks > var_block_map{};
      std::shared_ptr< monty::ndarray< int32_t,1 > > acon_elm_afe{};
      std::shared_ptr< monty::ndarray< int32_t,1 > > acon_elm_ofs{};
      std::shared_ptr< monty::ndarray< double,1 > > acon_elm_scale{};
      std::shared_ptr< monty::ndarray< int32_t,1 > > acon_elm_accid{};
      std::shared_ptr< monty::ndarray< int32_t,1 > > acon_afe{};
      std::shared_ptr< monty::ndarray< int32_t,1 > > acon_acc{};
      monty::rc_ptr< ::mosek::fusion::LinkedBlocks > acon_block_map{};
      monty::rc_ptr< ::mosek::fusion::LinkedBlocks > acc_block_map{};
      monty::rc_ptr< ::mosek::fusion::RowBlockManager > obj_blocks{};
      monty::rc_ptr< ::mosek::fusion::RowBlockManager > afe_blocks{};
      monty::rc_ptr< ::mosek::fusion::RowBlockManager > con_blocks{};
      int32_t num_task_acc{};
      int32_t num_task_afe{};
      int32_t num_task_con{};
      mosek::fusion::SolutionType solutionptr{};
      mosek::fusion::AccSolutionStatus acceptable_sol{};
      std::string model_name{};

      virtual void destroy();

      static Model::t _new_Model(monty::rc_ptr< ::mosek::fusion::Model > _697_m);
      void _initialize(monty::rc_ptr< ::mosek::fusion::Model > _697_m);
      static Model::t _new_Model(const std::string &  _704_name,int32_t _705_basesize);
      void _initialize(const std::string &  _704_name,int32_t _705_basesize);
      static Model::t _new_Model(int32_t _713_basesize);
      void _initialize(int32_t _713_basesize);
      static Model::t _new_Model(const std::string &  _714_name);
      void _initialize(const std::string &  _714_name);
      static Model::t _new_Model();
      void _initialize();
      virtual monty::rc_ptr< ::mosek::fusion::Disjunction > __mosek_2fusion_2Model__disjunction(const std::string &  _715_name,std::shared_ptr< monty::ndarray< monty::rc_ptr< ::mosek::fusion::Term >,1 > > _716_terms) ;
      virtual monty::rc_ptr< ::mosek::fusion::Disjunction > __mosek_2fusion_2Model__disjunction(std::shared_ptr< monty::ndarray< monty::rc_ptr< ::mosek::fusion::Term >,1 > > _762_terms) ;
      virtual monty::rc_ptr< ::mosek::fusion::Disjunction > __mosek_2fusion_2Model__disjunction(monty::rc_ptr< ::mosek::fusion::DisjunctionTerms > _763_djcterms) ;
      virtual monty::rc_ptr< ::mosek::fusion::Disjunction > __mosek_2fusion_2Model__disjunction(const std::string &  _764_name,monty::rc_ptr< ::mosek::fusion::DisjunctionTerms > _765_djcterms) ;
      virtual monty::rc_ptr< ::mosek::fusion::Disjunction > __mosek_2fusion_2Model__disjunction(const std::string &  _766_name,std::shared_ptr< monty::ndarray< monty::rc_ptr< ::mosek::fusion::ExprDomain >,1 > > _767_terms) ;
      virtual monty::rc_ptr< ::mosek::fusion::Disjunction > __mosek_2fusion_2Model__disjunction(std::shared_ptr< monty::ndarray< monty::rc_ptr< ::mosek::fusion::ExprDomain >,1 > > _769_terms) ;
      virtual monty::rc_ptr< ::mosek::fusion::Disjunction > __mosek_2fusion_2Model__disjunction(const std::string &  _771_name,monty::rc_ptr< ::mosek::fusion::ExprDomain > _772_term) ;
      virtual monty::rc_ptr< ::mosek::fusion::Disjunction > __mosek_2fusion_2Model__disjunction(monty::rc_ptr< ::mosek::fusion::ExprDomain > _773_term) ;
      virtual monty::rc_ptr< ::mosek::fusion::Disjunction > __mosek_2fusion_2Model__disjunction(monty::rc_ptr< ::mosek::fusion::Term > _774_t1,monty::rc_ptr< ::mosek::fusion::Term > _775_t2,monty::rc_ptr< ::mosek::fusion::Term > _776_t3) ;
      virtual monty::rc_ptr< ::mosek::fusion::Disjunction > __mosek_2fusion_2Model__disjunction(monty::rc_ptr< ::mosek::fusion::Term > _777_t1,monty::rc_ptr< ::mosek::fusion::Term > _778_t2) ;
      virtual monty::rc_ptr< ::mosek::fusion::Disjunction > __mosek_2fusion_2Model__disjunction(monty::rc_ptr< ::mosek::fusion::Term > _779_t1) ;
      virtual monty::rc_ptr< ::mosek::fusion::Disjunction > __mosek_2fusion_2Model__disjunction(const std::string &  _780_name,monty::rc_ptr< ::mosek::fusion::Term > _781_t1,monty::rc_ptr< ::mosek::fusion::Term > _782_t2,monty::rc_ptr< ::mosek::fusion::Term > _783_t3) ;
      virtual monty::rc_ptr< ::mosek::fusion::Disjunction > __mosek_2fusion_2Model__disjunction(const std::string &  _784_name,monty::rc_ptr< ::mosek::fusion::Term > _785_t1,monty::rc_ptr< ::mosek::fusion::Term > _786_t2) ;
      virtual monty::rc_ptr< ::mosek::fusion::Disjunction > __mosek_2fusion_2Model__disjunction(const std::string &  _787_name,monty::rc_ptr< ::mosek::fusion::Term > _788_t1) ;
      virtual monty::rc_ptr< ::mosek::fusion::Utils::StringBuffer > __mosek_2fusion_2Model__formstConstr(monty::rc_ptr< ::mosek::fusion::Utils::StringBuffer > _789_sb,std::shared_ptr< monty::ndarray< int32_t,1 > > _790_shape,std::shared_ptr< monty::ndarray< int32_t,1 > > _791_idxs) ;
      virtual void acon_release(int32_t _792_id) ;
      virtual int32_t acon_allocate(int64_t _800_domidx,int32_t _801_conesize,int32_t _802_numcone,std::shared_ptr< monty::ndarray< double,1 > > _803_g,std::shared_ptr< monty::ndarray< int32_t,1 > > _804_afeidxs,std::shared_ptr< monty::ndarray< int32_t,1 > > _805_accidxs) ;
      virtual void afe_release(int32_t _831_id) ;
      virtual int32_t afe_allocate(std::shared_ptr< monty::ndarray< int32_t,1 > > _834_nativeidxs) ;
      virtual void con_release(int32_t _840_id) ;
      virtual int32_t con_allocate(std::shared_ptr< monty::ndarray< int32_t,1 > > _843_nativeidxs) ;
      virtual int32_t barvar_alloc(int32_t _851_conedim,int32_t _852_numcone,std::shared_ptr< monty::ndarray< int32_t,1 > > _853_barvaridxs,std::shared_ptr< monty::ndarray< int64_t,1 > > _854_nativeidxs) ;
      virtual int32_t conicvar_alloc(int64_t _885_domidx,int32_t _886_conesize,int32_t _887_numcone,std::shared_ptr< monty::ndarray< int32_t,1 > > _888_accidxs,std::shared_ptr< monty::ndarray< int32_t,1 > > _889_nativeidxs) ;
      virtual int32_t linearvar_alloc(int32_t _901_n,std::shared_ptr< monty::ndarray< int32_t,1 > > _902_nativeidxs) ;
      virtual void make_continuous(std::shared_ptr< monty::ndarray< int64_t,1 > > _914_idxs) ;
      virtual void make_integer(std::shared_ptr< monty::ndarray< int64_t,1 > > _920_idxs) ;
      static  void putlicensewait(bool _926_wait);
      static  void putlicensepath(const std::string &  _927_licfile);
      static  void putlicensecode(std::shared_ptr< monty::ndarray< int32_t,1 > > _928_code);
      virtual /* override */ void dispose() ;
      virtual MSKtask_t __mosek_2fusion_2Model__getTask() ;
      virtual void getConstraintDuals(bool _934_lower,std::shared_ptr< monty::ndarray< int32_t,1 > > _935_nativeidxs,std::shared_ptr< monty::ndarray< double,1 > > _936_res,int32_t _937_offset) ;
      virtual void getConstraintValues(bool _942_primal,std::shared_ptr< monty::ndarray< int32_t,1 > > _943_nativeidxs,std::shared_ptr< monty::ndarray< double,1 > > _944_res,int32_t _945_offset) ;
      virtual void getVariableDuals(bool _957_lower,std::shared_ptr< monty::ndarray< int64_t,1 > > _958_nativeidxs,std::shared_ptr< monty::ndarray< double,1 > > _959_res,int32_t _960_offset) ;
      virtual void getVariableValues(bool _966_primal,std::shared_ptr< monty::ndarray< int64_t,1 > > _967_nativeidxs,std::shared_ptr< monty::ndarray< double,1 > > _968_res,int32_t _969_offset) ;
      virtual void setVariableValues(bool _981_primal,std::shared_ptr< monty::ndarray< int64_t,1 > > _982_nativeidxs,std::shared_ptr< monty::ndarray< double,1 > > _983_values) ;
      virtual void flushNames() ;
      virtual void writeTaskNoFlush(const std::string &  _994_filename) ;
      virtual void writeTaskStream(const std::string &  _995_ext,std::ostream&  _996_stream) ;
      virtual void dataReport() ;
      virtual void dataReport(int32_t _997_detail) ;
      virtual void writeTask(const std::string &  _998_filename) ;
      virtual int64_t getSolverLIntInfo(const std::string &  _999_name) ;
      virtual int32_t getSolverIntInfo(const std::string &  _1000_name) ;
      virtual double getSolverDoubleInfo(const std::string &  _1001_name) ;
      virtual void setCallbackHandler(mosek::cbhandler_t  _1002_h) ;
      virtual void setDataCallbackHandler(mosek::datacbhandler_t  _1003_h) ;
      virtual void setLogHandler(mosek::msghandler_t  _1004_h) ;
      virtual void setSolverParam(const std::string &  _1005_name,double _1006_floatval) ;
      virtual void setSolverParam(const std::string &  _1007_name,int32_t _1008_intval) ;
      virtual void setSolverParam(const std::string &  _1009_name,const std::string &  _1010_strval) ;
      virtual void breakSolver() ;
      virtual void optserverHost(const std::string &  _1011_addr) ;
      virtual /* override */ void report_solution(mosek::fusion::SolutionType _1012_soltype,mosek::fusion::ProblemStatus _1013_prosta,mosek::fusion::SolutionStatus _1014_psolsta,mosek::fusion::SolutionStatus _1015_dsolsta,double _1016_pobj,double _1017_dobj,int32_t _1018_numvar,int32_t _1019_numcon,int32_t _1020_numbarelm,int32_t _1021_numacc,int32_t _1022_numaccelm,bool _1023_hasprimal,bool _1024_hasdual) ;
      virtual /* override */ void clear_solutions() ;
      static  std::shared_ptr< monty::ndarray< mosek::fusion::SolverStatus,1 > > solveBatch(bool _1034_israce,double _1035_maxtime,int32_t _1036_numthreads,std::shared_ptr< monty::ndarray< monty::rc_ptr< ::mosek::fusion::Model >,1 > > _1037_models);
      virtual void solve(const std::string &  _1043_addr,const std::string &  _1044_accesstoken) ;
      virtual void solve() ;
      virtual void flush_parameters() ;
      virtual void flushParameters() ;
      virtual void evaluate_parameterized(monty::rc_ptr< ::mosek::fusion::WorkStack > _1057_xs,int32_t _1058_numrow,std::shared_ptr< monty::ndarray< int32_t,1 > > _1059_rowptrb,std::shared_ptr< monty::ndarray< int32_t,1 > > _1060_rowptre,std::shared_ptr< monty::ndarray< int64_t,1 > > _1061_codenidx,std::shared_ptr< monty::ndarray< int32_t,1 > > _1062_codeptr,std::shared_ptr< monty::ndarray< int32_t,1 > > _1063_codesizes,std::shared_ptr< monty::ndarray< int32_t,1 > > _1064_code,std::shared_ptr< monty::ndarray< double,1 > > _1065_cconst,std::shared_ptr< monty::ndarray< int32_t,1 > > _1066_subj,std::shared_ptr< monty::ndarray< double,1 > > _1067_val) ;
      virtual void flushSolutions() ;
      virtual void flush_initsol(mosek::fusion::SolutionType _1078_which) ;
      virtual mosek::fusion::SolutionStatus getDualSolutionStatus() ;
      virtual mosek::fusion::ProblemStatus getProblemStatus() ;
      virtual mosek::fusion::SolutionStatus getPrimalSolutionStatus() ;
      virtual double dualObjValue() ;
      virtual double primalObjValue() ;
      virtual monty::rc_ptr< ::mosek::fusion::SolutionStruct > __mosek_2fusion_2Model__get_sol_cache(mosek::fusion::SolutionType _1085_which_,bool _1086_primal,bool _1087_nothrow) ;
      virtual monty::rc_ptr< ::mosek::fusion::SolutionStruct > __mosek_2fusion_2Model__get_sol_cache(mosek::fusion::SolutionType _1093_which_,bool _1094_primal) ;
      virtual void setSolution_xx(std::shared_ptr< monty::ndarray< int32_t,1 > > _1095_subj,std::shared_ptr< monty::ndarray< double,1 > > _1096_val) ;
      virtual void ensure_initsol_xx() ;
      virtual std::shared_ptr< monty::ndarray< int32_t,1 > > getSolution_accptr(mosek::fusion::SolutionType _1103_which) ;
      virtual std::shared_ptr< monty::ndarray< double,1 > > getSolution_accy(mosek::fusion::SolutionType _1104_which) ;
      virtual std::shared_ptr< monty::ndarray< double,1 > > getSolution_accx(mosek::fusion::SolutionType _1105_which) ;
      virtual std::shared_ptr< monty::ndarray< double,1 > > getSolution_bars(mosek::fusion::SolutionType _1106_which) ;
      virtual std::shared_ptr< monty::ndarray< double,1 > > getSolution_barx(mosek::fusion::SolutionType _1107_which) ;
      virtual std::shared_ptr< monty::ndarray< double,1 > > getSolution_y(mosek::fusion::SolutionType _1108_which) ;
      virtual std::shared_ptr< monty::ndarray< double,1 > > getSolution_xc(mosek::fusion::SolutionType _1109_which) ;
      virtual std::shared_ptr< monty::ndarray< double,1 > > getSolution_suc(mosek::fusion::SolutionType _1110_which) ;
      virtual std::shared_ptr< monty::ndarray< double,1 > > getSolution_slc(mosek::fusion::SolutionType _1111_which) ;
      virtual std::shared_ptr< monty::ndarray< double,1 > > getSolution_sux(mosek::fusion::SolutionType _1112_which) ;
      virtual std::shared_ptr< monty::ndarray< double,1 > > getSolution_slx(mosek::fusion::SolutionType _1113_which) ;
      virtual std::shared_ptr< monty::ndarray< double,1 > > getSolution_yx(mosek::fusion::SolutionType _1114_which) ;
      virtual std::shared_ptr< monty::ndarray< double,1 > > getSolution_xx(mosek::fusion::SolutionType _1115_which) ;
      virtual void selectedSolution(mosek::fusion::SolutionType _1116_soltype) ;
      virtual mosek::fusion::AccSolutionStatus getAcceptedSolutionStatus() ;
      virtual void acceptedSolutionStatus(mosek::fusion::AccSolutionStatus _1117_what) ;
      virtual mosek::fusion::ProblemStatus getProblemStatus(mosek::fusion::SolutionType _1118_which) ;
      virtual mosek::fusion::SolutionStatus getDualSolutionStatus(mosek::fusion::SolutionType _1120_which) ;
      virtual mosek::fusion::SolutionStatus getPrimalSolutionStatus(mosek::fusion::SolutionType _1121_which) ;
      virtual mosek::fusion::SolutionStatus getSolutionStatus(mosek::fusion::SolutionType _1122_which,bool _1123_primal) ;
      virtual void update(std::shared_ptr< monty::ndarray< int32_t,1 > > _1126_conidxs,monty::rc_ptr< ::mosek::fusion::Expression > _1127_expr) ;
      virtual void update(std::shared_ptr< monty::ndarray< int32_t,1 > > _1194_conidxs,monty::rc_ptr< ::mosek::fusion::Expression > _1195_expr,std::shared_ptr< monty::ndarray< int32_t,1 > > _1196_varidxs) ;
      virtual void updateObjective(monty::rc_ptr< ::mosek::fusion::Expression > _1298_expr,monty::rc_ptr< ::mosek::fusion::Variable > _1299_x) ;
      virtual monty::rc_ptr< ::mosek::fusion::Parameter > __mosek_2fusion_2Model__parameter_unchecked(const std::string &  _1336_name,std::shared_ptr< monty::ndarray< int32_t,1 > > _1337_shape,std::shared_ptr< monty::ndarray< int64_t,1 > > _1338_sp) ;
      virtual monty::rc_ptr< ::mosek::fusion::Parameter > __mosek_2fusion_2Model__parameter_(const std::string &  _1348_name,std::shared_ptr< monty::ndarray< int32_t,1 > > _1349_shape,std::shared_ptr< monty::ndarray< int64_t,1 > > _1350_sp) ;
      virtual monty::rc_ptr< ::mosek::fusion::Parameter > __mosek_2fusion_2Model__parameter_(const std::string &  _1355_name,std::shared_ptr< monty::ndarray< int32_t,1 > > _1356_shape,std::shared_ptr< monty::ndarray< int32_t,2 > > _1357_sparsity) ;
      virtual monty::rc_ptr< ::mosek::fusion::Parameter > __mosek_2fusion_2Model__parameter(const std::string &  _1365_name) ;
      virtual monty::rc_ptr< ::mosek::fusion::Parameter > __mosek_2fusion_2Model__parameter(const std::string &  _1367_name,int32_t _1368_d1,int32_t _1369_d2,int32_t _1370_d3) ;
      virtual monty::rc_ptr< ::mosek::fusion::Parameter > __mosek_2fusion_2Model__parameter(const std::string &  _1372_name,int32_t _1373_d1,int32_t _1374_d2) ;
      virtual monty::rc_ptr< ::mosek::fusion::Parameter > __mosek_2fusion_2Model__parameter(const std::string &  _1376_name,int32_t _1377_d1) ;
      virtual monty::rc_ptr< ::mosek::fusion::Parameter > __mosek_2fusion_2Model__parameter(const std::string &  _1379_name,std::shared_ptr< monty::ndarray< int32_t,1 > > _1380_shape) ;
      virtual monty::rc_ptr< ::mosek::fusion::Parameter > __mosek_2fusion_2Model__parameter(const std::string &  _1382_name,std::shared_ptr< monty::ndarray< int32_t,1 > > _1383_shape,std::shared_ptr< monty::ndarray< int64_t,1 > > _1384_sp) ;
      virtual monty::rc_ptr< ::mosek::fusion::Parameter > __mosek_2fusion_2Model__parameter(const std::string &  _1385_name,std::shared_ptr< monty::ndarray< int32_t,1 > > _1386_shape,std::shared_ptr< monty::ndarray< int32_t,2 > > _1387_sparsity) ;
      virtual monty::rc_ptr< ::mosek::fusion::Parameter > __mosek_2fusion_2Model__parameter() ;
      virtual monty::rc_ptr< ::mosek::fusion::Parameter > __mosek_2fusion_2Model__parameter(int32_t _1389_d1,int32_t _1390_d2,int32_t _1391_d3) ;
      virtual monty::rc_ptr< ::mosek::fusion::Parameter > __mosek_2fusion_2Model__parameter(int32_t _1393_d1,int32_t _1394_d2) ;
      virtual monty::rc_ptr< ::mosek::fusion::Parameter > __mosek_2fusion_2Model__parameter(int32_t _1396_d1) ;
      virtual monty::rc_ptr< ::mosek::fusion::Parameter > __mosek_2fusion_2Model__parameter(std::shared_ptr< monty::ndarray< int32_t,1 > > _1398_shape) ;
      virtual monty::rc_ptr< ::mosek::fusion::Parameter > __mosek_2fusion_2Model__parameter(std::shared_ptr< monty::ndarray< int32_t,1 > > _1400_shape,std::shared_ptr< monty::ndarray< int64_t,1 > > _1401_sp) ;
      virtual monty::rc_ptr< ::mosek::fusion::Parameter > __mosek_2fusion_2Model__parameter(std::shared_ptr< monty::ndarray< int32_t,1 > > _1402_shape,std::shared_ptr< monty::ndarray< int32_t,2 > > _1403_sparsity) ;
      virtual void objective_(const std::string &  _1404_name,mosek::fusion::ObjectiveSense _1405_sense,monty::rc_ptr< ::mosek::fusion::Expression > _1406_expr) ;
      virtual void objective(double _1440_c) ;
      virtual void objective(mosek::fusion::ObjectiveSense _1441_sense,double _1442_c) ;
      virtual void objective(mosek::fusion::ObjectiveSense _1443_sense,monty::rc_ptr< ::mosek::fusion::Expression > _1444_expr) ;
      virtual void objective(const std::string &  _1445_name,double _1446_c) ;
      virtual void objective(const std::string &  _1447_name,mosek::fusion::ObjectiveSense _1448_sense,double _1449_c) ;
      virtual void objective(const std::string &  _1450_name,mosek::fusion::ObjectiveSense _1451_sense,monty::rc_ptr< ::mosek::fusion::Expression > _1452_expr) ;
      virtual monty::rc_ptr< ::mosek::fusion::RangedConstraint > __mosek_2fusion_2Model__constraint(monty::rc_ptr< ::mosek::fusion::ExprRangeDomain > _1453_exprdom) ;
      virtual monty::rc_ptr< ::mosek::fusion::RangedConstraint > __mosek_2fusion_2Model__constraint(const std::string &  _1454_name,monty::rc_ptr< ::mosek::fusion::ExprRangeDomain > _1455_exprdom) ;
      virtual monty::rc_ptr< ::mosek::fusion::RangedConstraint > __mosek_2fusion_2Model__constraint(monty::rc_ptr< ::mosek::fusion::Expression > _1456_expr,monty::rc_ptr< ::mosek::fusion::RangeDomain > _1457_rdom) ;
      virtual monty::rc_ptr< ::mosek::fusion::RangedConstraint > __mosek_2fusion_2Model__constraint(const std::string &  _1458_name,monty::rc_ptr< ::mosek::fusion::Expression > _1459_expr,monty::rc_ptr< ::mosek::fusion::RangeDomain > _1460_rdom) ;
      virtual monty::rc_ptr< ::mosek::fusion::Constraint > __mosek_2fusion_2Model__constraint(monty::rc_ptr< ::mosek::fusion::ExprConicDomain > _1461_exprdom) ;
      virtual monty::rc_ptr< ::mosek::fusion::Constraint > __mosek_2fusion_2Model__constraint(const std::string &  _1462_name,monty::rc_ptr< ::mosek::fusion::ExprConicDomain > _1463_exprdom) ;
      virtual monty::rc_ptr< ::mosek::fusion::Constraint > __mosek_2fusion_2Model__constraint(monty::rc_ptr< ::mosek::fusion::Expression > _1464_expr,monty::rc_ptr< ::mosek::fusion::ConeDomain > _1465_qdom) ;
      virtual monty::rc_ptr< ::mosek::fusion::Constraint > __mosek_2fusion_2Model__constraint(const std::string &  _1466_name,monty::rc_ptr< ::mosek::fusion::Expression > _1467_expr,monty::rc_ptr< ::mosek::fusion::ConeDomain > _1468_qdom) ;
      virtual monty::rc_ptr< ::mosek::fusion::Constraint > __mosek_2fusion_2Model__constraint(monty::rc_ptr< ::mosek::fusion::ExprLinearDomain > _1469_exprdom) ;
      virtual monty::rc_ptr< ::mosek::fusion::Constraint > __mosek_2fusion_2Model__constraint(const std::string &  _1470_name,monty::rc_ptr< ::mosek::fusion::ExprLinearDomain > _1471_exprdom) ;
      virtual monty::rc_ptr< ::mosek::fusion::Constraint > __mosek_2fusion_2Model__constraint(monty::rc_ptr< ::mosek::fusion::Expression > _1472_expr,monty::rc_ptr< ::mosek::fusion::LinearDomain > _1473_ldom) ;
      virtual monty::rc_ptr< ::mosek::fusion::Constraint > __mosek_2fusion_2Model__constraint(const std::string &  _1474_name,monty::rc_ptr< ::mosek::fusion::Expression > _1475_expr,monty::rc_ptr< ::mosek::fusion::LinearDomain > _1476_ldom) ;
      virtual monty::rc_ptr< ::mosek::fusion::Constraint > __mosek_2fusion_2Model__constraint(monty::rc_ptr< ::mosek::fusion::ExprPSDDomain > _1477_exprdom) ;
      virtual monty::rc_ptr< ::mosek::fusion::Constraint > __mosek_2fusion_2Model__constraint(const std::string &  _1478_name,monty::rc_ptr< ::mosek::fusion::ExprPSDDomain > _1479_exprdom) ;
      virtual monty::rc_ptr< ::mosek::fusion::Constraint > __mosek_2fusion_2Model__constraint(monty::rc_ptr< ::mosek::fusion::Expression > _1480_expr,monty::rc_ptr< ::mosek::fusion::PSDDomain > _1481_psddom) ;
      virtual monty::rc_ptr< ::mosek::fusion::Constraint > __mosek_2fusion_2Model__constraint(const std::string &  _1482_name,monty::rc_ptr< ::mosek::fusion::Expression > _1483_expr,monty::rc_ptr< ::mosek::fusion::PSDDomain > _1484_psddom) ;
      virtual monty::rc_ptr< ::mosek::fusion::Variable > __mosek_2fusion_2Model__variable(monty::rc_ptr< ::mosek::fusion::PSDDomain > _1485_psddom) ;
      virtual monty::rc_ptr< ::mosek::fusion::Variable > __mosek_2fusion_2Model__variable(int32_t _1486_n,int32_t _1487_m,monty::rc_ptr< ::mosek::fusion::PSDDomain > _1488_psddom) ;
      virtual monty::rc_ptr< ::mosek::fusion::Variable > __mosek_2fusion_2Model__variable(int32_t _1489_n,monty::rc_ptr< ::mosek::fusion::PSDDomain > _1490_psddom) ;
      virtual monty::rc_ptr< ::mosek::fusion::Variable > __mosek_2fusion_2Model__variable(const std::string &  _1491_name,monty::rc_ptr< ::mosek::fusion::PSDDomain > _1492_psddom) ;
      virtual monty::rc_ptr< ::mosek::fusion::Variable > __mosek_2fusion_2Model__variable(const std::string &  _1493_name,int32_t _1494_n,int32_t _1495_m,monty::rc_ptr< ::mosek::fusion::PSDDomain > _1496_psddom) ;
      virtual monty::rc_ptr< ::mosek::fusion::Variable > __mosek_2fusion_2Model__variable(const std::string &  _1497_name,int32_t _1498_n,monty::rc_ptr< ::mosek::fusion::PSDDomain > _1499_psddom) ;
      virtual monty::rc_ptr< ::mosek::fusion::Variable > __mosek_2fusion_2Model__variable(const std::string &  _1500_name,std::shared_ptr< monty::ndarray< int32_t,1 > > _1501_shp,monty::rc_ptr< ::mosek::fusion::PSDDomain > _1502_psddom) ;
      virtual monty::rc_ptr< ::mosek::fusion::Variable > __mosek_2fusion_2Model__variable(monty::rc_ptr< ::mosek::fusion::ConeDomain > _1503_qdom) ;
      virtual monty::rc_ptr< ::mosek::fusion::RangedVariable > __mosek_2fusion_2Model__variable(monty::rc_ptr< ::mosek::fusion::RangeDomain > _1504_rdom) ;
      virtual monty::rc_ptr< ::mosek::fusion::Variable > __mosek_2fusion_2Model__variable(monty::rc_ptr< ::mosek::fusion::LinearDomain > _1505_ldom) ;
      virtual monty::rc_ptr< ::mosek::fusion::Variable > __mosek_2fusion_2Model__variable(std::shared_ptr< monty::ndarray< int32_t,1 > > _1506_shp,monty::rc_ptr< ::mosek::fusion::ConeDomain > _1507_qdom) ;
      virtual monty::rc_ptr< ::mosek::fusion::RangedVariable > __mosek_2fusion_2Model__variable(std::shared_ptr< monty::ndarray< int32_t,1 > > _1508_shp,monty::rc_ptr< ::mosek::fusion::RangeDomain > _1509_rdom) ;
      virtual monty::rc_ptr< ::mosek::fusion::Variable > __mosek_2fusion_2Model__variable(std::shared_ptr< monty::ndarray< int32_t,1 > > _1510_shp,monty::rc_ptr< ::mosek::fusion::LinearDomain > _1511_ldom) ;
      virtual monty::rc_ptr< ::mosek::fusion::Variable > __mosek_2fusion_2Model__variable(std::shared_ptr< monty::ndarray< int32_t,1 > > _1512_shp) ;
      virtual monty::rc_ptr< ::mosek::fusion::Variable > __mosek_2fusion_2Model__variable(int32_t _1513_size,monty::rc_ptr< ::mosek::fusion::ConeDomain > _1514_qdom) ;
      virtual monty::rc_ptr< ::mosek::fusion::RangedVariable > __mosek_2fusion_2Model__variable(int32_t _1515_size,monty::rc_ptr< ::mosek::fusion::RangeDomain > _1516_rdom) ;
      virtual monty::rc_ptr< ::mosek::fusion::Variable > __mosek_2fusion_2Model__variable(int32_t _1517_size,monty::rc_ptr< ::mosek::fusion::LinearDomain > _1518_ldom) ;
      virtual monty::rc_ptr< ::mosek::fusion::Variable > __mosek_2fusion_2Model__variable(int32_t _1519_size) ;
      virtual monty::rc_ptr< ::mosek::fusion::Variable > __mosek_2fusion_2Model__variable() ;
      virtual monty::rc_ptr< ::mosek::fusion::Variable > __mosek_2fusion_2Model__variable(const std::string &  _1520_name,monty::rc_ptr< ::mosek::fusion::ConeDomain > _1521_qdom) ;
      virtual monty::rc_ptr< ::mosek::fusion::RangedVariable > __mosek_2fusion_2Model__variable(const std::string &  _1522_name,monty::rc_ptr< ::mosek::fusion::RangeDomain > _1523_rdom) ;
      virtual monty::rc_ptr< ::mosek::fusion::Variable > __mosek_2fusion_2Model__variable(const std::string &  _1524_name,monty::rc_ptr< ::mosek::fusion::LinearDomain > _1525_ldom) ;
      virtual monty::rc_ptr< ::mosek::fusion::Variable > __mosek_2fusion_2Model__variable(const std::string &  _1526_name,std::shared_ptr< monty::ndarray< int32_t,1 > > _1527_shp,monty::rc_ptr< ::mosek::fusion::ConeDomain > _1528_qdom) ;
      virtual monty::rc_ptr< ::mosek::fusion::RangedVariable > __mosek_2fusion_2Model__variable(const std::string &  _1529_name,std::shared_ptr< monty::ndarray< int32_t,1 > > _1530_shp,monty::rc_ptr< ::mosek::fusion::RangeDomain > _1531_rdom) ;
      virtual monty::rc_ptr< ::mosek::fusion::Variable > __mosek_2fusion_2Model__variable(const std::string &  _1532_name,std::shared_ptr< monty::ndarray< int32_t,1 > > _1533_shp,monty::rc_ptr< ::mosek::fusion::LinearDomain > _1534_ldom) ;
      virtual monty::rc_ptr< ::mosek::fusion::Variable > __mosek_2fusion_2Model__variable(const std::string &  _1535_name,std::shared_ptr< monty::ndarray< int32_t,1 > > _1536_shp) ;
      virtual monty::rc_ptr< ::mosek::fusion::Variable > __mosek_2fusion_2Model__variable(const std::string &  _1537_name,int32_t _1538_size,monty::rc_ptr< ::mosek::fusion::ConeDomain > _1539_qdom) ;
      virtual monty::rc_ptr< ::mosek::fusion::RangedVariable > __mosek_2fusion_2Model__variable(const std::string &  _1540_name,int32_t _1541_size,monty::rc_ptr< ::mosek::fusion::RangeDomain > _1542_rdom) ;
      virtual monty::rc_ptr< ::mosek::fusion::Variable > __mosek_2fusion_2Model__variable(const std::string &  _1543_name,int32_t _1544_size,monty::rc_ptr< ::mosek::fusion::LinearDomain > _1545_ldom) ;
      virtual monty::rc_ptr< ::mosek::fusion::Variable > __mosek_2fusion_2Model__variable(const std::string &  _1546_name,int32_t _1547_size) ;
      virtual monty::rc_ptr< ::mosek::fusion::Variable > __mosek_2fusion_2Model__variable(const std::string &  _1548_name) ;
      virtual void removeConstraintBlock(int32_t _1549_conid) ;
      virtual void removeVariableBlock(int64_t _1550_varid64) ;
      virtual monty::rc_ptr< ::mosek::fusion::RangedVariable > __mosek_2fusion_2Model__ranged_variable(const std::string &  _1555_name,std::shared_ptr< monty::ndarray< int32_t,1 > > _1556_shp,monty::rc_ptr< ::mosek::fusion::RangeDomain > _1557_dom_) ;
      virtual monty::rc_ptr< ::mosek::fusion::Variable > __mosek_2fusion_2Model__variable_(const std::string &  _1618_name,std::shared_ptr< monty::ndarray< int32_t,1 > > _1619_shp,monty::rc_ptr< ::mosek::fusion::ConeDomain > _1620_dom_) ;
      virtual monty::rc_ptr< ::mosek::fusion::Variable > __mosek_2fusion_2Model__variable_(const std::string &  _1658_name,std::shared_ptr< monty::ndarray< int32_t,1 > > _1659_shp,monty::rc_ptr< ::mosek::fusion::LinearDomain > _1660_dom_) ;
      virtual monty::rc_ptr< ::mosek::fusion::Variable > __mosek_2fusion_2Model__variable_(const std::string &  _1719_name,std::shared_ptr< monty::ndarray< int32_t,1 > > _1720_shp,monty::rc_ptr< ::mosek::fusion::PSDDomain > _1721_dom_) ;
      virtual void putfrows(std::shared_ptr< monty::ndarray< int32_t,1 > > _1750_nativeidxs,int32_t _1751_nativebaseptr,monty::rc_ptr< ::mosek::fusion::WorkStack > _1752_rs,int32_t _1753_nelem,int32_t _1754_nnz,int32_t _1755_ptr_base,int32_t _1756_nidxs_base,int32_t _1757_cof_base) ;
      virtual void putarows(std::shared_ptr< monty::ndarray< int32_t,1 > > _1797_nativeidxs,monty::rc_ptr< ::mosek::fusion::WorkStack > _1798_rs,int32_t _1799_nelem,int32_t _1800_nnz,int32_t _1801_ptr_base,int32_t _1802_nidxs_base,int32_t _1803_cof_base) ;
      virtual monty::rc_ptr< ::mosek::fusion::RangedConstraint > __mosek_2fusion_2Model__constraint_(const std::string &  _1840_name,monty::rc_ptr< ::mosek::fusion::Expression > _1841_expr,monty::rc_ptr< ::mosek::fusion::RangeDomain > _1842_dom_) ;
      virtual monty::rc_ptr< ::mosek::fusion::Constraint > __mosek_2fusion_2Model__constraint_(const std::string &  _1885_name,monty::rc_ptr< ::mosek::fusion::Expression > _1886_expr,monty::rc_ptr< ::mosek::fusion::PSDDomain > _1887_dom_) ;
      virtual monty::rc_ptr< ::mosek::fusion::Constraint > __mosek_2fusion_2Model__constraint_(const std::string &  _1981_name,monty::rc_ptr< ::mosek::fusion::Expression > _1982_expr,monty::rc_ptr< ::mosek::fusion::ConeDomain > _1983_dom_) ;
      virtual monty::rc_ptr< ::mosek::fusion::Constraint > __mosek_2fusion_2Model__constraint_(const std::string &  _2040_name,monty::rc_ptr< ::mosek::fusion::Expression > _2041_expr,monty::rc_ptr< ::mosek::fusion::LinearDomain > _2042_dom_) ;
      static  std::string getVersion();
      virtual bool hasParameter(const std::string &  _2081_name) ;
      virtual bool hasConstraint(const std::string &  _2082_name) ;
      virtual bool hasVariable(const std::string &  _2083_name) ;
      virtual monty::rc_ptr< ::mosek::fusion::Parameter > __mosek_2fusion_2Model__getParameter(const std::string &  _2084_name) ;
      virtual monty::rc_ptr< ::mosek::fusion::Constraint > __mosek_2fusion_2Model__getConstraint(int32_t _2085_index) ;
      virtual monty::rc_ptr< ::mosek::fusion::Constraint > __mosek_2fusion_2Model__getConstraint(const std::string &  _2087_name) ;
      virtual monty::rc_ptr< ::mosek::fusion::Variable > __mosek_2fusion_2Model__getVariable(int32_t _2090_index) ;
      virtual monty::rc_ptr< ::mosek::fusion::Variable > __mosek_2fusion_2Model__getVariable(const std::string &  _2091_name) ;
      virtual std::string getName() ;
      virtual std::shared_ptr< monty::ndarray< double,1 > > getParameterValue(std::shared_ptr< monty::ndarray< int32_t,1 > > _2093_idxs) ;
      virtual void setParameterValue(std::shared_ptr< monty::ndarray< int32_t,1 > > _2096_idxs,std::shared_ptr< monty::ndarray< double,1 > > _2097_vals) ;
      virtual monty::rc_ptr< ::mosek::fusion::Model > __mosek_2fusion_2Model__clone() ;
    }; // struct Model;

    // mosek.fusion.Debug from file 'src/fusion/cxx/Debug_p.h'
    // namespace mosek::fusion
    struct p_Debug
    {
      Debug * _pubthis;
    
      p_Debug(Debug * _pubthis) : _pubthis(_pubthis) {}
    
      static Debug::t o ()                 { return monty::rc_ptr<Debug>(new Debug()); }
      Debug::t p (const std::string & val) { std::cout << val; return Debug::t(_pubthis); }
      Debug::t p (      int val)           { std::cout << val; return Debug::t(_pubthis); }
      Debug::t p (int64_t val)           { std::cout << val; return Debug::t(_pubthis); }
      Debug::t p (   double val)           { std::cout << val; return Debug::t(_pubthis); }
      Debug::t p (     bool val)           { std::cout << val; return Debug::t(_pubthis); }
      Debug::t lf()                        { std::cout << std::endl; return Debug::t(_pubthis); }
    
      static p_Debug * _get_impl(Debug * _inst) { return _inst->ptr; }
    
      template<typename T>
      Debug::t p(const std::shared_ptr<monty::ndarray<T,1>> & val)
      {
        if (val->size() > 0)
        {
          std::cout << (*val)[0];
          for (int i = 1; i < val->size(); ++i)
            std::cout << "," << (*val)[i];
        }
        return Debug::t(_pubthis);
      }
    
      Debug::t __mosek_2fusion_2Debug__p (const std::string & val) { std::cout << val; return Debug::t(_pubthis); }
    
      template<class C>
      Debug::t __mosek_2fusion_2Debug__p (C val) { std::cout << val; return Debug::t(_pubthis); }
      Debug::t __mosek_2fusion_2Debug__lf()      { std::cout << std::endl; return Debug::t(_pubthis); }
    
    };
    // End of file 'src/fusion/cxx/Debug_p.h'
    struct p_Sort
    {
      Sort * _pubthis;
      static mosek::fusion::p_Sort* _get_impl(mosek::fusion::Sort * _inst){ assert(_inst); assert(_inst->_impl); return _inst->_impl; }
      static mosek::fusion::p_Sort * _get_impl(mosek::fusion::Sort::t _inst) { return _get_impl(_inst.get()); }
      p_Sort(Sort * _pubthis);
      virtual ~p_Sort() { /* std::cout << "~p_Sort" << std::endl;*/ };

      virtual void destroy();

      static  void argTransposeSort(std::shared_ptr< monty::ndarray< int64_t,1 > > _303_perm,std::shared_ptr< monty::ndarray< int64_t,1 > > _304_ptrb,int32_t _305_m,int32_t _306_n,int32_t _307_p,std::shared_ptr< monty::ndarray< int64_t,1 > > _308_val);
      static  void argsort(std::shared_ptr< monty::ndarray< int64_t,1 > > _316_idx,std::shared_ptr< monty::ndarray< int64_t,1 > > _317_vals1);
      static  void argsort(std::shared_ptr< monty::ndarray< int64_t,1 > > _318_idx,std::shared_ptr< monty::ndarray< int32_t,1 > > _319_vals1);
      static  void argsort(std::shared_ptr< monty::ndarray< int64_t,1 > > _320_idx,std::shared_ptr< monty::ndarray< int64_t,1 > > _321_vals1,std::shared_ptr< monty::ndarray< int64_t,1 > > _322_vals2);
      static  void argsort(std::shared_ptr< monty::ndarray< int64_t,1 > > _323_idx,std::shared_ptr< monty::ndarray< int32_t,1 > > _324_vals1,std::shared_ptr< monty::ndarray< int32_t,1 > > _325_vals2);
      static  void argsort(std::shared_ptr< monty::ndarray< int64_t,1 > > _326_idx,std::shared_ptr< monty::ndarray< int64_t,1 > > _327_vals1,int64_t _328_first,int64_t _329_last);
      static  void argsort(std::shared_ptr< monty::ndarray< int64_t,1 > > _330_idx,std::shared_ptr< monty::ndarray< int32_t,1 > > _331_vals1,int64_t _332_first,int64_t _333_last);
      static  void argsort(std::shared_ptr< monty::ndarray< int64_t,1 > > _334_idx,std::shared_ptr< monty::ndarray< int64_t,1 > > _335_vals1,std::shared_ptr< monty::ndarray< int64_t,1 > > _336_vals2,int64_t _337_first,int64_t _338_last);
      static  void argsort(std::shared_ptr< monty::ndarray< int64_t,1 > > _339_idx,std::shared_ptr< monty::ndarray< int32_t,1 > > _340_vals1,std::shared_ptr< monty::ndarray< int32_t,1 > > _341_vals2,int64_t _342_first,int64_t _343_last);
      static  void argsort(std::shared_ptr< monty::ndarray< int64_t,1 > > _344_idx,std::shared_ptr< monty::ndarray< int64_t,1 > > _345_vals1,int64_t _346_first,int64_t _347_last,bool _348_check);
      static  void argsort(std::shared_ptr< monty::ndarray< int64_t,1 > > _351_idx,std::shared_ptr< monty::ndarray< int32_t,1 > > _352_vals1,int64_t _353_first,int64_t _354_last,bool _355_check);
      static  void argsort(std::shared_ptr< monty::ndarray< int64_t,1 > > _358_idx,std::shared_ptr< monty::ndarray< int64_t,1 > > _359_vals1,std::shared_ptr< monty::ndarray< int64_t,1 > > _360_vals2,int64_t _361_first,int64_t _362_last,bool _363_check);
      static  void argsort(std::shared_ptr< monty::ndarray< int64_t,1 > > _366_idx,std::shared_ptr< monty::ndarray< int32_t,1 > > _367_vals1,std::shared_ptr< monty::ndarray< int32_t,1 > > _368_vals2,int64_t _369_first,int64_t _370_last,bool _371_check);
      static  void argbucketsort(std::shared_ptr< monty::ndarray< int64_t,1 > > _374_idx,std::shared_ptr< monty::ndarray< int64_t,1 > > _375_vals,int64_t _376_first,int64_t _377_last,int64_t _378_minv,int64_t _379_maxv);
      static  void argbucketsort(std::shared_ptr< monty::ndarray< int64_t,1 > > _380_idx,std::shared_ptr< monty::ndarray< int32_t,1 > > _381_vals,int64_t _382_first,int64_t _383_last,int32_t _384_minv,int32_t _385_maxv);
      static  void getminmax(std::shared_ptr< monty::ndarray< int64_t,1 > > _386_idx,std::shared_ptr< monty::ndarray< int64_t,1 > > _387_vals1,std::shared_ptr< monty::ndarray< int64_t,1 > > _388_vals2,int64_t _389_first,int64_t _390_last,std::shared_ptr< monty::ndarray< int64_t,1 > > _391_res);
      static  void getminmax(std::shared_ptr< monty::ndarray< int64_t,1 > > _394_idx,std::shared_ptr< monty::ndarray< int32_t,1 > > _395_vals1,std::shared_ptr< monty::ndarray< int32_t,1 > > _396_vals2,int64_t _397_first,int64_t _398_last,std::shared_ptr< monty::ndarray< int32_t,1 > > _399_res);
      static  bool issorted(std::shared_ptr< monty::ndarray< int64_t,1 > > _402_idx,std::shared_ptr< monty::ndarray< int64_t,1 > > _403_vals1,int64_t _404_first,int64_t _405_last,bool _406_check);
      static  bool issorted(std::shared_ptr< monty::ndarray< int64_t,1 > > _408_idx,std::shared_ptr< monty::ndarray< int32_t,1 > > _409_vals1,int64_t _410_first,int64_t _411_last,bool _412_check);
      static  bool issorted(std::shared_ptr< monty::ndarray< int64_t,1 > > _414_idx,std::shared_ptr< monty::ndarray< int64_t,1 > > _415_vals1,std::shared_ptr< monty::ndarray< int64_t,1 > > _416_vals2,int64_t _417_first,int64_t _418_last,bool _419_check);
      static  bool issorted(std::shared_ptr< monty::ndarray< int64_t,1 > > _421_idx,std::shared_ptr< monty::ndarray< int32_t,1 > > _422_vals1,std::shared_ptr< monty::ndarray< int32_t,1 > > _423_vals2,int64_t _424_first,int64_t _425_last,bool _426_check);
    }; // struct Sort;

    struct p_IndexCounter
    {
      IndexCounter * _pubthis;
      static mosek::fusion::p_IndexCounter* _get_impl(mosek::fusion::IndexCounter * _inst){ assert(_inst); assert(_inst->_impl); return _inst->_impl; }
      static mosek::fusion::p_IndexCounter * _get_impl(mosek::fusion::IndexCounter::t _inst) { return _get_impl(_inst.get()); }
      p_IndexCounter(IndexCounter * _pubthis);
      virtual ~p_IndexCounter() { /* std::cout << "~p_IndexCounter" << std::endl;*/ };
      int64_t start{};
      std::shared_ptr< monty::ndarray< int32_t,1 > > dims{};
      std::shared_ptr< monty::ndarray< int64_t,1 > > strides{};
      std::shared_ptr< monty::ndarray< int64_t,1 > > st{};
      std::shared_ptr< monty::ndarray< int32_t,1 > > ii{};
      int32_t n{};

      virtual void destroy();

      static IndexCounter::t _new_IndexCounter(std::shared_ptr< monty::ndarray< int32_t,1 > > _428_shape);
      void _initialize(std::shared_ptr< monty::ndarray< int32_t,1 > > _428_shape);
      static IndexCounter::t _new_IndexCounter(int64_t _430_start_,std::shared_ptr< monty::ndarray< int32_t,1 > > _431_dims_,std::shared_ptr< monty::ndarray< int32_t,1 > > _432_shape);
      void _initialize(int64_t _430_start_,std::shared_ptr< monty::ndarray< int32_t,1 > > _431_dims_,std::shared_ptr< monty::ndarray< int32_t,1 > > _432_shape);
      static IndexCounter::t _new_IndexCounter(int64_t _435_start_,std::shared_ptr< monty::ndarray< int32_t,1 > > _436_dims_,std::shared_ptr< monty::ndarray< int64_t,1 > > _437_strides_);
      void _initialize(int64_t _435_start_,std::shared_ptr< monty::ndarray< int32_t,1 > > _436_dims_,std::shared_ptr< monty::ndarray< int64_t,1 > > _437_strides_);
      virtual bool atEnd() ;
      virtual std::shared_ptr< monty::ndarray< int32_t,1 > > getIndex() ;
      virtual int64_t next() ;
      virtual int64_t get() ;
      virtual void inc() ;
      virtual void reset() ;
    }; // struct IndexCounter;

    struct p_CommonTools
    {
      CommonTools * _pubthis;
      static mosek::fusion::p_CommonTools* _get_impl(mosek::fusion::CommonTools * _inst){ assert(_inst); assert(_inst->_impl); return _inst->_impl; }
      static mosek::fusion::p_CommonTools * _get_impl(mosek::fusion::CommonTools::t _inst) { return _get_impl(_inst.get()); }
      p_CommonTools(CommonTools * _pubthis);
      virtual ~p_CommonTools() { /* std::cout << "~p_CommonTools" << std::endl;*/ };

      virtual void destroy();

      static  std::shared_ptr< monty::ndarray< int64_t,1 > > resize(std::shared_ptr< monty::ndarray< int64_t,1 > > _443_values,int32_t _444_newsize);
      static  std::shared_ptr< monty::ndarray< int32_t,1 > > resize(std::shared_ptr< monty::ndarray< int32_t,1 > > _446_values,int32_t _447_newsize);
      static  std::shared_ptr< monty::ndarray< double,1 > > resize(std::shared_ptr< monty::ndarray< double,1 > > _449_values,int32_t _450_newsize);
      static  int32_t binarySearch(std::shared_ptr< monty::ndarray< int32_t,1 > > _452_values,int32_t _453_target);
      static  int32_t binarySearch(std::shared_ptr< monty::ndarray< int64_t,1 > > _457_values,int64_t _458_target);
      static  int32_t binarySearchR(std::shared_ptr< monty::ndarray< int64_t,1 > > _460_values,int64_t _461_target);
      static  int32_t binarySearchL(std::shared_ptr< monty::ndarray< int64_t,1 > > _465_values,int64_t _466_target);
      static  void ndIncr(std::shared_ptr< monty::ndarray< int32_t,1 > > _470_ndidx,std::shared_ptr< monty::ndarray< int32_t,1 > > _471_first,std::shared_ptr< monty::ndarray< int32_t,1 > > _472_last);
      static  void transposeTriplets(std::shared_ptr< monty::ndarray< int32_t,1 > > _474_subi,std::shared_ptr< monty::ndarray< int32_t,1 > > _475_subj,std::shared_ptr< monty::ndarray< double,1 > > _476_val,std::shared_ptr< monty::ndarray< std::shared_ptr< monty::ndarray< int64_t,1 > >,1 > > _477_tsubi_,std::shared_ptr< monty::ndarray< std::shared_ptr< monty::ndarray< int64_t,1 > >,1 > > _478_tsubj_,std::shared_ptr< monty::ndarray< std::shared_ptr< monty::ndarray< double,1 > >,1 > > _479_tval_,int64_t _480_nelm,int32_t _481_dimi,int32_t _482_dimj);
      static  void transposeTriplets(std::shared_ptr< monty::ndarray< int32_t,1 > > _495_subi,std::shared_ptr< monty::ndarray< int32_t,1 > > _496_subj,std::shared_ptr< monty::ndarray< double,1 > > _497_val,std::shared_ptr< monty::ndarray< std::shared_ptr< monty::ndarray< int32_t,1 > >,1 > > _498_tsubi_,std::shared_ptr< monty::ndarray< std::shared_ptr< monty::ndarray< int32_t,1 > >,1 > > _499_tsubj_,std::shared_ptr< monty::ndarray< std::shared_ptr< monty::ndarray< double,1 > >,1 > > _500_tval_,int64_t _501_nelm,int32_t _502_dimi,int32_t _503_dimj);
      static  void tripletSort(std::shared_ptr< monty::ndarray< int32_t,1 > > _516_subi,std::shared_ptr< monty::ndarray< int32_t,1 > > _517_subj,std::shared_ptr< monty::ndarray< double,1 > > _518_val,std::shared_ptr< monty::ndarray< std::shared_ptr< monty::ndarray< int32_t,1 > >,1 > > _519_tsubi_,std::shared_ptr< monty::ndarray< std::shared_ptr< monty::ndarray< int32_t,1 > >,1 > > _520_tsubj_,std::shared_ptr< monty::ndarray< std::shared_ptr< monty::ndarray< double,1 > >,1 > > _521_tval_,int64_t _522_nelm,int32_t _523_dimi,int32_t _524_dimj);
      static  void argMSort(std::shared_ptr< monty::ndarray< int32_t,1 > > _550_idx,std::shared_ptr< monty::ndarray< int32_t,1 > > _551_vals);
      static  void mergeInto(std::shared_ptr< monty::ndarray< int32_t,1 > > _556_src,std::shared_ptr< monty::ndarray< int32_t,1 > > _557_tgt,std::shared_ptr< monty::ndarray< int32_t,1 > > _558_vals,int32_t _559_si0,int32_t _560_si1_,int32_t _561_si2_);
      static  void argQsort(std::shared_ptr< monty::ndarray< int64_t,1 > > _567_idx,std::shared_ptr< monty::ndarray< int64_t,1 > > _568_vals1,std::shared_ptr< monty::ndarray< int64_t,1 > > _569_vals2,int64_t _570_first,int64_t _571_last);
      static  void argQsort(std::shared_ptr< monty::ndarray< int64_t,1 > > _572_idx,std::shared_ptr< monty::ndarray< int32_t,1 > > _573_vals1,std::shared_ptr< monty::ndarray< int32_t,1 > > _574_vals2,int64_t _575_first,int64_t _576_last);
    }; // struct CommonTools;

    struct p_SolutionStruct
    {
      SolutionStruct * _pubthis;
      static mosek::fusion::p_SolutionStruct* _get_impl(mosek::fusion::SolutionStruct * _inst){ assert(_inst); assert(_inst->_impl); return _inst->_impl; }
      static mosek::fusion::p_SolutionStruct * _get_impl(mosek::fusion::SolutionStruct::t _inst) { return _get_impl(_inst.get()); }
      p_SolutionStruct(SolutionStruct * _pubthis);
      virtual ~p_SolutionStruct() { /* std::cout << "~p_SolutionStruct" << std::endl;*/ };
      std::shared_ptr< monty::ndarray< double,1 > > accy{};
      std::shared_ptr< monty::ndarray< double,1 > > accx{};
      std::shared_ptr< monty::ndarray< int32_t,1 > > accptr{};
      std::shared_ptr< monty::ndarray< double,1 > > yx{};
      std::shared_ptr< monty::ndarray< double,1 > > sux{};
      std::shared_ptr< monty::ndarray< double,1 > > slx{};
      std::shared_ptr< monty::ndarray< double,1 > > bars{};
      std::shared_ptr< monty::ndarray< double,1 > > barx{};
      std::shared_ptr< monty::ndarray< double,1 > > y{};
      std::shared_ptr< monty::ndarray< double,1 > > suc{};
      std::shared_ptr< monty::ndarray< double,1 > > slc{};
      std::shared_ptr< monty::ndarray< double,1 > > xx{};
      std::shared_ptr< monty::ndarray< double,1 > > xc{};
      double dobj{};
      double pobj{};
      mosek::fusion::ProblemStatus probstatus{};
      mosek::fusion::SolutionStatus dstatus{};
      mosek::fusion::SolutionStatus pstatus{};
      int32_t sol_numbarvar{};
      int32_t sol_numaccelm{};
      int32_t sol_numacc{};
      int32_t sol_numvar{};
      int32_t sol_numcon{};

      virtual void destroy();

      static SolutionStruct::t _new_SolutionStruct(int32_t _577_numvar,int32_t _578_numcon,int32_t _579_numbarvar,int32_t _580_numacc,int32_t _581_numaccelm);
      void _initialize(int32_t _577_numvar,int32_t _578_numcon,int32_t _579_numbarvar,int32_t _580_numacc,int32_t _581_numaccelm);
      static SolutionStruct::t _new_SolutionStruct(monty::rc_ptr< ::mosek::fusion::SolutionStruct > _582_that);
      void _initialize(monty::rc_ptr< ::mosek::fusion::SolutionStruct > _582_that);
      virtual monty::rc_ptr< ::mosek::fusion::SolutionStruct > __mosek_2fusion_2SolutionStruct__clone() ;
      virtual void resize(int32_t _583_numvar,int32_t _584_numcon,int32_t _585_numbarvar,int32_t _586_numacc,int32_t _587_numaccelm) ;
      virtual bool isDualAcceptable(mosek::fusion::AccSolutionStatus _608_acceptable_sol) ;
      virtual bool isPrimalAcceptable(mosek::fusion::AccSolutionStatus _609_acceptable_sol) ;
      virtual bool isAcceptable(mosek::fusion::SolutionStatus _610_stat,mosek::fusion::AccSolutionStatus _611_accstat) ;
    }; // struct SolutionStruct;

    struct p_RowBlockManager
    {
      RowBlockManager * _pubthis;
      static mosek::fusion::p_RowBlockManager* _get_impl(mosek::fusion::RowBlockManager * _inst){ assert(_inst); assert(_inst->_impl); return _inst->_impl; }
      static mosek::fusion::p_RowBlockManager * _get_impl(mosek::fusion::RowBlockManager::t _inst) { return _get_impl(_inst.get()); }
      p_RowBlockManager(RowBlockManager * _pubthis);
      virtual ~p_RowBlockManager() { /* std::cout << "~p_RowBlockManager" << std::endl;*/ };
      int32_t varidx_used{};
      int32_t code_used{};
      std::shared_ptr< monty::ndarray< double,1 > > cconst{};
      std::shared_ptr< monty::ndarray< int32_t,1 > > code{};
      int32_t first_free_codeitem{};
      std::shared_ptr< monty::ndarray< int32_t,1 > > param_code_sizes{};
      std::shared_ptr< monty::ndarray< int64_t,1 > > param_varidx{};
      int32_t first_free_entry{};
      std::shared_ptr< monty::ndarray< int32_t,1 > > row_code_ptr{};
      std::shared_ptr< monty::ndarray< int32_t,1 > > row_param_ptre{};
      std::shared_ptr< monty::ndarray< int32_t,1 > > row_param_ptrb{};
      monty::rc_ptr< ::mosek::fusion::LinkedBlocks > blocks{};

      virtual void destroy();

      static RowBlockManager::t _new_RowBlockManager(monty::rc_ptr< ::mosek::fusion::RowBlockManager > _612_that);
      void _initialize(monty::rc_ptr< ::mosek::fusion::RowBlockManager > _612_that);
      static RowBlockManager::t _new_RowBlockManager();
      void _initialize();
      virtual int32_t num_parameterized() ;
      virtual bool is_parameterized() ;
      virtual int32_t blocksize(int32_t _613_id) ;
      virtual int32_t block_capacity() ;
      virtual int32_t capacity() ;
      virtual void get(int32_t _614_id,std::shared_ptr< monty::ndarray< int32_t,1 > > _615_target,int32_t _616_offset) ;
      virtual void evaluate(monty::rc_ptr< ::mosek::fusion::WorkStack > _617_xs,std::shared_ptr< monty::ndarray< double,1 > > _618_param_value,std::shared_ptr< monty::ndarray< int32_t,1 > > _619_subi,std::shared_ptr< monty::ndarray< int32_t,1 > > _620_subj,std::shared_ptr< monty::ndarray< double,1 > > _621_val) ;
      virtual void replace_row_code(monty::rc_ptr< ::mosek::fusion::WorkStack > _632_rs,std::shared_ptr< monty::ndarray< int32_t,1 > > _633_nativeidxs,int32_t _634_ptr,int32_t _635_nidxs,int32_t _636_codeptr,int32_t _637_code_p,int32_t _638_cconst_p) ;
      virtual void clear_row_code(std::shared_ptr< monty::ndarray< int32_t,1 > > _661_nativeidxs) ;
      virtual void compress() ;
      virtual void ensure_code_cap(int32_t _674_nentry,int32_t _675_codesize) ;
      virtual void release(int32_t _685_id,std::shared_ptr< monty::ndarray< int32_t,1 > > _686_nativeidxs) ;
      virtual int32_t allocate(std::shared_ptr< monty::ndarray< int32_t,1 > > _690_nativeidxs) ;
      virtual bool row_is_parameterized(int32_t _696_i) ;
    }; // struct RowBlockManager;

    struct p_BaseVariable : public /*implements*/ virtual ::mosek::fusion::Variable
    {
      BaseVariable * _pubthis;
      static mosek::fusion::p_BaseVariable* _get_impl(mosek::fusion::BaseVariable * _inst){ assert(_inst); assert(_inst->_impl); return _inst->_impl; }
      static mosek::fusion::p_BaseVariable * _get_impl(mosek::fusion::BaseVariable::t _inst) { return _get_impl(_inst.get()); }
      p_BaseVariable(BaseVariable * _pubthis);
      virtual ~p_BaseVariable() { /* std::cout << "~p_BaseVariable" << std::endl;*/ };
      std::shared_ptr< monty::ndarray< int64_t,1 > > sparsity{};
      std::shared_ptr< monty::ndarray< int64_t,1 > > basevar_nativeidxs{};
      monty::rc_ptr< ::mosek::fusion::Model > model{};
      std::shared_ptr< monty::ndarray< int32_t,1 > > shape{};

      virtual void destroy();

      static BaseVariable::t _new_BaseVariable(monty::rc_ptr< ::mosek::fusion::BaseVariable > _2272_v,monty::rc_ptr< ::mosek::fusion::Model > _2273_m);
      void _initialize(monty::rc_ptr< ::mosek::fusion::BaseVariable > _2272_v,monty::rc_ptr< ::mosek::fusion::Model > _2273_m);
      static BaseVariable::t _new_BaseVariable(monty::rc_ptr< ::mosek::fusion::Model > _2274_m,std::shared_ptr< monty::ndarray< int32_t,1 > > _2275_shape,std::shared_ptr< monty::ndarray< int64_t,1 > > _2276_sparsity,std::shared_ptr< monty::ndarray< int64_t,1 > > _2277_basevar_nativeidxs);
      void _initialize(monty::rc_ptr< ::mosek::fusion::Model > _2274_m,std::shared_ptr< monty::ndarray< int32_t,1 > > _2275_shape,std::shared_ptr< monty::ndarray< int64_t,1 > > _2276_sparsity,std::shared_ptr< monty::ndarray< int64_t,1 > > _2277_basevar_nativeidxs);
      virtual /* override */ std::string toString() ;
      virtual void eval(monty::rc_ptr< ::mosek::fusion::WorkStack > _2280_rs,monty::rc_ptr< ::mosek::fusion::WorkStack > _2281_ws,monty::rc_ptr< ::mosek::fusion::WorkStack > _2282_xs) ;
      virtual void remove() ;
      virtual monty::rc_ptr< ::mosek::fusion::Variable > __mosek_2fusion_2BaseVariable__fromTril(int32_t _2300_dim0,int32_t _2301_d) ;
      virtual monty::rc_ptr< ::mosek::fusion::Variable > __mosek_2fusion_2BaseVariable__fromTril(int32_t _2334_d) ;
      virtual monty::rc_ptr< ::mosek::fusion::Variable > __mosek_2fusion_2Variable__fromTril(int32_t _2334_d) { return __mosek_2fusion_2BaseVariable__fromTril(_2334_d); }
      virtual monty::rc_ptr< ::mosek::fusion::Variable > __mosek_2fusion_2BaseVariable__tril(int32_t _2335_dim1,int32_t _2336_dim2) ;
      virtual monty::rc_ptr< ::mosek::fusion::Variable > __mosek_2fusion_2BaseVariable__tril() ;
      virtual monty::rc_ptr< ::mosek::fusion::Variable > __mosek_2fusion_2Variable__tril() { return __mosek_2fusion_2BaseVariable__tril(); }
      virtual monty::rc_ptr< ::mosek::fusion::Variable > __mosek_2fusion_2BaseVariable__reshape(int32_t _2390_dim0,int32_t _2391_dim1,int32_t _2392_dim2) ;
      virtual monty::rc_ptr< ::mosek::fusion::Variable > __mosek_2fusion_2Variable__reshape(int32_t _2390_dim0,int32_t _2391_dim1,int32_t _2392_dim2) { return __mosek_2fusion_2BaseVariable__reshape(_2390_dim0,_2391_dim1,_2392_dim2); }
      virtual monty::rc_ptr< ::mosek::fusion::Variable > __mosek_2fusion_2BaseVariable__reshape(int32_t _2393_dim0,int32_t _2394_dim1) ;
      virtual monty::rc_ptr< ::mosek::fusion::Variable > __mosek_2fusion_2Variable__reshape(int32_t _2393_dim0,int32_t _2394_dim1) { return __mosek_2fusion_2BaseVariable__reshape(_2393_dim0,_2394_dim1); }
      virtual monty::rc_ptr< ::mosek::fusion::Variable > __mosek_2fusion_2BaseVariable__reshape(int32_t _2395_dim0) ;
      virtual monty::rc_ptr< ::mosek::fusion::Variable > __mosek_2fusion_2Variable__reshape(int32_t _2395_dim0) { return __mosek_2fusion_2BaseVariable__reshape(_2395_dim0); }
      virtual monty::rc_ptr< ::mosek::fusion::Variable > __mosek_2fusion_2BaseVariable__reshape(std::shared_ptr< monty::ndarray< int32_t,1 > > _2396_shape) ;
      virtual monty::rc_ptr< ::mosek::fusion::Variable > __mosek_2fusion_2Variable__reshape(std::shared_ptr< monty::ndarray< int32_t,1 > > _2396_shape) { return __mosek_2fusion_2BaseVariable__reshape(_2396_shape); }
      virtual void setLevel(std::shared_ptr< monty::ndarray< double,1 > > _2400_v) ;
      virtual monty::rc_ptr< ::mosek::fusion::Model > __mosek_2fusion_2BaseVariable__getModel() ;
      virtual monty::rc_ptr< ::mosek::fusion::Model > __mosek_2fusion_2Variable__getModel() { return __mosek_2fusion_2BaseVariable__getModel(); }
      virtual int32_t getND() ;
      virtual int32_t getDim(int32_t _2403_d) ;
      virtual std::shared_ptr< monty::ndarray< int32_t,1 > > getShape() ;
      virtual int64_t getSize() ;
      virtual std::shared_ptr< monty::ndarray< double,1 > > dual() ;
      virtual std::shared_ptr< monty::ndarray< double,1 > > level() ;
      virtual void makeContinuous() ;
      virtual void makeInteger() ;
      virtual monty::rc_ptr< ::mosek::fusion::Variable > __mosek_2fusion_2BaseVariable__transpose() ;
      virtual monty::rc_ptr< ::mosek::fusion::Variable > __mosek_2fusion_2Variable__transpose() { return __mosek_2fusion_2BaseVariable__transpose(); }
      virtual monty::rc_ptr< ::mosek::fusion::Variable > __mosek_2fusion_2BaseVariable__index(int32_t _2424_i0,int32_t _2425_i1,int32_t _2426_i2) ;
      virtual monty::rc_ptr< ::mosek::fusion::Variable > __mosek_2fusion_2Variable__index(int32_t _2424_i0,int32_t _2425_i1,int32_t _2426_i2) { return __mosek_2fusion_2BaseVariable__index(_2424_i0,_2425_i1,_2426_i2); }
      virtual monty::rc_ptr< ::mosek::fusion::Variable > __mosek_2fusion_2BaseVariable__index(int32_t _2427_i0,int32_t _2428_i1) ;
      virtual monty::rc_ptr< ::mosek::fusion::Variable > __mosek_2fusion_2Variable__index(int32_t _2427_i0,int32_t _2428_i1) { return __mosek_2fusion_2BaseVariable__index(_2427_i0,_2428_i1); }
      virtual monty::rc_ptr< ::mosek::fusion::Variable > __mosek_2fusion_2BaseVariable__index(std::shared_ptr< monty::ndarray< int32_t,1 > > _2429_index) ;
      virtual monty::rc_ptr< ::mosek::fusion::Variable > __mosek_2fusion_2Variable__index(std::shared_ptr< monty::ndarray< int32_t,1 > > _2429_index) { return __mosek_2fusion_2BaseVariable__index(_2429_index); }
      virtual monty::rc_ptr< ::mosek::fusion::Variable > __mosek_2fusion_2BaseVariable__index(int32_t _2432_index) ;
      virtual monty::rc_ptr< ::mosek::fusion::Variable > __mosek_2fusion_2Variable__index(int32_t _2432_index) { return __mosek_2fusion_2BaseVariable__index(_2432_index); }
      virtual monty::rc_ptr< ::mosek::fusion::Variable > __mosek_2fusion_2BaseVariable__pick(std::shared_ptr< monty::ndarray< int32_t,1 > > _2433_i0,std::shared_ptr< monty::ndarray< int32_t,1 > > _2434_i1,std::shared_ptr< monty::ndarray< int32_t,1 > > _2435_i2) ;
      virtual monty::rc_ptr< ::mosek::fusion::Variable > __mosek_2fusion_2Variable__pick(std::shared_ptr< monty::ndarray< int32_t,1 > > _2433_i0,std::shared_ptr< monty::ndarray< int32_t,1 > > _2434_i1,std::shared_ptr< monty::ndarray< int32_t,1 > > _2435_i2) { return __mosek_2fusion_2BaseVariable__pick(_2433_i0,_2434_i1,_2435_i2); }
      virtual monty::rc_ptr< ::mosek::fusion::Variable > __mosek_2fusion_2BaseVariable__pick(std::shared_ptr< monty::ndarray< int32_t,1 > > _2438_i0,std::shared_ptr< monty::ndarray< int32_t,1 > > _2439_i1) ;
      virtual monty::rc_ptr< ::mosek::fusion::Variable > __mosek_2fusion_2Variable__pick(std::shared_ptr< monty::ndarray< int32_t,1 > > _2438_i0,std::shared_ptr< monty::ndarray< int32_t,1 > > _2439_i1) { return __mosek_2fusion_2BaseVariable__pick(_2438_i0,_2439_i1); }
      virtual monty::rc_ptr< ::mosek::fusion::Variable > __mosek_2fusion_2BaseVariable__pick(std::shared_ptr< monty::ndarray< int32_t,2 > > _2442_midxs) ;
      virtual monty::rc_ptr< ::mosek::fusion::Variable > __mosek_2fusion_2Variable__pick(std::shared_ptr< monty::ndarray< int32_t,2 > > _2442_midxs) { return __mosek_2fusion_2BaseVariable__pick(_2442_midxs); }
      virtual monty::rc_ptr< ::mosek::fusion::Variable > __mosek_2fusion_2BaseVariable__pick(std::shared_ptr< monty::ndarray< int32_t,1 > > _2464_idxs) ;
      virtual monty::rc_ptr< ::mosek::fusion::Variable > __mosek_2fusion_2Variable__pick(std::shared_ptr< monty::ndarray< int32_t,1 > > _2464_idxs) { return __mosek_2fusion_2BaseVariable__pick(_2464_idxs); }
      virtual monty::rc_ptr< ::mosek::fusion::Variable > __mosek_2fusion_2BaseVariable__antidiag(int32_t _2475_index) ;
      virtual monty::rc_ptr< ::mosek::fusion::Variable > __mosek_2fusion_2Variable__antidiag(int32_t _2475_index) { return __mosek_2fusion_2BaseVariable__antidiag(_2475_index); }
      virtual monty::rc_ptr< ::mosek::fusion::Variable > __mosek_2fusion_2BaseVariable__antidiag() ;
      virtual monty::rc_ptr< ::mosek::fusion::Variable > __mosek_2fusion_2Variable__antidiag() { return __mosek_2fusion_2BaseVariable__antidiag(); }
      virtual monty::rc_ptr< ::mosek::fusion::Variable > __mosek_2fusion_2BaseVariable__diag(int32_t _2476_index) ;
      virtual monty::rc_ptr< ::mosek::fusion::Variable > __mosek_2fusion_2Variable__diag(int32_t _2476_index) { return __mosek_2fusion_2BaseVariable__diag(_2476_index); }
      virtual monty::rc_ptr< ::mosek::fusion::Variable > __mosek_2fusion_2BaseVariable__diag() ;
      virtual monty::rc_ptr< ::mosek::fusion::Variable > __mosek_2fusion_2Variable__diag() { return __mosek_2fusion_2BaseVariable__diag(); }
      virtual monty::rc_ptr< ::mosek::fusion::Variable > __mosek_2fusion_2BaseVariable__general_diag(std::shared_ptr< monty::ndarray< int32_t,1 > > _2477_firstidx,std::shared_ptr< monty::ndarray< int32_t,1 > > _2478_step,int32_t _2479_num) ;
      virtual monty::rc_ptr< ::mosek::fusion::Variable > __mosek_2fusion_2BaseVariable__slice(std::shared_ptr< monty::ndarray< int32_t,1 > > _2500_first,std::shared_ptr< monty::ndarray< int32_t,1 > > _2501_last) ;
      virtual monty::rc_ptr< ::mosek::fusion::Variable > __mosek_2fusion_2Variable__slice(std::shared_ptr< monty::ndarray< int32_t,1 > > _2500_first,std::shared_ptr< monty::ndarray< int32_t,1 > > _2501_last) { return __mosek_2fusion_2BaseVariable__slice(_2500_first,_2501_last); }
      virtual monty::rc_ptr< ::mosek::fusion::Variable > __mosek_2fusion_2BaseVariable__slice(int32_t _2535_first,int32_t _2536_last) ;
      virtual monty::rc_ptr< ::mosek::fusion::Variable > __mosek_2fusion_2Variable__slice(int32_t _2535_first,int32_t _2536_last) { return __mosek_2fusion_2BaseVariable__slice(_2535_first,_2536_last); }
      virtual monty::rc_ptr< ::mosek::fusion::Expression > __mosek_2fusion_2BaseVariable__asExpr() ;
      virtual monty::rc_ptr< ::mosek::fusion::Expression > __mosek_2fusion_2Variable__asExpr() { return __mosek_2fusion_2BaseVariable__asExpr(); }
      virtual int32_t inst(int32_t _2545_spoffset,std::shared_ptr< monty::ndarray< int64_t,1 > > _2546_sparsity,int32_t _2547_nioffset,std::shared_ptr< monty::ndarray< int64_t,1 > > _2548_basevar_nativeidxs) ;
      virtual int32_t numInst() ;
      virtual void inst(int32_t _2553_offset,std::shared_ptr< monty::ndarray< int64_t,1 > > _2554_nindex) ;
      virtual void set_values(std::shared_ptr< monty::ndarray< double,1 > > _2561_values,bool _2562_primal) ;
      virtual void dual_lu(int32_t _2567_offset,std::shared_ptr< monty::ndarray< double,1 > > _2568_target,bool _2569_lower) ;
      virtual void values(int32_t _2572_offset,std::shared_ptr< monty::ndarray< double,1 > > _2573_target,bool _2574_primal) ;
      virtual void make_continuous() ;
      virtual void make_integer() ;
    }; // struct BaseVariable;

    struct p_SliceVariable : public ::mosek::fusion::p_BaseVariable
    {
      SliceVariable * _pubthis;
      static mosek::fusion::p_SliceVariable* _get_impl(mosek::fusion::SliceVariable * _inst){ return static_cast< mosek::fusion::p_SliceVariable* >(mosek::fusion::p_BaseVariable::_get_impl(_inst)); }
      static mosek::fusion::p_SliceVariable * _get_impl(mosek::fusion::SliceVariable::t _inst) { return _get_impl(_inst.get()); }
      p_SliceVariable(SliceVariable * _pubthis);
      virtual ~p_SliceVariable() { /* std::cout << "~p_SliceVariable" << std::endl;*/ };
      std::shared_ptr< monty::ndarray< int32_t,1 > > shape{};
      std::shared_ptr< monty::ndarray< int64_t,1 > > sparsity{};
      std::shared_ptr< monty::ndarray< int64_t,1 > > nativeidxs{};

      virtual void destroy();

      static SliceVariable::t _new_SliceVariable(monty::rc_ptr< ::mosek::fusion::Model > _2125_m,std::shared_ptr< monty::ndarray< int32_t,1 > > _2126_shape,std::shared_ptr< monty::ndarray< int64_t,1 > > _2127_sparsity,std::shared_ptr< monty::ndarray< int64_t,1 > > _2128_nativeidxs);
      void _initialize(monty::rc_ptr< ::mosek::fusion::Model > _2125_m,std::shared_ptr< monty::ndarray< int32_t,1 > > _2126_shape,std::shared_ptr< monty::ndarray< int64_t,1 > > _2127_sparsity,std::shared_ptr< monty::ndarray< int64_t,1 > > _2128_nativeidxs);
      static SliceVariable::t _new_SliceVariable(monty::rc_ptr< ::mosek::fusion::SliceVariable > _2129_v);
      void _initialize(monty::rc_ptr< ::mosek::fusion::SliceVariable > _2129_v);
    }; // struct SliceVariable;

    struct p_BoundInterfaceVariable : public ::mosek::fusion::p_SliceVariable
    {
      BoundInterfaceVariable * _pubthis;
      static mosek::fusion::p_BoundInterfaceVariable* _get_impl(mosek::fusion::BoundInterfaceVariable * _inst){ return static_cast< mosek::fusion::p_BoundInterfaceVariable* >(mosek::fusion::p_SliceVariable::_get_impl(_inst)); }
      static mosek::fusion::p_BoundInterfaceVariable * _get_impl(mosek::fusion::BoundInterfaceVariable::t _inst) { return _get_impl(_inst.get()); }
      p_BoundInterfaceVariable(BoundInterfaceVariable * _pubthis);
      virtual ~p_BoundInterfaceVariable() { /* std::cout << "~p_BoundInterfaceVariable" << std::endl;*/ };
      bool islower{};

      virtual void destroy();

      static BoundInterfaceVariable::t _new_BoundInterfaceVariable(monty::rc_ptr< ::mosek::fusion::Model > _2099_m,std::shared_ptr< monty::ndarray< int32_t,1 > > _2100_shape,std::shared_ptr< monty::ndarray< int64_t,1 > > _2101_sparsity,std::shared_ptr< monty::ndarray< int64_t,1 > > _2102_nativeidxs,bool _2103_islower);
      void _initialize(monty::rc_ptr< ::mosek::fusion::Model > _2099_m,std::shared_ptr< monty::ndarray< int32_t,1 > > _2100_shape,std::shared_ptr< monty::ndarray< int64_t,1 > > _2101_sparsity,std::shared_ptr< monty::ndarray< int64_t,1 > > _2102_nativeidxs,bool _2103_islower);
      static BoundInterfaceVariable::t _new_BoundInterfaceVariable(monty::rc_ptr< ::mosek::fusion::SliceVariable > _2104_v,bool _2105_islower);
      void _initialize(monty::rc_ptr< ::mosek::fusion::SliceVariable > _2104_v,bool _2105_islower);
      virtual /* override */ std::shared_ptr< monty::ndarray< double,1 > > dual() ;
      virtual /* override */ monty::rc_ptr< ::mosek::fusion::Variable > __mosek_2fusion_2BoundInterfaceVariable__transpose() ;
      virtual monty::rc_ptr< ::mosek::fusion::Variable > __mosek_2fusion_2BaseVariable__transpose() { return __mosek_2fusion_2BoundInterfaceVariable__transpose(); }
      virtual /* override */ monty::rc_ptr< ::mosek::fusion::Variable > __mosek_2fusion_2BoundInterfaceVariable__pick(std::shared_ptr< monty::ndarray< int32_t,1 > > _2107_i0,std::shared_ptr< monty::ndarray< int32_t,1 > > _2108_i1,std::shared_ptr< monty::ndarray< int32_t,1 > > _2109_i2) ;
      virtual monty::rc_ptr< ::mosek::fusion::Variable > __mosek_2fusion_2BaseVariable__pick(std::shared_ptr< monty::ndarray< int32_t,1 > > _2107_i0,std::shared_ptr< monty::ndarray< int32_t,1 > > _2108_i1,std::shared_ptr< monty::ndarray< int32_t,1 > > _2109_i2) { return __mosek_2fusion_2BoundInterfaceVariable__pick(_2107_i0,_2108_i1,_2109_i2); }
      virtual /* override */ monty::rc_ptr< ::mosek::fusion::Variable > __mosek_2fusion_2BoundInterfaceVariable__pick(std::shared_ptr< monty::ndarray< int32_t,1 > > _2110_i0,std::shared_ptr< monty::ndarray< int32_t,1 > > _2111_i1) ;
      virtual monty::rc_ptr< ::mosek::fusion::Variable > __mosek_2fusion_2BaseVariable__pick(std::shared_ptr< monty::ndarray< int32_t,1 > > _2110_i0,std::shared_ptr< monty::ndarray< int32_t,1 > > _2111_i1) { return __mosek_2fusion_2BoundInterfaceVariable__pick(_2110_i0,_2111_i1); }
      virtual /* override */ monty::rc_ptr< ::mosek::fusion::Variable > __mosek_2fusion_2BoundInterfaceVariable__pick(std::shared_ptr< monty::ndarray< int32_t,2 > > _2112_midxs) ;
      virtual monty::rc_ptr< ::mosek::fusion::Variable > __mosek_2fusion_2BaseVariable__pick(std::shared_ptr< monty::ndarray< int32_t,2 > > _2112_midxs) { return __mosek_2fusion_2BoundInterfaceVariable__pick(_2112_midxs); }
      virtual /* override */ monty::rc_ptr< ::mosek::fusion::Variable > __mosek_2fusion_2BoundInterfaceVariable__pick(std::shared_ptr< monty::ndarray< int32_t,1 > > _2113_idxs) ;
      virtual monty::rc_ptr< ::mosek::fusion::Variable > __mosek_2fusion_2BaseVariable__pick(std::shared_ptr< monty::ndarray< int32_t,1 > > _2113_idxs) { return __mosek_2fusion_2BoundInterfaceVariable__pick(_2113_idxs); }
      virtual /* override */ monty::rc_ptr< ::mosek::fusion::Variable > __mosek_2fusion_2BoundInterfaceVariable__antidiag(int32_t _2114_index) ;
      virtual monty::rc_ptr< ::mosek::fusion::Variable > __mosek_2fusion_2BaseVariable__antidiag(int32_t _2114_index) { return __mosek_2fusion_2BoundInterfaceVariable__antidiag(_2114_index); }
      virtual /* override */ monty::rc_ptr< ::mosek::fusion::Variable > __mosek_2fusion_2BoundInterfaceVariable__antidiag() ;
      virtual monty::rc_ptr< ::mosek::fusion::Variable > __mosek_2fusion_2BaseVariable__antidiag() { return __mosek_2fusion_2BoundInterfaceVariable__antidiag(); }
      virtual /* override */ monty::rc_ptr< ::mosek::fusion::Variable > __mosek_2fusion_2BoundInterfaceVariable__diag(int32_t _2115_index) ;
      virtual monty::rc_ptr< ::mosek::fusion::Variable > __mosek_2fusion_2BaseVariable__diag(int32_t _2115_index) { return __mosek_2fusion_2BoundInterfaceVariable__diag(_2115_index); }
      virtual /* override */ monty::rc_ptr< ::mosek::fusion::Variable > __mosek_2fusion_2BoundInterfaceVariable__diag() ;
      virtual monty::rc_ptr< ::mosek::fusion::Variable > __mosek_2fusion_2BaseVariable__diag() { return __mosek_2fusion_2BoundInterfaceVariable__diag(); }
      virtual /* override */ monty::rc_ptr< ::mosek::fusion::Variable > __mosek_2fusion_2BoundInterfaceVariable__slice(std::shared_ptr< monty::ndarray< int32_t,1 > > _2116_firsta,std::shared_ptr< monty::ndarray< int32_t,1 > > _2117_lasta) ;
      virtual monty::rc_ptr< ::mosek::fusion::Variable > __mosek_2fusion_2BaseVariable__slice(std::shared_ptr< monty::ndarray< int32_t,1 > > _2116_firsta,std::shared_ptr< monty::ndarray< int32_t,1 > > _2117_lasta) { return __mosek_2fusion_2BoundInterfaceVariable__slice(_2116_firsta,_2117_lasta); }
      virtual /* override */ monty::rc_ptr< ::mosek::fusion::Variable > __mosek_2fusion_2BoundInterfaceVariable__slice(int32_t _2118_first,int32_t _2119_last) ;
      virtual monty::rc_ptr< ::mosek::fusion::Variable > __mosek_2fusion_2BaseVariable__slice(int32_t _2118_first,int32_t _2119_last) { return __mosek_2fusion_2BoundInterfaceVariable__slice(_2118_first,_2119_last); }
      virtual monty::rc_ptr< ::mosek::fusion::BoundInterfaceVariable > __mosek_2fusion_2BoundInterfaceVariable__from_(monty::rc_ptr< ::mosek::fusion::Variable > _2120_v) ;
    }; // struct BoundInterfaceVariable;

    struct p_ModelVariable : public ::mosek::fusion::p_BaseVariable
    {
      ModelVariable * _pubthis;
      static mosek::fusion::p_ModelVariable* _get_impl(mosek::fusion::ModelVariable * _inst){ return static_cast< mosek::fusion::p_ModelVariable* >(mosek::fusion::p_BaseVariable::_get_impl(_inst)); }
      static mosek::fusion::p_ModelVariable * _get_impl(mosek::fusion::ModelVariable::t _inst) { return _get_impl(_inst.get()); }
      p_ModelVariable(ModelVariable * _pubthis);
      virtual ~p_ModelVariable() { /* std::cout << "~p_ModelVariable" << std::endl;*/ };
      std::shared_ptr< monty::ndarray< int64_t,1 > > sparsity{};
      std::shared_ptr< monty::ndarray< int32_t,1 > > shape{};
      std::shared_ptr< monty::ndarray< int64_t,1 > > modelvar_nativeidxs{};
      int64_t varid{};
      std::string name{};

      virtual void destroy();

      static ModelVariable::t _new_ModelVariable(monty::rc_ptr< ::mosek::fusion::ModelVariable > _2235_v,monty::rc_ptr< ::mosek::fusion::Model > _2236_m);
      void _initialize(monty::rc_ptr< ::mosek::fusion::ModelVariable > _2235_v,monty::rc_ptr< ::mosek::fusion::Model > _2236_m);
      static ModelVariable::t _new_ModelVariable(monty::rc_ptr< ::mosek::fusion::Model > _2237_model,const std::string &  _2238_name,std::shared_ptr< monty::ndarray< int32_t,1 > > _2239_shape,int64_t _2240_varid,std::shared_ptr< monty::ndarray< int64_t,1 > > _2241_sparsity,std::shared_ptr< monty::ndarray< int64_t,1 > > _2242_modelvar_nativeidxs);
      void _initialize(monty::rc_ptr< ::mosek::fusion::Model > _2237_model,const std::string &  _2238_name,std::shared_ptr< monty::ndarray< int32_t,1 > > _2239_shape,int64_t _2240_varid,std::shared_ptr< monty::ndarray< int64_t,1 > > _2241_sparsity,std::shared_ptr< monty::ndarray< int64_t,1 > > _2242_modelvar_nativeidxs);
      virtual void flushNames() { throw monty::AbstractClassError("Call to abstract method"); }
      virtual void elementName(int64_t _2243_index,monty::rc_ptr< ::mosek::fusion::Utils::StringBuffer > _2244_sb) ;
      virtual /* override */ void remove() ;
      virtual monty::rc_ptr< ::mosek::fusion::ModelVariable > __mosek_2fusion_2ModelVariable__clone(monty::rc_ptr< ::mosek::fusion::Model > _2245_m) { throw monty::AbstractClassError("Call to abstract method"); }
    }; // struct ModelVariable;

    struct p_RangedVariable : public ::mosek::fusion::p_ModelVariable
    {
      RangedVariable * _pubthis;
      static mosek::fusion::p_RangedVariable* _get_impl(mosek::fusion::RangedVariable * _inst){ return static_cast< mosek::fusion::p_RangedVariable* >(mosek::fusion::p_ModelVariable::_get_impl(_inst)); }
      static mosek::fusion::p_RangedVariable * _get_impl(mosek::fusion::RangedVariable::t _inst) { return _get_impl(_inst.get()); }
      p_RangedVariable(RangedVariable * _pubthis);
      virtual ~p_RangedVariable() { /* std::cout << "~p_RangedVariable" << std::endl;*/ };
      std::shared_ptr< monty::ndarray< int32_t,1 > > shape{};
      std::string name{};
      bool names_flushed{};
      std::shared_ptr< monty::ndarray< int32_t,1 > > nativeidxs{};
      std::shared_ptr< monty::ndarray< int64_t,1 > > sparsity{};

      virtual void destroy();

      static RangedVariable::t _new_RangedVariable(monty::rc_ptr< ::mosek::fusion::RangedVariable > _2130_v,monty::rc_ptr< ::mosek::fusion::Model > _2131_m);
      void _initialize(monty::rc_ptr< ::mosek::fusion::RangedVariable > _2130_v,monty::rc_ptr< ::mosek::fusion::Model > _2131_m);
      static RangedVariable::t _new_RangedVariable(monty::rc_ptr< ::mosek::fusion::Model > _2132_model,const std::string &  _2133_name,int64_t _2134_varid,std::shared_ptr< monty::ndarray< int32_t,1 > > _2135_shape,std::shared_ptr< monty::ndarray< int64_t,1 > > _2136_sparsity,std::shared_ptr< monty::ndarray< int32_t,1 > > _2137_nativeidxs);
      void _initialize(monty::rc_ptr< ::mosek::fusion::Model > _2132_model,const std::string &  _2133_name,int64_t _2134_varid,std::shared_ptr< monty::ndarray< int32_t,1 > > _2135_shape,std::shared_ptr< monty::ndarray< int64_t,1 > > _2136_sparsity,std::shared_ptr< monty::ndarray< int32_t,1 > > _2137_nativeidxs);
      virtual monty::rc_ptr< ::mosek::fusion::Utils::StringBuffer > __mosek_2fusion_2RangedVariable__elementDesc(int64_t _2138_index,monty::rc_ptr< ::mosek::fusion::Utils::StringBuffer > _2139_sb) ;
      virtual /* override */ void flushNames() ;
      virtual void dual_u(int32_t _2140_offset,std::shared_ptr< monty::ndarray< double,1 > > _2141_target) ;
      virtual void dual_l(int32_t _2142_offset,std::shared_ptr< monty::ndarray< double,1 > > _2143_target) ;
      virtual monty::rc_ptr< ::mosek::fusion::BoundInterfaceVariable > __mosek_2fusion_2RangedVariable__upperBoundVar() ;
      virtual monty::rc_ptr< ::mosek::fusion::BoundInterfaceVariable > __mosek_2fusion_2RangedVariable__lowerBoundVar() ;
      virtual /* override */ monty::rc_ptr< ::mosek::fusion::ModelVariable > __mosek_2fusion_2RangedVariable__clone(monty::rc_ptr< ::mosek::fusion::Model > _2146_m) ;
      virtual monty::rc_ptr< ::mosek::fusion::ModelVariable > __mosek_2fusion_2ModelVariable__clone(monty::rc_ptr< ::mosek::fusion::Model > _2146_m) { return __mosek_2fusion_2RangedVariable__clone(_2146_m); }
      static  std::shared_ptr< monty::ndarray< int64_t,1 > > globalNativeIndexes(std::shared_ptr< monty::ndarray< int32_t,1 > > _2147_nativeidxs);
    }; // struct RangedVariable;

    struct p_LinearPSDVariable : public ::mosek::fusion::p_ModelVariable
    {
      LinearPSDVariable * _pubthis;
      static mosek::fusion::p_LinearPSDVariable* _get_impl(mosek::fusion::LinearPSDVariable * _inst){ return static_cast< mosek::fusion::p_LinearPSDVariable* >(mosek::fusion::p_ModelVariable::_get_impl(_inst)); }
      static mosek::fusion::p_LinearPSDVariable * _get_impl(mosek::fusion::LinearPSDVariable::t _inst) { return _get_impl(_inst.get()); }
      p_LinearPSDVariable(LinearPSDVariable * _pubthis);
      virtual ~p_LinearPSDVariable() { /* std::cout << "~p_LinearPSDVariable" << std::endl;*/ };
      std::shared_ptr< monty::ndarray< int32_t,1 > > shape{};
      std::string name{};
      int32_t varid{};
      std::shared_ptr< monty::ndarray< int64_t,1 > > nativeidxs{};
      int32_t conedim{};

      virtual void destroy();

      static LinearPSDVariable::t _new_LinearPSDVariable(monty::rc_ptr< ::mosek::fusion::LinearPSDVariable > _2150_v,monty::rc_ptr< ::mosek::fusion::Model > _2151_m);
      void _initialize(monty::rc_ptr< ::mosek::fusion::LinearPSDVariable > _2150_v,monty::rc_ptr< ::mosek::fusion::Model > _2151_m);
      static LinearPSDVariable::t _new_LinearPSDVariable(monty::rc_ptr< ::mosek::fusion::Model > _2152_model,const std::string &  _2153_name,int32_t _2154_varid,std::shared_ptr< monty::ndarray< int32_t,1 > > _2155_shape,int32_t _2156_conedim,std::shared_ptr< monty::ndarray< int64_t,1 > > _2157_nativeidxs);
      void _initialize(monty::rc_ptr< ::mosek::fusion::Model > _2152_model,const std::string &  _2153_name,int32_t _2154_varid,std::shared_ptr< monty::ndarray< int32_t,1 > > _2155_shape,int32_t _2156_conedim,std::shared_ptr< monty::ndarray< int64_t,1 > > _2157_nativeidxs);
      virtual /* override */ void flushNames() ;
      virtual /* override */ std::string toString() ;
      virtual void make_continuous(std::shared_ptr< monty::ndarray< int64_t,1 > > _2160_idxs) ;
      virtual void make_integer(std::shared_ptr< monty::ndarray< int64_t,1 > > _2161_idxs) ;
      virtual /* override */ monty::rc_ptr< ::mosek::fusion::ModelVariable > __mosek_2fusion_2LinearPSDVariable__clone(monty::rc_ptr< ::mosek::fusion::Model > _2162_m) ;
      virtual monty::rc_ptr< ::mosek::fusion::ModelVariable > __mosek_2fusion_2ModelVariable__clone(monty::rc_ptr< ::mosek::fusion::Model > _2162_m) { return __mosek_2fusion_2LinearPSDVariable__clone(_2162_m); }
      static  std::shared_ptr< monty::ndarray< int64_t,1 > > globalNativeIndexes(std::shared_ptr< monty::ndarray< int64_t,1 > > _2163_nativeidxs);
    }; // struct LinearPSDVariable;

    struct p_PSDVariable : public ::mosek::fusion::p_ModelVariable
    {
      PSDVariable * _pubthis;
      static mosek::fusion::p_PSDVariable* _get_impl(mosek::fusion::PSDVariable * _inst){ return static_cast< mosek::fusion::p_PSDVariable* >(mosek::fusion::p_ModelVariable::_get_impl(_inst)); }
      static mosek::fusion::p_PSDVariable * _get_impl(mosek::fusion::PSDVariable::t _inst) { return _get_impl(_inst.get()); }
      p_PSDVariable(PSDVariable * _pubthis);
      virtual ~p_PSDVariable() { /* std::cout << "~p_PSDVariable" << std::endl;*/ };
      monty::rc_ptr< ::mosek::fusion::Model > model{};
      bool names_flushed{};
      std::shared_ptr< monty::ndarray< int32_t,1 > > barvaridxs{};
      int32_t conedim2{};
      int32_t conedim1{};
      std::shared_ptr< monty::ndarray< int32_t,1 > > shape{};
      std::string name{};
      std::shared_ptr< monty::ndarray< int64_t,1 > > nativeidxs{};
      int32_t varid{};

      virtual void destroy();

      static PSDVariable::t _new_PSDVariable(monty::rc_ptr< ::mosek::fusion::PSDVariable > _2165_v,monty::rc_ptr< ::mosek::fusion::Model > _2166_m);
      void _initialize(monty::rc_ptr< ::mosek::fusion::PSDVariable > _2165_v,monty::rc_ptr< ::mosek::fusion::Model > _2166_m);
      static PSDVariable::t _new_PSDVariable(monty::rc_ptr< ::mosek::fusion::Model > _2167_model,const std::string &  _2168_name,int32_t _2169_varid,std::shared_ptr< monty::ndarray< int32_t,1 > > _2170_shape,int32_t _2171_conedim1,int32_t _2172_conedim2,std::shared_ptr< monty::ndarray< int32_t,1 > > _2173_barvaridxs,std::shared_ptr< monty::ndarray< int64_t,1 > > _2174_nativeidxs);
      void _initialize(monty::rc_ptr< ::mosek::fusion::Model > _2167_model,const std::string &  _2168_name,int32_t _2169_varid,std::shared_ptr< monty::ndarray< int32_t,1 > > _2170_shape,int32_t _2171_conedim1,int32_t _2172_conedim2,std::shared_ptr< monty::ndarray< int32_t,1 > > _2173_barvaridxs,std::shared_ptr< monty::ndarray< int64_t,1 > > _2174_nativeidxs);
      virtual /* override */ void flushNames() ;
      virtual /* override */ std::string toString() ;
      virtual monty::rc_ptr< ::mosek::fusion::Utils::StringBuffer > __mosek_2fusion_2PSDVariable__elementDesc(int64_t _2177_index,monty::rc_ptr< ::mosek::fusion::Utils::StringBuffer > _2178_sb) ;
      virtual /* override */ monty::rc_ptr< ::mosek::fusion::ModelVariable > __mosek_2fusion_2PSDVariable__clone(monty::rc_ptr< ::mosek::fusion::Model > _2179_m) ;
      virtual monty::rc_ptr< ::mosek::fusion::ModelVariable > __mosek_2fusion_2ModelVariable__clone(monty::rc_ptr< ::mosek::fusion::Model > _2179_m) { return __mosek_2fusion_2PSDVariable__clone(_2179_m); }
      static  std::shared_ptr< monty::ndarray< int64_t,1 > > fullnativeidxs(std::shared_ptr< monty::ndarray< int32_t,1 > > _2180_shape,int32_t _2181_conedim1,int32_t _2182_conedim2,std::shared_ptr< monty::ndarray< int64_t,1 > > _2183_nativeidxs);
    }; // struct PSDVariable;

    struct p_LinearVariable : public ::mosek::fusion::p_ModelVariable
    {
      LinearVariable * _pubthis;
      static mosek::fusion::p_LinearVariable* _get_impl(mosek::fusion::LinearVariable * _inst){ return static_cast< mosek::fusion::p_LinearVariable* >(mosek::fusion::p_ModelVariable::_get_impl(_inst)); }
      static mosek::fusion::p_LinearVariable * _get_impl(mosek::fusion::LinearVariable::t _inst) { return _get_impl(_inst.get()); }
      p_LinearVariable(LinearVariable * _pubthis);
      virtual ~p_LinearVariable() { /* std::cout << "~p_LinearVariable" << std::endl;*/ };
      std::shared_ptr< monty::ndarray< int32_t,1 > > shape{};
      std::shared_ptr< monty::ndarray< int64_t,1 > > sparsity{};
      std::shared_ptr< monty::ndarray< int32_t,1 > > nativeidxs{};
      bool names_flushed{};
      std::string name{};

      virtual void destroy();

      static LinearVariable::t _new_LinearVariable(monty::rc_ptr< ::mosek::fusion::LinearVariable > _2208_v,monty::rc_ptr< ::mosek::fusion::Model > _2209_m);
      void _initialize(monty::rc_ptr< ::mosek::fusion::LinearVariable > _2208_v,monty::rc_ptr< ::mosek::fusion::Model > _2209_m);
      static LinearVariable::t _new_LinearVariable(monty::rc_ptr< ::mosek::fusion::Model > _2210_model,const std::string &  _2211_name,int64_t _2212_varid,std::shared_ptr< monty::ndarray< int32_t,1 > > _2213_shape,std::shared_ptr< monty::ndarray< int64_t,1 > > _2214_sparsity,std::shared_ptr< monty::ndarray< int32_t,1 > > _2215_nativeidxs);
      void _initialize(monty::rc_ptr< ::mosek::fusion::Model > _2210_model,const std::string &  _2211_name,int64_t _2212_varid,std::shared_ptr< monty::ndarray< int32_t,1 > > _2213_shape,std::shared_ptr< monty::ndarray< int64_t,1 > > _2214_sparsity,std::shared_ptr< monty::ndarray< int32_t,1 > > _2215_nativeidxs);
      virtual /* override */ std::string toString() ;
      virtual /* override */ void flushNames() ;
      virtual /* override */ monty::rc_ptr< ::mosek::fusion::ModelVariable > __mosek_2fusion_2LinearVariable__clone(monty::rc_ptr< ::mosek::fusion::Model > _2218_m) ;
      virtual monty::rc_ptr< ::mosek::fusion::ModelVariable > __mosek_2fusion_2ModelVariable__clone(monty::rc_ptr< ::mosek::fusion::Model > _2218_m) { return __mosek_2fusion_2LinearVariable__clone(_2218_m); }
      static  std::shared_ptr< monty::ndarray< int64_t,1 > > globalNativeIndexes(std::shared_ptr< monty::ndarray< int32_t,1 > > _2219_nativeidxs);
    }; // struct LinearVariable;

    struct p_ConicVariable : public ::mosek::fusion::p_ModelVariable
    {
      ConicVariable * _pubthis;
      static mosek::fusion::p_ConicVariable* _get_impl(mosek::fusion::ConicVariable * _inst){ return static_cast< mosek::fusion::p_ConicVariable* >(mosek::fusion::p_ModelVariable::_get_impl(_inst)); }
      static mosek::fusion::p_ConicVariable * _get_impl(mosek::fusion::ConicVariable::t _inst) { return _get_impl(_inst.get()); }
      p_ConicVariable(ConicVariable * _pubthis);
      virtual ~p_ConicVariable() { /* std::cout << "~p_ConicVariable" << std::endl;*/ };
      std::shared_ptr< monty::ndarray< int32_t,1 > > nativeidxs{};
      std::shared_ptr< monty::ndarray< int32_t,1 > > shape{};
      std::string name{};
      bool names_flushed{};
      int32_t varid{};

      virtual void destroy();

      static ConicVariable::t _new_ConicVariable(monty::rc_ptr< ::mosek::fusion::ConicVariable > _2222_v,monty::rc_ptr< ::mosek::fusion::Model > _2223_m);
      void _initialize(monty::rc_ptr< ::mosek::fusion::ConicVariable > _2222_v,monty::rc_ptr< ::mosek::fusion::Model > _2223_m);
      static ConicVariable::t _new_ConicVariable(monty::rc_ptr< ::mosek::fusion::Model > _2224_model,const std::string &  _2225_name,int32_t _2226_varid,std::shared_ptr< monty::ndarray< int32_t,1 > > _2227_shape,std::shared_ptr< monty::ndarray< int32_t,1 > > _2228_nativeidxs);
      void _initialize(monty::rc_ptr< ::mosek::fusion::Model > _2224_model,const std::string &  _2225_name,int32_t _2226_varid,std::shared_ptr< monty::ndarray< int32_t,1 > > _2227_shape,std::shared_ptr< monty::ndarray< int32_t,1 > > _2228_nativeidxs);
      virtual /* override */ std::string toString() ;
      virtual /* override */ void flushNames() ;
      virtual /* override */ monty::rc_ptr< ::mosek::fusion::ModelVariable > __mosek_2fusion_2ConicVariable__clone(monty::rc_ptr< ::mosek::fusion::Model > _2231_m) ;
      virtual monty::rc_ptr< ::mosek::fusion::ModelVariable > __mosek_2fusion_2ModelVariable__clone(monty::rc_ptr< ::mosek::fusion::Model > _2231_m) { return __mosek_2fusion_2ConicVariable__clone(_2231_m); }
      static  std::shared_ptr< monty::ndarray< int64_t,1 > > globalNativeIndexes(std::shared_ptr< monty::ndarray< int32_t,1 > > _2232_nativeidxs);
    }; // struct ConicVariable;

    struct p_NilVariable : public ::mosek::fusion::p_BaseVariable
    {
      NilVariable * _pubthis;
      static mosek::fusion::p_NilVariable* _get_impl(mosek::fusion::NilVariable * _inst){ return static_cast< mosek::fusion::p_NilVariable* >(mosek::fusion::p_BaseVariable::_get_impl(_inst)); }
      static mosek::fusion::p_NilVariable * _get_impl(mosek::fusion::NilVariable::t _inst) { return _get_impl(_inst.get()); }
      p_NilVariable(NilVariable * _pubthis);
      virtual ~p_NilVariable() { /* std::cout << "~p_NilVariable" << std::endl;*/ };
      std::shared_ptr< monty::ndarray< int32_t,1 > > shape{};

      virtual void destroy();

      static NilVariable::t _new_NilVariable(std::shared_ptr< monty::ndarray< int32_t,1 > > _2246_shape);
      void _initialize(std::shared_ptr< monty::ndarray< int32_t,1 > > _2246_shape);
      static NilVariable::t _new_NilVariable();
      void _initialize();
      virtual void flushNames() ;
      virtual monty::rc_ptr< ::mosek::fusion::Utils::StringBuffer > __mosek_2fusion_2NilVariable__elementDesc(int64_t _2248_index,monty::rc_ptr< ::mosek::fusion::Utils::StringBuffer > _2249_sb) ;
      virtual void elementName(int64_t _2250_index,monty::rc_ptr< ::mosek::fusion::Utils::StringBuffer > _2251_sb) ;
      virtual /* override */ int32_t numInst() ;
      virtual int32_t inst(int32_t _2252_offset,std::shared_ptr< monty::ndarray< int64_t,1 > > _2253_sparsity,std::shared_ptr< monty::ndarray< int64_t,1 > > _2254_basevar_nativeidxs) ;
      virtual /* override */ void inst(int32_t _2255_offset,std::shared_ptr< monty::ndarray< int64_t,1 > > _2256_nindex) ;
      virtual /* override */ void set_values(std::shared_ptr< monty::ndarray< double,1 > > _2257_target,bool _2258_primal) ;
      virtual /* override */ void values(int32_t _2259_offset,std::shared_ptr< monty::ndarray< double,1 > > _2260_target,bool _2261_primal) ;
      virtual /* override */ void make_continuous() ;
      virtual /* override */ void make_integer() ;
      virtual /* override */ std::string toString() ;
      virtual /* override */ monty::rc_ptr< ::mosek::fusion::Variable > __mosek_2fusion_2NilVariable__index(std::shared_ptr< monty::ndarray< int32_t,1 > > _2262_first) ;
      virtual monty::rc_ptr< ::mosek::fusion::Variable > __mosek_2fusion_2BaseVariable__index(std::shared_ptr< monty::ndarray< int32_t,1 > > _2262_first) { return __mosek_2fusion_2NilVariable__index(_2262_first); }
      virtual /* override */ monty::rc_ptr< ::mosek::fusion::Variable > __mosek_2fusion_2NilVariable__index(int32_t _2264_first) ;
      virtual monty::rc_ptr< ::mosek::fusion::Variable > __mosek_2fusion_2BaseVariable__index(int32_t _2264_first) { return __mosek_2fusion_2NilVariable__index(_2264_first); }
      virtual /* override */ monty::rc_ptr< ::mosek::fusion::Variable > __mosek_2fusion_2NilVariable__slice(std::shared_ptr< monty::ndarray< int32_t,1 > > _2266_first,std::shared_ptr< monty::ndarray< int32_t,1 > > _2267_last) ;
      virtual monty::rc_ptr< ::mosek::fusion::Variable > __mosek_2fusion_2BaseVariable__slice(std::shared_ptr< monty::ndarray< int32_t,1 > > _2266_first,std::shared_ptr< monty::ndarray< int32_t,1 > > _2267_last) { return __mosek_2fusion_2NilVariable__slice(_2266_first,_2267_last); }
      virtual /* override */ monty::rc_ptr< ::mosek::fusion::Variable > __mosek_2fusion_2NilVariable__slice(int32_t _2270_first,int32_t _2271_last) ;
      virtual monty::rc_ptr< ::mosek::fusion::Variable > __mosek_2fusion_2BaseVariable__slice(int32_t _2270_first,int32_t _2271_last) { return __mosek_2fusion_2NilVariable__slice(_2270_first,_2271_last); }
    }; // struct NilVariable;

    struct p_Var
    {
      Var * _pubthis;
      static mosek::fusion::p_Var* _get_impl(mosek::fusion::Var * _inst){ assert(_inst); assert(_inst->_impl); return _inst->_impl; }
      static mosek::fusion::p_Var * _get_impl(mosek::fusion::Var::t _inst) { return _get_impl(_inst.get()); }
      p_Var(Var * _pubthis);
      virtual ~p_Var() { /* std::cout << "~p_Var" << std::endl;*/ };

      virtual void destroy();

      static  monty::rc_ptr< ::mosek::fusion::Variable > empty(std::shared_ptr< monty::ndarray< int32_t,1 > > _2617_shape);
      static  monty::rc_ptr< ::mosek::fusion::Variable > compress(monty::rc_ptr< ::mosek::fusion::Variable > _2619_v);
      static  monty::rc_ptr< ::mosek::fusion::Variable > reshape(monty::rc_ptr< ::mosek::fusion::Variable > _2627_v,int32_t _2628_d1);
      static  monty::rc_ptr< ::mosek::fusion::Variable > reshape(monty::rc_ptr< ::mosek::fusion::Variable > _2629_v,int32_t _2630_d1,int32_t _2631_d2);
      static  monty::rc_ptr< ::mosek::fusion::Variable > flatten(monty::rc_ptr< ::mosek::fusion::Variable > _2632_v);
      static  monty::rc_ptr< ::mosek::fusion::Variable > reshape(monty::rc_ptr< ::mosek::fusion::Variable > _2633_v,std::shared_ptr< monty::ndarray< int32_t,1 > > _2634_shape);
      static  monty::rc_ptr< ::mosek::fusion::Variable > index_permute_(monty::rc_ptr< ::mosek::fusion::Variable > _2635_v,std::shared_ptr< monty::ndarray< int32_t,1 > > _2636_perm);
      static  monty::rc_ptr< ::mosek::fusion::Variable > hrepeat(monty::rc_ptr< ::mosek::fusion::Variable > _2665_v,int32_t _2666_n);
      static  monty::rc_ptr< ::mosek::fusion::Variable > vrepeat(monty::rc_ptr< ::mosek::fusion::Variable > _2667_v,int32_t _2668_n);
      static  monty::rc_ptr< ::mosek::fusion::Variable > repeat(monty::rc_ptr< ::mosek::fusion::Variable > _2669_v,int32_t _2670_n);
      static  monty::rc_ptr< ::mosek::fusion::Variable > repeat(monty::rc_ptr< ::mosek::fusion::Variable > _2671_v,int32_t _2672_dim,int32_t _2673_n);
      static  monty::rc_ptr< ::mosek::fusion::Variable > drepeat(monty::rc_ptr< ::mosek::fusion::Variable > _2674_v,int32_t _2675_dim,int32_t _2676_n);
      static  monty::rc_ptr< ::mosek::fusion::Variable > stack(std::shared_ptr< monty::ndarray< std::shared_ptr< monty::ndarray< monty::rc_ptr< ::mosek::fusion::Variable >,1 > >,1 > > _2740_vlist);
      static  monty::rc_ptr< ::mosek::fusion::Variable > vstack(monty::rc_ptr< ::mosek::fusion::Variable > _2742_v1,monty::rc_ptr< ::mosek::fusion::Variable > _2743_v2,monty::rc_ptr< ::mosek::fusion::Variable > _2744_v3);
      static  monty::rc_ptr< ::mosek::fusion::Variable > vstack(monty::rc_ptr< ::mosek::fusion::Variable > _2745_v1,monty::rc_ptr< ::mosek::fusion::Variable > _2746_v2);
      static  monty::rc_ptr< ::mosek::fusion::Variable > vstack(std::shared_ptr< monty::ndarray< monty::rc_ptr< ::mosek::fusion::Variable >,1 > > _2747_v);
      static  monty::rc_ptr< ::mosek::fusion::Variable > hstack(monty::rc_ptr< ::mosek::fusion::Variable > _2748_v1,monty::rc_ptr< ::mosek::fusion::Variable > _2749_v2,monty::rc_ptr< ::mosek::fusion::Variable > _2750_v3);
      static  monty::rc_ptr< ::mosek::fusion::Variable > hstack(monty::rc_ptr< ::mosek::fusion::Variable > _2751_v1,monty::rc_ptr< ::mosek::fusion::Variable > _2752_v2);
      static  monty::rc_ptr< ::mosek::fusion::Variable > hstack(std::shared_ptr< monty::ndarray< monty::rc_ptr< ::mosek::fusion::Variable >,1 > > _2753_v);
      static  monty::rc_ptr< ::mosek::fusion::Variable > stack(std::shared_ptr< monty::ndarray< monty::rc_ptr< ::mosek::fusion::Variable >,1 > > _2754_v,int32_t _2755_dim);
      static  monty::rc_ptr< ::mosek::fusion::Variable > stack(monty::rc_ptr< ::mosek::fusion::Variable > _2756_v1,monty::rc_ptr< ::mosek::fusion::Variable > _2757_v2,monty::rc_ptr< ::mosek::fusion::Variable > _2758_v3,int32_t _2759_dim);
      static  monty::rc_ptr< ::mosek::fusion::Variable > stack(monty::rc_ptr< ::mosek::fusion::Variable > _2760_v1,monty::rc_ptr< ::mosek::fusion::Variable > _2761_v2,int32_t _2762_dim);
      static  monty::rc_ptr< ::mosek::fusion::Variable > stack(int32_t _2763_dim,std::shared_ptr< monty::ndarray< monty::rc_ptr< ::mosek::fusion::Variable >,1 > > _2764_v);
      static  monty::rc_ptr< ::mosek::fusion::Variable > stack(int32_t _2767_dim,monty::rc_ptr< ::mosek::fusion::Variable > _2768_v1,monty::rc_ptr< ::mosek::fusion::Variable > _2769_v2,monty::rc_ptr< ::mosek::fusion::Variable > _2770_v3);
      static  monty::rc_ptr< ::mosek::fusion::Variable > stack(int32_t _2771_dim,monty::rc_ptr< ::mosek::fusion::Variable > _2772_v1,monty::rc_ptr< ::mosek::fusion::Variable > _2773_v2);
      static  monty::rc_ptr< ::mosek::fusion::Variable > promote(monty::rc_ptr< ::mosek::fusion::Variable > _2774_v,int32_t _2775_nd);
      static  monty::rc_ptr< ::mosek::fusion::Variable > dstack(std::shared_ptr< monty::ndarray< monty::rc_ptr< ::mosek::fusion::Variable >,1 > > _2780_v,int32_t _2781_dim);
    }; // struct Var;

    struct p_Constraint
    {
      Constraint * _pubthis;
      static mosek::fusion::p_Constraint* _get_impl(mosek::fusion::Constraint * _inst){ assert(_inst); assert(_inst->_impl); return _inst->_impl; }
      static mosek::fusion::p_Constraint * _get_impl(mosek::fusion::Constraint::t _inst) { return _get_impl(_inst.get()); }
      p_Constraint(Constraint * _pubthis);
      virtual ~p_Constraint() { /* std::cout << "~p_Constraint" << std::endl;*/ };
      std::shared_ptr< monty::ndarray< int32_t,1 > > con_nativeidxs{};
      std::shared_ptr< monty::ndarray< int32_t,1 > > shape{};
      monty::rc_ptr< ::mosek::fusion::Model > model{};

      virtual void destroy();

      static Constraint::t _new_Constraint(monty::rc_ptr< ::mosek::fusion::Constraint > _2972_c,monty::rc_ptr< ::mosek::fusion::Model > _2973_m);
      void _initialize(monty::rc_ptr< ::mosek::fusion::Constraint > _2972_c,monty::rc_ptr< ::mosek::fusion::Model > _2973_m);
      static Constraint::t _new_Constraint(monty::rc_ptr< ::mosek::fusion::Model > _2974_model,std::shared_ptr< monty::ndarray< int32_t,1 > > _2975_shape,std::shared_ptr< monty::ndarray< int32_t,1 > > _2976_con_nativeidxs);
      void _initialize(monty::rc_ptr< ::mosek::fusion::Model > _2974_model,std::shared_ptr< monty::ndarray< int32_t,1 > > _2975_shape,std::shared_ptr< monty::ndarray< int32_t,1 > > _2976_con_nativeidxs);
      virtual /* override */ std::string toString() ;
      virtual void toStringArray(std::shared_ptr< monty::ndarray< int64_t,1 > > _2977_subi,int64_t _2978_dstidx,std::shared_ptr< monty::ndarray< std::string,1 > > _2979_result) ;
      virtual void dual_lu(int32_t _2980_offset,std::shared_ptr< monty::ndarray< double,1 > > _2981_target,bool _2982_islower) ;
      virtual std::shared_ptr< monty::ndarray< double,1 > > dual() ;
      virtual std::shared_ptr< monty::ndarray< double,1 > > level() ;
      virtual void values(bool _2985_primal,int32_t _2986_offset,std::shared_ptr< monty::ndarray< double,1 > > _2987_target) ;
      virtual void remove() ;
      virtual void update(std::shared_ptr< monty::ndarray< double,1 > > _2988_bfix) ;
      virtual void update(monty::rc_ptr< ::mosek::fusion::Expression > _2989_expr) ;
      virtual void update(monty::rc_ptr< ::mosek::fusion::Expression > _2993_expr,monty::rc_ptr< ::mosek::fusion::Variable > _2994_x,bool _2995_bfixupdate) ;
      virtual void update(monty::rc_ptr< ::mosek::fusion::Expression > _3015_expr,monty::rc_ptr< ::mosek::fusion::Variable > _3016_x) ;
      virtual monty::rc_ptr< ::mosek::fusion::Model > __mosek_2fusion_2Constraint__get_model() ;
      virtual int32_t get_nd() ;
      virtual int64_t size() ;
      static  monty::rc_ptr< ::mosek::fusion::Constraint > stack(std::shared_ptr< monty::ndarray< monty::rc_ptr< ::mosek::fusion::Constraint >,1 > > _3019_clist,int32_t _3020_dim);
      static  monty::rc_ptr< ::mosek::fusion::Constraint > stack(monty::rc_ptr< ::mosek::fusion::Constraint > _3021_v1,monty::rc_ptr< ::mosek::fusion::Constraint > _3022_v2,monty::rc_ptr< ::mosek::fusion::Constraint > _3023_v3,int32_t _3024_dim);
      static  monty::rc_ptr< ::mosek::fusion::Constraint > stack(monty::rc_ptr< ::mosek::fusion::Constraint > _3025_v1,monty::rc_ptr< ::mosek::fusion::Constraint > _3026_v2,int32_t _3027_dim);
      static  monty::rc_ptr< ::mosek::fusion::Constraint > hstack(std::shared_ptr< monty::ndarray< monty::rc_ptr< ::mosek::fusion::Constraint >,1 > > _3028_clist);
      static  monty::rc_ptr< ::mosek::fusion::Constraint > vstack(std::shared_ptr< monty::ndarray< monty::rc_ptr< ::mosek::fusion::Constraint >,1 > > _3029_clist);
      static  monty::rc_ptr< ::mosek::fusion::Constraint > hstack(monty::rc_ptr< ::mosek::fusion::Constraint > _3030_v1,monty::rc_ptr< ::mosek::fusion::Constraint > _3031_v2,monty::rc_ptr< ::mosek::fusion::Constraint > _3032_v3);
      static  monty::rc_ptr< ::mosek::fusion::Constraint > vstack(monty::rc_ptr< ::mosek::fusion::Constraint > _3033_v1,monty::rc_ptr< ::mosek::fusion::Constraint > _3034_v2,monty::rc_ptr< ::mosek::fusion::Constraint > _3035_v3);
      static  monty::rc_ptr< ::mosek::fusion::Constraint > hstack(monty::rc_ptr< ::mosek::fusion::Constraint > _3036_v1,monty::rc_ptr< ::mosek::fusion::Constraint > _3037_v2);
      static  monty::rc_ptr< ::mosek::fusion::Constraint > vstack(monty::rc_ptr< ::mosek::fusion::Constraint > _3038_v1,monty::rc_ptr< ::mosek::fusion::Constraint > _3039_v2);
      static  monty::rc_ptr< ::mosek::fusion::Constraint > dstack(std::shared_ptr< monty::ndarray< monty::rc_ptr< ::mosek::fusion::Constraint >,1 > > _3040_c,int32_t _3041_dim);
      virtual monty::rc_ptr< ::mosek::fusion::Constraint > __mosek_2fusion_2Constraint__index(std::shared_ptr< monty::ndarray< int32_t,1 > > _3092_idxa) ;
      virtual monty::rc_ptr< ::mosek::fusion::Constraint > __mosek_2fusion_2Constraint__index(int32_t _3099_idx) ;
      virtual monty::rc_ptr< ::mosek::fusion::Constraint > __mosek_2fusion_2Constraint__slice(std::shared_ptr< monty::ndarray< int32_t,1 > > _3100_firsta,std::shared_ptr< monty::ndarray< int32_t,1 > > _3101_lasta) ;
      virtual monty::rc_ptr< ::mosek::fusion::Constraint > __mosek_2fusion_2Constraint__slice(int32_t _3120_first,int32_t _3121_last) ;
      virtual int32_t getND() ;
      virtual int32_t getSize() ;
      virtual monty::rc_ptr< ::mosek::fusion::Model > __mosek_2fusion_2Constraint__getModel() ;
      virtual std::shared_ptr< monty::ndarray< int32_t,1 > > getShape() ;
      virtual std::shared_ptr< monty::ndarray< int32_t,1 > > getNativeidxs() ;
    }; // struct Constraint;

    struct p_SliceConstraint : public ::mosek::fusion::p_Constraint
    {
      SliceConstraint * _pubthis;
      static mosek::fusion::p_SliceConstraint* _get_impl(mosek::fusion::SliceConstraint * _inst){ return static_cast< mosek::fusion::p_SliceConstraint* >(mosek::fusion::p_Constraint::_get_impl(_inst)); }
      static mosek::fusion::p_SliceConstraint * _get_impl(mosek::fusion::SliceConstraint::t _inst) { return _get_impl(_inst.get()); }
      p_SliceConstraint(SliceConstraint * _pubthis);
      virtual ~p_SliceConstraint() { /* std::cout << "~p_SliceConstraint" << std::endl;*/ };

      virtual void destroy();

      static SliceConstraint::t _new_SliceConstraint(monty::rc_ptr< ::mosek::fusion::SliceConstraint > _2924_c);
      void _initialize(monty::rc_ptr< ::mosek::fusion::SliceConstraint > _2924_c);
      static SliceConstraint::t _new_SliceConstraint(monty::rc_ptr< ::mosek::fusion::Model > _2925_model,std::shared_ptr< monty::ndarray< int32_t,1 > > _2926_shape,std::shared_ptr< monty::ndarray< int32_t,1 > > _2927_nativeidxs);
      void _initialize(monty::rc_ptr< ::mosek::fusion::Model > _2925_model,std::shared_ptr< monty::ndarray< int32_t,1 > > _2926_shape,std::shared_ptr< monty::ndarray< int32_t,1 > > _2927_nativeidxs);
      virtual /* override */ std::string toString() ;
    }; // struct SliceConstraint;

    struct p_BoundInterfaceConstraint : public ::mosek::fusion::p_SliceConstraint
    {
      BoundInterfaceConstraint * _pubthis;
      static mosek::fusion::p_BoundInterfaceConstraint* _get_impl(mosek::fusion::BoundInterfaceConstraint * _inst){ return static_cast< mosek::fusion::p_BoundInterfaceConstraint* >(mosek::fusion::p_SliceConstraint::_get_impl(_inst)); }
      static mosek::fusion::p_BoundInterfaceConstraint * _get_impl(mosek::fusion::BoundInterfaceConstraint::t _inst) { return _get_impl(_inst.get()); }
      p_BoundInterfaceConstraint(BoundInterfaceConstraint * _pubthis);
      virtual ~p_BoundInterfaceConstraint() { /* std::cout << "~p_BoundInterfaceConstraint" << std::endl;*/ };
      bool islower{};

      virtual void destroy();

      static BoundInterfaceConstraint::t _new_BoundInterfaceConstraint(monty::rc_ptr< ::mosek::fusion::Model > _2850_m,std::shared_ptr< monty::ndarray< int32_t,1 > > _2851_shape,std::shared_ptr< monty::ndarray< int32_t,1 > > _2852_nativeidxs,bool _2853_islower);
      void _initialize(monty::rc_ptr< ::mosek::fusion::Model > _2850_m,std::shared_ptr< monty::ndarray< int32_t,1 > > _2851_shape,std::shared_ptr< monty::ndarray< int32_t,1 > > _2852_nativeidxs,bool _2853_islower);
      static BoundInterfaceConstraint::t _new_BoundInterfaceConstraint(monty::rc_ptr< ::mosek::fusion::SliceConstraint > _2854_c,bool _2855_islower);
      void _initialize(monty::rc_ptr< ::mosek::fusion::SliceConstraint > _2854_c,bool _2855_islower);
      virtual /* override */ std::shared_ptr< monty::ndarray< double,1 > > dual() ;
      virtual /* override */ monty::rc_ptr< ::mosek::fusion::Constraint > __mosek_2fusion_2BoundInterfaceConstraint__slice(std::shared_ptr< monty::ndarray< int32_t,1 > > _2857_firsta,std::shared_ptr< monty::ndarray< int32_t,1 > > _2858_lasta) ;
      virtual monty::rc_ptr< ::mosek::fusion::Constraint > __mosek_2fusion_2Constraint__slice(std::shared_ptr< monty::ndarray< int32_t,1 > > _2857_firsta,std::shared_ptr< monty::ndarray< int32_t,1 > > _2858_lasta) { return __mosek_2fusion_2BoundInterfaceConstraint__slice(_2857_firsta,_2858_lasta); }
      virtual /* override */ monty::rc_ptr< ::mosek::fusion::Constraint > __mosek_2fusion_2BoundInterfaceConstraint__slice(int32_t _2860_first,int32_t _2861_last) ;
      virtual monty::rc_ptr< ::mosek::fusion::Constraint > __mosek_2fusion_2Constraint__slice(int32_t _2860_first,int32_t _2861_last) { return __mosek_2fusion_2BoundInterfaceConstraint__slice(_2860_first,_2861_last); }
      virtual /* override */ monty::rc_ptr< ::mosek::fusion::Constraint > __mosek_2fusion_2BoundInterfaceConstraint__index(std::shared_ptr< monty::ndarray< int32_t,1 > > _2863_idxa) ;
      virtual monty::rc_ptr< ::mosek::fusion::Constraint > __mosek_2fusion_2Constraint__index(std::shared_ptr< monty::ndarray< int32_t,1 > > _2863_idxa) { return __mosek_2fusion_2BoundInterfaceConstraint__index(_2863_idxa); }
      virtual /* override */ monty::rc_ptr< ::mosek::fusion::Constraint > __mosek_2fusion_2BoundInterfaceConstraint__index(int32_t _2865_idx) ;
      virtual monty::rc_ptr< ::mosek::fusion::Constraint > __mosek_2fusion_2Constraint__index(int32_t _2865_idx) { return __mosek_2fusion_2BoundInterfaceConstraint__index(_2865_idx); }
      virtual monty::rc_ptr< ::mosek::fusion::BoundInterfaceConstraint > __mosek_2fusion_2BoundInterfaceConstraint__from_(monty::rc_ptr< ::mosek::fusion::Constraint > _2867_c) ;
    }; // struct BoundInterfaceConstraint;

    struct p_ModelConstraint : public ::mosek::fusion::p_Constraint
    {
      ModelConstraint * _pubthis;
      static mosek::fusion::p_ModelConstraint* _get_impl(mosek::fusion::ModelConstraint * _inst){ return static_cast< mosek::fusion::p_ModelConstraint* >(mosek::fusion::p_Constraint::_get_impl(_inst)); }
      static mosek::fusion::p_ModelConstraint * _get_impl(mosek::fusion::ModelConstraint::t _inst) { return _get_impl(_inst.get()); }
      p_ModelConstraint(ModelConstraint * _pubthis);
      virtual ~p_ModelConstraint() { /* std::cout << "~p_ModelConstraint" << std::endl;*/ };
      int32_t conid{};
      std::shared_ptr< monty::ndarray< int32_t,1 > > shape{};
      std::shared_ptr< monty::ndarray< int32_t,1 > > modelcon_nativeidxs{};
      std::string name{};

      virtual void destroy();

      static ModelConstraint::t _new_ModelConstraint(monty::rc_ptr< ::mosek::fusion::ModelConstraint > _2963_c,monty::rc_ptr< ::mosek::fusion::Model > _2964_m);
      void _initialize(monty::rc_ptr< ::mosek::fusion::ModelConstraint > _2963_c,monty::rc_ptr< ::mosek::fusion::Model > _2964_m);
      static ModelConstraint::t _new_ModelConstraint(monty::rc_ptr< ::mosek::fusion::Model > _2965_model,const std::string &  _2966_name,std::shared_ptr< monty::ndarray< int32_t,1 > > _2967_shape,std::shared_ptr< monty::ndarray< int32_t,1 > > _2968_nidxs,int32_t _2969_conid);
      void _initialize(monty::rc_ptr< ::mosek::fusion::Model > _2965_model,const std::string &  _2966_name,std::shared_ptr< monty::ndarray< int32_t,1 > > _2967_shape,std::shared_ptr< monty::ndarray< int32_t,1 > > _2968_nidxs,int32_t _2969_conid);
      virtual /* override */ std::string toString() ;
      virtual void flushNames() ;
      virtual monty::rc_ptr< ::mosek::fusion::ModelConstraint > __mosek_2fusion_2ModelConstraint__clone(monty::rc_ptr< ::mosek::fusion::Model > _2971_m) { throw monty::AbstractClassError("Call to abstract method"); }
      virtual /* override */ void remove() ;
    }; // struct ModelConstraint;

    struct p_LinearPSDConstraint : public ::mosek::fusion::p_ModelConstraint
    {
      LinearPSDConstraint * _pubthis;
      static mosek::fusion::p_LinearPSDConstraint* _get_impl(mosek::fusion::LinearPSDConstraint * _inst){ return static_cast< mosek::fusion::p_LinearPSDConstraint* >(mosek::fusion::p_ModelConstraint::_get_impl(_inst)); }
      static mosek::fusion::p_LinearPSDConstraint * _get_impl(mosek::fusion::LinearPSDConstraint::t _inst) { return _get_impl(_inst.get()); }
      p_LinearPSDConstraint(LinearPSDConstraint * _pubthis);
      virtual ~p_LinearPSDConstraint() { /* std::cout << "~p_LinearPSDConstraint" << std::endl;*/ };
      int32_t conedim{};
      std::shared_ptr< monty::ndarray< int32_t,1 > > shape{};
      int32_t conid{};
      std::shared_ptr< monty::ndarray< int64_t,1 > > slackidxs{};
      std::shared_ptr< monty::ndarray< int32_t,1 > > nativeidxs{};

      virtual void destroy();

      static LinearPSDConstraint::t _new_LinearPSDConstraint(monty::rc_ptr< ::mosek::fusion::LinearPSDConstraint > _2870_c,monty::rc_ptr< ::mosek::fusion::Model > _2871_m);
      void _initialize(monty::rc_ptr< ::mosek::fusion::LinearPSDConstraint > _2870_c,monty::rc_ptr< ::mosek::fusion::Model > _2871_m);
      static LinearPSDConstraint::t _new_LinearPSDConstraint(monty::rc_ptr< ::mosek::fusion::Model > _2872_model,const std::string &  _2873_name,int32_t _2874_conid,std::shared_ptr< monty::ndarray< int32_t,1 > > _2875_shape,int32_t _2876_conedim,std::shared_ptr< monty::ndarray< int32_t,1 > > _2877_nativeidxs,std::shared_ptr< monty::ndarray< int64_t,1 > > _2878_slackidxs);
      void _initialize(monty::rc_ptr< ::mosek::fusion::Model > _2872_model,const std::string &  _2873_name,int32_t _2874_conid,std::shared_ptr< monty::ndarray< int32_t,1 > > _2875_shape,int32_t _2876_conedim,std::shared_ptr< monty::ndarray< int32_t,1 > > _2877_nativeidxs,std::shared_ptr< monty::ndarray< int64_t,1 > > _2878_slackidxs);
      virtual void domainToString(int64_t _2879_i,monty::rc_ptr< ::mosek::fusion::Utils::StringBuffer > _2880_sb) ;
      virtual /* override */ monty::rc_ptr< ::mosek::fusion::ModelConstraint > __mosek_2fusion_2LinearPSDConstraint__clone(monty::rc_ptr< ::mosek::fusion::Model > _2884_m) ;
      virtual monty::rc_ptr< ::mosek::fusion::ModelConstraint > __mosek_2fusion_2ModelConstraint__clone(monty::rc_ptr< ::mosek::fusion::Model > _2884_m) { return __mosek_2fusion_2LinearPSDConstraint__clone(_2884_m); }
    }; // struct LinearPSDConstraint;

    struct p_PSDConstraint : public ::mosek::fusion::p_ModelConstraint
    {
      PSDConstraint * _pubthis;
      static mosek::fusion::p_PSDConstraint* _get_impl(mosek::fusion::PSDConstraint * _inst){ return static_cast< mosek::fusion::p_PSDConstraint* >(mosek::fusion::p_ModelConstraint::_get_impl(_inst)); }
      static mosek::fusion::p_PSDConstraint * _get_impl(mosek::fusion::PSDConstraint::t _inst) { return _get_impl(_inst.get()); }
      p_PSDConstraint(PSDConstraint * _pubthis);
      virtual ~p_PSDConstraint() { /* std::cout << "~p_PSDConstraint" << std::endl;*/ };
      bool names_flushed{};
      int32_t conedim1{};
      int32_t conedim0{};
      std::shared_ptr< monty::ndarray< int32_t,1 > > shape{};
      std::string name{};
      std::shared_ptr< monty::ndarray< int64_t,1 > > slackidxs{};
      std::shared_ptr< monty::ndarray< int32_t,1 > > nativeidxs{};
      int32_t conid{};

      virtual void destroy();

      static PSDConstraint::t _new_PSDConstraint(monty::rc_ptr< ::mosek::fusion::PSDConstraint > _2885_c,monty::rc_ptr< ::mosek::fusion::Model > _2886_m);
      void _initialize(monty::rc_ptr< ::mosek::fusion::PSDConstraint > _2885_c,monty::rc_ptr< ::mosek::fusion::Model > _2886_m);
      static PSDConstraint::t _new_PSDConstraint(monty::rc_ptr< ::mosek::fusion::Model > _2887_model,const std::string &  _2888_name,int32_t _2889_conid,std::shared_ptr< monty::ndarray< int32_t,1 > > _2890_shape,int32_t _2891_conedim0,int32_t _2892_conedim1,std::shared_ptr< monty::ndarray< int64_t,1 > > _2893_slackidxs,std::shared_ptr< monty::ndarray< int32_t,1 > > _2894_nativeidxs);
      void _initialize(monty::rc_ptr< ::mosek::fusion::Model > _2887_model,const std::string &  _2888_name,int32_t _2889_conid,std::shared_ptr< monty::ndarray< int32_t,1 > > _2890_shape,int32_t _2891_conedim0,int32_t _2892_conedim1,std::shared_ptr< monty::ndarray< int64_t,1 > > _2893_slackidxs,std::shared_ptr< monty::ndarray< int32_t,1 > > _2894_nativeidxs);
      virtual /* override */ std::string toString() ;
      virtual /* override */ monty::rc_ptr< ::mosek::fusion::ModelConstraint > __mosek_2fusion_2PSDConstraint__clone(monty::rc_ptr< ::mosek::fusion::Model > _2895_m) ;
      virtual monty::rc_ptr< ::mosek::fusion::ModelConstraint > __mosek_2fusion_2ModelConstraint__clone(monty::rc_ptr< ::mosek::fusion::Model > _2895_m) { return __mosek_2fusion_2PSDConstraint__clone(_2895_m); }
      static  std::shared_ptr< monty::ndarray< int32_t,1 > > computenidxs(std::shared_ptr< monty::ndarray< int32_t,1 > > _2896_shape,int32_t _2897_d0,int32_t _2898_d1,std::shared_ptr< monty::ndarray< int32_t,1 > > _2899_nativeidxs);
    }; // struct PSDConstraint;

    struct p_RangedConstraint : public ::mosek::fusion::p_ModelConstraint
    {
      RangedConstraint * _pubthis;
      static mosek::fusion::p_RangedConstraint* _get_impl(mosek::fusion::RangedConstraint * _inst){ return static_cast< mosek::fusion::p_RangedConstraint* >(mosek::fusion::p_ModelConstraint::_get_impl(_inst)); }
      static mosek::fusion::p_RangedConstraint * _get_impl(mosek::fusion::RangedConstraint::t _inst) { return _get_impl(_inst.get()); }
      p_RangedConstraint(RangedConstraint * _pubthis);
      virtual ~p_RangedConstraint() { /* std::cout << "~p_RangedConstraint" << std::endl;*/ };
      std::shared_ptr< monty::ndarray< int32_t,1 > > nativeidxs{};
      std::shared_ptr< monty::ndarray< int32_t,1 > > shape{};

      virtual void destroy();

      static RangedConstraint::t _new_RangedConstraint(monty::rc_ptr< ::mosek::fusion::RangedConstraint > _2929_c,monty::rc_ptr< ::mosek::fusion::Model > _2930_m);
      void _initialize(monty::rc_ptr< ::mosek::fusion::RangedConstraint > _2929_c,monty::rc_ptr< ::mosek::fusion::Model > _2930_m);
      static RangedConstraint::t _new_RangedConstraint(monty::rc_ptr< ::mosek::fusion::Model > _2931_model,const std::string &  _2932_name,std::shared_ptr< monty::ndarray< int32_t,1 > > _2933_shape,std::shared_ptr< monty::ndarray< int32_t,1 > > _2934_nativeidxs,int32_t _2935_conid);
      void _initialize(monty::rc_ptr< ::mosek::fusion::Model > _2931_model,const std::string &  _2932_name,std::shared_ptr< monty::ndarray< int32_t,1 > > _2933_shape,std::shared_ptr< monty::ndarray< int32_t,1 > > _2934_nativeidxs,int32_t _2935_conid);
      virtual monty::rc_ptr< ::mosek::fusion::BoundInterfaceConstraint > __mosek_2fusion_2RangedConstraint__upperBoundCon() ;
      virtual monty::rc_ptr< ::mosek::fusion::BoundInterfaceConstraint > __mosek_2fusion_2RangedConstraint__lowerBoundCon() ;
      virtual /* override */ monty::rc_ptr< ::mosek::fusion::ModelConstraint > __mosek_2fusion_2RangedConstraint__clone(monty::rc_ptr< ::mosek::fusion::Model > _2936_m) ;
      virtual monty::rc_ptr< ::mosek::fusion::ModelConstraint > __mosek_2fusion_2ModelConstraint__clone(monty::rc_ptr< ::mosek::fusion::Model > _2936_m) { return __mosek_2fusion_2RangedConstraint__clone(_2936_m); }
    }; // struct RangedConstraint;

    struct p_ConicConstraint : public ::mosek::fusion::p_ModelConstraint
    {
      ConicConstraint * _pubthis;
      static mosek::fusion::p_ConicConstraint* _get_impl(mosek::fusion::ConicConstraint * _inst){ return static_cast< mosek::fusion::p_ConicConstraint* >(mosek::fusion::p_ModelConstraint::_get_impl(_inst)); }
      static mosek::fusion::p_ConicConstraint * _get_impl(mosek::fusion::ConicConstraint::t _inst) { return _get_impl(_inst.get()); }
      p_ConicConstraint(ConicConstraint * _pubthis);
      virtual ~p_ConicConstraint() { /* std::cout << "~p_ConicConstraint" << std::endl;*/ };
      std::shared_ptr< monty::ndarray< std::shared_ptr< monty::ndarray< std::string,1 > >,1 > > indexnames{};
      std::shared_ptr< monty::ndarray< int32_t,1 > > nativeidxs{};
      bool names_flushed{};
      std::string name{};
      std::shared_ptr< monty::ndarray< int32_t,1 > > shape{};
      monty::rc_ptr< ::mosek::fusion::ConeDomain > dom{};
      int32_t conid{};

      virtual void destroy();

      static ConicConstraint::t _new_ConicConstraint(monty::rc_ptr< ::mosek::fusion::ConicConstraint > _2937_c,monty::rc_ptr< ::mosek::fusion::Model > _2938_m);
      void _initialize(monty::rc_ptr< ::mosek::fusion::ConicConstraint > _2937_c,monty::rc_ptr< ::mosek::fusion::Model > _2938_m);
      static ConicConstraint::t _new_ConicConstraint(monty::rc_ptr< ::mosek::fusion::Model > _2939_model,const std::string &  _2940_name,monty::rc_ptr< ::mosek::fusion::ConeDomain > _2941_dom,std::shared_ptr< monty::ndarray< int32_t,1 > > _2942_shape,int32_t _2943_conid,std::shared_ptr< monty::ndarray< int32_t,1 > > _2944_nativeidxs,std::shared_ptr< monty::ndarray< std::shared_ptr< monty::ndarray< std::string,1 > >,1 > > _2945_indexnames);
      void _initialize(monty::rc_ptr< ::mosek::fusion::Model > _2939_model,const std::string &  _2940_name,monty::rc_ptr< ::mosek::fusion::ConeDomain > _2941_dom,std::shared_ptr< monty::ndarray< int32_t,1 > > _2942_shape,int32_t _2943_conid,std::shared_ptr< monty::ndarray< int32_t,1 > > _2944_nativeidxs,std::shared_ptr< monty::ndarray< std::shared_ptr< monty::ndarray< std::string,1 > >,1 > > _2945_indexnames);
      virtual /* override */ std::string toString() ;
      virtual void domainToString(int64_t _2948_i,monty::rc_ptr< ::mosek::fusion::Utils::StringBuffer > _2949_sb) ;
      virtual /* override */ monty::rc_ptr< ::mosek::fusion::ModelConstraint > __mosek_2fusion_2ConicConstraint__clone(monty::rc_ptr< ::mosek::fusion::Model > _2950_m) ;
      virtual monty::rc_ptr< ::mosek::fusion::ModelConstraint > __mosek_2fusion_2ModelConstraint__clone(monty::rc_ptr< ::mosek::fusion::Model > _2950_m) { return __mosek_2fusion_2ConicConstraint__clone(_2950_m); }
    }; // struct ConicConstraint;

    struct p_LinearConstraint : public ::mosek::fusion::p_ModelConstraint
    {
      LinearConstraint * _pubthis;
      static mosek::fusion::p_LinearConstraint* _get_impl(mosek::fusion::LinearConstraint * _inst){ return static_cast< mosek::fusion::p_LinearConstraint* >(mosek::fusion::p_ModelConstraint::_get_impl(_inst)); }
      static mosek::fusion::p_LinearConstraint * _get_impl(mosek::fusion::LinearConstraint::t _inst) { return _get_impl(_inst.get()); }
      p_LinearConstraint(LinearConstraint * _pubthis);
      virtual ~p_LinearConstraint() { /* std::cout << "~p_LinearConstraint" << std::endl;*/ };
      std::shared_ptr< monty::ndarray< std::shared_ptr< monty::ndarray< std::string,1 > >,1 > > indexnames{};
      bool names_flushed{};
      std::shared_ptr< monty::ndarray< int32_t,1 > > nidxs{};
      std::string name{};
      int32_t conid{};

      virtual void destroy();

      static LinearConstraint::t _new_LinearConstraint(monty::rc_ptr< ::mosek::fusion::LinearConstraint > _2951_c,monty::rc_ptr< ::mosek::fusion::Model > _2952_m);
      void _initialize(monty::rc_ptr< ::mosek::fusion::LinearConstraint > _2951_c,monty::rc_ptr< ::mosek::fusion::Model > _2952_m);
      static LinearConstraint::t _new_LinearConstraint(monty::rc_ptr< ::mosek::fusion::Model > _2953_model,const std::string &  _2954_name,int32_t _2955_conid,std::shared_ptr< monty::ndarray< int32_t,1 > > _2956_shape,std::shared_ptr< monty::ndarray< int32_t,1 > > _2957_nidxs,std::shared_ptr< monty::ndarray< std::shared_ptr< monty::ndarray< std::string,1 > >,1 > > _2958_indexnames);
      void _initialize(monty::rc_ptr< ::mosek::fusion::Model > _2953_model,const std::string &  _2954_name,int32_t _2955_conid,std::shared_ptr< monty::ndarray< int32_t,1 > > _2956_shape,std::shared_ptr< monty::ndarray< int32_t,1 > > _2957_nidxs,std::shared_ptr< monty::ndarray< std::shared_ptr< monty::ndarray< std::string,1 > >,1 > > _2958_indexnames);
      virtual /* override */ std::string toString() ;
      virtual /* override */ void flushNames() ;
      virtual void domainToString(int64_t _2960_i,monty::rc_ptr< ::mosek::fusion::Utils::StringBuffer > _2961_sb) ;
      virtual /* override */ monty::rc_ptr< ::mosek::fusion::ModelConstraint > __mosek_2fusion_2LinearConstraint__clone(monty::rc_ptr< ::mosek::fusion::Model > _2962_m) ;
      virtual monty::rc_ptr< ::mosek::fusion::ModelConstraint > __mosek_2fusion_2ModelConstraint__clone(monty::rc_ptr< ::mosek::fusion::Model > _2962_m) { return __mosek_2fusion_2LinearConstraint__clone(_2962_m); }
    }; // struct LinearConstraint;

    struct p_Set
    {
      Set * _pubthis;
      static mosek::fusion::p_Set* _get_impl(mosek::fusion::Set * _inst){ assert(_inst); assert(_inst->_impl); return _inst->_impl; }
      static mosek::fusion::p_Set * _get_impl(mosek::fusion::Set::t _inst) { return _get_impl(_inst.get()); }
      p_Set(Set * _pubthis);
      virtual ~p_Set() { /* std::cout << "~p_Set" << std::endl;*/ };

      virtual void destroy();

      static  int64_t size(std::shared_ptr< monty::ndarray< int32_t,1 > > _3126_shape);
      static  bool match(std::shared_ptr< monty::ndarray< int32_t,1 > > _3129_s1,std::shared_ptr< monty::ndarray< int32_t,1 > > _3130_s2);
      static  int64_t linearidx(std::shared_ptr< monty::ndarray< int32_t,1 > > _3132_shape,std::shared_ptr< monty::ndarray< int32_t,1 > > _3133_key);
      static  std::shared_ptr< monty::ndarray< int32_t,1 > > idxtokey(std::shared_ptr< monty::ndarray< int32_t,1 > > _3136_shape,int64_t _3137_idx);
      static  void idxtokey(std::shared_ptr< monty::ndarray< int32_t,1 > > _3139_shape,int64_t _3140_idx,std::shared_ptr< monty::ndarray< int32_t,1 > > _3141_dest);
      static  std::string indexToString(std::shared_ptr< monty::ndarray< int32_t,1 > > _3145_shape,int64_t _3146_key);
      static  std::string keyToString(std::shared_ptr< monty::ndarray< int32_t,1 > > _3153_key);
      static  void indexToKey(std::shared_ptr< monty::ndarray< int32_t,1 > > _3156_shape,int64_t _3157_key,std::shared_ptr< monty::ndarray< int32_t,1 > > _3158_res);
      static  std::shared_ptr< monty::ndarray< int64_t,1 > > strides(std::shared_ptr< monty::ndarray< int32_t,1 > > _3162_shape);
      static  std::shared_ptr< monty::ndarray< int32_t,1 > > make(std::shared_ptr< monty::ndarray< int32_t,1 > > _3166_set1,std::shared_ptr< monty::ndarray< int32_t,1 > > _3167_set2);
      static  std::shared_ptr< monty::ndarray< int32_t,1 > > make(std::shared_ptr< monty::ndarray< int32_t,1 > > _3171_sizes);
      static  std::shared_ptr< monty::ndarray< int32_t,1 > > make(int32_t _3173_s1,int32_t _3174_s2,int32_t _3175_s3);
      static  std::shared_ptr< monty::ndarray< int32_t,1 > > make(int32_t _3176_s1,int32_t _3177_s2);
      static  std::shared_ptr< monty::ndarray< int32_t,1 > > make(int32_t _3178_sz);
      static  std::shared_ptr< monty::ndarray< int32_t,1 > > scalar();
      static  std::shared_ptr< monty::ndarray< int32_t,1 > > make(std::shared_ptr< monty::ndarray< std::string,1 > > _3179_names);
    }; // struct Set;

    struct p_ConeDomain
    {
      ConeDomain * _pubthis;
      static mosek::fusion::p_ConeDomain* _get_impl(mosek::fusion::ConeDomain * _inst){ assert(_inst); assert(_inst->_impl); return _inst->_impl; }
      static mosek::fusion::p_ConeDomain * _get_impl(mosek::fusion::ConeDomain::t _inst) { return _get_impl(_inst.get()); }
      p_ConeDomain(ConeDomain * _pubthis);
      virtual ~p_ConeDomain() { /* std::cout << "~p_ConeDomain" << std::endl;*/ };
      std::shared_ptr< monty::ndarray< std::shared_ptr< monty::ndarray< std::string,1 > >,1 > > indexnames{};
      int64_t domsize{};
      std::shared_ptr< monty::ndarray< double,1 > > domofs{};
      std::shared_ptr< monty::ndarray< double,1 > > alpha{};
      std::shared_ptr< monty::ndarray< int32_t,1 > > shape{};
      bool int_flag{};
      bool axisset{};
      int32_t axisidx{};
      mosek::fusion::QConeKey key{};

      virtual void destroy();

      static ConeDomain::t _new_ConeDomain(mosek::fusion::QConeKey _3180_k,std::shared_ptr< monty::ndarray< double,1 > > _3181_alpha,std::shared_ptr< monty::ndarray< int32_t,1 > > _3182_d);
      void _initialize(mosek::fusion::QConeKey _3180_k,std::shared_ptr< monty::ndarray< double,1 > > _3181_alpha,std::shared_ptr< monty::ndarray< int32_t,1 > > _3182_d);
      static ConeDomain::t _new_ConeDomain(mosek::fusion::QConeKey _3183_k,std::shared_ptr< monty::ndarray< int32_t,1 > > _3184_d);
      void _initialize(mosek::fusion::QConeKey _3183_k,std::shared_ptr< monty::ndarray< int32_t,1 > > _3184_d);
      static ConeDomain::t _new_ConeDomain(monty::rc_ptr< ::mosek::fusion::ConeDomain > _3185_other);
      void _initialize(monty::rc_ptr< ::mosek::fusion::ConeDomain > _3185_other);
      virtual bool match_shape(std::shared_ptr< monty::ndarray< int32_t,1 > > _3186_shp) ;
      virtual monty::rc_ptr< ::mosek::fusion::ConeDomain > __mosek_2fusion_2ConeDomain__integral() ;
      virtual bool axisIsSet() ;
      virtual int32_t getAxis() ;
      virtual monty::rc_ptr< ::mosek::fusion::ConeDomain > __mosek_2fusion_2ConeDomain__axis(int32_t _3187_a) ;
      virtual monty::rc_ptr< ::mosek::fusion::ConeDomain > __mosek_2fusion_2ConeDomain__withShape(int32_t _3188_dim0,int32_t _3189_dim1,int32_t _3190_dim2) ;
      virtual monty::rc_ptr< ::mosek::fusion::ConeDomain > __mosek_2fusion_2ConeDomain__withShape(int32_t _3191_dim0,int32_t _3192_dim1) ;
      virtual monty::rc_ptr< ::mosek::fusion::ConeDomain > __mosek_2fusion_2ConeDomain__withShape(int32_t _3193_dim0) ;
      virtual monty::rc_ptr< ::mosek::fusion::ConeDomain > __mosek_2fusion_2ConeDomain__withShape(std::shared_ptr< monty::ndarray< int32_t,1 > > _3194_shp) ;
      virtual monty::rc_ptr< ::mosek::fusion::ConeDomain > __mosek_2fusion_2ConeDomain__withShape_(std::shared_ptr< monty::ndarray< int32_t,1 > > _3195_shp) ;
      virtual monty::rc_ptr< ::mosek::fusion::ConeDomain > __mosek_2fusion_2ConeDomain__withNamesOnAxis(std::shared_ptr< monty::ndarray< std::string,1 > > _3196_names,int32_t _3197_axis) ;
      virtual void finalize_and_validate_inplace(std::shared_ptr< monty::ndarray< int32_t,1 > > _3202_shp) ;
      virtual monty::rc_ptr< ::mosek::fusion::ConeDomain > __mosek_2fusion_2ConeDomain__finalize_and_validate(std::shared_ptr< monty::ndarray< int32_t,1 > > _3206_shp) ;
    }; // struct ConeDomain;

    struct p_PSDDomain
    {
      PSDDomain * _pubthis;
      static mosek::fusion::p_PSDDomain* _get_impl(mosek::fusion::PSDDomain * _inst){ assert(_inst); assert(_inst->_impl); return _inst->_impl; }
      static mosek::fusion::p_PSDDomain * _get_impl(mosek::fusion::PSDDomain::t _inst) { return _get_impl(_inst.get()); }
      p_PSDDomain(PSDDomain * _pubthis);
      virtual ~p_PSDDomain() { /* std::cout << "~p_PSDDomain" << std::endl;*/ };
      std::shared_ptr< monty::ndarray< std::shared_ptr< monty::ndarray< std::string,1 > >,1 > > indexnames{};
      bool axisIsSet{};
      int32_t conedim2{};
      int32_t conedim1{};
      mosek::fusion::PSDKey key{};
      std::shared_ptr< monty::ndarray< int32_t,1 > > shape{};

      virtual void destroy();

      static PSDDomain::t _new_PSDDomain(mosek::fusion::PSDKey _3208_k,std::shared_ptr< monty::ndarray< int32_t,1 > > _3209_shp,int32_t _3210_conedim1,int32_t _3211_conedim2);
      void _initialize(mosek::fusion::PSDKey _3208_k,std::shared_ptr< monty::ndarray< int32_t,1 > > _3209_shp,int32_t _3210_conedim1,int32_t _3211_conedim2);
      static PSDDomain::t _new_PSDDomain(mosek::fusion::PSDKey _3213_k,std::shared_ptr< monty::ndarray< int32_t,1 > > _3214_shp);
      void _initialize(mosek::fusion::PSDKey _3213_k,std::shared_ptr< monty::ndarray< int32_t,1 > > _3214_shp);
      static PSDDomain::t _new_PSDDomain(mosek::fusion::PSDKey _3215_k);
      void _initialize(mosek::fusion::PSDKey _3215_k);
      static PSDDomain::t _new_PSDDomain(monty::rc_ptr< ::mosek::fusion::PSDDomain > _3216_other);
      void _initialize(monty::rc_ptr< ::mosek::fusion::PSDDomain > _3216_other);
      virtual monty::rc_ptr< ::mosek::fusion::PSDDomain > __mosek_2fusion_2PSDDomain__axis(int32_t _3217_conedim1,int32_t _3218_conedim2) ;
      virtual monty::rc_ptr< ::mosek::fusion::PSDDomain > __mosek_2fusion_2PSDDomain__withNamesOnAxis(std::shared_ptr< monty::ndarray< std::string,1 > > _3219_names,int32_t _3220_axis) ;
      virtual void finalize_and_validate_inplace(std::shared_ptr< monty::ndarray< int32_t,1 > > _3227_shp) ;
      virtual monty::rc_ptr< ::mosek::fusion::PSDDomain > __mosek_2fusion_2PSDDomain__finalize_and_validate(std::shared_ptr< monty::ndarray< int32_t,1 > > _3230_shp) ;
    }; // struct PSDDomain;

    struct p_RangeDomain
    {
      RangeDomain * _pubthis;
      static mosek::fusion::p_RangeDomain* _get_impl(mosek::fusion::RangeDomain * _inst){ assert(_inst); assert(_inst->_impl); return _inst->_impl; }
      static mosek::fusion::p_RangeDomain * _get_impl(mosek::fusion::RangeDomain::t _inst) { return _get_impl(_inst.get()); }
      p_RangeDomain(RangeDomain * _pubthis);
      virtual ~p_RangeDomain() { /* std::cout << "~p_RangeDomain" << std::endl;*/ };
      int64_t domsize{};
      int64_t nelements{};
      std::shared_ptr< monty::ndarray< std::shared_ptr< monty::ndarray< std::string,1 > >,1 > > indexnames{};
      bool cardinal_flag{};
      bool scalable{};
      std::shared_ptr< monty::ndarray< double,1 > > ub{};
      std::shared_ptr< monty::ndarray< double,1 > > lb{};
      std::shared_ptr< monty::ndarray< int32_t,2 > > sparsity{};
      bool empty{};
      std::shared_ptr< monty::ndarray< int32_t,1 > > shape{};

      virtual void destroy();

      static RangeDomain::t _new_RangeDomain(bool _3233_scalable,std::shared_ptr< monty::ndarray< double,1 > > _3234_lb,std::shared_ptr< monty::ndarray< double,1 > > _3235_ub,std::shared_ptr< monty::ndarray< int32_t,1 > > _3236_dims);
      void _initialize(bool _3233_scalable,std::shared_ptr< monty::ndarray< double,1 > > _3234_lb,std::shared_ptr< monty::ndarray< double,1 > > _3235_ub,std::shared_ptr< monty::ndarray< int32_t,1 > > _3236_dims);
      static RangeDomain::t _new_RangeDomain(bool _3237_scalable,double _3238_lb,double _3239_ub,std::shared_ptr< monty::ndarray< int32_t,1 > > _3240_dims);
      void _initialize(bool _3237_scalable,double _3238_lb,double _3239_ub,std::shared_ptr< monty::ndarray< int32_t,1 > > _3240_dims);
      static RangeDomain::t _new_RangeDomain(bool _3245_scalable,std::shared_ptr< monty::ndarray< double,1 > > _3246_lb,std::shared_ptr< monty::ndarray< double,1 > > _3247_ub,std::shared_ptr< monty::ndarray< int32_t,1 > > _3248_dims,std::shared_ptr< monty::ndarray< int32_t,2 > > _3249_sp);
      void _initialize(bool _3245_scalable,std::shared_ptr< monty::ndarray< double,1 > > _3246_lb,std::shared_ptr< monty::ndarray< double,1 > > _3247_ub,std::shared_ptr< monty::ndarray< int32_t,1 > > _3248_dims,std::shared_ptr< monty::ndarray< int32_t,2 > > _3249_sp);
      static RangeDomain::t _new_RangeDomain(bool _3250_scalable,std::shared_ptr< monty::ndarray< double,1 > > _3251_lb,std::shared_ptr< monty::ndarray< double,1 > > _3252_ub,std::shared_ptr< monty::ndarray< int32_t,1 > > _3253_dims,std::shared_ptr< monty::ndarray< int32_t,2 > > _3254_sp,int32_t _3255_steal);
      void _initialize(bool _3250_scalable,std::shared_ptr< monty::ndarray< double,1 > > _3251_lb,std::shared_ptr< monty::ndarray< double,1 > > _3252_ub,std::shared_ptr< monty::ndarray< int32_t,1 > > _3253_dims,std::shared_ptr< monty::ndarray< int32_t,2 > > _3254_sp,int32_t _3255_steal);
      static RangeDomain::t _new_RangeDomain(monty::rc_ptr< ::mosek::fusion::RangeDomain > _3256_other);
      void _initialize(monty::rc_ptr< ::mosek::fusion::RangeDomain > _3256_other);
      virtual monty::rc_ptr< ::mosek::fusion::SymmetricRangeDomain > __mosek_2fusion_2RangeDomain__symmetric() ;
      virtual monty::rc_ptr< ::mosek::fusion::RangeDomain > __mosek_2fusion_2RangeDomain__sparse(std::shared_ptr< monty::ndarray< int32_t,2 > > _3258_sparsity) ;
      virtual monty::rc_ptr< ::mosek::fusion::RangeDomain > __mosek_2fusion_2RangeDomain__sparse(std::shared_ptr< monty::ndarray< int32_t,1 > > _3261_sparsity) ;
      virtual monty::rc_ptr< ::mosek::fusion::RangeDomain > __mosek_2fusion_2RangeDomain__sparse() ;
      virtual monty::rc_ptr< ::mosek::fusion::RangeDomain > __mosek_2fusion_2RangeDomain__integral() ;
      virtual monty::rc_ptr< ::mosek::fusion::RangeDomain > __mosek_2fusion_2RangeDomain__withShape(int32_t _3263_dim0,int32_t _3264_dim1,int32_t _3265_dim2) ;
      virtual monty::rc_ptr< ::mosek::fusion::RangeDomain > __mosek_2fusion_2RangeDomain__withShape(int32_t _3266_dim0,int32_t _3267_dim1) ;
      virtual monty::rc_ptr< ::mosek::fusion::RangeDomain > __mosek_2fusion_2RangeDomain__withShape(int32_t _3268_dim0) ;
      virtual monty::rc_ptr< ::mosek::fusion::RangeDomain > __mosek_2fusion_2RangeDomain__withShape(std::shared_ptr< monty::ndarray< int32_t,1 > > _3269_shp) ;
      virtual monty::rc_ptr< ::mosek::fusion::RangeDomain > __mosek_2fusion_2RangeDomain__withNamesOnAxis(std::shared_ptr< monty::ndarray< std::string,1 > > _3270_names,int32_t _3271_axis) ;
      virtual bool match_shape(std::shared_ptr< monty::ndarray< int32_t,1 > > _3278_shp) ;
      virtual void finalize_and_validate_inplace(std::shared_ptr< monty::ndarray< int32_t,1 > > _3280_shp) ;
      virtual monty::rc_ptr< ::mosek::fusion::RangeDomain > __mosek_2fusion_2RangeDomain__finalize_and_validate(std::shared_ptr< monty::ndarray< int32_t,1 > > _3287_shp) ;
    }; // struct RangeDomain;

    struct p_SymmetricRangeDomain : public ::mosek::fusion::p_RangeDomain
    {
      SymmetricRangeDomain * _pubthis;
      static mosek::fusion::p_SymmetricRangeDomain* _get_impl(mosek::fusion::SymmetricRangeDomain * _inst){ return static_cast< mosek::fusion::p_SymmetricRangeDomain* >(mosek::fusion::p_RangeDomain::_get_impl(_inst)); }
      static mosek::fusion::p_SymmetricRangeDomain * _get_impl(mosek::fusion::SymmetricRangeDomain::t _inst) { return _get_impl(_inst.get()); }
      p_SymmetricRangeDomain(SymmetricRangeDomain * _pubthis);
      virtual ~p_SymmetricRangeDomain() { /* std::cout << "~p_SymmetricRangeDomain" << std::endl;*/ };
      int32_t dim{};

      virtual void destroy();

      static SymmetricRangeDomain::t _new_SymmetricRangeDomain(monty::rc_ptr< ::mosek::fusion::RangeDomain > _3232_other);
      void _initialize(monty::rc_ptr< ::mosek::fusion::RangeDomain > _3232_other);
    }; // struct SymmetricRangeDomain;

    struct p_SymmetricLinearDomain
    {
      SymmetricLinearDomain * _pubthis;
      static mosek::fusion::p_SymmetricLinearDomain* _get_impl(mosek::fusion::SymmetricLinearDomain * _inst){ assert(_inst); assert(_inst->_impl); return _inst->_impl; }
      static mosek::fusion::p_SymmetricLinearDomain * _get_impl(mosek::fusion::SymmetricLinearDomain::t _inst) { return _get_impl(_inst.get()); }
      p_SymmetricLinearDomain(SymmetricLinearDomain * _pubthis);
      virtual ~p_SymmetricLinearDomain() { /* std::cout << "~p_SymmetricLinearDomain" << std::endl;*/ };
      std::shared_ptr< monty::ndarray< int32_t,2 > > sparsity{};
      bool cardinal_flag{};
      mosek::fusion::RelationKey key{};
      std::shared_ptr< monty::ndarray< int32_t,1 > > shape{};
      monty::rc_ptr< ::mosek::fusion::LinearDomain > dom{};
      int32_t dim{};

      virtual void destroy();

      static SymmetricLinearDomain::t _new_SymmetricLinearDomain(monty::rc_ptr< ::mosek::fusion::LinearDomain > _3289_other);
      void _initialize(monty::rc_ptr< ::mosek::fusion::LinearDomain > _3289_other);
      virtual monty::rc_ptr< ::mosek::fusion::SymmetricLinearDomain > __mosek_2fusion_2SymmetricLinearDomain__sparse(std::shared_ptr< monty::ndarray< int32_t,2 > > _3290_sparsity) ;
      virtual monty::rc_ptr< ::mosek::fusion::SymmetricLinearDomain > __mosek_2fusion_2SymmetricLinearDomain__sparse(std::shared_ptr< monty::ndarray< int32_t,1 > > _3293_sparsity) ;
      virtual monty::rc_ptr< ::mosek::fusion::SymmetricLinearDomain > __mosek_2fusion_2SymmetricLinearDomain__integral() ;
      virtual bool match_shape(std::shared_ptr< monty::ndarray< int32_t,1 > > _3295_shp) ;
    }; // struct SymmetricLinearDomain;

    struct p_LinearDomain
    {
      LinearDomain * _pubthis;
      static mosek::fusion::p_LinearDomain* _get_impl(mosek::fusion::LinearDomain * _inst){ assert(_inst); assert(_inst->_impl); return _inst->_impl; }
      static mosek::fusion::p_LinearDomain * _get_impl(mosek::fusion::LinearDomain::t _inst) { return _get_impl(_inst.get()); }
      p_LinearDomain(LinearDomain * _pubthis);
      virtual ~p_LinearDomain() { /* std::cout << "~p_LinearDomain" << std::endl;*/ };
      int64_t nelements{};
      int64_t domsize{};
      std::shared_ptr< monty::ndarray< std::shared_ptr< monty::ndarray< std::string,1 > >,1 > > indexnames{};
      bool empty{};
      bool scalable{};
      std::shared_ptr< monty::ndarray< int32_t,2 > > sparsity{};
      bool cardinal_flag{};
      mosek::fusion::RelationKey key{};
      std::shared_ptr< monty::ndarray< double,1 > > bnd{};
      std::shared_ptr< monty::ndarray< int32_t,1 > > shape{};

      virtual void destroy();

      static LinearDomain::t _new_LinearDomain(mosek::fusion::RelationKey _3297_k,bool _3298_scalable,std::shared_ptr< monty::ndarray< double,1 > > _3299_rhs,std::shared_ptr< monty::ndarray< int32_t,1 > > _3300_dims);
      void _initialize(mosek::fusion::RelationKey _3297_k,bool _3298_scalable,std::shared_ptr< monty::ndarray< double,1 > > _3299_rhs,std::shared_ptr< monty::ndarray< int32_t,1 > > _3300_dims);
      static LinearDomain::t _new_LinearDomain(mosek::fusion::RelationKey _3301_k,bool _3302_scalable,std::shared_ptr< monty::ndarray< double,1 > > _3303_rhs,std::shared_ptr< monty::ndarray< int32_t,1 > > _3304_dims,std::shared_ptr< monty::ndarray< int32_t,2 > > _3305_sp,int32_t _3306_steal);
      void _initialize(mosek::fusion::RelationKey _3301_k,bool _3302_scalable,std::shared_ptr< monty::ndarray< double,1 > > _3303_rhs,std::shared_ptr< monty::ndarray< int32_t,1 > > _3304_dims,std::shared_ptr< monty::ndarray< int32_t,2 > > _3305_sp,int32_t _3306_steal);
      static LinearDomain::t _new_LinearDomain(monty::rc_ptr< ::mosek::fusion::LinearDomain > _3307_other);
      void _initialize(monty::rc_ptr< ::mosek::fusion::LinearDomain > _3307_other);
      virtual monty::rc_ptr< ::mosek::fusion::SymmetricLinearDomain > __mosek_2fusion_2LinearDomain__symmetric() ;
      virtual monty::rc_ptr< ::mosek::fusion::LinearDomain > __mosek_2fusion_2LinearDomain__sparse(std::shared_ptr< monty::ndarray< int32_t,2 > > _3308_sparsity) ;
      virtual monty::rc_ptr< ::mosek::fusion::LinearDomain > __mosek_2fusion_2LinearDomain__sparse(std::shared_ptr< monty::ndarray< int32_t,1 > > _3311_sparsity) ;
      virtual monty::rc_ptr< ::mosek::fusion::LinearDomain > __mosek_2fusion_2LinearDomain__sparse() ;
      virtual monty::rc_ptr< ::mosek::fusion::LinearDomain > __mosek_2fusion_2LinearDomain__integral() ;
      virtual monty::rc_ptr< ::mosek::fusion::LinearDomain > __mosek_2fusion_2LinearDomain__withShape(int32_t _3313_dim0,int32_t _3314_dim1,int32_t _3315_dim2) ;
      virtual monty::rc_ptr< ::mosek::fusion::LinearDomain > __mosek_2fusion_2LinearDomain__withShape(int32_t _3316_dim0,int32_t _3317_dim1) ;
      virtual monty::rc_ptr< ::mosek::fusion::LinearDomain > __mosek_2fusion_2LinearDomain__withShape(int32_t _3318_dim0) ;
      virtual monty::rc_ptr< ::mosek::fusion::LinearDomain > __mosek_2fusion_2LinearDomain__withShape(std::shared_ptr< monty::ndarray< int32_t,1 > > _3319_shp) ;
      virtual monty::rc_ptr< ::mosek::fusion::LinearDomain > __mosek_2fusion_2LinearDomain__withNamesOnAxis(std::shared_ptr< monty::ndarray< std::string,1 > > _3320_names,int32_t _3321_axis) ;
      virtual bool match_shape(std::shared_ptr< monty::ndarray< int32_t,1 > > _3328_shp) ;
      virtual void finalize_and_validate_inplace(std::shared_ptr< monty::ndarray< int32_t,1 > > _3330_shp) ;
      virtual monty::rc_ptr< ::mosek::fusion::LinearDomain > __mosek_2fusion_2LinearDomain__finalize_and_validate(std::shared_ptr< monty::ndarray< int32_t,1 > > _3338_shp) ;
    }; // struct LinearDomain;

    struct p_Domain
    {
      Domain * _pubthis;
      static mosek::fusion::p_Domain* _get_impl(mosek::fusion::Domain * _inst){ assert(_inst); assert(_inst->_impl); return _inst->_impl; }
      static mosek::fusion::p_Domain * _get_impl(mosek::fusion::Domain::t _inst) { return _get_impl(_inst.get()); }
      p_Domain(Domain * _pubthis);
      virtual ~p_Domain() { /* std::cout << "~p_Domain" << std::endl;*/ };

      virtual void destroy();

      static  int64_t dimsize(std::shared_ptr< monty::ndarray< int32_t,1 > > _3340_dims);
      static  monty::rc_ptr< ::mosek::fusion::RangeDomain > mkRangedDomain(monty::rc_ptr< ::mosek::fusion::Matrix > _3343_lb,monty::rc_ptr< ::mosek::fusion::Matrix > _3344_ub);
      static  monty::rc_ptr< ::mosek::fusion::RangeDomain > mkRangedDomain(std::shared_ptr< monty::ndarray< double,2 > > _3373_lb,std::shared_ptr< monty::ndarray< double,2 > > _3374_ub);
      static  monty::rc_ptr< ::mosek::fusion::LinearDomain > mkLinearDomain(mosek::fusion::RelationKey _3383_k,monty::rc_ptr< ::mosek::fusion::Matrix > _3384_mx);
      static  int64_t prod(std::shared_ptr< monty::ndarray< int32_t,1 > > _3390_dim);
      static  monty::rc_ptr< ::mosek::fusion::RangeDomain > inRange(bool _3393_scalable,std::shared_ptr< monty::ndarray< double,1 > > _3394_lb,std::shared_ptr< monty::ndarray< double,1 > > _3395_ub,std::shared_ptr< monty::ndarray< int32_t,2 > > _3396_sp,std::shared_ptr< monty::ndarray< int32_t,1 > > _3397_dims);
      static  monty::rc_ptr< ::mosek::fusion::SymmetricRangeDomain > symmetric(monty::rc_ptr< ::mosek::fusion::RangeDomain > _3399_rd);
      static  monty::rc_ptr< ::mosek::fusion::SymmetricLinearDomain > symmetric(monty::rc_ptr< ::mosek::fusion::LinearDomain > _3400_ld);
      static  monty::rc_ptr< ::mosek::fusion::RangeDomain > sparse(monty::rc_ptr< ::mosek::fusion::RangeDomain > _3401_rd,std::shared_ptr< monty::ndarray< int32_t,2 > > _3402_sparsity);
      static  monty::rc_ptr< ::mosek::fusion::RangeDomain > sparse(monty::rc_ptr< ::mosek::fusion::RangeDomain > _3403_rd,std::shared_ptr< monty::ndarray< int32_t,1 > > _3404_sparsity);
      static  monty::rc_ptr< ::mosek::fusion::LinearDomain > sparse(monty::rc_ptr< ::mosek::fusion::LinearDomain > _3405_ld,std::shared_ptr< monty::ndarray< int32_t,2 > > _3406_sparsity);
      static  monty::rc_ptr< ::mosek::fusion::LinearDomain > sparse(monty::rc_ptr< ::mosek::fusion::LinearDomain > _3407_ld,std::shared_ptr< monty::ndarray< int32_t,1 > > _3408_sparsity);
      static  monty::rc_ptr< ::mosek::fusion::RangeDomain > integral(monty::rc_ptr< ::mosek::fusion::RangeDomain > _3409_rd);
      static  monty::rc_ptr< ::mosek::fusion::LinearDomain > integral(monty::rc_ptr< ::mosek::fusion::LinearDomain > _3410_ld);
      static  monty::rc_ptr< ::mosek::fusion::ConeDomain > integral(monty::rc_ptr< ::mosek::fusion::ConeDomain > _3411_c);
      static  monty::rc_ptr< ::mosek::fusion::ConeDomain > axis(monty::rc_ptr< ::mosek::fusion::ConeDomain > _3412_c,int32_t _3413_a);
      static  monty::rc_ptr< ::mosek::fusion::ConeDomain > inDPowerCone(std::shared_ptr< monty::ndarray< double,1 > > _3414_alphas,std::shared_ptr< monty::ndarray< int32_t,1 > > _3415_dims);
      static  monty::rc_ptr< ::mosek::fusion::ConeDomain > inDPowerCone(std::shared_ptr< monty::ndarray< double,1 > > _3417_alphas,int32_t _3418_m);
      static  monty::rc_ptr< ::mosek::fusion::ConeDomain > inDPowerCone(std::shared_ptr< monty::ndarray< double,1 > > _3419_alphas);
      static  monty::rc_ptr< ::mosek::fusion::ConeDomain > inDPowerCone(double _3420_alpha,std::shared_ptr< monty::ndarray< int32_t,1 > > _3421_dims);
      static  monty::rc_ptr< ::mosek::fusion::ConeDomain > inDPowerCone(double _3423_alpha,int32_t _3424_m);
      static  monty::rc_ptr< ::mosek::fusion::ConeDomain > inDPowerCone(double _3425_alpha);
      static  monty::rc_ptr< ::mosek::fusion::ConeDomain > inPPowerCone(std::shared_ptr< monty::ndarray< double,1 > > _3426_alphas,std::shared_ptr< monty::ndarray< int32_t,1 > > _3427_dims);
      static  monty::rc_ptr< ::mosek::fusion::ConeDomain > inPPowerCone(std::shared_ptr< monty::ndarray< double,1 > > _3429_alphas,int32_t _3430_m);
      static  monty::rc_ptr< ::mosek::fusion::ConeDomain > inPPowerCone(std::shared_ptr< monty::ndarray< double,1 > > _3431_alphas);
      static  monty::rc_ptr< ::mosek::fusion::ConeDomain > inPPowerCone(double _3432_alpha,std::shared_ptr< monty::ndarray< int32_t,1 > > _3433_dims);
      static  monty::rc_ptr< ::mosek::fusion::ConeDomain > inPPowerCone(double _3435_alpha,int32_t _3436_m);
      static  monty::rc_ptr< ::mosek::fusion::ConeDomain > inPPowerCone(double _3437_alpha);
      static  monty::rc_ptr< ::mosek::fusion::ConeDomain > inDExpCone(std::shared_ptr< monty::ndarray< int32_t,1 > > _3438_dims);
      static  monty::rc_ptr< ::mosek::fusion::ConeDomain > inDExpCone(int32_t _3440_m);
      static  monty::rc_ptr< ::mosek::fusion::ConeDomain > inDExpCone();
      static  monty::rc_ptr< ::mosek::fusion::ConeDomain > inPExpCone(std::shared_ptr< monty::ndarray< int32_t,1 > > _3441_dims);
      static  monty::rc_ptr< ::mosek::fusion::ConeDomain > inPExpCone(int32_t _3443_m);
      static  monty::rc_ptr< ::mosek::fusion::ConeDomain > inPExpCone();
      static  monty::rc_ptr< ::mosek::fusion::ConeDomain > inDGeoMeanCone(std::shared_ptr< monty::ndarray< int32_t,1 > > _3444_dims);
      static  monty::rc_ptr< ::mosek::fusion::ConeDomain > inDGeoMeanCone(int32_t _3446_m,int32_t _3447_n);
      static  monty::rc_ptr< ::mosek::fusion::ConeDomain > inDGeoMeanCone(int32_t _3448_n);
      static  monty::rc_ptr< ::mosek::fusion::ConeDomain > inDGeoMeanCone();
      static  monty::rc_ptr< ::mosek::fusion::ConeDomain > inPGeoMeanCone(std::shared_ptr< monty::ndarray< int32_t,1 > > _3449_dims);
      static  monty::rc_ptr< ::mosek::fusion::ConeDomain > inPGeoMeanCone(int32_t _3451_m,int32_t _3452_n);
      static  monty::rc_ptr< ::mosek::fusion::ConeDomain > inPGeoMeanCone(int32_t _3453_n);
      static  monty::rc_ptr< ::mosek::fusion::ConeDomain > inPGeoMeanCone();
      static  monty::rc_ptr< ::mosek::fusion::ConeDomain > inRotatedQCone(std::shared_ptr< monty::ndarray< int32_t,1 > > _3454_dims);
      static  monty::rc_ptr< ::mosek::fusion::ConeDomain > inRotatedQCone(int32_t _3456_m,int32_t _3457_n);
      static  monty::rc_ptr< ::mosek::fusion::ConeDomain > inRotatedQCone(int32_t _3458_n);
      static  monty::rc_ptr< ::mosek::fusion::ConeDomain > inRotatedQCone();
      static  monty::rc_ptr< ::mosek::fusion::ConeDomain > inQCone(std::shared_ptr< monty::ndarray< int32_t,1 > > _3459_dims);
      static  monty::rc_ptr< ::mosek::fusion::ConeDomain > inQCone(int32_t _3461_m,int32_t _3462_n);
      static  monty::rc_ptr< ::mosek::fusion::ConeDomain > inQCone(int32_t _3463_n);
      static  monty::rc_ptr< ::mosek::fusion::ConeDomain > inQCone();
      static  monty::rc_ptr< ::mosek::fusion::ConeDomain > inSVecPSDCone(std::shared_ptr< monty::ndarray< int32_t,1 > > _3464_dims);
      static  monty::rc_ptr< ::mosek::fusion::ConeDomain > inSVecPSDCone(int32_t _3465_d1,int32_t _3466_d2);
      static  monty::rc_ptr< ::mosek::fusion::ConeDomain > inSVecPSDCone(int32_t _3467_n);
      static  monty::rc_ptr< ::mosek::fusion::ConeDomain > inSVecPSDCone();
      static  monty::rc_ptr< ::mosek::fusion::PSDDomain > isTrilPSD(int32_t _3468_n,int32_t _3469_m);
      static  monty::rc_ptr< ::mosek::fusion::PSDDomain > isTrilPSD(int32_t _3470_n);
      static  monty::rc_ptr< ::mosek::fusion::PSDDomain > isTrilPSD();
      static  monty::rc_ptr< ::mosek::fusion::PSDDomain > inPSDCone(int32_t _3471_n,int32_t _3472_m);
      static  monty::rc_ptr< ::mosek::fusion::PSDDomain > inPSDCone(int32_t _3473_n);
      static  monty::rc_ptr< ::mosek::fusion::PSDDomain > inPSDCone();
      static  monty::rc_ptr< ::mosek::fusion::RangeDomain > binary();
      static  monty::rc_ptr< ::mosek::fusion::RangeDomain > binary(std::shared_ptr< monty::ndarray< int32_t,1 > > _3474_dims);
      static  monty::rc_ptr< ::mosek::fusion::RangeDomain > binary(int32_t _3475_m,int32_t _3476_n);
      static  monty::rc_ptr< ::mosek::fusion::RangeDomain > binary(int32_t _3477_n);
      static  monty::rc_ptr< ::mosek::fusion::RangeDomain > inRange(monty::rc_ptr< ::mosek::fusion::Matrix > _3478_lbm,monty::rc_ptr< ::mosek::fusion::Matrix > _3479_ubm);
      static  monty::rc_ptr< ::mosek::fusion::RangeDomain > inRange(std::shared_ptr< monty::ndarray< double,2 > > _3480_lba,std::shared_ptr< monty::ndarray< double,2 > > _3481_uba);
      static  monty::rc_ptr< ::mosek::fusion::RangeDomain > inRange(std::shared_ptr< monty::ndarray< double,1 > > _3482_lba,std::shared_ptr< monty::ndarray< double,1 > > _3483_uba,std::shared_ptr< monty::ndarray< int32_t,1 > > _3484_dims);
      static  monty::rc_ptr< ::mosek::fusion::RangeDomain > inRange(std::shared_ptr< monty::ndarray< double,1 > > _3485_lba,double _3486_ub,std::shared_ptr< monty::ndarray< int32_t,1 > > _3487_dims);
      static  monty::rc_ptr< ::mosek::fusion::RangeDomain > inRange(double _3489_lb,std::shared_ptr< monty::ndarray< double,1 > > _3490_uba,std::shared_ptr< monty::ndarray< int32_t,1 > > _3491_dims);
      static  monty::rc_ptr< ::mosek::fusion::RangeDomain > inRange(double _3493_lb,double _3494_ub,std::shared_ptr< monty::ndarray< int32_t,1 > > _3495_dims);
      static  monty::rc_ptr< ::mosek::fusion::RangeDomain > inRange(std::shared_ptr< monty::ndarray< double,1 > > _3496_lba,std::shared_ptr< monty::ndarray< double,1 > > _3497_uba);
      static  monty::rc_ptr< ::mosek::fusion::RangeDomain > inRange(std::shared_ptr< monty::ndarray< double,1 > > _3498_lba,double _3499_ub);
      static  monty::rc_ptr< ::mosek::fusion::RangeDomain > inRange(double _3501_lb,std::shared_ptr< monty::ndarray< double,1 > > _3502_uba);
      static  monty::rc_ptr< ::mosek::fusion::RangeDomain > inRange(double _3504_lb,double _3505_ub);
      static  monty::rc_ptr< ::mosek::fusion::LinearDomain > greaterThan(monty::rc_ptr< ::mosek::fusion::Matrix > _3506_mx);
      static  monty::rc_ptr< ::mosek::fusion::LinearDomain > greaterThan(std::shared_ptr< monty::ndarray< double,1 > > _3507_a1,std::shared_ptr< monty::ndarray< int32_t,1 > > _3508_dims);
      static  monty::rc_ptr< ::mosek::fusion::LinearDomain > greaterThan(std::shared_ptr< monty::ndarray< double,2 > > _3509_a2);
      static  monty::rc_ptr< ::mosek::fusion::LinearDomain > greaterThan(std::shared_ptr< monty::ndarray< double,1 > > _3512_a1);
      static  monty::rc_ptr< ::mosek::fusion::LinearDomain > greaterThan(double _3513_b,std::shared_ptr< monty::ndarray< int32_t,1 > > _3514_dims);
      static  monty::rc_ptr< ::mosek::fusion::LinearDomain > greaterThan(double _3516_b,int32_t _3517_m,int32_t _3518_n);
      static  monty::rc_ptr< ::mosek::fusion::LinearDomain > greaterThan(double _3520_b,int32_t _3521_n);
      static  monty::rc_ptr< ::mosek::fusion::LinearDomain > greaterThan(double _3523_b);
      static  monty::rc_ptr< ::mosek::fusion::LinearDomain > lessThan(monty::rc_ptr< ::mosek::fusion::Matrix > _3524_mx);
      static  monty::rc_ptr< ::mosek::fusion::LinearDomain > lessThan(std::shared_ptr< monty::ndarray< double,1 > > _3525_a1,std::shared_ptr< monty::ndarray< int32_t,1 > > _3526_dims);
      static  monty::rc_ptr< ::mosek::fusion::LinearDomain > lessThan(std::shared_ptr< monty::ndarray< double,2 > > _3527_a2);
      static  monty::rc_ptr< ::mosek::fusion::LinearDomain > lessThan(std::shared_ptr< monty::ndarray< double,1 > > _3530_a1);
      static  monty::rc_ptr< ::mosek::fusion::LinearDomain > lessThan(double _3531_b,std::shared_ptr< monty::ndarray< int32_t,1 > > _3532_dims);
      static  monty::rc_ptr< ::mosek::fusion::LinearDomain > lessThan(double _3533_b,int32_t _3534_m,int32_t _3535_n);
      static  monty::rc_ptr< ::mosek::fusion::LinearDomain > lessThan(double _3536_b,int32_t _3537_n);
      static  monty::rc_ptr< ::mosek::fusion::LinearDomain > lessThan(double _3538_b);
      static  monty::rc_ptr< ::mosek::fusion::LinearDomain > equalsTo(monty::rc_ptr< ::mosek::fusion::Matrix > _3539_mx);
      static  monty::rc_ptr< ::mosek::fusion::LinearDomain > equalsTo(std::shared_ptr< monty::ndarray< double,1 > > _3540_a1,std::shared_ptr< monty::ndarray< int32_t,1 > > _3541_dims);
      static  monty::rc_ptr< ::mosek::fusion::LinearDomain > equalsTo(std::shared_ptr< monty::ndarray< double,2 > > _3542_a2);
      static  monty::rc_ptr< ::mosek::fusion::LinearDomain > equalsTo(std::shared_ptr< monty::ndarray< double,1 > > _3545_a1);
      static  monty::rc_ptr< ::mosek::fusion::LinearDomain > equalsTo(double _3546_b,std::shared_ptr< monty::ndarray< int32_t,1 > > _3547_dims);
      static  monty::rc_ptr< ::mosek::fusion::LinearDomain > equalsTo(double _3548_b,int32_t _3549_m,int32_t _3550_n);
      static  monty::rc_ptr< ::mosek::fusion::LinearDomain > equalsTo(double _3551_b,int32_t _3552_n);
      static  monty::rc_ptr< ::mosek::fusion::LinearDomain > equalsTo(double _3553_b);
      static  monty::rc_ptr< ::mosek::fusion::LinearDomain > unbounded(std::shared_ptr< monty::ndarray< int32_t,1 > > _3554_dims);
      static  monty::rc_ptr< ::mosek::fusion::LinearDomain > unbounded(int32_t _3556_m,int32_t _3557_n);
      static  monty::rc_ptr< ::mosek::fusion::LinearDomain > unbounded(int32_t _3558_n);
      static  monty::rc_ptr< ::mosek::fusion::LinearDomain > unbounded();
    }; // struct Domain;

    struct p_ExprCode
    {
      ExprCode * _pubthis;
      static mosek::fusion::p_ExprCode* _get_impl(mosek::fusion::ExprCode * _inst){ assert(_inst); assert(_inst->_impl); return _inst->_impl; }
      static mosek::fusion::p_ExprCode * _get_impl(mosek::fusion::ExprCode::t _inst) { return _get_impl(_inst.get()); }
      p_ExprCode(ExprCode * _pubthis);
      virtual ~p_ExprCode() { /* std::cout << "~p_ExprCode" << std::endl;*/ };

      virtual void destroy();

      static  void inplace_relocate(std::shared_ptr< monty::ndarray< int32_t,1 > > _3559_code,int32_t _3560_from_offset,int32_t _3561_num,int32_t _3562_const_base);
      static  std::string op2str(int32_t _3564_op);
      static  void eval_add_list(std::shared_ptr< monty::ndarray< int32_t,1 > > _3565_code,std::shared_ptr< monty::ndarray< int32_t,1 > > _3566_ptr,std::shared_ptr< monty::ndarray< double,1 > > _3567_consts,int32_t _3568_offset,std::shared_ptr< monty::ndarray< double,1 > > _3569_target,std::shared_ptr< monty::ndarray< double,1 > > _3570_P,monty::rc_ptr< ::mosek::fusion::WorkStack > _3571_xs);
      static  void eval_add_list(std::shared_ptr< monty::ndarray< int32_t,1 > > _3579_code,std::shared_ptr< monty::ndarray< int32_t,1 > > _3580_ptr,std::shared_ptr< monty::ndarray< double,1 > > _3581_consts,std::shared_ptr< monty::ndarray< double,1 > > _3582_target,std::shared_ptr< monty::ndarray< double,1 > > _3583_P,monty::rc_ptr< ::mosek::fusion::WorkStack > _3584_xs);
      static  int32_t emit_sum(std::shared_ptr< monty::ndarray< int32_t,1 > > _3585_tgt,int32_t _3586_ofs,int32_t _3587_num);
      static  int32_t emit_inv(std::shared_ptr< monty::ndarray< int32_t,1 > > _3588_tgt,int32_t _3589_ofs);
      static  int32_t emit_mul(std::shared_ptr< monty::ndarray< int32_t,1 > > _3590_tgt,int32_t _3591_ofs);
      static  int32_t emit_neg(std::shared_ptr< monty::ndarray< int32_t,1 > > _3592_tgt,int32_t _3593_ofs);
      static  int32_t emit_add(std::shared_ptr< monty::ndarray< int32_t,1 > > _3594_tgt,int32_t _3595_ofs);
      static  int32_t emit_constref(std::shared_ptr< monty::ndarray< int32_t,1 > > _3596_tgt,int32_t _3597_ofs,int32_t _3598_i);
      static  int32_t emit_paramref(std::shared_ptr< monty::ndarray< int32_t,1 > > _3599_tgt,int32_t _3600_ofs,int32_t _3601_i);
      static  int32_t emit_nop(std::shared_ptr< monty::ndarray< int32_t,1 > > _3602_tgt,int32_t _3603_ofs);
    }; // struct ExprCode;

    struct p_Param
    {
      Param * _pubthis;
      static mosek::fusion::p_Param* _get_impl(mosek::fusion::Param * _inst){ assert(_inst); assert(_inst->_impl); return _inst->_impl; }
      static mosek::fusion::p_Param * _get_impl(mosek::fusion::Param::t _inst) { return _get_impl(_inst.get()); }
      p_Param(Param * _pubthis);
      virtual ~p_Param() { /* std::cout << "~p_Param" << std::endl;*/ };

      virtual void destroy();

      static  monty::rc_ptr< ::mosek::fusion::Parameter > repeat(monty::rc_ptr< ::mosek::fusion::Parameter > _3612_p,int32_t _3613_n,int32_t _3614_dim);
      static  monty::rc_ptr< ::mosek::fusion::Parameter > stack(int32_t _3616_dim,monty::rc_ptr< ::mosek::fusion::Parameter > _3617_p1,monty::rc_ptr< ::mosek::fusion::Parameter > _3618_p2,monty::rc_ptr< ::mosek::fusion::Parameter > _3619_p3);
      static  monty::rc_ptr< ::mosek::fusion::Parameter > stack(int32_t _3620_dim,monty::rc_ptr< ::mosek::fusion::Parameter > _3621_p1,monty::rc_ptr< ::mosek::fusion::Parameter > _3622_p2);
      static  monty::rc_ptr< ::mosek::fusion::Parameter > stack(int32_t _3623_dim,std::shared_ptr< monty::ndarray< monty::rc_ptr< ::mosek::fusion::Parameter >,1 > > _3624_p);
      static  monty::rc_ptr< ::mosek::fusion::Parameter > stack(std::shared_ptr< monty::ndarray< std::shared_ptr< monty::ndarray< monty::rc_ptr< ::mosek::fusion::Parameter >,1 > >,1 > > _3625_p);
      static  monty::rc_ptr< ::mosek::fusion::Parameter > hstack(monty::rc_ptr< ::mosek::fusion::Parameter > _3627_p1,monty::rc_ptr< ::mosek::fusion::Parameter > _3628_p2,monty::rc_ptr< ::mosek::fusion::Parameter > _3629_p3);
      static  monty::rc_ptr< ::mosek::fusion::Parameter > hstack(monty::rc_ptr< ::mosek::fusion::Parameter > _3630_p1,monty::rc_ptr< ::mosek::fusion::Parameter > _3631_p2);
      static  monty::rc_ptr< ::mosek::fusion::Parameter > hstack(std::shared_ptr< monty::ndarray< monty::rc_ptr< ::mosek::fusion::Parameter >,1 > > _3632_p);
      static  monty::rc_ptr< ::mosek::fusion::Parameter > vstack(monty::rc_ptr< ::mosek::fusion::Parameter > _3633_p1,monty::rc_ptr< ::mosek::fusion::Parameter > _3634_p2,monty::rc_ptr< ::mosek::fusion::Parameter > _3635_p3);
      static  monty::rc_ptr< ::mosek::fusion::Parameter > vstack(monty::rc_ptr< ::mosek::fusion::Parameter > _3636_p1,monty::rc_ptr< ::mosek::fusion::Parameter > _3637_p2);
      static  monty::rc_ptr< ::mosek::fusion::Parameter > vstack(std::shared_ptr< monty::ndarray< monty::rc_ptr< ::mosek::fusion::Parameter >,1 > > _3638_p);
      static  monty::rc_ptr< ::mosek::fusion::Parameter > dstack(std::shared_ptr< monty::ndarray< monty::rc_ptr< ::mosek::fusion::Parameter >,1 > > _3639_p,int32_t _3640_dim);
    }; // struct Param;

    struct p_ParameterImpl : public /*implements*/ virtual ::mosek::fusion::Parameter
    {
      ParameterImpl * _pubthis;
      static mosek::fusion::p_ParameterImpl* _get_impl(mosek::fusion::ParameterImpl * _inst){ assert(_inst); assert(_inst->_impl); return _inst->_impl; }
      static mosek::fusion::p_ParameterImpl * _get_impl(mosek::fusion::ParameterImpl::t _inst) { return _get_impl(_inst.get()); }
      p_ParameterImpl(ParameterImpl * _pubthis);
      virtual ~p_ParameterImpl() { /* std::cout << "~p_ParameterImpl" << std::endl;*/ };
      int64_t size{};
      std::shared_ptr< monty::ndarray< int32_t,1 > > nidxs{};
      std::shared_ptr< monty::ndarray< int64_t,1 > > sp{};
      std::shared_ptr< monty::ndarray< int32_t,1 > > shape{};
      monty::rc_ptr< ::mosek::fusion::Model > model{};

      virtual void destroy();

      static ParameterImpl::t _new_ParameterImpl(monty::rc_ptr< ::mosek::fusion::ParameterImpl > _4429_other,monty::rc_ptr< ::mosek::fusion::Model > _4430_model);
      void _initialize(monty::rc_ptr< ::mosek::fusion::ParameterImpl > _4429_other,monty::rc_ptr< ::mosek::fusion::Model > _4430_model);
      static ParameterImpl::t _new_ParameterImpl(monty::rc_ptr< ::mosek::fusion::Model > _4431_model,std::shared_ptr< monty::ndarray< int32_t,1 > > _4432_shape,std::shared_ptr< monty::ndarray< int64_t,1 > > _4433_sp,std::shared_ptr< monty::ndarray< int32_t,1 > > _4434_nidxs);
      void _initialize(monty::rc_ptr< ::mosek::fusion::Model > _4431_model,std::shared_ptr< monty::ndarray< int32_t,1 > > _4432_shape,std::shared_ptr< monty::ndarray< int64_t,1 > > _4433_sp,std::shared_ptr< monty::ndarray< int32_t,1 > > _4434_nidxs);
      virtual monty::rc_ptr< ::mosek::fusion::Parameter > __mosek_2fusion_2ParameterImpl__clone(monty::rc_ptr< ::mosek::fusion::Model > _4435_m) ;
      virtual monty::rc_ptr< ::mosek::fusion::Parameter > __mosek_2fusion_2Parameter__clone(monty::rc_ptr< ::mosek::fusion::Model > _4435_m) { return __mosek_2fusion_2ParameterImpl__clone(_4435_m); }
      virtual /* override */ std::string toString() ;
      virtual monty::rc_ptr< ::mosek::fusion::Parameter > __mosek_2fusion_2ParameterImpl__pick(std::shared_ptr< monty::ndarray< int32_t,2 > > _4438_midxs) ;
      virtual monty::rc_ptr< ::mosek::fusion::Parameter > __mosek_2fusion_2Parameter__pick(std::shared_ptr< monty::ndarray< int32_t,2 > > _4438_midxs) { return __mosek_2fusion_2ParameterImpl__pick(_4438_midxs); }
      virtual monty::rc_ptr< ::mosek::fusion::Parameter > __mosek_2fusion_2ParameterImpl__pick(std::shared_ptr< monty::ndarray< int32_t,1 > > _4460_idxs) ;
      virtual monty::rc_ptr< ::mosek::fusion::Parameter > __mosek_2fusion_2Parameter__pick(std::shared_ptr< monty::ndarray< int32_t,1 > > _4460_idxs) { return __mosek_2fusion_2ParameterImpl__pick(_4460_idxs); }
      virtual monty::rc_ptr< ::mosek::fusion::Expression > __mosek_2fusion_2ParameterImpl__index(std::shared_ptr< monty::ndarray< int32_t,1 > > _4471_indexes) ;
      virtual monty::rc_ptr< ::mosek::fusion::Expression > __mosek_2fusion_2Expression__index(std::shared_ptr< monty::ndarray< int32_t,1 > > _4471_indexes) { return __mosek_2fusion_2ParameterImpl__index(_4471_indexes); }
      virtual monty::rc_ptr< ::mosek::fusion::Expression > __mosek_2fusion_2ParameterImpl__index(int32_t _4480_i) ;
      virtual monty::rc_ptr< ::mosek::fusion::Expression > __mosek_2fusion_2Expression__index(int32_t _4480_i) { return __mosek_2fusion_2ParameterImpl__index(_4480_i); }
      virtual void eval(monty::rc_ptr< ::mosek::fusion::WorkStack > _4482_rs,monty::rc_ptr< ::mosek::fusion::WorkStack > _4483_ws,monty::rc_ptr< ::mosek::fusion::WorkStack > _4484_xs) ;
      virtual void getSp(std::shared_ptr< monty::ndarray< int64_t,1 > > _4506_dest,int32_t _4507_offset) ;
      virtual bool isSparse() ;
      virtual monty::rc_ptr< ::mosek::fusion::Parameter > __mosek_2fusion_2ParameterImpl__slice(std::shared_ptr< monty::ndarray< int32_t,1 > > _4510_astart,std::shared_ptr< monty::ndarray< int32_t,1 > > _4511_astop) ;
      virtual monty::rc_ptr< ::mosek::fusion::Parameter > __mosek_2fusion_2Parameter__slice(std::shared_ptr< monty::ndarray< int32_t,1 > > _4510_astart,std::shared_ptr< monty::ndarray< int32_t,1 > > _4511_astop) { return __mosek_2fusion_2ParameterImpl__slice(_4510_astart,_4511_astop); }
      virtual monty::rc_ptr< ::mosek::fusion::Parameter > __mosek_2fusion_2ParameterImpl__slice(int32_t _4543_start,int32_t _4544_stop) ;
      virtual monty::rc_ptr< ::mosek::fusion::Parameter > __mosek_2fusion_2Parameter__slice(int32_t _4543_start,int32_t _4544_stop) { return __mosek_2fusion_2ParameterImpl__slice(_4543_start,_4544_stop); }
      virtual monty::rc_ptr< ::mosek::fusion::Parameter > __mosek_2fusion_2ParameterImpl__reshape(std::shared_ptr< monty::ndarray< int32_t,1 > > _4552_dims) ;
      virtual monty::rc_ptr< ::mosek::fusion::Parameter > __mosek_2fusion_2Parameter__reshape(std::shared_ptr< monty::ndarray< int32_t,1 > > _4552_dims) { return __mosek_2fusion_2ParameterImpl__reshape(_4552_dims); }
      virtual monty::rc_ptr< ::mosek::fusion::Expression > __mosek_2fusion_2ParameterImpl__asExpr() ;
      virtual monty::rc_ptr< ::mosek::fusion::Expression > __mosek_2fusion_2Parameter__asExpr() { return __mosek_2fusion_2ParameterImpl__asExpr(); }
      virtual int64_t getSize() ;
      virtual int32_t getNumNonzero() ;
      virtual int32_t getND() ;
      virtual std::shared_ptr< monty::ndarray< int32_t,1 > > getShape() ;
      virtual int32_t getDim(int32_t _4553_i) ;
      virtual void getAllIndexes(std::shared_ptr< monty::ndarray< int32_t,1 > > _4554_dst,int32_t _4555_ofs) ;
      virtual int32_t getIndex(int32_t _4557_i) ;
      virtual std::shared_ptr< monty::ndarray< double,1 > > getValue() ;
      virtual void setValue(std::shared_ptr< monty::ndarray< double,2 > > _4558_values2) ;
      virtual void setValue(std::shared_ptr< monty::ndarray< double,1 > > _4564_values) ;
      virtual void setValue(double _4567_value) ;
      virtual monty::rc_ptr< ::mosek::fusion::Model > __mosek_2fusion_2ParameterImpl__getModel() ;
      virtual monty::rc_ptr< ::mosek::fusion::Model > __mosek_2fusion_2Parameter__getModel() { return __mosek_2fusion_2ParameterImpl__getModel(); }
    }; // struct ParameterImpl;

    struct p_ExprRangeDomain : public /*implements*/ virtual ::mosek::fusion::ExprDomain
    {
      ExprRangeDomain * _pubthis;
      static mosek::fusion::p_ExprRangeDomain* _get_impl(mosek::fusion::ExprRangeDomain * _inst){ assert(_inst); assert(_inst->_impl); return _inst->_impl; }
      static mosek::fusion::p_ExprRangeDomain * _get_impl(mosek::fusion::ExprRangeDomain::t _inst) { return _get_impl(_inst.get()); }
      p_ExprRangeDomain(ExprRangeDomain * _pubthis);
      virtual ~p_ExprRangeDomain() { /* std::cout << "~p_ExprRangeDomain" << std::endl;*/ };
      monty::rc_ptr< ::mosek::fusion::RangeDomain > dom{};
      monty::rc_ptr< ::mosek::fusion::Expression > expr{};

      virtual void destroy();

      static ExprRangeDomain::t _new_ExprRangeDomain(monty::rc_ptr< ::mosek::fusion::Expression > _7271_expr,monty::rc_ptr< ::mosek::fusion::RangeDomain > _7272_dom);
      void _initialize(monty::rc_ptr< ::mosek::fusion::Expression > _7271_expr,monty::rc_ptr< ::mosek::fusion::RangeDomain > _7272_dom);
      virtual monty::rc_ptr< ::mosek::fusion::Term > __mosek_2fusion_2ExprRangeDomain__toDJCTerm() ;
      virtual monty::rc_ptr< ::mosek::fusion::Term > __mosek_2fusion_2ExprDomain__toDJCTerm() { return __mosek_2fusion_2ExprRangeDomain__toDJCTerm(); }
    }; // struct ExprRangeDomain;

    struct p_ExprPSDDomain
    {
      ExprPSDDomain * _pubthis;
      static mosek::fusion::p_ExprPSDDomain* _get_impl(mosek::fusion::ExprPSDDomain * _inst){ assert(_inst); assert(_inst->_impl); return _inst->_impl; }
      static mosek::fusion::p_ExprPSDDomain * _get_impl(mosek::fusion::ExprPSDDomain::t _inst) { return _get_impl(_inst.get()); }
      p_ExprPSDDomain(ExprPSDDomain * _pubthis);
      virtual ~p_ExprPSDDomain() { /* std::cout << "~p_ExprPSDDomain" << std::endl;*/ };
      monty::rc_ptr< ::mosek::fusion::PSDDomain > dom{};
      monty::rc_ptr< ::mosek::fusion::Expression > expr{};

      virtual void destroy();

      static ExprPSDDomain::t _new_ExprPSDDomain(monty::rc_ptr< ::mosek::fusion::Expression > _7273_expr,monty::rc_ptr< ::mosek::fusion::PSDDomain > _7274_dom);
      void _initialize(monty::rc_ptr< ::mosek::fusion::Expression > _7273_expr,monty::rc_ptr< ::mosek::fusion::PSDDomain > _7274_dom);
    }; // struct ExprPSDDomain;

    struct p_ExprConicDomain
    {
      ExprConicDomain * _pubthis;
      static mosek::fusion::p_ExprConicDomain* _get_impl(mosek::fusion::ExprConicDomain * _inst){ assert(_inst); assert(_inst->_impl); return _inst->_impl; }
      static mosek::fusion::p_ExprConicDomain * _get_impl(mosek::fusion::ExprConicDomain::t _inst) { return _get_impl(_inst.get()); }
      p_ExprConicDomain(ExprConicDomain * _pubthis);
      virtual ~p_ExprConicDomain() { /* std::cout << "~p_ExprConicDomain" << std::endl;*/ };
      monty::rc_ptr< ::mosek::fusion::ConeDomain > dom{};
      monty::rc_ptr< ::mosek::fusion::Expression > expr{};

      virtual void destroy();

      static ExprConicDomain::t _new_ExprConicDomain(monty::rc_ptr< ::mosek::fusion::Expression > _7275_expr,monty::rc_ptr< ::mosek::fusion::ConeDomain > _7276_dom);
      void _initialize(monty::rc_ptr< ::mosek::fusion::Expression > _7275_expr,monty::rc_ptr< ::mosek::fusion::ConeDomain > _7276_dom);
    }; // struct ExprConicDomain;

    struct p_ExprLinearDomain : public /*implements*/ virtual ::mosek::fusion::ExprDomain
    {
      ExprLinearDomain * _pubthis;
      static mosek::fusion::p_ExprLinearDomain* _get_impl(mosek::fusion::ExprLinearDomain * _inst){ assert(_inst); assert(_inst->_impl); return _inst->_impl; }
      static mosek::fusion::p_ExprLinearDomain * _get_impl(mosek::fusion::ExprLinearDomain::t _inst) { return _get_impl(_inst.get()); }
      p_ExprLinearDomain(ExprLinearDomain * _pubthis);
      virtual ~p_ExprLinearDomain() { /* std::cout << "~p_ExprLinearDomain" << std::endl;*/ };
      monty::rc_ptr< ::mosek::fusion::LinearDomain > dom{};
      monty::rc_ptr< ::mosek::fusion::Expression > expr{};

      virtual void destroy();

      static ExprLinearDomain::t _new_ExprLinearDomain(monty::rc_ptr< ::mosek::fusion::Expression > _7277_expr,monty::rc_ptr< ::mosek::fusion::LinearDomain > _7278_dom);
      void _initialize(monty::rc_ptr< ::mosek::fusion::Expression > _7277_expr,monty::rc_ptr< ::mosek::fusion::LinearDomain > _7278_dom);
      virtual monty::rc_ptr< ::mosek::fusion::Term > __mosek_2fusion_2ExprLinearDomain__toDJCTerm() ;
      virtual monty::rc_ptr< ::mosek::fusion::Term > __mosek_2fusion_2ExprDomain__toDJCTerm() { return __mosek_2fusion_2ExprLinearDomain__toDJCTerm(); }
    }; // struct ExprLinearDomain;

    struct p_BaseExpression : public /*implements*/ virtual ::mosek::fusion::Expression
    {
      BaseExpression * _pubthis;
      static mosek::fusion::p_BaseExpression* _get_impl(mosek::fusion::BaseExpression * _inst){ assert(_inst); assert(_inst->_impl); return _inst->_impl; }
      static mosek::fusion::p_BaseExpression * _get_impl(mosek::fusion::BaseExpression::t _inst) { return _get_impl(_inst.get()); }
      p_BaseExpression(BaseExpression * _pubthis);
      virtual ~p_BaseExpression() { /* std::cout << "~p_BaseExpression" << std::endl;*/ };
      std::shared_ptr< monty::ndarray< int32_t,1 > > shape{};

      virtual void destroy();

      static BaseExpression::t _new_BaseExpression(std::shared_ptr< monty::ndarray< int32_t,1 > > _7279_shape);
      void _initialize(std::shared_ptr< monty::ndarray< int32_t,1 > > _7279_shape);
      virtual /* override */ std::string toString() ;
      virtual void printStack(monty::rc_ptr< ::mosek::fusion::WorkStack > _7280_rs) ;
      virtual void eval(monty::rc_ptr< ::mosek::fusion::WorkStack > _7306_rs,monty::rc_ptr< ::mosek::fusion::WorkStack > _7307_ws,monty::rc_ptr< ::mosek::fusion::WorkStack > _7308_xs) { throw monty::AbstractClassError("Call to abstract method"); }
      virtual monty::rc_ptr< ::mosek::fusion::Expression > __mosek_2fusion_2BaseExpression__pick(std::shared_ptr< monty::ndarray< int32_t,2 > > _7309_indexrows) ;
      virtual monty::rc_ptr< ::mosek::fusion::Expression > __mosek_2fusion_2Expression__pick(std::shared_ptr< monty::ndarray< int32_t,2 > > _7309_indexrows) { return __mosek_2fusion_2BaseExpression__pick(_7309_indexrows); }
      virtual monty::rc_ptr< ::mosek::fusion::Expression > __mosek_2fusion_2BaseExpression__pick(std::shared_ptr< monty::ndarray< int32_t,1 > > _7310_indexes) ;
      virtual monty::rc_ptr< ::mosek::fusion::Expression > __mosek_2fusion_2Expression__pick(std::shared_ptr< monty::ndarray< int32_t,1 > > _7310_indexes) { return __mosek_2fusion_2BaseExpression__pick(_7310_indexes); }
      virtual monty::rc_ptr< ::mosek::fusion::Expression > __mosek_2fusion_2BaseExpression__index(std::shared_ptr< monty::ndarray< int32_t,1 > > _7313_indexes) ;
      virtual monty::rc_ptr< ::mosek::fusion::Expression > __mosek_2fusion_2Expression__index(std::shared_ptr< monty::ndarray< int32_t,1 > > _7313_indexes) { return __mosek_2fusion_2BaseExpression__index(_7313_indexes); }
      virtual monty::rc_ptr< ::mosek::fusion::Expression > __mosek_2fusion_2BaseExpression__index(int32_t _7316_i) ;
      virtual monty::rc_ptr< ::mosek::fusion::Expression > __mosek_2fusion_2Expression__index(int32_t _7316_i) { return __mosek_2fusion_2BaseExpression__index(_7316_i); }
      virtual monty::rc_ptr< ::mosek::fusion::Expression > __mosek_2fusion_2BaseExpression__slice(std::shared_ptr< monty::ndarray< int32_t,1 > > _7318_firsta,std::shared_ptr< monty::ndarray< int32_t,1 > > _7319_lasta) ;
      virtual monty::rc_ptr< ::mosek::fusion::Expression > __mosek_2fusion_2Expression__slice(std::shared_ptr< monty::ndarray< int32_t,1 > > _7318_firsta,std::shared_ptr< monty::ndarray< int32_t,1 > > _7319_lasta) { return __mosek_2fusion_2BaseExpression__slice(_7318_firsta,_7319_lasta); }
      virtual monty::rc_ptr< ::mosek::fusion::Expression > __mosek_2fusion_2BaseExpression__slice(int32_t _7320_first,int32_t _7321_last) ;
      virtual monty::rc_ptr< ::mosek::fusion::Expression > __mosek_2fusion_2Expression__slice(int32_t _7320_first,int32_t _7321_last) { return __mosek_2fusion_2BaseExpression__slice(_7320_first,_7321_last); }
      virtual int64_t getSize() ;
      virtual int32_t getND() ;
      virtual int32_t getDim(int32_t _7322_d) ;
      virtual std::shared_ptr< monty::ndarray< int32_t,1 > > getShape() ;
    }; // struct BaseExpression;

    struct p_ExprParameter : public ::mosek::fusion::p_BaseExpression
    {
      ExprParameter * _pubthis;
      static mosek::fusion::p_ExprParameter* _get_impl(mosek::fusion::ExprParameter * _inst){ return static_cast< mosek::fusion::p_ExprParameter* >(mosek::fusion::p_BaseExpression::_get_impl(_inst)); }
      static mosek::fusion::p_ExprParameter * _get_impl(mosek::fusion::ExprParameter::t _inst) { return _get_impl(_inst.get()); }
      p_ExprParameter(ExprParameter * _pubthis);
      virtual ~p_ExprParameter() { /* std::cout << "~p_ExprParameter" << std::endl;*/ };
      monty::rc_ptr< ::mosek::fusion::Parameter > p{};

      virtual void destroy();

      static ExprParameter::t _new_ExprParameter(monty::rc_ptr< ::mosek::fusion::Parameter > _3604_p);
      void _initialize(monty::rc_ptr< ::mosek::fusion::Parameter > _3604_p);
      virtual /* override */ void eval(monty::rc_ptr< ::mosek::fusion::WorkStack > _3605_rs,monty::rc_ptr< ::mosek::fusion::WorkStack > _3606_ws,monty::rc_ptr< ::mosek::fusion::WorkStack > _3607_xs) ;
      virtual /* override */ monty::rc_ptr< ::mosek::fusion::Expression > __mosek_2fusion_2ExprParameter__slice(std::shared_ptr< monty::ndarray< int32_t,1 > > _3608_start,std::shared_ptr< monty::ndarray< int32_t,1 > > _3609_stop) ;
      virtual monty::rc_ptr< ::mosek::fusion::Expression > __mosek_2fusion_2BaseExpression__slice(std::shared_ptr< monty::ndarray< int32_t,1 > > _3608_start,std::shared_ptr< monty::ndarray< int32_t,1 > > _3609_stop) { return __mosek_2fusion_2ExprParameter__slice(_3608_start,_3609_stop); }
      virtual /* override */ monty::rc_ptr< ::mosek::fusion::Expression > __mosek_2fusion_2ExprParameter__slice(int32_t _3610_start,int32_t _3611_stop) ;
      virtual monty::rc_ptr< ::mosek::fusion::Expression > __mosek_2fusion_2BaseExpression__slice(int32_t _3610_start,int32_t _3611_stop) { return __mosek_2fusion_2ExprParameter__slice(_3610_start,_3611_stop); }
      virtual /* override */ std::string toString() ;
    }; // struct ExprParameter;

    struct p_ExprMulParamScalarExpr : public ::mosek::fusion::p_BaseExpression
    {
      ExprMulParamScalarExpr * _pubthis;
      static mosek::fusion::p_ExprMulParamScalarExpr* _get_impl(mosek::fusion::ExprMulParamScalarExpr * _inst){ return static_cast< mosek::fusion::p_ExprMulParamScalarExpr* >(mosek::fusion::p_BaseExpression::_get_impl(_inst)); }
      static mosek::fusion::p_ExprMulParamScalarExpr * _get_impl(mosek::fusion::ExprMulParamScalarExpr::t _inst) { return _get_impl(_inst.get()); }
      p_ExprMulParamScalarExpr(ExprMulParamScalarExpr * _pubthis);
      virtual ~p_ExprMulParamScalarExpr() { /* std::cout << "~p_ExprMulParamScalarExpr" << std::endl;*/ };
      monty::rc_ptr< ::mosek::fusion::Expression > e{};
      monty::rc_ptr< ::mosek::fusion::Parameter > p{};

      virtual void destroy();

      static ExprMulParamScalarExpr::t _new_ExprMulParamScalarExpr(monty::rc_ptr< ::mosek::fusion::Parameter > _3699_p,monty::rc_ptr< ::mosek::fusion::Expression > _3700_e);
      void _initialize(monty::rc_ptr< ::mosek::fusion::Parameter > _3699_p,monty::rc_ptr< ::mosek::fusion::Expression > _3700_e);
      virtual /* override */ void eval(monty::rc_ptr< ::mosek::fusion::WorkStack > _3701_rs,monty::rc_ptr< ::mosek::fusion::WorkStack > _3702_ws,monty::rc_ptr< ::mosek::fusion::WorkStack > _3703_xs) ;
      virtual /* override */ std::string toString() ;
    }; // struct ExprMulParamScalarExpr;

    struct p_ExprMulParamScalar : public ::mosek::fusion::p_BaseExpression
    {
      ExprMulParamScalar * _pubthis;
      static mosek::fusion::p_ExprMulParamScalar* _get_impl(mosek::fusion::ExprMulParamScalar * _inst){ return static_cast< mosek::fusion::p_ExprMulParamScalar* >(mosek::fusion::p_BaseExpression::_get_impl(_inst)); }
      static mosek::fusion::p_ExprMulParamScalar * _get_impl(mosek::fusion::ExprMulParamScalar::t _inst) { return _get_impl(_inst.get()); }
      p_ExprMulParamScalar(ExprMulParamScalar * _pubthis);
      virtual ~p_ExprMulParamScalar() { /* std::cout << "~p_ExprMulParamScalar" << std::endl;*/ };
      monty::rc_ptr< ::mosek::fusion::Expression > e{};
      monty::rc_ptr< ::mosek::fusion::Parameter > p{};

      virtual void destroy();

      static ExprMulParamScalar::t _new_ExprMulParamScalar(monty::rc_ptr< ::mosek::fusion::Parameter > _3754_p,monty::rc_ptr< ::mosek::fusion::Expression > _3755_e);
      void _initialize(monty::rc_ptr< ::mosek::fusion::Parameter > _3754_p,monty::rc_ptr< ::mosek::fusion::Expression > _3755_e);
      virtual /* override */ void eval(monty::rc_ptr< ::mosek::fusion::WorkStack > _3756_rs,monty::rc_ptr< ::mosek::fusion::WorkStack > _3757_ws,monty::rc_ptr< ::mosek::fusion::WorkStack > _3758_xs) ;
      virtual /* override */ std::string toString() ;
    }; // struct ExprMulParamScalar;

    struct p_ExprMulParamDiagLeft : public ::mosek::fusion::p_BaseExpression
    {
      ExprMulParamDiagLeft * _pubthis;
      static mosek::fusion::p_ExprMulParamDiagLeft* _get_impl(mosek::fusion::ExprMulParamDiagLeft * _inst){ return static_cast< mosek::fusion::p_ExprMulParamDiagLeft* >(mosek::fusion::p_BaseExpression::_get_impl(_inst)); }
      static mosek::fusion::p_ExprMulParamDiagLeft * _get_impl(mosek::fusion::ExprMulParamDiagLeft::t _inst) { return _get_impl(_inst.get()); }
      p_ExprMulParamDiagLeft(ExprMulParamDiagLeft * _pubthis);
      virtual ~p_ExprMulParamDiagLeft() { /* std::cout << "~p_ExprMulParamDiagLeft" << std::endl;*/ };
      monty::rc_ptr< ::mosek::fusion::Expression > e{};
      monty::rc_ptr< ::mosek::fusion::Parameter > p{};

      virtual void destroy();

      static ExprMulParamDiagLeft::t _new_ExprMulParamDiagLeft(monty::rc_ptr< ::mosek::fusion::Parameter > _3801_p,monty::rc_ptr< ::mosek::fusion::Expression > _3802_e);
      void _initialize(monty::rc_ptr< ::mosek::fusion::Parameter > _3801_p,monty::rc_ptr< ::mosek::fusion::Expression > _3802_e);
      virtual /* override */ void eval(monty::rc_ptr< ::mosek::fusion::WorkStack > _3803_rs,monty::rc_ptr< ::mosek::fusion::WorkStack > _3804_ws,monty::rc_ptr< ::mosek::fusion::WorkStack > _3805_xs) ;
      virtual /* override */ std::string toString() ;
    }; // struct ExprMulParamDiagLeft;

    struct p_ExprMulParamDiagRight : public ::mosek::fusion::p_BaseExpression
    {
      ExprMulParamDiagRight * _pubthis;
      static mosek::fusion::p_ExprMulParamDiagRight* _get_impl(mosek::fusion::ExprMulParamDiagRight * _inst){ return static_cast< mosek::fusion::p_ExprMulParamDiagRight* >(mosek::fusion::p_BaseExpression::_get_impl(_inst)); }
      static mosek::fusion::p_ExprMulParamDiagRight * _get_impl(mosek::fusion::ExprMulParamDiagRight::t _inst) { return _get_impl(_inst.get()); }
      p_ExprMulParamDiagRight(ExprMulParamDiagRight * _pubthis);
      virtual ~p_ExprMulParamDiagRight() { /* std::cout << "~p_ExprMulParamDiagRight" << std::endl;*/ };
      monty::rc_ptr< ::mosek::fusion::Expression > e{};
      monty::rc_ptr< ::mosek::fusion::Parameter > p{};

      virtual void destroy();

      static ExprMulParamDiagRight::t _new_ExprMulParamDiagRight(monty::rc_ptr< ::mosek::fusion::Expression > _3920_e,monty::rc_ptr< ::mosek::fusion::Parameter > _3921_p);
      void _initialize(monty::rc_ptr< ::mosek::fusion::Expression > _3920_e,monty::rc_ptr< ::mosek::fusion::Parameter > _3921_p);
      virtual /* override */ void eval(monty::rc_ptr< ::mosek::fusion::WorkStack > _3922_rs,monty::rc_ptr< ::mosek::fusion::WorkStack > _3923_ws,monty::rc_ptr< ::mosek::fusion::WorkStack > _3924_xs) ;
      virtual /* override */ std::string toString() ;
    }; // struct ExprMulParamDiagRight;

    struct p_ExprDotParam : public ::mosek::fusion::p_BaseExpression
    {
      ExprDotParam * _pubthis;
      static mosek::fusion::p_ExprDotParam* _get_impl(mosek::fusion::ExprDotParam * _inst){ return static_cast< mosek::fusion::p_ExprDotParam* >(mosek::fusion::p_BaseExpression::_get_impl(_inst)); }
      static mosek::fusion::p_ExprDotParam * _get_impl(mosek::fusion::ExprDotParam::t _inst) { return _get_impl(_inst.get()); }
      p_ExprDotParam(ExprDotParam * _pubthis);
      virtual ~p_ExprDotParam() { /* std::cout << "~p_ExprDotParam" << std::endl;*/ };
      monty::rc_ptr< ::mosek::fusion::Expression > e{};
      monty::rc_ptr< ::mosek::fusion::Parameter > p{};

      virtual void destroy();

      static ExprDotParam::t _new_ExprDotParam(monty::rc_ptr< ::mosek::fusion::Parameter > _4038_p,monty::rc_ptr< ::mosek::fusion::Expression > _4039_e);
      void _initialize(monty::rc_ptr< ::mosek::fusion::Parameter > _4038_p,monty::rc_ptr< ::mosek::fusion::Expression > _4039_e);
      virtual /* override */ void eval(monty::rc_ptr< ::mosek::fusion::WorkStack > _4041_rs,monty::rc_ptr< ::mosek::fusion::WorkStack > _4042_ws,monty::rc_ptr< ::mosek::fusion::WorkStack > _4043_xs) ;
      virtual /* override */ std::string toString() ;
    }; // struct ExprDotParam;

    struct p_ExprMulParamElem : public ::mosek::fusion::p_BaseExpression
    {
      ExprMulParamElem * _pubthis;
      static mosek::fusion::p_ExprMulParamElem* _get_impl(mosek::fusion::ExprMulParamElem * _inst){ return static_cast< mosek::fusion::p_ExprMulParamElem* >(mosek::fusion::p_BaseExpression::_get_impl(_inst)); }
      static mosek::fusion::p_ExprMulParamElem * _get_impl(mosek::fusion::ExprMulParamElem::t _inst) { return _get_impl(_inst.get()); }
      p_ExprMulParamElem(ExprMulParamElem * _pubthis);
      virtual ~p_ExprMulParamElem() { /* std::cout << "~p_ExprMulParamElem" << std::endl;*/ };
      monty::rc_ptr< ::mosek::fusion::Expression > e{};
      monty::rc_ptr< ::mosek::fusion::Parameter > p{};

      virtual void destroy();

      static ExprMulParamElem::t _new_ExprMulParamElem(monty::rc_ptr< ::mosek::fusion::Parameter > _4101_p,monty::rc_ptr< ::mosek::fusion::Expression > _4102_e);
      void _initialize(monty::rc_ptr< ::mosek::fusion::Parameter > _4101_p,monty::rc_ptr< ::mosek::fusion::Expression > _4102_e);
      virtual /* override */ void eval(monty::rc_ptr< ::mosek::fusion::WorkStack > _4104_rs,monty::rc_ptr< ::mosek::fusion::WorkStack > _4105_ws,monty::rc_ptr< ::mosek::fusion::WorkStack > _4106_xs) ;
      virtual /* override */ std::string toString() ;
    }; // struct ExprMulParamElem;

    struct p_ExprMulParamRight : public ::mosek::fusion::p_BaseExpression
    {
      ExprMulParamRight * _pubthis;
      static mosek::fusion::p_ExprMulParamRight* _get_impl(mosek::fusion::ExprMulParamRight * _inst){ return static_cast< mosek::fusion::p_ExprMulParamRight* >(mosek::fusion::p_BaseExpression::_get_impl(_inst)); }
      static mosek::fusion::p_ExprMulParamRight * _get_impl(mosek::fusion::ExprMulParamRight::t _inst) { return _get_impl(_inst.get()); }
      p_ExprMulParamRight(ExprMulParamRight * _pubthis);
      virtual ~p_ExprMulParamRight() { /* std::cout << "~p_ExprMulParamRight" << std::endl;*/ };
      monty::rc_ptr< ::mosek::fusion::Expression > e{};
      monty::rc_ptr< ::mosek::fusion::Parameter > p{};

      virtual void destroy();

      static ExprMulParamRight::t _new_ExprMulParamRight(monty::rc_ptr< ::mosek::fusion::Expression > _4168_e,monty::rc_ptr< ::mosek::fusion::Parameter > _4169_p);
      void _initialize(monty::rc_ptr< ::mosek::fusion::Expression > _4168_e,monty::rc_ptr< ::mosek::fusion::Parameter > _4169_p);
      virtual /* override */ void eval(monty::rc_ptr< ::mosek::fusion::WorkStack > _4170_rs,monty::rc_ptr< ::mosek::fusion::WorkStack > _4171_ws,monty::rc_ptr< ::mosek::fusion::WorkStack > _4172_xs) ;
      virtual /* override */ std::string toString() ;
    }; // struct ExprMulParamRight;

    struct p_ExprMulParamLeft : public ::mosek::fusion::p_BaseExpression
    {
      ExprMulParamLeft * _pubthis;
      static mosek::fusion::p_ExprMulParamLeft* _get_impl(mosek::fusion::ExprMulParamLeft * _inst){ return static_cast< mosek::fusion::p_ExprMulParamLeft* >(mosek::fusion::p_BaseExpression::_get_impl(_inst)); }
      static mosek::fusion::p_ExprMulParamLeft * _get_impl(mosek::fusion::ExprMulParamLeft::t _inst) { return _get_impl(_inst.get()); }
      p_ExprMulParamLeft(ExprMulParamLeft * _pubthis);
      virtual ~p_ExprMulParamLeft() { /* std::cout << "~p_ExprMulParamLeft" << std::endl;*/ };
      monty::rc_ptr< ::mosek::fusion::Expression > e{};
      monty::rc_ptr< ::mosek::fusion::Parameter > p{};

      virtual void destroy();

      static ExprMulParamLeft::t _new_ExprMulParamLeft(monty::rc_ptr< ::mosek::fusion::Parameter > _4272_p,monty::rc_ptr< ::mosek::fusion::Expression > _4273_e);
      void _initialize(monty::rc_ptr< ::mosek::fusion::Parameter > _4272_p,monty::rc_ptr< ::mosek::fusion::Expression > _4273_e);
      virtual /* override */ void eval(monty::rc_ptr< ::mosek::fusion::WorkStack > _4274_rs,monty::rc_ptr< ::mosek::fusion::WorkStack > _4275_ws,monty::rc_ptr< ::mosek::fusion::WorkStack > _4276_xs) ;
      virtual /* override */ std::string toString() ;
    }; // struct ExprMulParamLeft;

    struct p_ExprOptimizeCode : public ::mosek::fusion::p_BaseExpression
    {
      ExprOptimizeCode * _pubthis;
      static mosek::fusion::p_ExprOptimizeCode* _get_impl(mosek::fusion::ExprOptimizeCode * _inst){ return static_cast< mosek::fusion::p_ExprOptimizeCode* >(mosek::fusion::p_BaseExpression::_get_impl(_inst)); }
      static mosek::fusion::p_ExprOptimizeCode * _get_impl(mosek::fusion::ExprOptimizeCode::t _inst) { return _get_impl(_inst.get()); }
      p_ExprOptimizeCode(ExprOptimizeCode * _pubthis);
      virtual ~p_ExprOptimizeCode() { /* std::cout << "~p_ExprOptimizeCode" << std::endl;*/ };
      monty::rc_ptr< ::mosek::fusion::Expression > expr{};

      virtual void destroy();

      static ExprOptimizeCode::t _new_ExprOptimizeCode(monty::rc_ptr< ::mosek::fusion::Expression > _4587_expr);
      void _initialize(monty::rc_ptr< ::mosek::fusion::Expression > _4587_expr);
      static  void compress_code(monty::rc_ptr< ::mosek::fusion::WorkStack > _4588_xs,int32_t _4589_n,std::shared_ptr< monty::ndarray< int32_t,1 > > _4590_code,int32_t _4591_code_base,std::shared_ptr< monty::ndarray< int32_t,1 > > _4592_ptr,int32_t _4593_ptr_base,std::shared_ptr< monty::ndarray< double,1 > > _4594_fixterm,int32_t _4595_fixterm_base,std::shared_ptr< monty::ndarray< double,1 > > _4596_code_consts,int32_t _4597_code_consts_base,int32_t _4598_target_code_base,int32_t _4599_target_const_base,int32_t _4600_target_ptr_base);
      virtual /* override */ void eval(monty::rc_ptr< ::mosek::fusion::WorkStack > _4653_rs,monty::rc_ptr< ::mosek::fusion::WorkStack > _4654_ws,monty::rc_ptr< ::mosek::fusion::WorkStack > _4655_xs) ;
      virtual /* override */ std::string toString() ;
    }; // struct ExprOptimizeCode;

    struct p_ExprCompress : public ::mosek::fusion::p_BaseExpression
    {
      ExprCompress * _pubthis;
      static mosek::fusion::p_ExprCompress* _get_impl(mosek::fusion::ExprCompress * _inst){ return static_cast< mosek::fusion::p_ExprCompress* >(mosek::fusion::p_BaseExpression::_get_impl(_inst)); }
      static mosek::fusion::p_ExprCompress * _get_impl(mosek::fusion::ExprCompress::t _inst) { return _get_impl(_inst.get()); }
      p_ExprCompress(ExprCompress * _pubthis);
      virtual ~p_ExprCompress() { /* std::cout << "~p_ExprCompress" << std::endl;*/ };
      monty::rc_ptr< ::mosek::fusion::Expression > expr{};

      virtual void destroy();

      static ExprCompress::t _new_ExprCompress(monty::rc_ptr< ::mosek::fusion::Expression > _4722_expr);
      void _initialize(monty::rc_ptr< ::mosek::fusion::Expression > _4722_expr);
      static  void arg_sort(monty::rc_ptr< ::mosek::fusion::WorkStack > _4723_ws,monty::rc_ptr< ::mosek::fusion::WorkStack > _4724_xs,int32_t _4725_perm,int32_t _4726_nelem,int32_t _4727_nnz,int32_t _4728_ptr,int32_t _4729_nidxs);
      static  void merge_sort(int32_t _4765_origperm1,int32_t _4766_origperm2,int32_t _4767_nelem,int32_t _4768_nnz,int32_t _4769_ptr_base,int32_t _4770_nidxs_base,std::shared_ptr< monty::ndarray< int32_t,1 > > _4771_wi32,std::shared_ptr< monty::ndarray< int64_t,1 > > _4772_wi64);
      virtual /* override */ void eval(monty::rc_ptr< ::mosek::fusion::WorkStack > _4795_rs,monty::rc_ptr< ::mosek::fusion::WorkStack > _4796_ws,monty::rc_ptr< ::mosek::fusion::WorkStack > _4797_xs) ;
      virtual /* override */ std::string toString() ;
    }; // struct ExprCompress;

    struct p_ExprConst : public ::mosek::fusion::p_BaseExpression
    {
      ExprConst * _pubthis;
      static mosek::fusion::p_ExprConst* _get_impl(mosek::fusion::ExprConst * _inst){ return static_cast< mosek::fusion::p_ExprConst* >(mosek::fusion::p_BaseExpression::_get_impl(_inst)); }
      static mosek::fusion::p_ExprConst * _get_impl(mosek::fusion::ExprConst::t _inst) { return _get_impl(_inst.get()); }
      p_ExprConst(ExprConst * _pubthis);
      virtual ~p_ExprConst() { /* std::cout << "~p_ExprConst" << std::endl;*/ };
      std::shared_ptr< monty::ndarray< int64_t,1 > > sparsity{};
      std::shared_ptr< monty::ndarray< double,1 > > bfix{};

      virtual void destroy();

      static ExprConst::t _new_ExprConst(std::shared_ptr< monty::ndarray< int32_t,1 > > _4883_shape,std::shared_ptr< monty::ndarray< int64_t,1 > > _4884_sparsity,std::shared_ptr< monty::ndarray< double,1 > > _4885_bfix);
      void _initialize(std::shared_ptr< monty::ndarray< int32_t,1 > > _4883_shape,std::shared_ptr< monty::ndarray< int64_t,1 > > _4884_sparsity,std::shared_ptr< monty::ndarray< double,1 > > _4885_bfix);
      static ExprConst::t _new_ExprConst(std::shared_ptr< monty::ndarray< int32_t,1 > > _4886_shape,std::shared_ptr< monty::ndarray< int64_t,1 > > _4887_sparsity,double _4888_bfix);
      void _initialize(std::shared_ptr< monty::ndarray< int32_t,1 > > _4886_shape,std::shared_ptr< monty::ndarray< int64_t,1 > > _4887_sparsity,double _4888_bfix);
      virtual /* override */ void eval(monty::rc_ptr< ::mosek::fusion::WorkStack > _4891_rs,monty::rc_ptr< ::mosek::fusion::WorkStack > _4892_ws,monty::rc_ptr< ::mosek::fusion::WorkStack > _4893_xs) ;
      static  void validate(std::shared_ptr< monty::ndarray< int32_t,1 > > _4912_shape,std::shared_ptr< monty::ndarray< double,1 > > _4913_bfix,std::shared_ptr< monty::ndarray< int64_t,1 > > _4914_sparsity);
      virtual /* override */ std::string toString() ;
    }; // struct ExprConst;

    struct p_ExprPick : public ::mosek::fusion::p_BaseExpression
    {
      ExprPick * _pubthis;
      static mosek::fusion::p_ExprPick* _get_impl(mosek::fusion::ExprPick * _inst){ return static_cast< mosek::fusion::p_ExprPick* >(mosek::fusion::p_BaseExpression::_get_impl(_inst)); }
      static mosek::fusion::p_ExprPick * _get_impl(mosek::fusion::ExprPick::t _inst) { return _get_impl(_inst.get()); }
      p_ExprPick(ExprPick * _pubthis);
      virtual ~p_ExprPick() { /* std::cout << "~p_ExprPick" << std::endl;*/ };
      std::shared_ptr< monty::ndarray< int64_t,1 > > idxs{};
      monty::rc_ptr< ::mosek::fusion::Expression > expr{};

      virtual void destroy();

      static ExprPick::t _new_ExprPick(monty::rc_ptr< ::mosek::fusion::Expression > _4918_expr,std::shared_ptr< monty::ndarray< int32_t,2 > > _4919_idxs);
      void _initialize(monty::rc_ptr< ::mosek::fusion::Expression > _4918_expr,std::shared_ptr< monty::ndarray< int32_t,2 > > _4919_idxs);
      static ExprPick::t _new_ExprPick(monty::rc_ptr< ::mosek::fusion::Expression > _4931_expr,std::shared_ptr< monty::ndarray< int64_t,1 > > _4932_idxs);
      void _initialize(monty::rc_ptr< ::mosek::fusion::Expression > _4931_expr,std::shared_ptr< monty::ndarray< int64_t,1 > > _4932_idxs);
      virtual /* override */ void eval(monty::rc_ptr< ::mosek::fusion::WorkStack > _4937_rs,monty::rc_ptr< ::mosek::fusion::WorkStack > _4938_ws,monty::rc_ptr< ::mosek::fusion::WorkStack > _4939_xs) ;
      virtual /* override */ std::string toString() ;
    }; // struct ExprPick;

    struct p_ExprSlice : public ::mosek::fusion::p_BaseExpression
    {
      ExprSlice * _pubthis;
      static mosek::fusion::p_ExprSlice* _get_impl(mosek::fusion::ExprSlice * _inst){ return static_cast< mosek::fusion::p_ExprSlice* >(mosek::fusion::p_BaseExpression::_get_impl(_inst)); }
      static mosek::fusion::p_ExprSlice * _get_impl(mosek::fusion::ExprSlice::t _inst) { return _get_impl(_inst.get()); }
      p_ExprSlice(ExprSlice * _pubthis);
      virtual ~p_ExprSlice() { /* std::cout << "~p_ExprSlice" << std::endl;*/ };
      std::shared_ptr< monty::ndarray< int32_t,1 > > last{};
      std::shared_ptr< monty::ndarray< int32_t,1 > > first{};
      monty::rc_ptr< ::mosek::fusion::Expression > expr{};

      virtual void destroy();

      static ExprSlice::t _new_ExprSlice(monty::rc_ptr< ::mosek::fusion::Expression > _5004_expr,std::shared_ptr< monty::ndarray< int32_t,1 > > _5005_first,std::shared_ptr< monty::ndarray< int32_t,1 > > _5006_last);
      void _initialize(monty::rc_ptr< ::mosek::fusion::Expression > _5004_expr,std::shared_ptr< monty::ndarray< int32_t,1 > > _5005_first,std::shared_ptr< monty::ndarray< int32_t,1 > > _5006_last);
      virtual /* override */ void eval(monty::rc_ptr< ::mosek::fusion::WorkStack > _5007_rs,monty::rc_ptr< ::mosek::fusion::WorkStack > _5008_ws,monty::rc_ptr< ::mosek::fusion::WorkStack > _5009_xs) ;
      static  std::shared_ptr< monty::ndarray< int32_t,1 > > makeShape(std::shared_ptr< monty::ndarray< int32_t,1 > > _5074_shape,std::shared_ptr< monty::ndarray< int32_t,1 > > _5075_first,std::shared_ptr< monty::ndarray< int32_t,1 > > _5076_last);
      virtual /* override */ std::string toString() ;
    }; // struct ExprSlice;

    struct p_ExprPermuteDims : public ::mosek::fusion::p_BaseExpression
    {
      ExprPermuteDims * _pubthis;
      static mosek::fusion::p_ExprPermuteDims* _get_impl(mosek::fusion::ExprPermuteDims * _inst){ return static_cast< mosek::fusion::p_ExprPermuteDims* >(mosek::fusion::p_BaseExpression::_get_impl(_inst)); }
      static mosek::fusion::p_ExprPermuteDims * _get_impl(mosek::fusion::ExprPermuteDims::t _inst) { return _get_impl(_inst.get()); }
      p_ExprPermuteDims(ExprPermuteDims * _pubthis);
      virtual ~p_ExprPermuteDims() { /* std::cout << "~p_ExprPermuteDims" << std::endl;*/ };
      std::shared_ptr< monty::ndarray< int32_t,1 > > dperm{};
      monty::rc_ptr< ::mosek::fusion::Expression > expr{};

      virtual void destroy();

      static ExprPermuteDims::t _new_ExprPermuteDims(std::shared_ptr< monty::ndarray< int32_t,1 > > _5081_perm,monty::rc_ptr< ::mosek::fusion::Expression > _5082_expr);
      void _initialize(std::shared_ptr< monty::ndarray< int32_t,1 > > _5081_perm,monty::rc_ptr< ::mosek::fusion::Expression > _5082_expr);
      static ExprPermuteDims::t _new_ExprPermuteDims(std::shared_ptr< monty::ndarray< int32_t,1 > > _5088_perm,monty::rc_ptr< ::mosek::fusion::Expression > _5089_expr,int32_t _5090_validated);
      void _initialize(std::shared_ptr< monty::ndarray< int32_t,1 > > _5088_perm,monty::rc_ptr< ::mosek::fusion::Expression > _5089_expr,int32_t _5090_validated);
      virtual /* override */ void eval(monty::rc_ptr< ::mosek::fusion::WorkStack > _5091_rs,monty::rc_ptr< ::mosek::fusion::WorkStack > _5092_ws,monty::rc_ptr< ::mosek::fusion::WorkStack > _5093_xs) ;
      static  std::shared_ptr< monty::ndarray< int32_t,1 > > computeshape(std::shared_ptr< monty::ndarray< int32_t,1 > > _5147_perm,std::shared_ptr< monty::ndarray< int32_t,1 > > _5148_shape);
    }; // struct ExprPermuteDims;

    struct p_ExprTranspose : public ::mosek::fusion::p_BaseExpression
    {
      ExprTranspose * _pubthis;
      static mosek::fusion::p_ExprTranspose* _get_impl(mosek::fusion::ExprTranspose * _inst){ return static_cast< mosek::fusion::p_ExprTranspose* >(mosek::fusion::p_BaseExpression::_get_impl(_inst)); }
      static mosek::fusion::p_ExprTranspose * _get_impl(mosek::fusion::ExprTranspose::t _inst) { return _get_impl(_inst.get()); }
      p_ExprTranspose(ExprTranspose * _pubthis);
      virtual ~p_ExprTranspose() { /* std::cout << "~p_ExprTranspose" << std::endl;*/ };
      monty::rc_ptr< ::mosek::fusion::Expression > expr{};

      virtual void destroy();

      static ExprTranspose::t _new_ExprTranspose(monty::rc_ptr< ::mosek::fusion::Expression > _5150_expr);
      void _initialize(monty::rc_ptr< ::mosek::fusion::Expression > _5150_expr);
      virtual /* override */ void eval(monty::rc_ptr< ::mosek::fusion::WorkStack > _5151_rs,monty::rc_ptr< ::mosek::fusion::WorkStack > _5152_ws,monty::rc_ptr< ::mosek::fusion::WorkStack > _5153_xs) ;
      virtual /* override */ std::string toString() ;
      static  std::shared_ptr< monty::ndarray< int32_t,1 > > transposeShape(std::shared_ptr< monty::ndarray< int32_t,1 > > _5206_shape);
    }; // struct ExprTranspose;

    struct p_ExprRepeat : public ::mosek::fusion::p_BaseExpression
    {
      ExprRepeat * _pubthis;
      static mosek::fusion::p_ExprRepeat* _get_impl(mosek::fusion::ExprRepeat * _inst){ return static_cast< mosek::fusion::p_ExprRepeat* >(mosek::fusion::p_BaseExpression::_get_impl(_inst)); }
      static mosek::fusion::p_ExprRepeat * _get_impl(mosek::fusion::ExprRepeat::t _inst) { return _get_impl(_inst.get()); }
      p_ExprRepeat(ExprRepeat * _pubthis);
      virtual ~p_ExprRepeat() { /* std::cout << "~p_ExprRepeat" << std::endl;*/ };
      int32_t n{};
      int32_t dim{};
      monty::rc_ptr< ::mosek::fusion::Expression > expr{};

      virtual void destroy();

      static ExprRepeat::t _new_ExprRepeat(monty::rc_ptr< ::mosek::fusion::Expression > _5207_expr,int32_t _5208_dim,int32_t _5209_n);
      void _initialize(monty::rc_ptr< ::mosek::fusion::Expression > _5207_expr,int32_t _5208_dim,int32_t _5209_n);
      virtual /* override */ void eval(monty::rc_ptr< ::mosek::fusion::WorkStack > _5210_rs,monty::rc_ptr< ::mosek::fusion::WorkStack > _5211_ws,monty::rc_ptr< ::mosek::fusion::WorkStack > _5212_xs) ;
      static  std::shared_ptr< monty::ndarray< int32_t,1 > > getshape(monty::rc_ptr< ::mosek::fusion::Expression > _5277_e,int32_t _5278_dim,int32_t _5279_n);
      virtual /* override */ std::string toString() ;
    }; // struct ExprRepeat;

    struct p_ExprStack : public ::mosek::fusion::p_BaseExpression
    {
      ExprStack * _pubthis;
      static mosek::fusion::p_ExprStack* _get_impl(mosek::fusion::ExprStack * _inst){ return static_cast< mosek::fusion::p_ExprStack* >(mosek::fusion::p_BaseExpression::_get_impl(_inst)); }
      static mosek::fusion::p_ExprStack * _get_impl(mosek::fusion::ExprStack::t _inst) { return _get_impl(_inst.get()); }
      p_ExprStack(ExprStack * _pubthis);
      virtual ~p_ExprStack() { /* std::cout << "~p_ExprStack" << std::endl;*/ };
      int32_t dim{};
      std::shared_ptr< monty::ndarray< monty::rc_ptr< ::mosek::fusion::Expression >,1 > > exprs{};

      virtual void destroy();

      static ExprStack::t _new_ExprStack(std::shared_ptr< monty::ndarray< monty::rc_ptr< ::mosek::fusion::Expression >,1 > > _5284_exprs,int32_t _5285_dim);
      void _initialize(std::shared_ptr< monty::ndarray< monty::rc_ptr< ::mosek::fusion::Expression >,1 > > _5284_exprs,int32_t _5285_dim);
      virtual /* override */ void eval(monty::rc_ptr< ::mosek::fusion::WorkStack > _5287_rs,monty::rc_ptr< ::mosek::fusion::WorkStack > _5288_ws,monty::rc_ptr< ::mosek::fusion::WorkStack > _5289_xs) ;
      static  std::shared_ptr< monty::ndarray< int32_t,1 > > getshape(std::shared_ptr< monty::ndarray< monty::rc_ptr< ::mosek::fusion::Expression >,1 > > _5433_es,int32_t _5434_dim);
      virtual /* override */ std::string toString() ;
    }; // struct ExprStack;

    struct p_ExprInner : public ::mosek::fusion::p_BaseExpression
    {
      ExprInner * _pubthis;
      static mosek::fusion::p_ExprInner* _get_impl(mosek::fusion::ExprInner * _inst){ return static_cast< mosek::fusion::p_ExprInner* >(mosek::fusion::p_BaseExpression::_get_impl(_inst)); }
      static mosek::fusion::p_ExprInner * _get_impl(mosek::fusion::ExprInner::t _inst) { return _get_impl(_inst.get()); }
      p_ExprInner(ExprInner * _pubthis);
      virtual ~p_ExprInner() { /* std::cout << "~p_ExprInner" << std::endl;*/ };
      std::shared_ptr< monty::ndarray< double,1 > > vcof{};
      std::shared_ptr< monty::ndarray< int64_t,1 > > vsub{};
      monty::rc_ptr< ::mosek::fusion::Expression > expr{};

      virtual void destroy();

      static ExprInner::t _new_ExprInner(monty::rc_ptr< ::mosek::fusion::Expression > _5448_expr3,std::shared_ptr< monty::ndarray< int64_t,1 > > _5449_vsub3,std::shared_ptr< monty::ndarray< double,1 > > _5450_vcof3);
      void _initialize(monty::rc_ptr< ::mosek::fusion::Expression > _5448_expr3,std::shared_ptr< monty::ndarray< int64_t,1 > > _5449_vsub3,std::shared_ptr< monty::ndarray< double,1 > > _5450_vcof3);
      static ExprInner::t _new_ExprInner(monty::rc_ptr< ::mosek::fusion::Expression > _5456_expr2,std::shared_ptr< monty::ndarray< double,1 > > _5457_vcof2);
      void _initialize(monty::rc_ptr< ::mosek::fusion::Expression > _5456_expr2,std::shared_ptr< monty::ndarray< double,1 > > _5457_vcof2);
      static ExprInner::t _new_ExprInner(monty::rc_ptr< ::mosek::fusion::Expression > _5459_expr1,std::shared_ptr< monty::ndarray< int32_t,2 > > _5460_vsub1,std::shared_ptr< monty::ndarray< double,1 > > _5461_vcof1);
      void _initialize(monty::rc_ptr< ::mosek::fusion::Expression > _5459_expr1,std::shared_ptr< monty::ndarray< int32_t,2 > > _5460_vsub1,std::shared_ptr< monty::ndarray< double,1 > > _5461_vcof1);
      virtual /* override */ void eval(monty::rc_ptr< ::mosek::fusion::WorkStack > _5462_rs,monty::rc_ptr< ::mosek::fusion::WorkStack > _5463_ws,monty::rc_ptr< ::mosek::fusion::WorkStack > _5464_xs) ;
      static  std::shared_ptr< monty::ndarray< int64_t,1 > > range(int32_t _5508_n);
      static  std::shared_ptr< monty::ndarray< int64_t,1 > > convert(std::shared_ptr< monty::ndarray< int32_t,1 > > _5510_shape,std::shared_ptr< monty::ndarray< int32_t,2 > > _5511_vsub);
      virtual /* override */ std::string toString() ;
    }; // struct ExprInner;

    struct p_ExprMulDiagRight : public ::mosek::fusion::p_BaseExpression
    {
      ExprMulDiagRight * _pubthis;
      static mosek::fusion::p_ExprMulDiagRight* _get_impl(mosek::fusion::ExprMulDiagRight * _inst){ return static_cast< mosek::fusion::p_ExprMulDiagRight* >(mosek::fusion::p_BaseExpression::_get_impl(_inst)); }
      static mosek::fusion::p_ExprMulDiagRight * _get_impl(mosek::fusion::ExprMulDiagRight::t _inst) { return _get_impl(_inst.get()); }
      p_ExprMulDiagRight(ExprMulDiagRight * _pubthis);
      virtual ~p_ExprMulDiagRight() { /* std::cout << "~p_ExprMulDiagRight" << std::endl;*/ };
      monty::rc_ptr< ::mosek::fusion::Expression > expr{};
      std::shared_ptr< monty::ndarray< double,1 > > mval{};
      std::shared_ptr< monty::ndarray< int32_t,1 > > msubj{};
      std::shared_ptr< monty::ndarray< int32_t,1 > > msubi{};
      int32_t mdim1{};
      int32_t mdim0{};

      virtual void destroy();

      static ExprMulDiagRight::t _new_ExprMulDiagRight(int32_t _5518_mdim0,int32_t _5519_mdim1,std::shared_ptr< monty::ndarray< int32_t,1 > > _5520_msubi,std::shared_ptr< monty::ndarray< int32_t,1 > > _5521_msubj,std::shared_ptr< monty::ndarray< double,1 > > _5522_mval,monty::rc_ptr< ::mosek::fusion::Expression > _5523_expr,int32_t _5524_validated);
      void _initialize(int32_t _5518_mdim0,int32_t _5519_mdim1,std::shared_ptr< monty::ndarray< int32_t,1 > > _5520_msubi,std::shared_ptr< monty::ndarray< int32_t,1 > > _5521_msubj,std::shared_ptr< monty::ndarray< double,1 > > _5522_mval,monty::rc_ptr< ::mosek::fusion::Expression > _5523_expr,int32_t _5524_validated);
      static ExprMulDiagRight::t _new_ExprMulDiagRight(int32_t _5525_mdim0,int32_t _5526_mdim1,std::shared_ptr< monty::ndarray< int32_t,1 > > _5527_msubi,std::shared_ptr< monty::ndarray< int32_t,1 > > _5528_msubj,std::shared_ptr< monty::ndarray< double,1 > > _5529_mval,monty::rc_ptr< ::mosek::fusion::Expression > _5530_expr);
      void _initialize(int32_t _5525_mdim0,int32_t _5526_mdim1,std::shared_ptr< monty::ndarray< int32_t,1 > > _5527_msubi,std::shared_ptr< monty::ndarray< int32_t,1 > > _5528_msubj,std::shared_ptr< monty::ndarray< double,1 > > _5529_mval,monty::rc_ptr< ::mosek::fusion::Expression > _5530_expr);
      virtual /* override */ void eval(monty::rc_ptr< ::mosek::fusion::WorkStack > _5531_rs,monty::rc_ptr< ::mosek::fusion::WorkStack > _5532_ws,monty::rc_ptr< ::mosek::fusion::WorkStack > _5533_xs) ;
      static  int32_t validate(int32_t _5612_mdim0,int32_t _5613_mdim1,std::shared_ptr< monty::ndarray< int32_t,1 > > _5614_msubi,std::shared_ptr< monty::ndarray< int32_t,1 > > _5615_msubj,std::shared_ptr< monty::ndarray< double,1 > > _5616_mval,monty::rc_ptr< ::mosek::fusion::Expression > _5617_expr);
      virtual /* override */ std::string toString() ;
    }; // struct ExprMulDiagRight;

    struct p_ExprMulDiagLeft : public ::mosek::fusion::p_BaseExpression
    {
      ExprMulDiagLeft * _pubthis;
      static mosek::fusion::p_ExprMulDiagLeft* _get_impl(mosek::fusion::ExprMulDiagLeft * _inst){ return static_cast< mosek::fusion::p_ExprMulDiagLeft* >(mosek::fusion::p_BaseExpression::_get_impl(_inst)); }
      static mosek::fusion::p_ExprMulDiagLeft * _get_impl(mosek::fusion::ExprMulDiagLeft::t _inst) { return _get_impl(_inst.get()); }
      p_ExprMulDiagLeft(ExprMulDiagLeft * _pubthis);
      virtual ~p_ExprMulDiagLeft() { /* std::cout << "~p_ExprMulDiagLeft" << std::endl;*/ };
      monty::rc_ptr< ::mosek::fusion::Expression > expr{};
      std::shared_ptr< monty::ndarray< double,1 > > mval{};
      std::shared_ptr< monty::ndarray< int32_t,1 > > msubj{};
      std::shared_ptr< monty::ndarray< int32_t,1 > > msubi{};
      int32_t mdim1{};
      int32_t mdim0{};

      virtual void destroy();

      static ExprMulDiagLeft::t _new_ExprMulDiagLeft(int32_t _5626_mdim0,int32_t _5627_mdim1,std::shared_ptr< monty::ndarray< int32_t,1 > > _5628_msubi,std::shared_ptr< monty::ndarray< int32_t,1 > > _5629_msubj,std::shared_ptr< monty::ndarray< double,1 > > _5630_mval,monty::rc_ptr< ::mosek::fusion::Expression > _5631_expr,int32_t _5632_validated);
      void _initialize(int32_t _5626_mdim0,int32_t _5627_mdim1,std::shared_ptr< monty::ndarray< int32_t,1 > > _5628_msubi,std::shared_ptr< monty::ndarray< int32_t,1 > > _5629_msubj,std::shared_ptr< monty::ndarray< double,1 > > _5630_mval,monty::rc_ptr< ::mosek::fusion::Expression > _5631_expr,int32_t _5632_validated);
      static ExprMulDiagLeft::t _new_ExprMulDiagLeft(int32_t _5633_mdim0,int32_t _5634_mdim1,std::shared_ptr< monty::ndarray< int32_t,1 > > _5635_msubi,std::shared_ptr< monty::ndarray< int32_t,1 > > _5636_msubj,std::shared_ptr< monty::ndarray< double,1 > > _5637_mval,monty::rc_ptr< ::mosek::fusion::Expression > _5638_expr);
      void _initialize(int32_t _5633_mdim0,int32_t _5634_mdim1,std::shared_ptr< monty::ndarray< int32_t,1 > > _5635_msubi,std::shared_ptr< monty::ndarray< int32_t,1 > > _5636_msubj,std::shared_ptr< monty::ndarray< double,1 > > _5637_mval,monty::rc_ptr< ::mosek::fusion::Expression > _5638_expr);
      virtual /* override */ void eval(monty::rc_ptr< ::mosek::fusion::WorkStack > _5639_rs,monty::rc_ptr< ::mosek::fusion::WorkStack > _5640_ws,monty::rc_ptr< ::mosek::fusion::WorkStack > _5641_xs) ;
      static  int32_t validate(int32_t _5739_mdim0,int32_t _5740_mdim1,std::shared_ptr< monty::ndarray< int32_t,1 > > _5741_msubi,std::shared_ptr< monty::ndarray< int32_t,1 > > _5742_msubj,std::shared_ptr< monty::ndarray< double,1 > > _5743_mval,monty::rc_ptr< ::mosek::fusion::Expression > _5744_expr);
      virtual /* override */ std::string toString() ;
    }; // struct ExprMulDiagLeft;

    struct p_ExprMulElement : public ::mosek::fusion::p_BaseExpression
    {
      ExprMulElement * _pubthis;
      static mosek::fusion::p_ExprMulElement* _get_impl(mosek::fusion::ExprMulElement * _inst){ return static_cast< mosek::fusion::p_ExprMulElement* >(mosek::fusion::p_BaseExpression::_get_impl(_inst)); }
      static mosek::fusion::p_ExprMulElement * _get_impl(mosek::fusion::ExprMulElement::t _inst) { return _get_impl(_inst.get()); }
      p_ExprMulElement(ExprMulElement * _pubthis);
      virtual ~p_ExprMulElement() { /* std::cout << "~p_ExprMulElement" << std::endl;*/ };
      monty::rc_ptr< ::mosek::fusion::Expression > expr{};
      std::shared_ptr< monty::ndarray< int64_t,1 > > msp{};
      std::shared_ptr< monty::ndarray< double,1 > > mcof{};

      virtual void destroy();

      static ExprMulElement::t _new_ExprMulElement(std::shared_ptr< monty::ndarray< double,1 > > _5753_mcof,std::shared_ptr< monty::ndarray< int64_t,1 > > _5754_msp,monty::rc_ptr< ::mosek::fusion::Expression > _5755_expr);
      void _initialize(std::shared_ptr< monty::ndarray< double,1 > > _5753_mcof,std::shared_ptr< monty::ndarray< int64_t,1 > > _5754_msp,monty::rc_ptr< ::mosek::fusion::Expression > _5755_expr);
      static ExprMulElement::t _new_ExprMulElement(std::shared_ptr< monty::ndarray< double,1 > > _5762_cof,std::shared_ptr< monty::ndarray< int64_t,1 > > _5763_msp,monty::rc_ptr< ::mosek::fusion::Expression > _5764_expr,int32_t _5765_validated);
      void _initialize(std::shared_ptr< monty::ndarray< double,1 > > _5762_cof,std::shared_ptr< monty::ndarray< int64_t,1 > > _5763_msp,monty::rc_ptr< ::mosek::fusion::Expression > _5764_expr,int32_t _5765_validated);
      virtual /* override */ void eval(monty::rc_ptr< ::mosek::fusion::WorkStack > _5766_rs,monty::rc_ptr< ::mosek::fusion::WorkStack > _5767_ws,monty::rc_ptr< ::mosek::fusion::WorkStack > _5768_xs) ;
      virtual /* override */ std::string toString() ;
    }; // struct ExprMulElement;

    struct p_ExprMulScalarConst : public ::mosek::fusion::p_BaseExpression
    {
      ExprMulScalarConst * _pubthis;
      static mosek::fusion::p_ExprMulScalarConst* _get_impl(mosek::fusion::ExprMulScalarConst * _inst){ return static_cast< mosek::fusion::p_ExprMulScalarConst* >(mosek::fusion::p_BaseExpression::_get_impl(_inst)); }
      static mosek::fusion::p_ExprMulScalarConst * _get_impl(mosek::fusion::ExprMulScalarConst::t _inst) { return _get_impl(_inst.get()); }
      p_ExprMulScalarConst(ExprMulScalarConst * _pubthis);
      virtual ~p_ExprMulScalarConst() { /* std::cout << "~p_ExprMulScalarConst" << std::endl;*/ };
      monty::rc_ptr< ::mosek::fusion::Expression > expr{};
      double c{};

      virtual void destroy();

      static ExprMulScalarConst::t _new_ExprMulScalarConst(double _5826_c,monty::rc_ptr< ::mosek::fusion::Expression > _5827_expr);
      void _initialize(double _5826_c,monty::rc_ptr< ::mosek::fusion::Expression > _5827_expr);
      virtual /* override */ void eval(monty::rc_ptr< ::mosek::fusion::WorkStack > _5828_rs,monty::rc_ptr< ::mosek::fusion::WorkStack > _5829_ws,monty::rc_ptr< ::mosek::fusion::WorkStack > _5830_xs) ;
      virtual /* override */ std::string toString() ;
    }; // struct ExprMulScalarConst;

    struct p_ExprScalarMul : public ::mosek::fusion::p_BaseExpression
    {
      ExprScalarMul * _pubthis;
      static mosek::fusion::p_ExprScalarMul* _get_impl(mosek::fusion::ExprScalarMul * _inst){ return static_cast< mosek::fusion::p_ExprScalarMul* >(mosek::fusion::p_BaseExpression::_get_impl(_inst)); }
      static mosek::fusion::p_ExprScalarMul * _get_impl(mosek::fusion::ExprScalarMul::t _inst) { return _get_impl(_inst.get()); }
      p_ExprScalarMul(ExprScalarMul * _pubthis);
      virtual ~p_ExprScalarMul() { /* std::cout << "~p_ExprScalarMul" << std::endl;*/ };
      monty::rc_ptr< ::mosek::fusion::Expression > expr{};
      std::shared_ptr< monty::ndarray< double,1 > > mval{};
      std::shared_ptr< monty::ndarray< int32_t,1 > > msubj{};
      std::shared_ptr< monty::ndarray< int32_t,1 > > msubi{};
      int32_t mdim1{};
      int32_t mdim0{};

      virtual void destroy();

      static ExprScalarMul::t _new_ExprScalarMul(int32_t _5868_mdim0,int32_t _5869_mdim1,std::shared_ptr< monty::ndarray< int32_t,1 > > _5870_msubi,std::shared_ptr< monty::ndarray< int32_t,1 > > _5871_msubj,std::shared_ptr< monty::ndarray< double,1 > > _5872_mval,monty::rc_ptr< ::mosek::fusion::Expression > _5873_expr,int32_t _5874_validated);
      void _initialize(int32_t _5868_mdim0,int32_t _5869_mdim1,std::shared_ptr< monty::ndarray< int32_t,1 > > _5870_msubi,std::shared_ptr< monty::ndarray< int32_t,1 > > _5871_msubj,std::shared_ptr< monty::ndarray< double,1 > > _5872_mval,monty::rc_ptr< ::mosek::fusion::Expression > _5873_expr,int32_t _5874_validated);
      static ExprScalarMul::t _new_ExprScalarMul(int32_t _5875_mdim0,int32_t _5876_mdim1,std::shared_ptr< monty::ndarray< int32_t,1 > > _5877_msubi,std::shared_ptr< monty::ndarray< int32_t,1 > > _5878_msubj,std::shared_ptr< monty::ndarray< double,1 > > _5879_mval,monty::rc_ptr< ::mosek::fusion::Expression > _5880_expr);
      void _initialize(int32_t _5875_mdim0,int32_t _5876_mdim1,std::shared_ptr< monty::ndarray< int32_t,1 > > _5877_msubi,std::shared_ptr< monty::ndarray< int32_t,1 > > _5878_msubj,std::shared_ptr< monty::ndarray< double,1 > > _5879_mval,monty::rc_ptr< ::mosek::fusion::Expression > _5880_expr);
      virtual /* override */ void eval(monty::rc_ptr< ::mosek::fusion::WorkStack > _5881_rs,monty::rc_ptr< ::mosek::fusion::WorkStack > _5882_ws,monty::rc_ptr< ::mosek::fusion::WorkStack > _5883_xs) ;
      static  int32_t validate(int32_t _5919_mdim0,int32_t _5920_mdim1,std::shared_ptr< monty::ndarray< int32_t,1 > > _5921_msubi,std::shared_ptr< monty::ndarray< int32_t,1 > > _5922_msubj,std::shared_ptr< monty::ndarray< double,1 > > _5923_mval,monty::rc_ptr< ::mosek::fusion::Expression > _5924_expr);
      virtual /* override */ std::string toString() ;
    }; // struct ExprScalarMul;

    struct p_ExprMulRight : public ::mosek::fusion::p_BaseExpression
    {
      ExprMulRight * _pubthis;
      static mosek::fusion::p_ExprMulRight* _get_impl(mosek::fusion::ExprMulRight * _inst){ return static_cast< mosek::fusion::p_ExprMulRight* >(mosek::fusion::p_BaseExpression::_get_impl(_inst)); }
      static mosek::fusion::p_ExprMulRight * _get_impl(mosek::fusion::ExprMulRight::t _inst) { return _get_impl(_inst.get()); }
      p_ExprMulRight(ExprMulRight * _pubthis);
      virtual ~p_ExprMulRight() { /* std::cout << "~p_ExprMulRight" << std::endl;*/ };
      monty::rc_ptr< ::mosek::fusion::Expression > expr{};
      std::shared_ptr< monty::ndarray< double,1 > > mval{};
      std::shared_ptr< monty::ndarray< int32_t,1 > > msubj{};
      std::shared_ptr< monty::ndarray< int32_t,1 > > msubi{};
      int32_t mdim1{};
      int32_t mdim0{};

      virtual void destroy();

      static ExprMulRight::t _new_ExprMulRight(int32_t _5931_mdim0,int32_t _5932_mdim1,std::shared_ptr< monty::ndarray< int32_t,1 > > _5933_msubi,std::shared_ptr< monty::ndarray< int32_t,1 > > _5934_msubj,std::shared_ptr< monty::ndarray< double,1 > > _5935_mval,monty::rc_ptr< ::mosek::fusion::Expression > _5936_expr,int32_t _5937_validated);
      void _initialize(int32_t _5931_mdim0,int32_t _5932_mdim1,std::shared_ptr< monty::ndarray< int32_t,1 > > _5933_msubi,std::shared_ptr< monty::ndarray< int32_t,1 > > _5934_msubj,std::shared_ptr< monty::ndarray< double,1 > > _5935_mval,monty::rc_ptr< ::mosek::fusion::Expression > _5936_expr,int32_t _5937_validated);
      static ExprMulRight::t _new_ExprMulRight(int32_t _5938_mdim0,int32_t _5939_mdim1,std::shared_ptr< monty::ndarray< int32_t,1 > > _5940_msubi,std::shared_ptr< monty::ndarray< int32_t,1 > > _5941_msubj,std::shared_ptr< monty::ndarray< double,1 > > _5942_mval,monty::rc_ptr< ::mosek::fusion::Expression > _5943_expr);
      void _initialize(int32_t _5938_mdim0,int32_t _5939_mdim1,std::shared_ptr< monty::ndarray< int32_t,1 > > _5940_msubi,std::shared_ptr< monty::ndarray< int32_t,1 > > _5941_msubj,std::shared_ptr< monty::ndarray< double,1 > > _5942_mval,monty::rc_ptr< ::mosek::fusion::Expression > _5943_expr);
      virtual /* override */ void eval(monty::rc_ptr< ::mosek::fusion::WorkStack > _5944_rs,monty::rc_ptr< ::mosek::fusion::WorkStack > _5945_ws,monty::rc_ptr< ::mosek::fusion::WorkStack > _5946_xs) ;
      static  std::shared_ptr< monty::ndarray< int32_t,1 > > computeshape(int32_t _6090_d0,std::shared_ptr< monty::ndarray< int32_t,1 > > _6091_ds);
      static  int32_t validate(int32_t _6092_mdim0,int32_t _6093_mdim1,std::shared_ptr< monty::ndarray< int32_t,1 > > _6094_msubi,std::shared_ptr< monty::ndarray< int32_t,1 > > _6095_msubj,std::shared_ptr< monty::ndarray< double,1 > > _6096_mval,monty::rc_ptr< ::mosek::fusion::Expression > _6097_expr);
      virtual /* override */ std::string toString() ;
    }; // struct ExprMulRight;

    struct p_ExprMulLeft : public ::mosek::fusion::p_BaseExpression
    {
      ExprMulLeft * _pubthis;
      static mosek::fusion::p_ExprMulLeft* _get_impl(mosek::fusion::ExprMulLeft * _inst){ return static_cast< mosek::fusion::p_ExprMulLeft* >(mosek::fusion::p_BaseExpression::_get_impl(_inst)); }
      static mosek::fusion::p_ExprMulLeft * _get_impl(mosek::fusion::ExprMulLeft::t _inst) { return _get_impl(_inst.get()); }
      p_ExprMulLeft(ExprMulLeft * _pubthis);
      virtual ~p_ExprMulLeft() { /* std::cout << "~p_ExprMulLeft" << std::endl;*/ };
      monty::rc_ptr< ::mosek::fusion::Expression > expr{};
      std::shared_ptr< monty::ndarray< double,1 > > mval{};
      std::shared_ptr< monty::ndarray< int32_t,1 > > msubj{};
      std::shared_ptr< monty::ndarray< int32_t,1 > > msubi{};
      int32_t mdim1{};
      int32_t mdim0{};

      virtual void destroy();

      static ExprMulLeft::t _new_ExprMulLeft(int32_t _6106_mdim0,int32_t _6107_mdim1,std::shared_ptr< monty::ndarray< int32_t,1 > > _6108_msubi,std::shared_ptr< monty::ndarray< int32_t,1 > > _6109_msubj,std::shared_ptr< monty::ndarray< double,1 > > _6110_mval,monty::rc_ptr< ::mosek::fusion::Expression > _6111_expr,int32_t _6112_validated);
      void _initialize(int32_t _6106_mdim0,int32_t _6107_mdim1,std::shared_ptr< monty::ndarray< int32_t,1 > > _6108_msubi,std::shared_ptr< monty::ndarray< int32_t,1 > > _6109_msubj,std::shared_ptr< monty::ndarray< double,1 > > _6110_mval,monty::rc_ptr< ::mosek::fusion::Expression > _6111_expr,int32_t _6112_validated);
      static ExprMulLeft::t _new_ExprMulLeft(int32_t _6113_mdim0,int32_t _6114_mdim1,std::shared_ptr< monty::ndarray< int32_t,1 > > _6115_msubi,std::shared_ptr< monty::ndarray< int32_t,1 > > _6116_msubj,std::shared_ptr< monty::ndarray< double,1 > > _6117_mval,monty::rc_ptr< ::mosek::fusion::Expression > _6118_expr);
      void _initialize(int32_t _6113_mdim0,int32_t _6114_mdim1,std::shared_ptr< monty::ndarray< int32_t,1 > > _6115_msubi,std::shared_ptr< monty::ndarray< int32_t,1 > > _6116_msubj,std::shared_ptr< monty::ndarray< double,1 > > _6117_mval,monty::rc_ptr< ::mosek::fusion::Expression > _6118_expr);
      virtual /* override */ void eval(monty::rc_ptr< ::mosek::fusion::WorkStack > _6119_rs,monty::rc_ptr< ::mosek::fusion::WorkStack > _6120_ws,monty::rc_ptr< ::mosek::fusion::WorkStack > _6121_xs) ;
      static  std::shared_ptr< monty::ndarray< int32_t,1 > > computeshape(int32_t _6223_d0,int32_t _6224_d1,std::shared_ptr< monty::ndarray< int32_t,1 > > _6225_ds);
      static  int32_t validate(int32_t _6226_mdim0,int32_t _6227_mdim1,std::shared_ptr< monty::ndarray< int32_t,1 > > _6228_msubi,std::shared_ptr< monty::ndarray< int32_t,1 > > _6229_msubj,std::shared_ptr< monty::ndarray< double,1 > > _6230_mval,monty::rc_ptr< ::mosek::fusion::Expression > _6231_expr);
      virtual /* override */ std::string toString() ;
    }; // struct ExprMulLeft;

    struct p_ExprMulVar : public ::mosek::fusion::p_BaseExpression
    {
      ExprMulVar * _pubthis;
      static mosek::fusion::p_ExprMulVar* _get_impl(mosek::fusion::ExprMulVar * _inst){ return static_cast< mosek::fusion::p_ExprMulVar* >(mosek::fusion::p_BaseExpression::_get_impl(_inst)); }
      static mosek::fusion::p_ExprMulVar * _get_impl(mosek::fusion::ExprMulVar::t _inst) { return _get_impl(_inst.get()); }
      p_ExprMulVar(ExprMulVar * _pubthis);
      virtual ~p_ExprMulVar() { /* std::cout << "~p_ExprMulVar" << std::endl;*/ };
      bool left{};
      monty::rc_ptr< ::mosek::fusion::Variable > x{};
      std::shared_ptr< monty::ndarray< double,1 > > mcof{};
      std::shared_ptr< monty::ndarray< int32_t,1 > > msubj{};
      std::shared_ptr< monty::ndarray< int32_t,1 > > msubi{};
      int32_t mdimj{};
      int32_t mdimi{};

      virtual void destroy();

      static ExprMulVar::t _new_ExprMulVar(bool _6239_left,int32_t _6240_mdimi,int32_t _6241_mdimj,std::shared_ptr< monty::ndarray< int32_t,1 > > _6242_msubi,std::shared_ptr< monty::ndarray< int32_t,1 > > _6243_msubj,std::shared_ptr< monty::ndarray< double,1 > > _6244_mcof,monty::rc_ptr< ::mosek::fusion::Variable > _6245_x);
      void _initialize(bool _6239_left,int32_t _6240_mdimi,int32_t _6241_mdimj,std::shared_ptr< monty::ndarray< int32_t,1 > > _6242_msubi,std::shared_ptr< monty::ndarray< int32_t,1 > > _6243_msubj,std::shared_ptr< monty::ndarray< double,1 > > _6244_mcof,monty::rc_ptr< ::mosek::fusion::Variable > _6245_x);
      static ExprMulVar::t _new_ExprMulVar(bool _6248_left,int32_t _6249_mdimi,int32_t _6250_mdimj,std::shared_ptr< monty::ndarray< int32_t,1 > > _6251_msubi,std::shared_ptr< monty::ndarray< int32_t,1 > > _6252_msubj,std::shared_ptr< monty::ndarray< double,1 > > _6253_mcof,monty::rc_ptr< ::mosek::fusion::Variable > _6254_x,int32_t _6255_unchecked_);
      void _initialize(bool _6248_left,int32_t _6249_mdimi,int32_t _6250_mdimj,std::shared_ptr< monty::ndarray< int32_t,1 > > _6251_msubi,std::shared_ptr< monty::ndarray< int32_t,1 > > _6252_msubj,std::shared_ptr< monty::ndarray< double,1 > > _6253_mcof,monty::rc_ptr< ::mosek::fusion::Variable > _6254_x,int32_t _6255_unchecked_);
      virtual /* override */ void eval(monty::rc_ptr< ::mosek::fusion::WorkStack > _6256_rs,monty::rc_ptr< ::mosek::fusion::WorkStack > _6257_ws,monty::rc_ptr< ::mosek::fusion::WorkStack > _6258_xs) ;
      virtual void eval_right(monty::rc_ptr< ::mosek::fusion::WorkStack > _6259_rs,monty::rc_ptr< ::mosek::fusion::WorkStack > _6260_ws,monty::rc_ptr< ::mosek::fusion::WorkStack > _6261_xs) ;
      virtual void eval_left(monty::rc_ptr< ::mosek::fusion::WorkStack > _6366_rs,monty::rc_ptr< ::mosek::fusion::WorkStack > _6367_ws,monty::rc_ptr< ::mosek::fusion::WorkStack > _6368_xs) ;
      virtual void validate(int32_t _6441_mdimi,int32_t _6442_mdimj,std::shared_ptr< monty::ndarray< int32_t,1 > > _6443_msubi,std::shared_ptr< monty::ndarray< int32_t,1 > > _6444_msubj,std::shared_ptr< monty::ndarray< double,1 > > _6445_mcof) ;
      static  std::shared_ptr< monty::ndarray< int32_t,1 > > resshape(int32_t _6449_mdimi,int32_t _6450_mdimj,std::shared_ptr< monty::ndarray< int32_t,1 > > _6451_xshape,bool _6452_left);
      virtual /* override */ std::string toString() ;
    }; // struct ExprMulVar;

    struct p_ExprMulScalarVar : public ::mosek::fusion::p_BaseExpression
    {
      ExprMulScalarVar * _pubthis;
      static mosek::fusion::p_ExprMulScalarVar* _get_impl(mosek::fusion::ExprMulScalarVar * _inst){ return static_cast< mosek::fusion::p_ExprMulScalarVar* >(mosek::fusion::p_BaseExpression::_get_impl(_inst)); }
      static mosek::fusion::p_ExprMulScalarVar * _get_impl(mosek::fusion::ExprMulScalarVar::t _inst) { return _get_impl(_inst.get()); }
      p_ExprMulScalarVar(ExprMulScalarVar * _pubthis);
      virtual ~p_ExprMulScalarVar() { /* std::cout << "~p_ExprMulScalarVar" << std::endl;*/ };
      monty::rc_ptr< ::mosek::fusion::Variable > x{};
      std::shared_ptr< monty::ndarray< double,1 > > mcof{};
      std::shared_ptr< monty::ndarray< int32_t,1 > > msubj{};
      std::shared_ptr< monty::ndarray< int32_t,1 > > msubi{};
      int32_t mdimj{};
      int32_t mdimi{};

      virtual void destroy();

      static ExprMulScalarVar::t _new_ExprMulScalarVar(int32_t _6453_mdimi,int32_t _6454_mdimj,std::shared_ptr< monty::ndarray< int32_t,1 > > _6455_msubi,std::shared_ptr< monty::ndarray< int32_t,1 > > _6456_msubj,std::shared_ptr< monty::ndarray< double,1 > > _6457_mcof,monty::rc_ptr< ::mosek::fusion::Variable > _6458_x);
      void _initialize(int32_t _6453_mdimi,int32_t _6454_mdimj,std::shared_ptr< monty::ndarray< int32_t,1 > > _6455_msubi,std::shared_ptr< monty::ndarray< int32_t,1 > > _6456_msubj,std::shared_ptr< monty::ndarray< double,1 > > _6457_mcof,monty::rc_ptr< ::mosek::fusion::Variable > _6458_x);
      static ExprMulScalarVar::t _new_ExprMulScalarVar(int32_t _6463_mdimi,int32_t _6464_mdimj,std::shared_ptr< monty::ndarray< int32_t,1 > > _6465_msubi,std::shared_ptr< monty::ndarray< int32_t,1 > > _6466_msubj,std::shared_ptr< monty::ndarray< double,1 > > _6467_mcof,monty::rc_ptr< ::mosek::fusion::Variable > _6468_x,int32_t _6469_unchecked_);
      void _initialize(int32_t _6463_mdimi,int32_t _6464_mdimj,std::shared_ptr< monty::ndarray< int32_t,1 > > _6465_msubi,std::shared_ptr< monty::ndarray< int32_t,1 > > _6466_msubj,std::shared_ptr< monty::ndarray< double,1 > > _6467_mcof,monty::rc_ptr< ::mosek::fusion::Variable > _6468_x,int32_t _6469_unchecked_);
      virtual /* override */ void eval(monty::rc_ptr< ::mosek::fusion::WorkStack > _6470_rs,monty::rc_ptr< ::mosek::fusion::WorkStack > _6471_ws,monty::rc_ptr< ::mosek::fusion::WorkStack > _6472_xs) ;
      virtual /* override */ std::string toString() ;
    }; // struct ExprMulScalarVar;

    struct p_ExprMulVarScalarConst : public ::mosek::fusion::p_BaseExpression
    {
      ExprMulVarScalarConst * _pubthis;
      static mosek::fusion::p_ExprMulVarScalarConst* _get_impl(mosek::fusion::ExprMulVarScalarConst * _inst){ return static_cast< mosek::fusion::p_ExprMulVarScalarConst* >(mosek::fusion::p_BaseExpression::_get_impl(_inst)); }
      static mosek::fusion::p_ExprMulVarScalarConst * _get_impl(mosek::fusion::ExprMulVarScalarConst::t _inst) { return _get_impl(_inst.get()); }
      p_ExprMulVarScalarConst(ExprMulVarScalarConst * _pubthis);
      virtual ~p_ExprMulVarScalarConst() { /* std::cout << "~p_ExprMulVarScalarConst" << std::endl;*/ };
      double c{};
      monty::rc_ptr< ::mosek::fusion::Variable > x{};

      virtual void destroy();

      static ExprMulVarScalarConst::t _new_ExprMulVarScalarConst(monty::rc_ptr< ::mosek::fusion::Variable > _6489_x,double _6490_c);
      void _initialize(monty::rc_ptr< ::mosek::fusion::Variable > _6489_x,double _6490_c);
      virtual /* override */ void eval(monty::rc_ptr< ::mosek::fusion::WorkStack > _6491_rs,monty::rc_ptr< ::mosek::fusion::WorkStack > _6492_ws,monty::rc_ptr< ::mosek::fusion::WorkStack > _6493_xs) ;
      virtual /* override */ std::string toString() ;
    }; // struct ExprMulVarScalarConst;

    struct p_ExprAdd : public ::mosek::fusion::p_BaseExpression
    {
      ExprAdd * _pubthis;
      static mosek::fusion::p_ExprAdd* _get_impl(mosek::fusion::ExprAdd * _inst){ return static_cast< mosek::fusion::p_ExprAdd* >(mosek::fusion::p_BaseExpression::_get_impl(_inst)); }
      static mosek::fusion::p_ExprAdd * _get_impl(mosek::fusion::ExprAdd::t _inst) { return _get_impl(_inst.get()); }
      p_ExprAdd(ExprAdd * _pubthis);
      virtual ~p_ExprAdd() { /* std::cout << "~p_ExprAdd" << std::endl;*/ };
      double m2{};
      double m1{};
      monty::rc_ptr< ::mosek::fusion::Expression > e2{};
      monty::rc_ptr< ::mosek::fusion::Expression > e1{};

      virtual void destroy();

      static ExprAdd::t _new_ExprAdd(monty::rc_ptr< ::mosek::fusion::Expression > _6510_e1,monty::rc_ptr< ::mosek::fusion::Expression > _6511_e2,double _6512_m1,double _6513_m2);
      void _initialize(monty::rc_ptr< ::mosek::fusion::Expression > _6510_e1,monty::rc_ptr< ::mosek::fusion::Expression > _6511_e2,double _6512_m1,double _6513_m2);
      virtual /* override */ void eval(monty::rc_ptr< ::mosek::fusion::WorkStack > _6515_rs,monty::rc_ptr< ::mosek::fusion::WorkStack > _6516_ws,monty::rc_ptr< ::mosek::fusion::WorkStack > _6517_xs) ;
      virtual /* override */ std::string toString() ;
    }; // struct ExprAdd;

    struct p_ExprWSum : public ::mosek::fusion::p_BaseExpression
    {
      ExprWSum * _pubthis;
      static mosek::fusion::p_ExprWSum* _get_impl(mosek::fusion::ExprWSum * _inst){ return static_cast< mosek::fusion::p_ExprWSum* >(mosek::fusion::p_BaseExpression::_get_impl(_inst)); }
      static mosek::fusion::p_ExprWSum * _get_impl(mosek::fusion::ExprWSum::t _inst) { return _get_impl(_inst.get()); }
      p_ExprWSum(ExprWSum * _pubthis);
      virtual ~p_ExprWSum() { /* std::cout << "~p_ExprWSum" << std::endl;*/ };
      std::shared_ptr< monty::ndarray< double,1 > > w{};
      std::shared_ptr< monty::ndarray< monty::rc_ptr< ::mosek::fusion::Expression >,1 > > es{};

      virtual void destroy();

      static ExprWSum::t _new_ExprWSum(std::shared_ptr< monty::ndarray< monty::rc_ptr< ::mosek::fusion::Expression >,1 > > _6651_es,std::shared_ptr< monty::ndarray< double,1 > > _6652_w);
      void _initialize(std::shared_ptr< monty::ndarray< monty::rc_ptr< ::mosek::fusion::Expression >,1 > > _6651_es,std::shared_ptr< monty::ndarray< double,1 > > _6652_w);
      virtual /* override */ void eval(monty::rc_ptr< ::mosek::fusion::WorkStack > _6659_rs,monty::rc_ptr< ::mosek::fusion::WorkStack > _6660_ws,monty::rc_ptr< ::mosek::fusion::WorkStack > _6661_xs) ;
      virtual /* override */ std::string toString() ;
    }; // struct ExprWSum;

    struct p_ExprSumReduce : public ::mosek::fusion::p_BaseExpression
    {
      ExprSumReduce * _pubthis;
      static mosek::fusion::p_ExprSumReduce* _get_impl(mosek::fusion::ExprSumReduce * _inst){ return static_cast< mosek::fusion::p_ExprSumReduce* >(mosek::fusion::p_BaseExpression::_get_impl(_inst)); }
      static mosek::fusion::p_ExprSumReduce * _get_impl(mosek::fusion::ExprSumReduce::t _inst) { return _get_impl(_inst.get()); }
      p_ExprSumReduce(ExprSumReduce * _pubthis);
      virtual ~p_ExprSumReduce() { /* std::cout << "~p_ExprSumReduce" << std::endl;*/ };
      int32_t dim{};
      monty::rc_ptr< ::mosek::fusion::Expression > expr{};

      virtual void destroy();

      static ExprSumReduce::t _new_ExprSumReduce(int32_t _6755_dim,monty::rc_ptr< ::mosek::fusion::Expression > _6756_expr);
      void _initialize(int32_t _6755_dim,monty::rc_ptr< ::mosek::fusion::Expression > _6756_expr);
      virtual /* override */ void eval(monty::rc_ptr< ::mosek::fusion::WorkStack > _6758_rs,monty::rc_ptr< ::mosek::fusion::WorkStack > _6759_ws,monty::rc_ptr< ::mosek::fusion::WorkStack > _6760_xs) ;
      static  std::shared_ptr< monty::ndarray< int32_t,1 > > computeShape(int32_t _6876_dim,std::shared_ptr< monty::ndarray< int32_t,1 > > _6877_shape);
      virtual /* override */ std::string toString() ;
    }; // struct ExprSumReduce;

    struct p_ExprScaleVecPSD : public ::mosek::fusion::p_BaseExpression
    {
      ExprScaleVecPSD * _pubthis;
      static mosek::fusion::p_ExprScaleVecPSD* _get_impl(mosek::fusion::ExprScaleVecPSD * _inst){ return static_cast< mosek::fusion::p_ExprScaleVecPSD* >(mosek::fusion::p_BaseExpression::_get_impl(_inst)); }
      static mosek::fusion::p_ExprScaleVecPSD * _get_impl(mosek::fusion::ExprScaleVecPSD::t _inst) { return _get_impl(_inst.get()); }
      p_ExprScaleVecPSD(ExprScaleVecPSD * _pubthis);
      virtual ~p_ExprScaleVecPSD() { /* std::cout << "~p_ExprScaleVecPSD" << std::endl;*/ };
      int32_t dim1{};
      int32_t dim0{};
      monty::rc_ptr< ::mosek::fusion::Expression > expr{};

      virtual void destroy();

      static ExprScaleVecPSD::t _new_ExprScaleVecPSD(int32_t _6881_dim0,int32_t _6882_dim1,monty::rc_ptr< ::mosek::fusion::BaseExpression > _6883_expr);
      void _initialize(int32_t _6881_dim0,int32_t _6882_dim1,monty::rc_ptr< ::mosek::fusion::BaseExpression > _6883_expr);
      virtual /* override */ void eval(monty::rc_ptr< ::mosek::fusion::WorkStack > _6884_rs,monty::rc_ptr< ::mosek::fusion::WorkStack > _6885_ws,monty::rc_ptr< ::mosek::fusion::WorkStack > _6886_xs) ;
    }; // struct ExprScaleVecPSD;

    struct p_ExprDenseTril : public ::mosek::fusion::p_BaseExpression
    {
      ExprDenseTril * _pubthis;
      static mosek::fusion::p_ExprDenseTril* _get_impl(mosek::fusion::ExprDenseTril * _inst){ return static_cast< mosek::fusion::p_ExprDenseTril* >(mosek::fusion::p_BaseExpression::_get_impl(_inst)); }
      static mosek::fusion::p_ExprDenseTril * _get_impl(mosek::fusion::ExprDenseTril::t _inst) { return _get_impl(_inst.get()); }
      p_ExprDenseTril(ExprDenseTril * _pubthis);
      virtual ~p_ExprDenseTril() { /* std::cout << "~p_ExprDenseTril" << std::endl;*/ };
      int32_t dim1{};
      int32_t dim0{};
      monty::rc_ptr< ::mosek::fusion::Expression > expr{};

      virtual void destroy();

      static ExprDenseTril::t _new_ExprDenseTril(int32_t _6959_dim0,int32_t _6960_dim1,monty::rc_ptr< ::mosek::fusion::Expression > _6961_expr,int32_t _6962_unchecked_);
      void _initialize(int32_t _6959_dim0,int32_t _6960_dim1,monty::rc_ptr< ::mosek::fusion::Expression > _6961_expr,int32_t _6962_unchecked_);
      static ExprDenseTril::t _new_ExprDenseTril(int32_t _6963_dim0_,int32_t _6964_dim1_,monty::rc_ptr< ::mosek::fusion::Expression > _6965_expr);
      void _initialize(int32_t _6963_dim0_,int32_t _6964_dim1_,monty::rc_ptr< ::mosek::fusion::Expression > _6965_expr);
      virtual /* override */ void eval(monty::rc_ptr< ::mosek::fusion::WorkStack > _6967_rs,monty::rc_ptr< ::mosek::fusion::WorkStack > _6968_ws,monty::rc_ptr< ::mosek::fusion::WorkStack > _6969_xs) ;
      virtual /* override */ std::string toString() ;
    }; // struct ExprDenseTril;

    struct p_ExprDense : public ::mosek::fusion::p_BaseExpression
    {
      ExprDense * _pubthis;
      static mosek::fusion::p_ExprDense* _get_impl(mosek::fusion::ExprDense * _inst){ return static_cast< mosek::fusion::p_ExprDense* >(mosek::fusion::p_BaseExpression::_get_impl(_inst)); }
      static mosek::fusion::p_ExprDense * _get_impl(mosek::fusion::ExprDense::t _inst) { return _get_impl(_inst.get()); }
      p_ExprDense(ExprDense * _pubthis);
      virtual ~p_ExprDense() { /* std::cout << "~p_ExprDense" << std::endl;*/ };
      monty::rc_ptr< ::mosek::fusion::Expression > expr{};

      virtual void destroy();

      static ExprDense::t _new_ExprDense(monty::rc_ptr< ::mosek::fusion::Expression > _7053_expr);
      void _initialize(monty::rc_ptr< ::mosek::fusion::Expression > _7053_expr);
      virtual /* override */ void eval(monty::rc_ptr< ::mosek::fusion::WorkStack > _7054_rs,monty::rc_ptr< ::mosek::fusion::WorkStack > _7055_ws,monty::rc_ptr< ::mosek::fusion::WorkStack > _7056_xs) ;
      virtual /* override */ std::string toString() ;
    }; // struct ExprDense;

    struct p_ExprSymmetrize : public ::mosek::fusion::p_BaseExpression
    {
      ExprSymmetrize * _pubthis;
      static mosek::fusion::p_ExprSymmetrize* _get_impl(mosek::fusion::ExprSymmetrize * _inst){ return static_cast< mosek::fusion::p_ExprSymmetrize* >(mosek::fusion::p_BaseExpression::_get_impl(_inst)); }
      static mosek::fusion::p_ExprSymmetrize * _get_impl(mosek::fusion::ExprSymmetrize::t _inst) { return _get_impl(_inst.get()); }
      p_ExprSymmetrize(ExprSymmetrize * _pubthis);
      virtual ~p_ExprSymmetrize() { /* std::cout << "~p_ExprSymmetrize" << std::endl;*/ };
      int32_t dim1{};
      int32_t dim0{};
      monty::rc_ptr< ::mosek::fusion::Expression > expr{};

      virtual void destroy();

      static ExprSymmetrize::t _new_ExprSymmetrize(int32_t _7097_dim0,int32_t _7098_dim1,monty::rc_ptr< ::mosek::fusion::Expression > _7099_expr,int32_t _7100_unchecked_);
      void _initialize(int32_t _7097_dim0,int32_t _7098_dim1,monty::rc_ptr< ::mosek::fusion::Expression > _7099_expr,int32_t _7100_unchecked_);
      static ExprSymmetrize::t _new_ExprSymmetrize(int32_t _7101_dim0_,int32_t _7102_dim1_,monty::rc_ptr< ::mosek::fusion::Expression > _7103_expr);
      void _initialize(int32_t _7101_dim0_,int32_t _7102_dim1_,monty::rc_ptr< ::mosek::fusion::Expression > _7103_expr);
      virtual /* override */ void eval(monty::rc_ptr< ::mosek::fusion::WorkStack > _7105_rs,monty::rc_ptr< ::mosek::fusion::WorkStack > _7106_ws,monty::rc_ptr< ::mosek::fusion::WorkStack > _7107_xs) ;
      virtual /* override */ std::string toString() ;
    }; // struct ExprSymmetrize;

    struct p_ExprCondense : public ::mosek::fusion::p_BaseExpression
    {
      ExprCondense * _pubthis;
      static mosek::fusion::p_ExprCondense* _get_impl(mosek::fusion::ExprCondense * _inst){ return static_cast< mosek::fusion::p_ExprCondense* >(mosek::fusion::p_BaseExpression::_get_impl(_inst)); }
      static mosek::fusion::p_ExprCondense * _get_impl(mosek::fusion::ExprCondense::t _inst) { return _get_impl(_inst.get()); }
      p_ExprCondense(ExprCondense * _pubthis);
      virtual ~p_ExprCondense() { /* std::cout << "~p_ExprCondense" << std::endl;*/ };
      monty::rc_ptr< ::mosek::fusion::Expression > expr{};

      virtual void destroy();

      static ExprCondense::t _new_ExprCondense(monty::rc_ptr< ::mosek::fusion::Expression > _7231_expr);
      void _initialize(monty::rc_ptr< ::mosek::fusion::Expression > _7231_expr);
      virtual /* override */ void eval(monty::rc_ptr< ::mosek::fusion::WorkStack > _7232_rs,monty::rc_ptr< ::mosek::fusion::WorkStack > _7233_ws,monty::rc_ptr< ::mosek::fusion::WorkStack > _7234_xs) ;
      virtual /* override */ std::string toString() ;
    }; // struct ExprCondense;

    struct p_ExprFromVar : public ::mosek::fusion::p_BaseExpression
    {
      ExprFromVar * _pubthis;
      static mosek::fusion::p_ExprFromVar* _get_impl(mosek::fusion::ExprFromVar * _inst){ return static_cast< mosek::fusion::p_ExprFromVar* >(mosek::fusion::p_BaseExpression::_get_impl(_inst)); }
      static mosek::fusion::p_ExprFromVar * _get_impl(mosek::fusion::ExprFromVar::t _inst) { return _get_impl(_inst.get()); }
      p_ExprFromVar(ExprFromVar * _pubthis);
      virtual ~p_ExprFromVar() { /* std::cout << "~p_ExprFromVar" << std::endl;*/ };
      monty::rc_ptr< ::mosek::fusion::Variable > x{};

      virtual void destroy();

      static ExprFromVar::t _new_ExprFromVar(monty::rc_ptr< ::mosek::fusion::Variable > _7238_x);
      void _initialize(monty::rc_ptr< ::mosek::fusion::Variable > _7238_x);
      virtual /* override */ void eval(monty::rc_ptr< ::mosek::fusion::WorkStack > _7239_rs,monty::rc_ptr< ::mosek::fusion::WorkStack > _7240_ws,monty::rc_ptr< ::mosek::fusion::WorkStack > _7241_xs) ;
      virtual /* override */ std::string toString() ;
    }; // struct ExprFromVar;

    struct p_ExprReshape : public ::mosek::fusion::p_BaseExpression
    {
      ExprReshape * _pubthis;
      static mosek::fusion::p_ExprReshape* _get_impl(mosek::fusion::ExprReshape * _inst){ return static_cast< mosek::fusion::p_ExprReshape* >(mosek::fusion::p_BaseExpression::_get_impl(_inst)); }
      static mosek::fusion::p_ExprReshape * _get_impl(mosek::fusion::ExprReshape::t _inst) { return _get_impl(_inst.get()); }
      p_ExprReshape(ExprReshape * _pubthis);
      virtual ~p_ExprReshape() { /* std::cout << "~p_ExprReshape" << std::endl;*/ };
      monty::rc_ptr< ::mosek::fusion::Expression > e{};

      virtual void destroy();

      static ExprReshape::t _new_ExprReshape(std::shared_ptr< monty::ndarray< int32_t,1 > > _7258_shape,monty::rc_ptr< ::mosek::fusion::Expression > _7259_e);
      void _initialize(std::shared_ptr< monty::ndarray< int32_t,1 > > _7258_shape,monty::rc_ptr< ::mosek::fusion::Expression > _7259_e);
      virtual /* override */ void eval(monty::rc_ptr< ::mosek::fusion::WorkStack > _7261_rs,monty::rc_ptr< ::mosek::fusion::WorkStack > _7262_ws,monty::rc_ptr< ::mosek::fusion::WorkStack > _7263_xs) ;
      virtual /* override */ std::string toString() ;
    }; // struct ExprReshape;

    struct p_Expr : public ::mosek::fusion::p_BaseExpression
    {
      Expr * _pubthis;
      static mosek::fusion::p_Expr* _get_impl(mosek::fusion::Expr * _inst){ return static_cast< mosek::fusion::p_Expr* >(mosek::fusion::p_BaseExpression::_get_impl(_inst)); }
      static mosek::fusion::p_Expr * _get_impl(mosek::fusion::Expr::t _inst) { return _get_impl(_inst.get()); }
      p_Expr(Expr * _pubthis);
      virtual ~p_Expr() { /* std::cout << "~p_Expr" << std::endl;*/ };
      std::shared_ptr< monty::ndarray< int64_t,1 > > inst{};
      std::shared_ptr< monty::ndarray< double,1 > > cof_v{};
      std::shared_ptr< monty::ndarray< int64_t,1 > > subj{};
      std::shared_ptr< monty::ndarray< int64_t,1 > > ptrb{};
      std::shared_ptr< monty::ndarray< double,1 > > bfix{};
      std::shared_ptr< monty::ndarray< int32_t,1 > > shape{};

      virtual void destroy();

      static Expr::t _new_Expr(std::shared_ptr< monty::ndarray< int64_t,1 > > _7395_ptrb,std::shared_ptr< monty::ndarray< int64_t,1 > > _7396_subj,std::shared_ptr< monty::ndarray< double,1 > > _7397_cof,std::shared_ptr< monty::ndarray< double,1 > > _7398_bfix,std::shared_ptr< monty::ndarray< int32_t,1 > > _7399_shape,std::shared_ptr< monty::ndarray< int64_t,1 > > _7400_inst);
      void _initialize(std::shared_ptr< monty::ndarray< int64_t,1 > > _7395_ptrb,std::shared_ptr< monty::ndarray< int64_t,1 > > _7396_subj,std::shared_ptr< monty::ndarray< double,1 > > _7397_cof,std::shared_ptr< monty::ndarray< double,1 > > _7398_bfix,std::shared_ptr< monty::ndarray< int32_t,1 > > _7399_shape,std::shared_ptr< monty::ndarray< int64_t,1 > > _7400_inst);
      static Expr::t _new_Expr(std::shared_ptr< monty::ndarray< int64_t,1 > > _7411_ptrb,std::shared_ptr< monty::ndarray< int64_t,1 > > _7412_subj,std::shared_ptr< monty::ndarray< double,1 > > _7413_cof,std::shared_ptr< monty::ndarray< double,1 > > _7414_bfix,std::shared_ptr< monty::ndarray< int32_t,1 > > _7415_shp,std::shared_ptr< monty::ndarray< int64_t,1 > > _7416_inst,int32_t _7417_unchecked_);
      void _initialize(std::shared_ptr< monty::ndarray< int64_t,1 > > _7411_ptrb,std::shared_ptr< monty::ndarray< int64_t,1 > > _7412_subj,std::shared_ptr< monty::ndarray< double,1 > > _7413_cof,std::shared_ptr< monty::ndarray< double,1 > > _7414_bfix,std::shared_ptr< monty::ndarray< int32_t,1 > > _7415_shp,std::shared_ptr< monty::ndarray< int64_t,1 > > _7416_inst,int32_t _7417_unchecked_);
      static Expr::t _new_Expr(monty::rc_ptr< ::mosek::fusion::Expression > _7418_e);
      void _initialize(monty::rc_ptr< ::mosek::fusion::Expression > _7418_e);
      virtual int64_t prod(std::shared_ptr< monty::ndarray< int32_t,1 > > _7443_vals) ;
      static  std::shared_ptr< monty::ndarray< monty::rc_ptr< ::mosek::fusion::Variable >,1 > > varstack(std::shared_ptr< monty::ndarray< std::shared_ptr< monty::ndarray< monty::rc_ptr< ::mosek::fusion::Variable >,1 > >,1 > > _7446_vs);
      static  std::shared_ptr< monty::ndarray< monty::rc_ptr< ::mosek::fusion::Variable >,1 > > varstack(std::shared_ptr< monty::ndarray< monty::rc_ptr< ::mosek::fusion::Variable >,1 > > _7449_v1,std::shared_ptr< monty::ndarray< monty::rc_ptr< ::mosek::fusion::Variable >,1 > > _7450_v2);
      static  monty::rc_ptr< ::mosek::fusion::Expression > condense(monty::rc_ptr< ::mosek::fusion::Expression > _7454_e);
      static  monty::rc_ptr< ::mosek::fusion::Expression > flatten(monty::rc_ptr< ::mosek::fusion::Expression > _7455_e);
      static  monty::rc_ptr< ::mosek::fusion::Expression > reshape(monty::rc_ptr< ::mosek::fusion::Expression > _7457_e,int32_t _7458_dimi,int32_t _7459_dimj);
      static  monty::rc_ptr< ::mosek::fusion::Expression > reshape(monty::rc_ptr< ::mosek::fusion::Expression > _7460_e,int32_t _7461_size);
      static  monty::rc_ptr< ::mosek::fusion::Expression > reshape(monty::rc_ptr< ::mosek::fusion::Expression > _7462_e,std::shared_ptr< monty::ndarray< int32_t,1 > > _7463_newshape);
      static  monty::rc_ptr< ::mosek::fusion::Expression > zeros(std::shared_ptr< monty::ndarray< int32_t,1 > > _7464_shp);
      static  monty::rc_ptr< ::mosek::fusion::Expression > zeros(int32_t _7465_size);
      static  monty::rc_ptr< ::mosek::fusion::Expression > ones();
      static  monty::rc_ptr< ::mosek::fusion::Expression > ones(std::shared_ptr< monty::ndarray< int32_t,1 > > _7466_shp,std::shared_ptr< monty::ndarray< int32_t,2 > > _7467_sparsity);
      static  monty::rc_ptr< ::mosek::fusion::Expression > ones(std::shared_ptr< monty::ndarray< int32_t,1 > > _7468_shp);
      static  monty::rc_ptr< ::mosek::fusion::Expression > ones(int32_t _7469_size);
      static  monty::rc_ptr< ::mosek::fusion::Expression > constTerm(monty::rc_ptr< ::mosek::fusion::NDSparseArray > _7470_nda);
      static  monty::rc_ptr< ::mosek::fusion::Expression > constTerm(monty::rc_ptr< ::mosek::fusion::Matrix > _7471_m);
      static  monty::rc_ptr< ::mosek::fusion::Expression > constTerm(double _7480_val);
      static  monty::rc_ptr< ::mosek::fusion::Expression > constTerm(std::shared_ptr< monty::ndarray< int32_t,1 > > _7481_shp,std::shared_ptr< monty::ndarray< int32_t,2 > > _7482_sparsity,double _7483_val);
      static  monty::rc_ptr< ::mosek::fusion::Expression > constTerm(std::shared_ptr< monty::ndarray< int32_t,1 > > _7491_shp,std::shared_ptr< monty::ndarray< int32_t,2 > > _7492_sparsity,std::shared_ptr< monty::ndarray< double,1 > > _7493_vals1);
      static  monty::rc_ptr< ::mosek::fusion::Expression > constTerm(std::shared_ptr< monty::ndarray< int32_t,1 > > _7501_shp,double _7502_val);
      static  monty::rc_ptr< ::mosek::fusion::Expression > constTerm(int32_t _7503_size,double _7504_val);
      static  monty::rc_ptr< ::mosek::fusion::Expression > constTerm(std::shared_ptr< monty::ndarray< double,2 > > _7506_vals2);
      static  monty::rc_ptr< ::mosek::fusion::Expression > constTerm(std::shared_ptr< monty::ndarray< double,1 > > _7509_vals1);
      static  monty::rc_ptr< ::mosek::fusion::Expression > sum(monty::rc_ptr< ::mosek::fusion::Expression > _7510_expr,int32_t _7511_dim);
      static  monty::rc_ptr< ::mosek::fusion::Expression > sum(monty::rc_ptr< ::mosek::fusion::Expression > _7512_expr);
      static  monty::rc_ptr< ::mosek::fusion::Expression > neg(monty::rc_ptr< ::mosek::fusion::Expression > _7513_e);
      static  monty::rc_ptr< ::mosek::fusion::Expression > mulDiag(bool _7514_left,monty::rc_ptr< ::mosek::fusion::Matrix > _7515_mx,monty::rc_ptr< ::mosek::fusion::Expression > _7516_expr);
      static  monty::rc_ptr< ::mosek::fusion::Expression > mulDiag(monty::rc_ptr< ::mosek::fusion::Variable > _7523_v,monty::rc_ptr< ::mosek::fusion::Parameter > _7524_p);
      static  monty::rc_ptr< ::mosek::fusion::Expression > mulDiag(monty::rc_ptr< ::mosek::fusion::Parameter > _7525_p,monty::rc_ptr< ::mosek::fusion::Variable > _7526_v);
      static  monty::rc_ptr< ::mosek::fusion::Expression > mulDiag(monty::rc_ptr< ::mosek::fusion::Expression > _7527_expr,monty::rc_ptr< ::mosek::fusion::Parameter > _7528_p);
      static  monty::rc_ptr< ::mosek::fusion::Expression > mulDiag(monty::rc_ptr< ::mosek::fusion::Parameter > _7529_p,monty::rc_ptr< ::mosek::fusion::Expression > _7530_expr);
      static  monty::rc_ptr< ::mosek::fusion::Expression > mulDiag(monty::rc_ptr< ::mosek::fusion::Variable > _7531_v,monty::rc_ptr< ::mosek::fusion::Matrix > _7532_mx);
      static  monty::rc_ptr< ::mosek::fusion::Expression > mulDiag(monty::rc_ptr< ::mosek::fusion::Matrix > _7533_mx,monty::rc_ptr< ::mosek::fusion::Variable > _7534_v);
      static  monty::rc_ptr< ::mosek::fusion::Expression > mulDiag(monty::rc_ptr< ::mosek::fusion::Expression > _7535_expr,monty::rc_ptr< ::mosek::fusion::Matrix > _7536_mx);
      static  monty::rc_ptr< ::mosek::fusion::Expression > mulDiag(monty::rc_ptr< ::mosek::fusion::Matrix > _7537_mx,monty::rc_ptr< ::mosek::fusion::Expression > _7538_expr);
      static  monty::rc_ptr< ::mosek::fusion::Expression > mulDiag(monty::rc_ptr< ::mosek::fusion::Variable > _7539_v,std::shared_ptr< monty::ndarray< double,2 > > _7540_a);
      static  monty::rc_ptr< ::mosek::fusion::Expression > mulDiag(monty::rc_ptr< ::mosek::fusion::Expression > _7547_expr,std::shared_ptr< monty::ndarray< double,2 > > _7548_a);
      static  monty::rc_ptr< ::mosek::fusion::Expression > mulDiag(std::shared_ptr< monty::ndarray< double,2 > > _7555_a,monty::rc_ptr< ::mosek::fusion::Variable > _7556_v);
      static  monty::rc_ptr< ::mosek::fusion::Expression > mulDiag(std::shared_ptr< monty::ndarray< double,2 > > _7563_a,monty::rc_ptr< ::mosek::fusion::Expression > _7564_expr);
      static  monty::rc_ptr< ::mosek::fusion::Expression > mulElm_(monty::rc_ptr< ::mosek::fusion::Matrix > _7571_m,monty::rc_ptr< ::mosek::fusion::Expression > _7572_e);
      static  monty::rc_ptr< ::mosek::fusion::Expression > mulElm_(std::shared_ptr< monty::ndarray< double,1 > > _7581_a,monty::rc_ptr< ::mosek::fusion::Expression > _7582_expr);
      static  monty::rc_ptr< ::mosek::fusion::Expression > mulElm_(monty::rc_ptr< ::mosek::fusion::NDSparseArray > _7584_spm,monty::rc_ptr< ::mosek::fusion::Expression > _7585_expr);
      static  monty::rc_ptr< ::mosek::fusion::Expression > mul(monty::rc_ptr< ::mosek::fusion::Expression > _7588_expr,double _7589_c);
      static  monty::rc_ptr< ::mosek::fusion::Expression > mul(double _7590_c,monty::rc_ptr< ::mosek::fusion::Expression > _7591_expr);
      static  monty::rc_ptr< ::mosek::fusion::Expression > mul(monty::rc_ptr< ::mosek::fusion::Expression > _7592_expr,std::shared_ptr< monty::ndarray< double,1 > > _7593_a);
      static  monty::rc_ptr< ::mosek::fusion::Expression > mul(std::shared_ptr< monty::ndarray< double,1 > > _7594_a,monty::rc_ptr< ::mosek::fusion::Expression > _7595_expr);
      static  monty::rc_ptr< ::mosek::fusion::Expression > mul(monty::rc_ptr< ::mosek::fusion::Expression > _7596_expr,std::shared_ptr< monty::ndarray< double,2 > > _7597_a);
      static  monty::rc_ptr< ::mosek::fusion::Expression > mul(std::shared_ptr< monty::ndarray< double,2 > > _7598_a,monty::rc_ptr< ::mosek::fusion::Expression > _7599_expr);
      static  monty::rc_ptr< ::mosek::fusion::Expression > mul(monty::rc_ptr< ::mosek::fusion::Expression > _7600_expr,monty::rc_ptr< ::mosek::fusion::Matrix > _7601_mx);
      static  monty::rc_ptr< ::mosek::fusion::Expression > mul(monty::rc_ptr< ::mosek::fusion::Matrix > _7602_mx,monty::rc_ptr< ::mosek::fusion::Expression > _7603_expr);
      static  monty::rc_ptr< ::mosek::fusion::Expression > mul(bool _7604_left,std::shared_ptr< monty::ndarray< double,1 > > _7605_mx,monty::rc_ptr< ::mosek::fusion::Expression > _7606_e);
      static  monty::rc_ptr< ::mosek::fusion::Expression > mul(bool _7621_left,std::shared_ptr< monty::ndarray< double,2 > > _7622_mx,monty::rc_ptr< ::mosek::fusion::Expression > _7623_e);
      static  monty::rc_ptr< ::mosek::fusion::Expression > mul(bool _7635_left,monty::rc_ptr< ::mosek::fusion::Matrix > _7636_mx,monty::rc_ptr< ::mosek::fusion::Expression > _7637_e);
      static  monty::rc_ptr< ::mosek::fusion::Expression > mul(monty::rc_ptr< ::mosek::fusion::Variable > _7646_v,std::shared_ptr< monty::ndarray< double,2 > > _7647_mx);
      static  monty::rc_ptr< ::mosek::fusion::Expression > mul(std::shared_ptr< monty::ndarray< double,2 > > _7657_mx,monty::rc_ptr< ::mosek::fusion::Variable > _7658_v);
      static  monty::rc_ptr< ::mosek::fusion::Expression > mul(monty::rc_ptr< ::mosek::fusion::Variable > _7668_v,monty::rc_ptr< ::mosek::fusion::Matrix > _7669_mx);
      static  monty::rc_ptr< ::mosek::fusion::Expression > mul(monty::rc_ptr< ::mosek::fusion::Matrix > _7675_mx,monty::rc_ptr< ::mosek::fusion::Variable > _7676_v);
      static  monty::rc_ptr< ::mosek::fusion::Expression > mul(bool _7682_left,int32_t _7683_mdimi,int32_t _7684_mdimj,std::shared_ptr< monty::ndarray< int32_t,1 > > _7685_msubi,std::shared_ptr< monty::ndarray< int32_t,1 > > _7686_msubj,std::shared_ptr< monty::ndarray< double,1 > > _7687_mcof,monty::rc_ptr< ::mosek::fusion::Variable > _7688_v);
      static  monty::rc_ptr< ::mosek::fusion::Expression > mul(monty::rc_ptr< ::mosek::fusion::Expression > _7690_expr,monty::rc_ptr< ::mosek::fusion::Parameter > _7691_p);
      static  monty::rc_ptr< ::mosek::fusion::Expression > mul(monty::rc_ptr< ::mosek::fusion::Parameter > _7692_p,monty::rc_ptr< ::mosek::fusion::Expression > _7693_expr);
      static  monty::rc_ptr< ::mosek::fusion::Expression > dot(monty::rc_ptr< ::mosek::fusion::Expression > _7694_e,monty::rc_ptr< ::mosek::fusion::Matrix > _7695_m);
      static  monty::rc_ptr< ::mosek::fusion::Expression > dot(monty::rc_ptr< ::mosek::fusion::Expression > _7703_e,std::shared_ptr< monty::ndarray< double,2 > > _7704_c2);
      static  monty::rc_ptr< ::mosek::fusion::Expression > dot(monty::rc_ptr< ::mosek::fusion::Expression > _7708_e,monty::rc_ptr< ::mosek::fusion::NDSparseArray > _7709_nda);
      static  monty::rc_ptr< ::mosek::fusion::Expression > dot(monty::rc_ptr< ::mosek::fusion::Expression > _7710_e,std::shared_ptr< monty::ndarray< double,1 > > _7711_c1);
      static  monty::rc_ptr< ::mosek::fusion::Expression > dot(monty::rc_ptr< ::mosek::fusion::Matrix > _7716_m,monty::rc_ptr< ::mosek::fusion::Expression > _7717_e);
      static  monty::rc_ptr< ::mosek::fusion::Expression > dot(monty::rc_ptr< ::mosek::fusion::NDSparseArray > _7718_nda,monty::rc_ptr< ::mosek::fusion::Expression > _7719_e);
      static  monty::rc_ptr< ::mosek::fusion::Expression > dot(std::shared_ptr< monty::ndarray< double,2 > > _7720_c2,monty::rc_ptr< ::mosek::fusion::Expression > _7721_e);
      static  monty::rc_ptr< ::mosek::fusion::Expression > dot(std::shared_ptr< monty::ndarray< double,1 > > _7722_c1,monty::rc_ptr< ::mosek::fusion::Expression > _7723_e);
      static  monty::rc_ptr< ::mosek::fusion::Expression > dot(monty::rc_ptr< ::mosek::fusion::Expression > _7724_e,monty::rc_ptr< ::mosek::fusion::Parameter > _7725_p);
      static  monty::rc_ptr< ::mosek::fusion::Expression > dot(monty::rc_ptr< ::mosek::fusion::Parameter > _7726_p,monty::rc_ptr< ::mosek::fusion::Expression > _7727_e);
      static  monty::rc_ptr< ::mosek::fusion::Expression > outer(monty::rc_ptr< ::mosek::fusion::Parameter > _7728_p,monty::rc_ptr< ::mosek::fusion::Expression > _7729_e);
      static  monty::rc_ptr< ::mosek::fusion::Expression > outer(monty::rc_ptr< ::mosek::fusion::Expression > _7732_e,monty::rc_ptr< ::mosek::fusion::Parameter > _7733_p);
      static  monty::rc_ptr< ::mosek::fusion::Expression > outer(monty::rc_ptr< ::mosek::fusion::Matrix > _7736_m,monty::rc_ptr< ::mosek::fusion::Expression > _7737_e);
      static  monty::rc_ptr< ::mosek::fusion::Expression > outer(monty::rc_ptr< ::mosek::fusion::Expression > _7739_e,monty::rc_ptr< ::mosek::fusion::Matrix > _7740_m);
      static  monty::rc_ptr< ::mosek::fusion::Expression > outer(std::shared_ptr< monty::ndarray< double,1 > > _7742_a,monty::rc_ptr< ::mosek::fusion::Expression > _7743_e);
      static  monty::rc_ptr< ::mosek::fusion::Expression > outer(monty::rc_ptr< ::mosek::fusion::Expression > _7745_e,std::shared_ptr< monty::ndarray< double,1 > > _7746_a);
      static  monty::rc_ptr< ::mosek::fusion::Expression > outer_(int32_t _7748_edim,std::shared_ptr< monty::ndarray< int64_t,1 > > _7749_eptrb,std::shared_ptr< monty::ndarray< int64_t,1 > > _7750_esubj,std::shared_ptr< monty::ndarray< double,1 > > _7751_ecof,std::shared_ptr< monty::ndarray< double,1 > > _7752_ebfix,std::shared_ptr< monty::ndarray< int64_t,1 > > _7753_einst,std::shared_ptr< monty::ndarray< double,1 > > _7754_a,std::shared_ptr< monty::ndarray< int32_t,1 > > _7755_sub,int32_t _7756_dim,bool _7757_transpose);
      static  monty::rc_ptr< ::mosek::fusion::Expression > outer_(monty::rc_ptr< ::mosek::fusion::Variable > _7787_v,int32_t _7788_vdim,std::shared_ptr< monty::ndarray< double,1 > > _7789_a,std::shared_ptr< monty::ndarray< int32_t,1 > > _7790_sub,int32_t _7791_dim,bool _7792_transpose);
      static  monty::rc_ptr< ::mosek::fusion::Expression > stack(std::shared_ptr< monty::ndarray< std::shared_ptr< monty::ndarray< monty::rc_ptr< ::mosek::fusion::Expression >,1 > >,1 > > _7809_exprs);
      static  monty::rc_ptr< ::mosek::fusion::Expression > vstack(double _7815_a1,double _7816_a2,double _7817_a3);
      static  monty::rc_ptr< ::mosek::fusion::Expression > vstack(double _7818_a1,double _7819_a2,monty::rc_ptr< ::mosek::fusion::Expression > _7820_e3);
      static  monty::rc_ptr< ::mosek::fusion::Expression > vstack(double _7821_a1,monty::rc_ptr< ::mosek::fusion::Expression > _7822_e2,double _7823_a3);
      static  monty::rc_ptr< ::mosek::fusion::Expression > vstack(double _7824_a1,monty::rc_ptr< ::mosek::fusion::Expression > _7825_e2,monty::rc_ptr< ::mosek::fusion::Expression > _7826_e3);
      static  monty::rc_ptr< ::mosek::fusion::Expression > vstack(monty::rc_ptr< ::mosek::fusion::Expression > _7827_e1,double _7828_a2,double _7829_a3);
      static  monty::rc_ptr< ::mosek::fusion::Expression > vstack(monty::rc_ptr< ::mosek::fusion::Expression > _7830_e1,double _7831_a2,monty::rc_ptr< ::mosek::fusion::Expression > _7832_e3);
      static  monty::rc_ptr< ::mosek::fusion::Expression > vstack(monty::rc_ptr< ::mosek::fusion::Expression > _7833_e1,monty::rc_ptr< ::mosek::fusion::Expression > _7834_e2,double _7835_a3);
      static  monty::rc_ptr< ::mosek::fusion::Expression > vstack(monty::rc_ptr< ::mosek::fusion::Expression > _7836_e1,monty::rc_ptr< ::mosek::fusion::Expression > _7837_e2,monty::rc_ptr< ::mosek::fusion::Expression > _7838_e3);
      static  monty::rc_ptr< ::mosek::fusion::Expression > vstack(double _7839_a1,monty::rc_ptr< ::mosek::fusion::Expression > _7840_e2);
      static  monty::rc_ptr< ::mosek::fusion::Expression > vstack(monty::rc_ptr< ::mosek::fusion::Expression > _7841_e1,double _7842_a2);
      static  monty::rc_ptr< ::mosek::fusion::Expression > vstack(monty::rc_ptr< ::mosek::fusion::Expression > _7843_e1,monty::rc_ptr< ::mosek::fusion::Expression > _7844_e2);
      static  monty::rc_ptr< ::mosek::fusion::Expression > vstack(std::shared_ptr< monty::ndarray< monty::rc_ptr< ::mosek::fusion::Expression >,1 > > _7845_exprs);
      static  monty::rc_ptr< ::mosek::fusion::Expression > hstack(monty::rc_ptr< ::mosek::fusion::Expression > _7847_e1,monty::rc_ptr< ::mosek::fusion::Expression > _7848_e2,monty::rc_ptr< ::mosek::fusion::Expression > _7849_e3);
      static  monty::rc_ptr< ::mosek::fusion::Expression > hstack(monty::rc_ptr< ::mosek::fusion::Expression > _7850_e1,monty::rc_ptr< ::mosek::fusion::Expression > _7851_e2,double _7852_a3);
      static  monty::rc_ptr< ::mosek::fusion::Expression > hstack(monty::rc_ptr< ::mosek::fusion::Expression > _7853_e1,double _7854_a2,monty::rc_ptr< ::mosek::fusion::Expression > _7855_e3);
      static  monty::rc_ptr< ::mosek::fusion::Expression > hstack(monty::rc_ptr< ::mosek::fusion::Expression > _7856_e1,double _7857_a2,double _7858_a3);
      static  monty::rc_ptr< ::mosek::fusion::Expression > hstack(double _7859_a1,monty::rc_ptr< ::mosek::fusion::Expression > _7860_e2,monty::rc_ptr< ::mosek::fusion::Expression > _7861_e3);
      static  monty::rc_ptr< ::mosek::fusion::Expression > hstack(double _7862_a1,monty::rc_ptr< ::mosek::fusion::Expression > _7863_e2,double _7864_a3);
      static  monty::rc_ptr< ::mosek::fusion::Expression > hstack(double _7865_a1,double _7866_a2,monty::rc_ptr< ::mosek::fusion::Expression > _7867_e3);
      static  monty::rc_ptr< ::mosek::fusion::Expression > hstack(double _7868_a1,monty::rc_ptr< ::mosek::fusion::Expression > _7869_e2);
      static  monty::rc_ptr< ::mosek::fusion::Expression > hstack(monty::rc_ptr< ::mosek::fusion::Expression > _7870_e1,double _7871_a2);
      static  monty::rc_ptr< ::mosek::fusion::Expression > hstack(monty::rc_ptr< ::mosek::fusion::Expression > _7872_e1,monty::rc_ptr< ::mosek::fusion::Expression > _7873_e2);
      static  monty::rc_ptr< ::mosek::fusion::Expression > hstack(std::shared_ptr< monty::ndarray< monty::rc_ptr< ::mosek::fusion::Expression >,1 > > _7874_exprs);
      static  monty::rc_ptr< ::mosek::fusion::Expression > stack(int32_t _7876_dim,monty::rc_ptr< ::mosek::fusion::Expression > _7877_e1,monty::rc_ptr< ::mosek::fusion::Expression > _7878_e2,monty::rc_ptr< ::mosek::fusion::Expression > _7879_e3);
      static  monty::rc_ptr< ::mosek::fusion::Expression > stack(int32_t _7880_dim,monty::rc_ptr< ::mosek::fusion::Expression > _7881_e1,monty::rc_ptr< ::mosek::fusion::Expression > _7882_e2,double _7883_a3);
      static  monty::rc_ptr< ::mosek::fusion::Expression > stack(int32_t _7884_dim,monty::rc_ptr< ::mosek::fusion::Expression > _7885_e1,double _7886_a2,monty::rc_ptr< ::mosek::fusion::Expression > _7887_e3);
      static  monty::rc_ptr< ::mosek::fusion::Expression > stack(int32_t _7888_dim,monty::rc_ptr< ::mosek::fusion::Expression > _7889_e1,double _7890_a2,double _7891_a3);
      static  monty::rc_ptr< ::mosek::fusion::Expression > stack(int32_t _7892_dim,double _7893_a1,monty::rc_ptr< ::mosek::fusion::Expression > _7894_e2,monty::rc_ptr< ::mosek::fusion::Expression > _7895_e3);
      static  monty::rc_ptr< ::mosek::fusion::Expression > stack(int32_t _7896_dim,double _7897_a1,monty::rc_ptr< ::mosek::fusion::Expression > _7898_e2,double _7899_a3);
      static  monty::rc_ptr< ::mosek::fusion::Expression > stack(int32_t _7900_dim,double _7901_a1,double _7902_a2,monty::rc_ptr< ::mosek::fusion::Expression > _7903_e1);
      static  monty::rc_ptr< ::mosek::fusion::Expression > stack(int32_t _7904_dim,double _7905_a1,monty::rc_ptr< ::mosek::fusion::Expression > _7906_e2);
      static  monty::rc_ptr< ::mosek::fusion::Expression > stack(int32_t _7907_dim,monty::rc_ptr< ::mosek::fusion::Expression > _7908_e1,double _7909_a2);
      static  monty::rc_ptr< ::mosek::fusion::Expression > stack(int32_t _7910_dim,monty::rc_ptr< ::mosek::fusion::Expression > _7911_e1,monty::rc_ptr< ::mosek::fusion::Expression > _7912_e2);
      static  monty::rc_ptr< ::mosek::fusion::Expression > stack(int32_t _7913_dim,std::shared_ptr< monty::ndarray< monty::rc_ptr< ::mosek::fusion::Expression >,1 > > _7914_exprs);
      static  monty::rc_ptr< ::mosek::fusion::Expression > stack_(std::shared_ptr< monty::ndarray< monty::rc_ptr< ::mosek::fusion::Expression >,1 > > _7915_exprs,int32_t _7916_dim);
      static  std::shared_ptr< monty::ndarray< monty::rc_ptr< ::mosek::fusion::Expression >,1 > > promote(std::shared_ptr< monty::ndarray< monty::rc_ptr< ::mosek::fusion::Expression >,1 > > _7917_exprs,int32_t _7918_dim);
      static  monty::rc_ptr< ::mosek::fusion::Expression > repeat(monty::rc_ptr< ::mosek::fusion::Variable > _7931_x,int32_t _7932_n,int32_t _7933_d);
      static  monty::rc_ptr< ::mosek::fusion::Expression > repeat(monty::rc_ptr< ::mosek::fusion::Expression > _7934_e,int32_t _7935_n,int32_t _7936_d);
      static  monty::rc_ptr< ::mosek::fusion::Expression > add(std::shared_ptr< monty::ndarray< monty::rc_ptr< ::mosek::fusion::Expression >,1 > > _7939_exps);
      static  monty::rc_ptr< ::mosek::fusion::Expression > add(std::shared_ptr< monty::ndarray< monty::rc_ptr< ::mosek::fusion::Variable >,1 > > _7941_vs);
      static  monty::rc_ptr< ::mosek::fusion::Expression > add_(monty::rc_ptr< ::mosek::fusion::Expression > _7974_e1,double _7975_m1,monty::rc_ptr< ::mosek::fusion::Expression > _7976_e2,double _7977_m2);
      static  monty::rc_ptr< ::mosek::fusion::Expression > transpose(monty::rc_ptr< ::mosek::fusion::Expression > _7988_e);
      static  monty::rc_ptr< ::mosek::fusion::Expression > mulElm(monty::rc_ptr< ::mosek::fusion::Matrix > _7989_m,monty::rc_ptr< ::mosek::fusion::Expression > _7990_expr);
      static  monty::rc_ptr< ::mosek::fusion::Expression > mulElm(monty::rc_ptr< ::mosek::fusion::NDSparseArray > _7991_spm,monty::rc_ptr< ::mosek::fusion::Expression > _7992_expr);
      static  monty::rc_ptr< ::mosek::fusion::Expression > mulElm(std::shared_ptr< monty::ndarray< double,2 > > _7993_a2,monty::rc_ptr< ::mosek::fusion::Expression > _7994_expr);
      static  monty::rc_ptr< ::mosek::fusion::Expression > mulElm(std::shared_ptr< monty::ndarray< double,1 > > _7995_a1,monty::rc_ptr< ::mosek::fusion::Expression > _7996_expr);
      static  monty::rc_ptr< ::mosek::fusion::Expression > mulElm(monty::rc_ptr< ::mosek::fusion::Expression > _7997_expr,monty::rc_ptr< ::mosek::fusion::Matrix > _7998_m);
      static  monty::rc_ptr< ::mosek::fusion::Expression > mulElm(monty::rc_ptr< ::mosek::fusion::Expression > _7999_expr,std::shared_ptr< monty::ndarray< double,2 > > _8000_a2);
      static  monty::rc_ptr< ::mosek::fusion::Expression > mulElm(monty::rc_ptr< ::mosek::fusion::Expression > _8001_expr,std::shared_ptr< monty::ndarray< double,1 > > _8002_a1);
      static  monty::rc_ptr< ::mosek::fusion::Expression > mulElm(monty::rc_ptr< ::mosek::fusion::Expression > _8003_expr,monty::rc_ptr< ::mosek::fusion::NDSparseArray > _8004_spm);
      static  monty::rc_ptr< ::mosek::fusion::Expression > mulElm(monty::rc_ptr< ::mosek::fusion::Parameter > _8005_p,monty::rc_ptr< ::mosek::fusion::Expression > _8006_expr);
      static  monty::rc_ptr< ::mosek::fusion::Expression > mulElm(monty::rc_ptr< ::mosek::fusion::Expression > _8007_expr,monty::rc_ptr< ::mosek::fusion::Parameter > _8008_p);
      static  monty::rc_ptr< ::mosek::fusion::Expression > sub(monty::rc_ptr< ::mosek::fusion::NDSparseArray > _8009_n,monty::rc_ptr< ::mosek::fusion::Expression > _8010_e2);
      static  monty::rc_ptr< ::mosek::fusion::Expression > sub(monty::rc_ptr< ::mosek::fusion::Expression > _8011_e1,monty::rc_ptr< ::mosek::fusion::NDSparseArray > _8012_n);
      static  monty::rc_ptr< ::mosek::fusion::Expression > sub(monty::rc_ptr< ::mosek::fusion::Matrix > _8013_m,monty::rc_ptr< ::mosek::fusion::Expression > _8014_e2);
      static  monty::rc_ptr< ::mosek::fusion::Expression > sub(monty::rc_ptr< ::mosek::fusion::Expression > _8015_e1,monty::rc_ptr< ::mosek::fusion::Matrix > _8016_m);
      static  monty::rc_ptr< ::mosek::fusion::Expression > sub(double _8017_c,monty::rc_ptr< ::mosek::fusion::Expression > _8018_e2);
      static  monty::rc_ptr< ::mosek::fusion::Expression > sub(monty::rc_ptr< ::mosek::fusion::Expression > _8019_e1,double _8020_c);
      static  monty::rc_ptr< ::mosek::fusion::Expression > sub(std::shared_ptr< monty::ndarray< double,2 > > _8021_a2,monty::rc_ptr< ::mosek::fusion::Expression > _8022_e2);
      static  monty::rc_ptr< ::mosek::fusion::Expression > sub(std::shared_ptr< monty::ndarray< double,1 > > _8023_a1,monty::rc_ptr< ::mosek::fusion::Expression > _8024_e2);
      static  monty::rc_ptr< ::mosek::fusion::Expression > sub(monty::rc_ptr< ::mosek::fusion::Expression > _8025_e1,std::shared_ptr< monty::ndarray< double,2 > > _8026_a2);
      static  monty::rc_ptr< ::mosek::fusion::Expression > sub(monty::rc_ptr< ::mosek::fusion::Expression > _8027_e1,std::shared_ptr< monty::ndarray< double,1 > > _8028_a1);
      static  monty::rc_ptr< ::mosek::fusion::Expression > sub(monty::rc_ptr< ::mosek::fusion::Expression > _8029_e1,monty::rc_ptr< ::mosek::fusion::Expression > _8030_e2);
      static  monty::rc_ptr< ::mosek::fusion::Expression > add(monty::rc_ptr< ::mosek::fusion::NDSparseArray > _8031_n,monty::rc_ptr< ::mosek::fusion::Expression > _8032_e2);
      static  monty::rc_ptr< ::mosek::fusion::Expression > add(monty::rc_ptr< ::mosek::fusion::Expression > _8033_e1,monty::rc_ptr< ::mosek::fusion::NDSparseArray > _8034_n);
      static  monty::rc_ptr< ::mosek::fusion::Expression > add(monty::rc_ptr< ::mosek::fusion::Matrix > _8035_m,monty::rc_ptr< ::mosek::fusion::Expression > _8036_e2);
      static  monty::rc_ptr< ::mosek::fusion::Expression > add(monty::rc_ptr< ::mosek::fusion::Expression > _8037_e1,monty::rc_ptr< ::mosek::fusion::Matrix > _8038_m);
      static  monty::rc_ptr< ::mosek::fusion::Expression > add(double _8039_c,monty::rc_ptr< ::mosek::fusion::Expression > _8040_e2);
      static  monty::rc_ptr< ::mosek::fusion::Expression > add(monty::rc_ptr< ::mosek::fusion::Expression > _8041_e1,double _8042_c);
      static  monty::rc_ptr< ::mosek::fusion::Expression > add(std::shared_ptr< monty::ndarray< double,2 > > _8043_a2,monty::rc_ptr< ::mosek::fusion::Expression > _8044_e2);
      static  monty::rc_ptr< ::mosek::fusion::Expression > add(std::shared_ptr< monty::ndarray< double,1 > > _8045_a1,monty::rc_ptr< ::mosek::fusion::Expression > _8046_e2);
      static  monty::rc_ptr< ::mosek::fusion::Expression > add(monty::rc_ptr< ::mosek::fusion::Expression > _8047_e1,std::shared_ptr< monty::ndarray< double,2 > > _8048_a2);
      static  monty::rc_ptr< ::mosek::fusion::Expression > add(monty::rc_ptr< ::mosek::fusion::Expression > _8049_e1,std::shared_ptr< monty::ndarray< double,1 > > _8050_a1);
      static  monty::rc_ptr< ::mosek::fusion::Expression > add(monty::rc_ptr< ::mosek::fusion::Expression > _8051_e1,monty::rc_ptr< ::mosek::fusion::Expression > _8052_e2);
      virtual /* override */ void eval(monty::rc_ptr< ::mosek::fusion::WorkStack > _8053_rs,monty::rc_ptr< ::mosek::fusion::WorkStack > _8054_ws,monty::rc_ptr< ::mosek::fusion::WorkStack > _8055_xs) ;
      static  void validateData(std::shared_ptr< monty::ndarray< int64_t,1 > > _8071_ptrb,std::shared_ptr< monty::ndarray< int64_t,1 > > _8072_subj,std::shared_ptr< monty::ndarray< double,1 > > _8073_cof,std::shared_ptr< monty::ndarray< double,1 > > _8074_bfix,std::shared_ptr< monty::ndarray< int32_t,1 > > _8075_shape,std::shared_ptr< monty::ndarray< int64_t,1 > > _8076_inst);
      static  monty::rc_ptr< ::mosek::fusion::Model > extractModel(std::shared_ptr< monty::ndarray< monty::rc_ptr< ::mosek::fusion::Variable >,1 > > _8089_v);
    }; // struct Expr;

    struct p_WorkStack
    {
      WorkStack * _pubthis;
      static mosek::fusion::p_WorkStack* _get_impl(mosek::fusion::WorkStack * _inst){ assert(_inst); assert(_inst->_impl); return _inst->_impl; }
      static mosek::fusion::p_WorkStack * _get_impl(mosek::fusion::WorkStack::t _inst) { return _get_impl(_inst.get()); }
      p_WorkStack(WorkStack * _pubthis);
      virtual ~p_WorkStack() { /* std::cout << "~p_WorkStack" << std::endl;*/ };
      int32_t code_base{};
      int32_t cconst_base{};
      int32_t codeptr_base{};
      int32_t cof_base{};
      int32_t nidxs_base{};
      int32_t sp_base{};
      int32_t shape_base{};
      int32_t ptr_base{};
      bool hassp{};
      int32_t ncodeatom{};
      int32_t nelem{};
      int32_t nnz{};
      int32_t nd{};
      int32_t pf64{};
      int32_t pi64{};
      int32_t pi32{};
      std::shared_ptr< monty::ndarray< double,1 > > f64{};
      std::shared_ptr< monty::ndarray< int64_t,1 > > i64{};
      std::shared_ptr< monty::ndarray< int32_t,1 > > i32{};

      virtual void destroy();

      static WorkStack::t _new_WorkStack();
      void _initialize();
      virtual std::string formatCurrent() ;
      virtual bool peek_hassp() ;
      virtual int32_t peek_nnz() ;
      virtual int32_t peek_nelem() ;
      virtual int32_t peek_dim(int32_t _7334_i) ;
      virtual int32_t peek_nd() ;
      virtual void alloc_expr(int32_t _7335_nd,int32_t _7336_nelem,int32_t _7337_nnz,bool _7338_hassp) ;
      virtual void alloc_expr(int32_t _7339_nd,int32_t _7340_nelem,int32_t _7341_nnz,bool _7342_hassp,int32_t _7343_ncodeatom) ;
      virtual void pop_expr() ;
      virtual void move_expr(monty::rc_ptr< ::mosek::fusion::WorkStack > _7344_to) ;
      virtual void peek_expr() ;
      virtual void ensure_sparsity() ;
      virtual void clear() ;
      virtual int32_t allocf64(int32_t _7359_n) ;
      virtual int32_t alloci64(int32_t _7362_n) ;
      virtual int32_t alloci32(int32_t _7365_n) ;
      virtual void pushf64(double _7368_v) ;
      virtual void pushi64(int64_t _7369_v) ;
      virtual void pushi32(int32_t _7370_v) ;
      virtual void ensuref64(int32_t _7371_n) ;
      virtual void ensurei64(int32_t _7374_n) ;
      virtual void ensurei32(int32_t _7377_n) ;
      virtual int32_t popf64(int32_t _7380_n) ;
      virtual int32_t popi64(int32_t _7381_n) ;
      virtual int32_t popi32(int32_t _7382_n) ;
      virtual void popf64(int32_t _7383_n,std::shared_ptr< monty::ndarray< double,1 > > _7384_r,int32_t _7385_ofs) ;
      virtual void popi64(int32_t _7386_n,std::shared_ptr< monty::ndarray< int64_t,1 > > _7387_r,int32_t _7388_ofs) ;
      virtual void popi32(int32_t _7389_n,std::shared_ptr< monty::ndarray< int32_t,1 > > _7390_r,int32_t _7391_ofs) ;
      virtual double popf64() ;
      virtual int64_t popi64() ;
      virtual int32_t popi32() ;
      virtual double peekf64() ;
      virtual int64_t peeki64() ;
      virtual int32_t peeki32() ;
      virtual double peekf64(int32_t _7392_i) ;
      virtual int64_t peeki64(int32_t _7393_i) ;
      virtual int32_t peeki32(int32_t _7394_i) ;
    }; // struct WorkStack;

    struct p_SymmetricMatrix
    {
      SymmetricMatrix * _pubthis;
      static mosek::fusion::p_SymmetricMatrix* _get_impl(mosek::fusion::SymmetricMatrix * _inst){ assert(_inst); assert(_inst->_impl); return _inst->_impl; }
      static mosek::fusion::p_SymmetricMatrix * _get_impl(mosek::fusion::SymmetricMatrix::t _inst) { return _get_impl(_inst.get()); }
      p_SymmetricMatrix(SymmetricMatrix * _pubthis);
      virtual ~p_SymmetricMatrix() { /* std::cout << "~p_SymmetricMatrix" << std::endl;*/ };
      int32_t nnz{};
      double scale{};
      std::shared_ptr< monty::ndarray< double,1 > > vval{};
      std::shared_ptr< monty::ndarray< int32_t,1 > > vsubj{};
      std::shared_ptr< monty::ndarray< int32_t,1 > > vsubi{};
      std::shared_ptr< monty::ndarray< double,1 > > uval{};
      std::shared_ptr< monty::ndarray< int32_t,1 > > usubj{};
      std::shared_ptr< monty::ndarray< int32_t,1 > > usubi{};
      int32_t d1{};
      int32_t d0{};

      virtual void destroy();

      static SymmetricMatrix::t _new_SymmetricMatrix(int32_t _8102_dim0,int32_t _8103_dim1,std::shared_ptr< monty::ndarray< int32_t,1 > > _8104_usubi,std::shared_ptr< monty::ndarray< int32_t,1 > > _8105_usubj,std::shared_ptr< monty::ndarray< double,1 > > _8106_uval,std::shared_ptr< monty::ndarray< int32_t,1 > > _8107_vsubi,std::shared_ptr< monty::ndarray< int32_t,1 > > _8108_vsubj,std::shared_ptr< monty::ndarray< double,1 > > _8109_vval,double _8110_scale);
      void _initialize(int32_t _8102_dim0,int32_t _8103_dim1,std::shared_ptr< monty::ndarray< int32_t,1 > > _8104_usubi,std::shared_ptr< monty::ndarray< int32_t,1 > > _8105_usubj,std::shared_ptr< monty::ndarray< double,1 > > _8106_uval,std::shared_ptr< monty::ndarray< int32_t,1 > > _8107_vsubi,std::shared_ptr< monty::ndarray< int32_t,1 > > _8108_vsubj,std::shared_ptr< monty::ndarray< double,1 > > _8109_vval,double _8110_scale);
      static  monty::rc_ptr< ::mosek::fusion::SymmetricMatrix > rankOne(int32_t _8111_n,std::shared_ptr< monty::ndarray< int32_t,1 > > _8112_sub,std::shared_ptr< monty::ndarray< double,1 > > _8113_v);
      static  monty::rc_ptr< ::mosek::fusion::SymmetricMatrix > rankOne(std::shared_ptr< monty::ndarray< double,1 > > _8121_v);
      static  monty::rc_ptr< ::mosek::fusion::SymmetricMatrix > antiDiag(std::shared_ptr< monty::ndarray< double,1 > > _8129_vals);
      static  monty::rc_ptr< ::mosek::fusion::SymmetricMatrix > diag(std::shared_ptr< monty::ndarray< double,1 > > _8136_vals);
      virtual monty::rc_ptr< ::mosek::fusion::SymmetricMatrix > __mosek_2fusion_2SymmetricMatrix__add(monty::rc_ptr< ::mosek::fusion::SymmetricMatrix > _8142_m) ;
      virtual monty::rc_ptr< ::mosek::fusion::SymmetricMatrix > __mosek_2fusion_2SymmetricMatrix__sub(monty::rc_ptr< ::mosek::fusion::SymmetricMatrix > _8162_m) ;
      virtual monty::rc_ptr< ::mosek::fusion::SymmetricMatrix > __mosek_2fusion_2SymmetricMatrix__mul(double _8163_v) ;
      virtual int32_t getdim() ;
    }; // struct SymmetricMatrix;

    struct p_NDSparseArray
    {
      NDSparseArray * _pubthis;
      static mosek::fusion::p_NDSparseArray* _get_impl(mosek::fusion::NDSparseArray * _inst){ assert(_inst); assert(_inst->_impl); return _inst->_impl; }
      static mosek::fusion::p_NDSparseArray * _get_impl(mosek::fusion::NDSparseArray::t _inst) { return _get_impl(_inst.get()); }
      p_NDSparseArray(NDSparseArray * _pubthis);
      virtual ~p_NDSparseArray() { /* std::cout << "~p_NDSparseArray" << std::endl;*/ };
      std::shared_ptr< monty::ndarray< double,1 > > cof{};
      std::shared_ptr< monty::ndarray< int64_t,1 > > inst{};
      std::shared_ptr< monty::ndarray< int32_t,1 > > dims{};
      int64_t size{};

      virtual void destroy();

      static NDSparseArray::t _new_NDSparseArray(std::shared_ptr< monty::ndarray< int32_t,1 > > _8164_dims_,std::shared_ptr< monty::ndarray< int32_t,2 > > _8165_sub,std::shared_ptr< monty::ndarray< double,1 > > _8166_cof_);
      void _initialize(std::shared_ptr< monty::ndarray< int32_t,1 > > _8164_dims_,std::shared_ptr< monty::ndarray< int32_t,2 > > _8165_sub,std::shared_ptr< monty::ndarray< double,1 > > _8166_cof_);
      static NDSparseArray::t _new_NDSparseArray(std::shared_ptr< monty::ndarray< int32_t,1 > > _8187_dims_,std::shared_ptr< monty::ndarray< int64_t,1 > > _8188_inst_,std::shared_ptr< monty::ndarray< double,1 > > _8189_cof_);
      void _initialize(std::shared_ptr< monty::ndarray< int32_t,1 > > _8187_dims_,std::shared_ptr< monty::ndarray< int64_t,1 > > _8188_inst_,std::shared_ptr< monty::ndarray< double,1 > > _8189_cof_);
      static NDSparseArray::t _new_NDSparseArray(monty::rc_ptr< ::mosek::fusion::Matrix > _8205_m);
      void _initialize(monty::rc_ptr< ::mosek::fusion::Matrix > _8205_m);
      static  monty::rc_ptr< ::mosek::fusion::NDSparseArray > make(monty::rc_ptr< ::mosek::fusion::Matrix > _8213_m);
      static  monty::rc_ptr< ::mosek::fusion::NDSparseArray > make(std::shared_ptr< monty::ndarray< int32_t,1 > > _8214_dims,std::shared_ptr< monty::ndarray< int64_t,1 > > _8215_inst,std::shared_ptr< monty::ndarray< double,1 > > _8216_cof);
      static  monty::rc_ptr< ::mosek::fusion::NDSparseArray > make(std::shared_ptr< monty::ndarray< int32_t,1 > > _8217_dims,std::shared_ptr< monty::ndarray< int32_t,2 > > _8218_sub,std::shared_ptr< monty::ndarray< double,1 > > _8219_cof);
    }; // struct NDSparseArray;

    struct p_Matrix
    {
      Matrix * _pubthis;
      static mosek::fusion::p_Matrix* _get_impl(mosek::fusion::Matrix * _inst){ assert(_inst); assert(_inst->_impl); return _inst->_impl; }
      static mosek::fusion::p_Matrix * _get_impl(mosek::fusion::Matrix::t _inst) { return _get_impl(_inst.get()); }
      p_Matrix(Matrix * _pubthis);
      virtual ~p_Matrix() { /* std::cout << "~p_Matrix" << std::endl;*/ };
      int32_t dimj{};
      int32_t dimi{};

      virtual void destroy();

      static Matrix::t _new_Matrix(int32_t _8289_di,int32_t _8290_dj);
      void _initialize(int32_t _8289_di,int32_t _8290_dj);
      virtual /* override */ std::string toString() ;
      virtual void switchDims() ;
      static  monty::rc_ptr< ::mosek::fusion::Matrix > diag(int32_t _8292_num,monty::rc_ptr< ::mosek::fusion::Matrix > _8293_mv);
      static  monty::rc_ptr< ::mosek::fusion::Matrix > diag(std::shared_ptr< monty::ndarray< monty::rc_ptr< ::mosek::fusion::Matrix >,1 > > _8295_md);
      static  monty::rc_ptr< ::mosek::fusion::Matrix > antidiag(int32_t _8313_n,double _8314_val,int32_t _8315_k);
      static  monty::rc_ptr< ::mosek::fusion::Matrix > antidiag(int32_t _8316_n,double _8317_val);
      static  monty::rc_ptr< ::mosek::fusion::Matrix > diag(int32_t _8318_n,double _8319_val,int32_t _8320_k);
      static  monty::rc_ptr< ::mosek::fusion::Matrix > diag(int32_t _8321_n,double _8322_val);
      static  monty::rc_ptr< ::mosek::fusion::Matrix > antidiag(std::shared_ptr< monty::ndarray< double,1 > > _8323_d,int32_t _8324_k);
      static  monty::rc_ptr< ::mosek::fusion::Matrix > antidiag(std::shared_ptr< monty::ndarray< double,1 > > _8334_d);
      static  monty::rc_ptr< ::mosek::fusion::Matrix > diag(std::shared_ptr< monty::ndarray< double,1 > > _8335_d,int32_t _8336_k);
      static  monty::rc_ptr< ::mosek::fusion::Matrix > diag(std::shared_ptr< monty::ndarray< double,1 > > _8344_d);
      static  monty::rc_ptr< ::mosek::fusion::Matrix > ones(int32_t _8345_n,int32_t _8346_m);
      static  monty::rc_ptr< ::mosek::fusion::Matrix > eye(int32_t _8347_n);
      static  monty::rc_ptr< ::mosek::fusion::Matrix > dense(monty::rc_ptr< ::mosek::fusion::Matrix > _8349_other);
      static  monty::rc_ptr< ::mosek::fusion::Matrix > dense(int32_t _8350_dimi,int32_t _8351_dimj,double _8352_value);
      static  monty::rc_ptr< ::mosek::fusion::Matrix > dense(int32_t _8353_dimi,int32_t _8354_dimj,std::shared_ptr< monty::ndarray< double,1 > > _8355_data);
      static  monty::rc_ptr< ::mosek::fusion::Matrix > dense(std::shared_ptr< monty::ndarray< double,2 > > _8356_data);
      static  monty::rc_ptr< ::mosek::fusion::Matrix > sparse(monty::rc_ptr< ::mosek::fusion::Matrix > _8357_mx);
      static  monty::rc_ptr< ::mosek::fusion::Matrix > sparse(std::shared_ptr< monty::ndarray< std::shared_ptr< monty::ndarray< monty::rc_ptr< ::mosek::fusion::Matrix >,1 > >,1 > > _8361_blocks);
      static  monty::rc_ptr< ::mosek::fusion::Matrix > sparse(std::shared_ptr< monty::ndarray< double,2 > > _8392_data);
      static  monty::rc_ptr< ::mosek::fusion::Matrix > sparse(int32_t _8402_nrow,int32_t _8403_ncol);
      static  monty::rc_ptr< ::mosek::fusion::Matrix > sparse(int32_t _8404_nrow,int32_t _8405_ncol,std::shared_ptr< monty::ndarray< int32_t,1 > > _8406_subi,std::shared_ptr< monty::ndarray< int32_t,1 > > _8407_subj,double _8408_val);
      static  monty::rc_ptr< ::mosek::fusion::Matrix > sparse(std::shared_ptr< monty::ndarray< int32_t,1 > > _8410_subi,std::shared_ptr< monty::ndarray< int32_t,1 > > _8411_subj,double _8412_val);
      static  monty::rc_ptr< ::mosek::fusion::Matrix > sparse(std::shared_ptr< monty::ndarray< int32_t,1 > > _8417_subi,std::shared_ptr< monty::ndarray< int32_t,1 > > _8418_subj,std::shared_ptr< monty::ndarray< double,1 > > _8419_val);
      static  monty::rc_ptr< ::mosek::fusion::Matrix > sparse(int32_t _8424_nrow,int32_t _8425_ncol,std::shared_ptr< monty::ndarray< int32_t,1 > > _8426_subi,std::shared_ptr< monty::ndarray< int32_t,1 > > _8427_subj,std::shared_ptr< monty::ndarray< double,1 > > _8428_val);
      virtual double get(int32_t _8433_i,int32_t _8434_j) { throw monty::AbstractClassError("Call to abstract method"); }
      virtual bool isSparse() { throw monty::AbstractClassError("Call to abstract method"); }
      virtual std::shared_ptr< monty::ndarray< double,1 > > getDataAsArray() { throw monty::AbstractClassError("Call to abstract method"); }
      virtual void getDataAsTriplets(std::shared_ptr< monty::ndarray< int32_t,1 > > _8435_subi,std::shared_ptr< monty::ndarray< int32_t,1 > > _8436_subj,std::shared_ptr< monty::ndarray< double,1 > > _8437_val) { throw monty::AbstractClassError("Call to abstract method"); }
      virtual monty::rc_ptr< ::mosek::fusion::Matrix > __mosek_2fusion_2Matrix__transpose() { throw monty::AbstractClassError("Call to abstract method"); }
      virtual int64_t numNonzeros() { throw monty::AbstractClassError("Call to abstract method"); }
      virtual int32_t numColumns() ;
      virtual int32_t numRows() ;
    }; // struct Matrix;

    struct p_DenseMatrix : public ::mosek::fusion::p_Matrix
    {
      DenseMatrix * _pubthis;
      static mosek::fusion::p_DenseMatrix* _get_impl(mosek::fusion::DenseMatrix * _inst){ return static_cast< mosek::fusion::p_DenseMatrix* >(mosek::fusion::p_Matrix::_get_impl(_inst)); }
      static mosek::fusion::p_DenseMatrix * _get_impl(mosek::fusion::DenseMatrix::t _inst) { return _get_impl(_inst.get()); }
      p_DenseMatrix(DenseMatrix * _pubthis);
      virtual ~p_DenseMatrix() { /* std::cout << "~p_DenseMatrix" << std::endl;*/ };
      int64_t nnz{};
      std::shared_ptr< monty::ndarray< double,1 > > data{};

      virtual void destroy();

      static DenseMatrix::t _new_DenseMatrix(int32_t _8220_dimi_,int32_t _8221_dimj_,std::shared_ptr< monty::ndarray< double,1 > > _8222_cof);
      void _initialize(int32_t _8220_dimi_,int32_t _8221_dimj_,std::shared_ptr< monty::ndarray< double,1 > > _8222_cof);
      static DenseMatrix::t _new_DenseMatrix(monty::rc_ptr< ::mosek::fusion::Matrix > _8223_m_);
      void _initialize(monty::rc_ptr< ::mosek::fusion::Matrix > _8223_m_);
      static DenseMatrix::t _new_DenseMatrix(std::shared_ptr< monty::ndarray< double,2 > > _8228_d);
      void _initialize(std::shared_ptr< monty::ndarray< double,2 > > _8228_d);
      static DenseMatrix::t _new_DenseMatrix(int32_t _8231_dimi_,int32_t _8232_dimj_,double _8233_value_);
      void _initialize(int32_t _8231_dimi_,int32_t _8232_dimj_,double _8233_value_);
      virtual /* override */ std::string toString() ;
      virtual /* override */ monty::rc_ptr< ::mosek::fusion::Matrix > __mosek_2fusion_2DenseMatrix__transpose() ;
      virtual monty::rc_ptr< ::mosek::fusion::Matrix > __mosek_2fusion_2Matrix__transpose() { return __mosek_2fusion_2DenseMatrix__transpose(); }
      virtual /* override */ bool isSparse() ;
      virtual /* override */ std::shared_ptr< monty::ndarray< double,1 > > getDataAsArray() ;
      virtual /* override */ void getDataAsTriplets(std::shared_ptr< monty::ndarray< int32_t,1 > > _8246_subi,std::shared_ptr< monty::ndarray< int32_t,1 > > _8247_subj,std::shared_ptr< monty::ndarray< double,1 > > _8248_cof) ;
      virtual /* override */ double get(int32_t _8252_i,int32_t _8253_j) ;
      virtual /* override */ int64_t numNonzeros() ;
    }; // struct DenseMatrix;

    struct p_SparseMatrix : public ::mosek::fusion::p_Matrix
    {
      SparseMatrix * _pubthis;
      static mosek::fusion::p_SparseMatrix* _get_impl(mosek::fusion::SparseMatrix * _inst){ return static_cast< mosek::fusion::p_SparseMatrix* >(mosek::fusion::p_Matrix::_get_impl(_inst)); }
      static mosek::fusion::p_SparseMatrix * _get_impl(mosek::fusion::SparseMatrix::t _inst) { return _get_impl(_inst.get()); }
      p_SparseMatrix(SparseMatrix * _pubthis);
      virtual ~p_SparseMatrix() { /* std::cout << "~p_SparseMatrix" << std::endl;*/ };
      int64_t nnz{};
      std::shared_ptr< monty::ndarray< double,1 > > val{};
      std::shared_ptr< monty::ndarray< int32_t,1 > > subj{};
      std::shared_ptr< monty::ndarray< int32_t,1 > > subi{};

      virtual void destroy();

      static SparseMatrix::t _new_SparseMatrix(int32_t _8254_dimi_,int32_t _8255_dimj_,std::shared_ptr< monty::ndarray< int32_t,1 > > _8256_subi_,std::shared_ptr< monty::ndarray< int32_t,1 > > _8257_subj_,std::shared_ptr< monty::ndarray< double,1 > > _8258_val_,int64_t _8259_nelm);
      void _initialize(int32_t _8254_dimi_,int32_t _8255_dimj_,std::shared_ptr< monty::ndarray< int32_t,1 > > _8256_subi_,std::shared_ptr< monty::ndarray< int32_t,1 > > _8257_subj_,std::shared_ptr< monty::ndarray< double,1 > > _8258_val_,int64_t _8259_nelm);
      static SparseMatrix::t _new_SparseMatrix(int32_t _8265_dimi_,int32_t _8266_dimj_,std::shared_ptr< monty::ndarray< int32_t,1 > > _8267_subi_,std::shared_ptr< monty::ndarray< int32_t,1 > > _8268_subj_,std::shared_ptr< monty::ndarray< double,1 > > _8269_val_);
      void _initialize(int32_t _8265_dimi_,int32_t _8266_dimj_,std::shared_ptr< monty::ndarray< int32_t,1 > > _8267_subi_,std::shared_ptr< monty::ndarray< int32_t,1 > > _8268_subj_,std::shared_ptr< monty::ndarray< double,1 > > _8269_val_);
      virtual std::shared_ptr< monty::ndarray< int64_t,1 > > formPtrb() ;
      virtual /* override */ std::string toString() ;
      virtual /* override */ int64_t numNonzeros() ;
      virtual /* override */ monty::rc_ptr< ::mosek::fusion::Matrix > __mosek_2fusion_2SparseMatrix__transpose() ;
      virtual monty::rc_ptr< ::mosek::fusion::Matrix > __mosek_2fusion_2Matrix__transpose() { return __mosek_2fusion_2SparseMatrix__transpose(); }
      virtual /* override */ bool isSparse() ;
      virtual /* override */ std::shared_ptr< monty::ndarray< double,1 > > getDataAsArray() ;
      virtual /* override */ void getDataAsTriplets(std::shared_ptr< monty::ndarray< int32_t,1 > > _8281_subi_,std::shared_ptr< monty::ndarray< int32_t,1 > > _8282_subj_,std::shared_ptr< monty::ndarray< double,1 > > _8283_cof_) ;
      virtual /* override */ double get(int32_t _8284_i,int32_t _8285_j) ;
    }; // struct SparseMatrix;

    struct p_LinkedBlocks
    {
      LinkedBlocks * _pubthis;
      static mosek::fusion::p_LinkedBlocks* _get_impl(mosek::fusion::LinkedBlocks * _inst){ assert(_inst); assert(_inst->_impl); return _inst->_impl; }
      static mosek::fusion::p_LinkedBlocks * _get_impl(mosek::fusion::LinkedBlocks::t _inst) { return _get_impl(_inst.get()); }
      p_LinkedBlocks(LinkedBlocks * _pubthis);
      virtual ~p_LinkedBlocks() { /* std::cout << "~p_LinkedBlocks" << std::endl;*/ };
      std::shared_ptr< monty::ndarray< int32_t,1 > > bfirst{};
      std::shared_ptr< monty::ndarray< int32_t,1 > > bsize{};
      monty::rc_ptr< ::mosek::fusion::LinkedInts > blocks{};
      monty::rc_ptr< ::mosek::fusion::LinkedInts > ints{};

      virtual void destroy();

      static LinkedBlocks::t _new_LinkedBlocks();
      void _initialize();
      static LinkedBlocks::t _new_LinkedBlocks(int32_t _8463_n);
      void _initialize(int32_t _8463_n);
      static LinkedBlocks::t _new_LinkedBlocks(monty::rc_ptr< ::mosek::fusion::LinkedBlocks > _8465_other);
      void _initialize(monty::rc_ptr< ::mosek::fusion::LinkedBlocks > _8465_other);
      virtual void free(int32_t _8466_bkey) ;
      virtual int32_t alloc(int32_t _8468_size) ;
      virtual int32_t maxidx(int32_t _8473_bkey) ;
      virtual int32_t numallocated() ;
      virtual void get(int32_t _8474_bkey,std::shared_ptr< monty::ndarray< int32_t,1 > > _8475_target,int32_t _8476_offset) ;
      virtual int32_t numblocks() ;
      virtual int32_t blocksize(int32_t _8477_bkey) ;
      virtual int32_t block_capacity() ;
      virtual int32_t capacity() ;
      virtual bool validate() ;
    }; // struct LinkedBlocks;

    struct p_LinkedInts
    {
      LinkedInts * _pubthis;
      static mosek::fusion::p_LinkedInts* _get_impl(mosek::fusion::LinkedInts * _inst){ assert(_inst); assert(_inst->_impl); return _inst->_impl; }
      static mosek::fusion::p_LinkedInts * _get_impl(mosek::fusion::LinkedInts::t _inst) { return _get_impl(_inst.get()); }
      p_LinkedInts(LinkedInts * _pubthis);
      virtual ~p_LinkedInts() { /* std::cout << "~p_LinkedInts" << std::endl;*/ };
      int32_t nfree{};
      int32_t last_free{};
      int32_t first_free{};
      int32_t first_used{};
      std::shared_ptr< monty::ndarray< int32_t,1 > > prev{};
      std::shared_ptr< monty::ndarray< int32_t,1 > > next{};

      virtual void destroy();

      static LinkedInts::t _new_LinkedInts(int32_t _8478_cap_);
      void _initialize(int32_t _8478_cap_);
      static LinkedInts::t _new_LinkedInts();
      void _initialize();
      static LinkedInts::t _new_LinkedInts(monty::rc_ptr< ::mosek::fusion::LinkedInts > _8481_other);
      void _initialize(monty::rc_ptr< ::mosek::fusion::LinkedInts > _8481_other);
      virtual void free(int32_t _8482_i,int32_t _8483_num) ;
      virtual int32_t alloc() ;
      virtual int32_t alloc(int32_t _8489_n) ;
      virtual void alloc(int32_t _8490_num,std::shared_ptr< monty::ndarray< int32_t,1 > > _8491_target,int32_t _8492_offset) ;
      virtual void get(int32_t _8495_i,int32_t _8496_num,std::shared_ptr< monty::ndarray< int32_t,1 > > _8497_target,int32_t _8498_offset) ;
      virtual int32_t numallocated() ;
      virtual int32_t maxidx(int32_t _8501_i,int32_t _8502_num) ;
      virtual int32_t allocblock(int32_t _8506_num) ;
      virtual void recap(int32_t _8512_ncap) ;
      virtual int32_t capacity() ;
      virtual bool validate() ;
    }; // struct LinkedInts;

    struct p_Parameters
    {
      Parameters * _pubthis;
      static mosek::fusion::p_Parameters* _get_impl(mosek::fusion::Parameters * _inst){ assert(_inst); assert(_inst->_impl); return _inst->_impl; }
      static mosek::fusion::p_Parameters * _get_impl(mosek::fusion::Parameters::t _inst) { return _get_impl(_inst.get()); }
      p_Parameters(Parameters * _pubthis);
      virtual ~p_Parameters() { /* std::cout << "~p_Parameters" << std::endl;*/ };

      virtual void destroy();

      static  void setParameter(monty::rc_ptr< ::mosek::fusion::Model > _8521_M,const std::string &  _8522_name,double _8523_value);
      static  void setParameter(monty::rc_ptr< ::mosek::fusion::Model > _8625_M,const std::string &  _8626_name,int32_t _8627_value);
      static  void setParameter(monty::rc_ptr< ::mosek::fusion::Model > _8729_M,const std::string &  _8730_name,const std::string &  _8731_value);
      static  int32_t string_to_variabletype_value(const std::string &  _8993_v);
      static  int32_t string_to_value_value(const std::string &  _8994_v);
      static  int32_t string_to_streamtype_value(const std::string &  _8995_v);
      static  int32_t string_to_startpointtype_value(const std::string &  _8996_v);
      static  int32_t string_to_stakey_value(const std::string &  _8997_v);
      static  int32_t string_to_sparam_value(const std::string &  _8998_v);
      static  int32_t string_to_solveform_value(const std::string &  _8999_v);
      static  int32_t string_to_soltype_value(const std::string &  _9000_v);
      static  int32_t string_to_solsta_value(const std::string &  _9001_v);
      static  int32_t string_to_solitem_value(const std::string &  _9002_v);
      static  int32_t string_to_simseltype_value(const std::string &  _9003_v);
      static  int32_t string_to_sensitivitytype_value(const std::string &  _9004_v);
      static  int32_t string_to_scalingmethod_value(const std::string &  _9005_v);
      static  int32_t string_to_scalingtype_value(const std::string &  _9006_v);
      static  int32_t string_to_rescodetype_value(const std::string &  _9007_v);
      static  int32_t string_to_rescode_value(const std::string &  _9008_v);
      static  int32_t string_to_xmlwriteroutputtype_value(const std::string &  _9009_v);
      static  int32_t string_to_prosta_value(const std::string &  _9010_v);
      static  int32_t string_to_problemtype_value(const std::string &  _9011_v);
      static  int32_t string_to_problemitem_value(const std::string &  _9012_v);
      static  int32_t string_to_parametertype_value(const std::string &  _9013_v);
      static  int32_t string_to_presolvemode_value(const std::string &  _9014_v);
      static  int32_t string_to_orderingtype_value(const std::string &  _9015_v);
      static  int32_t string_to_optimizertype_value(const std::string &  _9016_v);
      static  int32_t string_to_onoffkey_value(const std::string &  _9017_v);
      static  int32_t string_to_objsense_value(const std::string &  _9018_v);
      static  int32_t string_to_mpsformat_value(const std::string &  _9019_v);
      static  int32_t string_to_miovarseltype_value(const std::string &  _9020_v);
      static  int32_t string_to_mionodeseltype_value(const std::string &  _9021_v);
      static  int32_t string_to_miomode_value(const std::string &  _9022_v);
      static  int32_t string_to_miocontsoltype_value(const std::string &  _9023_v);
      static  int32_t string_to_miodatapermmethod_value(const std::string &  _9024_v);
      static  int32_t string_to_miqcqoreformmethod_value(const std::string &  _9025_v);
      static  int32_t string_to_branchdir_value(const std::string &  _9026_v);
      static  int32_t string_to_iparam_value(const std::string &  _9027_v);
      static  int32_t string_to_iomode_value(const std::string &  _9028_v);
      static  int32_t string_to_internal_iinf_value(const std::string &  _9029_v);
      static  int32_t string_to_internal_dinf_value(const std::string &  _9030_v);
      static  int32_t string_to_inftype_value(const std::string &  _9031_v);
      static  int32_t string_to_iinfitem_value(const std::string &  _9032_v);
      static  int32_t string_to_internal_liinf_value(const std::string &  _9033_v);
      static  int32_t string_to_liinfitem_value(const std::string &  _9034_v);
      static  int32_t string_to_dparam_value(const std::string &  _9035_v);
      static  int32_t string_to_feature_value(const std::string &  _9036_v);
      static  int32_t string_to_dinfitem_value(const std::string &  _9037_v);
      static  int32_t string_to_solformat_value(const std::string &  _9038_v);
      static  int32_t string_to_dataformat_value(const std::string &  _9039_v);
      static  int32_t string_to_symmattype_value(const std::string &  _9040_v);
      static  int32_t string_to_nametype_value(const std::string &  _9041_v);
      static  int32_t string_to_domaintype_value(const std::string &  _9042_v);
      static  int32_t string_to_conetype_value(const std::string &  _9043_v);
      static  int32_t string_to_compresstype_value(const std::string &  _9044_v);
      static  int32_t string_to_callbackcode_value(const std::string &  _9045_v);
      static  int32_t string_to_purify_value(const std::string &  _9046_v);
      static  int32_t string_to_intpnthotstart_value(const std::string &  _9047_v);
      static  int32_t string_to_simhotstart_value(const std::string &  _9048_v);
      static  int32_t string_to_simdupvec_value(const std::string &  _9049_v);
      static  int32_t string_to_simreform_value(const std::string &  _9050_v);
      static  int32_t string_to_uplo_value(const std::string &  _9051_v);
      static  int32_t string_to_transpose_value(const std::string &  _9052_v);
      static  int32_t string_to_simdegen_value(const std::string &  _9053_v);
      static  int32_t string_to_mark_value(const std::string &  _9054_v);
      static  int32_t string_to_boundkey_value(const std::string &  _9055_v);
      static  int32_t string_to_basindtype_value(const std::string &  _9056_v);
      static  int32_t string_to_language_value(const std::string &  _9057_v);
    }; // struct Parameters;

  }
}
namespace mosek
{
  namespace fusion
  {
    namespace Utils
    {
      // mosek.fusion.Utils.IntMap from file 'src/fusion/cxx/IntMap_p.h'
      struct p_IntMap
      {
        IntMap * _pubself;
      
        static p_IntMap * _get_impl(IntMap * _inst) { return _inst->_impl.get(); }
      
        p_IntMap(IntMap * _pubself) : _pubself(_pubself) {}
      
        static IntMap::t _new_IntMap() { return new IntMap(); }
      
        ::std::unordered_map<int64_t,int> m;
      
        bool hasItem (int64_t key) { return m.find(key) != m.end(); }
        int  getItem (int64_t key) { return m.find(key)->second; } // will probably throw something or crash of no such key
        void setItem (int64_t key, int val) { m[key] = val; }
      
        std::shared_ptr<monty::ndarray<int64_t,1>> keys()
        {
          size_t size = m.size();
          auto res = std::shared_ptr<monty::ndarray<int64_t,1>>(new monty::ndarray<int64_t,1>(monty::shape((int)size)));
      
          ptrdiff_t i = 0;
          for (auto it = m.begin(); it != m.end(); ++it)
            (*res)[i++] = it->first;
      
          return res;
        }
      
        std::shared_ptr<monty::ndarray<int,1>> values()
        {
          size_t size = m.size();
          auto res = std::shared_ptr<monty::ndarray<int,1>>(new monty::ndarray<int,1>(monty::shape((int)size)));
      
          ptrdiff_t i = 0;
          for (auto it = m.begin(); it != m.end(); ++it)
            (*res)[i++] = it->second;
      
          return res;
        }
      
        IntMap::t clone();
        IntMap::t __mosek_2fusion_2Utils_2IntMap__clone();
      };
      
      
      
      struct p_StringIntMap
      {
        StringIntMap * _pubself;
      
        static p_StringIntMap * _get_impl(StringIntMap * _inst) { return _inst->_impl.get(); }
      
        p_StringIntMap(StringIntMap * _pubself) : _pubself(_pubself) {}
      
        static StringIntMap::t _new_StringIntMap() { return new StringIntMap(); }
      
        ::std::unordered_map<std::string,int> m;
      
        bool hasItem (const std::string & key) { return m.find(key) != m.end(); }
        int  getItem (const std::string & key) { return m.find(key)->second; } // will probably throw something or crash of no such key
        void setItem (const std::string & key, int val) { m[key] = val; }
      
        std::shared_ptr<monty::ndarray<std::string,1>> keys()
        {
          size_t size = m.size();
          auto res = std::shared_ptr<monty::ndarray<std::string,1>>(new monty::ndarray<std::string,1>(monty::shape((int)size)));
      
          ptrdiff_t i = 0;
          for (auto it = m.begin(); it != m.end(); ++it)
            (*res)[i++] = it->first;
      
          return res;
        }
      
        std::shared_ptr<monty::ndarray<int,1>> values()
        {
          size_t size = m.size();
          auto res = std::shared_ptr<monty::ndarray<int,1>>(new monty::ndarray<int,1>(monty::shape((int)size)));
      
          ptrdiff_t i = 0;
          for (auto it = m.begin(); it != m.end(); ++it)
            (*res)[i++] = it->second;
      
          return res;
        }
      
        StringIntMap::t clone();
        StringIntMap::t __mosek_2fusion_2Utils_2StringIntMap__clone();
      };
      // End of file 'src/fusion/cxx/IntMap_p.h'
      // mosek.fusion.Utils.StringBuffer from file 'src/fusion/cxx/StringBuffer_p.h'
      // namespace mosek::fusion::Utils
      struct p_StringBuffer
      {
        StringBuffer * _pubthis;
        std::stringstream ss;
      
        p_StringBuffer(StringBuffer * _pubthis) : _pubthis(_pubthis) {}
      
        static p_StringBuffer * _get_impl(StringBuffer::t ptr) { return ptr->_impl.get(); }
        static p_StringBuffer * _get_impl(StringBuffer * ptr) { return ptr->_impl.get(); }
      
        static StringBuffer::t _new_StringBuffer() { return new StringBuffer(); }
      
        StringBuffer::t clear ();
      
        StringBuffer::t a (const monty::ndarray<std::string,1> & val);
        StringBuffer::t a (const monty::ndarray<int,1> & val);
        StringBuffer::t a (const monty::ndarray<int64_t,1> & val);
        StringBuffer::t a (const monty::ndarray<double,1> & val);
      
      
        StringBuffer::t a (const int & val);
        StringBuffer::t a (const int64_t & val);
        StringBuffer::t a (const double & val);
        StringBuffer::t a (const bool & val);
        StringBuffer::t a (const std::string & val);
      
        StringBuffer::t lf ();
        StringBuffer::t __mosek_2fusion_2Utils_2StringBuffer__clear ();
      
        StringBuffer::t __mosek_2fusion_2Utils_2StringBuffer__a (const monty::ndarray<std::string,1> & val);
        StringBuffer::t __mosek_2fusion_2Utils_2StringBuffer__a (const monty::ndarray<int,1> & val);
        StringBuffer::t __mosek_2fusion_2Utils_2StringBuffer__a (const monty::ndarray<int64_t,1> & val);
        StringBuffer::t __mosek_2fusion_2Utils_2StringBuffer__a (const monty::ndarray<double,1> & val);
      
        StringBuffer::t __mosek_2fusion_2Utils_2StringBuffer__a (const int & val);
        StringBuffer::t __mosek_2fusion_2Utils_2StringBuffer__a (const int64_t & val);
        StringBuffer::t __mosek_2fusion_2Utils_2StringBuffer__a (const double & val);
        StringBuffer::t __mosek_2fusion_2Utils_2StringBuffer__a (const bool & val);
        StringBuffer::t __mosek_2fusion_2Utils_2StringBuffer__a (const std::string & val);
      
        StringBuffer::t __mosek_2fusion_2Utils_2StringBuffer__lf ();
      
        std::string toString () const;
      };
      // End of file 'src/fusion/cxx/StringBuffer_p.h'
    }
  }
}
#endif
