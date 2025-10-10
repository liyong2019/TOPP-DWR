#ifndef __FUSION_H__
#define __FUSION_H__
#include "cmath"
#include "fusion_fwd.h"
#include "monty.h"
#include "mosek.h"
#include "mosektask.h"
#include "stdlib.h"
namespace mosek {
namespace fusion {
enum class Opcode { NOP, PARAM, CONSTVAL, ADD, NEG, MUL, INV, SUM, ZERO };
std::ostream &operator<<(std::ostream &os, Opcode val);
enum class RelationKey { EqualsTo, LessThan, GreaterThan, IsFree, InRange };
std::ostream &operator<<(std::ostream &os, RelationKey val);
enum class PSDKey { IsSymPSD, IsTrilPSD };
std::ostream &operator<<(std::ostream &os, PSDKey val);
enum class QConeKey {
  InQCone,
  InRotatedQCone,
  InPExpCone,
  InPPowCone,
  InDExpCone,
  InDPowCone,
  InPGeoMeanCone,
  InDGeoMeanCone,
  Positive,
  Negative,
  Unbounded,
  Zero,
  InPSDCone,
  InSVecPSDCone
};
std::ostream &operator<<(std::ostream &os, QConeKey val);
enum class ObjectiveSense { Undefined, Minimize, Maximize };
std::ostream &operator<<(std::ostream &os, ObjectiveSense val);
enum class SolutionStatus {
  Undefined,
  Unknown,
  Optimal,
  Feasible,
  Certificate,
  IllposedCert
};
std::ostream &operator<<(std::ostream &os, SolutionStatus val);
enum class AccSolutionStatus { Anything, Optimal, Feasible, Certificate };
std::ostream &operator<<(std::ostream &os, AccSolutionStatus val);
enum class ProblemStatus {
  Unknown,
  PrimalAndDualFeasible,
  PrimalFeasible,
  DualFeasible,
  PrimalInfeasible,
  DualInfeasible,
  PrimalAndDualInfeasible,
  IllPosed,
  PrimalInfeasibleOrUnbounded
};
std::ostream &operator<<(std::ostream &os, ProblemStatus val);
enum class SolverStatus { OK, Error, LostRace };
std::ostream &operator<<(std::ostream &os, SolverStatus val);
enum class SolutionType { Default, Basic, Interior, Integer };
std::ostream &operator<<(std::ostream &os, SolutionType val);
enum class StatusKey { Unknown, Basic, SuperBasic, OnBound, Infinity };
std::ostream &operator<<(std::ostream &os, StatusKey val);
enum class DJCDomainType {
  EqualTo,
  LessThan,
  GreaterThan,
  IsFree,
  InRange,
  InQCone,
  InRotatedQCone,
  InPExpCone,
  InPPowCone,
  InDExpCone,
  InDPowCone,
  InOneNormCone,
  InInfNormCone,
  InPGeoMeanCone,
  InDGeoMeanCone,
  InPSDCone
};
std::ostream &operator<<(std::ostream &os, DJCDomainType val);
} // namespace fusion
} // namespace mosek
namespace mosek {
namespace fusion {
class /*interface*/ ExprDomain : public virtual monty::RefCounted {
public:
  typedef monty::rc_ptr<ExprDomain> t;

  virtual void destroy() = 0;
  virtual ~ExprDomain(){};
  virtual monty::rc_ptr<::mosek::fusion::Term>
  __mosek_2fusion_2ExprDomain__toDJCTerm() = 0;
  monty::rc_ptr<::mosek::fusion::Term> toDJCTerm();
}; // interface ExprDomain;

class /*interface*/ Expression : public virtual monty::RefCounted {
public:
  typedef monty::rc_ptr<Expression> t;

  virtual void destroy() = 0;
  virtual ~Expression(){};
  virtual std::string toString() = 0;
  virtual monty::rc_ptr<::mosek::fusion::Expression>
  __mosek_2fusion_2Expression__pick(
      std::shared_ptr<monty::ndarray<int32_t, 2>> indexrows) = 0;
  monty::rc_ptr<::mosek::fusion::Expression>
  pick(std::shared_ptr<monty::ndarray<int32_t, 2>> indexrows);
  virtual monty::rc_ptr<::mosek::fusion::Expression>
  __mosek_2fusion_2Expression__pick(
      std::shared_ptr<monty::ndarray<int32_t, 1>> indexes) = 0;
  monty::rc_ptr<::mosek::fusion::Expression>
  pick(std::shared_ptr<monty::ndarray<int32_t, 1>> indexes);
  virtual monty::rc_ptr<::mosek::fusion::Expression>
  __mosek_2fusion_2Expression__index(
      std::shared_ptr<monty::ndarray<int32_t, 1>> indexes) = 0;
  monty::rc_ptr<::mosek::fusion::Expression>
  index(std::shared_ptr<monty::ndarray<int32_t, 1>> indexes);
  virtual monty::rc_ptr<::mosek::fusion::Expression>
  __mosek_2fusion_2Expression__index(int32_t i) = 0;
  monty::rc_ptr<::mosek::fusion::Expression> index(int32_t i);
  virtual monty::rc_ptr<::mosek::fusion::Expression>
  __mosek_2fusion_2Expression__slice(
      std::shared_ptr<monty::ndarray<int32_t, 1>> firsta,
      std::shared_ptr<monty::ndarray<int32_t, 1>> lasta) = 0;
  monty::rc_ptr<::mosek::fusion::Expression>
  slice(std::shared_ptr<monty::ndarray<int32_t, 1>> firsta,
        std::shared_ptr<monty::ndarray<int32_t, 1>> lasta);
  virtual monty::rc_ptr<::mosek::fusion::Expression>
  __mosek_2fusion_2Expression__slice(int32_t first, int32_t last) = 0;
  monty::rc_ptr<::mosek::fusion::Expression> slice(int32_t first, int32_t last);
  virtual void eval(monty::rc_ptr<::mosek::fusion::WorkStack> rs,
                    monty::rc_ptr<::mosek::fusion::WorkStack> ws,
                    monty::rc_ptr<::mosek::fusion::WorkStack> xs) = 0;
  virtual int32_t getND() = 0;
  virtual int32_t getDim(int32_t d) = 0;
  virtual int64_t getSize() = 0;
  virtual std::shared_ptr<monty::ndarray<int32_t, 1>> getShape() = 0;
}; // interface Expression;

class /*interface*/ Variable : public virtual ::mosek::fusion::Expression {
public:
  typedef monty::rc_ptr<Variable> t;

  virtual void destroy() = 0;
  virtual ~Variable(){};
  virtual std::string toString() = 0;
  virtual int32_t numInst() = 0;
  virtual int32_t
  inst(int32_t spoffset, std::shared_ptr<monty::ndarray<int64_t, 1>> sparsity,
       int32_t nioffset,
       std::shared_ptr<monty::ndarray<int64_t, 1>> nativeidxs) = 0;
  virtual void inst(int32_t offset,
                    std::shared_ptr<monty::ndarray<int64_t, 1>> nativeidxs) = 0;
  virtual void remove() = 0;
  virtual int32_t getND() = 0;
  virtual std::shared_ptr<monty::ndarray<int32_t, 1>> getShape() = 0;
  virtual monty::rc_ptr<::mosek::fusion::Model>
  __mosek_2fusion_2Variable__getModel() = 0;
  monty::rc_ptr<::mosek::fusion::Model> getModel();
  virtual int64_t getSize() = 0;
  virtual void setLevel(std::shared_ptr<monty::ndarray<double, 1>> v) = 0;
  virtual monty::rc_ptr<::mosek::fusion::Variable>
  __mosek_2fusion_2Variable__fromTril(int32_t dim) = 0;
  monty::rc_ptr<::mosek::fusion::Variable> fromTril(int32_t dim);
  virtual monty::rc_ptr<::mosek::fusion::Variable>
  __mosek_2fusion_2Variable__tril() = 0;
  monty::rc_ptr<::mosek::fusion::Variable> tril();
  virtual monty::rc_ptr<::mosek::fusion::Variable>
  __mosek_2fusion_2Variable__reshape(int32_t dim0, int32_t dim1,
                                     int32_t dim2) = 0;
  monty::rc_ptr<::mosek::fusion::Variable> reshape(int32_t dim0, int32_t dim1,
                                                   int32_t dim2);
  virtual monty::rc_ptr<::mosek::fusion::Variable>
  __mosek_2fusion_2Variable__reshape(int32_t dim0, int32_t dim1) = 0;
  monty::rc_ptr<::mosek::fusion::Variable> reshape(int32_t dim0, int32_t dim1);
  virtual monty::rc_ptr<::mosek::fusion::Variable>
  __mosek_2fusion_2Variable__reshape(int32_t dim0) = 0;
  monty::rc_ptr<::mosek::fusion::Variable> reshape(int32_t dim0);
  virtual monty::rc_ptr<::mosek::fusion::Variable>
  __mosek_2fusion_2Variable__reshape(
      std::shared_ptr<monty::ndarray<int32_t, 1>> shape) = 0;
  monty::rc_ptr<::mosek::fusion::Variable>
  reshape(std::shared_ptr<monty::ndarray<int32_t, 1>> shape);
  virtual void set_values(std::shared_ptr<monty::ndarray<double, 1>> value,
                          bool primal) = 0;
  virtual std::shared_ptr<monty::ndarray<double, 1>> dual() = 0;
  virtual std::shared_ptr<monty::ndarray<double, 1>> level() = 0;
  virtual void values(int32_t offset,
                      std::shared_ptr<monty::ndarray<double, 1>> target,
                      bool primal) = 0;
  virtual void make_continuous() = 0;
  virtual void make_integer() = 0;
  virtual void makeContinuous() = 0;
  virtual void makeInteger() = 0;
  virtual monty::rc_ptr<::mosek::fusion::Variable>
  __mosek_2fusion_2Variable__transpose() = 0;
  monty::rc_ptr<::mosek::fusion::Variable> transpose();
  virtual monty::rc_ptr<::mosek::fusion::Variable>
  __mosek_2fusion_2Variable__pick(
      std::shared_ptr<monty::ndarray<int32_t, 1>> i1,
      std::shared_ptr<monty::ndarray<int32_t, 1>> i2,
      std::shared_ptr<monty::ndarray<int32_t, 1>> i3) = 0;
  monty::rc_ptr<::mosek::fusion::Variable>
  pick(std::shared_ptr<monty::ndarray<int32_t, 1>> i1,
       std::shared_ptr<monty::ndarray<int32_t, 1>> i2,
       std::shared_ptr<monty::ndarray<int32_t, 1>> i3);
  virtual monty::rc_ptr<::mosek::fusion::Variable>
  __mosek_2fusion_2Variable__pick(
      std::shared_ptr<monty::ndarray<int32_t, 1>> i1,
      std::shared_ptr<monty::ndarray<int32_t, 1>> i2) = 0;
  monty::rc_ptr<::mosek::fusion::Variable>
  pick(std::shared_ptr<monty::ndarray<int32_t, 1>> i1,
       std::shared_ptr<monty::ndarray<int32_t, 1>> i2);
  virtual monty::rc_ptr<::mosek::fusion::Variable>
  __mosek_2fusion_2Variable__pick(
      std::shared_ptr<monty::ndarray<int32_t, 2>> midxs) = 0;
  monty::rc_ptr<::mosek::fusion::Variable>
  pick(std::shared_ptr<monty::ndarray<int32_t, 2>> midxs);
  virtual monty::rc_ptr<::mosek::fusion::Expression>
  __mosek_2fusion_2Expression__pick(
      std::shared_ptr<monty::ndarray<int32_t, 2>> indexrows);
  virtual monty::rc_ptr<::mosek::fusion::Variable>
  __mosek_2fusion_2Variable__pick(
      std::shared_ptr<monty::ndarray<int32_t, 1>> idxs) = 0;
  monty::rc_ptr<::mosek::fusion::Variable>
  pick(std::shared_ptr<monty::ndarray<int32_t, 1>> idxs);
  virtual monty::rc_ptr<::mosek::fusion::Expression>
  __mosek_2fusion_2Expression__pick(
      std::shared_ptr<monty::ndarray<int32_t, 1>> indexes);
  virtual monty::rc_ptr<::mosek::fusion::Variable>
  __mosek_2fusion_2Variable__antidiag() = 0;
  monty::rc_ptr<::mosek::fusion::Variable> antidiag();
  virtual monty::rc_ptr<::mosek::fusion::Variable>
  __mosek_2fusion_2Variable__antidiag(int32_t index) = 0;
  monty::rc_ptr<::mosek::fusion::Variable> antidiag(int32_t index);
  virtual monty::rc_ptr<::mosek::fusion::Variable>
  __mosek_2fusion_2Variable__diag() = 0;
  monty::rc_ptr<::mosek::fusion::Variable> diag();
  virtual monty::rc_ptr<::mosek::fusion::Variable>
  __mosek_2fusion_2Variable__diag(int32_t index) = 0;
  monty::rc_ptr<::mosek::fusion::Variable> diag(int32_t index);
  virtual monty::rc_ptr<::mosek::fusion::Variable>
  __mosek_2fusion_2Variable__index(
      std::shared_ptr<monty::ndarray<int32_t, 1>> idx) = 0;
  monty::rc_ptr<::mosek::fusion::Variable>
  index(std::shared_ptr<monty::ndarray<int32_t, 1>> idx);
  virtual monty::rc_ptr<::mosek::fusion::Expression>
  __mosek_2fusion_2Expression__index(
      std::shared_ptr<monty::ndarray<int32_t, 1>> indexes);
  virtual monty::rc_ptr<::mosek::fusion::Variable>
  __mosek_2fusion_2Variable__index(int32_t i1, int32_t i2, int32_t i3) = 0;
  monty::rc_ptr<::mosek::fusion::Variable> index(int32_t i1, int32_t i2,
                                                 int32_t i3);
  virtual monty::rc_ptr<::mosek::fusion::Variable>
  __mosek_2fusion_2Variable__index(int32_t i1, int32_t i2) = 0;
  monty::rc_ptr<::mosek::fusion::Variable> index(int32_t i1, int32_t i2);
  virtual monty::rc_ptr<::mosek::fusion::Variable>
  __mosek_2fusion_2Variable__index(int32_t i1) = 0;
  monty::rc_ptr<::mosek::fusion::Variable> index(int32_t i1);
  virtual monty::rc_ptr<::mosek::fusion::Expression>
  __mosek_2fusion_2Expression__index(int32_t i);
  virtual monty::rc_ptr<::mosek::fusion::Variable>
  __mosek_2fusion_2Variable__slice(
      std::shared_ptr<monty::ndarray<int32_t, 1>> firsta,
      std::shared_ptr<monty::ndarray<int32_t, 1>> lasta) = 0;
  monty::rc_ptr<::mosek::fusion::Variable>
  slice(std::shared_ptr<monty::ndarray<int32_t, 1>> firsta,
        std::shared_ptr<monty::ndarray<int32_t, 1>> lasta);
  virtual monty::rc_ptr<::mosek::fusion::Expression>
  __mosek_2fusion_2Expression__slice(
      std::shared_ptr<monty::ndarray<int32_t, 1>> firsta,
      std::shared_ptr<monty::ndarray<int32_t, 1>> lasta);
  virtual monty::rc_ptr<::mosek::fusion::Variable>
  __mosek_2fusion_2Variable__slice(int32_t first, int32_t last) = 0;
  monty::rc_ptr<::mosek::fusion::Variable> slice(int32_t first, int32_t last);
  virtual monty::rc_ptr<::mosek::fusion::Expression>
  __mosek_2fusion_2Expression__slice(int32_t first, int32_t last);
  virtual monty::rc_ptr<::mosek::fusion::Expression>
  __mosek_2fusion_2Variable__asExpr() = 0;
  monty::rc_ptr<::mosek::fusion::Expression> asExpr();
}; // interface Variable;

class /*interface*/ Parameter : public virtual ::mosek::fusion::Expression {
public:
  typedef monty::rc_ptr<Parameter> t;

  virtual void destroy() = 0;
  virtual ~Parameter(){};
  virtual monty::rc_ptr<::mosek::fusion::Parameter>
  __mosek_2fusion_2Parameter__clone(
      monty::rc_ptr<::mosek::fusion::Model> m) = 0;
  monty::rc_ptr<::mosek::fusion::Parameter>
  clone(monty::rc_ptr<::mosek::fusion::Model> m);
  virtual monty::rc_ptr<::mosek::fusion::Model>
  __mosek_2fusion_2Parameter__getModel() = 0;
  monty::rc_ptr<::mosek::fusion::Model> getModel();
  virtual int64_t getSize() = 0;
  virtual void getAllIndexes(std::shared_ptr<monty::ndarray<int32_t, 1>> dst,
                             int32_t ofs) = 0;
  virtual int32_t getIndex(int32_t i) = 0;
  virtual void getSp(std::shared_ptr<monty::ndarray<int64_t, 1>> dest,
                     int32_t offset) = 0;
  virtual std::shared_ptr<monty::ndarray<int32_t, 1>> getShape() = 0;
  virtual int32_t getND() = 0;
  virtual int32_t getDim(int32_t i) = 0;
  virtual monty::rc_ptr<::mosek::fusion::Parameter>
  __mosek_2fusion_2Parameter__pick(
      std::shared_ptr<monty::ndarray<int32_t, 2>> midxs) = 0;
  monty::rc_ptr<::mosek::fusion::Parameter>
  pick(std::shared_ptr<monty::ndarray<int32_t, 2>> midxs);
  virtual monty::rc_ptr<::mosek::fusion::Expression>
  __mosek_2fusion_2Expression__pick(
      std::shared_ptr<monty::ndarray<int32_t, 2>> indexrows);
  virtual monty::rc_ptr<::mosek::fusion::Parameter>
  __mosek_2fusion_2Parameter__pick(
      std::shared_ptr<monty::ndarray<int32_t, 1>> idxs) = 0;
  monty::rc_ptr<::mosek::fusion::Parameter>
  pick(std::shared_ptr<monty::ndarray<int32_t, 1>> idxs);
  virtual monty::rc_ptr<::mosek::fusion::Expression>
  __mosek_2fusion_2Expression__pick(
      std::shared_ptr<monty::ndarray<int32_t, 1>> indexes);
  virtual monty::rc_ptr<::mosek::fusion::Parameter>
  __mosek_2fusion_2Parameter__slice(
      std::shared_ptr<monty::ndarray<int32_t, 1>> astart,
      std::shared_ptr<monty::ndarray<int32_t, 1>> astop) = 0;
  monty::rc_ptr<::mosek::fusion::Parameter>
  slice(std::shared_ptr<monty::ndarray<int32_t, 1>> astart,
        std::shared_ptr<monty::ndarray<int32_t, 1>> astop);
  virtual monty::rc_ptr<::mosek::fusion::Expression>
  __mosek_2fusion_2Expression__slice(
      std::shared_ptr<monty::ndarray<int32_t, 1>> firsta,
      std::shared_ptr<monty::ndarray<int32_t, 1>> lasta);
  virtual monty::rc_ptr<::mosek::fusion::Parameter>
  __mosek_2fusion_2Parameter__slice(int32_t start, int32_t stop) = 0;
  monty::rc_ptr<::mosek::fusion::Parameter> slice(int32_t start, int32_t stop);
  virtual monty::rc_ptr<::mosek::fusion::Expression>
  __mosek_2fusion_2Expression__slice(int32_t first, int32_t last);
  virtual bool isSparse() = 0;
  virtual monty::rc_ptr<::mosek::fusion::Parameter>
  __mosek_2fusion_2Parameter__reshape(
      std::shared_ptr<monty::ndarray<int32_t, 1>> dims) = 0;
  monty::rc_ptr<::mosek::fusion::Parameter>
  reshape(std::shared_ptr<monty::ndarray<int32_t, 1>> dims);
  virtual int32_t getNumNonzero() = 0;
  virtual monty::rc_ptr<::mosek::fusion::Expression>
  __mosek_2fusion_2Parameter__asExpr() = 0;
  monty::rc_ptr<::mosek::fusion::Expression> asExpr();
  virtual std::shared_ptr<monty::ndarray<double, 1>> getValue() = 0;
  virtual void setValue(std::shared_ptr<monty::ndarray<double, 2>> values2) = 0;
  virtual void setValue(std::shared_ptr<monty::ndarray<double, 1>> values) = 0;
  virtual void setValue(double value) = 0;
}; // interface Parameter;

class FusionException : public ::monty::Exception {
private:
  std::string msg;

protected:
public:
  typedef monty::rc_ptr<FusionException> t;

  FusionException(const std::string &msg_);
  virtual /* override */ std::string toString();
}; // class FusionException;

class SolutionError : public ::mosek::fusion::FusionException {
private:
protected:
public:
  typedef monty::rc_ptr<SolutionError> t;

  SolutionError();
  SolutionError(const std::string &msg);
}; // class SolutionError;

class UnimplementedError : public ::monty::RuntimeException {
private:
protected:
public:
  typedef monty::rc_ptr<UnimplementedError> t;

  UnimplementedError(const std::string &msg);
}; // class UnimplementedError;

class FatalError : public ::monty::RuntimeException {
private:
protected:
public:
  typedef monty::rc_ptr<FatalError> t;

  FatalError(const std::string &msg);
}; // class FatalError;

class UnexpectedError : public ::monty::RuntimeException {
private:
protected:
public:
  typedef monty::rc_ptr<UnexpectedError> t;

  UnexpectedError(::mosek::fusion::FusionException e);
  UnexpectedError(const std::string &msg);
}; // class UnexpectedError;

class FusionRuntimeException : public ::monty::RuntimeException {
private:
  std::string msg;

protected:
public:
  typedef monty::rc_ptr<FusionRuntimeException> t;

  FusionRuntimeException(const std::string &msg_);
  virtual /* override */ std::string toString();
}; // class FusionRuntimeException;

class SparseFormatError : public ::mosek::fusion::FusionRuntimeException {
private:
protected:
public:
  typedef monty::rc_ptr<SparseFormatError> t;

  SparseFormatError(const std::string &msg);
}; // class SparseFormatError;

class SliceError : public ::mosek::fusion::FusionRuntimeException {
private:
protected:
public:
  typedef monty::rc_ptr<SliceError> t;

  SliceError();
  SliceError(const std::string &msg);
}; // class SliceError;

class UpdateError : public ::mosek::fusion::FusionRuntimeException {
private:
protected:
public:
  typedef monty::rc_ptr<UpdateError> t;

  UpdateError();
  UpdateError(const std::string &msg);
}; // class UpdateError;

class SetDefinitionError : public ::mosek::fusion::FusionRuntimeException {
private:
protected:
public:
  typedef monty::rc_ptr<SetDefinitionError> t;

  SetDefinitionError(const std::string &msg);
}; // class SetDefinitionError;

class OptimizeError : public ::mosek::fusion::FusionRuntimeException {
private:
protected:
public:
  typedef monty::rc_ptr<OptimizeError> t;

  OptimizeError(const std::string &msg);
}; // class OptimizeError;

class NameError : public ::mosek::fusion::FusionRuntimeException {
private:
protected:
public:
  typedef monty::rc_ptr<NameError> t;

  NameError(const std::string &msg);
}; // class NameError;

class DeletionError : public ::mosek::fusion::FusionRuntimeException {
private:
protected:
public:
  typedef monty::rc_ptr<DeletionError> t;

  DeletionError(const std::string &msg);
}; // class DeletionError;

class ModelError : public ::mosek::fusion::FusionRuntimeException {
private:
protected:
public:
  typedef monty::rc_ptr<ModelError> t;

  ModelError(const std::string &msg);
}; // class ModelError;

class MatrixError : public ::mosek::fusion::FusionRuntimeException {
private:
protected:
public:
  typedef monty::rc_ptr<MatrixError> t;

  MatrixError(const std::string &msg);
}; // class MatrixError;

class DimensionError : public ::mosek::fusion::FusionRuntimeException {
private:
protected:
public:
  typedef monty::rc_ptr<DimensionError> t;

  DimensionError(const std::string &msg);
}; // class DimensionError;

class LengthError : public ::mosek::fusion::FusionRuntimeException {
private:
protected:
public:
  typedef monty::rc_ptr<LengthError> t;

  LengthError(const std::string &msg);
}; // class LengthError;

class RangeError : public ::mosek::fusion::FusionRuntimeException {
private:
protected:
public:
  typedef monty::rc_ptr<RangeError> t;

  RangeError(const std::string &msg);
}; // class RangeError;

class IndexError : public ::mosek::fusion::FusionRuntimeException {
private:
protected:
public:
  typedef monty::rc_ptr<IndexError> t;

  IndexError(const std::string &msg);
}; // class IndexError;

class DomainError : public ::mosek::fusion::FusionRuntimeException {
private:
protected:
public:
  typedef monty::rc_ptr<DomainError> t;

  DomainError(const std::string &msg);
}; // class DomainError;

class ValueConversionError : public ::mosek::fusion::FusionRuntimeException {
private:
protected:
public:
  typedef monty::rc_ptr<ValueConversionError> t;

  ValueConversionError(const std::string &msg);
}; // class ValueConversionError;

class ParameterError : public ::mosek::fusion::FusionRuntimeException {
private:
protected:
public:
  typedef monty::rc_ptr<ParameterError> t;

  ParameterError(const std::string &msg);
}; // class ParameterError;

class ExpressionError : public ::mosek::fusion::FusionRuntimeException {
private:
protected:
public:
  typedef monty::rc_ptr<ExpressionError> t;

  ExpressionError(const std::string &msg);
}; // class ExpressionError;

class IOError : public ::mosek::fusion::FusionRuntimeException {
private:
protected:
public:
  typedef monty::rc_ptr<IOError> t;

  IOError(const std::string &msg);
}; // class IOError;

class Disjunction : public virtual monty::RefCounted {
public:
  p_Disjunction *_impl;
  Disjunction(int64_t id);

protected:
  Disjunction(p_Disjunction *_impl);

public:
  Disjunction(const Disjunction &) = delete;
  const Disjunction &operator=(const Disjunction &) = delete;
  friend class p_Disjunction;
  virtual ~Disjunction();
  virtual void destroy();
  typedef monty::rc_ptr<Disjunction> t;

}; // class Disjunction;

class DisjunctionTerms : public virtual monty::RefCounted {
public:
  p_DisjunctionTerms *_impl;

protected:
  DisjunctionTerms(p_DisjunctionTerms *_impl);

public:
  DisjunctionTerms(const DisjunctionTerms &) = delete;
  const DisjunctionTerms &operator=(const DisjunctionTerms &) = delete;
  friend class p_DisjunctionTerms;
  virtual ~DisjunctionTerms();
  virtual void destroy();
  typedef monty::rc_ptr<DisjunctionTerms> t;

  DisjunctionTerms(monty::rc_ptr<::mosek::fusion::DisjunctionTerms> terms1,
                   monty::rc_ptr<::mosek::fusion::ExprRangeDomain> term);
  DisjunctionTerms(monty::rc_ptr<::mosek::fusion::DisjunctionTerms> terms1,
                   monty::rc_ptr<::mosek::fusion::ExprLinearDomain> term);
  DisjunctionTerms(monty::rc_ptr<::mosek::fusion::DisjunctionTerms> terms1,
                   monty::rc_ptr<::mosek::fusion::DisjunctionTerms> term2);
  DisjunctionTerms(
      monty::rc_ptr<::mosek::fusion::DisjunctionTerms> term1,
      std::shared_ptr<monty::ndarray<monty::rc_ptr<::mosek::fusion::Term>, 1>>
          term2);
  DisjunctionTerms(monty::rc_ptr<::mosek::fusion::DisjunctionTerms> term1,
                   monty::rc_ptr<::mosek::fusion::Term> term2);
  DisjunctionTerms(
      std::shared_ptr<monty::ndarray<monty::rc_ptr<::mosek::fusion::Term>, 1>>
          terms);
  DisjunctionTerms(monty::rc_ptr<::mosek::fusion::ExprRangeDomain> term);
  DisjunctionTerms(monty::rc_ptr<::mosek::fusion::ExprLinearDomain> term);
}; // class DisjunctionTerms;

class Term : public virtual ::mosek::fusion::ExprDomain {
public:
  p_Term *_impl;
  Term(std::shared_ptr<
           monty::ndarray<monty::rc_ptr<::mosek::fusion::Expression>, 1>>
           elist,
       std::shared_ptr<
           monty::ndarray<monty::rc_ptr<::mosek::fusion::DJCDomain>, 1>>
           dlist);
  Term(monty::rc_ptr<::mosek::fusion::Expression> e,
       monty::rc_ptr<::mosek::fusion::DJCDomain> d);
  Term(std::shared_ptr<
       monty::ndarray<monty::rc_ptr<::mosek::fusion::SimpleTerm>, 1>>
           t);

protected:
  Term(p_Term *_impl);

public:
  Term(const Term &) = delete;
  const Term &operator=(const Term &) = delete;
  friend class p_Term;
  virtual ~Term();
  virtual void destroy();
  typedef monty::rc_ptr<Term> t;

  virtual int32_t size();
  virtual monty::rc_ptr<::mosek::fusion::Term>
  __mosek_2fusion_2Term__toDJCTerm();
  monty::rc_ptr<::mosek::fusion::Term> toDJCTerm();
  /* override: mosek.fusion.ExprDomain.toDJCTerm*/
  virtual monty::rc_ptr<::mosek::fusion::Term>
  __mosek_2fusion_2ExprDomain__toDJCTerm();
}; // class Term;

class SimpleTerm : public ::mosek::fusion::Term {
  SimpleTerm(monty::rc_ptr<::mosek::fusion::Expression> e,
             monty::rc_ptr<::mosek::fusion::DJCDomain> d);

protected:
  SimpleTerm(p_SimpleTerm *_impl);

public:
  SimpleTerm(const SimpleTerm &) = delete;
  const SimpleTerm &operator=(const SimpleTerm &) = delete;
  friend class p_SimpleTerm;
  virtual ~SimpleTerm();
  virtual void destroy();
  typedef monty::rc_ptr<SimpleTerm> t;

}; // class SimpleTerm;

class DJCDomain : public virtual monty::RefCounted {
public:
  p_DJCDomain *_impl;
  DJCDomain(std::shared_ptr<monty::ndarray<double, 1>> b_,
            std::shared_ptr<monty::ndarray<double, 1>> par_,
            std::shared_ptr<monty::ndarray<int32_t, 1>> shape_,
            mosek::fusion::DJCDomainType dom_);
  DJCDomain(std::shared_ptr<monty::ndarray<double, 1>> b_,
            std::shared_ptr<monty::ndarray<double, 1>> par_,
            std::shared_ptr<monty::ndarray<int32_t, 1>> shape_,
            int32_t conedim_, mosek::fusion::DJCDomainType dom_);

protected:
  DJCDomain(p_DJCDomain *_impl);

public:
  DJCDomain(const DJCDomain &) = delete;
  const DJCDomain &operator=(const DJCDomain &) = delete;
  friend class p_DJCDomain;
  virtual ~DJCDomain();
  virtual void destroy();
  typedef monty::rc_ptr<DJCDomain> t;
  mosek::fusion::DJCDomainType get_dom();
  void set_dom(mosek::fusion::DJCDomainType val);
  int32_t get_conedim();
  void set_conedim(int32_t val);
  std::shared_ptr<monty::ndarray<int32_t, 1>> get_shape();
  void set_shape(std::shared_ptr<monty::ndarray<int32_t, 1>> val);
  std::shared_ptr<monty::ndarray<double, 1>> get_par();
  void set_par(std::shared_ptr<monty::ndarray<double, 1>> val);
  std::shared_ptr<monty::ndarray<double, 1>> get_b();
  void set_b(std::shared_ptr<monty::ndarray<double, 1>> val);

  virtual int32_t size();
}; // class DJCDomain;

class DJC : public virtual monty::RefCounted {
public:
  p_DJC *_impl;

protected:
  DJC(p_DJC *_impl);

public:
  DJC(const DJC &) = delete;
  const DJC &operator=(const DJC &) = delete;
  friend class p_DJC;
  virtual ~DJC();
  virtual void destroy();
  typedef monty::rc_ptr<DJC> t;

  static monty::rc_ptr<::mosek::fusion::Term> ANDFromTerms(
      std::shared_ptr<monty::ndarray<monty::rc_ptr<::mosek::fusion::Term>, 1>>
          tlist);
  static monty::rc_ptr<::mosek::fusion::Term>
  AND(std::shared_ptr<
      monty::ndarray<monty::rc_ptr<::mosek::fusion::ExprDomain>, 1>>
          elist);
  static monty::rc_ptr<::mosek::fusion::Term>
  AND(monty::rc_ptr<::mosek::fusion::SimpleTerm> s1,
      monty::rc_ptr<::mosek::fusion::SimpleTerm> s2,
      monty::rc_ptr<::mosek::fusion::SimpleTerm> s3);
  static monty::rc_ptr<::mosek::fusion::Term>
  AND(monty::rc_ptr<::mosek::fusion::SimpleTerm> s1,
      monty::rc_ptr<::mosek::fusion::SimpleTerm> s2);
  static monty::rc_ptr<::mosek::fusion::Term>
  AND(monty::rc_ptr<::mosek::fusion::SimpleTerm> s1);
  static monty::rc_ptr<::mosek::fusion::Term>
  AND(std::shared_ptr<
      monty::ndarray<monty::rc_ptr<::mosek::fusion::SimpleTerm>, 1>>
          slist);
  static monty::rc_ptr<::mosek::fusion::SimpleTerm>
  term(monty::rc_ptr<::mosek::fusion::Expression> expr,
       monty::rc_ptr<::mosek::fusion::RangeDomain> dom);
  static monty::rc_ptr<::mosek::fusion::SimpleTerm>
  term(monty::rc_ptr<::mosek::fusion::Variable> x,
       monty::rc_ptr<::mosek::fusion::RangeDomain> dom);
  static monty::rc_ptr<::mosek::fusion::SimpleTerm>
  term(monty::rc_ptr<::mosek::fusion::ExprRangeDomain> exprdom);
  static monty::rc_ptr<::mosek::fusion::SimpleTerm>
  term(monty::rc_ptr<::mosek::fusion::Expression> expr,
       monty::rc_ptr<::mosek::fusion::LinearDomain> dom);
  static monty::rc_ptr<::mosek::fusion::SimpleTerm>
  term(monty::rc_ptr<::mosek::fusion::ExprLinearDomain> exprdom);
  static monty::rc_ptr<::mosek::fusion::SimpleTerm>
  term(monty::rc_ptr<::mosek::fusion::Variable> x,
       monty::rc_ptr<::mosek::fusion::LinearDomain> dom);
}; // class DJC;

// mosek.fusion.BaseModel from file 'src/fusion/cxx/BaseModel.h'

class BaseModel : public monty::RefCounted {
private:
  p_BaseModel *_impl;

protected:
  BaseModel(p_BaseModel *ptr);

public:
  friend class p_BaseModel;

  virtual void dispose();

  ~BaseModel();
};
// End of file 'src/fusion/cxx/BaseModel.h'
class Model : public ::mosek::fusion::BaseModel {
  Model(monty::rc_ptr<::mosek::fusion::Model> m);

protected:
  Model(p_Model *_impl);

public:
  Model(const Model &) = delete;
  const Model &operator=(const Model &) = delete;
  friend class p_Model;
  virtual ~Model();
  virtual void destroy();
  typedef monty::rc_ptr<Model> t;

  Model(const std::string &name, int32_t basesize);
  Model(int32_t basesize);
  Model(const std::string &name);
  Model();
  virtual monty::rc_ptr<::mosek::fusion::Disjunction>
  __mosek_2fusion_2Model__disjunction(
      const std::string &name,
      std::shared_ptr<monty::ndarray<monty::rc_ptr<::mosek::fusion::Term>, 1>>
          terms);
  monty::rc_ptr<::mosek::fusion::Disjunction> disjunction(
      const std::string &name,
      std::shared_ptr<monty::ndarray<monty::rc_ptr<::mosek::fusion::Term>, 1>>
          terms);
  virtual monty::rc_ptr<::mosek::fusion::Disjunction>
  __mosek_2fusion_2Model__disjunction(
      std::shared_ptr<monty::ndarray<monty::rc_ptr<::mosek::fusion::Term>, 1>>
          terms);
  monty::rc_ptr<::mosek::fusion::Disjunction> disjunction(
      std::shared_ptr<monty::ndarray<monty::rc_ptr<::mosek::fusion::Term>, 1>>
          terms);
  virtual monty::rc_ptr<::mosek::fusion::Disjunction>
  __mosek_2fusion_2Model__disjunction(
      monty::rc_ptr<::mosek::fusion::DisjunctionTerms> djcterms);
  monty::rc_ptr<::mosek::fusion::Disjunction>
  disjunction(monty::rc_ptr<::mosek::fusion::DisjunctionTerms> djcterms);
  virtual monty::rc_ptr<::mosek::fusion::Disjunction>
  __mosek_2fusion_2Model__disjunction(
      const std::string &name,
      monty::rc_ptr<::mosek::fusion::DisjunctionTerms> djcterms);
  monty::rc_ptr<::mosek::fusion::Disjunction>
  disjunction(const std::string &name,
              monty::rc_ptr<::mosek::fusion::DisjunctionTerms> djcterms);
  virtual monty::rc_ptr<::mosek::fusion::Disjunction>
  __mosek_2fusion_2Model__disjunction(
      const std::string &name,
      std::shared_ptr<
          monty::ndarray<monty::rc_ptr<::mosek::fusion::ExprDomain>, 1>>
          terms);
  monty::rc_ptr<::mosek::fusion::Disjunction>
  disjunction(const std::string &name,
              std::shared_ptr<
                  monty::ndarray<monty::rc_ptr<::mosek::fusion::ExprDomain>, 1>>
                  terms);
  virtual monty::rc_ptr<::mosek::fusion::Disjunction>
  __mosek_2fusion_2Model__disjunction(
      std::shared_ptr<
          monty::ndarray<monty::rc_ptr<::mosek::fusion::ExprDomain>, 1>>
          terms);
  monty::rc_ptr<::mosek::fusion::Disjunction>
  disjunction(std::shared_ptr<
              monty::ndarray<monty::rc_ptr<::mosek::fusion::ExprDomain>, 1>>
                  terms);
  virtual monty::rc_ptr<::mosek::fusion::Disjunction>
  __mosek_2fusion_2Model__disjunction(
      const std::string &name, monty::rc_ptr<::mosek::fusion::ExprDomain> term);
  monty::rc_ptr<::mosek::fusion::Disjunction>
  disjunction(const std::string &name,
              monty::rc_ptr<::mosek::fusion::ExprDomain> term);
  virtual monty::rc_ptr<::mosek::fusion::Disjunction>
  __mosek_2fusion_2Model__disjunction(
      monty::rc_ptr<::mosek::fusion::ExprDomain> term);
  monty::rc_ptr<::mosek::fusion::Disjunction>
  disjunction(monty::rc_ptr<::mosek::fusion::ExprDomain> term);
  virtual monty::rc_ptr<::mosek::fusion::Disjunction>
  __mosek_2fusion_2Model__disjunction(monty::rc_ptr<::mosek::fusion::Term> t1,
                                      monty::rc_ptr<::mosek::fusion::Term> t2,
                                      monty::rc_ptr<::mosek::fusion::Term> t3);
  monty::rc_ptr<::mosek::fusion::Disjunction>
  disjunction(monty::rc_ptr<::mosek::fusion::Term> t1,
              monty::rc_ptr<::mosek::fusion::Term> t2,
              monty::rc_ptr<::mosek::fusion::Term> t3);
  virtual monty::rc_ptr<::mosek::fusion::Disjunction>
  __mosek_2fusion_2Model__disjunction(monty::rc_ptr<::mosek::fusion::Term> t1,
                                      monty::rc_ptr<::mosek::fusion::Term> t2);
  monty::rc_ptr<::mosek::fusion::Disjunction>
  disjunction(monty::rc_ptr<::mosek::fusion::Term> t1,
              monty::rc_ptr<::mosek::fusion::Term> t2);
  virtual monty::rc_ptr<::mosek::fusion::Disjunction>
  __mosek_2fusion_2Model__disjunction(monty::rc_ptr<::mosek::fusion::Term> t1);
  monty::rc_ptr<::mosek::fusion::Disjunction>
  disjunction(monty::rc_ptr<::mosek::fusion::Term> t1);
  virtual monty::rc_ptr<::mosek::fusion::Disjunction>
  __mosek_2fusion_2Model__disjunction(const std::string &name,
                                      monty::rc_ptr<::mosek::fusion::Term> t1,
                                      monty::rc_ptr<::mosek::fusion::Term> t2,
                                      monty::rc_ptr<::mosek::fusion::Term> t3);
  monty::rc_ptr<::mosek::fusion::Disjunction>
  disjunction(const std::string &name, monty::rc_ptr<::mosek::fusion::Term> t1,
              monty::rc_ptr<::mosek::fusion::Term> t2,
              monty::rc_ptr<::mosek::fusion::Term> t3);
  virtual monty::rc_ptr<::mosek::fusion::Disjunction>
  __mosek_2fusion_2Model__disjunction(const std::string &name,
                                      monty::rc_ptr<::mosek::fusion::Term> t1,
                                      monty::rc_ptr<::mosek::fusion::Term> t2);
  monty::rc_ptr<::mosek::fusion::Disjunction>
  disjunction(const std::string &name, monty::rc_ptr<::mosek::fusion::Term> t1,
              monty::rc_ptr<::mosek::fusion::Term> t2);
  virtual monty::rc_ptr<::mosek::fusion::Disjunction>
  __mosek_2fusion_2Model__disjunction(const std::string &name,
                                      monty::rc_ptr<::mosek::fusion::Term> t1);
  monty::rc_ptr<::mosek::fusion::Disjunction>
  disjunction(const std::string &name, monty::rc_ptr<::mosek::fusion::Term> t1);
  static void putlicensewait(bool wait);
  static void putlicensepath(const std::string &licfile);
  static void putlicensecode(std::shared_ptr<monty::ndarray<int32_t, 1>> code);
  virtual /* override */ void dispose();
  virtual MSKtask_t __mosek_2fusion_2Model__getTask();
  MSKtask_t getTask();
  virtual void getConstraintDuals(
      bool lower, std::shared_ptr<monty::ndarray<int32_t, 1>> nativeidxs,
      std::shared_ptr<monty::ndarray<double, 1>> res, int32_t offset);
  virtual void getConstraintValues(
      bool primal, std::shared_ptr<monty::ndarray<int32_t, 1>> nativeidxs,
      std::shared_ptr<monty::ndarray<double, 1>> res, int32_t offset);
  virtual void getVariableDuals(
      bool lower, std::shared_ptr<monty::ndarray<int64_t, 1>> nativeidxs,
      std::shared_ptr<monty::ndarray<double, 1>> res, int32_t offset);
  virtual void getVariableValues(
      bool primal, std::shared_ptr<monty::ndarray<int64_t, 1>> nativeidxs,
      std::shared_ptr<monty::ndarray<double, 1>> res, int32_t offset);
  virtual void
  setVariableValues(bool primal,
                    std::shared_ptr<monty::ndarray<int64_t, 1>> nativeidxs,
                    std::shared_ptr<monty::ndarray<double, 1>> values);
  virtual void flushNames();
  virtual void writeTaskNoFlush(const std::string &filename);
  virtual void writeTaskStream(const std::string &ext, std::ostream &stream);
  virtual void dataReport();
  virtual void dataReport(int32_t detail);
  virtual void writeTask(const std::string &filename);
  virtual int64_t getSolverLIntInfo(const std::string &name);
  virtual int32_t getSolverIntInfo(const std::string &name);
  virtual double getSolverDoubleInfo(const std::string &name);
  virtual void setCallbackHandler(mosek::cbhandler_t h);
  virtual void setDataCallbackHandler(mosek::datacbhandler_t h);
  virtual void setLogHandler(mosek::msghandler_t h);
  virtual void setSolverParam(const std::string &name, double floatval);
  virtual void setSolverParam(const std::string &name, int32_t intval);
  virtual void setSolverParam(const std::string &name,
                              const std::string &strval);
  virtual void breakSolver();
  virtual void optserverHost(const std::string &addr);
  static std::shared_ptr<monty::ndarray<mosek::fusion::SolverStatus, 1>>
  solveBatch(
      bool israce, double maxtime, int32_t numthreads,
      std::shared_ptr<monty::ndarray<monty::rc_ptr<::mosek::fusion::Model>, 1>>
          models);
  virtual void solve(const std::string &addr, const std::string &accesstoken);
  virtual void solve();
  virtual void flushParameters();
  virtual void flushSolutions();
  virtual mosek::fusion::SolutionStatus getDualSolutionStatus();
  virtual mosek::fusion::ProblemStatus getProblemStatus();
  virtual mosek::fusion::SolutionStatus getPrimalSolutionStatus();
  virtual double dualObjValue();
  virtual double primalObjValue();
  virtual void selectedSolution(mosek::fusion::SolutionType soltype);
  virtual mosek::fusion::AccSolutionStatus getAcceptedSolutionStatus();
  virtual void acceptedSolutionStatus(mosek::fusion::AccSolutionStatus what);
  virtual mosek::fusion::ProblemStatus
  getProblemStatus(mosek::fusion::SolutionType which);
  virtual mosek::fusion::SolutionStatus
  getDualSolutionStatus(mosek::fusion::SolutionType which);
  virtual mosek::fusion::SolutionStatus
  getPrimalSolutionStatus(mosek::fusion::SolutionType which);
  virtual void updateObjective(monty::rc_ptr<::mosek::fusion::Expression> expr,
                               monty::rc_ptr<::mosek::fusion::Variable> x);
  virtual monty::rc_ptr<::mosek::fusion::Parameter>
  __mosek_2fusion_2Model__parameter(const std::string &name);
  monty::rc_ptr<::mosek::fusion::Parameter> parameter(const std::string &name);
  virtual monty::rc_ptr<::mosek::fusion::Parameter>
  __mosek_2fusion_2Model__parameter(const std::string &name, int32_t d1,
                                    int32_t d2, int32_t d3);
  monty::rc_ptr<::mosek::fusion::Parameter>
  parameter(const std::string &name, int32_t d1, int32_t d2, int32_t d3);
  virtual monty::rc_ptr<::mosek::fusion::Parameter>
  __mosek_2fusion_2Model__parameter(const std::string &name, int32_t d1,
                                    int32_t d2);
  monty::rc_ptr<::mosek::fusion::Parameter> parameter(const std::string &name,
                                                      int32_t d1, int32_t d2);
  virtual monty::rc_ptr<::mosek::fusion::Parameter>
  __mosek_2fusion_2Model__parameter(const std::string &name, int32_t d1);
  monty::rc_ptr<::mosek::fusion::Parameter> parameter(const std::string &name,
                                                      int32_t d1);
  virtual monty::rc_ptr<::mosek::fusion::Parameter>
  __mosek_2fusion_2Model__parameter(
      const std::string &name,
      std::shared_ptr<monty::ndarray<int32_t, 1>> shape);
  monty::rc_ptr<::mosek::fusion::Parameter>
  parameter(const std::string &name,
            std::shared_ptr<monty::ndarray<int32_t, 1>> shape);
  virtual monty::rc_ptr<::mosek::fusion::Parameter>
  __mosek_2fusion_2Model__parameter(
      const std::string &name,
      std::shared_ptr<monty::ndarray<int32_t, 1>> shape,
      std::shared_ptr<monty::ndarray<int64_t, 1>> sp);
  monty::rc_ptr<::mosek::fusion::Parameter>
  parameter(const std::string &name,
            std::shared_ptr<monty::ndarray<int32_t, 1>> shape,
            std::shared_ptr<monty::ndarray<int64_t, 1>> sp);
  virtual monty::rc_ptr<::mosek::fusion::Parameter>
  __mosek_2fusion_2Model__parameter(
      const std::string &name,
      std::shared_ptr<monty::ndarray<int32_t, 1>> shape,
      std::shared_ptr<monty::ndarray<int32_t, 2>> sparsity);
  monty::rc_ptr<::mosek::fusion::Parameter>
  parameter(const std::string &name,
            std::shared_ptr<monty::ndarray<int32_t, 1>> shape,
            std::shared_ptr<monty::ndarray<int32_t, 2>> sparsity);
  virtual monty::rc_ptr<::mosek::fusion::Parameter>
  __mosek_2fusion_2Model__parameter();
  monty::rc_ptr<::mosek::fusion::Parameter> parameter();
  virtual monty::rc_ptr<::mosek::fusion::Parameter>
  __mosek_2fusion_2Model__parameter(int32_t d1, int32_t d2, int32_t d3);
  monty::rc_ptr<::mosek::fusion::Parameter> parameter(int32_t d1, int32_t d2,
                                                      int32_t d3);
  virtual monty::rc_ptr<::mosek::fusion::Parameter>
  __mosek_2fusion_2Model__parameter(int32_t d1, int32_t d2);
  monty::rc_ptr<::mosek::fusion::Parameter> parameter(int32_t d1, int32_t d2);
  virtual monty::rc_ptr<::mosek::fusion::Parameter>
  __mosek_2fusion_2Model__parameter(int32_t d1);
  monty::rc_ptr<::mosek::fusion::Parameter> parameter(int32_t d1);
  virtual monty::rc_ptr<::mosek::fusion::Parameter>
  __mosek_2fusion_2Model__parameter(
      std::shared_ptr<monty::ndarray<int32_t, 1>> shape);
  monty::rc_ptr<::mosek::fusion::Parameter>
  parameter(std::shared_ptr<monty::ndarray<int32_t, 1>> shape);
  virtual monty::rc_ptr<::mosek::fusion::Parameter>
  __mosek_2fusion_2Model__parameter(
      std::shared_ptr<monty::ndarray<int32_t, 1>> shape,
      std::shared_ptr<monty::ndarray<int64_t, 1>> sp);
  monty::rc_ptr<::mosek::fusion::Parameter>
  parameter(std::shared_ptr<monty::ndarray<int32_t, 1>> shape,
            std::shared_ptr<monty::ndarray<int64_t, 1>> sp);
  virtual monty::rc_ptr<::mosek::fusion::Parameter>
  __mosek_2fusion_2Model__parameter(
      std::shared_ptr<monty::ndarray<int32_t, 1>> shape,
      std::shared_ptr<monty::ndarray<int32_t, 2>> sparsity);
  monty::rc_ptr<::mosek::fusion::Parameter>
  parameter(std::shared_ptr<monty::ndarray<int32_t, 1>> shape,
            std::shared_ptr<monty::ndarray<int32_t, 2>> sparsity);
  virtual void objective(double c);
  virtual void objective(mosek::fusion::ObjectiveSense sense, double c);
  virtual void objective(mosek::fusion::ObjectiveSense sense,
                         monty::rc_ptr<::mosek::fusion::Expression> expr);
  virtual void objective(const std::string &name, double c);
  virtual void objective(const std::string &name,
                         mosek::fusion::ObjectiveSense sense, double c);
  virtual void objective(const std::string &name,
                         mosek::fusion::ObjectiveSense sense,
                         monty::rc_ptr<::mosek::fusion::Expression> expr);
  virtual monty::rc_ptr<::mosek::fusion::RangedConstraint>
  __mosek_2fusion_2Model__constraint(
      monty::rc_ptr<::mosek::fusion::ExprRangeDomain> exprdom);
  monty::rc_ptr<::mosek::fusion::RangedConstraint>
  constraint(monty::rc_ptr<::mosek::fusion::ExprRangeDomain> exprdom);
  virtual monty::rc_ptr<::mosek::fusion::RangedConstraint>
  __mosek_2fusion_2Model__constraint(
      const std::string &name,
      monty::rc_ptr<::mosek::fusion::ExprRangeDomain> exprdom);
  monty::rc_ptr<::mosek::fusion::RangedConstraint>
  constraint(const std::string &name,
             monty::rc_ptr<::mosek::fusion::ExprRangeDomain> exprdom);
  virtual monty::rc_ptr<::mosek::fusion::RangedConstraint>
  __mosek_2fusion_2Model__constraint(
      monty::rc_ptr<::mosek::fusion::Expression> expr,
      monty::rc_ptr<::mosek::fusion::RangeDomain> rdom);
  monty::rc_ptr<::mosek::fusion::RangedConstraint>
  constraint(monty::rc_ptr<::mosek::fusion::Expression> expr,
             monty::rc_ptr<::mosek::fusion::RangeDomain> rdom);
  virtual monty::rc_ptr<::mosek::fusion::RangedConstraint>
  __mosek_2fusion_2Model__constraint(
      const std::string &name, monty::rc_ptr<::mosek::fusion::Expression> expr,
      monty::rc_ptr<::mosek::fusion::RangeDomain> rdom);
  monty::rc_ptr<::mosek::fusion::RangedConstraint>
  constraint(const std::string &name,
             monty::rc_ptr<::mosek::fusion::Expression> expr,
             monty::rc_ptr<::mosek::fusion::RangeDomain> rdom);
  virtual monty::rc_ptr<::mosek::fusion::Constraint>
  __mosek_2fusion_2Model__constraint(
      monty::rc_ptr<::mosek::fusion::ExprConicDomain> exprdom);
  monty::rc_ptr<::mosek::fusion::Constraint>
  constraint(monty::rc_ptr<::mosek::fusion::ExprConicDomain> exprdom);
  virtual monty::rc_ptr<::mosek::fusion::Constraint>
  __mosek_2fusion_2Model__constraint(
      const std::string &name,
      monty::rc_ptr<::mosek::fusion::ExprConicDomain> exprdom);
  monty::rc_ptr<::mosek::fusion::Constraint>
  constraint(const std::string &name,
             monty::rc_ptr<::mosek::fusion::ExprConicDomain> exprdom);
  virtual monty::rc_ptr<::mosek::fusion::Constraint>
  __mosek_2fusion_2Model__constraint(
      monty::rc_ptr<::mosek::fusion::Expression> expr,
      monty::rc_ptr<::mosek::fusion::ConeDomain> qdom);
  monty::rc_ptr<::mosek::fusion::Constraint>
  constraint(monty::rc_ptr<::mosek::fusion::Expression> expr,
             monty::rc_ptr<::mosek::fusion::ConeDomain> qdom);
  virtual monty::rc_ptr<::mosek::fusion::Constraint>
  __mosek_2fusion_2Model__constraint(
      const std::string &name, monty::rc_ptr<::mosek::fusion::Expression> expr,
      monty::rc_ptr<::mosek::fusion::ConeDomain> qdom);
  monty::rc_ptr<::mosek::fusion::Constraint>
  constraint(const std::string &name,
             monty::rc_ptr<::mosek::fusion::Expression> expr,
             monty::rc_ptr<::mosek::fusion::ConeDomain> qdom);
  virtual monty::rc_ptr<::mosek::fusion::Constraint>
  __mosek_2fusion_2Model__constraint(
      monty::rc_ptr<::mosek::fusion::ExprLinearDomain> exprdom);
  monty::rc_ptr<::mosek::fusion::Constraint>
  constraint(monty::rc_ptr<::mosek::fusion::ExprLinearDomain> exprdom);
  virtual monty::rc_ptr<::mosek::fusion::Constraint>
  __mosek_2fusion_2Model__constraint(
      const std::string &name,
      monty::rc_ptr<::mosek::fusion::ExprLinearDomain> exprdom);
  monty::rc_ptr<::mosek::fusion::Constraint>
  constraint(const std::string &name,
             monty::rc_ptr<::mosek::fusion::ExprLinearDomain> exprdom);
  virtual monty::rc_ptr<::mosek::fusion::Constraint>
  __mosek_2fusion_2Model__constraint(
      monty::rc_ptr<::mosek::fusion::Expression> expr,
      monty::rc_ptr<::mosek::fusion::LinearDomain> ldom);
  monty::rc_ptr<::mosek::fusion::Constraint>
  constraint(monty::rc_ptr<::mosek::fusion::Expression> expr,
             monty::rc_ptr<::mosek::fusion::LinearDomain> ldom);
  virtual monty::rc_ptr<::mosek::fusion::Constraint>
  __mosek_2fusion_2Model__constraint(
      const std::string &name, monty::rc_ptr<::mosek::fusion::Expression> expr,
      monty::rc_ptr<::mosek::fusion::LinearDomain> ldom);
  monty::rc_ptr<::mosek::fusion::Constraint>
  constraint(const std::string &name,
             monty::rc_ptr<::mosek::fusion::Expression> expr,
             monty::rc_ptr<::mosek::fusion::LinearDomain> ldom);
  virtual monty::rc_ptr<::mosek::fusion::Constraint>
  __mosek_2fusion_2Model__constraint(
      monty::rc_ptr<::mosek::fusion::ExprPSDDomain> exprdom);
  monty::rc_ptr<::mosek::fusion::Constraint>
  constraint(monty::rc_ptr<::mosek::fusion::ExprPSDDomain> exprdom);
  virtual monty::rc_ptr<::mosek::fusion::Constraint>
  __mosek_2fusion_2Model__constraint(
      const std::string &name,
      monty::rc_ptr<::mosek::fusion::ExprPSDDomain> exprdom);
  monty::rc_ptr<::mosek::fusion::Constraint>
  constraint(const std::string &name,
             monty::rc_ptr<::mosek::fusion::ExprPSDDomain> exprdom);
  virtual monty::rc_ptr<::mosek::fusion::Constraint>
  __mosek_2fusion_2Model__constraint(
      monty::rc_ptr<::mosek::fusion::Expression> expr,
      monty::rc_ptr<::mosek::fusion::PSDDomain> psddom);
  monty::rc_ptr<::mosek::fusion::Constraint>
  constraint(monty::rc_ptr<::mosek::fusion::Expression> expr,
             monty::rc_ptr<::mosek::fusion::PSDDomain> psddom);
  virtual monty::rc_ptr<::mosek::fusion::Constraint>
  __mosek_2fusion_2Model__constraint(
      const std::string &name, monty::rc_ptr<::mosek::fusion::Expression> expr,
      monty::rc_ptr<::mosek::fusion::PSDDomain> psddom);
  monty::rc_ptr<::mosek::fusion::Constraint>
  constraint(const std::string &name,
             monty::rc_ptr<::mosek::fusion::Expression> expr,
             monty::rc_ptr<::mosek::fusion::PSDDomain> psddom);
  virtual monty::rc_ptr<::mosek::fusion::Variable>
  __mosek_2fusion_2Model__variable(
      monty::rc_ptr<::mosek::fusion::PSDDomain> psddom);
  monty::rc_ptr<::mosek::fusion::Variable>
  variable(monty::rc_ptr<::mosek::fusion::PSDDomain> psddom);
  virtual monty::rc_ptr<::mosek::fusion::Variable>
  __mosek_2fusion_2Model__variable(
      int32_t n, int32_t m, monty::rc_ptr<::mosek::fusion::PSDDomain> psddom);
  monty::rc_ptr<::mosek::fusion::Variable>
  variable(int32_t n, int32_t m,
           monty::rc_ptr<::mosek::fusion::PSDDomain> psddom);
  virtual monty::rc_ptr<::mosek::fusion::Variable>
  __mosek_2fusion_2Model__variable(
      int32_t n, monty::rc_ptr<::mosek::fusion::PSDDomain> psddom);
  monty::rc_ptr<::mosek::fusion::Variable>
  variable(int32_t n, monty::rc_ptr<::mosek::fusion::PSDDomain> psddom);
  virtual monty::rc_ptr<::mosek::fusion::Variable>
  __mosek_2fusion_2Model__variable(
      const std::string &name,
      monty::rc_ptr<::mosek::fusion::PSDDomain> psddom);
  monty::rc_ptr<::mosek::fusion::Variable>
  variable(const std::string &name,
           monty::rc_ptr<::mosek::fusion::PSDDomain> psddom);
  virtual monty::rc_ptr<::mosek::fusion::Variable>
  __mosek_2fusion_2Model__variable(
      const std::string &name, int32_t n, int32_t m,
      monty::rc_ptr<::mosek::fusion::PSDDomain> psddom);
  monty::rc_ptr<::mosek::fusion::Variable>
  variable(const std::string &name, int32_t n, int32_t m,
           monty::rc_ptr<::mosek::fusion::PSDDomain> psddom);
  virtual monty::rc_ptr<::mosek::fusion::Variable>
  __mosek_2fusion_2Model__variable(
      const std::string &name, int32_t n,
      monty::rc_ptr<::mosek::fusion::PSDDomain> psddom);
  monty::rc_ptr<::mosek::fusion::Variable>
  variable(const std::string &name, int32_t n,
           monty::rc_ptr<::mosek::fusion::PSDDomain> psddom);
  virtual monty::rc_ptr<::mosek::fusion::Variable>
  __mosek_2fusion_2Model__variable(
      const std::string &name, std::shared_ptr<monty::ndarray<int32_t, 1>> shp,
      monty::rc_ptr<::mosek::fusion::PSDDomain> psddom);
  monty::rc_ptr<::mosek::fusion::Variable>
  variable(const std::string &name,
           std::shared_ptr<monty::ndarray<int32_t, 1>> shp,
           monty::rc_ptr<::mosek::fusion::PSDDomain> psddom);
  virtual monty::rc_ptr<::mosek::fusion::Variable>
  __mosek_2fusion_2Model__variable(
      monty::rc_ptr<::mosek::fusion::ConeDomain> qdom);
  monty::rc_ptr<::mosek::fusion::Variable>
  variable(monty::rc_ptr<::mosek::fusion::ConeDomain> qdom);
  virtual monty::rc_ptr<::mosek::fusion::RangedVariable>
  __mosek_2fusion_2Model__variable(
      monty::rc_ptr<::mosek::fusion::RangeDomain> rdom);
  monty::rc_ptr<::mosek::fusion::RangedVariable>
  variable(monty::rc_ptr<::mosek::fusion::RangeDomain> rdom);
  virtual monty::rc_ptr<::mosek::fusion::Variable>
  __mosek_2fusion_2Model__variable(
      monty::rc_ptr<::mosek::fusion::LinearDomain> ldom);
  monty::rc_ptr<::mosek::fusion::Variable>
  variable(monty::rc_ptr<::mosek::fusion::LinearDomain> ldom);
  virtual monty::rc_ptr<::mosek::fusion::Variable>
  __mosek_2fusion_2Model__variable(
      std::shared_ptr<monty::ndarray<int32_t, 1>> shp,
      monty::rc_ptr<::mosek::fusion::ConeDomain> qdom);
  monty::rc_ptr<::mosek::fusion::Variable>
  variable(std::shared_ptr<monty::ndarray<int32_t, 1>> shp,
           monty::rc_ptr<::mosek::fusion::ConeDomain> qdom);
  virtual monty::rc_ptr<::mosek::fusion::RangedVariable>
  __mosek_2fusion_2Model__variable(
      std::shared_ptr<monty::ndarray<int32_t, 1>> shp,
      monty::rc_ptr<::mosek::fusion::RangeDomain> rdom);
  monty::rc_ptr<::mosek::fusion::RangedVariable>
  variable(std::shared_ptr<monty::ndarray<int32_t, 1>> shp,
           monty::rc_ptr<::mosek::fusion::RangeDomain> rdom);
  virtual monty::rc_ptr<::mosek::fusion::Variable>
  __mosek_2fusion_2Model__variable(
      std::shared_ptr<monty::ndarray<int32_t, 1>> shp,
      monty::rc_ptr<::mosek::fusion::LinearDomain> ldom);
  monty::rc_ptr<::mosek::fusion::Variable>
  variable(std::shared_ptr<monty::ndarray<int32_t, 1>> shp,
           monty::rc_ptr<::mosek::fusion::LinearDomain> ldom);
  virtual monty::rc_ptr<::mosek::fusion::Variable>
  __mosek_2fusion_2Model__variable(
      std::shared_ptr<monty::ndarray<int32_t, 1>> shp);
  monty::rc_ptr<::mosek::fusion::Variable>
  variable(std::shared_ptr<monty::ndarray<int32_t, 1>> shp);
  virtual monty::rc_ptr<::mosek::fusion::Variable>
  __mosek_2fusion_2Model__variable(
      int32_t size, monty::rc_ptr<::mosek::fusion::ConeDomain> qdom);
  monty::rc_ptr<::mosek::fusion::Variable>
  variable(int32_t size, monty::rc_ptr<::mosek::fusion::ConeDomain> qdom);
  virtual monty::rc_ptr<::mosek::fusion::RangedVariable>
  __mosek_2fusion_2Model__variable(
      int32_t size, monty::rc_ptr<::mosek::fusion::RangeDomain> rdom);
  monty::rc_ptr<::mosek::fusion::RangedVariable>
  variable(int32_t size, monty::rc_ptr<::mosek::fusion::RangeDomain> rdom);
  virtual monty::rc_ptr<::mosek::fusion::Variable>
  __mosek_2fusion_2Model__variable(
      int32_t size, monty::rc_ptr<::mosek::fusion::LinearDomain> ldom);
  monty::rc_ptr<::mosek::fusion::Variable>
  variable(int32_t size, monty::rc_ptr<::mosek::fusion::LinearDomain> ldom);
  virtual monty::rc_ptr<::mosek::fusion::Variable>
  __mosek_2fusion_2Model__variable(int32_t size);
  monty::rc_ptr<::mosek::fusion::Variable> variable(int32_t size);
  virtual monty::rc_ptr<::mosek::fusion::Variable>
  __mosek_2fusion_2Model__variable();
  monty::rc_ptr<::mosek::fusion::Variable> variable();
  virtual monty::rc_ptr<::mosek::fusion::Variable>
  __mosek_2fusion_2Model__variable(
      const std::string &name, monty::rc_ptr<::mosek::fusion::ConeDomain> qdom);
  monty::rc_ptr<::mosek::fusion::Variable>
  variable(const std::string &name,
           monty::rc_ptr<::mosek::fusion::ConeDomain> qdom);
  virtual monty::rc_ptr<::mosek::fusion::RangedVariable>
  __mosek_2fusion_2Model__variable(
      const std::string &name,
      monty::rc_ptr<::mosek::fusion::RangeDomain> rdom);
  monty::rc_ptr<::mosek::fusion::RangedVariable>
  variable(const std::string &name,
           monty::rc_ptr<::mosek::fusion::RangeDomain> rdom);
  virtual monty::rc_ptr<::mosek::fusion::Variable>
  __mosek_2fusion_2Model__variable(
      const std::string &name,
      monty::rc_ptr<::mosek::fusion::LinearDomain> ldom);
  monty::rc_ptr<::mosek::fusion::Variable>
  variable(const std::string &name,
           monty::rc_ptr<::mosek::fusion::LinearDomain> ldom);
  virtual monty::rc_ptr<::mosek::fusion::Variable>
  __mosek_2fusion_2Model__variable(
      const std::string &name, std::shared_ptr<monty::ndarray<int32_t, 1>> shp,
      monty::rc_ptr<::mosek::fusion::ConeDomain> qdom);
  monty::rc_ptr<::mosek::fusion::Variable>
  variable(const std::string &name,
           std::shared_ptr<monty::ndarray<int32_t, 1>> shp,
           monty::rc_ptr<::mosek::fusion::ConeDomain> qdom);
  virtual monty::rc_ptr<::mosek::fusion::RangedVariable>
  __mosek_2fusion_2Model__variable(
      const std::string &name, std::shared_ptr<monty::ndarray<int32_t, 1>> shp,
      monty::rc_ptr<::mosek::fusion::RangeDomain> rdom);
  monty::rc_ptr<::mosek::fusion::RangedVariable>
  variable(const std::string &name,
           std::shared_ptr<monty::ndarray<int32_t, 1>> shp,
           monty::rc_ptr<::mosek::fusion::RangeDomain> rdom);
  virtual monty::rc_ptr<::mosek::fusion::Variable>
  __mosek_2fusion_2Model__variable(
      const std::string &name, std::shared_ptr<monty::ndarray<int32_t, 1>> shp,
      monty::rc_ptr<::mosek::fusion::LinearDomain> ldom);
  monty::rc_ptr<::mosek::fusion::Variable>
  variable(const std::string &name,
           std::shared_ptr<monty::ndarray<int32_t, 1>> shp,
           monty::rc_ptr<::mosek::fusion::LinearDomain> ldom);
  virtual monty::rc_ptr<::mosek::fusion::Variable>
  __mosek_2fusion_2Model__variable(
      const std::string &name, std::shared_ptr<monty::ndarray<int32_t, 1>> shp);
  monty::rc_ptr<::mosek::fusion::Variable>
  variable(const std::string &name,
           std::shared_ptr<monty::ndarray<int32_t, 1>> shp);
  virtual monty::rc_ptr<::mosek::fusion::Variable>
  __mosek_2fusion_2Model__variable(
      const std::string &name, int32_t size,
      monty::rc_ptr<::mosek::fusion::ConeDomain> qdom);
  monty::rc_ptr<::mosek::fusion::Variable>
  variable(const std::string &name, int32_t size,
           monty::rc_ptr<::mosek::fusion::ConeDomain> qdom);
  virtual monty::rc_ptr<::mosek::fusion::RangedVariable>
  __mosek_2fusion_2Model__variable(
      const std::string &name, int32_t size,
      monty::rc_ptr<::mosek::fusion::RangeDomain> rdom);
  monty::rc_ptr<::mosek::fusion::RangedVariable>
  variable(const std::string &name, int32_t size,
           monty::rc_ptr<::mosek::fusion::RangeDomain> rdom);
  virtual monty::rc_ptr<::mosek::fusion::Variable>
  __mosek_2fusion_2Model__variable(
      const std::string &name, int32_t size,
      monty::rc_ptr<::mosek::fusion::LinearDomain> ldom);
  monty::rc_ptr<::mosek::fusion::Variable>
  variable(const std::string &name, int32_t size,
           monty::rc_ptr<::mosek::fusion::LinearDomain> ldom);
  virtual monty::rc_ptr<::mosek::fusion::Variable>
  __mosek_2fusion_2Model__variable(const std::string &name, int32_t size);
  monty::rc_ptr<::mosek::fusion::Variable> variable(const std::string &name,
                                                    int32_t size);
  virtual monty::rc_ptr<::mosek::fusion::Variable>
  __mosek_2fusion_2Model__variable(const std::string &name);
  monty::rc_ptr<::mosek::fusion::Variable> variable(const std::string &name);
  static std::string getVersion();
  virtual bool hasParameter(const std::string &name);
  virtual bool hasConstraint(const std::string &name);
  virtual bool hasVariable(const std::string &name);
  virtual monty::rc_ptr<::mosek::fusion::Parameter>
  __mosek_2fusion_2Model__getParameter(const std::string &name);
  monty::rc_ptr<::mosek::fusion::Parameter>
  getParameter(const std::string &name);
  virtual monty::rc_ptr<::mosek::fusion::Constraint>
  __mosek_2fusion_2Model__getConstraint(int32_t index);
  monty::rc_ptr<::mosek::fusion::Constraint> getConstraint(int32_t index);
  virtual monty::rc_ptr<::mosek::fusion::Constraint>
  __mosek_2fusion_2Model__getConstraint(const std::string &name);
  monty::rc_ptr<::mosek::fusion::Constraint>
  getConstraint(const std::string &name);
  virtual monty::rc_ptr<::mosek::fusion::Variable>
  __mosek_2fusion_2Model__getVariable(int32_t index);
  monty::rc_ptr<::mosek::fusion::Variable> getVariable(int32_t index);
  virtual monty::rc_ptr<::mosek::fusion::Variable>
  __mosek_2fusion_2Model__getVariable(const std::string &name);
  monty::rc_ptr<::mosek::fusion::Variable> getVariable(const std::string &name);
  virtual std::string getName();
  virtual monty::rc_ptr<::mosek::fusion::Model> __mosek_2fusion_2Model__clone();
  monty::rc_ptr<::mosek::fusion::Model> clone();
}; // class Model;

// mosek.fusion.Debug from file 'src/fusion/cxx/Debug.h'
class Debug : public monty::RefCounted {
  p_Debug *ptr;

public:
  friend struct p_Debug;
  typedef monty::rc_ptr<Debug> t;

  virtual ~Debug();

  Debug();

  static t o();
  t p(const std::string &val);
  t p(int val);
  t p(int64_t val);
  t p(double val);
  t p(bool val);
  t lf();

  t p(const std::shared_ptr<monty::ndarray<double, 1>> &val);
  t p(const std::shared_ptr<monty::ndarray<int, 1>> &val);
  t p(const std::shared_ptr<monty::ndarray<int64_t, 1>> &val);

  t __mosek_2fusion_2Debug__p(const std::string &val) { return p(val); }
  t __mosek_2fusion_2Debug__p(int val) { return p(val); }
  t __mosek_2fusion_2Debug__p(int64_t val) { return p(val); }
  t __mosek_2fusion_2Debug__p(double val) { return p(val); }
  t __mosek_2fusion_2Debug__p(bool val) { return p(val); }

  t __mosek_2fusion_2Debug__p(
      const std::shared_ptr<monty::ndarray<double, 1>> &val) {
    return p(val);
  }
  t __mosek_2fusion_2Debug__p(
      const std::shared_ptr<monty::ndarray<int, 1>> &val) {
    return p(val);
  }
  t __mosek_2fusion_2Debug__p(
      const std::shared_ptr<monty::ndarray<int64_t, 1>> &val) {
    return p(val);
  }

  t __mosek_2fusion_2Debug__lf() { return lf(); }
};
// End of file 'src/fusion/cxx/Debug.h'
class Sort : public virtual monty::RefCounted {
public:
  p_Sort *_impl;

protected:
  Sort(p_Sort *_impl);

public:
  Sort(const Sort &) = delete;
  const Sort &operator=(const Sort &) = delete;
  friend class p_Sort;
  virtual ~Sort();
  virtual void destroy();
  typedef monty::rc_ptr<Sort> t;

  static void argTransposeSort(std::shared_ptr<monty::ndarray<int64_t, 1>> perm,
                               std::shared_ptr<monty::ndarray<int64_t, 1>> ptrb,
                               int32_t m, int32_t n, int32_t p,
                               std::shared_ptr<monty::ndarray<int64_t, 1>> val);
  static void argsort(std::shared_ptr<monty::ndarray<int64_t, 1>> idx,
                      std::shared_ptr<monty::ndarray<int64_t, 1>> vals1);
  static void argsort(std::shared_ptr<monty::ndarray<int64_t, 1>> idx,
                      std::shared_ptr<monty::ndarray<int32_t, 1>> vals1);
  static void argsort(std::shared_ptr<monty::ndarray<int64_t, 1>> idx,
                      std::shared_ptr<monty::ndarray<int64_t, 1>> vals1,
                      std::shared_ptr<monty::ndarray<int64_t, 1>> vals2);
  static void argsort(std::shared_ptr<monty::ndarray<int64_t, 1>> idx,
                      std::shared_ptr<monty::ndarray<int32_t, 1>> vals1,
                      std::shared_ptr<monty::ndarray<int32_t, 1>> vals2);
  static void argsort(std::shared_ptr<monty::ndarray<int64_t, 1>> idx,
                      std::shared_ptr<monty::ndarray<int64_t, 1>> vals1,
                      int64_t first, int64_t last);
  static void argsort(std::shared_ptr<monty::ndarray<int64_t, 1>> idx,
                      std::shared_ptr<monty::ndarray<int32_t, 1>> vals1,
                      int64_t first, int64_t last);
  static void argsort(std::shared_ptr<monty::ndarray<int64_t, 1>> idx,
                      std::shared_ptr<monty::ndarray<int64_t, 1>> vals1,
                      std::shared_ptr<monty::ndarray<int64_t, 1>> vals2,
                      int64_t first, int64_t last);
  static void argsort(std::shared_ptr<monty::ndarray<int64_t, 1>> idx,
                      std::shared_ptr<monty::ndarray<int32_t, 1>> vals1,
                      std::shared_ptr<monty::ndarray<int32_t, 1>> vals2,
                      int64_t first, int64_t last);
  static void argsort(std::shared_ptr<monty::ndarray<int64_t, 1>> idx,
                      std::shared_ptr<monty::ndarray<int64_t, 1>> vals1,
                      int64_t first, int64_t last, bool check);
  static void argsort(std::shared_ptr<monty::ndarray<int64_t, 1>> idx,
                      std::shared_ptr<monty::ndarray<int32_t, 1>> vals1,
                      int64_t first, int64_t last, bool check);
  static void argsort(std::shared_ptr<monty::ndarray<int64_t, 1>> idx,
                      std::shared_ptr<monty::ndarray<int64_t, 1>> vals1,
                      std::shared_ptr<monty::ndarray<int64_t, 1>> vals2,
                      int64_t first, int64_t last, bool check);
  static void argsort(std::shared_ptr<monty::ndarray<int64_t, 1>> idx,
                      std::shared_ptr<monty::ndarray<int32_t, 1>> vals1,
                      std::shared_ptr<monty::ndarray<int32_t, 1>> vals2,
                      int64_t first, int64_t last, bool check);
  static void argbucketsort(std::shared_ptr<monty::ndarray<int64_t, 1>> idx,
                            std::shared_ptr<monty::ndarray<int64_t, 1>> vals,
                            int64_t first, int64_t last, int64_t minv,
                            int64_t maxv);
  static void argbucketsort(std::shared_ptr<monty::ndarray<int64_t, 1>> idx,
                            std::shared_ptr<monty::ndarray<int32_t, 1>> vals,
                            int64_t first, int64_t last, int32_t minv,
                            int32_t maxv);
  static void getminmax(std::shared_ptr<monty::ndarray<int64_t, 1>> idx,
                        std::shared_ptr<monty::ndarray<int64_t, 1>> vals1,
                        std::shared_ptr<monty::ndarray<int64_t, 1>> vals2,
                        int64_t first, int64_t last,
                        std::shared_ptr<monty::ndarray<int64_t, 1>> res);
  static void getminmax(std::shared_ptr<monty::ndarray<int64_t, 1>> idx,
                        std::shared_ptr<monty::ndarray<int32_t, 1>> vals1,
                        std::shared_ptr<monty::ndarray<int32_t, 1>> vals2,
                        int64_t first, int64_t last,
                        std::shared_ptr<monty::ndarray<int32_t, 1>> res);
  static bool issorted(std::shared_ptr<monty::ndarray<int64_t, 1>> idx,
                       std::shared_ptr<monty::ndarray<int64_t, 1>> vals1,
                       int64_t first, int64_t last, bool check);
  static bool issorted(std::shared_ptr<monty::ndarray<int64_t, 1>> idx,
                       std::shared_ptr<monty::ndarray<int32_t, 1>> vals1,
                       int64_t first, int64_t last, bool check);
  static bool issorted(std::shared_ptr<monty::ndarray<int64_t, 1>> idx,
                       std::shared_ptr<monty::ndarray<int64_t, 1>> vals1,
                       std::shared_ptr<monty::ndarray<int64_t, 1>> vals2,
                       int64_t first, int64_t last, bool check);
  static bool issorted(std::shared_ptr<monty::ndarray<int64_t, 1>> idx,
                       std::shared_ptr<monty::ndarray<int32_t, 1>> vals1,
                       std::shared_ptr<monty::ndarray<int32_t, 1>> vals2,
                       int64_t first, int64_t last, bool check);
}; // class Sort;

class IndexCounter : public virtual monty::RefCounted {
public:
  p_IndexCounter *_impl;

protected:
  IndexCounter(p_IndexCounter *_impl);

public:
  IndexCounter(const IndexCounter &) = delete;
  const IndexCounter &operator=(const IndexCounter &) = delete;
  friend class p_IndexCounter;
  virtual ~IndexCounter();
  virtual void destroy();
  typedef monty::rc_ptr<IndexCounter> t;

  IndexCounter(std::shared_ptr<monty::ndarray<int32_t, 1>> shape);
  IndexCounter(int64_t start_,
               std::shared_ptr<monty::ndarray<int32_t, 1>> dims_,
               std::shared_ptr<monty::ndarray<int32_t, 1>> shape);
  IndexCounter(int64_t start_,
               std::shared_ptr<monty::ndarray<int32_t, 1>> dims_,
               std::shared_ptr<monty::ndarray<int64_t, 1>> strides_);
  virtual bool atEnd();
  virtual std::shared_ptr<monty::ndarray<int32_t, 1>> getIndex();
  virtual int64_t next();
  virtual int64_t get();
  virtual void inc();
  virtual void reset();
}; // class IndexCounter;

class CommonTools : public virtual monty::RefCounted {
public:
  p_CommonTools *_impl;

protected:
  CommonTools(p_CommonTools *_impl);

public:
  CommonTools(const CommonTools &) = delete;
  const CommonTools &operator=(const CommonTools &) = delete;
  friend class p_CommonTools;
  virtual ~CommonTools();
  virtual void destroy();
  typedef monty::rc_ptr<CommonTools> t;

  static std::shared_ptr<monty::ndarray<int64_t, 1>>
  resize(std::shared_ptr<monty::ndarray<int64_t, 1>> values, int32_t newsize);
  static std::shared_ptr<monty::ndarray<int32_t, 1>>
  resize(std::shared_ptr<monty::ndarray<int32_t, 1>> values, int32_t newsize);
  static std::shared_ptr<monty::ndarray<double, 1>>
  resize(std::shared_ptr<monty::ndarray<double, 1>> values, int32_t newsize);
  static int32_t
  binarySearch(std::shared_ptr<monty::ndarray<int32_t, 1>> values,
               int32_t target);
  static int32_t
  binarySearch(std::shared_ptr<monty::ndarray<int64_t, 1>> values,
               int64_t target);
  static int32_t
  binarySearchR(std::shared_ptr<monty::ndarray<int64_t, 1>> values,
                int64_t target);
  static int32_t
  binarySearchL(std::shared_ptr<monty::ndarray<int64_t, 1>> values,
                int64_t target);
  static void ndIncr(std::shared_ptr<monty::ndarray<int32_t, 1>> ndidx,
                     std::shared_ptr<monty::ndarray<int32_t, 1>> first,
                     std::shared_ptr<monty::ndarray<int32_t, 1>> last);
  static void transposeTriplets(
      std::shared_ptr<monty::ndarray<int32_t, 1>> subi,
      std::shared_ptr<monty::ndarray<int32_t, 1>> subj,
      std::shared_ptr<monty::ndarray<double, 1>> val,
      std::shared_ptr<
          monty::ndarray<std::shared_ptr<monty::ndarray<int64_t, 1>>, 1>>
          tsubi_,
      std::shared_ptr<
          monty::ndarray<std::shared_ptr<monty::ndarray<int64_t, 1>>, 1>>
          tsubj_,
      std::shared_ptr<
          monty::ndarray<std::shared_ptr<monty::ndarray<double, 1>>, 1>>
          tval_,
      int64_t nelm, int32_t dimi, int32_t dimj);
  static void transposeTriplets(
      std::shared_ptr<monty::ndarray<int32_t, 1>> subi,
      std::shared_ptr<monty::ndarray<int32_t, 1>> subj,
      std::shared_ptr<monty::ndarray<double, 1>> val,
      std::shared_ptr<
          monty::ndarray<std::shared_ptr<monty::ndarray<int32_t, 1>>, 1>>
          tsubi_,
      std::shared_ptr<
          monty::ndarray<std::shared_ptr<monty::ndarray<int32_t, 1>>, 1>>
          tsubj_,
      std::shared_ptr<
          monty::ndarray<std::shared_ptr<monty::ndarray<double, 1>>, 1>>
          tval_,
      int64_t nelm, int32_t dimi, int32_t dimj);
  static void tripletSort(
      std::shared_ptr<monty::ndarray<int32_t, 1>> subi,
      std::shared_ptr<monty::ndarray<int32_t, 1>> subj,
      std::shared_ptr<monty::ndarray<double, 1>> val,
      std::shared_ptr<
          monty::ndarray<std::shared_ptr<monty::ndarray<int32_t, 1>>, 1>>
          tsubi_,
      std::shared_ptr<
          monty::ndarray<std::shared_ptr<monty::ndarray<int32_t, 1>>, 1>>
          tsubj_,
      std::shared_ptr<
          monty::ndarray<std::shared_ptr<monty::ndarray<double, 1>>, 1>>
          tval_,
      int64_t nelm, int32_t dimi, int32_t dimj);
  static void argMSort(std::shared_ptr<monty::ndarray<int32_t, 1>> idx,
                       std::shared_ptr<monty::ndarray<int32_t, 1>> vals);
  static void argQsort(std::shared_ptr<monty::ndarray<int64_t, 1>> idx,
                       std::shared_ptr<monty::ndarray<int64_t, 1>> vals1,
                       std::shared_ptr<monty::ndarray<int64_t, 1>> vals2,
                       int64_t first, int64_t last);
  static void argQsort(std::shared_ptr<monty::ndarray<int64_t, 1>> idx,
                       std::shared_ptr<monty::ndarray<int32_t, 1>> vals1,
                       std::shared_ptr<monty::ndarray<int32_t, 1>> vals2,
                       int64_t first, int64_t last);
}; // class CommonTools;

class SolutionStruct : public virtual monty::RefCounted {
public:
  p_SolutionStruct *_impl;

protected:
  SolutionStruct(p_SolutionStruct *_impl);

public:
  SolutionStruct(const SolutionStruct &) = delete;
  const SolutionStruct &operator=(const SolutionStruct &) = delete;
  friend class p_SolutionStruct;
  virtual ~SolutionStruct();
  virtual void destroy();
  typedef monty::rc_ptr<SolutionStruct> t;
  std::shared_ptr<monty::ndarray<double, 1>> get_accy();
  void set_accy(std::shared_ptr<monty::ndarray<double, 1>> val);
  std::shared_ptr<monty::ndarray<double, 1>> get_accx();
  void set_accx(std::shared_ptr<monty::ndarray<double, 1>> val);
  std::shared_ptr<monty::ndarray<int32_t, 1>> get_accptr();
  void set_accptr(std::shared_ptr<monty::ndarray<int32_t, 1>> val);
  std::shared_ptr<monty::ndarray<double, 1>> get_yx();
  void set_yx(std::shared_ptr<monty::ndarray<double, 1>> val);
  std::shared_ptr<monty::ndarray<double, 1>> get_sux();
  void set_sux(std::shared_ptr<monty::ndarray<double, 1>> val);
  std::shared_ptr<monty::ndarray<double, 1>> get_slx();
  void set_slx(std::shared_ptr<monty::ndarray<double, 1>> val);
  std::shared_ptr<monty::ndarray<double, 1>> get_bars();
  void set_bars(std::shared_ptr<monty::ndarray<double, 1>> val);
  std::shared_ptr<monty::ndarray<double, 1>> get_barx();
  void set_barx(std::shared_ptr<monty::ndarray<double, 1>> val);
  std::shared_ptr<monty::ndarray<double, 1>> get_y();
  void set_y(std::shared_ptr<monty::ndarray<double, 1>> val);
  std::shared_ptr<monty::ndarray<double, 1>> get_suc();
  void set_suc(std::shared_ptr<monty::ndarray<double, 1>> val);
  std::shared_ptr<monty::ndarray<double, 1>> get_slc();
  void set_slc(std::shared_ptr<monty::ndarray<double, 1>> val);
  std::shared_ptr<monty::ndarray<double, 1>> get_xx();
  void set_xx(std::shared_ptr<monty::ndarray<double, 1>> val);
  std::shared_ptr<monty::ndarray<double, 1>> get_xc();
  void set_xc(std::shared_ptr<monty::ndarray<double, 1>> val);
  double get_dobj();
  void set_dobj(double val);
  double get_pobj();
  void set_pobj(double val);
  mosek::fusion::ProblemStatus get_probstatus();
  void set_probstatus(mosek::fusion::ProblemStatus val);
  mosek::fusion::SolutionStatus get_dstatus();
  void set_dstatus(mosek::fusion::SolutionStatus val);
  mosek::fusion::SolutionStatus get_pstatus();
  void set_pstatus(mosek::fusion::SolutionStatus val);
  int32_t get_sol_numbarvar();
  void set_sol_numbarvar(int32_t val);
  int32_t get_sol_numaccelm();
  void set_sol_numaccelm(int32_t val);
  int32_t get_sol_numacc();
  void set_sol_numacc(int32_t val);
  int32_t get_sol_numvar();
  void set_sol_numvar(int32_t val);
  int32_t get_sol_numcon();
  void set_sol_numcon(int32_t val);

  SolutionStruct(int32_t numvar, int32_t numcon, int32_t numbarvar,
                 int32_t numacc, int32_t numaccelm);
  SolutionStruct(monty::rc_ptr<::mosek::fusion::SolutionStruct> that);
  virtual monty::rc_ptr<::mosek::fusion::SolutionStruct>
  __mosek_2fusion_2SolutionStruct__clone();
  monty::rc_ptr<::mosek::fusion::SolutionStruct> clone();
  virtual void resize(int32_t numvar, int32_t numcon, int32_t numbarvar,
                      int32_t numacc, int32_t numaccelm);
  virtual bool
  isDualAcceptable(mosek::fusion::AccSolutionStatus acceptable_sol);
  virtual bool
  isPrimalAcceptable(mosek::fusion::AccSolutionStatus acceptable_sol);
}; // class SolutionStruct;

class RowBlockManager : public virtual monty::RefCounted {
public:
  p_RowBlockManager *_impl;

protected:
  RowBlockManager(p_RowBlockManager *_impl);

public:
  RowBlockManager(const RowBlockManager &) = delete;
  const RowBlockManager &operator=(const RowBlockManager &) = delete;
  friend class p_RowBlockManager;
  virtual ~RowBlockManager();
  virtual void destroy();
  typedef monty::rc_ptr<RowBlockManager> t;
  int32_t get_varidx_used();
  void set_varidx_used(int32_t val);
  int32_t get_code_used();
  void set_code_used(int32_t val);
  std::shared_ptr<monty::ndarray<double, 1>> get_cconst();
  void set_cconst(std::shared_ptr<monty::ndarray<double, 1>> val);
  std::shared_ptr<monty::ndarray<int32_t, 1>> get_code();
  void set_code(std::shared_ptr<monty::ndarray<int32_t, 1>> val);
  int32_t get_first_free_codeitem();
  void set_first_free_codeitem(int32_t val);
  std::shared_ptr<monty::ndarray<int32_t, 1>> get_param_code_sizes();
  void set_param_code_sizes(std::shared_ptr<monty::ndarray<int32_t, 1>> val);
  std::shared_ptr<monty::ndarray<int64_t, 1>> get_param_varidx();
  void set_param_varidx(std::shared_ptr<monty::ndarray<int64_t, 1>> val);
  int32_t get_first_free_entry();
  void set_first_free_entry(int32_t val);
  std::shared_ptr<monty::ndarray<int32_t, 1>> get_row_code_ptr();
  void set_row_code_ptr(std::shared_ptr<monty::ndarray<int32_t, 1>> val);
  std::shared_ptr<monty::ndarray<int32_t, 1>> get_row_param_ptre();
  void set_row_param_ptre(std::shared_ptr<monty::ndarray<int32_t, 1>> val);
  std::shared_ptr<monty::ndarray<int32_t, 1>> get_row_param_ptrb();
  void set_row_param_ptrb(std::shared_ptr<monty::ndarray<int32_t, 1>> val);
  monty::rc_ptr<::mosek::fusion::LinkedBlocks> get_blocks();
  void set_blocks(monty::rc_ptr<::mosek::fusion::LinkedBlocks> val);

  RowBlockManager(monty::rc_ptr<::mosek::fusion::RowBlockManager> that);
  RowBlockManager();
  virtual int32_t num_parameterized();
  virtual bool is_parameterized();
  virtual int32_t blocksize(int32_t id);
  virtual int32_t block_capacity();
  virtual int32_t capacity();
  virtual void get(int32_t id,
                   std::shared_ptr<monty::ndarray<int32_t, 1>> target,
                   int32_t offset);
  virtual void evaluate(monty::rc_ptr<::mosek::fusion::WorkStack> xs,
                        std::shared_ptr<monty::ndarray<double, 1>> param_value,
                        std::shared_ptr<monty::ndarray<int32_t, 1>> subi,
                        std::shared_ptr<monty::ndarray<int32_t, 1>> subj,
                        std::shared_ptr<monty::ndarray<double, 1>> val);
  virtual void
  replace_row_code(monty::rc_ptr<::mosek::fusion::WorkStack> rs,
                   std::shared_ptr<monty::ndarray<int32_t, 1>> nativeidxs,
                   int32_t ptr, int32_t nidxs, int32_t codeptr, int32_t code_p,
                   int32_t cconst_p);
  virtual void
  clear_row_code(std::shared_ptr<monty::ndarray<int32_t, 1>> nativeidxs);
  virtual void release(int32_t id,
                       std::shared_ptr<monty::ndarray<int32_t, 1>> nativeidxs);
  virtual int32_t
  allocate(std::shared_ptr<monty::ndarray<int32_t, 1>> nativeidxs);
  virtual bool row_is_parameterized(int32_t i);
}; // class RowBlockManager;

class BaseVariable : public virtual ::mosek::fusion::Variable {
public:
  p_BaseVariable *_impl;

protected:
  BaseVariable(p_BaseVariable *_impl);

public:
  BaseVariable(const BaseVariable &) = delete;
  const BaseVariable &operator=(const BaseVariable &) = delete;
  friend class p_BaseVariable;
  virtual ~BaseVariable();
  virtual void destroy();
  typedef monty::rc_ptr<BaseVariable> t;

  BaseVariable(monty::rc_ptr<::mosek::fusion::Model> m,
               std::shared_ptr<monty::ndarray<int32_t, 1>> shape,
               std::shared_ptr<monty::ndarray<int64_t, 1>> sparsity,
               std::shared_ptr<monty::ndarray<int64_t, 1>> basevar_nativeidxs);
  virtual /* override */ std::string toString();
  virtual void eval(monty::rc_ptr<::mosek::fusion::WorkStack> rs,
                    monty::rc_ptr<::mosek::fusion::WorkStack> ws,
                    monty::rc_ptr<::mosek::fusion::WorkStack> xs);
  virtual void remove();
  virtual monty::rc_ptr<::mosek::fusion::Variable>
  __mosek_2fusion_2BaseVariable__fromTril(int32_t dim0, int32_t d);
  monty::rc_ptr<::mosek::fusion::Variable> fromTril(int32_t dim0, int32_t d);
  virtual monty::rc_ptr<::mosek::fusion::Variable>
  __mosek_2fusion_2BaseVariable__fromTril(int32_t d);
  monty::rc_ptr<::mosek::fusion::Variable> fromTril(int32_t d);
  /* override: mosek.fusion.Variable.fromTril*/
  virtual monty::rc_ptr<::mosek::fusion::Variable>
  __mosek_2fusion_2Variable__fromTril(int32_t d);
  virtual monty::rc_ptr<::mosek::fusion::Variable>
  __mosek_2fusion_2BaseVariable__tril(int32_t dim1, int32_t dim2);
  monty::rc_ptr<::mosek::fusion::Variable> tril(int32_t dim1, int32_t dim2);
  virtual monty::rc_ptr<::mosek::fusion::Variable>
  __mosek_2fusion_2BaseVariable__tril();
  monty::rc_ptr<::mosek::fusion::Variable> tril();
  /* override: mosek.fusion.Variable.tril*/
  virtual monty::rc_ptr<::mosek::fusion::Variable>
  __mosek_2fusion_2Variable__tril();
  virtual monty::rc_ptr<::mosek::fusion::Variable>
  __mosek_2fusion_2BaseVariable__reshape(int32_t dim0, int32_t dim1,
                                         int32_t dim2);
  monty::rc_ptr<::mosek::fusion::Variable> reshape(int32_t dim0, int32_t dim1,
                                                   int32_t dim2);
  /* override: mosek.fusion.Variable.reshape*/
  virtual monty::rc_ptr<::mosek::fusion::Variable>
  __mosek_2fusion_2Variable__reshape(int32_t dim0, int32_t dim1, int32_t dim2);
  virtual monty::rc_ptr<::mosek::fusion::Variable>
  __mosek_2fusion_2BaseVariable__reshape(int32_t dim0, int32_t dim1);
  monty::rc_ptr<::mosek::fusion::Variable> reshape(int32_t dim0, int32_t dim1);
  /* override: mosek.fusion.Variable.reshape*/
  virtual monty::rc_ptr<::mosek::fusion::Variable>
  __mosek_2fusion_2Variable__reshape(int32_t dim0, int32_t dim1);
  virtual monty::rc_ptr<::mosek::fusion::Variable>
  __mosek_2fusion_2BaseVariable__reshape(int32_t dim0);
  monty::rc_ptr<::mosek::fusion::Variable> reshape(int32_t dim0);
  /* override: mosek.fusion.Variable.reshape*/
  virtual monty::rc_ptr<::mosek::fusion::Variable>
  __mosek_2fusion_2Variable__reshape(int32_t dim0);
  virtual monty::rc_ptr<::mosek::fusion::Variable>
  __mosek_2fusion_2BaseVariable__reshape(
      std::shared_ptr<monty::ndarray<int32_t, 1>> shape);
  monty::rc_ptr<::mosek::fusion::Variable>
  reshape(std::shared_ptr<monty::ndarray<int32_t, 1>> shape);
  /* override: mosek.fusion.Variable.reshape*/
  virtual monty::rc_ptr<::mosek::fusion::Variable>
  __mosek_2fusion_2Variable__reshape(
      std::shared_ptr<monty::ndarray<int32_t, 1>> shape);
  virtual void setLevel(std::shared_ptr<monty::ndarray<double, 1>> v);
  virtual monty::rc_ptr<::mosek::fusion::Model>
  __mosek_2fusion_2BaseVariable__getModel();
  monty::rc_ptr<::mosek::fusion::Model> getModel();
  /* override: mosek.fusion.Variable.getModel*/
  virtual monty::rc_ptr<::mosek::fusion::Model>
  __mosek_2fusion_2Variable__getModel();
  virtual int32_t getND();
  virtual int32_t getDim(int32_t d);
  virtual std::shared_ptr<monty::ndarray<int32_t, 1>> getShape();
  virtual int64_t getSize();
  virtual std::shared_ptr<monty::ndarray<double, 1>> dual();
  virtual std::shared_ptr<monty::ndarray<double, 1>> level();
  virtual void makeContinuous();
  virtual void makeInteger();
  virtual monty::rc_ptr<::mosek::fusion::Variable>
  __mosek_2fusion_2BaseVariable__transpose();
  monty::rc_ptr<::mosek::fusion::Variable> transpose();
  /* override: mosek.fusion.Variable.transpose*/
  virtual monty::rc_ptr<::mosek::fusion::Variable>
  __mosek_2fusion_2Variable__transpose();
  virtual monty::rc_ptr<::mosek::fusion::Variable>
  __mosek_2fusion_2BaseVariable__index(int32_t i0, int32_t i1, int32_t i2);
  monty::rc_ptr<::mosek::fusion::Variable> index(int32_t i0, int32_t i1,
                                                 int32_t i2);
  /* override: mosek.fusion.Variable.index*/
  virtual monty::rc_ptr<::mosek::fusion::Variable>
  __mosek_2fusion_2Variable__index(int32_t i0, int32_t i1, int32_t i2);
  virtual monty::rc_ptr<::mosek::fusion::Variable>
  __mosek_2fusion_2BaseVariable__index(int32_t i0, int32_t i1);
  monty::rc_ptr<::mosek::fusion::Variable> index(int32_t i0, int32_t i1);
  /* override: mosek.fusion.Variable.index*/
  virtual monty::rc_ptr<::mosek::fusion::Variable>
  __mosek_2fusion_2Variable__index(int32_t i0, int32_t i1);
  virtual monty::rc_ptr<::mosek::fusion::Variable>
  __mosek_2fusion_2BaseVariable__index(
      std::shared_ptr<monty::ndarray<int32_t, 1>> index);
  monty::rc_ptr<::mosek::fusion::Variable>
  index(std::shared_ptr<monty::ndarray<int32_t, 1>> index);
  /* override: mosek.fusion.Variable.index*/
  virtual monty::rc_ptr<::mosek::fusion::Variable>
  __mosek_2fusion_2Variable__index(
      std::shared_ptr<monty::ndarray<int32_t, 1>> index);
  virtual monty::rc_ptr<::mosek::fusion::Variable>
  __mosek_2fusion_2BaseVariable__index(int32_t index);
  monty::rc_ptr<::mosek::fusion::Variable> index(int32_t index);
  /* override: mosek.fusion.Variable.index*/
  virtual monty::rc_ptr<::mosek::fusion::Variable>
  __mosek_2fusion_2Variable__index(int32_t index);
  virtual monty::rc_ptr<::mosek::fusion::Variable>
  __mosek_2fusion_2BaseVariable__pick(
      std::shared_ptr<monty::ndarray<int32_t, 1>> i0,
      std::shared_ptr<monty::ndarray<int32_t, 1>> i1,
      std::shared_ptr<monty::ndarray<int32_t, 1>> i2);
  monty::rc_ptr<::mosek::fusion::Variable>
  pick(std::shared_ptr<monty::ndarray<int32_t, 1>> i0,
       std::shared_ptr<monty::ndarray<int32_t, 1>> i1,
       std::shared_ptr<monty::ndarray<int32_t, 1>> i2);
  /* override: mosek.fusion.Variable.pick*/
  virtual monty::rc_ptr<::mosek::fusion::Variable>
  __mosek_2fusion_2Variable__pick(
      std::shared_ptr<monty::ndarray<int32_t, 1>> i0,
      std::shared_ptr<monty::ndarray<int32_t, 1>> i1,
      std::shared_ptr<monty::ndarray<int32_t, 1>> i2);
  virtual monty::rc_ptr<::mosek::fusion::Variable>
  __mosek_2fusion_2BaseVariable__pick(
      std::shared_ptr<monty::ndarray<int32_t, 1>> i0,
      std::shared_ptr<monty::ndarray<int32_t, 1>> i1);
  monty::rc_ptr<::mosek::fusion::Variable>
  pick(std::shared_ptr<monty::ndarray<int32_t, 1>> i0,
       std::shared_ptr<monty::ndarray<int32_t, 1>> i1);
  /* override: mosek.fusion.Variable.pick*/
  virtual monty::rc_ptr<::mosek::fusion::Variable>
  __mosek_2fusion_2Variable__pick(
      std::shared_ptr<monty::ndarray<int32_t, 1>> i0,
      std::shared_ptr<monty::ndarray<int32_t, 1>> i1);
  virtual monty::rc_ptr<::mosek::fusion::Variable>
  __mosek_2fusion_2BaseVariable__pick(
      std::shared_ptr<monty::ndarray<int32_t, 2>> midxs);
  monty::rc_ptr<::mosek::fusion::Variable>
  pick(std::shared_ptr<monty::ndarray<int32_t, 2>> midxs);
  /* override: mosek.fusion.Variable.pick*/
  virtual monty::rc_ptr<::mosek::fusion::Variable>
  __mosek_2fusion_2Variable__pick(
      std::shared_ptr<monty::ndarray<int32_t, 2>> midxs);
  virtual monty::rc_ptr<::mosek::fusion::Variable>
  __mosek_2fusion_2BaseVariable__pick(
      std::shared_ptr<monty::ndarray<int32_t, 1>> idxs);
  monty::rc_ptr<::mosek::fusion::Variable>
  pick(std::shared_ptr<monty::ndarray<int32_t, 1>> idxs);
  /* override: mosek.fusion.Variable.pick*/
  virtual monty::rc_ptr<::mosek::fusion::Variable>
  __mosek_2fusion_2Variable__pick(
      std::shared_ptr<monty::ndarray<int32_t, 1>> idxs);
  virtual monty::rc_ptr<::mosek::fusion::Variable>
  __mosek_2fusion_2BaseVariable__antidiag(int32_t index);
  monty::rc_ptr<::mosek::fusion::Variable> antidiag(int32_t index);
  /* override: mosek.fusion.Variable.antidiag*/
  virtual monty::rc_ptr<::mosek::fusion::Variable>
  __mosek_2fusion_2Variable__antidiag(int32_t index);
  virtual monty::rc_ptr<::mosek::fusion::Variable>
  __mosek_2fusion_2BaseVariable__antidiag();
  monty::rc_ptr<::mosek::fusion::Variable> antidiag();
  /* override: mosek.fusion.Variable.antidiag*/
  virtual monty::rc_ptr<::mosek::fusion::Variable>
  __mosek_2fusion_2Variable__antidiag();
  virtual monty::rc_ptr<::mosek::fusion::Variable>
  __mosek_2fusion_2BaseVariable__diag(int32_t index);
  monty::rc_ptr<::mosek::fusion::Variable> diag(int32_t index);
  /* override: mosek.fusion.Variable.diag*/
  virtual monty::rc_ptr<::mosek::fusion::Variable>
  __mosek_2fusion_2Variable__diag(int32_t index);
  virtual monty::rc_ptr<::mosek::fusion::Variable>
  __mosek_2fusion_2BaseVariable__diag();
  monty::rc_ptr<::mosek::fusion::Variable> diag();
  /* override: mosek.fusion.Variable.diag*/
  virtual monty::rc_ptr<::mosek::fusion::Variable>
  __mosek_2fusion_2Variable__diag();
  virtual monty::rc_ptr<::mosek::fusion::Variable>
  __mosek_2fusion_2BaseVariable__slice(
      std::shared_ptr<monty::ndarray<int32_t, 1>> first,
      std::shared_ptr<monty::ndarray<int32_t, 1>> last);
  monty::rc_ptr<::mosek::fusion::Variable>
  slice(std::shared_ptr<monty::ndarray<int32_t, 1>> first,
        std::shared_ptr<monty::ndarray<int32_t, 1>> last);
  /* override: mosek.fusion.Variable.slice*/
  virtual monty::rc_ptr<::mosek::fusion::Variable>
  __mosek_2fusion_2Variable__slice(
      std::shared_ptr<monty::ndarray<int32_t, 1>> first,
      std::shared_ptr<monty::ndarray<int32_t, 1>> last);
  virtual monty::rc_ptr<::mosek::fusion::Variable>
  __mosek_2fusion_2BaseVariable__slice(int32_t first, int32_t last);
  monty::rc_ptr<::mosek::fusion::Variable> slice(int32_t first, int32_t last);
  /* override: mosek.fusion.Variable.slice*/
  virtual monty::rc_ptr<::mosek::fusion::Variable>
  __mosek_2fusion_2Variable__slice(int32_t first, int32_t last);
  virtual monty::rc_ptr<::mosek::fusion::Expression>
  __mosek_2fusion_2BaseVariable__asExpr();
  monty::rc_ptr<::mosek::fusion::Expression> asExpr();
  /* override: mosek.fusion.Variable.asExpr*/
  virtual monty::rc_ptr<::mosek::fusion::Expression>
  __mosek_2fusion_2Variable__asExpr();
  virtual int32_t
  inst(int32_t spoffset, std::shared_ptr<monty::ndarray<int64_t, 1>> sparsity,
       int32_t nioffset,
       std::shared_ptr<monty::ndarray<int64_t, 1>> basevar_nativeidxs);
  virtual int32_t numInst();
  virtual void inst(int32_t offset,
                    std::shared_ptr<monty::ndarray<int64_t, 1>> nindex);
  virtual void set_values(std::shared_ptr<monty::ndarray<double, 1>> values,
                          bool primal);
  virtual void values(int32_t offset,
                      std::shared_ptr<monty::ndarray<double, 1>> target,
                      bool primal);
  virtual void make_continuous();
  virtual void make_integer();
}; // class BaseVariable;

class SliceVariable : public ::mosek::fusion::BaseVariable {
  SliceVariable(monty::rc_ptr<::mosek::fusion::Model> m,
                std::shared_ptr<monty::ndarray<int32_t, 1>> shape,
                std::shared_ptr<monty::ndarray<int64_t, 1>> sparsity,
                std::shared_ptr<monty::ndarray<int64_t, 1>> nativeidxs);
  SliceVariable(monty::rc_ptr<::mosek::fusion::SliceVariable> v);

protected:
  SliceVariable(p_SliceVariable *_impl);

public:
  SliceVariable(const SliceVariable &) = delete;
  const SliceVariable &operator=(const SliceVariable &) = delete;
  friend class p_SliceVariable;
  virtual ~SliceVariable();
  virtual void destroy();
  typedef monty::rc_ptr<SliceVariable> t;

}; // class SliceVariable;

class BoundInterfaceVariable : public ::mosek::fusion::SliceVariable {
  BoundInterfaceVariable(monty::rc_ptr<::mosek::fusion::Model> m,
                         std::shared_ptr<monty::ndarray<int32_t, 1>> shape,
                         std::shared_ptr<monty::ndarray<int64_t, 1>> sparsity,
                         std::shared_ptr<monty::ndarray<int64_t, 1>> nativeidxs,
                         bool islower);
  BoundInterfaceVariable(monty::rc_ptr<::mosek::fusion::SliceVariable> v,
                         bool islower);

protected:
  BoundInterfaceVariable(p_BoundInterfaceVariable *_impl);

public:
  BoundInterfaceVariable(const BoundInterfaceVariable &) = delete;
  const BoundInterfaceVariable &
  operator=(const BoundInterfaceVariable &) = delete;
  friend class p_BoundInterfaceVariable;
  virtual ~BoundInterfaceVariable();
  virtual void destroy();
  typedef monty::rc_ptr<BoundInterfaceVariable> t;

  virtual /* override */ std::shared_ptr<monty::ndarray<double, 1>> dual();
  virtual /* override */ monty::rc_ptr<::mosek::fusion::Variable>
  __mosek_2fusion_2BoundInterfaceVariable__transpose();
  monty::rc_ptr<::mosek::fusion::Variable> transpose();
  /* override: mosek.fusion.BaseVariable.transpose*/
  virtual monty::rc_ptr<::mosek::fusion::Variable>
  __mosek_2fusion_2BaseVariable__transpose();
  virtual /* override */ monty::rc_ptr<::mosek::fusion::Variable>
  __mosek_2fusion_2BoundInterfaceVariable__pick(
      std::shared_ptr<monty::ndarray<int32_t, 1>> i0,
      std::shared_ptr<monty::ndarray<int32_t, 1>> i1,
      std::shared_ptr<monty::ndarray<int32_t, 1>> i2);
  monty::rc_ptr<::mosek::fusion::Variable>
  pick(std::shared_ptr<monty::ndarray<int32_t, 1>> i0,
       std::shared_ptr<monty::ndarray<int32_t, 1>> i1,
       std::shared_ptr<monty::ndarray<int32_t, 1>> i2);
  /* override: mosek.fusion.BaseVariable.pick*/
  virtual monty::rc_ptr<::mosek::fusion::Variable>
  __mosek_2fusion_2BaseVariable__pick(
      std::shared_ptr<monty::ndarray<int32_t, 1>> i0,
      std::shared_ptr<monty::ndarray<int32_t, 1>> i1,
      std::shared_ptr<monty::ndarray<int32_t, 1>> i2);
  virtual /* override */ monty::rc_ptr<::mosek::fusion::Variable>
  __mosek_2fusion_2BoundInterfaceVariable__pick(
      std::shared_ptr<monty::ndarray<int32_t, 1>> i0,
      std::shared_ptr<monty::ndarray<int32_t, 1>> i1);
  monty::rc_ptr<::mosek::fusion::Variable>
  pick(std::shared_ptr<monty::ndarray<int32_t, 1>> i0,
       std::shared_ptr<monty::ndarray<int32_t, 1>> i1);
  /* override: mosek.fusion.BaseVariable.pick*/
  virtual monty::rc_ptr<::mosek::fusion::Variable>
  __mosek_2fusion_2BaseVariable__pick(
      std::shared_ptr<monty::ndarray<int32_t, 1>> i0,
      std::shared_ptr<monty::ndarray<int32_t, 1>> i1);
  virtual /* override */ monty::rc_ptr<::mosek::fusion::Variable>
  __mosek_2fusion_2BoundInterfaceVariable__pick(
      std::shared_ptr<monty::ndarray<int32_t, 2>> midxs);
  monty::rc_ptr<::mosek::fusion::Variable>
  pick(std::shared_ptr<monty::ndarray<int32_t, 2>> midxs);
  /* override: mosek.fusion.BaseVariable.pick*/
  virtual monty::rc_ptr<::mosek::fusion::Variable>
  __mosek_2fusion_2BaseVariable__pick(
      std::shared_ptr<monty::ndarray<int32_t, 2>> midxs);
  virtual /* override */ monty::rc_ptr<::mosek::fusion::Variable>
  __mosek_2fusion_2BoundInterfaceVariable__pick(
      std::shared_ptr<monty::ndarray<int32_t, 1>> idxs);
  monty::rc_ptr<::mosek::fusion::Variable>
  pick(std::shared_ptr<monty::ndarray<int32_t, 1>> idxs);
  /* override: mosek.fusion.BaseVariable.pick*/
  virtual monty::rc_ptr<::mosek::fusion::Variable>
  __mosek_2fusion_2BaseVariable__pick(
      std::shared_ptr<monty::ndarray<int32_t, 1>> idxs);
  virtual /* override */ monty::rc_ptr<::mosek::fusion::Variable>
  __mosek_2fusion_2BoundInterfaceVariable__antidiag(int32_t index);
  monty::rc_ptr<::mosek::fusion::Variable> antidiag(int32_t index);
  /* override: mosek.fusion.BaseVariable.antidiag*/
  virtual monty::rc_ptr<::mosek::fusion::Variable>
  __mosek_2fusion_2BaseVariable__antidiag(int32_t index);
  virtual /* override */ monty::rc_ptr<::mosek::fusion::Variable>
  __mosek_2fusion_2BoundInterfaceVariable__antidiag();
  monty::rc_ptr<::mosek::fusion::Variable> antidiag();
  /* override: mosek.fusion.BaseVariable.antidiag*/
  virtual monty::rc_ptr<::mosek::fusion::Variable>
  __mosek_2fusion_2BaseVariable__antidiag();
  virtual /* override */ monty::rc_ptr<::mosek::fusion::Variable>
  __mosek_2fusion_2BoundInterfaceVariable__diag(int32_t index);
  monty::rc_ptr<::mosek::fusion::Variable> diag(int32_t index);
  /* override: mosek.fusion.BaseVariable.diag*/
  virtual monty::rc_ptr<::mosek::fusion::Variable>
  __mosek_2fusion_2BaseVariable__diag(int32_t index);
  virtual /* override */ monty::rc_ptr<::mosek::fusion::Variable>
  __mosek_2fusion_2BoundInterfaceVariable__diag();
  monty::rc_ptr<::mosek::fusion::Variable> diag();
  /* override: mosek.fusion.BaseVariable.diag*/
  virtual monty::rc_ptr<::mosek::fusion::Variable>
  __mosek_2fusion_2BaseVariable__diag();
  virtual /* override */ monty::rc_ptr<::mosek::fusion::Variable>
  __mosek_2fusion_2BoundInterfaceVariable__slice(
      std::shared_ptr<monty::ndarray<int32_t, 1>> firsta,
      std::shared_ptr<monty::ndarray<int32_t, 1>> lasta);
  monty::rc_ptr<::mosek::fusion::Variable>
  slice(std::shared_ptr<monty::ndarray<int32_t, 1>> firsta,
        std::shared_ptr<monty::ndarray<int32_t, 1>> lasta);
  /* override: mosek.fusion.BaseVariable.slice*/
  virtual monty::rc_ptr<::mosek::fusion::Variable>
  __mosek_2fusion_2BaseVariable__slice(
      std::shared_ptr<monty::ndarray<int32_t, 1>> firsta,
      std::shared_ptr<monty::ndarray<int32_t, 1>> lasta);
  virtual /* override */ monty::rc_ptr<::mosek::fusion::Variable>
  __mosek_2fusion_2BoundInterfaceVariable__slice(int32_t first, int32_t last);
  monty::rc_ptr<::mosek::fusion::Variable> slice(int32_t first, int32_t last);
  /* override: mosek.fusion.BaseVariable.slice*/
  virtual monty::rc_ptr<::mosek::fusion::Variable>
  __mosek_2fusion_2BaseVariable__slice(int32_t first, int32_t last);
}; // class BoundInterfaceVariable;

class ModelVariable : public ::mosek::fusion::BaseVariable {
protected:
  ModelVariable(p_ModelVariable *_impl);

public:
  ModelVariable(const ModelVariable &) = delete;
  const ModelVariable &operator=(const ModelVariable &) = delete;
  friend class p_ModelVariable;
  virtual ~ModelVariable();
  virtual void destroy();
  typedef monty::rc_ptr<ModelVariable> t;

  virtual void
  elementName(int64_t index,
              monty::rc_ptr<::mosek::fusion::Utils::StringBuffer> sb);
  virtual /* override */ void remove();
}; // class ModelVariable;

class RangedVariable : public ::mosek::fusion::ModelVariable {
  RangedVariable(monty::rc_ptr<::mosek::fusion::RangedVariable> v,
                 monty::rc_ptr<::mosek::fusion::Model> m);
  RangedVariable(monty::rc_ptr<::mosek::fusion::Model> model,
                 const std::string &name, int64_t varid,
                 std::shared_ptr<monty::ndarray<int32_t, 1>> shape,
                 std::shared_ptr<monty::ndarray<int64_t, 1>> sparsity,
                 std::shared_ptr<monty::ndarray<int32_t, 1>> nativeidxs);

protected:
  RangedVariable(p_RangedVariable *_impl);

public:
  RangedVariable(const RangedVariable &) = delete;
  const RangedVariable &operator=(const RangedVariable &) = delete;
  friend class p_RangedVariable;
  virtual ~RangedVariable();
  virtual void destroy();
  typedef monty::rc_ptr<RangedVariable> t;

  virtual monty::rc_ptr<::mosek::fusion::Utils::StringBuffer>
  __mosek_2fusion_2RangedVariable__elementDesc(
      int64_t index, monty::rc_ptr<::mosek::fusion::Utils::StringBuffer> sb);
  monty::rc_ptr<::mosek::fusion::Utils::StringBuffer>
  elementDesc(int64_t index,
              monty::rc_ptr<::mosek::fusion::Utils::StringBuffer> sb);
  virtual monty::rc_ptr<::mosek::fusion::BoundInterfaceVariable>
  __mosek_2fusion_2RangedVariable__upperBoundVar();
  monty::rc_ptr<::mosek::fusion::BoundInterfaceVariable> upperBoundVar();
  virtual monty::rc_ptr<::mosek::fusion::BoundInterfaceVariable>
  __mosek_2fusion_2RangedVariable__lowerBoundVar();
  monty::rc_ptr<::mosek::fusion::BoundInterfaceVariable> lowerBoundVar();
}; // class RangedVariable;

class LinearPSDVariable : public ::mosek::fusion::ModelVariable {
  LinearPSDVariable(monty::rc_ptr<::mosek::fusion::LinearPSDVariable> v,
                    monty::rc_ptr<::mosek::fusion::Model> m);
  LinearPSDVariable(monty::rc_ptr<::mosek::fusion::Model> model,
                    const std::string &name, int32_t varid,
                    std::shared_ptr<monty::ndarray<int32_t, 1>> shape,
                    int32_t conedim,
                    std::shared_ptr<monty::ndarray<int64_t, 1>> nativeidxs);

protected:
  LinearPSDVariable(p_LinearPSDVariable *_impl);

public:
  LinearPSDVariable(const LinearPSDVariable &) = delete;
  const LinearPSDVariable &operator=(const LinearPSDVariable &) = delete;
  friend class p_LinearPSDVariable;
  virtual ~LinearPSDVariable();
  virtual void destroy();
  typedef monty::rc_ptr<LinearPSDVariable> t;

  virtual /* override */ std::string toString();
  virtual void
  make_continuous(std::shared_ptr<monty::ndarray<int64_t, 1>> idxs);
  virtual void make_integer(std::shared_ptr<monty::ndarray<int64_t, 1>> idxs);
}; // class LinearPSDVariable;

class PSDVariable : public ::mosek::fusion::ModelVariable {
  PSDVariable(monty::rc_ptr<::mosek::fusion::PSDVariable> v,
              monty::rc_ptr<::mosek::fusion::Model> m);
  PSDVariable(monty::rc_ptr<::mosek::fusion::Model> model,
              const std::string &name, int32_t varid,
              std::shared_ptr<monty::ndarray<int32_t, 1>> shape,
              int32_t conedim1, int32_t conedim2,
              std::shared_ptr<monty::ndarray<int32_t, 1>> barvaridxs,
              std::shared_ptr<monty::ndarray<int64_t, 1>> nativeidxs);

protected:
  PSDVariable(p_PSDVariable *_impl);

public:
  PSDVariable(const PSDVariable &) = delete;
  const PSDVariable &operator=(const PSDVariable &) = delete;
  friend class p_PSDVariable;
  virtual ~PSDVariable();
  virtual void destroy();
  typedef monty::rc_ptr<PSDVariable> t;

  virtual /* override */ std::string toString();
  virtual monty::rc_ptr<::mosek::fusion::Utils::StringBuffer>
  __mosek_2fusion_2PSDVariable__elementDesc(
      int64_t index, monty::rc_ptr<::mosek::fusion::Utils::StringBuffer> sb);
  monty::rc_ptr<::mosek::fusion::Utils::StringBuffer>
  elementDesc(int64_t index,
              monty::rc_ptr<::mosek::fusion::Utils::StringBuffer> sb);
}; // class PSDVariable;

class LinearVariable : public ::mosek::fusion::ModelVariable {
  LinearVariable(monty::rc_ptr<::mosek::fusion::LinearVariable> v,
                 monty::rc_ptr<::mosek::fusion::Model> m);
  LinearVariable(monty::rc_ptr<::mosek::fusion::Model> model,
                 const std::string &name, int64_t varid,
                 std::shared_ptr<monty::ndarray<int32_t, 1>> shape,
                 std::shared_ptr<monty::ndarray<int64_t, 1>> sparsity,
                 std::shared_ptr<monty::ndarray<int32_t, 1>> nativeidxs);

protected:
  LinearVariable(p_LinearVariable *_impl);

public:
  LinearVariable(const LinearVariable &) = delete;
  const LinearVariable &operator=(const LinearVariable &) = delete;
  friend class p_LinearVariable;
  virtual ~LinearVariable();
  virtual void destroy();
  typedef monty::rc_ptr<LinearVariable> t;

  virtual /* override */ std::string toString();
}; // class LinearVariable;

class ConicVariable : public ::mosek::fusion::ModelVariable {
  ConicVariable(monty::rc_ptr<::mosek::fusion::ConicVariable> v,
                monty::rc_ptr<::mosek::fusion::Model> m);
  ConicVariable(monty::rc_ptr<::mosek::fusion::Model> model,
                const std::string &name, int32_t varid,
                std::shared_ptr<monty::ndarray<int32_t, 1>> shape,
                std::shared_ptr<monty::ndarray<int32_t, 1>> nativeidxs);

protected:
  ConicVariable(p_ConicVariable *_impl);

public:
  ConicVariable(const ConicVariable &) = delete;
  const ConicVariable &operator=(const ConicVariable &) = delete;
  friend class p_ConicVariable;
  virtual ~ConicVariable();
  virtual void destroy();
  typedef monty::rc_ptr<ConicVariable> t;

  virtual /* override */ std::string toString();
}; // class ConicVariable;

class NilVariable : public ::mosek::fusion::BaseVariable {
  NilVariable(std::shared_ptr<monty::ndarray<int32_t, 1>> shape);
  NilVariable();

protected:
  NilVariable(p_NilVariable *_impl);

public:
  NilVariable(const NilVariable &) = delete;
  const NilVariable &operator=(const NilVariable &) = delete;
  friend class p_NilVariable;
  virtual ~NilVariable();
  virtual void destroy();
  typedef monty::rc_ptr<NilVariable> t;

  virtual monty::rc_ptr<::mosek::fusion::Utils::StringBuffer>
  __mosek_2fusion_2NilVariable__elementDesc(
      int64_t index, monty::rc_ptr<::mosek::fusion::Utils::StringBuffer> sb);
  monty::rc_ptr<::mosek::fusion::Utils::StringBuffer>
  elementDesc(int64_t index,
              monty::rc_ptr<::mosek::fusion::Utils::StringBuffer> sb);
  virtual void
  elementName(int64_t index,
              monty::rc_ptr<::mosek::fusion::Utils::StringBuffer> sb);
  virtual /* override */ int32_t numInst();
  virtual int32_t
  inst(int32_t offset, std::shared_ptr<monty::ndarray<int64_t, 1>> sparsity,
       std::shared_ptr<monty::ndarray<int64_t, 1>> basevar_nativeidxs);
  virtual /* override */ void
  inst(int32_t offset, std::shared_ptr<monty::ndarray<int64_t, 1>> nindex);
  virtual /* override */ void
  set_values(std::shared_ptr<monty::ndarray<double, 1>> target, bool primal);
  virtual /* override */ void
  values(int32_t offset, std::shared_ptr<monty::ndarray<double, 1>> target,
         bool primal);
  virtual /* override */ void make_continuous();
  virtual /* override */ void make_integer();
  virtual /* override */ std::string toString();
  virtual /* override */ monty::rc_ptr<::mosek::fusion::Variable>
  __mosek_2fusion_2NilVariable__index(
      std::shared_ptr<monty::ndarray<int32_t, 1>> first);
  monty::rc_ptr<::mosek::fusion::Variable>
  index(std::shared_ptr<monty::ndarray<int32_t, 1>> first);
  /* override: mosek.fusion.BaseVariable.index*/
  virtual monty::rc_ptr<::mosek::fusion::Variable>
  __mosek_2fusion_2BaseVariable__index(
      std::shared_ptr<monty::ndarray<int32_t, 1>> first);
  virtual /* override */ monty::rc_ptr<::mosek::fusion::Variable>
  __mosek_2fusion_2NilVariable__index(int32_t first);
  monty::rc_ptr<::mosek::fusion::Variable> index(int32_t first);
  /* override: mosek.fusion.BaseVariable.index*/
  virtual monty::rc_ptr<::mosek::fusion::Variable>
  __mosek_2fusion_2BaseVariable__index(int32_t first);
  virtual /* override */ monty::rc_ptr<::mosek::fusion::Variable>
  __mosek_2fusion_2NilVariable__slice(
      std::shared_ptr<monty::ndarray<int32_t, 1>> first,
      std::shared_ptr<monty::ndarray<int32_t, 1>> last);
  monty::rc_ptr<::mosek::fusion::Variable>
  slice(std::shared_ptr<monty::ndarray<int32_t, 1>> first,
        std::shared_ptr<monty::ndarray<int32_t, 1>> last);
  /* override: mosek.fusion.BaseVariable.slice*/
  virtual monty::rc_ptr<::mosek::fusion::Variable>
  __mosek_2fusion_2BaseVariable__slice(
      std::shared_ptr<monty::ndarray<int32_t, 1>> first,
      std::shared_ptr<monty::ndarray<int32_t, 1>> last);
  virtual /* override */ monty::rc_ptr<::mosek::fusion::Variable>
  __mosek_2fusion_2NilVariable__slice(int32_t first, int32_t last);
  monty::rc_ptr<::mosek::fusion::Variable> slice(int32_t first, int32_t last);
  /* override: mosek.fusion.BaseVariable.slice*/
  virtual monty::rc_ptr<::mosek::fusion::Variable>
  __mosek_2fusion_2BaseVariable__slice(int32_t first, int32_t last);
}; // class NilVariable;

class Var : public virtual monty::RefCounted {
public:
  p_Var *_impl;

protected:
  Var(p_Var *_impl);

public:
  Var(const Var &) = delete;
  const Var &operator=(const Var &) = delete;
  friend class p_Var;
  virtual ~Var();
  virtual void destroy();
  typedef monty::rc_ptr<Var> t;

  static monty::rc_ptr<::mosek::fusion::Variable>
  empty(std::shared_ptr<monty::ndarray<int32_t, 1>> shape);
  static monty::rc_ptr<::mosek::fusion::Variable>
  compress(monty::rc_ptr<::mosek::fusion::Variable> v);
  static monty::rc_ptr<::mosek::fusion::Variable>
  reshape(monty::rc_ptr<::mosek::fusion::Variable> v, int32_t d1);
  static monty::rc_ptr<::mosek::fusion::Variable>
  reshape(monty::rc_ptr<::mosek::fusion::Variable> v, int32_t d1, int32_t d2);
  static monty::rc_ptr<::mosek::fusion::Variable>
  flatten(monty::rc_ptr<::mosek::fusion::Variable> v);
  static monty::rc_ptr<::mosek::fusion::Variable>
  reshape(monty::rc_ptr<::mosek::fusion::Variable> v,
          std::shared_ptr<monty::ndarray<int32_t, 1>> shape);
  static monty::rc_ptr<::mosek::fusion::Variable>
  hrepeat(monty::rc_ptr<::mosek::fusion::Variable> v, int32_t n);
  static monty::rc_ptr<::mosek::fusion::Variable>
  vrepeat(monty::rc_ptr<::mosek::fusion::Variable> v, int32_t n);
  static monty::rc_ptr<::mosek::fusion::Variable>
  repeat(monty::rc_ptr<::mosek::fusion::Variable> v, int32_t n);
  static monty::rc_ptr<::mosek::fusion::Variable>
  repeat(monty::rc_ptr<::mosek::fusion::Variable> v, int32_t dim, int32_t n);
  static monty::rc_ptr<::mosek::fusion::Variable>
  stack(std::shared_ptr<
        monty::ndarray<std::shared_ptr<monty::ndarray<
                           monty::rc_ptr<::mosek::fusion::Variable>, 1>>,
                       1>>
            vlist);
  static monty::rc_ptr<::mosek::fusion::Variable>
  vstack(monty::rc_ptr<::mosek::fusion::Variable> v1,
         monty::rc_ptr<::mosek::fusion::Variable> v2,
         monty::rc_ptr<::mosek::fusion::Variable> v3);
  static monty::rc_ptr<::mosek::fusion::Variable>
  vstack(monty::rc_ptr<::mosek::fusion::Variable> v1,
         monty::rc_ptr<::mosek::fusion::Variable> v2);
  static monty::rc_ptr<::mosek::fusion::Variable>
  vstack(std::shared_ptr<
         monty::ndarray<monty::rc_ptr<::mosek::fusion::Variable>, 1>>
             v);
  static monty::rc_ptr<::mosek::fusion::Variable>
  hstack(monty::rc_ptr<::mosek::fusion::Variable> v1,
         monty::rc_ptr<::mosek::fusion::Variable> v2,
         monty::rc_ptr<::mosek::fusion::Variable> v3);
  static monty::rc_ptr<::mosek::fusion::Variable>
  hstack(monty::rc_ptr<::mosek::fusion::Variable> v1,
         monty::rc_ptr<::mosek::fusion::Variable> v2);
  static monty::rc_ptr<::mosek::fusion::Variable>
  hstack(std::shared_ptr<
         monty::ndarray<monty::rc_ptr<::mosek::fusion::Variable>, 1>>
             v);
  static monty::rc_ptr<::mosek::fusion::Variable>
  stack(std::shared_ptr<
            monty::ndarray<monty::rc_ptr<::mosek::fusion::Variable>, 1>>
            v,
        int32_t dim);
  static monty::rc_ptr<::mosek::fusion::Variable>
  stack(monty::rc_ptr<::mosek::fusion::Variable> v1,
        monty::rc_ptr<::mosek::fusion::Variable> v2,
        monty::rc_ptr<::mosek::fusion::Variable> v3, int32_t dim);
  static monty::rc_ptr<::mosek::fusion::Variable>
  stack(monty::rc_ptr<::mosek::fusion::Variable> v1,
        monty::rc_ptr<::mosek::fusion::Variable> v2, int32_t dim);
  static monty::rc_ptr<::mosek::fusion::Variable>
  stack(int32_t dim,
        std::shared_ptr<
            monty::ndarray<monty::rc_ptr<::mosek::fusion::Variable>, 1>>
            v);
  static monty::rc_ptr<::mosek::fusion::Variable>
  stack(int32_t dim, monty::rc_ptr<::mosek::fusion::Variable> v1,
        monty::rc_ptr<::mosek::fusion::Variable> v2,
        monty::rc_ptr<::mosek::fusion::Variable> v3);
  static monty::rc_ptr<::mosek::fusion::Variable>
  stack(int32_t dim, monty::rc_ptr<::mosek::fusion::Variable> v1,
        monty::rc_ptr<::mosek::fusion::Variable> v2);
  static monty::rc_ptr<::mosek::fusion::Variable>
  promote(monty::rc_ptr<::mosek::fusion::Variable> v, int32_t nd);
}; // class Var;

class Constraint : public virtual monty::RefCounted {
public:
  p_Constraint *_impl;

protected:
  Constraint(p_Constraint *_impl);

public:
  Constraint(const Constraint &) = delete;
  const Constraint &operator=(const Constraint &) = delete;
  friend class p_Constraint;
  virtual ~Constraint();
  virtual void destroy();
  typedef monty::rc_ptr<Constraint> t;

  Constraint(monty::rc_ptr<::mosek::fusion::Constraint> c,
             monty::rc_ptr<::mosek::fusion::Model> m);
  Constraint(monty::rc_ptr<::mosek::fusion::Model> model,
             std::shared_ptr<monty::ndarray<int32_t, 1>> shape,
             std::shared_ptr<monty::ndarray<int32_t, 1>> con_nativeidxs);
  virtual /* override */ std::string toString();
  virtual std::shared_ptr<monty::ndarray<double, 1>> dual();
  virtual std::shared_ptr<monty::ndarray<double, 1>> level();
  virtual void remove();
  virtual void update(std::shared_ptr<monty::ndarray<double, 1>> bfix);
  virtual void update(monty::rc_ptr<::mosek::fusion::Expression> expr);
  virtual void update(monty::rc_ptr<::mosek::fusion::Expression> expr,
                      monty::rc_ptr<::mosek::fusion::Variable> x,
                      bool bfixupdate);
  virtual void update(monty::rc_ptr<::mosek::fusion::Expression> expr,
                      monty::rc_ptr<::mosek::fusion::Variable> x);
  virtual monty::rc_ptr<::mosek::fusion::Model>
  __mosek_2fusion_2Constraint__get_model();
  monty::rc_ptr<::mosek::fusion::Model> get_model();
  virtual int32_t get_nd();
  virtual int64_t size();
  static monty::rc_ptr<::mosek::fusion::Constraint>
  stack(std::shared_ptr<
            monty::ndarray<monty::rc_ptr<::mosek::fusion::Constraint>, 1>>
            clist,
        int32_t dim);
  static monty::rc_ptr<::mosek::fusion::Constraint>
  stack(monty::rc_ptr<::mosek::fusion::Constraint> v1,
        monty::rc_ptr<::mosek::fusion::Constraint> v2,
        monty::rc_ptr<::mosek::fusion::Constraint> v3, int32_t dim);
  static monty::rc_ptr<::mosek::fusion::Constraint>
  stack(monty::rc_ptr<::mosek::fusion::Constraint> v1,
        monty::rc_ptr<::mosek::fusion::Constraint> v2, int32_t dim);
  static monty::rc_ptr<::mosek::fusion::Constraint>
  hstack(std::shared_ptr<
         monty::ndarray<monty::rc_ptr<::mosek::fusion::Constraint>, 1>>
             clist);
  static monty::rc_ptr<::mosek::fusion::Constraint>
  vstack(std::shared_ptr<
         monty::ndarray<monty::rc_ptr<::mosek::fusion::Constraint>, 1>>
             clist);
  static monty::rc_ptr<::mosek::fusion::Constraint>
  hstack(monty::rc_ptr<::mosek::fusion::Constraint> v1,
         monty::rc_ptr<::mosek::fusion::Constraint> v2,
         monty::rc_ptr<::mosek::fusion::Constraint> v3);
  static monty::rc_ptr<::mosek::fusion::Constraint>
  vstack(monty::rc_ptr<::mosek::fusion::Constraint> v1,
         monty::rc_ptr<::mosek::fusion::Constraint> v2,
         monty::rc_ptr<::mosek::fusion::Constraint> v3);
  static monty::rc_ptr<::mosek::fusion::Constraint>
  hstack(monty::rc_ptr<::mosek::fusion::Constraint> v1,
         monty::rc_ptr<::mosek::fusion::Constraint> v2);
  static monty::rc_ptr<::mosek::fusion::Constraint>
  vstack(monty::rc_ptr<::mosek::fusion::Constraint> v1,
         monty::rc_ptr<::mosek::fusion::Constraint> v2);
  virtual monty::rc_ptr<::mosek::fusion::Constraint>
  __mosek_2fusion_2Constraint__index(
      std::shared_ptr<monty::ndarray<int32_t, 1>> idxa);
  monty::rc_ptr<::mosek::fusion::Constraint>
  index(std::shared_ptr<monty::ndarray<int32_t, 1>> idxa);
  virtual monty::rc_ptr<::mosek::fusion::Constraint>
  __mosek_2fusion_2Constraint__index(int32_t idx);
  monty::rc_ptr<::mosek::fusion::Constraint> index(int32_t idx);
  virtual monty::rc_ptr<::mosek::fusion::Constraint>
  __mosek_2fusion_2Constraint__slice(
      std::shared_ptr<monty::ndarray<int32_t, 1>> firsta,
      std::shared_ptr<monty::ndarray<int32_t, 1>> lasta);
  monty::rc_ptr<::mosek::fusion::Constraint>
  slice(std::shared_ptr<monty::ndarray<int32_t, 1>> firsta,
        std::shared_ptr<monty::ndarray<int32_t, 1>> lasta);
  virtual monty::rc_ptr<::mosek::fusion::Constraint>
  __mosek_2fusion_2Constraint__slice(int32_t first, int32_t last);
  monty::rc_ptr<::mosek::fusion::Constraint> slice(int32_t first, int32_t last);
  virtual int32_t getND();
  virtual int32_t getSize();
  virtual monty::rc_ptr<::mosek::fusion::Model>
  __mosek_2fusion_2Constraint__getModel();
  monty::rc_ptr<::mosek::fusion::Model> getModel();
  virtual std::shared_ptr<monty::ndarray<int32_t, 1>> getShape();
}; // class Constraint;

class SliceConstraint : public ::mosek::fusion::Constraint {
  SliceConstraint(monty::rc_ptr<::mosek::fusion::SliceConstraint> c);
  SliceConstraint(monty::rc_ptr<::mosek::fusion::Model> model,
                  std::shared_ptr<monty::ndarray<int32_t, 1>> shape,
                  std::shared_ptr<monty::ndarray<int32_t, 1>> nativeidxs);

protected:
  SliceConstraint(p_SliceConstraint *_impl);

public:
  SliceConstraint(const SliceConstraint &) = delete;
  const SliceConstraint &operator=(const SliceConstraint &) = delete;
  friend class p_SliceConstraint;
  virtual ~SliceConstraint();
  virtual void destroy();
  typedef monty::rc_ptr<SliceConstraint> t;

  virtual /* override */ std::string toString();
}; // class SliceConstraint;

class BoundInterfaceConstraint : public ::mosek::fusion::SliceConstraint {
  BoundInterfaceConstraint(
      monty::rc_ptr<::mosek::fusion::Model> m,
      std::shared_ptr<monty::ndarray<int32_t, 1>> shape,
      std::shared_ptr<monty::ndarray<int32_t, 1>> nativeidxs, bool islower);
  BoundInterfaceConstraint(monty::rc_ptr<::mosek::fusion::SliceConstraint> c,
                           bool islower);

protected:
  BoundInterfaceConstraint(p_BoundInterfaceConstraint *_impl);

public:
  BoundInterfaceConstraint(const BoundInterfaceConstraint &) = delete;
  const BoundInterfaceConstraint &
  operator=(const BoundInterfaceConstraint &) = delete;
  friend class p_BoundInterfaceConstraint;
  virtual ~BoundInterfaceConstraint();
  virtual void destroy();
  typedef monty::rc_ptr<BoundInterfaceConstraint> t;

  virtual /* override */ std::shared_ptr<monty::ndarray<double, 1>> dual();
  virtual /* override */ monty::rc_ptr<::mosek::fusion::Constraint>
  __mosek_2fusion_2BoundInterfaceConstraint__slice(
      std::shared_ptr<monty::ndarray<int32_t, 1>> firsta,
      std::shared_ptr<monty::ndarray<int32_t, 1>> lasta);
  monty::rc_ptr<::mosek::fusion::Constraint>
  slice(std::shared_ptr<monty::ndarray<int32_t, 1>> firsta,
        std::shared_ptr<monty::ndarray<int32_t, 1>> lasta);
  /* override: mosek.fusion.Constraint.slice*/
  virtual monty::rc_ptr<::mosek::fusion::Constraint>
  __mosek_2fusion_2Constraint__slice(
      std::shared_ptr<monty::ndarray<int32_t, 1>> firsta,
      std::shared_ptr<monty::ndarray<int32_t, 1>> lasta);
  virtual /* override */ monty::rc_ptr<::mosek::fusion::Constraint>
  __mosek_2fusion_2BoundInterfaceConstraint__slice(int32_t first, int32_t last);
  monty::rc_ptr<::mosek::fusion::Constraint> slice(int32_t first, int32_t last);
  /* override: mosek.fusion.Constraint.slice*/
  virtual monty::rc_ptr<::mosek::fusion::Constraint>
  __mosek_2fusion_2Constraint__slice(int32_t first, int32_t last);
  virtual /* override */ monty::rc_ptr<::mosek::fusion::Constraint>
  __mosek_2fusion_2BoundInterfaceConstraint__index(
      std::shared_ptr<monty::ndarray<int32_t, 1>> idxa);
  monty::rc_ptr<::mosek::fusion::Constraint>
  index(std::shared_ptr<monty::ndarray<int32_t, 1>> idxa);
  /* override: mosek.fusion.Constraint.index*/
  virtual monty::rc_ptr<::mosek::fusion::Constraint>
  __mosek_2fusion_2Constraint__index(
      std::shared_ptr<monty::ndarray<int32_t, 1>> idxa);
  virtual /* override */ monty::rc_ptr<::mosek::fusion::Constraint>
  __mosek_2fusion_2BoundInterfaceConstraint__index(int32_t idx);
  monty::rc_ptr<::mosek::fusion::Constraint> index(int32_t idx);
  /* override: mosek.fusion.Constraint.index*/
  virtual monty::rc_ptr<::mosek::fusion::Constraint>
  __mosek_2fusion_2Constraint__index(int32_t idx);
}; // class BoundInterfaceConstraint;

class ModelConstraint : public ::mosek::fusion::Constraint {
protected:
  ModelConstraint(p_ModelConstraint *_impl);

public:
  ModelConstraint(const ModelConstraint &) = delete;
  const ModelConstraint &operator=(const ModelConstraint &) = delete;
  friend class p_ModelConstraint;
  virtual ~ModelConstraint();
  virtual void destroy();
  typedef monty::rc_ptr<ModelConstraint> t;

  virtual /* override */ std::string toString();
  virtual /* override */ void remove();
}; // class ModelConstraint;

class LinearPSDConstraint : public ::mosek::fusion::ModelConstraint {
  LinearPSDConstraint(monty::rc_ptr<::mosek::fusion::LinearPSDConstraint> c,
                      monty::rc_ptr<::mosek::fusion::Model> m);
  LinearPSDConstraint(monty::rc_ptr<::mosek::fusion::Model> model,
                      const std::string &name, int32_t conid,
                      std::shared_ptr<monty::ndarray<int32_t, 1>> shape,
                      int32_t conedim,
                      std::shared_ptr<monty::ndarray<int32_t, 1>> nativeidxs,
                      std::shared_ptr<monty::ndarray<int64_t, 1>> slackidxs);

protected:
  LinearPSDConstraint(p_LinearPSDConstraint *_impl);

public:
  LinearPSDConstraint(const LinearPSDConstraint &) = delete;
  const LinearPSDConstraint &operator=(const LinearPSDConstraint &) = delete;
  friend class p_LinearPSDConstraint;
  virtual ~LinearPSDConstraint();
  virtual void destroy();
  typedef monty::rc_ptr<LinearPSDConstraint> t;

}; // class LinearPSDConstraint;

class PSDConstraint : public ::mosek::fusion::ModelConstraint {
  PSDConstraint(monty::rc_ptr<::mosek::fusion::PSDConstraint> c,
                monty::rc_ptr<::mosek::fusion::Model> m);
  PSDConstraint(monty::rc_ptr<::mosek::fusion::Model> model,
                const std::string &name, int32_t conid,
                std::shared_ptr<monty::ndarray<int32_t, 1>> shape,
                int32_t conedim0, int32_t conedim1,
                std::shared_ptr<monty::ndarray<int64_t, 1>> slackidxs,
                std::shared_ptr<monty::ndarray<int32_t, 1>> nativeidxs);

protected:
  PSDConstraint(p_PSDConstraint *_impl);

public:
  PSDConstraint(const PSDConstraint &) = delete;
  const PSDConstraint &operator=(const PSDConstraint &) = delete;
  friend class p_PSDConstraint;
  virtual ~PSDConstraint();
  virtual void destroy();
  typedef monty::rc_ptr<PSDConstraint> t;

  virtual /* override */ std::string toString();
}; // class PSDConstraint;

class RangedConstraint : public ::mosek::fusion::ModelConstraint {
  RangedConstraint(monty::rc_ptr<::mosek::fusion::RangedConstraint> c,
                   monty::rc_ptr<::mosek::fusion::Model> m);
  RangedConstraint(monty::rc_ptr<::mosek::fusion::Model> model,
                   const std::string &name,
                   std::shared_ptr<monty::ndarray<int32_t, 1>> shape,
                   std::shared_ptr<monty::ndarray<int32_t, 1>> nativeidxs,
                   int32_t conid);

protected:
  RangedConstraint(p_RangedConstraint *_impl);

public:
  RangedConstraint(const RangedConstraint &) = delete;
  const RangedConstraint &operator=(const RangedConstraint &) = delete;
  friend class p_RangedConstraint;
  virtual ~RangedConstraint();
  virtual void destroy();
  typedef monty::rc_ptr<RangedConstraint> t;

  virtual monty::rc_ptr<::mosek::fusion::BoundInterfaceConstraint>
  __mosek_2fusion_2RangedConstraint__upperBoundCon();
  monty::rc_ptr<::mosek::fusion::BoundInterfaceConstraint> upperBoundCon();
  virtual monty::rc_ptr<::mosek::fusion::BoundInterfaceConstraint>
  __mosek_2fusion_2RangedConstraint__lowerBoundCon();
  monty::rc_ptr<::mosek::fusion::BoundInterfaceConstraint> lowerBoundCon();
}; // class RangedConstraint;

class ConicConstraint : public ::mosek::fusion::ModelConstraint {
  ConicConstraint(monty::rc_ptr<::mosek::fusion::ConicConstraint> c,
                  monty::rc_ptr<::mosek::fusion::Model> m);
  ConicConstraint(
      monty::rc_ptr<::mosek::fusion::Model> model, const std::string &name,
      monty::rc_ptr<::mosek::fusion::ConeDomain> dom,
      std::shared_ptr<monty::ndarray<int32_t, 1>> shape, int32_t conid,
      std::shared_ptr<monty::ndarray<int32_t, 1>> nativeidxs,
      std::shared_ptr<
          monty::ndarray<std::shared_ptr<monty::ndarray<std::string, 1>>, 1>>
          indexnames);

protected:
  ConicConstraint(p_ConicConstraint *_impl);

public:
  ConicConstraint(const ConicConstraint &) = delete;
  const ConicConstraint &operator=(const ConicConstraint &) = delete;
  friend class p_ConicConstraint;
  virtual ~ConicConstraint();
  virtual void destroy();
  typedef monty::rc_ptr<ConicConstraint> t;

  virtual /* override */ std::string toString();
}; // class ConicConstraint;

class LinearConstraint : public ::mosek::fusion::ModelConstraint {
  LinearConstraint(monty::rc_ptr<::mosek::fusion::LinearConstraint> c,
                   monty::rc_ptr<::mosek::fusion::Model> m);
  LinearConstraint(
      monty::rc_ptr<::mosek::fusion::Model> model, const std::string &name,
      int32_t conid, std::shared_ptr<monty::ndarray<int32_t, 1>> shape,
      std::shared_ptr<monty::ndarray<int32_t, 1>> nidxs,
      std::shared_ptr<
          monty::ndarray<std::shared_ptr<monty::ndarray<std::string, 1>>, 1>>
          indexnames);

protected:
  LinearConstraint(p_LinearConstraint *_impl);

public:
  LinearConstraint(const LinearConstraint &) = delete;
  const LinearConstraint &operator=(const LinearConstraint &) = delete;
  friend class p_LinearConstraint;
  virtual ~LinearConstraint();
  virtual void destroy();
  typedef monty::rc_ptr<LinearConstraint> t;

  virtual /* override */ std::string toString();
}; // class LinearConstraint;

class Set : public virtual monty::RefCounted {
public:
  p_Set *_impl;

protected:
  Set(p_Set *_impl);

public:
  Set(const Set &) = delete;
  const Set &operator=(const Set &) = delete;
  friend class p_Set;
  virtual ~Set();
  virtual void destroy();
  typedef monty::rc_ptr<Set> t;

  static int64_t size(std::shared_ptr<monty::ndarray<int32_t, 1>> shape);
  static bool match(std::shared_ptr<monty::ndarray<int32_t, 1>> s1,
                    std::shared_ptr<monty::ndarray<int32_t, 1>> s2);
  static int64_t linearidx(std::shared_ptr<monty::ndarray<int32_t, 1>> shape,
                           std::shared_ptr<monty::ndarray<int32_t, 1>> key);
  static std::shared_ptr<monty::ndarray<int32_t, 1>>
  idxtokey(std::shared_ptr<monty::ndarray<int32_t, 1>> shape, int64_t idx);
  static void idxtokey(std::shared_ptr<monty::ndarray<int32_t, 1>> shape,
                       int64_t idx,
                       std::shared_ptr<monty::ndarray<int32_t, 1>> dest);
  static std::string
  indexToString(std::shared_ptr<monty::ndarray<int32_t, 1>> shape, int64_t key);
  static std::string
  keyToString(std::shared_ptr<monty::ndarray<int32_t, 1>> key);
  static void indexToKey(std::shared_ptr<monty::ndarray<int32_t, 1>> shape,
                         int64_t key,
                         std::shared_ptr<monty::ndarray<int32_t, 1>> res);
  static std::shared_ptr<monty::ndarray<int64_t, 1>>
  strides(std::shared_ptr<monty::ndarray<int32_t, 1>> shape);
  static std::shared_ptr<monty::ndarray<int32_t, 1>>
  make(std::shared_ptr<monty::ndarray<int32_t, 1>> set1,
       std::shared_ptr<monty::ndarray<int32_t, 1>> set2);
  static std::shared_ptr<monty::ndarray<int32_t, 1>>
  make(std::shared_ptr<monty::ndarray<int32_t, 1>> sizes);
  static std::shared_ptr<monty::ndarray<int32_t, 1>>
  make(int32_t s1, int32_t s2, int32_t s3);
  static std::shared_ptr<monty::ndarray<int32_t, 1>> make(int32_t s1,
                                                          int32_t s2);
  static std::shared_ptr<monty::ndarray<int32_t, 1>> make(int32_t sz);
  static std::shared_ptr<monty::ndarray<int32_t, 1>> scalar();
  static std::shared_ptr<monty::ndarray<int32_t, 1>>
  make(std::shared_ptr<monty::ndarray<std::string, 1>> names);
}; // class Set;

class ConeDomain : public virtual monty::RefCounted {
public:
  p_ConeDomain *_impl;
  ConeDomain(mosek::fusion::QConeKey k,
             std::shared_ptr<monty::ndarray<double, 1>> alpha,
             std::shared_ptr<monty::ndarray<int32_t, 1>> d);
  ConeDomain(mosek::fusion::QConeKey k,
             std::shared_ptr<monty::ndarray<int32_t, 1>> d);
  ConeDomain(monty::rc_ptr<::mosek::fusion::ConeDomain> other);

protected:
  ConeDomain(p_ConeDomain *_impl);

public:
  ConeDomain(const ConeDomain &) = delete;
  const ConeDomain &operator=(const ConeDomain &) = delete;
  friend class p_ConeDomain;
  virtual ~ConeDomain();
  virtual void destroy();
  typedef monty::rc_ptr<ConeDomain> t;

  virtual monty::rc_ptr<::mosek::fusion::ConeDomain>
  __mosek_2fusion_2ConeDomain__integral();
  monty::rc_ptr<::mosek::fusion::ConeDomain> integral();
  virtual bool axisIsSet();
  virtual int32_t getAxis();
  virtual monty::rc_ptr<::mosek::fusion::ConeDomain>
  __mosek_2fusion_2ConeDomain__axis(int32_t a);
  monty::rc_ptr<::mosek::fusion::ConeDomain> axis(int32_t a);
  virtual monty::rc_ptr<::mosek::fusion::ConeDomain>
  __mosek_2fusion_2ConeDomain__withShape(int32_t dim0, int32_t dim1,
                                         int32_t dim2);
  monty::rc_ptr<::mosek::fusion::ConeDomain>
  withShape(int32_t dim0, int32_t dim1, int32_t dim2);
  virtual monty::rc_ptr<::mosek::fusion::ConeDomain>
  __mosek_2fusion_2ConeDomain__withShape(int32_t dim0, int32_t dim1);
  monty::rc_ptr<::mosek::fusion::ConeDomain> withShape(int32_t dim0,
                                                       int32_t dim1);
  virtual monty::rc_ptr<::mosek::fusion::ConeDomain>
  __mosek_2fusion_2ConeDomain__withShape(int32_t dim0);
  monty::rc_ptr<::mosek::fusion::ConeDomain> withShape(int32_t dim0);
  virtual monty::rc_ptr<::mosek::fusion::ConeDomain>
  __mosek_2fusion_2ConeDomain__withShape(
      std::shared_ptr<monty::ndarray<int32_t, 1>> shp);
  monty::rc_ptr<::mosek::fusion::ConeDomain>
  withShape(std::shared_ptr<monty::ndarray<int32_t, 1>> shp);
  virtual monty::rc_ptr<::mosek::fusion::ConeDomain>
  __mosek_2fusion_2ConeDomain__withNamesOnAxis(
      std::shared_ptr<monty::ndarray<std::string, 1>> names, int32_t axis);
  monty::rc_ptr<::mosek::fusion::ConeDomain>
  withNamesOnAxis(std::shared_ptr<monty::ndarray<std::string, 1>> names,
                  int32_t axis);
}; // class ConeDomain;

class PSDDomain : public virtual monty::RefCounted {
public:
  p_PSDDomain *_impl;
  PSDDomain(mosek::fusion::PSDKey k,
            std::shared_ptr<monty::ndarray<int32_t, 1>> shp, int32_t conedim1,
            int32_t conedim2);
  PSDDomain(mosek::fusion::PSDKey k,
            std::shared_ptr<monty::ndarray<int32_t, 1>> shp);
  PSDDomain(mosek::fusion::PSDKey k);
  PSDDomain(monty::rc_ptr<::mosek::fusion::PSDDomain> other);

protected:
  PSDDomain(p_PSDDomain *_impl);

public:
  PSDDomain(const PSDDomain &) = delete;
  const PSDDomain &operator=(const PSDDomain &) = delete;
  friend class p_PSDDomain;
  virtual ~PSDDomain();
  virtual void destroy();
  typedef monty::rc_ptr<PSDDomain> t;

  virtual monty::rc_ptr<::mosek::fusion::PSDDomain>
  __mosek_2fusion_2PSDDomain__axis(int32_t conedim1, int32_t conedim2);
  monty::rc_ptr<::mosek::fusion::PSDDomain> axis(int32_t conedim1,
                                                 int32_t conedim2);
  virtual monty::rc_ptr<::mosek::fusion::PSDDomain>
  __mosek_2fusion_2PSDDomain__withNamesOnAxis(
      std::shared_ptr<monty::ndarray<std::string, 1>> names, int32_t axis);
  monty::rc_ptr<::mosek::fusion::PSDDomain>
  withNamesOnAxis(std::shared_ptr<monty::ndarray<std::string, 1>> names,
                  int32_t axis);
}; // class PSDDomain;

class RangeDomain : public virtual monty::RefCounted {
public:
  p_RangeDomain *_impl;
  RangeDomain(bool scalable, std::shared_ptr<monty::ndarray<double, 1>> lb,
              std::shared_ptr<monty::ndarray<double, 1>> ub,
              std::shared_ptr<monty::ndarray<int32_t, 1>> dims);
  RangeDomain(bool scalable, double lb, double ub,
              std::shared_ptr<monty::ndarray<int32_t, 1>> dims);
  RangeDomain(bool scalable, std::shared_ptr<monty::ndarray<double, 1>> lb,
              std::shared_ptr<monty::ndarray<double, 1>> ub,
              std::shared_ptr<monty::ndarray<int32_t, 1>> dims,
              std::shared_ptr<monty::ndarray<int32_t, 2>> sp);
  RangeDomain(bool scalable, std::shared_ptr<monty::ndarray<double, 1>> lb,
              std::shared_ptr<monty::ndarray<double, 1>> ub,
              std::shared_ptr<monty::ndarray<int32_t, 1>> dims,
              std::shared_ptr<monty::ndarray<int32_t, 2>> sp, int32_t steal);
  RangeDomain(monty::rc_ptr<::mosek::fusion::RangeDomain> other);

protected:
  RangeDomain(p_RangeDomain *_impl);

public:
  RangeDomain(const RangeDomain &) = delete;
  const RangeDomain &operator=(const RangeDomain &) = delete;
  friend class p_RangeDomain;
  virtual ~RangeDomain();
  virtual void destroy();
  typedef monty::rc_ptr<RangeDomain> t;

  virtual monty::rc_ptr<::mosek::fusion::SymmetricRangeDomain>
  __mosek_2fusion_2RangeDomain__symmetric();
  monty::rc_ptr<::mosek::fusion::SymmetricRangeDomain> symmetric();
  virtual monty::rc_ptr<::mosek::fusion::RangeDomain>
  __mosek_2fusion_2RangeDomain__sparse(
      std::shared_ptr<monty::ndarray<int32_t, 2>> sparsity);
  monty::rc_ptr<::mosek::fusion::RangeDomain>
  sparse(std::shared_ptr<monty::ndarray<int32_t, 2>> sparsity);
  virtual monty::rc_ptr<::mosek::fusion::RangeDomain>
  __mosek_2fusion_2RangeDomain__sparse(
      std::shared_ptr<monty::ndarray<int32_t, 1>> sparsity);
  monty::rc_ptr<::mosek::fusion::RangeDomain>
  sparse(std::shared_ptr<monty::ndarray<int32_t, 1>> sparsity);
  virtual monty::rc_ptr<::mosek::fusion::RangeDomain>
  __mosek_2fusion_2RangeDomain__sparse();
  monty::rc_ptr<::mosek::fusion::RangeDomain> sparse();
  virtual monty::rc_ptr<::mosek::fusion::RangeDomain>
  __mosek_2fusion_2RangeDomain__integral();
  monty::rc_ptr<::mosek::fusion::RangeDomain> integral();
  virtual monty::rc_ptr<::mosek::fusion::RangeDomain>
  __mosek_2fusion_2RangeDomain__withShape(int32_t dim0, int32_t dim1,
                                          int32_t dim2);
  monty::rc_ptr<::mosek::fusion::RangeDomain>
  withShape(int32_t dim0, int32_t dim1, int32_t dim2);
  virtual monty::rc_ptr<::mosek::fusion::RangeDomain>
  __mosek_2fusion_2RangeDomain__withShape(int32_t dim0, int32_t dim1);
  monty::rc_ptr<::mosek::fusion::RangeDomain> withShape(int32_t dim0,
                                                        int32_t dim1);
  virtual monty::rc_ptr<::mosek::fusion::RangeDomain>
  __mosek_2fusion_2RangeDomain__withShape(int32_t dim0);
  monty::rc_ptr<::mosek::fusion::RangeDomain> withShape(int32_t dim0);
  virtual monty::rc_ptr<::mosek::fusion::RangeDomain>
  __mosek_2fusion_2RangeDomain__withShape(
      std::shared_ptr<monty::ndarray<int32_t, 1>> shp);
  monty::rc_ptr<::mosek::fusion::RangeDomain>
  withShape(std::shared_ptr<monty::ndarray<int32_t, 1>> shp);
  virtual monty::rc_ptr<::mosek::fusion::RangeDomain>
  __mosek_2fusion_2RangeDomain__withNamesOnAxis(
      std::shared_ptr<monty::ndarray<std::string, 1>> names, int32_t axis);
  monty::rc_ptr<::mosek::fusion::RangeDomain>
  withNamesOnAxis(std::shared_ptr<monty::ndarray<std::string, 1>> names,
                  int32_t axis);
}; // class RangeDomain;

class SymmetricRangeDomain : public ::mosek::fusion::RangeDomain {
  SymmetricRangeDomain(monty::rc_ptr<::mosek::fusion::RangeDomain> other);

protected:
  SymmetricRangeDomain(p_SymmetricRangeDomain *_impl);

public:
  SymmetricRangeDomain(const SymmetricRangeDomain &) = delete;
  const SymmetricRangeDomain &operator=(const SymmetricRangeDomain &) = delete;
  friend class p_SymmetricRangeDomain;
  virtual ~SymmetricRangeDomain();
  virtual void destroy();
  typedef monty::rc_ptr<SymmetricRangeDomain> t;

}; // class SymmetricRangeDomain;

class SymmetricLinearDomain : public virtual monty::RefCounted {
public:
  p_SymmetricLinearDomain *_impl;
  SymmetricLinearDomain(monty::rc_ptr<::mosek::fusion::LinearDomain> other);

protected:
  SymmetricLinearDomain(p_SymmetricLinearDomain *_impl);

public:
  SymmetricLinearDomain(const SymmetricLinearDomain &) = delete;
  const SymmetricLinearDomain &
  operator=(const SymmetricLinearDomain &) = delete;
  friend class p_SymmetricLinearDomain;
  virtual ~SymmetricLinearDomain();
  virtual void destroy();
  typedef monty::rc_ptr<SymmetricLinearDomain> t;

  virtual monty::rc_ptr<::mosek::fusion::SymmetricLinearDomain>
  __mosek_2fusion_2SymmetricLinearDomain__sparse(
      std::shared_ptr<monty::ndarray<int32_t, 2>> sparsity);
  monty::rc_ptr<::mosek::fusion::SymmetricLinearDomain>
  sparse(std::shared_ptr<monty::ndarray<int32_t, 2>> sparsity);
  virtual monty::rc_ptr<::mosek::fusion::SymmetricLinearDomain>
  __mosek_2fusion_2SymmetricLinearDomain__sparse(
      std::shared_ptr<monty::ndarray<int32_t, 1>> sparsity);
  monty::rc_ptr<::mosek::fusion::SymmetricLinearDomain>
  sparse(std::shared_ptr<monty::ndarray<int32_t, 1>> sparsity);
  virtual monty::rc_ptr<::mosek::fusion::SymmetricLinearDomain>
  __mosek_2fusion_2SymmetricLinearDomain__integral();
  monty::rc_ptr<::mosek::fusion::SymmetricLinearDomain> integral();
}; // class SymmetricLinearDomain;

class LinearDomain : public virtual monty::RefCounted {
public:
  p_LinearDomain *_impl;
  LinearDomain(mosek::fusion::RelationKey k, bool scalable,
               std::shared_ptr<monty::ndarray<double, 1>> rhs,
               std::shared_ptr<monty::ndarray<int32_t, 1>> dims);
  LinearDomain(mosek::fusion::RelationKey k, bool scalable,
               std::shared_ptr<monty::ndarray<double, 1>> rhs,
               std::shared_ptr<monty::ndarray<int32_t, 1>> dims,
               std::shared_ptr<monty::ndarray<int32_t, 2>> sp, int32_t steal);
  LinearDomain(monty::rc_ptr<::mosek::fusion::LinearDomain> other);

protected:
  LinearDomain(p_LinearDomain *_impl);

public:
  LinearDomain(const LinearDomain &) = delete;
  const LinearDomain &operator=(const LinearDomain &) = delete;
  friend class p_LinearDomain;
  virtual ~LinearDomain();
  virtual void destroy();
  typedef monty::rc_ptr<LinearDomain> t;

  virtual monty::rc_ptr<::mosek::fusion::SymmetricLinearDomain>
  __mosek_2fusion_2LinearDomain__symmetric();
  monty::rc_ptr<::mosek::fusion::SymmetricLinearDomain> symmetric();
  virtual monty::rc_ptr<::mosek::fusion::LinearDomain>
  __mosek_2fusion_2LinearDomain__sparse(
      std::shared_ptr<monty::ndarray<int32_t, 2>> sparsity);
  monty::rc_ptr<::mosek::fusion::LinearDomain>
  sparse(std::shared_ptr<monty::ndarray<int32_t, 2>> sparsity);
  virtual monty::rc_ptr<::mosek::fusion::LinearDomain>
  __mosek_2fusion_2LinearDomain__sparse(
      std::shared_ptr<monty::ndarray<int32_t, 1>> sparsity);
  monty::rc_ptr<::mosek::fusion::LinearDomain>
  sparse(std::shared_ptr<monty::ndarray<int32_t, 1>> sparsity);
  virtual monty::rc_ptr<::mosek::fusion::LinearDomain>
  __mosek_2fusion_2LinearDomain__sparse();
  monty::rc_ptr<::mosek::fusion::LinearDomain> sparse();
  virtual monty::rc_ptr<::mosek::fusion::LinearDomain>
  __mosek_2fusion_2LinearDomain__integral();
  monty::rc_ptr<::mosek::fusion::LinearDomain> integral();
  virtual monty::rc_ptr<::mosek::fusion::LinearDomain>
  __mosek_2fusion_2LinearDomain__withShape(int32_t dim0, int32_t dim1,
                                           int32_t dim2);
  monty::rc_ptr<::mosek::fusion::LinearDomain>
  withShape(int32_t dim0, int32_t dim1, int32_t dim2);
  virtual monty::rc_ptr<::mosek::fusion::LinearDomain>
  __mosek_2fusion_2LinearDomain__withShape(int32_t dim0, int32_t dim1);
  monty::rc_ptr<::mosek::fusion::LinearDomain> withShape(int32_t dim0,
                                                         int32_t dim1);
  virtual monty::rc_ptr<::mosek::fusion::LinearDomain>
  __mosek_2fusion_2LinearDomain__withShape(int32_t dim0);
  monty::rc_ptr<::mosek::fusion::LinearDomain> withShape(int32_t dim0);
  virtual monty::rc_ptr<::mosek::fusion::LinearDomain>
  __mosek_2fusion_2LinearDomain__withShape(
      std::shared_ptr<monty::ndarray<int32_t, 1>> shp);
  monty::rc_ptr<::mosek::fusion::LinearDomain>
  withShape(std::shared_ptr<monty::ndarray<int32_t, 1>> shp);
  virtual monty::rc_ptr<::mosek::fusion::LinearDomain>
  __mosek_2fusion_2LinearDomain__withNamesOnAxis(
      std::shared_ptr<monty::ndarray<std::string, 1>> names, int32_t axis);
  monty::rc_ptr<::mosek::fusion::LinearDomain>
  withNamesOnAxis(std::shared_ptr<monty::ndarray<std::string, 1>> names,
                  int32_t axis);
}; // class LinearDomain;

class Domain : public virtual monty::RefCounted {
public:
  p_Domain *_impl;

protected:
  Domain(p_Domain *_impl);

public:
  Domain(const Domain &) = delete;
  const Domain &operator=(const Domain &) = delete;
  friend class p_Domain;
  virtual ~Domain();
  virtual void destroy();
  typedef monty::rc_ptr<Domain> t;

  static monty::rc_ptr<::mosek::fusion::SymmetricRangeDomain>
  symmetric(monty::rc_ptr<::mosek::fusion::RangeDomain> rd);
  static monty::rc_ptr<::mosek::fusion::SymmetricLinearDomain>
  symmetric(monty::rc_ptr<::mosek::fusion::LinearDomain> ld);
  static monty::rc_ptr<::mosek::fusion::RangeDomain>
  sparse(monty::rc_ptr<::mosek::fusion::RangeDomain> rd,
         std::shared_ptr<monty::ndarray<int32_t, 2>> sparsity);
  static monty::rc_ptr<::mosek::fusion::RangeDomain>
  sparse(monty::rc_ptr<::mosek::fusion::RangeDomain> rd,
         std::shared_ptr<monty::ndarray<int32_t, 1>> sparsity);
  static monty::rc_ptr<::mosek::fusion::LinearDomain>
  sparse(monty::rc_ptr<::mosek::fusion::LinearDomain> ld,
         std::shared_ptr<monty::ndarray<int32_t, 2>> sparsity);
  static monty::rc_ptr<::mosek::fusion::LinearDomain>
  sparse(monty::rc_ptr<::mosek::fusion::LinearDomain> ld,
         std::shared_ptr<monty::ndarray<int32_t, 1>> sparsity);
  static monty::rc_ptr<::mosek::fusion::RangeDomain>
  integral(monty::rc_ptr<::mosek::fusion::RangeDomain> rd);
  static monty::rc_ptr<::mosek::fusion::LinearDomain>
  integral(monty::rc_ptr<::mosek::fusion::LinearDomain> ld);
  static monty::rc_ptr<::mosek::fusion::ConeDomain>
  integral(monty::rc_ptr<::mosek::fusion::ConeDomain> c);
  static monty::rc_ptr<::mosek::fusion::ConeDomain>
  axis(monty::rc_ptr<::mosek::fusion::ConeDomain> c, int32_t a);
  static monty::rc_ptr<::mosek::fusion::ConeDomain>
  inDPowerCone(std::shared_ptr<monty::ndarray<double, 1>> alphas,
               std::shared_ptr<monty::ndarray<int32_t, 1>> dims);
  static monty::rc_ptr<::mosek::fusion::ConeDomain>
  inDPowerCone(std::shared_ptr<monty::ndarray<double, 1>> alphas, int32_t m);
  static monty::rc_ptr<::mosek::fusion::ConeDomain>
  inDPowerCone(std::shared_ptr<monty::ndarray<double, 1>> alphas);
  static monty::rc_ptr<::mosek::fusion::ConeDomain>
  inDPowerCone(double alpha, std::shared_ptr<monty::ndarray<int32_t, 1>> dims);
  static monty::rc_ptr<::mosek::fusion::ConeDomain> inDPowerCone(double alpha,
                                                                 int32_t m);
  static monty::rc_ptr<::mosek::fusion::ConeDomain> inDPowerCone(double alpha);
  static monty::rc_ptr<::mosek::fusion::ConeDomain>
  inPPowerCone(std::shared_ptr<monty::ndarray<double, 1>> alphas,
               std::shared_ptr<monty::ndarray<int32_t, 1>> dims);
  static monty::rc_ptr<::mosek::fusion::ConeDomain>
  inPPowerCone(std::shared_ptr<monty::ndarray<double, 1>> alphas, int32_t m);
  static monty::rc_ptr<::mosek::fusion::ConeDomain>
  inPPowerCone(std::shared_ptr<monty::ndarray<double, 1>> alphas);
  static monty::rc_ptr<::mosek::fusion::ConeDomain>
  inPPowerCone(double alpha, std::shared_ptr<monty::ndarray<int32_t, 1>> dims);
  static monty::rc_ptr<::mosek::fusion::ConeDomain> inPPowerCone(double alpha,
                                                                 int32_t m);
  static monty::rc_ptr<::mosek::fusion::ConeDomain> inPPowerCone(double alpha);
  static monty::rc_ptr<::mosek::fusion::ConeDomain>
  inDExpCone(std::shared_ptr<monty::ndarray<int32_t, 1>> dims);
  static monty::rc_ptr<::mosek::fusion::ConeDomain> inDExpCone(int32_t m);
  static monty::rc_ptr<::mosek::fusion::ConeDomain> inDExpCone();
  static monty::rc_ptr<::mosek::fusion::ConeDomain>
  inPExpCone(std::shared_ptr<monty::ndarray<int32_t, 1>> dims);
  static monty::rc_ptr<::mosek::fusion::ConeDomain> inPExpCone(int32_t m);
  static monty::rc_ptr<::mosek::fusion::ConeDomain> inPExpCone();
  static monty::rc_ptr<::mosek::fusion::ConeDomain>
  inDGeoMeanCone(std::shared_ptr<monty::ndarray<int32_t, 1>> dims);
  static monty::rc_ptr<::mosek::fusion::ConeDomain> inDGeoMeanCone(int32_t m,
                                                                   int32_t n);
  static monty::rc_ptr<::mosek::fusion::ConeDomain> inDGeoMeanCone(int32_t n);
  static monty::rc_ptr<::mosek::fusion::ConeDomain> inDGeoMeanCone();
  static monty::rc_ptr<::mosek::fusion::ConeDomain>
  inPGeoMeanCone(std::shared_ptr<monty::ndarray<int32_t, 1>> dims);
  static monty::rc_ptr<::mosek::fusion::ConeDomain> inPGeoMeanCone(int32_t m,
                                                                   int32_t n);
  static monty::rc_ptr<::mosek::fusion::ConeDomain> inPGeoMeanCone(int32_t n);
  static monty::rc_ptr<::mosek::fusion::ConeDomain> inPGeoMeanCone();
  static monty::rc_ptr<::mosek::fusion::ConeDomain>
  inRotatedQCone(std::shared_ptr<monty::ndarray<int32_t, 1>> dims);
  static monty::rc_ptr<::mosek::fusion::ConeDomain> inRotatedQCone(int32_t m,
                                                                   int32_t n);
  static monty::rc_ptr<::mosek::fusion::ConeDomain> inRotatedQCone(int32_t n);
  static monty::rc_ptr<::mosek::fusion::ConeDomain> inRotatedQCone();
  static monty::rc_ptr<::mosek::fusion::ConeDomain>
  inQCone(std::shared_ptr<monty::ndarray<int32_t, 1>> dims);
  static monty::rc_ptr<::mosek::fusion::ConeDomain> inQCone(int32_t m,
                                                            int32_t n);
  static monty::rc_ptr<::mosek::fusion::ConeDomain> inQCone(int32_t n);
  static monty::rc_ptr<::mosek::fusion::ConeDomain> inQCone();
  static monty::rc_ptr<::mosek::fusion::ConeDomain>
  inSVecPSDCone(std::shared_ptr<monty::ndarray<int32_t, 1>> dims);
  static monty::rc_ptr<::mosek::fusion::ConeDomain> inSVecPSDCone(int32_t d1,
                                                                  int32_t d2);
  static monty::rc_ptr<::mosek::fusion::ConeDomain> inSVecPSDCone(int32_t n);
  static monty::rc_ptr<::mosek::fusion::ConeDomain> inSVecPSDCone();
  static monty::rc_ptr<::mosek::fusion::PSDDomain> isTrilPSD(int32_t n,
                                                             int32_t m);
  static monty::rc_ptr<::mosek::fusion::PSDDomain> isTrilPSD(int32_t n);
  static monty::rc_ptr<::mosek::fusion::PSDDomain> isTrilPSD();
  static monty::rc_ptr<::mosek::fusion::PSDDomain> inPSDCone(int32_t n,
                                                             int32_t m);
  static monty::rc_ptr<::mosek::fusion::PSDDomain> inPSDCone(int32_t n);
  static monty::rc_ptr<::mosek::fusion::PSDDomain> inPSDCone();
  static monty::rc_ptr<::mosek::fusion::RangeDomain> binary();
  static monty::rc_ptr<::mosek::fusion::RangeDomain>
  binary(std::shared_ptr<monty::ndarray<int32_t, 1>> dims);
  static monty::rc_ptr<::mosek::fusion::RangeDomain> binary(int32_t m,
                                                            int32_t n);
  static monty::rc_ptr<::mosek::fusion::RangeDomain> binary(int32_t n);
  static monty::rc_ptr<::mosek::fusion::RangeDomain>
  inRange(monty::rc_ptr<::mosek::fusion::Matrix> lbm,
          monty::rc_ptr<::mosek::fusion::Matrix> ubm);
  static monty::rc_ptr<::mosek::fusion::RangeDomain>
  inRange(std::shared_ptr<monty::ndarray<double, 2>> lba,
          std::shared_ptr<monty::ndarray<double, 2>> uba);
  static monty::rc_ptr<::mosek::fusion::RangeDomain>
  inRange(std::shared_ptr<monty::ndarray<double, 1>> lba,
          std::shared_ptr<monty::ndarray<double, 1>> uba,
          std::shared_ptr<monty::ndarray<int32_t, 1>> dims);
  static monty::rc_ptr<::mosek::fusion::RangeDomain>
  inRange(std::shared_ptr<monty::ndarray<double, 1>> lba, double ub,
          std::shared_ptr<monty::ndarray<int32_t, 1>> dims);
  static monty::rc_ptr<::mosek::fusion::RangeDomain>
  inRange(double lb, std::shared_ptr<monty::ndarray<double, 1>> uba,
          std::shared_ptr<monty::ndarray<int32_t, 1>> dims);
  static monty::rc_ptr<::mosek::fusion::RangeDomain>
  inRange(double lb, double ub,
          std::shared_ptr<monty::ndarray<int32_t, 1>> dims);
  static monty::rc_ptr<::mosek::fusion::RangeDomain>
  inRange(std::shared_ptr<monty::ndarray<double, 1>> lba,
          std::shared_ptr<monty::ndarray<double, 1>> uba);
  static monty::rc_ptr<::mosek::fusion::RangeDomain>
  inRange(std::shared_ptr<monty::ndarray<double, 1>> lba, double ub);
  static monty::rc_ptr<::mosek::fusion::RangeDomain>
  inRange(double lb, std::shared_ptr<monty::ndarray<double, 1>> uba);
  static monty::rc_ptr<::mosek::fusion::RangeDomain> inRange(double lb,
                                                             double ub);
  static monty::rc_ptr<::mosek::fusion::LinearDomain>
  greaterThan(monty::rc_ptr<::mosek::fusion::Matrix> mx);
  static monty::rc_ptr<::mosek::fusion::LinearDomain>
  greaterThan(std::shared_ptr<monty::ndarray<double, 1>> a1,
              std::shared_ptr<monty::ndarray<int32_t, 1>> dims);
  static monty::rc_ptr<::mosek::fusion::LinearDomain>
  greaterThan(std::shared_ptr<monty::ndarray<double, 2>> a2);
  static monty::rc_ptr<::mosek::fusion::LinearDomain>
  greaterThan(std::shared_ptr<monty::ndarray<double, 1>> a1);
  static monty::rc_ptr<::mosek::fusion::LinearDomain>
  greaterThan(double b, std::shared_ptr<monty::ndarray<int32_t, 1>> dims);
  static monty::rc_ptr<::mosek::fusion::LinearDomain>
  greaterThan(double b, int32_t m, int32_t n);
  static monty::rc_ptr<::mosek::fusion::LinearDomain> greaterThan(double b,
                                                                  int32_t n);
  static monty::rc_ptr<::mosek::fusion::LinearDomain> greaterThan(double b);
  static monty::rc_ptr<::mosek::fusion::LinearDomain>
  lessThan(monty::rc_ptr<::mosek::fusion::Matrix> mx);
  static monty::rc_ptr<::mosek::fusion::LinearDomain>
  lessThan(std::shared_ptr<monty::ndarray<double, 1>> a1,
           std::shared_ptr<monty::ndarray<int32_t, 1>> dims);
  static monty::rc_ptr<::mosek::fusion::LinearDomain>
  lessThan(std::shared_ptr<monty::ndarray<double, 2>> a2);
  static monty::rc_ptr<::mosek::fusion::LinearDomain>
  lessThan(std::shared_ptr<monty::ndarray<double, 1>> a1);
  static monty::rc_ptr<::mosek::fusion::LinearDomain>
  lessThan(double b, std::shared_ptr<monty::ndarray<int32_t, 1>> dims);
  static monty::rc_ptr<::mosek::fusion::LinearDomain>
  lessThan(double b, int32_t m, int32_t n);
  static monty::rc_ptr<::mosek::fusion::LinearDomain> lessThan(double b,
                                                               int32_t n);
  static monty::rc_ptr<::mosek::fusion::LinearDomain> lessThan(double b);
  static monty::rc_ptr<::mosek::fusion::LinearDomain>
  equalsTo(monty::rc_ptr<::mosek::fusion::Matrix> mx);
  static monty::rc_ptr<::mosek::fusion::LinearDomain>
  equalsTo(std::shared_ptr<monty::ndarray<double, 1>> a1,
           std::shared_ptr<monty::ndarray<int32_t, 1>> dims);
  static monty::rc_ptr<::mosek::fusion::LinearDomain>
  equalsTo(std::shared_ptr<monty::ndarray<double, 2>> a2);
  static monty::rc_ptr<::mosek::fusion::LinearDomain>
  equalsTo(std::shared_ptr<monty::ndarray<double, 1>> a1);
  static monty::rc_ptr<::mosek::fusion::LinearDomain>
  equalsTo(double b, std::shared_ptr<monty::ndarray<int32_t, 1>> dims);
  static monty::rc_ptr<::mosek::fusion::LinearDomain>
  equalsTo(double b, int32_t m, int32_t n);
  static monty::rc_ptr<::mosek::fusion::LinearDomain> equalsTo(double b,
                                                               int32_t n);
  static monty::rc_ptr<::mosek::fusion::LinearDomain> equalsTo(double b);
  static monty::rc_ptr<::mosek::fusion::LinearDomain>
  unbounded(std::shared_ptr<monty::ndarray<int32_t, 1>> dims);
  static monty::rc_ptr<::mosek::fusion::LinearDomain> unbounded(int32_t m,
                                                                int32_t n);
  static monty::rc_ptr<::mosek::fusion::LinearDomain> unbounded(int32_t n);
  static monty::rc_ptr<::mosek::fusion::LinearDomain> unbounded();
}; // class Domain;

class ExprCode : public virtual monty::RefCounted {
public:
  p_ExprCode *_impl;

protected:
  ExprCode(p_ExprCode *_impl);

public:
  ExprCode(const ExprCode &) = delete;
  const ExprCode &operator=(const ExprCode &) = delete;
  friend class p_ExprCode;
  virtual ~ExprCode();
  virtual void destroy();
  typedef monty::rc_ptr<ExprCode> t;

  static void inplace_relocate(std::shared_ptr<monty::ndarray<int32_t, 1>> code,
                               int32_t from_offset, int32_t num,
                               int32_t const_base);
  static std::string op2str(int32_t op);
  static void eval_add_list(std::shared_ptr<monty::ndarray<int32_t, 1>> code,
                            std::shared_ptr<monty::ndarray<int32_t, 1>> ptr,
                            std::shared_ptr<monty::ndarray<double, 1>> consts,
                            int32_t offset,
                            std::shared_ptr<monty::ndarray<double, 1>> target,
                            std::shared_ptr<monty::ndarray<double, 1>> P,
                            monty::rc_ptr<::mosek::fusion::WorkStack> xs);
  static void eval_add_list(std::shared_ptr<monty::ndarray<int32_t, 1>> code,
                            std::shared_ptr<monty::ndarray<int32_t, 1>> ptr,
                            std::shared_ptr<monty::ndarray<double, 1>> consts,
                            std::shared_ptr<monty::ndarray<double, 1>> target,
                            std::shared_ptr<monty::ndarray<double, 1>> P,
                            monty::rc_ptr<::mosek::fusion::WorkStack> xs);
  static int32_t emit_sum(std::shared_ptr<monty::ndarray<int32_t, 1>> tgt,
                          int32_t ofs, int32_t num);
  static int32_t emit_inv(std::shared_ptr<monty::ndarray<int32_t, 1>> tgt,
                          int32_t ofs);
  static int32_t emit_mul(std::shared_ptr<monty::ndarray<int32_t, 1>> tgt,
                          int32_t ofs);
  static int32_t emit_neg(std::shared_ptr<monty::ndarray<int32_t, 1>> tgt,
                          int32_t ofs);
  static int32_t emit_add(std::shared_ptr<monty::ndarray<int32_t, 1>> tgt,
                          int32_t ofs);
  static int32_t emit_constref(std::shared_ptr<monty::ndarray<int32_t, 1>> tgt,
                               int32_t ofs, int32_t i);
  static int32_t emit_paramref(std::shared_ptr<monty::ndarray<int32_t, 1>> tgt,
                               int32_t ofs, int32_t i);
  static int32_t emit_nop(std::shared_ptr<monty::ndarray<int32_t, 1>> tgt,
                          int32_t ofs);
}; // class ExprCode;

class Param : public virtual monty::RefCounted {
public:
  p_Param *_impl;

protected:
  Param(p_Param *_impl);

public:
  Param(const Param &) = delete;
  const Param &operator=(const Param &) = delete;
  friend class p_Param;
  virtual ~Param();
  virtual void destroy();
  typedef monty::rc_ptr<Param> t;

  static monty::rc_ptr<::mosek::fusion::Parameter>
  repeat(monty::rc_ptr<::mosek::fusion::Parameter> p, int32_t n, int32_t dim);
  static monty::rc_ptr<::mosek::fusion::Parameter>
  stack(int32_t dim, monty::rc_ptr<::mosek::fusion::Parameter> p1,
        monty::rc_ptr<::mosek::fusion::Parameter> p2,
        monty::rc_ptr<::mosek::fusion::Parameter> p3);
  static monty::rc_ptr<::mosek::fusion::Parameter>
  stack(int32_t dim, monty::rc_ptr<::mosek::fusion::Parameter> p1,
        monty::rc_ptr<::mosek::fusion::Parameter> p2);
  static monty::rc_ptr<::mosek::fusion::Parameter>
  stack(int32_t dim,
        std::shared_ptr<
            monty::ndarray<monty::rc_ptr<::mosek::fusion::Parameter>, 1>>
            p);
  static monty::rc_ptr<::mosek::fusion::Parameter>
  stack(std::shared_ptr<
        monty::ndarray<std::shared_ptr<monty::ndarray<
                           monty::rc_ptr<::mosek::fusion::Parameter>, 1>>,
                       1>>
            p);
  static monty::rc_ptr<::mosek::fusion::Parameter>
  hstack(monty::rc_ptr<::mosek::fusion::Parameter> p1,
         monty::rc_ptr<::mosek::fusion::Parameter> p2,
         monty::rc_ptr<::mosek::fusion::Parameter> p3);
  static monty::rc_ptr<::mosek::fusion::Parameter>
  hstack(monty::rc_ptr<::mosek::fusion::Parameter> p1,
         monty::rc_ptr<::mosek::fusion::Parameter> p2);
  static monty::rc_ptr<::mosek::fusion::Parameter>
  hstack(std::shared_ptr<
         monty::ndarray<monty::rc_ptr<::mosek::fusion::Parameter>, 1>>
             p);
  static monty::rc_ptr<::mosek::fusion::Parameter>
  vstack(monty::rc_ptr<::mosek::fusion::Parameter> p1,
         monty::rc_ptr<::mosek::fusion::Parameter> p2,
         monty::rc_ptr<::mosek::fusion::Parameter> p3);
  static monty::rc_ptr<::mosek::fusion::Parameter>
  vstack(monty::rc_ptr<::mosek::fusion::Parameter> p1,
         monty::rc_ptr<::mosek::fusion::Parameter> p2);
  static monty::rc_ptr<::mosek::fusion::Parameter>
  vstack(std::shared_ptr<
         monty::ndarray<monty::rc_ptr<::mosek::fusion::Parameter>, 1>>
             p);
}; // class Param;

class ParameterImpl : public virtual ::mosek::fusion::Parameter {
public:
  p_ParameterImpl *_impl;
  ParameterImpl(monty::rc_ptr<::mosek::fusion::ParameterImpl> other,
                monty::rc_ptr<::mosek::fusion::Model> model);
  ParameterImpl(monty::rc_ptr<::mosek::fusion::Model> model,
                std::shared_ptr<monty::ndarray<int32_t, 1>> shape,
                std::shared_ptr<monty::ndarray<int64_t, 1>> sp,
                std::shared_ptr<monty::ndarray<int32_t, 1>> nidxs);

protected:
  ParameterImpl(p_ParameterImpl *_impl);

public:
  ParameterImpl(const ParameterImpl &) = delete;
  const ParameterImpl &operator=(const ParameterImpl &) = delete;
  friend class p_ParameterImpl;
  virtual ~ParameterImpl();
  virtual void destroy();
  typedef monty::rc_ptr<ParameterImpl> t;

  virtual monty::rc_ptr<::mosek::fusion::Parameter>
  __mosek_2fusion_2ParameterImpl__clone(
      monty::rc_ptr<::mosek::fusion::Model> m);
  monty::rc_ptr<::mosek::fusion::Parameter>
  clone(monty::rc_ptr<::mosek::fusion::Model> m);
  /* override: mosek.fusion.Parameter.clone*/
  virtual monty::rc_ptr<::mosek::fusion::Parameter>
  __mosek_2fusion_2Parameter__clone(monty::rc_ptr<::mosek::fusion::Model> m);
  virtual /* override */ std::string toString();
  virtual monty::rc_ptr<::mosek::fusion::Parameter>
  __mosek_2fusion_2ParameterImpl__pick(
      std::shared_ptr<monty::ndarray<int32_t, 2>> midxs);
  monty::rc_ptr<::mosek::fusion::Parameter>
  pick(std::shared_ptr<monty::ndarray<int32_t, 2>> midxs);
  /* override: mosek.fusion.Parameter.pick*/
  virtual monty::rc_ptr<::mosek::fusion::Parameter>
  __mosek_2fusion_2Parameter__pick(
      std::shared_ptr<monty::ndarray<int32_t, 2>> midxs);
  virtual monty::rc_ptr<::mosek::fusion::Parameter>
  __mosek_2fusion_2ParameterImpl__pick(
      std::shared_ptr<monty::ndarray<int32_t, 1>> idxs);
  monty::rc_ptr<::mosek::fusion::Parameter>
  pick(std::shared_ptr<monty::ndarray<int32_t, 1>> idxs);
  /* override: mosek.fusion.Parameter.pick*/
  virtual monty::rc_ptr<::mosek::fusion::Parameter>
  __mosek_2fusion_2Parameter__pick(
      std::shared_ptr<monty::ndarray<int32_t, 1>> idxs);
  virtual monty::rc_ptr<::mosek::fusion::Expression>
  __mosek_2fusion_2ParameterImpl__index(
      std::shared_ptr<monty::ndarray<int32_t, 1>> indexes);
  monty::rc_ptr<::mosek::fusion::Expression>
  index(std::shared_ptr<monty::ndarray<int32_t, 1>> indexes);
  /* override: mosek.fusion.Expression.index*/
  virtual monty::rc_ptr<::mosek::fusion::Expression>
  __mosek_2fusion_2Expression__index(
      std::shared_ptr<monty::ndarray<int32_t, 1>> indexes);
  virtual monty::rc_ptr<::mosek::fusion::Expression>
  __mosek_2fusion_2ParameterImpl__index(int32_t i);
  monty::rc_ptr<::mosek::fusion::Expression> index(int32_t i);
  /* override: mosek.fusion.Expression.index*/
  virtual monty::rc_ptr<::mosek::fusion::Expression>
  __mosek_2fusion_2Expression__index(int32_t i);
  virtual void eval(monty::rc_ptr<::mosek::fusion::WorkStack> rs,
                    monty::rc_ptr<::mosek::fusion::WorkStack> ws,
                    monty::rc_ptr<::mosek::fusion::WorkStack> xs);
  virtual void getSp(std::shared_ptr<monty::ndarray<int64_t, 1>> dest,
                     int32_t offset);
  virtual bool isSparse();
  virtual monty::rc_ptr<::mosek::fusion::Parameter>
  __mosek_2fusion_2ParameterImpl__slice(
      std::shared_ptr<monty::ndarray<int32_t, 1>> astart,
      std::shared_ptr<monty::ndarray<int32_t, 1>> astop);
  monty::rc_ptr<::mosek::fusion::Parameter>
  slice(std::shared_ptr<monty::ndarray<int32_t, 1>> astart,
        std::shared_ptr<monty::ndarray<int32_t, 1>> astop);
  /* override: mosek.fusion.Parameter.slice*/
  virtual monty::rc_ptr<::mosek::fusion::Parameter>
  __mosek_2fusion_2Parameter__slice(
      std::shared_ptr<monty::ndarray<int32_t, 1>> astart,
      std::shared_ptr<monty::ndarray<int32_t, 1>> astop);
  virtual monty::rc_ptr<::mosek::fusion::Parameter>
  __mosek_2fusion_2ParameterImpl__slice(int32_t start, int32_t stop);
  monty::rc_ptr<::mosek::fusion::Parameter> slice(int32_t start, int32_t stop);
  /* override: mosek.fusion.Parameter.slice*/
  virtual monty::rc_ptr<::mosek::fusion::Parameter>
  __mosek_2fusion_2Parameter__slice(int32_t start, int32_t stop);
  virtual monty::rc_ptr<::mosek::fusion::Parameter>
  __mosek_2fusion_2ParameterImpl__reshape(
      std::shared_ptr<monty::ndarray<int32_t, 1>> dims);
  monty::rc_ptr<::mosek::fusion::Parameter>
  reshape(std::shared_ptr<monty::ndarray<int32_t, 1>> dims);
  /* override: mosek.fusion.Parameter.reshape*/
  virtual monty::rc_ptr<::mosek::fusion::Parameter>
  __mosek_2fusion_2Parameter__reshape(
      std::shared_ptr<monty::ndarray<int32_t, 1>> dims);
  virtual monty::rc_ptr<::mosek::fusion::Expression>
  __mosek_2fusion_2ParameterImpl__asExpr();
  monty::rc_ptr<::mosek::fusion::Expression> asExpr();
  /* override: mosek.fusion.Parameter.asExpr*/
  virtual monty::rc_ptr<::mosek::fusion::Expression>
  __mosek_2fusion_2Parameter__asExpr();
  virtual int64_t getSize();
  virtual int32_t getNumNonzero();
  virtual int32_t getND();
  virtual std::shared_ptr<monty::ndarray<int32_t, 1>> getShape();
  virtual int32_t getDim(int32_t i);
  virtual void getAllIndexes(std::shared_ptr<monty::ndarray<int32_t, 1>> dst,
                             int32_t ofs);
  virtual int32_t getIndex(int32_t i);
  virtual std::shared_ptr<monty::ndarray<double, 1>> getValue();
  virtual void setValue(std::shared_ptr<monty::ndarray<double, 2>> values2);
  virtual void setValue(std::shared_ptr<monty::ndarray<double, 1>> values);
  virtual void setValue(double value);
  virtual monty::rc_ptr<::mosek::fusion::Model>
  __mosek_2fusion_2ParameterImpl__getModel();
  monty::rc_ptr<::mosek::fusion::Model> getModel();
  /* override: mosek.fusion.Parameter.getModel*/
  virtual monty::rc_ptr<::mosek::fusion::Model>
  __mosek_2fusion_2Parameter__getModel();
}; // class ParameterImpl;

class ExprRangeDomain : public virtual ::mosek::fusion::ExprDomain {
public:
  p_ExprRangeDomain *_impl;

protected:
  ExprRangeDomain(p_ExprRangeDomain *_impl);

public:
  ExprRangeDomain(const ExprRangeDomain &) = delete;
  const ExprRangeDomain &operator=(const ExprRangeDomain &) = delete;
  friend class p_ExprRangeDomain;
  virtual ~ExprRangeDomain();
  virtual void destroy();
  typedef monty::rc_ptr<ExprRangeDomain> t;
  monty::rc_ptr<::mosek::fusion::RangeDomain> get_dom();
  void set_dom(monty::rc_ptr<::mosek::fusion::RangeDomain> val);
  monty::rc_ptr<::mosek::fusion::Expression> get_expr();
  void set_expr(monty::rc_ptr<::mosek::fusion::Expression> val);

  ExprRangeDomain(monty::rc_ptr<::mosek::fusion::Expression> expr,
                  monty::rc_ptr<::mosek::fusion::RangeDomain> dom);
  virtual monty::rc_ptr<::mosek::fusion::Term>
  __mosek_2fusion_2ExprRangeDomain__toDJCTerm();
  monty::rc_ptr<::mosek::fusion::Term> toDJCTerm();
  /* override: mosek.fusion.ExprDomain.toDJCTerm*/
  virtual monty::rc_ptr<::mosek::fusion::Term>
  __mosek_2fusion_2ExprDomain__toDJCTerm();
}; // class ExprRangeDomain;

class ExprPSDDomain : public virtual monty::RefCounted {
public:
  p_ExprPSDDomain *_impl;

protected:
  ExprPSDDomain(p_ExprPSDDomain *_impl);

public:
  ExprPSDDomain(const ExprPSDDomain &) = delete;
  const ExprPSDDomain &operator=(const ExprPSDDomain &) = delete;
  friend class p_ExprPSDDomain;
  virtual ~ExprPSDDomain();
  virtual void destroy();
  typedef monty::rc_ptr<ExprPSDDomain> t;
  monty::rc_ptr<::mosek::fusion::PSDDomain> get_dom();
  void set_dom(monty::rc_ptr<::mosek::fusion::PSDDomain> val);
  monty::rc_ptr<::mosek::fusion::Expression> get_expr();
  void set_expr(monty::rc_ptr<::mosek::fusion::Expression> val);

  ExprPSDDomain(monty::rc_ptr<::mosek::fusion::Expression> expr,
                monty::rc_ptr<::mosek::fusion::PSDDomain> dom);
}; // class ExprPSDDomain;

class ExprConicDomain : public virtual monty::RefCounted {
public:
  p_ExprConicDomain *_impl;

protected:
  ExprConicDomain(p_ExprConicDomain *_impl);

public:
  ExprConicDomain(const ExprConicDomain &) = delete;
  const ExprConicDomain &operator=(const ExprConicDomain &) = delete;
  friend class p_ExprConicDomain;
  virtual ~ExprConicDomain();
  virtual void destroy();
  typedef monty::rc_ptr<ExprConicDomain> t;
  monty::rc_ptr<::mosek::fusion::ConeDomain> get_dom();
  void set_dom(monty::rc_ptr<::mosek::fusion::ConeDomain> val);
  monty::rc_ptr<::mosek::fusion::Expression> get_expr();
  void set_expr(monty::rc_ptr<::mosek::fusion::Expression> val);

  ExprConicDomain(monty::rc_ptr<::mosek::fusion::Expression> expr,
                  monty::rc_ptr<::mosek::fusion::ConeDomain> dom);
}; // class ExprConicDomain;

class ExprLinearDomain : public virtual ::mosek::fusion::ExprDomain {
public:
  p_ExprLinearDomain *_impl;

protected:
  ExprLinearDomain(p_ExprLinearDomain *_impl);

public:
  ExprLinearDomain(const ExprLinearDomain &) = delete;
  const ExprLinearDomain &operator=(const ExprLinearDomain &) = delete;
  friend class p_ExprLinearDomain;
  virtual ~ExprLinearDomain();
  virtual void destroy();
  typedef monty::rc_ptr<ExprLinearDomain> t;
  monty::rc_ptr<::mosek::fusion::LinearDomain> get_dom();
  void set_dom(monty::rc_ptr<::mosek::fusion::LinearDomain> val);
  monty::rc_ptr<::mosek::fusion::Expression> get_expr();
  void set_expr(monty::rc_ptr<::mosek::fusion::Expression> val);

  ExprLinearDomain(monty::rc_ptr<::mosek::fusion::Expression> expr,
                   monty::rc_ptr<::mosek::fusion::LinearDomain> dom);
  virtual monty::rc_ptr<::mosek::fusion::Term>
  __mosek_2fusion_2ExprLinearDomain__toDJCTerm();
  monty::rc_ptr<::mosek::fusion::Term> toDJCTerm();
  /* override: mosek.fusion.ExprDomain.toDJCTerm*/
  virtual monty::rc_ptr<::mosek::fusion::Term>
  __mosek_2fusion_2ExprDomain__toDJCTerm();
}; // class ExprLinearDomain;

class BaseExpression : public virtual ::mosek::fusion::Expression {
public:
  p_BaseExpression *_impl;

protected:
  BaseExpression(p_BaseExpression *_impl);

public:
  BaseExpression(const BaseExpression &) = delete;
  const BaseExpression &operator=(const BaseExpression &) = delete;
  friend class p_BaseExpression;
  virtual ~BaseExpression();
  virtual void destroy();
  typedef monty::rc_ptr<BaseExpression> t;

  BaseExpression(std::shared_ptr<monty::ndarray<int32_t, 1>> shape);
  virtual /* override */ std::string toString();
  virtual void eval(monty::rc_ptr<::mosek::fusion::WorkStack> rs,
                    monty::rc_ptr<::mosek::fusion::WorkStack> ws,
                    monty::rc_ptr<::mosek::fusion::WorkStack> xs) = 0;
  virtual monty::rc_ptr<::mosek::fusion::Expression>
  __mosek_2fusion_2BaseExpression__pick(
      std::shared_ptr<monty::ndarray<int32_t, 2>> indexrows);
  monty::rc_ptr<::mosek::fusion::Expression>
  pick(std::shared_ptr<monty::ndarray<int32_t, 2>> indexrows);
  /* override: mosek.fusion.Expression.pick*/
  virtual monty::rc_ptr<::mosek::fusion::Expression>
  __mosek_2fusion_2Expression__pick(
      std::shared_ptr<monty::ndarray<int32_t, 2>> indexrows);
  virtual monty::rc_ptr<::mosek::fusion::Expression>
  __mosek_2fusion_2BaseExpression__pick(
      std::shared_ptr<monty::ndarray<int32_t, 1>> indexes);
  monty::rc_ptr<::mosek::fusion::Expression>
  pick(std::shared_ptr<monty::ndarray<int32_t, 1>> indexes);
  /* override: mosek.fusion.Expression.pick*/
  virtual monty::rc_ptr<::mosek::fusion::Expression>
  __mosek_2fusion_2Expression__pick(
      std::shared_ptr<monty::ndarray<int32_t, 1>> indexes);
  virtual monty::rc_ptr<::mosek::fusion::Expression>
  __mosek_2fusion_2BaseExpression__index(
      std::shared_ptr<monty::ndarray<int32_t, 1>> indexes);
  monty::rc_ptr<::mosek::fusion::Expression>
  index(std::shared_ptr<monty::ndarray<int32_t, 1>> indexes);
  /* override: mosek.fusion.Expression.index*/
  virtual monty::rc_ptr<::mosek::fusion::Expression>
  __mosek_2fusion_2Expression__index(
      std::shared_ptr<monty::ndarray<int32_t, 1>> indexes);
  virtual monty::rc_ptr<::mosek::fusion::Expression>
  __mosek_2fusion_2BaseExpression__index(int32_t i);
  monty::rc_ptr<::mosek::fusion::Expression> index(int32_t i);
  /* override: mosek.fusion.Expression.index*/
  virtual monty::rc_ptr<::mosek::fusion::Expression>
  __mosek_2fusion_2Expression__index(int32_t i);
  virtual monty::rc_ptr<::mosek::fusion::Expression>
  __mosek_2fusion_2BaseExpression__slice(
      std::shared_ptr<monty::ndarray<int32_t, 1>> firsta,
      std::shared_ptr<monty::ndarray<int32_t, 1>> lasta);
  monty::rc_ptr<::mosek::fusion::Expression>
  slice(std::shared_ptr<monty::ndarray<int32_t, 1>> firsta,
        std::shared_ptr<monty::ndarray<int32_t, 1>> lasta);
  /* override: mosek.fusion.Expression.slice*/
  virtual monty::rc_ptr<::mosek::fusion::Expression>
  __mosek_2fusion_2Expression__slice(
      std::shared_ptr<monty::ndarray<int32_t, 1>> firsta,
      std::shared_ptr<monty::ndarray<int32_t, 1>> lasta);
  virtual monty::rc_ptr<::mosek::fusion::Expression>
  __mosek_2fusion_2BaseExpression__slice(int32_t first, int32_t last);
  monty::rc_ptr<::mosek::fusion::Expression> slice(int32_t first, int32_t last);
  /* override: mosek.fusion.Expression.slice*/
  virtual monty::rc_ptr<::mosek::fusion::Expression>
  __mosek_2fusion_2Expression__slice(int32_t first, int32_t last);
  virtual int64_t getSize();
  virtual int32_t getND();
  virtual int32_t getDim(int32_t d);
  virtual std::shared_ptr<monty::ndarray<int32_t, 1>> getShape();
}; // class BaseExpression;

class ExprParameter : public ::mosek::fusion::BaseExpression {
protected:
  ExprParameter(p_ExprParameter *_impl);

public:
  ExprParameter(const ExprParameter &) = delete;
  const ExprParameter &operator=(const ExprParameter &) = delete;
  friend class p_ExprParameter;
  virtual ~ExprParameter();
  virtual void destroy();
  typedef monty::rc_ptr<ExprParameter> t;

  ExprParameter(monty::rc_ptr<::mosek::fusion::Parameter> p);
  virtual /* override */ void
  eval(monty::rc_ptr<::mosek::fusion::WorkStack> rs,
       monty::rc_ptr<::mosek::fusion::WorkStack> ws,
       monty::rc_ptr<::mosek::fusion::WorkStack> xs);
  virtual /* override */ monty::rc_ptr<::mosek::fusion::Expression>
  __mosek_2fusion_2ExprParameter__slice(
      std::shared_ptr<monty::ndarray<int32_t, 1>> start,
      std::shared_ptr<monty::ndarray<int32_t, 1>> stop);
  monty::rc_ptr<::mosek::fusion::Expression>
  slice(std::shared_ptr<monty::ndarray<int32_t, 1>> start,
        std::shared_ptr<monty::ndarray<int32_t, 1>> stop);
  /* override: mosek.fusion.BaseExpression.slice*/
  virtual monty::rc_ptr<::mosek::fusion::Expression>
  __mosek_2fusion_2BaseExpression__slice(
      std::shared_ptr<monty::ndarray<int32_t, 1>> start,
      std::shared_ptr<monty::ndarray<int32_t, 1>> stop);
  virtual /* override */ monty::rc_ptr<::mosek::fusion::Expression>
  __mosek_2fusion_2ExprParameter__slice(int32_t start, int32_t stop);
  monty::rc_ptr<::mosek::fusion::Expression> slice(int32_t start, int32_t stop);
  /* override: mosek.fusion.BaseExpression.slice*/
  virtual monty::rc_ptr<::mosek::fusion::Expression>
  __mosek_2fusion_2BaseExpression__slice(int32_t start, int32_t stop);
  virtual /* override */ std::string toString();
}; // class ExprParameter;

class ExprMulParamScalarExpr : public ::mosek::fusion::BaseExpression {
protected:
  ExprMulParamScalarExpr(p_ExprMulParamScalarExpr *_impl);

public:
  ExprMulParamScalarExpr(const ExprMulParamScalarExpr &) = delete;
  const ExprMulParamScalarExpr &
  operator=(const ExprMulParamScalarExpr &) = delete;
  friend class p_ExprMulParamScalarExpr;
  virtual ~ExprMulParamScalarExpr();
  virtual void destroy();
  typedef monty::rc_ptr<ExprMulParamScalarExpr> t;

  ExprMulParamScalarExpr(monty::rc_ptr<::mosek::fusion::Parameter> p,
                         monty::rc_ptr<::mosek::fusion::Expression> e);
  virtual /* override */ void
  eval(monty::rc_ptr<::mosek::fusion::WorkStack> rs,
       monty::rc_ptr<::mosek::fusion::WorkStack> ws,
       monty::rc_ptr<::mosek::fusion::WorkStack> xs);
  virtual /* override */ std::string toString();
}; // class ExprMulParamScalarExpr;

class ExprMulParamScalar : public ::mosek::fusion::BaseExpression {
protected:
  ExprMulParamScalar(p_ExprMulParamScalar *_impl);

public:
  ExprMulParamScalar(const ExprMulParamScalar &) = delete;
  const ExprMulParamScalar &operator=(const ExprMulParamScalar &) = delete;
  friend class p_ExprMulParamScalar;
  virtual ~ExprMulParamScalar();
  virtual void destroy();
  typedef monty::rc_ptr<ExprMulParamScalar> t;

  ExprMulParamScalar(monty::rc_ptr<::mosek::fusion::Parameter> p,
                     monty::rc_ptr<::mosek::fusion::Expression> e);
  virtual /* override */ void
  eval(monty::rc_ptr<::mosek::fusion::WorkStack> rs,
       monty::rc_ptr<::mosek::fusion::WorkStack> ws,
       monty::rc_ptr<::mosek::fusion::WorkStack> xs);
  virtual /* override */ std::string toString();
}; // class ExprMulParamScalar;

class ExprMulParamDiagLeft : public ::mosek::fusion::BaseExpression {
protected:
  ExprMulParamDiagLeft(p_ExprMulParamDiagLeft *_impl);

public:
  ExprMulParamDiagLeft(const ExprMulParamDiagLeft &) = delete;
  const ExprMulParamDiagLeft &operator=(const ExprMulParamDiagLeft &) = delete;
  friend class p_ExprMulParamDiagLeft;
  virtual ~ExprMulParamDiagLeft();
  virtual void destroy();
  typedef monty::rc_ptr<ExprMulParamDiagLeft> t;

  ExprMulParamDiagLeft(monty::rc_ptr<::mosek::fusion::Parameter> p,
                       monty::rc_ptr<::mosek::fusion::Expression> e);
  virtual /* override */ void
  eval(monty::rc_ptr<::mosek::fusion::WorkStack> rs,
       monty::rc_ptr<::mosek::fusion::WorkStack> ws,
       monty::rc_ptr<::mosek::fusion::WorkStack> xs);
  virtual /* override */ std::string toString();
}; // class ExprMulParamDiagLeft;

class ExprMulParamDiagRight : public ::mosek::fusion::BaseExpression {
protected:
  ExprMulParamDiagRight(p_ExprMulParamDiagRight *_impl);

public:
  ExprMulParamDiagRight(const ExprMulParamDiagRight &) = delete;
  const ExprMulParamDiagRight &
  operator=(const ExprMulParamDiagRight &) = delete;
  friend class p_ExprMulParamDiagRight;
  virtual ~ExprMulParamDiagRight();
  virtual void destroy();
  typedef monty::rc_ptr<ExprMulParamDiagRight> t;

  ExprMulParamDiagRight(monty::rc_ptr<::mosek::fusion::Expression> e,
                        monty::rc_ptr<::mosek::fusion::Parameter> p);
  virtual /* override */ void
  eval(monty::rc_ptr<::mosek::fusion::WorkStack> rs,
       monty::rc_ptr<::mosek::fusion::WorkStack> ws,
       monty::rc_ptr<::mosek::fusion::WorkStack> xs);
  virtual /* override */ std::string toString();
}; // class ExprMulParamDiagRight;

class ExprDotParam : public ::mosek::fusion::BaseExpression {
protected:
  ExprDotParam(p_ExprDotParam *_impl);

public:
  ExprDotParam(const ExprDotParam &) = delete;
  const ExprDotParam &operator=(const ExprDotParam &) = delete;
  friend class p_ExprDotParam;
  virtual ~ExprDotParam();
  virtual void destroy();
  typedef monty::rc_ptr<ExprDotParam> t;

  ExprDotParam(monty::rc_ptr<::mosek::fusion::Parameter> p,
               monty::rc_ptr<::mosek::fusion::Expression> e);
  virtual /* override */ void
  eval(monty::rc_ptr<::mosek::fusion::WorkStack> rs,
       monty::rc_ptr<::mosek::fusion::WorkStack> ws,
       monty::rc_ptr<::mosek::fusion::WorkStack> xs);
  virtual /* override */ std::string toString();
}; // class ExprDotParam;

class ExprMulParamElem : public ::mosek::fusion::BaseExpression {
protected:
  ExprMulParamElem(p_ExprMulParamElem *_impl);

public:
  ExprMulParamElem(const ExprMulParamElem &) = delete;
  const ExprMulParamElem &operator=(const ExprMulParamElem &) = delete;
  friend class p_ExprMulParamElem;
  virtual ~ExprMulParamElem();
  virtual void destroy();
  typedef monty::rc_ptr<ExprMulParamElem> t;

  ExprMulParamElem(monty::rc_ptr<::mosek::fusion::Parameter> p,
                   monty::rc_ptr<::mosek::fusion::Expression> e);
  virtual /* override */ void
  eval(monty::rc_ptr<::mosek::fusion::WorkStack> rs,
       monty::rc_ptr<::mosek::fusion::WorkStack> ws,
       monty::rc_ptr<::mosek::fusion::WorkStack> xs);
  virtual /* override */ std::string toString();
}; // class ExprMulParamElem;

class ExprMulParamRight : public ::mosek::fusion::BaseExpression {
protected:
  ExprMulParamRight(p_ExprMulParamRight *_impl);

public:
  ExprMulParamRight(const ExprMulParamRight &) = delete;
  const ExprMulParamRight &operator=(const ExprMulParamRight &) = delete;
  friend class p_ExprMulParamRight;
  virtual ~ExprMulParamRight();
  virtual void destroy();
  typedef monty::rc_ptr<ExprMulParamRight> t;

  ExprMulParamRight(monty::rc_ptr<::mosek::fusion::Expression> e,
                    monty::rc_ptr<::mosek::fusion::Parameter> p);
  virtual /* override */ void
  eval(monty::rc_ptr<::mosek::fusion::WorkStack> rs,
       monty::rc_ptr<::mosek::fusion::WorkStack> ws,
       monty::rc_ptr<::mosek::fusion::WorkStack> xs);
  virtual /* override */ std::string toString();
}; // class ExprMulParamRight;

class ExprMulParamLeft : public ::mosek::fusion::BaseExpression {
protected:
  ExprMulParamLeft(p_ExprMulParamLeft *_impl);

public:
  ExprMulParamLeft(const ExprMulParamLeft &) = delete;
  const ExprMulParamLeft &operator=(const ExprMulParamLeft &) = delete;
  friend class p_ExprMulParamLeft;
  virtual ~ExprMulParamLeft();
  virtual void destroy();
  typedef monty::rc_ptr<ExprMulParamLeft> t;

  ExprMulParamLeft(monty::rc_ptr<::mosek::fusion::Parameter> p,
                   monty::rc_ptr<::mosek::fusion::Expression> e);
  virtual /* override */ void
  eval(monty::rc_ptr<::mosek::fusion::WorkStack> rs,
       monty::rc_ptr<::mosek::fusion::WorkStack> ws,
       monty::rc_ptr<::mosek::fusion::WorkStack> xs);
  virtual /* override */ std::string toString();
}; // class ExprMulParamLeft;

class ExprOptimizeCode : public ::mosek::fusion::BaseExpression {
protected:
  ExprOptimizeCode(p_ExprOptimizeCode *_impl);

public:
  ExprOptimizeCode(const ExprOptimizeCode &) = delete;
  const ExprOptimizeCode &operator=(const ExprOptimizeCode &) = delete;
  friend class p_ExprOptimizeCode;
  virtual ~ExprOptimizeCode();
  virtual void destroy();
  typedef monty::rc_ptr<ExprOptimizeCode> t;

  ExprOptimizeCode(monty::rc_ptr<::mosek::fusion::Expression> expr);
  virtual /* override */ void
  eval(monty::rc_ptr<::mosek::fusion::WorkStack> rs,
       monty::rc_ptr<::mosek::fusion::WorkStack> ws,
       monty::rc_ptr<::mosek::fusion::WorkStack> xs);
  virtual /* override */ std::string toString();
}; // class ExprOptimizeCode;

class ExprCompress : public ::mosek::fusion::BaseExpression {
protected:
  ExprCompress(p_ExprCompress *_impl);

public:
  ExprCompress(const ExprCompress &) = delete;
  const ExprCompress &operator=(const ExprCompress &) = delete;
  friend class p_ExprCompress;
  virtual ~ExprCompress();
  virtual void destroy();
  typedef monty::rc_ptr<ExprCompress> t;

  ExprCompress(monty::rc_ptr<::mosek::fusion::Expression> expr);
  static void arg_sort(monty::rc_ptr<::mosek::fusion::WorkStack> ws,
                       monty::rc_ptr<::mosek::fusion::WorkStack> xs,
                       int32_t perm, int32_t nelem, int32_t nnz, int32_t ptr,
                       int32_t nidxs);
  static void merge_sort(int32_t origperm1, int32_t origperm2, int32_t nelem,
                         int32_t nnz, int32_t ptr_base, int32_t nidxs_base,
                         std::shared_ptr<monty::ndarray<int32_t, 1>> wi32,
                         std::shared_ptr<monty::ndarray<int64_t, 1>> wi64);
  virtual /* override */ void
  eval(monty::rc_ptr<::mosek::fusion::WorkStack> rs,
       monty::rc_ptr<::mosek::fusion::WorkStack> ws,
       monty::rc_ptr<::mosek::fusion::WorkStack> xs);
  virtual /* override */ std::string toString();
}; // class ExprCompress;

class ExprConst : public ::mosek::fusion::BaseExpression {
protected:
  ExprConst(p_ExprConst *_impl);

public:
  ExprConst(const ExprConst &) = delete;
  const ExprConst &operator=(const ExprConst &) = delete;
  friend class p_ExprConst;
  virtual ~ExprConst();
  virtual void destroy();
  typedef monty::rc_ptr<ExprConst> t;

  ExprConst(std::shared_ptr<monty::ndarray<int32_t, 1>> shape,
            std::shared_ptr<monty::ndarray<int64_t, 1>> sparsity,
            std::shared_ptr<monty::ndarray<double, 1>> bfix);
  ExprConst(std::shared_ptr<monty::ndarray<int32_t, 1>> shape,
            std::shared_ptr<monty::ndarray<int64_t, 1>> sparsity, double bfix);
  virtual /* override */ void
  eval(monty::rc_ptr<::mosek::fusion::WorkStack> rs,
       monty::rc_ptr<::mosek::fusion::WorkStack> ws,
       monty::rc_ptr<::mosek::fusion::WorkStack> xs);
  virtual /* override */ std::string toString();
}; // class ExprConst;

class ExprPick : public ::mosek::fusion::BaseExpression {
protected:
  ExprPick(p_ExprPick *_impl);

public:
  ExprPick(const ExprPick &) = delete;
  const ExprPick &operator=(const ExprPick &) = delete;
  friend class p_ExprPick;
  virtual ~ExprPick();
  virtual void destroy();
  typedef monty::rc_ptr<ExprPick> t;

  ExprPick(monty::rc_ptr<::mosek::fusion::Expression> expr,
           std::shared_ptr<monty::ndarray<int32_t, 2>> idxs);
  ExprPick(monty::rc_ptr<::mosek::fusion::Expression> expr,
           std::shared_ptr<monty::ndarray<int64_t, 1>> idxs);
  virtual /* override */ void
  eval(monty::rc_ptr<::mosek::fusion::WorkStack> rs,
       monty::rc_ptr<::mosek::fusion::WorkStack> ws,
       monty::rc_ptr<::mosek::fusion::WorkStack> xs);
  virtual /* override */ std::string toString();
}; // class ExprPick;

class ExprSlice : public ::mosek::fusion::BaseExpression {
protected:
  ExprSlice(p_ExprSlice *_impl);

public:
  ExprSlice(const ExprSlice &) = delete;
  const ExprSlice &operator=(const ExprSlice &) = delete;
  friend class p_ExprSlice;
  virtual ~ExprSlice();
  virtual void destroy();
  typedef monty::rc_ptr<ExprSlice> t;

  ExprSlice(monty::rc_ptr<::mosek::fusion::Expression> expr,
            std::shared_ptr<monty::ndarray<int32_t, 1>> first,
            std::shared_ptr<monty::ndarray<int32_t, 1>> last);
  virtual /* override */ void
  eval(monty::rc_ptr<::mosek::fusion::WorkStack> rs,
       monty::rc_ptr<::mosek::fusion::WorkStack> ws,
       monty::rc_ptr<::mosek::fusion::WorkStack> xs);
  virtual /* override */ std::string toString();
}; // class ExprSlice;

class ExprPermuteDims : public ::mosek::fusion::BaseExpression {
  ExprPermuteDims(std::shared_ptr<monty::ndarray<int32_t, 1>> perm,
                  monty::rc_ptr<::mosek::fusion::Expression> expr,
                  int32_t validated);

protected:
  ExprPermuteDims(p_ExprPermuteDims *_impl);

public:
  ExprPermuteDims(const ExprPermuteDims &) = delete;
  const ExprPermuteDims &operator=(const ExprPermuteDims &) = delete;
  friend class p_ExprPermuteDims;
  virtual ~ExprPermuteDims();
  virtual void destroy();
  typedef monty::rc_ptr<ExprPermuteDims> t;

  ExprPermuteDims(std::shared_ptr<monty::ndarray<int32_t, 1>> perm,
                  monty::rc_ptr<::mosek::fusion::Expression> expr);
  virtual /* override */ void
  eval(monty::rc_ptr<::mosek::fusion::WorkStack> rs,
       monty::rc_ptr<::mosek::fusion::WorkStack> ws,
       monty::rc_ptr<::mosek::fusion::WorkStack> xs);
}; // class ExprPermuteDims;

class ExprTranspose : public ::mosek::fusion::BaseExpression {
protected:
  ExprTranspose(p_ExprTranspose *_impl);

public:
  ExprTranspose(const ExprTranspose &) = delete;
  const ExprTranspose &operator=(const ExprTranspose &) = delete;
  friend class p_ExprTranspose;
  virtual ~ExprTranspose();
  virtual void destroy();
  typedef monty::rc_ptr<ExprTranspose> t;

  ExprTranspose(monty::rc_ptr<::mosek::fusion::Expression> expr);
  virtual /* override */ void
  eval(monty::rc_ptr<::mosek::fusion::WorkStack> rs,
       monty::rc_ptr<::mosek::fusion::WorkStack> ws,
       monty::rc_ptr<::mosek::fusion::WorkStack> xs);
  virtual /* override */ std::string toString();
}; // class ExprTranspose;

class ExprRepeat : public ::mosek::fusion::BaseExpression {
protected:
  ExprRepeat(p_ExprRepeat *_impl);

public:
  ExprRepeat(const ExprRepeat &) = delete;
  const ExprRepeat &operator=(const ExprRepeat &) = delete;
  friend class p_ExprRepeat;
  virtual ~ExprRepeat();
  virtual void destroy();
  typedef monty::rc_ptr<ExprRepeat> t;

  ExprRepeat(monty::rc_ptr<::mosek::fusion::Expression> expr, int32_t dim,
             int32_t n);
  virtual /* override */ void
  eval(monty::rc_ptr<::mosek::fusion::WorkStack> rs,
       monty::rc_ptr<::mosek::fusion::WorkStack> ws,
       monty::rc_ptr<::mosek::fusion::WorkStack> xs);
  virtual /* override */ std::string toString();
}; // class ExprRepeat;

class ExprStack : public ::mosek::fusion::BaseExpression {
protected:
  ExprStack(p_ExprStack *_impl);

public:
  ExprStack(const ExprStack &) = delete;
  const ExprStack &operator=(const ExprStack &) = delete;
  friend class p_ExprStack;
  virtual ~ExprStack();
  virtual void destroy();
  typedef monty::rc_ptr<ExprStack> t;

  ExprStack(std::shared_ptr<
                monty::ndarray<monty::rc_ptr<::mosek::fusion::Expression>, 1>>
                exprs,
            int32_t dim);
  virtual /* override */ void
  eval(monty::rc_ptr<::mosek::fusion::WorkStack> rs,
       monty::rc_ptr<::mosek::fusion::WorkStack> ws,
       monty::rc_ptr<::mosek::fusion::WorkStack> xs);
  virtual /* override */ std::string toString();
}; // class ExprStack;

class ExprInner : public ::mosek::fusion::BaseExpression {
protected:
  ExprInner(p_ExprInner *_impl);

public:
  ExprInner(const ExprInner &) = delete;
  const ExprInner &operator=(const ExprInner &) = delete;
  friend class p_ExprInner;
  virtual ~ExprInner();
  virtual void destroy();
  typedef monty::rc_ptr<ExprInner> t;

  ExprInner(monty::rc_ptr<::mosek::fusion::Expression> expr3,
            std::shared_ptr<monty::ndarray<int64_t, 1>> vsub3,
            std::shared_ptr<monty::ndarray<double, 1>> vcof3);
  ExprInner(monty::rc_ptr<::mosek::fusion::Expression> expr2,
            std::shared_ptr<monty::ndarray<double, 1>> vcof2);
  ExprInner(monty::rc_ptr<::mosek::fusion::Expression> expr1,
            std::shared_ptr<monty::ndarray<int32_t, 2>> vsub1,
            std::shared_ptr<monty::ndarray<double, 1>> vcof1);
  virtual /* override */ void
  eval(monty::rc_ptr<::mosek::fusion::WorkStack> rs,
       monty::rc_ptr<::mosek::fusion::WorkStack> ws,
       monty::rc_ptr<::mosek::fusion::WorkStack> xs);
  virtual /* override */ std::string toString();
}; // class ExprInner;

class ExprMulDiagRight : public ::mosek::fusion::BaseExpression {
  ExprMulDiagRight(int32_t mdim0, int32_t mdim1,
                   std::shared_ptr<monty::ndarray<int32_t, 1>> msubi,
                   std::shared_ptr<monty::ndarray<int32_t, 1>> msubj,
                   std::shared_ptr<monty::ndarray<double, 1>> mval,
                   monty::rc_ptr<::mosek::fusion::Expression> expr,
                   int32_t validated);

protected:
  ExprMulDiagRight(p_ExprMulDiagRight *_impl);

public:
  ExprMulDiagRight(const ExprMulDiagRight &) = delete;
  const ExprMulDiagRight &operator=(const ExprMulDiagRight &) = delete;
  friend class p_ExprMulDiagRight;
  virtual ~ExprMulDiagRight();
  virtual void destroy();
  typedef monty::rc_ptr<ExprMulDiagRight> t;

  ExprMulDiagRight(int32_t mdim0, int32_t mdim1,
                   std::shared_ptr<monty::ndarray<int32_t, 1>> msubi,
                   std::shared_ptr<monty::ndarray<int32_t, 1>> msubj,
                   std::shared_ptr<monty::ndarray<double, 1>> mval,
                   monty::rc_ptr<::mosek::fusion::Expression> expr);
  virtual /* override */ void
  eval(monty::rc_ptr<::mosek::fusion::WorkStack> rs,
       monty::rc_ptr<::mosek::fusion::WorkStack> ws,
       monty::rc_ptr<::mosek::fusion::WorkStack> xs);
  virtual /* override */ std::string toString();
}; // class ExprMulDiagRight;

class ExprMulDiagLeft : public ::mosek::fusion::BaseExpression {
  ExprMulDiagLeft(int32_t mdim0, int32_t mdim1,
                  std::shared_ptr<monty::ndarray<int32_t, 1>> msubi,
                  std::shared_ptr<monty::ndarray<int32_t, 1>> msubj,
                  std::shared_ptr<monty::ndarray<double, 1>> mval,
                  monty::rc_ptr<::mosek::fusion::Expression> expr,
                  int32_t validated);

protected:
  ExprMulDiagLeft(p_ExprMulDiagLeft *_impl);

public:
  ExprMulDiagLeft(const ExprMulDiagLeft &) = delete;
  const ExprMulDiagLeft &operator=(const ExprMulDiagLeft &) = delete;
  friend class p_ExprMulDiagLeft;
  virtual ~ExprMulDiagLeft();
  virtual void destroy();
  typedef monty::rc_ptr<ExprMulDiagLeft> t;

  ExprMulDiagLeft(int32_t mdim0, int32_t mdim1,
                  std::shared_ptr<monty::ndarray<int32_t, 1>> msubi,
                  std::shared_ptr<monty::ndarray<int32_t, 1>> msubj,
                  std::shared_ptr<monty::ndarray<double, 1>> mval,
                  monty::rc_ptr<::mosek::fusion::Expression> expr);
  virtual /* override */ void
  eval(monty::rc_ptr<::mosek::fusion::WorkStack> rs,
       monty::rc_ptr<::mosek::fusion::WorkStack> ws,
       monty::rc_ptr<::mosek::fusion::WorkStack> xs);
  virtual /* override */ std::string toString();
}; // class ExprMulDiagLeft;

class ExprMulElement : public ::mosek::fusion::BaseExpression {
  ExprMulElement(std::shared_ptr<monty::ndarray<double, 1>> cof,
                 std::shared_ptr<monty::ndarray<int64_t, 1>> msp,
                 monty::rc_ptr<::mosek::fusion::Expression> expr,
                 int32_t validated);

protected:
  ExprMulElement(p_ExprMulElement *_impl);

public:
  ExprMulElement(const ExprMulElement &) = delete;
  const ExprMulElement &operator=(const ExprMulElement &) = delete;
  friend class p_ExprMulElement;
  virtual ~ExprMulElement();
  virtual void destroy();
  typedef monty::rc_ptr<ExprMulElement> t;

  ExprMulElement(std::shared_ptr<monty::ndarray<double, 1>> mcof,
                 std::shared_ptr<monty::ndarray<int64_t, 1>> msp,
                 monty::rc_ptr<::mosek::fusion::Expression> expr);
  virtual /* override */ void
  eval(monty::rc_ptr<::mosek::fusion::WorkStack> rs,
       monty::rc_ptr<::mosek::fusion::WorkStack> ws,
       monty::rc_ptr<::mosek::fusion::WorkStack> xs);
  virtual /* override */ std::string toString();
}; // class ExprMulElement;

class ExprMulScalarConst : public ::mosek::fusion::BaseExpression {
protected:
  ExprMulScalarConst(p_ExprMulScalarConst *_impl);

public:
  ExprMulScalarConst(const ExprMulScalarConst &) = delete;
  const ExprMulScalarConst &operator=(const ExprMulScalarConst &) = delete;
  friend class p_ExprMulScalarConst;
  virtual ~ExprMulScalarConst();
  virtual void destroy();
  typedef monty::rc_ptr<ExprMulScalarConst> t;

  ExprMulScalarConst(double c, monty::rc_ptr<::mosek::fusion::Expression> expr);
  virtual /* override */ void
  eval(monty::rc_ptr<::mosek::fusion::WorkStack> rs,
       monty::rc_ptr<::mosek::fusion::WorkStack> ws,
       monty::rc_ptr<::mosek::fusion::WorkStack> xs);
  virtual /* override */ std::string toString();
}; // class ExprMulScalarConst;

class ExprScalarMul : public ::mosek::fusion::BaseExpression {
  ExprScalarMul(int32_t mdim0, int32_t mdim1,
                std::shared_ptr<monty::ndarray<int32_t, 1>> msubi,
                std::shared_ptr<monty::ndarray<int32_t, 1>> msubj,
                std::shared_ptr<monty::ndarray<double, 1>> mval,
                monty::rc_ptr<::mosek::fusion::Expression> expr,
                int32_t validated);

protected:
  ExprScalarMul(p_ExprScalarMul *_impl);

public:
  ExprScalarMul(const ExprScalarMul &) = delete;
  const ExprScalarMul &operator=(const ExprScalarMul &) = delete;
  friend class p_ExprScalarMul;
  virtual ~ExprScalarMul();
  virtual void destroy();
  typedef monty::rc_ptr<ExprScalarMul> t;

  ExprScalarMul(int32_t mdim0, int32_t mdim1,
                std::shared_ptr<monty::ndarray<int32_t, 1>> msubi,
                std::shared_ptr<monty::ndarray<int32_t, 1>> msubj,
                std::shared_ptr<monty::ndarray<double, 1>> mval,
                monty::rc_ptr<::mosek::fusion::Expression> expr);
  virtual /* override */ void
  eval(monty::rc_ptr<::mosek::fusion::WorkStack> rs,
       monty::rc_ptr<::mosek::fusion::WorkStack> ws,
       monty::rc_ptr<::mosek::fusion::WorkStack> xs);
  virtual /* override */ std::string toString();
}; // class ExprScalarMul;

class ExprMulRight : public ::mosek::fusion::BaseExpression {
  ExprMulRight(int32_t mdim0, int32_t mdim1,
               std::shared_ptr<monty::ndarray<int32_t, 1>> msubi,
               std::shared_ptr<monty::ndarray<int32_t, 1>> msubj,
               std::shared_ptr<monty::ndarray<double, 1>> mval,
               monty::rc_ptr<::mosek::fusion::Expression> expr,
               int32_t validated);

protected:
  ExprMulRight(p_ExprMulRight *_impl);

public:
  ExprMulRight(const ExprMulRight &) = delete;
  const ExprMulRight &operator=(const ExprMulRight &) = delete;
  friend class p_ExprMulRight;
  virtual ~ExprMulRight();
  virtual void destroy();
  typedef monty::rc_ptr<ExprMulRight> t;

  ExprMulRight(int32_t mdim0, int32_t mdim1,
               std::shared_ptr<monty::ndarray<int32_t, 1>> msubi,
               std::shared_ptr<monty::ndarray<int32_t, 1>> msubj,
               std::shared_ptr<monty::ndarray<double, 1>> mval,
               monty::rc_ptr<::mosek::fusion::Expression> expr);
  virtual /* override */ void
  eval(monty::rc_ptr<::mosek::fusion::WorkStack> rs,
       monty::rc_ptr<::mosek::fusion::WorkStack> ws,
       monty::rc_ptr<::mosek::fusion::WorkStack> xs);
  virtual /* override */ std::string toString();
}; // class ExprMulRight;

class ExprMulLeft : public ::mosek::fusion::BaseExpression {
  ExprMulLeft(int32_t mdim0, int32_t mdim1,
              std::shared_ptr<monty::ndarray<int32_t, 1>> msubi,
              std::shared_ptr<monty::ndarray<int32_t, 1>> msubj,
              std::shared_ptr<monty::ndarray<double, 1>> mval,
              monty::rc_ptr<::mosek::fusion::Expression> expr,
              int32_t validated);

protected:
  ExprMulLeft(p_ExprMulLeft *_impl);

public:
  ExprMulLeft(const ExprMulLeft &) = delete;
  const ExprMulLeft &operator=(const ExprMulLeft &) = delete;
  friend class p_ExprMulLeft;
  virtual ~ExprMulLeft();
  virtual void destroy();
  typedef monty::rc_ptr<ExprMulLeft> t;

  ExprMulLeft(int32_t mdim0, int32_t mdim1,
              std::shared_ptr<monty::ndarray<int32_t, 1>> msubi,
              std::shared_ptr<monty::ndarray<int32_t, 1>> msubj,
              std::shared_ptr<monty::ndarray<double, 1>> mval,
              monty::rc_ptr<::mosek::fusion::Expression> expr);
  virtual /* override */ void
  eval(monty::rc_ptr<::mosek::fusion::WorkStack> rs,
       monty::rc_ptr<::mosek::fusion::WorkStack> ws,
       monty::rc_ptr<::mosek::fusion::WorkStack> xs);
  virtual /* override */ std::string toString();
}; // class ExprMulLeft;

class ExprMulVar : public ::mosek::fusion::BaseExpression {
  ExprMulVar(bool left, int32_t mdimi, int32_t mdimj,
             std::shared_ptr<monty::ndarray<int32_t, 1>> msubi,
             std::shared_ptr<monty::ndarray<int32_t, 1>> msubj,
             std::shared_ptr<monty::ndarray<double, 1>> mcof,
             monty::rc_ptr<::mosek::fusion::Variable> x, int32_t unchecked_);

protected:
  ExprMulVar(p_ExprMulVar *_impl);

public:
  ExprMulVar(const ExprMulVar &) = delete;
  const ExprMulVar &operator=(const ExprMulVar &) = delete;
  friend class p_ExprMulVar;
  virtual ~ExprMulVar();
  virtual void destroy();
  typedef monty::rc_ptr<ExprMulVar> t;

  ExprMulVar(bool left, int32_t mdimi, int32_t mdimj,
             std::shared_ptr<monty::ndarray<int32_t, 1>> msubi,
             std::shared_ptr<monty::ndarray<int32_t, 1>> msubj,
             std::shared_ptr<monty::ndarray<double, 1>> mcof,
             monty::rc_ptr<::mosek::fusion::Variable> x);
  virtual /* override */ void
  eval(monty::rc_ptr<::mosek::fusion::WorkStack> rs,
       monty::rc_ptr<::mosek::fusion::WorkStack> ws,
       monty::rc_ptr<::mosek::fusion::WorkStack> xs);
  virtual void eval_right(monty::rc_ptr<::mosek::fusion::WorkStack> rs,
                          monty::rc_ptr<::mosek::fusion::WorkStack> ws,
                          monty::rc_ptr<::mosek::fusion::WorkStack> xs);
  virtual void eval_left(monty::rc_ptr<::mosek::fusion::WorkStack> rs,
                         monty::rc_ptr<::mosek::fusion::WorkStack> ws,
                         monty::rc_ptr<::mosek::fusion::WorkStack> xs);
  virtual /* override */ std::string toString();
}; // class ExprMulVar;

class ExprMulScalarVar : public ::mosek::fusion::BaseExpression {
  ExprMulScalarVar(int32_t mdimi, int32_t mdimj,
                   std::shared_ptr<monty::ndarray<int32_t, 1>> msubi,
                   std::shared_ptr<monty::ndarray<int32_t, 1>> msubj,
                   std::shared_ptr<monty::ndarray<double, 1>> mcof,
                   monty::rc_ptr<::mosek::fusion::Variable> x,
                   int32_t unchecked_);

protected:
  ExprMulScalarVar(p_ExprMulScalarVar *_impl);

public:
  ExprMulScalarVar(const ExprMulScalarVar &) = delete;
  const ExprMulScalarVar &operator=(const ExprMulScalarVar &) = delete;
  friend class p_ExprMulScalarVar;
  virtual ~ExprMulScalarVar();
  virtual void destroy();
  typedef monty::rc_ptr<ExprMulScalarVar> t;

  ExprMulScalarVar(int32_t mdimi, int32_t mdimj,
                   std::shared_ptr<monty::ndarray<int32_t, 1>> msubi,
                   std::shared_ptr<monty::ndarray<int32_t, 1>> msubj,
                   std::shared_ptr<monty::ndarray<double, 1>> mcof,
                   monty::rc_ptr<::mosek::fusion::Variable> x);
  virtual /* override */ void
  eval(monty::rc_ptr<::mosek::fusion::WorkStack> rs,
       monty::rc_ptr<::mosek::fusion::WorkStack> ws,
       monty::rc_ptr<::mosek::fusion::WorkStack> xs);
  virtual /* override */ std::string toString();
}; // class ExprMulScalarVar;

class ExprMulVarScalarConst : public ::mosek::fusion::BaseExpression {
protected:
  ExprMulVarScalarConst(p_ExprMulVarScalarConst *_impl);

public:
  ExprMulVarScalarConst(const ExprMulVarScalarConst &) = delete;
  const ExprMulVarScalarConst &
  operator=(const ExprMulVarScalarConst &) = delete;
  friend class p_ExprMulVarScalarConst;
  virtual ~ExprMulVarScalarConst();
  virtual void destroy();
  typedef monty::rc_ptr<ExprMulVarScalarConst> t;

  ExprMulVarScalarConst(monty::rc_ptr<::mosek::fusion::Variable> x, double c);
  virtual /* override */ void
  eval(monty::rc_ptr<::mosek::fusion::WorkStack> rs,
       monty::rc_ptr<::mosek::fusion::WorkStack> ws,
       monty::rc_ptr<::mosek::fusion::WorkStack> xs);
  virtual /* override */ std::string toString();
}; // class ExprMulVarScalarConst;

class ExprAdd : public ::mosek::fusion::BaseExpression {
protected:
  ExprAdd(p_ExprAdd *_impl);

public:
  ExprAdd(const ExprAdd &) = delete;
  const ExprAdd &operator=(const ExprAdd &) = delete;
  friend class p_ExprAdd;
  virtual ~ExprAdd();
  virtual void destroy();
  typedef monty::rc_ptr<ExprAdd> t;

  ExprAdd(monty::rc_ptr<::mosek::fusion::Expression> e1,
          monty::rc_ptr<::mosek::fusion::Expression> e2, double m1, double m2);
  virtual /* override */ void
  eval(monty::rc_ptr<::mosek::fusion::WorkStack> rs,
       monty::rc_ptr<::mosek::fusion::WorkStack> ws,
       monty::rc_ptr<::mosek::fusion::WorkStack> xs);
  virtual /* override */ std::string toString();
}; // class ExprAdd;

class ExprWSum : public ::mosek::fusion::BaseExpression {
protected:
  ExprWSum(p_ExprWSum *_impl);

public:
  ExprWSum(const ExprWSum &) = delete;
  const ExprWSum &operator=(const ExprWSum &) = delete;
  friend class p_ExprWSum;
  virtual ~ExprWSum();
  virtual void destroy();
  typedef monty::rc_ptr<ExprWSum> t;

  ExprWSum(std::shared_ptr<
               monty::ndarray<monty::rc_ptr<::mosek::fusion::Expression>, 1>>
               es,
           std::shared_ptr<monty::ndarray<double, 1>> w);
  virtual /* override */ void
  eval(monty::rc_ptr<::mosek::fusion::WorkStack> rs,
       monty::rc_ptr<::mosek::fusion::WorkStack> ws,
       monty::rc_ptr<::mosek::fusion::WorkStack> xs);
  virtual /* override */ std::string toString();
}; // class ExprWSum;

class ExprSumReduce : public ::mosek::fusion::BaseExpression {
protected:
  ExprSumReduce(p_ExprSumReduce *_impl);

public:
  ExprSumReduce(const ExprSumReduce &) = delete;
  const ExprSumReduce &operator=(const ExprSumReduce &) = delete;
  friend class p_ExprSumReduce;
  virtual ~ExprSumReduce();
  virtual void destroy();
  typedef monty::rc_ptr<ExprSumReduce> t;

  ExprSumReduce(int32_t dim, monty::rc_ptr<::mosek::fusion::Expression> expr);
  virtual /* override */ void
  eval(monty::rc_ptr<::mosek::fusion::WorkStack> rs,
       monty::rc_ptr<::mosek::fusion::WorkStack> ws,
       monty::rc_ptr<::mosek::fusion::WorkStack> xs);
  virtual /* override */ std::string toString();
}; // class ExprSumReduce;

class ExprScaleVecPSD : public ::mosek::fusion::BaseExpression {
protected:
  ExprScaleVecPSD(p_ExprScaleVecPSD *_impl);

public:
  ExprScaleVecPSD(const ExprScaleVecPSD &) = delete;
  const ExprScaleVecPSD &operator=(const ExprScaleVecPSD &) = delete;
  friend class p_ExprScaleVecPSD;
  virtual ~ExprScaleVecPSD();
  virtual void destroy();
  typedef monty::rc_ptr<ExprScaleVecPSD> t;

  ExprScaleVecPSD(int32_t dim0, int32_t dim1,
                  monty::rc_ptr<::mosek::fusion::BaseExpression> expr);
  virtual /* override */ void
  eval(monty::rc_ptr<::mosek::fusion::WorkStack> rs,
       monty::rc_ptr<::mosek::fusion::WorkStack> ws,
       monty::rc_ptr<::mosek::fusion::WorkStack> xs);
}; // class ExprScaleVecPSD;

class ExprDenseTril : public ::mosek::fusion::BaseExpression {
  ExprDenseTril(int32_t dim0, int32_t dim1,
                monty::rc_ptr<::mosek::fusion::Expression> expr,
                int32_t unchecked_);

protected:
  ExprDenseTril(p_ExprDenseTril *_impl);

public:
  ExprDenseTril(const ExprDenseTril &) = delete;
  const ExprDenseTril &operator=(const ExprDenseTril &) = delete;
  friend class p_ExprDenseTril;
  virtual ~ExprDenseTril();
  virtual void destroy();
  typedef monty::rc_ptr<ExprDenseTril> t;

  ExprDenseTril(int32_t dim0_, int32_t dim1_,
                monty::rc_ptr<::mosek::fusion::Expression> expr);
  virtual /* override */ void
  eval(monty::rc_ptr<::mosek::fusion::WorkStack> rs,
       monty::rc_ptr<::mosek::fusion::WorkStack> ws,
       monty::rc_ptr<::mosek::fusion::WorkStack> xs);
  virtual /* override */ std::string toString();
}; // class ExprDenseTril;

class ExprDense : public ::mosek::fusion::BaseExpression {
protected:
  ExprDense(p_ExprDense *_impl);

public:
  ExprDense(const ExprDense &) = delete;
  const ExprDense &operator=(const ExprDense &) = delete;
  friend class p_ExprDense;
  virtual ~ExprDense();
  virtual void destroy();
  typedef monty::rc_ptr<ExprDense> t;

  ExprDense(monty::rc_ptr<::mosek::fusion::Expression> expr);
  virtual /* override */ void
  eval(monty::rc_ptr<::mosek::fusion::WorkStack> rs,
       monty::rc_ptr<::mosek::fusion::WorkStack> ws,
       monty::rc_ptr<::mosek::fusion::WorkStack> xs);
  virtual /* override */ std::string toString();
}; // class ExprDense;

class ExprSymmetrize : public ::mosek::fusion::BaseExpression {
protected:
  ExprSymmetrize(p_ExprSymmetrize *_impl);

public:
  ExprSymmetrize(const ExprSymmetrize &) = delete;
  const ExprSymmetrize &operator=(const ExprSymmetrize &) = delete;
  friend class p_ExprSymmetrize;
  virtual ~ExprSymmetrize();
  virtual void destroy();
  typedef monty::rc_ptr<ExprSymmetrize> t;

  ExprSymmetrize(int32_t dim0, int32_t dim1,
                 monty::rc_ptr<::mosek::fusion::Expression> expr,
                 int32_t unchecked_);
  ExprSymmetrize(int32_t dim0_, int32_t dim1_,
                 monty::rc_ptr<::mosek::fusion::Expression> expr);
  virtual /* override */ void
  eval(monty::rc_ptr<::mosek::fusion::WorkStack> rs,
       monty::rc_ptr<::mosek::fusion::WorkStack> ws,
       monty::rc_ptr<::mosek::fusion::WorkStack> xs);
  virtual /* override */ std::string toString();
}; // class ExprSymmetrize;

class ExprCondense : public ::mosek::fusion::BaseExpression {
protected:
  ExprCondense(p_ExprCondense *_impl);

public:
  ExprCondense(const ExprCondense &) = delete;
  const ExprCondense &operator=(const ExprCondense &) = delete;
  friend class p_ExprCondense;
  virtual ~ExprCondense();
  virtual void destroy();
  typedef monty::rc_ptr<ExprCondense> t;

  ExprCondense(monty::rc_ptr<::mosek::fusion::Expression> expr);
  virtual /* override */ void
  eval(monty::rc_ptr<::mosek::fusion::WorkStack> rs,
       monty::rc_ptr<::mosek::fusion::WorkStack> ws,
       monty::rc_ptr<::mosek::fusion::WorkStack> xs);
  virtual /* override */ std::string toString();
}; // class ExprCondense;

class ExprFromVar : public ::mosek::fusion::BaseExpression {
protected:
  ExprFromVar(p_ExprFromVar *_impl);

public:
  ExprFromVar(const ExprFromVar &) = delete;
  const ExprFromVar &operator=(const ExprFromVar &) = delete;
  friend class p_ExprFromVar;
  virtual ~ExprFromVar();
  virtual void destroy();
  typedef monty::rc_ptr<ExprFromVar> t;

  ExprFromVar(monty::rc_ptr<::mosek::fusion::Variable> x);
  virtual /* override */ void
  eval(monty::rc_ptr<::mosek::fusion::WorkStack> rs,
       monty::rc_ptr<::mosek::fusion::WorkStack> ws,
       monty::rc_ptr<::mosek::fusion::WorkStack> xs);
  virtual /* override */ std::string toString();
}; // class ExprFromVar;

class ExprReshape : public ::mosek::fusion::BaseExpression {
protected:
  ExprReshape(p_ExprReshape *_impl);

public:
  ExprReshape(const ExprReshape &) = delete;
  const ExprReshape &operator=(const ExprReshape &) = delete;
  friend class p_ExprReshape;
  virtual ~ExprReshape();
  virtual void destroy();
  typedef monty::rc_ptr<ExprReshape> t;

  ExprReshape(std::shared_ptr<monty::ndarray<int32_t, 1>> shape,
              monty::rc_ptr<::mosek::fusion::Expression> e);
  virtual /* override */ void
  eval(monty::rc_ptr<::mosek::fusion::WorkStack> rs,
       monty::rc_ptr<::mosek::fusion::WorkStack> ws,
       monty::rc_ptr<::mosek::fusion::WorkStack> xs);
  virtual /* override */ std::string toString();
}; // class ExprReshape;

class Expr : public ::mosek::fusion::BaseExpression {
  Expr(std::shared_ptr<monty::ndarray<int64_t, 1>> ptrb,
       std::shared_ptr<monty::ndarray<int64_t, 1>> subj,
       std::shared_ptr<monty::ndarray<double, 1>> cof,
       std::shared_ptr<monty::ndarray<double, 1>> bfix,
       std::shared_ptr<monty::ndarray<int32_t, 1>> shp,
       std::shared_ptr<monty::ndarray<int64_t, 1>> inst, int32_t unchecked_);
  Expr(monty::rc_ptr<::mosek::fusion::Expression> e);

protected:
  Expr(p_Expr *_impl);

public:
  Expr(const Expr &) = delete;
  const Expr &operator=(const Expr &) = delete;
  friend class p_Expr;
  virtual ~Expr();
  virtual void destroy();
  typedef monty::rc_ptr<Expr> t;

  Expr(std::shared_ptr<monty::ndarray<int64_t, 1>> ptrb,
       std::shared_ptr<monty::ndarray<int64_t, 1>> subj,
       std::shared_ptr<monty::ndarray<double, 1>> cof,
       std::shared_ptr<monty::ndarray<double, 1>> bfix,
       std::shared_ptr<monty::ndarray<int32_t, 1>> shape,
       std::shared_ptr<monty::ndarray<int64_t, 1>> inst);
  static monty::rc_ptr<::mosek::fusion::Expression>
  condense(monty::rc_ptr<::mosek::fusion::Expression> e);
  static monty::rc_ptr<::mosek::fusion::Expression>
  flatten(monty::rc_ptr<::mosek::fusion::Expression> e);
  static monty::rc_ptr<::mosek::fusion::Expression>
  reshape(monty::rc_ptr<::mosek::fusion::Expression> e, int32_t dimi,
          int32_t dimj);
  static monty::rc_ptr<::mosek::fusion::Expression>
  reshape(monty::rc_ptr<::mosek::fusion::Expression> e, int32_t size);
  static monty::rc_ptr<::mosek::fusion::Expression>
  reshape(monty::rc_ptr<::mosek::fusion::Expression> e,
          std::shared_ptr<monty::ndarray<int32_t, 1>> newshape);
  static monty::rc_ptr<::mosek::fusion::Expression>
  zeros(std::shared_ptr<monty::ndarray<int32_t, 1>> shp);
  static monty::rc_ptr<::mosek::fusion::Expression> zeros(int32_t size);
  static monty::rc_ptr<::mosek::fusion::Expression> ones();
  static monty::rc_ptr<::mosek::fusion::Expression>
  ones(std::shared_ptr<monty::ndarray<int32_t, 1>> shp,
       std::shared_ptr<monty::ndarray<int32_t, 2>> sparsity);
  static monty::rc_ptr<::mosek::fusion::Expression>
  ones(std::shared_ptr<monty::ndarray<int32_t, 1>> shp);
  static monty::rc_ptr<::mosek::fusion::Expression> ones(int32_t size);
  static monty::rc_ptr<::mosek::fusion::Expression>
  constTerm(monty::rc_ptr<::mosek::fusion::NDSparseArray> nda);
  static monty::rc_ptr<::mosek::fusion::Expression>
  constTerm(monty::rc_ptr<::mosek::fusion::Matrix> m);
  static monty::rc_ptr<::mosek::fusion::Expression> constTerm(double val);
  static monty::rc_ptr<::mosek::fusion::Expression>
  constTerm(std::shared_ptr<monty::ndarray<int32_t, 1>> shp,
            std::shared_ptr<monty::ndarray<int32_t, 2>> sparsity, double val);
  static monty::rc_ptr<::mosek::fusion::Expression>
  constTerm(std::shared_ptr<monty::ndarray<int32_t, 1>> shp,
            std::shared_ptr<monty::ndarray<int32_t, 2>> sparsity,
            std::shared_ptr<monty::ndarray<double, 1>> vals1);
  static monty::rc_ptr<::mosek::fusion::Expression>
  constTerm(std::shared_ptr<monty::ndarray<int32_t, 1>> shp, double val);
  static monty::rc_ptr<::mosek::fusion::Expression> constTerm(int32_t size,
                                                              double val);
  static monty::rc_ptr<::mosek::fusion::Expression>
  constTerm(std::shared_ptr<monty::ndarray<double, 2>> vals2);
  static monty::rc_ptr<::mosek::fusion::Expression>
  constTerm(std::shared_ptr<monty::ndarray<double, 1>> vals1);
  static monty::rc_ptr<::mosek::fusion::Expression>
  sum(monty::rc_ptr<::mosek::fusion::Expression> expr, int32_t dim);
  static monty::rc_ptr<::mosek::fusion::Expression>
  sum(monty::rc_ptr<::mosek::fusion::Expression> expr);
  static monty::rc_ptr<::mosek::fusion::Expression>
  neg(monty::rc_ptr<::mosek::fusion::Expression> e);
  static monty::rc_ptr<::mosek::fusion::Expression>
  mulDiag(monty::rc_ptr<::mosek::fusion::Variable> v,
          monty::rc_ptr<::mosek::fusion::Parameter> p);
  static monty::rc_ptr<::mosek::fusion::Expression>
  mulDiag(monty::rc_ptr<::mosek::fusion::Parameter> p,
          monty::rc_ptr<::mosek::fusion::Variable> v);
  static monty::rc_ptr<::mosek::fusion::Expression>
  mulDiag(monty::rc_ptr<::mosek::fusion::Expression> expr,
          monty::rc_ptr<::mosek::fusion::Parameter> p);
  static monty::rc_ptr<::mosek::fusion::Expression>
  mulDiag(monty::rc_ptr<::mosek::fusion::Parameter> p,
          monty::rc_ptr<::mosek::fusion::Expression> expr);
  static monty::rc_ptr<::mosek::fusion::Expression>
  mulDiag(monty::rc_ptr<::mosek::fusion::Variable> v,
          monty::rc_ptr<::mosek::fusion::Matrix> mx);
  static monty::rc_ptr<::mosek::fusion::Expression>
  mulDiag(monty::rc_ptr<::mosek::fusion::Matrix> mx,
          monty::rc_ptr<::mosek::fusion::Variable> v);
  static monty::rc_ptr<::mosek::fusion::Expression>
  mulDiag(monty::rc_ptr<::mosek::fusion::Expression> expr,
          monty::rc_ptr<::mosek::fusion::Matrix> mx);
  static monty::rc_ptr<::mosek::fusion::Expression>
  mulDiag(monty::rc_ptr<::mosek::fusion::Matrix> mx,
          monty::rc_ptr<::mosek::fusion::Expression> expr);
  static monty::rc_ptr<::mosek::fusion::Expression>
  mulDiag(monty::rc_ptr<::mosek::fusion::Variable> v,
          std::shared_ptr<monty::ndarray<double, 2>> a);
  static monty::rc_ptr<::mosek::fusion::Expression>
  mulDiag(monty::rc_ptr<::mosek::fusion::Expression> expr,
          std::shared_ptr<monty::ndarray<double, 2>> a);
  static monty::rc_ptr<::mosek::fusion::Expression>
  mulDiag(std::shared_ptr<monty::ndarray<double, 2>> a,
          monty::rc_ptr<::mosek::fusion::Variable> v);
  static monty::rc_ptr<::mosek::fusion::Expression>
  mulDiag(std::shared_ptr<monty::ndarray<double, 2>> a,
          monty::rc_ptr<::mosek::fusion::Expression> expr);
  static monty::rc_ptr<::mosek::fusion::Expression>
  mul(monty::rc_ptr<::mosek::fusion::Expression> expr, double c);
  static monty::rc_ptr<::mosek::fusion::Expression>
  mul(double c, monty::rc_ptr<::mosek::fusion::Expression> expr);
  static monty::rc_ptr<::mosek::fusion::Expression>
  mul(monty::rc_ptr<::mosek::fusion::Expression> expr,
      std::shared_ptr<monty::ndarray<double, 1>> a);
  static monty::rc_ptr<::mosek::fusion::Expression>
  mul(std::shared_ptr<monty::ndarray<double, 1>> a,
      monty::rc_ptr<::mosek::fusion::Expression> expr);
  static monty::rc_ptr<::mosek::fusion::Expression>
  mul(monty::rc_ptr<::mosek::fusion::Expression> expr,
      std::shared_ptr<monty::ndarray<double, 2>> a);
  static monty::rc_ptr<::mosek::fusion::Expression>
  mul(std::shared_ptr<monty::ndarray<double, 2>> a,
      monty::rc_ptr<::mosek::fusion::Expression> expr);
  static monty::rc_ptr<::mosek::fusion::Expression>
  mul(monty::rc_ptr<::mosek::fusion::Expression> expr,
      monty::rc_ptr<::mosek::fusion::Matrix> mx);
  static monty::rc_ptr<::mosek::fusion::Expression>
  mul(monty::rc_ptr<::mosek::fusion::Matrix> mx,
      monty::rc_ptr<::mosek::fusion::Expression> expr);
  static monty::rc_ptr<::mosek::fusion::Expression>
  mul(monty::rc_ptr<::mosek::fusion::Variable> v,
      std::shared_ptr<monty::ndarray<double, 2>> mx);
  static monty::rc_ptr<::mosek::fusion::Expression>
  mul(std::shared_ptr<monty::ndarray<double, 2>> mx,
      monty::rc_ptr<::mosek::fusion::Variable> v);
  static monty::rc_ptr<::mosek::fusion::Expression>
  mul(monty::rc_ptr<::mosek::fusion::Variable> v,
      monty::rc_ptr<::mosek::fusion::Matrix> mx);
  static monty::rc_ptr<::mosek::fusion::Expression>
  mul(monty::rc_ptr<::mosek::fusion::Matrix> mx,
      monty::rc_ptr<::mosek::fusion::Variable> v);
  static monty::rc_ptr<::mosek::fusion::Expression>
  mul(monty::rc_ptr<::mosek::fusion::Expression> expr,
      monty::rc_ptr<::mosek::fusion::Parameter> p);
  static monty::rc_ptr<::mosek::fusion::Expression>
  mul(monty::rc_ptr<::mosek::fusion::Parameter> p,
      monty::rc_ptr<::mosek::fusion::Expression> expr);
  static monty::rc_ptr<::mosek::fusion::Expression>
  dot(monty::rc_ptr<::mosek::fusion::Expression> e,
      monty::rc_ptr<::mosek::fusion::Matrix> m);
  static monty::rc_ptr<::mosek::fusion::Expression>
  dot(monty::rc_ptr<::mosek::fusion::Expression> e,
      std::shared_ptr<monty::ndarray<double, 2>> c2);
  static monty::rc_ptr<::mosek::fusion::Expression>
  dot(monty::rc_ptr<::mosek::fusion::Expression> e,
      monty::rc_ptr<::mosek::fusion::NDSparseArray> nda);
  static monty::rc_ptr<::mosek::fusion::Expression>
  dot(monty::rc_ptr<::mosek::fusion::Expression> e,
      std::shared_ptr<monty::ndarray<double, 1>> c1);
  static monty::rc_ptr<::mosek::fusion::Expression>
  dot(monty::rc_ptr<::mosek::fusion::Matrix> m,
      monty::rc_ptr<::mosek::fusion::Expression> e);
  static monty::rc_ptr<::mosek::fusion::Expression>
  dot(monty::rc_ptr<::mosek::fusion::NDSparseArray> nda,
      monty::rc_ptr<::mosek::fusion::Expression> e);
  static monty::rc_ptr<::mosek::fusion::Expression>
  dot(std::shared_ptr<monty::ndarray<double, 2>> c2,
      monty::rc_ptr<::mosek::fusion::Expression> e);
  static monty::rc_ptr<::mosek::fusion::Expression>
  dot(std::shared_ptr<monty::ndarray<double, 1>> c1,
      monty::rc_ptr<::mosek::fusion::Expression> e);
  static monty::rc_ptr<::mosek::fusion::Expression>
  dot(monty::rc_ptr<::mosek::fusion::Expression> e,
      monty::rc_ptr<::mosek::fusion::Parameter> p);
  static monty::rc_ptr<::mosek::fusion::Expression>
  dot(monty::rc_ptr<::mosek::fusion::Parameter> p,
      monty::rc_ptr<::mosek::fusion::Expression> e);
  static monty::rc_ptr<::mosek::fusion::Expression>
  outer(monty::rc_ptr<::mosek::fusion::Parameter> p,
        monty::rc_ptr<::mosek::fusion::Expression> e);
  static monty::rc_ptr<::mosek::fusion::Expression>
  outer(monty::rc_ptr<::mosek::fusion::Expression> e,
        monty::rc_ptr<::mosek::fusion::Parameter> p);
  static monty::rc_ptr<::mosek::fusion::Expression>
  outer(monty::rc_ptr<::mosek::fusion::Matrix> m,
        monty::rc_ptr<::mosek::fusion::Expression> e);
  static monty::rc_ptr<::mosek::fusion::Expression>
  outer(monty::rc_ptr<::mosek::fusion::Expression> e,
        monty::rc_ptr<::mosek::fusion::Matrix> m);
  static monty::rc_ptr<::mosek::fusion::Expression>
  outer(std::shared_ptr<monty::ndarray<double, 1>> a,
        monty::rc_ptr<::mosek::fusion::Expression> e);
  static monty::rc_ptr<::mosek::fusion::Expression>
  outer(monty::rc_ptr<::mosek::fusion::Expression> e,
        std::shared_ptr<monty::ndarray<double, 1>> a);
  static monty::rc_ptr<::mosek::fusion::Expression>
  stack(std::shared_ptr<
        monty::ndarray<std::shared_ptr<monty::ndarray<
                           monty::rc_ptr<::mosek::fusion::Expression>, 1>>,
                       1>>
            exprs);
  static monty::rc_ptr<::mosek::fusion::Expression> vstack(double a1, double a2,
                                                           double a3);
  static monty::rc_ptr<::mosek::fusion::Expression>
  vstack(double a1, double a2, monty::rc_ptr<::mosek::fusion::Expression> e3);
  static monty::rc_ptr<::mosek::fusion::Expression>
  vstack(double a1, monty::rc_ptr<::mosek::fusion::Expression> e2, double a3);
  static monty::rc_ptr<::mosek::fusion::Expression>
  vstack(double a1, monty::rc_ptr<::mosek::fusion::Expression> e2,
         monty::rc_ptr<::mosek::fusion::Expression> e3);
  static monty::rc_ptr<::mosek::fusion::Expression>
  vstack(monty::rc_ptr<::mosek::fusion::Expression> e1, double a2, double a3);
  static monty::rc_ptr<::mosek::fusion::Expression>
  vstack(monty::rc_ptr<::mosek::fusion::Expression> e1, double a2,
         monty::rc_ptr<::mosek::fusion::Expression> e3);
  static monty::rc_ptr<::mosek::fusion::Expression>
  vstack(monty::rc_ptr<::mosek::fusion::Expression> e1,
         monty::rc_ptr<::mosek::fusion::Expression> e2, double a3);
  static monty::rc_ptr<::mosek::fusion::Expression>
  vstack(monty::rc_ptr<::mosek::fusion::Expression> e1,
         monty::rc_ptr<::mosek::fusion::Expression> e2,
         monty::rc_ptr<::mosek::fusion::Expression> e3);
  static monty::rc_ptr<::mosek::fusion::Expression>
  vstack(double a1, monty::rc_ptr<::mosek::fusion::Expression> e2);
  static monty::rc_ptr<::mosek::fusion::Expression>
  vstack(monty::rc_ptr<::mosek::fusion::Expression> e1, double a2);
  static monty::rc_ptr<::mosek::fusion::Expression>
  vstack(monty::rc_ptr<::mosek::fusion::Expression> e1,
         monty::rc_ptr<::mosek::fusion::Expression> e2);
  static monty::rc_ptr<::mosek::fusion::Expression>
  vstack(std::shared_ptr<
         monty::ndarray<monty::rc_ptr<::mosek::fusion::Expression>, 1>>
             exprs);
  static monty::rc_ptr<::mosek::fusion::Expression>
  hstack(monty::rc_ptr<::mosek::fusion::Expression> e1,
         monty::rc_ptr<::mosek::fusion::Expression> e2,
         monty::rc_ptr<::mosek::fusion::Expression> e3);
  static monty::rc_ptr<::mosek::fusion::Expression>
  hstack(monty::rc_ptr<::mosek::fusion::Expression> e1,
         monty::rc_ptr<::mosek::fusion::Expression> e2, double a3);
  static monty::rc_ptr<::mosek::fusion::Expression>
  hstack(monty::rc_ptr<::mosek::fusion::Expression> e1, double a2,
         monty::rc_ptr<::mosek::fusion::Expression> e3);
  static monty::rc_ptr<::mosek::fusion::Expression>
  hstack(monty::rc_ptr<::mosek::fusion::Expression> e1, double a2, double a3);
  static monty::rc_ptr<::mosek::fusion::Expression>
  hstack(double a1, monty::rc_ptr<::mosek::fusion::Expression> e2,
         monty::rc_ptr<::mosek::fusion::Expression> e3);
  static monty::rc_ptr<::mosek::fusion::Expression>
  hstack(double a1, monty::rc_ptr<::mosek::fusion::Expression> e2, double a3);
  static monty::rc_ptr<::mosek::fusion::Expression>
  hstack(double a1, double a2, monty::rc_ptr<::mosek::fusion::Expression> e3);
  static monty::rc_ptr<::mosek::fusion::Expression>
  hstack(double a1, monty::rc_ptr<::mosek::fusion::Expression> e2);
  static monty::rc_ptr<::mosek::fusion::Expression>
  hstack(monty::rc_ptr<::mosek::fusion::Expression> e1, double a2);
  static monty::rc_ptr<::mosek::fusion::Expression>
  hstack(monty::rc_ptr<::mosek::fusion::Expression> e1,
         monty::rc_ptr<::mosek::fusion::Expression> e2);
  static monty::rc_ptr<::mosek::fusion::Expression>
  hstack(std::shared_ptr<
         monty::ndarray<monty::rc_ptr<::mosek::fusion::Expression>, 1>>
             exprs);
  static monty::rc_ptr<::mosek::fusion::Expression>
  stack(int32_t dim, monty::rc_ptr<::mosek::fusion::Expression> e1,
        monty::rc_ptr<::mosek::fusion::Expression> e2,
        monty::rc_ptr<::mosek::fusion::Expression> e3);
  static monty::rc_ptr<::mosek::fusion::Expression>
  stack(int32_t dim, monty::rc_ptr<::mosek::fusion::Expression> e1,
        monty::rc_ptr<::mosek::fusion::Expression> e2, double a3);
  static monty::rc_ptr<::mosek::fusion::Expression>
  stack(int32_t dim, monty::rc_ptr<::mosek::fusion::Expression> e1, double a2,
        monty::rc_ptr<::mosek::fusion::Expression> e3);
  static monty::rc_ptr<::mosek::fusion::Expression>
  stack(int32_t dim, monty::rc_ptr<::mosek::fusion::Expression> e1, double a2,
        double a3);
  static monty::rc_ptr<::mosek::fusion::Expression>
  stack(int32_t dim, double a1, monty::rc_ptr<::mosek::fusion::Expression> e2,
        monty::rc_ptr<::mosek::fusion::Expression> e3);
  static monty::rc_ptr<::mosek::fusion::Expression>
  stack(int32_t dim, double a1, monty::rc_ptr<::mosek::fusion::Expression> e2,
        double a3);
  static monty::rc_ptr<::mosek::fusion::Expression>
  stack(int32_t dim, double a1, double a2,
        monty::rc_ptr<::mosek::fusion::Expression> e1);
  static monty::rc_ptr<::mosek::fusion::Expression>
  stack(int32_t dim, double a1, monty::rc_ptr<::mosek::fusion::Expression> e2);
  static monty::rc_ptr<::mosek::fusion::Expression>
  stack(int32_t dim, monty::rc_ptr<::mosek::fusion::Expression> e1, double a2);
  static monty::rc_ptr<::mosek::fusion::Expression>
  stack(int32_t dim, monty::rc_ptr<::mosek::fusion::Expression> e1,
        monty::rc_ptr<::mosek::fusion::Expression> e2);
  static monty::rc_ptr<::mosek::fusion::Expression>
  stack(int32_t dim,
        std::shared_ptr<
            monty::ndarray<monty::rc_ptr<::mosek::fusion::Expression>, 1>>
            exprs);
  static monty::rc_ptr<::mosek::fusion::Expression>
  repeat(monty::rc_ptr<::mosek::fusion::Variable> x, int32_t n, int32_t d);
  static monty::rc_ptr<::mosek::fusion::Expression>
  repeat(monty::rc_ptr<::mosek::fusion::Expression> e, int32_t n, int32_t d);
  static monty::rc_ptr<::mosek::fusion::Expression>
  add(std::shared_ptr<
      monty::ndarray<monty::rc_ptr<::mosek::fusion::Expression>, 1>>
          exps);
  static monty::rc_ptr<::mosek::fusion::Expression>
  add(std::shared_ptr<
      monty::ndarray<monty::rc_ptr<::mosek::fusion::Variable>, 1>>
          vs);
  static monty::rc_ptr<::mosek::fusion::Expression>
  transpose(monty::rc_ptr<::mosek::fusion::Expression> e);
  static monty::rc_ptr<::mosek::fusion::Expression>
  mulElm(monty::rc_ptr<::mosek::fusion::Matrix> m,
         monty::rc_ptr<::mosek::fusion::Expression> expr);
  static monty::rc_ptr<::mosek::fusion::Expression>
  mulElm(monty::rc_ptr<::mosek::fusion::NDSparseArray> spm,
         monty::rc_ptr<::mosek::fusion::Expression> expr);
  static monty::rc_ptr<::mosek::fusion::Expression>
  mulElm(std::shared_ptr<monty::ndarray<double, 2>> a2,
         monty::rc_ptr<::mosek::fusion::Expression> expr);
  static monty::rc_ptr<::mosek::fusion::Expression>
  mulElm(std::shared_ptr<monty::ndarray<double, 1>> a1,
         monty::rc_ptr<::mosek::fusion::Expression> expr);
  static monty::rc_ptr<::mosek::fusion::Expression>
  mulElm(monty::rc_ptr<::mosek::fusion::Expression> expr,
         monty::rc_ptr<::mosek::fusion::Matrix> m);
  static monty::rc_ptr<::mosek::fusion::Expression>
  mulElm(monty::rc_ptr<::mosek::fusion::Expression> expr,
         std::shared_ptr<monty::ndarray<double, 2>> a2);
  static monty::rc_ptr<::mosek::fusion::Expression>
  mulElm(monty::rc_ptr<::mosek::fusion::Expression> expr,
         std::shared_ptr<monty::ndarray<double, 1>> a1);
  static monty::rc_ptr<::mosek::fusion::Expression>
  mulElm(monty::rc_ptr<::mosek::fusion::Expression> expr,
         monty::rc_ptr<::mosek::fusion::NDSparseArray> spm);
  static monty::rc_ptr<::mosek::fusion::Expression>
  mulElm(monty::rc_ptr<::mosek::fusion::Parameter> p,
         monty::rc_ptr<::mosek::fusion::Expression> expr);
  static monty::rc_ptr<::mosek::fusion::Expression>
  mulElm(monty::rc_ptr<::mosek::fusion::Expression> expr,
         monty::rc_ptr<::mosek::fusion::Parameter> p);
  static monty::rc_ptr<::mosek::fusion::Expression>
  sub(monty::rc_ptr<::mosek::fusion::NDSparseArray> n,
      monty::rc_ptr<::mosek::fusion::Expression> e2);
  static monty::rc_ptr<::mosek::fusion::Expression>
  sub(monty::rc_ptr<::mosek::fusion::Expression> e1,
      monty::rc_ptr<::mosek::fusion::NDSparseArray> n);
  static monty::rc_ptr<::mosek::fusion::Expression>
  sub(monty::rc_ptr<::mosek::fusion::Matrix> m,
      monty::rc_ptr<::mosek::fusion::Expression> e2);
  static monty::rc_ptr<::mosek::fusion::Expression>
  sub(monty::rc_ptr<::mosek::fusion::Expression> e1,
      monty::rc_ptr<::mosek::fusion::Matrix> m);
  static monty::rc_ptr<::mosek::fusion::Expression>
  sub(double c, monty::rc_ptr<::mosek::fusion::Expression> e2);
  static monty::rc_ptr<::mosek::fusion::Expression>
  sub(monty::rc_ptr<::mosek::fusion::Expression> e1, double c);
  static monty::rc_ptr<::mosek::fusion::Expression>
  sub(std::shared_ptr<monty::ndarray<double, 2>> a2,
      monty::rc_ptr<::mosek::fusion::Expression> e2);
  static monty::rc_ptr<::mosek::fusion::Expression>
  sub(std::shared_ptr<monty::ndarray<double, 1>> a1,
      monty::rc_ptr<::mosek::fusion::Expression> e2);
  static monty::rc_ptr<::mosek::fusion::Expression>
  sub(monty::rc_ptr<::mosek::fusion::Expression> e1,
      std::shared_ptr<monty::ndarray<double, 2>> a2);
  static monty::rc_ptr<::mosek::fusion::Expression>
  sub(monty::rc_ptr<::mosek::fusion::Expression> e1,
      std::shared_ptr<monty::ndarray<double, 1>> a1);
  static monty::rc_ptr<::mosek::fusion::Expression>
  sub(monty::rc_ptr<::mosek::fusion::Expression> e1,
      monty::rc_ptr<::mosek::fusion::Expression> e2);
  static monty::rc_ptr<::mosek::fusion::Expression>
  add(monty::rc_ptr<::mosek::fusion::NDSparseArray> n,
      monty::rc_ptr<::mosek::fusion::Expression> e2);
  static monty::rc_ptr<::mosek::fusion::Expression>
  add(monty::rc_ptr<::mosek::fusion::Expression> e1,
      monty::rc_ptr<::mosek::fusion::NDSparseArray> n);
  static monty::rc_ptr<::mosek::fusion::Expression>
  add(monty::rc_ptr<::mosek::fusion::Matrix> m,
      monty::rc_ptr<::mosek::fusion::Expression> e2);
  static monty::rc_ptr<::mosek::fusion::Expression>
  add(monty::rc_ptr<::mosek::fusion::Expression> e1,
      monty::rc_ptr<::mosek::fusion::Matrix> m);
  static monty::rc_ptr<::mosek::fusion::Expression>
  add(double c, monty::rc_ptr<::mosek::fusion::Expression> e2);
  static monty::rc_ptr<::mosek::fusion::Expression>
  add(monty::rc_ptr<::mosek::fusion::Expression> e1, double c);
  static monty::rc_ptr<::mosek::fusion::Expression>
  add(std::shared_ptr<monty::ndarray<double, 2>> a2,
      monty::rc_ptr<::mosek::fusion::Expression> e2);
  static monty::rc_ptr<::mosek::fusion::Expression>
  add(std::shared_ptr<monty::ndarray<double, 1>> a1,
      monty::rc_ptr<::mosek::fusion::Expression> e2);
  static monty::rc_ptr<::mosek::fusion::Expression>
  add(monty::rc_ptr<::mosek::fusion::Expression> e1,
      std::shared_ptr<monty::ndarray<double, 2>> a2);
  static monty::rc_ptr<::mosek::fusion::Expression>
  add(monty::rc_ptr<::mosek::fusion::Expression> e1,
      std::shared_ptr<monty::ndarray<double, 1>> a1);
  static monty::rc_ptr<::mosek::fusion::Expression>
  add(monty::rc_ptr<::mosek::fusion::Expression> e1,
      monty::rc_ptr<::mosek::fusion::Expression> e2);
  virtual /* override */ void
  eval(monty::rc_ptr<::mosek::fusion::WorkStack> rs,
       monty::rc_ptr<::mosek::fusion::WorkStack> ws,
       monty::rc_ptr<::mosek::fusion::WorkStack> xs);
}; // class Expr;

class WorkStack : public virtual monty::RefCounted {
public:
  p_WorkStack *_impl;

protected:
  WorkStack(p_WorkStack *_impl);

public:
  WorkStack(const WorkStack &) = delete;
  const WorkStack &operator=(const WorkStack &) = delete;
  friend class p_WorkStack;
  virtual ~WorkStack();
  virtual void destroy();
  typedef monty::rc_ptr<WorkStack> t;
  int32_t get_code_base();
  void set_code_base(int32_t val);
  int32_t get_cconst_base();
  void set_cconst_base(int32_t val);
  int32_t get_codeptr_base();
  void set_codeptr_base(int32_t val);
  int32_t get_cof_base();
  void set_cof_base(int32_t val);
  int32_t get_nidxs_base();
  void set_nidxs_base(int32_t val);
  int32_t get_sp_base();
  void set_sp_base(int32_t val);
  int32_t get_shape_base();
  void set_shape_base(int32_t val);
  int32_t get_ptr_base();
  void set_ptr_base(int32_t val);
  bool get_hassp();
  void set_hassp(bool val);
  int32_t get_ncodeatom();
  void set_ncodeatom(int32_t val);
  int32_t get_nelem();
  void set_nelem(int32_t val);
  int32_t get_nnz();
  void set_nnz(int32_t val);
  int32_t get_nd();
  void set_nd(int32_t val);
  int32_t get_pf64();
  void set_pf64(int32_t val);
  int32_t get_pi64();
  void set_pi64(int32_t val);
  int32_t get_pi32();
  void set_pi32(int32_t val);
  std::shared_ptr<monty::ndarray<double, 1>> get_f64();
  void set_f64(std::shared_ptr<monty::ndarray<double, 1>> val);
  std::shared_ptr<monty::ndarray<int64_t, 1>> get_i64();
  void set_i64(std::shared_ptr<monty::ndarray<int64_t, 1>> val);
  std::shared_ptr<monty::ndarray<int32_t, 1>> get_i32();
  void set_i32(std::shared_ptr<monty::ndarray<int32_t, 1>> val);

  WorkStack();
  virtual std::string formatCurrent();
  virtual bool peek_hassp();
  virtual int32_t peek_nnz();
  virtual int32_t peek_nelem();
  virtual int32_t peek_dim(int32_t i);
  virtual int32_t peek_nd();
  virtual void alloc_expr(int32_t nd, int32_t nelem, int32_t nnz, bool hassp);
  virtual void alloc_expr(int32_t nd, int32_t nelem, int32_t nnz, bool hassp,
                          int32_t ncodeatom);
  virtual void pop_expr();
  virtual void move_expr(monty::rc_ptr<::mosek::fusion::WorkStack> to);
  virtual void peek_expr();
  virtual void ensure_sparsity();
  virtual void clear();
  virtual int32_t allocf64(int32_t n);
  virtual int32_t alloci64(int32_t n);
  virtual int32_t alloci32(int32_t n);
  virtual void pushf64(double v);
  virtual void pushi64(int64_t v);
  virtual void pushi32(int32_t v);
  virtual void ensuref64(int32_t n);
  virtual void ensurei64(int32_t n);
  virtual void ensurei32(int32_t n);
  virtual int32_t popf64(int32_t n);
  virtual int32_t popi64(int32_t n);
  virtual int32_t popi32(int32_t n);
  virtual void popf64(int32_t n, std::shared_ptr<monty::ndarray<double, 1>> r,
                      int32_t ofs);
  virtual void popi64(int32_t n, std::shared_ptr<monty::ndarray<int64_t, 1>> r,
                      int32_t ofs);
  virtual void popi32(int32_t n, std::shared_ptr<monty::ndarray<int32_t, 1>> r,
                      int32_t ofs);
  virtual double popf64();
  virtual int64_t popi64();
  virtual int32_t popi32();
  virtual double peekf64();
  virtual int64_t peeki64();
  virtual int32_t peeki32();
  virtual double peekf64(int32_t i);
  virtual int64_t peeki64(int32_t i);
  virtual int32_t peeki32(int32_t i);
}; // class WorkStack;

class SymmetricMatrix : public virtual monty::RefCounted {
public:
  p_SymmetricMatrix *_impl;
  SymmetricMatrix(int32_t dim0, int32_t dim1,
                  std::shared_ptr<monty::ndarray<int32_t, 1>> usubi,
                  std::shared_ptr<monty::ndarray<int32_t, 1>> usubj,
                  std::shared_ptr<monty::ndarray<double, 1>> uval,
                  std::shared_ptr<monty::ndarray<int32_t, 1>> vsubi,
                  std::shared_ptr<monty::ndarray<int32_t, 1>> vsubj,
                  std::shared_ptr<monty::ndarray<double, 1>> vval,
                  double scale);

protected:
  SymmetricMatrix(p_SymmetricMatrix *_impl);

public:
  SymmetricMatrix(const SymmetricMatrix &) = delete;
  const SymmetricMatrix &operator=(const SymmetricMatrix &) = delete;
  friend class p_SymmetricMatrix;
  virtual ~SymmetricMatrix();
  virtual void destroy();
  typedef monty::rc_ptr<SymmetricMatrix> t;

  static monty::rc_ptr<::mosek::fusion::SymmetricMatrix>
  rankOne(int32_t n, std::shared_ptr<monty::ndarray<int32_t, 1>> sub,
          std::shared_ptr<monty::ndarray<double, 1>> v);
  static monty::rc_ptr<::mosek::fusion::SymmetricMatrix>
  rankOne(std::shared_ptr<monty::ndarray<double, 1>> v);
  static monty::rc_ptr<::mosek::fusion::SymmetricMatrix>
  antiDiag(std::shared_ptr<monty::ndarray<double, 1>> vals);
  static monty::rc_ptr<::mosek::fusion::SymmetricMatrix>
  diag(std::shared_ptr<monty::ndarray<double, 1>> vals);
  virtual monty::rc_ptr<::mosek::fusion::SymmetricMatrix>
  __mosek_2fusion_2SymmetricMatrix__add(
      monty::rc_ptr<::mosek::fusion::SymmetricMatrix> m);
  monty::rc_ptr<::mosek::fusion::SymmetricMatrix>
  add(monty::rc_ptr<::mosek::fusion::SymmetricMatrix> m);
  virtual monty::rc_ptr<::mosek::fusion::SymmetricMatrix>
  __mosek_2fusion_2SymmetricMatrix__sub(
      monty::rc_ptr<::mosek::fusion::SymmetricMatrix> m);
  monty::rc_ptr<::mosek::fusion::SymmetricMatrix>
  sub(monty::rc_ptr<::mosek::fusion::SymmetricMatrix> m);
  virtual monty::rc_ptr<::mosek::fusion::SymmetricMatrix>
  __mosek_2fusion_2SymmetricMatrix__mul(double v);
  monty::rc_ptr<::mosek::fusion::SymmetricMatrix> mul(double v);
  virtual int32_t getdim();
}; // class SymmetricMatrix;

class NDSparseArray : public virtual monty::RefCounted {
public:
  p_NDSparseArray *_impl;
  NDSparseArray(std::shared_ptr<monty::ndarray<int32_t, 1>> dims_,
                std::shared_ptr<monty::ndarray<int32_t, 2>> sub,
                std::shared_ptr<monty::ndarray<double, 1>> cof_);
  NDSparseArray(std::shared_ptr<monty::ndarray<int32_t, 1>> dims_,
                std::shared_ptr<monty::ndarray<int64_t, 1>> inst_,
                std::shared_ptr<monty::ndarray<double, 1>> cof_);
  NDSparseArray(monty::rc_ptr<::mosek::fusion::Matrix> m);

protected:
  NDSparseArray(p_NDSparseArray *_impl);

public:
  NDSparseArray(const NDSparseArray &) = delete;
  const NDSparseArray &operator=(const NDSparseArray &) = delete;
  friend class p_NDSparseArray;
  virtual ~NDSparseArray();
  virtual void destroy();
  typedef monty::rc_ptr<NDSparseArray> t;

  static monty::rc_ptr<::mosek::fusion::NDSparseArray>
  make(monty::rc_ptr<::mosek::fusion::Matrix> m);
  static monty::rc_ptr<::mosek::fusion::NDSparseArray>
  make(std::shared_ptr<monty::ndarray<int32_t, 1>> dims,
       std::shared_ptr<monty::ndarray<int64_t, 1>> inst,
       std::shared_ptr<monty::ndarray<double, 1>> cof);
  static monty::rc_ptr<::mosek::fusion::NDSparseArray>
  make(std::shared_ptr<monty::ndarray<int32_t, 1>> dims,
       std::shared_ptr<monty::ndarray<int32_t, 2>> sub,
       std::shared_ptr<monty::ndarray<double, 1>> cof);
}; // class NDSparseArray;

class Matrix : public virtual monty::RefCounted {
public:
  p_Matrix *_impl;

protected:
  Matrix(p_Matrix *_impl);

public:
  Matrix(const Matrix &) = delete;
  const Matrix &operator=(const Matrix &) = delete;
  friend class p_Matrix;
  virtual ~Matrix();
  virtual void destroy();
  typedef monty::rc_ptr<Matrix> t;

  virtual /* override */ std::string toString();
  static monty::rc_ptr<::mosek::fusion::Matrix>
  diag(int32_t num, monty::rc_ptr<::mosek::fusion::Matrix> mv);
  static monty::rc_ptr<::mosek::fusion::Matrix> diag(
      std::shared_ptr<monty::ndarray<monty::rc_ptr<::mosek::fusion::Matrix>, 1>>
          md);
  static monty::rc_ptr<::mosek::fusion::Matrix> antidiag(int32_t n, double val,
                                                         int32_t k);
  static monty::rc_ptr<::mosek::fusion::Matrix> antidiag(int32_t n, double val);
  static monty::rc_ptr<::mosek::fusion::Matrix> diag(int32_t n, double val,
                                                     int32_t k);
  static monty::rc_ptr<::mosek::fusion::Matrix> diag(int32_t n, double val);
  static monty::rc_ptr<::mosek::fusion::Matrix>
  antidiag(std::shared_ptr<monty::ndarray<double, 1>> d, int32_t k);
  static monty::rc_ptr<::mosek::fusion::Matrix>
  antidiag(std::shared_ptr<monty::ndarray<double, 1>> d);
  static monty::rc_ptr<::mosek::fusion::Matrix>
  diag(std::shared_ptr<monty::ndarray<double, 1>> d, int32_t k);
  static monty::rc_ptr<::mosek::fusion::Matrix>
  diag(std::shared_ptr<monty::ndarray<double, 1>> d);
  static monty::rc_ptr<::mosek::fusion::Matrix> ones(int32_t n, int32_t m);
  static monty::rc_ptr<::mosek::fusion::Matrix> eye(int32_t n);
  static monty::rc_ptr<::mosek::fusion::Matrix>
  dense(monty::rc_ptr<::mosek::fusion::Matrix> other);
  static monty::rc_ptr<::mosek::fusion::Matrix>
  dense(int32_t dimi, int32_t dimj, double value);
  static monty::rc_ptr<::mosek::fusion::Matrix>
  dense(int32_t dimi, int32_t dimj,
        std::shared_ptr<monty::ndarray<double, 1>> data);
  static monty::rc_ptr<::mosek::fusion::Matrix>
  dense(std::shared_ptr<monty::ndarray<double, 2>> data);
  static monty::rc_ptr<::mosek::fusion::Matrix>
  sparse(monty::rc_ptr<::mosek::fusion::Matrix> mx);
  static monty::rc_ptr<::mosek::fusion::Matrix>
  sparse(std::shared_ptr<
         monty::ndarray<std::shared_ptr<monty::ndarray<
                            monty::rc_ptr<::mosek::fusion::Matrix>, 1>>,
                        1>>
             blocks);
  static monty::rc_ptr<::mosek::fusion::Matrix>
  sparse(std::shared_ptr<monty::ndarray<double, 2>> data);
  static monty::rc_ptr<::mosek::fusion::Matrix> sparse(int32_t nrow,
                                                       int32_t ncol);
  static monty::rc_ptr<::mosek::fusion::Matrix>
  sparse(int32_t nrow, int32_t ncol,
         std::shared_ptr<monty::ndarray<int32_t, 1>> subi,
         std::shared_ptr<monty::ndarray<int32_t, 1>> subj, double val);
  static monty::rc_ptr<::mosek::fusion::Matrix>
  sparse(std::shared_ptr<monty::ndarray<int32_t, 1>> subi,
         std::shared_ptr<monty::ndarray<int32_t, 1>> subj, double val);
  static monty::rc_ptr<::mosek::fusion::Matrix>
  sparse(std::shared_ptr<monty::ndarray<int32_t, 1>> subi,
         std::shared_ptr<monty::ndarray<int32_t, 1>> subj,
         std::shared_ptr<monty::ndarray<double, 1>> val);
  static monty::rc_ptr<::mosek::fusion::Matrix>
  sparse(int32_t nrow, int32_t ncol,
         std::shared_ptr<monty::ndarray<int32_t, 1>> subi,
         std::shared_ptr<monty::ndarray<int32_t, 1>> subj,
         std::shared_ptr<monty::ndarray<double, 1>> val);
  virtual double get(int32_t i, int32_t j) = 0;
  virtual bool isSparse() = 0;
  virtual std::shared_ptr<monty::ndarray<double, 1>> getDataAsArray() = 0;
  virtual void
  getDataAsTriplets(std::shared_ptr<monty::ndarray<int32_t, 1>> subi,
                    std::shared_ptr<monty::ndarray<int32_t, 1>> subj,
                    std::shared_ptr<monty::ndarray<double, 1>> val) = 0;
  virtual monty::rc_ptr<::mosek::fusion::Matrix>
  __mosek_2fusion_2Matrix__transpose() = 0;
  monty::rc_ptr<::mosek::fusion::Matrix> transpose();
  virtual int64_t numNonzeros() = 0;
  virtual int32_t numColumns();
  virtual int32_t numRows();
}; // class Matrix;

class DenseMatrix : public ::mosek::fusion::Matrix {
  DenseMatrix(int32_t dimi_, int32_t dimj_,
              std::shared_ptr<monty::ndarray<double, 1>> cof);
  DenseMatrix(monty::rc_ptr<::mosek::fusion::Matrix> m_);
  DenseMatrix(std::shared_ptr<monty::ndarray<double, 2>> d);
  DenseMatrix(int32_t dimi_, int32_t dimj_, double value_);

protected:
  DenseMatrix(p_DenseMatrix *_impl);

public:
  DenseMatrix(const DenseMatrix &) = delete;
  const DenseMatrix &operator=(const DenseMatrix &) = delete;
  friend class p_DenseMatrix;
  virtual ~DenseMatrix();
  virtual void destroy();
  typedef monty::rc_ptr<DenseMatrix> t;

  virtual /* override */ std::string toString();
  virtual /* override */ monty::rc_ptr<::mosek::fusion::Matrix>
  __mosek_2fusion_2DenseMatrix__transpose();
  monty::rc_ptr<::mosek::fusion::Matrix> transpose();
  /* override: mosek.fusion.Matrix.transpose*/
  virtual monty::rc_ptr<::mosek::fusion::Matrix>
  __mosek_2fusion_2Matrix__transpose();
  virtual /* override */ bool isSparse();
  virtual /* override */ std::shared_ptr<monty::ndarray<double, 1>>
  getDataAsArray();
  virtual /* override */ void
  getDataAsTriplets(std::shared_ptr<monty::ndarray<int32_t, 1>> subi,
                    std::shared_ptr<monty::ndarray<int32_t, 1>> subj,
                    std::shared_ptr<monty::ndarray<double, 1>> cof);
  virtual /* override */ double get(int32_t i, int32_t j);
  virtual /* override */ int64_t numNonzeros();
}; // class DenseMatrix;

class SparseMatrix : public ::mosek::fusion::Matrix {
  SparseMatrix(int32_t dimi_, int32_t dimj_,
               std::shared_ptr<monty::ndarray<int32_t, 1>> subi_,
               std::shared_ptr<monty::ndarray<int32_t, 1>> subj_,
               std::shared_ptr<monty::ndarray<double, 1>> val_, int64_t nelm);
  SparseMatrix(int32_t dimi_, int32_t dimj_,
               std::shared_ptr<monty::ndarray<int32_t, 1>> subi_,
               std::shared_ptr<monty::ndarray<int32_t, 1>> subj_,
               std::shared_ptr<monty::ndarray<double, 1>> val_);

protected:
  SparseMatrix(p_SparseMatrix *_impl);

public:
  SparseMatrix(const SparseMatrix &) = delete;
  const SparseMatrix &operator=(const SparseMatrix &) = delete;
  friend class p_SparseMatrix;
  virtual ~SparseMatrix();
  virtual void destroy();
  typedef monty::rc_ptr<SparseMatrix> t;

  virtual /* override */ std::string toString();
  virtual /* override */ int64_t numNonzeros();
  virtual /* override */ monty::rc_ptr<::mosek::fusion::Matrix>
  __mosek_2fusion_2SparseMatrix__transpose();
  monty::rc_ptr<::mosek::fusion::Matrix> transpose();
  /* override: mosek.fusion.Matrix.transpose*/
  virtual monty::rc_ptr<::mosek::fusion::Matrix>
  __mosek_2fusion_2Matrix__transpose();
  virtual /* override */ bool isSparse();
  virtual /* override */ std::shared_ptr<monty::ndarray<double, 1>>
  getDataAsArray();
  virtual /* override */ void
  getDataAsTriplets(std::shared_ptr<monty::ndarray<int32_t, 1>> subi_,
                    std::shared_ptr<monty::ndarray<int32_t, 1>> subj_,
                    std::shared_ptr<monty::ndarray<double, 1>> cof_);
  virtual /* override */ double get(int32_t i, int32_t j);
}; // class SparseMatrix;

class LinkedBlocks : public virtual monty::RefCounted {
public:
  p_LinkedBlocks *_impl;

protected:
  LinkedBlocks(p_LinkedBlocks *_impl);

public:
  LinkedBlocks(const LinkedBlocks &) = delete;
  const LinkedBlocks &operator=(const LinkedBlocks &) = delete;
  friend class p_LinkedBlocks;
  virtual ~LinkedBlocks();
  virtual void destroy();
  typedef monty::rc_ptr<LinkedBlocks> t;

  LinkedBlocks();
  LinkedBlocks(int32_t n);
  LinkedBlocks(monty::rc_ptr<::mosek::fusion::LinkedBlocks> other);
  virtual void free(int32_t bkey);
  virtual int32_t alloc(int32_t size);
  virtual int32_t maxidx(int32_t bkey);
  virtual int32_t numallocated();
  virtual void get(int32_t bkey,
                   std::shared_ptr<monty::ndarray<int32_t, 1>> target,
                   int32_t offset);
  virtual int32_t numblocks();
  virtual int32_t blocksize(int32_t bkey);
  virtual int32_t block_capacity();
  virtual int32_t capacity();
  virtual bool validate();
}; // class LinkedBlocks;

class LinkedInts : public virtual monty::RefCounted {
public:
  p_LinkedInts *_impl;

protected:
  LinkedInts(p_LinkedInts *_impl);

public:
  LinkedInts(const LinkedInts &) = delete;
  const LinkedInts &operator=(const LinkedInts &) = delete;
  friend class p_LinkedInts;
  virtual ~LinkedInts();
  virtual void destroy();
  typedef monty::rc_ptr<LinkedInts> t;

  LinkedInts(int32_t cap_);
  LinkedInts();
  LinkedInts(monty::rc_ptr<::mosek::fusion::LinkedInts> other);
  virtual void free(int32_t i, int32_t num);
  virtual int32_t alloc();
  virtual int32_t alloc(int32_t n);
  virtual void alloc(int32_t num,
                     std::shared_ptr<monty::ndarray<int32_t, 1>> target,
                     int32_t offset);
  virtual void get(int32_t i, int32_t num,
                   std::shared_ptr<monty::ndarray<int32_t, 1>> target,
                   int32_t offset);
  virtual int32_t numallocated();
  virtual int32_t maxidx(int32_t i, int32_t num);
  virtual int32_t capacity();
  virtual bool validate();
}; // class LinkedInts;

class Parameters : public virtual monty::RefCounted {
public:
  p_Parameters *_impl;

protected:
  Parameters(p_Parameters *_impl);

public:
  Parameters(const Parameters &) = delete;
  const Parameters &operator=(const Parameters &) = delete;
  friend class p_Parameters;
  virtual ~Parameters();
  virtual void destroy();
  typedef monty::rc_ptr<Parameters> t;

  static void setParameter(monty::rc_ptr<::mosek::fusion::Model> M,
                           const std::string &name, double value);
  static void setParameter(monty::rc_ptr<::mosek::fusion::Model> M,
                           const std::string &name, int32_t value);
  static void setParameter(monty::rc_ptr<::mosek::fusion::Model> M,
                           const std::string &name, const std::string &value);
}; // class Parameters;

} // namespace fusion
} // namespace mosek
namespace mosek {
namespace fusion {
namespace Utils {
// class mosek.fusion.Utils.StringIntMap
// mosek.fusion.Utils.IntMap from file 'src/fusion/cxx/IntMap.h'
// namespace mosek::fusion::Utils
class IntMap : public monty::RefCounted {
  std::unique_ptr<p_IntMap> _impl;

public:
  friend class p_IntMap;
  typedef monty::rc_ptr<IntMap> t;

  IntMap();
  bool hasItem(int64_t key);
  int getItem(int64_t key);
  void setItem(int64_t key, int val);
  std::shared_ptr<monty::ndarray<int64_t, 1>> keys();
  std::shared_ptr<monty::ndarray<int, 1>> values();

  t clone();
  t __mosek_2fusion_2Utils_2IntMap__clone();
};

class StringIntMap : public monty::RefCounted {
  std::unique_ptr<p_StringIntMap> _impl;

public:
  friend class p_StringIntMap;
  typedef monty::rc_ptr<StringIntMap> t;

  StringIntMap();
  bool hasItem(const std::string &key);
  int getItem(const std::string &key);
  void setItem(const std::string &key, int val);
  std::shared_ptr<monty::ndarray<std::string, 1>> keys();
  std::shared_ptr<monty::ndarray<int, 1>> values();

  t clone();
  t __mosek_2fusion_2Utils_2StringIntMap__clone();
};
// End of file 'src/fusion/cxx/IntMap.h'
// mosek.fusion.Utils.StringBuffer from file 'src/fusion/cxx/StringBuffer.h'
// namespace mosek::fusion::Utils
class StringBuffer : public monty::RefCounted {
private:
  std::unique_ptr<p_StringBuffer> _impl;

public:
  friend class p_StringBuffer;

  typedef monty::rc_ptr<StringBuffer> t;

  StringBuffer();
  t clear();
  t a(int value);
  t a(int64_t value);
  t a(double value);
  t a(const std::string &value);
  t a(bool value);
  t a(std::shared_ptr<monty::ndarray<std::string, 1>> value);
  t a(std::shared_ptr<monty::ndarray<int, 1>> value);
  t a(std::shared_ptr<monty::ndarray<int64_t, 1>> value);
  t a(std::shared_ptr<monty::ndarray<double, 1>> value);
  t lf();

  t __mosek_2fusion_2Utils_2StringBuffer__clear();
  t __mosek_2fusion_2Utils_2StringBuffer__a(int value);
  t __mosek_2fusion_2Utils_2StringBuffer__a(int64_t value);
  t __mosek_2fusion_2Utils_2StringBuffer__a(double value);
  t __mosek_2fusion_2Utils_2StringBuffer__a(const std::string &value);
  t __mosek_2fusion_2Utils_2StringBuffer__a(bool value);
  t __mosek_2fusion_2Utils_2StringBuffer__a(
      std::shared_ptr<monty::ndarray<std::string, 1>> value);
  t __mosek_2fusion_2Utils_2StringBuffer__a(
      std::shared_ptr<monty::ndarray<int, 1>> value);
  t __mosek_2fusion_2Utils_2StringBuffer__a(
      std::shared_ptr<monty::ndarray<int64_t, 1>> value);
  t __mosek_2fusion_2Utils_2StringBuffer__a(
      std::shared_ptr<monty::ndarray<double, 1>> value);
  t __mosek_2fusion_2Utils_2StringBuffer__lf();

  std::string toString() const;
};
// End of file 'src/fusion/cxx/StringBuffer.h'
// mosek.fusion.Utils.Tools from file 'src/fusion/cxx/Tools.h'
namespace Tools {
template <typename T, int N>
void arraycopy(const std::shared_ptr<monty::ndarray<T, N>> &src, int srcoffset,
               const std::shared_ptr<monty::ndarray<T, N>> &tgt, int tgtoffset,
               int size) {
  std::copy(src->flat_begin() + srcoffset, src->flat_begin() + srcoffset + size,
            tgt->flat_begin() + tgtoffset);
}

template <typename T, int N>
void arraycopy(const std::shared_ptr<monty::ndarray<T, N>> &src,
               int64_t srcoffset,
               const std::shared_ptr<monty::ndarray<T, N>> &tgt,
               int64_t tgtoffset, int64_t size) {
  std::copy(src->flat_begin() + srcoffset, src->flat_begin() + srcoffset + size,
            tgt->flat_begin() + tgtoffset);
}

template <typename T, int N>
std::shared_ptr<monty::ndarray<T, N>>
arraycopy(const std::shared_ptr<monty::ndarray<T, N>> &a) {
  return std::shared_ptr<monty::ndarray<T, N>>(
      new monty::ndarray<T, N>(a->shape, a->flat_begin(), a->flat_end()));
}

template <typename T> std::shared_ptr<monty::ndarray<T, 1>> range(T last) {
  return std::shared_ptr<monty::ndarray<T, 1>>(new monty::ndarray<T, 1>(
      monty::shape(last), monty::iterable(monty::range_t<T>(0, last))));
}

template <typename T>
std::shared_ptr<monty::ndarray<T, 1>> range(T first, T last) {
  return std::shared_ptr<monty::ndarray<T, 1>>(new monty::ndarray<T, 1>(
      monty::shape(last - first),
      monty::iterable(monty::range_t<T>(first, last))));
}

template <typename T>
std::shared_ptr<monty::ndarray<T, 1>> range(T first, T last, T step) {
  size_t num = last > first && step > 0 ? (last - first - 1) / step + 1 : 0;
  if (num > 0)
    return std::shared_ptr<monty::ndarray<T, 1>>(new monty::ndarray<T, 1>(
        monty::shape(num),
        monty::iterable(monty::range_t<T>(first, last, step))));
  else
    return std::shared_ptr<monty::ndarray<T, 1>>(
        new monty::ndarray<T, 1>(monty::shape(0)));
}

static std::shared_ptr<monty::ndarray<double, 1>> zeros(int num) {
  return std::shared_ptr<monty::ndarray<double, 1>>(
      new monty::ndarray<double, 1>(monty::shape(num), 0.0));
}
static std::shared_ptr<monty::ndarray<double, 2>> zeros(int dimi, int dimj) {
  return std::shared_ptr<monty::ndarray<double, 2>>(
      new monty::ndarray<double, 2>(monty::shape(dimi, dimj), 0.0));
}
static std::shared_ptr<monty::ndarray<double, 1>> ones(int num) {
  return std::shared_ptr<monty::ndarray<double, 1>>(
      new monty::ndarray<double, 1>(monty::shape(num), 1.0));
}

template <typename T>
std::shared_ptr<monty::ndarray<T, 1>> makevector(T v, int num) {
  return std::shared_ptr<monty::ndarray<T, 1>>(
      new monty::ndarray<T, 1>(monty::shape(num), v));
}

template <typename T>
std::shared_ptr<monty::ndarray<T, 1>> repeatrange(T first, T last, T num) {
  return std::shared_ptr<monty::ndarray<T, 1>>(new monty::ndarray<T, 1>(
      monty::shape((last - first) * num),
      [=](ptrdiff_t i) { return (T)(i % num + first); }));
}

template <typename T> std::string stringvalue(T val) {
  std::stringstream os;
  os << val;
  return os.str();
}

static int toInt(const std::string &v) { return atoi(v.c_str()); }
static double toDouble(const std::string &v) { return atof(v.c_str()); }
static double sqrt(double v) { return std::sqrt(v); }

template <typename T>
void sort(const std::shared_ptr<monty::ndarray<T, 1>> &vals, int first,
          int last) {
  std::sort(vals->flat_begin() + first, vals->flat_begin() + last,
            [&](T lhs, T rhs) { return lhs < rhs; });
}

#if (defined(_WIN32) && _WIN32) || (defined(_WIN64) && _WIN64)
static int randInt(int max) {
  int64_t lo = rand(), hi = rand();
  return (int)(((double)(((hi << 31) | lo)) / (double)(LLONG_MAX)) * max);
}
#else
static int randInt(int max) {
  return (int)(((double)(random()) / (double)(RAND_MAX)) * max);
}
#endif

template <typename T>
static void argsort(const std::shared_ptr<monty::ndarray<int64_t, 1>> &perm,
                    const std::shared_ptr<monty::ndarray<T, 1>> &v,
                    int64_t first, int64_t last) {
  std::sort(perm->begin() + first, perm->begin() + last,
            [&v](int64_t lhs, int64_t rhs) { return (*v)[lhs] < (*v)[rhs]; });
}

template <typename T>
static void argsort(const std::shared_ptr<monty::ndarray<int64_t, 1>> &perm,
                    const std::shared_ptr<monty::ndarray<T, 1>> &v0,
                    const std::shared_ptr<monty::ndarray<T, 1>> &v1,
                    int64_t first, int64_t last) {
  std::sort(perm->begin() + first, perm->begin() + last,
            [&v0, &v1](int64_t lhs, int64_t rhs) {
              return (*v0)[lhs] < (*v0)[rhs] ||
                     ((*v0)[lhs] == (*v0)[rhs] && (*v1)[lhs] < (*v1)[rhs]);
            });
}

template <typename T>
static void bucketsort(const std::shared_ptr<monty::ndarray<int64_t, 1>> &perm,
                       int64_t first, int64_t last,
                       const std::shared_ptr<monty::ndarray<T, 1>> &v, T minval,
                       T maxval) {
  T N = maxval - minval + 1;
  int64_t M = last - first;
  std::vector<ptrdiff_t> ptrb(N + 1);
  std::vector<int64_t> nperm(M);
  for (ptrdiff_t i = first; i < last; ++i)
    ++ptrb[(*v)[(*perm)[i]] - minval + 1];
  for (ptrdiff_t i = 1; i < N; ++i)
    ptrb[i] += ptrb[i - 1];
  for (ptrdiff_t i = first; i < last; ++i)
    nperm[ptrb[(*v)[(*perm)[i]] - minval]++] = (*perm)[i];

  std::copy(nperm.begin(), nperm.end(), perm->begin() + first);
}
} // namespace Tools
// End of file 'src/fusion/cxx/Tools.h'
} // namespace Utils
} // namespace fusion
} // namespace mosek
#endif
