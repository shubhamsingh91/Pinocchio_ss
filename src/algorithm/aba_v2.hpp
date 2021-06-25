//
// Copyright (c) 2016-2018 CNRS
//

#ifndef __pinocchio_aba_v2_hpp__
#define __pinocchio_aba_v2_hpp__

#include "pinocchio/multibody/model.hpp"
#include "pinocchio/multibody/data.hpp"
#include "pinocchio/algorithm/check.hpp"

namespace pinocchio
{
 
  /// \param[in] tau_mat The joint torque vector (dim model.nv,2*model.nv).

  template<typename Scalar, int Options, template<typename,int> class JointCollectionTpl, typename ConfigVectorType, typename MatrixType1>
  inline const typename DataTpl<Scalar,Options,JointCollectionTpl>::RowMatrixXs &
  computeMinverse_v2(const ModelTpl<Scalar,Options,JointCollectionTpl> & model,
                  DataTpl<Scalar,Options,JointCollectionTpl> & data,
                  const Eigen::MatrixBase<ConfigVectorType> & q,
                  const Eigen::MatrixBase<MatrixType1> & tau_mat);


 // PINOCCHIO_DEFINE_ALGO_CHECKER(ABA);

} // namespace pinocchio

/* --- Details -------------------------------------------------------------------- */
#include "pinocchio/algorithm/aba_v2.hxx"

#endif // ifndef __pinocchio_aba_v2_hpp__
