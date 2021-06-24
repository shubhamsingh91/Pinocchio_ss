//
// Copyright (c) 2016-2018 CNRS
//

#ifndef __pinocchio_azamat_v4_hpp__
#define __pinocchio_azamat_v4_hpp__

#include "pinocchio/multibody/model.hpp"
#include "pinocchio/multibody/data.hpp"
#include "pinocchio/algorithm/check.hpp"

namespace pinocchio
{
  ///
  /// \brief The Articulated-Body algorithm. It computes the forward dynamics, aka the joint accelerations given the current state and actuation of the model.
  ///
  /// \tparam JointCollection Collection of Joint types.
  /// \tparam ConfigVectorType Type of the joint configuration vector.
  /// \tparam TangentVectorType1 Type of the joint velocity vector.
  /// \tparam TangentVectorType2 Type of the joint torque vector.
  ///
  /// \param[in] model The model structure of the rigid body system.
  /// \param[in] data The data structure of the rigid body system.
  /// \param[in] q The joint configuration vector (dim model.nq).
  /// \param[in] tau_mat The joint torque vector (dim model.nv,2*model.nv).
  ///
  /// \note This also overwrites data.f, possibly leaving it in an inconsistent state
  ///
  /// \return The current joint acceleration stored in data.Minv_mat_prod.
  ///
  template<typename Scalar, int Options, template<typename,int> class JointCollectionTpl, typename ConfigVectorType, typename MatrixType1>
  inline const typename DataTpl<Scalar,Options,JointCollectionTpl>::RowMatrixXs &
  azamat_v4(const ModelTpl<Scalar,Options,JointCollectionTpl> & model,
                  DataTpl<Scalar,Options,JointCollectionTpl> & data,
                  const Eigen::MatrixBase<ConfigVectorType> & q,
                  const Eigen::MatrixBase<MatrixType1> & tau_mat);

 // PINOCCHIO_DEFINE_ALGO_CHECKER(AZA);

} // namespace pinocchio

/* --- Details -------------------------------------------------------------------- */
#include "pinocchio/algorithm/azamat_v4.hxx"

#endif // ifndef __pinocchio_azamat_v4_hpp__
