//
// Copyright (c) 2016-2020 CNRS, INRIA
// THIS VERSION IS NOT CORRECT HERE, NOT TO BE USED

#ifndef __pinocchio_azamat_v2_hxx__
#define __pinocchio_azamat_v2_hxx__

#include "pinocchio/spatial/act-on-set.hpp"
#include "pinocchio/multibody/visitor.hpp"
#include "pinocchio/algorithm/check.hpp"
#include <iostream>
#include <stdlib.h>

/// @cond DEV

namespace pinocchio
{
  template<typename Scalar, int Options, template<typename,int> class JointCollectionTpl, typename ConfigVectorType, typename TangentVectorType>
  struct Azamatv2ForwardStep1
  : public fusion::JointUnaryVisitorBase< Azamatv2ForwardStep1<Scalar,Options,JointCollectionTpl,ConfigVectorType,TangentVectorType> >
  {
    typedef ModelTpl<Scalar,Options,JointCollectionTpl> Model;
    typedef DataTpl<Scalar,Options,JointCollectionTpl> Data;
    
    typedef boost::fusion::vector<const Model &,
                                  Data &,
                                  const ConfigVectorType &,
                                  const TangentVectorType &
                                  > ArgsType;
    
    template<typename JointModel>
    static void algo(const pinocchio::JointModelBase<JointModel> & jmodel,
                     pinocchio::JointDataBase<typename JointModel::JointDataDerived> & jdata,
                     const Model & model,
                     Data & data,
                     const Eigen::MatrixBase<ConfigVectorType> & q,
                     const Eigen::MatrixBase<TangentVectorType> & v)
    {
      typedef typename Model::JointIndex JointIndex;
      
      const JointIndex & i = jmodel.id();
      data.liMi[i] = model.jointPlacements[i] * jdata.M(); // Calculating transformation matrices here

        data.a_gf_v2[i].setZero();
        data.Yaba_v1[i] = model.inertias[i].matrix(); // Ii^A quantity 
        data.f_v2[i].setZero(); // pi^A, nx2n

    }
    
  };


  namespace internal
  {
    
    template<typename Scalar>
    struct SE3actOnazamat_v2
    {
      template<int Options, typename Matrix6Type>
      static typename PINOCCHIO_EIGEN_PLAIN_TYPE(Matrix6Type)
      run(const SE3Tpl<Scalar,Options> & M,
          const Eigen::MatrixBase<Matrix6Type> & I)
      {
        typedef SE3Tpl<Scalar,Options> SE3;
        typedef typename SE3::Matrix3 Matrix3;
        typedef typename SE3::Vector3 Vector3;
        
        typedef const Eigen::Block<Matrix6Type,3,3> constBlock3;
        
        typedef typename PINOCCHIO_EIGEN_PLAIN_TYPE(Matrix6Type) ReturnType;
        typedef Eigen::Block<ReturnType,3,3> Block3;
        
        Matrix6Type & I_ = PINOCCHIO_EIGEN_CONST_CAST(Matrix6Type,I);
        const constBlock3 & Ai = I_.template block<3,3>(Inertia::LINEAR, Inertia::LINEAR);
        const constBlock3 & Bi = I_.template block<3,3>(Inertia::LINEAR, Inertia::ANGULAR);
        const constBlock3 & Di = I_.template block<3,3>(Inertia::ANGULAR, Inertia::ANGULAR);
        
        const Matrix3 & R = M.rotation();
        const Vector3 & t = M.translation();
        
        ReturnType res;
        Block3 Ao = res.template block<3,3>(Inertia::LINEAR, Inertia::LINEAR);
        Block3 Bo = res.template block<3,3>(Inertia::LINEAR, Inertia::ANGULAR);
        Block3 Co = res.template block<3,3>(Inertia::ANGULAR, Inertia::LINEAR);
        Block3 Do = res.template block<3,3>(Inertia::ANGULAR, Inertia::ANGULAR);
        
        Do.noalias() = R*Ai; // tmp variable
        Ao.noalias() = Do*R.transpose();
        
        Do.noalias() = R*Bi; // tmp variable
        Bo.noalias() = Do*R.transpose();
        
        Co.noalias() = R*Di; // tmp variable
        Do.noalias() = Co*R.transpose();
        
        Do.row(0) += t.cross(Bo.col(0));
        Do.row(1) += t.cross(Bo.col(1));
        Do.row(2) += t.cross(Bo.col(2));
        
        Co.col(0) = t.cross(Ao.col(0));
        Co.col(1) = t.cross(Ao.col(1));
        Co.col(2) = t.cross(Ao.col(2));
        Co += Bo.transpose();
        
        Bo = Co.transpose();
        Do.col(0) += t.cross(Bo.col(0));
        Do.col(1) += t.cross(Bo.col(1));
        Do.col(2) += t.cross(Bo.col(2));
        
        return res;
      }
    };
    
  }

    
  template<typename Scalar, int Options, template<typename,int> class JointCollectionTpl>
  struct Azamatv2BackwardStep
  : public fusion::JointUnaryVisitorBase< Azamatv2BackwardStep<Scalar,Options,JointCollectionTpl> >
  {
    typedef ModelTpl<Scalar,Options,JointCollectionTpl> Model;
    typedef DataTpl<Scalar,Options,JointCollectionTpl> Data;
    
    typedef boost::fusion::vector<const Model &,
                                  Data &> ArgsType;
    
    template<typename JointModel>
    static void algo(const JointModelBase<JointModel> & jmodel,
                     JointDataBase<typename JointModel::JointDataDerived> & jdata,
                     const Model & model,
                     Data & data)
    {
      typedef typename Model::JointIndex JointIndex;
      typedef typename Data::Inertia Inertia;
      typedef typename Data::Force Force;
      
      const JointIndex & i = jmodel.id();
      const JointIndex & parent  = model.parents[i];
      typename Inertia::Matrix6 & Ia = data.Yaba_v1[i];

       jmodel.jointVelocitySelector_mat(data.u_v2) -= jdata.S().transpose()*data.f_v2[i];
     
      //   std:: cout << "--------------------------------------------------"  << std::endl;
      //   std:: cout << "size of data.u_v2 is" << data.u_v2.rows() <<"x" << data.u_v2.cols() << std::endl;
      //   std:: cout << "size of data.f_v2[i] is" << data.f_v2[i].rows() << "x" << data.f_v2[i].cols() << std::endl;
      //   std:: cout << "i here is" <<  i << std::endl;
      // // std:: cout << "data.u here is for i" <<  jmodel.jointVelocitySelector_mat(data.u_v2) << std::endl;
      //   std:: cout << "data.u is" <<  data.u_v2 << std:: endl;
      //   std:: cout << "--------------------------------------------------"  << std::endl;

       jmodel.calc_aba(jdata.derived(), Ia, parent > 0);

       if (parent > 0)
       {

         // first line
         data.pa_v2 = data.f_v2[i];
         
         // second line
         data.pa_v2 +=  jdata.UDinv() * jmodel.jointVelocitySelector_mat(data.u_v2);

       // std::cout << "jdata.UDinv() is" << jdata.UDinv()<< std::endl;
       // std::cout << "data.u_v2 is" << data.u_v2 << std::endl;
       // std::cout << "pa_v2 is" << data.pa_v2 << std::endl;

         // third line

         data.Yaba_v1[parent] += internal::SE3actOnazamat_v2<Scalar>::run(data.liMi[i], Ia);

          // fourth line

         forceSet::se3Action(data.liMi[i],data.pa_v2 ,data.pa_v3);
        
         data.f_v2[parent] += data.pa_v3;

          // std::cout << "data.limi[i] here is" << data.liMi[i] << std::endl;
         // std::cout << "pa_v3 is" << data.pa_v3 << std::endl;


       }

    }
    
  };
  
  template<typename Scalar, int Options, template<typename,int> class JointCollectionTpl>
  struct Azamatv2ForwardStep2
  : public fusion::JointUnaryVisitorBase< Azamatv2ForwardStep2<Scalar,Options,JointCollectionTpl> >
  {
    typedef ModelTpl<Scalar,Options,JointCollectionTpl> Model;
    typedef DataTpl<Scalar,Options,JointCollectionTpl> Data;
    
    typedef boost::fusion::vector<const Model &,
                                  Data &> ArgsType;
    
    template<typename JointModel>
    static void algo(const pinocchio::JointModelBase<JointModel> & jmodel,
                     pinocchio::JointDataBase<typename JointModel::JointDataDerived> & jdata,
                     const Model & model,
                     Data & data)
    {
        typedef typename Model::JointIndex JointIndex;
        
        const JointIndex & i = jmodel.id();
        const JointIndex & parent = model.parents[i];
        
        //--- line-1

        //  std::cout <<"Testing loop 3 here ---------------------" << std:: endl;
         motionSet::se3ActionInverse(data.liMi[i],data.a_gf_v2[parent] , data.a_gf_v2[i]);

        //  data.a_gf_v2[i] += data.liMi[i].actInv(data.a_gf_v2[parent]);

        // std::cout <<"data.a_gf[i] here is " << data.a_gf_v2[i] << std:: endl;
        // std::cout <<"data.a_gf[parent] here is " << data.a_gf_v2[parent] << std:: endl;

        //--- line-2

        // jmodel.jointVelocitySelector(data.ddq_new_v2).noalias() =
        // jdata.Dinv() * jmodel.jointVelocitySelector(data.u) - jdata.UDinv().transpose() * data.a_gf_v2[i].toVector();
      
         jmodel.jointVelocitySelector_mat(data.Minv_mat_prod_v2).noalias() =
         jdata.Dinv() * jmodel.jointVelocitySelector_mat(data.u_v2) - jdata.UDinv().transpose() * data.a_gf_v2[i];
       
       // std::cout << "jdata.Dinv() here is" << jdata.Dinv() << std::endl;
       // std::cout << "Minv matrix here is" << data.Minv_mat_prod_v2 << std::endl;

        //--- line-3

       // std::cout << "a_gf_v2[i] here is" << data.a_gf_v2[i] << std::endl;
        
       // Eigen::MatrixXd temp1(6,1);
        // temp1.row(0) = 0.1;
        // temp1.row(1) = 0.1;
        // temp1.row(2) = 0.1;
        // temp1.row(3) = 0.1;    
        // temp1.row(4) = 0.1;
        // temp1.row(5) = 0.1;

       // std::cout << "temp1 vector is" << temp1 << std::endl;

       // data.a_gf_v2[i] = temp1*jmodel.jointVelocitySelector_mat(data.Minv_mat_prod_v2);

       
        data.a_gf_v2[i] = jdata.S().matrix()*jmodel.jointVelocitySelector_mat(data.Minv_mat_prod_v2);
       
      //  std::cout << "data.a_gf_v2[i] is" << data.a_gf_v2[i] << std::endl;
     
      //  std::cout << "jdata.S here is " << jdata.S().matrix() << std:: endl;

    }
    
  };
  
  template<typename Scalar, int Options, template<typename,int> class JointCollectionTpl, typename ConfigVectorType, typename TangentVectorType1, typename MatrixType1>
  inline const typename DataTpl<Scalar,Options,JointCollectionTpl>::TangentVectorType &
  azamat_v2(const ModelTpl<Scalar,Options,JointCollectionTpl> & model,
      DataTpl<Scalar,Options,JointCollectionTpl> & data,
      const Eigen::MatrixBase<ConfigVectorType> & q,
      const Eigen::MatrixBase<TangentVectorType1> & v,
      const Eigen::MatrixBase<MatrixType1> & tau_mat)
  {
    assert(model.check(data) && "data is not consistent with model.");
    PINOCCHIO_CHECK_ARGUMENT_SIZE(q.size(), model.nq, "The joint configuration vector is not of right size");
    PINOCCHIO_CHECK_ARGUMENT_SIZE(v.size(), model.nv, "The joint velocity vector is not of right size");
    //PINOCCHIO_CHECK_ARGUMENT_SIZE(tau_mat.cols(),model.nv);
    //PINOCCHIO_CHECK_ARGUMENT_SIZE(tau_mat.rows(),model.nv);

    typedef typename ModelTpl<Scalar,Options,JointCollectionTpl>::JointIndex JointIndex;
    
    data.a_gf_v2[0].setZero();
    data.u_v2 = tau_mat;

      // std::cout <<"u matrix  is" << data.u_v2  << std::endl;

    typedef Azamatv2ForwardStep1<Scalar,Options,JointCollectionTpl,ConfigVectorType,TangentVectorType1> Pass1;
    
    typedef Azamatv2BackwardStep<Scalar,Options,JointCollectionTpl> Pass2;

    typedef Azamatv2ForwardStep2<Scalar,Options,JointCollectionTpl> Pass3;

     //------------- Pass 1--------------------------------//

      for(JointIndex i=1; i<(JointIndex)model.njoints; ++i)
      {
        Pass1::run(model.joints[i],data.joints[i],
                  typename Pass1::ArgsType(model,data,q.derived(),v.derived()));
      }   
     
       //------------- Pass 2--------------------------------//

      //  std::cout << "data.u is" << data.u << std::endl;

        for(JointIndex i=(JointIndex)model.njoints-1;i>0; --i)
        {

       //  std::cout <<"i index in the second loop is" << i << std::endl;

        Pass2::run(model.joints[i],data.joints[i],
                    typename Pass2::ArgsType(model,data));
        }

      //  //------------- Pass 3--------------------------------//

        for(JointIndex i=1; i<(JointIndex)model.njoints; ++i)
        {
       //  std::cout <<"i index in the third loop is" << i << std::endl;

        Pass3::run(model.joints[i],data.joints[i],
                    typename Pass3::ArgsType(model,data));

        }

    return data.Minv_mat_prod_v2;
  }

  template<typename Scalar, int Options, template<typename,int> class JointCollectionTpl, typename ConfigVectorType, typename TangentVectorType1, typename MatrixType1, typename ForceDerived>
  inline const typename DataTpl<Scalar,Options,JointCollectionTpl>::TangentVectorType &
  azamat_v2(const ModelTpl<Scalar,Options,JointCollectionTpl> & model,
      DataTpl<Scalar,Options,JointCollectionTpl> & data,
      const Eigen::MatrixBase<ConfigVectorType> & q,
      const Eigen::MatrixBase<TangentVectorType1> & v,
      const Eigen::MatrixBase<MatrixType1> & tau_mat,
      const container::aligned_vector<ForceDerived> & fext)

  {
    assert(model.check(data) && "data is not consistent with model.");
    PINOCCHIO_CHECK_ARGUMENT_SIZE(q.size(), model.nq, "The joint configuration vector is not of right size");
    PINOCCHIO_CHECK_ARGUMENT_SIZE(v.size(), model.nv, "The joint velocity vector is not of right size");
    //PINOCCHIO_CHECK_ARGUMENT_SIZE(tau_mat.cols(),model.nv);
    //PINOCCHIO_CHECK_ARGUMENT_SIZE(tau_mat.rows(),model.nv);
    
    typedef typename ModelTpl<Scalar,Options,JointCollectionTpl>::JointIndex JointIndex;
    
    // data.v[0].setZero();
    // data.a_gf[0].setZero();
    
    // typedef AzamatForwardStep1<Scalar,Options,JointCollectionTpl,ConfigVectorType,TangentVectorType1> Pass1;
    // for(JointIndex i=1;i<(JointIndex)model.njoints;++i)
    // {
    //   Pass1::run(model.joints[i],data.joints[i],
    //              typename Pass1::ArgsType(model,data,q.derived(),v.derived()));
    //   data.f[i] -= fext[i];
    // }
    
    // typedef AzamatBackwardStep<Scalar,Options,JointCollectionTpl> Pass2;
 
    // typedef AzamatForwardStep2<Scalar,Options,JointCollectionTpl> Pass3;

    //   for (int ii=0; ii=(JointIndex)model.njoints;++ii)
    // {
    //     data.u = tau_mat.row(ii);

    //     for(JointIndex i=(JointIndex)model.njoints-1;i>0; --i)
    //     {
    //     Pass2::run(model.joints[i],data.joints[i],
    //                 typename Pass2::ArgsType(model,data));
    //     }

    //     for(JointIndex i=1; i<(JointIndex)model.njoints; ++i)
    //     {
    //     Pass3::run(model.joints[i],data.joints[i],
    //                 typename Pass3::ArgsType(model,data));
    //     }

    // }

    return data.Minv_mat_prod_v2;
  }
  


  // --- CHECKER ---------------------------------------------------------------
  // --- CHECKER ---------------------------------------------------------------
  // --- CHECKER ---------------------------------------------------------------

  // Check whether all masses are nonzero and diagonal of inertia is nonzero
  // The second test is overconstraining.
//   template<typename Scalar, int Options, template<typename,int> class JointCollectionTpl>
//   inline bool AZAmatChecker::checkModel_impl(const ModelTpl<Scalar,Options,JointCollectionTpl> & model) const
//   {
//     typedef ModelTpl<Scalar,Options,JointCollectionTpl> Model;
//     typedef typename Model::JointIndex JointIndex;
    
//     for(JointIndex j=1;j<(JointIndex)model.njoints;j++)
//       if(    (model.inertias[j].mass   ()           < 1e-5) 
//           || (model.inertias[j].inertia().data()[0] < 1e-5)
//           || (model.inertias[j].inertia().data()[3] < 1e-5)
//           || (model.inertias[j].inertia().data()[5] < 1e-5) )
//         return false;
//     return true;
//   }

} // namespace pinocchio

/// @endcond

#endif //