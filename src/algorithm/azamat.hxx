//
// Copyright (c) 2016-2020 CNRS, INRIA
//

#ifndef __pinocchio_azamat_hxx__
#define __pinocchio_azamat_hxx__

#include "pinocchio/spatial/act-on-set.hpp"
#include "pinocchio/multibody/visitor.hpp"
#include "pinocchio/algorithm/check.hpp"
#include <iostream>
#include <stdlib.h>

/// @cond DEV

namespace pinocchio
{
  template<typename Scalar, int Options, template<typename,int> class JointCollectionTpl, typename ConfigVectorType, typename TangentVectorType>
  struct AzamatForwardStep1
  : public fusion::JointUnaryVisitorBase< AzamatForwardStep1<Scalar,Options,JointCollectionTpl,ConfigVectorType,TangentVectorType> >
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
      //jmodel.calc(jdata.derived(),q.derived(),v.derived()); // joint calcs, not needed at all
      
      // const JointIndex & parent = model.parents[i];
      // data.liMi[i] = model.jointPlacements[i] * jdata.M(); // not needed here, working in the "new" function
      
      // Not needed here, since vj is zero

      // data.v[i] = jdata.v();
      // data.v[i].setZero();
       
      // std:: cout << "AZA data here----------------------------------" << std::endl;
      
      // std:: cout << "v(i) is" << data.v[i] <<  std::endl;

      // if (parent>0)
      // data.v[i] += data.liMi[i].actInv(data.v[parent]);

      // std::cout << "v(i) after if loop is" << data.v[i] << std::endl;

      // data.a_gf[i] = jdata.c() + (data.v[i] ^ jdata.v());

       data.a_gf_v1[i].setZero();

      //  std:: cout << "jdata.c() is" << jdata.c() <<  std::endl;
  
      //  std:: cout << "a_gf(i) is" << data.a_gf[i] <<  std::endl;
  
       data.Yaba_v1[i] = model.inertias[i].matrix();
       data.f_v1[i].setZero(); // -f_ext

     //  data.f[i] = model.inertias[i].vxiv(data.v[i]); // -f_ext

     //   std:: cout << "f(i) is" << data.f[i] <<  std::endl;
     //   std:: cout << "size of f(i) is" << sizeof(data.f[i]) <<  std::endl;
     //   std:: cout << "type of f(i) is" << typeof(data.f[i]) <<  std::endl;

    }
    
  };


  namespace internal
  {
    
    template<typename Scalar>
    struct SE3actOnazamat
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

  template<typename Scalar, int Options, template<typename,int> class JointCollectionTpl, typename ConfigVectorType, typename TangentVectorType>
  struct AzamatForwardStep1_new
  : public fusion::JointUnaryVisitorBase< AzamatForwardStep1_new<Scalar,Options,JointCollectionTpl,ConfigVectorType,TangentVectorType> >
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
      const JointIndex & parent  = model.parents[i];
      typename Inertia::Matrix6 & Ia = data.Yaba[i];

       data.liMi[i] = model.jointPlacements[i] * jdata.M();
     
     if (parent > 0)
      {
        data.Yaba[parent] += internal::SE3actOnazamat<Scalar>::run(data.liMi[i], Ia);

        
      }

    }
    
  };
  
  template<typename Scalar, int Options, template<typename,int> class JointCollectionTpl>
  struct AzamatBackwardStep
  : public fusion::JointUnaryVisitorBase< AzamatBackwardStep<Scalar,Options,JointCollectionTpl> >
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
    
    // data.u is already initialized as tau 
    
      jmodel.jointVelocitySelector(data.u) -= jdata.S().transpose()*data.f[i];

     // jmodel.calc_aba(jdata.derived(), Ia, parent > 0);

      if (parent > 0)
      {
        Force & pa = data.f[i]; 

        //pa.toVector() += Ia * data.a_gf[i].toVector() + jdata.UDinv() * jmodel.jointVelocitySelector(data.u);

        pa.toVector() +=  jdata.UDinv() * jmodel.jointVelocitySelector(data.u);
    
        //data.Yaba[parent] += internal::SE3actOnazamat<Scalar>::run(data.liMi[i], Ia);
        data.f[parent] += data.liMi[i].act(pa);

       // std:: cout << "Testing data.f[parent] here" << std::endl;
       // std:: cout << "f(parent) is" << data.f[parent] <<  std::endl;
        
      }
    }
    
  }; 
  
  template<typename Scalar, int Options, template<typename,int> class JointCollectionTpl>
  struct AzamatForwardStep2
  : public fusion::JointUnaryVisitorBase< AzamatForwardStep2<Scalar,Options,JointCollectionTpl> >
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
        
        //   std::cout <<"Testing loop 3 here ---------------------" << std:: endl;

        data.a_gf[i] += data.liMi[i].actInv(data.a_gf[parent]);

        //   std::cout <<"data.a_gf[i] here is " << data.a_gf[i] << std:: endl;

        jmodel.jointVelocitySelector(data.ddq_new).noalias() =
        jdata.Dinv() * jmodel.jointVelocitySelector(data.u) - jdata.UDinv().transpose() * data.a_gf[i].toVector();
       
        data.a_gf[i] += jdata.S() * jmodel.jointVelocitySelector(data.ddq_new);
        

      //   std::cout <<"Minv_mat_prod here is " << data.ddq_new << std:: endl;


    }
    
  };
  
  template<typename Scalar, int Options, template<typename,int> class JointCollectionTpl, typename ConfigVectorType, typename TangentVectorType1, typename MatrixType1>
  inline const typename DataTpl<Scalar,Options,JointCollectionTpl>::TangentVectorType &
  azamat(const ModelTpl<Scalar,Options,JointCollectionTpl> & model,
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
    
    data.v[0].setZero();
    data.a_gf[0].setZero();
    
    // data.a_gf[0] = -model.gravity;
    // data.u = tau;
    // data.f.setZero();    // didn't work
    // data.a_gf.setZero(); // didn't work

    typedef AzamatForwardStep1_new<Scalar,Options,JointCollectionTpl,ConfigVectorType,TangentVectorType1> Pass1_new;
  
    typedef AzamatForwardStep1<Scalar,Options,JointCollectionTpl,ConfigVectorType,TangentVectorType1> Pass1;
    
    typedef AzamatBackwardStep<Scalar,Options,JointCollectionTpl> Pass2;

    typedef AzamatForwardStep2<Scalar,Options,JointCollectionTpl> Pass3;

     //------------- Pass 1--------------------------------//

      for(JointIndex i=1; i<(JointIndex)model.njoints; ++i)
      {
        Pass1::run(model.joints[i],data.joints[i],
                  typename Pass1::ArgsType(model,data,q.derived(),v.derived()));
      }
        data.Yaba = data.Yaba_v1;
   
  // Need this only once

    for(JointIndex i=1; i<(JointIndex)model.njoints; ++i)
      {
        Pass1_new::run(model.joints[i],data.joints[i],
                  typename Pass1_new::ArgsType(model,data,q.derived(),v.derived()));
      }



    for (JointIndex ii=1; ii<2*model.nv+1; ++ii)
     {
      //   std::cout <<"model.nv  is" << model.nv  << std::endl;
      //   std::cout <<"(JointIndex)model.njoints  is" << (JointIndex)model.njoints  << std::endl;

     // setting variables reqd to zero again

     data.f = data.f_v1;
     data.a_gf = data.a_gf_v1;
     data.Yaba = data.Yaba_v1;

      //   std::cout <<"ii index is " << ii << std::endl;

        data.u = tau_mat.col(ii-1);

       //------------- Pass 2--------------------------------//

      //  std::cout << "data.u is" << data.u << std::endl;

        for(JointIndex i=(JointIndex)model.njoints-1;i>0; --i)
        {

       //  std::cout <<"i index in the second loop is" << i << std::endl;

        Pass2::run(model.joints[i],data.joints[i],
                    typename Pass2::ArgsType(model,data));
        }

       //------------- Pass 3--------------------------------//

        for(JointIndex i=1; i<(JointIndex)model.njoints; ++i)
        {
       //  std::cout <<"i index in the third loop is" << i << std::endl;

        Pass3::run(model.joints[i],data.joints[i],
                    typename Pass3::ArgsType(model,data));

        }

       data.Minv_mat_prod.col(ii-1) = data.ddq_new;

       //std::cout <<"data.ddq_new here is" << data.ddq_new << std::endl;
       //std::cout <<"Minv matrix here is" << data.Minv_mat_prod.col(ii-1) << std::endl;

     }

    
    return data.Minv_mat_prod;
  }

  template<typename Scalar, int Options, template<typename,int> class JointCollectionTpl, typename ConfigVectorType, typename TangentVectorType1, typename MatrixType1, typename ForceDerived>
  inline const typename DataTpl<Scalar,Options,JointCollectionTpl>::TangentVectorType &
  azamat(const ModelTpl<Scalar,Options,JointCollectionTpl> & model,
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
    
    data.v[0].setZero();
//    data.a_gf[0] = -model.gravity;
    data.a_gf[0].setZero();
    
    typedef AzamatForwardStep1<Scalar,Options,JointCollectionTpl,ConfigVectorType,TangentVectorType1> Pass1;
    for(JointIndex i=1;i<(JointIndex)model.njoints;++i)
    {
      Pass1::run(model.joints[i],data.joints[i],
                 typename Pass1::ArgsType(model,data,q.derived(),v.derived()));
      data.f[i] -= fext[i];
    }
    
    typedef AzamatBackwardStep<Scalar,Options,JointCollectionTpl> Pass2;
 
    typedef AzamatForwardStep2<Scalar,Options,JointCollectionTpl> Pass3;

      for (int ii=0; ii=(JointIndex)model.njoints;++ii)
    {
        data.u = tau_mat.row(ii);

        for(JointIndex i=(JointIndex)model.njoints-1;i>0; --i)
        {
        Pass2::run(model.joints[i],data.joints[i],
                    typename Pass2::ArgsType(model,data));
        }

        for(JointIndex i=1; i<(JointIndex)model.njoints; ++i)
        {
        Pass3::run(model.joints[i],data.joints[i],
                    typename Pass3::ArgsType(model,data));
        }

    }

    return data.Minv_mat_prod;
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

#endif // ifndef __pinocchio_aza_hxx__
