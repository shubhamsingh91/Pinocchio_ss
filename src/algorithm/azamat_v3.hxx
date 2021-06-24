//
// Copyright (c) 2016-2020 CNRS, INRIA
//

#ifndef __pinocchio_azamat_v3_hxx__
#define __pinocchio_azamat_v3_hxx__

#include "pinocchio/spatial/act-on-set.hpp"
#include "pinocchio/multibody/visitor.hpp"
#include "pinocchio/algorithm/check.hpp"
#include <iostream>
#include <stdlib.h>
#include <string>

/// @cond DEV

namespace pinocchio
{

  namespace internal
  
  {
    
    template<typename Scalar>
    struct SE3actOn_azamat_v3
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
  
  template<typename Scalar, int Options, template<typename,int> class JointCollectionTpl, typename ConfigVectorType>
  struct azamat_v3ForwardStep1
  : public fusion::JointUnaryVisitorBase< azamat_v3ForwardStep1<Scalar,Options,JointCollectionTpl,ConfigVectorType> >
  {
    typedef ModelTpl<Scalar,Options,JointCollectionTpl> Model;
    typedef DataTpl<Scalar,Options,JointCollectionTpl> Data;
    
    typedef boost::fusion::vector<const Model &,
                                  Data &,
                                  const ConfigVectorType &
                                  > ArgsType;
    
    template<typename JointModel>
    static void algo(const pinocchio::JointModelBase<JointModel> & jmodel,
                     pinocchio::JointDataBase<typename JointModel::JointDataDerived> & jdata,
                     const Model & model,
                     Data & data,
                     const Eigen::MatrixBase<ConfigVectorType> & q)
    {
      typedef typename Model::JointIndex JointIndex;
      
      const JointIndex & i = jmodel.id();
      jmodel.calc(jdata.derived(),q.derived());
      
      const JointIndex & parent = model.parents[i];
      data.liMi[i] = model.jointPlacements[i] * jdata.M();
      
      if (parent>0)
        data.oMi[i] = data.oMi[parent] * data.liMi[i];
      else
        data.oMi[i] = data.liMi[i];
      
      typedef typename SizeDepType<JointModel::NV>::template ColsReturn<typename Data::Matrix6x>::Type ColsBlock;
      ColsBlock J_cols = jmodel.jointCols(data.J);
      J_cols = data.oMi[i].act(jdata.S());
      
      data.Yaba[i] = model.inertias[i].matrix();
    }
    
  };
  
  template<typename Scalar, int Options, template<typename,int> class JointCollectionTpl>
  struct azamat_v3BackwardStep
  : public fusion::JointUnaryVisitorBase< azamat_v3BackwardStep<Scalar,Options,JointCollectionTpl> >
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
      
      const JointIndex & i = jmodel.id();
      const JointIndex & parent  = model.parents[i];
      
      typename Inertia::Matrix6 & Ia = data.Yaba[i];
      typename Data::RowMatrixXs & qdd_mat = data.Minv_mat_prod_v3;
      typename Data::Matrix6x & Fcrb = data.Fcrb_v2[0]; // first element of the array
      typename Data::Matrix6x & FcrbTmp = data.Fcrb_v2.back(); // last element of the array

      
      std::cout << "---------- Loop 2 variables here--------------" << std::endl;
      std::cout << "jmodel.nv is" << jmodel.nv() << std::endl;
      std::cout << "i is" << i << std::endl;
      std::cout << "jmodel.idx_v() is" << jmodel.idx_v() << std::endl;

      std::cout << "parent is" << parent << std::endl;

      std::cout <<"Fcrb here in second loop is" << Fcrb << std:: endl;
      std::cout <<"FcrbTmp here in second loop is" << FcrbTmp << std:: endl;

      jmodel.calc_aba(jdata.derived(), Ia, parent > 0);
      
      typedef typename SizeDepType<JointModel::NV>::template ColsReturn<typename Data::Matrix6x>::Type ColsBlock;

     // U{i} calculation
     
      ColsBlock U_cols = jmodel.jointCols(data.IS);
      forceSet::se3Action(data.oMi[i],jdata.U(),U_cols); //U{i} is expressed in the world frame
     
     // SD_inv calculation
       
        ColsBlock J_cols = jmodel.jointCols(data.J);
        ColsBlock SDinv_cols = jmodel.jointCols(data.SDinv);
        SDinv_cols.noalias() = J_cols * jdata.Dinv();

       // std::cout << "SDinv_cols is" << SDinv_cols << std::endl;
       // std::cout << "SDinv_cols transpose is" << SDinv_cols.transpose() << std::endl;   
       // std::cout << "J_cols here is" << J_cols << std::endl;
       // std::cout << "Fcrb_v2 here is" << Fcrb << std::endl;
       // std::cout << "qdd_mat here is" << qdd_mat << std::endl;

     //-- First line    

       // qdd_mat.row(i)  = -SDinv_cols.transpose()*Fcrb; // didn't work
       
       // std::cout  << "parent is is" << parent << std:: endl;
       // std::cout  << "i is" << i << std:: endl;
       
       // std:: cout <<  "Fcrb is" << Fcrb.middleCols(0,2*model.nv) << std::endl;

       // std::cout << "jdata.Dinv() is" << jdata.Dinv() << std::endl;

       // i-1 is replaced by jmodel.idx_v()

       //-- original expression--------- 

       //  qdd_mat.block(i-1,0,1,2*model.nv).noalias() = jdata.Dinv()*data.tau_mat_v2.block(i-1,0,1,2*model.nv)
       //     - SDinv_cols.transpose() * Fcrb.middleCols(0,2*model.nv);

       //-- alternate expression--------- 


          qdd_mat.block(jmodel.idx_v(),0,jmodel.nv(),2*model.nv).noalias() = 
          jdata.Dinv()*data.tau_mat_v2.block(jmodel.idx_v(),0,jmodel.nv(),2*model.nv)
            - SDinv_cols.transpose() * Fcrb.middleCols(0,2*model.nv);

          std::cout << "jdata.Dinv() here is" << jdata.Dinv() << std::endl;
          std::cout << "SDinv_cols here is" << SDinv_cols.transpose() << std::endl;
          std::cout << "qdd_mat here is" << qdd_mat << std::endl;

       // qdd_mat.block(i-1,0,1,2*model.nv).noalias() = jdata.Dinv()*data.tau_mat_v2.block(i-1,0,1,2*model.nv);

       // std::cout << "i here is" << i << std::endl;
       // std::cout << "qdd block matrix here is" <<   qdd_mat.block(i-1,0,1,2*model.nv) << std::endl;
       // std::cout << "Fcrb middle cols here is"<<Fcrb.middleCols(0,2*model.nv) << std::endl;

         if(parent > 0)
          {
          
            FcrbTmp.noalias() = U_cols* qdd_mat.block(jmodel.idx_v(),0,jmodel.nv(),2*model.nv);

            std::cout << "FcrbTmp in parent loop is" << FcrbTmp << std::endl;

              // std::cout << "FcrbTmp here is " << FcrbTmp << std::endl;

              // data.Fcrb_v2[parent]  += FcrbTmp; // this also works here

              // Fcrb  += FcrbTmp; // this also works here
      
            Fcrb.middleCols(0,2*model.nv)  += FcrbTmp.middleCols(0,2*model.nv);

            std::cout << "Fcrb in parent loop is" << Fcrb << std::endl;

          }
        
        if(parent > 0)
            data.Yaba[parent] += internal::SE3actOn_azamat_v3<Scalar>::run(data.liMi[i], Ia);
    
     }
    
  };
  
  template<typename Scalar, int Options, template<typename,int> class JointCollectionTpl>
  struct azamat_v3ForwardStep2
  : public fusion::JointUnaryVisitorBase< azamat_v3ForwardStep2<Scalar,Options,JointCollectionTpl> >
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
      typename Data::RowMatrixXs & qdd_mat = data.Minv_mat_prod_v3;
      typename Data::Matrix6x & FcrbTmp = data.Fcrb_v2.back();
      
      typedef typename SizeDepType<JointModel::NV>::template ColsReturn<typename Data::Matrix6x>::Type ColsBlock;
      ColsBlock UDinv_cols = jmodel.jointCols(data.UDinv);
      forceSet::se3Action(data.oMi[i],jdata.UDinv(),UDinv_cols); // expressed in the world frame
      ColsBlock J_cols = jmodel.jointCols(data.J);

      // std::cout << "FcrbTmp is" << FcrbTmp << std::endl;
      // std::cout << "data.Fcrb_v2[parent] is" << data.Fcrb_v2[parent] << std::endl;

      // std::cout << "i here in third loop is" << i << std::endl;
      // std::cout << "jmodel.idx_v() here in third loop is" << jmodel.idx_v() << std::endl;
      
      std::cout << "---------- Loop 3 variables here--------------" << std::endl;
      std::cout << "jmodel.nv is" << jmodel.nv() << std::endl;
      std::cout << "i is" << i << std::endl;
      std::cout << "jmodel.idx_v() is" << jmodel.idx_v() << std::endl;
      std::cout << "parent is" << parent << std::endl;

      if(parent > 0)
      {

       // std::cout << "data.Fcrb[parent] in third loop is" << data.Fcrb_v2[parent] << std::endl;
       // std::cout << "data.Fcrb_v2 is" << data.Fcrb_v2 << std::endl;

        std::cout << "data.Fcrb_v2[parent] in third loop is" << data.Fcrb_v2[parent] << std::endl;

        FcrbTmp.block(jmodel.idx_v(),0,jmodel.nv(),2*model.nv).noalias() = UDinv_cols.transpose() * data.Fcrb_v2[parent];
        qdd_mat.block(jmodel.idx_v(),0,jmodel.nv(),2*model.nv).noalias() -=  FcrbTmp.block(jmodel.idx_v(),0,jmodel.nv(),2*model.nv);

      }
          data.Fcrb_v2[i] = J_cols * qdd_mat.block(jmodel.idx_v(),0,jmodel.nv(),2*model.nv);

          std::cout << " data.Fcrb_v2[i] in third loop is" << data.Fcrb_v2[i] << std::endl;
    
          if(parent > 0)
            data.Fcrb_v2[i] += data.Fcrb_v2[parent];

     }
    
  };

  template<typename Scalar, int Options, template<typename,int> class JointCollectionTpl, typename ConfigVectorType,  typename MatrixType1>
  inline const typename DataTpl<Scalar,Options,JointCollectionTpl>::RowMatrixXs &
  azamat_v3(const ModelTpl<Scalar,Options,JointCollectionTpl> & model,
                  DataTpl<Scalar,Options,JointCollectionTpl> & data,
                  const Eigen::MatrixBase<ConfigVectorType> & q,
                  const Eigen::MatrixBase<MatrixType1> & tau_mat)


  {
    assert(model.check(data) && "data is not consistent with model.");
    PINOCCHIO_CHECK_ARGUMENT_SIZE(q.size(), model.nq, "The joint configuration vector is not of right size");
    
    typedef typename ModelTpl<Scalar,Options,JointCollectionTpl>::JointIndex JointIndex;
    
    typedef azamat_v3ForwardStep1<Scalar,Options,JointCollectionTpl,ConfigVectorType> Pass1;
    for(JointIndex i=1; i<(JointIndex)model.njoints; ++i)
    {
      Pass1::run(model.joints[i],data.joints[i],
                 typename Pass1::ArgsType(model,data,q.derived()));
    }
    
    data.Fcrb_v2[0].setZero();
    
    //std::cout << "--------------------------------------------" << std::endl;
    //std::cout << "data.Fcrb_v2 [0] is " <<  data.Fcrb_v2[0] << std::endl;

    data.tau_mat_v2 = tau_mat;

    typedef azamat_v3BackwardStep<Scalar,Options,JointCollectionTpl> Pass2;
    for(JointIndex i=(JointIndex)model.njoints-1; i>0; --i)
    {
      Pass2::run(model.joints[i],data.joints[i],
                 typename Pass2::ArgsType(model,data));
    }   

    typedef azamat_v3ForwardStep2<Scalar,Options,JointCollectionTpl> Pass3;
    for(JointIndex i=1; i<(JointIndex)model.njoints; ++i)
    {
      Pass3::run(model.joints[i],data.joints[i],
                 typename Pass3::ArgsType(model,data));
    }
    
    return data.Minv_mat_prod_v3;
  }


    // --- CHECKER ---------------------------------------------------------------
    // --- CHECKER ---------------------------------------------------------------
    // --- CHECKER ---------------------------------------------------------------

    // Check whether all masses are nonzero and diagonal of inertia is nonzero
    // The second test is overconstraining.
    //   template<typename Scalar, int Options, template<typename,int> class JointCollectionTpl>
    //   inline bool ABAChecker::checkModel_impl(const ModelTpl<Scalar,Options,JointCollectionTpl> & model) const
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

#endif 
