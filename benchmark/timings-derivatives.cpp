//
// Copyright (c) 2018-2020 CNRS INRIA
//

#include "pinocchio/algorithm/joint-configuration.hpp"
#include "pinocchio/algorithm/rnea-derivatives.hpp"
#include "pinocchio/algorithm/rnea-derivatives-faster.hpp"
#include "pinocchio/algorithm/cholesky.hpp"
#include "pinocchio/parsers/urdf.hpp"
#include "pinocchio/parsers/sample-models.hpp"
#include "pinocchio/container/aligned-vector.hpp"

#include <iostream>

#include "pinocchio/utils/timer.hpp"


int main(int argc, const char ** argv)
{
  using namespace Eigen;
  using namespace pinocchio;

  PinocchioTicToc timer(PinocchioTicToc::US);
  #ifdef NDEBUG
  const int NBT = 100;
  #else
    const int NBT = 1;
    std::cout << "(the time score in debug mode is not relevant) " << std::endl;
  #endif
    
  Model model;

  std::string filename = std::string("/home/pwensing/Source/pinocchio/models/simple_humanoid.urdf");
  if(argc>1) filename = argv[1];
  bool with_ff = true;
  
  if(argc>2)
  {
    const std::string ff_option = argv[2];
    if(ff_option == "-no-ff")
      with_ff = false;
  }
    
  if( filename == "HS") 
    buildModels::humanoidRandom(model,true);
  else
    if(with_ff)
      pinocchio::urdf::buildModel(filename,JointModelFreeFlyer(),model);
//      pinocchio::urdf::buildModel(filename,JointModelRX(),model);
    else
      pinocchio::urdf::buildModel(filename,model);
  std::cout << "nq = " << model.nq << std::endl;
  std::cout << "nv = " << model.nv << std::endl;

  Data data(model);
  VectorXd qmax = Eigen::VectorXd::Ones(model.nq);

  PINOCCHIO_ALIGNED_STD_VECTOR(VectorXd) qs     (NBT);
  PINOCCHIO_ALIGNED_STD_VECTOR(VectorXd) qdots  (NBT);
  PINOCCHIO_ALIGNED_STD_VECTOR(VectorXd) qddots (NBT);
  PINOCCHIO_ALIGNED_STD_VECTOR(VectorXd) taus (NBT);
  
  for(size_t i=0;i<NBT;++i)
  {
    qs[i]     = randomConfiguration(model,-qmax,qmax);
    qdots[i]  = Eigen::VectorXd::Random(model.nv);
    qddots[i] = Eigen::VectorXd::Random(model.nv);
    taus[i] = Eigen::VectorXd::Random(model.nv);
  }
  
  PINOCCHIO_EIGEN_PLAIN_ROW_MAJOR_TYPE(MatrixXd) drnea_dq(MatrixXd::Zero(model.nv,model.nv));
  PINOCCHIO_EIGEN_PLAIN_ROW_MAJOR_TYPE(MatrixXd) drnea_dv(MatrixXd::Zero(model.nv,model.nv));
  MatrixXd drnea_da(MatrixXd::Zero(model.nv,model.nv));
 
  MatrixXd drnea2_dq(MatrixXd::Zero(model.nv,model.nv)), eq(MatrixXd::Zero(model.nv,model.nv));
  MatrixXd drnea2_dv(MatrixXd::Zero(model.nv,model.nv)), ev(MatrixXd::Zero(model.nv,model.nv));
  Data::RowMatrixXs daba_dtau(Data::RowMatrixXs::Zero(model.nv,model.nv));
  
  timer.tic();
  SMOOTH(NBT)
  {
    computeRNEADerivatives(model,data,qs[_smooth],qdots[_smooth],qddots[_smooth],
                           drnea_dq,drnea_dv,drnea_da);
  }
  std::cout << "RNEA derivatives= \t\t"; timer.toc(std::cout,NBT);

  timer.tic();
  SMOOTH(NBT)
  {
    computeRNEADerivativesFaster(model,data,qs[_smooth],qdots[_smooth],qddots[_smooth],
                           drnea2_dq,drnea2_dv,drnea_da);
  }
  std::cout << "RNEA derivativeF= \t\t"; timer.toc(std::cout,NBT);

  eq = drnea_dq - drnea2_dq;
  std::cout << eq.squaredNorm() << std::endl;

  ev = drnea_dv - drnea2_dv;

  std::cout << ev.squaredNorm() << std::endl;




  std::cout << "--" << std::endl;
  return 0;
}
