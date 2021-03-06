//
// Copyright (c) 2018-2020 CNRS INRIA
//
// Modified by Shubham Singh for testing different robot configurations

#include "pinocchio/algorithm/joint-configuration.hpp"
#include "pinocchio/algorithm/rnea-derivatives.hpp"
#include "pinocchio/algorithm/rnea-derivatives-faster.hpp"
#include "pinocchio/algorithm/cholesky.hpp"
#include "pinocchio/parsers/urdf.hpp"
#include "pinocchio/parsers/sample-models.hpp"
#include "pinocchio/container/aligned-vector.hpp"

#include <iostream>

#include "pinocchio/utils/timer.hpp"
#include <string>
#include <iostream>
#include <ctime>

using namespace std;

int main(int argc, const char ** argv)
{
  using namespace Eigen;
  using namespace pinocchio;

  PinocchioTicToc timer(PinocchioTicToc::US);
  #ifdef NDEBUG
  const int NBT = 50000;
  #else
    const int NBT = 1;
    std::cout << "(the time score in debug mode is not relevant) " << std::endl;
  #endif
    
    cout << "value of NBT is" << NBT << "\n \n";

  Model model;

// Saving data to a file

  string txt_file_str;

  string str_urdf[10];

  str_urdf[0] = "/home/ss86299/Desktop/test_pinocchio/pinocchio/models/double_pendulum.urdf"; // double pendulum
  str_urdf[1] = "/home/ss86299/Desktop/test_pinocchio/pinocchio/models/ur3_robot.urdf"; // ur_3
  str_urdf[2] = "/home/ss86299/Desktop/test_pinocchio/pinocchio/models/ur5_robot.urdf"; // ur_5
  str_urdf[3] = "/home/ss86299/Desktop/test_pinocchio/pinocchio/models/kuka_peg_lwr.urdf"; // KUKA_lWR
  str_urdf[4] = "/home/ss86299/Desktop/test_pinocchio/pinocchio/models/ur10_robot.urdf"; // ur_10
  str_urdf[5] = "/home/ss86299/Desktop/test_pinocchio/pinocchio/models/anymal.urdf"; // anyMal
  str_urdf[6] = "/home/ss86299/Desktop/test_pinocchio/pinocchio/models/hyq.urdf"; // HYQ
  str_urdf[7] = "/home/ss86299/Desktop/test_pinocchio/pinocchio/models/simple_humanoid.urdf"; // simple humanoid
  str_urdf[8] = "/home/ss86299/Desktop/test_pinocchio/pinocchio/models/atlas.urdf"; //atlas 
  str_urdf[9] = "/home/ss86299/Desktop/test_pinocchio/pinocchio/models/baxter_simple.urdf"; //baxter 

  std :: string filename = str_urdf[0];


   cout << " model is " << filename << endl;

   txt_file_str = "atlas";
   ofstream myfile;
   myfile.open(txt_file_str);

  //-------------------------------------------------------------------

    if(argc>1) filename = argv[1];
    bool with_ff = false; // true originally
    
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

    cout << "with ff = " << with_ff <<endl;    
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


    // Writing to myfile

    //myfile << "test here" << "\n \n";
    myfile.close();

    //cout << "print timing here by new method" << difftime(tend, tstart)/NBT;
    
    timer.tic();
    SMOOTH(NBT)
    {
      computeRNEADerivativesFaster(model,data,qs[_smooth],qdots[_smooth],qddots[_smooth],
                            drnea2_dq,drnea2_dv,drnea_da);
    }
    std::cout << "RNEA derivativeF= \t\t"; timer.toc(std::cout,NBT);

    // Writing to myfile
    // myfile << timer.toc(std::cout,NBT) << "\n \n";

    eq = drnea_dq - drnea2_dq;
    std::cout << eq.topRows(6).squaredNorm() << std::endl;
    std::cout << eq.squaredNorm() << std::endl;

    double upErr= 0;
    double downErr = 0;
    for (int ii = 0 ; ii < model.nv ; ii+=1)
    {
      for(int jj = 0 ; jj <= ii ; jj+=1)
      {
        upErr+=eq(ii,jj)*eq(ii,jj);
        downErr+=eq(jj,ii)*eq(jj,ii);

      }
    }

    std::cout <<"up " << upErr << std::endl;
    std::cout <<"down " << downErr << std::endl;


    ev = drnea_dv - drnea2_dv;

    std::cout << ev.topRows(6).squaredNorm() << std::endl;
    std::cout << ev.squaredNorm() << std::endl;

    std::cout << "--" << std::endl;



  return 0;

}
