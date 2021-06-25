
// checking some boost options, ABA forward dynamics etc.
#include "pinocchio/algorithm/joint-configuration.hpp"
#include "pinocchio/algorithm/cholesky.hpp"
#include "pinocchio/parsers/urdf.hpp"
#include "pinocchio/parsers/sample-models.hpp"
#include "pinocchio/container/aligned-vector.hpp"
#include "pinocchio/algorithm/azamat.hpp"
#include "pinocchio/algorithm/azamat_v2.hpp"
#include "pinocchio/algorithm/azamat_v3.hpp"
#include "pinocchio/algorithm/azamat_v4.hpp"
#include "pinocchio/utils/timer.hpp"
#include "pinocchio/algorithm/aba-derivatives.hpp"
#include "pinocchio/algorithm/rnea-derivatives.hpp"
#include "pinocchio/algorithm/azamat_v4.hpp"
#include "pinocchio/algorithm/rnea-derivatives-faster.hpp"

#include <iostream>
#include <fstream>

#include <string>
#include <iostream>
#include <ctime>

//-----------------Running ABA  first---- //

using namespace std;

int main(int argc, const char ** argv)
{
  using namespace Eigen;
  using namespace pinocchio;

  PinocchioTicToc timer(PinocchioTicToc::US);

  #ifdef NDEBUG
  const int NBT = 1; // 50000 initially
  #else
    const int NBT = 1;
    std::cout << "(the time score in debug mode is not relevant) " << std::endl;
  #endif

    Model model;

    string str_urdf[5];
    int str_int;

   str_int = 0; // integer to change the str_urdf

   // cout << "Enter the str_int here" << endl;
   // cin >> str_int;

    str_urdf[0] = "/home/ss86299/Desktop/test_pinocchio/pinocchio/models/double_pendulum_v1.urdf"; // double pendulum
    str_urdf[1] = "/home/ss86299/Desktop/test_pinocchio/pinocchio/models/ur3_robot.urdf"; // UR3
    str_urdf[2] = "/home/ss86299/Desktop/test_pinocchio/pinocchio/models/ur5_robot.urdf"; // UR5
    str_urdf[3] = "/home/ss86299/Desktop/test_pinocchio/pinocchio/models/ur10_robot.urdf"; // UR10
    str_urdf[4] = "/home/ss86299/Desktop/test_pinocchio/pinocchio/models/kuka_peg_lwr.urdf"; // kuka   
    str_urdf[5] = "/home/ss86299/Desktop/test_pinocchio/pinocchio/models/anymal.urdf"; // anymal
    str_urdf[6] = "/home/ss86299/Desktop/test_pinocchio/pinocchio/models/hyq.urdf"; // hyq   
    str_urdf[7] = "/home/ss86299/Desktop/test_pinocchio/pinocchio/models/baxter_simple.urdf"; // baxter
    str_urdf[8] = "/home/ss86299/Desktop/test_pinocchio/pinocchio/models/atlas.urdf"; //atlas
    str_urdf[9] = "/home/ss86299/Desktop/test_pinocchio/pinocchio/models/simple_humanoid.urdf"; //simple humanoid

    std :: string filename = str_urdf[str_int];

    if(argc>1) filename = argv[1];
    bool with_ff = true; // true originally
    
    if ((str_int >4)&&(str_int<11))
    {
      with_ff = true; // True for anymal and atlas
    }

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

    cout << "with ff = " << with_ff << "\n \n";    
    std::cout << "nq = " << model.nq << std::endl;
    std::cout << "nv = " << model.nv << std::endl;
    cout << "Model is" << str_urdf[str_int] << endl;

    Data data(model);
    VectorXd qmax = Eigen::VectorXd::Ones(model.nq);

    PINOCCHIO_ALIGNED_STD_VECTOR(VectorXd) qs     (NBT);
    PINOCCHIO_ALIGNED_STD_VECTOR(VectorXd) qdots  (NBT);
    PINOCCHIO_ALIGNED_STD_VECTOR(VectorXd) taus (NBT); 
    PINOCCHIO_ALIGNED_STD_VECTOR(VectorXd) qddots (NBT); 
    PINOCCHIO_ALIGNED_STD_VECTOR(MatrixXd) tau_mat (NBT); 
    PINOCCHIO_ALIGNED_STD_VECTOR(MatrixXd) tau_mat_n2n (NBT); 
    PINOCCHIO_ALIGNED_STD_VECTOR(MatrixXd) tau_mat_n2n_v2 (NBT); 
    PINOCCHIO_ALIGNED_STD_VECTOR(MatrixXd) tau_mat_n2n_v3 (NBT); 

  // MatrixXd tau_mat(MatrixXd::Identity(model.nv,model.nv));

  // randomizing input data here

    std:: cout << "NBT variable is" << NBT << endl;

    for(size_t i=0;i<NBT;++i)
    {
      qs[i]     = randomConfiguration(model,-qmax,qmax);
      qdots[i]  = Eigen::VectorXd::Random(model.nv);
      taus[i] = Eigen::VectorXd::Random(model.nv);
      qddots[i] =  Eigen::VectorXd::Random(model.nv);
      tau_mat_n2n[i] = Eigen::MatrixXd::Identity(model.nv,2*model.nv);
      tau_mat_n2n_v2[i] = Eigen::MatrixXd::Identity(model.nv,2*model.nv);
      tau_mat_n2n_v3[i] = Eigen::MatrixXd::Identity(model.nv,2*model.nv);

    }   

 // std::cout << "tau_mat input here is" << tau_mat_n2n[0] << std::endl;

    //----------------------------------------------------//
    // Compute RNEA derivatives here ---------------------//
    //----------------------------------------------------//

    PINOCCHIO_EIGEN_PLAIN_ROW_MAJOR_TYPE(MatrixXd) drnea_dq(MatrixXd::Zero(model.nv,model.nv));
    PINOCCHIO_EIGEN_PLAIN_ROW_MAJOR_TYPE(MatrixXd) drnea_dv(MatrixXd::Zero(model.nv,model.nv));
    MatrixXd drnea_da(MatrixXd::Zero(model.nv,model.nv));

    timer.tic();
    SMOOTH(NBT)
    {
        computeRNEADerivatives(model,data,qs[_smooth],qdots[_smooth],qddots[_smooth],
                            drnea_dq,drnea_dv,drnea_da);
    }

     std::cout << "RNEA derivatives= \t\t"; timer.toc(std::cout,NBT);

    //----------------------------------------------------//
    // Compute RNEA derivatives faster--------------------//
    //----------------------------------------------------//

    PINOCCHIO_EIGEN_PLAIN_ROW_MAJOR_TYPE(MatrixXd) drnea2_dq(MatrixXd::Zero(model.nv,model.nv));
    PINOCCHIO_EIGEN_PLAIN_ROW_MAJOR_TYPE(MatrixXd) drnea2_dv(MatrixXd::Zero(model.nv,model.nv));
    MatrixXd drnea2_da(MatrixXd::Zero(model.nv,model.nv));

    timer.tic();
    SMOOTH(NBT)
    {
        computeRNEADerivativesFaster(model,data,qs[_smooth],qdots[_smooth],qddots[_smooth],
                            drnea2_dq,drnea2_dv,drnea2_da);
    }
    std::cout << "RNEA derivativeF= \t\t" ; timer.toc(std::cout,NBT);

    //-----------------------------------------------------------------//
    // FD partials using AZAmat_v3/4 function here---------------------//
    //-----------------------------------------------------------------//

    tau_mat_n2n_v3[0] << -drnea_dq,-drnea_dv; // concatenating partial wrt q and qdot

   // std::cout<< "tau_mat_n2n_v3 is" << tau_mat_n2n_v3[0] << std::endl;

    timer.tic();
    SMOOTH(NBT)
    azamat_v4(model,data,qs[_smooth],tau_mat_n2n_v3[_smooth]);
    std::cout << "IPR using AZA_mat_v4 method is = \t\t" ;timer.toc(std::cout,NBT);

    std::cout <<"---------------------------------------" << endl;

    //-----------------------------------------------------//
    //------------- Minv here------------------------------//
    //-----------------------------------------------------//

    timer.tic();
    SMOOTH(NBT)
    {
       computeMinverse(model,data,qs[_smooth]);
    }
    std::cout << "Minv =\t\t"; timer.toc(std::cout,NBT);

    // Filling the Minv matrix here

    for (int ii=0; ii < model.nv ;ii++)
    {
      for (int jj=0; jj<ii; jj++)
          {
              if (ii!=jj)
              {
                  data.Minv.coeffRef(ii,jj) = data.Minv.coeffRef(jj,ii);
              }
          }
    }

    //----------------------------------------------------//
    // Calculate ABA derivatives here
    //----------------------------------------------------//

    PINOCCHIO_EIGEN_PLAIN_ROW_MAJOR_TYPE(MatrixXd) daba_dq(MatrixXd::Zero(model.nv,model.nv));
    PINOCCHIO_EIGEN_PLAIN_ROW_MAJOR_TYPE(MatrixXd) daba_dv(MatrixXd::Zero(model.nv,model.nv));
    MatrixXd daba_dtau(MatrixXd::Zero(model.nv,model.nv));

    taus[0] = data.tau;

    timer.tic();
    SMOOTH(NBT)
    {
        computeABADerivatives(model,data,qs[_smooth],qdots[_smooth],taus[_smooth],
                            daba_dq,daba_dv,daba_dtau);
    }
     std::cout << "ABA derivatives= \t\t"; timer.toc(std::cout,NBT);

    //-----------------------------------------------------------------//
    // FD partials using DMM here -------------------------------------//
    //-----------------------------------------------------------------//
   
    MatrixXd aba_partial_dq(MatrixXd::Zero(model.nv,model.nv));
    MatrixXd aba_partial_dv(MatrixXd::Zero(model.nv,model.nv));
 
    timer.tic();
    SMOOTH(NBT)

    aba_partial_dq = -data.Minv*drnea_dq;
    aba_partial_dv = -data.Minv*drnea_dv;
     std::cout << "DMM= \t\t"; timer.toc(std::cout,NBT);

    //-----------------------------------------------------------------//
    // FD partials using AZAmat function here--------------------------//
    //-----------------------------------------------------------------//

    tau_mat_n2n[0] << -drnea_dq,-drnea_dv; // concatenating partial wrt q and qdot

    timer.tic();
    SMOOTH(NBT)
    azamat(model,data,qs[_smooth],qdots[_smooth],tau_mat_n2n[_smooth]);
    std::cout << "IPR using AZA_mat method is = \t\t" ;timer.toc(std::cout,NBT);
    
   // std::cout << "Minvmat_v1 is" << data.Minv_mat_prod << std::endl;
    std::cout << "---------------------------------------" << endl;



// Comparison of results here


    MatrixXd diff_daba_dq2(MatrixXd::Zero(model.nv,model.nv));
    MatrixXd diff_daba_dqd2(MatrixXd::Zero(model.nv,model.nv));
    MatrixXd diff_mat1(MatrixXd::Zero(model.nv,2*model.nv));

    MatrixXd diff_daba_dq3(MatrixXd::Zero(model.nv,model.nv));
    MatrixXd diff_daba_dqd3(MatrixXd::Zero(model.nv,model.nv));

    // comparison of ABA derivs with DMM result
    diff_daba_dq3 = daba_dq-aba_partial_dq;
    diff_daba_dqd3 = daba_dv-aba_partial_dv;


    // comparison of ABA derivs with AZAMAT_v4

    diff_daba_dq2 = daba_dq-data.Minv_mat_prod_v3.middleCols(0,model.nv);
    diff_daba_dqd2 = daba_dv-data.Minv_mat_prod_v3.middleCols(model.nv,model.nv);

    diff_mat1 = data.Minv_mat_prod - data.Minv_mat_prod_v3;



    std::cout << "------------------------------" << std::endl;

    std::cout << "Norm of the difference matrix for AZAmat_v4 FD partial wrt q from orig FD partial wrt q is " << diff_daba_dq2.squaredNorm() << std::endl;
    std::cout << "Norm of the difference matrix for AZAmat_v4 FD partial wrt qd from orig FD partial wrt qd is " << diff_daba_dqd2.squaredNorm() << std::endl;

    std::cout << "\n" << endl;

    std::cout << "Norm of difference between mat_v1 and mat_v4 is" << diff_mat1.squaredNorm() << std::endl;

    std::cout << "------------------------------" << std::endl;

    std::cout << "Norm of the difference matrix for DMM FD partial wrt q from orig FD partial wrt q is " << diff_daba_dq3.squaredNorm() << std::endl;
    std::cout << "Norm of the difference matrix for DMM FD partial wrt qd from orig FD partial wrt qd is " << diff_daba_dqd3.squaredNorm() << std::endl;

    std::cout << "\n" << endl;



  return 0;

}