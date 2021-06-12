
// checking some boost options, ABA forward dynamics etc.
#include "pinocchio/algorithm/joint-configuration.hpp"
#include "pinocchio/algorithm/rnea-derivatives.hpp"
#include "pinocchio/algorithm/aba-derivatives.hpp"
#include "pinocchio/algorithm/rnea-derivatives-faster.hpp"
#include "pinocchio/algorithm/cholesky.hpp"
#include "pinocchio/parsers/urdf.hpp"
#include "pinocchio/parsers/sample-models.hpp"
#include "pinocchio/container/aligned-vector.hpp"
#include "pinocchio/algorithm/aba.hpp"
#include "pinocchio/algorithm/aza.hpp"
#include "pinocchio/algorithm/azamat.hpp"
#include "pinocchio/algorithm/azamat_v2.hpp"
#include "pinocchio/algorithm/azamat_v3.hpp"

#include <iostream>
#include <fstream>

#include "pinocchio/utils/timer.hpp"
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

    string str_urdf[15];
    int str_int;

   str_int = 9; // integer to change the str_urdf

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
    str_urdf[10] = "/home/ss86299/Desktop/test_pinocchio/pinocchio/models/simple_humanoid.urdf"; //simple humanoid
    str_urdf[11] = "/home/ss86299/Desktop/test_pinocchio/pinocchio/models/test_robot.urdf"; //test_robot
    str_urdf[12] = "/home/ss86299/Desktop/test_pinocchio/pinocchio/models/2link.urdf"; //2link

    std :: string filename = str_urdf[str_int];

    // Opening a file to write to it

    ofstream file1;
    file1.open("example.txt");

    // Using serial chain-link models here


    // Using usual models here


    if(argc>1) filename = argv[1];
    bool with_ff = false; // true originally
    
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
      tau_mat[i] = Eigen::MatrixXd::Identity(model.nv,model.nv);
      tau_mat_n2n[i] = Eigen::MatrixXd::Identity(model.nv,2*model.nv);
      tau_mat_n2n_v3[i] = Eigen::MatrixXd::Identity(model.nv,2*model.nv);

    }   

 double time_ABA[7];

  timer.tic();
  SMOOTH(NBT)
  aba(model,data,qs[_smooth],qdots[_smooth],taus[_smooth]);
  // time_ABA[0] = timer.toc()/NBT; // ABA timing

  std::cout << "ABA= \t\t"; timer.toc(std::cout,NBT);
  //std::cout << "ABA= \t\t" <<  time_ABA[0] << endl;




  timer.tic();
  SMOOTH(NBT)
  {
    computeMinverse(model,data,qs[_smooth]);
  }
   time_ABA[1] = timer.toc()/NBT; // Minv timing

  //std::cout << "Minv =\t\t"; timer.toc(std::cout,NBT);
  std::cout << "Minv =\t\t" <<  time_ABA[1]<< endl;

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
   time_ABA[2] = timer.toc()/NBT; // RNEA timing

//  std::cout << "RNEA derivatives= \t\t"; timer.toc(std::cout,NBT);
  std::cout << "RNEA derivatives= \t\t" <<  time_ABA[2] << endl;

//----------------------------------------------------//
// Compute RNEA derivatives faster--------------------//
//----------------------------------------------------//

  MatrixXd drnea2_dq(MatrixXd::Zero(model.nv,model.nv));
  MatrixXd drnea2_dv(MatrixXd::Zero(model.nv,model.nv));

timer.tic();
  SMOOTH(NBT)
  {
    computeRNEADerivativesFaster(model,data,qs[_smooth],qdots[_smooth],qddots[_smooth],
                           drnea2_dq,drnea2_dv,drnea_da);
  }
   time_ABA[3] = timer.toc()/NBT; // RNEAF timing

  std::cout << "RNEA derivativeF= \t\t" << time_ABA[3] << endl;

  //std::cout << "RNEA derivativeF= \t\t"; timer.toc(std::cout,NBT);

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
   time_ABA[4] = timer.toc()/NBT; // ABA derivatives timing

  std::cout << "ABA derivatives= \t\t" << time_ABA[4] << endl;
  // std::cout << "ABA derivatives= \t\t"; timer.toc(std::cout,NBT);
  
//-------------------------------------------------------//
//----- Forward Dynamics partials using the AZA method---//
//-------------------------------------------------------//
//IPR using AZA method here------------------------------//

  // PINOCCHIO_EIGEN_PLAIN_ROW_MAJOR_TYPE(MatrixXd) daba_dq_2(MatrixXd::Zero(model.nv,model.nv));
  // PINOCCHIO_EIGEN_PLAIN_ROW_MAJOR_TYPE(MatrixXd) daba_dv_2(MatrixXd::Zero(model.nv,model.nv));

// timer.tic();

// for (int ii=0; ii < model.nv ;ii++){

//   taus[0] = -drnea_dq.col(ii);

//  SMOOTH(NBT)
//  aza(model,data,qs[_smooth],qdots[_smooth],taus[_smooth]);

//  daba_dq_2.col(ii) = data.ddq_new;

//   // cout << "ID partial original is " << drnea_dq.col(ii) << endl;
//   // cout << "\n" << endl;
//   // cout << "FD partial using AZA method is " << daba_dq_2.col(ii) << endl;
//   // cout << "\n" << endl;
//   // cout << "FD partial using original method is " << daba_dq.col(ii) << endl;
//   // cout << "-------------------------" << endl;

// }

// std::cout << "IPR using AZA method is = \t\t"; timer.toc(std::cout,NBT);

//-------------------------------------------------------//
// ---------- IPR using the regular DMM method-----------//
//IPR using AZA method here------------------------------//
//-------------------------------------------------------//
// timer.tic();

// std::cout << "IPR using DMM method is = \t\t"; timer.toc(std::cout,NBT);



//  // Difference matrix calculations here
// MatrixXd diff_daba_dq(MatrixXd::Zero(model.nv,model.nv));

// diff_daba_dq = daba_dq_2-daba_dq;

// std::cout << "\n " << endl;
// //std::cout << "Difference matrix of FD partials with AZA vs originial FD partials is" << diff_daba_dq << endl;
// std::cout << "Norm of the difference matrix of FD partials with AZA FD partials is " << diff_daba_dq.squaredNorm() << std::endl;
// std::cout << "---------------------------------------------" << endl;

//-------------------------------------------------//
// running Minv using ABA (AZA) -------------------//
//-------------------------------------------------//

// Zero inputs
  
// PINOCCHIO_EIGEN_PLAIN_ROW_MAJOR_TYPE(MatrixXd) Minv_aza(MatrixXd::Zero(model.nv,model.nv));
// MatrixXd eq(MatrixXd::Zero(model.nv,model.nv));
// //MatrixXd eq2(MatrixXd::Zero(model.nv,model.nv));

// // Inputs for AZA

// qdots[0].setZero();

// timer.tic();
// for (int ii=0; ii < model.nv ;ii++){

//   taus[0] = tau_mat[0].col(ii);

//  SMOOTH(NBT)
//  aza(model,data,qs[_smooth],qdots[_smooth],taus[_smooth]);
  
//    Minv_aza.row(ii) = data.ddq_new;

// }

// std::cout << "Minv using AZA timing is = \t\t"; timer.toc(std::cout,NBT);
  // eq = Minv_aza - data.Minv;

//-----------------------------------------------------------------//
// FD partials using AZAmat function here--------------------------//
//-----------------------------------------------------------------//

   cout <<"\n " << endl;
   cout <<"--------------------------------------------------------" << endl;
 //  cout << "ID partial wrt q using original method is " << drnea_dq << endl;
 //  cout << "ID partial wrt qd using original method is " << drnea_dv << endl;
 // cout << "FD partial wrt q using original method is " << daba_dq << endl;
 // cout << "FD partial wrt qd using original method is " << daba_dv << endl;
   cout <<"--------------------------------------------------------" << endl;

  tau_mat_n2n[0] << -drnea_dq,-drnea_dv; // concatenating partial wrt q and qdot

 timer.tic();
  //cout << "concatenated matrix is" << tau_mat_n2n[0] << endl;
  SMOOTH(NBT)
  azamat(model,data,qs[_smooth],qdots[_smooth],tau_mat_n2n[_smooth]);
  time_ABA[5] = timer.toc()/NBT;

  std::cout << "IPR using AZA_mat method is = \t\t" << time_ABA[5] << endl;

  //std::cout << "IPR using AZA_mat method is = \t\t"; timer.toc(std::cout,NBT);
 
 //-----------------------------------------------------------------//
    // FD partials using AZAmat_v3 function here-----------------------//
    //-----------------------------------------------------------------//

    tau_mat_n2n_v3[0] << -drnea_dq,-drnea_dv; // concatenating partial wrt q and qdot

    timer.tic();
    SMOOTH(NBT)
    azamat_v3(model,data,qs[_smooth],tau_mat_n2n_v3[_smooth]);
    time_ABA[6] = timer.toc()/NBT;
    std::cout << "IPR using AZA_mat_v3 method is = \t\t" << time_ABA[6] << endl;

 // Difference matrix calculations here
  
  MatrixXd diff_daba_dq2(MatrixXd::Zero(model.nv,model.nv));
  MatrixXd diff_daba_dqd2(MatrixXd::Zero(model.nv,model.nv));
  MatrixXd diff_mat1(MatrixXd::Zero(model.nv,2*model.nv));

  //eq2 = data.Minv_mat_prod - data.Minv;

  diff_daba_dq2 = daba_dq-data.Minv_mat_prod.middleCols(0,model.nv);
  diff_daba_dqd2 = daba_dv-data.Minv_mat_prod.middleCols(model.nv,model.nv);

  diff_mat1 = data.Minv_mat_prod - data.Minv_mat_prod_v3;

  std::cout << "---------------------------------------------" << endl;

  std::cout << "Norm of the difference matrix for AZAmat FD partial wrt q from orig FD partial wrt q is " << diff_daba_dq2.squaredNorm() << std::endl;
  std::cout << "Norm of the difference matrix for AZAmat FD partial wrt qd from orig FD partial wrt qd is " << diff_daba_dqd2.squaredNorm() << std::endl;

  std::cout << "\n" << endl;

  std::cout << "Norm of difference between mat_v1 and mat_v3 is" << diff_mat1.squaredNorm() << std::endl;

  //-------------------------------------------------------------
  // Just the IPR equation here using DMM
  //-------------------------------------------------------------



  //---------------------------------------------------------------

  // Writing all the timings to the file
  for (int ii=0; ii<8 ; ii++)
  {
    file1 << time_ABA[ii] << "\n" << endl;
  }
    file1.close();

  // std::cout<< "Minv original is " << data.Minv << endl;
  // std::cout << "---------------------------------------------" << endl;
  // std::cout << "Minv using AZA value is " << Minv_aza << endl;
  // std::cout << "---------------------------------------------" << endl;
  // std::cout << "Minv using AZA_mat value is " << data.Minv_mat_prod << endl;
  // std::cout << "---------------------------------------------" << endl;

  // Testing some eigen matrices here

//   MatrixXd temp1(MatrixXd::Random(5,5));

//  cout <<"temp1 matrix is" << temp1 << endl;
//  cout << "\n " << endl;

//  cout << "temp1 fourth row is" << temp1.row(3) << endl;
//  cout << "temp1 fourth row is" << temp1.middleRows(3,1) << endl;



  return 0;

}