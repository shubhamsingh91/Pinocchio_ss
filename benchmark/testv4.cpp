
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

   str_int = 2; // integer to change the str_urdf

   // cout << "Enter the str_int here" << endl;
   // cin >> str_int;

    str_urdf[0] = "/home/ss86299/Desktop/test_pinocchio/pinocchio/models/double_pendulum_v1.urdf"; // double pendulum
    str_urdf[1] = "/home/ss86299/Desktop/test_pinocchio/pinocchio/models/40link.urdf"; // 40_link
    str_urdf[2] = "/home/ss86299/Desktop/test_pinocchio/pinocchio/models/simple_humanoid.urdf"; //simple humanoid

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

  // MatrixXd tau_mat(MatrixXd::Identity(model.nv,model.nv));

  // randomizing input data here

    std:: cout << "NBT variable is" << NBT << endl;

    for(size_t i=0;i<NBT;++i)
    {
      qs[i]     = randomConfiguration(model,-qmax,qmax);
      qdots[i]  = Eigen::VectorXd::Random(model.nv);
      taus[i] = Eigen::VectorXd::Random(model.nv);
      qddots[i] =  Eigen::VectorXd::Random(model.nv);
      tau_mat_n2n[i] = Eigen::MatrixXd::Random(model.nv,2*model.nv);
      tau_mat_n2n_v2[i] = tau_mat_n2n[i];

    }   

  std::cout << "tau_mat input here is" << tau_mat_n2n[0] << std::endl;


 timer.tic();
  SMOOTH(NBT)
   azamat(model,data,qs[_smooth],qdots[_smooth],tau_mat_n2n[_smooth]); // this one is correct
  std::cout << "IPR using AZA_mat method is = \t\t"; timer.toc(std::cout,NBT);
  std::cout <<"------------------------------" << std::endl;
  std::cout << "Mat v1 is" << data.Minv_mat_prod << std::endl;


//  timer.tic();
//   SMOOTH(NBT)
//    azamat_v2(model,data,qs[_smooth],qdots[_smooth],tau_mat_n2n_v2[_smooth]); // this is not correct
//   std::cout << "IPR using AZA_mat_v2 method is = \t\t"; timer.toc(std::cout,NBT);

// Method-3

  // timer.tic();
  // SMOOTH(NBT)
  //   azamat_v3(model,data,qs[_smooth],tau_mat_n2n_v2[_smooth]);
  // std::cout << "IPR using AZA_mat_v3 method is = \t\t"; timer.toc(std::cout,NBT);

// Method-4

  timer.tic();
  SMOOTH(NBT)
    azamat_v4(model,data,qs[_smooth],tau_mat_n2n_v2[_smooth]);
  std::cout << "IPR using AZA_mat_v4 method is = \t\t"; timer.toc(std::cout,NBT);


// Comparison of results here

  MatrixXd diff_mat(MatrixXd::Zero(model.nv,2*model.nv));

  diff_mat = data.Minv_mat_prod - data.Minv_mat_prod_v3;

  std::cout <<"------------------------------" << std::endl;
  std::cout << "\n " << std::endl; 
 // std::cout << "Mat v1 is" << data.Minv_mat_prod << std::endl;
  std::cout <<"------------------------------" << std::endl;
  std::cout << "\n " << std::endl; 
  //std::cout << "Mat v3 is" << data.Minv_mat_prod_v3 << std::endl;
  std::cout <<"------------------------------" << std::endl;
 // std::cout << "Difference matrix here is" << diff_mat << std::endl;

  std::cout << "Norm of the difference matrix is " << diff_mat.squaredNorm() << std::endl;

//   MatrixXd temp1(MatrixXd::Random(5,5));

//  cout <<"temp1 matrix is" << temp1 << endl;
//  cout << "\n " << endl;

//  cout << "temp1 second row is" << temp1.row(2) << endl;
//  cout << "temp1 second row is" << temp1.middleRows(2,1) << endl;




  return 0;

}