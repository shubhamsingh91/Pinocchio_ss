
// checking some boost options, ABA forward dynamics etc.
// This script is only for testing the timings for serial chains

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
   int NBT = 10000; // 50000 initially
  #else
     int NBT = 1;
    std::cout << "(the time score in debug mode is not relevant) " << std::endl;
  #endif

 //   int n_vec[]={2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,25,30,35,40,45,50,60,70,80,90,100,120,150,180,200,220,250,280,300,320,350,380,400,420,450,480,500};
   // int n_vec[]={220,250,280,300,320,350,380,400,420,450,480,500};
    int n_vec[]={320,350,380,400,420,450,480,500};

    // for (int jj=0; jj<19; jj++)
    // {
    //     n_vec[jj]=jj+2;
    //   // cout << "n_vec is " << n_vec[jj] << endl;
    // }



  for(int jj=0; jj<46; jj++){


        Model model;

        string str_urdf="", str_file_ext ,n_str;
        string robot_name="";

        int n_links;

         //cout << "Enter the n here" << endl;
        // cin >> n_links;

        n_links = n_vec[jj];

        if(n_links>100)
        {
            NBT = 1000;
        }
        if (n_links>200)
        {
            NBT = 100;
        }
         if(n_links>300)
        {
            NBT = 1;
        }

        cout << "n = " << n_links << endl;

        n_str = to_string(n_links);

        robot_name.append(n_str);
        robot_name.append("link");

    // cout <<"robot_name is " << robot_name << std::endl;

        str_file_ext = "/home/ss86299/Desktop/test_pinocchio/pinocchio/models/";
        str_urdf.append(str_file_ext);
        str_urdf.append(robot_name);
        str_urdf.append(".urdf"); 
    
        cout <<"str_urdf is "<< str_urdf << endl;

    //  str_urdf= "/home/ss86299/Desktop/test_pinocchio/pinocchio/models/500link.urdf"; //2link

        std :: string filename = str_urdf;

        // Opening a file to write to it
        
        string filewrite="";

        ofstream file1;

        filewrite.append(robot_name);
        filewrite.append(".txt");

        file1.open(filewrite);


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

        cout << "with ff = " << with_ff << "\n \n";    
        std::cout << "nq = " << model.nq << std::endl;
        std::cout << "nv = " << model.nv << std::endl;
        cout << "Model is " << robot_name << endl;

        Data data(model);
        VectorXd qmax = Eigen::VectorXd::Ones(model.nq);

        PINOCCHIO_ALIGNED_STD_VECTOR(VectorXd) qs     (NBT);
        PINOCCHIO_ALIGNED_STD_VECTOR(VectorXd) qdots  (NBT);
        PINOCCHIO_ALIGNED_STD_VECTOR(VectorXd) taus (NBT); 
        PINOCCHIO_ALIGNED_STD_VECTOR(VectorXd) qddots (NBT); 
        PINOCCHIO_ALIGNED_STD_VECTOR(MatrixXd) tau_mat (NBT); 
        PINOCCHIO_ALIGNED_STD_VECTOR(MatrixXd) tau_mat_n2n (NBT); 
    
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
    
        }   

    double time_ABA[6];

    timer.tic();
    SMOOTH(NBT)
    aba(model,data,qs[_smooth],qdots[_smooth],taus[_smooth]);
    time_ABA[0] = timer.toc()/NBT; // ABA timing

   // std::cout << "ABA= \t\t"; timer.toc(std::cout,NBT);
    std::cout << "ABA= \t\t" <<  time_ABA[0] << endl;


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
    
    //-----------------------------------------------------------------//
    // FD partials using AZAmat function here--------------------------//
    //-----------------------------------------------------------------//

    //  cout << "ID partial wrt q using original method is " << drnea_dq << endl;
    //  cout << "ID partial wrt qd using original method is " << drnea_dv << endl;
    // cout << "FD partial wrt q using original method is " << daba_dq << endl;
    // cout << "FD partial wrt qd using original method is " << daba_dv << endl;

    tau_mat_n2n[0] << -drnea_dq,-drnea_dv; // concatenating partial wrt q and qdot

    timer.tic();
    SMOOTH(NBT)
    azamat(model,data,qs[_smooth],qdots[_smooth],tau_mat_n2n[_smooth]);
    time_ABA[5] = timer.toc()/NBT;

    std::cout << "IPR using AZA_mat method is = \t\t" << time_ABA[5] << endl;
    
    std::cout <<"---------------------------------------" << endl;

    // Difference matrix calculations here
    
    MatrixXd diff_daba_dq2(MatrixXd::Zero(model.nv,model.nv));
    MatrixXd diff_daba_dqd2(MatrixXd::Zero(model.nv,model.nv));

    //eq2 = data.Minv_mat_prod - data.Minv;

    diff_daba_dq2 = daba_dq-data.Minv_mat_prod.middleCols(0,model.nv);
    diff_daba_dqd2 = daba_dv-data.Minv_mat_prod.middleCols(model.nv,model.nv);


    std::cout << "Norm of the difference matrix for AZAmat FD partial wrt q from orig FD partial wrt q is " << diff_daba_dq2.squaredNorm() << std::endl;
    std::cout << "Norm of the difference matrix for AZAmat FD partial wrt qd from orig FD partial wrt qd is " << diff_daba_dqd2.squaredNorm() << std::endl;

    std::cout << "\n" << endl;

    // Writing all the timings to the file
    for (int ii=0; ii<6 ; ii++)
    {
        file1 << time_ABA[ii] << "\n" << endl;
    }
        file1.close();

  }

  return 0;

}