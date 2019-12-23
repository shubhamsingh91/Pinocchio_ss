//
// Copyright (c) 2016-2019 CNRS INRIA
//

#include "pinocchio/multibody/data.hpp"
#include "pinocchio/multibody/model.hpp"

#include "pinocchio/algorithm/check.hpp"
#include "pinocchio/algorithm/model.hpp"
#include "pinocchio/algorithm/kinematics.hpp"
#include "pinocchio/algorithm/frames.hpp"
#include "pinocchio/algorithm/joint-configuration.hpp"

#include "pinocchio/parsers/sample-models.hpp"

#include <boost/test/unit_test.hpp>
#include <boost/utility/binary.hpp>

using namespace pinocchio;

BOOST_AUTO_TEST_SUITE ( BOOST_TEST_MODULE )

  BOOST_AUTO_TEST_CASE(test_model_subtree)
  {
    Model model;
    buildModels::humanoidRandom(model);
    
    Model::JointIndex idx_larm1 = model.getJointId("larm1_joint");
    BOOST_CHECK(idx_larm1<(Model::JointIndex)model.njoints);
    Model::IndexVector subtree = model.subtrees[idx_larm1];
    BOOST_CHECK(subtree.size()==6);
    
    for(size_t i=1; i<subtree.size();++i)
    {
      BOOST_CHECK(model.parents[subtree[i]]==subtree[i-1]);
    }
      
    // Check that i starts subtree[i]
    for(JointIndex joint_id = 1; joint_id < (JointIndex)model.njoints; ++joint_id)
    {
      BOOST_CHECK(model.subtrees[joint_id][0] == joint_id);
    }
    
    // Check that subtree[0] contains all joint ids
    for(JointIndex joint_id = 1; joint_id < (JointIndex)model.njoints; ++joint_id)
    {
      BOOST_CHECK(model.subtrees[0][joint_id-1] == joint_id);
    }
  }

  BOOST_AUTO_TEST_CASE(test_model_support)
  {
    Model model;
    buildModels::humanoidRandom(model);
    const Model::IndexVector support0_ref(1,0);
    BOOST_CHECK(model.supports[0] == support0_ref);

    // Check that i ends supports[i]
    for(JointIndex joint_id = 1; joint_id < (JointIndex)model.njoints; ++joint_id)
    {
      BOOST_CHECK(model.supports[joint_id].back() == joint_id);
      Model::IndexVector & support = model.supports[joint_id];
      
      size_t current_id = support.size()-2;
      for(JointIndex parent_id = model.parents[joint_id];
          parent_id > 0;
          parent_id = model.parents[parent_id],
          current_id--)
      {
        BOOST_CHECK(parent_id == support[current_id]);
      }
    }
  }

  BOOST_AUTO_TEST_CASE(test_model_subspace_dimensions)
  {
    Model model;
    buildModels::humanoidRandom(model);

    // Check that i ends supports[i]
    for(JointIndex joint_id = 1; joint_id < (JointIndex)model.njoints; ++joint_id)
    {
      const Model::JointModel & jmodel = model.joints[joint_id];
      
      BOOST_CHECK(model.nqs[joint_id] == jmodel.nq());
      BOOST_CHECK(model.idx_qs[joint_id] == jmodel.idx_q());
      
      BOOST_CHECK(model.nvs[joint_id] == jmodel.nv());
      BOOST_CHECK(model.idx_vs[joint_id] == jmodel.idx_v());
    }
  }

  BOOST_AUTO_TEST_CASE(comparison)
  {
    Model model;
    buildModels::humanoidRandom(model);
    
    BOOST_CHECK(model == model);
  }
  
  BOOST_AUTO_TEST_CASE(cast)
  {
    Model model;
    buildModels::humanoidRandom(model);
    
    BOOST_CHECK(model.cast<double>() == model);
    BOOST_CHECK(model.cast<long double>().cast<double>() == model);
  }

#ifdef PINOCCHIO_WITH_HPP_FCL
  struct AddPrefix {
    std::string p;
    std::string operator() (const std::string& n) { return p + n; }
    Frame operator() (const Frame& _f) { Frame f (_f); f.name = p + f.name; return f; }
    AddPrefix (const char* _p) : p(_p) {}
  };

  BOOST_AUTO_TEST_CASE(append)
  {
    Model manipulator, humanoid;
    GeometryModel geomManipulator, geomHumanoid;

    buildModels::manipulator(manipulator);
    buildModels::manipulatorGeometries(manipulator, geomManipulator);
    // Add prefix to joint and frame names
    AddPrefix addManipulatorPrefix ("manipulator/");
    std::transform (++manipulator.names.begin(), manipulator.names.end(),
        ++manipulator.names.begin(), addManipulatorPrefix);
    std::transform (++manipulator.frames.begin(), manipulator.frames.end(),
        ++manipulator.frames.begin(), addManipulatorPrefix);

    BOOST_TEST_MESSAGE(manipulator);

    buildModels::humanoid(humanoid);
    buildModels::humanoidGeometries(humanoid, geomHumanoid);
    // Add prefix to joint and frame names
    AddPrefix addHumanoidPrefix ("humanoid/");
    std::transform (++humanoid.names.begin(), humanoid.names.end(),
        ++humanoid.names.begin(), addHumanoidPrefix);
    std::transform (++humanoid.frames.begin(), humanoid.frames.end(),
        ++humanoid.frames.begin(), addHumanoidPrefix);

    BOOST_TEST_MESSAGE(humanoid);
    
    //TODO fix inertia of the base
    manipulator.inertias[0].setRandom();
    SE3 aMb = SE3::Random();

    // First append a model to the universe frame.
    Model model1;
    GeometryModel geomModel1;
    int fid = 0;
    appendModel (humanoid, manipulator, geomHumanoid, geomManipulator, (FrameIndex)fid,
        SE3::Identity(), model1, geomModel1);

    Data data1 (model1);
    BOOST_CHECK(model1.check(data1));

    BOOST_TEST_MESSAGE(model1);

    // Second, append a model to a moving frame.
    Model model;
    
    GeometryModel geomModel;
    fid = humanoid.addFrame (Frame ("humanoid/add_manipulator", 
          humanoid.getJointId("humanoid/chest2_joint"),
          humanoid.getFrameId("humanoid/chest2_joint"), aMb,
          OP_FRAME));

    BOOST_CHECK(fid >= 0);

    appendModel (humanoid, manipulator, geomHumanoid, geomManipulator, (FrameIndex)fid,
        SE3::Identity(), model, geomModel);

    BOOST_TEST_MESSAGE(model);

    Data data (model);
    BOOST_CHECK(model.check(data));

    // Check the model
    BOOST_CHECK_EQUAL(model.getJointId("humanoid/chest2_joint"),
        model.parents[model.getJointId("manipulator/shoulder1_joint")]);

    // check the joint order and the inertias
    JointIndex chest2 = model.getJointId("humanoid/chest2_joint");
    for (JointIndex jid = 1; jid < chest2; ++jid) {
      BOOST_TEST_MESSAGE("Checking joint " << jid << " " << model.names[jid]);
      BOOST_CHECK_EQUAL(model.names[jid], humanoid.names[jid]);
      BOOST_CHECK_EQUAL(model.inertias[jid], humanoid.inertias[jid]);
      BOOST_CHECK_EQUAL(model.jointPlacements[jid], humanoid.jointPlacements[jid]);
    }
    BOOST_TEST_MESSAGE("Checking joint " << chest2 << " " << model.names[chest2]);
    BOOST_CHECK_EQUAL(model.names[chest2], humanoid.names[chest2]);
    BOOST_CHECK_MESSAGE(model.inertias[chest2].isApprox(manipulator.inertias[0].se3Action(aMb) + humanoid.inertias[chest2]),
        model.inertias[chest2] << " != " << manipulator.inertias[0].se3Action(aMb) + humanoid.inertias[chest2]);
    BOOST_CHECK_EQUAL(model.jointPlacements[chest2], humanoid.jointPlacements[chest2]);

    for (JointIndex jid = 1; jid < manipulator.joints.size(); ++jid) {
      BOOST_TEST_MESSAGE("Checking joint " << chest2+jid << " " << model.names[chest2+jid]);
      BOOST_CHECK_EQUAL(model.names[chest2+jid], manipulator.names[jid]);
      BOOST_CHECK_EQUAL(model.inertias[chest2+jid], manipulator.inertias[jid]);
      if (jid==1)
        BOOST_CHECK_EQUAL(model.jointPlacements[chest2+jid], aMb*manipulator.jointPlacements[jid]);
      else
        BOOST_CHECK_EQUAL(model.jointPlacements[chest2+jid], manipulator.jointPlacements[jid]);
    }
    for (JointIndex jid = chest2+1; jid < humanoid.joints.size(); ++jid) {
      BOOST_TEST_MESSAGE("Checking joint " << jid+manipulator.joints.size()-1 << " " << model.names[jid+manipulator.joints.size()-1]);
      BOOST_CHECK_EQUAL(model.names[jid+manipulator.joints.size()-1], humanoid.names[jid]);
      BOOST_CHECK_EQUAL(model.inertias[jid+manipulator.joints.size()-1], humanoid.inertias[jid]);
      BOOST_CHECK_EQUAL(model.jointPlacements[jid+manipulator.joints.size()-1], humanoid.jointPlacements[jid]);
    }

    // Check the frames
    for (FrameIndex fid = 1; fid < humanoid.frames.size(); ++fid) {
      const Frame& frame  = humanoid.frames[fid],
                   parent = humanoid.frames[frame.previousFrame];
      BOOST_CHECK(model.existFrame (frame.name, frame.type));
      const Frame& nframe  = model.frames[model.getFrameId(frame.name, frame.type)],
                   nparent = model.frames[nframe.previousFrame];
      BOOST_CHECK_EQUAL(parent.name, nparent.name);
      BOOST_CHECK_EQUAL(frame.placement, nframe.placement);
    }
    for (FrameIndex fid = 1; fid < manipulator.frames.size(); ++fid) {
      const Frame& frame  = manipulator.frames[fid],
                   parent = manipulator.frames[frame.previousFrame];
      BOOST_CHECK(model.existFrame (frame.name, frame.type));
      const Frame& nframe  = model.frames[model.getFrameId(frame.name, frame.type)],
                   nparent = model.frames[nframe.previousFrame];
      if (frame.previousFrame > 0) {
        BOOST_CHECK_EQUAL(parent.name, nparent.name);
        BOOST_CHECK_EQUAL(frame.placement, nframe.placement);
      }
    }
  }
#endif

BOOST_AUTO_TEST_CASE(test_buildReducedModel)
{
  Model humanoid_model;
  buildModels::humanoid(humanoid_model);
  
  static const std::vector<JointIndex> empty_index_vector;
  
  humanoid_model.lowerPositionLimit.head<3>().fill(-1.);
  humanoid_model.upperPositionLimit.head<3>().fill( 1.);
  
  Eigen::VectorXd reference_config_humanoid = randomConfiguration(humanoid_model);
  Model humanoid_copy_model = buildReducedModel(humanoid_model,empty_index_vector,reference_config_humanoid);
  
  BOOST_CHECK(humanoid_copy_model.names == humanoid_model.names);
  BOOST_CHECK(humanoid_copy_model.joints == humanoid_model.joints);
  BOOST_CHECK(humanoid_copy_model == humanoid_model);
  
  std::vector<JointIndex> joints_to_lock;
  const std::string joint1_to_lock = "rarm_shoulder2_joint";
  joints_to_lock.push_back(humanoid_model.getJointId(joint1_to_lock));
  const std::string joint2_to_lock = "larm_shoulder2_joint";
  joints_to_lock.push_back(humanoid_model.getJointId(joint2_to_lock));

  Model reduced_humanoid_model = buildReducedModel(humanoid_model,joints_to_lock,reference_config_humanoid);
  BOOST_CHECK(reduced_humanoid_model.njoints == humanoid_model.njoints-(int)joints_to_lock.size());
  BOOST_CHECK(reduced_humanoid_model != humanoid_model);
  
  Model reduced_humanoid_model_other_signature;
  buildReducedModel(humanoid_model,joints_to_lock,reference_config_humanoid,reduced_humanoid_model_other_signature);
  BOOST_CHECK(reduced_humanoid_model == reduced_humanoid_model_other_signature);
  
  // Check that forward kinematics give same results
  Data data(humanoid_model), reduced_data(reduced_humanoid_model);
  Eigen::VectorXd q = reference_config_humanoid;
  Eigen::VectorXd reduced_q(reduced_humanoid_model.nq);
  
  for(JointIndex joint_id = 1; joint_id < (JointIndex)reduced_humanoid_model.njoints; ++joint_id)
  {
    const JointIndex reference_joint_id = humanoid_model.getJointId(reduced_humanoid_model.names[joint_id]);
    
    reduced_humanoid_model.joints[joint_id].jointConfigSelector(reduced_q) =
    humanoid_model.joints[reference_joint_id].jointConfigSelector(q);
  }
  
  framesForwardKinematics(humanoid_model,data,q);
  framesForwardKinematics(reduced_humanoid_model,reduced_data,reduced_q);
  
  for(size_t frame_id = 0; frame_id < reduced_humanoid_model.frames.size(); ++frame_id)
  {
    const Frame & reduced_frame = reduced_humanoid_model.frames[frame_id];
    switch(reduced_frame.type)
    {
      case JOINT:
      case FIXED_JOINT:
      {
        // May not be present in the original model
        if(humanoid_model.existJointName(reduced_frame.name))
        {
          const JointIndex joint_id = humanoid_model.getJointId(reduced_frame.name);
          BOOST_CHECK(data.oMi[joint_id].isApprox(reduced_data.oMf[frame_id]));
        }
        else if(humanoid_model.existFrame(reduced_frame.name))
        {
          const FrameIndex humanoid_frame_id = humanoid_model.getFrameId(reduced_frame.name);
          BOOST_CHECK(data.oMf[humanoid_frame_id].isApprox(reduced_data.oMf[frame_id]));
        }
        else
        {
          BOOST_CHECK_MESSAGE(false, "The frame " << reduced_frame.name << " is not presend in the humanoid_model");
        }
        break;
      }
      default:
      {
        BOOST_CHECK(humanoid_model.existFrame(reduced_frame.name));
        const FrameIndex humanoid_frame_id = humanoid_model.getFrameId(reduced_frame.name);
        BOOST_CHECK(data.oMf[humanoid_frame_id].isApprox(reduced_data.oMf[frame_id]));
        break;
      }
        
    }
  }
}

BOOST_AUTO_TEST_SUITE_END()
