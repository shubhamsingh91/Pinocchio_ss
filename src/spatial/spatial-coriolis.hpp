//
// Copyright (c) 2015-2019 CNRS INRIA
// Copyright (c) 2016 Wandercraft, 86 rue de Paris 91400 Orsay, France.
//

#ifndef __pinocchio_coriolis_hpp__
#define __pinocchio_coriolis_hpp__

#include <iostream>

#include "pinocchio/math/fwd.hpp"
#include "pinocchio/spatial/symmetric3.hpp"
#include "pinocchio/spatial/force.hpp"
#include "pinocchio/spatial/motion.hpp"
#include "pinocchio/spatial/skew.hpp"

namespace pinocchio
{

  template< class Derived>
  class CoriolisBase
  {
  protected:

    typedef Derived  Derived_t;
    SPATIAL_TYPEDEF_TEMPLATE(Derived_t);

  public:
    Derived_t & derived() { return *static_cast<Derived_t*>(this); }
    const Derived_t & derived() const { return *static_cast<const Derived_t*>(this); }

    const Vector3 &    lever()   const { return static_cast<const Derived_t*>(this)->lever(); }
    Vector3 &          lever() { return static_cast<const Derived_t*>(this)->lever(); }

    const Matrix3 & inertia() const { return static_cast<const Derived_t*>(this)->inertia(); }
    Matrix3 &       inertia() { return static_cast<const Derived_t*>(this)->inertia(); }

    Matrix6 matrix() const { return derived().matrix_impl(); }
    operator Matrix6 () const { return matrix(); }

    Derived_t& operator= (const Derived_t& clone){return derived().__equl__(clone);}
    bool operator==(const Derived_t & other) const {return derived().isEqual(other);}
    bool operator!=(const Derived_t & other) const { return !(*this == other); }
    
    Derived_t& operator+= (const Derived_t & Yb) { return derived().__pequ__(Yb); }
    Derived_t operator+(const Derived_t & Yb) const { return derived().__plus__(Yb); }
    
    template<typename MotionDerived>
    ForceTpl<typename traits<MotionDerived>::Scalar,traits<MotionDerived>::Options>
    operator*(const MotionDense<MotionDerived> & v) const
    { return derived().__mult__(v); }

    void setZero() { derived().setZero(); }
    void setRandom() { derived().setRandom(); }
    
    bool isApprox(const Derived & other, const Scalar & prec = Eigen::NumTraits<Scalar>::dummy_precision()) const
    { return derived().isApprox_impl(other, prec); }
    
    bool isZero(const Scalar & prec = Eigen::NumTraits<Scalar>::dummy_precision()) const
    { return derived().isZero_impl(prec); }

    /// aI = aXb.act(bI)
    //Derived_t se3Action(const SE3 & M) const { return derived().se3Action_impl(M); }

    /// bI = aXb.actInv(aI)
    //Derived_t se3ActionInverse(const SE3 & M) const { return derived().se3ActionInverse_impl(M); }

    void disp(std::ostream & os) const { static_cast<const Derived_t*>(this)->disp_impl(os); }
    friend std::ostream & operator << (std::ostream & os,const InertiaBase<Derived_t> & X)
    { 
      X.disp(os);
      return os;
    }

  }; // class InertiaBase


  template<typename T, int U>
  struct traits< CoriolisTpl<T, U> >
  {
    typedef T Scalar;
    typedef Eigen::Matrix<T,3,1,U> Vector3;
    typedef Eigen::Matrix<T,4,1,U> Vector4;
    typedef Eigen::Matrix<T,6,1,U> Vector6;
    typedef Eigen::Matrix<T,3,3,U> Matrix3;
    typedef Eigen::Matrix<T,4,4,U> Matrix4;
    typedef Eigen::Matrix<T,6,6,U> Matrix6;
    typedef Matrix6 ActionMatrix_t;
    typedef Vector3 Angular_t;
    typedef Vector3 Linear_t;
    typedef const Vector3 ConstAngular_t;
    typedef const Vector3 ConstLinear_t;
    typedef Eigen::Quaternion<T,U> Quaternion_t;
    typedef SE3Tpl<T,U> SE3;
    typedef ForceTpl<T,U> Force;
    typedef MotionTpl<T,U> Motion;
    typedef Symmetric3Tpl<T,U> Symmetric3;
    enum {
      LINEAR = 0,
      ANGULAR = 3
    };
  }; // traits CoriolisTpl

  template<typename _Scalar, int _Options>
  class CoriolisTpl : public CoriolisBase< CoriolisTpl< _Scalar, _Options > >
  {
  public:
    friend class CoriolisBase< CoriolisTpl< _Scalar, _Options > >;
    SPATIAL_TYPEDEF_TEMPLATE(CoriolisTpl);
    enum { Options = _Options };
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    
    typedef typename Eigen::Matrix<_Scalar, 10, 1, _Options> Vector10;

  public:
    // Constructors
    CoriolisTpl()
    {}

    CoriolisTpl(const Vector3 & com, const Matrix3 & rotational_inertia)
    : m_com(com), m_inertia(rotational_inertia)
    {}
    
    CoriolisTpl(const Matrix6 & I6)
    {
      const Matrix3 & mc_cross = I6.template block <3,3>(LINEAR,ANGULAR);
      lever() = unSkew(mc_cross);
      inertia() = I6.template block <3,3>(ANGULAR,ANGULAR);
    }


    CoriolisTpl(const Inertia & I, const Motion & v)
    {
      const Matrix3 & Ic = I.inertia();
      Matrix3 tmp;

      cross( v.angular() , Ic, tmp);
      const Vector3 hw = Ic*v.angular();
      inertia() = tmp + tmp.transpose();

      addSkew(-hw , inertia());

      lever()= -2*I.mass()*(v.linear() + v.angular().cross(I.lever()));
      skewSquare(I.lever(), lever() , tmp );
      inertia() += tmp;

      //inertia() = Bang;
      //CoriolisTpl(mv2, Bang);
    }
    
    CoriolisTpl(const CoriolisTpl & clone)  // Copy constructor
    : m_com(clone.lever()), m_inertia(clone.inertia())
    {}

    CoriolisTpl& operator=(const CoriolisTpl & clone)  // Copy assignment operator
    {
      m_com = clone.lever();
      m_inertia = clone.inertia();
      return *this;
    }

    template<int O2>
    CoriolisTpl(const CoriolisTpl<Scalar,O2> & clone)
    : m_com(clone.lever())
    , m_inertia(clone.inertia().matrix())
    {}

    // Initializers
    static CoriolisTpl Zero() 
    {
      return CoriolisTpl(Vector3::Zero(), 
                        Matrix3::Zero());
    }
    
    void setZero() { lever().setZero(); inertia().setZero(); }
    

    static CoriolisTpl Random()
    {
        // We have to shoot "I" definite positive and not only symmetric.
      return CoriolisTpl(Vector3::Random(),
                        Matrix3::Random());
    }
    

    void setRandom()
    {
      lever().setRandom(); inertia().setRandom();
    }

    Matrix6 matrix_impl() const
    {
      Matrix6 M;
      
      M.template block<6,3>(0, LINEAR ).setZero();
      M.template block<3,3>(LINEAR,ANGULAR ) = skew(lever());
      M.template block<3,3>(ANGULAR, ANGULAR) = inertia();
      return M;
    }



    // Arithmetic operators
    CoriolisTpl & __equl__(const CoriolisTpl & clone)
    {
      lever() = clone.lever(); inertia() = clone.inertia();
      return *this;
    }

    // Required by std::vector boost::python bindings.
    bool isEqual( const CoriolisTpl& Y2 ) const
    { 
      return (lever()==Y2.lever()) && (inertia()==Y2.inertia());
    }
    
    bool isApprox_impl(const CoriolisTpl & other,
                       const Scalar & prec = Eigen::NumTraits<Scalar>::dummy_precision()) const
    {
      using math::fabs;
      return lever().isApprox(other.lever(),prec)
      && inertia().isApprox(other.inertia(),prec);
    }
    
    bool isZero_impl(const Scalar & prec = Eigen::NumTraits<Scalar>::dummy_precision()) const
    {
      using math::fabs;
      return lever().isZero(prec)
          && inertia().isZero(prec);
    }
    
    CoriolisTpl __plus__(const CoriolisTpl & Yb) const
    {
      return CoriolisTpl(lever() + Yb.lever(), inertia()+Yb.inertia());
    }

    CoriolisTpl& __pequ__(const CoriolisTpl & Yb)
    {
      lever() += Yb.lever();
      inertia() += Yb.inertia();
      return *this;
    }

    template<typename MotionDerived>
    ForceTpl<typename traits<MotionDerived>::Scalar,traits<MotionDerived>::Options>
    __mult__(const MotionDense<MotionDerived> & v) const
    {
      typedef ForceTpl<typename traits<MotionDerived>::Scalar,traits<MotionDerived>::Options> ReturnType;
      ReturnType f;
      __mult__(v,f);
      return f;
    }
    
    template<typename MotionDerived, typename ForceDerived>
    void __mult__(const MotionDense<MotionDerived> & v, ForceDense<ForceDerived> & f) const
    {
      f.linear().noalias() = lever().cross(v.angular() );
      f.angular().noalias() = inertia()*v.angular();
    }

    template<typename MotionDerived, typename ForceDerived>
    void __transpose_mult__(const MotionDense<MotionDerived> & v, ForceDense<ForceDerived> & f) const
    {
      f.linear().setZero();
      f.angular() = inertia().transpose()*v.angular() + v.linear().cross( lever() );
    }
    
    template<typename MotionDerived, typename ForceDerived>
    void __transpose_mult_add__(const MotionDense<MotionDerived> & v, ForceDense<ForceDerived> & f) const
    {
      f.angular() += inertia().transpose()*v.angular() + v.linear().cross( lever() );
    }

    template<typename MotionDerived, typename ForceDerived>
    void __transpose_mult_sub__(const MotionDense<MotionDerived> & v, ForceDense<ForceDerived> & f) const
    {
      f.angular() -= inertia().transpose()*v.angular() + v.linear().cross( lever() );
    }



    // Getters
    const Vector3 &    lever()   const { return m_com; }
    const Matrix3 & inertia() const { return m_inertia; }
    
    Vector3 &    lever()   { return m_com; }
    Matrix3 & inertia() { return m_inertia; }



    void disp_impl(std::ostream & os) const
    {
      os
      << "  c = " << lever().transpose() << "\n"
      << "  I = \n" << inertia().matrix() << "";
    }
    
    /// \returns An expression of *this with the Scalar type casted to NewScalar.
    template<typename NewScalar>
    CoriolisTpl<NewScalar,Options> cast() const
    {
      return CoriolisTpl<NewScalar,Options>(lever().template cast<NewScalar>(),
                                           inertia().template cast<NewScalar>());
    }
    

  protected:
    Vector3 m_com;
    Matrix3 m_inertia;

  }; // class CoriolisTpl
    
} // namespace pinocchio

#endif // ifndef __pinocchio_inertia_hpp__
