#ifndef COOLFluiD_Physics_GReKO_NavierStokes2DGReKLogOPuvt_hh
#define COOLFluiD_Physics_GReKO_NavierStokes2DGReKLogOPuvt_hh

//////////////////////////////////////////////////////////////////////////////

#include "GReKO/TransitionOmegaSource.hh"
#include "GReKO/LM2009Source.hh"
#include "KOmega/NavierStokesKLogOmegaVarSet.hh"
#include "KOmega/NavierStokesKLogOmegaPvt.hh"
#include "NavierStokes/NavierStokes2DVarSet.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace GReKO {

//////////////////////////////////////////////////////////////////////////////

  /**
   * This class represents a NavierStokes physical model 2D for primitive
   * variables and Gamma-Re trasition model coupled with K-LogOmega turbulence model
   *
   * @author Khalil Bensassi
   * @author Ray Vandenhoeck
   */

class NavierStokes2DGReKLogOPuvt : public 
KOmega::NavierStokesKLogOmegaPvt<KOmega::NavierStokesKLogOmegaVarSet<NavierStokes::NavierStokes2DVarSet, 0> > {
public: // classes

  /// Point-local constitutive inputs for the shared LM2009 residual/Jacobian.
  lm2009::Inputs getLM2009Inputs(const RealVector& state, const std::vector<RealVector*>& gradients);

  /// Checked transition transport; overrides the existing virtual slot.
  virtual CFreal getLaminarDynViscosityFromGradientVars(const RealVector& state);

  virtual CFreal getDensity(const RealVector& state);

  /// rho*max(k,0)/mu_t evaluated after cancelling k, including k=0.
  /// Nonvirtual dispatch in the primitive base preserves the existing ABI.
  CFreal getRhoKOverMuT(const RealVector& state, const std::vector<RealVector*>& gradients);


  /// Local k/log-omega derivatives of the implemented volume log-omega source.
  /// Gradients, p/T and other state components are fixed; no cache is changed.
  omega_detail::Partial getLogOmegaSourcePartial(const RealVector& state,
      const std::vector<RealVector*>& gradients, CFreal shearProduction,
      CFreal dilatation, omega_detail::Partial destructionFactor, bool limitProduction);

  /**
   * Constructor
   * @see NavierStokes2D
   */
  NavierStokes2DGReKLogOPuvt(const std::string& name,
			  Common::SafePtr<Framework::PhysicalModelImpl> model);

  /**
   * Default destructor
   */
  ~NavierStokes2DGReKLogOPuvt();
  
  virtual void setComposition(const RealVector& state,
			      const bool isPerturb,
			      const CFuint iVar)
  {
    _isPerturb = isPerturb;
    _iPerturbVar = iVar;
  }
  
  /**
   * Get number of turbulent variables
   */
  virtual CFuint getNbTurbVars() const 
  {
    const CFuint nbTurbVars = _eulerModel->getNbScalarVars(0);
    cf_assert(nbTurbVars == 4);
    return nbTurbVars;
  }
  
  /**
   * Get the diffusive flux
   */
  virtual RealVector& getFlux(const RealVector& values,
                              const std::vector<RealVector*>& gradients,
                              const RealVector& normal,
                              const CFreal& radius);

private :
  
  CFreal _unperturbedFluxGa;
  CFreal _unperturbedFluxRe;

}; // end of class NavierStokes2DGReKLogOPuvt

//////////////////////////////////////////////////////////////////////////////

    } // namespace GReKO

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_GReKO_NavierStokes2DGReKLogOPuvt_hh
