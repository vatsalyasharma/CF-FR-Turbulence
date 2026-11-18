#ifndef COOLFluiD_Physics_GReKO_NavierStokes3DGReKLogOPuvt_hh
#define COOLFluiD_Physics_GReKO_NavierStokes3DGReKLogOPuvt_hh

//////////////////////////////////////////////////////////////////////////////

#include "GReKO/LM2009Source.hh"
#include "KOmega/NavierStokesKLogOmegaVarSet.hh"
#include "KOmega/NavierStokesKLogOmegaPvt.hh"
#include "NavierStokes/NavierStokes3DVarSet.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace GReKO {

//////////////////////////////////////////////////////////////////////////////

  /**
   * This class represents a NavierStokes physical model 3D for primitive
   * variables and Gamma-Re trasition model coupled with K-LogOmega turbulence model
   *
   * @author Khalil Bensassi
   * @author Ray Vandenhoeck
   */

class NavierStokes3DGReKLogOPuvt : public 
KOmega::NavierStokesKLogOmegaPvt<KOmega::NavierStokesKLogOmegaVarSet<NavierStokes::NavierStokes3DVarSet, 0> > {
public: // classes

  /// Point-local constitutive inputs for the shared LM2009 residual/Jacobian.
  lm2009::Inputs getLM2009Inputs(const RealVector& state, const std::vector<RealVector*>& gradients);

  /// Checked transition transport; overrides the existing virtual slot.
  virtual CFreal getLaminarDynViscosityFromGradientVars(const RealVector& state);

  virtual CFreal getDensity(const RealVector& state);

  /// rho*max(k,0)/mu_t evaluated after cancelling k, including k=0.
  /// Nonvirtual dispatch in the primitive base preserves the existing ABI.
  CFreal getRhoKOverMuT(const RealVector& state, const std::vector<RealVector*>& gradients);

  
  /**
   * Constructor
   * @see NavierStokes3D
   */
  NavierStokes3DGReKLogOPuvt(const std::string& name,
			  Common::SafePtr<Framework::PhysicalModelImpl> model);

  /**
   * Default destructor
   */
  ~NavierStokes3DGReKLogOPuvt();
  
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

}; // end of class NavierStokes3DGReKLogOPuvt

//////////////////////////////////////////////////////////////////////////////

    } // namespace GReKO

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_GReKO_NavierStokes3DGReKLogOPuvt_hh
