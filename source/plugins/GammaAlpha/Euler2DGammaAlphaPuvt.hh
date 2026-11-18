#ifndef COOLFluiD_Physics_GammaAlpha_Euler2DGammaAlphaPuvt_hh
#define COOLFluiD_Physics_GammaAlpha_Euler2DGammaAlphaPuvt_hh

//////////////////////////////////////////////////////////////////////////////

#include "NavierStokes/MultiScalarVarSet.hh"
#include "NavierStokes/Euler2DVarSet.hh"
#include "NavierStokes/EulerTerm.hh"
#include "NavierStokes/Euler2DPuvt.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace GammaAlpha {

//////////////////////////////////////////////////////////////////////////////

  /**
   * This class represents a Euler physical model 2D for primitive
   * variables with Gamma-Alpha transition Model coupled with k-LogOmega turbulence model
   *
   * @author Khalil Bensassi
   * @author Ray Vandenhoeck
   */

class Euler2DGammaAlphaPuvt : public Physics::NavierStokes::MultiScalarVarSet<Physics::NavierStokes::Euler2DPuvt> {
public: // classes

  /**
   * Puvttructor
   * @see Euler2D
   */
  Euler2DGammaAlphaPuvt(Common::SafePtr<Framework::BaseTerm> term);

  /**
   * Default destructor
   */
  ~Euler2DGammaAlphaPuvt();

  /**
   * Set up the private data and give the maximum size of states physical
   * data to store
   */
  void setup();

  /**
   * Get extra variable names
   */
  std::vector<std::string> getExtraVarNames() const;

  /**
   * Gets the block separator for this variable set
   */
  CFuint getBlockSeparator() const;

  /**
   * Split the jacobian
   */
  void splitJacobian(RealMatrix& jacobPlus,
		     RealMatrix& jacobMin,
		     RealVector& eValues,
		     const RealVector& normal);
  
  /**
   * Set the matrix of the right eigenvectors and the matrix of the eigenvalues
   */
  void computeEigenValuesVectors(RealMatrix& rightEv,
				 RealMatrix& leftEv,
				 RealVector& eValues,
				 const RealVector& normal);
  
  /**
   * Get the speed
   */
  CFreal getSpeed(const Framework::State& state) const;

  /**
   * Give dimensional values to the adimensional state variables
   */
  void setDimensionalValues(const Framework::State& state,
                            RealVector& result);

  /**
   * Give adimensional values to the dimensional state variables
   */
  virtual void setAdimensionalValues(const Framework::State& state,
                                    RealVector& result);


  /**
   * Set other adimensional values for useful physical quantities
   */
  void setDimensionalValuesPlusExtraValues
  (const Framework::State& state, RealVector& result,
   RealVector& extra);
  
  /// Compute the perturbed physical data
  virtual void computePerturbedPhysicalData(const Framework::State& state,
					    const RealVector& pdataBkp,
					    RealVector& pdata,
					    CFuint iVar);
  
  /**
   * Set the PhysicalData corresponding to the given State
   * @see EulerPhysicalModel
   */
  /// Physical flux; adds the wall-normal alpha transport rho*alpha*(w.n), w = c*(v,-u), when the
  /// N_alpha flux form is active (nAlphaFluxCoefficient() > 0, set by GammaAlpha2DSourceTerm), and the
  /// transport of gamma away from the wall gamma*(rho w_E).n of E_gamma (GammaAlphaEGammaFlux.hh) when the
  /// E_gamma flux form is active (eGammaFluxCoefficient() > 0).
  using Physics::NavierStokes::MultiScalarVarSet<Physics::NavierStokes::Euler2DPuvt>::computeFlux;
  virtual void computeFlux(const RealVector& pdata, const RealVector& normals);

  virtual void computePhysicalData(const Framework::State& state,
			       RealVector& data);
  
  /**
   * Set a State starting from the given PhysicalData
   * @see EulerPhysicalModel
   */
  virtual void computeStateFromPhysicalData(const RealVector& data,
				       Framework::State& state);
     
  /// Set the IDs corresponding to the velocity components in a State
  virtual void setStateVelocityIDs (std::vector<CFuint>& velIDs)
  {
    velIDs.resize(2); velIDs[XX] = 1; velIDs[YY] = 2;
  }
  
}; // end of class Euler2DGammaAlphaPuvt

/// Coefficient C_Nalpha of the gamma-alpha N_alpha term when it is evaluated in flux form (0 = pointwise
/// source form, the original implementation). One definition in libGammaAlpha, shared by the variable set,
/// the GA Riemann flux and the GA source term (which sets it from its NAlphaFlux option).
CFreal& nAlphaFluxCoefficient();
/// Coefficient C_Egamma of the 2025 E_gamma (eq. 37) when it is evaluated as a transport of gamma away from the
/// wall in flux form (GammaAlphaEGammaFlux.hh); 0 = the pointwise source of GammaAlpha2DSourceTerm. Shared like
/// nAlphaFluxCoefficient() (set by GammaAlpha2DSourceTerm from its EGammaFlux option).
CFreal& eGammaFluxCoefficient();
/// k and omega diffuse with the SST mu_t (mu_t*/f(gamma), A4) instead of mu_t*: with the SST k production of
/// GammaAlpha2DSourceTerm this is SST not changed (Jatayu v1.5 gamma_alpha_k_production = sst). Set by
/// GammaAlpha2DSourceTerm (KOmegaDiffusionSST); default false.
bool& gammaAlphaKOmegaDiffusionSST();
/// Laminar flow share of the mean-flow eddy viscosity (Jatayu v1.5 gamma_alpha_laminar_flow_share, not in the paper):
/// the share of mu_t* (A4) the mean flow gets at gamma = 0.01, all of it from gamma = 0.05 on (GammaAlphaMuTFactor.hh).
/// 1 = A4 as printed (default). Set by GammaAlpha2DSourceTerm (LaminarFlowShare).
CFreal& gammaAlphaLaminarFlowShare();

/// Gamma-alpha formulation: 2018 = Van den Eynde & Steelant / Vandenhoeck thesis (legacy default),
/// 2025 = Jacobs, Van den Eynde & Steelant, CEAS Space J. 2025. Set by GammaAlpha2DSourceTerm (ModelVersion).
CFuint& gammaAlphaModelVersion();
      
//////////////////////////////////////////////////////////////////////////////

    } // namespace GammaAlpha

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_GammaAlpha_Euler2DGammaAlphaPuvt_hh
