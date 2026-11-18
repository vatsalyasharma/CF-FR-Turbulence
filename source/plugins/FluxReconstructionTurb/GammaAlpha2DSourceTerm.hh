#ifndef COOLFluiD_FluxReconstructionMethod_GammaAlpha2DSourceTerm_hh
#define COOLFluiD_FluxReconstructionMethod_GammaAlpha2DSourceTerm_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/VectorialFunction.hh"
#include "FluxReconstructionTurb/KLogOmega2DSourceTerm.hh"
#include "Common/SafePtr.hh"
#include "NavierStokes/MultiScalarVarSet.hh"
#include "NavierStokes/Euler2DVarSet.hh"
#include "KOmega/NavierStokesKLogOmegaVarSetTypes.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {
    class GeometricEntity;
  }
  
    namespace FluxReconstructionMethod {
      
//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents the 2D Gamma-Re_theta-K-Omega Source Term for FR
 *
 * @author Ray Vandenhoeck
 *
 */
class GammaAlpha2DSourceTerm : 
	public KLogOmega2DSourceTerm {
  
public:
  
  /**
   * Constructor
   */
  GammaAlpha2DSourceTerm(const std::string& name);

  /**
   * Default destructor
   */
  ~GammaAlpha2DSourceTerm();
  
  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Set up the member data
   */
  virtual void setup();

  /**
   * UnSet up private data and data of the aggregated classes
   * in this command after processing phase
   */
  virtual void unsetup();
  
  /**
   * add the source term
   */
  void addSourceTerm(RealVector& resUpdates);

  /**
   * With OnsetFunction: gather the wall onset table over the ranks, then the usual source-term loop
   */
  virtual void execute();
  
  /// Returns the DataSocket's that this command provides as sources
  /// @return a vector of SafePtr with the DataSockets
  virtual std::vector< Common::SafePtr< Framework::BaseDataSocketSource > >
    providesSockets();
  
  /// Returns the DataSocket's that this command needs as sinks
  /// @return a vector of SafePtr with the DataSockets
  std::vector< Common::SafePtr< Framework::BaseDataSocketSink > >
    needsSockets();
  
  virtual void getSToStateJacobian(const CFuint iState);
  
  virtual void getSToGradJacobian(const CFuint iState){};
  
  virtual bool isGradDependent(){return true;};
  
private: // helper functions

  /// F_onset at a state (nearest wall flux point of the gathered table), cached per residual evaluation
  CFreal onsetFunction(const Framework::State& state);

  /// wall-law value at the projection of the state onto the wall (WallLaw; 0.01 before the first wall data)
  CFreal wallLawCap(const Framework::State& state);

  /// WallLaw: march of the wall intermittency along the wall over the gathered wall points (all ranks)
  void marchWallLaw();

  /// EGammaArtificialDiffusion: fill the artificial diffusivity of the gamma equation of every state (ghosts
  /// included) from the current solution (GammaAlpha/GammaAlphaArtificialDiffusion.hh)
  void updateGammaArtificialDiffusion();
    
  /// storage for local freestream M
  Framework::DataSocketSource<CFreal> socket_MInfLocal;
  
  /// storage for local freestream u
  Framework::DataSocketSource<CFreal> socket_uInfLocal;
  
  /// storage for local freestream T
  Framework::DataSocketSource<CFreal> socket_TInfLocal;
  
  /// storage for local freestream Tu
  Framework::DataSocketSource<CFreal> socket_TuInfLocal;
  
  /// storage for alpha_c-alpha
  Framework::DataSocketSource<CFreal> socket_alphaDiff;
  
  /// denominators of the coefficients for the update
  Framework::DataSocketSink<CFreal > socket_updateCoeff;
  
  /// corresponding diffusive variable set
  void getVorticity(const CFuint iState);

  void getStrain(const CFuint iState);

  CFreal getRethetat(const CFreal Tu, bool prime);  
  
  CFreal getRethetatwithPressureGradient(const CFreal avAlpha, const CFreal rhoInfLocal, const CFreal uInfLocal, const CFreal muInfLocal, const CFreal mu, const CFreal rho, const CFreal dpx, const CFreal dpy);

  CFreal getTuInfLocal(const CFreal Rethetat, const CFreal MInfLocal, const CFreal TuLocal);

private: // data
  
  CFreal   m_Rethetat; 
  CFreal   m_Rethetac; 
  CFreal   m_Flength; 
  CFreal   m_vorticity;
  CFreal  m_strain;
  bool     m_PGrad; 
  CFreal m_FLambda;
  bool m_decouple;
  bool m_limPRe;
  bool m_addUpdateCoeff;
  bool m_addDGDA;

  /// use the unclamped gamma wherever it carries the equilibrium at gamma = 1, keeping the
  /// [0.01,0.99] clamp only where it guards log(1-gamma). Default false = old behaviour
  bool m_gammaEquilibrium;
  /// accept k<0 at solution points (the residual already uses max(k,0))
  bool m_signedK;
  /// N_alpha in flux form (AUSMPlusUpFluxGammaAlpha2D + Euler2DGammaAlphaPuvt) instead of a pointwise source
  bool m_nAlphaFlux;
  /// 2018 (thesis/Van den Eynde & Steelant) or 2025 (Jacobs, Van den Eynde & Steelant) formulation
  CFuint m_modelVersion;
  CFreal  m_lambdaLim;
  CFreal  m_ceg;
  
  CFuint m_order;

  /// P_gamma multiplied by a smoothed F_onset (2025 paper eq. A7) from the nearest wall flux point's
  /// tau_w* <= alpha criterion; use with a zero-flux wall gamma (Wall.GammaThreshold very negative)
  bool m_onsetFunction;
  /// relative width w of F_onset = (1 + tanh(((alpha - tau_w*)/alpha)/w))/2
  CFreal m_onsetWidth;
  /// wall table gathered over all ranks: x, y, margin per wall flux point
  std::vector<CFreal> m_onsetPoints;
  /// F_onset per state local ID, valid while m_onsetCacheStamp equals m_onsetStamp
  std::vector<CFreal> m_onsetCache;
  std::vector<CFuint> m_onsetCacheStamp;
  CFuint m_onsetStamp;
  /// numerical source Jacobian with the cell-local gradient dependence: in perturbed evaluations the stored
  /// gradients get the FR volume-term change dq_j grad(L_j) (makes E_gamma, N_alpha... implicit within the cell)
  bool m_gradientJacobian;
  /// 2025: E_gamma fades linearly to zero between GammaFloor + 0.01 and GammaFloor (default 0)
  CFreal m_gammaFloor;
  /// SENSITIVITY (not in the paper): E_gamma multiplied by (1 - F_onset), i.e. the sink acts only before onset
  bool m_sinkPreOnsetOnly;
  /// wall law (Jatayu v1.5 gamma_alpha_wall_law; requires OnsetFunction): the intermittency at the wall is 0.01 up to
  /// the place where the onset criterion (A7) is first met along the wall and follows from there the gamma equation
  /// at the wall, dgamma/ds = 2 f_gamma (1-gamma) sqrt(-ln(1-gamma)) beta0; BCNoSlipWallGammaAlpha2D (WallLaw) imposes
  /// it, and the production of a cell is switched on behind that place (not by the local criterion)
  bool m_wallLaw;
  /// k production of SST unchanged (Jatayu v1.5 gamma_alpha_k_production = sst): P_k with the SST mu_t = mu_t*/f(gamma)
  /// (A4 factor) and the plain (2/3) rho k div(u); the mean flow keeps mu_t*. Default false: P_k with mu_t* (effective)
  bool m_kProductionSST;
  /// k and omega diffuse with the SST mu_t instead of mu_t* (sets Physics::GammaAlpha::gammaAlphaKOmegaDiffusionSST());
  /// with KProductionSST: SST not changed (Jatayu v1.5 gamma_alpha_k_production = sst). Default false
  bool m_kOmegaDiffusionSST;
  /// laminar flow share of the mean-flow eddy viscosity at gamma = 0.01 (Jatayu v1.5, fitted, not in the paper);
  /// 1 = A4 as printed (sets Physics::GammaAlpha::gammaAlphaLaminarFlowShare())
  CFreal m_laminarFlowShare;
  /// numbers per gathered wall point: x, y, margin (3); with WallLaw also alpha and p (5)
  CFuint m_onsetStride;
  /// per gathered wall point: production switch behind the onset place, and the wall value of the law
  std::vector<CFreal> m_wallBehind;
  std::vector<CFreal> m_wallGamma;
  /// WallLaw: wall points with x at or below this are laminar (gamma_w = 0.01, no production) and the march starts
  /// anew behind them (Jatayu v1.5 gamma_alpha_laminar_until_x; a very large value gives a laminar wall everywhere)
  CFreal m_wallLawLaminarUntilX;
  /// x where the march starts (-1: criterion met nowhere)
  CFreal m_wallLawStart;
  /// WallLaw + 2025: keep the part of E_gamma (eq. 37) that carries the wall intermittency into the layer
  /// (dgamma/dn < 0 makes -E_gamma > 0), behind the onset, faded to zero as gamma approaches the wall-law value along
  /// the wall; Jatayu v1.5 transports the wall gamma away from the wall with an unclipped upwind E_gamma. The sink
  /// part (dgamma/dn > 0) is unchanged. Default false
  bool m_eGammaWallSource;
  /// 2025: in the numerical source Jacobian (perturbed evaluations only; the residual is unchanged) the gamma equation
  /// gets the upwind diagonal -rho |w_E| / dn of the E_gamma advection (w_E = 0.25 mu_gamma 2c sqrt(1-c/c_e)/mu, dn the
  /// solution-point spacing along the wall normal), so the gradient-dependent E_gamma terms, explicit otherwise, do
  /// not limit the CFL number. Default false
  bool m_eGammaImplicitDiagonal;
  /// 2025: factor C of the artificial diffusivity C rho |w_E| dn/2 of the gamma equation (0 = off)
  CFreal m_eGammaArtificialDiffusion;
  /// EGammaArtificialDiffusion weighted by the wall-law intermittency at the state's projection onto the wall:
  /// 0 up to Start, 1 from End on (Start < 0: weight 1 everywhere)
  CFreal m_eGammaADWallGammaStart;
  CFreal m_eGammaADWallGammaEnd;
  /// 2025 paper eqs. 10-16 (Karsch et al.): the spot growth factor f_sigma in place of f_M,n_sigma (eq. 9), in the
  /// gamma production and the wall law; the paper uses it for supersonic cases only. Default false
  bool m_karschSpotGrowth;
  /// wall temperature T_w of eq. 16 [K] (isothermal wall) and the recovery factor r of eq. 16 (default sqrt(0.72))
  CFreal m_karschWallTemperature;
  CFreal m_karschRecoveryFactor;
  /// per state local ID (valid while m_onsetCacheStamp equals m_onsetStamp): the wall-law value interpolated along
  /// the wall at the projection of the state
  std::vector<CFreal> m_wallCapCache;
  /// the gathered wall points in their order along the wall, and the place of each point in that order
  std::vector<CFuint> m_wallOrder;
  std::vector<CFuint> m_wallRank;
  /// 2025: E_gamma (eq. 37) as a transport of gamma away from the wall in flux form (GammaAlphaEGammaFlux.hh):
  /// gamma rho w_E in the convective flux (AUSMPlusUpFluxGammaAlpha2D, Euler2DGammaAlphaPuvt), gamma div(rho w_E)
  /// here, with the FR volume divergence of the cell; replaces the pointwise (clipped) E_gamma
  bool m_eGammaFlux;
  /// 2025: E_gamma (eq. 37) in upwind non-conservative form: -rho w_E.grad(gamma) here with the central gradient,
  /// and the upwind dissipation -|rho w_E.n|/2 jump(gamma) in the gamma interface flux (AUSMPlusUpFluxGammaAlpha2D);
  /// not clipped; replaces the pointwise clipped E_gamma
  bool m_eGammaUpwind;
  /// 2025, with the pointwise E_gamma (default form or EGammaWallSource): factor of the upwind dissipation
  /// -|rho w_E.n|/2 jump(gamma) added to the gamma interface flux (AUSMPlusUpFluxGammaAlpha2D); 0 = off
  CFreal m_eGammaFaceDissipation;
  /// rho w_E at the solution points of the current cell, and gamma's local term div(rho w_E) there
  std::vector<CFreal> m_eGammaFx;
  std::vector<CFreal> m_eGammaFy;
  std::vector<CFreal> m_eGammaDivergence;
  std::vector<RealVector> m_unpertStates;
  Common::SafePtr< std::vector< std::vector< std::vector< CFreal > > > > m_solPolyDerivAtSolPnts;
  /// metric (projection vectors) at the solution points
  Framework::DataSocketSink< CFreal > socket_solPntNormals;
  /// all states and the solution point Jacobian determinants (EGammaArtificialDiffusion)
  Framework::DataSocketSink< Framework::State*, Framework::GLOBAL > socket_states;
  Framework::DataSocketSink< CFreal > socket_volumes;
  /// diagnostic: log the gamma source pieces at near-wall states within 0.02 m of this x (default off)
  CFreal m_debugProbeX;
  std::vector<CFuint> m_debugLogged;

}; // end of class GammaAlpha2DSourceTerm

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluxReconstructionMethod

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_FluxReconstructionMethod_NavierStokesGammaAlpha2DSourceTerm_hh
