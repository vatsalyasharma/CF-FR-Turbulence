#include "NavierStokes/TransitionTransport.hh"
#include "FluxReconstructionTurb/TransitionSourceScales.hh"
#include "FluxReconstructionTurb/GammaAlpha2DSourceTerm.hh"
#include "FluxReconstructionTurb/GammaAlphaOnsetTable.hh"
#include "Common/PE.hh"
#include <limits>
#include <sstream>
#include <algorithm>
#include "Common/CFLog.hh"
#include "Framework/GeometricEntity.hh"
#include "Framework/MeshData.hh"
#include "FluxReconstructionTurb/FluxReconstructionKOmega.hh"
#include "Framework/SubSystemStatus.hh"
#include "FluxReconstructionTurb/KLogOmega2DSourceTerm.hh"
#include "GammaAlpha/NavierStokes2DGammaAlphaPuvt.hh"
#include "GammaAlpha/Euler2DGammaAlphaPuvt.hh"
#include "GammaAlpha/GammaAlphaEGammaFlux.hh"
#include "GammaAlpha/GammaAlphaMuTFactor.hh"
#include "GammaAlpha/GammaAlphaArtificialDiffusion.hh"
#include "FluxReconstructionTurb/AUSMPlusUpFluxGammaAlpha2D.hh"

#include "Framework/MethodCommandProvider.hh"
#include "MathTools/MathConsts.hh"
#include "KOmega/NavierStokesKLogOmegaVarSetTypes.hh"

#include "FluxReconstructionMethod/FluxReconstructionElementData.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Physics::NavierStokes;
using namespace COOLFluiD::Physics::KOmega;
using namespace COOLFluiD::Physics::GammaAlpha;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace FluxReconstructionMethod {

////////////////////////////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<GammaAlpha2DSourceTerm, FluxReconstructionSolverData, FluxReconstructionKOmegaModule>
GammaAlpha2DSourceTermFRProvider("GammaAlpha2DSourceTerm");

///////////////////////////////////////////////////////////////////////////////////////////////////

void GammaAlpha2DSourceTerm::defineConfigOptions(Config::OptionList& options)
{
 options.addConfigOption< bool >("PGrad","pressure Gradient");
 options.addConfigOption< bool >("Decouple","Decouple y-ReTheta from k and log(omega), simply solving as fully turbulent.");
 options.addConfigOption< bool >("LimPAlpha","Limit P_alpha.");
 options.addConfigOption< bool >("AddUpdateCoeff","Add the ST time step restriction.");
 options.addConfigOption< bool,Config::DynamicOption<> >("AddDGamma","Add destruction terms for gamma and alpha.");
 options.addConfigOption< CFreal >("LimLambda","Limit Lambda pressure term.");
 options.addConfigOption< CFreal >("CEg","CEg term to be used for gamma-alpha.");
 options.addConfigOption< bool >("GammaEquilibrium","Use the unclamped gamma in the (1-gamma) production prefactor and in the axisymmetric transport, so gamma = 1 is an equilibrium.");
 options.addConfigOption< bool >("SignedK","Accept k<0 at solution points; the source uses max(k,0).");
 options.addConfigOption< bool >("NAlphaFlux","Evaluate N_alpha = 0.4 rho |U| dalpha/dn in flux form (requires RiemannFlux AUSMPlusUpFluxGammaAlpha2D and the numerical source Jacobian); default false = pointwise source.");
 options.addConfigOption< bool >("OnsetFunction","Multiply P_gamma by a smoothed F_onset (Jacobs et al. 2025, eq. A7) from the nearest wall flux point's tau_w* <= alpha criterion. Use with a zero-flux wall gamma (Wall.GammaThreshold very negative). Default false.");
 options.addConfigOption< CFreal >("OnsetWidth","Relative width w of F_onset = (1 + tanh(((alpha - tau_w*)/alpha)/w))/2 (default 0.05).");
 options.addConfigOption< bool >("EGammaFlux","ModelVersion 2025: E_gamma (eq. 37) as a transport of gamma away from the wall in flux form, as Jatayu v1.5 evaluates it (requires RiemannFlux AUSMPlusUpFluxGammaAlpha2D, AnalyticalJacob = false and AddDGamma = true): gamma rho w_E is added to the convective flux, gamma div(rho w_E) here; the pointwise E_gamma and its clip are not used. Default false.");
 options.addConfigOption< CFreal >("LaminarFlowShare","Share of mu_t* (A4) that the mean flow gets at gamma = 0.01, all of it from gamma = 0.05 on (Jatayu v1.5 gamma_alpha_laminar_flow_share; fitted, not in the paper). Default 1 = A4 as printed.");
 options.addConfigOption< bool >("KOmegaDiffusionSST","k and omega diffuse with the SST mu_t = mu_t*/f(gamma) (A4) instead of mu_t*; with KProductionSST the SST is not changed (Jatayu v1.5 gamma_alpha_k_production = sst). Default false.");
 options.addConfigOption< bool >("KProductionSST","k production of SST unchanged (Jatayu v1.5 gamma_alpha_k_production = sst): P_k with mu_t = mu_t*/f(gamma) (A4) and (2/3) rho k div(u); the mean flow keeps mu_t*. Default false (P_k with mu_t*).");
 options.addConfigOption< CFreal >("WallLawLaminarUntilX","WallLaw: wall points with x at or below this value are laminar (gamma_w = 0.01, production off) and the march starts behind them (Jatayu v1.5 gamma_alpha_laminar_until_x); a very large value gives a laminar wall everywhere. Default -1e30 (none).");
 options.addConfigOption< bool >("WallLaw","With OnsetFunction: the intermittency at the wall by the wall law (Jatayu v1.5 gamma_alpha_wall_law): 0.01 up to the first wall point where tau_w* <= alpha, then dgamma/ds = 2 f_gamma (1-gamma) sqrt(-ln(1-gamma)) beta0 along the wall (flat wall ordered by x); imposed by BCNoSlipWallGammaAlpha2D (WallLaw); the production is switched on behind that place. Default false.");
 options.addConfigOption< bool >("EGammaUpwind","ModelVersion 2025: E_gamma (eq. 37) as an upwind transport of gamma away from the wall in non-conservative form (requires RiemannFlux AUSMPlusUpFluxGammaAlpha2D, AnalyticalJacob = false, AddDGamma = true; not with EGammaFlux): -rho w_E.grad(gamma) here with the central gradient, the upwind dissipation in the gamma interface flux; not clipped. Default false.");
 options.addConfigOption< CFreal >("EGammaFaceDissipation","ModelVersion 2025, with the pointwise E_gamma (default form or EGammaWallSource; not with EGammaUpwind or EGammaFlux; requires RiemannFlux AUSMPlusUpFluxGammaAlpha2D and AnalyticalJacob = false): add this factor times the upwind dissipation -|rho w_E.n|/2 (gamma_R - gamma_L) of the E_gamma advection to the gamma interface flux (1 = the dissipation of the upwind face flux, as in EGammaUpwind and the upwind transport of Jatayu v1.5). Zero for a continuous gamma, so consistent; it damps inter-element gamma oscillations that the centred pointwise E_gamma leaves undamped. Default 0 (off).");
 options.addConfigOption< CFreal >("EGammaArtificialDiffusion","ModelVersion 2025 (requires the turbulence diffusion commands of FluxReconstructionTurb and the NavierStokes2DGammaAlpha variable sets): artificial diffusivity C rho |w_E| dn/2 added to the gamma equation only, with this factor C; rho w_E = 0.25 mu_gamma 2c sqrt(1 - c/c_e) rho/mu n is the E_gamma advection (eq. 37) and dn the solution-point spacing along n. C = 1 is the numerical diffusion of first-order upwinding of the E_gamma advection over one solution-point spacing; O(h), so consistent. It is zero at walls (c = 0) and outside the layers (c = c_e), and no artificial flux crosses the boundaries. Evaluated from the solution at the start of each residual evaluation, so frozen in the Jacobian. Default 0 (off).");
 options.addConfigOption< CFreal >("EGammaADWallGammaStart","EGammaArtificialDiffusion with WallLaw: weight the artificial diffusivity by the wall-law intermittency gamma_w interpolated along the wall at the state's projection, linearly from 0 at gamma_w = EGammaADWallGammaStart to 1 at EGammaADWallGammaEnd, so that it acts only where the layer is (nearly) turbulent and leaves the transition front alone. Negative: weight 1 everywhere (default -1).");
 options.addConfigOption< CFreal >("EGammaADWallGammaEnd","End of the wall-law intermittency ramp of EGammaADWallGammaStart (default 0.9).");
 options.addConfigOption< bool >("EGammaWallSource","WallLaw + ModelVersion 2025: behind the onset keep the part of E_gamma (eq. 37) that carries the wall intermittency into the layer (dgamma/dn < 0 makes -E_gamma a source), faded to zero as gamma approaches the wall-law value interpolated along the wall; Jatayu v1.5 transports the wall gamma away from the wall with an unclipped upwind E_gamma. The sink part is unchanged. Default false.");
 options.addConfigOption< bool >("EGammaImplicitDiagonal","ModelVersion 2025: in the numerical source Jacobian only (residual unchanged), add to the gamma equation the upwind diagonal -rho |w_E|/dn of the E_gamma advection (dn: solution-point spacing along the wall normal), so the explicit gradient-dependent E_gamma terms do not limit the CFL number. Default false.");
 options.addConfigOption< bool >("KarschSpotGrowth","ModelVersion 2025: spot growth factor f_sigma of Karsch et al. (paper eqs. 10-16) in place of f_M,n_sigma (eq. 9), in the gamma production and the wall law; the paper uses it for supersonic cases. Needs KarschWallTemperature. Default false.");
 options.addConfigOption< CFreal >("KarschWallTemperature","Wall temperature T_w [K] of eq. 16 (isothermal wall) for KarschSpotGrowth.");
 options.addConfigOption< CFreal >("KarschRecoveryFactor","Recovery factor r of eq. 16 for KarschSpotGrowth (default sqrt(0.72)).");
 options.addConfigOption< bool >("SinkPreOnsetOnly","SENSITIVITY OPTION, not in Jacobs et al. 2025: multiply E_gamma by (1 - F_onset) so the sink acts only before onset (requires OnsetFunction). Default false.");
 options.addConfigOption< CFreal >("GammaFloor","ModelVersion 2025: E_gamma fades linearly to zero between GammaFloor + 0.01 and GammaFloor, so the sink cannot push gamma below GammaFloor (the paper keeps gamma >= 0.01). Default 0.");
 options.addConfigOption< bool >("GradientJacobian","Numerical source Jacobian including the cell-local dependence of the gradients on the states (FR volume term); the residual is unchanged (default false).");
 options.addConfigOption< CFreal >("DebugProbeX","Diagnostic: log P_gamma, E_gamma, F_onset, f_gamma, beta, Tu at near-wall states (wall distance < 2 mm) within 0.02 m of this x, every 50 residual evaluations (OnsetFunction only). Default off.");
 options.addConfigOption< CFuint >("ModelVersion","2018 = thesis formulation (default); 2025 = Jacobs, Van den Eynde & Steelant (CEAS Space J. 2025): mu_gamma eq. 31, E_gamma eq. 37, F_Palpha (gamma-0.02)/0.98, alpha diffusion with mu_t, local freestream from inlet totals (App. B). Requires AnalyticalJacob = false.");
}

/////////////////////////////////////////////////////////////////////////////////////////////////////

GammaAlpha2DSourceTerm::GammaAlpha2DSourceTerm(const std::string& name) :
  KLogOmega2DSourceTerm(name),
  socket_MInfLocal("MInfLocal"),
  socket_uInfLocal("uInfLocal"),
  socket_TInfLocal("TInfLocal"),
  socket_TuInfLocal("TuInfLocal"),
  socket_alphaDiff("alphaDiff"),
  socket_updateCoeff("updateCoeff"),
  m_Rethetat(),
  m_Rethetac(),
  m_Flength(),
  m_vorticity(),
  m_strain(),
  m_order(),
  m_onsetPoints(),
  m_onsetCache(),
  m_onsetCacheStamp(),
  m_onsetStamp(0),
  m_onsetStride(3),
  m_wallBehind(),
  m_wallGamma(),
  m_wallLawStart(-1.),
  m_unpertStates(),
  m_solPolyDerivAtSolPnts(),
  socket_solPntNormals("solPntNormals"),
  socket_states("states"),
  socket_volumes("volumes"),
  m_debugLogged()
{ 
  addConfigOptionsTo(this);
  
  m_PGrad = true;
  setParameter("PGrad",&m_PGrad);
  
  m_decouple = false;
  setParameter("Decouple",&m_decouple);
  
  m_limPRe = false;
  setParameter("LimPAlpha",&m_limPRe);
  
  m_addUpdateCoeff = false;
  setParameter("AddUpdateCoeff",&m_addUpdateCoeff);
  
  m_addDGDA = true;
  setParameter("AddDGamma",&m_addDGDA);

  m_gammaEquilibrium = false;
  setParameter("GammaEquilibrium",&m_gammaEquilibrium);

  m_signedK = false;
  setParameter("SignedK",&m_signedK);

  m_nAlphaFlux = false;
  setParameter("NAlphaFlux",&m_nAlphaFlux);

  m_modelVersion = 2018;
  setParameter("ModelVersion",&m_modelVersion);
  
  m_lambdaLim = 0.04;
  setParameter("LimLambda",&m_lambdaLim);
  
  m_ceg = 20.0;
  setParameter("CEg",&m_ceg);

  m_onsetFunction = false;
  setParameter("OnsetFunction",&m_onsetFunction);

  m_onsetWidth = 0.05;
  setParameter("OnsetWidth",&m_onsetWidth);

  m_eGammaFlux = false;
  setParameter("EGammaFlux",&m_eGammaFlux);

  m_eGammaUpwind = false;
  setParameter("EGammaUpwind",&m_eGammaUpwind);

  m_eGammaFaceDissipation = 0.0;
  setParameter("EGammaFaceDissipation",&m_eGammaFaceDissipation);

  m_wallLaw = false;
  setParameter("WallLaw",&m_wallLaw);

  m_wallLawLaminarUntilX = -1.0e30;
  setParameter("WallLawLaminarUntilX",&m_wallLawLaminarUntilX);

  m_kProductionSST = false;
  setParameter("KProductionSST",&m_kProductionSST);

  m_kOmegaDiffusionSST = false;
  setParameter("KOmegaDiffusionSST",&m_kOmegaDiffusionSST);

  m_laminarFlowShare = 1.0;
  setParameter("LaminarFlowShare",&m_laminarFlowShare);

  m_sinkPreOnsetOnly = false;
  setParameter("SinkPreOnsetOnly",&m_sinkPreOnsetOnly);

  m_eGammaWallSource = false;
  setParameter("EGammaWallSource",&m_eGammaWallSource);

  m_eGammaImplicitDiagonal = false;
  setParameter("EGammaImplicitDiagonal",&m_eGammaImplicitDiagonal);

  m_eGammaArtificialDiffusion = 0.0;
  setParameter("EGammaArtificialDiffusion",&m_eGammaArtificialDiffusion);

  m_eGammaADWallGammaStart = -1.0;
  setParameter("EGammaADWallGammaStart",&m_eGammaADWallGammaStart);

  m_eGammaADWallGammaEnd = 0.9;
  setParameter("EGammaADWallGammaEnd",&m_eGammaADWallGammaEnd);

  m_karschSpotGrowth = false;
  setParameter("KarschSpotGrowth",&m_karschSpotGrowth);

  m_karschWallTemperature = -1.0;
  setParameter("KarschWallTemperature",&m_karschWallTemperature);

  m_karschRecoveryFactor = std::sqrt(0.72);
  setParameter("KarschRecoveryFactor",&m_karschRecoveryFactor);

  m_gammaFloor = 0.0;
  setParameter("GammaFloor",&m_gammaFloor);

  m_gradientJacobian = false;
  setParameter("GradientJacobian",&m_gradientJacobian);

  m_debugProbeX = -1.0e30;
  setParameter("DebugProbeX",&m_debugProbeX);
}

/////////////////////////////////////////////////////////////////////////////////////////////////

GammaAlpha2DSourceTerm::~GammaAlpha2DSourceTerm()
{
}

//////////////////////////////////////////////////////////////////////////////////////////////////

void GammaAlpha2DSourceTerm::setup()
{
  KLogOmega2DSourceTerm::setup();

  // N_alpha flux form: the coefficient is shared with the variable set and the GA Riemann flux.
  if (m_nAlphaFlux && m_useAnaJacob)
    throw Common::BadValueException(FromHere(), "GammaAlpha2DSourceTerm: NAlphaFlux requires AnalyticalJacob = false (the analytical Jacobian differentiates the pointwise N_alpha)");
  Physics::GammaAlpha::nAlphaFluxCoefficient() = m_nAlphaFlux ? 0.4 : 0.;
  if (m_modelVersion != 2018 && m_modelVersion != 2025)
    throw Common::BadValueException(FromHere(), "GammaAlpha2DSourceTerm: ModelVersion must be 2018 or 2025");
  if (m_modelVersion == 2025 && m_useAnaJacob)
    throw Common::BadValueException(FromHere(), "GammaAlpha2DSourceTerm: ModelVersion 2025 requires AnalyticalJacob = false");
  Physics::GammaAlpha::gammaAlphaModelVersion() = m_modelVersion;
  if (m_eGammaFlux && (m_modelVersion != 2025 || m_useAnaJacob || !m_addDGDA))
    throw Common::BadValueException(FromHere(), "GammaAlpha2DSourceTerm: EGammaFlux requires ModelVersion = 2025, AnalyticalJacob = false and AddDGamma = true");
  // E_gamma flux form: C_Egamma (eq. 37) shared with the variable set and the GA Riemann flux
  Physics::GammaAlpha::eGammaFluxCoefficient() = m_eGammaFlux ? 0.25 : 0.;
  // k and omega diffusion with the SST mu_t (shared with the GA diffusive variable set)
  Physics::GammaAlpha::gammaAlphaKOmegaDiffusionSST() = m_kOmegaDiffusionSST;
  if (!(m_laminarFlowShare > 0.0 && m_laminarFlowShare <= 1.0))
    throw Common::BadValueException(FromHere(), "GammaAlpha2DSourceTerm: LaminarFlowShare must be in (0, 1]");
  Physics::GammaAlpha::gammaAlphaLaminarFlowShare() = m_laminarFlowShare;
  if (m_eGammaUpwind && (m_modelVersion != 2025 || m_useAnaJacob || !m_addDGDA || m_eGammaFlux))
    throw Common::BadValueException(FromHere(), "GammaAlpha2DSourceTerm: EGammaUpwind requires ModelVersion = 2025, AnalyticalJacob = false, AddDGamma = true and EGammaFlux = false");
  if (m_eGammaFaceDissipation < 0.0 || (m_eGammaFaceDissipation > 0.0 && (m_modelVersion != 2025 || m_useAnaJacob || !m_addDGDA || m_eGammaUpwind || m_eGammaFlux)))
    throw Common::BadValueException(FromHere(), "GammaAlpha2DSourceTerm: EGammaFaceDissipation must be >= 0; > 0 requires ModelVersion = 2025, AnalyticalJacob = false, AddDGamma = true, EGammaUpwind = false and EGammaFlux = false");
  // E_gamma upwind non-conservative form (or the face dissipation alone, with the pointwise E_gamma): C_Egamma for
  // the interface dissipation of the GA Riemann flux
  AUSMPlusUpFluxGammaAlpha2D::eGammaUpwindCoefficient() = m_eGammaUpwind ? 0.25 : 0.25*m_eGammaFaceDissipation;
  if (m_wallLaw && !m_onsetFunction)
    throw Common::BadValueException(FromHere(), "GammaAlpha2DSourceTerm: WallLaw requires OnsetFunction = true");
  if (m_karschSpotGrowth && (m_modelVersion != 2025 || !(m_karschWallTemperature > 0.0)))
    throw Common::BadValueException(FromHere(), "GammaAlpha2DSourceTerm: KarschSpotGrowth requires ModelVersion = 2025 and KarschWallTemperature > 0");
  if (m_eGammaWallSource && (!m_wallLaw || m_modelVersion != 2025 || m_eGammaUpwind || m_eGammaFlux || !m_addDGDA))
    throw Common::BadValueException(FromHere(), "GammaAlpha2DSourceTerm: EGammaWallSource requires WallLaw = true, ModelVersion = 2025, AddDGamma = true, EGammaUpwind = false and EGammaFlux = false");
  if (m_eGammaArtificialDiffusion < 0.0 || (m_eGammaArtificialDiffusion > 0.0 && m_modelVersion != 2025))
    throw Common::BadValueException(FromHere(), "GammaAlpha2DSourceTerm: EGammaArtificialDiffusion must be >= 0; > 0 requires ModelVersion = 2025");
  if (m_eGammaADWallGammaStart >= 0.0 && (!m_wallLaw || !(m_eGammaADWallGammaEnd > m_eGammaADWallGammaStart)))
    throw Common::BadValueException(FromHere(), "GammaAlpha2DSourceTerm: EGammaADWallGammaStart >= 0 requires WallLaw = true and EGammaADWallGammaEnd > EGammaADWallGammaStart");
  Physics::GammaAlpha::gammaArtificialDiffusion().active = m_eGammaArtificialDiffusion > 0.0;
  m_onsetStride = m_wallLaw ? 5 : 3;
  
  DataHandle< CFreal > MInfLocal = socket_MInfLocal.getDataHandle();
  
  DataHandle< CFreal > uInfLocal = socket_uInfLocal.getDataHandle();
  
  DataHandle< CFreal > TInfLocal = socket_TInfLocal.getDataHandle();
  
  DataHandle< CFreal > TuInfLocal = socket_TuInfLocal.getDataHandle();
  
  DataHandle< CFreal > alphaDiff = socket_alphaDiff.getDataHandle();
  
  DataHandle< CFreal > wallDistance = socket_wallDistance.getDataHandle();
  
  // resize socket
  MInfLocal.resize(wallDistance.size());
  uInfLocal.resize(wallDistance.size());
  TInfLocal.resize(wallDistance.size());
  TuInfLocal.resize(wallDistance.size());
  alphaDiff.resize(wallDistance.size());
  
  // get the local FR data
  vector< FluxReconstructionElementData* >& frLocalData = getMethodData().getFRLocalData();
  cf_assert(frLocalData.size() > 0);
  // for now, there should be only one type of element
  cf_assert(frLocalData.size() == 1);
  
  const CFPolyOrder::Type order = frLocalData[0]->getPolyOrder();
  
  m_order = static_cast<CFuint>(order);

  m_solPolyDerivAtSolPnts = frLocalData[0]->getCoefSolPolyDerivInSolPnts();
  m_eGammaFx.assign(frLocalData[0]->getNbrOfSolPnts(), 0.);
  m_eGammaFy.assign(frLocalData[0]->getNbrOfSolPnts(), 0.);
  m_eGammaDivergence.assign(frLocalData[0]->getNbrOfSolPnts(), 0.);
  m_unpertStates.resize(frLocalData[0]->getNbrOfSolPnts());
  for (CFuint i = 0; i < m_unpertStates.size(); ++i) m_unpertStates[i].resize(m_nbrEqs);
}

//////////////////////////////////////////////////////////////////////////////////////////////////

void GammaAlpha2DSourceTerm::execute()
{
  if (m_onsetFunction)
  {
    // the boundary pass of this residual evaluation has filled the local wall table; gather it
    const std::map< std::pair<CFreal,CFreal>, CFreal >& table = gammaAlphaOnsetTable();
    const std::map< std::pair<CFreal,CFreal>, std::pair<CFreal,CFreal> >& points = gammaAlphaWallPointTable();
    std::vector<CFreal> local;
    local.reserve(m_onsetStride*table.size());
    for (std::map< std::pair<CFreal,CFreal>, CFreal >::const_iterator it = table.begin(); it != table.end(); ++it)
    {
      local.push_back(it->first.first);
      local.push_back(it->first.second);
      local.push_back(it->second);
      if (m_wallLaw)
      {
        const std::map< std::pair<CFreal,CFreal>, std::pair<CFreal,CFreal> >::const_iterator pt = points.find(it->first);
        local.push_back(pt != points.end() ? pt->second.first : 0.);
        local.push_back(pt != points.end() ? pt->second.second : 0.);
      }
    }
#ifdef CF_HAVE_MPI
    const std::string nsp = getMethodData().getNamespace();
    MPI_Comm comm = Common::PE::GetPE().GetCommunicator(nsp);
    const int nbProc = Common::PE::GetPE().GetProcessorCount(nsp);
    int localSize = static_cast<int>(local.size());
    std::vector<int> sizes(nbProc), offsets(nbProc,0);
    MPI_Allgather(&localSize, 1, MPI_INT, &sizes[0], 1, MPI_INT, comm);
    for (int i = 1; i < nbProc; ++i) offsets[i] = offsets[i-1] + sizes[i-1];
    m_onsetPoints.resize(offsets[nbProc-1] + sizes[nbProc-1]);
    MPI_Allgatherv(local.empty() ? CFNULL : &local[0], localSize, MPI_DOUBLE,
                   m_onsetPoints.empty() ? CFNULL : &m_onsetPoints[0], &sizes[0], &offsets[0], MPI_DOUBLE, comm);
#else
    m_onsetPoints = local;
#endif
    ++m_onsetStamp;
    if (m_wallLaw) marchWallLaw();
    if (m_onsetStamp % 50 == 1)
    {
      // diagnostic: wall points with tau_w* <= alpha and their most upstream x
      CFuint nOn = 0;
      CFreal xOn = std::numeric_limits<CFreal>::max();
      for (CFuint i = 0; i + 2 < m_onsetPoints.size(); i += m_onsetStride)
      {
        if (m_onsetPoints[i+2] > 0.0) { ++nOn; xOn = std::min(xOn, m_onsetPoints[i]); }
      }
      CFLog(INFO, "GammaAlpha onset: evaluation " << m_onsetStamp << ", " << nOn << " of " << m_onsetPoints.size()/m_onsetStride
            << " wall points with tau_w* <= alpha, most upstream x = " << (nOn > 0 ? xOn : -1.0) << "\n");
      if (m_wallLaw)
      {
        // diagnostic: the wall law at a few places along the wall
        std::ostringstream law;
        law << "GammaAlpha wall law: start x = " << m_wallLawStart << "; gamma_w";
        const CFreal stations[8] = {0.25, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0, 5.0};
        for (CFuint j = 0; j < 8; ++j)
        {
          CFreal best = std::numeric_limits<CFreal>::max();
          CFuint kBest = 0;
          for (CFuint k = 0; k < m_wallGamma.size(); ++k)
          {
            const CFreal d = std::fabs(m_onsetPoints[k*m_onsetStride] - stations[j]);
            if (d < best) { best = d; kBest = k; }
          }
          if (!m_wallGamma.empty()) law << " " << m_onsetPoints[kBest*m_onsetStride] << ":" << m_wallGamma[kBest];
        }
        CFLog(INFO, law.str() << "\n");
      }
    }
  }

  if (m_eGammaArtificialDiffusion > 0.0) updateGammaArtificialDiffusion();

  KLogOmega2DSourceTerm::execute();
}

//////////////////////////////////////////////////////////////////////////////////////////////////

void GammaAlpha2DSourceTerm::updateGammaArtificialDiffusion()
{
  // C rho |w_E| dn/2 at every state: rho w_E of the E_gamma advection (GammaAlphaEGammaFlux.hh, C_Eg = 0.25, along
  // n = (-v, u)/c), dn = 1/invDn the solution-point spacing along n, invDn = (p+1)/2 sum_d |grad xi_d . n|
  // (as EGammaImplicitDiagonal); all states, so that both cells of a partition face see the same values
  DataHandle< State*, GLOBAL > states = socket_states.getDataHandle();
  DataHandle< CFreal > volumes = socket_volumes.getDataHandle();
  DataHandle< CFreal > solPntNormals = socket_solPntNormals.getDataHandle();
  Physics::GammaAlpha::GammaArtificialDiffusion& artDiff = Physics::GammaAlpha::gammaArtificialDiffusion();
  artDiff.active = true;
  artDiff.perState.assign(states.size(), 0.0);
  const CFreal R = m_eulerVarSet->getModel()->getR();
  const CFreal kappa = m_eulerVarSet->getModel()->getGamma();
  const CFreal machInf = m_eulerVarSet->getModel()->getMachInf();
  const CFreal uInf = m_eulerVarSet->getModel()->getVelInf();
  const CFreal pInf = m_eulerVarSet->getModel()->getStaticPressInf();
  for (CFuint iState = 0; iState < states.size(); ++iState)
  {
    const State& state = *states[iState];
    const CFuint sID = state.getLocalID();
    const CFuint gammaID = state.size() - 2; // p, u, v, T, k, omega, gamma, alpha
    const CFreal p = state[0];
    const CFreal u = state[1];
    const CFreal v = state[2];
    const CFreal T = state[3];
    const CFreal c = std::sqrt(u*u + v*v);
    if (!(c > 0.0) || !(p > 0.0) || !(T > 0.0) || !(volumes[sID] > 0.0)) continue;
    CFreal fx = 0., fy = 0.;
    Physics::GammaAlpha::eGammaMassFlux(0.25, p/(R*T), u, v, p, state[gammaID], R, kappa, machInf, uInf, pInf, fx, fy);
    const CFreal nx = -v/c;
    const CFreal ny = u/c;
    CFreal invDn = 0.0;
    for (CFuint iDir = 0; iDir < 2; ++iDir)
    {
      invDn += std::fabs(solPntNormals[sID*4+iDir*2]*nx + solPntNormals[sID*4+iDir*2+1]*ny);
    }
    invDn *= 0.5*(m_order + 1.0)/volumes[sID];
    if (!(invDn > 0.0)) continue;
    CFreal weight = 1.0;
    if (m_eGammaADWallGammaStart >= 0.0)
    {
      weight = std::min(1.0, std::max(0.0, (wallLawCap(state) - m_eGammaADWallGammaStart)/(m_eGammaADWallGammaEnd - m_eGammaADWallGammaStart)));
      if (weight == 0.0) continue;
    }
    artDiff.perState[sID] = weight*m_eGammaArtificialDiffusion*std::sqrt(fx*fx + fy*fy)/(2.0*invDn);
  }
}

//////////////////////////////////////////////////////////////////////////////////////////////////

CFreal GammaAlpha2DSourceTerm::onsetFunction(const State& state)
{
  // no wall data yet (first evaluation): no onset control
  if (m_onsetPoints.empty()) return 1.0;

  const CFuint id = state.getLocalID();
  if (id >= m_onsetCache.size())
  {
    m_onsetCache.resize(id+1,0.0);
    m_onsetCacheStamp.resize(id+1,0);
    m_wallCapCache.resize(id+1,0.01);
  }
  if (m_onsetCacheStamp[id] != m_onsetStamp)
  {
    const CFreal x = state.getCoordinates()[XX];
    const CFreal y = state.getCoordinates()[YY];
    CFreal best = std::numeric_limits<CFreal>::max();
    CFreal margin = 0.0;
    CFuint nearest = 0;
    for (CFuint i = 0; i + 2 < m_onsetPoints.size(); i += m_onsetStride)
    {
      const CFreal dx = m_onsetPoints[i] - x;
      const CFreal dy = m_onsetPoints[i+1] - y;
      const CFreal d2 = dx*dx + dy*dy;
      if (d2 < best)
      {
        best = d2;
        margin = m_onsetPoints[i+2];
        nearest = i/m_onsetStride;
      }
    }
    // WallLaw: production behind the place where the march starts (the turbulent wall shear behind the transition
    // is above alpha again, so the local criterion would switch it off there)
    m_onsetCache[id] = (m_wallLaw && nearest < m_wallBehind.size()) ? m_wallBehind[nearest]
                                                                     : 0.5*(1.0 + std::tanh(margin/m_onsetWidth));
    // WallLaw: the wall-law value at the projection of the state onto the segment from the nearest wall point to
    // one of its neighbours along the wall (the nearest point's value if the projection falls on neither)
    CFreal cap = 0.01;
    if (m_wallLaw && nearest < m_wallGamma.size() && nearest < m_wallRank.size())
    {
      cap = m_wallGamma[nearest];
      const CFreal xk = m_onsetPoints[nearest*m_onsetStride];
      const CFreal yk = m_onsetPoints[nearest*m_onsetStride+1];
      const CFuint place = m_wallRank[nearest];
      for (int side = 1; side >= -1; side -= 2)
      {
        if ((side < 0 && place == 0) || (side > 0 && place + 1 >= m_wallOrder.size())) continue;
        const CFuint kn = m_wallOrder[side > 0 ? place + 1 : place - 1];
        const CFreal ex = m_onsetPoints[kn*m_onsetStride] - xk;
        const CFreal ey = m_onsetPoints[kn*m_onsetStride+1] - yk;
        const CFreal len2 = ex*ex + ey*ey;
        if (!(len2 > 0.0)) continue;
        const CFreal t = ((x - xk)*ex + (y - yk)*ey)/len2;
        if (t > 0.0 && t <= 1.0)
        {
          cap = m_wallGamma[nearest] + t*(m_wallGamma[kn] - m_wallGamma[nearest]);
          break;
        }
      }
    }
    m_wallCapCache[id] = cap;
    m_onsetCacheStamp[id] = m_onsetStamp;
  }
  return m_onsetCache[id];
}

//////////////////////////////////////////////////////////////////////////////////////////////////

CFreal GammaAlpha2DSourceTerm::wallLawCap(const State& state)
{
  if (m_onsetPoints.empty()) return 0.01;
  onsetFunction(state);
  return m_wallCapCache[state.getLocalID()];
}

//////////////////////////////////////////////////////////////////////////////////////////////////

namespace {
// Spot growth factor f_sigma of Karsch et al. (Jacobs et al. 2025, eqs. 11-16): edge Mach number Me, edge temperature
// Te, wall temperature Tw, Re_theta,cr, ratio of specific heats kappa, recovery factor r
inline CFreal karschSpotFactor(const CFreal Me, const CFreal Te, const CFreal Tw, const CFreal reThetaCr,
                               const CFreal kappa, const CFreal r)
{
  const CFreal e = std::exp(-0.004*reThetaCr);
  const CFreal cleOverCe = 0.15*e + 0.85;                                               // eq. 12
  const CFreal cteOverCe = 0.61*e + 0.39;                                               // eq. 13
  const CFreal tjetOverTe = 0.45 + 0.55*Tw/Te + 0.25*r*0.5*(kappa - 1.0)*Me*Me;         // eq. 16
  const CFreal mc = Me*(1.0 - 0.45)/(1.0 + std::sqrt(std::max(tjetOverTe, 1.0e-12)));  // eq. 15, c_jet = 0.45 c_e
  const CFreal betaStar = (10.0*M_PI/180.0)*std::pow(1.0 + 7.06*std::pow(mc, 2.86), -0.5)
                          *(1.0 - 0.29*std::exp(-0.0035*reThetaCr));                    // eq. 14, beta0* = 10 deg
  return (1.0/0.25)*(1.0/cteOverCe - 1.0/cleOverCe)*std::tan(betaStar);                // eq. 11, sigma0 = 0.25
}

// Ramp f_gamma of Eq. 3 and the rate of the gamma equation at the wall (Jatayu v1.5 wallLawRate, f_K = 1)
inline CFreal wallLawRate(const CFreal gamma, const CFreal beta0)
{
  const CFreal g = std::min(std::max(gamma, 0.01), 0.99);
  const CFreal fg = (g < 0.45) ? 1.0 - std::exp(-1.735*std::tan(5.45*g - 0.95375) - 2.2) : 1.0;
  return 2.0*fg*(1.0 - g)*std::sqrt(-std::log(1.0 - g))*beta0;
}
// gamma at the end of a stretch ds that starts with gammaIn, beta0 linear from betaA to betaB over it: RK4 in equal
// steps with at most 1e-3 change of gamma per step (Jatayu v1.5 wallLawAdvance); gamma kept in [0.01, 0.99]
inline CFreal wallLawAdvance(const CFreal gammaIn, const CFreal betaA, const CFreal betaB, const CFreal ds)
{
  CFreal g = std::min(std::max(gammaIn, 0.01), 0.99);
  if (!(ds > 0.0) || !(std::max(betaA, betaB) > 0.0)) return g;
  const CFreal fastest = 0.86*std::max(betaA, betaB);
  const int steps = static_cast<int>(std::min(std::ceil(fastest*ds/1.0e-3) + 1.0, 200000.0));
  const CFreal h = ds/steps;
  for (int i = 0; i < steps && g < 0.99; ++i)
  {
    const CFreal t0 = CFreal(i)/steps, tm = (i + 0.5)/steps, t1 = CFreal(i + 1)/steps;
    const CFreal k1 = wallLawRate(g, betaA + t0*(betaB - betaA));
    const CFreal k2 = wallLawRate(g + 0.5*h*k1, betaA + tm*(betaB - betaA));
    const CFreal k3 = wallLawRate(g + 0.5*h*k2, betaA + tm*(betaB - betaA));
    const CFreal k4 = wallLawRate(g + h*k3, betaA + t1*(betaB - betaA));
    g = std::min(g + h*(k1 + 2.0*k2 + 2.0*k3 + k4)/6.0, 0.99);
  }
  return g;
}
} // namespace

void GammaAlpha2DSourceTerm::marchWallLaw()
{
  const CFuint stride = m_onsetStride;
  const CFuint n = m_onsetPoints.size()/stride;
  m_wallBehind.assign(n, 0.0);
  m_wallGamma.assign(n, 0.01);
  m_wallLawStart = -1.0;
  if (n == 0) return;

  // the wall points in their order along the wall (one flat wall: increasing x, then y)
  std::vector<CFuint> order(n);
  for (CFuint k = 0; k < n; ++k) order[k] = k;
  const std::vector<CFreal>& pts = m_onsetPoints;
  std::sort(order.begin(), order.end(), [&pts, stride](const CFuint a, const CFuint b) {
    return pts[a*stride] < pts[b*stride] || (pts[a*stride] == pts[b*stride] && pts[a*stride+1] < pts[b*stride+1]);
  });
  m_wallOrder = order;
  m_wallRank.assign(n, 0);
  for (CFuint idx = 0; idx < n; ++idx) m_wallRank[order[idx]] = idx;

  // beta0 of each wall point as beta of the cells (Eqs. 4-5 without f_K): Tu_inf,loc from its alpha
  SafePtr< NavierStokes2DGammaAlphaPuvt > nsVarSet = m_diffVarSet.d_castTo< NavierStokes2DGammaAlphaPuvt >();
  const CFreal Rgas = m_eulerVarSet->getModel()->getR();
  const CFreal kappa = m_eulerVarSet->getModel()->getGamma();
  const CFreal machInf = m_eulerVarSet->getModel()->getMachInf();
  const CFreal uInf = m_eulerVarSet->getModel()->getVelInf();
  const CFreal pInf = m_eulerVarSet->getModel()->getStaticPressInf();
  RealVector fsState(m_nbrEqs);
  fsState = 0.0;
  std::vector<CFreal> beta0(n, 0.0);
  for (CFuint k = 0; k < n; ++k)
  {
    const CFreal alphaW = std::max(pts[k*stride+3], 1.0e-12);
    const CFreal pW = std::max(pts[k*stride+4], 1.0e-3);
    const Physics::NavierStokes::transition_detail::GammaAlphaFreestream fs =
      Physics::NavierStokes::transition_detail::gammaAlphaFreestreamInlet(pW, Rgas, kappa, machInf, uInf, pInf);
    fsState[0] = pW;
    fsState[3] = fs.temperature;
    const CFreal muInf = nsVarSet->getLaminarDynViscosityFromGradientVars(fsState);
    const CFreal reTheta = fs.density*fs.speed*fs.speed*0.2247141166/(alphaW*std::sqrt(fs.density*muInf));
    const CFreal tuLoc = getTuInfLocal(reTheta, fs.mach, 2.0);
    const CFreal fMn = m_karschSpotGrowth ? karschSpotFactor(fs.mach, fs.temperature, m_karschWallTemperature, reTheta, kappa, m_karschRecoveryFactor)
                                          : 1.0/std::pow(1.0 + 0.58*std::pow(fs.mach, 0.6), 2.0);
    beta0[k] = std::sqrt(1.25e-11*std::pow(tuLoc, 1.75)*fMn)*fs.speed*fs.density/muInf;
  }

  // the march: 0.01 up to the place where the criterion (margin > 0) is first met, interpolated between the last
  // point without it and the first with it; from there each point from the one before (Jatayu v1.5 wallLawMarch)
  bool started = false;
  for (CFuint idx = 0; idx < n; ++idx)
  {
    const CFuint k = order[idx];
    // laminar stretch (WallLawLaminarUntilX): 0.01, no production, and the law starts anew behind it
    if (pts[k*stride] <= m_wallLawLaminarUntilX)
    {
      m_wallGamma[k] = 0.01;
      m_wallBehind[k] = 0.0;
      started = false;
      continue;
    }
    if (!started)
    {
      if (pts[k*stride+2] > 0.0)
      {
        started = true;
        m_wallBehind[k] = 1.0;
        if (idx > 0 && pts[order[idx-1]*stride] > m_wallLawLaminarUntilX)
        {
          const CFuint kp = order[idx-1];
          const CFreal mp = pts[kp*stride+2];
          const CFreal mk = pts[k*stride+2];
          const CFreal theta = (mp < 0.0) ? -mp/(mk - mp) : 0.0;
          const CFreal ds = std::sqrt(std::pow(pts[k*stride] - pts[kp*stride], 2) + std::pow(pts[k*stride+1] - pts[kp*stride+1], 2));
          m_wallLawStart = pts[kp*stride] + theta*(pts[k*stride] - pts[kp*stride]);
          m_wallGamma[k] = wallLawAdvance(0.01, beta0[kp] + theta*(beta0[k] - beta0[kp]), beta0[k], (1.0 - theta)*ds);
        }
        else
        {
          m_wallLawStart = pts[k*stride];
        }
      }
    }
    else
    {
      const CFuint kp = order[idx-1];
      const CFreal ds = std::sqrt(std::pow(pts[k*stride] - pts[kp*stride], 2) + std::pow(pts[k*stride+1] - pts[kp*stride+1], 2));
      m_wallGamma[k] = wallLawAdvance(m_wallGamma[kp], beta0[kp], beta0[k], ds);
      m_wallBehind[k] = 1.0;
    }
  }

  // the wall values for the wall condition (read in the next residual evaluation)
  std::map< std::pair<CFreal,CFreal>, CFreal >& law = gammaAlphaWallLawTable();
  for (CFuint k = 0; k < n; ++k)
  {
    law[std::make_pair(pts[k*stride], pts[k*stride+1])] = m_wallGamma[k];
  }
}

//////////////////////////////////////////////////////////////////////////////////////////////////

void GammaAlpha2DSourceTerm::unsetup()
{
  KLogOmega2DSourceTerm::unsetup();
}

//////////////////////////////////////////////////////////////////////////////

std::vector< Common::SafePtr< BaseDataSocketSource > >
  GammaAlpha2DSourceTerm::providesSockets()
{
  std::vector< Common::SafePtr< BaseDataSocketSource > > result = KLogOmega2DSourceTerm::providesSockets();
  result.push_back(&socket_MInfLocal);
  result.push_back(&socket_uInfLocal);
  result.push_back(&socket_TInfLocal);
  result.push_back(&socket_TuInfLocal);
  result.push_back(&socket_alphaDiff);
  return result;
}

//////////////////////////////////////////////////////////////////////////////

std::vector< Common::SafePtr< BaseDataSocketSink > >
GammaAlpha2DSourceTerm::needsSockets()
{
  std::vector< Common::SafePtr< BaseDataSocketSink > > result = KLogOmega2DSourceTerm::needsSockets();
  result.push_back(&socket_updateCoeff);
  result.push_back(&socket_solPntNormals);
  result.push_back(&socket_states);
  result.push_back(&socket_volumes);
  return result;
}

//////////////////////////////////////////////////////////////////////////////////////////////////

void  GammaAlpha2DSourceTerm::getStrain(const CFuint iState)
{
  const CFreal gradU_X = (*(m_cellGrads[iState][1]))[XX];
  const CFreal gradU_Y = (*(m_cellGrads[iState][1]))[YY];
  const CFreal gradV_X = (*(m_cellGrads[iState][2]))[XX];
  const CFreal gradV_Y = (*(m_cellGrads[iState][2]))[YY];
  const CFreal gradSum = (gradU_Y + gradV_X);
  const CFreal strain = gradU_X*gradU_X + 0.5*gradSum*gradSum + gradV_Y*gradV_Y + m_vOverRadius*m_vOverRadius ;
  m_strain = std::sqrt(2.*strain);
}

///////////////////////////////////////////////////////////////////////////////////////////////////

void  GammaAlpha2DSourceTerm::getVorticity(const CFuint iState)
{
  const CFreal Vorticity1 = (*(m_cellGrads[iState][2]))[XX] - (*(m_cellGrads[iState][1]))[YY];
  m_vorticity = fabs(Vorticity1);
}

//////////////////////////////////////////////////////////////////////////////

void GammaAlpha2DSourceTerm::addSourceTerm(RealVector& resUpdates)
{       
  std::vector<Physics::NavierStokes::transition_detail::GammaAlphaFreestream> freestreamStates(m_cellStates->size());
  // Preflight the entire cell: a bad later point must not leave partial RHS writes.
  for (CFuint i = 0; i < m_cellStates->size(); ++i) {
    const State& state = *(*m_cellStates)[i];
    detail::validateTransitionSourceState(state, 2, state.getLocalID(), !m_signedK);
    Physics::NavierStokes::transition_detail::validateThermodynamics(state[0], state[3], m_eulerVarSet->getModel()->getR(), m_eulerVarSet->getModel()->getGamma());
    freestreamStates[i] = (m_modelVersion == 2025 ? Physics::NavierStokes::transition_detail::gammaAlphaFreestreamInlet : Physics::NavierStokes::transition_detail::gammaAlphaFreestream)(std::max(state[0],1.0e-3), m_eulerVarSet->getModel()->getR(), m_eulerVarSet->getModel()->getGamma(), m_eulerVarSet->getModel()->getMachInf(), m_eulerVarSet->getModel()->getVelInf(), m_eulerVarSet->getModel()->getStaticPressInf());
    detail::validateTransitionSourceEnvironment(socket_gradients.getDataHandle()[state.getLocalID()], socket_wallDistance.getDataHandle()[state.getLocalID()], 2, m_nbrEqs);
    CFreal velocityGradients[3][3] = {};
    for (CFuint e=0; e<2; ++e) for (CFuint d=0; d<2; ++d)
      velocityGradients[e][d] = socket_gradients.getDataHandle()[state.getLocalID()][e+1][d];
    const CFreal axisGradient = m_isAxisymmetric ? Physics::NavierStokes::transition_detail::finiteProduct(1./std::max(state.getCoordinates()[YY],1.e-4),state[2]) : 0.;
    Physics::NavierStokes::transition_detail::validateVelocityGradientQuadratics(velocityGradients,2,axisGradient);
    auto transport = m_diffVarSet.d_castTo<NavierStokes2DGammaAlphaPuvt>();
    transport->getLaminarDynViscosityFromGradientVars(state);
    Physics::NavierStokes::transition_detail::laminarViscosity(transport->getModel(), state[0], freestreamStates[i].temperature,
        m_eulerVarSet->getModel()->getPressRef(), m_eulerVarSet->getModel()->getTempRef());
  }

  const CFuint kID = (*m_cellStates)[0]->size() - 4;
  const CFuint omegaID = kID + 1;
  const CFuint gammaID = kID + 2;
  const CFuint alphaID = kID + 3;

  const CFuint iKPD = m_eulerVarSet->getModel()->getFirstScalarVar(0);
      
  SafePtr< NavierStokes2DGammaAlphaPuvt > navierStokesVarSet = m_diffVarSet.d_castTo< NavierStokes2DGammaAlphaPuvt >();
   
  // get the gradients datahandle
  DataHandle< vector< RealVector > > gradients = socket_gradients.getDataHandle();

  // set gradients
  const CFuint nbrStates = m_cellStates->size();

  for (CFuint iState = 0; iState < nbrStates; ++iState)
  {
    const CFuint stateID = (*m_cellStates)[iState]->getLocalID();
    
    for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
    {
      *(m_cellGrads[iState][iEq]) = gradients[stateID][iEq];
    }
    
    // Get the wall distance
    DataHandle< CFreal > wallDist = socket_wallDistance.getDataHandle();

    m_currWallDist[iState] = wallDist[stateID];//((*m_cellStates)[iState]->getCoordinates())[YY];
  }

  if (m_eGammaImplicitDiagonal && !m_gradientJacobian && !m_isPerturbed)
  {
    // unperturbed evaluation: remember the states (EGammaImplicitDiagonal measures the perturbation of gamma)
    for (CFuint iState = 0; iState < nbrStates; ++iState)
    {
      for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq) m_unpertStates[iState][iEq] = (*(*m_cellStates)[iState])[iEq];
    }
  }

  if (m_gradientJacobian)
  {
    if (!m_isPerturbed)
    {
      // unperturbed evaluation: remember the states the stored gradients belong to
      for (CFuint iState = 0; iState < nbrStates; ++iState)
      {
        for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq) m_unpertStates[iState][iEq] = (*(*m_cellStates)[iState])[iEq];
      }
    }
    else
    {
      // perturbed evaluation: grad q_e += sum_d D[e][d][p] proj_d(e) dq_p / J_e (FR volume term, gradient vars = states)
      DataHandle< CFreal > solPntNormals = socket_solPntNormals.getDataHandle();
      for (CFuint pSol = 0; pSol < nbrStates; ++pSol)
      {
        for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
        {
          const CFreal dq = (*(*m_cellStates)[pSol])[iEq] - m_unpertStates[pSol][iEq];
          if (dq == 0.0) continue;
          for (CFuint eSol = 0; eSol < nbrStates; ++eSol)
          {
            const CFuint eID = (*m_cellStates)[eSol]->getLocalID();
            for (CFuint iDir = 0; iDir < 2; ++iDir)
            {
              const CFreal coef = (*m_solPolyDerivAtSolPnts)[eSol][iDir][pSol]*dq/m_solPntJacobDets[eSol];
              if (coef == 0.0) continue;
              for (CFuint jDim = 0; jDim < 2; ++jDim)
              {
                (*(m_cellGrads[eSol][iEq]))[jDim] += coef*solPntNormals[eID*4+iDir*2+jDim];
              }
            }
          }
        }
      }
    }
  }
  
  if (m_eGammaFlux)
  {
    // rho w_E of E_gamma at the solution points, and its divergence with the FR volume term of the cell
    // (the transformed flux rho w_E . proj_d at each point, differentiated as the convective flux is)
    const CFreal Rgas = m_eulerVarSet->getModel()->getR();
    const CFreal kappa = m_eulerVarSet->getModel()->getGamma();
    const CFreal machInf = m_eulerVarSet->getModel()->getMachInf();
    const CFreal uInf = m_eulerVarSet->getModel()->getVelInf();
    const CFreal pInf = m_eulerVarSet->getModel()->getStaticPressInf();
    for (CFuint jSol = 0; jSol < nbrStates; ++jSol)
    {
      const State& sj = *(*m_cellStates)[jSol];
      Physics::GammaAlpha::eGammaMassFlux(0.25, sj[0]/(Rgas*sj[3]), sj[1], sj[2], sj[0], sj[6], Rgas, kappa, machInf,
                                          uInf, pInf, m_eGammaFx[jSol], m_eGammaFy[jSol]);
    }
    DataHandle< CFreal > solPntNormals = socket_solPntNormals.getDataHandle();
    for (CFuint eSol = 0; eSol < nbrStates; ++eSol)
    {
      CFreal divergence = 0.;
      for (CFuint iDir = 0; iDir < 2; ++iDir)
      {
        for (CFuint jSol = 0; jSol < nbrStates; ++jSol)
        {
          const CFreal coef = (*m_solPolyDerivAtSolPnts)[eSol][iDir][jSol];
          if (coef == 0.) continue;
          const CFuint jID = (*m_cellStates)[jSol]->getLocalID();
          divergence += coef*(m_eGammaFx[jSol]*solPntNormals[jID*4+iDir*2] + m_eGammaFy[jSol]*solPntNormals[jID*4+iDir*2+1]);
        }
      }
      m_eGammaDivergence[eSol] = divergence/m_solPntJacobDets[eSol];
    }
  }

  const EquationSubSysDescriptor& eqData = PhysicalModelStack::getActive()->getEquationSubSysDescriptor();
  const CFuint nbEqs = eqData.getNbEqsSS();
  const CFuint totalNbEqs = PhysicalModelStack::getActive()->getNbEq(); 
  
  for (CFuint iSol = 0; iSol < nbrStates; ++iSol)
  {
    const State& sourceState = *(*m_cellStates)[iSol];
    const CFreal avV = detail::transitionResidualSpeed(sourceState[1], sourceState[2], 0.0, m_PGrad);
    m_eulerVarSet->computePhysicalData(*((*m_cellStates)[iSol]), m_solPhysData);
    // The physical-data norm squares components; retain tiny positive speeds.
    m_solPhysData[EulerTerm::V] = avV;
      
    // Set the wall distance before computing the turbulent viscosity
    navierStokesVarSet->setWallDistance(m_currWallDist[iSol]);
    
    const CFreal mut = navierStokesVarSet->getTurbDynViscosityFromGradientVars(*((*m_cellStates)[iSol]), m_cellGrads[iSol]);
    const CFreal rhoKOverMut = m_diffVarSet.d_castTo< NavierStokes2DGammaAlphaPuvt >()->getRhoKOverMuT(*((*m_cellStates)[iSol]), m_cellGrads[iSol]);
    const CFreal mu = navierStokesVarSet->getLaminarDynViscosityFromGradientVars(*((*m_cellStates)[iSol]));
    
    navierStokesVarSet->computeBlendingCoefFromGradientVars(*((*m_cellStates)[iSol]), *(m_cellGrads[iSol][kID]), *(m_cellGrads[iSol][omegaID]));
    
    if (m_isAxisymmetric)
    {
      m_overRadius = 1.0/max(((*m_cellStates)[iSol]->getCoordinates())[YY],1.0e-4);
      m_vOverRadius = m_overRadius*m_solPhysData[EulerTerm::VY];
    }
    
    // Get Vorticity
    getVorticity(iSol);
    
    // get current state values
    const CFreal avP     = std::max(m_solPhysData[EulerTerm::P],1.0e-3);
    const CFreal avK     = std::max(m_solPhysData[iKPD],0.0);
    const CFreal avOmega = std::exp(m_solPhysData[iKPD+1]);
    
    const CFreal avGa    = std::min(std::max(m_solPhysData[iKPD+2],0.01),0.99);
    // gamma as it is, used where the (1-gamma) factor has to vanish at gamma = 1; avGa stays
    // the guarded value for log(1-gamma), which is what the clamp is there for
    const CFreal gaRaw   = m_gammaEquilibrium ? m_solPhysData[iKPD+2] : avGa;
    //const CFreal avGa    = std::min(std::max(m_solPhysData[iKPD+2],0.0),1.0);

    const CFreal rho = navierStokesVarSet->getDensity(*((*m_cellStates)[iSol]));
    const CFreal u = m_solPhysData[EulerTerm::VX];
    const CFreal v = m_solPhysData[EulerTerm::VY];
    const CFreal dpx = (*(m_cellGrads[iSol][0]))[XX];
    const CFreal dpy = (*(m_cellGrads[iSol][0]))[YY];
    
    const auto& freestream = freestreamStates[iSol];
    const CFreal gammaIsentropic = m_eulerVarSet->getModel()->getGamma();
    const CFreal R = m_eulerVarSet->getModel()->getR();
    const CFreal rhoInfLocal = freestream.density;
    const CFreal MInfLocal = freestream.mach;
    const CFreal uInfLocal = freestream.speed;
    const CFreal TInfLocal = freestream.temperature;
    
    const CFreal Tu = detail::transitionIntensity(avK, avV, 13.0);
    
    if (!(Tu >= 0.0)) CFLog(INFO, "Tu problem, k: " << avK << ", V: " << avV << "\n");
    
    // Owned temporary also preserves the state if viscosity evaluation throws.
    RealVector freestreamState((*m_cellStates)[iSol]->size());
    freestreamState = *((*m_cellStates)[iSol]);
    freestreamState[3] = TInfLocal;
    const CFreal muInfLocal = navierStokesVarSet->getLaminarDynViscosityFromGradientVars(freestreamState);
    
    const CFreal alphaMin = rhoInfLocal*uInfLocal*uInfLocal*0.2247/(std::sqrt(rhoInfLocal*muInfLocal*(1.0+0.38*pow(MInfLocal,0.6)))*4900.0);
    
    const CFreal avAlpha = std::max(m_solPhysData[iKPD+3],alphaMin);    
    
    CFreal ReThetat;

    if (!m_PGrad)
    {
      ReThetat = rhoInfLocal*uInfLocal*uInfLocal*0.2247141166/(avAlpha*sqrt(rhoInfLocal*muInfLocal));
    }
    else
    {
      ReThetat = getRethetatwithPressureGradient(avAlpha,rhoInfLocal,uInfLocal,muInfLocal,mu,rho,dpx,dpy);
    }

    const CFreal TuInfLocal = getTuInfLocal(ReThetat,MInfLocal,Tu);

    // F_onset (2025 eq. A7, smoothed) from the nearest wall flux point; 1 without OnsetFunction
    const CFreal fOnset = m_onsetFunction ? onsetFunction(*(*m_cellStates)[iSol]) : 1.0;
    
    // compute destruction terms for k and logOmega
    // Destruction term: k
    m_destructionTerm_k = (-1.) * rho * avOmega * avK * navierStokesVarSet->getBetaStar(*((*m_cellStates)[iSol]));
      
    // Destruction term: Omega
    m_destructionTerm_Omega = (-1.) * rho * avOmega * navierStokesVarSet->getBeta(*((*m_cellStates)[iSol]));//(-1.) * rho * avOmega * avOmega * navierStokesVarSet->getBeta(*((*m_cellStates)[iState]));

  
    // Make sure negative values dont propagate...
    m_destructionTerm_k     = std::min(0., m_destructionTerm_k );
    m_destructionTerm_Omega = std::min(0., m_destructionTerm_Omega);
        
    ///The Modified Production  terms for k and logOmega
    const CFuint uID = 1;//getStateVelocityIDs()[XX];
    const CFuint vID = 2;//getStateVelocityIDs()[YY];
    const CFreal dux = (*(m_cellGrads[iSol][uID]))[XX];
    const CFreal duy = (*(m_cellGrads[iSol][uID]))[YY]; 
    const CFreal dvx = (*(m_cellGrads[iSol][vID]))[XX]; 
    const CFreal dvy = (*(m_cellGrads[iSol][vID]))[YY]; 
    const CFreal dgammax = (*(m_cellGrads[iSol][6]))[XX];
    const CFreal dgammay = (*(m_cellGrads[iSol][6]))[YY];
    const CFreal dax = (*(m_cellGrads[iSol][7]))[XX];
    const CFreal day = (*(m_cellGrads[iSol][7]))[YY];
    const CFreal dkx = (*(m_cellGrads[iSol][4]))[XX];
    const CFreal dky = (*(m_cellGrads[iSol][4]))[YY];

    const CFreal coeffTauMu = navierStokesVarSet->getModel().getCoeffTau();
    // KProductionSST: SST's own k production, with mu_t = mu_t*/f(gamma) (the A4 factor of the variable set) and
    // without the k_t* factor; otherwise the effective form (mu_t* and k_t*)
    const CFreal muTStarFactor = Physics::GammaAlpha::gammaAlphaMeanFlowMuTFactor(avGa);
    const CFreal mutK = m_kProductionSST ? mut/muTStarFactor : mut;
    const CFreal gammaTerm = m_kProductionSST ? 1.0 : avGa+avGa*(1.0-avGa);
    const CFreal rhoKOverMutK = m_kProductionSST ? rhoKOverMut*muTStarFactor : rhoKOverMut;
    CFreal productionOverMut;
    const CFreal twoThirdRhoK = (2./3.)*(avK * rho * gammaTerm);
    
    if (!m_isSSTV)
    {
      m_prodTerm_k = coeffTauMu*(mutK*((4./3.)*((dux-dvy)*(dux-dvy)+(dux*dvy)-(dux+dvy-m_vOverRadius)*m_vOverRadius)
			         +(duy+dvx)*(duy+dvx)))
                                 -twoThirdRhoK*(dux+dvy+m_vOverRadius);
      productionOverMut = coeffTauMu*(((4./3.)*((dux-dvy)*(dux-dvy)+(dux*dvy)-(dux+dvy-m_vOverRadius)*m_vOverRadius)
			         +(duy+dvx)*(duy+dvx)))
                                 -(2./3.)*rhoKOverMutK*gammaTerm*(dux+dvy+m_vOverRadius);
    }
    else if (m_neglectSSTVTerm)
    {
      const CFreal vorticity2 = (duy-dvx)*(duy-dvx);

      m_prodTerm_k = coeffTauMu*mutK*vorticity2;
      productionOverMut = coeffTauMu*vorticity2;
    }
    else
    {
      const CFreal vorticity2 = (duy-dvx)*(duy-dvx);

      m_prodTerm_k = coeffTauMu*mutK*vorticity2 - twoThirdRhoK*(dux+dvy+m_vOverRadius);
      productionOverMut = coeffTauMu*vorticity2 - (2./3.)*rhoKOverMutK*gammaTerm*(dux+dvy+m_vOverRadius);
    }
  
//    m_prodTerm_k = coeffTauMu*(mut*((4./3.)*((dux-dvy)*(dux-dvy)+(dux*dvy)-(dux+dvy-m_vOverRadius)*m_vOverRadius)
//			         +(duy+dvx)*(duy+dvx)))
//                                 -twoThirdRhoK*(dux+dvy+m_vOverRadius);
  
    ///Production term: Omega
    const CFreal blendingCoefF1 = navierStokesVarSet->getBlendingCoefficientF1();
    const CFreal sigmaOmega2 = navierStokesVarSet->getSigmaOmega2();
     
    const CFreal overOmega = 1./avOmega;
    m_prodTerm_Omega  = navierStokesVarSet->getGammaCoef()*rho * productionOverMut * overOmega;
  
    if (m_limitP)
    {
      m_prodTerm_k     = std::min(10.*fabs(m_destructionTerm_k), m_prodTerm_k);
      // the omega counterpart of the k limit, (gamma rho/mu_t) 10 beta* rho k omega: with KProductionSST the SST mu_t
      // (rhoKOverMutK), as the production itself (v21; before, rho k/mu_t*, i.e. 1/A4 times too loose)
      m_prodTerm_Omega = std::min(10.*rho*rhoKOverMutK*navierStokesVarSet->getBetaStar(*((*m_cellStates)[iSol]))*navierStokesVarSet->getGammaCoef(), m_prodTerm_Omega);
    }  
  
    ///This is used in (BSL,SST), not for normal kOmeg
    m_prodTerm_Omega += (1. - blendingCoefF1) * 2. * rho * overOmega * sigmaOmega2* ((*(m_cellGrads[iSol][kID]))[XX]*(*(m_cellGrads[iSol][omegaID]))[XX] + (*(m_cellGrads[iSol][kID]))[YY]*(*(m_cellGrads[iSol][omegaID]))[YY]);
  //OmegaProdTerm += (1. - blendingCoefF1) * 2. * rho * overOmega * sigmaOmega2* ((*(m_cellGrads[iState][kID]))[XX]*(*(m_cellGrads[iState][omegaID]))[XX] + (*(m_cellGrads[iState][kID]))[YY]*(*(m_cellGrads[iState][omegaID]))[YY]);
    //MathFunctions::innerProd(*(m_cellGrads[iState][kID]), *(m_cellGrads[iState][omegaID]));
//  OmegaProdTerm *= _Radius; 
     
    const CFreal coeffTauMu1 = coeffTauMu*mu;
    // log(omega) form: (mu + sigma_omega mu_t)|grad log omega|^2 with the mu_t of the omega diffusion, i.e. the SST
    // mu_t = mu_t*/f(gamma) with KOmegaDiffusionSST (v21; before, mu_t* although the diffusion had the SST mu_t)
    const CFreal mutOmega = m_kOmegaDiffusionSST ? mut/muTStarFactor : mut;
    const CFreal coeffTauMu3 = coeffTauMu*mutOmega*navierStokesVarSet->getSigmaOmega();//(blendingCoefF1 * navierStokesVarSet->getSigmaOmega1() + (1.0-blendingCoefF1)*sigmaOmega2);
  
    m_prodTerm_Omega += (coeffTauMu1 + coeffTauMu3)*((*(m_cellGrads[iSol][omegaID]))[XX]*(*(m_cellGrads[iSol][omegaID]))[XX] + (*(m_cellGrads[iSol][omegaID]))[YY]*(*(m_cellGrads[iSol][omegaID]))[YY]);
  
    //Make sure negative values dont propagate...
    m_prodTerm_k            = std::max(0., m_prodTerm_k);
    m_prodTerm_Omega        = std::max(0., m_prodTerm_Omega);
    
    // compute production term of gamma
    const CFreal cfg1 = 1.735;
    const CFreal cfg2 = 5.45;
    const CFreal cfg3 = 0.95375;
    const CFreal cfg4 = 2.2;
    
    CFreal fg = 1.0;
    
    if (avGa < 0.45) fg -= exp(-cfg1*tan(cfg2*avGa-cfg3)-cfg4);
    
    const CFreal fMnsigmaTerm = 1.0+0.58*pow(MInfLocal,0.6);
    // KarschSpotGrowth: f_sigma of eqs. 11-16 in place of f_M,n_sigma (eq. 9)
    const CFreal fMnsigma = m_karschSpotGrowth ?
      karschSpotFactor(MInfLocal, TInfLocal, m_karschWallTemperature, ReThetat, m_eulerVarSet->getModel()->getGamma(), m_karschRecoveryFactor)
      : 1.0/(fMnsigmaTerm*fMnsigmaTerm);
    
    CFreal fk = 1.0;
    
    CFreal Rethetac = getRethetat(Tu, false);
    
    ///@todo check if this is local M or local MInf
    Rethetac *= sqrt(1.0+0.38*pow(MInfLocal,0.6));
    
    CFreal lambda = 0.0;

    const CFreal rawDpds = avV > 0.0 ? (u/avV)*dpx + (v/avV)*dpy : std::numeric_limits<CFreal>::quiet_NaN();
    CFreal dpds = rawDpds;
    
    if (m_PGrad)
    {
      const CFreal lambdaTerm1 = -mu/(rho*rho*uInfLocal*uInfLocal*uInfLocal);

      const CFreal dpdsLimit = -m_lambdaLim/(lambdaTerm1*Rethetac*Rethetac);
      dpds = min(max(dpds,-dpdsLimit),dpdsLimit);
    
      const CFreal kPGrad = -muInfLocal/(rhoInfLocal*rhoInfLocal*uInfLocal*uInfLocal*uInfLocal)*fabs(1.0-MInfLocal*MInfLocal)*dpds;
    
      CFreal PRC;
      if (kPGrad < 0)
      {
        const CFreal PRCfactor1 = 474.0*pow(TuInfLocal,-2.9);
        const CFreal PRCexponent = 1.0-exp(2.0e6*kPGrad);
        PRC = pow(PRCfactor1,PRCexponent);
      }
      else
      {
        const CFreal PRCexponent = -3227.0*pow(kPGrad,0.5985);
        PRC = pow(10.0,PRCexponent);
      }
    
      fk += fg*(PRC-1.0); 
      
      lambda = -Rethetac*Rethetac*mu/(rho*rho*uInfLocal*uInfLocal*uInfLocal)*dpds;
    }
    
    const CFreal nsigma = 1.25e-11*pow(TuInfLocal,7.0/4.0)*fk*fMnsigma;///8.0;
    
    const CFreal beta = sqrt(nsigma)*uInfLocal*rhoInfLocal/muInfLocal;
    
    //const CFreal limGa =  std::min(avGa,0.9999);
    CFreal prodTerm_Ga = fOnset*2.0*fg*(1.0-gaRaw)*sqrt(-log(1.0-avGa))*beta*rho*avV;
    
    // compute dissipation term of gamma
    const CFreal cEg = m_ceg;//20.0;//10.0;//20.0;
    const CFreal Cmg2 = 256.0;//1.0e-7;//
  
    const CFreal fMuGamma = 1.0-exp(-Cmg2*(m_currWallDist[iSol]*uInfLocal*rhoInfLocal/muInfLocal)*(m_currWallDist[iSol]*uInfLocal*rhoInfLocal/muInfLocal));
    const CFreal fMMuGamma = (1.0+0.26*(gammaIsentropic-1.0)/2.0*MInfLocal*MInfLocal)*sqrt(1+0.38*pow(MInfLocal,0.6));
  
    //const CFreal gammaLim = std::min(std::max(0.01,avGa),0.99);
    // 2018: thesis eq. 7.65; 2025: Jacobs et al. eq. 31
    const CFreal muGamma = (m_modelVersion == 2025 ? 0.076*pow(avGa/pow(1.0-avGa,7.0),1.0/6.0)*(1.0-pow(avGa,5.0))
                                                   : 0.57*pow(-log(1.0-avGa),-5.0/6.0*(1.0-avGa)))*fMuGamma*fMMuGamma*mu;
    
    // The complete transition gradient/relaxation terms vanish at rest.
    // Do not evaluate an undefined direction and then multiply it by zero.
    CFreal destructionTerm_Ga = 0.0;
    CFreal prodTerm_alpha = 0.0;
    CFreal destructionTerm_alpha = 0.0;
    CFreal dkdnFactor = std::numeric_limits<CFreal>::quiet_NaN();
    if (avV > 0.0)
    {
      const CFreal dudn = -((u/avV)*(u/avV)*duy-(v/avV)*(v/avV)*dvx+(u/avV)*(v/avV)*(dvy-dux));//-duy; //
      const CFreal dgammadn = (v/avV)*dgammax - (u/avV)*dgammay;//-dgammay; //
    
      if (m_modelVersion == 2025) {
        // Jacobs et al. eq. 37, boundary-layer branch: E = C_Eg mu_g dgamma/dn 2c sqrt(1-c/c_e) rho/mu (C_Eg = 0.25),
        // entering as -E. dgammadn above is -dgamma/dn for n pointing away from the wall. The stagnation branch,
        // C_Eg2 c/c_e^2 dc/dn max(-c/|c|.n, 0), vanishes with this flow-aligned n (c.n = 0); it needs a geometric
        // wall normal and matters only in stagnation regions (none on the sharp T3A plate).
        const CFreal blTerm = 2.0*avV*std::sqrt(std::max(0.0,1.0-avV/uInfLocal))*rho/mu;
        // E_gamma is a sink: for dgamma/dn < 0 eq. 37 would turn into an unbounded gamma source (gamma > 1 runaway,
        // vv041). The paper keeps gamma in [0.01, 0.99] and its App. C Jacobian uses max(E_gamma, 0).
        destructionTerm_Ga = std::min(0.0, -0.25*muGamma*(-dgammadn)*blTerm);
        // Floor: the paper clips gamma at 0.01; mu_gamma here is evaluated at max(gamma, 0.01), so the sink would
        // otherwise drive gamma below zero (vv042). Fade it linearly to zero between gamma = 0.01 and gamma = 0.
        destructionTerm_Ga *= std::min(1.0, std::max(m_solPhysData[iKPD+2] - m_gammaFloor, 0.0)/0.01);
        // sensitivity option (not in the paper): sink only before onset
        if (m_sinkPreOnsetOnly && m_onsetFunction) destructionTerm_Ga *= (1.0 - fOnset);
        // EGammaWallSource: behind the onset the wall law raises the wall gamma above the layer's, and eq. 37 with
        // dgamma/dn < 0 (dgammadn > 0) carries it away from the wall. Kept up to the local wall-law value, so gamma
        // cannot rise above the wall's this way (the transport of Jatayu v1.5 brings no new extrema either).
        if (m_eGammaWallSource && fOnset > 0.0)
        {
          const CFreal cap = wallLawCap(*(*m_cellStates)[iSol]);
          const CFreal feed = std::max(0.0, 0.25*muGamma*dgammadn*blTerm);
          destructionTerm_Ga += fOnset*feed*std::min(1.0, std::max(cap - m_solPhysData[iKPD+2], 0.0)/0.01);
        }
        // EGammaImplicitDiagonal: Jacobian only. The E_gamma terms advect gamma along n with rho w = 0.25 mu_gamma blTerm;
        // their dependence on gamma goes through the gradient, which the cell-local perturbation does not see. The
        // upwind diagonal rho |w|/dn, with 1/dn = (p+1)/2 sum_d |grad xi_d . n|, restores it; zero in the residual.
        if (m_eGammaImplicitDiagonal && m_isPerturbed)
        {
          const CFreal dGamma = m_solPhysData[iKPD+2] - m_unpertStates[iSol][gammaID];
          if (dGamma != 0.0)
          {
            DataHandle< CFreal > solPntNormals = socket_solPntNormals.getDataHandle();
            const CFuint sID = (*m_cellStates)[iSol]->getLocalID();
            const CFreal nx = -v/avV;
            const CFreal ny = u/avV;
            CFreal invDn = 0.0;
            for (CFuint iDir = 0; iDir < 2; ++iDir)
            {
              invDn += std::fabs(solPntNormals[sID*4+iDir*2]*nx + solPntNormals[sID*4+iDir*2+1]*ny);
            }
            invDn *= 0.5*(m_order + 1.0)/m_solPntJacobDets[iSol];
            destructionTerm_Ga -= 0.25*muGamma*blTerm*invDn*dGamma;
          }
        }
      }
      else {
        destructionTerm_Ga = (-1.0) * cEg/0.57 * muGamma * avV/(uInfLocal*uInfLocal) * dudn * dgammadn;
      }
    
      //destructionTerm_Ga = min(max(-10.0*fabs(prodTerm_Ga),destructionTerm_Ga),10.0*fabs(prodTerm_Ga));
      //if (avGa<0.1) destructionTerm_Ga = min(max(-fabs(prodTerm_Ga),destructionTerm_Ga),fabs(prodTerm_Ga));//if (avGa<0.4) destructionTerm_Ga = max(-fabs(prodTerm_Ga),destructionTerm_Ga);
      //destructionTerm_Ga = min(0.0,destructionTerm_Ga);
    
      // dissipation limit of the thesis FR implementation (eq. 7.90 states -20 P_gamma; the legacy code used -17).
      // The limit uses the clamped-gamma P_gamma: with GammaEquilibrium the raw (1-gamma) makes P_gamma = 0 at
      // gamma = 1, which would clip every E_gamma to zero and lock gamma = 1 (no laminar boundary layer can form).
      const CFreal prodTerm_GaLimit = fOnset*2.0*fg*(1.0-avGa)*sqrt(-log(1.0-avGa))*beta*rho*avV;
      // 2025: the paper has no such limiter; near the leading edge (gamma ~ 0.99, small U) it throttled E_gamma so
      // much that gamma stayed ~0.99 at the wall up to x ~ 0.15 m (vv045)
      if (m_modelVersion != 2025) destructionTerm_Ga = max(-17.0*fabs(prodTerm_GaLimit),destructionTerm_Ga);
    
      // compute production term of alpha
      const CFreal cpa1 = 0.03;
      //const CFreal cpa2 = 100.0;//50.0;
    
      const CFreal alphac = rhoInfLocal*uInfLocal*uInfLocal*pow(0.09+lambda,0.62)/(Rethetac*sqrt(rhoInfLocal*muInfLocal));
    
      const CFreal invT = detail::transitionRelaxationRate(rho, mu, avV);
    
      const CFreal dkdn = (v/avV)*dkx - (u/avV)*dky;//-dky; //
    
      const CFreal  Rew         = (rho * m_currWallDist[iSol] * m_currWallDist[iSol] * avOmega)/(mu);
      const CFreal  Fwake = Physics::NavierStokes::transition_detail::expNegativeSquare(1.0e-5 * Rew);
      const CFreal  coefFtheta1 = detail::transitionWakeExponential(rho, mu, ReThetat, m_vorticity, avV);
      const CFreal  Ftheta1     = Fwake * coefFtheta1;
      // F_Palpha gamma term: 2025 eq. A12 (gamma - 1/C_Pa2)/(1 - 1/C_Pa2), C_Pa2 = 50; legacy code (0.01, 0.98)
      const CFreal gF = (m_modelVersion == 2025) ? (avGa-0.02)/(1.0-0.02) : (avGa-0.01)/(0.99-0.01);
      const CFreal  Ftheta3     = 1.0-gF*gF;
      const CFreal  Ftheta4     = std::max(Ftheta1,Ftheta3);
      CFreal  Fthetat     = std::min(Ftheta4,1.0);
    
      dkdnFactor = fabs(dkdn)/max(avK,1.0e-6);

      if (dkdnFactor>100.0) Fthetat = 1.0;
    
      prodTerm_alpha = cpa1 * (rho*invT) * (alphac - avAlpha) * (1.0 - Fthetat);

      // compute destruction term of alpha
      const CFreal cna = 0.4;
      const CFreal speedDadn = detail::gammaAlphaCrossGradient(u, v, dax, day);
    
      if (m_nAlphaFlux) {
        // Flux form: rho*alpha*w (w = cna (v,-u)) is transported by the convective operator; the remaining
        // part of -rho w.grad(alpha) = -div(rho w alpha) + alpha div(rho w) is this local term.
        const CFreal Tl = m_solPhysData[EulerTerm::T];
        const CFreal dTx = (*(m_cellGrads[iSol][3]))[XX];
        const CFreal dTy = (*(m_cellGrads[iSol][3]))[YY];
        const CFreal drhox = (dpx - rho*R*dTx)/(R*Tl);
        const CFreal drhoy = (dpy - rho*R*dTy)/(R*Tl);
        destructionTerm_alpha = cna*m_solPhysData[iKPD+3]*(rho*(dvx - duy) + v*drhox - u*drhoy);
      }
      else {
        destructionTerm_alpha = -cna*rho*speedDadn;
      }
    
    }

    //destructionTerm_alpha = max(-10.0*fabs(prodTerm_alpha),destructionTerm_alpha);

    ///Make sure negative values dont propagate
    // Ceiling (2025 with GammaEquilibrium): keep P_gamma ~ (1-gamma) negative above gamma = 1, so gamma = 1 attracts
    // from both sides; with the clip an implicit overshoot above 1 has no restoring term (vv042).
    if (!(m_gammaEquilibrium && m_modelVersion == 2025)) prodTerm_Ga = max(0., prodTerm_Ga);

    if (m_limPRe) prodTerm_alpha     = max(0., prodTerm_alpha);
    
    if (m_limPRe) destructionTerm_Ga = min(0., destructionTerm_Ga);
    
    if (m_limPRe) destructionTerm_alpha = min(0., destructionTerm_alpha);

    // E_gamma flux form: the local part gamma div(rho w_E) of -E_gamma; the transport part gamma rho w_E is in the
    // convective flux (Euler2DGammaAlphaPuvt, AUSMPlusUpFluxGammaAlpha2D). Replaces the pointwise E_gamma.
    if (m_eGammaFlux) destructionTerm_Ga = m_solPhysData[iKPD+2]*m_eGammaDivergence[iSol];
    if (m_eGammaUpwind)
    {
      // E_gamma upwind non-conservative form: -rho w_E.grad(gamma) with the central (BR2) gradient of the cell; the
      // upwind dissipation is in the gamma interface flux (AUSMPlusUpFluxGammaAlpha2D). Not clipped.
      CFreal fx = 0., fy = 0.;
      Physics::GammaAlpha::eGammaMassFlux(0.25, rho, u, v, m_solPhysData[EulerTerm::P], m_solPhysData[iKPD+2],
                                          m_eulerVarSet->getModel()->getR(), m_eulerVarSet->getModel()->getGamma(),
                                          m_eulerVarSet->getModel()->getMachInf(), m_eulerVarSet->getModel()->getVelInf(),
                                          m_eulerVarSet->getModel()->getStaticPressInf(), fx, fy);
      destructionTerm_Ga = -(fx*(*(m_cellGrads[iSol][gammaID]))[XX] + fy*(*(m_cellGrads[iSol][gammaID]))[YY]);
    }
      
    /// Compute the rhs contribution
    // and Store the unperturbed source terms
    resUpdates[m_nbrEqs*iSol + kID] = m_prodTerm_k + m_destructionTerm_k;
    resUpdates[m_nbrEqs*iSol + omegaID] = m_prodTerm_Omega + m_destructionTerm_Omega;
    resUpdates[m_nbrEqs*iSol + gammaID] = prodTerm_Ga;// + destructionTerm_Ga;
    resUpdates[m_nbrEqs*iSol + alphaID] = prodTerm_alpha + destructionTerm_alpha;
    
    if (m_addDGDA)
    {
      resUpdates[m_nbrEqs*iSol + gammaID] += destructionTerm_Ga;
    }

    if (m_onsetFunction && m_debugProbeX > -1.0e29 && m_onsetStamp % 50 == 1)
    {
      const State& probeState = *(*m_cellStates)[iSol];
      const CFuint pid = probeState.getLocalID();
      if (fabs(probeState.getCoordinates()[XX] - m_debugProbeX) < 0.02 && m_currWallDist[iSol] < 2.0e-3 && probeState.isParUpdatable())
      {
        if (pid >= m_debugLogged.size()) m_debugLogged.resize(pid+1,0);
        if (m_debugLogged[pid] != m_onsetStamp)
        {
          m_debugLogged[pid] = m_onsetStamp;
          CFLog(INFO, "GammaAlpha probe: x " << probeState.getCoordinates()[XX] << " d " << m_currWallDist[iSol]
                << " gamma " << m_solPhysData[iKPD+2] << " c " << avV << " Fonset " << fOnset << " fg " << fg
                << " beta " << beta << " TuInfLocal " << TuInfLocal << " Tu " << Tu << " alpha " << m_solPhysData[iKPD+3]
                << " Pg " << prodTerm_Ga << " Eg " << destructionTerm_Ga << "\n");
        }
      }
    }
    
    resUpdates[m_nbrEqs*iSol + 3] = -m_prodTerm_k - m_destructionTerm_k;
    
    if (m_isAxisymmetric)
    { 
      const CFreal rhovr = rho*m_vOverRadius;
      
      const CFreal divV = dux + dvy;
      const CFreal lambdaQ = navierStokesVarSet->getThermConductivity(*((*m_cellStates)[iSol]), mu + mut * navierStokesVarSet->getModel().getPrOverPrT());
      const CFreal coeffTauMuAS = coeffTauMu*(mu+mut);
      const CFreal tauRX = coeffTauMuAS*(duy + dvx);
      const CFreal tauRR = coeffTauMuAS*(2.*dvy - 2.0/3.0*(divV + m_vOverRadius));
      const CFreal tauTT = -coeffTauMuAS*2.0/3.0*(divV - 2.*m_vOverRadius);
      const CFreal qr = -navierStokesVarSet->getModel().getCoeffQ()*lambdaQ*(*(m_cellGrads[iSol][3]))[YY];
      
      //rho
      resUpdates[m_nbrEqs*iSol + 0] = -rhovr;
      
      //rho u
      resUpdates[m_nbrEqs*iSol + 1] = -rhovr*u + m_overRadius*tauRX;
      
      //rho v
      resUpdates[m_nbrEqs*iSol + 2] = -rhovr*v + m_overRadius*(tauRR - tauTT);
      
      //rho e
      resUpdates[m_nbrEqs*iSol + 3] += -rhovr*m_solPhysData[EulerTerm::H] + m_overRadius*(-qr + tauRX*u + tauRR*v);
      
      // the radial parts of the k, omega and alpha diffusion take the mu_t of the planar diffusive fluxes
      // (NavierStokes2DGammaAlphaPuvt): the SST mu_t = mu_t*/f(gamma) for k and omega with KOmegaDiffusionSST, and for
      // alpha with ModelVersion 2025 (eq. A15); mu_t* otherwise (v23; before, always mu_t*)
      const CFreal mutKOmegaAxi = m_kOmegaDiffusionSST ? mut/muTStarFactor : mut;
      const CFreal mutAlphaAxi = (m_modelVersion == 2025) ? mut/muTStarFactor : mut;

      //rho k
      resUpdates[m_nbrEqs*iSol + 4] += -rhovr*avK + m_overRadius*(mu+navierStokesVarSet->getSigmaK()*mutKOmegaAxi)*(*(m_cellGrads[iSol][kID]))[YY];
      
      //rho logOmega
      resUpdates[m_nbrEqs*iSol + 5] += -rhovr*m_solPhysData[iKPD+1] + m_overRadius*(mu+navierStokesVarSet->getSigmaOmega()*mutKOmegaAxi)*(*(m_cellGrads[iSol][omegaID]))[YY];
      
      //rho gamma
      resUpdates[m_nbrEqs*iSol + 6] += -rhovr*gaRaw + m_overRadius*muGamma*dgammay;
      
      //rho alpha
      resUpdates[m_nbrEqs*iSol + 7] += -rhovr*avAlpha + m_overRadius*2.0*(mu+mutAlphaAxi)*day;
    }
    
    //if (m_solPhysData[iKPD+2]<0) CFLog(INFO, "ST update: " << resUpdates << "\n");
    
    if (!m_isPerturbed)
    {
      // store the shear stress velocity
      DataHandle< CFreal > wallShearStressVelocity = socket_wallShearStressVelocity.getDataHandle();
      
      DataHandle< CFreal > MInfLocalSocket = socket_MInfLocal.getDataHandle();
      
      DataHandle< CFreal > uInfLocalSocket = socket_uInfLocal.getDataHandle();
      
      DataHandle< CFreal > TInfLocalSocket = socket_TInfLocal.getDataHandle();
      
      DataHandle< CFreal > TuInfLocalSocket = socket_TuInfLocal.getDataHandle();
      
      DataHandle< CFreal > alphaDiffSocket = socket_alphaDiff.getDataHandle();
      
      const CFreal mu = navierStokesVarSet->getLaminarDynViscosityFromGradientVars(*((*m_cellStates)[iSol]));
      
      const CFreal rho = navierStokesVarSet->getDensity(*((*m_cellStates)[iSol]));
    
      const CFreal nuTot = (mu + mut)/rho;
      
      const CFreal dUdY = (*(m_cellGrads[iSol][1]))[YY];
            
      // take the absolute value of dUdY to avoid nan which causes tecplot to be unable to load the file
      wallShearStressVelocity[(((*m_cellStates)[iSol]))->getLocalID()] = sqrt(nuTot*fabs(dUdY));//m_prodTerm_Omega;//std::min(m_prodTerm_Omega,200.0);//m_prodTerm_Omega;//
      
      MInfLocalSocket[(((*m_cellStates)[iSol]))->getLocalID()] = dkdnFactor;//MInfLocal;//destructionTerm_Ga;//prodTerm_Ga;//dudn;//
      
      uInfLocalSocket[(((*m_cellStates)[iSol]))->getLocalID()] = rawDpds;//uInfLocal;//muGamma;//destructionTerm_Ga;//dgammadn;//
      
      TInfLocalSocket[(((*m_cellStates)[iSol]))->getLocalID()] = dpds;//TInfLocal;//dgammadn;//dudn;//prodTerm_alpha;//dadn;//
      
      TuInfLocalSocket[(((*m_cellStates)[iSol]))->getLocalID()] = TuInfLocal;//avV/(uInfLocal*uInfLocal);//destructionTerm_alpha;//destructionTerm_Ga;//
      
      alphaDiffSocket[(((*m_cellStates)[iSol]))->getLocalID()] = ((*m_cellStates)[iSol]->getCoordinates())[YY];//m_currWallDist[iSol];//alphac - avAlpha;//dudn;//muGamma;//
    }
  }
}
 
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

CFreal GammaAlpha2DSourceTerm::getRethetat(const CFreal Tu, bool prime) 
{
   //if (!(Tu >= 0.0)) CFLog(INFO, "Tu problem: " << Tu << "\n");
  //cf_assert(Tu >= 0.0);   
    
//  const CFreal overTu = 1.0/Tu;
//           
//  if (Tu<=1.3) 
//  {
//    m_Rethetat = (1173.51-589.428*Tu + 0.2196*overTu*overTu);
//  }
//  else 
//  {
//    const CFreal lamco5   = Tu - 0.5658;
//    const CFreal pwtu   = -0.671;
//    
//    m_Rethetat = 331.5*std::pow(lamco5,pwtu);
//  }
  CFreal Rethetat;
  
  if (!prime)
  {
    Rethetat = 0.664*sqrt(400094.0*pow(Tu,-1.38) - 105254.0*pow(Tu,-7.0/8.0));
  }
  else
  {
    Rethetat = 0.332/sqrt(400094.0*pow(Tu,-1.38) - 105254.0*pow(Tu,-7.0/8.0))*(-1.38*400094.0*pow(Tu,-2.38)+105254.0*7.0/8.0*pow(Tu,-15.0/8.0));
  }
  
  return Rethetat;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

CFreal GammaAlpha2DSourceTerm::getRethetatwithPressureGradient(const CFreal avAlpha, const CFreal rhoInfLocal, const CFreal uInfLocal, const CFreal muInfLocal, const CFreal mu, const CFreal rho, const CFreal dpx, const CFreal dpy)
{  
  const CFreal alphaTerm = rhoInfLocal*uInfLocal*uInfLocal/(avAlpha*sqrt(rhoInfLocal*muInfLocal));
  CFreal ReThetat0 = alphaTerm*0.2247141166;
  CFreal ReThetat1 = 0.0;
  const CFreal avV = m_solPhysData[EulerTerm::V]; 
  const CFreal u = m_solPhysData[EulerTerm::VX];
  const CFreal v = m_solPhysData[EulerTerm::VY]; 
  CFreal dpds = (u/avV)*dpx + (v/avV)*dpy;
  const CFreal lambdaTerm1 = -mu/(rho*rho*uInfLocal*uInfLocal*uInfLocal);

  const CFreal dpdsLimit = -m_lambdaLim/(lambdaTerm1*ReThetat0*ReThetat0);
  dpds = min(max(dpds,-dpdsLimit),dpdsLimit);

  const CFreal lambdaTerm = lambdaTerm1*dpds;
  
  const CFuint MAXITER   = 10;
  const CFreal TOL       = 1.0e-6;
  
  for (CFuint iter = 0; iter < MAXITER; ++iter)
  {
    CFreal resReThetat = std::fabs(ReThetat0*TOL);
    
    const CFreal lambda = max(ReThetat0*ReThetat0*lambdaTerm,-0.089);
    const CFreal SLambda = pow(0.09+lambda,0.62);
    //ReThetat1 = alphaTerm*SLambda;

    const CFreal derivDenom = 1.0/pow(0.09+lambda,0.38);

    const CFreal mainF = alphaTerm*SLambda - ReThetat0;
    const CFreal mainFprime = 1.24*alphaTerm*lambdaTerm*ReThetat0*derivDenom - 1.0; 
//CFLog(INFO, "lambdaTerm: " << lambdaTerm << ", alphaTerm: " << alphaTerm << ", Ret0: " << ReThetat0 << ", dpds: " << dpds << ", avV: " << avV << ", lambda: " << lambda << "\n");
    ReThetat1 = max(ReThetat0 - mainF/mainFprime,10.0);
    
    if ( fabs(ReThetat0-ReThetat1) < resReThetat || iter == MAXITER-1) 
    { 
      if (iter == MAXITER-1) CFLog(INFO, "Max iter thetat reached!, Re0: " << ReThetat0 << ", Re1: " << ReThetat1 << ", mainF: " << mainF << ", derivF: " << mainFprime << ", alpha: " << avAlpha << "\n");
          
      break;
    }
    
    if (ReThetat0 != ReThetat0 || ReThetat1 != ReThetat1) 
    { 
      CFLog(INFO, "Nan detected!, Re0: " << ReThetat0 << ", Re1: " << ReThetat1 << ", mainF: " << mainF << ", derivF: " << mainFprime << ", alpha: " << avAlpha << "\n");
          
      break;
    }
    
    ReThetat0 = ReThetat1;
  }
  
  return ReThetat1;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

CFreal GammaAlpha2DSourceTerm::getTuInfLocal(const CFreal Rethetat, const CFreal MInfLocal, const CFreal TuLocal)
{  
  const CFreal fM = sqrt(1.0+0.38*pow(MInfLocal,0.6));
  CFreal Tu0 = TuLocal;
  CFreal Tu1 = 0.0;
  
  const CFuint MAXITER   = 50;
  const CFreal TOL       = 1.0e-5;
  
  for (CFuint iter = 0; iter < MAXITER; ++iter)
  {
    CFreal resTu = fabs(Tu0*TOL);
    
    const CFreal mainF = fM*getRethetat(Tu0,false) - Rethetat;
    const CFreal mainFprime = fM*getRethetat(Tu0,true);   
    
    Tu1 = min(max(Tu0 - mainF/mainFprime,0.027),13.0);
    
    //if (!(Tu1>=0.0)) CFLog(INFO, "Tu1: " << Tu1 << ", Tu0: " << Tu0 << ", fM: " << fM << ", mainFPrime: " << mainFprime << "\n");
    
    if ( fabs(Tu0-Tu1) < resTu || iter == MAXITER-1) 
    { 
      if (iter == MAXITER-1) CFLog(INFO, "Max iter TuInfLocal reached! Tu0: " << Tu0 << ", Tu1: " << Tu1 << ", Re: " << Rethetat << ", MInf: " << MInfLocal << ", TuLocal: " << TuLocal << "\n");
          
      break;
    }
    
    Tu0 = Tu1;
  }
  
  return Tu1;
}

//////////////////////////////////////////////////////////////////////////////////////////////////

void  GammaAlpha2DSourceTerm::getSToStateJacobian(const CFuint iState)
{
  // Check the analytical domain before resetting or filling the matrix.
  const State& sourceState = *(*m_cellStates)[iState];
  detail::validateTransitionSourceState(sourceState, 2, sourceState.getLocalID(), !m_signedK);
  Physics::NavierStokes::transition_detail::validateThermodynamics(sourceState[0], sourceState[3], m_eulerVarSet->getModel()->getR(), m_eulerVarSet->getModel()->getGamma());
  const auto freestream = Physics::NavierStokes::transition_detail::gammaAlphaFreestream(std::max(sourceState[0],1.0e-3), m_eulerVarSet->getModel()->getR(), m_eulerVarSet->getModel()->getGamma(), m_eulerVarSet->getModel()->getMachInf(), m_eulerVarSet->getModel()->getVelInf(), m_eulerVarSet->getModel()->getStaticPressInf());
  if (!(std::isfinite(m_currWallDist[iState]) && m_currWallDist[iState] >= 0.))
    throw Common::BadValueException(FromHere(), "Transition source requires finite nonnegative wall distance");
  for (CFuint e = 0; e < m_nbrEqs; ++e)
    for (CFuint d = 0; d < 2; ++d)
      if (!std::isfinite((*m_cellGrads[iState][e])[d]))
        throw Common::BadValueException(FromHere(), "Transition source requires finite gradients");
  CFreal velocityGradients[3][3] = {};
  for (CFuint e=0; e<2; ++e) for (CFuint d=0; d<2; ++d)
    velocityGradients[e][d] = (*m_cellGrads[iState][e+1])[d];
  const CFreal axisGradient = m_isAxisymmetric ? Physics::NavierStokes::transition_detail::finiteProduct(1./std::max(sourceState.getCoordinates()[YY],1.e-4),sourceState[2]) : 0.;
  Physics::NavierStokes::transition_detail::validateVelocityGradientQuadratics(velocityGradients,2,axisGradient);
  auto transport = m_diffVarSet.d_castTo<NavierStokes2DGammaAlphaPuvt>();
  transport->getLaminarDynViscosityFromGradientVars(sourceState);
  Physics::NavierStokes::transition_detail::laminarViscosity(transport->getModel(), sourceState[0], freestream.temperature,
      m_eulerVarSet->getModel()->getPressRef(), m_eulerVarSet->getModel()->getTempRef());
  const CFreal avV = detail::transitionAnalyticalSpeed(sourceState[1], sourceState[2], 0.0);

  // reset the jacobian
  for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
  {
    m_stateJacobian[iEq] = 0.0;
  }
  
  SafePtr< NavierStokes2DKLogOmega > navierStokesVarSet = m_diffVarSet.d_castTo< NavierStokes2DKLogOmega >();
  
  // Set the wall distance before computing the turbulent viscosity
  navierStokesVarSet->setWallDistance(m_currWallDist[iState]);
    
  /// destruction term of k
    
  const CFreal avK     = max((*((*m_cellStates)[iState]))[4],0.0);
  const CFreal avOmega = exp((*((*m_cellStates)[iState]))[5]);
  const CFreal T = max(1.0e-3,(*((*m_cellStates)[iState]))[3]);
  const CFreal p = max(1.0e-3,(*((*m_cellStates)[iState]))[0]);
  const CFreal rho = max(navierStokesVarSet->getDensity(*((*m_cellStates)[iState])),1.0e-8);
  const CFreal R = m_eulerVarSet->getModel()->getR();
  const CFreal overRT = 1.0/(R*T);
  const CFreal pOverRTT = p/(R*T*T);
  const CFreal overOmega = 1.0/avOmega;
  
  const CFreal avGa    = min(max((*((*m_cellStates)[iState]))[6],0.01),0.99);
    
  const CFreal avAlpha    = max((*((*m_cellStates)[iState]))[7],0.0);
  
  const CFreal mu = navierStokesVarSet->getLaminarDynViscosityFromGradientVars(*((*m_cellStates)[iState]));
  const CFreal mut = navierStokesVarSet->getTurbDynViscosityFromGradientVars(*((*m_cellStates)[iState]), m_cellGrads[iState]);
  const CFreal rhoKOverMut = m_diffVarSet.d_castTo< NavierStokes2DGammaAlphaPuvt >()->getRhoKOverMuT(*((*m_cellStates)[iState]), m_cellGrads[iState]);
  
  if (m_isAxisymmetric)
  {
    m_overRadius = 1.0/max(((*m_cellStates)[iState]->getCoordinates())[YY],1.0e-4);
    m_vOverRadius = m_overRadius*(*((*m_cellStates)[iState]))[2];
  }
  
  //Compute Strain 
  getStrain(iState);
  
  // Get Vorticity
  getVorticity(iState);
  
  const CFreal u = (*((*m_cellStates)[iState]))[1];
  const CFreal v = (*((*m_cellStates)[iState]))[2];
  const CFreal dpx = (*(m_cellGrads[iState][0]))[XX];
  const CFreal dpy = (*(m_cellGrads[iState][0]))[YY];
    
  const CFreal gammaIsentropic = m_eulerVarSet->getModel()->getGamma();
  const CFreal rhoInfLocal = freestream.density;
  const CFreal MInfLocal = freestream.mach;
  const CFreal uInfLocal = freestream.speed;
  const CFreal TInfLocal = freestream.temperature;
    
  const CFreal Tu = detail::transitionIntensity(avK, avV, 13.0);
      
  // Allocate owned storage: copy construction can alias a non-owning state.
  RealVector freestreamState((*m_cellStates)[iState]->size());
  freestreamState = *((*m_cellStates)[iState]);
  freestreamState[3] = TInfLocal;
  const CFreal muInfLocal = navierStokesVarSet->getLaminarDynViscosityFromGradientVars(freestreamState);
      
  CFreal ReThetat;

  if (!m_PGrad)
  {
    ReThetat = rhoInfLocal*uInfLocal*uInfLocal*0.2247141166/(avAlpha*sqrt(rhoInfLocal*muInfLocal));
  }
  else
  {
    ReThetat = getRethetatwithPressureGradient(avAlpha,rhoInfLocal,uInfLocal,muInfLocal,mu,rho,dpx,dpy);
  }

  const CFreal TuInfLocal = getTuInfLocal(ReThetat,MInfLocal,Tu);
  
  const CFreal dux = (*(m_cellGrads[iState][1]))[XX];
  const CFreal duy = (*(m_cellGrads[iState][1]))[YY]; 
  const CFreal dvx = (*(m_cellGrads[iState][2]))[XX]; 
  const CFreal dvy = (*(m_cellGrads[iState][2]))[YY]; 
  const CFreal dgammax = (*(m_cellGrads[iState][6]))[XX];
  const CFreal dgammay = (*(m_cellGrads[iState][6]))[YY];
  const CFreal dkx = (*(m_cellGrads[iState][4]))[XX];
  const CFreal dky = (*(m_cellGrads[iState][4]))[YY];
  const CFreal dax = (*(m_cellGrads[iState][7]))[XX];
  const CFreal day = (*(m_cellGrads[iState][7]))[YY];
  
  const CFreal betaStar = navierStokesVarSet->getBetaStar(*((*m_cellStates)[iState]));
  
  const CFreal DkTerm = avOmega * betaStar;
  
  const CFreal Dk = -rho * avK * DkTerm;
  
  const CFreal dmudg = 0.5*(1-2.0*avGa)*(tanh((avGa-0.25)/0.1)+1) + 0.5/0.1*(1-avGa)*avGa/(pow(cosh((avGa-0.25)/0.1),2)) + 1.0;

  //const CFreal limGamma = max(0.01,avGa);
  const CFreal mutGaMod = (avGa+avGa*(1.0-avGa)*0.5*(1.0+tanh((avGa-0.25)/(0.1))));
    
  if (Dk < 0.0)
  {
    CFreal tempSTTerm = 0.0;
  
    if (!m_blockDecoupled)
    {
      //p
      tempSTTerm = -avK * overRT * DkTerm;
      m_stateJacobian[0][4] += tempSTTerm;
      m_stateJacobian[0][3] -= tempSTTerm;
  
      //T
      tempSTTerm = DkTerm * avK * pOverRTT;
      m_stateJacobian[3][4] += tempSTTerm;
      m_stateJacobian[3][3] -= tempSTTerm;
    }
  
    //k
    tempSTTerm = -DkTerm * rho;
    m_stateJacobian[4][4] += tempSTTerm;
    m_stateJacobian[4][3] -= tempSTTerm;
  
    //logOmega
    tempSTTerm = -DkTerm * avK * rho;
    m_stateJacobian[5][4] += tempSTTerm;
    m_stateJacobian[5][3] -= tempSTTerm;
  }
  
  /// destruction term of logOmega
  
  const CFreal beta = navierStokesVarSet->getBeta(*((*m_cellStates)[iState]));
  
  const CFreal DomegaTerm = avOmega * beta;
  
  const CFreal Domega = -DomegaTerm * rho;
  
  if (Domega < 0.0)
  {
    if (!m_blockDecoupled)
    {
      //p
      m_stateJacobian[0][5] = -DomegaTerm * overRT;
  
      //T
      m_stateJacobian[3][5] = DomegaTerm * pOverRTT;
    }
  
    //logOmega
    m_stateJacobian[5][5] = -DomegaTerm * rho;
  }
  
  /// production term of k

  const CFreal coeffTauMu = navierStokesVarSet->getModel().getCoeffTau();
  //const CFreal mutTerm = coeffTauMu*((4./3.)*((dux-dvy)*(dux-dvy)+(dux*dvy)-(dux+dvy-m_vOverRadius)*m_vOverRadius)+(duy+dvx)*(duy+dvx));
  const CFreal gammaTerm = avGa + avGa*(1.0-avGa);
  
  CFreal mutTerm;
  CFreal Pk;
  CFreal productionOverMut;
  
  if (!m_isSSTV)
  {
    mutTerm = coeffTauMu*((4./3.)*((dux-dvy)*(dux-dvy)+(dux*dvy)-(dux+dvy-m_vOverRadius)*m_vOverRadius)+(duy+dvx)*(duy+dvx));

    Pk = (mutTerm*mut - (2./3.)*(avK * rho * gammaTerm)*(dux+dvy+m_vOverRadius));
    productionOverMut = (mutTerm - (2./3.)*(rhoKOverMut * gammaTerm)*(dux+dvy+m_vOverRadius));
  }
  else if (m_neglectSSTVTerm)
  {
    mutTerm = coeffTauMu*(duy-dvx)*(duy-dvx);

    Pk = mutTerm*mut;
    productionOverMut = mutTerm;
  }
  else
  {
    mutTerm = coeffTauMu*(duy-dvx)*(duy-dvx);

    Pk = (mutTerm*mut - (2./3.)*(avK * rho)*(dux+dvy+m_vOverRadius));
    productionOverMut = (mutTerm - (2./3.)*rhoKOverMut*(dux+dvy+m_vOverRadius));
  }
  
  //const CFreal Pk = (mutTerm*mut - (2./3.)*(avK * rho * gammaTerm)*(dux+dvy+m_vOverRadius));
  
  const CFreal twoThirdduxduy = (2./3.)*(dux+dvy+m_vOverRadius);
  
  if (Pk > 0.0)
  {
    CFreal tempSTTerm = 0.0;
       
    if (!m_blockDecoupled)
    {
      //p
      tempSTTerm = mutTerm*avK*overRT*overOmega*mutGaMod;
      if (!m_neglectSSTVTerm) tempSTTerm += -twoThirdduxduy*avK*overRT * gammaTerm;
      
      m_stateJacobian[0][4] += tempSTTerm;
      m_stateJacobian[0][3] -= tempSTTerm;
  
      //T
      tempSTTerm = -mutTerm*avK*overOmega*mutGaMod*pOverRTT;
      if (!m_neglectSSTVTerm) tempSTTerm += twoThirdduxduy*avK*pOverRTT * gammaTerm;
      
      m_stateJacobian[3][4] += tempSTTerm;
      m_stateJacobian[3][3] -= tempSTTerm;
 
      //gamma
      tempSTTerm = mutTerm*mut*dmudg/mutGaMod;
      if (!m_neglectSSTVTerm) tempSTTerm += -twoThirdduxduy * avK * rho*2.0*(1.0-avGa);
      
      m_stateJacobian[6][4] += tempSTTerm;
      m_stateJacobian[6][3] -= tempSTTerm;
      
      if (m_isAxisymmetric && !m_isSSTV)
      {
        //v
        tempSTTerm = mut*coeffTauMu*4./3.*m_overRadius*(-(dux+dvy)+v*2.0)- (2./3.)*(avK * rho * gammaTerm)*m_overRadius;
        m_stateJacobian[2][4] += tempSTTerm;
        m_stateJacobian[2][3] -= tempSTTerm;
      }
    }
  
    //k
    tempSTTerm = mutTerm*rho*overOmega*mutGaMod;
    if (!m_neglectSSTVTerm) tempSTTerm += -twoThirdduxduy*rho * gammaTerm;
    
    m_stateJacobian[4][4] += tempSTTerm;
    m_stateJacobian[4][3] -= tempSTTerm;
  
    //logOmega
    tempSTTerm = -mutTerm*rho*avK*overOmega*mutGaMod;
    m_stateJacobian[5][4] += tempSTTerm;
    m_stateJacobian[5][3] -= tempSTTerm;
  }
  
  /// production of logOmega
  
  const CFreal gamma = navierStokesVarSet->getGammaCoef();
  const CFreal blendingCoefF1 = navierStokesVarSet->getBlendingCoefficientF1();
  const CFreal sigmaOmega2 = navierStokesVarSet->getSigmaOmega2();
  
  const CFreal pOmegaFactor = (1. - blendingCoefF1) * 2. * sigmaOmega2 * ((*(m_cellGrads[iState][4]))[XX]*(*(m_cellGrads[iState][5]))[XX] + (*(m_cellGrads[iState][4]))[YY]*(*(m_cellGrads[iState][5]))[YY]);
  
  const CFreal Pomega = gamma*rho * productionOverMut * overOmega + rho * overOmega * pOmegaFactor;
  
  if (Pomega > 0.0)
  {
    if (!m_blockDecoupled)
    {
      //p
      m_stateJacobian[0][5] += gamma*(mutTerm*overRT*overOmega - twoThirdduxduy*overRT/mutGaMod*gammaTerm) + pOmegaFactor*overRT*overOmega;
  
      //T
      m_stateJacobian[3][5] += gamma*pOverRTT*(-mutTerm*overOmega + twoThirdduxduy/mutGaMod*gammaTerm) - pOmegaFactor*pOverRTT*overOmega;
    }
 
    //logOmega
    m_stateJacobian[5][5] +=  -gamma*rho*mutTerm*overOmega - pOmegaFactor*rho*overOmega;
  }
  
  // P gamma
  const CFreal cfg1 = 1.735;
  const CFreal cfg2 = 5.45;
  const CFreal cfg3 = 0.95375;
  const CFreal cfg4 = 2.2;
    
  CFreal fg = 1.0;
    
  if (avGa < 0.45) fg -= exp(-cfg1*tan(cfg2*avGa-cfg3)-cfg4);
    
  const CFreal fMnsigmaTerm = 1.0+0.58*pow(MInfLocal,0.6);
  const CFreal fMnsigma = 1.0/(fMnsigmaTerm*fMnsigmaTerm);
    
  CFreal dpds = (u/avV)*dpx + (v/avV)*dpy;

  CFreal PRC = 1.0;
    
  if (m_PGrad)
  {
    const CFreal lambdaTerm1 = -mu/(rho*rho*uInfLocal*uInfLocal*uInfLocal);

    const CFreal dpdsLimit = -m_lambdaLim/(lambdaTerm1*ReThetat*ReThetat);
    dpds = min(max(dpds,-dpdsLimit),dpdsLimit);
    
    const CFreal kPGrad = -muInfLocal/(rhoInfLocal*rhoInfLocal*uInfLocal*uInfLocal*uInfLocal)*fabs(1.0-MInfLocal*MInfLocal)*dpds;
    
      
    if (kPGrad < 0)
    {
      const CFreal PRCfactor1 = 474.0*pow(TuInfLocal,-2.9);
      const CFreal PRCexponent = 1.0-exp(2.0e6*kPGrad);
      PRC = pow(PRCfactor1,PRCexponent);
    }
    else
    {
      const CFreal PRCexponent = -3227.0*pow(kPGrad,0.5985);
      PRC = pow(10.0,PRCexponent);
    }
  }
    
  const CFreal fk = 1.0+fg*(PRC-1.0); 
   
  const CFreal nsigma = 1.25e-11*pow(TuInfLocal,7.0/4.0)*fk*fMnsigma;
  
  const CFreal betaGA = sqrt(nsigma)*uInfLocal*rhoInfLocal/muInfLocal;
  
  const CFreal gTerm = (1-avGa)*sqrt(-log(1.0-avGa));
  
  const CFreal brc2 = betaGA*rho*avV*2.0;
  
  CFreal prodTerm_Ga = fg*gTerm*brc2;
  
  if (prodTerm_Ga > 0.0)
  {    
    const CFreal dgTermdg = (2.0*log(1-avGa)+1.0)/(2.0*sqrt(-log(1-avGa))); 
    
    CFreal dfgdg = 0.0;
    
    if (avGa<0.45) dfgdg = cfg1*cfg2/pow(cos(cfg2*avGa-cfg3),2) * exp(cfg1*tan(cfg3-cfg2*avGa)-cfg4);
    
    //const CFreal dbdg = 0.5*uInfLocal*rhoInfLocal/muInfLocal/sqrt(nsigma) * 1.25e-11*pow(TuInfLocal,7.0/4.0) * fMnsigma * (PRC-1.0) * dfgdg;
     
    //gamma
    //m_stateJacobian[6][6] += brc2 * fg * dgTermdg + brc2 * gTerm * dfgdg;// + gTerm * fg * 2.0 * rho * avV * dbdg;
   
    if (!m_blockDecoupled)
    {
      const CFreal gTermFgBeta2 = gTerm * fg * betaGA *2.0;
      
      //p
      m_stateJacobian[0][6] +=  gTermFgBeta2 * overRT * avV;
  
      //T
      m_stateJacobian[3][6] +=  -gTermFgBeta2 * pOverRTT * avV;
      
      //u
      m_stateJacobian[1][6] +=  gTermFgBeta2 * rho * (u/avV);
              
      //v
      m_stateJacobian[2][6] +=  gTermFgBeta2 * rho * (v/avV);
    }
  }
  
  // D gamma
  const CFreal cEg = m_ceg/0.57;
  
  const CFreal fMuGamma = 1.0-exp(-256.0*(m_currWallDist[iState]*uInfLocal*rhoInfLocal/muInfLocal)*(m_currWallDist[iState]*uInfLocal*rhoInfLocal/muInfLocal));
  const CFreal fMMuGamma = (1.0+0.26*(gammaIsentropic-1.0)/2.0*MInfLocal*MInfLocal)*sqrt(1+0.38*pow(MInfLocal,0.6));
  
  //const CFreal gammaLim = std::min(std::max(0.01,avGa),0.99);
  const CFreal muGamma = 0.57*pow(-log(1.0-avGa),-5.0/6.0*(1.0-avGa))*fMuGamma*fMMuGamma*mu;
    
  const CFreal dudn = -((u/avV)*(u/avV)*duy-(v/avV)*(v/avV)*dvx+(u/avV)*(v/avV)*(dvy-dux));//-duy; //
  const CFreal dgammadn = (v/avV)*dgammax - (u/avV)*dgammay;//-dgammay; //
  
  // Factor out speed before the velocity derivative; do not form V^2.
  const CFreal destructionTerm_GaOverSpeed = -cEg/(uInfLocal*uInfLocal) * dudn * dgammadn * muGamma;
  const CFreal dgTerm = -cEg * avV/(uInfLocal*uInfLocal) * dudn * dgammadn;
    
  CFreal  destructionTerm_Ga  = dgTerm * muGamma;
    
  //destructionTerm_Ga = max(-10.0*fabs(prodTerm_Ga),destructionTerm_Ga);
  //if (avGa<0.4) destructionTerm_Ga = max(-fabs(prodTerm_Ga),destructionTerm_Ga);
  
  if (destructionTerm_Ga < 0.0 && m_addDGDA)
  {
    const CFreal mgdg = muGamma * 5.0/6.0 * (log(-log(1-avGa)) + 1.0/log(1-avGa));
    
    // gamma
    //m_stateJacobian[6][6] +=  dgTerm * mgdg;
  
    if (!m_blockDecoupled)
    {
      // u
      m_stateJacobian[1][6] +=  destructionTerm_GaOverSpeed * (u/avV);
      
      // v
      m_stateJacobian[2][6] +=  destructionTerm_GaOverSpeed * (v/avV);
    }
  }
  
  // production  Alpha
  
  const CFreal cpa1 = 0.03;
  //const CFreal cpa2 = 50.0;

  CFreal Rethetac = getRethetat(Tu, false);
    
  ///@todo check if this is local M or local MInf
  Rethetac *= sqrt(1.0+0.38*pow(MInfLocal,0.6));
    
  CFreal lambda = 0.0;

  if (m_PGrad)
  {
    lambda = -Rethetac*Rethetac*mu/(rho*rho*uInfLocal*uInfLocal*uInfLocal)*dpds;
  }
    
  const CFreal alphac = rhoInfLocal*uInfLocal*uInfLocal*pow(0.09+lambda,0.62)/(Rethetac*sqrt(rhoInfLocal*muInfLocal));
    
  const CFreal invT = detail::transitionRelaxationRate(rho, mu, avV);
  
  const CFreal dkdn = (v/avV)*dkx - (u/avV)*dky;//-dky; //
    
  const CFreal  Rew         = (rho * m_currWallDist[iState] * m_currWallDist[iState] * avOmega)/(mu);   
  const CFreal  Fwake = Physics::NavierStokes::transition_detail::expNegativeSquare(1.0e-5 * Rew);
  const CFreal  coefFtheta1 = detail::transitionWakeExponential(rho, mu, ReThetat, m_vorticity, avV);
  const CFreal  Ftheta1     = Fwake * coefFtheta1;
  const CFreal  Ftheta3     = 1.0-(((avGa-0.01)/(0.99-0.01))*((avGa-0.01)/(0.99-0.01)));
  const CFreal  Ftheta4     = std::max(Ftheta1,Ftheta3);
  CFreal  Fthetat     = std::min(Ftheta4,1.0);
  
  const CFreal dkdnFactor = fabs(dkdn)/max(avK,1.0e-6);

  if (dkdnFactor>100.0) Fthetat = 1.0;
    
  CFreal prodTerm_alpha = cpa1 * (rho*invT) * (alphac - avAlpha) * (1.0 - Fthetat);

  // compute destruction term of alpha
  const CFreal cna = 0.4;
  const CFreal dadn = (v/avV)*dax - (u/avV)*day;
    
  CFreal destructionTerm_alpha = -cna*rho*avV*dadn;
  
  ///@todo Fthetat and alpha_c comtribution
  if (prodTerm_alpha > 0.0 || !m_limPRe)
  {
    // alpha
    m_stateJacobian[7][7] +=  -cpa1 * (rho*invT) * (1.0 - Fthetat);
  
    if (!m_blockDecoupled)
    {
      const CFreal tTerm = 2.0*rho*avV*avV/(500.0*mu);
        
      // p
      m_stateJacobian[0][7] +=  cpa1 * tTerm * overRT * (alphac - avAlpha) * (1.0 - Fthetat);
  
      // T
      m_stateJacobian[3][7] +=  -cpa1 * tTerm * pOverRTT * (alphac - avAlpha) * (1.0 - Fthetat);
  
      // u
      m_stateJacobian[1][7] +=  cpa1 * 2.0*rho*rho*u/(500.0*mu) * (alphac - avAlpha) * (1.0 - Fthetat);
  
      // v
      m_stateJacobian[2][7] +=  cpa1 * 2.0*rho*rho*v/(500.0*mu) * (alphac - avAlpha) * (1.0 - Fthetat);
    }
  }
  
  if (destructionTerm_alpha < 0.0 && !m_blockDecoupled)
  {
      const CFreal daTerm = -cna*dadn;
        
      // p
      m_stateJacobian[0][7] +=  daTerm*avV*overRT;
  
      // T
      m_stateJacobian[3][7] +=  -daTerm*avV*pOverRTT;
  
      // u
      m_stateJacobian[1][7] +=  daTerm*rho*(u/avV);
  
      // v
      m_stateJacobian[2][7] +=  daTerm*rho*(v/avV);
  }
  
  if (m_isAxisymmetric)
  { 
      const CFreal rhovr = rho*v*m_overRadius;
      const CFreal vrOverRT = v*m_overRadius*overRT;
      const CFreal vrPOverRTT = v*m_overRadius*pOverRTT;
      const CFreal rhor = m_overRadius*rho;
      
      const CFreal coeffTauMuAS = coeffTauMu*(mu+mut);
      
      //rho
      m_stateJacobian[0][0] += -vrOverRT;
      m_stateJacobian[2][0] += -rhor;
      m_stateJacobian[3][0] += vrPOverRTT;
      
      //rho u
      m_stateJacobian[0][1] += -vrOverRT*u;
      m_stateJacobian[2][1] += -rhor*u;
      m_stateJacobian[3][1] += vrPOverRTT*u;
      m_stateJacobian[1][1] += -rhovr;
      
      //rho v
      m_stateJacobian[0][2] += -vrOverRT*v;
      m_stateJacobian[2][2] += -2.0*rhovr-m_overRadius*m_overRadius*2.0*coeffTauMuAS;
      m_stateJacobian[3][2] += vrPOverRTT*v;
      
      //rho e
      const CFreal isentropicTerm = gammaIsentropic/(gammaIsentropic-1.0);
      m_stateJacobian[0][3] += -isentropicTerm*m_vOverRadius-0.5*avV*avV*vrOverRT;
      m_stateJacobian[2][3] += -isentropicTerm*p*m_overRadius-0.5*(avV*avV+2.0*v*v)*rho*m_overRadius - 4.0/3.0*m_overRadius*m_overRadius*coeffTauMuAS*v;
      m_stateJacobian[3][3] += 0.5*avV*avV*vrPOverRTT;
      m_stateJacobian[1][3] += -rhovr*u;
      
      //rho k
      if (!m_blockDecoupled)
      {
        m_stateJacobian[0][4] += -vrOverRT*avK;
        m_stateJacobian[2][4] += -rhor*avK;
        m_stateJacobian[3][4] += vrPOverRTT*avK;
      }
      m_stateJacobian[4][4] += -rhovr;
      
      //rho logOmega
      const CFreal logOmega = (*((*m_cellStates)[iState]))[5];
      if (!m_blockDecoupled)
      {
        m_stateJacobian[0][5] += -vrOverRT*logOmega;
        m_stateJacobian[2][5] += -rhor*logOmega;
        m_stateJacobian[3][5] += vrPOverRTT*logOmega;
      }
      m_stateJacobian[5][5] += -rhovr;
      
      //rho gamma
      if (!m_blockDecoupled)
      {
        m_stateJacobian[0][6] += -vrOverRT*avGa;
        m_stateJacobian[2][6] += -rhor*avGa;
        m_stateJacobian[3][6] += vrPOverRTT*avGa;
      }
      m_stateJacobian[6][6] += -rhovr;
      
      //rho alpha
      if (!m_blockDecoupled)
      {
        m_stateJacobian[0][7] += -vrOverRT*avAlpha;
        m_stateJacobian[2][7] += -rhor*avAlpha;
        m_stateJacobian[3][7] += vrPOverRTT*avAlpha;
      }
      m_stateJacobian[7][7] += -rhovr;
  }

  if (m_addUpdateCoeff)
  {
    DataHandle<CFreal> updateCoeff = socket_updateCoeff.getDataHandle();
  
    CFreal update = 0.0;
  
    const CFreal PkUpdate = max(Pk,0.0);
  
    const CFreal DkUpdate = min(Dk,-1.0e-10);
  
    const CFreal G = sqrt(-PkUpdate*rho/DkUpdate);
  
    const CFreal Gk = pow(G,1.0/gamma);
  
    const CFreal kUpdateTerm = (PkUpdate + DkUpdate)/max(avK,0.001)/(Gk-1) + m_stateJacobian[4][4];
  
    update = max(0.0,kUpdateTerm);
  
    //if (update>1000) CFLog(INFO,"kUpdate: " << update << "\n");
  
    const CFreal logOmega = (*((*m_cellStates)[iState]))[5];
  
    const CFreal Gomega = pow(G,1.0/logOmega);
  
    const CFreal PomegaUpdate = max(Pomega,0.0);
  
    const CFreal DomegaUpdate = min(Domega,0.0);
  
    const CFreal omegaUpdateTerm = (PomegaUpdate + DomegaUpdate)/max(logOmega,0.1)/(Gomega-1) + m_stateJacobian[5][5];
  
    update = max(update,omegaUpdateTerm);
  
    // get the local ID of the current sol pnt
    const CFuint solID = (*m_cellStates)[iState]->getLocalID();
  
    //if (update>1000) CFLog(INFO,"old: " << updateCoeff[solID] << ", new: " << update << ", J: " << m_solPntJacobDets[iState] << "\n\n");
  
    update *= (2.0*m_order+1)*m_solPntJacobDets[iState];//*m_solPntJacobDets[iState];
  
    if (update>0.5*updateCoeff[solID]) CFLog(INFO,"Large influence ST Dt: oldCoeff: " << updateCoeff[solID] << ", ST addition: " << update << ", J: " << m_solPntJacobDets[iState] << "\n\n");
    
    // add the wave speed update previously computed
    updateCoeff[solID] += update;
  }
  
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluxReconstructionMethod

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
