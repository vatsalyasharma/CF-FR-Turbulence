#include "Framework/MethodStrategyProvider.hh"
#include "FluxReconstructionMethod/FluxReconstructionSolverData.hh"
#include "FluxReconstructionNavierStokes/FluxReconstructionNavierStokes.hh"
#include "FluxReconstructionTurb/AUSMPlusUpFluxGammaAlpha2D.hh"
#include "GammaAlpha/Euler2DGammaAlphaPuvt.hh"
#include "GammaAlpha/GammaAlphaEGammaFlux.hh"

using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Physics::NavierStokes;

namespace COOLFluiD {
  namespace FluxReconstructionMethod {

Framework::MethodStrategyProvider<AUSMPlusUpFluxGammaAlpha2D, FluxReconstructionSolverData, RiemannFlux,
                                  FluxReconstructionNavierStokesModule>
AUSMPlusUpFluxGammaAlpha2DProvider("AUSMPlusUpFluxGammaAlpha2D");

AUSMPlusUpFluxGammaAlpha2D::AUSMPlusUpFluxGammaAlpha2D(const std::string& name) : Base(name) {}

AUSMPlusUpFluxGammaAlpha2D::~AUSMPlusUpFluxGammaAlpha2D() {}

CFreal& AUSMPlusUpFluxGammaAlpha2D::eGammaUpwindCoefficient()
{
  static CFreal coefficient = 0.;
  return coefficient;
}

RealVector& AUSMPlusUpFluxGammaAlpha2D::computeFlux(State& lState, State& rState, const RealVector& normal)
{
  Base::computeFlux(lState, rState, normal); // fills m_pData[LEFT/RIGHT] and m_rFlux
  const CFreal c = Physics::GammaAlpha::nAlphaFluxCoefficient();
  if (c > 0.) {
    const RealVector& L = m_pData[LEFT];
    const RealVector& R = m_pData[RIGHT];
    const CFuint a = m_updateVarSet->getModel()->getFirstScalarVar(0) + 3;
    const CFreal wnL = c*(L[EulerTerm::VY]*normal[XX] - L[EulerTerm::VX]*normal[YY]);
    const CFreal wnR = c*(R[EulerTerm::VY]*normal[XX] - R[EulerTerm::VX]*normal[YY]);
    const CFreal wn = 0.5*(wnL + wnR);
    m_rFlux[7] += (wn > 0.) ? L[EulerTerm::RHO]*L[a]*wn : R[EulerTerm::RHO]*R[a]*wn;
  }
  const CFreal cE = Physics::GammaAlpha::eGammaFluxCoefficient();
  if (cE > 0.) {
    // E_gamma as a transport of gamma away from the wall (GammaAlphaEGammaFlux.hh): the face mass flux is the mean
    // of the two sides' rho w_E.n, carrying gamma of the upwind side (Jatayu v1.5, wallNormalTransport)
    const RealVector& L = m_pData[LEFT];
    const RealVector& R = m_pData[RIGHT];
    const CFuint g = m_updateVarSet->getModel()->getFirstScalarVar(0) + 2;
    const CFreal Rgas = m_updateVarSet->getModel()->getR();
    const CFreal kappa = m_updateVarSet->getModel()->getGamma();
    const CFreal machInf = m_updateVarSet->getModel()->getMachInf();
    const CFreal uInf = m_updateVarSet->getModel()->getVelInf();
    const CFreal pInf = m_updateVarSet->getModel()->getStaticPressInf();
    CFreal fxL = 0., fyL = 0., fxR = 0., fyR = 0.;
    Physics::GammaAlpha::eGammaMassFlux(cE, L[EulerTerm::RHO], L[EulerTerm::VX], L[EulerTerm::VY], L[EulerTerm::P],
                                        L[g], Rgas, kappa, machInf, uInf, pInf, fxL, fyL);
    Physics::GammaAlpha::eGammaMassFlux(cE, R[EulerTerm::RHO], R[EulerTerm::VX], R[EulerTerm::VY], R[EulerTerm::P],
                                        R[g], Rgas, kappa, machInf, uInf, pInf, fxR, fyR);
    const CFreal mn = 0.5*((fxL + fxR)*normal[XX] + (fyL + fyR)*normal[YY]);
    m_rFlux[6] += mn*((mn > 0.) ? L[g] : R[g]);
  }
  const CFreal cU = eGammaUpwindCoefficient();
  if (cU > 0.) {
    // upwind non-conservative E_gamma: only the dissipation of the upwind flux, -|rho w_E.n|/2 (gamma_R - gamma_L),
    // with the mean of the two sides' rho w_E.n; the central part is the source -rho w_E.grad(gamma)
    const RealVector& L = m_pData[LEFT];
    const RealVector& R = m_pData[RIGHT];
    const CFuint g = m_updateVarSet->getModel()->getFirstScalarVar(0) + 2;
    const CFreal Rgas = m_updateVarSet->getModel()->getR();
    const CFreal kappa = m_updateVarSet->getModel()->getGamma();
    const CFreal machInf = m_updateVarSet->getModel()->getMachInf();
    const CFreal uInf = m_updateVarSet->getModel()->getVelInf();
    const CFreal pInf = m_updateVarSet->getModel()->getStaticPressInf();
    CFreal fxL = 0., fyL = 0., fxR = 0., fyR = 0.;
    Physics::GammaAlpha::eGammaMassFlux(cU, L[EulerTerm::RHO], L[EulerTerm::VX], L[EulerTerm::VY], L[EulerTerm::P],
                                        L[g], Rgas, kappa, machInf, uInf, pInf, fxL, fyL);
    Physics::GammaAlpha::eGammaMassFlux(cU, R[EulerTerm::RHO], R[EulerTerm::VX], R[EulerTerm::VY], R[EulerTerm::P],
                                        R[g], Rgas, kappa, machInf, uInf, pInf, fxR, fyR);
    const CFreal mn = 0.5*((fxL + fxR)*normal[XX] + (fyL + fyR)*normal[YY]);
    m_rFlux[6] -= 0.5*std::fabs(mn)*(R[g] - L[g]);
  }
  return m_rFlux;
}

  } // namespace FluxReconstructionMethod
} // namespace COOLFluiD
