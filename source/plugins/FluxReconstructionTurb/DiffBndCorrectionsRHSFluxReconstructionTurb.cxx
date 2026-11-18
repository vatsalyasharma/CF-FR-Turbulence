// Copyright (C) 2019 KU Leuven, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Framework/MethodCommandProvider.hh"
#include "FluxReconstructionTurb/GReKODiffusiveUpdate.hh"

#include "FluxReconstructionTurb/DiffBndCorrectionsRHSFluxReconstructionTurb.hh"
#include "FluxReconstructionTurb/TurbWallDistance.hh"
#include "FluxReconstructionNavierStokes/FluxReconstructionNavierStokes.hh"
#include "FluxReconstructionMethod/FluxReconstructionElementData.hh"
#include "GammaAlpha/GammaAlphaArtificialDiffusion.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Physics::NavierStokes;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider< DiffBndCorrectionsRHSFluxReconstructionTurb,
		       FluxReconstructionSolverData,
		       FluxReconstructionNavierStokesModule >
diffBndCorrectionsRHSTurbFluxReconstructionProvider("DiffBndCorrectionsRHSTurb");

//////////////////////////////////////////////////////////////////////////////

DiffBndCorrectionsRHSFluxReconstructionTurb::DiffBndCorrectionsRHSFluxReconstructionTurb(const std::string& name) :
  DiffBndCorrectionsRHSFluxReconstructionNS(name),
  socket_wallDistance("wallDistance"),
  m_closestSolToFlxIdx(CFNULL),
  m_navierStokesVarSet(CFNULL)
{
}

//////////////////////////////////////////////////////////////////////////////

DiffBndCorrectionsRHSFluxReconstructionTurb::~DiffBndCorrectionsRHSFluxReconstructionTurb()
{
}

//////////////////////////////////////////////////////////////////////////////

void DiffBndCorrectionsRHSFluxReconstructionTurb::setup()
{
  CFAUTOTRACE;

  // setup parent class
  DiffBndCorrectionsRHSFluxReconstructionNS::setup();

  // get the local FR data
  vector< FluxReconstructionElementData* >& frLocalData = getMethodData().getFRLocalData();
  cf_assert(frLocalData.size() > 0);
  // for now, there should be only one type of element
  cf_assert(frLocalData.size() == 1);

  m_closestSolToFlxIdx = frLocalData[0]->getClosestSolToFlxIdx();

  m_navierStokesVarSet = m_diffusiveVarSet.d_castTo< NavierStokesVarSet >();
}

//////////////////////////////////////////////////////////////////////////////

void DiffBndCorrectionsRHSFluxReconstructionTurb::unsetup()
{
  CFAUTOTRACE;

  // unsetup parent class
  DiffBndCorrectionsRHSFluxReconstructionNS::unsetup();
}

//////////////////////////////////////////////////////////////////////////////

std::vector< Common::SafePtr< BaseDataSocketSink > >
    DiffBndCorrectionsRHSFluxReconstructionTurb::needsSockets()
{
  std::vector< Common::SafePtr< BaseDataSocketSink > > result = DiffBndCorrectionsRHSFluxReconstructionNS::needsSockets();

  result.push_back(&socket_wallDistance);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

void DiffBndCorrectionsRHSFluxReconstructionTurb::computeWaveSpeedUpdates(CFreal& waveSpeedUpd)
{
  if(!detail::isGReKLogO()){DiffBndCorrectionsRHSFluxReconstructionNS::computeWaveSpeedUpdates(waveSpeedUpd);return;}
  waveSpeedUpd=0.;auto wall=socket_wallDistance.getDataHandle();
  for(CFuint f=0;f<m_nbrFaceFlxPnts;++f){
    m_navierStokesVarSet->setWallDistance(wallDistanceAtFlxPnt(wall,*m_cellStates,*m_closestSolToFlxIdx,(*m_faceFlxPntConn)[m_orient][f]));
    const CFreal nu=detail::greDiffusivity(*m_navierStokesVarSet,*m_cellStatesFlxPnt[f],m_cellGradFlxPnt[f]);
    waveSpeedUpd+=nu*m_faceJacobVecAbsSizeFlxPnts[f]*m_faceJacobVecAbsSizeFlxPnts[f]*(*m_faceIntegrationCoefs)[f]*m_cflConvDiffRatio/m_cellVolume;
  }
}

void DiffBndCorrectionsRHSFluxReconstructionTurb::prepareFlxPntFluxComputation(const CFuint iFlx)
{
  DataHandle< CFreal > wallDist = socket_wallDistance.getDataHandle();

  // the interior cell, at its solution point closest to the flux point
  m_navierStokesVarSet->setWallDistance(wallDistanceAtFlxPnt(wallDist,*m_cellStates,*m_closestSolToFlxIdx,
                                                             (*m_faceFlxPntConn)[m_orient][iFlx]));

  // EGammaArtificialDiffusion of the gamma-alpha source term: no artificial gamma flux through boundaries
  Physics::GammaAlpha::setGammaArtificialDiffusionZero();

  prepareFluxComputation();
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluxReconstructionMethod

} // namespace COOLFluiD
