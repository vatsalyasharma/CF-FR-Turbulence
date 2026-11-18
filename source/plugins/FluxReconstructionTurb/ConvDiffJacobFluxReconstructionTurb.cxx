// Copyright (C) 2019 KU Leuven, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Framework/MethodCommandProvider.hh"
#include "FluxReconstructionTurb/GReKODiffusiveUpdate.hh"
#include "Framework/CFSide.hh"

#include "FluxReconstructionTurb/ConvDiffJacobFluxReconstructionTurb.hh"
#include "FluxReconstructionTurb/TurbWallDistance.hh"
#include "FluxReconstructionNavierStokes/FluxReconstructionNavierStokes.hh"
#include "FluxReconstructionMethod/FluxReconstructionElementData.hh"
#include "GammaAlpha/GammaAlphaArtificialDiffusion.hh"
#include "NavierStokes/NavierStokesVarSet.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Physics::NavierStokes;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider< ConvDiffJacobFluxReconstructionTurb,
		       FluxReconstructionSolverData,
		       FluxReconstructionNavierStokesModule >
convDiffRHSJacobTurbFluxReconstructionProvider("ConvDiffRHSJacobTurb");

//////////////////////////////////////////////////////////////////////////////

ConvDiffJacobFluxReconstructionTurb::ConvDiffJacobFluxReconstructionTurb(const std::string& name) :
  ConvDiffJacobFluxReconstructionNS(name),
  socket_wallDistance("wallDistance"),
  m_closestSolToFlxIdx(CFNULL)
{
}

//////////////////////////////////////////////////////////////////////////////

void ConvDiffJacobFluxReconstructionTurb::configure ( Config::ConfigArgs& args )
{
  ConvDiffJacobFluxReconstructionNS::configure(args);
}

//////////////////////////////////////////////////////////////////////////////

std::vector< Common::SafePtr< BaseDataSocketSink > >
    ConvDiffJacobFluxReconstructionTurb::needsSockets()
{
  std::vector< Common::SafePtr< BaseDataSocketSink > > result = ConvDiffJacobFluxReconstructionNS::needsSockets();

  result.push_back(&socket_wallDistance);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

void ConvDiffJacobFluxReconstructionTurb::computeWaveSpeedUpdates(std::vector<CFreal>& waveSpeedUpd)
{
  if(!detail::isGReKLogO()){ConvDiffJacobFluxReconstructionNS::computeWaveSpeedUpdates(waveSpeedUpd);return;}
  ConvDiffJacobFluxReconstruction::computeWaveSpeedUpdates(waveSpeedUpd);
  auto wall=socket_wallDistance.getDataHandle();
  for(CFuint side=0;side<2;++side)for(CFuint f=0;f<m_nbrFaceFlxPnts;++f){
    m_diffusiveVarSetNS->setWallDistance(wallDistanceAtFlxPnt(wall,*m_states[side],*m_closestSolToFlxIdx,(*m_faceFlxPntConnPerOrient)[m_orient][side][f]));
    const CFreal nu=detail::greDiffusivity(*m_diffusiveVarSetNS,*m_cellStatesFlxPnt[side][f],m_cellGradFlxPnt[side][f]);
    waveSpeedUpd[side]+=nu*m_faceJacobVecAbsSizeFlxPnts[f]*m_faceJacobVecAbsSizeFlxPnts[f]*(*m_faceIntegrationCoefs)[f]*m_cflConvDiffRatio/m_cellVolume[side];
  }
}

void ConvDiffJacobFluxReconstructionTurb::prepareSolPntFluxComputation(const CFuint stateID)
{
  DataHandle< CFreal > wallDist = socket_wallDistance.getDataHandle();

  m_diffusiveVarSetNS->setWallDistance(wallDist[stateID]);
  Physics::GammaAlpha::setGammaArtificialDiffusionAtState(stateID);

  prepareFluxComputation();
}

//////////////////////////////////////////////////////////////////////////////

void ConvDiffJacobFluxReconstructionTurb::prepareFlxPntFluxComputation(const CFuint iFlx)
{
  DataHandle< CFreal > wallDist = socket_wallDistance.getDataHandle();

  // average of the two cells, each taken at its solution point closest to the flux point
  const CFreal wallDistL = wallDistanceAtFlxPnt(wallDist,*m_states[LEFT],*m_closestSolToFlxIdx,
                                                (*m_faceFlxPntConnPerOrient)[m_orient][LEFT][iFlx]);
  const CFreal wallDistR = wallDistanceAtFlxPnt(wallDist,*m_states[RIGHT],*m_closestSolToFlxIdx,
                                                (*m_faceFlxPntConnPerOrient)[m_orient][RIGHT][iFlx]);

  m_diffusiveVarSetNS->setWallDistance(0.5*(wallDistL + wallDistR));

  // EGammaArtificialDiffusion of the gamma-alpha source term (zero otherwise): mean of the same two solution points
  Physics::GammaAlpha::setGammaArtificialDiffusionAtFace(
    (*m_states[LEFT])[(*m_closestSolToFlxIdx)[(*m_faceFlxPntConnPerOrient)[m_orient][LEFT][iFlx]]]->getLocalID(),
    (*m_states[RIGHT])[(*m_closestSolToFlxIdx)[(*m_faceFlxPntConnPerOrient)[m_orient][RIGHT][iFlx]]]->getLocalID());

  prepareFluxComputation();
}

//////////////////////////////////////////////////////////////////////////////

void ConvDiffJacobFluxReconstructionTurb::setup()
{
  CFAUTOTRACE;

  // setup parent class
  ConvDiffJacobFluxReconstructionNS::setup();

  // get the local FR data
  vector< FluxReconstructionElementData* >& frLocalData = getMethodData().getFRLocalData();
  cf_assert(frLocalData.size() > 0);
  // for now, there should be only one type of element
  cf_assert(frLocalData.size() == 1);

  m_closestSolToFlxIdx = frLocalData[0]->getClosestSolToFlxIdx();
}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod

}  // namespace COOLFluiD
