//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//
//        This file is part of E-Cell Simulation Environment package
//
//                Copyright (C) 2006-2009 Keio University
//
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//
//
// E-Cell is free software; you can redistribute it and/or
// modify it under the terms of the GNU General Public
// License as published by the Free Software Foundation; either
// version 2 of the License, or (at your option) any later version.
// 
// E-Cell is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// See the GNU General Public License for more details.
// 
// You should have received a copy of the GNU General Public
// License along with E-Cell -- see the file COPYING.
// If not, write to the Free Software Foundation, Inc.,
// 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
// 
//END_HEADER
//
// written by Satya Arjunan <satya.arjunan@gmail.com>
// E-Cell Project, Institute for Advanced Biosciences, Keio University.
//


#include "PolymerizationProcess.hpp"
#include "SpatiocyteSpecies.hpp"

LIBECS_DM_INIT(PolymerizationProcess, Process);

void PolymerizationProcess::pushNewBend(Subunit* aSubunit, double aBendAngle)
{
  Point& aRefPoint(aSubunit->subunitPoint);
  Point aPoint;
  Bend aBend;
  while(aBendAngle > M_PI)
    {
      aBendAngle -= 2*M_PI;
    }
  while(aBendAngle < -M_PI)
    {
      aBendAngle += 2*M_PI;
    }
  aBend.angle = aBendAngle;
  if(aRefPoint.z < theOriZ)
    {
      aBendAngle = aBendAngle-M_PI;
    }
  double x(aRefPoint.x);
  if(x > theMaxX)
    { 
      x = x-theMaxX;
    }
  else if(x < theMinX)
    {
      x = x-theMinX;
    }
  else
    {
      x = 0;
    }     
  double currX[3];
  currX[0] = aRefPoint.x-theMinX-((theMaxX-theMinX)/2);
  currX[1] = aRefPoint.y-theOriY;
  currX[2] = aRefPoint.z-theOriZ;
  double tmpDcm[9];
  getOneDcm(tmpDcm);
  rotXrotY(tmpDcm, atan2(currX[1],sqrt(x*x+currX[2]*currX[2])),
           -atan2(x,-currX[2]));
  getCylinderDcm(-aBendAngle, -aBendAngle-CylinderYaw, aBend.cylinderDcm);
  getSphereDcm(0, -SphereYaw, aBend.sphereDcm);
  double aDcm[9];
  if(getLocation(aRefPoint.x) == CYLINDER)
    {
      getCylinderDcm(0, M_PI-aBendAngle, aDcm);
      dcmXdcm(aDcm, tmpDcm, aDcm);
      pinStep(currX, aBend.cylinderDcm, aDcm, aBend.dcm);
    }
  else
    {
      getSphereDcm(0, M_PI-aBendAngle, aDcm);
      dcmXdcm(aDcm, theInitSphereDcm, aDcm);
      dcmXdcm(aDcm, tmpDcm, aDcm);
      pinStep(currX, aBend.sphereDcm, aDcm, aBend.dcm);
    }
  dcmXdcm(aBend.dcm, aDcm, aBend.dcm);
  aPoint.x = aRefPoint.x + theMonomerLength*aBend.dcm[0];
  aPoint.y = aRefPoint.y + theMonomerLength*aBend.dcm[1];
  aPoint.z = aRefPoint.z + theMonomerLength*aBend.dcm[2];
  aSubunit->targetPoints.push_back(aPoint);
  aSubunit->targetBends.push_back(aBend);
}

void PolymerizationProcess::pushJoinBend(Subunit* aSubunit, Subunit* refSubunit,
                                         unsigned int aBendIndex)
{
  Point& aRefPoint(refSubunit->targetPoints[aBendIndex]);
  Bend& aRefBend(refSubunit->targetBends[aBendIndex]);
  Bend aBend(aRefBend);
  Point aPoint;
  double currX[3];
  currX[0] = aRefPoint.x-theMinX-((theMaxX-theMinX)/2);
  currX[1] = aRefPoint.y-theOriY;
  currX[2] = aRefPoint.z-theOriZ;
  if(getLocation(aRefPoint.x) == CYLINDER)
    {
      pinStep(currX, aBend.cylinderDcm, aRefBend.dcm, aBend.dcm);
    }
  else
    {
      pinStep(currX, aBend.sphereDcm, aRefBend.dcm, aBend.dcm);
    }
  dcmXdcm(aBend.dcm, aRefBend.dcm, aBend.dcm);
  aPoint.x = aRefPoint.x + theMonomerLength*aBend.dcm[0];
  aPoint.y = aRefPoint.y + theMonomerLength*aBend.dcm[1];
  aPoint.z = aRefPoint.z + theMonomerLength*aBend.dcm[2];
  aBend.angle = aRefBend.angle; //maybe can remove this line
  aSubunit->targetPoints.push_back(aPoint);
  aSubunit->targetBends.push_back(aBend);
}

Bend* PolymerizationProcess::getNewReverseBend(Point* aRefPoint,
                                               Point* aPoint, Bend* aRefBend)
{ 
  Bend* aBend(getReverseBend(aRefBend));
  double currX[3];
  currX[0] = aRefPoint->x-theMinX-((theMaxX-theMinX)/2);
  currX[1] = aRefPoint->y-theOriY;
  currX[2] = aRefPoint->z-theOriZ;
  double aDcm[9];
  if( getLocation( aRefPoint->x ) == CYLINDER )
    {
      pinStep( currX, aBend->cylinderDcm, aBend->dcm, aDcm );
    }
  else
    {
      pinStep( currX, aBend->sphereDcm, aBend->dcm, aDcm );
    }
  dcmXdcm( aDcm, aBend->dcm, aBend->dcm );
  aPoint->x = aRefPoint->x + theMonomerLength*aBend->dcm[0];
  aPoint->y = aRefPoint->y + theMonomerLength*aBend->dcm[1];
  aPoint->z = aRefPoint->z + theMonomerLength*aBend->dcm[2];
  return aBend;
}

Bend* PolymerizationProcess::getReverseBend( Bend* aRefBend )
{ 
  Bend* aBend( new Bend );
  reverseDcm(aRefBend->dcm, aBend->dcm);
  reverseYpr(aRefBend->sphereDcm, aBend->sphereDcm);
  reverseYpr(aRefBend->cylinderDcm, aBend->cylinderDcm);
  aBend->angle = aRefBend->angle-M_PI;
  while( aBend->angle > M_PI )
    {
      aBend->angle -= 2*M_PI;
    }
  while( aBend->angle < -M_PI )
    {
      aBend->angle += 2*M_PI;
    }
  return aBend;
}

Point PolymerizationProcess::getNextPoint(Point* aRefPoint,
                                          Bend* aRefBend ) 
{
  double currX[3];
  double aDcm[9];
  currX[0] = aRefPoint->x-theMinX-((theMaxX-theMinX)/2);
  currX[1] = aRefPoint->y-theOriY;
  currX[2] = aRefPoint->z-theOriZ;
  if( getLocation( aRefPoint->x ) == CYLINDER )
    {
      pinStep( currX, aRefBend->cylinderDcm,
               aRefBend->dcm, aDcm );
    }
  else
    {
      pinStep( currX, aRefBend->sphereDcm, aRefBend->dcm,
               aDcm );
    }
  dcmXdcm( aDcm, aRefBend->dcm, aDcm );
  Point aPoint;
  aPoint.x = aRefPoint->x + theMonomerLength*aDcm[0];
  aPoint.y = aRefPoint->y + theMonomerLength*aDcm[1];
  aPoint.z = aRefPoint->z + theMonomerLength*aDcm[2];
  return aPoint;
}

void PolymerizationProcess::getCylinderDcm(double b1, double b2, double *dcm)
{
  double sb1(sin(b1));
  double cb1(cos(b1));
  double sb2(sin(b2));
  double cb2(cos(b2));
  double x(asin(theMonomerLength*sb1/(2*theRadius))+
           asin(theMonomerLength*sb2/(2*theRadius)));
  double cx(cos(x));
  double sx(sin(x));
  double ypr[3];
  ypr[0] = atan2(-sb1*cb1+cb1*cx*sb2, cb1*cb2+sb1*cx*sb2);
  ypr[1] = asin(-sx*sb2);
  ypr[2] = atan2(sx*cb2, cx);
  ypr2dcm( ypr, dcm );
}  

void PolymerizationProcess::getSphereDcm(double b1, double b2, double *dcm)
{
  double x(asin(theMonomerLength/(2*theRadius)));
  double cx(cos(x));
  double sx(sin(x));
  double a(b2/cx);
  double ca(cos(a));
  double sa(sin(a));
  double ypr[3];
  ypr[0] = atan2(sa*cx,ca*cx*cx-sx*sx);
  ypr[1] = asin(-ca*cx*sx-sx*cx);
  ypr[2] = atan2(-sa*sx,-ca*sx*sx+cx*cx);
  ypr2dcm( ypr, dcm );
}  

void PolymerizationProcess::initSphereDcm()
{
  double rot[9];
  double x(0);
  double y(0);
  double z(-theRadius);
  getSphereDcm(0, 0, theInitSphereDcm); //dcm is const because angle is const
  rotY(rot, -atan2(theMonomerLength*theInitSphereDcm[0],
                   theRadius-theMonomerLength*theInitSphereDcm[2]));
  rotate(rot,&x,&y,&z);
  double thetaY(-atan2(x,-z)+atan2((theMonomerLength-1)/2,theRadius));
  // makes dcm start at correct position (-5, 0, -79) so that the first
  // coord after that is (0, 0, -r) given x,y,z = (-5,0,-79)
  getOneDcm(rot);
  rotY(rot, thetaY);
  dcmXdcm(theInitSphereDcm, rot, theInitSphereDcm);
}

void PolymerizationProcess::pinStep(double* currX, double *fixedDcm,
                                    double* currDcm, double* resDcm)
{
  double nextDcm[9];
  double nextX[3];
  dcmXdcm( fixedDcm, currDcm, nextDcm );
  nextX[0] = currX[0] + theMonomerLength*nextDcm[0];
  nextX[1] = currX[1] + theMonomerLength*nextDcm[1];
  nextX[2] = currX[2] + theMonomerLength*nextDcm[2];
  pinPoint( currX, nextX, 'e', theRadius, (theMaxX-theMinX)/2, nextDcm);
  dcmXdcmt( nextDcm, currDcm, resDcm );
  double ypr[3];
  dcm2ypr(resDcm, ypr);
  ypr2dcm(ypr, resDcm);
}

void PolymerizationProcess::initializeThird()
{
  theMinX = C->getWestPoint().x;
  theMaxX = C->getEastPoint().x;
  theOriY = C->getWestPoint().y;
  theOriZ = C->getWestPoint().z;
  theRadius = C->getRadius();
  theBendIndexA = A->getBendIndex(BendAngle);
  theBendIndexB = B->getBendIndex(BendAngle);
  initSphereDcm();  
  //Create polymer bends for each subunit of every polymer species:
  initSubunits(A);
  if(A != B)
    {
      initSubunits(B);
    }
  if(A != C && B != C)
    {
      initSubunits(C);
    }
  if(D && A != D && B != D && C != D)
    {
      initSubunits(D);
    }
}

void PolymerizationProcess::initSubunits(Species* aSpecies)
{
  if(aSpecies->getIsPolymer() && !aSpecies->getIsSubunitInitialized())
    {
      std::vector<Voxel*>& molecules(aSpecies->getMolecules());
      for(std::vector<Voxel*>::iterator i(molecules.begin());
          i != molecules.end(); ++i)
        {
          initSubunit(*i, aSpecies);
        }
      aSpecies->setIsSubunitInitialized();
    }
}

void PolymerizationProcess::initSubunit(Voxel* aMolecule, Species* aSpecies)
{
  const std::vector<double>& bendAngles(aSpecies->getBendAngles());
  Subunit* aSubunit(aMolecule->subunit);
  aSubunit->bendSize = bendAngles.size();
  aSubunit->voxel = aMolecule;
  //Use the surfacePoint as the subunitPoint since the voxel is the
  //origin of the polymer:
  aSubunit->subunitPoint = aSubunit->surfacePoint;
  for(std::vector<double>::const_iterator j(bendAngles.begin());
      j != bendAngles.end(); ++j)
    {
      pushNewBend(aSubunit, *j);
    } 
  aSubunit->targetVoxels.resize(aSubunit->bendSize);
  aSubunit->sourceVoxels.resize(aSubunit->bendSize);
  aSubunit->sharedLipids.resize(aSubunit->bendSize);
  aSubunit->tmpVoxels.resize(aSubunit->bendSize);
  aSubunit->boundBends.resize(aSubunit->bendSize);
  for(unsigned int i(0); i != aSubunit->bendSize; ++i)
    {
      aSubunit->targetVoxels[i] = NULL;
      aSubunit->sourceVoxels[i] = NULL;
      aSubunit->sharedLipids[i] = NULL;
      aSubunit->tmpVoxels[i] = NULL;
      aSubunit->boundBends[i] = false;
    }
  //Keep track of the continuous points represented by the voxel.
  //Use its subunit structure to do this:
  addContPoint(aSubunit, &aSubunit->subunitPoint);
}


void PolymerizationProcess::finalizeReaction()
{
  DiffusionInfluencedReactionProcess::finalizeReaction();
}

bool PolymerizationProcess::isInterrupting(Process* aProcess)
{
  if(aProcess->getPropertyInterface().getClassName() ==
     "PolymerFragmentationProcess") 
    {
      PolymerFragmentationProcess* aDepolymerizeProcess(
           dynamic_cast<PolymerFragmentationProcess*>(aProcess));
      aDepolymerizeProcess->setPolymerizeProcess(this);
    }
  return ReactionProcess::isInterrupting(aProcess);
}

void PolymerizationProcess::initJoinSubunit(Voxel* aMolecule, Species* aSpecies,
                                            Subunit* refSubunit)
{
  const std::vector<double>& bendAngles(aSpecies->getBendAngles());
  Subunit* aSubunit(aMolecule->subunit);
  aSubunit->bendSize = bendAngles.size();
  aSubunit->voxel = aMolecule;
  aSubunit->subunitPoint = refSubunit->targetPoints[theBendIndexA];
  for(unsigned int i(0); i != bendAngles.size(); ++i)
    {
      if(i == theBendIndexB)
        {
          pushJoinBend(aSubunit, refSubunit, theBendIndexA);
        }
      else
        {
          pushNewBend(aSubunit, bendAngles[i]);
        }
    } 
  aSubunit->targetVoxels.resize(aSubunit->bendSize);
  aSubunit->sourceVoxels.resize(aSubunit->bendSize);
  aSubunit->sharedLipids.resize(aSubunit->bendSize);
  aSubunit->tmpVoxels.resize(aSubunit->bendSize);
  aSubunit->boundBends.resize(aSubunit->bendSize);
  for(unsigned int i(0); i != aSubunit->bendSize; ++i)
    {
      aSubunit->targetVoxels[i] = NULL;
      aSubunit->sourceVoxels[i] = NULL;
      aSubunit->sharedLipids[i] = NULL;
      aSubunit->tmpVoxels[i] = NULL;
      aSubunit->boundBends[i] = false;
    }
  //Keep track of the continuous points represented by the voxel.
  //Use its subunit structure to do this:
  addContPoint(aSubunit, &aSubunit->subunitPoint);
}



//A voxel may represent different continuous points to surrounding subunits.
//The contPoints contains the list of the continuous points.
//Although it is more intuitive to put the contPoints under the voxel
//structure, I have placed it under the subunit structure to save memory space,
//since the subunit structure is only created when required.
//Each continuous point of the subunit is used by one surroundings subunit, so
//there shouldn't be duplicates by the same surrounding subunit.
//Each continuous point of the subunit may be used by more than one surrounding
//subunit. As such, we keep track of their number using contPointSize.
void PolymerizationProcess::addContPoint(Subunit* aSubunit, Point* aPoint)
{
  //Add a duplicate contPoint:
  if(aSubunit->contPoints.size())
    {
      for(unsigned int i(0); i != aSubunit->contPoints.size(); ++i)
        {
          if(getDistance(&aSubunit->contPoints[i], aPoint) < 0.1)
            { 
              ++aSubunit->contPointSize[i];
              return;
            }
        }
    }
  //Add a new unique contPoint:
  aSubunit->contPoints.push_back(*aPoint);
  aSubunit->contPointSize.resize(aSubunit->contPointSize.size()+1);
  ++aSubunit->contPointSize.back();
}


void PolymerizationProcess::removeContPoint(Subunit* aSubunit,  Point* aPoint)
{
  if(aSubunit->contPoints.size() == 1)
    {
      --aSubunit->contPointSize[0];
      if(!aSubunit->contPointSize[0])
        {
          aSubunit->contPoints.clear();
          aSubunit->contPointSize.clear();
        }
      return;
    }
  //If the size of contPoints is more than 1:
  else
    { 
      for(unsigned int i(0); i != aSubunit->contPoints.size(); ++i)
        {
          if(getDistance(&aSubunit->contPoints[i], aPoint) < 0.1)
            { 
              --aSubunit->contPointSize[i];
              //If the size of the continuous point is zero, we need to remove
              //the point from the contPoints, and update the contPointSize
              //list size:
              if(!aSubunit->contPointSize[i])
                {
                  //Remove the continuous point:
                  aSubunit->contPoints[i] = aSubunit->contPoints.back();
                  aSubunit->contPoints.pop_back();
                  //Update the size:
                  aSubunit->contPointSize[i] = aSubunit->contPointSize.back();
                  aSubunit->contPointSize.pop_back();
                }
              return;
            }
        }
    }
  std::cout << "error in remove contPoint at time:" << getStepper()->getCurrentTime() << std::endl;
}

void PolymerizationProcess::resetSubunit(Subunit* aSubunit)
{
  //There are at least two identical contPoints if the subunit
  //is the target of other subunit:
  //1. the subunit's point
  //2. the point set when getTargetVoxel was called to get the voxel
  //for this subunit.
  removeContPoint(aSubunit, &aSubunit->subunitPoint);
  for(unsigned int i(0); i != aSubunit->bendSize; ++i)
    {
      if(aSubunit->boundBends[i])
        {
          Subunit* boundSubunit(aSubunit->targetVoxels[i]->subunit);
          std::vector<Voxel*>& sourceVoxels(boundSubunit->sourceVoxels);
          for(unsigned int j(0); j != sourceVoxels.size(); ++j)
            {
              if(sourceVoxels[j] == aSubunit->voxel)
                {
                  sourceVoxels[j] = NULL;
                }
            }
          removeContPoint(aSubunit->targetVoxels[i]->subunit,
                          &aSubunit->targetPoints[i]); 
        }
      removeLipid(aSubunit, i);
      if(aSubunit->sourceVoxels[i] != NULL)
        {
          Subunit* boundSubunit(aSubunit->sourceVoxels[i]->subunit);
          std::vector<Voxel*>& targetVoxels(boundSubunit->targetVoxels);
          std::vector<bool>& boundBends(boundSubunit->boundBends);
          for(unsigned int j(0); j != targetVoxels.size(); ++j)
            {
              if(targetVoxels[j] == aSubunit->voxel)
                {
                  boundBends[j] = false;
                  break; 
                }
            }
        }
    }
  aSubunit->targetPoints.resize(0);
  aSubunit->targetVoxels.resize(0);
  aSubunit->sourceVoxels.resize(0);
  aSubunit->targetBends.resize(0);
  aSubunit->tmpVoxels.resize(0);
  aSubunit->boundBends.resize(0);
  //Do not remove aSubunit->contPoints and aSubunit->contPointSize because
  //they hold persistent information of the voxel, not the subunit.
  //std::cout << "reset done" << std::endl;
}

void PolymerizationProcess::removeLipid(Subunit* aSubunit, 
                                        unsigned int aBendIndex)
{
  if(aSubunit->sharedLipids[aBendIndex] != NULL)
    { 
      for(unsigned int i(0); i!=aSubunit->bendSize; ++i)
        {
          //If the shared lipid is still used by other bends:
          if(i!=aBendIndex && aSubunit->sharedLipids[i] != NULL &&
             aSubunit->sharedLipids[i] == aSubunit->sharedLipids[aBendIndex])
            {
              //Remove the one of the contPoints in the subunit of the shared
              //lipid:
              removeContPoint(aSubunit->sharedLipids[aBendIndex]->subunit,
                              &aSubunit->subunitPoint);
              aSubunit->sharedLipids[aBendIndex] = NULL;
              return;
            }
        } 
      aSubunit->sharedLipids[aBendIndex]->subunit->voxel =
        aSubunit->sharedLipids[aBendIndex];
      removeContPoint(aSubunit->sharedLipids[aBendIndex]->subunit,
                      &aSubunit->subunitPoint);
      aSubunit->sharedLipids[aBendIndex]->id = 
        id2species(aSubunit->sharedLipids[aBendIndex]->id)->getVacantID(); 
      aSubunit->sharedLipids[aBendIndex] = NULL;
    }
}

//Do the reaction A + B -> C + D. So that A <- C and B <- D.
//We need to consider that the source molecule can be either A or B.
//If A and C belong to the same compartment, A <- C.
//Otherwise, find a vacant adjoining voxel of A, X which is the same compartment
//as C and X <- C.
//Similarly, if B and D belong to the same compartment, B <- D.
//Otherwise, find a vacant adjoining voxel of C, Y which is the same compartment
//as D and Y <- D.
//We need to consider that the source molecule can be either A or B.
bool PolymerizationProcess::react(Voxel* moleculeB, Voxel** target)
{
  Voxel* moleculeA(*target);
  //moleculeB is the source molecule. It will be soft-removed (id kept intact)
  //by the calling Species if this method returns true.
  //moleculeA is the target molecule,  it will also be soft-removed by the
  //calling Species.
  //Usually the source molecule, moleculeB is not a protomer since it is
  //diffusing and the target molecule, moleculeA is a static protomer.
  //First let us make sure moleculeA and moleculeB belong to the
  //correct species.
  if(moleculeA->id != A->getID())
    {
      Voxel* tempA(moleculeA);
      moleculeA = moleculeB;
      moleculeB = tempA;
    }
  //C && D must be protomers:
  if(C && D)
    {
      //Dimerization reaction:
      if(!A->getIsPolymer() && !B->getIsPolymer())
        {
          initSubunit(moleculeA, C);
          Subunit* subunitA(moleculeA->subunit);
          if(setTargetVoxels(subunitA))
            { 
              Voxel* moleculeD(subunitA->targetVoxels[theBendIndexA]);
              if(moleculeD == moleculeB || 
                 id2species(moleculeD->id)->getIsLipid())
                { 
                  subunitA->boundBends[theBendIndexA] = true;
                  moleculeB->id = B->getVacantID();
                  initJoinSubunit(moleculeD, D, subunitA); 
                  moleculeD->subunit->sourceVoxels[theBendIndexB] = moleculeA;
                  C->addMolecule(moleculeA);
                  updateSharedLipidsID(moleculeA);
                  D->addMolecule(moleculeD);
                  return true;
                }
            }
          resetSubunit(subunitA);
        }
      //Polymer elongation reaction:
      //A is the reference polymer subunit:
      else if(A->getIsPolymer() && !B->getIsPolymer())
        {
          //Other molecules can react with a shared molecule by collision.
          //Make sure the moleculeA is not a shared molecule by updating
          //it to the actual molecule pointed by the subunit:
          moleculeA = moleculeA->subunit->voxel;
          Subunit* subunitA(moleculeA->subunit);
          //Only proceed if the target bend is unbound:
          if(!subunitA->boundBends[theBendIndexA])
            {
              removeUnboundTargetVoxels(subunitA);
              //Set up all unbound target voxels. Proceed if at least the
              //targetVoxel of this reaction was found successfully:
              if(setTargetVoxels(subunitA))
                {
                  Voxel* moleculeD(subunitA->targetVoxels[theBendIndexA]);
                  if(moleculeD == moleculeB || 
                     id2species(moleculeD->id)->getIsLipid())
                    { 
                      subunitA->boundBends[theBendIndexA] = true;
                      moleculeB->id = B->getVacantID();
                      initJoinSubunit(moleculeD, D, subunitA); 
                      moleculeD->subunit->sourceVoxels[theBendIndexB] =
                        moleculeA;
                      C->addMolecule(moleculeA);
                      updateSharedLipidsID(moleculeA);
                      D->addMolecule(moleculeD);
                      *target = moleculeA;
                      return true;
                    }
                }
            }
        }
    }
  //Single product polymerization:
  else
    {
      //Depolymerize once react: 
      if(A->getIsPolymer() && !B->getIsPolymer() && !C->getIsPolymer())
        {
          Voxel* moleculeC;
          if(A->getVacantID() != C->getVacantID())
            {
              if(B->getVacantID() != C->getVacantID())
                {
                  moleculeC = C->getRandomAdjoiningVoxel(moleculeA);
                  if(moleculeC == NULL)
                    {
                      moleculeC = C-> getRandomAdjoiningVoxel(moleculeB);
                      if(moleculeC == NULL)
                        {
                          //std::cout << "unavailable" << std::endl;
                          return false;
                        }
                    }
                }
              else
                {
                  moleculeC = moleculeB;
                }
            }
          else
            {
              moleculeC = moleculeA;
            }
          moleculeA->id = A->getVacantID();
          moleculeB->id = B->getVacantID();
          C->addMolecule(moleculeC);
          resetSubunit(moleculeA->subunit);
          //std::cout << "left" << std::endl;
          return true;
        }
    }
  return false;
}

void PolymerizationProcess::updateSharedLipidsID(Voxel* aMolecule)
{
  std::vector<Voxel*>& sharedLipids(aMolecule->subunit->sharedLipids);
  for(unsigned int i(0); i != sharedLipids.size(); ++i)
    {
      if(sharedLipids[i])
        {
          sharedLipids[i]->id = aMolecule->id;
        }
    }
}

void PolymerizationProcess::removeUnboundTargetVoxels(Subunit* aSubunit)
{
  for(unsigned int i(0); i != aSubunit->bendSize; ++i)
    {
      if(!aSubunit->boundBends[i] && aSubunit->targetVoxels[i])
        {
          removeContPoint(aSubunit->targetVoxels[i]->subunit, 
                          &aSubunit->targetPoints[i]); 
          removeLipid(aSubunit, i);
          aSubunit->targetVoxels[i] = NULL;
        }
    }
}

Voxel* PolymerizationProcess::setNewTargetVoxel(Subunit* aSubunit, 
                                                unsigned int aBendIndex)
{
  double aDist(setImmediateTargetVoxel(aSubunit, aBendIndex));
  //If we have not found a targetVoxel from immediate voxels of the subunit:
  if(aDist)
    {
      //If we have not found a targetVoxel from all immediate and
      //extended voxels of the subunit:
      if(!setExtendedTargetVoxel(aSubunit, aBendIndex, aDist)) 
        {
          aSubunit->targetVoxels[aBendIndex] = NULL;
        }
    }
  return aSubunit->targetVoxels[aBendIndex];
}

//Returns false if the targetVoxel[theBendIndexA] is not successful 
//otherwise it returns true:
bool PolymerizationProcess::setTargetVoxels(Subunit* aSubunit)
{
  bool retValue(true);
  std::vector<unsigned int>  nearestBends;
  for(unsigned int i(0); i != aSubunit->bendSize; ++i)
    {
      int j(setExistingTargetVoxel(aSubunit, i));
      if(j == 0)
        {
          nearestBends.push_back(i);
        }
      else if(j == 2 && i == theBendIndexA)
        {
          retValue = false;
        }
    } 
  std::vector<std::vector<Voxel*> > targetVoxels; 
  std::vector<std::vector<double> > distances;
  std::vector<std::vector<int> > sharedIndices; 
  std::vector<unsigned int> sizes;
  std::vector<unsigned int> curr;
  targetVoxels.resize(nearestBends.size());
  distances.resize(nearestBends.size());
  sharedIndices.resize(nearestBends.size());
  curr.resize(nearestBends.size());
  for(unsigned int i(0); i != nearestBends.size(); ++i)
    {
      sizes.push_back(sortNearestTargetVoxels(aSubunit, nearestBends[i],
                                              targetVoxels[i],
                                              distances[i],
                                              sharedIndices[i]));
      curr[i] = 0;
      std::cout << "size:" << sizes[i] << std::endl;
    }
  for(unsigned int i(0); i != nearestBends.size(); ++i)
    {
      std::cout << "aBend:" << nearestBends[i] << std::endl;
      for(unsigned int j(0); j != sizes[i]; ++j)
        {
          std::cout << "dist:" << j << " " << distances[i][j] << std::endl;
        }
    }
  std::cout << "before" << std::endl;
  unsigned int j(1);
  while(j != nearestBends.size())
    {
      unsigned int i(0);
      while(i != nearestBends.size())
        {
          if(i+j < nearestBends.size())
            {
              if(sizes[i] > curr[i])
                {
                  if(sizes[i+j] > curr[i+j])
                    {
                      if(targetVoxels[i][curr[i]] ==
                         targetVoxels[i+j][curr[i+j]])
                        {
                          if(sizes[i+j] > curr[i+j]+1)
                            {
                              if(sizes[i] > curr[i]+1)
                                {
                                  if(distances[i][curr[i]] +
                                     distances[i][curr[i+j]+1] < 
                                     distances[i][curr[i]+1] + 
                                     distances[i][curr[i+j]])
                                    {
                                      ++curr[i+j];
                                    }
                                  else
                                    {
                                      ++curr[i];
                                    }
                                }
                              else
                                {
                                  ++curr[i+j];
                                }
                            }
                          else if(sizes[i] > curr[i]+1)
                            {
                              ++curr[i];
                            }
                        }
                    }
                }
            }
          ++i;
        }
      ++j;
    }
  /*
  std::cout << "end" << std::endl;
  for(unsigned int i(0); i != nearestBends.size(); ++i)
    {
      std::cout << "dist:" << distances[i][curr[i]] << std::endl;
      if(!i && nearestBends.size() > 1)
        {
          std::cout << "In:" << curr[i] << " out:" << curr[i+1] << std::endl;
          if(targetVoxels[i][curr[i]] == targetVoxels[i+1][curr[i+1]])
            {
              std::cout << "same" << std::endl;
            }
        }
    }
    */

  for(unsigned int k(0); k != nearestBends.size(); ++k)
    { 
      if(sizes[k])
        {
          //std::cout << "k:" << k << std::endl;
          unsigned int i(nearestBends[k]);
          aSubunit->targetVoxels[i] = targetVoxels[k][curr[k]]; 
          Voxel* aRefVoxel(aSubunit->voxel);
          Point* aRefPoint(&aSubunit->targetPoints[i]);
          if(sharedIndices[k][curr[k]] != -1)
            { 
              //std::cout << "here:" << sharedIndices[k][curr[k]] << std::endl;
              std::vector<Voxel*>& aSharedList((*aRefVoxel->surfaceVoxels)
                                          [SHARED+sharedIndices[k][curr[k]]]);
              double aSharedDist(LARGE_DISTANCE);
              Voxel* aSelectedSharedVoxel(NULL);
              for(unsigned int j(0); j != aSharedList.size(); ++j)
                { 
                  //std::cout << "sharedlist:" << j << std::endl;
                  Voxel* aSharedVoxel(aSharedList[j]);
                  Subunit* aSharedSubunit(aSharedVoxel->subunit);
                  //First find an existing shared voxel which connects
                  //the reference voxel to the extended voxel:
                  if(aSharedSubunit->voxel == aSubunit->voxel)
                    {
                      //std::cout << "break" << std::endl;
                      aSelectedSharedVoxel = aSharedVoxel;
                      break;
                    }
                  //Otherwise find a shared voxel that is unoccupied by a
                  //protomer:
                  else if(aSharedSubunit->contPoints.empty())
                    {
                      //std::cout << "elsse" << std::endl;
                      double aDist(getDistance(aRefPoint,
                                           &aSharedSubunit->surfacePoint));
                      if(aDist < aSharedDist)
                        {
                          //std::cout << "af" << std::endl;
                          aSharedDist = aDist;
                          aSelectedSharedVoxel = aSharedVoxel;
                        }
                    }
                  else
                    {
                      //std::cout << "failure" << std::endl;
                      return false;
                    }
                }
              //If the selected shared voxel is not an existing shared voxel nor
              //a lipid:
              if(!id2species(aSelectedSharedVoxel->id)->getIsLipid() &&
                 aSelectedSharedVoxel->subunit->voxel != aRefVoxel)
                {
                  //std::cout << "fail1 immed"  << std::endl;
                  return false;
                }
                          //std::cout << "afdaff" << std::endl;
              aSelectedSharedVoxel->subunit->voxel = aRefVoxel;
              //Add the subunitPoint of the reference protomer to the 
              //list of contPoints of the selected shared subunit:
              addContPoint(aSelectedSharedVoxel->subunit,
                           &aSubunit->subunitPoint);
              aSubunit->sharedLipids[i] = aSelectedSharedVoxel;
            }
                          //std::cout << "afllll" << std::endl;
          addContPoint(aSubunit->targetVoxels[i]->subunit, aRefPoint);
        }
    }
  //std::cout << "Herell" << std::endl;
  return retValue;
}


int PolymerizationProcess::setExistingTargetVoxel(Subunit* aRefSubunit,
                                                  unsigned int aBendIndex) 
{
  Voxel* aRefVoxel(aRefSubunit->voxel);
  Point* aRefPoint(&aRefSubunit->targetPoints[aBendIndex]);
  std::vector<Voxel*>& immediateSurface(aRefVoxel->surfaceVoxels[0][IMMEDIATE]);
  //Check the immediate 6 (usually) surface voxels adjoining the voxel of the
  //reference subunit:
  for(unsigned int i(0); i != immediateSurface.size(); ++i)
    {
      Voxel* aVoxel(immediateSurface[i]);
      Subunit* aSubunit(aVoxel->subunit);
      if(!aSubunit->contPoints.empty())
        {
          for(unsigned int j(0); j!=aSubunit->contPoints.size(); ++j)
            {
              double aDist(getDistance(aRefPoint, &aSubunit->contPoints[j]));
              //If there is an existing target point:
              if(aDist < 0.1)
                {
                  //If it is a shared voxel:
                  if(aSubunit->voxel != aVoxel)
                    {
                      //Lets use the actual voxel if it is an immediate voxel
                      //of the reference voxel:
                      if(std::find(immediateSurface.begin(), immediateSurface.end(),
                              aSubunit->voxel) != immediateSurface.end() )
                        {
                          aVoxel = aSubunit->voxel;
                        }
                    }
                  //We add the actual existing voxel if it is an immediate
                  //voxel, otherwise we add the shared voxel which is
                  //an immediate voxel:
                  aRefSubunit->targetVoxels[aBendIndex] = aVoxel;
                  addContPoint(aRefSubunit->targetVoxels[aBendIndex]->subunit,
                               aRefPoint);
                  //Successfully found an existing Voxel which is not blocked:
                  return 1;
                }
              //If another protomer is blocking the target position:
              else if(aDist < 0.8 && id2species(aVoxel->id)->getIsStatic())
                {
                  //There is existing Voxel but blocked by a protomer:
                  aRefSubunit->targetVoxels[aBendIndex] = NULL;
                  return 2;
                }
            }
        }

    }
  std::vector<Voxel*>& extendedSurface(aRefVoxel->surfaceVoxels[0][EXTENDED]);
  //Check the immediate surface voxels adjoining the immediate surface voxels
  //of the reference subunit, defined as the extended surface voxels:
  for(unsigned int i(0); i != extendedSurface.size(); ++i)
    { 
      Voxel* aVoxel(extendedSurface[i]);
      Subunit* aSubunit(aVoxel->subunit);
      if(!aSubunit->contPoints.empty())
        {
          for(unsigned int j(0); j !=aSubunit->contPoints.size(); ++j)
            {
              double aDist(getDistance(aRefPoint, &aSubunit->contPoints[j]));
              //If it is an existing target point:
              if(aDist < 0.1)
                {
                  //If it is a shared voxel:
                  if(aSubunit->voxel != aVoxel)
                    {
                      //Lets use the actual voxel if it is an extended voxel
                      //of the reference voxel:
                      unsigned int k(std::find(extendedSurface.begin(),
                        extendedSurface.end(), aSubunit->voxel)-
                                     extendedSurface.begin());
                      if(k != extendedSurface.size())
                        {
                          aVoxel = aSubunit->voxel;
                          i = k;
                        }
                    }
                  std::vector<Voxel*>& aSharedList(aRefVoxel->surfaceVoxels[0]
                                              [SHARED+i]);
                  //Find out if we can occupy at least one of the shared 
                  //voxels which is either already set to the reference
                  //molecule or vacant:
                  for(unsigned int l(0); l!=aSharedList.size(); ++l)
                    {
                      Voxel* aSharedVoxel(aSharedList[l]);
                      Subunit* aSharedSubunit(aSharedVoxel->subunit);
                      if(aSharedSubunit->voxel == aRefVoxel ||
                         id2species(aSharedVoxel->id)->getIsLipid())
                        {
                          aSharedVoxel->subunit->voxel = aRefVoxel;
                          aSharedVoxel->id = aRefVoxel->id;
                          addContPoint(aSharedVoxel->subunit,
                                       &aRefSubunit->subunitPoint);
                          aRefSubunit->sharedLipids[aBendIndex] = aSharedVoxel;
                          aRefSubunit->targetVoxels[aBendIndex] = aVoxel;
                          addContPoint(aRefSubunit->targetVoxels[
                                       aBendIndex]->subunit, aRefPoint); 
                          //Successfully found an existing Voxel which is not
                          //blocked:
                          return 1;
                        }
                    }
                  //Shared voxel has not been created previously by the
                  //reference voxel or it is blocked by a non-lipid molecule:
                  aRefSubunit->targetVoxels[aBendIndex] = NULL;
                  return 2;
                }
              else if(aDist < 0.8 && id2species(aVoxel->id)->getIsStatic())
                {
                  aRefSubunit->targetVoxels[aBendIndex] = NULL;
                  return 2;
                }
            }
        }
    }
  return 0;
}


unsigned int PolymerizationProcess::sortNearestTargetVoxels(
                                                  Subunit* aRefSubunit,
                                                  unsigned int aBendIndex,
                                                  std::vector<Voxel*>& targetVoxels, 
                                                  std::vector<double>& distances,
                                                  std::vector<int>& sharedIndices) 
{
  Voxel* aRefVoxel(aRefSubunit->voxel);
  Point* aRefPoint(&aRefSubunit->targetPoints[aBendIndex]);
  std::vector<Voxel*>& immediateSurface((*aRefVoxel->surfaceVoxels)[IMMEDIATE]);
  unsigned int bendSize(aRefSubunit->bendSize);
  unsigned int aSize(0);
  distances.resize(bendSize);
  targetVoxels.resize(bendSize);
  sharedIndices.resize(bendSize);
  //Get the top bendSize nearest Voxels from immediate Voxels:
  for(unsigned int i(0); i!=immediateSurface.size(); ++i)
    {
      Voxel* aVoxel(immediateSurface[i]);
      Subunit* aSubunit(aVoxel->subunit);
      if(aSubunit->contPoints.empty())
        {
          double aDist(getDistance(aRefPoint, &aSubunit->surfacePoint));
          //std::cout << "immediate:" << aDist << std::endl;
          unsigned int j(aSize);
          while(j)
            {
              if(aDist >= distances[j-1])
                {
                  break;
                }
              --j;
            }
          if(j < bendSize)
            {
              if(aSize < bendSize)
                {
                  ++aSize;
                }
              for(unsigned int k(aSize-1); k > j; --k)
                {
                  distances[k] = distances[k-1];
                  targetVoxels[k] = targetVoxels[k-1];
                  sharedIndices[k] = sharedIndices[k-1];
                }
              distances[j] = aDist;
              targetVoxels[j] = aVoxel;
              sharedIndices[j] = -1;
            }
        }
    }
  std::vector<Voxel*>& extendedSurface((*aRefVoxel->surfaceVoxels)[EXTENDED]);
  for(unsigned int i(0); i != extendedSurface.size(); ++i)
    { 
      Voxel* aVoxel(extendedSurface[i]);
      Subunit* aSubunit(aVoxel->subunit);
      if(aSubunit->contPoints.empty())
        {
          std::vector<Voxel*>& aSharedList((*aRefVoxel->surfaceVoxels)[SHARED+i]);
          for(unsigned int j(0); j!=aSharedList.size(); ++j)
            { 
              Voxel* aSharedVoxel(aSharedList[j]);
              if(aSharedVoxel->subunit->contPoints.empty() ||
                 aSharedVoxel->subunit->voxel == aRefVoxel)
                { 
                  double aDist(getDistance(aRefPoint, &aSubunit->surfacePoint));
                  //std::cout << "extended:" << aDist << std::endl;
                  unsigned int j(aSize);
                  while(j)
                    {
                      if(aDist >= distances[j-1])
                        {
                          break;
                        }
                      --j;
                    }
                  if(j < bendSize)
                    {
                      if(aSize < bendSize)
                        {
                          ++aSize;
                        }
                      for(unsigned int k(aSize-1); k > j; --k)
                        {
                          distances[k] = distances[k-1];
                          targetVoxels[k] = targetVoxels[k-1];
                          sharedIndices[k] = sharedIndices[k-1];
                        }
                      distances[j] = aDist;
                      targetVoxels[j] = aVoxel;
                      sharedIndices[j] = (int) i;
                    }
                  break;
                }
            }
        }
    }
  return aSize;
}



double PolymerizationProcess::setImmediateTargetVoxel(Subunit* aRefSubunit,
                                                      unsigned int aBendIndex) 
{
  Voxel* aRefVoxel(aRefSubunit->voxel);
  Point* aRefPoint(&aRefSubunit->targetPoints[aBendIndex]);
  double anImmediateDist(LARGE_DISTANCE);
  std::vector<Voxel*>& immediateSurface(aRefVoxel->surfaceVoxels[0][IMMEDIATE]);
  //Check the immediate 6 (usually) surface voxels adjoining the voxel of the
  //reference subunit:
  for(unsigned int i(0); i!=immediateSurface.size(); ++i)
    {
      Voxel* aVoxel(immediateSurface[i]);
      Subunit* aSubunit(aVoxel->subunit);
      //If the voxel does not already occupy a protomer
      //(A protomer has at least one contPoint -- its subunitPoint)
      //and it is not a shared voxel:
      if(aSubunit->contPoints.empty())
        {
          double aDist(getDistance(aRefPoint, &aSubunit->surfacePoint));
          if(aDist < anImmediateDist)
            {
              anImmediateDist = aDist;
              aRefSubunit->targetVoxels[aBendIndex] = aVoxel;
            }
        }
    }
  if(anImmediateDist < 0.7)
    {
      //We are definitely going to use the targetVoxel[aBendIndex] since
      //anImmediateDist is less than the cut off, so add the calculated
      //contPoint of the target protomer to its subunit:
      addContPoint(aRefSubunit->targetVoxels[aBendIndex]->subunit, aRefPoint);
      return 0;
    }
  //Note that at this point, aRefSubunit->targetVoxels[aBendIndex] is not set
  //to NULL, unless anImmediateDist == LARGE_DISTANCE:
  return anImmediateDist;
}

bool PolymerizationProcess::setExtendedTargetVoxel(Subunit* aRefSubunit,
                                                   unsigned int aBendIndex,
                                                   double extDist) 
{
  double immeDist(extDist);
  Voxel* aRefVoxel(aRefSubunit->voxel);
  Point* aRefPoint(&aRefSubunit->targetPoints[aBendIndex]);
  std::vector<Voxel*>& extendedSurface(aRefVoxel->surfaceVoxels[0][EXTENDED]);
  int extIndex(-1);
  //Check the immediate surface voxels adjoining the immediate surface voxels
  //of the reference subunit, defined as the extended surface voxels:
  for(unsigned int i(0); i != extendedSurface.size(); ++i)
    { 
      Voxel* aVoxel(extendedSurface[i]);
      Subunit* aSubunit(aVoxel->subunit);
      //If the voxel does not already occupy a protomer
      //(A protomer has at least one contPoint -- its subunitPoint)
      //and it is not a shared voxel:
      if(aSubunit->contPoints.empty())
        {
          double aDist(getDistance(aRefPoint, &aSubunit->surfacePoint));
          if(aDist < extDist)
            { 
              //Find the shared voxel which connects the reference voxel
              //to the extended voxel:
              std::vector<Voxel*>& 
                aSharedList(aRefVoxel->surfaceVoxels[0][SHARED+i]);
              //Check if there is an existing shared voxel or a voxel
              //that is unoccupied by a protomer, which connects
              //the reference voxel to the extended voxel:
              for(unsigned int j(0); j!=aSharedList.size(); ++j)
                { 
                  Voxel* aSharedVoxel(aSharedList[j]);
                  if(aSharedVoxel->subunit->contPoints.empty() ||
                     aSharedVoxel->subunit->voxel == aRefVoxel)
                    {
                      extDist = aDist;
                      extIndex = i;
                      break;
                    } 
                }
            }
        }
    }
  //If the distance is within cut off:
  //(Note that this distance could also be from an immediate voxel)
  if(extDist < 1.25)
    {
      //If we found an extended voxel, let us select the best shared voxel:
      if(extIndex != -1)
        {
          Voxel* aVoxel(extendedSurface[extIndex]);
          std::vector<Voxel*>& 
            aSharedList(aRefVoxel->surfaceVoxels[0][SHARED+extIndex]);
          double aSharedDist(LARGE_DISTANCE);
          Voxel* aSelectedSharedVoxel(NULL);
          for(unsigned int i(0); i!=aSharedList.size(); ++i)
            { 
              Voxel* aSharedVoxel(aSharedList[i]);
              Subunit* aSharedSubunit(aSharedVoxel->subunit);
              //First find an existing shared voxel which connects
              //the reference voxel to the extended voxel:
              if(aSharedSubunit->voxel == aRefVoxel)
                {
                  aSelectedSharedVoxel = aSharedVoxel;
                  break;
                }
              //Otherwise find a shared voxel that is unoccupied by a protomer:
              else if(aSharedSubunit->contPoints.empty())
                {
                  double aDist(getDistance(aRefPoint,
                                           &aSharedSubunit->surfacePoint));
                  if(aDist < aSharedDist)
                    {
                      aSharedDist = aDist;
                      aSelectedSharedVoxel = aSharedVoxel;
                    }
                }
            }
          //If the selected shared voxel is not an existing shared voxel nor
          //a lipid:
          if(!id2species(aSelectedSharedVoxel->id)->getIsLipid() &&
             aSelectedSharedVoxel->subunit->voxel != aRefVoxel)
            {
              //std::cout << "fail1 immed:" << immeDist << std::endl;
              return false;
            }
          aSelectedSharedVoxel->subunit->voxel = aRefVoxel;
          //Add the subunitPoint of the reference protomer to the 
          //list of contPoints of the selected shared subunit:
          addContPoint(aSelectedSharedVoxel->subunit,
                       &aRefSubunit->subunitPoint);
          aRefSubunit->sharedLipids[aBendIndex] = aSelectedSharedVoxel;
          aRefSubunit->targetVoxels[aBendIndex] = aVoxel;
        }
      addContPoint(aRefSubunit->targetVoxels[aBendIndex]->subunit, aRefPoint);
      return true;
    }
  return false;
}


/*


bool PolymerizationProcess::setTargetVoxels(Subunit* aSubunit)
{

  vector<int> anImmediateBend;
  vector<int> anExtendedBend;
  for(int i(0); i != aSubunit->targetBends.size(); ++i)
    {
      // j = 0; successfully found an existing Voxel which is not blocked
      // j = 1: there is existing Voxel but blocked
      // j = 1: there is another Voxel that overlaps the nextPoint
      // j = 2; no existing Voxel, but got an immediate Voxel with dist < 0.7
      //        targetVoxels[index] != NULL
      // j = 3; no existing Voxel, and no immediate Voxel with dist < 0.7 
      //        targetVoxels[index] == NULL
      int j(setExistingTargetVoxel(aSubunit, i));
      if(j == 1)
        {
          return false;
        }
      else if(j == 2)
        {
          anImmediateBend.push_back(i);
        }
      else if(j == 3)
        {
          anExtendedBend.push_back(i);
        }
    }
  for(vector<int>::iterator i(anImmediateBend.begin());
      i != anImmediateBend.end(); ++i)
    {
      if(aSubunit->tmpVoxels[*i]->subunit->contPoints.empty())
        {
          aSubunit->targetVoxels[*i] = aSubunit->tmpVoxels[*i];
          addContPoint(aSubunit->targetVoxels[*i], &aSubunit->targetPoints[*i]);
        }
      else
        {
          double aDist(setImmediateTargetVoxel(aSubunit, *i));
          if(aDist != 0)
            {
              if(!setExtendedTargetVoxel(aSubunit, *i, aDist))
                {
                  return false;
                }
            }
        }
    }
  for(vector<int>::iterator i(anExtendedBend.begin());
      i != anExtendedBend.end(); ++i)
    {
      if(!setExtendedTargetVoxel(aSubunit, *i, BIG_NUMBER))
        {
          return false;
        }
    }
  return true;
}
          double aDist(getDistance(aRefPoint, &aSubunit->surfacePoint));
          unsigned int j(aSize);
          while(j > 0)
            {
              if(aDist >= distances[--j] && aSize < bendSize)
                {
                  for(unsigned int k(aSize); k > j; --k)
                    {
                      distances[k] = distances[k-1];
                      targetVoxel[k] = targetVoxel[k-1];
                    }
                  break;
                }
            }
          if(j < bendSize)
            {
              sharedVoxels[j] = NULL;
              distances[j] = aDist;
              targetVoxel[j] = aVoxel;
              ++aSize;
            }






          while(j > 0)
            {
              if(aDist >= distances[--j] && aSize < bendSize)
                {
                  for(unsigned int k(aSize); k > j; --k)
                    {
                      distances[k] = distances[k-1];
                      targetVoxel[k] = targetVoxel[k-1];
                    }
                  break;
                }
            }
          if(j < bendSize)
            {
              sharedVoxels[j] = NULL;
              distances[j] = aDist;
              targetVoxel[j] = aVoxel;
              ++aSize;
            }


          unsigned int j(aSize);
          while(j > 0)
            {
              if(aDist >= distances[--j] && aSize-1 < bendSize)
                {
                  for(unsigned int k(aSize); k > j; --k)
                    {
                      distances[k] = distances[k-1];
                      targetVoxel[k] = targetVoxel[k-1];
                    }
                  distances[j] = aDist;
                  targetVoxel[j] = aVoxel;
                  break;
                }
            }
          if(j < bendSize)
            {
                  sharedVoxels[j] = NULL;

            }


*/
