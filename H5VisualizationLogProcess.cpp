//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//
//                This file is part of E-Cell Simulation Environment package
//
//                                Copyright (C) 2006-2009 Keio University
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

#include <H5Cpp.h>
#include "H5Support.hpp"
#include <MethodProxy.hpp>
#include <boost/scoped_array.hpp>
#include <boost/mpl/and.hpp>
#include <boost/utility/enable_if.hpp>
#include <boost/type_traits/is_integral.hpp>
#include <boost/type_traits/is_signed.hpp>
#include <boost/type_traits/is_unsigned.hpp>
#include <boost/type_traits/make_unsigned.hpp>
#include "SpatiocyteProcess.hpp"
#include "SpatiocyteSpecies.hpp"

class ParticleData
{
public:
    std::size_t getID() const
    {
        return theID;
    }

    unsigned short getSpeciesID() const
    {
        return theSpeciesID;
    }

    unsigned short getVoxelID() const
    {
        return theVoxelID;
    }

    ParticleData(std::size_t id, unsigned short speciesID, unsigned short voxelID)
        : theID(id), theSpeciesID(speciesID), theVoxelID(voxelID) {}

private:
    const unsigned int theID;
    const unsigned short theSpeciesID;
    const unsigned short theVoxelID;
};


template<typename Tarc>
struct PointDataPacker
{
    typedef Tarc archiver_type;

    void operator()(archiver_type& arc, Point const* data = 0) const
    {
        arc << field("x", &Point::x, data);
        arc << field("y", &Point::y, data);
        arc << field("z", &Point::z, data);
    }
};

template<typename Tarc>
struct ParticleDataPacker
{
    typedef Tarc archiver_type;

    void operator()(archiver_type& arc, ParticleData const* data = 0) const
    {
        arc << field<uint64_t>("id", &ParticleData::getID, data);
        arc << field<uint64_t>("species_id", &ParticleData::getSpeciesID, data);
        arc << field<uint64_t>("voxel_id", &ParticleData::getVoxelID, data);
    }
};

template<typename Tarc>
struct SpeciesDataPacker
{
    typedef Tarc archiver_type;

    static const char* Species_getName(Species const& self)
    {
        return self.getVariable()->getName().c_str();
    }

    void operator()(archiver_type& arc, Species const* data = 0) const
    {
        arc << field<uint64_t>("id", &Species::getID, data);
        arc << field<char[32]>("name", &Species_getName, data);
        arc << field<Species, double, double>("radius", 0.);
        arc << field<Species, double, double>("D", 0.);
    }
};

template<template<typename> class TTserialize_>
H5::CompType getH5Type()
{
    return get_h5_type<h5_le_traits, TTserialize_>();
}

template<typename T>
unsigned char* pack(unsigned char* buffer, T const& container)
{
    return pack<h5_le_traits>(buffer, container);
}

template<template<typename> class TTserialize_, typename T>
unsigned char* pack(unsigned char* buffer, T const& container)
{
    return pack<h5_le_traits, TTserialize_>(buffer, container);
}

LIBECS_DM_CLASS(H5VisualizationLogProcess, SpatiocyteProcess)
{ 
public:
    typedef h5_le_traits traits_type;
    typedef traits_type::packer_type packer_type;

public:
    LIBECS_DM_OBJECT(H5VisualizationLogProcess, Process)
    {
        INHERIT_PROPERTIES(Process);
        PROPERTYSLOT_SET_GET(Integer, Polymer);
        PROPERTYSLOT_SET_GET(Real, LogInterval);
        PROPERTYSLOT_SET_GET(String, FileName);
    }

    H5VisualizationLogProcess();

    virtual ~H5VisualizationLogProcess() {}
    SIMPLE_SET_GET_METHOD(Integer, Polymer);
    SIMPLE_SET_GET_METHOD(Real, LogInterval);
    SIMPLE_SET_GET_METHOD(String, FileName);

    virtual void initializeFourth()
    {
        if(LogInterval > 0)
        {
            theStepInterval = LogInterval;
        }
        else
        {
            //Use the smallest time step of all queued events for
            //the step interval:
            theTime = libecs::INF;
            thePriorityQueue->move(theQueueID);
            theStepInterval = thePriorityQueue->getTop()->getTime();
        }
        theTime = theStepInterval;
        thePriorityQueue->move(theQueueID);
    }

    virtual void initializeLastOnce()
    {
        theLogFile = H5::H5File(FileName, H5F_ACC_TRUNC);
        theDataGroup = theLogFile.createGroup("data");
        initializeLog();
        logSpecies();
    }

    virtual void fire()
    {
        logSpecies();
        if(LogInterval > 0)
        {
            theTime += LogInterval;
            thePriorityQueue->moveTop();
        }
        else
        {
            //get the next step interval of the SpatiocyteStepper:
            double aTime(theTime);
            theTime = libecs::INF;
            thePriorityQueue->moveTop();
            if(thePriorityQueue->getTop()->getTime() > aTime)
            {
                theStepInterval = thePriorityQueue->getTop()->getTime() -
                    theSpatiocyteStepper->getCurrentTime();
            }
            theTime = aTime + theStepInterval;
            thePriorityQueue->move(theQueueID);
        }
    }
protected:
    virtual void initializeLog();
    void logSpecies();
    void logMolecules(H5::DataSpace const& space, H5::DataSet const& dataSet, hsize_t (&dims)[1], Species *);

    template<typename T>
    void setH5Attribute(const char* name, T const& data);

    void setH5Attribute(const char* name, Point const& data);

protected:
    unsigned int Polymer;
    unsigned int theLogMarker;
    double LogInterval;
    String FileName;
    packer_type packer;
    H5::H5File theLogFile;
    H5::Group theDataGroup;
    H5::CompType pointDataType;
    H5::CompType particleDataType;
    H5::CompType speciesDataType;
};

H5VisualizationLogProcess::H5VisualizationLogProcess()
    :   SpatiocyteProcess(false),
        Polymer(1),
        theLogMarker(UINT_MAX),
        LogInterval(0),
        FileName("visualLog.h5"),
        pointDataType(getH5Type<PointDataPacker>()),
        particleDataType(getH5Type<ParticleDataPacker>()),
        speciesDataType(getH5Type<SpeciesDataPacker>())
{
}

template<typename T>
void H5VisualizationLogProcess::setH5Attribute(const char* name, T const& data)
{
    H5::DataType dataType = get_h5_scalar_data_type_le<T>()();
    unsigned char buf[sizeof(T)];
    packer(buf, data);
    theDataGroup.createAttribute(name, dataType, H5::DataSpace()).write(dataType, buf); 
}

void H5VisualizationLogProcess::setH5Attribute(const char* name, Point const& data)
{
    H5::Attribute attr(theDataGroup.createAttribute(name, pointDataType, H5::DataSpace()));
    unsigned char buf[24];
    BOOST_ASSERT(sizeof(buf) >= pointDataType.getSize());
    pack<PointDataPacker>(buf, data);
    attr.write(pointDataType, buf); 
}

void H5VisualizationLogProcess::initializeLog()
{
    setH5Attribute("start_coord", theSpatiocyteStepper->getStartCoord());
    setH5Attribute("row_size", theSpatiocyteStepper->getRowSize());
    setH5Attribute("layer_size", theSpatiocyteStepper->getLayerSize());
    setH5Attribute("column_size", theSpatiocyteStepper->getColSize());
    Point aCenterPoint(theSpatiocyteStepper->getCenterPoint());
    setH5Attribute("center_point", aCenterPoint);
    setH5Attribute("real_row_size", aCenterPoint.z * 2);
    setH5Attribute("real_layer_size", aCenterPoint.y * 2);
    setH5Attribute("real_col_size", aCenterPoint.x * 2);
    setH5Attribute("normalized_voxel_radius",  theSpatiocyteStepper->getNormalizedVoxelRadius());
}

void H5VisualizationLogProcess::logMolecules(H5::DataSpace const& space, H5::DataSet const& dataSet, hsize_t (&dims)[1], Species* aSpecies)
{
    //No need to log lipid or vacant molecules since the size is 0:
    if (aSpecies->getIsVacant())
    {
        return;
    }
    //theLogFile.write((char*)(&anIndex), sizeof(anIndex));
    //The species molecule size:
    const int aSize(aSpecies->size());

    if (!aSize)
        return;

    hsize_t wdim[] = { aSize };
    const hsize_t offset[] = { dims[0] };

    dims[0] += wdim[0];
    dataSet.extend(dims);

    H5::DataSpace mem(1, wdim);
    H5::DataSpace slab(dataSet.getSpace());
    slab.selectHyperslab(H5S_SELECT_SET, wdim, offset);
    boost::scoped_array<unsigned char> buf(new unsigned char[particleDataType.getSize() * aSize]);
    unsigned char* p = buf.get();
    for(int i(0); i != aSize; ++i)
    {
        p = pack<ParticleDataPacker>(p, ParticleData(i, aSpecies->getID(), aSpecies->getMolecule(i)->id));
    }
    dataSet.write(buf.get(), particleDataType, mem, slab);
}

void H5VisualizationLogProcess::logSpecies()
{
    const Time currentTime(theSpatiocyteStepper->getCurrentTime());

    H5::Group perTimeDataGroup(theDataGroup.createGroup(boost::lexical_cast<std::string>(currentTime).c_str()));
    H5::DataSpace space;
    H5::DataSet dataSet;
    {
        static const hsize_t initdims[] = { 0 };
        static const hsize_t maxdims[] = { H5S_UNLIMITED };
        static const hsize_t chunkdims[] = { 4096 };
        space = H5::DataSpace(1, initdims, maxdims);
        H5::DSetCreatPropList props;
        props.setChunk(1, chunkdims);
        dataSet = perTimeDataGroup.createDataSet("particles", particleDataType, space, props);
    }

    hsize_t dims[] = { 0 };

    for(unsigned int i(0); i != theProcessSpecies.size(); ++i)
    {
        logMolecules(space, dataSet, dims, theProcessSpecies[i]);
    }
}

LIBECS_DM_INIT(H5VisualizationLogProcess, Process); 
