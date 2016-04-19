// $Id$
//
//  Copyright (C) 2016 Paolo Tosco
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <RDBoost/python.h>

#include <RDBoost/Wrap.h>
#include <Geometry/Trajectory.h>

namespace python = boost::python;

namespace RDGeom {

class PySnapshot;
class PyTrajectory;
PySnapshot *getSnapshot_wrap(const PyTrajectory *pyTrajectory, unsigned int snapshotNum);

class PySnapshot {
 public:
  PySnapshot(Snapshot *s) : snapshot(s) {};
  PySnapshot(boost::shared_ptr<Snapshot> s) : snapshot(s) {};
  ~PySnapshot() {};
  Point2D getPoint2D(unsigned int pointNum) const {
    Point2D res = snapshot.get()->getPoint2D(pointNum);
    return res;
  }
  Point3D getPoint3D(unsigned int pointNum) const {
    Point3D res = snapshot.get()->getPoint3D(pointNum);
    return res;
  }
  double getEnergy() const {
    return snapshot.get()->getEnergy();
  }
  void setEnergy(double energy) {
    snapshot.get()->setEnergy(energy);
  }
  void freePos() {
    snapshot.get()->freePos();
  }
  boost::shared_ptr<Snapshot> snapshot;
};

class PyTrajectory {
 public:
  PyTrajectory(Trajectory *t) : traj(t) {};
  PyTrajectory(boost::shared_ptr<Trajectory> t) : traj(t) {};
  ~PyTrajectory() {};
  unsigned int dimension() const {
    return traj.get()->dimension();
  }
  unsigned int numPoints() const {
    return traj.get()->numPoints();
  }
  size_t size() const {
    return traj.get()->size();
  }
  unsigned int addSnapshot(Snapshot s) {
    return traj.get()->addSnapshot(s);
  }
  PySnapshot *getSnapshot(unsigned int snapshotNum) const {
    return getSnapshot_wrap(this, snapshotNum);
  }
  unsigned int insertSnapshot(unsigned int snapshotNum, Snapshot s) {
    return traj.get()->insertSnapshot(snapshotNum, s);
  }
  unsigned int removeSnapshot(unsigned int snapshotNum) {
    return traj.get()->removeSnapshot(snapshotNum);
  }
  unsigned int readAmber(const std::string &fName) {
    return traj.get()->readAmber(fName);
  }
  unsigned int readGromos(const std::string &fName) {
    return traj.get()->readGromos(fName);
  }
  boost::shared_ptr<Trajectory> traj;
};

PyTrajectory *constructTrajectory(unsigned int dimension, unsigned int numPoints) {
  Trajectory *traj;
  {
    NOGIL gil;
    traj = new Trajectory(dimension, numPoints);
  }
  PyTrajectory *pyTrajectory = new PyTrajectory(traj);

  return pyTrajectory;
}

PySnapshot *constructSnapshot(python::list &coordList, double energy) {
  Snapshot *snapshot;
  unsigned int l = python::len(coordList);
  double *coords = new double[l];
  for (unsigned int i = 0; i < l; ++i)
    coords[i] = python::extract<double>(coordList[i]);
  {
    NOGIL gil;
    snapshot = new Snapshot(coords, energy, true);
  }
  PySnapshot *pySnapshot = new PySnapshot(snapshot);
  
  return pySnapshot;
}

PySnapshot *getSnapshot_wrap(const PyTrajectory *pyTrajectory, unsigned int snapshotNum) {
  Snapshot *snapshot;
  {
    NOGIL gil;
    snapshot = new Snapshot(pyTrajectory->traj.get()->getSnapshot(snapshotNum));
  }
  PySnapshot *pySnapshot = new PySnapshot(snapshot);
  
  return pySnapshot;
}

std::string TrajectoryClassDoc =
    "A class which allows storing coordinates from a trajectory.\n";

struct Trajectory_wrapper {
  static void wrap() {
    python::class_<PyTrajectory, boost::shared_ptr<PyTrajectory> >(
        "Trajectory", "Trajectory object", python::no_init)
        .def("Dimension", &PyTrajectory::dimension, (python::arg("self")),
             "return the dimensionality of this Trajectory's coordinate tuples")
        .def("NumPoints", &PyTrajectory::numPoints, (python::arg("self")),
             "returns the number of coordinate tuples associated to each Snapshot")
        .def("__len__", &PyTrajectory::size)
        .def("AddSnapshot", &PyTrajectory::addSnapshot,
             (python::arg("self"), python::arg("s")),
             "appends Snapshot s to this Trajectory; returns the zero-based index "
             "position of the added snapshot\n")
        .def("GetSnapshot", &PyTrajectory::getSnapshot,
             (python::arg("self"), python::arg("snapshotNum")),
             "returns the Snapshot snapshotNum, where the latter is the zero-based "
             "index of the retrieved Snapshot\n",
             python::return_value_policy<python::manage_new_object>())
        .def("InsertSnapshot", &PyTrajectory::insertSnapshot,
             (python::arg("self"), python::arg("snapshotNum"), python::arg("s")),
             "inserts Snapshot s into the Trajectory at the position snapshotNum, "
             "where the latter is the zero-based index of the Trajectory's Snapshot "
             "before which the Snapshot s will be inserted; returns the zero-based "
             "index position of the inserted snapshot\n")
        .def("RemoveSnapshot", &PyTrajectory::removeSnapshot,
             (python::arg("self"), python::arg("snapshotNum")),
             "removes Snapshot snapshotNum from the Trajectory, where "
             "snapshotNum is the zero-based index of Snapshot to be removed\n")
        .def("ReadAmber", &PyTrajectory::readAmber,
             (python::arg("self"), python::arg("fName")),
             "reads coordinates from an AMBER trajectory file into the Trajectory object; "
             "returns the number of Snapshot objects read in\n")
        .def("ReadGromos", &PyTrajectory::readGromos,
             (python::arg("self"), python::arg("fName")),
             "reads coordinates from a GROMOS trajectory file into the Trajectory object; "
             "returns the number of Snapshot objects read in\n");
    std::string docString =
        "Get a Trajectory object\n\
       \n\
       ARGUMENTS\n\
        - dimension   dimensionality of this Trajectory's coordinate tuples\n\
        - numPoints   number of coordinate tuples associated to each Snapshot\n\
         \n\
        RETURNS\n\
        the Trajectory object\n\
      \n";
    python::def("Trajectory", constructTrajectory,
                (python::arg("dimension"), python::arg("numPoints")),
                 python::return_value_policy<python::manage_new_object>(),
                 docString.c_str());

    python::class_<PySnapshot, boost::shared_ptr<PySnapshot> >(
        "Snapshot", "Snapshot object", python::no_init)
        .def("GetPoint2D", &PySnapshot::getPoint2D, (python::arg("self"), python::arg("pointNum")),
             "return the coordinates at pointNum as a Point2D object; "
             "requires the Trajectory dimension to be == 2")
        .def("GetPoint3D", &PySnapshot::getPoint3D, (python::arg("self"), python::arg("pointNum")),
             "return the coordinates at pointNum as a Point3D object; "
             "requires the Trajectory dimension to be >= 2")
        .def("GetEnergy", &PySnapshot::getEnergy, (python::arg("self")),
             "returns the energy for this Snapshot")
        .def("SetEnergy", &PySnapshot::setEnergy, (python::arg("self"), python::arg("energy")),
             "sets the energy for this Snapshot")
        .def("FreePos", &PySnapshot::freePos, (python::arg("self")),
             "Frees the pointer to the C++ array of doubles where the "
             "coordinates for this snapshot are stored\n");
    docString =
        "Get a Snapshot object\n\
       \n\
       ARGUMENTS\n\
        - pos      array of floats containing the coordinates for this Snapshot\n\
        - energy   energy for this Snapshot\n\
         \n\
        RETURNS\n\
        the Snapshot object\n\
      \n";
    python::def("Snapshot", constructSnapshot,
                (python::arg("pos") = python::list(), python::arg("energy") = 0.0),
                 python::return_value_policy<python::manage_new_object>(),
                 docString.c_str());
  };
  
};
}

void wrap_trajectory() { RDGeom::Trajectory_wrapper::wrap(); }
