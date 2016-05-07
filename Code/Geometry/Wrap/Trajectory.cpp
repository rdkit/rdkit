// $Id$
//
//  Copyright (C) 2016 Sereina Riniker, Paolo Tosco
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

Snapshot *getSnapshot_wrap(Trajectory *traj, unsigned int snapshotNum) {
  return new Snapshot(traj->getSnapshot(snapshotNum));
}

Snapshot *constructSnapshot_wrap(python::list &coordList, double energy) {
  unsigned int l = python::len(coordList);
  boost::shared_array<double> c;
  if (l)
    c.reset(new double[l]);
  for (unsigned int i = 0; i < l; ++i)
    c[i] = python::extract<double>(coordList[i]);
  return new Snapshot(c, energy);
}

Snapshot *copyConstructSnapshot_wrap(Snapshot *other) {
  return new Snapshot(*other);
}

struct Trajectory_wrapper {
  static void wrap() {
    std::string docString =
        "A class which allows storing Snapshots from a trajectory.\n\n\
        Usage example:\n\
        traj = Trajectory(dimension, numPoints)\n\
       \n\
       ARGUMENTS\n\
        - dimension   dimensionality of this Trajectory's coordinate tuples\n\
        - numPoints   number of coordinate tuples associated to each Snapshot\n\
         \n\
        RETURNS\n\
        the Trajectory object\n\
      \n";
    python::class_<Trajectory>(
        "Trajectory", docString.c_str(), python::init<unsigned int, unsigned int>(
        (python::arg("dimension"), python::arg("numPoints"))))
        .def(python::init<const Trajectory&>(python::arg("other")))
        .def("Dimension", &Trajectory::dimension, (python::arg("self")),
             "return the dimensionality of this Trajectory's coordinate tuples")
        .def("NumPoints", &Trajectory::numPoints, (python::arg("self")),
             "returns the number of coordinate tuples associated to each Snapshot")
        .def("__len__", &Trajectory::size)
        .def("AddSnapshot", &Trajectory::addSnapshot,
             (python::arg("self"), python::arg("s")),
             "appends Snapshot s to this Trajectory; returns the zero-based index "
             "position of the added snapshot\n")
        .def("GetSnapshot", getSnapshot_wrap,
             (python::arg("self"), python::arg("snapshotNum")),
             "returns the Snapshot snapshotNum, where the latter is the zero-based "
             "index of the retrieved Snapshot\n",
             python::return_value_policy<python::manage_new_object>())
        .def("InsertSnapshot", &Trajectory::insertSnapshot,
             (python::arg("self"), python::arg("snapshotNum"), python::arg("s")),
             "inserts Snapshot s into the Trajectory at the position snapshotNum, "
             "where the latter is the zero-based index of the Trajectory's Snapshot "
             "before which the Snapshot s will be inserted; returns the zero-based "
             "index position of the inserted snapshot\n")
        .def("RemoveSnapshot", &Trajectory::removeSnapshot,
             (python::arg("self"), python::arg("snapshotNum")),
             "removes Snapshot snapshotNum from the Trajectory, where "
             "snapshotNum is the zero-based index of Snapshot to be removed\n");

    docString =
        "A class which allows storing coordinates from a trajectory.\n\n\
        Usage example:\n\
        snapshot = Snapshot(pos, energy)\n\n\
       \n\
       ARGUMENTS\n\
        - pos      list of floats containing the coordinates for this Snapshot\n\
        - energy   energy for this Snapshot\n\
         \n\
        RETURNS\n\
        the Snapshot object\n\
      \n";
    python::class_<Snapshot>(
        "Snapshot", docString.c_str(), python::no_init)
        .def("GetPoint2D", &Snapshot::getPoint2D, (python::arg("self"), python::arg("pointNum")),
             "return the coordinates at pointNum as a Point2D object; "
             "requires the Trajectory dimension to be == 2")
        .def("GetPoint3D", &Snapshot::getPoint3D, (python::arg("self"), python::arg("pointNum")),
             "return the coordinates at pointNum as a Point3D object; "
             "requires the Trajectory dimension to be >= 2")
        .def("GetEnergy", &Snapshot::getEnergy, (python::arg("self")),
             "returns the energy for this Snapshot")
        .def("SetEnergy", &Snapshot::setEnergy, (python::arg("self"), python::arg("energy")),
             "sets the energy for this Snapshot");
    python::def("Snapshot", constructSnapshot_wrap,
        (python::arg("coordList"), python::arg("energy") = 0.0),
        "Constructor", python::return_value_policy<python::manage_new_object>());
    python::def("Snapshot", copyConstructSnapshot_wrap, (python::arg("other")),
        "Copy constructor", python::return_value_policy<python::manage_new_object>());
    python::def("ReadAmberTrajectory", &readAmberTrajectory,
        (python::arg("fName"), python::arg("traj")),
        "reads coordinates from an AMBER trajectory file into the Trajectory object; "
        "returns the number of Snapshot objects read in\n");
    python::def("ReadGromosTrajectory", &readGromosTrajectory,
        (python::arg("fName"), python::arg("traj")),
        "reads coordinates from a GROMOS trajectory file into the Trajectory object; "
        "returns the number of Snapshot objects read in\n");
  };
  
};
}

void wrap_trajectory() { RDGeom::Trajectory_wrapper::wrap(); }
