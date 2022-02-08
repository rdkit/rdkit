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
#include "GraphMol/Trajectory/Trajectory.h"
#include "GraphMol/ROMol.h"

namespace python = boost::python;

namespace RDKit {

Snapshot *getSnapshot_wrap(Trajectory *traj, unsigned int snapshotNum) {
  return new Snapshot(traj->getSnapshot(snapshotNum));
}

Snapshot *constructSnapshot_wrap(python::list &coordList, double energy) {
  unsigned int l = python::len(coordList);
  boost::shared_array<double> c;
  if (l) {
    c.reset(new double[l]);
  }
  for (unsigned int i = 0; i < l; ++i) {
    c[i] = python::extract<double>(coordList[i]);
  }
  return new Snapshot(c, energy);
}

Snapshot *copyConstructSnapshot_wrap(Snapshot *other) {
  return new Snapshot(*other);
}

Trajectory *constructTrajectory_wrap(unsigned int dimension,
                                     unsigned int numPoints,
                                     python::list snapshotList) {
  unsigned int l = python::len(snapshotList);
  auto *traj = new Trajectory(dimension, numPoints);
  for (unsigned int i = 0; i < l; ++i) {
    Snapshot *s = python::extract<Snapshot *>(snapshotList[i]);
    traj->addSnapshot(*s);
  }
  return traj;
}

Trajectory *copyConstructTrajectory_wrap(Trajectory *other) {
  return new Trajectory(*other);
}

struct Trajectory_wrapper {
  static void wrap() {
    python::class_<Trajectory>(
        "Trajectory",
        "A class which allows storing Snapshots from a trajectory",
        python::no_init)
        .def("__init__",
             python::make_constructor(
                 &constructTrajectory_wrap, python::default_call_policies(),
                 (python::arg("dimension"), python::arg("numPoints"),
                  python::arg("snapshotList") = python::list())),
             "Constructor;\n"
             "dimension:    dimensionality of this Trajectory's coordinate "
             "tuples;\n"
             "numPoints:    number of coordinate tuples associated to each "
             "Snapshot;\n"
             "snapshotList: list of Snapshot objects used to initialize the "
             "Trajectory (optional; defaults to []).\n")
        .def("__init__",
             python::make_constructor(&copyConstructTrajectory_wrap,
                                      python::default_call_policies(),
                                      python::arg("other")),
             "Copy constructor")
        .def("Dimension", &Trajectory::dimension, (python::arg("self")),
             "returns the dimensionality of this Trajectory's coordinate "
             "tuples")
        .def("NumPoints", &Trajectory::numPoints, (python::arg("self")),
             "returns the number of coordinate tuples associated to each "
             "Snapshot")
        .def("__len__", &Trajectory::size)
        .def("AddSnapshot", &Trajectory::addSnapshot,
             (python::arg("self"), python::arg("s")),
             "appends Snapshot s to this Trajectory; returns the zero-based "
             "index position of the added snapshot\n")
        .def("GetSnapshot", getSnapshot_wrap,
             (python::arg("self"), python::arg("snapshotNum")),
             "returns the Snapshot snapshotNum, where the latter is the "
             "zero-based index of the retrieved Snapshot\n",
             python::return_value_policy<python::manage_new_object>())
        .def(
            "InsertSnapshot", &Trajectory::insertSnapshot,
            (python::arg("self"), python::arg("snapshotNum"), python::arg("s")),
            "inserts Snapshot s into the Trajectory at the position "
            "snapshotNum, where the latter is the zero-based index of the "
            "Trajectory's Snapshot before which the Snapshot s will be "
            "inserted; returns the zero-based index position of the inserted "
            "snapshot\n")
        .def("RemoveSnapshot", &Trajectory::removeSnapshot,
             (python::arg("self"), python::arg("snapshotNum")),
             "removes Snapshot snapshotNum from the Trajectory, where "
             "snapshotNum is the zero-based index of Snapshot to be removed\n")
        .def("Clear", &Trajectory::clear, (python::arg("self")),
             "removes all Snapshots from the Trajectory\n")
        .def("AddConformersToMol", &Trajectory::addConformersToMol,
             (python::arg("self"), python::arg("mol"), python::arg("from") = -1,
              python::arg("to") = -1),
             "adds conformations from the Trajectory to mol\n"
             "from is the first Snapshot that will be added as a Conformer; "
             "defaults to -1 (first available)\n"
             "to is the last Snapshot that will be added as a Conformer; "
             "defaults to -1 (all)\n");

    python::class_<Snapshot>(
        "Snapshot",
        "A class which allows storing coordinates from a trajectory",
        python::no_init)
        .def("__init__",
             python::make_constructor(
                 &constructSnapshot_wrap, python::default_call_policies(),
                 (python::arg("coordList"), python::arg("energy") = 0.0)),
             "Constructor;\n"
             "coordList: list of floats containing the coordinates for this "
             "Snapshot;\n"
             "energy:    the energy for this Snapshot.\n")
        .def("__init__",
             python::make_constructor(&copyConstructSnapshot_wrap,
                                      python::default_call_policies(),
                                      python::arg("other")),
             "Copy constructor")
        .def("GetPoint2D", &Snapshot::getPoint2D,
             (python::arg("self"), python::arg("pointNum")),
             "return the coordinates at pointNum as a Point2D object; "
             "requires the Trajectory dimension to be == 2")
        .def("GetPoint3D", &Snapshot::getPoint3D,
             (python::arg("self"), python::arg("pointNum")),
             "return the coordinates at pointNum as a Point3D object; "
             "requires the Trajectory dimension to be >= 2")
        .def("GetEnergy", &Snapshot::getEnergy, (python::arg("self")),
             "returns the energy for this Snapshot")
        .def("SetEnergy", &Snapshot::setEnergy,
             (python::arg("self"), python::arg("energy")),
             "sets the energy for this Snapshot");
    python::def("ReadAmberTrajectory", &readAmberTrajectory,
                (python::arg("fName"), python::arg("traj")),
                "reads coordinates from an AMBER trajectory file into the "
                "Trajectory object; "
                "returns the number of Snapshot objects read in\n");
    python::def("ReadGromosTrajectory", &readGromosTrajectory,
                (python::arg("fName"), python::arg("traj")),
                "reads coordinates from a GROMOS trajectory file into the "
                "Trajectory object; "
                "returns the number of Snapshot objects read in\n");
  };
};
}  // namespace RDKit

void wrap_trajectory() { RDKit::Trajectory_wrapper::wrap(); }
