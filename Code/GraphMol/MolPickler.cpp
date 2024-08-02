//
//  Copyright (C) 2001-2021 Greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <GraphMol/RDKitBase.h>
#include <GraphMol/RDKitQueries.h>
#include <GraphMol/MolPickler.h>
#include <GraphMol/QueryOps.h>
#include <GraphMol/MonomerInfo.h>
#include <GraphMol/StereoGroup.h>
#include <GraphMol/SubstanceGroup.h>
#include <RDGeneral/utils.h>
#include <RDGeneral/RDLog.h>
#include <RDGeneral/StreamOps.h>
#include <RDGeneral/types.h>
#include <DataStructs/DatastructsStreamOps.h>
#include <Query/QueryObjects.h>
#include <map>
#include <iostream>
#include <cstdint>
#include <boost/algorithm/string.hpp>

#ifdef RDK_BUILD_THREADSAFE_SSS
#include <mutex>
#endif

using std::int32_t;
using std::uint32_t;

namespace RDKit {

const int32_t MolPickler::versionMajor = 16;
const int32_t MolPickler::versionMinor = 2;
const int32_t MolPickler::versionPatch = 0;
const int32_t MolPickler::endianId = 0xDEADBEEF;

void streamWrite(std::ostream &ss, MolPickler::Tags tag) {
  auto tmp = static_cast<unsigned char>(tag);
  streamWrite(ss, tmp);
}
template <typename T>
void streamWrite(std::ostream &ss, MolPickler::Tags tag, const T &what) {
  streamWrite(ss, tag);
  streamWrite(ss, what);
};

void streamRead(std::istream &ss, MolPickler::Tags &tag, int version) {
  if (version < 7000) {
    int32_t tmp;
    streamRead(ss, tmp, version);
    if (tmp < 0 || tmp >= MolPickler::Tags::INVALID_TAG) {
      throw MolPicklerException("Invalid tag found.");
    }
    tag = static_cast<MolPickler::Tags>(tmp);
  } else {
    unsigned char tmp;
    streamRead(ss, tmp, version);
    if (tmp >= MolPickler::Tags::INVALID_TAG) {
      throw MolPicklerException("Invalid tag found.");
    }
    tag = static_cast<MolPickler::Tags>(tmp);
  }
}

void streamReadPositiveChar(std::istream &ss, char &res, int version) {
  streamRead(ss, res, version);
  if (res < 0) {
    throw MolPicklerException("invalid value in pickle");
  }
}

namespace {
static unsigned int defaultProperties = PicklerOps::NoProps;
static CustomPropHandlerVec defaultPropHandlers = {};

#ifdef RDK_BUILD_THREADSAFE_SSS
std::mutex &propmutex_get() {
  // create on demand
  static std::mutex _mutex;
  return _mutex;
}

void propmutex_create() {
  std::mutex &mutex = propmutex_get();
  std::lock_guard<std::mutex> test_lock(mutex);
}

std::mutex &GetPropMutex() {
  static std::once_flag flag;
  std::call_once(flag, propmutex_create);
  return propmutex_get();
}
#endif

void write_sstream_to_stream(std::ostream &outStream,
                             const std::stringstream &toWrite) {
  auto ts = toWrite.str();
  int32_t tmpInt = static_cast<int32_t>(ts.size());
  streamWrite(outStream, tmpInt);
  outStream.write(ts.c_str(), sizeof(char) * tmpInt);
}

}  // namespace

void MolPickler::_pickleProperties(std::ostream &ss, const RDProps &props,
                                   unsigned int pickleFlags) {
  if (!pickleFlags) {
    return;
  }

  streamWriteProps<std::uint16_t>(ss, props,
                                  pickleFlags & PicklerOps::PrivateProps,
                                  pickleFlags & PicklerOps::ComputedProps,
                                  MolPickler::getCustomPropHandlers());
}

namespace {

template <typename SAVEAS, typename STOREAS>
void unpickleExplicitProperties(
    std::istream &ss, RDProps &props, int version,
    const std::vector<std::pair<std::string, std::uint16_t>> &explicitProps) {
  if (version >= 14000) {
    std::uint8_t bprops;
    streamRead(ss, bprops, version);
    for (const auto &pr : explicitProps) {
      if (bprops & pr.second) {
        SAVEAS bv;
        streamRead(ss, bv, version);
        props.setProp(pr.first, static_cast<STOREAS>(bv));
      }
    }
  }
}

template <typename SAVEAS>
bool pickleExplicitProperties(
    std::ostream &ss, const RDProps &props,
    const std::vector<std::pair<std::string, std::uint16_t>> &explicitProps) {
  std::uint8_t bprops = 0;
  std::vector<SAVEAS> ps;
  SAVEAS bv;
  for (const auto &pr : explicitProps) {
    if (props.getPropIfPresent(pr.first, bv)) {
      bprops |= pr.second;
      ps.push_back(bv);
    }
  }

  streamWrite(ss, bprops);
  for (auto v : ps) {
    streamWrite(ss, v);
  }
  return !ps.empty();
}

// you're going to scratch your head and wonder why on earth this is using a
// class instead of a set of globals. This is a workaround for some very strange
// problem with globals and the cartridge
class PropTracker {
 public:
  // this is stored as bitflags in a byte, so don't exceed 8 entries or we need
  // to update the pickle format.
  // the properties themselves are stored as std::int8_t
  const std::vector<std::pair<std::string, std::uint16_t>> explicitBondProps = {
      {RDKit::common_properties::_MolFileBondType, 0x1},
      {RDKit::common_properties::_MolFileBondStereo, 0x2},
      {RDKit::common_properties::_MolFileBondCfg, 0x4},
      {RDKit::common_properties::_MolFileBondQuery, 0x8},
      {RDKit::common_properties::molStereoCare, 0x10},
  };
  // this is stored as bitflags in a byte, so don't exceed 8 entries or we need
  // to update the pickle format.
  // the properties themselves are stored as std::int16_t
  const std::vector<std::pair<std::string, std::uint16_t>> explicitAtomProps = {
      {common_properties::molStereoCare, 0x1},
      {common_properties::molParity, 0x2},
      {common_properties::molInversionFlag, 0x4},
      {common_properties::_ChiralityPossible, 0x8},

  };
  const std::vector<std::string> ignoreAtomProps = {
      common_properties::molAtomMapNumber,
      common_properties::dummyLabel,
  };
  std::unordered_set<std::string> ignoreBondProps;
  PropTracker() {
    for (const auto &pr : explicitBondProps) {
      ignoreBondProps.insert(pr.first);
    }
  };
};

bool pickleAtomProperties(std::ostream &ss, const RDProps &props,
                          unsigned int pickleFlags) {
  const static PropTracker aprops;
  static std::unordered_set<std::string> ignoreProps;
  if (ignoreProps.empty()) {
    for (const auto &pr : aprops.explicitAtomProps) {
      ignoreProps.insert(pr.first);
    }
    for (const auto &pn : aprops.ignoreAtomProps) {
      ignoreProps.insert(pn);
    }
  }

  if (!pickleFlags) {
    return false;
  }

  bool res = streamWriteProps<std::uint16_t>(
      ss, props, pickleFlags & PicklerOps::PrivateProps,
      pickleFlags & PicklerOps::ComputedProps,
      MolPickler::getCustomPropHandlers(), ignoreProps);

  res |= pickleExplicitProperties<std::int16_t>(ss, props,
                                                aprops.explicitAtomProps);
  return res;
}

void unpickleAtomProperties(std::istream &ss, RDProps &props, int version) {
  const static PropTracker aprops;
  if (version >= 14000) {
    streamReadProps<std::uint16_t>(ss, props,
                                   MolPickler::getCustomPropHandlers(), false);
  } else {
    streamReadProps<unsigned int>(ss, props,
                                  MolPickler::getCustomPropHandlers(), false);
  }
  unpickleExplicitProperties<std::int16_t, int>(ss, props, version,
                                                aprops.explicitAtomProps);
}

bool pickleBondProperties(std::ostream &ss, const RDProps &props,
                          unsigned int pickleFlags) {
  const static PropTracker bprops;
  if (!pickleFlags) {
    return false;
  }
  bool res = streamWriteProps<std::uint16_t>(
      ss, props, pickleFlags & PicklerOps::PrivateProps,
      pickleFlags & PicklerOps::ComputedProps,
      MolPickler::getCustomPropHandlers(), bprops.ignoreBondProps);
  res |= pickleExplicitProperties<std::int8_t>(ss, props,
                                               bprops.explicitBondProps);
  return res;
}

void unpickleBondProperties(std::istream &ss, RDProps &props, int version) {
  const static PropTracker bprops;
  if (version >= 14000) {
    streamReadProps<std::uint16_t>(ss, props,
                                   MolPickler::getCustomPropHandlers(), false);
  } else {
    streamReadProps<unsigned int>(ss, props,
                                  MolPickler::getCustomPropHandlers(), false);
  }
  unpickleExplicitProperties<std::int8_t, int>(ss, props, version,
                                               bprops.explicitBondProps);
}

}  // namespace

unsigned int MolPickler::getDefaultPickleProperties() {
#ifdef RDK_BUILD_THREADSAFE_SSS
  std::lock_guard<std::mutex> lock(GetPropMutex());
#endif
  unsigned int props = defaultProperties;
  return props;
}

void MolPickler::setDefaultPickleProperties(unsigned int props) {
#ifdef RDK_BUILD_THREADSAFE_SSS
  std::lock_guard<std::mutex> lock(GetPropMutex());
#endif
  defaultProperties = props;
}

const CustomPropHandlerVec &MolPickler::getCustomPropHandlers() {
#ifdef RDK_BUILD_THREADSAFE_SSS
  std::lock_guard<std::mutex> lock(GetPropMutex());
#endif
  if (defaultPropHandlers.size() == 0) {
    // initialize handlers
    defaultPropHandlers.push_back(
        std::make_shared<DataStructsExplicitBitVecPropHandler>());
  }
  return defaultPropHandlers;
}

void MolPickler::addCustomPropHandler(const CustomPropHandler &handler) {
#ifdef RDK_BUILD_THREADSAFE_SSS
  std::lock_guard<std::mutex> lock(GetPropMutex());
#endif
  if (defaultPropHandlers.size() == 0) {
    // initialize handlers
    defaultPropHandlers.push_back(
        std::make_shared<DataStructsExplicitBitVecPropHandler>());
  }
  defaultPropHandlers.push_back(
      std::shared_ptr<CustomPropHandler>(handler.clone()));
}

using namespace Queries;
namespace PicklerOps {

template <class T>
QueryDetails getQueryDetails(const Query<int, T const *, true> *query) {
  PRECONDITION(query, "no query");
  if (typeid(*query) == typeid(AndQuery<int, T const *, true>)) {
    return QueryDetails(MolPickler::QUERY_AND);
  } else if (typeid(*query) == typeid(OrQuery<int, T const *, true>)) {
    return QueryDetails(MolPickler::QUERY_OR);
  } else if (typeid(*query) == typeid(XOrQuery<int, T const *, true>)) {
    return QueryDetails(MolPickler::QUERY_XOR);
  } else if (typeid(*query) == typeid(EqualityQuery<int, T const *, true>)) {
    return QueryDetails(std::make_tuple(
        MolPickler::QUERY_EQUALS,
        static_cast<const EqualityQuery<int, T const *, true> *>(query)
            ->getVal(),
        static_cast<const EqualityQuery<int, T const *, true> *>(query)
            ->getTol()));
  } else if (typeid(*query) == typeid(GreaterQuery<int, T const *, true>)) {
    return QueryDetails(std::make_tuple(
        MolPickler::QUERY_GREATER,
        static_cast<const GreaterQuery<int, T const *, true> *>(query)
            ->getVal(),
        static_cast<const GreaterQuery<int, T const *, true> *>(query)
            ->getTol()));
  } else if (typeid(*query) ==
             typeid(GreaterEqualQuery<int, T const *, true>)) {
    return QueryDetails(std::make_tuple(
        MolPickler::QUERY_GREATEREQUAL,
        static_cast<const GreaterEqualQuery<int, T const *, true> *>(query)
            ->getVal(),
        static_cast<const GreaterEqualQuery<int, T const *, true> *>(query)
            ->getTol()));
  } else if (typeid(*query) == typeid(LessQuery<int, T const *, true>)) {
    return QueryDetails(std::make_tuple(
        MolPickler::QUERY_LESS,
        static_cast<const LessQuery<int, T const *, true> *>(query)->getVal(),
        static_cast<const LessQuery<int, T const *, true> *>(query)->getTol()));
  } else if (typeid(*query) == typeid(LessEqualQuery<int, T const *, true>)) {
    return QueryDetails(std::make_tuple(
        MolPickler::QUERY_LESSEQUAL,
        static_cast<const LessEqualQuery<int, T const *, true> *>(query)
            ->getVal(),
        static_cast<const LessEqualQuery<int, T const *, true> *>(query)
            ->getTol()));
  } else if (typeid(*query) == typeid(AtomRingQuery)) {
    return QueryDetails(std::make_tuple(
        MolPickler::QUERY_ATOMRING,
        static_cast<const EqualityQuery<int, T const *, true> *>(query)
            ->getVal(),
        static_cast<const EqualityQuery<int, T const *, true> *>(query)
            ->getTol()));
  } else if (typeid(*query) == typeid(HasPropQuery<T const *>)) {
    return QueryDetails(std::make_tuple(
        MolPickler::QUERY_PROPERTY,
        static_cast<const HasPropQuery<T const *> *>(query)->getPropName()));
  } else if (typeid(*query) == typeid(Query<int, T const *, true>)) {
    return QueryDetails(MolPickler::QUERY_NULL);
  } else if (typeid(*query) == typeid(RangeQuery<int, T const *, true>)) {
    char ends;
    bool lowerOpen, upperOpen;
    boost::tie(lowerOpen, upperOpen) =
        static_cast<const RangeQuery<int, T const *, true> *>(query)
            ->getEndsOpen();
    ends = 0 | (rdcast<int>(lowerOpen) << 1) | rdcast<int>(upperOpen);
    return QueryDetails(std::make_tuple(
        MolPickler::QUERY_RANGE,
        static_cast<const RangeQuery<int, T const *, true> *>(query)
            ->getLower(),
        static_cast<const RangeQuery<int, T const *, true> *>(query)
            ->getUpper(),
        static_cast<const RangeQuery<int, T const *, true> *>(query)->getTol(),
        ends));
  } else if (typeid(*query) == typeid(SetQuery<int, T const *, true>)) {
    std::set<int32_t> tset(
        static_cast<const SetQuery<int, T const *, true> *>(query)->beginSet(),
        static_cast<const SetQuery<int, T const *, true> *>(query)->endSet());
    return QueryDetails(
        std::make_tuple(MolPickler::QUERY_SET, std::move(tset)));
  } else {
    throw MolPicklerException("do not know how to pickle part of the query.");
  }
}
template RDKIT_GRAPHMOL_EXPORT QueryDetails getQueryDetails<RDKit::Atom>(
    const Queries::Query<int, RDKit::Atom const *, true> *query);
template RDKIT_GRAPHMOL_EXPORT QueryDetails getQueryDetails<RDKit::Bond>(
    const Queries::Query<int, RDKit::Bond const *, true> *query);

}  // namespace PicklerOps

namespace {
template <class T>
void pickleQuery(std::ostream &ss, const Query<int, T const *, true> *query) {
  PRECONDITION(query, "no query");
  streamWrite(ss, query->getDescription());
  if (!query->getTypeLabel().empty()) {
    streamWrite(ss, MolPickler::QUERY_TYPELABEL);
    streamWrite(ss, query->getTypeLabel());
  }
  if (query->getNegation()) {
    streamWrite(ss, MolPickler::QUERY_ISNEGATED);
  }
  if (typeid(*query) == typeid(RecursiveStructureQuery)) {
    streamWrite(ss, MolPickler::QUERY_RECURSIVE);
    streamWrite(ss, MolPickler::QUERY_VALUE);
    MolPickler::pickleMol(
        ((const RecursiveStructureQuery *)query)->getQueryMol(), ss);
  } else {
    auto qdetails = PicklerOps::getQueryDetails(query);
    switch (qdetails.which()) {
      case 0:
        streamWrite(ss, boost::get<MolPickler::Tags>(qdetails));
        break;
      case 1: {
        auto v = boost::get<std::tuple<MolPickler::Tags, int32_t>>(qdetails);
        streamWrite(ss, std::get<0>(v));
        streamWrite(ss, MolPickler::QUERY_VALUE, std::get<1>(v));
      } break;
      case 2: {
        auto v = boost::get<std::tuple<MolPickler::Tags, int32_t, int32_t>>(
            qdetails);
        streamWrite(ss, std::get<0>(v));
        streamWrite(ss, MolPickler::QUERY_VALUE, std::get<1>(v));
        streamWrite(ss, std::get<2>(v));
      } break;
      case 3: {
        auto v = boost::get<
            std::tuple<MolPickler::Tags, int32_t, int32_t, int32_t, char>>(
            qdetails);
        streamWrite(ss, std::get<0>(v));
        streamWrite(ss, MolPickler::QUERY_VALUE, std::get<1>(v));
        streamWrite(ss, std::get<2>(v));
        streamWrite(ss, std::get<3>(v));
        streamWrite(ss, std::get<4>(v));
      } break;
      case 4: {
        auto v = boost::get<std::tuple<MolPickler::Tags, std::set<int32_t>>>(
            qdetails);
        streamWrite(ss, std::get<0>(v));
        const auto &tset = std::get<1>(v);
        int32_t sz = tset.size();
        streamWrite(ss, MolPickler::QUERY_VALUE, sz);
        for (auto val : tset) {
          streamWrite(ss, val);
        }

      } break;
      case 5: {
        auto v =
            boost::get<std::tuple<MolPickler::Tags, std::string>>(qdetails);
        streamWrite(ss, std::get<0>(v));
        const auto &pval = std::get<1>(v);
        streamWrite(ss, MolPickler::QUERY_VALUE, pval);
      } break;
      default:
        throw MolPicklerException(
            "do not know how to pickle part of the query.");
    }
  }

  // now the children:
  streamWrite(ss, MolPickler::QUERY_NUMCHILDREN,
              static_cast<unsigned char>(query->endChildren() -
                                         query->beginChildren()));
  typename Query<int, T const *, true>::CHILD_VECT_CI cit;
  for (cit = query->beginChildren(); cit != query->endChildren(); ++cit) {
    pickleQuery(ss, cit->get());
  }
}

template <class T>
Query<int, T const *, true> *buildBaseQuery(std::istream &ss, T const *owner,
                                            MolPickler::Tags tag, int version) {
  PRECONDITION(owner, "no query");
  std::string descr;
  Query<int, T const *, true> *res = nullptr;
  int32_t val;
  int32_t nMembers;
  char cval;
  const unsigned int lowerOpen = 1 << 1;
  const unsigned int upperOpen = 1;
  switch (tag) {
    case MolPickler::QUERY_AND:
      res = new AndQuery<int, T const *, true>();
      break;
    case MolPickler::QUERY_OR:
      res = new OrQuery<int, T const *, true>();
      break;
    case MolPickler::QUERY_XOR:
      res = new XOrQuery<int, T const *, true>();
      break;
    case MolPickler::QUERY_EQUALS:
      res = new EqualityQuery<int, T const *, true>();
      streamRead(ss, tag, version);
      if (tag != MolPickler::QUERY_VALUE) {
        delete res;
        throw MolPicklerException(
            "Bad pickle format: QUERY_VALUE tag not found.");
      }
      streamRead(ss, val, version);
      static_cast<EqualityQuery<int, T const *, true> *>(res)->setVal(val);
      streamRead(ss, val, version);
      static_cast<EqualityQuery<int, T const *, true> *>(res)->setTol(val);
      break;
    case MolPickler::QUERY_GREATER:
      res = new GreaterQuery<int, T const *, true>();
      streamRead(ss, tag, version);
      if (tag != MolPickler::QUERY_VALUE) {
        delete res;
        throw MolPicklerException(
            "Bad pickle format: QUERY_VALUE tag not found.");
      }
      streamRead(ss, val, version);
      static_cast<GreaterQuery<int, T const *, true> *>(res)->setVal(val);
      streamRead(ss, val, version);
      static_cast<GreaterQuery<int, T const *, true> *>(res)->setTol(val);
      break;
    case MolPickler::QUERY_GREATEREQUAL:
      res = new GreaterEqualQuery<int, T const *, true>();
      streamRead(ss, tag, version);
      if (tag != MolPickler::QUERY_VALUE) {
        delete res;
        throw MolPicklerException(
            "Bad pickle format: QUERY_VALUE tag not found.");
      }
      streamRead(ss, val, version);
      static_cast<GreaterEqualQuery<int, T const *, true> *>(res)->setVal(val);
      streamRead(ss, val, version);
      static_cast<GreaterEqualQuery<int, T const *, true> *>(res)->setTol(val);
      break;
    case MolPickler::QUERY_LESS:
      res = new LessQuery<int, T const *, true>();
      streamRead(ss, tag, version);
      if (tag != MolPickler::QUERY_VALUE) {
        delete res;
        throw MolPicklerException(
            "Bad pickle format: QUERY_VALUE tag not found.");
      }
      streamRead(ss, val, version);
      static_cast<LessQuery<int, T const *, true> *>(res)->setVal(val);
      streamRead(ss, val, version);
      static_cast<LessQuery<int, T const *, true> *>(res)->setTol(val);
      break;
    case MolPickler::QUERY_LESSEQUAL:
      res = new LessEqualQuery<int, T const *, true>();
      streamRead(ss, tag, version);
      if (tag != MolPickler::QUERY_VALUE) {
        delete res;
        throw MolPicklerException(
            "Bad pickle format: QUERY_VALUE tag not found.");
      }
      streamRead(ss, val, version);
      static_cast<LessEqualQuery<int, T const *, true> *>(res)->setVal(val);
      streamRead(ss, val, version);
      static_cast<LessEqualQuery<int, T const *, true> *>(res)->setTol(val);
      break;
    case MolPickler::QUERY_RANGE:
      res = new RangeQuery<int, T const *, true>();
      streamRead(ss, tag, version);
      if (tag != MolPickler::QUERY_VALUE) {
        delete res;
        throw MolPicklerException(
            "Bad pickle format: QUERY_VALUE tag not found.");
      }
      streamRead(ss, val, version);
      static_cast<RangeQuery<int, T const *, true> *>(res)->setLower(val);
      streamRead(ss, val, version);
      static_cast<RangeQuery<int, T const *, true> *>(res)->setUpper(val);
      streamRead(ss, val, version);
      static_cast<RangeQuery<int, T const *, true> *>(res)->setTol(val);
      streamRead(ss, cval, version);
      static_cast<RangeQuery<int, T const *, true> *>(res)->setEndsOpen(
          cval & lowerOpen, cval & upperOpen);
      break;
    case MolPickler::QUERY_SET:
      res = new SetQuery<int, T const *, true>();
      streamRead(ss, tag, version);
      if (tag != MolPickler::QUERY_VALUE) {
        delete res;
        throw MolPicklerException(
            "Bad pickle format: QUERY_VALUE tag not found.");
      }
      streamRead(ss, nMembers);
      while (nMembers > 0) {
        streamRead(ss, val, version);
        static_cast<SetQuery<int, T const *, true> *>(res)->insert(val);
        --nMembers;
      }
      break;
    case MolPickler::QUERY_NULL:
      res = new Query<int, T const *, true>();
      break;
    case MolPickler::QUERY_PROPERTY: {
      streamRead(ss, tag, version);
      if (tag != MolPickler::QUERY_VALUE) {
        throw MolPicklerException(
            "Bad pickle format: QUERY_VALUE tag not found.");
      }
      std::string propName = "";
      streamRead(ss, propName, version);
      res = makeHasPropQuery<T>(propName);
    } break;
    default:
      throw MolPicklerException("unknown query-type tag encountered");
  }

  POSTCONDITION(res, "no match found");
  return res;
}

Query<int, Atom const *, true> *unpickleQuery(std::istream &ss,
                                              Atom const *owner, int version) {
  PRECONDITION(owner, "no query");
  std::string descr;
  std::string typeLabel = "";
  bool isNegated = false;
  Query<int, Atom const *, true> *res;
  streamRead(ss, descr, version);
  MolPickler::Tags tag;
  streamRead(ss, tag, version);
  if (tag == MolPickler::QUERY_TYPELABEL) {
    streamRead(ss, typeLabel, version);
    streamRead(ss, tag, version);
  }
  if (tag == MolPickler::QUERY_ISNEGATED) {
    isNegated = true;
    streamRead(ss, tag, version);
  }
  int32_t val;
  ROMol *tmpMol;
  switch (tag) {
    case MolPickler::QUERY_ATOMRING:
      streamRead(ss, tag, version);
      if (tag != MolPickler::QUERY_VALUE) {
        throw MolPicklerException(
            "Bad pickle format: QUERY_VALUE tag not found.");
      }
      res = new AtomRingQuery();
      streamRead(ss, val, version);
      static_cast<EqualityQuery<int, Atom const *, true> *>(res)->setVal(val);
      streamRead(ss, val, version);
      static_cast<EqualityQuery<int, Atom const *, true> *>(res)->setTol(val);
      break;
    case MolPickler::QUERY_RECURSIVE:
      streamRead(ss, tag, version);
      if (tag != MolPickler::QUERY_VALUE) {
        throw MolPicklerException(
            "Bad pickle format: QUERY_VALUE tag not found.");
      }
      tmpMol = new ROMol();
      MolPickler::molFromPickle(ss, tmpMol);
      res = new RecursiveStructureQuery(tmpMol);
      break;
    default:
      res = buildBaseQuery(ss, owner, tag, version);
      break;
  }
  CHECK_INVARIANT(res, "no query!");

  res->setNegation(isNegated);
  res->setDescription(descr);
  if (!typeLabel.empty()) {
    res->setTypeLabel(typeLabel);
  }

  QueryOps::finalizeQueryFromDescription(res, owner);

  // read in the children:
  streamRead(ss, tag, version);
  if (tag != MolPickler::QUERY_NUMCHILDREN) {
    throw MolPicklerException(
        "Bad pickle format: QUERY_NUMCHILDREN tag not found.");
  }
  unsigned char numChildren;
  streamRead(ss, numChildren, version);
  while (numChildren > 0) {
    Query<int, Atom const *, true> *child = unpickleQuery(ss, owner, version);
    res->addChild(Query<int, Atom const *, true>::CHILD_TYPE(child));
    --numChildren;
  }
  return res;
}
Query<int, Bond const *, true> *unpickleQuery(std::istream &ss,
                                              Bond const *owner, int version) {
  PRECONDITION(owner, "no query");
  std::string descr;
  std::string typeLabel = "";
  bool isNegated = false;
  Query<int, Bond const *, true> *res;
  streamRead(ss, descr, version);
  MolPickler::Tags tag;
  streamRead(ss, tag, version);
  if (tag == MolPickler::QUERY_TYPELABEL) {
    streamRead(ss, typeLabel, version);
    streamRead(ss, tag, version);
  }
  if (tag == MolPickler::QUERY_ISNEGATED) {
    isNegated = true;
    streamRead(ss, tag, version);
  }
  res = buildBaseQuery(ss, owner, tag, version);
  CHECK_INVARIANT(res, "no query!");

  res->setNegation(isNegated);
  res->setDescription(descr);

  QueryOps::finalizeQueryFromDescription(res, owner);

  // read in the children:
  streamRead(ss, tag, version);
  if (tag != MolPickler::QUERY_NUMCHILDREN) {
    throw MolPicklerException(
        "Bad pickle format: QUERY_NUMCHILDREN tag not found.");
  }
  unsigned char numChildren;
  streamRead(ss, numChildren, version);
  while (numChildren > 0) {
    Query<int, Bond const *, true> *child = unpickleQuery(ss, owner, version);
    res->addChild(Query<int, Bond const *, true>::CHILD_TYPE(child));
    --numChildren;
  }
  return res;
}

void pickleAtomPDBResidueInfo(std::ostream &ss,
                              const AtomPDBResidueInfo *info) {
  PRECONDITION(info, "no info");
  if (info->getSerialNumber()) {
    streamWrite(ss, MolPickler::ATOM_PDB_RESIDUE_SERIALNUMBER,
                info->getSerialNumber());
  }
  if (info->getAltLoc() != "") {
    streamWrite(ss, MolPickler::ATOM_PDB_RESIDUE_ALTLOC, info->getAltLoc());
  }
  if (info->getResidueName() != "") {
    streamWrite(ss, MolPickler::ATOM_PDB_RESIDUE_RESIDUENAME,
                info->getResidueName());
  }
  if (info->getResidueNumber()) {
    streamWrite(ss, MolPickler::ATOM_PDB_RESIDUE_RESIDUENUMBER,
                info->getResidueNumber());
  }
  if (info->getChainId() != "") {
    streamWrite(ss, MolPickler::ATOM_PDB_RESIDUE_CHAINID, info->getChainId());
  }
  if (info->getInsertionCode() != "") {
    streamWrite(ss, MolPickler::ATOM_PDB_RESIDUE_INSERTIONCODE,
                info->getInsertionCode());
  }
  if (info->getOccupancy()) {
    streamWrite(ss, MolPickler::ATOM_PDB_RESIDUE_OCCUPANCY,
                info->getOccupancy());
  }
  if (info->getTempFactor()) {
    streamWrite(ss, MolPickler::ATOM_PDB_RESIDUE_TEMPFACTOR,
                info->getTempFactor());
  }
  if (info->getIsHeteroAtom()) {
    streamWrite(ss, MolPickler::ATOM_PDB_RESIDUE_ISHETEROATOM,
                static_cast<char>(info->getIsHeteroAtom()));
  }
  if (info->getSecondaryStructure()) {
    streamWrite(ss, MolPickler::ATOM_PDB_RESIDUE_SECONDARYSTRUCTURE,
                info->getSecondaryStructure());
  }
  if (info->getSegmentNumber()) {
    streamWrite(ss, MolPickler::ATOM_PDB_RESIDUE_SEGMENTNUMBER,
                info->getSegmentNumber());
  }
}

void unpickleAtomPDBResidueInfo(std::istream &ss, AtomPDBResidueInfo *info,
                                int version) {
  PRECONDITION(info, "no info");
  std::string sval;
  double dval;
  char cval;
  unsigned int uival;
  int ival;
  MolPickler::Tags tag = MolPickler::BEGIN_ATOM_MONOMER;
  while (tag != MolPickler::END_ATOM_MONOMER) {
    streamRead(ss, tag, version);
    switch (tag) {
      case MolPickler::ATOM_PDB_RESIDUE_SERIALNUMBER:
        streamRead(ss, ival, version);
        info->setSerialNumber(ival);
        break;
      case MolPickler::ATOM_PDB_RESIDUE_ALTLOC:
        streamRead(ss, sval, version);
        info->setAltLoc(sval);
        break;
      case MolPickler::ATOM_PDB_RESIDUE_RESIDUENAME:
        streamRead(ss, sval, version);
        info->setResidueName(sval);
        break;
      case MolPickler::ATOM_PDB_RESIDUE_RESIDUENUMBER:
        streamRead(ss, ival, version);
        info->setResidueNumber(ival);
        break;
      case MolPickler::ATOM_PDB_RESIDUE_CHAINID:
        streamRead(ss, sval, version);
        info->setChainId(sval);
        break;
      case MolPickler::ATOM_PDB_RESIDUE_INSERTIONCODE:
        streamRead(ss, sval, version);
        info->setInsertionCode(sval);
        break;
      case MolPickler::ATOM_PDB_RESIDUE_OCCUPANCY:
        streamRead(ss, dval, version);
        info->setOccupancy(dval);
        break;
      case MolPickler::ATOM_PDB_RESIDUE_TEMPFACTOR:
        streamRead(ss, dval, version);
        info->setTempFactor(dval);
        break;
      case MolPickler::ATOM_PDB_RESIDUE_ISHETEROATOM:
        streamRead(ss, cval, version);
        info->setIsHeteroAtom(cval);
        break;
      case MolPickler::ATOM_PDB_RESIDUE_SECONDARYSTRUCTURE:
        streamRead(ss, uival, version);
        info->setSecondaryStructure(uival);
        break;
      case MolPickler::ATOM_PDB_RESIDUE_SEGMENTNUMBER:
        streamRead(ss, uival, version);
        info->setSegmentNumber(uival);
        break;
      case MolPickler::END_ATOM_MONOMER:
        break;
      default:
        throw MolPicklerException(
            "unrecognized tag while parsing atom peptide residue info");
    }
  }
}

void pickleAtomMonomerInfo(std::ostream &ss, const AtomMonomerInfo *info) {
  PRECONDITION(info, "no info");
  streamWrite(ss, info->getName());
  streamWrite(ss, static_cast<unsigned int>(info->getMonomerType()));
  switch (info->getMonomerType()) {
    case AtomMonomerInfo::UNKNOWN:
    case AtomMonomerInfo::OTHER:
      break;
    case AtomMonomerInfo::PDBRESIDUE:
      pickleAtomPDBResidueInfo(ss,
                               static_cast<const AtomPDBResidueInfo *>(info));
      break;
    default:
      throw MolPicklerException("unrecognized MonomerType");
  }
}
AtomMonomerInfo *unpickleAtomMonomerInfo(std::istream &ss, int version) {
  MolPickler::Tags tag;
  std::string nm;
  streamRead(ss, nm, version);
  unsigned int typ;
  streamRead(ss, typ, version);

  AtomMonomerInfo *res;
  switch (typ) {
    case AtomMonomerInfo::UNKNOWN:
    case AtomMonomerInfo::OTHER:
      streamRead(ss, tag, version);
      if (tag != MolPickler::END_ATOM_MONOMER) {
        throw MolPicklerException(
            "did not find expected end of atom monomer info");
      }
      res =
          new AtomMonomerInfo(RDKit::AtomMonomerInfo::AtomMonomerType(typ), nm);
      break;
    case AtomMonomerInfo::PDBRESIDUE:
      res = static_cast<AtomMonomerInfo *>(new AtomPDBResidueInfo(nm));
      unpickleAtomPDBResidueInfo(ss, static_cast<AtomPDBResidueInfo *>(res),
                                 version);
      break;
    default:
      throw MolPicklerException("unrecognized MonomerType");
  }
  return res;
}

}  // namespace

// Resets the `exceptionState` of the passed stream `ss` in the destructor to
// the `exceptionState` the stream ss was in, before setting
// `newExceptionState`.
struct IOStreamExceptionStateResetter {
  std::ios &originalStream;
  std::ios_base::iostate originalExceptionState;
  IOStreamExceptionStateResetter(std::ios &ss,
                                 std::ios_base::iostate newExceptionState)
      : originalStream(ss), originalExceptionState(ss.exceptions()) {
    ss.exceptions(newExceptionState);
  }

  ~IOStreamExceptionStateResetter() {
    if (originalStream) {
      originalStream.exceptions(originalExceptionState);
    }
  }
};

void MolPickler::pickleMol(const ROMol *mol, std::ostream &ss) {
  pickleMol(mol, ss, MolPickler::getDefaultPickleProperties());
}

void MolPickler::pickleMol(const ROMol *mol, std::ostream &ss,
                           unsigned int propertyFlags) {
  PRECONDITION(mol, "empty molecule");

  // Ensure that the exception state of the `ostream` is reset to the previous
  // state after we're done.
  // Also enable exceptions here, so we're notified when we've reached EOF or
  // any other problem.
  IOStreamExceptionStateResetter resetter(ss, std::ios_base::eofbit |
                                                  std::ios_base::failbit |
                                                  std::ios_base::badbit);

  try {
    streamWrite(ss, endianId);
    streamWrite(ss, static_cast<int>(VERSION));
    streamWrite(ss, versionMajor);
    streamWrite(ss, versionMinor);
    streamWrite(ss, versionPatch);
#ifndef OLD_PICKLE
    if (mol->getNumAtoms() > 255) {
      _pickle<int32_t>(mol, ss, propertyFlags);
    } else {
      _pickle<unsigned char>(mol, ss, propertyFlags);
    }
#else
    _pickleV1(mol, ss);
#endif
  } catch (const std::ios_base::failure &) {
    if (ss.eof()) {
      throw MolPicklerException(
          "Bad pickle format: unexpected End-of-File while writing");
    } else if (ss.bad()) {
      throw MolPicklerException("Bad pickle format: write error while writing");
    } else if (ss.fail()) {
      throw MolPicklerException(
          "Bad pickle format: logical error while writing");
    } else {
      throw MolPicklerException(
          "Bad pickle format: unexpected error while writing");
    }
  }
}

void MolPickler::pickleMol(const ROMol &mol, std::ostream &ss) {
  pickleMol(&mol, ss, MolPickler::getDefaultPickleProperties());
}

void MolPickler::pickleMol(const ROMol *mol, std::string &res) {
  pickleMol(mol, res, MolPickler::getDefaultPickleProperties());
}

void MolPickler::pickleMol(const ROMol *mol, std::string &res,
                           unsigned int pickleFlags) {
  PRECONDITION(mol, "empty molecule");
  std::stringstream ss(std::ios_base::binary | std::ios_base::out |
                       std::ios_base::in);
  MolPickler::pickleMol(mol, ss, pickleFlags);
  res = ss.str();
}

void MolPickler::pickleMol(const ROMol &mol, std::string &ss) {
  pickleMol(&mol, ss, MolPickler::getDefaultPickleProperties());
}

// NOTE: if the mol passed in here already has atoms and bonds, they will
// be left intact.  The side effect is that ALL atom and bond bookmarks
// will be blown out by the end of this process.
void MolPickler::molFromPickle(std::istream &ss, ROMol *mol,
                               unsigned int propertyFlags) {
  PRECONDITION(mol, "empty molecule");

  // Ensure that the exception state of the `istream` is reset to the previous
  // state after we're done.
  // Also enable exceptions here, so we're notified when we've reached EOF or
  // any other problem.
  IOStreamExceptionStateResetter resetter(ss, std::ios_base::eofbit |
                                                  std::ios_base::failbit |
                                                  std::ios_base::badbit);

  try {
    int32_t tmpInt;

    mol->clearAllAtomBookmarks();
    mol->clearAllBondBookmarks();

    streamRead(ss, tmpInt);
    if (tmpInt != endianId) {
      throw MolPicklerException(
          "Bad pickle format: bad endian ID or invalid file format");
    }

    streamRead(ss, tmpInt);
    if (static_cast<Tags>(tmpInt) != VERSION) {
      throw MolPicklerException("Bad pickle format: no version tag");
    }
    int32_t majorVersion, minorVersion, patchVersion;
    streamRead(ss, majorVersion);
    streamRead(ss, minorVersion);
    streamRead(ss, patchVersion);
    if (majorVersion > versionMajor ||
        (majorVersion == versionMajor && minorVersion > versionMinor)) {
      BOOST_LOG(rdWarningLog)
          << "Depickling from a version number (" << majorVersion << "."
          << minorVersion << ")"
          << "that is higher than our version (" << versionMajor << "."
          << versionMinor << ").\nThis probably won't work." << std::endl;
    }
    // version sanity checking
    if (majorVersion > 1000 || minorVersion > 100 || patchVersion > 100) {
      throw MolPicklerException("unreasonable version numbers");
    }
    majorVersion = 1000 * majorVersion + minorVersion * 10 + patchVersion;
    if (majorVersion == 1) {
      _depickleV1(ss, mol);
    } else {
      int32_t numAtoms;
      streamRead(ss, numAtoms, majorVersion);
      if (numAtoms > 255) {
        _depickle<int32_t>(ss, mol, majorVersion, numAtoms, propertyFlags);
      } else {
        _depickle<unsigned char>(ss, mol, majorVersion, numAtoms,
                                 propertyFlags);
      }
    }
    mol->clearAllAtomBookmarks();
    mol->clearAllBondBookmarks();
    if (majorVersion < 4000) {
      // FIX for issue 220 - probably better to change the pickle format later
      MolOps::assignStereochemistry(*mol, true);
    }
  } catch (const std::ios_base::failure &) {
    if (ss.eof()) {
      throw MolPicklerException(
          "Bad pickle format: unexpected End-of-File while reading");
    } else if (ss.bad()) {
      throw MolPicklerException("Bad pickle format: read error while reading");
    } else if (ss.fail()) {
      throw MolPicklerException(
          "Bad pickle format: logical error while reading");
    } else {
      throw MolPicklerException(
          "Bad pickle format: unexpected error while reading");
    }
  }
}
void MolPickler::molFromPickle(const std::string &pickle, ROMol *mol,
                               unsigned int propertyFlags) {
  PRECONDITION(mol, "empty molecule");
  std::stringstream ss(std::ios_base::binary | std::ios_base::out |
                       std::ios_base::in);
  ss.write(pickle.c_str(), pickle.length());
  MolPickler::molFromPickle(ss, mol, propertyFlags);
}

//--------------------------------------
//
//            Molecules
//
//--------------------------------------
template <typename T>
void MolPickler::_pickle(const ROMol *mol, std::ostream &ss,
                         unsigned int propertyFlags) {
  PRECONDITION(mol, "empty molecule");
  int32_t tmpInt;
  std::map<int, int> atomIdxMap;
  std::map<int, int> bondIdxMap;

  tmpInt = static_cast<int32_t>(mol->getNumAtoms());
  streamWrite(ss, tmpInt);
  tmpInt = static_cast<int32_t>(mol->getNumBonds());
  streamWrite(ss, tmpInt);

  unsigned char flag = 0x1 << 7;
  streamWrite(ss, flag);

  // -------------------
  //
  // Write Atoms
  //
  // -------------------
  streamWrite(ss, BEGINATOM);
  ROMol::ConstAtomIterator atIt;
  int nWritten = 0;
  for (atIt = mol->beginAtoms(); atIt != mol->endAtoms(); ++atIt) {
    _pickleAtom<T>(ss, *atIt);
    atomIdxMap[(*atIt)->getIdx()] = nWritten;
    nWritten++;
  }

  // -------------------
  //
  // Write Bonds
  //
  // -------------------
  streamWrite(ss, BEGINBOND);
  for (unsigned int i = 0; i < mol->getNumBonds(); i++) {
    auto bond = mol->getBondWithIdx(i);
    _pickleBond<T>(ss, bond, atomIdxMap);
    bondIdxMap[bond->getIdx()] = i;
  }

  // -------------------
  //
  // Write Rings (if present)
  //
  // -------------------
  const RingInfo *ringInfo = mol->getRingInfo();
  if (ringInfo && ringInfo->isInitialized()) {
    switch (ringInfo->getRingType()) {
      case RDKit::FIND_RING_TYPE::FIND_RING_TYPE_FAST:
        streamWrite(ss, BEGINFASTFIND);
        break;
      case RDKit::FIND_RING_TYPE::FIND_RING_TYPE_SSSR:
        streamWrite(ss, BEGINSSSR);
        break;
      case RDKit::FIND_RING_TYPE::FIND_RING_TYPE_SYMM_SSSR:
        streamWrite(ss, BEGINSYMMSSSR);
        break;
      default:
        streamWrite(ss, BEGINFINDOTHERORUNKNOWN);
        break;
    }
    _pickleSSSR<T>(ss, ringInfo, atomIdxMap);
  }

  // -------------------
  //
  // Write SubstanceGroups (if present)
  //
  // -------------------
  const auto &sgroups = getSubstanceGroups(*mol);
  if (!sgroups.empty()) {
    streamWrite(ss, BEGINSGROUP);

    tmpInt = static_cast<int32_t>(sgroups.size());
    streamWrite(ss, tmpInt);

    for (const auto &sgroup : sgroups) {
      _pickleSubstanceGroup<T>(ss, sgroup, atomIdxMap, bondIdxMap);
    }
  }
  // Write Stereo Groups
  {
    auto &stereo_groups = mol->getStereoGroups();
    if (stereo_groups.size() > 0u) {
      streamWrite(ss, BEGINSTEREOGROUP);
      _pickleStereo<T>(ss, stereo_groups, atomIdxMap, bondIdxMap);
    }
  }

  if (!(propertyFlags & PicklerOps::NoConformers)) {
    std::stringstream tss;
    if (propertyFlags & PicklerOps::CoordsAsDouble) {
      // pickle the conformations
      streamWrite(ss, BEGINCONFS_DOUBLE);
      tmpInt = static_cast<int32_t>(mol->getNumConformers());
      streamWrite(tss, tmpInt);
      for (auto ci = mol->beginConformers(); ci != mol->endConformers(); ++ci) {
        const Conformer *conf = ci->get();
        _pickleConformer<T, double>(tss, conf);
      }
    } else {
      // pickle the conformations
      streamWrite(ss, BEGINCONFS);
      tmpInt = static_cast<int32_t>(mol->getNumConformers());
      streamWrite(tss, tmpInt);
      for (auto ci = mol->beginConformers(); ci != mol->endConformers(); ++ci) {
        const Conformer *conf = ci->get();
        _pickleConformer<T, float>(tss, conf);
      }
    }

    if (propertyFlags & PicklerOps::MolProps) {
      streamWrite(tss, BEGINCONFPROPS);
      std::stringstream tss2;
      for (auto ci = mol->beginConformers(); ci != mol->endConformers(); ++ci) {
        const Conformer *conf = ci->get();
        _pickleProperties(tss2, *conf, propertyFlags);
      }
      write_sstream_to_stream(tss, tss2);
    }
    write_sstream_to_stream(ss, tss);
  }

  if (propertyFlags & PicklerOps::MolProps) {
    std::stringstream tss;
    _pickleProperties(tss, *mol, propertyFlags);
    if (!tss.str().empty()) {
      streamWrite(ss, BEGINPROPS);
      write_sstream_to_stream(ss, tss);
      streamWrite(ss, ENDPROPS);
    }
  }

  if (propertyFlags & PicklerOps::AtomProps) {
    std::stringstream tss;
    bool anyWritten = false;
    for (const auto atom : mol->atoms()) {
      anyWritten |= pickleAtomProperties(tss, *atom, propertyFlags);
    }
    if (anyWritten) {
      streamWrite(ss, BEGINATOMPROPS);
      write_sstream_to_stream(ss, tss);
      streamWrite(ss, ENDPROPS);
    }
  }

  if (propertyFlags & PicklerOps::BondProps) {
    std::stringstream tss;
    bool anyWritten = false;
    for (const auto bond : mol->bonds()) {
      anyWritten |= pickleBondProperties(tss, *bond, propertyFlags);
    }
    if (anyWritten) {
      streamWrite(ss, BEGINBONDPROPS);
      write_sstream_to_stream(ss, tss);
      streamWrite(ss, ENDPROPS);
    }
  }
  streamWrite(ss, ENDMOL);
}

template <typename T>
void MolPickler::_depickle(std::istream &ss, ROMol *mol, int version,
                           int numAtoms, unsigned int propertyFlags) {
  PRECONDITION(mol, "empty molecule");
  bool directMap = mol->getNumAtoms() == 0;
  Tags tag;
  int32_t tmpInt;
  // int numAtoms,numBonds;
  int numBonds;
  bool haveQuery = false;

  streamRead(ss, tmpInt, version);
  numBonds = tmpInt;

  // did we include coordinates
  bool includeCoords = false;
  if (version >= 3000) {
    unsigned char flag;
    streamRead(ss, flag, version);
    if (flag & 0x1 << 7) {
      includeCoords = true;
    }
  }
  // -------------------
  //
  // Read Atoms
  //
  // -------------------
  streamRead(ss, tag, version);
  if (tag != BEGINATOM) {
    throw MolPicklerException("Bad pickle format: BEGINATOM tag not found.");
  }
  Conformer *conf = nullptr;
  if ((version >= 2000 && version < 3000) && includeCoords) {
    // there can only one conformation - since the positions were stored on
    // the atoms themselves in this version
    conf = new Conformer(numAtoms);
    mol->addConformer(conf, true);
  }
  for (int i = 0; i < numAtoms; i++) {
    RDGeom::Point3D pos;
    Atom *atom = _addAtomFromPickle<T>(ss, mol, pos, version, directMap);
    if ((version >= 2000 && version < 3000) && includeCoords) {
      // this is a older pickle so we go the pos
      conf->setAtomPos(i, pos);
    }
    if (!directMap) {
      mol->setAtomBookmark(atom, i);
    }
    if (atom->hasQuery()) {
      haveQuery = true;
    }
  }

  // -------------------
  //
  // Read Bonds
  //
  // -------------------
  streamRead(ss, tag, version);
  if (tag != BEGINBOND) {
    throw MolPicklerException("Bad pickle format: BEGINBOND tag not found.");
  }
  for (int i = 0; i < numBonds; i++) {
    Bond *bond = _addBondFromPickle<T>(ss, mol, version, directMap);
    if (!directMap) {
      mol->setBondBookmark(bond, i);
    }
  }

  // -------------------
  //
  // Read Rings (if needed)
  //
  // -------------------
  streamRead(ss, tag, version);
  bool ringFound = false;
  FIND_RING_TYPE ringType =
      RDKit::FIND_RING_TYPE::FIND_RING_TYPE_OTHER_OR_UNKNOWN;
  if (tag == BEGINSSSR) {
    ringFound = true;
    ringType = FIND_RING_TYPE::FIND_RING_TYPE_SSSR;
  } else if (tag == BEGINSYMMSSSR) {
    ringFound = true;
    ringType = FIND_RING_TYPE::FIND_RING_TYPE_SYMM_SSSR;
  } else if (tag == BEGINFASTFIND) {
    ringFound = true;
    ringType = FIND_RING_TYPE::FIND_RING_TYPE_FAST;
  } else if (tag == BEGINFINDOTHERORUNKNOWN) {
    ringFound = true;
    ringType = FIND_RING_TYPE::FIND_RING_TYPE_OTHER_OR_UNKNOWN;
  }
  if (ringFound) {
    _addRingInfoFromPickle<T>(ss, mol, version, directMap, ringType);
    streamRead(ss, tag, version);
  }

  // -------------------
  //
  // Read SubstanceGroups (if needed)
  //
  // -------------------

  if (tag == BEGINSGROUP) {
    streamRead(ss, tmpInt, version);

    // Create SubstanceGroups
    for (int i = 0; i < tmpInt; ++i) {
      auto sgroup = _getSubstanceGroupFromPickle<T>(ss, mol, version);
      addSubstanceGroup(*mol, sgroup);
    }

    streamRead(ss, tag, version);
  }

  if (tag == BEGINSTEREOGROUP) {
    _depickleStereo<T>(ss, mol, version);
    streamRead(ss, tag, version);
  }

  if (tag == BEGINCONFS || tag == BEGINCONFS_DOUBLE) {
    int32_t blkSize = 0;
    if (version >= 13000) {
      streamRead(ss, blkSize, version);
    }
    if (version >= 13000 && (propertyFlags & PicklerOps::NoConformers)) {
      // not reading coordinates, so just skip over those bytes
      ss.seekg(blkSize, std::ios_base::cur);
      streamRead(ss, tag, version);
    } else {
      // read in the conformation
      streamRead(ss, tmpInt, version);
      std::vector<unsigned int> cids(tmpInt);
      for (auto i = 0; i < tmpInt; i++) {
        Conformer *conf;
        if (tag == BEGINCONFS) {
          conf = _conformerFromPickle<T, float>(ss, version);
        } else {
          conf = _conformerFromPickle<T, double>(ss, version);
        }
        mol->addConformer(conf);
        cids[i] = conf->getId();
      }
      streamRead(ss, tag, version);
      if (tag == BEGINCONFPROPS) {
        int32_t blkSize = 0;
        if (version >= 13000) {
          streamRead(ss, blkSize, version);
        }
        if (version >= 13000 && !(propertyFlags & PicklerOps::MolProps)) {
          ss.seekg(blkSize, std::ios_base::cur);
        } else {
          for (auto cid : cids) {
            _unpickleProperties(ss, mol->getConformer(cid), version);
          }
        }
        streamRead(ss, tag, version);
      }
    }
  }

  while (tag != ENDMOL) {
    if (tag == BEGINPROPS) {
      int32_t blkSize = 0;
      if (version >= 13000) {
        streamRead(ss, blkSize, version);
      }
      if (version >= 13000 && !(propertyFlags & PicklerOps::MolProps)) {
        ss.seekg(blkSize, std::ios_base::cur);
      } else {
        _unpickleProperties(ss, *mol, version);
      }
      streamRead(ss, tag, version);
    } else if (tag == BEGINATOMPROPS) {
      int32_t blkSize = 0;
      if (version >= 13000) {
        streamRead(ss, blkSize, version);
      }
      if (version >= 13000 && !(propertyFlags & PicklerOps::AtomProps)) {
        ss.seekg(blkSize, std::ios_base::cur);
      } else {
        for (const auto atom : mol->atoms()) {
          unpickleAtomProperties(ss, *atom, version);
        }
      }
      streamRead(ss, tag, version);
    } else if (tag == BEGINBONDPROPS) {
      int32_t blkSize = 0;
      if (version >= 13000) {
        streamRead(ss, blkSize, version);
      }
      if (version >= 13000 && !(propertyFlags & PicklerOps::BondProps)) {
        ss.seekg(blkSize, std::ios_base::cur);
      } else {
        for (const auto bond : mol->bonds()) {
          unpickleBondProperties(ss, *bond, version);
        }
      }
      streamRead(ss, tag, version);
    } else if (tag == BEGINQUERYATOMDATA) {
      for (const auto atom : mol->atoms()) {
        _unpickleAtomData(ss, atom, version);
      }
      streamRead(ss, tag, version);
    } else {
      break;  // break to tag != ENDMOL
    }
    if (tag != ENDPROPS) {
      throw MolPicklerException("Bad pickle format: ENDPROPS tag not found.");
    }
    streamRead(ss, tag, version);
  }
  if (tag != ENDMOL) {
    throw MolPicklerException("Bad pickle format: ENDMOL tag not found.");
  }

  if (haveQuery) {
    // we didn't read any property info for atoms with associated
    // queries. update their property caches
    // (was sf.net Issue 3316407)
    for (const auto atom : mol->atoms()) {
      if (atom->hasQuery()) {
        atom->updatePropertyCache(false);
      }
    }
  }
}

//--------------------------------------
//
//            Atoms
//
//--------------------------------------

namespace {
bool getAtomMapNumber(const Atom *atom, int &mapNum) {
  PRECONDITION(atom, "bad atom");
  if (!atom->hasProp(common_properties::molAtomMapNumber)) {
    return false;
  }
  bool res = true;
  int tmpInt;
  try {
    atom->getProp(common_properties::molAtomMapNumber, tmpInt);
  } catch (std::bad_any_cast &) {
    const std::string &tmpSVal =
        atom->getProp<std::string>(common_properties::molAtomMapNumber);
    try {
      tmpInt = boost::lexical_cast<int>(tmpSVal);
    } catch (boost::bad_lexical_cast &) {
      res = false;
    }
  }
  if (res) {
    mapNum = tmpInt;
  }
  return res;
}
}  // namespace

int32_t MolPickler::_pickleAtomData(std::ostream &tss, const Atom *atom) {
  int32_t propFlags = 0;
  char tmpChar;
  signed char tmpSchar;
  // tmpFloat=atom->getMass()-PeriodicTable::getTable()->getAtomicWeight(atom->getAtomicNum());
  // if(fabs(tmpFloat)>.0001){
  //   propFlags |= 1;
  //   streamWrite(tss,tmpFloat);
  // }
  tmpSchar = static_cast<signed char>(atom->getFormalCharge());
  if (tmpSchar != 0) {
    propFlags |= 1 << 1;
    streamWrite(tss, tmpSchar);
  }
  tmpChar = static_cast<char>(atom->getChiralTag());
  if (tmpChar != 0) {
    propFlags |= 1 << 2;
    streamWrite(tss, tmpChar);
  }
  tmpChar = static_cast<char>(atom->getHybridization());
  if (tmpChar != static_cast<char>(Atom::SP3)) {
    propFlags |= 1 << 3;
    streamWrite(tss, tmpChar);
  }

  tmpChar = static_cast<char>(atom->getNumExplicitHs());
  if (tmpChar != 0) {
    propFlags |= 1 << 4;
    streamWrite(tss, tmpChar);
  }
  if (atom->d_explicitValence > 0) {
    tmpChar = static_cast<char>(atom->d_explicitValence);
    propFlags |= 1 << 5;
    streamWrite(tss, tmpChar);
  }
  if (atom->d_implicitValence > 0) {
    tmpChar = static_cast<char>(atom->d_implicitValence);
    propFlags |= 1 << 6;
    streamWrite(tss, tmpChar);
  }
  tmpChar = static_cast<char>(atom->getNumRadicalElectrons());
  if (tmpChar != 0) {
    propFlags |= 1 << 7;
    streamWrite(tss, tmpChar);
  }

  unsigned int tmpuint = atom->getIsotope();
  if (tmpuint > 0) {
    propFlags |= 1 << 8;
    streamWrite(tss, tmpuint);
  }
  return propFlags;
}

void MolPickler::_unpickleAtomData(std::istream &ss, Atom *atom, int version) {
  int propFlags;
  char tmpChar;
  signed char tmpSchar;

  streamRead(ss, propFlags, version);
  if (propFlags & 1) {
    float tmpFloat;
    streamRead(ss, tmpFloat, version);
    int iso = static_cast<int>(floor(tmpFloat + atom->getMass() + .0001));
    atom->setIsotope(iso);
  }

  if (propFlags & (1 << 1)) {
    streamRead(ss, tmpSchar, version);
  } else {
    tmpSchar = 0;
  }
  atom->setFormalCharge(static_cast<int>(tmpSchar));

  if (propFlags & (1 << 2)) {
    streamReadPositiveChar(ss, tmpChar, version);
  } else {
    tmpChar = 0;
  }
  atom->setChiralTag(static_cast<Atom::ChiralType>(tmpChar));

  if (propFlags & (1 << 3)) {
    streamReadPositiveChar(ss, tmpChar, version);
  } else {
    tmpChar = Atom::SP3;
  }
  atom->setHybridization(static_cast<Atom::HybridizationType>(tmpChar));

  if (propFlags & (1 << 4)) {
    streamRead(ss, tmpChar, version);
  } else {
    tmpChar = 0;
  }
  atom->setNumExplicitHs(tmpChar);

  if (propFlags & (1 << 5)) {
    streamRead(ss, tmpChar, version);
  } else {
    tmpChar = 0;
  }
  atom->d_explicitValence = tmpChar;

  if (propFlags & (1 << 6)) {
    streamRead(ss, tmpChar, version);
  } else {
    tmpChar = 0;
  }
  atom->d_implicitValence = tmpChar;
  if (propFlags & (1 << 7)) {
    streamReadPositiveChar(ss, tmpChar, version);
  } else {
    tmpChar = 0;
  }
  atom->d_numRadicalElectrons = static_cast<unsigned int>(tmpChar);

  atom->d_isotope = 0;
  if (propFlags & (1 << 8)) {
    unsigned int tmpuint;
    streamRead(ss, tmpuint, version);
    atom->setIsotope(tmpuint);
  }
}

// T refers to the type of the atom indices written
template <typename T>
void MolPickler::_pickleAtom(std::ostream &ss, const Atom *atom) {
  PRECONDITION(atom, "empty atom");
  char tmpChar;
  unsigned char tmpUchar;
  int tmpInt;
  char flags;

  tmpUchar = atom->getAtomicNum() % 256;
  streamWrite(ss, tmpUchar);

  flags = 0;
  if (atom->getIsAromatic()) {
    flags |= 0x1 << 6;
  }
  if (atom->getNoImplicit()) {
    flags |= 0x1 << 5;
  }
  if (atom->hasQuery()) {
    flags |= 0x1 << 4;
  }
  if (getAtomMapNumber(atom, tmpInt)) {
    flags |= 0x1 << 3;
  }
  if (atom->hasProp(common_properties::dummyLabel)) {
    flags |= 0x1 << 2;
  }
  if (atom->getMonomerInfo()) {
    flags |= 0x1 << 1;
  }

  streamWrite(ss, flags);

  std::stringstream tss(std::ios_base::binary | std::ios_base::out |
                        std::ios_base::in);
  int32_t propFlags = _pickleAtomData(tss, atom);
  streamWrite(ss, propFlags);
  ss.write(tss.str().c_str(), tss.str().size());
  if (atom->hasQuery()) {
    streamWrite(ss, BEGINQUERY);
    pickleQuery(ss, static_cast<const QueryAtom *>(atom)->getQuery());
    streamWrite(ss, ENDQUERY);
  }
  if (getAtomMapNumber(atom, tmpInt)) {
    if (tmpInt >= 0 && tmpInt < 128) {
      tmpChar = static_cast<char>(tmpInt);
      streamWrite(ss, ATOM_MAPNUMBER, tmpChar);
    } else {
      tmpChar = static_cast<char>(255);
      streamWrite(ss, ATOM_MAPNUMBER, tmpChar);
      streamWrite(ss, tmpInt);
    }
  }
  if (atom->hasProp(common_properties::dummyLabel)) {
    streamWrite(ss, ATOM_DUMMYLABEL,
                atom->getProp<std::string>(common_properties::dummyLabel));
  }
  if (atom->getMonomerInfo()) {
    streamWrite(ss, BEGIN_ATOM_MONOMER);
    pickleAtomMonomerInfo(ss, atom->getMonomerInfo());
    streamWrite(ss, END_ATOM_MONOMER);
  }
}

template <typename T, typename C>
void MolPickler::_pickleConformer(std::ostream &ss, const Conformer *conf) {
  PRECONDITION(conf, "empty conformer");
  char tmpChr = static_cast<int>(conf->is3D());
  streamWrite(ss, tmpChr);
  auto tmpInt = static_cast<int32_t>(conf->getId());
  streamWrite(ss, tmpInt);
  T tmpT = static_cast<T>(conf->getNumAtoms());
  streamWrite(ss, tmpT);
  const RDGeom::POINT3D_VECT &pts = conf->getPositions();
  for (const auto &pt : pts) {
    C tmpFloat;
    tmpFloat = static_cast<C>(pt.x);
    streamWrite(ss, tmpFloat);
    tmpFloat = static_cast<C>(pt.y);
    streamWrite(ss, tmpFloat);
    tmpFloat = static_cast<C>(pt.z);
    streamWrite(ss, tmpFloat);
  }
}

template <typename T, typename C>
Conformer *MolPickler::_conformerFromPickle(std::istream &ss, int version) {
  C tmpFloat;
  bool is3D = true;
  if (version > 4000) {
    char tmpChr;
    streamRead(ss, tmpChr, version);
    is3D = static_cast<bool>(tmpChr);
  }
  int tmpInt;
  streamRead(ss, tmpInt, version);
  auto cid = static_cast<unsigned int>(tmpInt);
  T tmpT;
  streamRead(ss, tmpT, version);
  auto numAtoms = static_cast<unsigned int>(tmpT);
  auto *conf = new Conformer(numAtoms);
  conf->setId(cid);
  conf->set3D(is3D);
  try {
    for (unsigned int i = 0; i < numAtoms; i++) {
      streamRead(ss, tmpFloat, version);
      conf->getAtomPos(i).x = static_cast<double>(tmpFloat);
      streamRead(ss, tmpFloat, version);
      conf->getAtomPos(i).y = static_cast<double>(tmpFloat);
      streamRead(ss, tmpFloat, version);
      conf->getAtomPos(i).z = static_cast<double>(tmpFloat);
    }
  } catch (...) {
    delete conf;
    throw;
  }
  return conf;
}

template <typename T>
Atom *MolPickler::_addAtomFromPickle(std::istream &ss, ROMol *mol,
                                     RDGeom::Point3D &pos, int version, bool) {
  PRECONDITION(mol, "empty molecule");
  float x, y, z;
  char tmpChar;
  unsigned char tmpUchar;
  signed char tmpSchar;
  char flags;
  Tags tag;
  Atom *atom = nullptr;
  int atomicNum = 0;

  streamRead(ss, tmpUchar, version);
  atomicNum = tmpUchar;

  bool hasQuery = false;
  streamRead(ss, flags, version);
  if (version > 5000) {
    hasQuery = flags & 0x1 << 4;
  }
  if (!hasQuery) {
    atom = new Atom(atomicNum);
  } else {
    atom = new QueryAtom();
    if (atomicNum) {
      // can't set this in the constructor because that builds a
      // query and we're going to take care of that later:
      atom->setAtomicNum(atomicNum);
    }
  }
  atom->setIsAromatic(flags & 0x1 << 6);
  atom->setNoImplicit(flags & 0x1 << 5);

  bool hasAtomMap = 0, hasDummyLabel = 0;
  if (version >= 6020) {
    hasAtomMap = flags & 0x1 << 3;
    hasDummyLabel = flags & 0x1 << 2;
  }
  bool hasMonomerInfo = 0;
  if (version >= 7020) {
    hasMonomerInfo = flags & 0x1 << 1;
  }

  // are coordinates present?
  if (flags & 0x1 << 7) {
    streamRead(ss, x, version);
    pos.x = static_cast<double>(x);
    streamRead(ss, y, version);
    pos.y = static_cast<double>(y);
    streamRead(ss, z, version);
    pos.z = static_cast<double>(z);
  }

  if (version <= 5000 || !hasQuery) {
    if (version < 7000) {
      if (version < 6030) {
        streamRead(ss, tmpSchar, version);
        // FIX: technically should be handling this in order to maintain true
        // backwards compatibility
        // atom->setMass(PeriodicTable::getTable()->getAtomicWeight(atom->getAtomicNum())+
        //              static_cast<int>(tmpSchar));
      } else {
        float tmpFloat;
        streamRead(ss, tmpFloat, version);
        // FIX: technically should be handling this in order to maintain true
        // backwards compatibility
        // atom->setMass(tmpFloat);
      }

      streamRead(ss, tmpSchar, version);
      atom->setFormalCharge(static_cast<int>(tmpSchar));

      streamRead(ss, tmpChar, version);
      atom->setChiralTag(static_cast<Atom::ChiralType>(tmpChar));
      streamRead(ss, tmpChar, version);
      atom->setHybridization(static_cast<Atom::HybridizationType>(tmpChar));
      streamRead(ss, tmpChar, version);
      atom->setNumExplicitHs(static_cast<int>(tmpChar));
      streamRead(ss, tmpChar, version);
      atom->d_explicitValence = tmpChar;
      streamRead(ss, tmpChar, version);
      atom->d_implicitValence = tmpChar;
      if (version > 6000) {
        streamRead(ss, tmpChar, version);
        atom->d_numRadicalElectrons = static_cast<unsigned int>(tmpChar);
      }
    } else {
      _unpickleAtomData(ss, atom, version);
    }

  } else if (version > 5000) {
    // we have a query
    if (version >= 9000) {
      _unpickleAtomData(ss, atom, version);
    }
    streamRead(ss, tag, version);
    if (tag != BEGINQUERY) {
      throw MolPicklerException("Bad pickle format: BEGINQUERY tag not found.");
    }
    static_cast<QueryAtom *>(atom)->setQuery(unpickleQuery(ss, atom, version));
    streamRead(ss, tag, version);
    if (tag != ENDQUERY) {
      throw MolPicklerException("Bad pickle format: ENDQUERY tag not found.");
    }
    // atom->setNumExplicitHs(0);
  }

  if (version > 5000) {
    if (version < 6020) {
      unsigned int sPos = rdcast<unsigned int>(ss.tellg());
      Tags tag;
      streamRead(ss, tag, version);
      if (tag == ATOM_MAPNUMBER) {
        int32_t tmpInt;
        streamRead(ss, tmpChar, version);
        tmpInt = tmpChar;
        atom->setProp(common_properties::molAtomMapNumber, tmpInt);
      } else {
        ss.seekg(sPos);
      }
    } else {
      if (hasAtomMap) {
        Tags tag;
        streamRead(ss, tag, version);
        if (tag != ATOM_MAPNUMBER) {
          throw MolPicklerException(
              "Bad pickle format: ATOM_MAPNUMBER tag not found.");
        }
        int tmpInt;
        streamRead(ss, tmpChar, version);
        // the test for tmpChar below seems redundant, but on at least
        // the POWER8 architecture it seems that chars may be unsigned
        // by default.
#ifdef __GNUC__
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wtype-limits"
#endif
        if ((tmpChar < 0 || tmpChar > 127) && version > 9000) {
#ifdef __GNUC__
#pragma GCC diagnostic pop
#endif
          streamRead(ss, tmpInt, version);
        } else {
          tmpInt = tmpChar;
        }
        atom->setProp(common_properties::molAtomMapNumber, tmpInt);
      }
      if (hasDummyLabel) {
        streamRead(ss, tag, version);
        if (tag != ATOM_DUMMYLABEL) {
          throw MolPicklerException(
              "Bad pickle format: ATOM_DUMMYLABEL tag not found.");
        }
        std::string tmpStr;
        streamRead(ss, tmpStr, version);
        atom->setProp(common_properties::dummyLabel, tmpStr);
      }
    }
  }
  if (version >= 7020) {
    if (hasMonomerInfo) {
      streamRead(ss, tag, version);
      if (tag != BEGIN_ATOM_MONOMER) {
        throw MolPicklerException(
            "Bad pickle format: BEGIN_ATOM_MONOMER tag not found.");
      }
      atom->setMonomerInfo(unpickleAtomMonomerInfo(ss, version));
    }
  }
  mol->addAtom(atom, false, true);
  return atom;
}

//--------------------------------------
//
//            Bonds
//
//--------------------------------------

template <typename T>
void MolPickler::_pickleBond(std::ostream &ss, const Bond *bond,
                             std::map<int, int> &atomIdxMap) {
  PRECONDITION(bond, "empty bond");
  T tmpT;
  char tmpChar;
  char flags;

  tmpT = static_cast<T>(atomIdxMap[bond->getBeginAtomIdx()]);
  streamWrite(ss, tmpT);
  tmpT = static_cast<T>(atomIdxMap[bond->getEndAtomIdx()]);
  streamWrite(ss, tmpT);

  flags = 0;
  if (bond->getIsAromatic()) {
    flags |= 0x1 << 6;
  }
  if (bond->getIsConjugated()) {
    flags |= 0x1 << 5;
  }
  if (bond->hasQuery()) {
    flags |= 0x1 << 4;
  }
  if (bond->getBondType() != Bond::SINGLE) {
    flags |= 0x1 << 3;
  }
  if (bond->getBondDir() != Bond::NONE) {
    flags |= 0x1 << 2;
  }
  if (bond->getStereo() != Bond::STEREONONE) {
    flags |= 0x1 << 1;
  }
  std::string endpts;
  if (bond->getPropIfPresent(common_properties::_MolFileBondEndPts, endpts)) {
    flags |= 0x1;
  }
  streamWrite(ss, flags);

  if (bond->getBondType() != Bond::SINGLE) {
    tmpChar = static_cast<char>(bond->getBondType());
    streamWrite(ss, tmpChar);
  }
  if (bond->getBondDir() != Bond::NONE) {
    tmpChar = static_cast<char>(bond->getBondDir());
    streamWrite(ss, tmpChar);
  }

  // write info about the stereochemistry:
  if (bond->getStereo() != Bond::STEREONONE) {
    tmpChar = static_cast<char>(bond->getStereo());
    streamWrite(ss, tmpChar);
    const INT_VECT &stereoAts = bond->getStereoAtoms();
    tmpChar = rdcast<unsigned int>(stereoAts.size());
    streamWrite(ss, tmpChar);
    for (int stereoAt : stereoAts) {
      tmpT = static_cast<T>(stereoAt);
      streamWrite(ss, tmpT);
    }
  }
  if (bond->hasQuery()) {
    streamWrite(ss, BEGINQUERY);
    pickleQuery(ss, static_cast<const QueryBond *>(bond)->getQuery());
    streamWrite(ss, ENDQUERY);
  }
  if (!endpts.empty()) {
    streamWrite(ss, endpts);
    std::string attach = "ALL";
    bond->getPropIfPresent(common_properties::_MolFileBondAttach, attach);
    streamWrite(ss, attach);
  }
}

template <typename T>
Bond *MolPickler::_addBondFromPickle(std::istream &ss, ROMol *mol, int version,
                                     bool directMap) {
  PRECONDITION(mol, "empty molecule");
  char tmpChar;
  char flags;
  int begIdx, endIdx;
  T tmpT;

  Bond *bond = nullptr;
  streamRead(ss, tmpT, version);
  if (directMap) {
    begIdx = tmpT;
  } else {
    begIdx = mol->getAtomWithBookmark(static_cast<int>(tmpT))->getIdx();
  }
  streamRead(ss, tmpT, version);
  if (directMap) {
    endIdx = tmpT;

  } else {
    endIdx = mol->getAtomWithBookmark(static_cast<int>(tmpT))->getIdx();
  }
  streamRead(ss, flags, version);
  bool hasQuery = flags & 0x1 << 4;

  if (version <= 5000 || (version <= 7000 && !hasQuery) || version > 7000) {
    bond = new Bond();
    bond->setIsAromatic(flags & 0x1 << 6);
    bond->setIsConjugated(flags & 0x1 << 5);

    if (version < 7000) {
      streamReadPositiveChar(ss, tmpChar, version);
      bond->setBondType(static_cast<Bond::BondType>(tmpChar));
      streamReadPositiveChar(ss, tmpChar, version);
      bond->setBondDir(static_cast<Bond::BondDir>(tmpChar));

      if (version > 3000) {
        streamReadPositiveChar(ss, tmpChar, version);
        auto stereo = static_cast<Bond::BondStereo>(tmpChar);
        bond->setStereo(stereo);
        if (stereo != Bond::STEREONONE) {
          streamRead(ss, tmpChar, version);
          for (char i = 0; i < tmpChar; ++i) {
            streamRead(ss, tmpT, version);
            bond->getStereoAtoms().push_back(static_cast<int>(tmpT));
          }
        }
      }
    } else {
      if (flags & (0x1 << 3)) {
        streamReadPositiveChar(ss, tmpChar, version);
        bond->setBondType(static_cast<Bond::BondType>(tmpChar));
      } else {
        bond->setBondType(Bond::SINGLE);
      }

      if (flags & (0x1 << 2)) {
        streamReadPositiveChar(ss, tmpChar, version);
        bond->setBondDir(static_cast<Bond::BondDir>(tmpChar));
      } else {
        bond->setBondDir(Bond::NONE);
      }

      if (flags & (0x1 << 1)) {
        streamReadPositiveChar(ss, tmpChar, version);
        auto stereo = static_cast<Bond::BondStereo>(tmpChar);
        streamRead(ss, tmpChar, version);
        for (char i = 0; i < tmpChar; ++i) {
          streamRead(ss, tmpT, version);
          bond->getStereoAtoms().push_back(static_cast<int>(tmpT));
        }
        bond->setStereo(stereo);
      } else {
        bond->setStereo(Bond::STEREONONE);
      }
      if (flags & 0x1) {
        std::string tmpStr;
        streamRead(ss, tmpStr, version);
        bond->setProp(common_properties::_MolFileBondEndPts, tmpStr);
        streamRead(ss, tmpStr, version);
        bond->setProp(common_properties::_MolFileBondAttach, tmpStr);
      }
    }
  }
  if (version > 5000 && hasQuery) {
    Tags tag;
    if (bond) {
      Bond *tbond = bond;
      bond = new QueryBond(*bond);
      delete tbond;
    } else {
      bond = new QueryBond();
    }

    // we have a query:
    streamRead(ss, tag, version);
    if (tag != BEGINQUERY) {
      delete bond;
      throw MolPicklerException("Bad pickle format: BEGINQUERY tag not found.");
    }
    static_cast<QueryBond *>(bond)->setQuery(unpickleQuery(ss, bond, version));
    streamRead(ss, tag, version);
    if (tag != ENDQUERY) {
      delete bond;
      throw MolPicklerException("Bad pickle format: ENDQUERY tag not found.");
    }
  }
  if (bond) {
    bond->setBeginAtomIdx(begIdx);
    bond->setEndAtomIdx(endIdx);
    mol->addBond(bond, true);
  }
  return bond;
}

//--------------------------------------
//
//            Rings
//
//--------------------------------------
template <typename T>
void MolPickler::_pickleSSSR(std::ostream &ss, const RingInfo *ringInfo,
                             std::map<int, int> &atomIdxMap) {
  PRECONDITION(ringInfo, "missing ring info");
  std::uint32_t nrings = ringInfo->numRings();
  streamWrite(ss, nrings);
  for (unsigned int i = 0; i < ringInfo->numRings(); i++) {
    INT_VECT ring;
    ring = ringInfo->atomRings()[i];
    T tmpT = static_cast<T>(ring.size());
    streamWrite(ss, tmpT);
    for (int &j : ring) {
      tmpT = static_cast<T>(atomIdxMap[j]);
      streamWrite(ss, tmpT);
    }
  }
}

template <typename T>
void MolPickler::_addRingInfoFromPickle(std::istream &ss, ROMol *mol,
                                        int version, bool directMap,
                                        FIND_RING_TYPE ringType) {
  PRECONDITION(mol, "empty molecule");
  RingInfo *ringInfo = mol->getRingInfo();
  // if (!ringInfo->isInitialized()) {
  ringInfo->initialize(ringType);
  //}

  std::uint32_t numRings;
  if (version >= 13002) {
    streamRead(ss, numRings, version);
  } else {
    T tmpV;
    streamRead(ss, tmpV, version);
    numRings = tmpV;
  }

  if (numRings > 0) {
    ringInfo->preallocate(mol->getNumAtoms(), mol->getNumBonds());
    for (unsigned int i = 0; i < static_cast<unsigned int>(numRings); i++) {
      T tmpT;
      T ringSize;
      streamRead(ss, ringSize, version);
      if (ringSize < 0) {
        throw MolPicklerException("negative ring size");
      }
      INT_VECT atoms(static_cast<int>(ringSize));
      INT_VECT bonds(static_cast<int>(ringSize));
      for (unsigned int j = 0; j < static_cast<unsigned int>(ringSize); j++) {
        streamRead(ss, tmpT, version);
        if (directMap) {
          atoms[j] = static_cast<int>(tmpT);
          if (atoms[j] < 0 ||
              static_cast<unsigned int>(atoms[j]) >= mol->getNumAtoms()) {
            throw MolPicklerException("ring-atom index out of range");
          }
        } else {
          atoms[j] = mol->getAtomWithBookmark(static_cast<int>(tmpT))->getIdx();
        }
      }
      if (version < 7000) {
        for (unsigned int j = 0; j < static_cast<unsigned int>(ringSize); j++) {
          streamRead(ss, tmpT, version);
          if (directMap) {
            bonds[j] = static_cast<int>(tmpT);
            if (bonds[j] < 0 ||
                static_cast<unsigned int>(bonds[j]) >= mol->getNumBonds()) {
              throw MolPicklerException("ring-bond index out of range");
            }

          } else {
            bonds[j] =
                mol->getBondWithBookmark(static_cast<int>(tmpT))->getIdx();
          }
        }
      } else {
        for (unsigned int j = 1; j < static_cast<unsigned int>(ringSize); ++j) {
          bonds[j - 1] =
              mol->getBondBetweenAtoms(atoms[j - 1], atoms[j])->getIdx();
        }
        bonds[ringSize - 1] =
            mol->getBondBetweenAtoms(atoms[0], atoms[ringSize - 1])->getIdx();
      }
      ringInfo->addRing(atoms, bonds);
    }
  }
}

//--------------------------------------
//
//            SubstanceGroups
//
//--------------------------------------

template <typename T>
void MolPickler::_pickleSubstanceGroup(std::ostream &ss,
                                       const SubstanceGroup &sgroup,
                                       std::map<int, int> &atomIdxMap,
                                       std::map<int, int> &bondIdxMap) {
  T tmpT;

  streamWriteProps<std::uint16_t>(ss, sgroup);

  const auto &atoms = sgroup.getAtoms();
  streamWrite(ss, static_cast<T>(atoms.size()));
  for (const auto &atom : atoms) {
    tmpT = static_cast<T>(atomIdxMap[atom]);
    streamWrite(ss, tmpT);
  }

  const auto &p_atoms = sgroup.getParentAtoms();
  streamWrite(ss, static_cast<T>(p_atoms.size()));
  for (const auto &p_atom : p_atoms) {
    tmpT = static_cast<T>(atomIdxMap[p_atom]);
    streamWrite(ss, tmpT);
  }

  const auto &bonds = sgroup.getBonds();
  streamWrite(ss, static_cast<T>(bonds.size()));
  for (const auto &bond : bonds) {
    tmpT = static_cast<T>(bondIdxMap[bond]);
    streamWrite(ss, tmpT);
  }

  const auto &brackets = sgroup.getBrackets();
  streamWrite(ss, static_cast<T>(brackets.size()));
  for (const auto &bracket : brackets) {
    // 3 point per bracket; 3rd point and all z are zeros,
    // but this might change in the future.
    for (const auto &pt : bracket) {
      float tmpFloat;
      tmpFloat = static_cast<float>(pt.x);
      streamWrite(ss, tmpFloat);
      tmpFloat = static_cast<float>(pt.y);
      streamWrite(ss, tmpFloat);
      tmpFloat = static_cast<float>(pt.z);
      streamWrite(ss, tmpFloat);
    }
  }

  const auto &cstates = sgroup.getCStates();
  streamWrite(ss, static_cast<T>(cstates.size()));
  for (const auto &cstate : cstates) {
    // Bond
    tmpT = static_cast<T>(bondIdxMap[cstate.bondIdx]);
    streamWrite(ss, tmpT);

    // Vector -- existence depends on SubstanceGroup type
    if ("SUP" == sgroup.getProp<std::string>("TYPE")) {
      float tmpFloat;
      tmpFloat = static_cast<float>(cstate.vector.x);
      streamWrite(ss, tmpFloat);
      tmpFloat = static_cast<float>(cstate.vector.y);
      streamWrite(ss, tmpFloat);
      tmpFloat = static_cast<float>(cstate.vector.z);
      streamWrite(ss, tmpFloat);
    }
  }

  const auto &attachPoints = sgroup.getAttachPoints();
  streamWrite(ss, static_cast<T>(attachPoints.size()));
  for (const auto &attachPoint : attachPoints) {
    // aIdx -- always present
    tmpT = static_cast<T>(atomIdxMap[attachPoint.aIdx]);
    streamWrite(ss, tmpT);

    // lvIdx -- may be -1 if not used (0 in spec)
    auto tmpInt = -1;
    if (attachPoint.lvIdx != -1) {
      tmpInt = static_cast<signed int>(atomIdxMap[attachPoint.lvIdx]);
    }
    streamWrite(ss, tmpInt);

    // id -- may be blank
    auto tmpS = static_cast<const std::string>(attachPoint.id);
    streamWrite(ss, tmpS);
  }
}

template <typename T>
SubstanceGroup MolPickler::_getSubstanceGroupFromPickle(std::istream &ss,
                                                        ROMol *mol,
                                                        int version) {
  T tmpT;
  T numItems;
  float tmpFloat;
  int tmpInt = -1;
  std::string tmpS;

  // temporarily accept empty TYPE
  SubstanceGroup sgroup(mol, "");

  // Read RDProps, overriding ID, TYPE and COMPNO
  if (version >= 14000) {
    streamReadProps<std::uint16_t>(ss, sgroup,
                                   MolPickler::getCustomPropHandlers());
  } else {
    streamReadProps<unsigned int>(ss, sgroup,
                                  MolPickler::getCustomPropHandlers());
  }

  streamRead(ss, numItems, version);
  for (int i = 0; i < numItems; ++i) {
    streamRead(ss, tmpT, version);
    sgroup.addAtomWithIdx(tmpT);
  }

  streamRead(ss, numItems, version);
  for (int i = 0; i < numItems; ++i) {
    streamRead(ss, tmpT, version);
    sgroup.addParentAtomWithIdx(tmpT);
  }

  streamRead(ss, numItems, version);
  for (int i = 0; i < numItems; ++i) {
    streamRead(ss, tmpT, version);
    sgroup.addBondWithIdx(tmpT);
  }

  streamRead(ss, numItems, version);
  for (int i = 0; i < numItems; ++i) {
    SubstanceGroup::Bracket bracket;
    for (auto j = 0; j < 3; ++j) {
      streamRead(ss, tmpFloat, version);
      auto x = static_cast<double>(tmpFloat);
      streamRead(ss, tmpFloat, version);
      auto y = static_cast<double>(tmpFloat);
      streamRead(ss, tmpFloat, version);
      auto z = static_cast<double>(tmpFloat);
      bracket[j] = RDGeom::Point3D(x, y, z);
    }
    sgroup.addBracket(bracket);
  }

  streamRead(ss, numItems, version);
  for (int i = 0; i < numItems; ++i) {
    streamRead(ss, tmpT, version);
    RDGeom::Point3D vector;

    if ("SUP" == sgroup.getProp<std::string>("TYPE")) {
      streamRead(ss, tmpFloat, version);
      vector.x = static_cast<double>(tmpFloat);
      streamRead(ss, tmpFloat, version);
      vector.y = static_cast<double>(tmpFloat);
      streamRead(ss, tmpFloat, version);
      vector.z = static_cast<double>(tmpFloat);
    }

    sgroup.addCState(tmpT, vector);
  }

  streamRead(ss, numItems, version);
  for (int i = 0; i < numItems; ++i) {
    streamRead(ss, tmpT, version);
    unsigned int aIdx = tmpT;

    streamRead(ss, tmpInt, version);
    int lvIdx = tmpInt;

    std::string id;
    streamRead(ss, id, version);

    sgroup.addAttachPoint(aIdx, lvIdx, id);
  }

  return sgroup;
}

template <typename T>
void MolPickler::_pickleStereo(std::ostream &ss,
                               std::vector<StereoGroup> groups,
                               std::map<int, int> &atomIdxMap,
                               std::map<int, int> &bondIdxMap) {
  T tmpT = static_cast<T>(groups.size());
  streamWrite(ss, tmpT);
  assignStereoGroupIds(groups);
  for (auto &&group : groups) {
    streamWrite(ss, static_cast<T>(group.getGroupType()));
    if (group.getGroupType() != StereoGroupType::STEREO_ABSOLUTE) {
      streamWrite(ss, static_cast<T>(group.getWriteId()));
    }
    auto &atoms = group.getAtoms();
    streamWrite(ss, static_cast<T>(atoms.size()));
    for (auto &&atom : atoms) {
      tmpT = static_cast<T>(atomIdxMap[atom->getIdx()]);
      streamWrite(ss, tmpT);
    }

    auto &bonds = group.getBonds();
    streamWrite(ss, static_cast<T>(bonds.size()));
    for (auto &&bond : bonds) {
      tmpT = static_cast<T>(bondIdxMap[bond->getIdx()]);
      streamWrite(ss, tmpT);
    }
  }
}

template <typename T>
void MolPickler::_depickleStereo(std::istream &ss, ROMol *mol, int version) {
  T tmpT;
  streamRead(ss, tmpT, version);
  const auto numGroups = static_cast<unsigned>(tmpT);

  if (numGroups > 0u) {
    std::vector<StereoGroup> groups;
    for (unsigned group = 0u; group < numGroups; ++group) {
      T tmpT;
      streamRead(ss, tmpT, version);
      const auto groupType = static_cast<RDKit::StereoGroupType>(tmpT);

      unsigned gId = 0;
      if (version >= 14010 && groupType != StereoGroupType::STEREO_ABSOLUTE) {
        streamRead(ss, tmpT, version);
        gId = static_cast<unsigned>(tmpT);
      }

      streamRead(ss, tmpT, version);
      const auto numAtoms = static_cast<unsigned>(tmpT);

      std::vector<Atom *> atoms;
      atoms.reserve(numAtoms);
      for (unsigned i = 0u; i < numAtoms; ++i) {
        streamRead(ss, tmpT, version);
        atoms.push_back(mol->getAtomWithIdx(tmpT));
      }

      std::vector<Bond *> bonds;
      if (version > 16000) {
        streamRead(ss, tmpT, version);
        const auto numBonds = static_cast<unsigned>(tmpT);

        bonds.reserve(numBonds);
        for (unsigned i = 0u; i < numBonds; ++i) {
          streamRead(ss, tmpT, version);
          bonds.push_back(mol->getBondWithIdx(tmpT));
        }
      }
      groups.emplace_back(groupType, std::move(atoms), std::move(bonds), gId);
    }

    mol->setStereoGroups(std::move(groups));
  }
}

//! unpickle standard properties
void MolPickler::_unpickleProperties(std::istream &ss, RDProps &props,
                                     int version) {
  if (version >= 14000) {
    streamReadProps<std::uint16_t>(ss, props,
                                   MolPickler::getCustomPropHandlers());
  } else {
    streamReadProps<unsigned int>(ss, props,
                                  MolPickler::getCustomPropHandlers());
  }
}

//--------------------------------------
//
//            Version 1 Pickler:
//
//  NOTE: this is not 64bit clean, but it shouldn't be used anymore anyway
//
//--------------------------------------

void MolPickler::_pickleV1(const ROMol *mol, std::ostream &ss) {
  PRECONDITION(mol, "empty molecule");
  ROMol::ConstAtomIterator atIt;
  const Conformer *conf = nullptr;
  if (mol->getNumConformers() > 0) {
    conf = &(mol->getConformer());
  }
  for (atIt = mol->beginAtoms(); atIt != mol->endAtoms(); ++atIt) {
    const Atom *atom = *atIt;

    streamWrite(ss, BEGINATOM);
    streamWrite(ss, ATOM_NUMBER, atom->getAtomicNum());

    streamWrite(ss, ATOM_INDEX, atom->getIdx());

    streamWrite(ss, ATOM_POS);
    RDGeom::Point3D p;
    if (conf) {
      p = conf->getAtomPos(atom->getIdx());
    }
    streamWrite(ss, p.x);
    streamWrite(ss, p.y);
    streamWrite(ss, p.z);

    if (atom->getFormalCharge() != 0) {
      streamWrite(ss, ATOM_CHARGE, atom->getFormalCharge());
    }
    if (atom->getNumExplicitHs() != 0) {
      streamWrite(ss, ATOM_NEXPLICIT, atom->getNumExplicitHs());
    }
    if (atom->getChiralTag() != 0) {
      streamWrite(ss, ATOM_CHIRALTAG, atom->getChiralTag());
    }
    if (atom->getIsAromatic()) {
      streamWrite(ss, ATOM_ISAROMATIC,
                  static_cast<char>(atom->getIsAromatic()));
    }
    streamWrite(ss, ENDATOM);
  }

  ROMol::ConstBondIterator bondIt;
  for (bondIt = mol->beginBonds(); bondIt != mol->endBonds(); ++bondIt) {
    const Bond *bond = *bondIt;
    streamWrite(ss, BEGINBOND);
    streamWrite(ss, BOND_INDEX, bond->getIdx());
    streamWrite(ss, BOND_BEGATOMIDX, bond->getBeginAtomIdx());
    streamWrite(ss, BOND_ENDATOMIDX, bond->getEndAtomIdx());
    streamWrite(ss, BOND_TYPE, bond->getBondType());
    if (bond->getBondDir()) {
      streamWrite(ss, BOND_DIR, bond->getBondDir());
    }
    streamWrite(ss, ENDBOND);
  }
  streamWrite(ss, ENDMOL);
}

void MolPickler::_depickleV1(std::istream &ss, ROMol *mol) {
  PRECONDITION(mol, "empty molecule");
  Tags tag;

  auto *conf = new Conformer();
  mol->addConformer(conf);
  streamRead(ss, tag, 1);
  while (tag != ENDMOL) {
    switch (tag) {
      case BEGINATOM:
        _addAtomFromPickleV1(ss, mol);
        break;
      case BEGINBOND:
        _addBondFromPickleV1(ss, mol);
        break;
      default:
        UNDER_CONSTRUCTION("bad tag in pickle");
    }
    streamRead(ss, tag, 1);
  }
  mol->clearAllAtomBookmarks();
  mol->clearAllBondBookmarks();
}

void MolPickler::_addAtomFromPickleV1(std::istream &ss, ROMol *mol) {
  PRECONDITION(mol, "empty molecule");
  Tags tag;
  int intVar;
  double dblVar;
  char charVar;
  int version = 1;
  streamRead(ss, tag, version);
  auto *atom = new Atom();
  Conformer &conf = mol->getConformer();
  RDGeom::Point3D pos;
  while (tag != ENDATOM) {
    switch (tag) {
      case ATOM_INDEX:
        streamRead(ss, intVar, version);
        mol->setAtomBookmark(atom, intVar);
        break;
      case ATOM_NUMBER:
        streamRead(ss, intVar, version);
        atom->setAtomicNum(intVar);
        break;
      case ATOM_POS:
        streamRead(ss, pos.x, version);
        streamRead(ss, pos.y, version);
        streamRead(ss, pos.z, version);
        break;
      case ATOM_CHARGE:
        streamRead(ss, intVar, version);
        atom->setFormalCharge(intVar);
        break;
      case ATOM_NEXPLICIT:
        streamRead(ss, intVar, version);
        atom->setNumExplicitHs(intVar);
        break;
      case ATOM_CHIRALTAG:
        streamRead(ss, intVar, version);
        atom->setChiralTag(static_cast<Atom::ChiralType>(intVar));
        break;
      case ATOM_MASS:
        streamRead(ss, dblVar, version);
        // we don't need to set this anymore, but we do need to read it in
        // order to maintain backwards compatibility
        break;
      case ATOM_ISAROMATIC:
        streamRead(ss, charVar, version);
        atom->setIsAromatic(charVar);
        break;
      default:
        ASSERT_INVARIANT(0, "bad tag in atom block of pickle");
    }
    streamRead(ss, tag, version);
  }
  unsigned int id = mol->addAtom(atom, false, true);
  conf.setAtomPos(id, pos);
}
void MolPickler::_addBondFromPickleV1(std::istream &ss, ROMol *mol) {
  PRECONDITION(mol, "empty molecule");
  Tags tag;
  int intVar, idx = -1;
  int version = 1;
  Bond::BondType bt;
  Bond::BondDir bd;
  streamRead(ss, tag, version);
  auto *bond = new Bond();
  while (tag != ENDBOND) {
    switch (tag) {
      case BOND_INDEX:
        streamRead(ss, idx, version);
        break;
      case BOND_BEGATOMIDX:
        streamRead(ss, intVar, version);
        bond->setBeginAtomIdx(mol->getAtomWithBookmark(intVar)->getIdx());
        break;
      case BOND_ENDATOMIDX:
        streamRead(ss, intVar, version);
        bond->setEndAtomIdx(mol->getAtomWithBookmark(intVar)->getIdx());
        break;
      case BOND_TYPE:
        streamRead(ss, bt, version);
        bond->setBondType(bt);
        break;
      case BOND_DIR:
        streamRead(ss, bd, version);
        bond->setBondDir(bd);
        break;
      default:
        ASSERT_INVARIANT(0, "bad tag in bond block of pickle");
    }
    streamRead(ss, tag, version);
  }
  mol->addBond(bond, true);
}
};  // namespace RDKit
