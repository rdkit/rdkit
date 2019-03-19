#include <GraphMol/RDKitBase.h>
#include <GraphMol/RDKitQueries.h>
#include <GraphMol/SubstructLibrary/EditableSubstructLibraryTrustedSmilesWithPattern.h>

#include <GraphMol/Substruct/SubstructMatch.h>

#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/FileParsers/MolSupplier.h>

using namespace RDKit;

namespace {

template <typename K, typename M, typename F>
void checkMatches(EditableSubstructLibrary<K, M, F> &sssLib,
                  const ROMol &pattern, const std::vector<K> idMatches,
                  const unsigned int maxHits) {
  boost::dynamic_bitset<> hasMatch(sssLib.size());
  BOOST_FOREACH (const auto &id, idMatches) {
    unsigned int idx = sssLib.idToIndex(id);
    hasMatch[idx] = 1;
  }
  // std::cerr << "sss search got " << idMatches.size() << " hits " << std::endl;

  unsigned int noHits(0);
  for (auto i = 0u; i < sssLib.size(); i++) {
    auto id = sssLib.indexToId(i);
    auto mol = sssLib.getMol(id);
    MatchVectType match;
    bool matched = SubstructMatch(*mol, pattern, match);
    if (matched) noHits++;
    // If we have a limit on the numbrer of hits we can't validate misses
    if (noHits == 0 || hasMatch[i]) TEST_ASSERT(hasMatch[i] == matched);
    if (maxHits > 0 && noHits == maxHits) break;
  }
}

template <typename K, typename M, typename F>
void runTest(EditableSubstructLibrary<K, M, F> &sssLib, const ROMol &pattern,
             const int nThreads, const unsigned int nExpected) {
  auto idMatches = sssLib.getMatches(pattern, true, true, false, nThreads);
  TEST_ASSERT(idMatches.size() == nExpected);
  checkMatches(sssLib, pattern, idMatches, 0);
}

template <typename K, typename M, typename F>
void runChunkedTest(EditableSubstructLibrary<K, M, F> &sssLib,
                    const ROMol &pattern, int nThreads,
                    const unsigned int nExpected,
                    const unsigned int maxHits = 0) {
  auto hitlist = sssLib.getHitlistMatches(pattern, true, true, false, 10,
                                          nThreads, maxHits);
  std::vector<K> idMatches;
  while (true) {
    auto chunk = hitlist->next();
    if (chunk.empty()) break;
    //std::cerr << "received chunk of size " << chunk.size() << std::endl;
    idMatches.insert(idMatches.end(), chunk.cbegin(), chunk.cend());
  }
  TEST_ASSERT(idMatches.size() == nExpected);
  checkMatches(sssLib, pattern, idMatches, maxHits);
}

void loadFromSdf(EditableSubstructLibraryTrustedSmilesWithPattern &sssLib) {
  std::string fName = getenv("RDBASE");
  fName += "/Data/NCI/first_200.names.sdf";
  SDMolSupplier suppl(fName);
  while (!suppl.atEnd()) {
    ROMol *mol = nullptr;
    try {
      mol = suppl.next();
    } catch (...) {
      continue;
    }
    if (!mol) continue;
    sssLib.addMol(*mol);
    delete mol;
  }
}

void loadFromSmiles(EditableSubstructLibraryTrustedSmilesWithPattern &sssLib) {
  std::string fName = getenv("RDBASE");
  fName += "/Data/NCI/first_200.names.smi";
  std::ifstream in(fName);
  std::string smi, id;
  std::vector<std::string> smiles, fingerprints;
  while (in >> smi >> id) {
    smiles.push_back(smi + " " + id);
    auto fp = sssLib.makeStringFingerprint(smi);
    fingerprints.push_back(fp);
  }

  sssLib.addSmiles(smiles, fingerprints);
}

}  // namespace

void testSimple() {
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "    Simple tests" << std::endl;

  EditableSubstructLibraryTrustedSmilesWithPattern sssLib{};
  loadFromSdf(sssLib);

  {
    ROMol *query = SmartsToMol("[#6;$([#6]([#6])[!#6])]");
    runTest(sssLib, *query, 1, 185);
#ifdef RDK_TEST_MULTITHREADED
    runTest(sssLib, *query, -1, 185);
#endif

    std::vector<std::string> ids{"10", "20", "30", "40"};
    auto startSize = sssLib.size();
    auto smi200 = MolToSmiles(*sssLib.getMol("202"), true);
    auto mol40 = sssLib.getMol("40");
    sssLib.removeMols(ids);
    auto currentSize = sssLib.size();
    TEST_ASSERT(currentSize == startSize - 4);
    auto currentSmi200 = MolToSmiles(*sssLib.getMol("202"), true);
    TEST_ASSERT(currentSmi200.compare(smi200) == 0);

    runTest(sssLib, *query, 1, 181);
#ifdef RDK_TEST_MULTITHREADED
    runTest(sssLib, *query, -1, 181);
#endif

    sssLib.addMol("40", *mol40);

    runTest(sssLib, *query, 1, 182);
#ifdef RDK_TEST_MULTITHREADED
    runTest(sssLib, *query, -1, 182);
#endif

    delete query;
  }
}

void testHitlistOnSdfInput() {
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "    Hitlist search on library built from SDF"
                        << std::endl;

  EditableSubstructLibraryTrustedSmilesWithPattern sssLib{};
  loadFromSdf(sssLib);
  {
    ROMol *query = SmartsToMol("[#6;$([#6]([#6])[!#6])]");
    runChunkedTest(sssLib, *query, 1, 185);
#ifdef RDK_TEST_MULTITHREADED
    runChunkedTest(sssLib, *query, -1, 185);
#endif
  }
}

void testHitlistOnSmilesInput() {
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "    Hitlist search on library built from smiles"
                        << std::endl;

  EditableSubstructLibraryTrustedSmilesWithPattern sssLib{};
  loadFromSmiles(sssLib);
  {
    ROMol *query = SmartsToMol("[#6;$([#6]([#6])[!#6])]");
    runChunkedTest(sssLib, *query, 1, 185);
#ifdef RDK_TEST_MULTITHREADED
    runChunkedTest(sssLib, *query, -1, 185);
#endif
  }
}

void testMaxHits() {
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "    Hitlist search with max hits "
                        << std::endl;

  EditableSubstructLibraryTrustedSmilesWithPattern sssLib{};
  loadFromSmiles(sssLib);
  {
    ROMol *query = SmartsToMol("[#6;$([#6]([#6])[!#6])]");
    runChunkedTest(sssLib, *query, 1, 7, 7);
#ifdef RDK_TEST_MULTITHREADED
    runChunkedTest(sssLib, *query, -1, 7, 7);
#endif
  }
}

int main() {
  RDLog::InitLogs();
  boost::logging::enable_logs("rdApp.*");
  testSimple();
  testHitlistOnSdfInput();
  testHitlistOnSmilesInput();
  testMaxHits();
  return 0;
}