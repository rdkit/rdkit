/* MaxMinPicker.cpp */
#include <cstring>
#include <cstdio>
#include <sys/time.h>

#include <vector>
#include <string>

#include <GraphMol/GraphMol.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/Fingerprints/Fingerprints.h>
#include <GraphMol/Fingerprints/MorganFingerprints.h>
#include <DataStructs/BitOps.h>

#include "MaxMinPicker.h"
#include "LeaderPicker.h"

std::vector<RDKit::ROMol *> mols;
std::vector<ExplicitBitVect *> fps;
// ExplicitBitVect **fps;

const char *poolname;
const char *pickname;
unsigned int origsize;
unsigned int picksize;
unsigned int poolsize;
unsigned int newpicks;
int seed;

unsigned long ticks = 0;

static unsigned int LoadDatabase(FILE *fp) {
  char buffer[32768];
  unsigned int result = 0;

  while (fgets(buffer, sizeof(buffer), fp)) {
    if (buffer[0] == '#' || buffer[0] == ' ' || buffer[0] == '\t') continue;
    char *ptr = buffer;
    while (*ptr && *ptr != ' ' && *ptr != '\t') ptr++;
    if (*ptr) {
      *ptr++ = '\0';
      char *end = ptr;
      while (*end && *end != '\n' && *end != '\r') end++;
      *end = '\0';
    }

    RDKit::RWMol *mol;
    try {
      mol = RDKit::SmilesToMol(buffer);
    } catch (...) {
      mol = (RDKit::RWMol *)0;
    }

    if (mol) {
      ExplicitBitVect *fp;
      fp = RDKit::MorganFingerprints::getFingerprintAsBitVect(*mol, 2, 2048);
      fps.push_back(fp);
      delete mol;
      result++;
    }
  }
  fprintf(stderr, "%u molecules\n", result);
  return result;
}

static unsigned int LoadDatabase(const char *fname) {
  if (fname && strcmp(fname, "-")) {
    FILE *fp = fopen(fname, "rb");
    if (!fp) return 0;
    unsigned int result = LoadDatabase(fp);
    fclose(fp);
    return result;
  } else
    return LoadDatabase(stdin);
}

#ifdef UNUSED
static void GenerateFingerprints() {
  unsigned int count = (unsigned int)mols.size();
  unsigned int size = (unsigned int)(count * sizeof(void *));
  fps = (ExplicitBitVect **)malloc(size);

  for (unsigned int i = 0; i < count; i++)
    fps[i] =
        RDKit::MorganFingerprints::getFingerprintAsBitVect(*mols[i], 2, 2048);
  fprintf(stderr, "%u Fingerprints\n", count);
}

static void DestroyFingerprints() {
  unsigned int count = (unsigned int)mols.size();
  for (unsigned int i = 0; i < count; i++) delete fps[i];
  free(fps);
}
#endif

static void DisplayUsage() {
  fputs("usage:  MaxMinPicker [options] <poolfile> [<pickfile>]\n", stderr);
  exit(1);
}

static void ProcessCommandLine(int argc, char *argv[]) {
  poolname = (const char *)0;
  pickname = (const char *)0;
  newpicks = 1000;

  for (int i = 1; i < argc; i++) {
    const char *ptr = argv[i];
    if (ptr[0] == '-' && ptr[1]) {
      if (ptr[1] >= '0' && ptr[1] <= '9') {
        newpicks = atoi(ptr + 1);
      } else
        DisplayUsage();
    } else if (!poolname) {
      poolname = ptr;
    } else if (!pickname) {
      pickname = ptr;
    } else
      DisplayUsage();
  }

  if (!poolname) DisplayUsage();
}

double MyDist(int i, int j) {
  ticks++;
  return 1.0 - TanimotoSimilarity(*fps[i], *fps[j]);
}

int main(int argc, char *argv[]) {
  struct timeval beg, end;
  double elapsed;

  ProcessCommandLine(argc, argv);

  gettimeofday(&beg, (struct timezone *)0);
  origsize = LoadDatabase(poolname);
  gettimeofday(&end, (struct timezone *)0);
  elapsed = (end.tv_sec + 0.000001 * end.tv_usec) -
            (beg.tv_sec + 0.000001 * beg.tv_usec);
  fprintf(stderr, "Elapsed time: %g secs\n", elapsed);

  if (pickname) {
    gettimeofday(&beg, (struct timezone *)0);
    picksize = LoadDatabase(pickname);
    gettimeofday(&end, (struct timezone *)0);
    elapsed = (end.tv_sec + 0.000001 * end.tv_usec) -
              (beg.tv_sec + 0.000001 * beg.tv_usec);
    fprintf(stderr, "Elapsed time: %g secs\n", elapsed);
  } else
    picksize = 0;
  poolsize = origsize + picksize;
  if (newpicks == 0) newpicks = origsize;

#ifdef UNUSED
  gettimeofday(&beg, (struct timezone *)0);
  GenerateFingerprints();
  gettimeofday(&end, (struct timezone *)0);
  elapsed = (end.tv_sec + 0.000001 * end.tv_usec) -
            (beg.tv_sec + 0.000001 * beg.tv_usec);
  fprintf(stderr, "Elapsed time: %g secs\n", elapsed);
#endif

  RDKit::INT_VECT firstPicks;
  if (picksize) {
    firstPicks.reserve(picksize);
    for (unsigned int i = 0; i < picksize; i++)
      firstPicks.push_back((int)(i + origsize));
  }

  //  RDPickers::MaxMinPicker mmpicker;
  double threshold = 0.9;  // 0.675;
  RDPickers::LeaderPicker ldpicker(threshold, 4);

  gettimeofday(&beg, (struct timezone *)0);
  RDKit::INT_VECT iv = ldpicker.lazyPick(MyDist, poolsize, picksize + newpicks,
                                         firstPicks, threshold, 16);
  gettimeofday(&end, (struct timezone *)0);
  elapsed = (end.tv_sec + 0.000001 * end.tv_usec) -
            (beg.tv_sec + 0.000001 * beg.tv_usec);
  fprintf(stderr, "Elapsed time: %g secs\n", elapsed);

  unsigned int count = (unsigned int)iv.size();
  for (unsigned int i = 0; i < std::min(iv.size(), (size_t)10); ++i) {
    for (unsigned int j = 0; j < i; ++j) {
      std::cerr << iv[i] << " " << iv[j] << ": " << MyDist(iv[i], iv[j])
                << std::endl;
    }
  }
  fprintf(stderr, "%u picks\n", count);
  fprintf(stderr, "%g threshold\n", threshold);
  fprintf(stderr, "%lu comparisons\n", ticks);

  // DestroyFingerprints();
  return 0;
}
