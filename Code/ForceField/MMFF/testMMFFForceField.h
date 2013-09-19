class BondStretchInstance {
  public:
    unsigned int idx;
    unsigned int iAtomType;
    unsigned int jAtomType;
    unsigned int ffType;
    double kb;
};

class AngleBendInstance {
  public:
    unsigned int idx;
    unsigned int iAtomType;
    unsigned int jAtomType;
    unsigned int kAtomType;
    unsigned int ffType;
    double ka;
};

class StretchBendInstance {
  public:
    unsigned int idx;
    unsigned int iAtomType;
    unsigned int jAtomType;
    unsigned int kAtomType;
    unsigned int ffType;
    double kba;
};

class OopBendInstance {
  public:
    unsigned int idx;
    unsigned int iAtomType;
    unsigned int jAtomType;
    unsigned int kAtomType;
    unsigned int lAtomType;
    double koop;
};

class TorsionInstance {
  public:
    unsigned int idx;
    unsigned int iAtomType;
    unsigned int jAtomType;
    unsigned int kAtomType;
    unsigned int lAtomType;
    unsigned int ffType;
    double V1;
    double V2;
    double V3;
};

bool fexist(std::string filename);
bool fgrep(std::fstream &fStream, std::string key);
bool fgrep(std::fstream &fStream, std::string key, std::string &line);
#ifndef RDK_NOGZIP
void fgunzip(std::string filename);
#endif
void skipLines(std::istream& stream, unsigned int nLines);
bool sortAngleBendInstanceVec(AngleBendInstance *a, AngleBendInstance *b);
bool sortBondStretchInstanceVec(BondStretchInstance *a, BondStretchInstance *b);
bool sortOopBendInstanceVec(OopBendInstance *a, OopBendInstance *b);
bool sortStretchBendInstanceVec(StretchBendInstance *a, StretchBendInstance *b);
bool sortTorsionInstanceVec(TorsionInstance *a, TorsionInstance *b);
void fixAngleBendInstance(AngleBendInstance *angleBendInstance);
void fixBondStretchInstance(BondStretchInstance *bondStretchInstance);
void fixOopBendInstance(OopBendInstance *oopBendInstance);
void fixTorsionInstance(TorsionInstance *torsionInstance);
