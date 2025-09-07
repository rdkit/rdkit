/*
 * File: jni_JNIMatch.cpp
 *
 * Implements the wrappers for the JNIMatch class.
 */

#include "jni_JNIMatch.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#ifdef __WIN32__
#include <windows.h>
#endif

#ifdef __cplusplus
extern "C" {
#endif

// FORTIFY #define NOFORTIFY
#include "local.h"
#include "reaccs.h"
#include "reaccsio.h"
#include "utilities.h"
#include "perceive.h"
#include "smi2mol.h"
#include "ssmatch.h"
#include "patclean.h"
#include "canonizer.h"

/* Referred to in local.h but not used */
FILE *log_file = (FILE *)NULL;

/* Values between 0 and 999 enable heap debugging completely */
/* Values above 1000 enable the segment that has the same thousends digit */
#define DEBUG_HEAP 0000

/**
 * Returns a full path name for a file fname in the directory pointed
 * to by the environment variables TMP or TEMP in that order.
 */
static char *GetTempName(const char *fname)
{
   static char buffer[1000];
   char *varp;
   varp = getenv("TMP");
   if (!varp) varp = getenv("TEMP");
   if (varp)
   {
      strcpy(buffer, varp);
      strcat(buffer, "\\");
   }
   else      strcpy(buffer, "C:\\TEMP\\");
   strcat(buffer, fname);
   return (buffer);
}

void logString(char *str)
{
   FILE *fp;
   fp = fopen(GetTempName("tmp.log"), "a");
   fprintf(fp, "%s\n", str);
   fclose(fp);
}

// static int checkHeap(int where)
// {
//    int result = TRUE;
// #ifdef __WIN32__
//    HGLOBAL heap;
//    if (DEBUG_HEAP == 0) return TRUE;
//    if (DEBUG_HEAP/1000 != 0  &&
//        where / 1000 != DEBUG_HEAP/ 1000) return TRUE;
// 
//    heap = GetProcessHeap();
//    result = HeapValidate(heap,0,NULL);
//    if (!result)
//       fprintf(stderr, "HeapValidate(%p, %d, %d) at %d returns %d\n",
//          heap, 0, NULL, where, result);
// #else
//    if (DEBUG_HEAP == 0) return TRUE;
//    if (DEBUG_HEAP/1000 != 0  &&
//        where / 1000 != DEBUG_HEAP/ 1000) return TRUE;
// 
// #endif
//    return (result);
// }

/*
 * Class:     jni_JNIMatch
 * Method:    heapCheck
 * Signature: ()Z
 */
JNIEXPORT jboolean
JNICALL Java_jni_JNIMatch_heapCheck (JNIEnv *env, jobject thisObj)
{
   // if (checkHeap(0)) // disabled
   if (TRUE)
      return JNI_TRUE;
   else
      return JNI_FALSE;
}

/*
 * Class:     jni_JNIMatch
 * Method:    querySmilesToHandleNative
 * Signature: (Ljava/lang/String;)I
 */
JNIEXPORT jlong JNICALL Java_jni_JNIMatch_querySmilesToHandleNative
  (JNIEnv *env, jobject thisObject, jstring querySmiles)
{
   struct reaccs_molecule_t *qp;
   const char *const_query;
   jboolean isCopy = JNI_FALSE;
   neighbourhood_t *nbp;

// checkHeap(2002);
//
   // Read query
   const_query = env->GetStringUTFChars(querySmiles, &isCopy);

   qp = SMIToMOL(const_query, DO_LAYOUT);

   // if (isCopy == JNI_TRUE)
   env->ReleaseStringUTFChars(querySmiles, const_query);

// checkHeap(2012);
//
   if (!qp) return ((jlong)qp);

   /* Compute query_H_count to match the required explicit hydrogens */
   MakeHydrogensImplicit(qp);
   nbp     = TypeAlloc(qp->n_atoms, neighbourhood_t);
   SetupNeighbourhood(qp,nbp,qp->n_atoms);
   PerceiveMarkush(qp, nbp);
   SetRingSizeFlags(qp, 14, nbp);
   MyFree((char *)nbp);
   qp = ThreadQueryAtoms(qp);

// checkHeap(2032);
//
   return ((jlong)qp);
}


/*
 * Class:     jni_JNIMatch
 * Method:    queryCTToHandleNative
 * Signature: (Ljava/lang/String;)I
 */
JNIEXPORT jlong JNICALL Java_jni_JNIMatch_queryCTToHandleNative
  (JNIEnv *env, jobject thisObject, jstring queryCT)
{
   struct reaccs_molecule_t *qp;
   const char *const_query;
   jboolean isCopy = JNI_FALSE;
   Fortran_FILE *ffp;
   neighbourhood_t *nbp;

// checkHeap(1002);
//
   // Read query
   const_query = env->GetStringUTFChars(queryCT, &isCopy);

   ffp = FortranStringOpen((char *)const_query);
   qp = TypeAlloc(1,struct reaccs_molecule_t);
   if (FORTRAN_NORMAL != ReadREACCSMolecule(ffp,qp,""))
   {
      fprintf(stderr, "failed to read template\n");
      FreeMolecule(qp);
      qp = (struct reaccs_molecule_t *)NULL;
   }
   else if (qp->n_atoms <= 0)
   {
      fprintf(stderr, "template without atoms\n");
      FreeMolecule(qp);
      qp = (struct reaccs_molecule_t *)NULL;
   }
   FortranClose(ffp);

// if (isCopy == JNI_TRUE)
   env->ReleaseStringUTFChars(queryCT, const_query);

// checkHeap(1012);
//
   if (qp == (struct reaccs_molecule_t *)NULL) return 0;
   /* Compute query_H_count to match the required explicit hydrogens */
   MakeHydrogensImplicit(qp);
   nbp     = TypeAlloc(qp->n_atoms, neighbourhood_t);
   SetupNeighbourhood(qp,nbp,qp->n_atoms);
   PerceiveMarkush(qp, nbp);
   SetRingSizeFlags(qp, 14, nbp);
   MyFree((char *)nbp);
   qp = ThreadQueryAtoms(qp);

// checkHeap(1032);
//
   return ((jlong)qp);
}

/*
 * Class:     jni_JNIMatch
 * Method:    disposeCTHandleNative
 * Signature: (I)V
 */
JNIEXPORT void JNICALL Java_jni_JNIMatch_disposeCTHandleNative
  (JNIEnv *env, jobject thisObject, jlong ctHandle)
{
   struct reaccs_molecule_t *mp;

   mp = (struct reaccs_molecule_t *)ctHandle;
   if (mp) FreeMolecule(mp);
}

static struct reaccs_molecule_t *SMIToMOLNoStereo(const char *const_smiles)
{
   char *smiles, *cp;
   struct reaccs_molecule_t *mp;

   if (const_smiles)
   {
      smiles = TypeAlloc(strlen(const_smiles)+1, char);
      cp = smiles;
      do
      {
         if ((*const_smiles) != '@') (*cp++) = (*const_smiles);
         const_smiles++;
      }
      while (*const_smiles);
      mp = SMIToMOL(smiles, DROP_TRIVIAL_HYDROGENS);
      MyFree(smiles);
   }
   else
      mp = (struct reaccs_molecule_t *)NULL;
   return (mp);
}

/*
 * Class:     jni_JNIMatch
 * Method:    smilesMatchesHandleNative
 * Signature: (Ljava/lang/String;I)Z
 */
JNIEXPORT jboolean JNICALL Java_jni_JNIMatch_smilesMatchesHandleNative
  (JNIEnv *env, jobject thisObject, jstring smilesObj, jlong queryHandle)
{
   struct reaccs_molecule_t *mp;
   struct reaccs_molecule_t *qp;
   const char *const_smiles;
   jboolean isCopy = JNI_FALSE;
   jboolean result = JNI_FALSE;
   ssmatch_t * matches;
   struct reaccs_bond_t *bp;
   unsigned int i;
   int *H_count;
   neighbourhood_t *nbp;

// checkHeap(1002);
//
   const_smiles = env->GetStringUTFChars(smilesObj, &isCopy);
   mp = SMIToMOLNoStereo(const_smiles);
// checkHeap(1003);
// checkHeap(1004);
// if (isCopy == JNI_TRUE)
       env->ReleaseStringUTFChars(smilesObj, const_smiles);
   /* Illegal molecules don't match at all */
   if (mp == (struct reaccs_molecule_t *)NULL)
      return (result);
   nbp     = TypeAlloc(mp->n_atoms, neighbourhood_t);
   SetupNeighbourhood(mp,nbp,mp->n_atoms);
   SetRingSizeFlags(mp, 14, nbp);
   MyFree((char *)nbp);

   qp = (struct reaccs_molecule_t *)queryHandle;

// checkHeap(1012);
//
   if (qp == (struct reaccs_molecule_t *)NULL)  /* illegal query => no match */
   {
      FreeMolecule(mp);
      return (result);
   }

   /* Set up hydrogen count fields in structure for matching */
   H_count = TypeAlloc(mp->n_atoms+1, int);
   ComputeImplicitH(mp, H_count);
   /* Add the explicit hydrogens to the implicit counts */
   for (i=0, bp=mp->bond_array; i<mp->n_bonds; i++, bp++)
   {
      if (0 == strcmp("H", mp->atom_array[bp->atoms[0]-1].atom_symbol))
         H_count[bp->atoms[1]]++;
      else if (0 == strcmp("H", mp->atom_array[bp->atoms[1]-1].atom_symbol))
         H_count[bp->atoms[0]]++;
   }
   /* set the 'query_H_count' field to the correct value */
   for (i=0; i<mp->n_atoms; i++)
      if (H_count[i+1] >= 0) 
         mp->atom_array[i].query_H_count = ZERO_COUNT+H_count[i+1];

// PrintREACCSMolecule(stderr,mp,"STRUCTURE");
// PrintREACCSMolecule(stderr,qp,"QUERY");
// 
// checkHeap(1022);
//
   matches = SSMatch(mp, qp, TRUE);
   if (matches != (ssmatch_t *)NULL)
   {
      result = JNI_TRUE;
      FreeSSMatch(matches);
      FreeSSMatchHeap();
   }

   FreeMolecule(mp);

// checkHeap(1032);
//
   MyFree((char *)H_count);
   return (result);
}

/*
 * Class:     jni_JNIMatch
 * Method:    computeMatchColorsForCT
 * Signature: (Ljava/lang/String;Ljava/lang/String;[I)I
 */
JNIEXPORT jboolean
JNICALL Java_jni_JNIMatch_computeMatchColorsForCT(JNIEnv *env,
                                                  jobject thisObj,
                                                  jstring molCT,
                                                  jstring queryCT,
                                                  jintArray color_buffer)
{
   struct reaccs_molecule_t *mp;
   struct reaccs_molecule_t *qp;
   const char *const_mol;
   const char *const_query;
   Fortran_FILE *ffp;
   jint *colors;
   jboolean isCopy = JNI_FALSE;
   jboolean result = JNI_FALSE;
   ssmatch_t * matches;
   struct reaccs_bond_t *bp;
   int size;
   unsigned int i;
   int *H_count;
   neighbourhood_t *nbp;

   // Make molecule from molCT
   const_mol = env->GetStringUTFChars(molCT, &isCopy);
   ffp = FortranStringOpen((char *)const_mol);
   mp = TypeAlloc(1,struct reaccs_molecule_t);
   if (FORTRAN_NORMAL != ReadREACCSMolecule(ffp,mp,""))
   {
      fprintf(stderr, "failed to read molecule\n");
      FreeMolecule(mp);
      mp = (struct reaccs_molecule_t *)NULL;
   }
   FortranClose(ffp);
   env->ReleaseStringUTFChars(molCT, const_mol);
   if (mp == NULL) return result;

   // Make molecule from queryCT
   const_query = env->GetStringUTFChars(queryCT, &isCopy);
   ffp = FortranStringOpen((char *)const_query);
   qp = TypeAlloc(1,struct reaccs_molecule_t);
   if (FORTRAN_NORMAL != ReadREACCSMolecule(ffp,qp,""))
   {
      fprintf(stderr, "failed to read template\n");
      FreeMolecule(qp);
      qp = (struct reaccs_molecule_t *)NULL;
   }
   else if (qp->n_atoms <= 0)
   {
      fprintf(stderr, "template without atoms\n");
      FreeMolecule(qp);
      qp = (struct reaccs_molecule_t *)NULL;
   }
   FortranClose(ffp);
   if (qp == (struct reaccs_molecule_t *)NULL)  /* illegal query => no match */
   {
      FreeMolecule(mp);
   }
   env->ReleaseStringUTFChars(queryCT, const_query);
   if (qp == NULL) return result;

   // perceive rings to support matching
   nbp     = TypeAlloc(mp->n_atoms, neighbourhood_t);
   SetupNeighbourhood(mp,nbp,mp->n_atoms);
   SetRingSizeFlags(mp, 14, nbp);
   MyFree((char *)nbp);

   /* Compute query_H_count to match the required explicit hydrogens */
   MakeHydrogensImplicit(qp);
   nbp     = TypeAlloc(qp->n_atoms, neighbourhood_t);
   SetupNeighbourhood(qp,nbp,qp->n_atoms);
   PerceiveMarkush(qp, nbp);
   SetRingSizeFlags(qp, 14, nbp);
   MyFree((char *)nbp);
   qp = ThreadQueryAtoms(qp);

   /* Set up hydrogen count fields in structure for matching */
   H_count = TypeAlloc(mp->n_atoms+1, int);
   ComputeImplicitH(mp, H_count);
   /* Add the explicit hydrogens to the implicit counts */
   for (i=0, bp=mp->bond_array; i<mp->n_bonds; i++, bp++)
   {
      if (0 == strcmp("H", mp->atom_array[bp->atoms[0]-1].atom_symbol))
         H_count[bp->atoms[1]]++;
      else if (0 == strcmp("H", mp->atom_array[bp->atoms[1]-1].atom_symbol))
         H_count[bp->atoms[0]]++;
   }
   /* set the 'query_H_count' field to the correct value */
   for (i=0; i<mp->n_atoms; i++)
      if (H_count[i+1] >= 0) 
         mp->atom_array[i].query_H_count = ZERO_COUNT+H_count[i+1];

   // only look for first match since single_match == TRUE
   matches = SSMatch(mp, qp, TRUE);
   if (matches != (ssmatch_t *)NULL)
   {
      size = env->GetArrayLength(color_buffer);
      isCopy = JNI_FALSE;
      colors = env->GetIntArrayElements(color_buffer, &isCopy);
      for (i=0; (int)i<matches->n_match; i++)
         if (size > matches->match_atoms[i])
            colors[matches->match_atoms[i]] = 1;
         else
            fprintf(stderr, "match atom %d exceeds size of colors array %d\n", matches->match_atoms[i], size);
      env->ReleaseIntArrayElements(color_buffer, colors, 0);
      result = JNI_TRUE;
      FreeSSMatch(matches);
      FreeSSMatchHeap();
   }

   FreeMolecule(mp); FreeMolecule(qp);
   MyFree((char *)H_count);

   return (result);
}

/*
 * Class:     JNIMatch
 * Method:    smilesMatchesQueryCTNative
 * Signature: (Ljava/lang/String;Ljava/lang/String;)Z
 */
JNIEXPORT jboolean
JNICALL Java_jni_JNIMatch_smilesMatchesQueryCTNative(JNIEnv *env,
                                                     jobject thisObj,
                                                     jstring smilesObj,
                                                     jstring queryCT)
{
   struct reaccs_molecule_t *mp;
   struct reaccs_molecule_t *qp;
   const char *const_smiles;
   const char *const_query;
   jboolean isCopy = JNI_FALSE;
   jboolean result = JNI_FALSE;
   Fortran_FILE *ffp;
   ssmatch_t * matches;
   struct reaccs_bond_t *bp;
   unsigned int i;
   int *H_count;
   neighbourhood_t *nbp;

// checkHeap(1002);
//
   const_smiles = env->GetStringUTFChars(smilesObj, &isCopy);
   mp = SMIToMOLNoStereo(const_smiles);
   // if (isCopy == JNI_TRUE)
   env->ReleaseStringUTFChars(smilesObj, const_smiles);
   /* Illegal molecules don't match at all */
   if (mp == (struct reaccs_molecule_t *)NULL)
      return (result);

   // Read query
   const_query = env->GetStringUTFChars(queryCT, &isCopy);

   ffp = FortranStringOpen((char *)const_query);
   qp = TypeAlloc(1,struct reaccs_molecule_t);
   if (FORTRAN_NORMAL != ReadREACCSMolecule(ffp,qp,""))
   {
      fprintf(stderr, "failed to read template\n");
      FreeMolecule(qp);
      qp = (struct reaccs_molecule_t *)NULL;
   }
   FortranClose(ffp);

   // if (isCopy == JNI_TRUE)
   env->ReleaseStringUTFChars(queryCT, const_query);

// checkHeap(1012);
//
   if (qp == (struct reaccs_molecule_t *)NULL)  /* illegal query => no match */
   {
      FreeMolecule(mp);
      return (result);
   }
   if (qp->n_atoms <= 0)
   {
      FreeMolecule(qp);
      FreeMolecule(mp);
      return (result);
   }
   nbp     = TypeAlloc(mp->n_atoms, neighbourhood_t);
   SetupNeighbourhood(mp,nbp,mp->n_atoms);
   SetRingSizeFlags(mp, 14, nbp);
   MyFree((char *)nbp);

   /* Compute query_H_count to match the required explicit hydrogens */
   MakeHydrogensImplicit(qp);
   nbp     = TypeAlloc(qp->n_atoms, neighbourhood_t);
   SetupNeighbourhood(qp,nbp,qp->n_atoms);
   PerceiveMarkush(qp, nbp);
   SetRingSizeFlags(qp, 14, nbp);
   MyFree((char *)nbp);
   qp = ThreadQueryAtoms(qp);

   /* Set up hydrogen count fields in structure for matching */
   H_count = TypeAlloc(mp->n_atoms+1, int);
   ComputeImplicitH(mp, H_count);
   /* Add the explicit hydrogens to the implicit counts */
   for (i=0, bp=mp->bond_array; i<mp->n_bonds; i++, bp++)
   {
      if (0 == strcmp("H", mp->atom_array[bp->atoms[0]-1].atom_symbol))
         H_count[bp->atoms[1]]++;
      else if (0 == strcmp("H", mp->atom_array[bp->atoms[1]-1].atom_symbol))
         H_count[bp->atoms[0]]++;
   }
   /* set the 'query_H_count' field to the correct value */
   for (i=0; i<mp->n_atoms; i++)
      if (H_count[i+1] >= 0) 
         mp->atom_array[i].query_H_count = ZERO_COUNT+H_count[i+1];


// PrintREACCSMolecule(stderr,mp,"STRUCTURE");
// PrintREACCSMolecule(stderr,qp,"QUERY");
// 
// checkHeap(1022);
//
   matches = SSMatch(mp, qp, TRUE);
   if (matches != (ssmatch_t *)NULL)
   {
      result = JNI_TRUE;
      FreeSSMatch(matches);
      FreeSSMatchHeap();
   }

   FreeMolecule(mp); FreeMolecule(qp);

// checkHeap(1032);
//
   MyFree((char *)H_count);
   return (result);
}

/*
 * Class:     jni_JNIMatch
 * Method:    canonicalSmilesNative
 * Signature: (Ljava/lang/String;[BI)I
 */
JNIEXPORT jint JNICALL Java_jni_JNIMatch_canonicalSmilesNative
  (JNIEnv *env, jobject thisObj,
   jstring inSmiles, jbyteArray outBuffer, jint flags)
{
   const char *const_smiles;
   jboolean isCopy = JNI_FALSE;
   int size;
   jbyte *buffer_bytes;
   char *outsmiles;
   int result;

   const_smiles = env->GetStringUTFChars(inSmiles, &isCopy);
   if (const_smiles == NULL)
   {
// if (isCopy == JNI_TRUE)
      env->ReleaseStringUTFChars(inSmiles, const_smiles);
      return -1;
   }
   outsmiles = CanSmiles((char *)const_smiles, flags);
// if (isCopy == JNI_TRUE)
   env->ReleaseStringUTFChars(inSmiles, const_smiles);
   if (outsmiles == NULL) return -1;
   result = strlen(outsmiles);

   size = env->GetArrayLength(outBuffer);
   /* lock the outBuffer and return a pointer to it */
   buffer_bytes = env->GetByteArrayElements(outBuffer, &isCopy);
   /* copy the canonical SMILES into buffer */
   if (size > (int)strlen(outsmiles))
      strcpy((char *)buffer_bytes, outsmiles);
   else
      result = -2;
   MyFree((char *)outsmiles);
// Release buffer and possibly copy bytes back
// if (isCopy == JNI_TRUE)
   env->ReleaseByteArrayElements(outBuffer, buffer_bytes, 0);

   return (result);
}

/*
 * Class:     jni_JNIMatch
 * Method:    computeFingerprintColorsForCT
 * Signature: (Ljava/lang/String;[BZI[I)I
 */
JNIEXPORT jint JNICALL Java_jni_JNIMatch_computeFingerprintColorsForCT
   (JNIEnv *env, jobject thisObj,
    jstring CTObj, jbyteArray buffer, jboolean as_query, jint which_bits,
    jintArray color_buffer)
{
    const char *const_CT;
    Fortran_FILE *ffp;
    struct reaccs_molecule_t *mp;
    int result;
    jboolean isCopy = JNI_FALSE;
    int size;
    int nbits;
    int nmasks;
    int nzero, ndiff;
    jbyte *fingerprint;
    int *mol_flags;
    int *flag_counts;
    int *mol_counts;
    int *tmp_counts;
    int *core_bits;
    int *count_changed;
    int *tmp_colors;
    jint *colors;
    unsigned int i;
    int j;

    // Read CT
    const_CT = env->GetStringUTFChars(CTObj, &isCopy);
    ffp = FortranStringOpen((char *)const_CT);
    mp = TypeAlloc(1,struct reaccs_molecule_t);
    if (FORTRAN_NORMAL != ReadREACCSMolecule(ffp,mp,""))
    {
       fprintf(stderr, "failed to read template\n");
       FreeMolecule(mp);
       mp = (struct reaccs_molecule_t *)NULL;
    }
    FortranClose(ffp);
// if (isCopy == JNI_TRUE)
    env->ReleaseStringUTFChars(CTObj, const_CT);
    if (mp == NULL) return (0);

    size = env->GetArrayLength(buffer);
    isCopy = JNI_FALSE;
    fingerprint = env->GetByteArrayElements(buffer, &isCopy);
    mol_flags = TypeAlloc(size*8, int);          // indicates which pattern flags contributed to this bit position
    mol_counts = TypeAlloc(size*8, int);         // accumulated FP counts for full molecule
    flag_counts = TypeAlloc(size*8, int);        // counts how many pattern flags contributed to this bit position
    tmp_counts = TypeAlloc(size*8, int);
    core_bits = TypeAlloc(size*8, int);
    count_changed = TypeAlloc(mp->n_atoms, int);      // record those atoms that do show some count change
    tmp_colors = TypeAlloc(mp->n_atoms, int);
    result = 0;
    if (size > 0)
    {
        // count number of bits in input fingerprint
        nbits = 0;
        nmasks = 0;
        for (j=0; j<size*8; j++)
            if ((fingerprint[j/8] & (1<<(j%8))) != 0) nbits++;
        /* for speed reasons, find pattern flags that apply to each atom to be tested */
        for (i=0; i<32; i++)
        {
            if ((which_bits & (1<<i)) == 0) continue;
            nmasks++;
            SetFingerprintCountsWithFocus(mp, mol_counts, size*8,
                                          which_bits & (1<<i),
                                          as_query == JNI_TRUE,
                                          0, 0);
            SetFingerprintCountsWithFocus(mp, mol_counts, size*8,
                                          which_bits & (1<<i),
                                          as_query == JNI_TRUE,
                                          ACCUMULATE_BITS | USE_DY_AROMATICITY, 0);
            for (j=0; j<size*8; j++)
                if (mol_counts[j] != 0)
                {
                    mol_flags[j] |= (1<<i);     // bit position j is lit by FP class i
                    flag_counts[j]++;           // accumulate number of flags affecting bit position j
                }
        }
// for (j=0; j<size*8; j++)
//    if (flag_counts[j] >= 1)
//       fprintf(stderr, "bit %d has flags %X\n", j, mol_flags[j]);

        // re-initialize counts for all requested pattern flags
        SetFingerprintCountsWithFocus(mp, mol_counts, size*8,
                                      which_bits, as_query == JNI_TRUE,
                                      0, 0);
        SetFingerprintCountsWithFocus(mp, mol_counts, size*8,
                                      which_bits, as_query == JNI_TRUE,
                                      ACCUMULATE_BITS | USE_DY_AROMATICITY, 0);

        /* first, we identify core atoms, i.e. those atoms required to set some of the bits */
        for (i=0; i<mp->n_atoms; i++)
        {
// fprintf(stderr, "fingerprinting atom %d of %d with fpsize %d\n", i+1, mp->n_atoms, size);
            SetFingerprintCountsWithFocus(mp, tmp_counts, size*8,
                                          which_bits, as_query == JNI_TRUE,
                                          0, i+1);
// fprintf(stderr, "fingerprinting with USE_DY_AROMATICITY\n");
            SetFingerprintCountsWithFocus(mp, tmp_counts, size*8,
                                          which_bits, as_query == JNI_TRUE,
                                          ACCUMULATE_BITS | USE_DY_AROMATICITY, i+1);
// fprintf(stderr, "collecting color\n");
// for (j=0; j<size; j++) fprintf(stderr, "%02X",
// fingerprint[j] ^ (fingerprint[j]&tmpfp[j]));
// fprintf(stderr, "\n");
            nzero = 0;
            for (j=0; j<size*8; j++)
               if ((fingerprint[j/8] & (1<<(j%8))) != 0)        // only look at relevant bit positions
               {
                  if (mol_counts[j] == 0) continue;
                  if (tmp_counts[j] == 0) nzero++;
                  if (mol_counts[j] != tmp_counts[j]) count_changed[i] = TRUE;
                  if (mol_counts[j] != 0  &&  tmp_counts[j] == 0) core_bits[j]++;       // this input bit is changed by a core atom
               }
// if (nzero > nbits/2.0  ||  nzero > 3)
            if (nzero > 0)
            {
                tmp_colors[i] = 1;
// fprintf(stderr, "removing '%s' atom %d clears %d bits\n", mp->atom_array[i].atom_symbol, i, nzero);
            }
        }

        /* now, we find additional important atoms */
        for (i=0; i<mp->n_atoms; i++)
        {
            if (tmp_colors[i] > 0) continue; /* core atoms => don't need to reanalyse */
            if (!count_changed[i]) continue; /* no change => do not reanalyze */
            SetFingerprintCountsWithFocus(mp, tmp_counts, size*8,
                                          which_bits, as_query == JNI_TRUE,
                                          0, i+1);
            SetFingerprintCountsWithFocus(mp, tmp_counts, size*8,
                                          which_bits, as_query == JNI_TRUE,
                                          ACCUMULATE_BITS | USE_DY_AROMATICITY, i+1);
// fprintf(stderr, "collecting color\n");
// for (j=0; j<size; j++) fprintf(stderr, "%02X",
// fingerprint[j] ^ (fingerprint[j]&tmpfp[j]));
// fprintf(stderr, "\n");
            ndiff = 0;
            for (j=0; j<size*8; j++)
                if ((fingerprint[j/8] & (1<<(j%8))) != 0)
                {
// fprintf(stderr, "mol_counts[%d] = %d, tmp_counts[%d] = %d\n", j, mol_counts[j], j, tmp_counts[j]);
                    if (mol_counts[j] == 0) continue;
                    if (tmp_counts[j] < 0.2*mol_counts[j]  ||        // at least 1/5th of the support for this bit depends on this atom
                        (tmp_counts[j] < mol_counts[j]  &&  nmasks == 1))        // single mask bit case is special
                    {
                       ndiff++;
// fprintf(stderr, "atom %d: bit %d has flags %X\n", i, j, mol_flags[j]);
                    }
                }
// fprintf(stderr, "removing '%s' atom %d changed some (%d) relevant bits\n", mp->atom_array[i].atom_symbol, i, ndiff);
            if ((ndiff > 0  &&  nmasks == 1)  ||       // one FP mask at a time case
                ndiff > 6  ||                              // this atom is supported by more than six bit changes, or
                (ndiff > 0  &&  ndiff >= nbits/2.0))  // this atom is supported by more than half of the input fingerprint
            {
                tmp_colors[i] = 2;
// fprintf(stderr, "'%s' atom %d has %d core bits updated\n", mp->atom_array[i].atom_symbol, i, ndiff);
            }
            if (tmp_colors[i] != 0) result++;
        }
    }
    MyFree((char *)tmp_counts);
    MyFree((char *)mol_counts);
    MyFree((char *)mol_flags);
    MyFree((char *)core_bits);
    MyFree((char *)count_changed);
    MyFree((char *)flag_counts);
// if (isCopy == JNI_TRUE)
    env->ReleaseByteArrayElements(buffer, fingerprint, 0);

// fprintf(stderr, "copying colors\n");
    size = env->GetArrayLength(color_buffer);
    isCopy = JNI_FALSE;
    colors = env->GetIntArrayElements(color_buffer, &isCopy);
    for (i=0; i<mp->n_atoms  &&  (int)i<size; i++)
        colors[i] = tmp_colors[i];
    MyFree((char *)tmp_colors);
// if (isCopy == JNI_TRUE)
    env->ReleaseIntArrayElements(color_buffer, colors, 0);
// fprintf(stderr, "%d atoms of %d are colored\n", result, mp->n_atoms);
    FreeMolecule(mp);
    return (result);
}

/*
 * Class:     JNIMatch
 * Method:    computeFingerprintFromCT
 * Signature: (Ljava/lang/String;[B)I
 */
JNIEXPORT jint JNICALL Java_jni_JNIMatch_computeFingerprintFromCT
  (JNIEnv *env, jobject thisObj,
   jstring CTObj, jbyteArray buffer, jboolean as_query, jint which_bits)
{
   const char *const_CT;
   Fortran_FILE *ffp;
   struct reaccs_molecule_t *mp;
   int result;
   jboolean isCopy = JNI_FALSE;
   int size;
   jbyte *fingerprint;
   unsigned int i;

   // Read CT
   const_CT = env->GetStringUTFChars(CTObj, &isCopy);
   ffp = FortranStringOpen((char *)const_CT);
   mp = TypeAlloc(1,struct reaccs_molecule_t);
   if (FORTRAN_NORMAL != ReadREACCSMolecule(ffp,mp,""))
   {
      fprintf(stderr, "failed to read template\n");
      FreeMolecule(mp);
      mp = (struct reaccs_molecule_t *)NULL;
   }
   FortranClose(ffp);
// if (isCopy == JNI_TRUE)
   env->ReleaseStringUTFChars(CTObj, const_CT);
   if (mp == NULL) return (0);

   // convert R-atoms to substitution indicators
   if (as_query ==  JNI_TRUE)
   {
      int has_R = FALSE;
      // quick check to avoid somewhat expensive processing
      for (i=0; i<mp->n_atoms; i++)
         if (0 == strcmp(mp->atom_array[i].atom_symbol, "R")) has_R = TRUE;
      if (has_R)
      {
         MakeHydrogensImplicit(mp);
         neighbourhood_t *nbp;
         nbp     = TypeAlloc(mp->n_atoms, neighbourhood_t);
         SetupNeighbourhood(mp,nbp,mp->n_atoms);
         PerceiveMarkush(mp, nbp);
         SetRingSizeFlags(mp, 14, nbp);
         MyFree((char *)nbp);
      }
   }

   size = env->GetArrayLength(buffer);
   fingerprint = env->GetByteArrayElements(buffer, &isCopy);
   result = SetFingerprintBits(mp, (char *)fingerprint, size,
                               which_bits, as_query == JNI_TRUE,
                               0);
   if (as_query == JNI_FALSE)
       result = SetFingerprintBits(mp, (char *)fingerprint, size,
                                   which_bits, as_query == JNI_TRUE,
                                   ACCUMULATE_BITS | USE_DY_AROMATICITY);
// if (isCopy == JNI_TRUE)
   env->ReleaseByteArrayElements(buffer, fingerprint, 0);
   FreeMolecule(mp);
   return (result);
}

/*
 * Class:     JNIMatch
 * Method:    computeFingerprintFromSmiles
 * Signature: (Ljava/lang/String;[B)I
 */
JNIEXPORT jint JNICALL Java_jni_JNIMatch_computeFingerprintFromSmiles
  (JNIEnv *env, jobject thisObj,
   jstring smilesObj, jbyteArray buffer, jboolean as_query, jint which_bits)
{
   const char *const_smiles;
   struct reaccs_molecule_t *mp;
   jboolean isCopy = JNI_FALSE;
   jint result = 0;
   int size;
   jbyte *fingerprint;

   const_smiles = env->GetStringUTFChars(smilesObj, &isCopy);
   mp = SMIToMOL(const_smiles, DROP_TRIVIAL_HYDROGENS);
// if (isCopy == JNI_TRUE)
   env->ReleaseStringUTFChars(smilesObj, const_smiles);
   /* Illegal molecules don't have fingerprint bits */
   if (mp == (struct reaccs_molecule_t *)NULL)
      return (result);

   size = env->GetArrayLength(buffer);
   fingerprint = env->GetByteArrayElements(buffer, &isCopy);
   result = SetFingerprintBits(mp, (char *)fingerprint, size,
                               which_bits, as_query == JNI_TRUE,
                               0);
   if (as_query == JNI_FALSE)
      result = SetFingerprintBits(mp, (char *)fingerprint, size,
                                  which_bits, as_query == JNI_TRUE,
                                  ACCUMULATE_BITS | USE_DY_AROMATICITY);
// if (isCopy == JNI_TRUE)
   env->ReleaseByteArrayElements(buffer, fingerprint, 0);
   FreeMolecule(mp);
   return (result);
}

#ifdef __cplusplus
}
#endif

