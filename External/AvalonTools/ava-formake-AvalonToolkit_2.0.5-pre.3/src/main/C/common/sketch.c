//
//  Copyright (c) 2010, Novartis Institutes for BioMedical Research Inc.
//  All rights reserved.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met: 
//
//     * Redistributions of source code must retain the above copyright 
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following 
//       disclaimer in the documentation and/or other materials provided 
//       with the distribution.
//     * Neither the name of Novartis Institutes for BioMedical Research Inc. 
//       nor the names of its contributors may be used to endorse or promote 
//       products derived from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
/************************************************************************/
/*                                                                      */
/*  File:     sketch.c                                                  */
/*                                                                      */
/*  Purpose:  Module to implement writing (and later reading) of MDL    */
/*            .skc sketch format files/buffers.                         */
/*                                                                      */
/************************************************************************/

#include <stdio.h>
#include <string.h>

#include "reaccs.h"
#include "reaccsio.h"

/**
 * Helper function to add a record tbuffer[0..tsize-1] to the beginning of
 * (*buffer)[0..*bufsize-1].
 *
 * The function updates (*buffer), *bufsize, and *size t reflect the addition.
 * TRUE is returned if everything went OK and FALSE if the new record would
 * overrun the buffer. No data is written in this case but *size is updated
 * to inform the caller of the required size.
 */
static int addRecord(char **buffer, int *bufsize, int *size,
                     char *tbuffer, int tsize)
{
    if (tsize <= *bufsize)      /* record fits in */
    {
        memcpy(*buffer, tbuffer, tsize);
        (*size) += tsize; (*buffer) += tsize; (*bufsize) -= tsize;
        return (TRUE);
    }
    else                        /* record does not fit */
    {
        (*size) += tsize;
        return (FALSE);
    }
}

#define VERSION 1
/**
 * The function writes a version record into buffer[0..bsize-1].
 */
static int versionRecord(char *buffer, int bsize)
{
    if (bsize < 4) return (0);
    buffer[0] = VERSION;                /* $Version */
    buffer[1] = 3; buffer[2] = 0;       /* size of data */
    buffer[3] = 4;                      /* version number */
    return (4);
}

#define TOTOBJS 2
/**
 * The function writes a Totobjs record into buffer[0..bsize-1].
 */
static int totobjsRecord(char *buffer, int bsize, int totobjs)
{
    if (bsize < 5) return (0);
    buffer[0] = TOTOBJS;                /* $Totobjs */
    buffer[1] = 4; buffer[2] = 0;       /* size of data */
    buffer[3] = totobjs & 0xFF;         /* low order byte */
    buffer[4] = totobjs >> 8;           /* high order byte */
    return (5);
}

#define OBJ 3
#define TYPE_MOL        12
#define TYPE_CHIRAL     19
#define TYPE_ATOM       33
#define TYPE_BOND       32
/**
 * The function writes an Obj record into buffer[0..bsize-1].
 */
static int objRecord(char *buffer, int bsize, int type)
{
    if (bsize < 4) return (0);
    buffer[0] = OBJ;                    /* $Obj */
    buffer[1] = 3; buffer[2] = 0;       /* size of data */
    buffer[3] = type & 0xFF;            /* type of object */
    return (4);
}

#define BEGSKETCH       19
/**
 * The function writes an Begsketch record into buffer[0..bsize-1].
 */
static int begsketchRecord(char *buffer, int bsize)
{
    if (bsize < 3) return (0);
    buffer[0] = BEGSKETCH;              /* $Begsketch */
    buffer[1] = 2; buffer[2] = 0;       /* size of data */
    return (3);
}

#define ENDSKETCH       20
/**
 * The function writes an Begsketch record into buffer[0..bsize-1].
 */
static int endsketchRecord(char *buffer, int bsize)
{
    if (bsize < 3) return (0);
    buffer[0] = ENDSKETCH;              /* $Endsketch */
    buffer[1] = 2; buffer[2] = 0;       /* size of data */
    return (3);
}

#define FONT 11
/**
 * The function writes a Font record into buffer[0..bsize-1].
 *
 * Note: font_size is measured in 1/10th of a point, e.g. 100 for 10pt font.
 */
static int fontRecord(char *buffer, int bsize,
                      char *font_name, int font_size, int font_flags)
{
    int data_size;
    
    data_size = 2+1+strlen(font_name) + 2 + 1;
    if (bsize < 1 + data_size) return (0);
    buffer[0] = FONT;                   /* $Font */
    buffer++;

    buffer[0] = data_size & 0xFF;       /* low order byte of data_size */
    buffer[1] = data_size >> 8;         /* high order byte of data_size */
    buffer += 2;

    buffer[0] = strlen(font_name);
    buffer++;
    strcpy(buffer, font_name);
    buffer += strlen(font_name);

    buffer[0] = font_size & 0xFF;         /* low order byte of font_size */
    buffer[1] = font_size >> 8;           /* high order byte of font_size */
    buffer += 2;

    buffer[0] = font_flags;
    return (1 + data_size);
}

union number_t
   {
      unsigned char bytestring[4];
      float         single_num;
      long          long_num;
   };

#define ATOM_COORDS 22
/**
 * The function writes an atom coordinates record into buffer[0..bsize-1].
 */
static int atom_coordsRecord(char *buffer, int bsize,
                             float x, float y)
{
    int data_size;
    union number_t number_union;

    data_size = 2+4+4;
    if (bsize < 1 + data_size) return (0);
    buffer[0] = ATOM_COORDS;                   /* $Atom_coords */
    buffer++;

    buffer[0] = data_size & 0xFF;       /* low order byte of data_size */
    buffer[1] = data_size >> 8;         /* high order byte of data_size */
    buffer += 2;

    number_union.single_num = x*240+60;
    memcpy(buffer, number_union.bytestring, 4);
    buffer += 4;

    number_union.single_num = y*240+60;
    memcpy(buffer, number_union.bytestring, 4);
    buffer += 4;

    return (1 + data_size);
}

#define BOND_ATOMS 39
/**
 * The function writes a bond atoms record into buffer[0..bsize-1].
 */
static int bond_atomsRecord(char *buffer, int bsize,
                            int at1, int at2)
{
    int data_size;

    data_size = 2+4;
    if (bsize < 1 + data_size) return (0);
    buffer[0] = BOND_ATOMS;                   /* $Bond_atoms */
    buffer++;

    buffer[0] = data_size & 0xFF;       /* low order byte of data_size */
    buffer[1] = data_size >> 8;         /* high order byte of data_size */
    buffer += 2;

    buffer[0] = at1 & 0xFF;
    buffer[1] = (at1 >> 8) & 0xFF;
    buffer += 2;

    buffer[0] = at2 & 0xFF;
    buffer[1] = (at2 >> 8) & 0xFF;
    buffer += 2;

    return (1 + data_size);
}

static char* bond_type[] =
{
    "none",
    "single",
    "double",
    "triple",
    "aromatic",
    "sgl/dbl",
    "sgl/arom",
    "dbl/arom",
    "any"
};

#define BOND_TYPE 40
/**
 * The function writes a bond type record into buffer[0..bsize-1].
 */
static int bond_typeRecord(char *buffer, int bsize, int type)
{
    int data_size;

    data_size = 2+1;
    if (bsize < 1 + data_size) return (0);
    buffer[0] = BOND_TYPE;                   /* $Bond_type */
    buffer++;

    buffer[0] = data_size & 0xFF;       /* low order byte of data_size */
    buffer[1] = data_size >> 8;         /* high order byte of data_size */
    buffer += 2;

    buffer[0] = type;                   /* bond type is one byte */
    buffer++;

    return (1 + data_size);
}

#define BOND_STEREO_TYPE 41
/**
 * The function writes a bond stereo type record into buffer[0..bsize-1].
 */
static int bond_stereo_typeRecord(char *buffer, int bsize, int stereo_type)
{
    int data_size;

    data_size = 2+1;
    if (bsize < 1 + data_size) return (0);
    buffer[0] = BOND_STEREO_TYPE;                   /* $Bond_stereo_type */
    buffer++;

    buffer[0] = data_size & 0xFF;       /* low order byte of data_size */
    buffer[1] = data_size >> 8;         /* high order byte of data_size */
    buffer += 2;

    buffer[0] = stereo_type;            /* bond stereo type is one byte */
    buffer++;

    return (1 + data_size);
}

#define ATOM_SYMBOL 232
/**
 * The function writes an atom symbol record into buffer[0..bsize-1].
 *
 * Note: font_size is measured in 1/10th of a point, e.g. 100 for 10pt font.
 */
static int atom_symbolRecord(char *buffer, int bsize,
                             char *symbol)
{
    int data_size;
    
    data_size = 2+1+strlen(symbol);
    if (bsize < 1 + data_size) return (0);
    buffer[0] = ATOM_SYMBOL;                   /* $Atom_symbol */
    buffer++;

    buffer[0] = data_size & 0xFF;       /* low order byte of data_size */
    buffer[1] = data_size >> 8;         /* high order byte of data_size */
    buffer += 2;

    buffer[0] = strlen(symbol);
    buffer++;
    strcpy(buffer, symbol);
    buffer += strlen(symbol);

    return (1 + data_size);
}

int  MoleculeToSKCBuffer(struct reaccs_molecule_t *mp,
                         char *buffer, int bufsize)
/*
 * Convert chemical structure into an SKC byte stream.
 *
 * This byte stream can then be written to a binary .SKC file or
 * even be used to be put on the clipboard as MDLSK and/or Native
 * Format.
 *
 * The function returns the real size required to hold the sketch
 * or <=0 if the conversion was unsuccessful. If the return value is < 0
 * then then negative of it is the buffer size that would be required.
 * the buffer is, however, overwritten.
 */
{
    char tbuffer[255];
    int tsize = 255;
    int i, tmp, size;
    int overflow = FALSE;
    float minx, maxy;

    size = 0;

    tmp = versionRecord(tbuffer, tsize);
    if (!addRecord(&buffer, &bufsize, &size, tbuffer, tmp)) overflow = TRUE;

    tmp = totobjsRecord(tbuffer, tsize, 1);
    if (!addRecord(&buffer, &bufsize, &size, tbuffer, tmp)) overflow = TRUE;

    tmp = fontRecord(tbuffer, tsize, "Arial", 100, NONE);
    if (!addRecord(&buffer, &bufsize, &size, tbuffer, tmp)) overflow = TRUE;

    tmp = objRecord(tbuffer, tsize, TYPE_MOL);
    if (!addRecord(&buffer, &bufsize, &size, tbuffer, tmp)) overflow = TRUE;

    tmp = begsketchRecord(tbuffer, tsize);
    if (!addRecord(&buffer, &bufsize, &size, tbuffer, tmp)) overflow = TRUE;

    tmp = totobjsRecord(tbuffer, tsize, mp->n_atoms + mp->n_bonds);
    if (!addRecord(&buffer, &bufsize, &size, tbuffer, tmp)) overflow = TRUE;

    minx = 1.0e7;
    maxy = -1.0e7;
    for (i=0; i<mp->n_atoms; i++)
    {
        if (mp->atom_array[i].x < minx) minx = mp->atom_array[i].x;
        if (mp->atom_array[i].y > maxy) maxy = mp->atom_array[i].y;
    }
    for (i=0; i<mp->n_atoms; i++)
    {
        tmp = objRecord(tbuffer, tsize, TYPE_ATOM);
        if (!addRecord(&buffer, &bufsize, &size, tbuffer, tmp)) overflow = TRUE;
        tmp = atom_coordsRecord(tbuffer, tsize,
                                mp->atom_array[i].x - minx,
                                maxy - mp->atom_array[i].y);
        if (!addRecord(&buffer, &bufsize, &size, tbuffer, tmp)) overflow = TRUE;
        if (strcmp(mp->atom_array[i].atom_symbol, "C") != 0)
        {
            tmp = atom_symbolRecord(tbuffer, tsize,
                                    mp->atom_array[i].atom_symbol);
            if (!addRecord(&buffer, &bufsize, &size, tbuffer, tmp))
                overflow = TRUE;
        }
    }

    for (i=0; i<mp->n_bonds; i++)
    {
        tmp = objRecord(tbuffer, tsize, TYPE_BOND);
        if (!addRecord(&buffer, &bufsize, &size, tbuffer, tmp)) overflow = TRUE;
        tmp = bond_atomsRecord(tbuffer, tsize,
                               mp->bond_array[i].atoms[0]-1,
                               mp->bond_array[i].atoms[1]-1);
        if (!addRecord(&buffer, &bufsize, &size, tbuffer, tmp)) overflow = TRUE;
        if (mp->bond_array[i].bond_type != SINGLE)
        {
            tmp = bond_typeRecord(tbuffer, tsize,
                                  mp->bond_array[i].bond_type);
            if (!addRecord(&buffer, &bufsize, &size, tbuffer, tmp))
                overflow = TRUE;
        }
        if (mp->bond_array[i].stereo_symbol != NONE)
        {
            tmp = bond_stereo_typeRecord(tbuffer, tsize,
                                         mp->bond_array[i].stereo_symbol);
            if (!addRecord(&buffer, &bufsize, &size, tbuffer, tmp))
                overflow = TRUE;
        }
    }

    tmp = endsketchRecord(tbuffer, tsize);
    if (!addRecord(&buffer, &bufsize, &size, tbuffer, tmp)) overflow = TRUE;

    tmp = endsketchRecord(tbuffer, tsize);
    if (!addRecord(&buffer, &bufsize, &size, tbuffer, tmp)) overflow = TRUE;

    if (overflow) return (-size);
    else          return (size);
}

#ifdef MAIN
#define NBUFFER 4000
int main(int argc, char *argv[])
{
    struct reaccs_molecule_t *mp = NULL;
    FILE *fp;
    Fortran_FILE *ffp;
    char buffer[NBUFFER];
    int size;

    ffp = FortranOpen(argv[1], "r");
    fp = fopen(argv[2], "wb");
    mp = TypeAlloc(1,struct reaccs_molecule_t);
    ReadREACCSMolecule(ffp,mp,"");
    size = MoleculeToSKCBuffer(mp, buffer, NBUFFER);
    if (size > 0)
        fwrite(buffer, sizeof(char), size, fp);
    else
        fwrite(buffer, sizeof(char), NBUFFER, fp);
    fclose(fp);
}
#endif
