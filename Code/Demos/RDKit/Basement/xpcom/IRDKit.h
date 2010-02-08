/*
 * DO NOT EDIT.  THIS FILE IS GENERATED FROM IRDKit.idl
 */

#ifndef __gen_IRDKit_h__
#define __gen_IRDKit_h__


#ifndef __gen_nsISupports_h__
#include "nsISupports.h"
#endif

/* For IDL files that don't want to include root IDL files. */
#ifndef NS_NO_VTABLE
#define NS_NO_VTABLE
#endif

/* starting interface:    IRDMolecule */
#define IRDMOLECULE_IID_STR "2fd55049-0125-48be-88e6-270b1b83a8a8"

#define IRDMOLECULE_IID \
  {0x2fd55049, 0x0125, 0x48be, \
    { 0x88, 0xe6, 0x27, 0x0b, 0x1b, 0x83, 0xa8, 0xa8 }}

class NS_NO_VTABLE IRDMolecule : public nsISupports {
 public: 

  NS_DEFINE_STATIC_IID_ACCESSOR(IRDMOLECULE_IID)

  /* double GetMW (); */
  NS_IMETHOD GetMW(double *_retval) = 0;

  /* string GetSmiles (); */
  NS_IMETHOD GetSmiles(char **_retval) = 0;

  /* string GetMolBlock (); */
  NS_IMETHOD GetMolBlock(char **_retval) = 0;

  /* unsigned long GetSmartsMatchCount (in string smarts); */
  NS_IMETHOD GetSmartsMatchCount(const char *smarts, PRUint32 *_retval) = 0;

  /* double LogP (); */
  NS_IMETHOD LogP(double *_retval) = 0;

  /* double MR (); */
  NS_IMETHOD MR(double *_retval) = 0;

  /* void Generate3DCoords (); */
  NS_IMETHOD Generate3DCoords(void) = 0;

};

/* Use this macro when declaring classes that implement this interface. */
#define NS_DECL_IRDMOLECULE \
  NS_IMETHOD GetMW(double *_retval); \
  NS_IMETHOD GetSmiles(char **_retval); \
  NS_IMETHOD GetMolBlock(char **_retval); \
  NS_IMETHOD GetSmartsMatchCount(const char *smarts, PRUint32 *_retval); \
  NS_IMETHOD LogP(double *_retval); \
  NS_IMETHOD MR(double *_retval); \
  NS_IMETHOD Generate3DCoords(void); 

/* Use this macro to declare functions that forward the behavior of this interface to another object. */
#define NS_FORWARD_IRDMOLECULE(_to) \
  NS_IMETHOD GetMW(double *_retval) { return _to GetMW(_retval); } \
  NS_IMETHOD GetSmiles(char **_retval) { return _to GetSmiles(_retval); } \
  NS_IMETHOD GetMolBlock(char **_retval) { return _to GetMolBlock(_retval); } \
  NS_IMETHOD GetSmartsMatchCount(const char *smarts, PRUint32 *_retval) { return _to GetSmartsMatchCount(smarts, _retval); } \
  NS_IMETHOD LogP(double *_retval) { return _to LogP(_retval); } \
  NS_IMETHOD MR(double *_retval) { return _to MR(_retval); } \
  NS_IMETHOD Generate3DCoords(void) { return _to Generate3DCoords(); } 

/* Use this macro to declare functions that forward the behavior of this interface to another object in a safe way. */
#define NS_FORWARD_SAFE_IRDMOLECULE(_to) \
  NS_IMETHOD GetMW(double *_retval) { return !_to ? NS_ERROR_NULL_POINTER : _to->GetMW(_retval); } \
  NS_IMETHOD GetSmiles(char **_retval) { return !_to ? NS_ERROR_NULL_POINTER : _to->GetSmiles(_retval); } \
  NS_IMETHOD GetMolBlock(char **_retval) { return !_to ? NS_ERROR_NULL_POINTER : _to->GetMolBlock(_retval); } \
  NS_IMETHOD GetSmartsMatchCount(const char *smarts, PRUint32 *_retval) { return !_to ? NS_ERROR_NULL_POINTER : _to->GetSmartsMatchCount(smarts, _retval); } \
  NS_IMETHOD LogP(double *_retval) { return !_to ? NS_ERROR_NULL_POINTER : _to->LogP(_retval); } \
  NS_IMETHOD MR(double *_retval) { return !_to ? NS_ERROR_NULL_POINTER : _to->MR(_retval); } \
  NS_IMETHOD Generate3DCoords(void) { return !_to ? NS_ERROR_NULL_POINTER : _to->Generate3DCoords(); } 

#if 0
/* Use the code below as a template for the implementation class for this interface. */

/* Header file */
class _MYCLASS_ : public IRDMolecule
{
public:
  NS_DECL_ISUPPORTS
  NS_DECL_IRDMOLECULE

  _MYCLASS_();

private:
  ~_MYCLASS_();

protected:
  /* additional members */
};

/* Implementation file */
NS_IMPL_ISUPPORTS1(_MYCLASS_, IRDMolecule)

_MYCLASS_::_MYCLASS_()
{
  /* member initializers and constructor code */
}

_MYCLASS_::~_MYCLASS_()
{
  /* destructor code */
}

/* double GetMW (); */
NS_IMETHODIMP _MYCLASS_::GetMW(double *_retval)
{
    return NS_ERROR_NOT_IMPLEMENTED;
}

/* string GetSmiles (); */
NS_IMETHODIMP _MYCLASS_::GetSmiles(char **_retval)
{
    return NS_ERROR_NOT_IMPLEMENTED;
}

/* string GetMolBlock (); */
NS_IMETHODIMP _MYCLASS_::GetMolBlock(char **_retval)
{
    return NS_ERROR_NOT_IMPLEMENTED;
}

/* unsigned long GetSmartsMatchCount (in string smarts); */
NS_IMETHODIMP _MYCLASS_::GetSmartsMatchCount(const char *smarts, PRUint32 *_retval)
{
    return NS_ERROR_NOT_IMPLEMENTED;
}

/* double LogP (); */
NS_IMETHODIMP _MYCLASS_::LogP(double *_retval)
{
    return NS_ERROR_NOT_IMPLEMENTED;
}

/* double MR (); */
NS_IMETHODIMP _MYCLASS_::MR(double *_retval)
{
    return NS_ERROR_NOT_IMPLEMENTED;
}

/* void Generate3DCoords (); */
NS_IMETHODIMP _MYCLASS_::Generate3DCoords()
{
    return NS_ERROR_NOT_IMPLEMENTED;
}

/* End of implementation class template. */
#endif


/* starting interface:    IRDMolSupplier */
#define IRDMOLSUPPLIER_IID_STR "056a8da1-7820-41d7-b254-5ef7dd1693ce"

#define IRDMOLSUPPLIER_IID \
  {0x056a8da1, 0x7820, 0x41d7, \
    { 0xb2, 0x54, 0x5e, 0xf7, 0xdd, 0x16, 0x93, 0xce }}

class NS_NO_VTABLE IRDMolSupplier : public nsISupports {
 public: 

  NS_DEFINE_STATIC_IID_ACCESSOR(IRDMOLSUPPLIER_IID)

  /* boolean atEnd (); */
  NS_IMETHOD AtEnd(PRBool *_retval) = 0;

  /* IRDMolecule next (); */
  NS_IMETHOD Next(IRDMolecule **_retval) = 0;

};

/* Use this macro when declaring classes that implement this interface. */
#define NS_DECL_IRDMOLSUPPLIER \
  NS_IMETHOD AtEnd(PRBool *_retval); \
  NS_IMETHOD Next(IRDMolecule **_retval); 

/* Use this macro to declare functions that forward the behavior of this interface to another object. */
#define NS_FORWARD_IRDMOLSUPPLIER(_to) \
  NS_IMETHOD AtEnd(PRBool *_retval) { return _to AtEnd(_retval); } \
  NS_IMETHOD Next(IRDMolecule **_retval) { return _to Next(_retval); } 

/* Use this macro to declare functions that forward the behavior of this interface to another object in a safe way. */
#define NS_FORWARD_SAFE_IRDMOLSUPPLIER(_to) \
  NS_IMETHOD AtEnd(PRBool *_retval) { return !_to ? NS_ERROR_NULL_POINTER : _to->AtEnd(_retval); } \
  NS_IMETHOD Next(IRDMolecule **_retval) { return !_to ? NS_ERROR_NULL_POINTER : _to->Next(_retval); } 

#if 0
/* Use the code below as a template for the implementation class for this interface. */

/* Header file */
class _MYCLASS_ : public IRDMolSupplier
{
public:
  NS_DECL_ISUPPORTS
  NS_DECL_IRDMOLSUPPLIER

  _MYCLASS_();

private:
  ~_MYCLASS_();

protected:
  /* additional members */
};

/* Implementation file */
NS_IMPL_ISUPPORTS1(_MYCLASS_, IRDMolSupplier)

_MYCLASS_::_MYCLASS_()
{
  /* member initializers and constructor code */
}

_MYCLASS_::~_MYCLASS_()
{
  /* destructor code */
}

/* boolean atEnd (); */
NS_IMETHODIMP _MYCLASS_::AtEnd(PRBool *_retval)
{
    return NS_ERROR_NOT_IMPLEMENTED;
}

/* IRDMolecule next (); */
NS_IMETHODIMP _MYCLASS_::Next(IRDMolecule **_retval)
{
    return NS_ERROR_NOT_IMPLEMENTED;
}

/* End of implementation class template. */
#endif


/* starting interface:    IRDKit */
#define IRDKIT_IID_STR "bfb9acf3-9349-47ec-8984-f6f8e2f02f65"

#define IRDKIT_IID \
  {0xbfb9acf3, 0x9349, 0x47ec, \
    { 0x89, 0x84, 0xf6, 0xf8, 0xe2, 0xf0, 0x2f, 0x65 }}

class NS_NO_VTABLE IRDKit : public nsISupports {
 public: 

  NS_DEFINE_STATIC_IID_ACCESSOR(IRDKIT_IID)

  /* unsigned long strlen (in string arg); */
  NS_IMETHOD Strlen(const char *arg, PRUint32 *_retval) = 0;

  /* IRDMolecule MolFromSmiles (in string smiles); */
  NS_IMETHOD MolFromSmiles(const char *smiles, IRDMolecule **_retval) = 0;

  /* IRDMolecule MolFromMolBlock (in string molBlock); */
  NS_IMETHOD MolFromMolBlock(const char *molBlock, IRDMolecule **_retval) = 0;

  /* IRDMolSupplier SupplierFromSDFile (in string fileName); */
  NS_IMETHOD SupplierFromSDFile(const char *fileName, IRDMolSupplier **_retval) = 0;

};

/* Use this macro when declaring classes that implement this interface. */
#define NS_DECL_IRDKIT \
  NS_IMETHOD Strlen(const char *arg, PRUint32 *_retval); \
  NS_IMETHOD MolFromSmiles(const char *smiles, IRDMolecule **_retval); \
  NS_IMETHOD MolFromMolBlock(const char *molBlock, IRDMolecule **_retval); \
  NS_IMETHOD SupplierFromSDFile(const char *fileName, IRDMolSupplier **_retval); 

/* Use this macro to declare functions that forward the behavior of this interface to another object. */
#define NS_FORWARD_IRDKIT(_to) \
  NS_IMETHOD Strlen(const char *arg, PRUint32 *_retval) { return _to Strlen(arg, _retval); } \
  NS_IMETHOD MolFromSmiles(const char *smiles, IRDMolecule **_retval) { return _to MolFromSmiles(smiles, _retval); } \
  NS_IMETHOD MolFromMolBlock(const char *molBlock, IRDMolecule **_retval) { return _to MolFromMolBlock(molBlock, _retval); } \
  NS_IMETHOD SupplierFromSDFile(const char *fileName, IRDMolSupplier **_retval) { return _to SupplierFromSDFile(fileName, _retval); } 

/* Use this macro to declare functions that forward the behavior of this interface to another object in a safe way. */
#define NS_FORWARD_SAFE_IRDKIT(_to) \
  NS_IMETHOD Strlen(const char *arg, PRUint32 *_retval) { return !_to ? NS_ERROR_NULL_POINTER : _to->Strlen(arg, _retval); } \
  NS_IMETHOD MolFromSmiles(const char *smiles, IRDMolecule **_retval) { return !_to ? NS_ERROR_NULL_POINTER : _to->MolFromSmiles(smiles, _retval); } \
  NS_IMETHOD MolFromMolBlock(const char *molBlock, IRDMolecule **_retval) { return !_to ? NS_ERROR_NULL_POINTER : _to->MolFromMolBlock(molBlock, _retval); } \
  NS_IMETHOD SupplierFromSDFile(const char *fileName, IRDMolSupplier **_retval) { return !_to ? NS_ERROR_NULL_POINTER : _to->SupplierFromSDFile(fileName, _retval); } 

#if 0
/* Use the code below as a template for the implementation class for this interface. */

/* Header file */
class _MYCLASS_ : public IRDKit
{
public:
  NS_DECL_ISUPPORTS
  NS_DECL_IRDKIT

  _MYCLASS_();

private:
  ~_MYCLASS_();

protected:
  /* additional members */
};

/* Implementation file */
NS_IMPL_ISUPPORTS1(_MYCLASS_, IRDKit)

_MYCLASS_::_MYCLASS_()
{
  /* member initializers and constructor code */
}

_MYCLASS_::~_MYCLASS_()
{
  /* destructor code */
}

/* unsigned long strlen (in string arg); */
NS_IMETHODIMP _MYCLASS_::Strlen(const char *arg, PRUint32 *_retval)
{
    return NS_ERROR_NOT_IMPLEMENTED;
}

/* IRDMolecule MolFromSmiles (in string smiles); */
NS_IMETHODIMP _MYCLASS_::MolFromSmiles(const char *smiles, IRDMolecule **_retval)
{
    return NS_ERROR_NOT_IMPLEMENTED;
}

/* IRDMolecule MolFromMolBlock (in string molBlock); */
NS_IMETHODIMP _MYCLASS_::MolFromMolBlock(const char *molBlock, IRDMolecule **_retval)
{
    return NS_ERROR_NOT_IMPLEMENTED;
}

/* IRDMolSupplier SupplierFromSDFile (in string fileName); */
NS_IMETHODIMP _MYCLASS_::SupplierFromSDFile(const char *fileName, IRDMolSupplier **_retval)
{
    return NS_ERROR_NOT_IMPLEMENTED;
}

/* End of implementation class template. */
#endif


#endif /* __gen_IRDKit_h__ */
