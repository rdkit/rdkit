#ifndef _RDKITIMPL_H_
#define _RDKITIMPL_H_
#include "IRDKit.h"

/* Header file */
class RDKitImpl : public IRDKit
{
public:
  NS_DECL_ISUPPORTS
  NS_DECL_IRDKIT

  RDKitImpl();

private:
  ~RDKitImpl();

protected:
  /* additional members */
};
#endif
