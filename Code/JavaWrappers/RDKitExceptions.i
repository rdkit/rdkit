/* 
* $Id$
*
*  Copyright (c) 2010, Novartis Institutes for BioMedical Research Inc.
*  All rights reserved.
* 
* Redistribution and use in source and binary forms, with or without
* modification, are permitted provided that the following conditions are
* met: 
*
*     * Redistributions of source code must retain the above copyright 
*       notice, this list of conditions and the following disclaimer.
*     * Redistributions in binary form must reproduce the above
*       copyright notice, this list of conditions and the following 
*       disclaimer in the documentation and/or other materials provided 
*       with the distribution.
*     * Neither the name of Novartis Institutes for BioMedical Research Inc. 
*       nor the names of its contributors may be used to endorse or promote 
*       products derived from this software without specific prior written permission.
*
* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
* "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
* LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
* A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
* OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
* SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
* LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
* DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
* THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
* (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
* OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

// Defines all of the C++ -> Java exception handling

%{
#include <RDBoost/Exceptions.h>
#include <JavaWrappers/GenericRDKitException.h>
%}

// ===== ChemicalReactionException =====
%typemap(javabase) RDKit::ChemicalReactionException "java.lang.RuntimeException";
%typemap(throws, throws="org.RDKit.ChemicalReactionException") RDKit::ChemicalReactionException {
  jclass excep = jenv->FindClass("org/RDKit/ChemicalReactionException");
  if (excep)
    jenv->ThrowNew(excep, $1.message());
  return $null;
}

// ===== ChemicalReactionParserException =====
%typemap(javabase) RDKit::ChemicalReactionParserException "java.lang.RuntimeException";
%typemap(throws, throws="org.RDKit.ChemicalReactionParserException") RDKit::ChemicalReactionParserException {
  jclass excep = jenv->FindClass("org/RDKit/ChemicalReactionParserException");
  if (excep)
    jenv->ThrowNew(excep, $1.message());
  return $null;
}

// ===== ConformerException =====
%typemap(javabase) RDKit::ConformerException "java.lang.RuntimeException";
%typemap(throws, throws="org.RDKit.ConformerException") RDKit::ConformerException {
  jclass excep = jenv->FindClass("org/RDKit/ConformerException");
  if (excep)
    jenv->ThrowNew(excep, $1.message());
  return $null;
}

// ===== MolPicklerException =====
%typemap(javabase) RDKit::MolPicklerException "java.lang.RuntimeException";
%typemap(throws, throws="org.RDKit.MolPicklerException") RDKit::MolPicklerException {
  jclass excep = jenv->FindClass("org/RDKit/MolPicklerException");
  if (excep)
    jenv->ThrowNew(excep, $1.message());
  return $null;
}

// ===== MolSanitizeException =====
%typemap(javabase) RDKit::MolSanitizeException "java.lang.RuntimeException";
%typemap(throws, throws="org.RDKit.MolSanitizeException") RDKit::MolSanitizeException {
  jclass excep = jenv->FindClass("org/RDKit/MolSanitizeException");
  if (excep)
    jenv->ThrowNew(excep, $1.message());
  return $null;
}

// ===== SmilesParseException =====
%typemap(javabase) RDKit::SmilesParseException "java.lang.RuntimeException";
%typemap(throws, throws="org.RDKit.SmilesParseException") RDKit::SmilesParseException {
  jclass excep = jenv->FindClass("org/RDKit/SmilesParseException");
  if (excep)
    jenv->ThrowNew(excep, $1.message());
  return $null;
}

// ===== KeyErrorException =====
%typemap(javabase) KeyErrorException "java.lang.RuntimeException";
%typemap(throws, throws="org.RDKit.KeyErrorException") KeyErrorException {
  jclass excep = jenv->FindClass("org/RDKit/KeyErrorException");
  if (excep)
    jenv->ThrowNew(excep, $1.key());
  return $null;
}
%extend KeyErrorException {
	std::string message() {
		return "Unknown key: " + ($self)->key();
	}
}

// ===== GenericRDKitException =====
%typemap(javabase) RDKit::GenericRDKitException "java.lang.RuntimeException";
%typemap(throws, throws="org.RDKit.GenericRDKitException") RDKit::GenericRDKitException {
  jclass excep = jenv->FindClass("org/RDKit/GenericRDKitException");
  if (excep)
    jenv->ThrowNew(excep, $1.message());
  return $null;
}

// Note that these files must follow the typemap declarations
%include <RDBoost/Exceptions.h>
%include <JavaWrappers/GenericRDKitException.h>

%exception {
  try {
     $action
  } catch (RDKit::ChemicalReactionException &e) {
    jclass clazz = jenv->FindClass("org/RDKit/ChemicalReactionException");
    jenv->ThrowNew(clazz, e.message());
    return $null;
  } catch (RDKit::ChemicalReactionParserException &e) {
    jclass clazz = jenv->FindClass("org/RDKit/ChemicalReactionParserException");
    jenv->ThrowNew(clazz, e.message());
    return $null;
  } catch (RDKit::ConformerException &e) {
    jclass clazz = jenv->FindClass("org/RDKit/ConformerException");
    jenv->ThrowNew(clazz, e.message());
    return $null;
  } catch (RDKit::MolPicklerException &e) {
    jclass clazz = jenv->FindClass("org/RDKit/MolPicklerException");
    jenv->ThrowNew(clazz, e.message());
    return $null;
  } catch (RDKit::MolSanitizeException &e) {
    jclass clazz = jenv->FindClass("org/RDKit/MolSanitizeException");
    jenv->ThrowNew(clazz, e.message());
    return $null;
  } catch (RDKit::SmilesParseException &e) {
    jclass clazz = jenv->FindClass("org/RDKit/SmilesParseException");
    jenv->ThrowNew(clazz, e.message());
    return $null;
  } catch (KeyErrorException &e) {
    jclass clazz = jenv->FindClass("org/RDKit/KeyErrorException");
    jenv->ThrowNew(clazz, e.key().c_str());
    return $null;

  // Generic exception -- anything else
  } catch (std::exception &e) {
    jclass clazz = jenv->FindClass("org/RDKit/GenericRDKitException");
    jenv->ThrowNew(clazz, "Unknown exception");
    return $null;

  }

}
