/*
*
*  Copyright (c) 2015, Novartis Institutes for BioMedical Research Inc.
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
%include "boost_tuple.i"
%include "std_vector.i"
%include "std_map.i"
%include "std_pair.i"
%include "boost_tuple.i"
%{
#include <GraphMol/RDKitBase.h>
#include <GraphMol/MolDraw2D/MolDraw2DHelpers.h>
#include <GraphMol/MolDraw2D/MolDraw2D.h>
#include <GraphMol/MolDraw2D/MolDraw2DSVG.h>
#include <GraphMol/MolDraw2D/MolDraw2DUtils.h>
#include <map>
#include <vector>
#include <utility>
%}

#ifdef RDK_BUILD_CAIRO_SUPPORT
%{
#include <GraphMol/MolDraw2D/MolDraw2DCairo.h>
 %}
#endif

%template(Int_String_Map) std::map< int, std::string >;

%template(ColourPalette) std::map< int, RDKit::DrawColour >;

%template(Int_Double_Map) std::map< int, double >;
%template(Float_Pair) std::pair<float,float>;
%template(Float_Pair_Vect) std::vector< std::pair<float,float> >;
%template(ROMol_Ptr_Vect) std::vector<RDKit::ROMol*>;
%template(Point2D_Vect) std::vector<RDGeom::Point2D *>;
%template(ColourPalette_Vect) std::vector< std::map< int, RDKit::DrawColour > >;

%ignore RDKit::MolDraw2DSVG::MolDraw2DSVG(int,int,std::ostream &);
%ignore RDKit::MolDraw2DUtils::contourAndDrawGaussians(
    MolDraw2D &, const std::vector<Point2D> &,
    const std::vector<double> &, const std::vector<double> &,
    size_t, std::vector<double> &,
    const ContourParams &);
%ignore RDKit::MolDraw2DUtils::contourAndDrawGaussians(
    MolDraw2D &, const std::vector<Point2D> &,
    const std::vector<double> &, const std::vector<double> &,
    size_t, ContourParams &);
%ignore RDKit::MolDraw2DUtils::contourAndDrawGaussians;
%ignore RDKit::MolDraw2DUtils::contourAndDrawGrid;


%include <GraphMol/MolDraw2D/MolDraw2DHelpers.h>
%include <GraphMol/MolDraw2D/MolDraw2D.h>
%include <GraphMol/MolDraw2D/MolDraw2DSVG.h>
#ifdef RDK_BUILD_CAIRO_SUPPORT
    
#ifdef SWIGJAVA
%typemap(jni) std::string RDKit::MolDraw2DCairo::toByteArray "jbyteArray"
%typemap(jtype) std::string RDKit::MolDraw2DCairo::toByteArray "byte[]"
%typemap(jstype) std::string RDKit::MolDraw2DCairo::toByteArray "byte[]"
%typemap(javaout) std::string RDKit::MolDraw2DCairo::toByteArray {
  return $jnicall;
}
%typemap(out) std::string RDKit::MolDraw2DCairo::toByteArray {
  $result = JCALL1(NewByteArray, jenv, $1.size());
  JCALL4(SetByteArrayRegion, jenv, $result, 0, $1.size(), (const jbyte*)$1.c_str());
}
#endif

#ifdef SWIGCSHARP
%template(UChar_Vect) std::vector<unsigned char>;
#endif
          
%include <GraphMol/MolDraw2D/MolDraw2DCairo.h>
    
#ifdef SWIGJAVA
%extend RDKit::MolDraw2DCairo {
  const std::string toByteArray() {
     return ($self)->getDrawingText();
  }
}
#endif

#ifdef SWIGCSHARP
%extend RDKit::MolDraw2DCairo {
    const std::vector<unsigned char> getImage() {
        const auto text = ($self)->getDrawingText();
        const std::vector<unsigned char> image(text.begin(), text.end());
        return image;
    }
};
#endif
        
#endif // RDK_BUILD_CAIRO_SUPPORT
        
%include <GraphMol/MolDraw2D/MolDraw2DUtils.h>

%inline %{
    void ContourAndDrawGaussians(
    RDKit::MolDraw2D &drawer, const std::vector<RDGeom::Point2D *> &p_locs,
    const std::vector<double> &heights, const std::vector<double> &widths,
    size_t nContours, std::vector<double> &levels,
    const RDKit::MolDraw2DUtils::ContourParams &ps = RDKit::MolDraw2DUtils::ContourParams()){
        std::vector<Point2D> locs;
        locs.reserve(p_locs.size());
        for(const auto *p : p_locs){
            locs.push_back(*p);
        }
        contourAndDrawGaussians(drawer,locs,heights,widths,nContours,levels,ps);
    };
    void ContourAndDrawGaussians(
    RDKit::MolDraw2D &drawer, const std::vector<RDGeom::Point2D *> &p_locs,
    const std::vector<double> &heights, const std::vector<double> &widths,
    size_t nContours = 10, 
    const RDKit::MolDraw2DUtils::ContourParams &ps = RDKit::MolDraw2DUtils::ContourParams()){
        std::vector<Point2D> locs;
        locs.reserve(p_locs.size());
        for(const auto *p : p_locs){
            locs.push_back(*p);
        }
        contourAndDrawGaussians(drawer,locs,heights,widths,nContours,ps);
    };
    void ContourAndDrawGrid(
    RDKit::MolDraw2D &drawer, const std::vector<double> &p_grid,
    const std::vector<double> &xcoords, const std::vector<double> &ycoords,
    size_t nContours, std::vector<double> &levels,
    const RDKit::MolDraw2DUtils::ContourParams &ps = RDKit::MolDraw2DUtils::ContourParams()){
        contourAndDrawGrid(drawer,p_grid.data(),xcoords,ycoords,nContours,levels,ps);
    };
    void ContourAndDrawGrid(
    RDKit::MolDraw2D &drawer, const std::vector<double> &p_grid,
    const std::vector<double> &xcoords, const std::vector<double> &ycoords,
    size_t nContours=10,
    const RDKit::MolDraw2DUtils::ContourParams &ps = RDKit::MolDraw2DUtils::ContourParams()){
        contourAndDrawGrid(drawer,p_grid.data(),xcoords,ycoords,nContours,ps);
    };
%}
