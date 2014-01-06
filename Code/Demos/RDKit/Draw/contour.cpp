// $Id$
//
//  Copyright (C) 2014 Greg Landrum
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#define USE_CAIRO 1

#include <GraphMol/RDKitBase.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/Depictor/RDDepictor.h>
#include <GraphMol/MolDrawing/MolDrawing.h>
#include <GraphMol/MolDrawing/DrawingToSVG.h>
#ifdef USE_CAIRO
#include <GraphMol/MolDrawing/DrawingToCairo.h>
#endif

#include <Geometry/point.h>

#include <RDGeneral/RDLog.h>
#include <vector>
#include <algorithm>
#include <sstream>
#include <fstream>
#include <boost/foreach.hpp>
#include <boost/shared_ptr.hpp>

#include "conrec.h"
using namespace RDKit;
using namespace RDGeom;

struct ContourAccumulator {
  void operator()(double x1,double y1,double x2,double y2,double contour){
    contours[contour].push_back(std::make_pair(Point2D(x1,y1),Point2D(x2,y2)));
  };
  std::map<double,std::vector<std::pair<Point2D,Point2D> >  > contours;
};

class Grid {
public:
  unsigned int nX,nY;
  double **grid,*xCoords,*yCoords;
  Grid() : nX(0), nY(0), grid(NULL), xCoords(NULL), yCoords(NULL) {};
  Grid(unsigned int inX,unsigned inY) : nX(inX), nY(inY) {
    grid = new double * [nX];
    for(unsigned int i=0;i<nX;++i) grid[i]=new double [nY];
    xCoords = new double [nX];
    yCoords = new double [nY];
  }
  ~Grid(){
    if(grid){
      for(unsigned int i=0;i<nX;++i){
        if(grid[i]) delete [] grid[i];
      }
      delete [] grid;
    }
    if(xCoords) delete[] xCoords;
    if(yCoords) delete[] yCoords;
  }
private:
  Grid(const Grid &o);
};

void genMolGrid(const ROMol &mol,
                Point2D &minV,Point2D &maxV,
                Grid &grid,
                double &minGridVal,double &maxGridVal,
                double sigma=0.5){
  Point2D step=maxV-minV;
  step.x/=(grid.nX-1);
  step.y/=(grid.nY-1);    
  for(unsigned int i=0;i<grid.nX;++i){
    grid.xCoords[i] = minV.x + step.x*i;
  }
  for(unsigned int i=0;i<grid.nY;++i){
    grid.yCoords[i] = minV.y + step.y*i;
  }
  maxGridVal=-1000;
  minGridVal=1000;
  double denom = 2*sigma*sigma;
  double norm = 1./(sigma*sqrt(2*M_PI));
  for(unsigned int i=0;i<grid.nX;++i){
    for(unsigned int j=0;j<grid.nY;++j){
      Point3D pt(grid.xCoords[i],grid.yCoords[j],0);
      grid.grid[i][j] = 0.0;
      ROMol::VERTEX_ITER bAts,eAts;
      boost::tie(bAts,eAts)=mol.getVertices();
      while(bAts!=eAts){
        const Point3D &apt=mol.getConformer().getAtomPos(mol[*bAts]->getIdx());
        ++bAts;
        double d2=(apt-pt).lengthSq();
        double v = norm*exp(-d2/denom);
        grid.grid[i][j]+=v;
      }
      maxGridVal = max(maxGridVal,grid.grid[i][j]);
      minGridVal = min(minGridVal,grid.grid[i][j]);
      //std::cerr<<"      "<<pt.x<<" "<<pt.y<<" | "<<grid[i][j]<<std::endl;
    }
  }
}

void blah(const ROMol &mol,std::vector<int> &drawing,const Point2D &minV,const Point2D &maxV){
  static int nX=200,nY=200;
  Grid grid(nX,nY);
  Point2D gMinV(minV.x-2,minV.y-2);
  Point2D gMaxV(maxV.x+2,maxV.y+2);
  double maxGridVal,minGridVal;
  genMolGrid(mol,gMinV,gMaxV,grid,minGridVal,maxGridVal);
  std::cerr<<" max: "<<maxGridVal<<" min: "<<minGridVal<<std::endl;

  const int ncontours=5;
  double contours[ncontours];
  for(unsigned int i=0;i<ncontours;++i){
    contours[i] = minGridVal+(i+1)*(maxGridVal-minGridVal)/(ncontours+1);
  }
  ContourAccumulator accum;
  conrec(grid.grid,0,nX-1,0,nY-1,grid.xCoords,grid.yCoords,ncontours,contours,accum);

  int resolution=drawing[1];
  
  typedef std::pair<Point2D,Point2D> PtPair;
  for(unsigned int i=0;i<ncontours;++i){
    const std::vector<PtPair> &pv=accum.contours[contours[i]];
    BOOST_FOREACH(const PtPair &pt,pv){
      //std::cerr<<"   "<<pt.first<<" "<<pt.second;
      Point2D tpt1 = pt.first - minV;
      Point2D tpt2 = pt.second - minV;
      tpt1 *= resolution;
      tpt2 *= resolution;
      drawing.push_back(Drawing::PLINE);
      drawing.push_back(2);
      drawing.push_back(0);
      drawing.push_back(50);drawing.push_back(50);drawing.push_back(50);
      drawing.push_back(tpt1.x);drawing.push_back(tpt1.y);
      drawing.push_back(tpt2.x);drawing.push_back(tpt2.y);
      //std::cerr<<"   "<<tpt1<<" "<<tpt2<<std::endl;
    }
  }
}

std::string MolToSVG(ROMol &mol){
  if(!mol.getNumConformers()) {
    RDDepict::compute2DCoords(mol);
  }

  Point2D minV;
  Point2D maxV;
  std::vector<int> drawing=RDKit::Drawing::MolToDrawing(mol,NULL,true,&minV,&maxV);

  blah(mol,drawing,minV,maxV);

  std::string svg=RDKit::Drawing::DrawingToSVG(drawing);
  return svg;
}


#ifdef USE_CAIRO

void MolToCairo(ROMol &mol,cairo_t *cr,int width,int height){
  PRECONDITION(cr,"no context");
  PRECONDITION(width>0 && height>0,"bad dimensions");
  
  if(!mol.getNumConformers()) {
    RDDepict::compute2DCoords(mol);
  }

  Point2D minV;
  Point2D maxV;
  std::vector<int> drawing=RDKit::Drawing::MolToDrawing(mol,
                                                        NULL,true,
                                                        &minV,&maxV);
  blah(mol,drawing,minV,maxV);

  RDKit::Drawing::DrawingToCairo(drawing,cr,width,height);
#if 0
  double gridMin=0,gridMax=1.0;
  Point2D gstep=gMaxV-gMinV;
  gstep.x/=(nX-1);
  gstep.y/=(nY-1);
  gstep*=resolution;
  for(unsigned int i=0;i<nX;++i){
    for(unsigned int j=0;j<nY;++j){
      Point2D pti(xCoords[i],yCoords[j]);
      pti -= minV;
      pti *= resolution;
      double cv = (grid[i][j])/gridMax;
      cairo_set_source_rgba (cr, cv,cv,cv,0.5);
      cairo_rectangle(cr,pti.x,pti.y,gstep.x,gstep.y);
      cairo_fill(cr);
    }
  }
#endif

}
#endif

void DrawDemo(){

#if 1
  {
    RWMol *mol=SmilesToMol("Oc1c(C#N)cc(C(=O)NCc2sc([NH3+])c([NH3+])c2)cc1");
    std::string svg=MolToSVG(*mol);
    std::ofstream ostr("mol1.svg");
    ostr<<svg<<std::endl;
    delete mol;
  }
#endif

#ifdef USE_CAIRO
  {
    RWMol *mol=SmilesToMol("Oc1c(C#N)cc(C(=O)NCc2sc([NH3+])c([NH3+])c2)cc1");
    //RWMol *mol=SmilesToMol("Oc1cc(O)cc(O)c1");
    cairo_surface_t *surface =
      cairo_image_surface_create (CAIRO_FORMAT_ARGB32, 600, 600);
    cairo_t *cr = cairo_create (surface);
    MolToCairo(*mol,cr,600,600);

    cairo_destroy (cr);
    cairo_surface_write_to_png (surface, "mol1.png");
    cairo_surface_destroy (surface);
    delete mol;
  }
#endif
}

int
main(int argc, char *argv[])
{
  RDLog::InitLogs();
  DrawDemo();
}
