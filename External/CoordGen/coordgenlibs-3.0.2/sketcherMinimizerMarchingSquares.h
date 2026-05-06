/*
 *  sketcherMarchingSquares.h
 *
 *
 *  Created by Nicola Zonta on 19/11/2010.
 *  Copyright Schrodinger, LLC. All rights reserved
 *
 */

#ifndef sketcherMINIMIZERMARCHINGSQUARES_H
#define sketcherMINIMIZERMARCHINGSQUARES_H

#include "CoordgenConfig.hpp"
#include <cstddef>
#include <vector>

class sketcherMinimizerPointF;

struct sketcherMinimizerMarchingSquaresPoint;

struct sketcherMinimizerMarchingSquaresSide {
    sketcherMinimizerMarchingSquaresPoint* p1;
    sketcherMinimizerMarchingSquaresPoint* p2;
};

struct sketcherMinimizerMarchingSquaresPoint {
  public:
    sketcherMinimizerMarchingSquaresPoint(float ix, float iy)
    {
        x = ix;
        y = iy;
        side1 = nullptr;
        side2 = nullptr;
        visited = false;
    }
    float x, y;
    sketcherMinimizerMarchingSquaresSide *side1, *side2;
    bool visited;
};

/* implementation of a marching squares algorithm */
class EXPORT_COORDGEN sketcherMinimizerMarchingSquares
{
  public:
    sketcherMinimizerMarchingSquares();
    ~sketcherMinimizerMarchingSquares();
    //  inline void clearGrid ();
    void setValue(float v, unsigned int x, unsigned int y);
    void initialize(float minx, float maxx, float miny, float maxy,
                    float x_interval, float y_interval = 0.f);

    void clear();

    void setThreshold(float t);
    float getThreshold() const;

    float toRealx(float x) const;
    float toRealy(float y) const;

    unsigned int getXN() const { return m_XN; };
    unsigned int getYN() const { return m_YN; };

    void run(); // computes the isovalue points and segments

    /* call after run () is executed, returns the coordinates of all the
       isovalue line points [x1, y1, x2, y2 .. xn, yn] in the order they
       were created */
    std::vector<float> getCoordinatesPoints() const;

    /* call after run () is executed. Returns a vector of isovalue closed
       lines [x1, y1, x2, y2 .. xn, yn]. The points are ordered as they
       appear along the line. */
    std::vector<std::vector<float>> getOrderedCoordinatesPoints() const;

    inline std::vector<float> getRawData() const
    {
        return m_grid;
    }; // returns a vector of all the data set with setValue.

    float getNodeValue(unsigned int x, unsigned int y) const;

  private:
    void addSide(sketcherMinimizerMarchingSquaresPoint* p1,
                 sketcherMinimizerMarchingSquaresPoint* p2);
    float m_xinterval, m_yinterval, m_left, m_bottom;
    std::vector<float> m_grid;
    unsigned int m_XN, m_YN;
    float m_threshold;
    std::vector<sketcherMinimizerMarchingSquaresPoint*> m_lastRowPoints;
    sketcherMinimizerMarchingSquaresPoint* m_lastCellRightPoint;
    float interpolate(float v1, float v2) const;
    std::vector<sketcherMinimizerMarchingSquaresPoint*> m_points;
    std::vector<sketcherMinimizerMarchingSquaresSide*> m_sides;
};

#endif // sketcherMINIMIZERMARCHINGSQUARES_H
