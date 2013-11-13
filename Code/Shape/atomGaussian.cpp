/*******************************************************************************
atomGaussian.cpp - Shape-it
 
Copyright 2012 by Silicos-it, a division of Imacosi BVBA
 
This file is part of Shape-it.

	Shape-it is free software: you can redistribute it and/or modify
	it under the terms of the GNU Lesser General Public License as published 
	by the Free Software Foundation, either version 3 of the License, or
	(at your option) any later version.

	Shape-it is distributed in the hope that it will be useful,
	but WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
	GNU Lesser General Public License for more details.

	You should have received a copy of the GNU Lesser General Public License
	along with Shape-it.  If not, see <http://www.gnu.org/licenses/>.

Shape-it is linked against OpenBabel version 2.

	OpenBabel is free software; you can redistribute it and/or modify
	it under the terms of the GNU General Public License as published by
	the Free Software Foundation version 2 of the License.

***********************************************************************/



#include <Shape/atomGaussian.h>



AtomGaussian::AtomGaussian(void)
:center(0.0, 0.0, 0.0)
, alpha(0.0)
, volume(0.0)
, C(0.0)
, nbr(0)
{
}



AtomGaussian::~AtomGaussian(void)
{
}



AtomGaussian atomIntersection(AtomGaussian & a, AtomGaussian & b)
{
    AtomGaussian c;

    // new alpha 
    c.alpha = a.alpha + b.alpha;

    // new center
    c.center.x = (a.alpha * a.center.x + b.alpha * b.center.x) / c.alpha;
    c.center.y = (a.alpha * a.center.y + b.alpha * b.center.y) / c.alpha;
    c.center.z = (a.alpha * a.center.z + b.alpha * b.center.z) / c.alpha;

    // self-volume 
    double d = (a.center.x - b.center.x) * (a.center.x - b.center.x)
	+ (a.center.y - b.center.y) * (a.center.y - b.center.y)
	+ (a.center.z - b.center.z) * (a.center.z - b.center.z);

    c.C = a.C * b.C * exp(-a.alpha * b.alpha / c.alpha * d);

    double scale = PI / (c.alpha);

    c.volume = c.C * scale * sqrt(scale);

    // set the number of gaussians 
    c.nbr = a.nbr + b.nbr;

    return c;
}
