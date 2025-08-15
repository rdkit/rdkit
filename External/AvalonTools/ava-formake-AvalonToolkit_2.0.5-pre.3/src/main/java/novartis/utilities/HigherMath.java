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
package novartis.utilities;

/**
 * Provides math utilities exceeding the tools in java.lang.Math.
 */
public class HigherMath
{
   public static final double sqr(double x)
   {
      return (x*x);
   }



   public static double Distance(double x1, double y1, 
			         double x2, double y2)
   /* 
    * Returns the distance between the two space points (x1,y1) and (x2,y2)	
    */

  { 
     return (Distance(x1,y1,0,x2,y2,0));
  }


   public static double Distance(double x1, double y1, double z1,
			         double x2, double y2, double z2)
   /* 
    * Returns the distance between the two space points (x1,y1,z1) and (x2,y2,z2)	
    */

  {  double result;
     result = Math.sqrt(sqr(x1-x2) + sqr(y1-y2) + sqr(z1-z2));
     return (result);
  }


   public static double Angle(double x1, double y1, double x2, double y2)
   /*
    * Returns the angle between the two vectors (x1,y1,0) and (x2,y2,0) at
    * (0,0).
    */
    {return Angle(x1,y1,0,x2,y2,0);}   

    public static double Angle(double x1, double y1, double z1, double x2, double y2, double z2)
   /*
    * Returns the angle between the two vectors (x1,y1,z1) and (x2,y2,z2) at
    * (0,0).
    */
   {
      double l1, l2;
      double cos_alpha, sin_alpha;
      double X_diff, Y_diff, Z_diff; 
      double result;
      l1 = Math.sqrt(x1*x1+y1*y1+z1*z1); l2 = Math.sqrt(x2*x2+y2*y2+z2*z2);
      if (l1 < 0.00001  ||  l2 < 0.00001) return (0.0);
      cos_alpha = (x1*x2 + y1*y2 + z1*z2)/(l1*l2);
      if (cos_alpha > 1.0)          /* safeguard against round off erros */
         cos_alpha = 1.0;
      X_diff = y1*z2 - z1*y2;
      Y_diff = z1*x2 - z2*x1;
      Z_diff = x1*y2 - x2*y1;
      sin_alpha = (X_diff + Y_diff + Z_diff)/(l1*l2);
      result = Math.acos(cos_alpha);
      if (sin_alpha < 0.0) result = 2*Math.PI-result;
      return (result);
   }
   
   
   
     public static double AngleSmall(double x1, double y1, double x2, double y2)
   /*
    * Returns the angle between the two vectors (x1,y1,0) and (x2,y2,0) at
    * (0,0).
    */
    {return AngleSmall(x1,y1,0,x2,y2,0);}   

    public static double AngleSmall(double x1, double y1, double z1, double x2, double y2, double z2)
   /*
    * Returns the angle between the two vectors (x1,y1,z1) and (x2,y2,z2) at
    * (0,0).
    */
   {
      double l1, l2;
      double cos_alpha, sin_alpha;
      double X_diff, Y_diff, Z_diff; 
      double result;
      l1 = Math.sqrt(x1*x1+y1*y1+z1*z1); l2 = Math.sqrt(x2*x2+y2*y2+z2*z2);
      if (l1 < 0.00001  ||  l2 < 0.00001) return (0.0);
      cos_alpha = (x1*x2 + y1*y2 + z1*z2)/(l1*l2);
      if (cos_alpha > 1.0)          /* safeguard against round off erros */
         cos_alpha = 1.0;
      X_diff = y1*z2 - z1*y2;
      Y_diff = z1*x2 - z2*x1;
      Z_diff = x1*y2 - x2*y1;
      result = Math.acos(cos_alpha);
      if (result > Math.PI) result = 2*Math.PI-result;
      if (result < -Math.PI) result = 2*Math.PI + result;	
      return (result);
   }
   
   
   
   

    public static double Angle_3p(double x1, double y1, double z1,
				  double x2, double y2, double z2,
                                  double x3, double y3, double z3)	
   /* returns the angle 123 */
      {
        double xa, ya, za, xb, yb, zb; 
        double result;    
        xa = x1 - x2; 
 	ya = y1 - y2;
 	za = z1 - z2;
 
	xb = x3 - x2;
 	yb = y3 - y2;
 	zb = z3 - z2; 
        result = HigherMath.AngleSmall(xa, ya, za, xb, yb, zb);
        return (result);
       } 
          

      public static double[] Normale(double x1, double y1, double z1,
 				     double x2, double y2, double z2,
 				     double x3, double y3, double z3) 
 	/*
   	* Returns the perpendicular vector to the plane defined 
    	* by points 1,2,3 as a 3-element array
    	*/

	 {
 	double nx, ny, nz;        
 	double length;
 	double xa, ya, za, xb, yb, zb;
 	double norm_nx, norm_ny, norm_nz;
 
 	xa = x1 - x2; 
 	ya = y1 - y2;
 	za = z1 - z2;
 
	xb = x3 - x2;
 	yb = y3 - y2;
 	zb = z3 - z2;
	
 	nx = ya*zb - yb*za;	    // "Kreuzprodukt"
 	ny = za*xb - zb*xa;
 	nz = xa*yb - xb*ya;
	  
 	length = Math.sqrt(nx*nx + ny*ny + nz*nz);
 
 	norm_nx = nx/length;
 	norm_ny = ny/length;
 	norm_nz = nz/length;
	
	double[] N_vektor = new double[3];
 	N_vektor[0] = norm_nx;
 	N_vektor[1] = norm_ny;
 	N_vektor[2] = norm_nz;
	return (N_vektor);
 	}



	public static double Torsion_4p (double x1, double y1, double z1,
	        			double x2, double y2, double z2,
				       	double x3, double y3, double z3,
					double x4, double y4, double z4)	

	 /*
    	  * returns the torsion angle defined by four points in 3D-space
   	  * 1,2,3,4; absolute value is given 
   	  */

	{     
	double[] normale_123;
	double[] normale_234; 
	double TOR;
	normale_123 = HigherMath.Normale(x1,y1,z1,x2,y2,z2,x3,y3,z3);
	normale_234 = HigherMath.Normale(x2,y2,z2,x3,y3,z3,x4,y4,z4);
	TOR = HigherMath.AngleSmall(normale_123[0], normale_123[1], normale_123[2],
			               normale_234[0], normale_234[1], normale_234[2]);
	return (TOR);
	}




   /**
    * Transforms the 2D points in r[0...n-1][2] such that the points p1[] and
    * p2[] would be transformed to p1p[] and p2p[], resp. in a best fit
    * manner by just rotating and moving them without scaling.
    */
   public static void transformPoints(double r[][], int n,
                                      double p1[],  double p2[],
                                      double p1p[], double p2p[])
   {
      double rp[] = new double[2];       /* r prime                              */
      double dp[] = new double[2];       /* delta p                              */
      double dpp[] = new double[2];      /* delta p prime                        */
      double mp[] = new double[2];       /* mean p                               */
      double mpp[] = new double[2];      /* mean p prime                         */
      double a, b, c;                    /* parameters of rotation matrix        */
      double ls, lsp;                    /* distances between source points and  */
                                         /* destination points, resp.            */
      int i, j;

      ls = lsp = 0.0;
      for (i=0; i<2; i++)          /* compute transformation parameters    */
      {
         dp[i] = p1[i]-p2[i];         dpp[i] = p1p[i]-p2p[i];
         mp[i] = (p1[i]+p2[i])/2.0;   mpp[i] = (p1p[i]+p2p[i])/2.0;
         ls += dp[i]*dp[i];           lsp += dpp[i]*dpp[i];
      }
      if (ls < 0.00001  ||  lsp < 0.00001)
      {
         a = 1.0; b = 0.0; c = 1.0;
      }
      else
      {
         a = (dp[0]*dpp[0] + dp[1]*dpp[1])/ls;
         b = (dpp[0]*dp[1] - dp[0]*dpp[1])/ls;
         c = Math.sqrt(lsp/ls);
      }

      for (j=0; j<n; j++)              /* transform points                     */
      {
         for (i=0; i<2; i++) /* move center of gravity */
            rp[i] = r[j][i]-mp[i];

         r[j][0] =  a*rp[0] + b*rp[1];  /* rotate point */
         r[j][1] = -b*rp[0] + a*rp[1];
         for (i=0; i<2; i++) rp[i] = r[j][i]/c;

         for (i=0; i<2; i++)      /* restore new center of gravity */
         r[j][i] = rp[i]+mpp[i];
      }
   }

   /**
    * Finds the next best point to substitute the graph defined by
    * coords[][2] and edges[][2] at node index seed.
    * The 2D coordinates found are returned.
    *
    * The parameter use_inwards tells whether candidated positions may
    * also be put between existing substituents.
    */
   public static double[] nextSubstituentPoint(double coordinates[][],
                                               int edges[][],
                                               int seed,
			                                      boolean use_inwards)
   {
      double point[] = new double[2];

      int i, j;
      double value, dist, len, vtmp;
      double m[] = new double[2];                  /* coordinates of seed */
      double neigh[][]; /* coordinates of neighbours of seed */
      int nneigh;
      double cand_points[][];       /* coordinates of candidates for point */
      int ncand;
      double t[] = new double[2], th[] = new double[2];            /* temporary points */

      int    nphi=2;  /* directions of sprouting (120 degrees, 90 degrees are ignored) */
      double sinphi[] = {0.866025, -0.866025, 0.707107, -0.707107};
      double cosphi[] = {-0.5,     -0.5,      0.707107,  0.707107};

                                       /* get # of neighbours of seed */
      m[0] = coordinates[seed][0];
      m[1] = coordinates[seed][1];
      for (i=0, nneigh=0; i<edges.length; i++)
         if (edges[i][0] == seed  ||  edges[i][1] == seed) nneigh++;

      if (nneigh == 0)                        /* trivial boundary cases */
      {
         cand_points = new double[4*nphi+2][];
         ncand=0;
         th[0] = 1.514;
         th[1] = 0.0;
         cand_points[ncand] = new double[2];
         cand_points[ncand][0] = th[0] + m[0];
         cand_points[ncand][1] = th[1] + m[1];
         ncand++;
         cand_points[ncand] = new double[2];
         cand_points[ncand][0] = -th[0] + m[0];
         cand_points[ncand][1] = -th[1] + m[1];
         ncand++;

         for (i=0; i<2*nphi; i++)
         {
            cand_points[ncand] = new double[2];
            cand_points[ncand][0] = m[0] + th[0]*cosphi[i] - th[1]*sinphi[i];
	         cand_points[ncand][1] = m[1] + th[0]*sinphi[i] + th[1]*cosphi[i];
            ncand++;
         }

         for (i=0; i<2*nphi; i++)
         {
            cand_points[ncand] = new double[2];
            cand_points[ncand][0] = m[0] - th[0]*cosphi[i] + th[1]*sinphi[i];
	         cand_points[ncand][1] = m[1] - th[0]*sinphi[i] - th[1]*cosphi[i];
            ncand++;
         }
      }
      else
      {
         neigh = new double[nneigh][]; /* set-up neighbourhood */
         for (i=0, nneigh=0; i<edges.length; i++)
         {
           if (edges[i][0] == seed)
           {
              neigh[nneigh] = new double[2];
              neigh[nneigh][0] = coordinates[edges[i][1]][0];
              neigh[nneigh][1] = coordinates[edges[i][1]][1];
              nneigh++;
           }
           if (edges[i][1] == seed)
           {
              neigh[nneigh] = new double[2];
              neigh[nneigh][0] = coordinates[edges[i][0]][0];
              neigh[nneigh][1] = coordinates[edges[i][0]][1];
              nneigh++;
           }
         }

                                                  /* set-up candidates */
         cand_points = new double[nphi+(nneigh-1)*nneigh][];
         ncand=0;

         for (i=0, len=0.0; i<nneigh; i++)   /* compute average edge length */
         {
            len += distSquare(m[0],m[1], neigh[i][0], neigh[i][1]);
         }


         len = Math.sqrt(len/nneigh);          /* add angle candidates */
         dist = distSquare(m[0],m[1],neigh[0][0],neigh[0][1]);
         dist = Math.sqrt(dist);
         if (dist > 0.001)
         {
            th[0] = (neigh[0][0]-m[0])*len/dist;
            th[1] = (neigh[0][1]-m[1])*len/dist;
            for (i=0; i<nphi; i++)
            {
               cand_points[ncand] = new double[2];
               cand_points[ncand][0] = th[0]*cosphi[i] - th[1]*sinphi[i] + m[0];
               cand_points[ncand][1] = th[0]*sinphi[i] + th[1]*cosphi[i] + m[1];
               ncand++;
            }
         }

         for (i=0; i<nneigh; i++)   /* add vector sum and diff. candidates */
            for (j=i+1; j<nneigh; j++)
            {
               t[0] = neigh[i][0]-m[0]; t[1] = neigh[i][1]-m[1];
               dist = Math.sqrt(t[0]*t[0]+t[1]*t[1]);
               if (dist > 0.001) {t[0] /= dist; t[1] /= dist;}
               th[0] = neigh[j][0]-m[0]; th[1] = neigh[j][1]-m[1];
               dist = Math.sqrt(th[0]*th[0]+th[1]*th[1]);
               if (dist > 0.001) {th[0] /= dist; th[1] /= dist;}
               t[0] = t[0]+th[0];
               t[1] = t[1]+th[1];
               dist = Math.sqrt(t[0]*t[0]+t[1]*t[1]);
	            if (dist > 0.001)
               {
	               /* pointing outwards */
                  cand_points[ncand] = new double[2];
                  cand_points[ncand][0] = m[0]-t[0]*len/dist;
                  cand_points[ncand][1] = m[1]-t[1]*len/dist;
                  ncand++;
	               /* pointing inwards */
	               if (use_inwards && nneigh>2)
	               {
                     cand_points[ncand] = new double[2];
		               cand_points[ncand][0] = m[0]+t[0]*len/dist;
		               cand_points[ncand][1] = m[1]+t[1]*len/dist;
		               ncand++;
	               }
               }
            }
      }

      value = 1.0e7;             /* find best candidate */
      for (i=0; i<ncand; i++)
      {
         for (j=0, vtmp=0.0; j<coordinates.length; j++)
            vtmp +=
               1.0/(0.001+distSquare(cand_points[i][0],cand_points[i][1],
                               coordinates[j][0],coordinates[j][1]));

         if (i < nphi) vtmp *= 1.5;        /* prefer bisection candidates */

         if (vtmp < value)
         {
            value = vtmp;
            point[0] = cand_points[i][0];
            point[1] = cand_points[i][1];
         }
      }

      return (point);
   }

   /**
    * Computes the distance square of x1/y1 and x2/y2.
    */
   public static double distSquare(double x1, double y1, double x2, double y2)
   {
      return (((x1)-(x2))*((x1)-(x2)) + ((y1)-(y2))*((y1)-(y2)));
   }
}
