#ifndef FAUNUS_SPHEROCYLINDER_H
#define FAUNUS_SPHEROCYLINDER_H

#include <faunus/geometry.h>
#include <faunus/potentials.h>

namespace Faunus {

  namespace Geometry {

    /**
     * @brief Calculate minimum distance between two line segments
     *
     * Find closest distance between line segments and return its vector
     * gets orientations and lengths of line segments and the vector connecting
     * their center os masses (from vec1 to vec2)
     * Copyright 2001, softSurfer (www.softsurfer.com)
     * This code may be freely used and modified for any purpose
     * providing that this copyright notice is included with it.
     * SoftSurfer makes no warranty for this code, and cannot be held
     * liable for any real or imagined damage resulting from its use.
     * Users of this code must verify correctness for their application.
     *
     * \param dir1 Direction of first segment
     * \param halfl1 Half length of first segment
     * \param dir2 Direction of second segment
     * \param halfl2 Half length of second segment
     * \param r_cm Distance vector between the middle of the two segments
     */
    template<class T=double>
      Point mindist_segment2segment(const Point &dir1, T halfl1,
          const Point &dir2, T halfl2, const Point &r_cm)
      {
        Point u = dir1 * (halfl1*2); //S1.P1 - S1.P0;
        Point v = dir2 * (halfl2*2); //S2.P1 - S2.P0;
        Point w = dir2*halfl2 - dir1*halfl1 - r_cm; //S1.P0-S2.P0
        T a = u.dot(u);        // always >= 0
        T b = u.dot(v);
        T c = v.dot(v);        // always >= 0
        T d = u.dot(w);
        T e = v.dot(w);
        T D = a*c - b*b;       // always >= 0
        T sc = D;
        T sN = D;
        T sD = D;      // sc = sN / sD, default sD = D >= 0
        T tc = D;
        T tN = D;
        T tD = D;      // tc = tN / tD, default tD = D >= 0

        // compute the line parameters of the two closest points
        if (D < 0.00000001) { // the lines are almost parallel
          sN = 0.0;        // force using point P0 on segment S1
          sD = 1.0;        // to prevent possible division by 0.0 later
          tN = e;
          tD = c;
        }
        else {                // get the closest points on the infinite lines
          sN = (b*e - c*d);
          tN = (a*e - b*d);
          if (sN < 0.0) {       // sc < 0 => the s=0 edge is visible
            sN = 0.0;
            tN = e;
            tD = c;
          }
          else if (sN > sD) {  // sc > 1 => the s=1 edge is visible
            sN = sD;
            tN = e + b;
            tD = c;
          }
        }

        if (tN < 0.0) {           // tc < 0 => the t=0 edge is visible
          tN = 0.0;
          // recompute sc for this edge
          if (-d < 0.0)
            sN = 0.0;
          else if (-d > a)
            sN = sD;
          else {
            sN = -d;
            sD = a;
          }
        }
        else if (tN > tD) {      // tc > 1 => the t=1 edge is visible
          tN = tD;
          // recompute sc for this edge
          if ((-d + b) < 0.0)
            sN = 0;
          else if ((-d + b) > a)
            sN = sD;
          else {
            sN = (-d + b);
            sD = a;
          }
        }
        // finally do the division to get sc and tc
        if (std::abs(sN) < 0.00000001)
          sc = 0.0 ;
        else
          sc = sN / sD;
        if (std::abs(tN) < 0.00000001)
          tc = 0.0 ;
        else
          tc = tN / tD;

        // get the difference of the two closest points
        //Vector = w + (sc * u) - (tc * v);  // = S1(sc) - S2(tc)
        //vec.x = u.x*sc + w.x - v.x*tc;
        //vec.y = u.y*sc + w.y - v.y*tc;
        //vec.z = u.z*sc + w.z - v.z*tc;
        return u*sc + w - v*tc;
      }    

    /*!
     * \param dir Direction of segment
     * \param halfl Half length of segment
     * \param r_cm Distance vector between the middle segment to point
     */
    template<class T=double>
      Point mindist_segment2point(const Point &dir, T halfl, const Point &r_cm) {
        T d, c = dir.dot(r_cm);
        if (c > halfl)
          d = halfl;
        else {
          if (c > -halfl)
            d = c;
          else
            d = -halfl;
        }
        return -r_cm + (dir*d);
      }

    template<class Tpoint1, class Tpoint2>
      Point vec_perpproject(const Tpoint1 &a, const Tpoint2 &b) {
        return a-b*a.dot(b);
      }

    /**
     * Finds intersections of spherocylinder and plane defined by vector
     * "w_vec" and if they are in all-way patch then returns number of them (PSC)
     */
    template<class Tpsc>
      int find_intersect_plane(
          const Tpsc &part1,
          const Tpsc &part2,  
          const Point &r_cm,
          const Point &w_vec,
          double rcut,
          double cospatch,
          double intersections[5] )
      {
        int intrs;
        double disti;

        Point nplane=part1.dir.cross(w_vec);
        nplane.normalize();
        double a =  nplane.dot(part2.dir);
        if (a == 0.0)
          intrs=0; /* there is no intersection plane and sc are paralel*/
        else {
          double ti = nplane.dot(r_cm)/a;
          if ((ti  > part2.halfl ) || (ti < -part2.halfl))
            intrs=0; /* there is no intersection plane sc is too short*/
          else {
            Point d_vec = ti * part2.dir - r_cm; /*vector from intersection point to CM*/
            double c = d_vec.dot(w_vec);
            if ( c * cospatch < 0)
              intrs=0; /* the intersection in plane is on other side of patch */
            else {
              double d = fabs(d_vec.dot(part1.dir)) - part2.halfl;
              if (d <= 0)
                disti = c*c; /*is inside cylinder*/
              else
                disti = d*d + c*c; /*is inside patch*/
              if (disti > rcut*rcut)
                intrs=0; /* the intersection is outside sc */
              else {
                intrs=1;
                int i=0;
                while (intersections[i] !=0) {
                  if (ti == intersections[i])
                    intrs=0; /* found intersection we already have -it is at boundary*/
                  i++;
                }
                if (intrs > 0) intersections[i]=ti;
              }
            }
          }

        }
        return intrs;
      }

    /*!
     * @brief Finds if vector "vec" has angular intersection w. patch of part1 
     */
    template<class Tpsc, class T=double>
      int test_intrpatch(const Tpsc &part1, Point &vec, T cospatch, 
          T ti, T intersections[5])
      {
        /*test if we have intersection*/
        /* do projection to patch plane*/
        vec=vec_perpproject(vec,part1.dir);
        vec.normalize();
        /* test angle distance from patch*/
        T a = part1.patchdir.dot(vec);
        int intrs=0;
        if (a >= cospatch) {
          intrs=1;
          int i=0;
          while (intersections[i] !=0) {
            if (ti == intersections[i]) 
              intrs=0; /* found intersection we already have -it is at boundary*/
            i++;
          }
          if (intrs > 0) 
            intersections[i]=ti;
        }
        return intrs;
      }

    /**
     * @brief Intersect of plane
     *
     * Finds intersections of plane defined by vector "w_vec"
     * and if they are in cylindrical patch then returns number of them (CPSC)
     */
    template<class Tpsc>
      int find_intersect_planec(
          const Tpsc &part1, const Tpsc &part2, 
          const Point &r_cm, const Point &w_vec,
          double rcut, double cospatch, double intersections[5])
      {
        int i, intrs=0;
        double a, c, d, ti, disti;
        Point nplane=part1.dir.cross(w_vec);
        nplane.normalize();
        a =  nplane.dot( part2.dir);
        if (a == 0.0) intrs=0; /* there is no intersection plane and sc are paralel*/
        else {
          ti = nplane.dot(r_cm)/a;
          if ((ti  > part2.halfl ) || (ti < -part2.halfl))
            intrs=0; /* there is no intersection plane sc is too short*/
          else {
            Point d_vec = ti*part2.dir - r_cm; /*vector from intersection point to CM*/
            c = d_vec.dot( w_vec);
            if ( c *cospatch < 0)
              intrs=0; /* the intersection in plane is on other side of patch */
            else {
              d = fabs(d_vec.dot( part1.dir)) - part2.halfl;
              if (d <= 0) {
                disti= c*c; /*is inside cylinder*/
                if (disti > rcut*rcut)
                  intrs=0; /* the intersection is outside sc */
                else {
                  intrs=1;
                  i=0;
                  while (intersections[i] !=0) {
                    if (ti == intersections[i])
                      intrs=0; /* found intersection we already have -it is at boundary*/
                    i++;
                  }
                  if (intrs > 0)
                    intersections[i]=ti;
                }
              }
            }
          }
        }
        return intrs;
      }

    /*!
      \brief Calculates intersections of spherocylinder2 with a all-way patch of spherocylinder1 and return them (PSC)
      */
    template<class Tpsc, class T=double>
      int psc_intersect(
          const Tpsc &part1, const Tpsc &part2,
          const Point &r_cm, T intersections[5], T rcut)
      {
        T a, b, c, d, e, x1, x2;
        Point cm21, vec1, vec2, vec3, vec4;          

        T rcut2=rcut*rcut;
        /*1- do intersections of spherocylinder2 with patch of spherocylinder1 at 
          cut distance C*/
        /*1a- test intersection with half planes of patch and look how far they are 
          from spherocylinder. If closer then C  we got itersection*/

        /* plane1 */
        /* find intersections of part2 with plane by par1 and patchsides[0] */
        int intrs=0;
        intrs+=find_intersect_plane(part1,part2,r_cm,part1.patchsides[0],rcut,part1.pcanglsw,intersections);
        //	    printf("plane1 %d\n", intrs);
        /* plane2 */
        /* find intersections of part2 with plane by par1 and patchsides[1] */
        intrs+=find_intersect_plane(part1,part2,r_cm,part1.patchsides[1],rcut,part1.pcanglsw,intersections);


        if ( (intrs == 2 ) && (part1.pcanglsw <0) ) {
          assert("Patch>180 -> two segments.");
        }

        /*1b- test intersection with cylinder - it is at distance C*/
        if (intrs < 2 )  {
          cm21=-r_cm;
          vec1=cm21.cross(part1.dir);
          vec2=part2.dir.cross(part1.dir);
          a = vec2.dot(vec2);
          b = 2*vec1.dot(vec2);
          c = -rcut*rcut + vec1.dot(vec1);
          d = b*b - 4*a*c;
          if ( d >= 0) { /*there is intersection with infinite cylinder */
            x1 = (-b+sqrt(d))*0.5/a;/*parameter on line of SC2 determining intersection*/
            if ((x1 >=part2.halfl) || (x1 <= -part2.halfl))
              intrs+=0; /*intersection is outside sc2*/
            else {
              /* vectors from center os sc1 to intersection with infinite cylinder*/
              vec1=part2.dir*x1-r_cm;
              e = part1.dir.dot(vec1);
              if ((e >=part1.halfl) || (e <= -part1.halfl))
                intrs+=0; /*intersection is outside sc1*/
              else {
                intrs+=test_intrpatch(part1,vec1,part1.pcanglsw,x1,intersections);
              }
            }
            if ( d > 0 ){
              x2 = (-b-sqrt(d))*0.5/a;/*parameter on line of SC2 determining intersection*/
              if ((x2 >=part2.halfl) || (x2 <= -part2.halfl))
                intrs+=0; /*intersection is outside sc2*/
              else {
                vec2 = part2.dir*x2-r_cm;
                e = part1.dir.dot(vec2);
                if ((e >=part1.halfl) || (e <= -part1.halfl)) intrs+=0; /*intersection is outside sc1*/
                else {
                  intrs+=test_intrpatch(part1,vec2,part1.pcanglsw,x2,intersections);
                }
              }
            }
          }
        }
        //	    printf ("cylinder %d x1 %f x2 %f e %f\n", intrs, x1, x2, e);
        /*1c- test intersection with spheres at the end - it is at distace C*/
        if (intrs < 2 )  {
          /*centers of spheres*/
          /*relative to the CM of sc2*/
          vec1 =  part1.dir*part1.halfl - r_cm;
          vec2 = -part1.dir*part1.halfl - r_cm;

          /*sphere1*/
          a = part2.dir.dot(part2.dir);
          b = 2.0*vec1.dot(part2.dir);
          c = vec1.dot(vec1)-rcut*rcut;
          d = b*b-4*a*c;
          if (d >= 0) { /*if d<0 there are no intersections*/
            x1= (-b + sqrt(d))*0.5/a; /*parameter on line of SC2 determining intersection*/
            if ((x1 >=part2.halfl) || (x1 <= -part2.halfl)) intrs+=0; /*intersection is outside sc2*/
            else {
              vec3 = part2.dir*x1-r_cm;
              e = part1.dir.dot(vec3);
              if ((e >= part1.halfl) || (e <= -part1.halfl)) { /*if not intersection is inside sc1*/
                intrs+=test_intrpatch(part1,vec3,part1.pcanglsw,x1,intersections);
              }
            }
            if ( d > 0) {
              x2= (-b - sqrt(d))*0.5/a; /*parameter on line of SC2 determining intersection*/
              if ((x2 >=part2.halfl) || (x2 <= -part2.halfl)) intrs+=0; /*intersection is outside sc2*/
              else {
                vec4 = part2.dir*x2 - r_cm;
                e = part1.dir.dot(vec4);
                if ((e >=part1.halfl) || (e <= -part1.halfl)) { /*if not intersection is inside sc1*/
                  intrs+=test_intrpatch(part1,vec4,part1.pcanglsw,x2,intersections);
                }
              }
            }
          }
          //		printf ("sphere1 %d x1 %f x2 %f e %f\n", intrs, x1, x2, e);
          /*sphere2*/
          a = part2.dir.dot(part2.dir);
          b = 2.0*vec2.dot(part2.dir);
          c = vec2.dot(vec2)-rcut*rcut;
          d = b*b-4*a*c;
          if (d >= 0) { /*if d<0 there are no intersections*/
            x1= (-b + sqrt(d))*0.5/a; /*parameter on line of SC2 determining intersection*/
            if ((x1 >=part2.halfl) || (x1 <= -part2.halfl)) intrs+=0; /*intersection is outside sc2*/
            else {
              vec3 = part2.dir*x1 - r_cm;
              e = part1.dir.dot(vec3);
              if ((e >=part1.halfl) || (e <= -part1.halfl)) { /*if not intersection is inside sc1*/
                intrs+=test_intrpatch(part1,vec3,part1.pcanglsw,x1,intersections);
              }
            }
            if ( d > 0 ) {
              x2= (-b - sqrt(d))*0.5/a; /*parameter on line of SC2 determining intersection*/
              if ((x2 >=part2.halfl) || (x2 <= -part2.halfl)) intrs+=0; /*intersection is outside sc2*/
              else {
                vec4 = part2.dir*x2 - r_cm;
                e = part1.dir.dot(vec4);
                if ((e >=part1.halfl) || (e <= -part1.halfl)) { /*if not intersection is inside sc1*/
                  intrs+=test_intrpatch(part1,vec4,part1.pcanglsw,x2,intersections);
                }
              }
            }
          }
          //		printf ("sphere2 %d\n", intrs);
        }

        /*1d- if there is only one itersection shperocylinder ends within patch wedge 
          set as second intersection end inside patch*/
        if (intrs < 2 )  {
          /*whole spherocylinder is in or all out if intrs ==0*/
          vec1 = part2.dir*part2.halfl - r_cm;
          /*vector from CM of sc1 to end of sc2*/
          /*check is is inside sc1*/
          a=vec1.dot(part1.dir);
          vec3 = vec1 - part1.dir*a;
          b=vec3.dot(vec3);
          d = fabs(a)-part1.halfl;
          if ( d <= 0) 
            c = b; /*is inside cylindrical part*/
          else 
            c = d*d + b; /*is inside caps*/
          /*c is distance squared from line or end to test if is inside sc*/
          if (c < rcut2) 
            intrs+=test_intrpatch(part1,vec1,part1.pcanglsw,part2.halfl,intersections);
          if (intrs < 2 ) {
            vec2 = -part2.dir*part2.halfl - r_cm;
            /*check is is inside sc1*/
            a=vec2.dot(part1.dir);
            vec4 = vec2 - part1.dir*a;
            b=vec4.dot(vec4);
            d = fabs(a) -part1.halfl;
            if (d <= 0) 
              c = b; /*is inside cylindrical part*/
            else 
              c = d*d + b; /*is inside caps*/
            /*c is distance squared from line or end to test if is inside sc*/
            if (c < rcut2) 
              intrs+=test_intrpatch(part1,vec2,part1.pcanglsw,-1.0*part2.halfl,intersections);
          }
        }
        return intrs;
      }

    /**
      @brief Intersection of PSC2 with cylindrical patch of PSC1 and return them (CPSC)
      */
    template<class Tpsc>
      int cpsc_intersect(const Tpsc &part1, const Tpsc &part2,
          const Point &r_cm, double intersections[5], double rcut)
      {
        int intrs;
        double a, b, c, d, e, x1, x2, rcut2;
        Point cm21, vec1, vec2, vec3, vec4;

        intrs=0;
        rcut2=rcut*rcut;
        /*1- do intersections of spherocylinder2 with patch of spherocylinder1 at 
          cut distance C*/
        /*1a- test intersection with half planes of patch and look how far they are 
          from spherocylinder. If closer then C  we got itersection*/

        /* plane1 */
        /* find intersections of part2 with plane by par1 and part1.patchsides[0] */
        intrs+=find_intersect_planec(part1,part2,r_cm,part1.patchsides[0],rcut,part1.pcanglsw,intersections);
        //	    printf("plane1 %d\n", intrs);
        /* plane2 */
        /* find intersections of part2 with plane by par1 and part1.patchsides[1] */
        intrs+=find_intersect_planec(part1,part2,r_cm,part1.patchsides[1],rcut,part1.pcanglsw,intersections);

        if ( (intrs == 2 ) && (part1.pcanglsw < 0) ) {
          fprintf (stderr, "ERROR: Patch is larger than 180 degrees and we are getting two segments - this hasnot been programed yet.\n\n");
          exit (1);
        }
        //	    printf("plane2 %d\n", intrs);

        /*1b- test intersection with cylinder - it is at distance C*/
        if (intrs < 2 )  {
          cm21=-r_cm;
          vec1=cm21.cross(part1.dir);
          vec2=part2.dir.cross(part1.dir);
          a = vec2.dot(vec2);
          b = 2*vec1.dot(vec2);
          c = -rcut*rcut + vec1.dot(vec1);
          d = b*b - 4*a*c;
          if ( d >= 0) { /*there is intersection with infinite cylinder */
            x1 = (-b+sqrt(d))*0.5/a; /*parameter on line of SC2 determining intersection*/
            if ((x1 >=part2.halfl) || (x1 <= -part2.halfl)) intrs+=0; /*intersection is outside sc2*/
            else {
              /* vectors from center os sc1 to intersection with infinite cylinder*/
              vec1=part2.dir*x1-r_cm;
              e = part1.dir.dot(vec1);
              if ((e >=part1.halfl) || (e <= -part1.halfl)) intrs+=0; /*intersection is outside sc1*/
              else {
                intrs+=test_intrpatch(part1,vec1,part1.pcanglsw,x1,intersections);
              }
            }
            if ( d > 0 ){
              x2 = (-b-sqrt(d))*0.5/a; /*parameter on line of SC2 determining intersection*/
              if ((x2 >=part2.halfl) || (x2 <= -part2.halfl)) intrs+=0; /*intersection is outside sc2*/
              else {
                vec2 = part2.dir*x2-r_cm;
                e = part1.dir.dot(vec2);
                if ((e >=part1.halfl) || (e <= -part1.halfl)) intrs+=0; /*intersection is outside sc1*/
                else {
                  intrs+=test_intrpatch(part1,vec2,part1.pcanglsw,x2,intersections);
                }
              }
            }
          }
        }
        //	    printf ("cylinder %d x1 %f x2 %f e %f\n", intrs, x1, x2, e);
        /*1c- test intersection with plates at the end - it is at distace C and in wedge*/
        if (intrs < 2 )  {
          a =  part1.dir.dot( part2.dir);
          if (a == 0.0) intrs=0; /* there is no intersection plane and sc are paralel*/
          else {
            /*plane cap1*/
            vec1= r_cm + part1.halfl*part1.dir;
            x1 = part1.dir.dot(vec1)/a; /*parameter on line of SC2 determining intersection*/
            if ((x1 > part2.halfl ) || (x1 < -part2.halfl)) intrs+=0; /* there is no intersection plane sc is too short*/
            else {
              vec2 = x1*part2.dir - vec1; /*vector from ENDPOINT to intersection point */
              b = vec2.dot( vec2);
              if (b > rcut*rcut) intrs+=0;   /* the intersection is outside sc */
              else {
                intrs+=test_intrpatch(part1,vec2,part1.pcanglsw,x1,intersections);
              }
            }
            //		    printf ("plane cap1 %d %f\n", intrs, x1);
            /*plane cap2*/
            vec1= r_cm - part1.halfl*part1.dir;
            x2 = part1.dir.dot(vec1)/a; /*parameter on line of SC2 determining intersection*/
            if ((x2  > part2.halfl ) || (x2 < -part2.halfl)) intrs+=0; /* there is no intersection plane sc is too short*/
            else {
              vec2 = x2*part2.dir - vec1; /*vector from ENDPOINT to intersection point */
              b = vec2.dot( vec2);
              if (b > rcut*rcut) intrs+=0;   /* the intersection is outside sc */
              else {
                intrs+=test_intrpatch(part1,vec2,part1.pcanglsw,x2,intersections);
              }
            }
            //		    printf ("plane cap2 %d %f\n", intrs,x2);

          }
        }

        /*1d- if there is only one itersection shperocylinder ends within patch wedge 
          set as second intersection end inside patch*/
        if (intrs < 2 )  {
          /*whole spherocylinder is in or all out if intrs ==0*/
          vec1 = part2.dir*part2.halfl - r_cm;
          /*vector from CM of sc1 to end of sc2*/
          /*check is is inside sc1*/
          a=vec1.dot(part1.dir);
          vec3 = vec1 - part1.dir*a;
          b=vec3.dot(vec3);
          d = fabs(a)-part1.halfl;
          if ( d <= 0) { /*is in cylindrical part*/
            /*c is distance squared from line or end to test if is inside sc*/
            if (b < rcut2) intrs+=test_intrpatch(part1,vec1,part1.pcanglsw,part2.halfl,intersections);
          }
          if (intrs < 2 ) {
            vec2 = -part2.dir*part2.halfl - r_cm;
            /*check is is inside sc1*/
            a=vec2.dot(part1.dir);
            vec4 = vec2 - part1.dir*a;
            b=vec4.dot(vec4);
            d = fabs(a) -part1.halfl;
            if (d <= 0) {
              /*c is distance squared from line or end to test if is inside sc*/
              if (b < rcut2) intrs+=test_intrpatch(part1,vec2,part1.pcanglsw,-1.0*part2.halfl,intersections);
            }
          }
          //		    printf ("ends %d\n", intrs);
        }
        return intrs;
      }


  }//namespace Geometry

  namespace Potential {

    inline double fanglscale(double a, const CigarParticle &p) {
      // a = r_ij * n_i
      double f;
      if (a <= p.pcanglsw)
        f=0.0;
      else {
        if (a >= p.pcangl)
          f=1.0;
        else {
          f = 0.5 - ((p.pcanglsw + p.pcangl)*0.5 - a )/(p.pcangl - p.pcanglsw);
        }
      }
      return f;
    }

    /*!
     * \brief Hard pair potential for spherocylinders
     */
    class HardSpheroCylinder : public PairPotentialBase {
      private:
        string _brief() { return name; };
        Geometry::Geometrybase *geoPtr;
      public:
        struct prop {
          double halfl;
        };

        std::map<particle::Tid, prop> m;

        HardSpheroCylinder(InputMap &in) {
          name="HardspheroCylinder";
        }

        inline double operator() (const CigarParticle &p1, const CigarParticle &p2, double r2) {
          Point r_cm = geoPtr->vdist(p1,p2);
          Point distvec = Geometry::mindist_segment2segment(p1.dir, m[p1.id].halfl, p2.dir, m[p2.id].halfl, r_cm );
          double mindist=p1.radius+p2.radius;
          if ( distvec.dot(distvec) < mindist*mindist)
            return pc::infty;
          return 0;		
        }

        string info(char w) {
          using namespace Faunus::textio;
          return textio::indent(SUB)+name+"\n";
        }
    };

    template<typename Tcigarsphere>
      class PatchyCigarSphere : public PairPotentialBase {
        private:
          string _brief() {
            return pairpot.brief();
          }
        public:
          Tcigarsphere pairpot;
          vector< vector<double> >* rcutPtr;

          PatchyCigarSphere(InputMap &in) : pairpot(in) {
            rcutPtr=nullptr;
          }

          double operator() (const CigarParticle &a, const CigarParticle &b, const Point &r_cm) {
            assert(rcutPtr!=nullptr);
            //0- isotropic, 1-PSC all-way patch,2 -CPSC cylindrical patch
            //b is sphere, a is spherocylinder
            double s, t, f0, f1, contt;
            double rcut=(*rcutPtr)[a.id][b.id];

            assert( a.halfl < 1e-6 && "First should be cigar then sphere, not opposite!");
            double c = a.dir.dot(r_cm);
            if (c > a.halfl)
              contt = a.halfl;
            else {
              if (c > -a.halfl) contt = c;
              else contt = -a.halfl;
            }
            Point distvec = -r_cm + (a.dir*contt);

            if (atom.list[a.id].patchtype ==0 ) 
              if (atom.list[b.id].patchtype == 0) 
                return pairpot(a,b,distvec.dot(distvec));

            // scaling function: angular dependence of patch1
            Point vec1=Geometry::vec_perpproject(distvec, a.dir);
            vec1.normalize();
            s = vec1.dot(a.patchdir);
            f1 = fanglscale(s,a);

            // scaling function for the length of spherocylinder within cutoff
            double ndistsq = distvec.dot(distvec);
            t = sqrt(rcut*rcut- ndistsq);//TODO cutoff
            if ( contt + t > a.halfl )
              f0 = a.halfl;
            else
              f0 = contt + t;
            if ( contt - t < -a.halfl )
              f0 -= -a.halfl;
            else
              f0 -= contt - t;
            return pairpot(a,b,ndistsq)*f1*(f0+1.0);
          }

          string info(char w) { return pairpot.info(w); }
      };


    template<typename Tcigarcigar>
      class PatchyCigarCigar : public PairPotentialBase {
        private:
          string _brief() {
            return pairpot.brief();
          }
        public:
          vector< vector<double> >* rcutPtr;
          Tcigarcigar pairpot;

          PatchyCigarCigar(InputMap &in) : pairpot(in) {
            rcutPtr=nullptr;
          }

          double operator() (const CigarParticle &a, const CigarParticle &b, const Point &r_cm) {
            assert(rcutPtr!=nullptr);
            //0- isotropic, 1-PSC all-way patch,2 -CPSC cylindrical patch
            if (atom.list[a.id].patchtype >0 ) {
              if (atom.list[b.id].patchtype > 0) {
                //patchy sc with patchy sc
                int i, intrs;
                double rcut=(*rcutPtr)[a.id][b.id];
                double ndistsq;
                double v1, v2, f0, f1, f2, T1, T2, S1, S2,s;
                double intersections[5];
                Point vec1, vec2, vec_intrs, vec_mindist;

                for(i=0;i<5;i++)
                  intersections[i]=0;
                //1- do intersections of spherocylinder2 with patch of spherocylinder1 at.
                // cut distance C
                if (atom.list[a.id].patchtype == 1) {
                  intrs=Geometry::psc_intersect(a,b,r_cm, intersections, rcut);
                } else {
                  if (atom.list[a.id].patchtype == 2) {
                    intrs=Geometry::cpsc_intersect(a,b,r_cm, intersections, rcut);
                  } else {
                    //we dont have anything like this
                    assert(!"Patchtype not implemented!");
                  }
                }
                if (intrs ==0){
                  return 0.0; //sc is all outside patch, attractive energy is 0
                }
                T1=intersections[0]; //points on sc2
                T2=intersections[1];
                //2- now do the same oposite way psc1 in patch of psc2
                for(i=0;i<5;i++)
                  intersections[i]=0;
                if (atom.list[a.id].patchtype == 1) {
                  intrs=Geometry::psc_intersect(b,a,-r_cm, intersections, rcut);
                } else {
                  if (atom.list[a.id].patchtype == 2) {
                    intrs=Geometry::cpsc_intersect(b,a,-r_cm, intersections, rcut);
                  } else {
                    assert(!"Patchtype not implemented!");
                    //we dont have anything like this
                  }
                }
                if (intrs ==0) {
                  return 0.0; //sc is all outside patch, attractive energy is 0
                }
                S1=intersections[0]; //points on sc1
                S2=intersections[1];

                //3a- with two intersection pices calculate vector between their CM
                //-this is for angular orientation
                v1=fabs(S1-S2);
                v2=fabs(T1-T2);
                vec1=a.dir*(S1+S2)*0.5;
                vec2=b.dir*(T1+T2)*0.5;
                vec_intrs=vec2-vec1-r_cm;
                //vec_intrs should be from sc1 t sc2

                //3b - calculate closest distance attractive energy from it
                vec_mindist = Geometry::mindist_segment2segment(a.dir,v1,b.dir,v2,vec_intrs);
                ndistsq=vec_mindist.dot(vec_mindist);

                //4- scaling function1: dependence on the length of intersetions
                //F0=(V1+V2)*0.5/interact->param->sigma;
                f0=(v1+v2)*0.5+1.0;
                //5- scaling function2: angular dependence of patch1
                vec1=Geometry::vec_perpproject(vec_intrs, a.dir);
                vec1.normalize();
                s = vec1.dot(a.patchdir);
                f1 = fanglscale(s,a);

                //6- scaling function3: angular dependence of patch2
                vec1=Geometry::vec_perpproject(-vec_intrs, b.dir);
                vec1.normalize();
                s = vec1.dot(b.patchdir);
                f2 = fanglscale(s,b);

                //7- put it all together and output scale
                return f0*f1*f2*pairpot(a,b,ndistsq);
              }
              else {
                assert(!"PSC w. isotropic cigar not implemented!");
                //patchy sc with isotropic sc
                //we dont have at the moment
              }
            } else {
              if (atom.list[b.id].patchtype > 0) {
                assert(!"PSC w. isotropic cigar not implemented!");
                //isotropic sc with patchy sc
                //we dont have at the moment
              }
              else {
                //isotropic sc with isotropic sc
                Point rclose=Geometry::mindist_segment2segment(a.dir, a.halfl, b.dir, b.halfl, r_cm);
                return pairpot(a,b,rclose.squaredNorm());
              }

            }
            //something we have not implemented
            return 0.0;
          }
          string info(char w) { return pairpot.info(w); }
      };

    template<typename Tcigarcigar, typename Tspheresphere, typename Tcigarsphere>
      class CigarSphereSplit : public PairPotentialBase {
        private:
          string _brief() {
            return pairpot_cc.brief() + " "
              + pairpot_ss.brief() + " "
              + pairpot_cs.brief();
          }

        public:
          PatchyCigarCigar<Tcigarcigar> pairpot_cc;
          Tspheresphere pairpot_ss;
          PatchyCigarSphere<Tcigarsphere> pairpot_cs;

          CigarSphereSplit(InputMap &in) : pairpot_cc(in), pairpot_ss(in), pairpot_cs(in){
            name="CigarSphereSplit";
          }

          double operator() (const particle &a, const particle &b, double r2) const {
            assert(!"Undefined!");
            return 0;
          }

          double operator() (const CigarParticle &a, const CigarParticle &b, const Point &r_cm)
          {
            if (a.halfl<1e-6) {
              // a sphere - b sphere
              if (b.halfl<1e-6) {
                return pairpot_ss(a,b,r_cm.squaredNorm());
              }
              // a sphere - b cigar
              else {
                //PatchyCigarSphere(b,a)
                return pairpot_cs(b,a,r_cm);
              }
            } else {
              // a cigar - b sphere
              if (b.halfl<1e-6) {
                //PatchyCigarSphere(a,b)
                return pairpot_cs(a,b,r_cm);
              }
              // a cigar - b cigar
              else {
                //PatchyCigarCigar
                return pairpot_cc(a,b,r_cm);
              }
            }
            return 0;
          }

          string info(char w) {
            return pairpot_cc.info(w) + pairpot_ss.info(w) + pairpot_cs.info(w);
          }
      };

  }//namespace Potential
}//namespace Faunus

#endif
