#ifndef FAU_taulate_h
#define FAU_taulate_h

#ifndef SWIG

#include <iostream>
#include <fstream>
#include <vector>
#include <functional>
#include <map>
#include <algorithm>
#include <cmath>
#include <memory>
#include <assert.h>

#endif

namespace Faunus {
  
  namespace Tabulate {
    /* base class for all tabulators - no dependencies */
    template<typename T=double>
    class tabulatorbase {
    protected:
      T utol=0.01, ftol = -1, umaxtol = -1, fmaxtol = -1, rmin, rmax;
      
      T numdr = 0.0001; // dr for derivative evaluation
      
      // First derivative with respect to r2 numerically
      T f1(std::function<T(T)> &f, T r2) {//Tparticle a, Tparticle b, std::function<T(Tparticle&,Tparticle&,T)> f, T r2) {
        return (f(r2+numdr*0.5)-f(r2-numdr*0.5))/(numdr);
      }
      
      // Second derivative with respect to r2 numerically
      T f2(std::function<T(T)> &f, T r2) {//Tparticle a, Tparticle b, std::function<T(Tparticle&,Tparticle&,T)> f, T r2) {
        return (f1(f,r2+numdr*0.5)-f1(f,r2-numdr*0.5))/(numdr);
      }
      
    public:
      struct data {
        std::vector<T> r2;  // r2 for intervals
        std::vector<T> c;   // c for coefficents
        T rmin2, rmax2;     // useful to save these with table
      };
      void setTolerance(T _utol, T _ftol=-1, T _umaxtol=-1, T _fmaxtol=-1) {
        utol = _utol;
        ftol = _ftol;
        umaxtol = _umaxtol;
        fmaxtol = _fmaxtol;
      }
      void setRange(T _rmin, T _rmax) {
        rmin = _rmin;
        rmax = _rmax;
      }
      void setNumdr(T _numdr) {
        numdr = _numdr;
      }
      
      std::string info(char w=20) {
        using namespace Faunus::textio;
        std::ostringstream o( "" );
        o << pad(SUB,w,"Rmin") << rmin << std::endl;
        o << pad(SUB,w,"Rmax") << rmax << std::endl;
        o << pad(SUB,w,"Utol") << utol << std::endl;
        o << pad(SUB,w,"Ftol") << ftol << std::endl;
        o << pad(SUB,w,"Umaxtol") << umaxtol << std::endl;
        o << pad(SUB,w,"Fmaxtol") << fmaxtol << std::endl;
        return o.str();
      }
    };

    /*  Tabulator
     Code mainly from MolSim (Per Linse) with some upgrades
     Translated from Fortran into C++
     Ref of algorithm Andrea,Swope,Andersen JCP. 79 (1983)
     */
    template<typename T=double>
    struct tabulator : public tabulatorbase<T> {
      typedef tabulatorbase<T> base; // for convenience, only
      
      int mngrid = 1200;    // Max number of controlpoints
      int ndr = 100;        // Max number of trials to decr dr
      T drfrac = 0.9;  // Multiplicative factor to decr dr
      
      void set_mngrid(int _mngrid) {
        mngrid = _mngrid;
      }
      
      void set_ndr(int _ndr) {
        mngrid = _ndr;
      }
      
      // gets tabulated value at r2
      T eval(const typename base::data& d, T r2) const {
        const typename std::vector<T>::const_iterator low = std::lower_bound(d.r2.begin(), d.r2.end(), r2);
        int pos = (low-d.r2.begin()-1);
        T min = d.r2[pos];
        T dz = r2-min;
        int pos6 = 6.0*pos;
        T usum =  d.c[pos6+0]+
                  dz*(d.c[pos6+1]+
                  dz*(d.c[pos6+2]+
                  dz*(d.c[pos6+3]+
                  dz*(d.c[pos6+4]+
                  dz*(d.c[pos6+5])
                  ))));
        return usum;
        
      }
      
      std::vector<T> SetUBuffer(T rlow, 
                                T zlow, 
                                T rupp, 
                                T zupp,
                                T u0low,
                                T u1low,
                                T u2low,
                                T u0upp,
                                T u1upp,
                                T u2upp) {
        
        std::vector<T> ubuft;
        
        ubuft.push_back(zlow);
        // Zero potential and force return no coefficients
        if (u0low == 0.0 && u1low == 0.0) {
          for (int i = 0; i < 6; i++)
            ubuft.push_back(0.0);
          return ubuft;
        }
        
        T dz1 = zupp-zlow;
        T dz2 = dz1*dz1;
        T dz3 = dz2*dz1;
        
        T w0low = u0low;
        T w1low = u1low;
        T w2low = u2low;
        
        T w0upp = u0upp;
        T w1upp = u1upp;
        T w2upp = u2upp;
        
        T c0 = w0low;
        T c1 = w1low;
        T c2 = w2low*0.5;
        
        T a = 6.0*(w0upp-c0-c1*dz1-c2*dz2)/dz3;
        T b = 2.0*(w1upp-c1-2.0*c2*dz1)/dz2;
        T c = (w2upp-2*c2)/dz1;
        
        T c3 = (10.0*a-12.0*b+3.0*c)/6.0;
        T c4 = (-15.0*a+21.0*b-6.0*c)/(6.0*dz1);
        T c5 = (2.0*a-3.0*b+c)/(2.0*dz2);
        
        ubuft.push_back(c0);
        ubuft.push_back(c1);
        ubuft.push_back(c2);
        ubuft.push_back(c3);
        ubuft.push_back(c4);
        ubuft.push_back(c5);
        
    #ifdef D_SetUBuffer
        std::cout << "zlow="<<zlow<<std::endl;
        std::cout << "zupp="<<zupp<<std::endl;
        std::cout << "dz1="<<dz1<<std::endl;
        std::cout<< "dz2="<<dz2<<std::endl;
        std::cout<< "dz3="<<dz3<<std::endl;
        std::cout<< "a="<<a<<std::endl;
        std::cout<< "b="<<b<<std::endl;
        std::cout<< "c="<<c<<std::endl;
        std::cout<< "c0="<<c0<<" , ubuft.at(1)="<<ubuft.at(1)<< std::endl;
        std::cout<< "c1="<<c1<<" , ubuft.at(2)="<<ubuft.at(2)<< std::endl;
        std::cout<< "c2="<<c2<<" , ubuft.at(3)="<<ubuft.at(3)<< std::endl;
        std::cout<< "c3="<<c3<<" , ubuft.at(4)="<<ubuft.at(4)<< std::endl;
        std::cout<< "c4="<<c4<<" , ubuft.at(5)="<<ubuft.at(5)<< std::endl;
        std::cout<< "c5="<<c5<<" , ubuft.at(6)="<<ubuft.at(6)<< std::endl;
    #endif
        
        return ubuft;
      }
      
      std::vector<bool> CheckUBuffer(std::vector<T>& ubuft, 
                                     T rlow,
                                     T rupp,
                                     std::function<T(T)> &f) {/*Tparticle a, 
                                     Tparticle b,
                                     std::function<T(Tparticle&,Tparticle&,T)> f) {*/
        
        // vb[0]: Tolerance is approved
        // vb[1]: A repulsive part is found
        std::vector<bool> vb;
        vb.push_back(false);
        vb.push_back(false);
        
        // Number of points to control
        int ncheck = 11;
        double dr = (rupp-rlow)/(ncheck-1);
        
        for (int i = 0; i < ncheck; i++) {
          T r1 = rlow+dr*((T)i);
          T r2 = r1*r1;
          T u0 = f(r2);
          T u1 = base::f1(f,r2);
          T dz = r2-rlow*rlow;
          T usum =  ubuft.at(1)+
                    dz*(ubuft.at(2)+
                    dz*(ubuft.at(3)+
                    dz*(ubuft.at(4)+
                    dz*(ubuft.at(5)+
                    dz*ubuft.at(6)
                    ))));
          
          T fsum =  ubuft.at(2)+
                    dz*(2.0*ubuft.at(3)+
                    dz*(3.0*ubuft.at(4)+
                    dz*(4.0*ubuft.at(5)+
                    dz*(5.0*ubuft.at(6)
                    ))));
          
    #ifdef D_CheckUBuffer
          std::cout<< "u0=pot("<<sqrt(r2)<<","<<r2<<"): " << u0 << std::endl;
          std::cout<< "u1=pot("<<sqrt(r2)<<","<<r2<<"): " << u1 << std::endl;
          std::cout << "usum=" << usum << " abs(usum-u0) " << std::abs(usum-u0) << " u0 " << u0 << " base::utol " << base::utol << std::endl;
          std::cout << "fsum=" << fsum << " abs(fsum-u1) " << std::abs(fsum-u1) << " u1 " << u1 << " base::ftol " << base::ftol << std::endl;
    #endif
          
          if (std::abs(usum-u0) > base::utol) {
      #ifdef D_CheckUBuffer
            std::cout << "Failed by u" << std::endl;
      #endif
            return vb;
          }
          if (base::ftol != -1 && std::abs(fsum-u1) > base::ftol) {
      #ifdef D_CheckUBuffer
            std::cout << "Failed by f" << std::endl;
      #endif
            return vb;    
          }
          if (base::umaxtol != -1 && std::abs(usum) > base::umaxtol)
            vb.at(1) = true;
          if (base::fmaxtol != -1 && std::abs(usum) > base::fmaxtol)
            vb.at(1) = true;
          
    #ifdef D_CheckUBuffer
          std::cout << "Next point" << std::endl;
    #endif
          
          
        }
        
        vb.at(0) = true;
        
        return vb;
      }
      
      
      typename base::data
      generate(std::function<T(T)> &f) {//(Tparticle a, Tparticle b, std::function<T(Tparticle&,Tparticle&,T)> f) {
        
        assert(base::rmin >= 0.0);
        assert(base::rmax >= 0.0);
        // rmin larger than rmax
        assert(base::rmin < base::rmax);
        // Tolerance of potential not set or negative
        assert(base::utol > 0.0000000000001);
        // Tolerance of potential set and negative
        assert(base::ftol == -1 || base::ftol > 0.0);
        // Tolerance of repulsion of potential set and negative
        assert(base::umaxtol == -1 || base::umaxtol > 0.0);
        // Tolerance of repulsion of potential set and negative
        assert(base::fmaxtol == -1 || base::fmaxtol > 0.0);
        
        typename base::data td;
        td.rmax2 = base::rmax*base::rmax;
        td.rmin2 = base::rmin*base::rmin;
        
        T minv = base::rmin;
        T maxv = base::rmax;
        
        T rumin = minv;
        
        T maxv2 = maxv*maxv;
        
        T dr = maxv-minv;
        
        T rupp = maxv;
        T zupp = maxv2;
        bool repul = false; // Stop tabulation if repul is true
        
        td.r2.push_back(zupp);
        
        
        int i;
        for (i=0; i < mngrid; i++) {
          T rlow = rupp;
          T zlow;
          std::vector<T> ubuft;
          int j;
          
          dr=(rupp-minv);
          
          for (j=0; j < ndr; j++) {
            zupp=rupp*rupp;
            rlow = rupp-dr;
            if (rumin > rlow)
              rlow = rumin;
            
    #ifdef D_create_molsim_r2
            std::cout <<std::endl<<"D_create_molsim_r2 dr="<<dr<<std::endl;
            std::cout << "D_create_molsim_r2 rupp="<<rupp<<" minv " <<minv<<std::endl;
            std::cout << "D_create_molsim_r2 rlow="<<rlow<<" maxv " <<maxv<<std::endl<<std::endl;
    #endif
            
            zlow = rlow*rlow;
            
            T u0low = f(zlow);
            T u1low = base::f1(f,zlow);
            T u2low = base::f2(f,zlow);
            
            T u0upp = f(zupp);
            T u1upp = base::f1(f,zupp);
            T u2upp = base::f2(f,zupp);
            
            ubuft = SetUBuffer(rlow,zlow,rupp,zupp,u0low,u1low,u2low,u0upp,u1upp,u2upp);
            std::vector<bool> vb = CheckUBuffer(ubuft,rlow,rupp,f);
            repul = vb.at(1);
            if (vb.at(0) == true) {
    #ifdef D_MR2
              std::cout << std::endl << "SUCCESS" << std::endl << std::endl;
    #endif
              rupp=rlow;
              break;
            } else {
    #ifdef D_MR2
              std::cout << std::endl << "FAILED (decr dr)" << std::endl << std::endl;
    #endif          
            }
            dr*=drfrac;
          }
          // Error if j has gone over ndr
          assert(j < ndr);
          td.r2.push_back(zlow);
          // Wrong size of ubuft, minvalue + 6 coefficients
          assert(ubuft.size() == 7);
    #ifdef D_MR2
          std::cout<< "saving into td zlow="<<zlow<<" rlow "<<rlow<<" rupp "<<rupp<<std::endl;
          std::cout<< std::endl << "Printing ubuft for dr: " << dr << std::endl;
          for (int k = 0; k < ubuft.size(); k++)
            std::cout<< k << ": " << ubuft.at(k) << std::endl;
    #endif
          for (unsigned int k = 1; k < ubuft.size(); k++)
            td.c.push_back(ubuft.at(k));
          
          // Entered a highly repulsive part, stop tabulation
          if (repul == true) {
            rumin = rlow;
            td.rmin2 = rlow*rlow;
    #ifdef D_MR2
            std::cout<< "Repulsive region found" << std::endl;
    #endif
          }
    #ifdef D_MR2
          std::cout<< "rumin="<<rumin<<std::endl;
    #endif
          if (rlow <= rumin || repul == true)
            break;
        }
        // mngrid not enough (increase utol or ftol)
        assert(i < mngrid);
        
        // Sort td
        typename base::data tdsort;
        tdsort.rmax2 = td.rmax2;
        tdsort.rmin2 = td.rmin2;
        tdsort.r2.push_back(td.r2.at(td.r2.size()-1));
        for (int i = td.r2.size()-2; i >= 0; i--) {
          tdsort.r2.push_back(td.r2.at(i));
          for (int j = 0; j < 6; j++) {
            tdsort.c.push_back(td.c.at(6*i+j));
          }
        }
        
        return tdsort;
        
      }
      
      
      
      typename base::data
      generate_full(std::function<T(T)> &f) {
        
        typename base::data tg = generate(f);
        
        // Assume zero at from max to infinity
        tg.rmax2 = 100000000.0;
        tg.r2.push_back(pc::infty);
        tg.c.push_back(0.0);
        tg.c.push_back(0.0);
        tg.c.push_back(0.0);
        tg.c.push_back(0.0);
        tg.c.push_back(0.0);
        tg.c.push_back(0.0);
        
        // Assume infinity from min to zero
        typename std::vector<T>::iterator it = tg.r2.begin();
        tg.rmin2 = -0.1;
        tg.r2.insert(it, 0.0);
        
        it = tg.c.begin();
        tg.c.insert ( it , 0.0 );
        it = tg.c.begin();
        tg.c.insert ( it , 0.0 );
        it = tg.c.begin();
        tg.c.insert ( it , 0.0 );
        it = tg.c.begin();
        tg.c.insert ( it , 0.0 );
        it = tg.c.begin();
        tg.c.insert ( it , 0.0 );
        it = tg.c.begin();
        tg.c.insert ( it , 100000.0 );
        
        return tg;
        
      }
      
      
      
      typename base::data
      generate_empty() {
        
        typename base::data tg;
        
        // Assume zero at from zero to infinity
        tg.rmin2 = -0.1;
        tg.rmax2 = 100000000.0;
        tg.r2.push_back(0.0);
        tg.c.push_back(0.0);
        tg.c.push_back(0.0);
        tg.c.push_back(0.0);
        tg.c.push_back(0.0);
        tg.c.push_back(0.0);
        tg.c.push_back(0.0);
        
        return tg;
        
      }

      
      std::string print(typename base::data &d) {
        std::ostringstream o( "Size of r2: " + d.r2.size() );
        o << "rmax2 r2=" << d.rmax2 << " r=" << std::sqrt(d.rmax2) << std::endl;
        o << "rmin2 r2=" << d.rmin2 << " r=" << std::sqrt(d.rmin2) << std::endl;
        for (unsigned int i = 0; i < d.r2.size(); i++) {
          o << i << ": r2=" << d.r2.at(i) << " r=" << std::sqrt(d.r2.at(i)) << std::endl;
          if (i != (d.r2.size()-1)) {
            o << "coeffs:";
            for (unsigned int j = 0; j < 6; j++) {
              o << " "<< d.c.at(i*6+j) <<",";
            }
          }
          o << std::endl;
        }
        return o.str();
      }
      
    };
    
  }

}

#endif