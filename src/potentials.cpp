#ifndef FAU_POT_BASE_H
#define FAU_POT_BASE_H

#include <faunus/common.h>
#include <faunus/point.h>
#include <faunus/inputfile.h>
#include <faunus/physconst.h>
#include <faunus/geometry.h>
#include <faunus/potentials.h>
#include <faunus/textio.h>

namespace Faunus {

  namespace Potential {

    PairPotentialBase::PairPotentialBase() {
      setScale(1);
    }

    PairPotentialBase::~PairPotentialBase() { }
    
    /*!
     * This functions sets the energy scaling as returned by tokT().
     */
    void PairPotentialBase::setScale(double s)  {
      assert(s!=0 && "Energy scaling should be non-zero.");
      _tokT=s;
      _setScale(_tokT);
    }
    
    void PairPotentialBase::_setScale(double s) {}
   
    /*!
     * To optimize certain potentials, the returned energy is not required to be
     * in units of kT. Instead, unit scaling can be done after energy summations,
     * thereby reducing computational time. Therefore, to convert energy to kT,
     * multiply with the scaling factor returned by this function.
     */
    double PairPotentialBase::tokT() const { return _tokT; }

    /*!
     * This will reset the temperature to the specified value. By default this function
     * does nothing, although in Debug mode it will throw an exception if derived classes
     * does not implement it (and is called).
     */
    void PairPotentialBase::setTemperature(double) {
      assert(!"Not implemented.");
    }

    string PairPotentialBase::brief() {
      assert(!name.empty() && "Potential must have a name.");
      return _brief();
    }

    Harmonic::Harmonic(double forceconst, double eqdist) : k(forceconst), req(eqdist) {
      name="Harmonic";
    }

    Harmonic::Harmonic(InputMap &in, string pfx) {
      name="Harmonic";
      k  = in.get<double>( pfx+"forceconst", 0);
      req = in.get<double>( pfx+"eqdist", 0);
    }

    void Harmonic::_setScale(double s) {
      _tokT=s;
      k=k/_tokT;
    }

    double Harmonic::operator() (const particle &a, const particle &b, double r2) const {
      double d=sqrt(r2)-req;
      return k*d*d;
    }

    string Harmonic::_brief() {
      using namespace Faunus::textio;
      std::ostringstream o;
      o << name << ": k=" << k*tokT() << kT << "/" << angstrom << squared << " req=" << req << _angstrom; 
      return o.str();
    }

    HardSphere::HardSphere() {
      name="Hardsphere";
    }

    /*!
     * This is a compatibility constructor - no data is read from the InputMap.
     */
    HardSphere::HardSphere(InputMap& in) {
      name="Hardsphere";
    }
   
    string HardSphere::_brief() {
      return name;
    }

    string HardSphere::info(char w) {
      using namespace Faunus::textio;
      return textio::indent(SUB)+name+"\n";
    }

    LennardJones::LennardJones() {
      name="Lennard-Jones";
    }
   
    /*!
     * \param in InputMap is scanned for the keyword \c lj_eps and should be in units of kT
     */
    LennardJones::LennardJones(InputMap &in, string pfx) {
      name="Lennard-Jones";
      eps = 4*in.get<double>( pfx+"eps", 0.04, name+" epsilon (kT)" );
    }
    
    string LennardJones::_brief() {
      std::ostringstream o;
      o << name << ": " << textio::epsilon+"(LJ)=" << eps*tokT()/4 << textio::kT;
      return o.str();
    }
    
    void LennardJones::_setScale(double s) {
      _tokT=s;
      eps=eps/_tokT;
    }

    string LennardJones::info(char w) {
      using namespace Faunus::textio;
      std::ostringstream o;
      o << pad(SUB,w+1,epsilon+"(LJ)") << eps*tokT()/4 << kT << endl;
      return o.str();
    }

    /*!
     * \param in is scanned for the keywords \c prefix_threshold (angstrom) and \c prefix_depth (kT).
     * \param prefix InputMap keyword prefix. Default is "squareWell"
     */
    SquareWell::SquareWell(InputMap &in, string prefix) {
      name="Square Well";
      threshold = in.get<double>(prefix+"_threshold", 0, name+" upper threshold (AA)");
      depth     = in.get<double>(prefix+"_depth", 0, name+" depth (kT)");
    }
    
    void SquareWell::_setScale(double s) {
      _tokT=s;
      depth=depth/_tokT;
    }
    
    string SquareWell::_brief() {
      std::ostringstream o;
      o << name << ": u=" << depth*tokT() << textio::kT << " r=" << threshold;
      return o.str();
    }

    string SquareWell::info(char w) {
      using namespace Faunus::textio;
      std::ostringstream o;
      o << pad(SUB,w,"Threshold") << threshold << " "+angstrom+" (surface-surface)" << endl;
      o << pad(SUB,w,"Depth") << depth*tokT() << kT << endl;
      return o.str();
    }

    SquareWellHydrophobic::SquareWellHydrophobic(InputMap &in, string prefix) : SquareWell(in,prefix) {
      name="Hydrophobic " + name;
    }

    /*!
     * \param in InputMap is scanned for the keyword \c softrep_sigma which should be in angstrom
     */
    SoftRepulsion::SoftRepulsion(InputMap &in) {
      name="Repulsive r6";
      sigma6 = pow( in.get<double>( "softrep_sigma", 5 ), 6);
    }
    
    string SoftRepulsion::_brief() {
      std::ostringstream o;
      o << name << ": " << textio::sigma  << pow(sigma6*tokT(),1/6.) << textio::_angstrom;
      return o.str();
    }
    
    void SoftRepulsion::_setScale(double s) {
      _tokT=s;
      sigma6=sigma6/_tokT;
    }

    string SoftRepulsion::info(char w) {
      using namespace Faunus::textio;
      std::ostringstream o;
      o << textio::pad(SUB,w+1,textio::sigma) << pow(sigma6*tokT(),1/6.) << textio::_angstrom << endl;
      return o.str();
    }
 
    /*!
     * The following input keywords are searched searched:
     * \li \c temperature [Kelvin, default = 298.15]
     * \li \c epsilon_r - relative dielectric constant. Default is 80.
     * \li \c depsdt - temperature dependence of dielectric constant, \f$ \partial\epsilon_r/\partial T\approx-0.368\f$ for water.
     */
    Coulomb::Coulomb(InputMap &in) {
      name="Coulomb";
      pc::setT ( in.get<double>("temperature", 298.15, "Absolute temperature (K)") );
      epsilon_r = in.get<double>("epsilon_r",80., "Dielectric constant");
      depsdt = in.get<double>("depsdt", -0.368, "See documentation") * pc::T() / epsilon_r;
      lB=pc::lB( epsilon_r );
      setScale(lB);
    }

    void Coulomb::_setScale(double s) {
    }

    string Coulomb::_brief() {
      std::ostringstream o;
      o << name << ": lB=" << lB << " eps_r=" << epsilon_r << " T=" << pc::T();
      return o.str();
    }

    double Coulomb::bjerrumLength() const {
      return lB;
    }

    string Coulomb::info(char w) {
      using namespace textio;
      std::ostringstream o;
      o << pad(SUB,w,"Temperature") << pc::T() << " K" << endl
        << pad(SUB,w,"Dielectric constant") << epsilon_r << endl
        << pad(SUB,w+6,"T"+partial+epsilon+"/"+epsilon+partial+"T") << depsdt << endl
        << pad(SUB,w,"Bjerrum length") << lB << " "+angstrom << endl;

      return o.str();
    }

    /*!
     * In addition to the keywords from Potential::Coulomb, InputMap is searched for:
     * \li \c dh_ionicstrength [mol/l] 
     * \li \c dh_debyelength [angstrom] (only if I=0, default)
     */
    DebyeHuckel::DebyeHuckel(InputMap &in) : Coulomb(in) {
      double I;
      const double zero=1e-10;
      name="Debye-Huckel";
      c=8 * lB * pc::pi * pc::Nav / 1e27;
      I=in.get<double>("dh_ionicstrength",0, "Ionic strength (mol/l)");  // [mol/l]
      k=sqrt( I*c );
      if (k<zero)
        k=1/in.get<double>("dh_debyelength", 1/zero, "Debye length (AA)"); // [A]
    }

    string DebyeHuckel::_brief() {
      std::ostringstream o;
      o << Coulomb::_brief() << " I=" << ionicStrength();
      return o.str();
    }

    double DebyeHuckel::ionicStrength() const {
      return k*k/c;
    }

    double DebyeHuckel::debyeLength() const {
      return 1/k;
    }

    /*!
     * The Debye-Huckel potential is temperature dependent and contains entropy
     * contributions from both solvent and salt degrees of freedom.
     * This function return the entropy of interaction between a pair of
     * particles interacting with an effective Debye-Huckel potential. This is done by
     * taking the temperature derivate of w(R):
     * \f[
     * S(r_{ij})/k_B = -\frac{ \partial w(r_{ij},T) } {k_B \partial T} = \beta w_{ij}\left [ \alpha - \frac{\kappa r_{ij}(\alpha+1)}{2}\right ]
     * \f]
     * where \f$ \alpha=T \partial \epsilon_r/\epsilon_r\partial T\f$
     * is determined experimentally for pure water. To get the entropy from salt ions
     * only, set \f$\alpha=0\f$ via the InputMap.
     *
     * \param  betaw    Inter particle free energy, \f$\beta w\f$, in units of kT.
     * \param  r        Inter particle distance
     * \return Interaction entropy \f$ S(r_{ij})/k_B = \beta TS(r_{ij})\f$
     * \todo   Optimize
     */
    double DebyeHuckel::entropy(double betaw, double r) const {
      return betaw * (depsdt - 0.5*k*r*(depsdt+1));
    }

    string DebyeHuckel::info(char w) {
      using namespace Faunus::textio;
      std::ostringstream o;
      o << Coulomb::info(w);
      o << pad(SUB,w,"Ionic strength") << ionicStrength() << " mol/l" << endl;
      o << pad(SUB,w+1,"Debye length, 1/\u03BA") << debyeLength() << " "+angstrom << endl;
      return o.str();
    }

    /*!
     * \f$ \beta u(r) = l_B \frac{ z_1 z_2 }{r}\f$
     */
     double MultipoleEnergy::ionion(double z1, double z2, double r) {
      return lB*z1*z2/r;
    }

    /*!
     * \f$ \beta u(r) = -l_B \frac{ z a_z }{r^2}\f$
     */
     double MultipoleEnergy::iondip(double z, const Point &a, double r) {
      return -lB*z*a.z/(r*r);
    }

    /*!
     * \f$ \beta u(r) = l_B \frac{a_x b_x + a_y b_y - 2a_z b_z  }{r^3}\f$
     */
    double MultipoleEnergy::dipdip(const Point &a, const Point &b, double r) {
      return lB*( a.x*b.x + a.y*b.y - 2*a.z*b.z ) / (r*r*r);
    }

  } //Potential namespace

} //Faunus namespace
#endif
