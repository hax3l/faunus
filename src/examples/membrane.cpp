#include <faunus/faunus.h>

//#define tab
//#define tabopt
#define taboptintel
//#define tabhermopt
//#define tablinopt

//#define org
//#define orgopt

using namespace Faunus;

/**
 * @brief Setup interactions for coarse grained membrane
 *
 * More information: doi:10/chqzjk
 */
template<class Tbonded, class Tnonbonded, class Tlipid, class Tinput>
void MakeDesernoMembrane(const Tlipid &lipid, Tbonded &bond, Tnonbonded &nb, Tinput &in) {

  using namespace Potential;

  // non-bonded interactions
  auto hid=atom["HD"].id;
  auto tid=atom["TL"].id;
  double sigma   = in("lipid_sigma", 10);  // angstrom
  double epsilon = in("lipid_epsilon", 1); // kT

  CombinedPairPotential<WeeksChandlerAndersen,DebyeHuckel> headhead(in);
  CombinedPairPotential<WeeksChandlerAndersen,CosAttract> tailtail(in);
  CombinedPairPotential<WeeksChandlerAndersen,ChargeNonpolar> headtail(in);

  headhead.first.customSigma(hid, hid, 0.95*sigma);            
  headhead.first.customSigma(hid, tid, 0.95*sigma);            
  headhead.first.customSigma(tid, tid, sigma);                 
  headhead.first.customEpsilon(hid, hid, epsilon);             
  headhead.first.customEpsilon(hid, tid, epsilon);             
  headhead.first.customEpsilon(tid, tid, epsilon);   

  tailtail.first=headhead.first;      // copy initialized WCA to tail-tail and
  headtail.first=headhead.first;      // head-tail pair potential

  nb.pairpot.add(hid, hid, headhead); // Add to main pair potential
  nb.pairpot.add(tid, tid, tailtail);
  nb.pairpot.add(hid, tid, headtail);
#ifdef tab
  //nb.pairpot.print_tabulation(10000); // Print files to compare tabulation with real potential
#endif
#ifdef tabopt
  //nb.pairpot.print_tabulation(10000); // Print files to compare tabulation with real potential
#endif
#ifdef taboptintel
  nb.pairpot.print_tabulation(10000); // Print files to compare tabulation with real potential
#endif
#ifdef tabhermopt
  //nb.pairpot.print_tabulation(10000); // Print files to compare tabulation with real potential
#endif
#ifdef tablinopt
  //nb.pairpot.print_tabulation(10000); // Print files to compare tabulation with real potential
#endif

  // bonded interactions
  double headtail_k=0.5*10*epsilon/(sigma*sigma);
  double headtail_req=4*sigma;
  double fene_k=30*epsilon/(sigma*sigma);
  double fene_rmax=1.5*sigma;

  assert(lipid.size() % 3 == 0);

  auto fene = std::shared_ptr<PairPotentialBase>(new FENE(fene_k,fene_rmax));
  auto harm = std::shared_ptr<PairPotentialBase>(new Harmonic(headtail_k,headtail_req));

  for (int i=lipid.front(); i<lipid.back(); i=i+3) {
    bond.add(i,  i+1, fene ); // Add potentials as *pointers*
    bond.add(i+1,i+2, fene ); // to reduce memory usage
    bond.add(i,  i+2, harm );
  }
}

typedef Geometry::Cuboid Tgeometry;   // specify geometry - here cube w. periodic boundaries
#ifdef tab
typedef Potential::PotentialMapTabulated<Potential::DebyeHuckelLJ> Tpairpot;
#endif
#ifdef tabopt
typedef Potential::PotentialVecTabulated<Potential::DebyeHuckelLJ> Tpairpot;
#endif
#ifdef taboptintel
typedef Potential::PotentialVecTabulated<Potential::DebyeHuckelLJ,Tabulate::tabulatorintel<double>> Tpairpot;
#endif
#ifdef tabhermopt
typedef Potential::PotentialVecTabulated<Potential::DebyeHuckelLJ,Tabulate::tabulatorherm<double>> Tpairpot;
#endif
#ifdef tablinopt
typedef Potential::PotentialVecTabulated<Potential::DebyeHuckelLJ,Tabulate::tabulatorlin<double>> Tpairpot;
#endif
#ifdef org
typedef Potential::PotentialMap<Potential::DebyeHuckelLJ> Tpairpot;
#endif
#ifdef orgopt
typedef Potential::PotentialVec<Potential::DebyeHuckelLJ> Tpairpot;
#endif


int main() {

  cout << textio::splash();           // show faunus banner and credits
  InputMap mcp("membrane.input");   //read input file

  MCLoop loop(mcp);                   // class for handling mc loops
  FormatPQR pqr;                      // PQR structure file I/O
  FormatXTC xtc(1000);
  EnergyDrift sys;                    // class for tracking system energy drifts

  // Energy functions and space
  Energy::Hamiltonian pot;

  auto nonbonded = pot.create( Energy::Nonbonded<Tpairpot,Tgeometry>(mcp) );
  auto bonded    = pot.create( Energy::Bonded() );
  Space spc( pot.getGeometry() );

  // Markov moves and analysis
  Move::AtomicTranslation mv(mcp, pot, spc);
  Move::TranslateRotate gmv(mcp,pot,spc);
  Move::Pivot piv(mcp,pot,spc);
  Move::Isobaric iso(mcp,pot,spc);

  Analysis::BilayerStructure lipidstruct;

  // Load and add polymer to Space
  GroupArray lipids(3);
  FormatAAM aam;                                      // AAM structure file I/O
  string flipid = mcp.get<string>("lipid_file","");
  int Nlipid=mcp("lipid_N",1);
  for (int i=0; i<Nlipid; i++) {
    aam.load(flipid);                                 // Load polymer structure into aam class
    Geometry::FindSpace().find(*spc.geo, spc.p, aam.particles()); // find empty spot
    Group pol = spc.insert( aam.particles() );        // Insert into Space
    if (slp_global()>0.5)
      std::swap( spc.p[pol.front()].z(), spc.p[pol.back()].z() );
    if (slp_global()<mcp("lipid_chargefraction", 0.0))
      spc.p[ pol.front() ].charge = -1;
    lipids.setrange(0, pol.back());
  }
  spc.enroll(lipids);   // Enroll polymer in Space

  // Set up bonded and non-bonded interactions
  MakeDesernoMembrane(lipids, *bonded, *nonbonded, mcp);

  // Place all lipids in xy plane (z=0);
  for (int l=0; l<lipids.numMolecules(); l++) {
    auto g=lipids[l];
    double dz=spc.p[ g.back() ].z();
    for (auto i : g) {
      spc.p[i].z() -= dz;
      spc.geo->boundary(spc.p[i]);
    }
  }
  spc.trial=spc.p; // sync. particle trial vector

  spc.load("state");                                // load old config. from disk (if any)
  sys.init( Energy::systemEnergy(spc,pot,spc.p)  ); // store initial total system energy

  cout << atom.info() + spc.info() + pot.info();

  while ( loop.macroCnt() ) {  // Markov chain 
    while ( loop.microCnt() ) {
      int k=lipids.numMolecules(); //number of lipids
      int i=slp_global.rand() % 3;
      Group g;
      switch (i) {
        case 0:
          mv.setGroup(lipids);
          sys+=mv.move( lipids.size() ); // translate lipid monomers
          break;
        case 1:
          while (k-->0) {
            i=lipids.randomMol(); // pick random lipid molecule
            g=lipids[i];
            g.name="lipid";
            g.setMassCenter(spc);     // mass center needed for rotation
            if (slp_global()>0.5) {
              gmv.setGroup(g);          // tell what to move
              sys+=gmv.move();          // translate/rotate polymers
            } else {
              piv.setGroup(g);          // tell what to move
              sys+=piv.move();          // translate/rotate polymers
            }
          }
          break;
        case 2:
          sys+=iso.move();
          break;
      }
      double ran = slp_global();
      if (ran>0.99) {
        xtc.setbox( nonbonded->geometry.len );
        xtc.save("traj.xtc", spc.p);
      }
      if (ran>0.90)
        lipidstruct.sample(nonbonded->geometry, spc.p, lipids);

    } // end of micro loop

    sys.checkDrift(Energy::systemEnergy(spc,pot,spc.p)); // compare energy sum with current

    // save to disk
    pqr.save("confout.pqr", spc.p);
    spc.save("state");

    cout << loop.timing();
  } // end of macro loop

  // perform unit tests
  UnitTest test(mcp);
  iso.test(test);
  mv.test(test);
  gmv.test(test);
  piv.test(test);
  sys.test(test);
  lipidstruct.test(test);

  // print information
  cout << loop.info() + sys.info() + spc.info() + mv.info() + gmv.info() + iso.info()
    << piv.info() + lipidstruct.info() + test.info();

  return test.numFailed();
}

/** @page example_membrane Example: Membrane Bilayer
  
 This will simulate a 3-bead coarse grained membrane according to
 Cooke and Deserno (doi:10/chqzjk). Each bead interacts with a
 Weeks-Chandler-Andersen potential, while tail-tail interactions
 have an additional long range attractive potential. There is preliminary
 support for charged head groups (effect on elastic properties is unknown).

 The following moves are included:
 - Lipid rotation, translation and pivot
 - Monomer translation
 - Iso-tension move

 Run this example from the `examples` directory:

 ~~~~~~~~~~~~~~~~~~~
 $ make
 $ cd src/examples
 $ ./membrane.run
 ~~~~~~~~~~~~~~~~~~~

 ![Bilayer formed by 3-bead CG lipid model](membrane-3bead.jpg)

 membrane.cpp
 ============

 \includelineno examples/membrane.cpp
*/

