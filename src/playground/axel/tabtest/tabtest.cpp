#include <faunus/faunus.h>
#include <faunus/tabulate.h>

//#define tab
//#define tabopt

//#define tabherm
//#define taboptherm

//#define tablin
//#define taboptlin

#define org


//#define tabsingle
//#define taboptsingle
//#define orgsingle

using namespace Faunus;
using namespace std;


typedef Geometry::Cuboid Tgeometry;
typedef Potential::DebyeHuckelLJ Tpairpot;

int main(int argc, char** argv) {
  cout << textio::splash();
  InputMap mcp("cluster.input");
  MCLoop loop(mcp);                    // class for handling mc loops
  FormatPQR pqr;                       // PQR structure file I/O
  FormatAAM aam;                       // AAM structure file I/O
  FormatXTC xtc(1000);                 // XTC gromacs trajectory format
  EnergyDrift sys;                     // class for tracking system energy drifts
  UnitTest test(mcp);

  Energy::Hamiltonian pot;
#ifdef tab
  Energy::Nonbonded<Potential::PotentialMapTabulated<Tpairpot,Tabulate::tabulator<double>>,Tgeometry> nb(mcp); // Tabulation
  nb.pairpot.add(atom["Ca"].id,atom["Ca"].id,Tpairpot(mcp));
  nb.pairpot.add(atom["Ca"].id,atom["Cl"].id,Tpairpot(mcp));
  nb.pairpot.add(atom["Ca"].id,atom["Na"].id,Tpairpot(mcp));
  nb.pairpot.add(atom["Na"].id,atom["Na"].id,Tpairpot(mcp));
  nb.pairpot.add(atom["Na"].id,atom["Cl"].id,Tpairpot(mcp));
  nb.pairpot.add(atom["Cl"].id,atom["Cl"].id,Tpairpot(mcp));
  auto nonbonded = pot.create(nb);
#endif
#ifdef tablin
  Energy::Nonbonded<Potential::PotentialMapTabulated<Tpairpot,Tabulate::tabulatorlin<double>>,Tgeometry> nb(mcp); // Tabulation
  nb.pairpot.add(atom["Ca"].id,atom["Ca"].id,Tpairpot(mcp));
  nb.pairpot.add(atom["Ca"].id,atom["Cl"].id,Tpairpot(mcp));
  nb.pairpot.add(atom["Ca"].id,atom["Na"].id,Tpairpot(mcp));
  nb.pairpot.add(atom["Na"].id,atom["Na"].id,Tpairpot(mcp));
  nb.pairpot.add(atom["Na"].id,atom["Cl"].id,Tpairpot(mcp));
  nb.pairpot.add(atom["Cl"].id,atom["Cl"].id,Tpairpot(mcp));
  auto nonbonded = pot.create(nb);
#endif
#ifdef tabherm
  Energy::Nonbonded<Potential::PotentialMapTabulated<Tpairpot,Tabulate::tabulatorherm<double>>,Tgeometry> nb(mcp); // Tabulation
  nb.pairpot.add(atom["Ca"].id,atom["Ca"].id,Tpairpot(mcp));
  nb.pairpot.add(atom["Ca"].id,atom["Cl"].id,Tpairpot(mcp));
  nb.pairpot.add(atom["Ca"].id,atom["Na"].id,Tpairpot(mcp));
  nb.pairpot.add(atom["Na"].id,atom["Na"].id,Tpairpot(mcp));
  nb.pairpot.add(atom["Na"].id,atom["Cl"].id,Tpairpot(mcp));
  nb.pairpot.add(atom["Cl"].id,atom["Cl"].id,Tpairpot(mcp));
  auto nonbonded = pot.create(nb);
#endif
#ifdef tabopt
  Energy::Nonbonded<Potential::PotentialVecTabulated<Tpairpot,Tabulate::tabulator<double>>,Tgeometry> nb(mcp); // Tabulation optimized
  nb.pairpot.add(atom["Ca"].id,atom["Ca"].id,Tpairpot(mcp));
  nb.pairpot.add(atom["Ca"].id,atom["Cl"].id,Tpairpot(mcp));
  nb.pairpot.add(atom["Ca"].id,atom["Na"].id,Tpairpot(mcp));
  nb.pairpot.add(atom["Na"].id,atom["Na"].id,Tpairpot(mcp));
  nb.pairpot.add(atom["Na"].id,atom["Cl"].id,Tpairpot(mcp));
  nb.pairpot.add(atom["Cl"].id,atom["Cl"].id,Tpairpot(mcp));
  auto nonbonded = pot.create(nb);
#endif
#ifdef taboptlin
  Energy::Nonbonded<Potential::PotentialVecTabulated<Tpairpot,Tabulate::tabulatorlin<double>>,Tgeometry> nb(mcp); // Tabulation optimized
  nb.pairpot.add(atom["Ca"].id,atom["Ca"].id,Tpairpot(mcp));
  nb.pairpot.add(atom["Ca"].id,atom["Cl"].id,Tpairpot(mcp));
  nb.pairpot.add(atom["Ca"].id,atom["Na"].id,Tpairpot(mcp));
  nb.pairpot.add(atom["Na"].id,atom["Na"].id,Tpairpot(mcp));
  nb.pairpot.add(atom["Na"].id,atom["Cl"].id,Tpairpot(mcp));
  nb.pairpot.add(atom["Cl"].id,atom["Cl"].id,Tpairpot(mcp));
  auto nonbonded = pot.create(nb);
#endif
#ifdef taboptherm
  Energy::Nonbonded<Potential::PotentialVecTabulated<Tpairpot,Tabulate::tabulatorherm<double>>,Tgeometry> nb(mcp); // Tabulation optimized
  nb.pairpot.add(atom["Ca"].id,atom["Ca"].id,Tpairpot(mcp));
  nb.pairpot.add(atom["Ca"].id,atom["Cl"].id,Tpairpot(mcp));
  nb.pairpot.add(atom["Ca"].id,atom["Na"].id,Tpairpot(mcp));
  nb.pairpot.add(atom["Na"].id,atom["Na"].id,Tpairpot(mcp));
  nb.pairpot.add(atom["Na"].id,atom["Cl"].id,Tpairpot(mcp));
  nb.pairpot.add(atom["Cl"].id,atom["Cl"].id,Tpairpot(mcp));
  auto nonbonded = pot.create(nb);
#endif
#ifdef org
  Energy::Nonbonded<Potential::PotentialMap<Tpairpot>,Tgeometry> nb(mcp); //  Non-tabulation
  nb.pairpot.add(atom["Ca"].id,atom["Ca"].id,Tpairpot(mcp));
  nb.pairpot.add(atom["Ca"].id,atom["Cl"].id,Tpairpot(mcp));
  nb.pairpot.add(atom["Ca"].id,atom["Na"].id,Tpairpot(mcp));
  nb.pairpot.add(atom["Na"].id,atom["Na"].id,Tpairpot(mcp));
  nb.pairpot.add(atom["Na"].id,atom["Cl"].id,Tpairpot(mcp));
  nb.pairpot.add(atom["Cl"].id,atom["Cl"].id,Tpairpot(mcp));
  auto nonbonded = pot.create(nb);
#endif
  

#ifdef tabsingle
  auto nonbonded = pot.create( Energy::Nonbonded<Potential::PotentialTabulate<Tpairpot>,Tgeometry>(mcp) );
#endif 
#ifdef taboptsingle
  auto nonbonded = pot.create( Energy::Nonbonded<Potential::PotentialTabulateVec<Tpairpot>,Tgeometry>(mcp) );
#endif 
#ifdef orgsingle
  auto nonbonded = pot.create( Energy::Nonbonded<Tpairpot,Tgeometry>(mcp) );
#endif 
  
  //nb.pairpot.print_tabulation(); // Print files to compare tabulation with real potential
  
  
  Space spc( pot.getGeometry() );

  // Add salt
  GroupAtomic salt(spc, mcp);
  salt.name="Salt";
  spc.enroll(salt);

  spc.load("state");

  Move::AtomicTranslation mv(mcp, pot, spc);
  mv.setGroup(salt);   // specify atomic particles to be moved

  short na = atom["Na"].id;
  short ca = atom["Ca"].id;
  short cl = atom["Cl"].id;
  
  Analysis::RadialDistribution<float,unsigned int> rdf1(0.2); // 0.2 Å resolution
  Analysis::RadialDistribution<float,unsigned int> rdf2(0.2); // 0.2 Å resolution

  
  
  sys.init( Energy::systemEnergy(spc,pot,spc.p) );

  cout << atom.info() << spc.info() << pot.info() << textio::header("MC Simulation Begins!");
  
  while ( loop.macroCnt() ) {  // Markov chain 
    while ( loop.microCnt() ) {
      xtc.save("out.xtc", spc.p);
      int i=rand() % 1;
      switch (i) {
        case 0:
          mv.setGroup(salt);
          sys+=mv.move( salt.size()/2+1 );
          break;
      }
      if (slp_global() > 0.9) {
        rdf1.sample( spc, salt, na, cl );
        rdf2.sample( spc, salt, ca, cl );
      }
    } // end of micro loop

    sys.checkDrift( Energy::systemEnergy(spc,pot,spc.p) );
    spc.save("state");
    cout << loop.timing();
  } // end of macro loop
  
#ifdef tab
  rdf1.save("naclrdf_tab.dat");
  rdf2.save("caclrdf_tab.dat");
#endif
#ifdef tabopt
  rdf1.save("naclrdf_tabopt.dat");
  rdf2.save("caclrdf_tabopt.dat");
#endif
#ifdef org
  rdf1.save("naclrdf_notab.dat");
  rdf2.save("caclrdf_notab.dat");
#endif

  pqr.save("confout.pqr", spc.p);

  cout << loop.info() << spc.info() << sys.info() << mv.info();
}
