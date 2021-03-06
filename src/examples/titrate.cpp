#include <faunus/faunus.h>

using namespace Faunus;

int main() {
  cout << textio::splash();            // faunus splash info!
  ::atom.includefile("titrate.json");  // load atom properties
  InputMap mcp("titrate.input");
  EnergyDrift sys;                     // class for tracking system energy drifts
  UnitTest test(mcp);

  Energy::Hamiltonian pot;
  pot.create( Energy::Nonbonded<Potential::DebyeHuckel,Geometry::Sphere>(mcp) );
  Space spc( pot.getGeometry() );

  // Add molecule to middle of simulation container
  FormatAAM aam;
  aam.load( mcp.get<string>("molecule", string()) );
  Geometry::cm2origo(*spc.geo, aam.particles() ); // center molecule
  GroupMolecular g = spc.insert( aam.particles() ); // insert into space
  g.name="peptide"; // babtise
  spc.enroll(g); // enroll group in space

  Move::SwapMove tit(mcp,pot,spc);
  Analysis::ChargeMultipole mp;

  spc.load("state");
  sys.init( Energy::systemEnergy(spc,pot,spc.p) );

  cout << atom.info() + spc.info() + pot.info() + tit.info()
    + textio::header("MC Simulation Begins!");

  MCLoop loop(mcp);            // class for handling mc loops
  while ( loop.macroCnt() ) {  // Markov chain 
    while ( loop.microCnt() ) {
      sys+=tit.move();
      mp.sample(g,spc);
    } // end of micro loop

    sys.checkDrift( Energy::systemEnergy(spc,pot,spc.p) );
    cout << loop.timing();

  } // end of macro loop

  cout << loop.info() + sys.info() + tit.info() + g.info() + mp.info();

  FormatPQR().save("confout.pqr", spc.p);
  spc.save("state");
}
