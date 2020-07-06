// Martin Goettlich @ DESY
//
 
#ifndef EventAction_h
#define EventAction_h 1

#include "G4UserEventAction.hh"
#include "globals.hh"
#include <vector>

const G4int kDim = 2;


class G4Event;

class EventAction : public G4UserEventAction
{
public:
  EventAction(const G4int& modulo);
  ~EventAction();
  
  void BeginOfEventAction(const G4Event*);
  void EndOfEventAction(const G4Event*);
  int GetNPhotons(const G4Event*) ;
 
private:
  G4int printModulo;
//  G4int fDRHCID_f;
  G4int fDRHCID;
};


#endif
