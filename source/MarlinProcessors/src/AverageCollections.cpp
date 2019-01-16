#include "AverageCollections.hh"

#include <EVENT/LCCollection.h>
#include <EVENT/LCEvent.h>
#include <EVENT/LCIO.h>
#include <EVENT/SimCalorimeterHit.h>
#include <EVENT/SimTrackerHit.h>

#include <Exceptions.h>

// ----- include for verbosity dependend logging ---------
#include <streamlog/baselevels.h>
#include <streamlog/loglevels.h>
#include <streamlog/logstream.h>
#include <streamlog/streamlog.h>

#include <algorithm>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <vector>

using namespace lcio;
using namespace marlin;

AverageCollections aAverageCollections;

// #pragma GCC diagnostic push
// #pragma GCC diagnostic ignored "-Wshadow"

AverageCollections::AverageCollections() : Processor("AverageCollections") {
  _description = "AverageCollections reads the collections and returns an average energy deposit for each collection";

  registerInputCollections(LCIO::TRACK, "Collections", "Names of the Collections", m_collections, {});

  registerProcessorParameter("EventsPerBX", "Number of events per BX", m_eventsPerBX, m_eventsPerBX);
  registerProcessorParameter("EventsPerFile", "Number of events per File", m_eventsPerFile, m_eventsPerFile);
}
//#pragma GCC diagnostic pop

void AverageCollections::init() {
  printParameters();

  for (size_t index = 0; index < m_collections.size(); ++index) {
    m_totalEnergy.push_back(0);

    m_maxLength = std::max(int(m_collections[index].size()), m_maxLength);
  }

}  //init

void AverageCollections::processRunHeader(LCRunHeader*) { m_nRuns += 1; }

void AverageCollections::processEvent(LCEvent* evt) {
  m_nEvents += 1;

  for (size_t index = 0; index < m_collections.size(); ++index) {
    auto const colName = m_collections[index];
    LCCollection* collection;

    try {
      collection = evt->getCollection(colName);
    } catch (Exception& e) {
      collection = 0;
    }

    if (not collection)
      continue;

    for (int i = 0; i < collection->getNumberOfElements(); ++i) {
      double energy(0.0);
      SimCalorimeterHit* cHit = dynamic_cast<SimCalorimeterHit*>(collection->getElementAt(i));
      if (cHit) {
        energy = cHit->getEnergy();
        // for(int iCont = 0; iCont < cHit->getNMCContributions(); ++iCont){
        //   if(cHit->getTimeCont(iCont) < 300){
        //     energy += cHit->getEnergyCont(iCont);
        //   }
        // }
      } else {
        SimTrackerHit* tHit = dynamic_cast<SimTrackerHit*>(collection->getElementAt(i));
        if (tHit) {
          energy = tHit->getEDep();
        }
      }
      m_totalEnergy[index] += energy;
    }  //for all entries in the collection
  }    // for all collections
  return;
}  //processEvent

void AverageCollections::check(LCEvent*) {
  // nothing to check here - could be used to fill checkplots in reconstruction processor
}

void AverageCollections::end() {
  streamlog_out(DEBUG) << "\nNumber of events " << m_nEvents << "\nNumber of Files " << m_nRuns << "\nEvents Per BX "
                       << m_eventsPerBX << "\nEvents per File " << m_eventsPerFile << "\nScaled with "
                       << m_eventsPerBX / (m_nRuns * m_eventsPerFile) << std::endl;

  if (streamlog::out.write<DEBUG>()) {
    std::cout << "Debug only...." << std::endl;

  }  //Only draw all the things in DEBUG Mode

  for (size_t index = 0; index < m_collections.size(); ++index) {
    streamlog_out(MESSAGE) << std::setw(m_maxLength + 4) << m_collections[index] << std::setw(14)
                           << m_totalEnergy[index] / (m_nRuns * m_eventsPerFile) * m_eventsPerBX << std::endl;

  }  //for all collections
}
