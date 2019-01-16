#ifndef AverageCollections_h
#define AverageCollections_h 1

#include <marlin/Processor.h>

#include <string>

namespace EVENT {
  class LCEvent;
  class LCRunHeader;
}

class AverageCollections: public marlin::Processor {
  
 public:
  
  virtual Processor*  newProcessor() { return new AverageCollections ; }
  
  AverageCollections() ;
  AverageCollections(const AverageCollections&) = delete;
  AverageCollections& operator=(const AverageCollections&) = delete;

  /** Called at the begin of the job before anything is read.
   * Use to initialize the processor, e.g. book histograms.
   */
  virtual void init() ;
  
  /** Called for every run.
   */
  virtual void processRunHeader( LCRunHeader* run ) ;
  
  /** Called for every event - the working horse.
   */
  virtual void processEvent( LCEvent * evt ) ; 
  
  
  virtual void check( LCEvent * evt ) ; 
  
  
  /** Called after data processing for clean up.
   */
  virtual void end() ;
  
  
 protected:

  std::vector<std::string> m_collections{};
  std::vector<double> m_totalEnergy{};

  int m_nEvents = 0;
  int m_maxLength = 0;
  int m_nRuns = 0;
  double m_eventsPerBX = 0.0;
  double m_eventsPerFile = 0.0;
};


#endif



