#ifndef ReadBeamCal_h
#define ReadBeamCal_h 1

#include <string>

#include <marlin/Processor.h>

class BCPadEnergies;
class BeamCalGeo;

class TRandom3;
class TFile;
class TTree;

namespace EVENT {
  class LCEvent;
  class LCRunHeader;
}

class ReadBeamCal : public marlin::Processor {
  
 public:
  
  virtual Processor*  newProcessor() { return new ReadBeamCal ; }
  
  
  ReadBeamCal() ;

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

  /** Input collection name.
   */
  std::string m_colNameBCal, m_nameOutputFile, m_nameFinalOutputFile, m_nameInputFile;
  int m_startingLayer=1;
  std::string m_detectorName = "BeamCal";

  //std::string m_pdftitle;
  int m_nRun=0;
  int m_nEvt=0;
  double m_probFactor;

  TRandom3 *m_random3;
  TFile *m_rootfile=nullptr;
  TTree *m_tree=nullptr;

  BCPadEnergies* m_padEnergiesLeft;
  BCPadEnergies* m_padEnergiesRight;

  BeamCalGeo* m_bcg;
  bool m_usingDD4HEP;
  bool m_writeEachEvent=false;
  int m_rootFileCounter = 0;
  double m_eventsPerBX = -1;
  int m_eventsToCount = 0;
  double m_averageNumber = 0;
private://to shut the warnings up
  ReadBeamCal(const ReadBeamCal&);
  ReadBeamCal& operator=(const ReadBeamCal&);

} ;


#endif



