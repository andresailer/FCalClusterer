#ifndef FCalAna_hh
#define FCalAna_hh 1

#include <marlin/Processor.h>
#include <lcio.h>


#include <string>
#include <vector>


class TH1D;
class TH2D;
class TH3D;
class TFile;
class TRandom3;
class TTree;

class FCalAna : public marlin::Processor {
  
 public:
  
  virtual Processor*  newProcessor() { return new FCalAna ; }
  
  
  FCalAna() ;

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
  std::string m_nameOutputFile, m_mcParticleName, m_recoParticleName;

  int m_nRun ;
  int m_nEvt ;

  TFile* m_file;
  TTree* m_tree;

  double m_thetaMC, m_thetaReco, m_eMC, m_eReco, m_phiMC, m_phiReco;


private://to shut the warnings up
  FCalAna(const FCalAna&);
  FCalAna& operator=(const FCalAna&);

} ;

#endif // FCalAna_hh
