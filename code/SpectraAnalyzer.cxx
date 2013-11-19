#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/SYSTEM/File.h>

#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/MSSpectrum.h>

#include <OpenMS/FORMAT/MzMLFile.h>

#include <vector>
#include <map>

using namespace OpenMS;
using namespace std;

class SpectraAnalyzer :
public TOPPBase
{
public:
  SpectraAnalyzer() :
    TOPPBase("TMTExporter", "Analyze the type of spectra (HCID, CID, analyzer type).", false, true)
  {
  }
  
protected:
  void registerOptionsAndFlags_()
  {
    registerInputFile_("in", "<input-file>", "", "The input file.",true);
    setValidFormats_("in", StringList::create("mzML"));
        
    addEmptyLine_();
  }
  
  Param getSubsectionDefaults_(const String & /*section*/) const
  {
    Param tmp;
    return tmp;
  }
   
  ExitCodes main_(int, const char **)
  {
    String filename = getStringOption_("in");
    
    MzMLFile mzMLFile;
    MSExperiment<> exp;
    
    // we only want to load the MS2 spectra
    IntList msLevels = IntList::create("2");
    mzMLFile.getOptions().setMSLevels(msLevels);
    
    mzMLFile.load(filename, exp);
    
    for(MSExperiment<>::Iterator specIt = exp.begin(); specIt != exp.end(); ++specIt)
    {
      MSSpectrum<> & spec = *specIt;
      InstrumentSettings instSettings = spec.getInstrumentSettings();

    }
    
    
    
    return EXECUTION_OK;
  }
};

int main(int argc, const char ** argv)
{
  SpectraAnalyzer tool;
  return tool.main(argc, argv);
}