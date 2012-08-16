#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/FORMAT/FASTAFile.h>

#include <OpenMS/CHEMISTRY/EnzymaticDigestion.h>

#include <vector>
#include <map>

using namespace OpenMS;
using namespace std;

class TMTExporter :
  public TOPPBase
{
public:
  IsotopeCalculator() :
    TOPPBase("TMTExporter", "Exports ConsensusXML created by TMT pipeline to tsv.", false, true)
  {
  }

protected:
  void registerOptionsAndFlags_()
  {
    registerStringOption_("in", "<input-file>", "", "The input file.",true);
    addEmptyLine_();
  }

  Param getSubsectionDefaults_(const String & /*section*/) const
  {
    Param tmp;
    return tmp;
  }

  ExitCodes main_(int, const char **)
  {
    return EXECUTION_OK;
  }
};

int main(int argc, const char ** argv)
{
  IsotopeCalculator tool;
  return tool.main(argc, argv);
}
