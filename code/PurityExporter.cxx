#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/SYSTEM/File.h>

#include <OpenMS/FORMAT/ConsensusXMLFile.h>
#include <OpenMS/FORMAT/TextFile.h>

#include <OpenMS/KERNEL/ConsensusMap.h>
#include <OpenMS/KERNEL/ConsensusFeature.h>

#include <OpenMS/ANALYSIS/QUANTITATION/ItraqConstants.h>

#include <vector>
#include <map>

using namespace OpenMS;
using namespace std;

class PurityExporter :
  public TOPPBase
{
public:
  PurityExporter() :
    TOPPBase("PurityExporter", "Exports ConsensusXML created by TMT/Itraq pipeline to tsv.", false, true)
  {
  }

protected:
  void registerOptionsAndFlags_()
  {
    registerInputFile_("in", "<input-file>", "", "The input file.",true);
    setValidFormats_("in",  ListUtils::create<String>("ConsensusXML"));

    registerOutputFile_("out", "<output-file>", "", "The output file.",true);
    setValidFormats_("out",  ListUtils::create<String>("tsv"));

    registerStringOption_("quant-type", "<ms2-type>", "ITRAQ4PLEX", "Type of MS2 quantification technique.",false);
    setValidStrings_("quant-type",  ListUtils::create<String>("TMT,ITRAQ4PLEX,ITRAQ8PLEX"));

    addEmptyLine_();
  }

  Param getSubsectionDefaults_(const String & /*section*/) const
  {
    Param tmp;
    return tmp;
  }

  std::vector<Int> getChannels(const String & quant_type)
  {
    if (quant_type == "TMT")
    {
      std::vector<Int> channles(ItraqConstants::CHANNELS_TMT_SIXPLEX[0], ItraqConstants::CHANNELS_TMT_SIXPLEX[0] + 6);
      return channles;
    }
    else if(quant_type == "ITRAQ4PLEX")
    {
      std::vector<Int> channles(ItraqConstants::CHANNELS_FOURPLEX[0], ItraqConstants::CHANNELS_FOURPLEX[0] + 4);
      return channles;
    }
    else // its 8plex
    {
      std::vector<Int> channles(ItraqConstants::CHANNELS_EIGHTPLEX[0], ItraqConstants::CHANNELS_EIGHTPLEX[0] + 8);
      return channles;
    }
  }
  
  template<typename ReturnType>
  ReturnType getValueOrDefault(const MetaInfoInterface& container, const String& key, ReturnType default_value)
  {
    if(container.metaValueExists(key))
    {
      return (ReturnType) container.getMetaValue(key);
    }
    else
    {
      return default_value;
    }
  }

  ExitCodes main_(int, const char **)
  {
    String quant_type = getStringOption_("quant-type");
    std::vector<Int> channels = getChannels(quant_type);

    ConsensusXMLFile file;
    ConsensusMap consensusMap;
    String filename = getStringOption_("in");
    file.load(filename, consensusMap);

    String tsvFileName = getStringOption_("out");

    TextFile textFile;

    // construct tsv file header
    StringList header;
    header.push_back("#scan_id");
    header.push_back("prec.rt");
    header.push_back("prec.mz");
    header.push_back("prec.intensity");
    header.push_back("charge");
    header.push_back("prec.purity");

    for(Size i = 0; i < channels.size(); ++i)
    {
      header.push_back(String("i") + String(channels[i]));
    }

    header.push_back("cf.intensity");

    textFile.push_back(ListUtils::concatenate(header, "\t"));

    for (ConsensusMap::Iterator cf = consensusMap.begin(); cf != consensusMap.end(); ++cf)
    {
      StringList currentLine;

      ConsensusFeature & cFeature = *cf;

      currentLine.push_back(getValueOrDefault(cFeature, "scan_id", String("")));
      currentLine.push_back(cFeature.getRT());
      currentLine.push_back(cFeature.getMZ());
      currentLine.push_back(getValueOrDefault(cFeature, "precursor_intensity", -1.0));
      currentLine.push_back(getValueOrDefault(cFeature, "precursor_charge", 0));
      currentLine.push_back(getValueOrDefault(cFeature, "precursor_purity", -1.0));
      
      std::map<Int, DoubleReal> intensityMap;


      // skip features with 0 intensity
      if(cf->getIntensity() == 0)
      {
        continue;
      }

      ConsensusFeature::HandleSetType features = cFeature.getFeatures();
      for(ConsensusFeature::HandleSetType::const_iterator fIt = features.begin();
          fIt != features.end();
          ++fIt)
      {
        intensityMap[Int(fIt->getMZ())] = fIt->getIntensity();
      }

      // export
      for (Size i = 0; i < channels.size(); ++i)
      {
        currentLine.push_back(String(intensityMap[channels[i]]));
      }

      currentLine.push_back(cFeature.getIntensity());
      textFile.push_back(ListUtils::concatenate(currentLine, "\t"));
    }

    textFile.store(tsvFileName);

    return EXECUTION_OK;
  }
};

int main(int argc, const char ** argv)
{
  PurityExporter tool;
  return tool.main(argc, argv);
}
