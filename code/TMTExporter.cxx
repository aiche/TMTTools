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

class TMTExporter :
  public TOPPBase
{
public:
  TMTExporter() :
    TOPPBase("TMTExporter", "Exports ConsensusXML created by TMT/Itraq pipeline to tsv.", false, true)
  {
  }

protected:
  void registerOptionsAndFlags_()
  {
    registerInputFile_("in", "<input-file>", "", "The input file.",true);
    setValidFormats_("in", StringList::create("ConsensusXML"));

    registerOutputFile_("out", "<output-file>", "", "The output file.",true);
    setValidFormats_("out", StringList::create("tsv"));
    
    registerStringOption_("quant-type", "<ms2-type>", "TMT", "Type of MS2 quantification technique.",false);
    setValidStrings_("quant-type", StringList::create("TMT,ITRAQ4PLEX,ITRAQ8PLEX"));
    
    registerFlag_("no-prot", "Do not check for protein identifications.");
    
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

  ExitCodes main_(int, const char **)
  {
    String quant_type = getStringOption_("quant-type");
    std::vector<Int> channels = getChannels(quant_type);
    
    ConsensusXMLFile file;
    ConsensusMap consensusMap;
    String filename = getStringOption_("in");
    file.load(filename, consensusMap);

    String tsvFileName = getStringOption_("out");
    
    bool useProt = !getFlag_("no-prot");
    
    TextFile textFile;

    // construct tsv file header
    StringList header;
    header.push_back("#sequence");
    header.push_back("prot.name");
    header.push_back("prot.accession");
    
    header.push_back("spot.nr");
    
    for(Size i = 0; i < channels.size(); ++i)
    {
      header.push_back(String("i") + String(channels[i]));
    }
    
    header.push_back("cf_id");
    header.push_back("MS2tic");
    header.push_back("charge");
    header.push_back("RT");
    
    textFile.push_back(header.concatenate("\t"));
    
    //"sequence,prot.name,prot.accession,spot.nr,i114,i115,i116,i117,cf_id,
    //overlap,pre_mono_contribution,tic,MS2tic,charge,RT"
    
    typedef std::vector< ProteinHit >::iterator ProtHitIt;
    ProteinIdentification protIdent;
    
    if (useProt)
    {
      protIdent = consensusMap.getProteinIdentifications()[0];
    }
    
    Int noId(0);
    Int notUnique(0);
    Int noSignal(0);
    Int remaining(0);
    
    for (ConsensusMap::Iterator cf = consensusMap.begin(); cf != consensusMap.end(); ++cf)
    {
      StringList currentLine;
      
      ConsensusFeature & cFeature = *cf;
      
      Int charge(0);
      
      if (cFeature.getPeptideIdentifications().size() == 0 || !useProt)
      {
        // we store unidentified hits anyway, because the iTRAQ quant is still helpful for normalization
        ++noId;
        currentLine.push_back("UNIDENTIFIED_PEPTIDE");
        currentLine.push_back("UNIDENTIFIED_PROTEIN");
        currentLine.push_back("UNIDENTIFIED_PROTEIN");
      }
      else
      {
        charge = cFeature.getPeptideIdentifications()[0].getHits()[0].getCharge();
        currentLine.push_back( cFeature.getPeptideIdentifications()[0].getHits()[0].getSequence().toUnmodifiedString());
        
        // protein name:
        if (cFeature.getPeptideIdentifications()[0].getHits()[0].getProteinAccessions ().size() != 1)
        {
          ++notUnique;
          continue; // we only want unique peptides
        }

        ProtHitIt proteinHit = protIdent.findHit (cFeature.getPeptideIdentifications()[0].getHits()[0].getProteinAccessions()[0]);
        if (proteinHit == protIdent.getHits().end())
        {
          std::cerr << "Protein referenced in peptide not found...\n";
          continue; // protein not found
        }
        
        currentLine.push_back( proteinHit->getAccession() );
        
        // protein accession
        currentLine.push_back( proteinHit->getAccession() );
      }
      
      // spot number (always constant here)
      currentLine.push_back("1");
      
      std::map<Int, DoubleReal> intensityMap;
      
      
      // skip features with 0 intensity
      if(cf->getIntensity() == 0)
      {
        ++noSignal;
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

      currentLine.push_back(String(cFeature.getUniqueId()));
      
      currentLine.push_back((String(cFeature.getIntensity())));
      currentLine.push_back(String(charge));
      currentLine.push_back(String(cFeature.getRT()));
      textFile.push_back(currentLine.concatenate("\t"));
      ++remaining;
    }
    
    textFile.store(tsvFileName);
    

    std::cout << "Useless: " << noId << " peptides (no ID - but kept for normalization)\n"
      << "Removed: " << notUnique << " peptides (not unique)\n"
      << "         " << noSignal << " peptides (no signal)\n";
    std::cout << "Remaining " << remaining << " peptides\n";
    
    return EXECUTION_OK;
  }
};

int main(int argc, const char ** argv)
{
  TMTExporter tool;
  return tool.main(argc, argv);
}
