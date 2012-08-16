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
    TOPPBase("TMTExporter", "Exports ConsensusXML created by TMT pipeline to tsv.", false, true)
  {
  }

protected:
  void registerOptionsAndFlags_()
  {
    registerInputFile_("in", "<input-file>", "", "The input file.",true);
    setValidFormats_("in", StringList::create("ConsensusXML"));

    registerOutputFile_("out", "<output-file>", "", "The output file.",true);
    setValidFormats_("out", StringList::create("tsv"));
    
    addEmptyLine_();
  }

  Param getSubsectionDefaults_(const String & /*section*/) const
  {
    Param tmp;
    return tmp;
  }

  ExitCodes main_(int, const char **)
  {
    std::vector<Int> tmtChannels(ItraqConstants::CHANNELS_TMT_SIXPLEX[0], ItraqConstants::CHANNELS_TMT_SIXPLEX[0] + 6);
    
    ConsensusXMLFile file;
    ConsensusMap consensusMap;
    String filename = getStringOption_("in");
    file.load(filename, consensusMap);

    String tsvFileName = getStringOption_("out");
    
    TextFile textFile;

    // construct tsv file header
    StringList header;
    header.push_back("sequence");
    header.push_back("prot.name");
    header.push_back("prot.accession");
    
    header.push_back("spot.nr");
    
    for(Size i = 0; i < tmtChannels.size(); ++i)
    {
      header.push_back(String("i") + String(tmtChannels[i]));
    }
    
    header.push_back("cf_id");
    header.push_back("MS2tic");
    header.push_back("charge");
    header.push_back("RT");
    
    textFile.push_back(header.concatenate("\t"));
    
    //"sequence,prot.name,prot.accession,spot.nr,i114,i115,i116,i117,cf_id,
    //overlap,pre_mono_contribution,tic,MS2tic,charge,RT"
    
    ProteinIdentification protIdent = consensusMap.getProteinIdentifications()[0];
    typedef std::vector< ProteinHit >::iterator ProtHitIt;
    
    Int noId(0);
    Int notUnique(0);
    Int noSignal(0);
    Int remaining(0);
    
    for (ConsensusMap::Iterator cf = consensusMap.begin(); cf != consensusMap.end(); ++cf)
    {
      StringList currentLine;
      
      ConsensusFeature & cFeature = *cf;
      
      Int charge(0);
      
      if (cFeature.getPeptideIdentifications().size() == 0)
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
      for (Size i = 0; i < tmtChannels.size(); ++i)
      {
        currentLine.push_back(String(intensityMap[tmtChannels[i]]));
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
