#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/SYSTEM/File.h>

#include <OpenMS/FORMAT/ConsensusXMLFile.h>
#include <OpenMS/FORMAT/TextFile.h>

#include <OpenMS/KERNEL/ConsensusMap.h>
#include <OpenMS/KERNEL/ConsensusFeature.h>

#include <OpenMS/ANALYSIS/QUANTITATION/IsobaricQuantitationMethod.h>
#include <OpenMS/ANALYSIS/QUANTITATION/ItraqFourPlexQuantitationMethod.h>
#include <OpenMS/ANALYSIS/QUANTITATION/ItraqEightPlexQuantitationMethod.h>
#include <OpenMS/ANALYSIS/QUANTITATION/TMTSixPlexQuantitationMethod.h>

#include <vector>
#include <map>

using namespace OpenMS;
using namespace std;

class IBSpectraExporter :
  public TOPPBase
{
private:
  /// available quantitiation methods
  std::map<String, IsobaricQuantitationMethod*> quant_methods_;
  /// user selected quantitation method
  IsobaricQuantitationMethod * selected_quant_method_;
  /// only unique protein matches are allowed
  bool allow_non_unique_;
  ///
  bool use_prot_;
  
  /**
    @brief Holds all id information contained in a id csv line
  */
  struct IdCSV
  {
    String accession; // ￼ Protein AC
    String peptide; // ￼ Peptide sequence
    String modif; //￼ Peptide modification string
    Int charge; // ￼ Charge state
    DoubleReal theo_mass; // Theoretical peptide mass
    DoubleReal exp_mass; //￼ Experimentally observed mass
    DoubleReal parent_intens; //￼ Parent intensity
    DoubleReal retention_time; // Retention time
    String spectrum; // ￼ Spectrum identifier
    String search_engine; // ￼ Protein search engine and score
    
    void toStringList(StringList & target_list)
    {
      target_list.push_back(accession);
      target_list.push_back(peptide);
      target_list.push_back(modif);
      target_list.push_back(charge);
      target_list.push_back(theo_mass);
      target_list.push_back(exp_mass);
      target_list.push_back(parent_intens);
      target_list.push_back(retention_time);
      target_list.push_back(spectrum);
      target_list.push_back(search_engine);
    }
    
    IdCSV()
    : accession("UNIDENTIFIED_PROTEIN"),
      peptide("UNIDENTIFIED_PEPTIDE"),
      modif(""),
      charge(0),
      theo_mass(-1.),
      exp_mass(-1.),
      parent_intens(-1.),
      retention_time(-1.),
      spectrum(""),
      search_engine("open-ms-generic")
    {}
  };
  
public:
  IBSpectraExporter() :
    TOPPBase("IBSpectra", "Exports ConsensusXML created by TMT/Itraq pipeline to ibspectra format that can directly be loaded into isobar.", false, true)
  {
    ItraqFourPlexQuantitationMethod* itraq4plex = new ItraqFourPlexQuantitationMethod();
    ItraqEightPlexQuantitationMethod* itraq8plex = new ItraqEightPlexQuantitationMethod();
    TMTSixPlexQuantitationMethod* tmt6plex = new TMTSixPlexQuantitationMethod();
    quant_methods_[itraq4plex->getName()] = itraq4plex;
    quant_methods_[itraq8plex->getName()] = itraq8plex;
    quant_methods_[tmt6plex->getName()] = tmt6plex;
    
    selected_quant_method_ = NULL;
  }
  
protected:
  void registerOptionsAndFlags_()
  {
    registerInputFile_("in", "<input-file>", "", "The input files.",true);
    setValidFormats_("in", StringList::create("consensusXML"));

    registerOutputFile_("out", "<output-file>", "", "The output file.",true);
    // setValidFormats_("out", StringList::create("tsv"));
    
    registerStringOption_("type", "<mode>", quant_methods_.begin()->first, "Isobaric Quantitation method used in the experiment.", false);
    StringList valid_types;
    for (std::map<String, IsobaricQuantitationMethod*>::iterator it = quant_methods_.begin();
         it != quant_methods_.end();
         ++it)
    {
      valid_types.push_back(it->first);
    }
    setValidStrings_("type", valid_types);
    
    registerFlag_("no-prot", "Do not check for protein identifications.");
    registerFlag_("allow-non-unique", "Allow also non-unique peptides to be exported.");
    
    addEmptyLine_();     
  }
  
  StringList construct_header()
  {
    // construct tsv file header
    StringList header;
    
    
    header.push_back("accession"); // ￼ Protein AC 
    header.push_back("peptide"); // ￼ Peptide sequence
    header.push_back("modif"); // ￼ Peptide modification string 
    header.push_back("charge"); // ￼ Charge state
    header.push_back("theo.mass"); // ￼ Theoretical peptide mass 
    header.push_back("exp.mass"); // ￼ Experimentally observed mass
    header.push_back("parent.intens"); // ￼ Parent intensity
    header.push_back("retention.time"); // ￼ Retention time
    header.push_back("spectrum"); // ￼ Spectrum identifier 
    header.push_back("search.engine"); // ￼ Protein search engine and score
      
    for(IsobaricQuantitationMethod::IsobaricChannelList::const_iterator it = selected_quant_method_->getChannelInformation().begin();
        it != selected_quant_method_->getChannelInformation().end();
        ++it)
    {
      header.push_back("X" + String(int(it->center)) +  "_mass");
    }

    for(IsobaricQuantitationMethod::IsobaricChannelList::const_iterator it = selected_quant_method_->getChannelInformation().begin();
        it != selected_quant_method_->getChannelInformation().end();
        ++it)
    {
      header.push_back("X" + String(int(it->center)) +  "_ions");
    }

    return header;
  }
  
  std::vector<IdCSV> createLinesFromConsensusFeature(const ConsensusFeature & cf)
  {
    std::vector<IdCSV> entries;
    
    // no id information available or not requested
    if (cf.getPeptideIdentifications().size() == 0 || !use_prot_)
    {
      // just add a dummy entry
      entries.push_back(IdCSV());
    }
    else
    {
      
    }
    
    return entries;
  }
  
  /**  
    accession ￼ Protein AC 
    peptide ￼ Peptide sequence
    modif ￼ Peptide modification string 
    charge ￼ Charge state
    theo.mass ￼ Theoretical peptide mass 
    exp.mass ￼ Experimentally observed mass
    parent.intens ￼ Parent intensity
    retention.time ￼ Retention time
    spectrum ￼ Spectrum identifier 
    search.engine ￼ Protein search engine and score
    ---
    X114_mass ￼ reporter ion mass 
    X115_mass ￼ reporter ion mass 
    X116_mass ￼ reporter ion mass 
    X117_mass ￼ reporter ion mass
    X114_ions ￼ reporter ion intensity 
    X115_ions ￼ reporter ion intensity
    X116_ions ￼ reporter ion intensity 
    X117_ions ￼ reporter ion intensity
  */
  
  ExitCodes main_(int, const char **)
  {
    typedef std::vector< ProteinHit >::iterator ProtHitIt;

    //-------------------------------------------------------------
    // init quant method
    //-------------------------------------------------------------
    selected_quant_method_ = quant_methods_[getStringOption_("type")];
    
    // load input files
    ConsensusXMLFile file;
    ConsensusMap consensusMap;
    String filename = getStringOption_("in");
    file.load(filename, consensusMap);
    
    // get prot flag
    use_prot_ = !getFlag_("no-prot");
    allow_non_unique_ = getFlag_("allow-non-unique");
    
    // check if output file is valid
    String tsvFileName = getStringOption_("out");
    TextFile textFile;
    
    textFile.push_back(construct_header().concatenate("\t"));
    // std::cerr << *textFile.begin() << std::endl;
    
    // stats counter
    Int noId(0);
    Int notUnique(0);
    Int noSignal(0);
    Int remaining(0);
    
    // protein identifications
    ProteinIdentification protIdent;
    if (use_prot_ && consensusMap.getProteinIdentifications().size() > 0)
    {
      protIdent = consensusMap.getProteinIdentifications()[0];
    }
    
    for(ConsensusMap::Iterator cm_iter = consensusMap.begin();
        cm_iter != consensusMap.end();
        ++cm_iter)
    {
      ConsensusFeature & cFeature = *cm_iter;
      std::vector<IdCSV> entries;
      
      if (cFeature.getPeptideIdentifications().size() == 0 || !use_prot_)
      {
        // we store unidentified hits anyway, because the iTRAQ quant is still helpful for normalization
        ++noId;
        entries.push_back(IdCSV());
      }
      else
      {
        // protein name:
        if (cFeature.getPeptideIdentifications()[0].getHits()[0].getProteinAccessions().size() != 1)
        {
          ++notUnique;
          if(!allow_non_unique_)
            continue; // we only want unique peptides
        }
        
        for (std::vector<String>::const_iterator prot_ac = cFeature.getPeptideIdentifications()[0].getHits()[0].getProteinAccessions().begin();
             prot_ac != cFeature.getPeptideIdentifications()[0].getHits()[0].getProteinAccessions().end();
             ++prot_ac)
        {
          IdCSV entry;
          entry.charge = cFeature.getPeptideIdentifications()[0].getHits()[0].getCharge();
          entry.peptide = cFeature.getPeptideIdentifications()[0].getHits()[0].getSequence().toUnmodifiedString();
          entry.theo_mass = cFeature.getPeptideIdentifications()[0].getHits()[0].getSequence().getMonoWeight(Residue::Full, cFeature.getPeptideIdentifications()[0].getHits()[0].getCharge());
          
          // write modif
          entry.modif = cFeature.getPeptideIdentifications()[0].getHits()[0].getSequence().getNTerminalModification();
          for(AASequence::ConstIterator aa_it = cFeature.getPeptideIdentifications()[0].getHits()[0].getSequence().begin();
              aa_it != cFeature.getPeptideIdentifications()[0].getHits()[0].getSequence().end();
              ++aa_it)
          {
            entry.modif += ":" + aa_it->getModification();
          }
          entry.modif += ":" + cFeature.getPeptideIdentifications()[0].getHits()[0].getSequence().getCTerminalModification();
          
          
          ProtHitIt proteinHit = protIdent.findHit(*prot_ac);
          if (proteinHit == protIdent.getHits().end())
          {
            std::cerr << "Protein referenced in peptide not found...\n";
            continue; // protein not found
          }

          entry.accession = proteinHit->getAccession();
          entries.push_back(entry);
        }
      }
      
      // skip features with 0 intensity
      if(cFeature.getIntensity() == 0)
      {
        ++noSignal;
        continue;
      }

      for (std::vector<IdCSV>::iterator entry = entries.begin();
           entry != entries.end();
           ++entry)
      {
        // set parent intensity
        entry->parent_intens = cFeature.getIntensity();
        entry->retention_time = cFeature.getRT();
        entry->spectrum = cFeature.getUniqueId();
        entry->exp_mass = cFeature.getMZ();
        
        // create output line
        StringList currentLine;
        
        // add entry to currentLine
        entry->toStringList(currentLine);
        
        // extract channel intensities and positions
        std::map<Int, DoubleReal> intensityMap;
        ConsensusFeature::HandleSetType features = cFeature.getFeatures();
        
        for(ConsensusFeature::HandleSetType::const_iterator fIt = features.begin();
            fIt != features.end();
            ++fIt)
        {
          intensityMap[Int(fIt->getMZ())] = fIt->getIntensity();
        }
        for(IsobaricQuantitationMethod::IsobaricChannelList::const_iterator it = selected_quant_method_->getChannelInformation().begin();
            it != selected_quant_method_->getChannelInformation().end();
            ++it)
        {
          currentLine.push_back(String(it->center));
        }
        for(IsobaricQuantitationMethod::IsobaricChannelList::const_iterator it = selected_quant_method_->getChannelInformation().begin();
            it != selected_quant_method_->getChannelInformation().end();
            ++it)
        {
          currentLine.push_back(String(intensityMap[int(it->center)]));
        }
        
        textFile.push_back(currentLine.concatenate("\t"));
        ++remaining;
      }
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
  IBSpectraExporter tool;
  return tool.main(argc, argv);
}
