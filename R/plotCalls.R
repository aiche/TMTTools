# iTRAQ 4plex 1:1
plotRatios(filenames_=c("120628_001_HeLa_1-1-1-101.tsv", "SR030712_0016_HeLa_1-1-1-1_EL1.tsv", "SR030712_0026_HeLa_1-1-1-1_EL2.tsv"), title_="iTRAQ-4plex 1:1:1:1", ref_="i114", channels=c("i114", "i115", "i116", "i117"), ylim_=c(-2,2), expectedRatios_=c(1))

# iTRAQ 4plex mixed
plotRatios(filenames_=c("SR030712_0010.tsv", "SR030712_0014_HeLa_50-5-25-1_EL1.tsv", "SR030712_0020_HeLa_50-5-25-1_EL2.tsv"), title_="iTRAQ-4plex 50:5:25:1", ref_="i117", channels=c("i114", "i115", "i116", "i117"), ylim_=c(-2,6), expectedRatios_=c(50,5,25))

################################################################################################################

# iTRAQ 8plex 1:1
plotRatios(filenames_=c("SR030712_0028_HeLa_1-1-1-1-1-1_8-PLEX.tsv", "SR030712_0037_HeLa_1-1-1-1-1-1_8-PLEX_EL1.tsv", "SR030712_0045_HeLa_1-1-1-1-1-1_8-PLEX_EL2.tsv"), title_="iTRAQ-8plex 1:1:1:1:1:1", ref_="i113", channels=c("i113", "i114", "i115", "i116", "i117", "i118"), ylim_=c(-2,2), expectedRatios_=c(1))

# iTRAQ 8plex mixed
plotRatios(filenames_=c("SR030712_0030_HeLa_50-5-25-1-25-5_8-PLEX.tsv", "SR030712_0039_HeLa_50-5-25-1-25-5_8-PLEX_EL1.tsv", "SR030712_0047_HeLa_50-5-25-1-25-5_8-PLEX_EL2.tsv"), title_="iTRAQ-8plex 50:5:25:1:25:5", ref_="i116", channels=c("i113", "i114", "i115", "i116", "i117", "i118"), ylim_=c(-2,6), expectedRatios_=c(50,5,25))

################################################################################################################

# TMT 1:1
plotRatios(filenames_=c("SR030712_0022_HeLa_1-1-1-1-1-1.tsv", "SR030712_0033_HeLa_1-1-1-1-1-1_TMT_EL1_no2.tsv", "SR030712_0041_HeLa_1-1-1-1-1-1_TMT_EL2.tsv"), title_="TMT-6plex 1:1:1:1:1:1",  ref_="i126", channels=c("i126", "i127", "i128", "i129", "i130", "i131"), ylim_=c(-2,2), expectedRatios_=c(1))

# TMT mixed
plotRatios(filenames_=c("SR030712_0024_HeLa_50-5-25-1-25-5_TMT.tsv", "SR030712_0035_HeLa_50-5-25-1-25-5_TMT_EL1.tsv", "SR030712_0043_HeLa_50-5-25-1-25-5_TMT_EL2.tsv"), title_="TMT-6plex 50:5:25:1:25:5", ref_="i129", channels=c("i126", "i127", "i128", "i129", "i130", "i131"), ylim_=c(-2,6), expectedRatios_=c(50,5,25))
