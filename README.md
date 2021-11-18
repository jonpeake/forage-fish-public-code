# Community dynamics of estuarine forage fishes are associated with a latitudinal basal resource regime
This is the collective data and R code for Jonathan Peake's work on forage fish community dynamics in the eastern Gulf of Mexico.

The data files Forage_Permanova and Forage_AbioticComps include data on samples constructed from raw seine and trawl samples by Estuary, Zone, Subzone, Year, and Season. Sample codes are a concatenation of these elements. Estuaries are noted in the column "Bay" and are abbreviated as follows:
  1) AP - Apalachicola Bay (note that bay codes are different from their abbreviations in the paper)
  2) CK - Cedar Key
  3) TB - Tampa Bay
  4) CH - Charlotte Harbor


These two data files note the Estuary (Bay), Stratification Zone, Subzone (Re_grid), Season, Number of Seines (nseines), and number of trawls (ntrawls) -
  1) Bay - AP, CK, TB, CH
  2) Zone - A, B, C, D, E, dependent on Bay
  2) Re_grid - 1-30, dependent on Bay and Zone
  3) Year - 1998-2017
  4) Season - 1 (Winter), 2 (Spring), 3 (Summer), 4 (Fall)
  5) nseines - Number of seines combined into sample, used to correct for effort
  6) ntrawls - Number of trawls combined into sample, used to correct for effort

These two data files include columns denoting unitless Mean Standardized Catch values for 50 different taxa of forage fishes. Forage fishes are denoted using the first 3 letters of the genus followed by the species. In the case of multi-species generic or familial complexes, genus or family names are spelled out.

The third data file, Forage_BasalResource, includes the basal resource associations of these taxa.

The data file Forage_AbioticComps additionally includes columns with local habitat values averaged across surveys by sample -
1) bottomveg - Proportion of bottom vegetation (seagrass, SAV) (percent cover)
2) temp - Temperature (degrees Celsius)
3) pH (unitless)
4) sal - Salinity (PSU)
5) DO - Dissolved Oxygen (milligrams per Liter)
6) vis - Vertical visibility estimated from Secchi disk depth (meters)

The R files include all code necessary to conduct the analyses presented in the manuscript, commented to indicate the usage of functions and code sections.
  
 For further information, please reach out to Jonathan Peake - jpeake1@usf.edu
