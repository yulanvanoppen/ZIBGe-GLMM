# Aeshna viridis population data

_Aeshna viridis population measures and environmental factors collected from 2015–2018 at 17 locations (ditches or ponds) in the northern Netherlands by [1]._

## Contents
The data in `/populationdata.csv` contain 114 observations of 15 variables (yearly aggregates; totals and averages):

- `manager`: the manager responsible for maintaining and surveilling the location, labeled _GL_ (_Gemeente Leeuwarden_), _GV_ (_Gemeente Veendam_), _SBBF_ (_Staatsbosbeheer Twijzel_), _SBBG_ (_Staatsbosbeheer Groningen_), or _WHA_ (_Waterschap Hunze en Aa's_);
- `area`: the location with labels composed of the manager label and a counter;
- `transect`: one of two parts in which the location is divided, labeled 'a' or 'b';
- `treatment`: either the 'regular' or the 'checker' treatment, randomly allocated to each location's part 'a' or 'b' (same for each year);
- `year`: the year in which the (aggregate) observation was made;
- `exuviae`: the total number of _exuviae_ (exoskeletons shed by the _A. viridis_ upon transitioning from the larval to the adult stage) counted;
- `egglaying_females`: the total number of egg-laying females counted.

The remaining variables relate to the water quality (yearly averages):

- `emers_frac`: the surface coverage fraction of _A. viridis_' host plant _Stratiotes aloides_.
- `pH`: the acidity measured on the pH scale;
- `redox_V`: the redox measured in volts;
- `O2_frac`: the oxygen level measured in the dissolved oxygen fraction;
- `EC_mS_cm`: the electrical conductivity measured in millisiemens per centimeter;
- `temp_C`: the temperatere measured in degrees Celcius;
- `water_depth_m`: the water depth measured in meters;
- `sludge_thickness_m`: the thickness of the sludge layer formed by decayed _S. aloides_ measured in meters.


## References
[1] Milder-Mulderij, G. and Brochard, C. and Wiggers, R. and de Vries, S. (2020). Alternatief krabbenscheerbeheer in Fryslân, Groningen en Drenthe.  Interpretatie op basis van vier jaar onderzoek op diverse locaties. Bureau Biota.
