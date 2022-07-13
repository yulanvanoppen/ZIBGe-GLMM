# Aeshna viridis population data

_Aeshna viridis population measures and environmental factors at 17 locations (ditches or ponds) in the Northern Netherlands during 2015--2018, collected and owned by [Bureau Biota (Groningen, NL)](https://www.bureaubiota.com "https://www.bureaubiota.com") [1]_

<img src="https://github.com/yulanvanoppen/storage/blob/main/Aeshna_viridis_with_exuviae.jpg" width="400">

**Figure 1.** A female green hawker (left) and exuviae (right). Photo copyright by [Bureau Biota (Groningen, NL)](https://www.bureaubiota.com "https://www.bureaubiota.com").

&nbsp;

## Contents
The data in `/populationdata.csv` contain 114 observations of 15 variables (yearly aggregates; totals and averages):

- `manager`: the manager responsible for maintaining and surveilling the location, labeled _GL_ (_Groninger Landschap_), _GV_ (_Gemeente Veendam_), _SBBF_ (_Staatsbosbeheer Twijzel_), _SBBG_ (_Staatsbosbeheer Groningen_), or _WHA_ (_Waterschap Hunze en Aa's_);
- `area`: the location with labels composed of the manager label and a counter;
- `transect`: one of two parts in which the location is divided, labeled 'a' or 'b';
- `treatment`: either the 'regular' or the 'checker' treatment (host plant removal procedure), randomly allocated to each location's part 'a' or 'b' (same for each year);
- `year`: the year in which the (aggregate) observation was made;
- `exuviae`: the total number of _exuviae_ (exoskeletons shed by the _A. viridis_ upon transitioning from the larval to the adult stage) counted;
- `egglaying_females`: the total number of egg-laying females counted.

The remaining variables relate to the water quality (yearly averages):

- `emers_frac`: the surface coverage fraction of _A. viridis_' host plant _Stratiotes aloides_.
- `pH`: the acidity measured on the pH scale;
- `redox_V`: the redox potential measured in volts;
- `O2_frac`: the oxygen concentration measured in the dissolved oxygen fraction;
- `EC_mS_cm`: the electrical conductivity measured in millisiemens per centimeter;
- `temp_C`: the temperatere measured in degrees Celcius;
- `water_depth_m`: the water depth measured in meters;
- `sludge_thickness_m`: the thickness of the sludge layer formed by decayed _S. aloides_ measured in meters.

&nbsp;

<img src="https://github.com/yulanvanoppen/storage/blob/main/ritsbeheer_GL1.jpg" width="550">

**Figure 2.** The part of ditch 'GL1' where the 'checker' treatment was applied. Photo copyright by [Bureau Biota (Groningen, NL)](https://www.bureaubiota.com "https://www.bureaubiota.com").

&nbsp;

<img src="https://github.com/yulanvanoppen/storage/blob/main/map.png" width="700">

**Figure B.1 in [2].** Simplified map of the Northern Netherlands showing the 17 locations at which the data were collected. The labels are composed of a manager label and a number as in the data.

&nbsp;


## References
[1] Milder-Mulderij, G.,Brochard, C., Wiggers, R., & de Vries, S. (2020). Alternatief krabbenscheerbeheer in Fryslân, Groningen en Drenthe.  Interpretatie op basis van vier jaar onderzoek op diverse locaties. Bureau Biota.

[2] van Oppen, Y. B., Milder-Mulderij, G., Brochard, C., Wiggers, R., de Vries, S., Krijnen, W. P., & Grzegorczyk, M. A. (2022). Modeling dragonfly population data with a Bayesian bivariate geometric mixed-effects model. _Journal of Applied Statistics_, 1-23.
