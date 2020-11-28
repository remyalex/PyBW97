**PyBW97 v1.0.0 [Pyhon 2.7]**:

Bakun & Wenworth (1997) Method [(article)](http://usuarios.geofisica.unam.mx/cruz/Sismociones_Libres/Biblio_Sismocion/Bakun_and_Wentworth_BSSA_1997.pdf) for Earthquake Physical parameters estimation Program:
Computes Physical parameters from Earthquake Intensities file in sample output format:

* ** Intensities File Format **
```bash
    Head [Coordinates of Intensities Center and Grid Parameters]:
         <Long> <Lat>;<GridRadious(kms)>;<GridStep(kms)>
    Following lines [Coordinates and intensity value (EMS98 scale)]:
         <Long> <Lat> <Intensity(ems98)>\
  ** End (<filename>.int)**
```
*  --- Required Parameters ---
```bash
         -if Path of intensities file: /<path>/<to>/<filename>.int
         -m IPE Models (Intensity Prediction Equation):
         <1> - Palme et al (2005)       :  I = 2.2 - 1.6*M + 4e-2*(Dist) + 0*(LogDist)
         <2> - Gomez-Capera et al (2016):  I = 1.92 - 2.3*M + 2.1e-3*(Dist) + 3.68*(LogDist)
   --- End ---
```   
*  --- Optional Parameters ---
```bash
       -mn Name of magnitudes file: <magnitudesName>
       -rn Name of rms file: <rmsName>
       -b  Number of Bootstrap resamples
   --- End ---
```
Author: Remy Galan
Licence: GNU GPLv3
