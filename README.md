# LABS
Bioturbation model 

Bioturbation code from Choi et al. (2002) Conmuters &amp; Geosciences 28, 213-222
https://www.iamg.org/index.php/publisher/articleview/frmArticleID/115

Extended to include calculation of water flow and oxygen and organic matter concentration fields 
Modules of mod_o2_dif+adv.f90 and mod_NS_MAC_2D.f90 do the calculations of oxygen concentration and water flow fields, respectively. 
Module of pfit.f90 does the polynomial fit to the tracer activity profiles. 
