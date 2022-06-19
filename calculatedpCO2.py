def calculatedpco2(pH,temp,alkalinity,DIC,IS):
    from numpy import log10
    #output is log10(pCO2 in atm).
    #temp is the water temperature in Kelvin
    #alkalinity in mol/l
    #DIC in mol/l
    #IS in mol/l
    H=10.0**(-1.0*pH)
    #H2O(l) <--> H+ + OH-
    DHw_25=55.90661 # kJ/mol
    log10Kw_25=-14.0
    # 
    # #CO2(g) <--> CO2(aq)
    DHpco2_25=-19.983 # kJ/mol
    log10Kpco2_25=-1.468
    #CO2(aq) + H2O <--> HCO3- + H+
    DHa1_25=9.109 # kJ/mol
    log10Ka1_25=-6.352
    #HCO3- <--> CO32- + H+
    DHa2_25=14.90 # kJ/mol
    log10Ka2_25=-10.329
    R=0.00831439 # kJ/Kmol
    Kw=10.0**(((DHw_25/(2.303*R))*((1.0/temp)-(1.0/298.15)))+log10Kw_25)
    Kpco2=10.0**(((DHpco2_25/(2.303*R))*((1.0/temp)-(1.0/298.15)))+log10Kpco2_25)
    Ka1=10.0**(((DHa1_25/(2.303*R))*((1.0/temp)-(1.0/298.15)))+log10Ka1_25)
    Ka2=10.0**(((DHa2_25/(2.303*R))*((1.0/temp)-(1.0/298.15)))+log10Ka2_25)
    E=87.74-(0.4008*(temp-273.15))+(0.0009398*((temp-273.15)**2.0))-(0.00000141*((temp-273.15)**3.0))
    Ai=1824830*((E*temp)**(-1.5))
    alphaz1=10.0**(((-1.0*Ai*(IS**(1.0/2.0)))/(1.0+(IS**(1.0/2.0))))+(1.0*0.3*Ai*IS))
    alphaz2=10.0**(((-4.0*Ai*(IS**(1.0/2.0)))/(1.0+(IS**(1.0/2.0))))+(1.0*0.3*4.0*Ai*IS))
    #alphaz1=1;alphaz2=1
    H=H/alphaz1
    #Activity coefficient for CO2(aq) from Drummond et al. (1981)
    g1=1.0312
    g2=0.0012806
    g3=255.9
    g4=0.4445
    g5=0.001606
    alphaz0=2.71828**(((g1+(g2*temp)+(g3/temp))*IS)-((IS*(g4+(g5*temp)))/(IS+1.0)))
    #alphaz0=1
    #CO2 From alkalintiy
    pCO2a=(alkalinity-H+(Kw/(H*alphaz1*alphaz1)))/(((Ka1*Kpco2)/(H*alphaz1*alphaz1))+((2.0*Ka2*Ka1*Kpco2)/((H**2.0)*(alphaz1**2.0)*alphaz2)))
    pCO2a=log10(pCO2a)
    #CO2 from DIC
    pCO2c=DIC/(((Kpco2)/(alphaz0))+((Ka1*Kpco2)/(H*alphaz1*alphaz1))+((Ka2*Ka1*Kpco2)/((H**2.0)*(alphaz1**2.0)*alphaz2)))
    pCO2c=log10(pCO2c)
    return pCO2a,pCO2c # pCO2a (calculated from alkalinity); pCO2c (calculated from DIC)




