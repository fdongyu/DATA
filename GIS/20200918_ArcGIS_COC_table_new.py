# -*- coding: utf-8 -*-
"""
Created on Tue Apr 28 08:32:00 2020

@author: dfeng
"""

import pandas as pd
from datetime import datetime
import numpy as np
import math
import utm
import osgeo.ogr
from decimal import Decimal

import pdb


class ArcGIS_COC(object):
    
    """
    class to write COC data into ArcGIS shapefile
    """
    
    def __init__(self, **kwargs):
        
        self.readCOCtable()
        
    def readCOCtable(self):
        """
        read COC table, contaminant name, solubility, density, 
        treatability, degradation factor, volatility
        """
        
        #filename = r'\\aus1.aus.apai\share\Projects\0326\099-09\2-0 Wrk Prod\Dongyu_work\spill_modeling\treatability\20200603_COC_updates\tblCOC_06012020.xlsx'
        #filename = r'M:\Projects\0326\099-09\2-0 Wrk Prod\Dongyu_work\spill_modeling\treatability\COC_table.xlsx'
        #filename = r'M:\Projects\0326\099-09\2-0 Wrk Prod\Dongyu_work\spill_modeling\treatability\20200807_COC_updates\tblCOC_06012020_ET_edits_20200806_1953.xlsx'
        #filename = r'M:\Projects\0326\099-09\2-0 Wrk Prod\Dongyu_work\spill_modeling\treatability\20200825_COC_updates\tblCOC_20200825_0815.xlsx'
        filename = r'M:\Projects\0326\099-09\2-0 Wrk Prod\Dongyu_work\spill_modeling\treatability\DRAFT_COC_database_20200915_2056_SENT_TO_NTMWD.xlsx'
        
        xl = pd.ExcelFile(filename)
        df = xl.parse('COC', skiprows=0)
        #pdb.set_trace()
        #df = df.drop(df.index[0])
        
        self.COC = df['COC'].values
        self.CD_APAI = df['CD_APAI'].values.astype('i4').astype('U25')
        self.SYNONYM = df['SYNONYM'].values.astype('U50')   # Alternate name(s) for chemical contaminant 
        self.CHEM_FORM = df['CHEM_FORM'].values.astype('U50') 
        
        self.GR_EPA = df['GR_EPA'].values.astype('U25') 
        self.GR_TNRCC = df['GR_TNRCC'].values.astype('U50') 
        self.GR_USGS = df['GR_USGS'].values.astype('U25') 
        self.GR_TCEQEPA = df['GR_TCEQEPA'].values.astype('U25') 
        
        self.CAS = df['CAS'].values.astype('U25') 
        self.CD_TNRCC = df['CD_TNRCC'].values.astype('i4').astype('U25')
        self.CD_STORET = df['CD_STORET'].values.astype('i4').astype('U25')
        self.CD_TCEQEPA = df['CD_TCEQEPA'].values.astype('U25')
        self.ASSESS = df['ASSESS'].values.astype('U25')
        
        self.PWS_RULE = df['PWS_RULE'].values.astype('U25')
        self.REGULATED = df['REGULATED'].values.astype('U25')
        
        self.TH_VALUE = df['TH_VALUE'].values.astype('U25')
        self.TH_UNITS = df['TH_UNITS'].values.astype('U25')
        self.TH_REASON = df['TH_REASON'].values.astype('U50')
        self.TH_HRF = df['TH_HRF'].values.astype('U25')
        
        self.MCL = df['MCL'].values.astype('U25')
        self.MCL_UNITS = df['MCL_UNITS'].values.astype('U25')
        self.MCL_HRF = df['MCL_HRF'].values.astype('U25')
        
        self.DET_LIM = df['DET_LIM'].values.astype('U25')
        self.DET_UNITS = df['DET_UNITS'].values.astype('U25')
        
        self.WQS = df['WQS'].values
        self.WQS_UNITS = df['WQS_UNITS'].values.astype('U25')
        self.WQS_SOURCE = df['WQS_SOURCE'].values.astype('U25')
        
        self.RBEL = df['RBEL'].values
        self.RBEL_UNITS = df['RBEL_UNITS'].values.astype('U25')
        
        self.SOL_VALUE = df['SOL_VALUE'].values.astype('U25')
        self.SOL_EXT = df['SOL_EXT'].values.astype('U25')
        self.SOL_GROUP = df['SOL_GROUP'].values.astype('i4').astype('U25')
        
        self.DENSITY = df['DENSITY'].values.astype('U25')
        self.SINK_FLAG = df['SINK_FLAG'].values.astype('U25')
        self.KH = df['KH'].values
        self.KOC = df['KOC'].values
        self.LOG_KD = df['LOG_KD'].values
        self.DEG_RATE = df['DEG_RATE'].values
        self.HALFLIFE_D = df['HALFLIFE_D'].values.astype('U25')
        
        self.AIR_DECAY = df['AIR_DECAY'].values
        self.H2O_DECAY = df['H2O_DECAY'].values
        self.LOG_KOW = df['LOG_KOW'].values.astype('U25')
        
        self.TREATMENT = df['TREATMENT'].values.astype('U100')
        self.TR_COMMENT = df['TR_COMMENT'].values.astype('U300')
        self.PCT_REMOV = df['PCT_REMOV'].values.astype('U300')
        self.OTHER_STR = df['OTHER_STR'].values.astype('U300')
        

        #### convert nan string 
        self.SYNONYM[self.SYNONYM=='nan'] = ''
        self.CHEM_FORM[self.CHEM_FORM=='nan'] = ''
        self.GR_EPA[self.GR_EPA=='nan'] = ''
        self.GR_TNRCC[self.GR_TNRCC=='nan'] = ''
        self.GR_USGS[self.GR_USGS=='nan'] = ''
        self.GR_TCEQEPA[self.GR_TCEQEPA=='nan'] = ''
        
        self.CAS[self.CAS=='nan'] = ''
        self.CD_TNRCC[self.CD_TNRCC=='nan'] = ''
        self.CD_STORET[self.CD_STORET=='nan'] = ''
        self.CD_TCEQEPA[self.CD_TCEQEPA=='nan'] = ''
        self.ASSESS[self.ASSESS=='1.0'] = 'True'
        self.ASSESS[self.ASSESS=='0.0'] = 'False'
        
        self.PWS_RULE[self.PWS_RULE=='nan'] = ''
        self.REGULATED[self.REGULATED=='1.0'] = 'True'
        self.REGULATED[self.REGULATED=='0.0'] = 'False'
        
        def conversion(inarray):
            outarray = []
            for i in range(len(inarray)):
                if np.isnan(inarray[i]):
                    outarray.append('')
                elif inarray[i] == 0:
                    outarray.append("{:d}".format(int(inarray[i])))
                else:
                    outarray.append("{:.3E}".format(Decimal(inarray[i]))) 
            
            return outarray
    
        
        
        self.TH_VALUE[self.TH_VALUE=='nan'] = ''
        self.TH_UNITS[self.TH_UNITS=='nan'] = ''
        self.TH_REASON[self.TH_REASON=='nan'] = ''
        self.TH_HRF[self.TH_HRF=='nan'] = ''
        
        self.MCL[self.MCL=='nan'] = ''
        self.MCL_UNITS[self.MCL_UNITS=='nan'] = ''
        self.MCL_HRF[self.MCL_HRF=='nan'] = ''
        
        
        self.DET_LIM[self.DET_LIM=='nan'] = ''
        self.DET_UNITS[self.DET_UNITS=='nan'] = ''
        
        self.WQS = conversion(self.WQS)
        self.WQS_UNITS[self.WQS_UNITS=='nan'] = ''
        self.WQS_SOURCE[self.WQS_SOURCE=='nan'] = ''
        
        self.RBEL = conversion(self.RBEL)
        self.RBEL_UNITS[self.RBEL_UNITS=='nan'] = ''
        
        self.SOL_VALUE[self.SOL_VALUE=='nan'] = ''
        self.SOL_EXT[self.SOL_EXT=='nan'] = ''
        self.SOL_GROUP[(self.SOL_GROUP!='0')&(self.SOL_GROUP!='1')&(self.SOL_GROUP!='2')] = ''
        
        self.DENSITY[self.DENSITY=='nan'] = ''
        self.SINK_FLAG[self.SINK_FLAG=='nan'] = ''
        
        self.KH = conversion(self.KH)
        self.KOC = conversion(self.KOC)
        self.LOG_KD = conversion(self.LOG_KD)
        self.DEG_RATE = conversion(self.DEG_RATE)
        self.HALFLIFE_D[self.HALFLIFE_D=='nan'] = ''
        
        self.AIR_DECAY = conversion(self.AIR_DECAY)
        self.H2O_DECAY = conversion(self.H2O_DECAY)
        
        self.LOG_KOW[self.LOG_KOW=='nan'] = ''
        self.TREATMENT[self.TREATMENT=='nan'] = ''
        self.TR_COMMENT[self.TR_COMMENT=='nan'] = ''
        self.PCT_REMOV[self.PCT_REMOV=='nan'] = ''
        self.OTHER_STR[self.OTHER_STR=='nan'] = ''
        
        #### convert nan string 
        self.TYPE = []
        for i in range(len(self.SOL_EXT)):
            if self.SOL_EXT[i] == 'Soluble':
                self.TYPE.append('Soluble')
                
            elif self.SOL_EXT[i]=='Insoluble' or self.SOL_EXT[i] =='Slightly Soluble':
                
                if self.SINK_FLAG[i]=='Float' or self.SINK_FLAG[i]=='Neutral' or self.SINK_FLAG[i]=='NoData':
                    self.TYPE.append('Light Insoluble')
                    
                elif self.SINK_FLAG[i]=='Sink':
                    self.TYPE.append('Heavy Insoluble')
                    
                elif self.SINK_FLAG[i] == '':
                    self.TYPE.append('')
                    
            elif self.SOL_EXT[i] == '':
                self.TYPE.append('')
        
        #pdb.set_trace()
        if len(self.TYPE) != len(self.COC):
            raise IOError('Check the contaminant type variable!!!')
                
        self.FULL_NAME = [] ## treatment full name
        #AA = activated alumina; C/F = coagulation/filtration (not BAT for systems < 500 service connections); 
        #DDF = direct and diatomite filtration; GAC = granular activated carbon; IE = ion exchange; 
        #LS = lime softening (not BAT for systems < 500 service connections); RO = reverse osmosis; 
        #CC = corrosion control; ED = electrodialysis; Cl2 = chlorine; UV = ultraviolet; 
        #O/F = oxidation/filtration; AC = alkaline chlorination (pH ≥ 8.5); PTA = packed tower aeration; 
        #O3 = ozone; PAC = packed activated carbon; UV+H2O2 = ultraviolet + hydrogen peroxide; 
        #O3+H2O2 = ozone + hydrogen peroxide; H2O2 = hydrogen peroxide; AS = air stripping, 
        #AC = activated carbon
        
        #AA = 'AA = activated alumina; '; 
        #C_F = 'C/F = coagulation/filtration (not BAT for systems < 500 service connections); '; 
        #DDF = 'DDF = direct and diatomite filtration; '; GAC = 'GAC = granular activated carbon; '; IE = 'IE = ion exchange; '; 
        #LS = "LS = lime softening (not BAT for systems < 500 service connections); "; RO = 'RO = reverse osmosis; '; 
        #CC = 'CC = corrosion control; '; ED = 'ED = electrodialysis; '; Cl2 = 'Cl2 = chlorine; '; UV = 'UV = ultraviolet; '; 
        #O_F = 'O/F = oxidation/filtration; '; AC = 'AC = alkaline chlorination (pH ≥ 8.5); '; PTA = 'PTA = packed tower aeration; '; 
        #O3 = 'O3 = ozone; '; PAC = 'PAC = packed activated carbon; '; UV_H2O2 = 'UV+H2O2 = ultraviolet + hydrogen peroxide; '; 
        #O3_H2O2 = 'O3+H2O2 = ozone + hydrogen peroxide; '; H2O2 = 'H2O2 = hydrogen peroxide; '; AS = 'AS = air stripping; '; 
        MS = 'MS = Membrane Separation; '; GAC = 'GAC = Granular Activated Carbon; '; LS = 'LS = Lime Or Precipitative Softening; ';
        CHE = 'CHE = Chemical Oxidation; '; MF = 'MF = Micro- and Ultra-membrane filtration; '; AIR_AS = 'AIR+AS = Aeration And Air Stripping; ';
        ED = 'ED = Electrodialysis; '; PAC = 'PAC = Powdered Activated Carbon; '; CONV = 'CONV = Coagulation/Filtration; ';
        AA = 'AA = Activated Alumina; '; IE = 'IE = Ion Exchange; '; O3 = 'O3 = Ozonation; ';
        Ads = 'Ads = Adsorptive Media; '; CL2 = 'CL2 = Chlorination (with Chloramine Or Chlorine); '; UV_H2O2 = 'UV+H2O2 = Ultraviolet Irradiation + Hydrogen Peroxide; ';
        UV = 'UV = Ultraviolet Radiation; '; AOP = 'AOP = Advanced Oxidative Processes; '; ClO2 = 'ClO2 = Chlorine Dioxide; ';
        H2O2 = 'H2O2 = Hydrogen Peroxide; '; UV_O3 = 'UV+O3 = Ultraviolet Irradiation + Ozone; '; Bio = 'Bio = Biological Treatment; ';
        MnO4 = 'MnO4 = Permanganate; '; BAC = 'BAC = Biologically Active Filters; '; EC = 'EC = Electrocoagulation; '; 
        MCA = 'MCA = Monochloramination; '; Acid = 'Acid = Acidification; '; Clar = 'Clar = Clarification; '; 
        CP = 'CP = Chemical Precipitation; '; DAF = 'DAF = Dissolved Air Flotation; '; Dist = 'Dist = Distillation; ';
        EF = 'EF = Electroflotation; '; MBR = 'MBR = Membrane Bioreactor; '; O_F = 'O/F = Oxidation/filtration; ';
        PTA = 'PTA = Packed Tower Aeration; '; AcSl = 'AcSl = Activated Sludge; '; ANFF = 'ANFF = Anaerobic Fixed Film Biological Treatment; ';
        BNR = 'BNR = Biological Nutrient Removal; '; CC = 'CC = Corrosion Control; '; DDF = 'DDF = Direct And Diatomite Filtration; ';  
        ECO = 'ECO = Electrochemical Oxidation; '; EQ = 'EQ = Flow Equalization; '; Evap = 'Evap = Evaporation; ';
        FO = 'FO = Foward Osmosis; '; FR = 'FR = Fenton Reaction; '; GMF = 'GMF = Granular Media Filtration; '; 
        Hydro = 'Hydro = Hydrolysis; '; LGE = 'LGE = Liquified Gas Extraction; '; MBBR = 'MBBR = Moving Bed Bioreactor; '; 
        MD = 'MD = Membrane Distillation; '; Photo = 'Photo = Photodegradation; '; ZVI = 'ZVI = Zero Valent Iron; '
        
        
        #AIR = 'AIR = Aeration; '; 
        #AS = 'AS = Air Stripping; '; 
        #C_F = 'C/F = Coagulation/Filtration; ';
        #RO = 'RO = reverse osmosis; '; 
        
        for i in range(len(self.TREATMENT)):
            tem = ''
            
            if 'MS' in self.TREATMENT[i]:
                tem+=MS
            if 'GAC' in self.TREATMENT[i]:
                tem+=GAC
            if 'LS' in self.TREATMENT[i]:
                tem+=LS
            if 'CHE' in self.TREATMENT[i]:
                tem+=CHE
            if 'MF' in self.TREATMENT[i]:
                tem+=MF
            if 'AIR+AS' in self.TREATMENT[i]:
                tem+=AIR_AS
            if 'ED' in self.TREATMENT[i]:
                tem+=ED
            if 'PAC' in self.TREATMENT[i]:
                tem+=PAC
            if 'CONV' in self.TREATMENT[i]:
                tem+=CONV
            if 'AA' in self.TREATMENT[i]:
                tem+=AA
            if 'IE' in self.TREATMENT[i]:
                tem+=IE
            if 'O3' in self.TREATMENT[i] and 'UV+O3' not in self.TREATMENT[i]:
                tem+=O3
            if 'Ads' in self.TREATMENT[i]:
                tem+=Ads
            if 'CL2' in self.TREATMENT[i]:
                tem+=CL2
            if 'UV+H2O2' in self.TREATMENT[i]:
                tem+=UV_H2O2
            if 'UV' in self.TREATMENT[i] and 'UV+H2O2' not in self.TREATMENT[i] and 'UV+O3' not in self.TREATMENT[i]:
                tem+=UV
            if 'AOP' in self.TREATMENT[i]:
                tem+=AOP
            if 'ClO2' in self.TREATMENT[i]:
                tem+=ClO2
            if 'H2O2' in self.TREATMENT[i] and 'UV+H2O2' not in self.TREATMENT[i]:
                tem+=H2O2    
            if 'UV+O3' in self.TREATMENT[i]:
                tem+=UV_O3
            if 'Bio' in self.TREATMENT[i]:
                tem+=Bio
            if 'MnO4' in self.TREATMENT[i]:
                tem+=MnO4
            if 'BAC' in self.TREATMENT[i]:
                tem+=BAC
            if 'EC' in self.TREATMENT[i] and 'ECO' not in self.TREATMENT[i]:
                tem+=EC
            if 'MCA' in self.TREATMENT[i]:
                tem+=MCA
            if 'Acid' in self.TREATMENT[i]: 
                tem+=Acid
            if 'Clar' in self.TREATMENT[i]:
                tem+=Clar
            if 'CP' in self.TREATMENT[i]:
                tem+=CP
            if 'DAF' in self.TREATMENT[i]:
                tem+=DAF
            if 'Dist' in self.TREATMENT[i]:
                tem+=Dist
            if 'EF' in self.TREATMENT[i]:
                tem+=EF
            if 'MBR' in self.TREATMENT[i]:
                tem+=MBR
            if 'O/F' in self.TREATMENT[i]:
                tem+=O_F
            if 'PTA' in self.TREATMENT[i]:
                tem+=PTA
            if 'AcSl' in self.TREATMENT[i]:
                tem+=AcSl
            if 'ANFF' in self.TREATMENT[i]:
                tem+=ANFF
            if 'BNR' in self.TREATMENT[i]:
                tem+=BNR
            if 'CC' in self.TREATMENT[i]:
                tem+=CC 
            if 'DDF' in self.TREATMENT[i]:
                tem+=DDF  
            if 'ECO' in self.TREATMENT[i]:
                tem+=ECO 
            if 'EQ' in self.TREATMENT[i]:
                tem+=EQ
            if 'Evap' in self.TREATMENT[i]:
                tem+=Evap
            if 'FO' in self.TREATMENT[i]:
                tem+=FO    
            if 'FR' in self.TREATMENT[i]:
                tem+=FR  
            if 'GMF' in self.TREATMENT[i]:
                tem+=GMF
            if 'Hydro' in self.TREATMENT[i]:
                tem+=Hydro
            if 'LGE' in self.TREATMENT[i]:
                tem+=LGE
            if 'MBBR' in self.TREATMENT[i]:
                tem+=MBBR
            if 'MD' in self.TREATMENT[i]:
                tem+=MD   
            if 'Photo' in self.TREATMENT[i]:
                tem+=Photo
            if 'ZVI' in self.TREATMENT[i]:
                tem+=ZVI
                
                         
#            if 'AIR' in self.TREATMENT[i]:
#                tem+=AIR
#            if 'AS' in self.TREATMENT[i]:
#                tem+=AS        
#            if 'C/F' in self.TREATMENT[i]:
#                tem+=C_F                                                 
#            if 'RO' in self.TREATMENT[i]:
#                tem+=RO

            self.FULL_NAME.append(tem)
       
        
        
        self.NOTES = [] ## additional notes
        for i in range(len(self.TYPE)):
            if self.TYPE[i] == 'Heavy Insoluble':
                self.NOTES.append('Model shows generally very little current movement at the bottom of the lake. As such heavy insoluble contaminants that sink to the lake bottom are likely to spread very slowly, the model did not predict any usable time-of-travel for this contaminant. ')
            else:
                self.NOTES.append('')
        
        
#        #### QA/QC
#        for i in range(len(self.TREATMENT)):
#            print ('Number: ', str(i+1))
#            print ('Name: ', self.COC[i])
#            print ('Type of Contaminant: ', self.TYPE[i])
#            print ('Half life: ', self.HALFLIFE_D[i])
#            print ('Volatility: ', self.KH[i])
#            print ('Treatment Technology: ', self.TREATMENT[i])
#            print ('Abbreviation: ', self.FULL_NAME[i])
#            print ('Comments: ', self.TR_COMMENT[i])
#            print ('Other_STR: ', self.OTHER_STR[i])
#            print ('Additional notes: ', self.NOTES[i])
#            print ('------------------------------------------------------------------------')
            
            
        
        pdb.set_trace()
        
    
    def writeShp(self,shpname):
        """
        COC name
        Soluble in Water/Heavier than Water/Lighter than Water
        Effective and available treatment technologies
        degradation factor
        volatility
        """
        
        #### Create the shapefile
        # Create the projection
        spatialReference = osgeo.osr.SpatialReference()
        spatialReference.ImportFromProj4('+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs')
        
        # Create the shape file
        outfile = r'%s'%shpname
        driver = osgeo.ogr.GetDriverByName('ESRI Shapefile')
        shapeData = driver.CreateDataSource(outfile)
        
        # Create the layer
        layer = shapeData.CreateLayer('Contour', spatialReference, osgeo.ogr.wkbPoint)
        layerDefinition = layer.GetLayerDefn()
        
        
        # Create fields containing 
        field_def = osgeo.ogr.FieldDefn('COC', osgeo.ogr.OFTString)
        layer.CreateField(field_def)
        
        field_def = osgeo.ogr.FieldDefn('CD_APAI', osgeo.ogr.OFTString)
        layer.CreateField(field_def)
        
        field_def = osgeo.ogr.FieldDefn('SYNONYM', osgeo.ogr.OFTString)
        layer.CreateField(field_def)
        
        field_def = osgeo.ogr.FieldDefn('CHEM_FORM', osgeo.ogr.OFTString)
        layer.CreateField(field_def)
        
        field_def = osgeo.ogr.FieldDefn('GR_EPA', osgeo.ogr.OFTString)
        layer.CreateField(field_def)
        
        field_def = osgeo.ogr.FieldDefn('GR_TNRCC', osgeo.ogr.OFTString)
        layer.CreateField(field_def)
        
        field_def = osgeo.ogr.FieldDefn('GR_USGS', osgeo.ogr.OFTString)
        layer.CreateField(field_def)
        
        field_def = osgeo.ogr.FieldDefn('GR_TCEQEPA', osgeo.ogr.OFTString)
        layer.CreateField(field_def)
        
        field_def = osgeo.ogr.FieldDefn('CAS', osgeo.ogr.OFTString)
        layer.CreateField(field_def)
        
        field_def = osgeo.ogr.FieldDefn('CD_TNRCC', osgeo.ogr.OFTString)
        layer.CreateField(field_def)
        
        field_def = osgeo.ogr.FieldDefn('CD_STORET', osgeo.ogr.OFTString)
        layer.CreateField(field_def)
        
        field_def = osgeo.ogr.FieldDefn('CD_TCEQEPA', osgeo.ogr.OFTString)
        layer.CreateField(field_def)
        
        field_def = osgeo.ogr.FieldDefn('ASSESS', osgeo.ogr.OFTString) 
        layer.CreateField(field_def)
        
        field_def = osgeo.ogr.FieldDefn('PWS_RULE', osgeo.ogr.OFTString) 
        layer.CreateField(field_def)
        
        field_def = osgeo.ogr.FieldDefn('REGULATED', osgeo.ogr.OFTString) 
        layer.CreateField(field_def)
        
        field_def = osgeo.ogr.FieldDefn('TH_VALUE', osgeo.ogr.OFTString) 
        layer.CreateField(field_def)
        
        field_def = osgeo.ogr.FieldDefn('TH_UNITS', osgeo.ogr.OFTString) 
        layer.CreateField(field_def)
        
        field_def = osgeo.ogr.FieldDefn('TH_REASON', osgeo.ogr.OFTString) 
        layer.CreateField(field_def)
        
        field_def = osgeo.ogr.FieldDefn('TH_HRF', osgeo.ogr.OFTString) 
        layer.CreateField(field_def)
        
        field_def = osgeo.ogr.FieldDefn('MCL', osgeo.ogr.OFTString) 
        layer.CreateField(field_def)
        
        field_def = osgeo.ogr.FieldDefn('MCL_UNITS', osgeo.ogr.OFTString) 
        layer.CreateField(field_def)
        
        field_def = osgeo.ogr.FieldDefn('MCL_HRF', osgeo.ogr.OFTString) 
        layer.CreateField(field_def)
        
        field_def = osgeo.ogr.FieldDefn('DET_LIM', osgeo.ogr.OFTString) 
        layer.CreateField(field_def)
        
        field_def = osgeo.ogr.FieldDefn('DET_UNITS', osgeo.ogr.OFTString) 
        layer.CreateField(field_def)
        
        field_def = osgeo.ogr.FieldDefn('WQS', osgeo.ogr.OFTString) 
        layer.CreateField(field_def)
        
        field_def = osgeo.ogr.FieldDefn('WQS_UNITS', osgeo.ogr.OFTString) 
        layer.CreateField(field_def)
        
        field_def = osgeo.ogr.FieldDefn('WQS_SOURCE', osgeo.ogr.OFTString) 
        layer.CreateField(field_def)
        
        field_def = osgeo.ogr.FieldDefn('RBEL', osgeo.ogr.OFTString) 
        layer.CreateField(field_def)
        
        field_def = osgeo.ogr.FieldDefn('RBEL_UNITS', osgeo.ogr.OFTString) 
        layer.CreateField(field_def)
        
        field_def = osgeo.ogr.FieldDefn('SOL_VALUE', osgeo.ogr.OFTString) 
        layer.CreateField(field_def)
        
        field_def = osgeo.ogr.FieldDefn('SOL_EXT', osgeo.ogr.OFTString) 
        layer.CreateField(field_def)
        
        field_def = osgeo.ogr.FieldDefn('SOL_GROUP', osgeo.ogr.OFTString) 
        layer.CreateField(field_def)
        
        field_def = osgeo.ogr.FieldDefn('DENSITY', osgeo.ogr.OFTString) 
        layer.CreateField(field_def)
        
        field_def = osgeo.ogr.FieldDefn('SINK_FLAG', osgeo.ogr.OFTString) 
        layer.CreateField(field_def)
        
        field_def = osgeo.ogr.FieldDefn('KH', osgeo.ogr.OFTString) 
        layer.CreateField(field_def)
        
        field_def = osgeo.ogr.FieldDefn('KOC', osgeo.ogr.OFTString) 
        layer.CreateField(field_def)
        
        field_def = osgeo.ogr.FieldDefn('LOG_KD', osgeo.ogr.OFTString) 
        layer.CreateField(field_def)
        
        field_def = osgeo.ogr.FieldDefn('DEG_RATE', osgeo.ogr.OFTString) 
        layer.CreateField(field_def)
        
        field_def = osgeo.ogr.FieldDefn('HALFLIFE_D', osgeo.ogr.OFTString) 
        layer.CreateField(field_def)
        
        field_def = osgeo.ogr.FieldDefn('AIR_DECAY', osgeo.ogr.OFTString) 
        layer.CreateField(field_def)
        
        field_def = osgeo.ogr.FieldDefn('H2O_DECAY', osgeo.ogr.OFTString) 
        layer.CreateField(field_def)
        
        field_def = osgeo.ogr.FieldDefn('LOG_KOW', osgeo.ogr.OFTString) 
        layer.CreateField(field_def)
        
        field_def = osgeo.ogr.FieldDefn('TREATMENT', osgeo.ogr.OFTString) 
        layer.CreateField(field_def)
        
        field_def = osgeo.ogr.FieldDefn('TR_COMMENT', osgeo.ogr.OFTString) 
        layer.CreateField(field_def)
        
        field_def = osgeo.ogr.FieldDefn('PCT_REMOV', osgeo.ogr.OFTString) 
        layer.CreateField(field_def)
        
        field_def = osgeo.ogr.FieldDefn('OTHER_STR', osgeo.ogr.OFTString)
        layer.CreateField(field_def)
        
        field_def = osgeo.ogr.FieldDefn('TYPE', osgeo.ogr.OFTString)
        layer.CreateField(field_def)
        
        field_def = osgeo.ogr.FieldDefn('FULL_NAME', osgeo.ogr.OFTString) 
        layer.CreateField(field_def)
        
        field_def = osgeo.ogr.FieldDefn('NOTES', osgeo.ogr.OFTString) 
        layer.CreateField(field_def)
        
        
        def add_feature(layer, COC, CD_APAI, SYNONYM, CHEM_FORM, 
                        GR_EPA, GR_TNRCC, GR_USGS, GR_TCEQEPA, CAS, 
                        CD_TNRCC, CD_STORET, CD_TCEQEPA, ASSESS, 
                        PWS_RULE, REGULATED, 
                        TH_VALUE, TH_UNITS, TH_REASON, TH_HRF, 
                        MCL, MCL_UNITS, MCL_HRF, DET_LIM, DET_UNITS, 
                        WQS, WQS_UNITS, WQS_SOURCE, RBEL, RBEL_UNITS, 
                        SOL_VALUE, SOL_EXT, SOL_GROUP, DENSITY, SINK_FLAG, 
                        KH, KOC, LOG_KD, DEG_RATE, HALFLIFE_D, 
                        AIR_DECAY, H2O_DECAY, LOG_KOW, TREATMENT, TR_COMMENT, PCT_REMOV, OTHER_STR,
                        TYPE, FULL_NAME, NOTES):
            """
            function that adds feature to layer
            """    
            ctr=0
            for i in range(len(COC)):
                ctr+=1
#                point = osgeo.ogr.Geometry(osgeo.ogr.wkbPoint)
                # Add points individually to the line
                #point.AddPoint(lon[i], lat[i])
                # Update the feature with the line data
                featureIndex = ctr
                feature = osgeo.ogr.Feature(layerDefinition)
#                feature.SetGeometry(point)
                feature.SetFID(featureIndex)
#                feature.SetGeometryDirectly(point)
                
                # Set the attribute table
                feature.SetField('COC', COC[i])
                feature.SetField('CD_APAI', CD_APAI[i])
                feature.SetField('SYNONYM', SYNONYM[i])
                feature.SetField('CHEM_FORM', CHEM_FORM[i])
                feature.SetField('GR_EPA', GR_EPA[i])
                
                feature.SetField('GR_TNRCC', GR_TNRCC[i])
                feature.SetField('GR_USGS', GR_USGS[i])
                feature.SetField('GR_TCEQEPA', GR_TCEQEPA[i])
                feature.SetField('CAS', CAS[i])
                feature.SetField('CD_TNRCC', CD_TNRCC[i])
                feature.SetField('CD_STORET', CD_STORET[i])
                
                feature.SetField('CD_TCEQEPA', CD_TCEQEPA[i])
                feature.SetField('ASSESS', ASSESS[i])
                feature.SetField('PWS_RULE', PWS_RULE[i])
                feature.SetField('REGULATED', REGULATED[i])
                
                feature.SetField('TH_VALUE', TH_VALUE[i])
                feature.SetField('TH_UNITS', TH_UNITS[i])
                feature.SetField('TH_REASON', TH_REASON[i])
                feature.SetField('TH_HRF', TH_HRF[i])
                
                feature.SetField('MCL', MCL[i])
                feature.SetField('MCL_UNITS', MCL_UNITS[i])
                feature.SetField('MCL_HRF', MCL_HRF[i])
                
                feature.SetField('DET_LIM', DET_LIM[i])
                feature.SetField('DET_UNITS', DET_UNITS[i])
                feature.SetField('WQS', WQS[i])
                feature.SetField('WQS_UNITS', WQS_UNITS[i])
                feature.SetField('WQS_SOURCE', WQS_SOURCE[i])
                
                feature.SetField('RBEL', RBEL[i])
                feature.SetField('RBEL_UNITS', RBEL_UNITS[i])
                
                feature.SetField('SOL_VALUE', SOL_VALUE[i])
                feature.SetField('SOL_EXT', SOL_EXT[i])
                feature.SetField('SOL_GROUP', SOL_GROUP[i])
                
                feature.SetField('DENSITY', DENSITY[i])
                feature.SetField('SINK_FLAG', SINK_FLAG[i])
                
                feature.SetField('KH', KH[i])
                feature.SetField('KOC', KOC[i])
                feature.SetField('LOG_KD', LOG_KD[i])
                feature.SetField('DEG_RATE', DEG_RATE[i])
                feature.SetField('HALFLIFE_D', HALFLIFE_D[i])
                
                feature.SetField('AIR_DECAY', AIR_DECAY[i])
                feature.SetField('H2O_DECAY', H2O_DECAY[i])
                feature.SetField('LOG_KOW', LOG_KOW[i])
                
                feature.SetField('TREATMENT', TREATMENT[i])
                feature.SetField('TR_COMMENT', TR_COMMENT[i])
                feature.SetField('PCT_REMOV', PCT_REMOV[i])
                feature.SetField('OTHER_STR', OTHER_STR[i]) 
                
                feature.SetField('TYPE', TYPE[i])
                feature.SetField('FULL_NAME', FULL_NAME[i])
                feature.SetField('NOTES', NOTES[i])
                
                layer.CreateFeature(feature)
        

        #pdb.set_trace()
        add_feature(layer, self.COC, self.CD_APAI, self.SYNONYM, self.CHEM_FORM, \
        self.GR_EPA, self.GR_TNRCC, self.GR_USGS, self.GR_TCEQEPA, self.CAS, \
        self.CD_TNRCC, self.CD_STORET, self.CD_TCEQEPA, self.ASSESS, \
        self.PWS_RULE, self.REGULATED, \
        self.TH_VALUE, self.TH_UNITS, self.TH_REASON, self.TH_HRF, \
        self.MCL, self.MCL_UNITS, self.MCL_HRF, self.DET_LIM, self.DET_UNITS, \
        self.WQS, self.WQS_UNITS, self.WQS_SOURCE, self.RBEL, self.RBEL_UNITS, \
        self.SOL_VALUE, self.SOL_EXT, self.SOL_GROUP, self.DENSITY, self.SINK_FLAG, \
        self.KH, self.KOC, self.LOG_KD, self.DEG_RATE, self.HALFLIFE_D, \
        self.AIR_DECAY, self.H2O_DECAY, self.LOG_KOW, self.TREATMENT, self.TR_COMMENT, self.PCT_REMOV, self.OTHER_STR, \
        self.TYPE, self.FULL_NAME, self.NOTES)
        
        
#        add_feature(layer, self.CONTAMINANT_CD, self.contaminant_code, self.CAS, self.ASSESS, \
#        self.contaminant_names, self.analyte_group, self.analyte_group2, self.regulated, \
#        self.treatabilities, self.water_solubility, self.classes, self.density, \
#        self.volatilities, self.koc_string, self.degradations, \
#        self.threshold_value, self.threshold_units, self.threshold_reason, \
#        self.air_decay_string, self.H2O_decay_string, self.LOG_Kow_string, \
#        self.details)
        
        """
        long string names, equal or less than 10 letters
        CONTAMINANT
        GROUP_TNRCC
        GROUP_TCEQ_EPA
        CD_TCEQ_EPA
        THRESH_VALUE
        THRESH_UNITS
        THRESH_REASON
        """
        
        
        
        
    
        
if __name__ == "__main__": 
    
    ArcCOC = ArcGIS_COC()
    ArcCOC.writeShp(shpname='COC_table.shp')
        