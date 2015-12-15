#!/usr/bin/env python
"""
Perform a search on the Absynte web server, save the results as full HTML. Use this parser to parse the HTML to obtain a list of genbank ids, then retrieve a fasta file containing those genes, then perform a multiple sequence alignment, then use espript to color according to conservation. In the middle of the process it is necessary to look visually at the Absynte results to see what threshold to use to exclude less well-conserved genes. Also use the representative taxid table to help cull to only have one strain represented from each species. Ultimately it will be good to do a further taxonomic analysis. 
"""

#import sys,re
#from Bio import Entrez,AlignIO,SeqIO
from bs4 import BeautifulSoup

def checkCoords(coords):
    coordinates = [int(x) for x in coords.split(',')]
    return 400>=coordinates[0] and 400<=coordinates[2]

def Absynte2ids(infile,threshold,outfile):
    soup = BeautifulSoup(open(infile),"html.parser")
    maintable=soup.find('table',id="MainContent_Table1")
    rows = maintable.findAll('tr')[1:]
    Results = []
    for row in maintable.findAll('tr')[1:]:
        tds = row.findAll('td')
        organism = tds[0].i.string
        info = tds[0].text
        organism = info.split()[0]
        score = float(info.split()[2])
        #print score
        geneID = None
        if score>threshold:
            for area in tds[1].findAll('area'):
                if checkCoords(area['coords']):
                    if geneID == None:
                        geneID = area['class'][0].split('|')[1]
                        Results.append((organism,geneID,score))
                    else:
                        print "Error: %s has more than one candidate, %s and %s"%(organism,geneID,area['class'][0].split('|')[1])
    openfile = open(outfile,'w')
    for organism,ID,score in Results:
        openfile.write('%s\t%s\t%0.2f\n'%(organism,ID,score))
    openfile.close()

infile = '/Users/edmonds/Documents/GiedrocLab/CzrA/SoxRs/Absynte/Absynte Results.html'
threshold = 28
outfile = '/Users/edmonds/Documents/GiedrocLab/CzrA/SoxRs/Absynte/SoxR_ID_list.txt'

#Absynte2ids('/Users/edmonds/Documents/GiedrocLab/CzrA/SoxRs/Absynte/test.html',28,'/Users/edmonds/Documents/GiedrocLab/CzrA/SoxRs/Absynte/SoxR_ID_list.txt')

Absynte2ids(infile, threshold, outfile)

# Now take this list and look at the actual HTML file and remove lines that don't have the recognizable gene organization
# Here is the result:
## Paracoccus_denitrificans_PD1222_bp13020_C2	119377446	63.96
## Rhodobacter_sphaeroides_ATCC_17025_bp15755_C2	145558164	38.31
## Loktanella_hongkongensis_DSM_17492_bp183668_C3	598658841	38.31
## Pseudovibrio_sp_FO_BEG1_bp73563_C1	359342511	36.04
## Roseobacter_litoralis_Och_149_bp19357_C1	338758684	35.71
## Octadecabacter_temperatus_bp290421_C1	906413296	35.06
## Dinoroseobacter_shibae_DFL_12_DSM_16493_bp17417_C1	157913112	34.74
## Confluentimicrobium_sp_EMB200_NS6_bp274184_C1	933694109	34.09
## Phaeobacter_gallaeciensis_DSM_26640_bp227502_C1	568230024	34.09
## Leisingera_methylohalidivorans_DSM_14336_bp74371_C1	568224622	34.09
## Phaeobacter_inhibens_DSM_17395_bp19355_C1	398656989	34.09
## Phaeobacter_gallaeciensis_2_10_bp19353_C1	398653273	34.09
## Roseibacterium_elongatum_DSM_19469_bp189501	594549561	33.44
## Ruegeria_pomeroyi_DSS_3_bp281_C1	56677627	33.12
## Polymorphum_gilvum_SL003B_26A1_bp63293_C1	326414051	32.01
## Bradyrhizobium_sp_BTAi1_bp16137_C1	146407240	31.53
## Methylobacterium_populi_BJ001_bp19559_C1	179346856	31.14
## Methylobacterium_extorquens_AM1_bp20_C1	240009106	31.01
## Bradyrhizobium_diazoefficiens_USDA_110_bp57599	27378617	30.52
## Rhodopseudomonas_palustris_BisA53_bp15751	115519035	28.90
## Xanthobacter_autotrophicus_Py2_bp15756_C1	154159559	28.64
## Rhodopseudomonas_palustris_BisB5_bp15749	91685177	28.38
## Rhodopseudomonas_palustris_DX_1_bp38503	315602196	28.02
## Hyphomicrobium_denitrificans_ATCC_51888_bp33261	299524051	28.02
## Rhodopseudomonas_palustris_TIE_1_bp20167	192287075	28.02
## Rhodopseudomonas_palustris_HaA2_bp15747	86574508	28.02
