DNA1 = "TTGTCGAAAAATATCAGGTGTAAGTTGATAGGTCGATCCATTTATCTGTATATTATATAGAAATGTATGTATGATCTACCGTACCCATTTTATTCAAAAATTCGCCTCGGCTGCATCTTTGAATAAAGTGAGTGGGACTAAGTCATATTATATAGAATCGATGGATTCCGTAGTCAAAANATCCCTCCGTGATACATTACAACGACTGATAGAGGGATCAAATGGTATAGTTTGTTTGTTGGTAGCTTGGAGGATTATAAGTATGACTATAGCTTTCCAATTGGCTGTTTTTGCATTAATTGCTACTTCATCA"
name1 = "LITHOPS GEYERI"

DNA2 = "TTGTCGAAAAATATCAGGTGTAAGTTGATAGGTCGATCCATTTATCTGTATATTATATAGAAATGTATGTATGATCTACCGTACCCATTTTATTCAAAAATTCGCCTCGGCTGCATCTTTGAATAAAGTGAGTGGGACTAAGTCATATTATATAGAATCGATGGATTCCGTAGTCAAAANANCCCTCCGTGATACATTACAACGACTGATAGAGGGATCAAATGGTATAGTTTGTTTGTTGGTAGCTTGGAGGAATATAAGTATGACTATAGCTTTCCAATTGGCTGTTTTTGCATTAATTGCTACTTCA"
name2 = "LITHOPS HOOKERI VAR. HOOKERI"

DNA3 = "TTGTCGAAAAATATCAGGTGTAAGTTGATAGGTCGATCCATTTATCTGTATATTATATAGAAATGTATGTATGATCTACCGTACCCATTTTATTCAAAAATTCGCCTCGGCTGCATCTTTGAATAAAGTGAGTGGGACTAAGTCATATTATATAGAATCGATGGATTCCGGAGTCAAAANATCCCTCCGTGATACATTACAACGACTGATAGAGGGATCAAATGGTAKAGTTTGTTTGTTGGTAGCTTGGAGGATTATAAGTATGACTATAGCTTTCCAATTGGCTGTTTTTGCATTAATTGCGACTTCATC"
name3 = "LITHOPS DOROTHEAE"

aa = ['F','L','I','M','V','S','P','T','A','Y',
      '|','H','Q','N','K','D','E','C','W','R',
      'G']

codons = [['TTT', 'TTC'],
          ['TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG'],
          ['ATT', 'ATC', 'ATA'],
          ['ATG'],
          ['GTT', 'GTC', 'GTA', 'GTG'],
          ['TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC'],
          ['CCT', 'CCC', 'CCA', 'CCG'],
          ['ACT', 'ACC', 'ACA', 'ACG'],
          ['GCT', 'GCC', 'GCA', 'GCG'],
          ['TAT', 'TAC'],
          ['TAA', 'TAG', 'TGA'],
          ['CAT', 'CAC'],
          ['CAA', 'CAG'],
          ['AAT', 'AAC'],
          ['AAA', 'AAG'],
          ['GAT', 'GAC'],
          ['GAA', 'GAG'],
          ['TGT', 'TGC'],
          ['TGG'],
          ['CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'],
          ['GGT', 'GGC', 'GGA', 'GGG']]


#returns the amino acid that corresponds to the codon input
#used in codingStrandToAA(DNA)
def amino(codon):
    for i in range (0, len(codons)):
        for triple in codons[i]:
            if(codon == triple):
                return aa[i]

#returns the protein encoded by input DNA in the form of an array of amino acids
#invokes the amino(codon) method
def codingStrandToAA(DNA):
    finalSeq = ""
    codon = ""
    for i in range(0, len(DNA), 3):
        codon = DNA[i:i+3]
        if(amino(codon) != None):
            finalSeq += amino(codon)
    return finalSeq

#checks the condition that needs to be met in order for an addition mutation to have occured

#Checks to see if there is a difference in amino acids at one index in the protein sequence. If #there is, it checks if a specific amino acid would be better lined up if it were one position #forward or one position back. It changes all sequences accordingly so that amino acids line #up. 
def findMutations(): 
    global AA1
    global AA2
    global AA3
    for i in range(0,90):    
        if(((AA1[i] != AA2[i]) and (AA1[i] != AA3[i]))):
            if((AA1[i] == AA2[i+1]) and (AA1[i] == AA3[i+1])):
                AA1 = (AA1[:i]) + " " + (AA1[i:])
            elif((AA1[i+1] == AA2[i]) and (AA1[i+1] == AA3[i])):
                AA3 = (AA3[:i]) + " " + (AA3[i:])
                AA2 = (AA2[:i]) + " " + (AA2[i:])
        elif(((AA2[i] != AA1[i]) and (AA2[i] != AA3[i]))):
            if((AA2[i] == AA1[i+1]) and (AA2[i] == AA3[i+1])):
                AA2 = (AA2[:i]) + " " + (AA2[i:])
            elif((AA2[i+1] == AA1[i]) and (AA2[i+1] == AA3[i])):
                AA1 = (AA1[:i]) + " " + (AA1[i:])
                AA3 = (AA3[:i]) + " " + (AA3[i:])
        elif(((AA3[i] != AA1[i]) and (AA3[i] != AA2[i]))): 
            if((AA3[i] == AA2[i+1]) and (AA3[i] == AA1[i+1])):
                AA3 = (AA3[:i]) + " " + (AA3[i:])
            elif((AA3[i+1] == AA1[i]) and (AA3[i+1] == AA2[i])):
                AA1 = (AA1[:i]) + " " + (AA1[i:])
                AA2 = (AA2[:i]) + " " + (AA2[i:])

#Assigns each amino acid to one of 20 colors, makes a list of colors representing the protein
def aminoAcidToColor(aminoList):
    colorList = []
    for i in range(0,len(aminoList)-1):
        if(aminoList[i] == "A"):
            colorList.append("#800000")
        if(aminoList[i] == "R"):
            colorList.append("#42d4f4")
        if(aminoList[i] == "N"):
            colorList.append("#808000")
        if(aminoList[i] == "D"):
            colorList.append("#bfef45")
        if(aminoList[i] == "C"):
            colorList.append("#469990")
        if(aminoList[i] == "Q"):
            colorList.append("#9A6324")
        if(aminoList[i] == "E"):
            colorList.append("#ffe119")
        if(aminoList[i] == "G"):
            colorList.append("#aaffc3")
        if(aminoList[i] == "H"):
            colorList.append("#f032e6")
        if(aminoList[i] == "I"):
            colorList.append("#a9a9a9")
        if(aminoList[i] == "L"):
            colorList.append("#4363d8")
        if(aminoList[i] == "K"):
            colorList.append("#000075")
        if(aminoList[i] == "M"):
            colorList.append("#fffac8")
        if(aminoList[i] == "F"):
            colorList.append("#fabebe")
        if(aminoList[i] == "P"):
            colorList.append("#f58231")
        if(aminoList[i] == "S"):
            colorList.append("#e6beff")
        if(aminoList[i] == "T"):
            colorList.append("#911eb4")
        if(aminoList[i] == "W"):
            colorList.append("#e6194B")
        if(aminoList[i] == "Y"):
            colorList.append("#ffd8b1")
        if(aminoList[i] == "V"):
            colorList.append("#3cb44b")
        if(aminoList[i] == " "):
            colorList.append("#ffffff")
    return colorList

#draws a 10 x 30 pixel box for each amino acid and shades it in with the corresponding amino #acid color
def drawBoxes(colorList, speciesNum, name):
    #noStroke()
    counter = 0
    Xpoint = 275
    Ypoint = 30 * speciesNum
    for i in colorList:
        fill(i)
        rect(Xpoint, Ypoint, 10, 30)
        if(counter % 10 == 0):
            textSize(12)
            fill("#000000")
            text(counter, Xpoint, 150)
        fill("#000000")
        text(name,20, Ypoint + 25)
        Xpoint += 10
        counter += 1



#runner method
#calls findMutations() and drawBoxes() methods
def main():
    findMutations()
    colorList1 = aminoAcidToColor(AA1)
    colorList2 = aminoAcidToColor(AA2)
    colorList3 = aminoAcidToColor(AA3)
    drawBoxes(colorList1, 1, name1)
    drawBoxes(colorList2, 2, name2)
    drawBoxes(colorList3, 3, name3)
    
AA1 = codingStrandToAA(DNA1)
AA2 = codingStrandToAA(DNA2)
AA3 = codingStrandToAA(DNA3)

size(1500, 200)
background(255)
main()
