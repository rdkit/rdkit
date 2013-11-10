from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import SDWriter

def CSVToSDF(csvFile,outFile):
    csvFile_h = open(csvFile,'r')
    csvEntries = csvFile_h.readlines()
    csvFile_h.close()
    sdw = Chem.SDWriter(outFile)
    SDtags = csvEntries[0].split('\t')
    counter = 0
    for line in csvEntries[1:]:
        cpd = CSVToMOL(line,SDtags,counter)
        sdw.write(cpd)
        counter += 1  
    sdw.flush()
    sdw.close()

def _CSVToMOL(csvEntry,colNames,name=''):
    if type(csvEntry)==str:
        csvValues = csvEntry.split('\t')
        cpd = Chem.MolFromSmiles(csvValues[0])
        AllChem.Compute2DCoords(cpd)
        i = 1
        if name == '': name=str(i)
        cpd.SetProp("_Name",str(name))      
        for prop in csvValues[1:]:
            cpd.SetProp(colNames[i],prop)
            i += 1
            return cpd
    else:
        i = 1
        for key in colNames:
            if key=='Smiles':
#                print csvEntry[key]
                cpd = Chem.MolFromSmiles(str(csvEntry[key]))
                AllChem.Compute2DCoords(cpd)
                cpd.SetProp("_Name",str(i))
            else:
#                print key
 #               print csvEntry[key]
                cpd.SetProp(str(key),str(csvEntry[key]))
        return cpd

#def __CSVToMOL(csvEntry,colNames,name=''):
 #   import sys
  #  for key in colNames:
#	print csvEntry[key]
 #   sys.exit(1)

def CSVToMOL(csvData,colNames,name=''):
#    import sys
#    print type(csvData)
    if type(csvData) != str:
        result = []
#	print csvData
        for csvEntry in csvData:
#	    if type(csvEntry)==dict:
#		__CSVToMOL(csvEntry,colNames,name='')
#	    print type(csvEntry)
#	    print csvEntry
#	    sys.exit(1)
            result.append(_CSVToMOL(csvEntry,colNames,len(result)))
        return result
    else:
        return _CSVToMOL(csvData,colNames,name)
