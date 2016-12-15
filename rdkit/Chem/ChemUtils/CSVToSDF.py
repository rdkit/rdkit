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
    i = 1
    if name == '': name=str(i)
    if type(csvEntry)==str:
        csvValues = csvEntry.split('\t')
        cpd = Chem.MolFromSmiles(csvValues[0])
        AllChem.Compute2DCoords(cpd)
        cpd.SetProp("_Name",str(name))      
        for prop in csvValues[1:]:
            cpd.SetProp(colNames[i],prop)
            i += 1
            return cpd
    else:
        for key in colNames:
            if key=='Smiles':
                cpd = Chem.MolFromSmiles(str(csvEntry[key]))
                AllChem.Compute2DCoords(cpd)
                cpd.SetProp("_Name",str(name))
            else:
                cpd.SetProp(str(key),str(csvEntry[key]))
        return cpd

def CSVToMOL(csvData,colNames,name=''):
    if type(csvData) != str:
        result = []
        t = 0
        for csvEntry in csvData:
            result.append(_CSVToMOL(csvEntry,colNames,t))
            t += 1
        return result
    else:
        return _CSVToMOL(csvData,colNames,name)
