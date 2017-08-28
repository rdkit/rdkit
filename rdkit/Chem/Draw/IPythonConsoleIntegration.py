
"""RDKit Conformer Browser
Derived by Malitha Humayun Kabir from Greg Landrum's code as a part of GSoC 2017
Project : RDKit - 3Dmol.js integration
Mentors: Paul Czodrowski and Greg Landrum
Acknowledgement: Peter Gedeck reviewed, provided advice on restructuring, and wrote initial MolViewState class.
Date: 28rd August 2017
Email# malitha12345@gmail.com
"""

try:
    import py3Dmol
    _canUse3D = True
except ImportError:
    _canUse3D = False
    
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from ipywidgets import Layout, Label, Button, Box, HBox, VBox
from ipywidgets import Dropdown, HTML, Checkbox, Button
from IPython.display import display
from IPython.display import HTML as scriptHTML
import time
import sys
import copy

PY3 = sys.version_info[0] >= 3

BGCOLORS_3D = ('0x000000', '0xeeeeee', '0xffffff')

PROP_RDKIT = tuple(sorted(prop for prop, _ in Descriptors._descList))

DRAWING_LIGAND_3D=('line', 'cross', 'stick', 'sphere', 'surface', 'ballstick')

DRAWING_PROTEIN_3D=('line', 'cartoon', 'surface')

LIGAND_COLOR_SCHEME_3D=('default', 'greenCarbon', 'cyanCarbon', 'magentaCarbon',
                        'yellowCarbon', 'whiteCarbon', 'orangeCarbon', 'purpleCarbon', 
                        'blueCarbon', 'ssPyMOL', 'ssJmol', 'Jmol', 'amino', 
                        'shapely', 'nucleic', 'chain', 'chainHetatm', 'prop')

def Check3Dmolpromise():
    uniqueid = str(time.time()).replace('.','')
    input_form = """<div id="%s"></div>""" % uniqueid
    javascript = """
    <script type="text/Javascript">
        function checkLibraryLoad()
            {
                if(typeof $3Dmolpromise === 'undefined')
                {
                    res = '3DMol.js not loaded'
                } else {
                    res = '3DMol.js already loaded'
                }
                return res
            }
        document.getElementById('%s').innerHTML = checkLibraryLoad()
    </script>
    """ % uniqueid
    return(scriptHTML(input_form + javascript))
    
    
# I think this function is not required
# User should be able to supply appropriate input
# A notebook showing 3Dmol visualization should be present in rdkit doc
def ProcessLigandDict(ligandDict, keyForParentMol = 'parent'):
    """This function adds another key to the dictionary of the molecule."""
    
    if isinstance(ligandDict, dict) is False:
        raise TypeError("ligandDict must be a dictionary")
        
    keys = list(ligandDict.keys())
    firstKey = keys[0]
    
    if isinstance(ligandDict[firstKey], dict):
        raise TypeError("ProcessLigandDict doesn't support nested ligandDict")
        
    newLigandDict = {}
    
    for molId in keys:
        newLigandDict[molId] = {}
        newLigandDict[molId][keyForParentMol] = ligandDict[molId]
        
    return newLigandDict
    
def MinimizeLigand(ligandDict, 
                   keyForParentMol = 'parent', keyForMinimizedMol = 'minimized',
                   energyDataKey = 'energy',
                   molAndConfIds = 'allConfs', 
                   ff = 'UFF', 
                   maxIters = 50):
    
    """ This function takes a dictionary of ligand and does energy minimization to the ligands"""
    
    if ff not in ('MMFF', 'UFF'):
        raise TypeError("ff can be either MMFF or UFF")
        
    if isinstance(ligandDict, dict) is False:
        raise TypeError("ligandDict must be dict")
        
    if all(isinstance(key, str if PY3 else basestring) for key in ligandDict.keys()) is False:
        raise TypeError("keys of ligandDict must be str")
        
    
    if molAndConfIds is 'allConfs':
        molAndConfIdsTuple = ()
        for molId in list(ligandDict.keys()):
            mol = ligandDict[molId][keyForParentMol]
            confIds = list(range(mol.GetNumConformers()))
            for cid in confIds:
                molAndConfIdsTuple = molAndConfIdsTuple + ((molId, cid),)
    else:
        molAndConfIdsTuple = molAndConfIds
        
    
    if all(isinstance(molId, str if PY3 else basestring) for (molId, confId) in molAndConfIdsTuple) is False:
        raise TypeError("keys of ligandDict must be str")
        
    
    for Ids in molAndConfIdsTuple:
        
        molId = Ids[0]
        confId = Ids[1]
        
        energy = ligandDict[molId].setdefault(energyDataKey, {})
        
        if keyForMinimizedMol not in list(ligandDict[molId].keys()):
            oldMol = ligandDict[molId][keyForParentMol]
            newMol = copy.deepcopy(oldMol)
            ligandDict[molId][keyForMinimizedMol] = newMol
        else:
            newMol = ligandDict[molId][keyForMinimizedMol]
            # print('molecule already exists with given key' + keyForMinimizedMol + 'using existing molecule')
            
        if ff == 'MMFF':
            getFF = AllChem.MMFFGetMoleculeForceField(newMol, confId=confId)
            getFF.Minimize(maxIts = maxIters)
            energy[confId] = getFF.CalcEnergy()
            
        elif ff == 'UFF':
            getFF = AllChem.UFFGetMoleculeForceField(newMol, confId=confId)
            getFF.Minimize(maxIts = maxIters)
            energy[confId] = getFF.CalcEnergy()
        else:
            raise ValueError("invalid conformer id")
            
        
        
    return ligandDict
    
    
def AddPropToLigandDict(ligandDict, keyForParentMol = 'parent'):
    """ Add property to the mol """
    for molId in ligandDict:
        mol = ligandDict[molId][keyForParentMol]
        for prop_name in PROP_RDKIT:
            calculator = Descriptors.__dict__[prop_name]
            mol.SetProp(prop_name, str(calculator(mol)))
            
    return ligandDict
    
class MolViewState(object):
    def __init__(self, 
                 molecules, 
                 protein, 
                 keyForParentMol = 'parent', 
                 keyForMinimizedMol = 'minimized',
                 energyDataKey = 'energy'):
        """ Molecules is dictionary of molecules """
        
        self.ligandDict = molecules
        self.protein = protein
        self.keyForParentMol = keyForParentMol
        self.keyForMinimizedMol = keyForMinimizedMol
        self.energyDataKey = energyDataKey
        # These should have reasonable initial values
        self.rdkit_mol_select = set()
        self.rdkit_conf_select = set()
        self.allConfIDs = []
        
    def selectMolecules(self, selectAllMols, selectMultiMols, selectMol):
        """ Select either all moleculs or add selectMol or show only selectMol """
        if selectAllMols:
            self.rdkit_mol_select = set(self.ligandDict)
        elif selectMultiMols:
            self.rdkit_mol_select.add(selectMol)
        else:
            self.rdkit_mol_select = {selectMol}
            
    def selectConformations(self, selectAllConfs, selectMultiConfs, selectConf):
        """ For all selected molecules, select either all conformations or add selectConf or show only selectConf """
        
        for molIds in self.selectedMolNames:
            
            mol = self.ligandDict[molIds][self.keyForParentMol]
            nConformers = mol.GetNumConformers()
            
            if nConformers > 1:
                if selectAllConfs:
                    self.rdkit_conf_select = set(range(nConformers))
                elif selectMultiConfs:
                    self.rdkit_conf_select.add(selectConf)
                else:
                    self.rdkit_conf_select = {selectConf}
            elif mol.GetNumConformers() == 1:
                self.rdkit_conf_select = {0}
                
    def generateIds(self):
        """generate paired id tuple from rdkit_mol_select and rdkit_conf_select"""
        Ids = []
        for molId in self.rdkit_mol_select:
            for confId in self.rdkit_conf_select:
                Ids.append((molId,confId))
        self.idPaired = Ids
        
    def update_rdkit_select_unique_ids(self):
        """generate unique molId and confId from paired id tuple"""
        molIds = []
        confIds = []
        for (molId, confId) in self.idPaired:
            molIds.append(molId)
            confIds.append(confId)
        self.rdkit_mol_select = set(molIds)
        self.rdkit_conf_select = set(confIds)
        
    def parseIds(self, molAndConfIds):
        """ takes directly parsed paired id and store """
        self.idPaired = molAndConfIds
        self.update_rdkit_select_unique_ids()
        
        
    @property
    def selectedModels(self):
        """ Iterator over all selected models (molecules/conformations) """
        for (molId, confId) in self.idPaired:
            
            molData = self.ligandDict[molId]
            
            if self.keyForParentMol in list(molData.keys()) and self.keyForMinimizedMol in list(molData.keys()):
                
                yield {self.keyForParentMol: 
                       Chem.MolToMolBlock(molData[self.keyForParentMol], confId=confId),
                       
                       self.keyForMinimizedMol: 
                       Chem.MolToMolBlock(molData[self.keyForMinimizedMol], confId=confId)
                      }
            else:
                
                yield {self.keyForParentMol: 
                       Chem.MolToMolBlock(molData[self.keyForParentMol], confId=confId)
                      }
            
    @property
    def selectedMolNames(self):
        """ Return the names of all selected molecules """
        return self.rdkit_mol_select
        
    @property
    def selectedConfIds(self):
        """ Return the names of all selected confIds """
        return self.rdkit_conf_select
        
    @property
    def selectedMolecules(self):
        """ Return the selected molecules """
        return [self.ligandDict[name][self.keyForParentMol] for name in self.selectedMolNames]
        
    @property
    def allConfIds(self):
        """ Return the number of conformations - use the first selected molecule to determine """
        nconfIds = self.selectedMolecules[0].GetNumConformers()
        return list(range(nconfIds))
        
    @property
    def getPropPrecalculated(self):
        """ Return the precalculated properties """
        if len(self.selectedMolNames) != 1:
            return None
        else:
            mol = self.ligandDict[list(self.selectedMolNames)[0]][self.keyForParentMol]
             
            if len(mol.GetPropNames()) == 0:
                return 'Not found'
            else:
                return {prop: mol.GetProp(prop) for prop in mol.GetPropNames()}
                
    @property
    def getMinimizationEnergy(self):
        """ Return the minimization data """
        if len(self.selectedMolNames) != 1:
            return None
        else:
            molData = self.ligandDict[list(self.selectedMolNames)[0]]
            if self.energyDataKey in molData:
                return molData[self.energyDataKey]
            else:
                return 'energyDataKey not found'
                
    @property
    def atomLabel(self):
        """ Return the atomLabel """
        if len(self.selectedMolNames) == 1 and len(self.selectedConfIds) == 1:
            
            mol = self.ligandDict[list(self.selectedMolNames)[0]][self.keyForParentMol]
            sconf = mol.GetConformer(list(self.selectedConfIds)[0])
            xyz = sconf.GetPositions()
            atomLabel = []
            for i in list(range(sconf.GetNumAtoms())):
                label = ('{}{}'.format(mol.GetAtomWithIdx(i).GetSymbol(), i), 
                         xyz[i]
                        )
                atomLabel.append(label)
            return atomLabel
    
    @property
    def confLabel(self):
        """ Return the conformer label """
        if len(self.selectedMolNames) == 1 and len(self.selectedConfIds) == 1:
            confLabel = list(self.selectedMolNames)[0] + ':' + str(list(self.selectedConfIds)[0])
            return confLabel
            
            
        
        
class MolView3D(object):
    
    def __init__(self, 
                 ligandDict = None,
                 protein = None, 
                 keyForParentMol = 'parent', keyForMinimizedMol = 'minimized', energyDataKey = 'energy',
                 ligStyle = 'stick', protStyle='cartoon', emLigStyle = 'stick',
                 ligSelPanel = 'full',
                 stylePanel = None, 
                 labelPanel = False,
                 propertyPanel = False,
                 emPanel = False):
        """This function initiates required widgets and 3Dmol.js viewer"""
        
        if ligSelPanel not in ('full', 'minimal'):
            raise KeyError('ligSelPanel can be full or minimal')
        
        if stylePanel is not None:
            if stylePanel not in ('lig', 'em', 'emlig', 'prot', 'ligprot', 'emprot', 'emligprot'):
                raise KeyError('stylePanel can be: lig, em, emlig, prot, ligprot, emprot or emligprot')
            
        
        self.onStart = True
        
        self.ligStyle = ligStyle
        self.protStyle = protStyle
        self.emLigStyle = emLigStyle
        
        self.ligSelPanel = ligSelPanel # ligSelPanel can be full or minimal
        self.propertyPanel = propertyPanel
        self.stylePanel = stylePanel
        self.labelPanel = labelPanel
        self.emPanel = emPanel
        
        
        # Right hand panel (widgets)
        
        self.wgBox = []
        self.itemLayout=Layout(display='flex', flex_flow='row', justify_content='space-between')
        self.widgetsFor3DView = []
        
        # Molecule selection panel
        self.MoleculeSelectionPanelUI()
        
        # Molecule visibility (show/hide)
        self.MoleculeVisibilityPanelUI()
        
        # property show
        if propertyPanel:
            self.PropertyViewPanelUI()
        
        # minimization energy data show
        if emPanel:
            self.EnergyMinimizationDataViewPanelUI()
            
        # stylePanel
        if stylePanel is not None:
            if stylePanel:
                self.StylePanelUI()
        
        # add conf labels and atom labels 
        if labelPanel:
            self.LabelPanelUI()
        
        
        # background color
        self.BackgroundColorPanelUI()
        
        # zoomTo button
        self.ZoomToButtomUI()
        
        # Add observer to interactive widgets
        self.AddObserverToWidgets()
        
        # Rendering left (molecule viewer) and right (all widgets) panels
        self.RenderWidgetsWithViewer()
        
        # adding model to viewer
        self.SetMolData(ligandDict, protein, 
                        keyForParentMol, keyForMinimizedMol, energyDataKey, 
                        molAndConfIds = None)
        
                 
        
    def MoleculeSelectionPanelUI(self):
        """ligand selection widgets"""
        self.selectedMolsConfsView = HTML(description='', value='')
        self.wgBox.append(Box([Label(value='ids'),self.selectedMolsConfsView], layout=self.itemLayout))
        
        self.molId = Dropdown(description='', options=['select'],value='select')
        self.confId = Dropdown(description='', options=['select'], value='select')
        
        if self.ligSelPanel == 'full':
            self.selectMultiMols = Checkbox(description='selectMultiMols', value=False)
            self.selectAllMols = Checkbox(description='selectAllMols', value=False)
            self.selectMultiConfs = Checkbox(description='selectMultiConfs', value=False)
            self.selectAllConfs = Checkbox(description='selectAllConfs', value=False)
            
            box_for_multimol_allmol = HBox([self.selectMultiMols, self.selectAllMols])
            cbMolSelect=Box([Label(value=''), box_for_multimol_allmol], layout=self.itemLayout)
            
            box_for_multiconf_allconf = HBox([self.selectMultiConfs, self.selectAllConfs])
            cbConfSelect=Box([Label(value=''), box_for_multiconf_allconf], layout=self.itemLayout)
            
            self.wgBox.append(Box([Label(value='molId'),self.molId], layout=self.itemLayout))
            self.wgBox.append(cbMolSelect)
            self.wgBox.append(Box([Label(value='confId'),self.confId], layout=self.itemLayout))
            self.wgBox.append(cbConfSelect)
            
            self.widgetsFor3DView.extend(['molId', 'selectMultiMols', 'selectAllMols', 
                                          'confId', 'selectMultiConfs', 'selectAllConfs'])
        else:
            self.wgBox.append(Box([Label(value='molId'),self.molId], layout=self.itemLayout))
            self.wgBox.append(Box([Label(value='confId'),self.confId], layout=self.itemLayout))
            self.widgetsFor3DView.extend(['molId', 'confId'])
            
    def MoleculeVisibilityPanelUI(self):
        """ligand, protein, and energy minimized ligand visibility widgets"""
        self.ligandVisible = Checkbox(description='ligandVisible', value=True)
        self.wgBox.append(Box([Label(value=''),self.ligandVisible], layout=self.itemLayout))
        self.widgetsFor3DView.append('ligandVisible')
        
        self.proteinVisible = Checkbox(description='proteinVisible', value=True)
        self.wgBox.append(Box([Label(value=''),self.proteinVisible], layout=self.itemLayout))
        self.widgetsFor3DView.append('proteinVisible')
        
        self.emLigandVisible = Checkbox(description='emLigandVisible', value=True)
        self.wgBox.append(Box([Label(value=''),self.emLigandVisible], layout=self.itemLayout))
        self.widgetsFor3DView.append('emLigandVisible')
        
    def PropertyViewPanelUI(self):
        """widgets for ligand property selection and property show"""
        self.prop_view = HTML(description='', value='initializing...')
        self.prop_wg = Dropdown(description='',options=['select'], value='select')
        self.wgBox.append(Box([Label(value='prop'),self.prop_view], layout=self.itemLayout))
        self.wgBox.append(Box([Label(value='prop'),self.prop_wg], layout=self.itemLayout))
        self.widgetsFor3DView.append('prop_wg')
        
        
    def EnergyMinimizationDataViewPanelUI(self):
        """widgets for showing energy data for each conformer"""
        self.energy_wg = HTML(description='', value='initializing...')
        self.wgBox.append(Box([Label(value=''),self.energy_wg], layout=self.itemLayout))
        self.widgetsFor3DView.append('energy_wg')
        
        
        
    def LigStylePanelUI(self):
        """widgets for selecting ligand representation and coloring scheme"""
        self.ligStyle_wg = Dropdown(description='', options=DRAWING_LIGAND_3D, value=self.ligStyle)
        self.colorScheme = Dropdown(description='', options=LIGAND_COLOR_SCHEME_3D, value='default')
        self.wgBox.append(Box([Label(value='ligStyle'),self.ligStyle_wg], layout=self.itemLayout))
        self.widgetsFor3DView.append('ligStyle_wg')
        self.wgBox.append(Box([Label(value='ligColor'),self.colorScheme], layout=self.itemLayout))
        self.widgetsFor3DView.append('colorScheme')
        
    def EMLigStylePanelUI(self):
        """widgets for selecting energy minimized ligand representation and coloring scheme"""
        self.emLigStyle_wg = Dropdown(description='', options=DRAWING_LIGAND_3D, value=self.emLigStyle)
        self.emColorScheme = Dropdown(description='', options=LIGAND_COLOR_SCHEME_3D, value='whiteCarbon')
        self.wgBox.append(Box([Label(value='emLigStyle'),self.emLigStyle_wg], layout=self.itemLayout))
        self.widgetsFor3DView.append('emLigStyle_wg')
        self.wgBox.append(Box([Label(value='emLigColor'),self.emColorScheme], layout=self.itemLayout))
        self.widgetsFor3DView.append('emColorScheme')
        
    def ProteinStylePanelUI(self):
        """widgets for selecting protein representation"""
        self.protStyle_wg = Dropdown(description='', options=DRAWING_PROTEIN_3D, value=self.protStyle)
        self.wgBox.append(Box([Label(value='protStyle'),self.protStyle_wg], layout=self.itemLayout))
        self.protStyle_helicesAsTubes_wg = Checkbox(description='helicesAsTubes', value=False)
        self.wgBox.append(Box([Label(value=''),self.protStyle_helicesAsTubes_wg], layout=self.itemLayout))
        self.widgetsFor3DView.extend(['protStyle_wg', 'protStyle_helicesAsTubes_wg'])
        
    def StylePanelUI(self):
        """widgets for selecting which style panel will be appeared"""
        if self.stylePanel == 'lig':
            self.LigStylePanelUI()
        if self.stylePanel == 'em':
            self.EMLigStylePanelUI()
        if self.stylePanel == 'emlig':
            self.LigStylePanelUI()
            self.EMLigStylePanelUI()
        if self.stylePanel == 'prot':
            self.ProteinStylePanelUI()
        if self.stylePanel == 'ligprot':
            self.LigStylePanelUI()
            self.ProteinStylePanelUI()
        if self.stylePanel == 'emprot':
            self.EMLigStylePanelUI()
            self.ProteinStylePanelUI()
        if self.stylePanel == 'emligprot':
            self.LigStylePanelUI()
            self.ProteinStylePanelUI()
            self.EMLigStylePanelUI()
            
    def LabelPanelUI(self):
        """widgets for labeling ligand (not energy minimized ligand)"""
        self.confLabel = Checkbox(description='confLabel', value=False)
        self.atomLabel = Checkbox(description='atomLabel', value=False)
        box_for_labeling = HBox([self.confLabel,self.atomLabel])
        cbLabel=Box([Label(value=''), box_for_labeling], layout=self.itemLayout)
        self.wgBox.append(cbLabel)
        self.widgetsFor3DView.extend(['confLabel', 'atomLabel'])
        
    def BackgroundColorPanelUI(self):
        """widgets for selecting viewer background"""
        self.background = Dropdown(description='', options=BGCOLORS_3D, value=BGCOLORS_3D[0])
        self.wgBox.append(Box([Label(value='background'),self.background], layout=self.itemLayout))
        self.widgetsFor3DView.append('background')
        
    def ZoomToButtomUI(self):
        """widgets for zoomTo button"""
        self.zoomTo = Button(description="zoomTo", button_style='success')
        box_for_zoomTo = Box([Label(value=''),HBox([self.zoomTo])],layout=self.itemLayout)
        self.wgBox.append(box_for_zoomTo)
        
    def AddObserverToWidgets(self):
        """This function sets observer to all the interactive widgets except zoomTo"""
        self.zoomTo.on_click(self.handle_zoomTo_button)
        
        for i in self.widgetsFor3DView:
            selectedWidget =self.__dict__[i]
            selectedWidget.observe(self.handle_change, names='value')
            
    def RenderWidgetsWithViewer(self):
        """this creates a html table having an id as and then insert the 3DMol.js viewer in table"""
        tableUID = str(time.time()).replace('.','')
        
        # left panel (container holding table for viewer)
        size = (435, 485)
        viewerLayout=Layout(width=str(size[0]+5)+'px', height=str(size[1]+5)+'px', border='solid')
        wgLeftWidget = HTML('''<table><tr><td id="%s"></td></tr></table>'''% tableUID, layout=viewerLayout)
        
        # right panel
        wgRightBox=VBox(self.wgBox, layout=Layout(border='solid', width='45%'))
        
        # Combining left and right panel
        viewer = HBox([wgLeftWidget, wgRightBox])
        
        # displaying everything
        display(viewer)
        
        # inserting 3DMol.js viewer in existing container (table)
        self.view = py3Dmol.view(width=size[0],height=size[1])
        self.view.setBackgroundColor('0x000000')
        self.view.zoomTo()
        
        display(self.view.insert(tableUID))
        
        
    def SetMolData(self, 
                   ligandDict = None, 
                   protein = None, 
                   keyForParentMol = 'parent', 
                   keyForMinimizedMol = 'minimized',
                   energyDataKey = 'energy',
                   molAndConfIds = None):
        """This function sets ligand dictionary, protein, and dict keys and initiates MolViewState class"""
        
        if isinstance(ligandDict, dict) is False:
            raise TypeError("ligandDict must be dict")
            
        if all(isinstance(key, str if PY3 else basestring) for key in ligandDict.keys()) is False:
            raise TypeError("keys of ligandDict must be str")
            
        
        self.keyForParentMol = keyForParentMol
        self.keyForMinimizedMol = keyForMinimizedMol
        self.energyDataKey = energyDataKey
        
        self.molViewState = MolViewState(ligandDict, protein, 
                                         keyForParentMol, keyForMinimizedMol, energyDataKey)
        
        if self.molViewState.ligandDict is not None:
            
            keys = list(self.molViewState.ligandDict.keys())
            self.onStartIds = ((keys[0], 0),)
            
            self.molId.options = keys
            self.molId.value = keys[0]
            self.confId.options = list(range(9))
            self.confId.value=0
            
            self.selectedMolsConfsView.value = ', '.join([str(x) for x in self.onStartIds])
            
    
    def handle_change(self, change):
        """This function handles all the interactive widgets except buttons and 3Dmol.js viewer"""
        
        self.SelectMolAndConf()
        
    
    def handle_zoomTo_button(self, b):
        """This function handles zoomTo button"""
        self.view.zoomTo()
        display(self.view.update())
        
    
    def SelectMolAndConf(self, molAndConfIds = None):
        """ instantiates mol and conformer selection function of the MolViewState"""
        if self.molViewState.ligandDict is None:
            raise TypeError("please provide ligandDict before selection")
        
        if self.onStart:
            molAndConfIds = self.onStartIds
        
        if molAndConfIds is None and self.ligSelPanel == 'full':
            
            self.molViewState.selectMolecules(self.selectAllMols.value, 
                                              self.selectMultiMols.value,
                                              self.molId.value)
            
            self.molViewState.selectConformations(self.selectAllConfs.value, 
                                                  self.selectMultiConfs.value,
                                                  self.confId.value)
            
            self.molViewState.generateIds()
            
            self.confId.options = self.molViewState.allConfIds
            
            self.selectedMolsConfsView.value = ', '.join([str(x) for x in self.molViewState.idPaired])
            
        elif molAndConfIds is None and self.ligSelPanel == 'minimal':
            
            self.molViewState.selectMolecules(False, False, self.molId.value)
            self.molViewState.selectConformations(False, False, self.confId.value)
            
            self.molViewState.generateIds()
            
            self.confId.options = self.molViewState.allConfIds
            
            self.selectedMolsConfsView.value = ', '.join([str(x) for x in self.molViewState.idPaired])
            
        else:
            self.molViewState.parseIds(molAndConfIds = molAndConfIds)
            
            if self.ligSelPanel == 'full' and len(self.molViewState.selectedMolNames) > 1:
                self.selectMultiMols.value = True
            if self.ligSelPanel == 'full' and len(self.molViewState.selectedConfIds) > 1:
                self.selectMultiConfs.value = True
                
            #For multiple molecules with unequal number of conformers, it is hard to determine acceptable confIds
            self.confId.options = self.molViewState.allConfIds
            
            self.selectedMolsConfsView.value = ', '.join([str(x) for x in self.molViewState.idPaired])
            
        self.render3D()
        
        
    def ShowLigandProperty(self):
        """ Handles property in render3D function """
        preCalcProp = self.molViewState.getPropPrecalculated
        if isinstance(preCalcProp, dict):
            self.prop_wg.options = sorted(list(preCalcProp.keys()))
            prop=self.prop_wg.value
            self.prop_view.value = prop + ' : ' + str(preCalcProp[prop])
        elif preCalcProp == 'Not found':
            self.prop_view.value = 'data not found!'
        else:
            self.prop_view.value = 'single molecule selection required!'
            
    def ShowMinimizationEnergy(self):
        """ Handles energy minimization data in render3D function """
        energy = self.molViewState.getMinimizationEnergy
        energyDataKey = self.molViewState.energyDataKey
        confId = self.confId.value
        
        if isinstance(energy, dict) and confId in energy:
            self.energy_wg.value = 'minimization data ('+ energyDataKey +') : ' + str(energy[confId])
        elif isinstance(energy, dict) and confId not in energy:
            self.energy_wg.value = 'data not found'
        elif energy is None:
            self.energy_wg.value = 'single molecule selection required'
        else:
            self.energy_wg.value = 'data not found'
            
    def AddLigandStyle(self, modelCategory, modelId):
        """ Handles ligand and energy minimized ligand drawing style and color in AddLigandWithStyle function """
        
        if modelCategory == self.keyForParentMol:
            
            ligStyle = self.ligStyle_wg.value if 'ligStyle_wg' in self.__dict__  else self.ligStyle
            color = self.colorScheme.value if 'colorScheme' in self.__dict__ else 'default'
            self.ligStyle = ligStyle
        else:
            ligStyle = self.emLigStyle_wg.value if 'emLigStyle_wg' in self.__dict__ else self.emLigStyle
            color = self.emColorScheme.value if 'emColorScheme' in self.__dict__ else 'whiteCarbon'
            self.emLigStyle = ligStyle
            
            
        if ligStyle == 'surface':
            self.view.addSurface('SES', {'model':modelId});
        elif ligStyle == 'ballstick':
            self.view.setStyle({'model':modelId},{'stick':{'radius':'0.2','colorscheme': color},
                                                  'sphere':{'radius':'0.4', 'colorscheme': color}
                                                 });
        else:
            self.view.setStyle({'model': modelId},{ligStyle:{'colorscheme': color}})
            
        
    def AddLigandLabels(self):
        """ Handles ligand labeling (conf label and atom label) in render3D function """
        
        if self.confLabel.value:
            if self.molViewState.confLabel is not None:
                self.view.addLabel(self.molViewState.confLabel, {'backgroundColor':'gray', 
                                                                 'fontColor':'white',
                                                                 'showBackground':'true',
                                                                 'alignment':'bottomCenter'})
        if self.atomLabel.value:
            if self.molViewState.atomLabel is not None:
                for (label,xyz) in self.molViewState.atomLabel:
                    self.view.addLabel(label, {'inFront' : 'false', 
                                               'fontSize' : '12',
                                               'fontColor':'gray',
                                               'showBackground':'false',
                                               'position' : {'x':xyz[0], 'y':xyz[1], 'z':xyz[2]}
                                              })
    
    
    def AddProteinStyle(self):
        """ Handles protein drawing style in AddProteinWithStyle function """
        # 'protStyle_wg' and 'protStyle_helicesAsTubes_wg' rendered as pair
        self.protStyle = self.protStyle_wg.value if 'protStyle_wg' in self.__dict__ else self.protStyle
        
        if self.protStyle == 'surface':
            self.view.addSurface('SES', {'model': self.proteinModelId});
            
        elif self.protStyle == 'line':
            self.view.setStyle({'model': self.proteinModelId},{'line':{}});
            
        elif self.protStyle == 'cartoon' and 'protStyle_helicesAsTubes_wg' not in self.__dict__:
            self.view.setStyle({'model': self.proteinModelId},{'cartoon':{'color': 'spectrum', 
                                                                          'arrows': 'true'}
                                                              })
        else:
            if self.protStyle_helicesAsTubes_wg.value:
                self.view.setStyle({'model': self.proteinModelId},{'cartoon':{'color': 'spectrum',
                                                                              'arrows': 'true', 
                                                                              'tubes' : 'true'}
                                                                  })
            else:
                self.view.setStyle({'model': self.proteinModelId},{'cartoon':{'color': 'spectrum', 
                                                                              'arrows': 'true'}
                                                                  })
                
    
    def AddLigandWithStyle(self):
        """ add ligand and energy minimized ligand in viewer (called in render3D function) """
        
        self.modelId = -1
        
        # add models (ligands)
        if 'ligandVisible' in self.__dict__ and self.ligandVisible.value:
            
            for models in self.molViewState.selectedModels:
                
                if self.keyForParentMol in list(models.keys()):
                    self.view.addModel(models[self.keyForParentMol], 'sdf')
                    self.modelId = self.modelId + 1
                    self.AddLigandStyle(self.keyForParentMol, self.modelId)
                    
                
        if 'emLigandVisible' in self.__dict__ and self.emLigandVisible.value:
            
            for models in self.molViewState.selectedModels:
                
                if self.keyForMinimizedMol in list(models.keys()):
                    self.view.addModel(models[self.keyForMinimizedMol], 'sdf')
                    self.modelId = self.modelId + 1
                    self.AddLigandStyle(self.keyForMinimizedMol, self.modelId)
                    
        
        
    def AddProteinWithStyle(self):
        """ add protein in viewer (called in render3D function) """
        if self.molViewState.protein is not None:
            if 'proteinVisible' in self.__dict__ and self.proteinVisible.value:
                pdb = Chem.MolToPDBBlock(self.molViewState.protein)
                self.view.addModel(pdb,'pdb')
                self.proteinModelId = self.modelId + 1
                self.AddProteinStyle()
                
    
    def render3D(self):
        """ This function updates the 3DMol.js viewer"""
        self.view.removeAllLabels()
        self.view.removeAllModels()
        self.view.removeAllSurfaces()
        
        self.view.setBackgroundColor(self.background.value)
        
        self.AddLigandWithStyle()
        
        # add label if required
        if self.labelPanel:
            if 'ligandVisible' in self.__dict__ and self.ligandVisible.value:
                self.AddLigandLabels()
                
        self.AddProteinWithStyle()
        
        # zoomTo does not work well for surface and label... so, zoomTo should not be default settings
        if self.onStart:
            self.view.zoomTo()
            self.onStart = False
            
        if self.propertyPanel:
            self.ShowLigandProperty()
            
        if self.emPanel:
            self.ShowMinimizationEnergy()
            
        display(self.view.update())
        
        
        
        