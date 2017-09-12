
"""MolView3D and a few utility functions
Project : RDKit - 3Dmol.js integration
Student: Malitha Humayun Kabir
Mentors: Paul Czodrowski and Greg Landrum
Acknowledgement: Peter Gedeck reviewed, provided advice on restructuring, and wrote initial MolViewState class. The starting point of MolView3D was Greg Landrum's codes on conformer browser.
Date: 12th September 2017
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

DRAWING_PROTEIN_3D=('line', 'cartoon', 'surface', 'cartoonWithTube')

COLOR_SCHEME_3D=('default', 'greenCarbon', 'cyanCarbon', 'magentaCarbon',
                 'yellowCarbon', 'whiteCarbon', 'orangeCarbon', 'purpleCarbon', 
                 'blueCarbon', 'ssPyMOL', 'ssJmol', 'Jmol', 'amino', 
                 'shapely', 'nucleic', 'chain', 'chainHetatm', 'prop')
# TODO: implement COLOR_3D
COLOR_3D = ('white', 'silver', 'gray', 'black', 'red', 'maroon', 'yellow', 
            'orange', 'olive', 'lime', 'green', 'aqua', 'cyan', 'teal', 'blue', 
            'navy', 'fuchsia', 'magenta', 'purple')

def Check3Dmolpromise():
    """This function checks whether 3Dmol.js loaded or not."""
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
    
def ProcessLigandDict(ligandDict, createChildKey = 'lig'):
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
        newLigandDict[molId][createChildKey] = ligandDict[molId]
    return newLigandDict
    
def CopyMolWithChildKey(ligandDict, copyFrom = 'lig', newChildKey = 'emLig'):
    """This function adds another key to the dictionary of the molecule."""
    if isinstance(ligandDict, dict) is False:
        raise TypeError("ligandDict must be a dictionary")
        
    keys = list(ligandDict.keys())
    
    for molId in keys:
        oldMol = ligandDict[molId][copyFrom]
        newMol = copy.deepcopy(oldMol)
        ligandDict[molId][newChildKey] = newMol
        
def ConfGenInLigandDict(ligandDict, onChildKey = 'emLig', numConfs = 10, molIds = 'allMols'):
    """This function adds another key to the dictionary of the molecule."""
    if molIds is 'allMols':
        molIds = list(ligandDict.keys())
        
    params = AllChem.ETKDG()
    params.numThreads=3
    for molId in molIds:
        AllChem.EmbedMultipleConfs(ligandDict[molId][onChildKey], numConfs=numConfs, params=params)
        
def MinimizeLigand(ligandDict, onChildKey = 'emLig',
                   createPerConfDataKey = 'energy',
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
            mol = ligandDict[molId][onChildKey]
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
        
        energy = ligandDict[molId].setdefault(createPerConfDataKey, {})
        mol = ligandDict[molId][onChildKey]
        
        if ff == 'MMFF':
            getFF = AllChem.MMFFGetMoleculeForceField(mol, confId=confId)
            getFF.Minimize(maxIts = maxIters)
            energy[confId] = getFF.CalcEnergy()
            
        elif ff == 'UFF':
            getFF = AllChem.UFFGetMoleculeForceField(mol, confId=confId)
            getFF.Minimize(maxIts = maxIters)
            energy[confId] = getFF.CalcEnergy()
        else:
            raise ValueError("invalid conformer id")
            
def AddPropToLigandDict(ligandDict, onChildKey = 'emLig'):
    """ Add property to the mol """
    for molId in ligandDict:
        mol = ligandDict[molId][onChildKey]
        for prop_name in PROP_RDKIT:
            calculator = Descriptors.__dict__[prop_name]
            mol.SetProp(prop_name, str(calculator(mol)))
            
class MolViewState(object):
    def __init__(self, 
                 molecules, 
                 protein, 
                 additionalChildKeysForMolRender = ['lig'],
                 childKeyForConfSelection = 'emLig',
                 childKeyForDataPerConf = 'energy',
                 childKeyForPropSelection = 'emLig'):
        """ Molecules is dictionary of molecules """
        
        self.ligandDict = molecules
        self.protein = protein
        
        self.additionalChildKeysForMolRender = additionalChildKeysForMolRender
        self.childKeyForConfSelection = childKeyForConfSelection
        
        self.childKeyForDataPerConf = childKeyForDataPerConf
        self.childKeyForPropSelection = childKeyForPropSelection
        # These should have reasonable initial values
        self.rdkit_mol_select = set()
        self.rdkit_conf_select = set()
        
    def selectAllMols(self):
        """ Select all molecules"""
        self.rdkit_mol_select = set(self.ligandDict)
        
    def selectAllConfs(self):
        """ For all selected molecules, select all conformations"""
        for molIds in self.selectedMolNames:
            
            mol = self.ligandDict[molIds][self.childKeyForConfSelection]
            nConformers = mol.GetNumConformers()
            
            if nConformers > 1:
                self.rdkit_conf_select = set(range(nConformers))
            else:
                self.rdkit_conf_select = {0}
                
    def selectSingleMol(self, molId):
        """ Select molecule"""
        self.rdkit_mol_select = {molId}
        
    def selectSingleConf(self, confId):
        """ Select conformation"""
        self.rdkit_conf_select = {confId}
        
    def generateIds(self):
        """generate paired id tuple from rdkit_mol_select and rdkit_conf_select"""
        self.idPaired = set()
        for molId in self.rdkit_mol_select:
            for confId in self.rdkit_conf_select:
                self.idPaired.add((molId,confId))
                
    def update_rdkit_select_unique_ids(self):
        """generate unique molId and confId from paired id tuple"""
        self.rdkit_mol_select = set()
        self.rdkit_conf_select = set()
        for (molId, confId) in self.idPaired:
            self.rdkit_mol_select.add(molId)
            self.rdkit_conf_select.add(confId)
            
    def deleteIds(self, molAndConfIds):
        """ takes directly parsed paired id and store """
        for ids in molAndConfIds:
            self.idPaired.discard(ids)
        self.update_rdkit_select_unique_ids()
        
    def addIds(self, molAndConfIds):
        """ takes directly parsed paired id and store """
        for ids in molAndConfIds:
            self.idPaired.add(ids)
        self.update_rdkit_select_unique_ids()
        
    def selectSingleId(self, molAndConfIds):
        """ takes directly parsed paired id and store """
        self.idPaired = set()
        for ids in molAndConfIds:
            self.idPaired.add(ids)
        self.update_rdkit_select_unique_ids()
        
    @property
    def selectedModels(self):
        """ Iterator over all selected models (molecules/conformations) """
        k_for_conf = self.childKeyForConfSelection
        k_for_additional = self.additionalChildKeysForMolRender
        
        modelDict = {}
        modelId = -1
        
        for (molId, confId) in self.idPaired:
            molData = self.ligandDict[molId]
            confKey = k_for_conf in list(molData.keys())
            if confKey:
                mol_with_conf = molData[k_for_conf]
                modelId = modelId + 1
                modelDict[modelId] = {}
                modelDict[modelId]['geom'] = Chem.MolToMolBlock(mol_with_conf,confId=confId)
                modelDict[modelId]['categoryForDrawing'] = k_for_conf
                modelDict[modelId]['molAndConfId'] = (molId, confId)
                modelDict[modelId]['molType'] = 'ligand'
                for key in k_for_additional:
                    additionalkeyLogical = key in list(molData.keys())
                    if additionalkeyLogical:
                        mol_additional = molData[key]
                        if confId in list(range(mol_additional.GetNumConformers())):
                            modelId = modelId + 1
                            modelDict[modelId] = {}
                            modelDict[modelId]['geom'] = Chem.MolToMolBlock(mol_additional,confId=confId)
                            modelDict[modelId]['categoryForDrawing'] = key
                            modelDict[modelId]['molAndConfId'] = (molId, confId)
                            modelDict[modelId]['molType'] = 'ligand'
        if self.protein is not None:
            modelId = modelId + 1
            modelDict[modelId] = {}
            modelDict[modelId]['geom'] = Chem.MolToPDBBlock(self.protein)
            modelDict[modelId]['categoryForDrawing'] = 'protein'
            modelDict[modelId]['molAndConfId'] = 'notRequired'
            modelDict[modelId]['molType'] = 'protein'
        yield modelDict
        
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
        return [self.ligandDict[name][self.childKeyForConfSelection] for name in self.selectedMolNames]
        
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
            mol = self.ligandDict[list(self.selectedMolNames)[0]][self.childKeyForPropSelection]
            
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
            if self.childKeyForDataPerConf in molData:
                return molData[self.childKeyForDataPerConf]
            else:
                return 'childKeyForDataPerConf not found'
                
    @property
    def atomLabel(self):
        """ Return the atomLabel """
        if len(self.selectedMolNames) == 1 and len(self.selectedConfIds) == 1:
            
            mol = self.ligandDict[list(self.selectedMolNames)[0]][self.childKeyForConfSelection]
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
                 childKeysForMolRender = ['lig', 'emLig'], 
                 childKeyForConfSelection = 'emLig',
                 childKeyForDataPerConf = 'energy',
                 childKeyForPropSelection = 'emLig',
                 ligStyle = ['stick', 'stick'], protStyle='cartoon',
                 ligColor = ['default','orangeCarbon'], protColor='ssJmol',
                 molVisibilityPanel = True, 
                 stylePanel = True, 
                 labelPanel = False,
                 propertyPanel = False,
                 emPanel = False,
                 molAndConfIds = None):
        """This function initiates required widgets and 3Dmol.js viewer"""
        
        self.onStart = True
        self.lastMoldelId = -1
        self.model_prop = {}
        self.model_temp = {}
        
        if ligandDict is not None:
            keys = list(ligandDict.keys())
            keysChild = list(ligandDict[keys[0]].keys())
            if len(keysChild) == 1:
                childKeysForMolRender = keysChild
                
        if type(childKeysForMolRender) is not list:
            raise TypeError('childKeysForMolRender must be a list')
            
        if len(childKeysForMolRender) == 1:
            childKeyForConfSelection = childKeysForMolRender[0]
            childKeyForPropSelection = childKeysForMolRender[0]
            
        self.childKeysForMolRender = childKeysForMolRender
        self.childKeyForConfSelection = childKeyForConfSelection
        
        self.childKeyForDataPerConf = childKeyForDataPerConf
        self.childKeyForPropSelection = childKeyForPropSelection
        
        self.ligStyleScheme = {}
        self.ligColorScheme = {}
        for (ids,lig) in enumerate(childKeysForMolRender):
            self.ligStyleScheme[lig] = ligStyle[ids]
            self.ligColorScheme[lig] = ligColor[ids]
            
        self.protStyle = protStyle
        self.protColorScheme = protColor
        
        self.molVisibilityPanel = molVisibilityPanel
        self.propertyPanel = propertyPanel
        self.stylePanel = stylePanel
        self.labelPanel = labelPanel
        self.emPanel = emPanel
        
        # Right hand panel (widgets)
        self.style = {'description_width': 'initial'}
        self.rightBoxWidth = '45%'
        self.wg_height = '30px'
        self.wgBox = []
        self.itemLayout=Layout(display='flex', flex_flow='row', 
                               justify_content='flex-start', 
                               height=self.wg_height)
        
        # Molecule selection panel
        self.MoleculeSelectionPanelUI()
        
        # Molecule visibility (show/hide)
        if self.molVisibilityPanel:
            self.MoleculeVisibilityPanelUI(ligandDict, protein)
            
        # property show
        if propertyPanel:
            self.PropertyViewPanelUI()
        
        # minimization data show
        if self.emPanel:
            self.EnergyMinimizationDataViewPanelUI()
            
        # stylePanel
        if stylePanel:
            if ligandDict is not None and protein is not None:
                self.stylePanel = 'ligprot'
            elif ligandDict is not None:
                self.stylePanel = 'lig'
            else:
                self.stylePanel = 'prot'
            self.StylePanelUI()
            
        # add conf labels and atom labels 
        if labelPanel:
            self.LabelPanelUI()
            
        # background color
        self.BackgroundColorPanelUI()
        
        # Add observer to interactive widgets
        self.AddObserverToWidgets()
        
        # Rendering left (molecule viewer) and right (all widgets) panels
        self.RenderWidgetsWithViewer()
        
        # adding model to viewer
        self.SetMolData(ligandDict, protein, 
                        childKeysForMolRender, 
                        childKeyForConfSelection, 
                        childKeyForDataPerConf,
                        childKeyForPropSelection,
                        molAndConfIds)
        
    def MoleculeSelectionPanelUI(self):
        """ligand selection widgets"""
        self.wg_selectedMolsConfsView = HTML(description='', value='')
        b1 = Box([Label(value='ids'),self.wg_selectedMolsConfsView], layout=self.itemLayout)
        
        tempLayout = Layout(display='flex', justify_content='flex-start', width='150px', height=self.wg_height)
        
        self.wg_selectConfIdfrom = Dropdown(description='', options=self.childKeysForMolRender,
                                            value=self.childKeyForConfSelection, layout=tempLayout)
        self.wg_updateConfId = Button(description="updateConfId", button_style='success',
                                      layout=Layout(width='100px'))
        box_for_confId_selection = HBox([self.wg_selectConfIdfrom, self.wg_updateConfId])
        b2 = Box([Label(value='confId from'),box_for_confId_selection], layout=self.itemLayout)
        
        self.wg_selectAllMols = Checkbox(description='selectAllMols', value=False)
        self.wg_selectAllConfs = Checkbox(description='selectAllConfs', value=False)
        box_for_allmol_allconf = HBox([self.wg_selectAllMols, self.wg_selectAllConfs])
        b3 = Box([Label(value=''), box_for_allmol_allconf], layout=self.itemLayout)
        
        molIdLayout = Layout(display='flex', justify_content='flex-start', width='200px', height=self.wg_height)
        confIdLayout = Layout(display='flex', justify_content='flex-start',width='140px', height=self.wg_height)
        
        self.wg_molId = Dropdown(description='molId', options=['select'],value='select', layout=molIdLayout)
        self.wg_confId = Dropdown(description='confId', options=['select'], value='select', layout=confIdLayout)
        b4 = Box([HBox([self.wg_molId, self.wg_confId])],layout=self.itemLayout)
        
        self.wg_modelDelete = Button(description='delete', button_style='success', layout=Layout(width='100px'))
        self.wg_modelAdd = Button(description='add', button_style='success', layout=Layout(width='100px'))
        self.wg_modelSelect = Button(description='select', button_style='success', layout=Layout(width='100px'))
        b5 = Box([Label(value=''),HBox([self.wg_modelDelete, self.wg_modelAdd, self.wg_modelSelect])
                 ],layout=self.itemLayout)
        
        self.wg_zoomTo = Button(description="zoomTo", button_style='success', layout=Layout(width='150px'))
        self.wg_reDraw = Button(description='reDraw', button_style='success', layout=Layout(width='150px'))
        b6 = Box([Label(value=''),HBox([self.wg_zoomTo, self.wg_reDraw])],layout=self.itemLayout)
        
        molConfBox = VBox([b1, b2, b3, b4, b5, b6], layout=Layout(border='solid', border_color = 'blue'))
        
        self.wgBox.append(molConfBox)
        
    def MoleculeVisibilityPanelUI(self, ligandDict, protein):
        """ligand, protein, and energy minimized ligand visibility widgets"""
        protein_sh = True if protein is not None else False
        if protein_sh:
            self.wg_proteinVisible = Checkbox(description='proteinVisible', value=True)
            self.wgBox.append(Box([Label(value=''),self.wg_proteinVisible], layout=self.itemLayout))
            
        ligand_sh = True if ligandDict is not None else False    
        if ligand_sh:
            for name in self.childKeysForMolRender:
                wgName = 'wg_' + name + 'Visible'
                self.__dict__[wgName] = Checkbox(description = wgName[3:], value=True)
                self.wgBox.append(Box([Label(value=''), self.__dict__[wgName]], layout=self.itemLayout))
                    
    def PropertyViewPanelUI(self):
        """widgets for ligand property selection and property show"""
        self.wg_propView = HTML(description='prop', value='initializing...', style = self.style)
        self.wgBox.append(Box([Label(value=''),self.wg_propView], layout=self.itemLayout))
        
        self.wg_propSelect = Dropdown(description='prop',options=['select'], value='select', style = self.style)
        self.wgBox.append(Box([Label(value=''),self.wg_propSelect], layout=self.itemLayout))
        
    def EnergyMinimizationDataViewPanelUI(self):
        """widgets for showing energy data for each conformer"""
        self.wg_energyView = HTML(description='', value='initializing...')
        self.wgBox.append(Box([Label(value=''),self.wg_energyView], layout=self.itemLayout))
        
    def LigStylePanelUI(self):
        """widgets for selecting ligand representation and coloring scheme"""
        for (ids,name) in enumerate(self.childKeysForMolRender):
            
            wgNameForStyle = 'wg_' + name + 'Style'
            
            self.__dict__[wgNameForStyle] = Dropdown(description=wgNameForStyle[3:], 
                                                     options=DRAWING_LIGAND_3D,
                                                     value=self.ligStyleScheme[name], style = self.style)
            self.wgBox.append(Box([Label(value=''), 
                                   self.__dict__[wgNameForStyle]
                                  ], layout=self.itemLayout))
            
            wgNameForColor = 'wg_' + name + 'Color'
            
            self.__dict__[wgNameForColor] = Dropdown(description=wgNameForColor[3:], 
                                                     options=COLOR_SCHEME_3D,
                                                     value=self.ligColorScheme[name], style = self.style)
            self.wgBox.append(Box([Label(value=''), 
                                   self.__dict__[wgNameForColor]
                                  ], layout=self.itemLayout))
            
    def ProteinStylePanelUI(self):
        """widgets for selecting protein representation"""
        self.wg_protStyle = Dropdown(description='protStyle', 
                                     options=DRAWING_PROTEIN_3D, 
                                     value=self.protStyle,
                                     style = self.style)
        self.wgBox.append(Box([Label(value=''),self.wg_protStyle], layout=self.itemLayout))
        
        self.wg_protColor = Dropdown(description='protColor', 
                                     options=COLOR_SCHEME_3D, 
                                     value=self.protColorScheme, 
                                     style = self.style)
        self.wgBox.append(Box([Label(value=''),self.wg_protColor], layout=self.itemLayout))
        
    def StylePanelUI(self):
        """widgets for selecting which style panel will be appeared"""
        if self.stylePanel == 'lig':
            self.LigStylePanelUI()
        if self.stylePanel == 'prot':
            self.ProteinStylePanelUI()
        if self.stylePanel == 'ligprot':
            self.LigStylePanelUI()
            self.ProteinStylePanelUI()
            
    def LabelPanelUI(self):
        """widgets for labeling ligand (not energy minimized ligand)"""
        self.wg_confLabel = Checkbox(description='confLabel', value=False)
        self.wg_atomLabel = Checkbox(description='atomLabel', value=False)
        box_for_labeling = HBox([self.wg_confLabel,self.wg_atomLabel])
        self.wgBox.append(Box([Label(value=''), box_for_labeling], layout=self.itemLayout))
        
    def BackgroundColorPanelUI(self):
        """widgets for selecting viewer background"""
        self.wg_background = Dropdown(description='background', 
                                      options=BGCOLORS_3D, 
                                      value=BGCOLORS_3D[0], style = self.style)
        self.wgBox.append(Box([Label(value=''),self.wg_background], layout=self.itemLayout))
        
    def AddObserverToWidgets(self):
        """This function sets observer to all the interactive widgets except zoomTo"""
        self.wg_updateConfId.on_click(self.handle_updateConfId_button)
        self.wg_zoomTo.on_click(self.handle_zoomTo_button)
        self.wg_reDraw.on_click(self.handle_reDraw_button)
        self.wg_modelDelete.on_click(self.handle_modelDelete_button)
        self.wg_modelAdd.on_click(self.handle_modelAdd_button)
        self.wg_modelSelect.on_click(self.handle_modelSelect_button)
        self.wg_selectAllMols.observe(self.handle_selectAllMols_checkbox, names='value')
        self.wg_selectAllConfs.observe(self.handle_selectAllConfs_checkbox, names='value')
        
    def RenderWidgetsWithViewer(self):
        """this creates a html table having an id as and then insert the 3DMol.js viewer in table"""
        tableUID = str(time.time()).replace('.','')
        
        # left panel (container holding table for viewer)
        size = (435, 485)
        viewerLayout=Layout(width=str(size[0]+5)+'px', height=str(size[1]+5)+'px', border='solid')
        wg_leftWidget = HTML('''<table><tr><td id="%s"></td></tr></table>'''% tableUID, layout=viewerLayout)
        
        # right panel
        wg_rightBox=VBox(self.wgBox, layout=Layout(border='solid', width=self.rightBoxWidth))
        
        # Combining left and right panels
        viewer = HBox([wg_leftWidget, wg_rightBox])
        
        # displaying everything
        display(viewer)
        
        # inserting 3DMol.js viewer in existing container (table)
        self.view = py3Dmol.view(width=size[0],height=size[1])
        self.view.setBackgroundColor('0x000000')
        self.view.zoomTo()
        
        display(self.view.insert(tableUID))
        
    def handle_updateConfId_button(self, b):
        """This function handles modelDelete button"""
        self.UpdateConfId(childKeyForConfSelection = self.wg_selectConfIdfrom.value)
        
    def handle_reDraw_button(self, b):
        """This function handles reDraw button"""
        self.previous_state = self.model_prop
        self.Render3D()
        
    def handle_zoomTo_button(self, b):
        """This function handles zoomTo button"""
        self.view.zoomTo()
        display(self.view.update())
        
    def handle_modelDelete_button(self, b):
        """This function handles modelDelete button"""
        molAndConfIds = ((self.wg_molId.value, self.wg_confId.value),)
        self.DeleteMolAndConf(molAndConfIds = molAndConfIds)
        
    def handle_modelAdd_button(self, b):
        """This function handles modelAdd button"""
        molAndConfIds = ((self.wg_molId.value, self.wg_confId.value),)
        self.AddMolAndConf(molAndConfIds = molAndConfIds)
        
    def handle_modelSelect_button(self, b):
        """This function handles modelSelect button"""
        molAndConfIds = ((self.wg_molId.value, self.wg_confId.value),)
        self.SelectMolAndConf(molAndConfIds = molAndConfIds)
    
    def handle_selectAllMols_checkbox(self, change):
        """This function handles selectAllMols checkbox"""
        self.SelectAllMols(confId = self.wg_confId.value, callFrom = 'widget')
        
    def handle_selectAllConfs_checkbox(self, change):
        """This function handles selectAllConfs checkbox"""
        self.SelectAllConfs(molId = self.wg_molId.value, callFrom = 'widget')
        
    def updateWidgets(self, callFrom):
        """ update mol and conf id viewer"""
        if callFrom == 'SelectAllMols':
            self.molViewState.generateIds()
            
        if callFrom == 'SelectAllConfs':
            self.molViewState.generateIds()
            
        #For multiple molecules with unequal number of conformers, it is hard to determine acceptable confIds
        self.wg_confId.options = self.molViewState.allConfIds
        self.wg_selectedMolsConfsView.value = ', '.join([str(x) for x in self.molViewState.idPaired])
        if callFrom == 'SelectMolAndConf':
            if self.onStart:
                self.Render3D()
                
    def SelectAllMols(self, confId = None, callFrom = 'nonWidget'):
        """ instantiates mol and conformer selection function of the MolViewState"""
        self.molViewState.selectAllMols()
        if confId is not None:
            self.molViewState.selectSingleConf(confId)
        if callFrom == 'widget':
            self.updateWidgets(callFrom = 'SelectAllMols')
            
    def SelectAllConfs(self, molId = None, callFrom = 'nonWidget'):
        """ instantiates mol and conformer selection function of the MolViewState"""
        if molId is not None:
            self.molViewState.selectSingleMol(molId)
        self.molViewState.selectAllConfs()
        if callFrom == 'widget':
            self.updateWidgets(callFrom = 'SelectAllConfs')
            
    def DeleteMolAndConf(self, molAndConfIds):
        """ instantiates mol and conformer selection function of the MolViewState"""
        self.molViewState.deleteIds(molAndConfIds = molAndConfIds)
        self.updateWidgets(callFrom = 'DeleteMolAndConf')
        
    def AddMolAndConf(self, molAndConfIds):
        """ instantiates mol and conformer selection function of the MolViewState"""
        self.molViewState.addIds(molAndConfIds = molAndConfIds)
        self.updateWidgets(callFrom = 'AddMolAndConf')
        
    def SelectMolAndConf(self, molAndConfIds):
        """ instantiates mol and conformer selection function of the MolViewState"""
        self.molViewState.selectSingleId(molAndConfIds = molAndConfIds)
        self.updateWidgets(callFrom = 'SelectMolAndConf')
        
    def SetMolData(self, 
                   ligandDict = None, 
                   protein = None, 
                   childKeysForMolRender = ['lig', 'emLig'], 
                   childKeyForConfSelection = 'emLig',
                   childKeyForDataPerConf = 'energy',
                   childKeyForPropSelection = 'emLig',
                   molAndConfIds = None, reinstantiated = False):
        """This function sets ligand dictionary, protein, and dict keys and initiates MolViewState class"""
        
        if isinstance(ligandDict, dict) is False:
            raise TypeError("ligandDict must be dict")
            
        if all(isinstance(key, str if PY3 else basestring) for key in ligandDict.keys()) is False:
            raise TypeError("keys of ligandDict must be str")
            
        self.childKeysForMolRender = childKeysForMolRender
        self.childKeyForConfSelection = childKeyForConfSelection
        
        additionalChildKeysForMolRender = set(self.childKeysForMolRender)
        additionalChildKeysForMolRender.discard(self.childKeyForConfSelection)
        self.additionalChildKeysForMolRender = list(additionalChildKeysForMolRender)
        
        self.childKeyForDataPerConf = childKeyForDataPerConf
        self.childKeyForPropSelection = childKeyForPropSelection
        
        if reinstantiated:
            self.molViewState = None
            
        self.molViewState = MolViewState(ligandDict, protein, 
                                         self.additionalChildKeysForMolRender, 
                                         childKeyForConfSelection, 
                                         childKeyForDataPerConf,
                                         childKeyForPropSelection)
        
        if self.molViewState.ligandDict is not None:    
            keys = list(self.molViewState.ligandDict.keys())
            self.onStartIds = ((keys[0], 0),) if molAndConfIds is None else molAndConfIds
            self.wg_molId.options = keys
            self.wg_molId.value = keys[0]
            self.SelectMolAndConf(molAndConfIds = self.onStartIds)
            
    def reInstantiateMolState(self):
        """This function resets mol data and mol viewer"""
        self.SetMolData(ligandDict = self.molViewState.ligandDict,
                        protein = self.molViewState.protein, 
                        childKeysForMolRender  = self.childKeysForMolRender ,
                        childKeyForConfSelection = self.childKeyForConfSelection,
                        childKeyForDataPerConf = self.childKeyForDataPerConf,
                        childKeyForPropSelection = self.childKeyForPropSelection,
                        molAndConfIds = None, reinstantiated = True)
        
    def UpdateConfId(self, childKeyForConfSelection = 'emLig'):
        """This function update the childKeys from which conf ids rendered"""
        self.childKeyForConfSelection = childKeyForConfSelection
        self.reInstantiateMolState()
        
    def ShowLigandProperty(self):
        """ Handles property in Render3D function """
        preCalcProp = self.molViewState.getPropPrecalculated
        if isinstance(preCalcProp, dict):
            self.wg_propSelect.options = sorted(list(preCalcProp.keys()))
            prop=self.wg_propSelect.value
            self.wg_propView.value = prop + ' : ' + str(preCalcProp[prop])
        elif preCalcProp == 'data not found':
            self.wg_propView.value = 'data not found'
        else:
            self.wg_propView.value = 'select single molecule'
            
    def ShowMinimizationEnergy(self):
        """ Handles energy minimization data in Render3D function """
        energy = self.molViewState.getMinimizationEnergy
        energyDataKey = self.molViewState.childKeyForDataPerConf
        confId = self.wg_confId.value
        if isinstance(energy, dict) and confId in energy:
            self.wg_energyView.value = energyDataKey +' : ' + str(energy[confId])
        elif isinstance(energy, dict) and confId not in energy:
            self.wg_energyView.value = 'data not found'
        elif energy is None:
            self.wg_energyView.value = 'select single molecule'
        else:
            self.wg_energyView.value = 'data not found'
            
    def GetAppliedStateRed(self):
        """this captures values from all the widgets except wg_modelSelect, wg_modelAdd and wg_zoomTo"""
        for wg in self.allWidgets:
            if wg in self.__dict__:
                self.wg_values[wg] = self.__dict__[wg].value
        self.wg_values['selectedMolsConfs'] = self.molViewState.idPaired
        
    def CompareWithPreviousStateRed(self):
        """this compares widget values previous and current state"""
        self.previous_state = self.wg_values
        self.GetAppliedState()
        self.new_wg_values = self.wg_values
        
    def AddLigandLabels(self):
        """ Handles ligand labeling (conf label and atom label) in Render3D function """
        if self.wg_confLabel.value:
            if self.molViewState.confLabel is not None:
                self.view.addLabel(self.molViewState.confLabel, {'backgroundColor':'gray', 
                                                                 'fontColor':'white',
                                                                 'showBackground':'true',
                                                                 'alignment':'bottomCenter'})
                
        if self.wg_atomLabel.value:
            if self.molViewState.atomLabel is not None:
                for (label,xyz) in self.molViewState.atomLabel:
                    self.view.addLabel(label, {'inFront' : 'false', 
                                               'fontSize' : '12',
                                               'fontColor':'gray',
                                               'showBackground':'false',
                                               'position' : {'x':xyz[0], 'y':xyz[1], 'z':xyz[2]}
                                              })
                    
    def AddProteinStyle(self, modelId):
        """ Handles protein drawing style in AddProteinWithStyle function """
        selectedModel = self.model_prop[modelId]
        style = selectedModel['style']
        color = selectedModel['color']
        hidden = selectedModel['hidden']
        colorFunc = 'colorscheme'
        
        if hidden is False:
            if style == 'surface':
                self.view.addSurface('SES', {'model': modelId});
                
            elif style == 'line':
                self.view.setStyle({'model': modelId},{'line':{colorFunc: color}});

            elif style == 'cartoon':
                self.view.setStyle({'model': modelId},{'cartoon':{colorFunc: color,
                                                                  'arrows': 'true'}
                                                      })
            else:
                self.view.setStyle({'model': modelId},{'cartoon':{colorFunc: color,
                                                                  'arrows': 'true', 
                                                                  'tubes' : 'true'}
                                                      })
        if hidden and style != 'surface':
            if style == 'line':
                self.view.setStyle({'model': modelId},{'line':{'hidden' : 'true'}});
                
            elif style == 'cartoon':
                self.view.setStyle({'model': modelId},{'cartoon':{'hidden' : 'true'}
                                                      })
            else:
                self.view.setStyle({'model': modelId},{'cartoon':{'hidden' : 'true'}
                                                      })
                
    def AddLigandStyle(self, modelId):
        """ Handles ligand drawing style and color in AddLigandWithStyle function """
        colorFunc = 'colorscheme'
        selectedModel = self.model_prop[modelId]
        style = selectedModel['style']
        color = selectedModel['color']
        hidden = selectedModel['hidden']
        
        if hidden is False:
            if style == 'surface':
                self.view.addSurface('SES', {'model':modelId});
            elif style == 'ballstick':
                self.view.setStyle({'model':modelId},{'stick':{'radius':'0.2', 
                                                               colorFunc: color},
                                                      'sphere':{'radius':'0.4', 
                                                                colorFunc: color}
                                                     });
            else:
                self.view.setStyle({'model': modelId},{style:{colorFunc: color}})
                
        if hidden and style != 'surface':
            if style == 'ballstick':
                self.view.setStyle({'model':modelId},{'stick':{'hidden' : 'true'},
                                                      'sphere':{'hidden' : 'true'}
                                                     });
            else:
                self.view.setStyle({'model': modelId},{style:{'hidden' : 'true'}})
                
    def RenderModel(self, callFrom = 'onStart'):
        """ this adds models and style for ligand and protein in viewer (called in Render3D function) """
        modelIds = list(self.model_prop.keys())
        for modelId in modelIds:
            selectedModel = self.model_prop[modelId]
            if callFrom == 'onStart':
                if selectedModel['molType'] == 'ligand':
                    self.view.addModel(selectedModel['geom'], 'sdf')
                else:
                    self.view.addModel(selectedModel['geom'], 'pdb')
                    
            if callFrom == 'onUpdate' and modelId > self.lastMoldelId:
                if selectedModel['molType'] == 'ligand':
                    self.view.addModel(selectedModel['geom'], 'sdf')
                else:
                    self.view.addModel(selectedModel['geom'], 'pdb')
                    
            if modelId > self.lastMoldelId:
                if selectedModel['molType'] == 'ligand':
                    self.AddLigandStyle(modelId)
                else:
                    self.AddProteinStyle(modelId)
            else:
                if selectedModel['reDraw']:
                    if selectedModel['molType'] == 'ligand':
                        self.AddLigandStyle(modelId)
                    else:
                        self.AddProteinStyle(modelId)
                        
    def CreateModelData(self, modelId, molAndConfId, geom, category, molType, callFrom = 'onStart'):
        """ this creates data before adding and styling models"""
        if molType == 'ligand':
            wgStyle = 'wg_' + category + 'Style'
            wgColor = 'wg_' + category + 'Color'
            style = self.__dict__[wgStyle].value if wgStyle in self.__dict__  else self.ligStyleScheme[category]
            color = self.__dict__[wgColor].value if wgColor in self.__dict__ else self.ligColorScheme[category]
            
        if molType == 'protein':
            wgStyle = 'wg_protStyle'
            wgColor = 'wg_protColor'
            style=self.__dict__[wgStyle].value if wgStyle in self.__dict__  else self.protStyleScheme
            color = self.__dict__[wgColor].value if wgColor in self.__dict__ else self.protColorScheme
            
        if callFrom == 'onStart':
            self.model_prop[modelId] = {}
            self.model_prop[modelId]['style'] = style
            self.model_prop[modelId]['color'] = color
            self.model_prop[modelId]['molType'] = molType
            self.model_prop[modelId]['categoryForDrawing'] = category
            self.model_prop[modelId]['molAndConfId'] = molAndConfId
            self.model_prop[modelId]['geom'] = geom
            self.model_prop[modelId]['reDraw'] = 'notApplicableOption'
            self.model_prop[modelId]['hidden'] = False
            self.model_prop[modelId]['newModel'] = True
            self.model_prop[modelId]['uniqueid'] = 'onStart'
            
        if callFrom == 'onUpdate':
            self.model_temp[modelId] = {}
            self.model_temp[modelId]['style'] = style
            self.model_temp[modelId]['color'] = color
            self.model_temp[modelId]['molType'] = molType
            self.model_temp[modelId]['categoryForDrawing'] = category
            self.model_temp[modelId]['molAndConfId'] = molAndConfId
            self.model_temp[modelId]['geom'] = geom
            self.model_temp[modelId]['hidden'] = False
            
    def GetModelFromMolViewState(self, callFrom = 'onStart'):
        """ this gererates data from molStateView """
        data = next(self.molViewState.selectedModels)
        for modelId in list(data.keys()):
            selectedModel = data[modelId]
            #print(selectedModel)
            molAndConfId = selectedModel['molAndConfId']
            geom = selectedModel['geom']
            category = selectedModel['categoryForDrawing']
            molType = selectedModel['molType']
            self.CreateModelData(modelId = modelId, 
                                 molAndConfId = molAndConfId, 
                                 geom = geom, 
                                 category = category, 
                                 molType = molType, 
                                 callFrom = callFrom)
            
    def ToggleHidden(self, category, preferredValue):
        """ this accepts values from visibility widgets"""
        # Decide hidden as per widget
        wgVisible = 'wg_' + category + 'Visible'
        if wgVisible in self.__dict__:
            if self.__dict__[wgVisible].value:
                hidden = False
            else:
                hidden = True
        else:
            hidden = preferredValue
        return hidden
        
    def ModelUpdate(self):
        """this is the main update function for mol viewer"""
        uniqueid = str(time.time()).replace('.','')
        oldData = self.model_prop
        newData = self.model_temp
        self.lastMoldelId = max(list(oldData.keys()))
        self.modelIdsForSurface = []
        molAndConfId = []
        categoryForDrawing = []
        modelIdsAll = []
        
        for oldModelId in list(oldData.keys()):
            oldSelectedData = oldData[oldModelId]
            molAndConfId.append(oldSelectedData['molAndConfId'])
            categoryForDrawing.append(oldSelectedData['categoryForDrawing'])
            modelIdsAll.append(oldModelId)
            
            if oldSelectedData['style'] == 'surface' and oldSelectedData['hidden'] is False:
                self.modelIdsForSurface.append(oldModelId)
            
        for newModelId in list(newData.keys()):
            newSelectedData = newData[newModelId]
            newMolAndConfId = newSelectedData['molAndConfId']
            newCategoryForDrawing = newSelectedData['categoryForDrawing']
            
            if newMolAndConfId in molAndConfId:
                indices = [i for i,x in enumerate(molAndConfId) if x == newMolAndConfId]
                ind = [i for i in indices if categoryForDrawing[i] == newCategoryForDrawing]
                
                if len(ind) == 1:
                    modelId = modelIdsAll[ind[0]]
                    oldSelectedData = oldData[modelId]
                    
                    # decide whether redraw
                    if oldSelectedData['hidden'] is False:
                        # was active in view previously
                        newStyle = newSelectedData['style']
                        newColor = newSelectedData['color']
                        oldStyle = oldSelectedData['style']
                        oldColor = oldSelectedData['color']
                        
                        tempData = False if newStyle == oldStyle and newColor == oldColor else True
                        oldSelectedData['reDraw'] = tempData
                        
                        # For surface ignore color change since color change doesn't work for surface
                        # If color change works for surface then we can remove reDraw... hidden is sufficient
                        if newStyle == 'surface' and oldStyle == 'surface': 
                            oldSelectedData['reDraw'] = False
                    else:
                        oldSelectedData['reDraw'] = True
                        
                    category = oldSelectedData['categoryForDrawing']
                    oldSelectedData['hidden'] = self.ToggleHidden(category = category, preferredValue = False)
                    oldSelectedData['uniqueid'] = uniqueid
                    oldSelectedData['newModel'] = False
                    oldSelectedData['style'] = newSelectedData['style']
                    oldSelectedData['color'] = newSelectedData['color']
                else:
                    # if category doesn't match with existing model list
                    #len(ind) == 0
                    newMoldelId = max(list(oldData.keys())) + 1
                    oldData[newMoldelId] = newSelectedData
                    oldData[newMoldelId]['reDraw'] = True
                    category = oldData[newMoldelId]['categoryForDrawing']
                    oldData[newMoldelId]['hidden'] = self.ToggleHidden(category=category,preferredValue =False)
                    oldData[newMoldelId]['newModel'] = True
                    oldData[newMoldelId]['uniqueid'] = uniqueid
            else:
                # if molAndConfId doesn't match with existing model list
                newMoldelId = max(list(oldData.keys())) + 1
                oldData[newMoldelId] = newSelectedData
                oldData[newMoldelId]['reDraw'] = True
                category = oldData[newMoldelId]['categoryForDrawing']
                oldData[newMoldelId]['hidden'] = self.ToggleHidden(category = category, preferredValue = False)
                oldData[newMoldelId]['newModel'] = True
                oldData[newMoldelId]['uniqueid'] = uniqueid
                
        for oldModelId in list(oldData.keys()):
            oldSelectedData = oldData[oldModelId]
            if oldSelectedData['uniqueid'] != uniqueid:
                oldSelectedData['uniqueid'] = uniqueid
                oldSelectedData['reDraw'] = True
                oldSelectedData['hidden'] = True
                oldSelectedData['newModel'] = False
                
    def RemoveAllSurfaceIfToggledToHidden(self):
        """ this decides whether surface will be removed or not """
        surfaceRemove = False
        oldData = self.model_prop
        for modelId in self.modelIdsForSurface:
            oldSelectedData = oldData[modelId]
            if oldSelectedData['style'] != 'surface':
                surfaceRemove = True
            if oldSelectedData['hidden']:
                surfaceRemove = True
                
        if surfaceRemove:
            self.view.removeAllSurfaces()
            
    def Render3D(self):
        """ This function updates the 3DMol.js viewer"""
        self.model_temp = {}
        self.view.removeAllLabels()
        
        if self.onStart:
            self.GetModelFromMolViewState(callFrom = 'onStart')
            self.RenderModel(callFrom = 'onStart')
        else:
            self.GetModelFromMolViewState(callFrom = 'onUpdate')
            self.ModelUpdate()
            self.RemoveAllSurfaceIfToggledToHidden()
            self.RenderModel(callFrom = 'onUpdate')
            
        self.view.setBackgroundColor(self.wg_background.value)
        
        # add label if required
        if self.labelPanel:
            self.AddLigandLabels()
            
        # zoomTo does not work well for surface and label... so, zoomTo should not be default settings
        if self.onStart:
            self.view.zoomTo()
            self.onStart = False
            
        if self.propertyPanel:
            self.ShowLigandProperty()
            
        if self.emPanel:
            self.ShowMinimizationEnergy()
            
        display(self.view.update())
        