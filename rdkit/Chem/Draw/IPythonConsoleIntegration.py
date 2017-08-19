
"""RDKit Conformer Browser
Derived by Malitha Humayun Kabir from Greg Landrum's code as a part of GSoC 2017
Project : RDKit - 3Dmol.js integration
Mentors: Paul Czodrowski and Greg Landrum
Acknowledgement: Peter Gedeck reviewed and provided advice on restructuring codes and wrote MolViewState class.
Date: 18th August 2017
Email# malitha12345@gmail.com
"""

import py3Dmol
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from ipywidgets import Layout, Label, Button, Box, HBox, VBox
from ipywidgets import Dropdown, SelectMultiple, IntSlider, HTML, Checkbox, Button, Text
from IPython.display import display
from IPython.display import HTML as scriptHTML
import time
import sys


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
                    res = '3DMol.js loaded'
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
def CreateMolDict(mol):
    """This function checks whether the object type fits the requirements for rendering.
    If the oject doesn't have necessary attributes, it takes action to include that.
    """
    if isinstance(mol, dict):
        moldict = dict(mol)
    elif type(mol) is tuple:
        moldict = {mol[0]: mol[1]}
    elif hasattr(mol, '__iter__') is False:
        moldict = {'0': mol}
    elif type(mol) is list:
        mol_keys=[str(x) for x in range(len(mol))]
        mol_2=[(mol_keys[i],mol[i]) for i in range(len(mol))]
        moldict = dict(mol_2)
    else:
        moldict = None
    return moldict
    
    
class MolViewState(object):
    def __init__(self, uid, molecules, protein, calculateProperties = False):
        """ Molecules is dictionary of molecules """
        if calculateProperties:
            for i in molecules:
                mol = molecules[i]
                for prop_name in PROP_RDKIT:
                    calculator = Descriptors.__dict__[prop_name]
                    mol.SetProp(prop_name, str(calculator(mol)))
                    
        self.uid = uid
        self.moldict = molecules    
        self.protein = protein
        # These should have reasonable initial values
        self.rdkit_mol_select = set()
        self.rdkit_conf_select = set()
        self.allConfIDs = []
        
    def selectMolecules(self, selectAllMols, selectMultiMols, selectMol):
        """ Select either all moleculs or add selectMol or show only selectMol """
        if selectAllMols:
            self.rdkit_mol_select = set(self.moldict)
        elif selectMultiMols:
            self.rdkit_mol_select.add(selectMol)
        else:
            self.rdkit_mol_select = {selectMol}
            
    def selectConformations(self, selectAllConfs, selectMultiConfs, selectConf):
        """ For all selected molecules, select either all conformations or add selectConf or show only selectConf """
        for mol in self.selectedMolecules:
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
                
    def selectedIds(self, molAndConfIds = None):
        if molAndConfIds is not None:
            self.idPaired = molAndConfIds
            molIds = []
            confIds = []
            for (molId, confId) in molAndConfIds:
                molIds.append(molId)
                confIds.append(confId)
            self.rdkit_mol_select = set(molIds)
            self.rdkit_conf_select = set(confIds)
        else:
            Ids = []
            for molId in self.rdkit_mol_select:
                for confId in self.rdkit_conf_select:
                    Ids.append((molId,confId))
                    self.idPaired = Ids
            
    @property
    def selectedModels(self):
        """ Iterator over all selected models (molecules/conformations) """
        for (molId, confId) in self.idPaired:
            yield Chem.MolToMolBlock(self.moldict[molId], confId=confId)
            
    @property
    def selectedMolNames(self):
        """ Return the names of all selected molecules """
        return self.rdkit_mol_select
        
    @property
    def selectedConfIds(self):
        """ Return the names of all selected molecules """
        return self.rdkit_conf_select
        
    @property
    def selectedMolecules(self):
        """ Return the selected molecules """
        return [self.moldict[name] for name in self.selectedMolNames]
        
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
            mol = self.moldict[list(self.selectedMolNames)[0]]
            if len(mol.GetPropNames()) == 0:
                return 'Not found'
            else:
                return {prop: mol.GetProp(prop) for prop in mol.GetPropNames()}
                
    
    
class MolView3D(object):
    """This function initiates required widgets and 3Dmol.js viewer
    uid : an unique id for 3DMol.js viewer, each of the ipywidgets, and all the global objects
    moldict : it is a molecule dictionary where is keys are string
    protein : it is rdkit standard mol object
    molStyle : representation type for moldict or small molecule
    pStyle : protein representation style
    ** the rest of the arguments are panels (user may or may not need those)
    """
    def __init__(self, 
                 moldict = None, protein = None,
                 molStyle = 'stick', pStyle='cartoon',
                 molPanel = 'full',
                 stylePanel = False, 
                 labelPanel = False,
                 propPanel = False, calculateProperties = False
                 ):
        
        if molPanel not in ('full', 'minimal'):
            raise KeyError('molPanel can be full or minimal')
            
        uid = str(time.time()).replace('.','')
        
        self.uid = uid
        self.molStyle = molStyle
        self.pStyle = pStyle
        self.onStart = True
        self.molPanel = molPanel # molPanel can be full or minimal
        self.propPanel = propPanel
        self.stylePanel = stylePanel
        self.labelPanel = labelPanel
        
        wgData = globals().setdefault('wgData', {})
        wgData[uid] = {}
        
        # Right hand panel (widgets)
        
        wgBox = []
        
        itemLayout=Layout(display='flex', flex_flow='row', justify_content='space-between')
        
        wgData[uid]['molId'] = Dropdown(description='', options=['select'],value='select')
        wgData[uid]['confId'] = Dropdown(description='', options=['select'], value='select')
        
        if molPanel == 'full':
            
            wgData[uid]['selectedMolsConfsView'] = HTML(description='', value='')
            
            wgData[uid]['selectMultiMols'] = Checkbox(description='selectMultiMols', value=False)
            wgData[uid]['selectAllMols'] = Checkbox(description='selectAllMols', value=False)
            
            wgData[uid]['selectMultiConfs'] = Checkbox(description='selectMultiConfs', value=False)
            wgData[uid]['selectAllConfs'] = Checkbox(description='selectAllConfs', value=False)
            
            
            wgBox.append(Box([Label(value='ids'),wgData[uid]['selectedMolsConfsView']], layout=itemLayout))
            
            wgBox.append(Box([Label(value='molId'),wgData[uid]['molId']], layout=itemLayout))
            
            box_for_multimol_allmol = HBox([wgData[uid]['selectMultiMols'], wgData[uid]['selectAllMols']])
            cbMolSelect=Box([Label(value=''), box_for_multimol_allmol], layout=itemLayout)
            
            wgBox.append(cbMolSelect)
            
            wgBox.append(Box([Label(value='confId'),wgData[uid]['confId']], layout=itemLayout))
            
            box_for_multiconf_allconf = HBox([wgData[uid]['selectMultiConfs'], wgData[uid]['selectAllConfs']])
            cbConfSelect=Box([Label(value=''), box_for_multiconf_allconf], layout=itemLayout)
            
            wgBox.append(cbConfSelect)
            
            widgetsFor3DView = ['molId', 'selectMultiMols', 'selectAllMols', 
                                'confId', 'selectMultiConfs', 'selectAllConfs']
        else:
            wgBox.append(Box([Label(value='molId'),wgData[uid]['molId']], layout=itemLayout))
            wgBox.append(Box([Label(value='confId'),wgData[uid]['confId']], layout=itemLayout))
            widgetsFor3DView = ['molId', 'confId']
            
        if protein is not None:
            wgData[uid]['proteinVisible'] = Checkbox(description='proteinVisible', value=True)
            wgBox.append(Box([Label(value=''),wgData[uid]['proteinVisible']], layout=itemLayout))
            widgetsFor3DView.append('proteinVisible')
            
        # prop
        if propPanel:
            wgData[uid]['prop_view'] = HTML(description='', value='initializing...')
            wgData[uid]['prop_wg'] = Dropdown(description='',options=['select'], value='select')
            wgBox.append(Box([Label(value='prop'),wgData[uid]['prop_view']], layout=itemLayout))
            wgBox.append(Box([Label(value='prop'),wgData[uid]['prop_wg']], layout=itemLayout))
            widgetsFor3DView.append('prop_wg')
        
        # stylePanel
        if stylePanel:
            if moldict is not None:
                wgData[uid]['molStyle_wg'] = Dropdown(description='', options=DRAWING_LIGAND_3D, 
                                                      value=molStyle)
                wgData[uid]['colorScheme'] = Dropdown(description='', options=LIGAND_COLOR_SCHEME_3D,
                                                      value='default')
                wgBox.append(Box([Label(value='molStyle'),wgData[uid]['molStyle_wg']], layout=itemLayout))
                widgetsFor3DView.append('molStyle_wg')
                wgBox.append(Box([Label(value='ligand color'),wgData[uid]['colorScheme']], layout=itemLayout))
                widgetsFor3DView.append('colorScheme')
            if protein is not None:
                wgData[uid]['pStyle_wg'] = Dropdown(description='', options=DRAWING_PROTEIN_3D, value=pStyle)
                wgBox.append(Box([Label(value='pStyle'),wgData[uid]['pStyle_wg']], layout=itemLayout))
                wgData[uid]['pStyle_helicesAsTubes_wg'] = Checkbox(description='helicesAsTubes', value=False)
                wgBox.append(Box([Label(value=''),wgData[uid]['pStyle_helicesAsTubes_wg']], layout=itemLayout))
                widgetsFor3DView.extend(['pStyle_wg', 'pStyle_helicesAsTubes_wg'])
        
        # labels
        if labelPanel:
            wgData[uid]['confLabel'] = Checkbox(description='confLabel', value=False)
            wgData[uid]['atomLabel'] = Checkbox(description='atomLabel', value=False)
            box_for_labeling = HBox([wgData[uid]['confLabel'],wgData[uid]['atomLabel']])
            cbLabel=Box([Label(value=''), box_for_labeling], layout=itemLayout)
            wgBox.append(cbLabel)
            widgetsFor3DView.extend(['confLabel', 'atomLabel'])
        
        
        # background
        wgData[uid]['background'] = Dropdown(description='', options=BGCOLORS_3D, value=BGCOLORS_3D[0])
        wgBox.append(Box([Label(value='background'),wgData[uid]['background']], layout=itemLayout))
        widgetsFor3DView.append('background')
        
        # buttons
        wgData[uid]['zoomTo'] = Button(description="zoomTo", button_style='success')
        
        box_for_zoomTo = Box([Label(value=''),HBox([wgData[uid]['zoomTo']])],layout=itemLayout)
        wgBox.append(box_for_zoomTo)
        
        # Observe and update dict with model_id
        wgData[uid]['zoomTo'].on_click(self.handle_zoomTo_button)
        
        wgDataUid = wgData[uid]
        for i in widgetsFor3DView:
            selectedWidget = wgDataUid[i]
            selectedWidget.observe(self.handle_change, names='value')
            
            
        # right panel
        wgRightBox=VBox(wgBox, layout=Layout(border='solid', width='45%'))
        
        # left panel (container holding table for viewer)
        size = (435, 485)
        viewerLayout=Layout(width=str(size[0]+5)+'px', height=str(size[1]+5)+'px', border='solid')
        wgLeftWidget = HTML('''<table><tr><td id="%s"></td></tr></table>'''%uid, layout=viewerLayout)
        
        # Combining left and right panel
        viewer = HBox([wgLeftWidget, wgRightBox])
        
        # displaying everything
        display(viewer)
        
        # inserting 3DMol.js viewer in existing container (table)
        wgData[uid]['view'] = py3Dmol.view(width=size[0],height=size[1])
        wgData[uid]['view'].setBackgroundColor('0x000000')
        wgData[uid]['view'].zoomTo()
        
        globals()['wgData'] = wgData
        
        display(globals()['wgData'][uid]['view'].insert(uid))
        
        if moldict is not None:
            self.SetMolData(moldict, protein, calculateProperties)
            
    
    def SetMolData(self, moldict = None, protein = None, calculateProperties = False, molAndConfIds = None):
        """This function update molecules."""
        
        if isinstance(moldict, dict) is False:
            raise TypeError("moldict must be dict")
        
        if all(isinstance(key, str if PY3 else basestring) for key in moldict.keys()) is False:
            raise TypeError("keys of moldict must be str")
        
        
        self.molViewState = MolViewState(self.uid, moldict, protein, calculateProperties)
        
        keys = list(self.molViewState.moldict.keys())
        
        wgData = globals()['wgData'][self.uid]
        
        wgData['molId'].options = keys
        wgData['molId'].value = keys[0]
        
        wgData['confId'].options = list(range(9))
        wgData['confId'].value=0
        
        if self.molPanel == 'full':
            ids = ((keys[0], 0),)
            wgData['selectedMolsConfsView'].value = ', '.join([str(x) for x in ids])
            
        self.render3D()
        
    
    def handle_change(self, change):
        """This function handles all the interactive widgets except buttons and 3Dmol.js viewer"""
        self.SelectMolAndConf()
        
    
    def handle_zoomTo_button(self, b):
        """This function handles zoomTo button"""
        
        wgData = globals()['wgData'][self.uid]
        wgData['view'].zoomTo()
        display(wgData['view'].update())
        
    
    def SelectMolAndConf(self, molAndConfIds = None):
        """ Handles molecules and conformers in update3D function """
        
        wgData = globals()['wgData'][self.uid]
          
        if molAndConfIds is None and self.molPanel == 'full':
            
            self.molViewState.selectMolecules(wgData['selectAllMols'].value, 
                                              wgData['selectMultiMols'].value,
                                              wgData['molId'].value)
            
            self.molViewState.selectConformations(wgData['selectAllConfs'].value, 
                                                  wgData['selectMultiConfs'].value,
                                                  wgData['confId'].value)
            
            self.molViewState.selectedIds()
            
            wgData['confId'].options = self.molViewState.allConfIds
            
            wgData['selectedMolsConfsView'].value = ', '.join([str(x) for x in self.molViewState.idPaired])
            
        elif molAndConfIds is None and self.molPanel == 'minimal':
            self.molViewState.selectMolecules(False, False, wgData['molId'].value)
            self.molViewState.selectConformations(False, False, wgData['confId'].value)
            wgData['confId'].options = self.molViewState.allConfIds
        
        else:
            self.molViewState.selectedIds(molAndConfIds)
            
            wgData['confId'].options = self.molViewState.allConfIds
            
            sMolNames = self.molViewState.selectedMolNames
            sConfIds = self.molViewState.selectedConfIds
            
            wgData['selectedMolsConfsView'].value = ', '.join([str(x) for x in self.molViewState.idPaired])
        
        self.render3D()
        
    def handleProperty(self):
        """ Handles property in update3D function """
        wgData = globals()['wgData'][self.uid]
        preCalcProp = self.molViewState.getPropPrecalculated
        if 'prop_wg' in wgData:
            if isinstance(preCalcProp, dict):
                wgData['prop_wg'].options = sorted(list(preCalcProp.keys()))
                prop=wgData['prop_wg'].value
                wgData['prop_view'].value = prop + ' : ' + str(preCalcProp[prop])
            elif preCalcProp == 'Not found':
                wgData['prop_view'].value = 'No precalculated property found!'
            else:
                wgData['prop_view'].value = 'Single molecule selection required!'
                
    
    def handleLigandStyle(self, view):
        """ Handles ligand drawing style and color in update3D function """
        
        wgData = globals()['wgData'][self.uid]
        molStyle = wgData['molStyle_wg'].value if self.stylePanel else self.molStyle 
        color = wgData['colorScheme'].value if self.stylePanel else 'default'
        
        if molStyle == 'surface':
            view.addSurface('SES', {});
        elif molStyle == 'ballstick':
            view.setStyle({},{'stick':{'radius':'0.2','colorscheme': color},
                              'sphere':{'radius':'0.4', 'colorscheme': color}
                             }
                         );
        else:
            view.setStyle({},{molStyle:{'colorscheme': color}})
            
    
    def handleLigandLabeling(self,view):
        """ Handles ligand labeling (conf label and atom label) in update3D function """
        
        wgData = globals()['wgData'][self.uid]
        sMolNames = self.molViewState.selectedMolNames
        sConfIds = self.molViewState.selectedConfIds
        
        if not len(sMolNames)==1 and len(sConfIds)==1:
            return None
        else:
            if 'confLabel' in wgData:
                if wgData['confLabel'].value:
                    confLabel = list(sMolNames)[0] + ':' + str(wgData['confId'].value)
                    view.addLabel(confLabel, {'backgroundColor':'gray', 'fontColor':'white',
                                              'showBackground':'true', 'alignment':'bottomCenter'})
                
            if 'atomLabel' in wgData:
                
                mol = self.molViewState.moldict[list(sMolNames)[0]]
                confId = wgData['confId'].value
                sconf = mol.GetConformer(confId)
                xyz = sconf.GetPositions()
                
                if wgData['atomLabel'].value:
                    for i in list(range(sconf.GetNumAtoms())):
                        atomLabel = '{}{}'.format(mol.GetAtomWithIdx(i).GetSymbol(), i)
                        view.addLabel(atomLabel, {'inFront' : 'false', 
                                                  'fontSize' : '12',
                                                  'fontColor':'gray',
                                                  'showBackground':'false',
                                                  'position' : {'x':xyz[i][0], 'y':xyz[i][1], 'z':xyz[i][2]}
                                                 })
                        
    
    def handleProteinStyle(self, view, pModelNum):
        """ Handles protein drawing style in update3D function """
        
        wgData = globals()['wgData'][self.uid]
        
        # 'pStyle_wg' and 'pStyle_helicesAsTubes_wg' rendered as pair
        pStyle = wgData['pStyle_wg'].value if self.stylePanel else self.pStyle
        
        if pStyle == 'surface':
            view.addSurface('SES', {'model':pModelNum});
        elif pStyle == 'line':
            view.setStyle({'model':pModelNum},{'line':{}});
        elif pStyle == 'cartoon' and 'pStyle_helicesAsTubes_wg' not in wgData:
            view.setStyle({'model':pModelNum},{'cartoon':{'color': 'spectrum', 'arrows': 'true'}})
        else:
            if wgData['pStyle_helicesAsTubes_wg'].value:
                view.setStyle({'model':pModelNum},{'cartoon':{'color': 'spectrum',
                                                              'arrows': 'true', 'tubes' : 'true'}})
            else:
                view.setStyle({'model':pModelNum},{'cartoon':{'color': 'spectrum', 'arrows': 'true'}})
                
    
    
    def render3D(self):
        """ This function invoked whenever user interacts with widgets.
        It runs first time through handle_button() when the start button clicked.
        """
        wgData = globals()['wgData'][self.uid]
        
        view = wgData['view']
        view.removeAllModels()
        view.removeAllSurfaces()
        view.removeAllLabels()
        
        view.setBackgroundColor(wgData['background'].value)
        
        # add models
        for model in self.molViewState.selectedModels:
            view.addModel(model, 'sdf')
            
        # add style
        self.handleLigandStyle(view)
        
        # add label if required
        if self.labelPanel:
            self.handleLigandLabeling(view)
            
        # Add model (protein)
        if self.molViewState.protein is not None:
            # Show protein if visibitity is true
            if 'proteinVisible' in wgData:
                if wgData['proteinVisible'].value:
                    pModelNum = len(list(self.molViewState.selectedModels))
                    pdb = Chem.MolToPDBBlock(self.molViewState.protein)
                    view.addModel(pdb,'pdb')
                    self.handleProteinStyle(view, pModelNum)
                    
        # zoomTo does not work well for surface and label... so, zoomTo should not be default settings
        if self.onStart:
            view.zoomTo()
            self.onStart = False
            
        if self.propPanel:
            self.handleProperty()
            
        display(view.update())
        
        
        
        