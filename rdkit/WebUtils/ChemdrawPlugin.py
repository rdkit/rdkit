#
#Copyright (C) 2002 Greg Landrum and Rational Discovery LLC
# All rights are reserved.
#

""" Code for dealing with the ChemDraw Plugin on web pages

"""

cdpHeader = """
<!-- grab the chemdraw javascript file -->
<script language="JavaScript" src="/chemdraw/chemdraw.js"></script>
<script> cd_includeWrapperFile("/chemdraw/"); </script>
<script language="JavaScript">
<!--
        function GetDataFromCDP(element,pluginName,format)
	{
	  element.value = cd_getData(pluginName,format,false)
	}
        function SetDataInCDP(val,pluginName,format)
	{
	  cd_putData(pluginName,format,val)
	}
<!-- -->
</script>
"""
pluginID = 1

def IncrID():
  """ increments the global ID number
  used to ensure no two plugins on the same page have the same ID

  """
  global pluginID
  pluginID = pluginID + 1
  if pluginID == 5000:
    pluginID = 1
    
cdpEmbedCode= """
<script language="JavaScript">
cd_insertObject("chemical/x-cdx",300,300,"%s","%s")
</script>
"""

cdpEmbedWithDataCode= """
<script language="JavaScript">
cd_insertObject("chemical/x-cdx",300,300,"%s","%s",false,true,dataurl="data:%s,%s")
</script>
"""

cdpAppletCode= """
"""

def EmbedPlugin(molData=None,id=None,template="",fmt="chemical/smiles"):
  """ returns the text to embed a plugin allowing editing

   **Arguments**

     - molData: (optional) the molecular data

     - id: (optional) the id to use

     - template: (optional) URL to a file to be used as a template

     - fmt: (optional) MIME code for the data format

  """
  if id is None:
    global pluginID
    id = pluginID
    IncrID()    
  if molData is not None:
    if fmt.upper().index('SMILES') != -1:
      molData = molData.replace('\\','\\\\')
    res = cdpEmbedWithDataCode%(str(id),template,fmt,molData)
  else:
    res = cdpEmbedCode%(str(id),template)
  return res

cdpEmbedViewerCode= """
<script language="JavaScript">
cd_insertObject("chemical/x-cdx",%d,%d,"cdp%d","%s",true,true,dataurl="data:%s,%s")
</script>
"""

def EmbedPluginViewer(molData,size=(200,200),id=None,template="",fmt="chemical/smiles"):
  """ returns the text to embed a plugin only for viewing

   **Arguments**

     - molData: the molecular data

     - size: a 2-tuple with the size to be used for the plugin

     - id: (optional) the id to use

     - template: (optional) URL to a file to be used as a template

     - fmt: (optional) MIME code for the data format

  """
  if id is None:
    global pluginID
    id = pluginID
    IncrID()
  if fmt.upper().index('SMILES') != -1:
    molData = molData.replace('\\','\\\\')
  res = cdpEmbedViewerCode%(size[0],size[1],int(id),template,fmt,molData)
  return res





