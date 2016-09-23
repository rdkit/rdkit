from mod_python import apache
import sys, os, tempfile, urllib
from WebUtils import General

General._version = "1.0.0"


def page(req, smiles='', width=300, height=300, highlight='[]', numbers=0, **kwargs):
  req.content_type = 'text/html'
  page = General.ConstructHtmlHeader(title='RD Depict')

  if smiles:
    uSmiles = urllib.quote_plus(smiles)
    uHighlight = urllib.quote_plus(highlight)
    url = '/RDExtras/MolImage.py/svg?smiles=%s&height=%s&width=%s&highlight=%s&numbers=%s' % (
      uSmiles, height, width, uHighlight, numbers)
    imgUrl = '/RDExtras/MolImage.py/gif?smiles=%s&height=%s&width=%s&useCactvs=1' % (uSmiles,
                                                                                     height, width)

    page += """<center>
    <table>
    <tr>
      <td><embed src="%s" Name="SVGEmbed" width=%s height=%s type="image/svg-xml">
      <td><img src="%s" width=%s height=%s></td>
    </tr>
    </table>
    </center>""" % (url, width, height, imgUrl, width, height)
    page += '<center><b>SMILES:</b>%s</center>' % smiles
  if not smiles:
    smiles = '""'
  if not numbers or numbers == 'off':
    checked = ""
  else:
    checked = "checked"
  page += """<center>
  <form action="/RDExtras/MolDepict.py/page">
  <input type=text name="smiles" value=%s size=80>
  <input type=submit Value="Depict">
  <br><input type=checkbox name="numbers" %s> Numbers
  <input type=hidden name="width" value=%s>
  <input type=hidden name="height" value=%s>
  <input type=hidden name="highlight" value=%s>
  </center>
  """ % (smiles, checked, width, height, highlight)

  page += General.ConstructHtmlFooter(includeRestart=0, logoutText='')

  return page
