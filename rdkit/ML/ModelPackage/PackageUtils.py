#
# Copyright (C) 2003 Rational Discovery LLC
# All rights are reserved.
#

# from elementtree.ElementTree import ElementTree, Element, SubElement
import time
from xml.etree.ElementTree import ElementTree, Element, SubElement


def _ConvertModelPerformance(perf, modelPerf):
    if len(modelPerf) > 3:
        confMat = modelPerf[3]
        accum = 0
        for row in confMat:
            for entry in row:
                accum += entry
        accum = str(accum)
    else:
        confMat = None
        accum = 'N/A'

    if len(modelPerf) > 4:
        elem = SubElement(perf, "ScreenThreshold")
        elem.text = str(modelPerf[4])
    elem = SubElement(perf, "NumScreened")
    elem.text = accum
    if len(modelPerf) > 4:
        elem = SubElement(perf, "NumSkipped")
        elem.text = str(modelPerf[6])
    elem = SubElement(perf, "Accuracy")
    elem.text = str(modelPerf[0])
    elem = SubElement(perf, "AvgCorrectConf")
    elem.text = str(modelPerf[1])
    elem = SubElement(perf, "AvgIncorrectConf")
    elem.text = str(modelPerf[2])
    if len(modelPerf) > 4:
        elem = SubElement(perf, "AvgSkipConf")
        elem.text = str(modelPerf[5])
    if confMat:
        elem = SubElement(perf, "ConfusionMatrix")
        elem.text = str(confMat)


def PackageToXml(pkg, summary="N/A", trainingDataId='N/A', dataPerformance=[],
                 recommendedThreshold=None, classDescriptions=None, modelType=None,
                 modelOrganism=None):
    """ generates XML for a package that follows the RD_Model.dtd

    If provided, dataPerformance should be a sequence of 2-tuples:
      ( note, performance )
    where performance is of the form:
      ( accuracy, avgCorrectConf, avgIncorrectConf, confusionMatrix, thresh, avgSkipConf, nSkipped )
      the last four elements are optional

    """
    head = Element("RDModelInfo")
    name = SubElement(head, "ModelName")
    notes = pkg.GetNotes()
    if not notes:
        notes = "Unnamed model"
    name.text = notes
    summ = SubElement(head, "ModelSummary")
    summ.text = summary
    calc = pkg.GetCalculator()
    descrs = SubElement(head, "ModelDescriptors")
    for name, summary, func in zip(calc.GetDescriptorNames(), calc.GetDescriptorSummaries(),
                                   calc.GetDescriptorFuncs()):
        descr = SubElement(descrs, "Descriptor")
        elem = SubElement(descr, "DescriptorName")
        elem.text = name
        elem = SubElement(descr, "DescriptorDetail")
        elem.text = summary
        if hasattr(func, 'version'):
            vers = SubElement(descr, "DescriptorVersion")
            major, minor, patch = func.version.split('.')
            elem = SubElement(vers, "VersionMajor")
            elem.text = major
            elem = SubElement(vers, "VersionMinor")
            elem.text = minor
            elem = SubElement(vers, "VersionPatch")
            elem.text = patch

    elem = SubElement(head, "TrainingDataId")
    elem.text = trainingDataId

    for description, perfData in dataPerformance:
        dataNode = SubElement(head, "ValidationData")
        note = SubElement(dataNode, 'ScreenNote')
        note.text = description
        perf = SubElement(dataNode, "PerformanceData")
        _ConvertModelPerformance(perf, perfData)

    if recommendedThreshold:
        elem = SubElement(head, "RecommendedThreshold")
        elem.text = str(recommendedThreshold)

    if classDescriptions:
        elem = SubElement(head, "ClassDescriptions")
        for val, text in classDescriptions:
            descr = SubElement(elem, 'ClassDescription')
            valElem = SubElement(descr, 'ClassVal')
            valElem.text = str(val)
            valText = SubElement(descr, 'ClassText')
            valText.text = str(text)

    if modelType:
        elem = SubElement(head, "ModelType")
        elem.text = modelType
    if modelOrganism:
        elem = SubElement(head, "ModelOrganism")
        elem.text = modelOrganism

    hist = SubElement(head, "ModelHistory")
    revision = SubElement(hist, "Revision")
    tm = time.localtime()
    date = SubElement(revision, "RevisionDate")
    elem = SubElement(date, "Year")
    elem.text = str(tm[0])
    elem = SubElement(date, "Month")
    elem.text = str(tm[1])
    elem = SubElement(date, "Day")
    elem.text = str(tm[2])
    note = SubElement(revision, "RevisionNote")
    note.text = "Created"
    return ElementTree(head)


if __name__ == '__main__':  # pragma: nocover
    import sys
    import pickle
    from io import StringIO
    pkg = pickle.load(open(sys.argv[1], 'rb'))
    perf = (.80, .95, .70, [[4, 1], [1, 4]])
    tree = PackageToXml(pkg, dataPerformance=[('training data performance', perf)])
    io = StringIO()
    tree.write(io)
    txt = io.getvalue()
    header = """<?xml version="1.0"?>
<!DOCTYPE RDModelInfo PUBLIC "-//RD//DTD RDModelInfo //EN" "RD_Model.dtd">
"""
    print(header)
    print(txt.replace('><', '>\n<'))
