<%@ taglib uri="http://java.sun.com/jsp/jstl/core" prefix="c" %>
<script type="text/javascript">
<!--
function setSmiles()
{
    var cbv;
    if (form1.cb.checked) {
        cbv="true";
    }
    else {
        cbv="false";
    }
    form1.image.src="mol-renderer/test.png?h=400&w=600&flags="+form1.flags.value+form1.other.value+"&cb="+cbv+"&smiles="+encodeURIComponent(form1.smiles.value)+"&random="+Math.random();
    form1.smallimage.src="mol-renderer/test.png?h=140&w=210&flags="+form1.flags.value+form1.other.value+"&cb="+cbv+"&smiles="+encodeURIComponent(form1.smiles.value)+"&random="+Math.random();
    molref.href="mol-renderer/test.mol?flags="+form1.flags.value+form1.other.value+"&smiles="+encodeURIComponent(form1.smiles.value)+"&random="+Math.random();
    smiref.href="mol-renderer/test.smi?flags="+form1.flags.value+form1.other.value+"&smiles="+encodeURIComponent(form1.smiles.value)+"&random="+Math.random();
}
//-->
</script>
<html>
    <head>
        <title>Interactive SMILES Visualizer</title>
    </head>
    <body>
        <h1>Interactive SMILES Visualizer</h1>
        <c:url var="url_mol" value="mol-renderer/test.mol">
            <c:param name="smiles" value="C12(C)CCC1(C)OCC(=C3C(CCO)CCC(CCC=C)C(CCC=C)CCC3(CCO))CO2"/>
        </c:url>
        <c:url var="url_smi" value="mol-renderer/test.smi">
            <c:param name="smiles" value="C12(C)CCC1(C)OCC(=C3C(CCO)CCC(CCC=C)C(CCC=C)CCC3(CCO))CO2"/>
        </c:url>
        <c:url var="url" value="mol-renderer/test.png">
            <c:param name="smiles" value="C12(C)CCC1(C)OCC(=C3C(CCO)CCC(CCC=C)C(CCC=C)CCC3(CCO))CO2"/>
            <c:param name="w" value="600"/>
            <c:param name="h" value="400"/>
        </c:url>
        <c:url var="smallurl" value="mol-renderer/test.png">
            <c:param name="smiles" value="C12(C)CCC1(C)OCC(=C3C(CCO)CCC(CCC=C)C(CCC=C)CCC3(CCO))CO2"/>
            <c:param name="w" value="210"/>
            <c:param name="h" value="140"/>
        </c:url>
        <form name="form1">
            <input type="text"
                   name="smiles"
                   onKeyUp="javascript:setSmiles()"
                   value='<c:out value="C12(C)CCC1(C)OCC(=C3C(CCO)CCC(CCC=C)C(CCC=C)CCC3(CCO))CO2"/>'
		   size="50"/> <br />

                Shortcut Flags:
                <select name="flags" onChange="javascript:setSmiles()">
                    <option value="">None</option>
                    <option value="aminoacids" selected="selected">aminoacids</option>
                    <option value="aminoacids,extended">aminoacids,extended</option>
                    <option value="aminoacids,extended,non-standard">aminoacids,extended,non-standard</option>
                    <option value="aminoacids,extended,non-standard,catch-all">aminoacids,extended,non-standard,catch-all</option>
                </select> <br/>
                Color-blind support: <input name="cb" type="checkbox" onclick="javascript:setSmiles()"/> <br/>
                Other Flags:
                <select name="other" onChange="javascript:setSmiles()">
                    <option value="" selected="selected">None</option>
                    <option value=",chains">Abbreviate Chains</option>
                </select>
            <p/>
            <a name="molref" href='<c:out value="${url_mol}"/>'>
            <img name="image" src='<c:out value="${url}"/>'/>
            </a>
            <img name="smallimage" src='<c:out value="${smallurl}"/>'/>
            <br/>
            <a name="smiref" href='<c:out value="${url_smi}"/>'>Converted SMILES</a>
        </form>
    </body>
</html>
