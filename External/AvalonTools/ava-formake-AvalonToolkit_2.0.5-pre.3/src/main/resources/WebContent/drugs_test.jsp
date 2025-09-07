<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN" "http://www.w3.org/TR/html4/loose.dtd">
<%@ taglib prefix="c" uri="http://java.sun.com/jsp/jstl/core" %>
<%@ taglib prefix="x" uri="http://java.sun.com/jsp/jstl/xml" %>
<html><body>

<c:set var="renderer" value="../depicter/mol-renderer" />

<c:set var="size" value="100" />
<c:if test="${!empty param.size}">
    <c:set var="size" value="${param.size}" />
</c:if>
<c:import var="source" url="SampleDrugs.xml" />
<x:parse var="smiles_list" xml="${source}" />
<x:forEach select="$smiles_list//smiles">
<c:set var="irow" value="${irow+1}" />
<c:set var="x" value="${x+size+2}" />
    <c:set var="smiles"><x:out select="./@value" /></c:set>
    <c:set var="name"><x:out select="./@name" /></c:set>
    <c:url value="${renderer}/test.gif" var="url">
      <c:param name="smiles" value="${smiles}" />
      <c:param name="footer" value="${irow}: ${name}" />
      <c:param name="w" value="${size}" />
      <c:param name="h" value="${size}" />
    </c:url>
    
    <div style="float: left; border: 1px solid;">
        <img width=${size} height=${size} src="${url}"/>
    </div>
</x:forEach>
</body></html>
