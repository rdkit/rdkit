<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN" "http://www.w3.org/TR/html4/loose.dtd">
<%@ taglib prefix="c" uri="http://java.sun.com/jsp/jstl/core" %>
<%@ taglib prefix="x" uri="http://java.sun.com/jsp/jstl/xml" %>
<html><body>

<c:set var="renderer" value="../depicter/mol-renderer" />

<c:set var="size" value="100" />
<c:if test="${!empty param.size}">
    <c:set var="size" value="${param.size}" />
</c:if>
<c:import var="source" url="depicter_source.xml" />
<x:parse var="smiles_list" xml="${source}" />
<x:forEach select="$smiles_list//smiles">
<c:set var="irow" value="${irow+1}" />
<c:set var="x" value="${x+size+2}" />
    <c:set var="smiles"><x:out select="." /></c:set>
    <img width=${size} height=${size} src="${renderer}/test.png?smiles=${smiles}&w=${size}&h=${size}"/>
    <c:if test="${x > 1024}"><br /> <c:set var="x" value="0" /> </c:if>
</x:forEach>
</body></html>
