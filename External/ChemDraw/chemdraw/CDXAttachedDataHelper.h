#pragma once

#include "CoreChemistryAPI.h"
#include "cs_univDefs.h"
#include "CDXStdObjects.h"
#include "ffMDL.h"
#include "BasicTypes.h"

class CORE_CHEMISTRY_API CDXAttachedDataHelper
{
public:
	enum AttachedDataTagPosition
	{
		AttachedDataTagPositionAuto,
		AttachedDataTagPositionTopLeft,
		AttachedDataTagPositionTopCenter,
		AttachedDataTagPositionTopRight,
		AttachedDataTagPositionMiddleLeft,
		AttachedDataTagPositionCenter,
		AttachedDataTagPositionMiddleRight,
		AttachedDataTagPositionBottomLeft,
		AttachedDataTagPositionBottomCenter,
		AttachedDataTagPositionBottomRight
	};

	static std::vector<std::string> GetSearchOperators();
	static std::vector<std::string> GetTagPositions();

	static std::string GetSearchOperatorByIndex(int index);
	static int GetSearchOperatorIndex(const std::string& op);

	static AttachedDataTagPosition ReadTagPosition(const CDXText& text);
	static AttachedDataTagPosition ReadTagPosition(const CDXPage& page, CDXObjectID propertyId);
	static int ReadQueryOperatorIndex(const CDXText& text);

	static void AddAttachedDataInfoTags(CDXText& propText, const std::string& updatedQueryOperator, char tag, AttachedDataTagPosition tagPosition, CDXObjectID& objectId);
	static void AddObjectTagsToBasisObjects(const CDXObject& cdxParent, const CDXChemicalProperty& cdxProp, char tag, AttachedDataTagPosition tagPosition, CDXObjectID& objectId);

  /* XXX Glysade - Unused
	static cs::PointTemplate2D<InternCoord> GetTagPositionForObject(const InternRect& tagBounds, const InternRect& objectBounds, AttachedDataTagPosition tagPosition);
  */
	static std::string GetDataFromDisplayText(const std::string& displayText, int searchOperatorIndex, const std::string& tag);

	static bool IsTag(char tag) { return tag != ' '; }

	static std::string MakeDisplayText(const std::string& data, char tag, const std::string& searchOperator);

private:
	static void AddAttachedDataObjectTag(CDXObject& object, CDXObjectID propertyId, char tag, AttachedDataTagPosition tagPosition, CDXObjectID& objectId);
};

