#include "CDXAttachedDataHelper.h"
#include "cs_stringUtils.h"
#include "CDXUtils.h"

#if TARGET_OS_WIN32
   #define strcasecmp _strcmpi
#endif

const char *searchOperators[] = { "None", ">", ">=", "<", "<=", "<>", "between", "contains", "like" };
const char *tagPositions[] = { "Auto", "Top Left", "Top Center", "Top Right", "Middle Left", "Center", "Middle Right", "Bottom Left", "Bottom Center", "Bottom Right" };

std::vector<std::string> CDXAttachedDataHelper::GetSearchOperators()
{
	return std::vector<std::string>(searchOperators, searchOperators + (sizeof(searchOperators) / sizeof(searchOperators[0])));
}

std::vector<std::string> CDXAttachedDataHelper::GetTagPositions()
{
	return std::vector<std::string>(tagPositions, tagPositions + (sizeof(tagPositions) / sizeof(tagPositions[0])));
}

std::string CDXAttachedDataHelper::GetSearchOperatorByIndex(int index)
{
	if ((index > 0) && (index < (sizeof(searchOperators) / sizeof(searchOperators[0]))))
	{
		return searchOperators[index];
	}
	else
	{
		return "";
	}
}

int CDXAttachedDataHelper::GetSearchOperatorIndex(const std::string& op)
{
	const char* searchOp = op.c_str();
	int totalItems = (sizeof(searchOperators) / sizeof(searchOperators[0]));
	for (int index = 0; index < totalItems; index++)
	{
		if (strcasecmp(searchOp, searchOperators[index]) == 0)
		{
			return index;
		}
	}

	return -1;
}

CDXAttachedDataHelper::AttachedDataTagPosition CDXAttachedDataHelper::ReadTagPosition(const CDXText& text)
{
	// Get tag position index from the tag position object tag. If tag is 
	// not found we default to 5 (Center).
	CDXAttachedDataHelper::AttachedDataTagPosition tagPosition = AttachedDataTagPositionCenter;
	std::vector<CDXObjectTag*> tagPosTags = GetObjectTags(text, kCDXTagType_AttachedDataTagPosition);
	if (tagPosTags.size() == 1)
	{
		// Get the first item from list. Ideally the list will 
		// only have one tag of tag position type
		CDXObjectTag* tagPositionTag = (*tagPosTags.begin());
		tagPosition = (CDXAttachedDataHelper::AttachedDataTagPosition)tagPositionTag->m_Int32Val;
	}

	return tagPosition;
}

CDXAttachedDataHelper::AttachedDataTagPosition CDXAttachedDataHelper::ReadTagPosition(const CDXPage& page, CDXObjectID propertyId)
{
	CDXAttachedDataHelper::AttachedDataTagPosition tagPos = CDXAttachedDataHelper::AttachedDataTagPositionCenter;

	// Attached data tag's value points to the main attached data chemical property
	CDXChemicalProperty* prop = dynamic_cast<CDXChemicalProperty*>(page.FindByIDRecursive(propertyId));
	if (prop)
	{
		// Get the CDXText associated with this property
		CDXText* text = dynamic_cast<CDXText*>(page.FindByIDRecursive(prop->m_displayID));
		if (text)
		{
			tagPos = CDXAttachedDataHelper::ReadTagPosition(*text);
		}
	}

	return tagPos;
}

int CDXAttachedDataHelper::ReadQueryOperatorIndex(const CDXText& text)
{
	// Get search operator index (aka query operator) from the search operator object tag
	// If the operator is not found we default to 0 (None)
	int queryOperatorIndex = 0;
	std::vector<CDXObjectTag*> searchOptags = GetObjectTags(text, kCDXTagType_AttachedDataSearchOperator);
	if (searchOptags.size() == 1)
	{
		// Get the first item from list. Ideally the list will 
		// only have one tag of search operator type
		CDXObjectTag* queryOperatorTag = (*searchOptags.begin());
		queryOperatorIndex = queryOperatorTag->m_Int32Val;
	}

	return queryOperatorIndex;
}

/*
* AddAttachedDataTags
* The function makes the following hierarchy:
* <t>		- (propText) The text attached to <chemicalproperty> for attached data
*	<objecttag Name="attachedDataSearchOperator">
*	<objecttag Name="attachedDataTagPosition">
*/
void CDXAttachedDataHelper::AddAttachedDataInfoTags(CDXText& propText, const std::string& updatedQueryOperator, char tag, AttachedDataTagPosition tagPosition, CDXObjectID& objectId)
{
	// Add the search operator object tag under property text
	if (!updatedQueryOperator.empty())
	{
		CDXObjectTag* searchOperatorTag = new CDXObjectTag(++objectId);
		searchOperatorTag->m_Name = kCDXTagType_AttachedDataSearchOperator;
		searchOperatorTag->m_Type = kCDXObjectTagType_Int32;
		int searchOpIndex = CDXAttachedDataHelper::GetSearchOperatorIndex(updatedQueryOperator);
		if (searchOpIndex >= 0)
		{
			searchOperatorTag->m_Int32Val = searchOpIndex;
		}
		propText.AddChild(searchOperatorTag);
	}

	if (IsTag(tag))
	{
		// Add the tag position object tag under property text
		CDXObjectTag* tagPositionTag = new CDXObjectTag(++objectId);
		tagPositionTag->m_Name = kCDXTagType_AttachedDataTagPosition;
		tagPositionTag->m_Type = kCDXObjectTagType_Int32;
		tagPositionTag->m_Int32Val = tagPosition;
		propText.AddChild(tagPositionTag);
	}
}

/*
* AddAttachedDataObjectTag
* The function adds an <objecttag> to the specified CDXObject.
* <n>|<b>	- The target node or bond to which data is attached
*	<objecttag Name="attachedData">
*/
void CDXAttachedDataHelper::AddAttachedDataObjectTag(CDXObject& object, CDXObjectID propertyId, char tag, AttachedDataTagPosition tagPosition, CDXObjectID& objectId)
{
	// Create the main object tag that holds the tag string
	CDXObjectTag* objectTag = new CDXObjectTag(++objectId);
	objectTag->m_Name = kCDXTagType_AttachedData;
	objectTag->m_Type = kCDXObjectTagType_Int32;
	objectTag->m_Int32Val = propertyId;
	
	// Create tag string and add it to object tag
	CDXText* tagText = new CDXText(++objectId);
	tagText->SetText(CDXString(std::string(1, tag)));
	objectTag->AddChild(tagText);

	object.AddChild(objectTag);
}

/*
* AddObjectTagsToBasisObjects
* Adds attached data tags to atoms, bonds and brackets. Used by both MOL and SKC code.
*/
void CDXAttachedDataHelper::AddObjectTagsToBasisObjects(const CDXObject& cdxParent, const CDXChemicalProperty& cdxProp, char tag, AttachedDataTagPosition tagPosition, CDXObjectID& objectId)
{
	std::vector<CDXNode*> atoms;
	std::vector<CDXObjectID> usedObjects;
	CDXObjectID propObjectId = cdxProp.GetObjectID();
	const std::vector<CDXObjectID>& basisObjects = cdxProp.m_basisObjects;

	// Go through all the objects and add tag to bonds and brackets first
	for (std::vector<CDXObjectID>::const_iterator it = basisObjects.begin(); it != basisObjects.end(); ++it)
	{
		CDXObject* cdxObject = cdxParent.FindByIDRecursive(*it);
		if (dynamic_cast<CDXBond*>(cdxObject))
		{
			CDXBond* bond = (CDXBond*)cdxObject;

			usedObjects.push_back(bond->m_beginNodeID);
			usedObjects.push_back(bond->m_endNodeID);

			CDXAttachedDataHelper::AddAttachedDataObjectTag(*cdxObject, propObjectId, tag, tagPosition, objectId);
		}
		else if (dynamic_cast<CDXNode*>(cdxObject))
		{
			atoms.push_back((CDXNode*)cdxObject);
		}
		else if (dynamic_cast<CDXGraphic*>(cdxObject))
		{
			CDXGraphic* graphic = (CDXGraphic*)cdxObject;
			if (graphic->m_graphicType == kCDXGraphicType_Bracket)
			{
				// Check if this brackets has alredy been used up
				if (std::find(usedObjects.begin(), usedObjects.end(), graphic->GetObjectID()) == usedObjects.end())
				{
					CDXAttachedDataHelper::AddAttachedDataObjectTag(*cdxObject, propObjectId, tag, tagPosition, objectId);

					// Add the other bracket to used objects list to that tag is added to just one bracket in the group
					std::vector<CDXObjectID> otherBrackets = graphic->GetOtherBracketsInGroup();
					usedObjects.insert(usedObjects.end(), otherBrackets.begin(), otherBrackets.end());
				}
			}
		}
	}

	// Now add tag for any remaining atoms that were not used by other objects such as bonds
	for (std::vector<CDXNode*>::const_iterator it = atoms.begin(); it != atoms.end(); ++it)
	{
		if (std::find(usedObjects.begin(), usedObjects.end(), (*it)->GetObjectID()) == usedObjects.end())
		{
			CDXAttachedDataHelper::AddAttachedDataObjectTag(*(*it), propObjectId, tag, tagPosition, objectId);
		}
	}
}

std::string CDXAttachedDataHelper::GetDataFromDisplayText(const std::string& displayText, int searchOperatorIndex, const std::string& tag)
{
	std::string text = displayText;

	// If there is a search operator and there is a = sign before the search operator then
	// it can be assumed that tag exists even though the 'tag' property is not set.
	// This can happen when an object (atom, bond etc) containing an object tag was deleted
	// by the user, causing object tag for that object to be deleted silently. There is 
	// currently no way to catch that in code effectively for all objects.
	std::string::size_type tagLength = tag.length();
	bool tagPossiblyExists = false;
	if (tag.empty() && (searchOperatorIndex > 0))
	{
		std::string searchOp = CDXAttachedDataHelper::GetSearchOperatorByIndex(searchOperatorIndex);
		std::size_t searchOpPos = displayText.find_first_of(searchOp);
		std::size_t equalsToPos = displayText.find_first_of('=');
		if ((equalsToPos != std::string::npos) && (equalsToPos < searchOpPos))
		{
			// There is an equals to sign in the display string and it appears before the search operator
			tagLength = 1;
			tagPossiblyExists = true;
		}
	}


	if (!tag.empty() || tagPossiblyExists)
	{
		// tag is present. Remove the tag plus three characters ( = ) from front
		std::string::size_type subLength = tagLength + 3; 
		if (text.length() >= subLength)
		{
			text = text.substr(subLength);
		}
	}

	if (searchOperatorIndex > 0)
	{
		// search operator is followed by a space that we will skip
		std::string searchOp = CDXAttachedDataHelper::GetSearchOperatorByIndex(searchOperatorIndex);
		std::string::size_type subLength = searchOp.length() + 1;
		if (text.length() >= subLength)
		{
			text = text.substr(subLength);
		}
	}

	return text;
}

/* XXX GLYSADE - UNUSED
cs::PointTemplate2D<InternCoord> CDXAttachedDataHelper::GetTagPositionForObject(const InternRect& tagBounds, const InternRect& objectBounds, AttachedDataTagPosition tagPosition)
{
	InternCoord tagWidth = tagBounds.Width();
	InternCoord tagHeight = tagBounds.Height();
	cs::PointTemplate2D<InternCoord> tagCoordinate = objectBounds.Center();

	switch (tagPosition)
	{
		case AttachedDataTagPositionTopLeft:
			tagCoordinate = objectBounds.TopLeft();
			tagCoordinate.Offset(-tagWidth, 0);
			break;

		case AttachedDataTagPositionTopCenter:
			tagCoordinate = objectBounds.TopCenter();
			tagCoordinate.Offset(-tagWidth/2, 0);
			break;

		case AttachedDataTagPositionTopRight:
			tagCoordinate = objectBounds.TopRight();
			break;

		case AttachedDataTagPositionMiddleLeft:
			tagCoordinate = objectBounds.CenterLeft();
			tagCoordinate.Offset(-tagWidth, tagHeight/2);
			break;

		case AttachedDataTagPositionAuto:		// Auto defaults to Center
		case AttachedDataTagPositionCenter:
			tagCoordinate = objectBounds.Center();
			tagCoordinate.Offset(-tagWidth/2, tagHeight/2);
			break;

		case AttachedDataTagPositionMiddleRight:
			tagCoordinate = objectBounds.CenterRight();
			tagCoordinate.Offset(0, tagHeight/2);
			break;

		case AttachedDataTagPositionBottomLeft:
			tagCoordinate = objectBounds.BottomLeft();
			tagCoordinate.Offset(-tagWidth, tagHeight);
			break;

		case AttachedDataTagPositionBottomCenter:
			tagCoordinate = objectBounds.BottomCenter();
			tagCoordinate.Offset(-tagWidth/2, tagHeight);
			break;

		case AttachedDataTagPositionBottomRight:
			tagCoordinate = objectBounds.BottomRight();
			tagCoordinate.Offset(0, tagHeight);
			break;
	}

	return tagCoordinate;
}
*/
std::string CDXAttachedDataHelper::MakeDisplayText(const std::string& data, char tag, const std::string& searchOperator)
{
	std::string displayText = data;

	// Make display text as TAG = SEARCHOP DATA
	if (!searchOperator.empty())
	{
		displayText = searchOperator + " " + displayText;
	}

	if (IsTag(tag))
	{
		displayText = std::string(1, tag) + " = " + displayText;
	}

	return displayText;
}
