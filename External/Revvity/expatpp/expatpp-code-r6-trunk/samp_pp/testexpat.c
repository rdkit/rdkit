/* This is simple demonstration of how to use expat. This program
reads an XML document from standard input and writes a line with the
name of each element to standard output indenting child elements by
one tab stop more than their parent element. 

Copied from elements.c 98/08/31 AD to open specified file rather than stdin
and added a CharacterDataHandler
*/

#include <stdio.h>
#include "xmlparse.h"

void startElement(void *userData, const char *name, const char **atts)
{
  int i;
  int *depthPtr = userData;
  for (i = 0; i < *depthPtr; i++)
    putchar('\t');
  puts(name);
  *depthPtr += 1;
}

void endElement(void *userData, const char *name)
{
  int *depthPtr = userData;
  *depthPtr -= 1;
}


void charData(void *userData, const XML_Char *s, int len)
{
/* write indent level stored in userData */
  int i;
  int *depthPtr = userData;
  
  if (len==0)
  	return;  // lots of extra calls?
  	
  for (i = 0; i < *depthPtr; i++)
    putchar('\t');

/* write out the user data bracketed by ()*/
  putchar('(');
  fwrite(s, len, 1, stdout);
  puts(")");
}


int main()
{	
  char filename[80];
  FILE* xmlFile;
  char buf[BUFSIZ];
  for (;;) {	  
	  XML_Parser parser;
	  int done;
	  int depth = 0;

	  puts("\n\nXML test: enter filename");
	  gets(filename);
	  if (strlen(filename)==0)
	  	break;
	  	
	  xmlFile = fopen(filename, "r");
	  if (!xmlFile)
	  	break;

	  parser = XML_ParserCreate(NULL);
	  XML_SetUserData(parser, &depth);
	  XML_SetElementHandler(parser, startElement, endElement);
	  XML_SetCharacterDataHandler(parser, charData);
	  do {
	    size_t len = fread(buf, 1, sizeof(buf), xmlFile /*stdin*/);
	    done = len < sizeof(buf);
	    if (!XML_Parse(parser, buf, len, done)) {
	      fprintf(stderr,
		      "%s at line %d\n",
		      XML_ErrorString(XML_GetErrorCode(parser)),
		      XML_GetCurrentLineNumber(parser));
	      return 1;
	    }
	  } while (!done);
	  XML_ParserFree(parser);
/*	  return 0;  */
  }
  puts("\nfinished!");
}