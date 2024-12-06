#pragma once

bool  csIsNetVariant (const char *s, long len);
bool  csIsValidBase64 (const char *s, long len);
long  csBase64Encode_Ext (const char	*input,
			  long			inputLen,
			  char			*&output,
			  bool			bNetVariant,
			  unsigned int	lineBreakFrequency = 0);
long  csBase64Decode(const char *input, long inputLen, char *&output, bool bNetVariantOnly);
void csSetBase64DecodeAuthorization (bool bAuthorized);
