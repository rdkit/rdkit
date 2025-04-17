// CommonCS/LibCommon/Hdr/cs_specialchardefs.h
// Copyright © 1998-2004, CambridgeSoft Corp., All Rights Reserved

#pragma once

// Special codes used in atom labels
// These should be non-printing characters unused for any other purpose
#define kCharMoveUp      1
#define kCharMoveDown    2
#define kCharSymb        3
#define kCharBold        4
#define kCharItalic      5
#define kCharLarger      6
#define kCharSmaller     7
#define kCharFormula    11

// Names for non-printing characters used for various purposes
#define kCharEnter       3
#define kCharBackSpace   8
#define kCharTab         9
#define kCharLineFeed   10
#define kCharOption		11
#define kCharFormFeed   12
//#define kCharReturn     13
#define kCharOpenApple  16
//#define kCharCommand    17
#define kCharCheck      18
#define kCharDiamond    19
#define kCharApple      20
//#define kCharShift		21
#define kCharPaste		12
#define kCharReturn     13
#define kCharHelp		14
#define kCharUndo		15
#define kCharCut		16
#define kCharControl    18
#if TARGET_OS_MAC || TARGET_WEB
#define kCharCommand	17
#else // TARGET_OS_MAC
#define kCharCommand kCharControl
#endif // TARGET_OS_MAC
#define kCharShift      19
#define kCharCopy		20
#if TARGET_OS_MAC
#define kCharPrior		kPageUpCharCode
#define kCharNext		kPageDownCharCode
#define kCharEnd		kEndCharCode
#define kCharHome		kHomeCharCode
#else // TARGET_OS_MAC
#define kCharPrior		21
#define kCharNext		22
#define kCharEnd		23
#define kCharHome		24
#endif // TARGET_OS_MAC
#define kCharInsert		25
#define kCharEscape     26
#define kCharClear      27
#define kCharArrowLeft  28
#define kCharArrowRight 29
#define kCharArrowUp    30
#define kCharArrowDown  31
#define kCharSpace		32
#define kCharDelete    127

enum
{
    kUnicodeCodePointDegree               = 0x00B0, // °
    kUnicodeCodePointCent                 = 0x00A2, // ¢
    kUnicodeCodePointPound                = 0x00A3, // £ 
    kUnicodeCodePointCopyright            = 0x00A9, // ©
    kUnicodeCodePointReg                  = 0x00AE, // ®
    kUnicodeCodePointPlusMinus            = 0x00B1, // ±
    kUnicodeCodePointMicro                = 0x00B5, // µ
    kUnicodeCodePointEnDash               = 0x2013, // –
    kUnicodeCodePointEmDash               = 0x2014, // —
    kUnicodeCodePointBullet               = 0x2022, // •
    kUnicodeCodePointEllipsis             = 0x2026, // …
    kUnicodeCodePointAngstrom             = 0x212B, // Å
    kUnicodeCodePointCenterDot            = 0x2219, // ·
    kUnicodeCodePointReplacement          = 0xFFFD
};

