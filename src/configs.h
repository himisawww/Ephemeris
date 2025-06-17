#pragma once
#include"definitions.h"

namespace Configs{
constexpr size_t MAX_LINESIZE=1024;
constexpr size_t MAX_PATHSIZE=260;

extern const int_t CheckPointMagic;
extern const int_t CheckPointVersion;
extern const int_t DataFormatVersion;

extern const int ExportHeaderCount;
extern const char *SaveNameCheckpoint;
extern const char *SaveNameReadme;
extern const char *SaveNameIndex;
extern const char *SaveNameDirectory;
extern const char *SaveRotationalDataExtension;
extern const char *SaveOrbitalDataExtension;
extern const char *SaveSubstepDataExtension;
extern const char *SaveNameInitialDirectory;
extern const char *VersionString;
extern const char AuthorName[];

//deprecated
extern const char *SaveNameStructure;
extern const char *SaveNameTimestamps;
extern const char *SaveDataExtension;
}
