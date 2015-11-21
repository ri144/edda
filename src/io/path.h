#ifndef PATH_H_
#define PATH_H_

#include <string>
namespace edda{
std::string getPath(const std::string &file);
std::string getFilename(const std::string &filepath);
std::string removeFileExtension(const std::string &filename);
bool isFilenameOnly(const std::string &filename);
std::string getFileExtension(const std::string &filename);
}
#endif


