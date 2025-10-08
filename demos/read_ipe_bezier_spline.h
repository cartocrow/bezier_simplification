#ifndef CARTOCROW_READIPEBEZIERSPLINE_H
#define CARTOCROW_READIPEBEZIERSPLINE_H

#include <cartocrow/core/cubic_bezier.h>

#include <ipeattributes.h>
#include <ipedoc.h>
#include <ipeshape.h>

#include <filesystem>

using namespace cartocrow;

std::vector<CubicBezierSpline> ipeSplinesToIsolines(const std::filesystem::path& file);
std::vector<CubicBezierSpline> splinesInPage(ipe::Page* page);
CubicBezierSpline customConvertPathToSpline(const ipe::SubPath& path, const ipe::Matrix& matrix);

#endif //CARTOCROW_READIPEBEZIERSPLINE_H
