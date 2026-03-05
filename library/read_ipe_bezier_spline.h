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
CubicBezierSpline convertPathToSpline(const ipe::SubPath& path, const ipe::Matrix& matrix);
Point<Inexact> ipeVectorToPoint(const ipe::Vector& v);
/// Precondition: the CurveSegment is not of type EArc
void addCurveSegmentToSpline(const ipe::CurveSegment& seg, CubicBezierSpline& spline, const ipe::Matrix& matrix);
#endif //CARTOCROW_READIPEBEZIERSPLINE_H
