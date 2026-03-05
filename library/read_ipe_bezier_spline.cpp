#include "read_ipe_bezier_spline.h"

#include <cartocrow/reader/ipe_reader.h>

#include <ipebase.h>
#include <ipeshape.h>
#include <ipepath.h>

std::vector<CubicBezierSpline> ipeSplinesToIsolines(const std::filesystem::path& file) {
	std::shared_ptr<ipe::Document> document = IpeReader::loadIpeFile(file);

	if (document->countPages() == 0) {
		throw std::runtime_error("Cannot read isolines from an Ipe file with no pages");
	} else if (document->countPages() > 1) {
		throw std::runtime_error("Cannot read isolines from an Ipe file with more than one page");
	}

	ipe::Page* page = document->page(0);
	return splinesInPage(page);
}

std::vector<CubicBezierSpline> splinesInPage(ipe::Page* page) {
	auto splines = std::vector<CubicBezierSpline>();

	for (int i = 0; i < page->count(); i++) {
		auto object = page->object(i);
		if (object->type() != ipe::Object::Type::EPath) continue;
		auto path = object->asPath();
		auto matrix = object->matrix();
		auto shape = path->shape();
		for (int j = 0; j < shape.countSubPaths(); j++) {
			auto subpath = shape.subPath(j);
			if (subpath->type() != ipe::SubPath::Type::ECurve) continue;
			splines.emplace_back(convertPathToSpline(*subpath, matrix));
		}
	}

	return splines;
}


CubicBezierSpline convertPathToSpline(const ipe::SubPath& path, const ipe::Matrix& matrix) {
    CubicBezierSpline spline;
    if (path.type() == ipe::SubPath::EClosedSpline) {
        std::vector<ipe::Bezier> beziers;
        path.asClosedSpline()->beziers(beziers);
        for (auto bezier : beziers) {
            ipe::Vector v0 = matrix * bezier.iV[0];
            ipe::Vector v1 = matrix * bezier.iV[1];;
            ipe::Vector v2 = matrix * bezier.iV[2];
            ipe::Vector v3 = matrix * bezier.iV[3];
            spline.appendCurve(
                    Point<Inexact>(v0.x, v0.y), Point<Inexact>(v1.x, v1.y),
                    Point<Inexact>(v2.x, v2.y), Point<Inexact>(v3.x, v3.y));
        }
    } else if (path.type() == ipe::SubPath::ECurve) {
        const ipe::Curve* curve = path.asCurve();

        for (int i = 0; i < curve->countSegmentsClosing(); ++i) {
            ipe::CurveSegment seg = curve->segment(i);
            addCurveSegmentToSpline(seg, spline, matrix);
        }
    } else {
        throw std::runtime_error("Ipe SubPath is neither a ClosedSpline nor a Curve.");
    }
    return spline;
}

void addCurveSegmentToSpline(const ipe::CurveSegment& seg, CubicBezierSpline& spline, const ipe::Matrix& matrix) {
    if (seg.type() == ipe::CurveSegment::ESegment) {
        spline.appendCurve(ipeVectorToPoint(matrix * seg.cp(0)), ipeVectorToPoint(matrix * seg.cp(1)));
    } else {
        std::vector<ipe::Bezier> bzs;
        seg.beziers(bzs);
        for (const auto& bz : bzs) {
            spline.appendCurve(ipeVectorToPoint(matrix * bz.iV[0]),
                               ipeVectorToPoint(matrix * bz.iV[1]),
                               ipeVectorToPoint(matrix * bz.iV[2]),
                               ipeVectorToPoint(matrix * bz.iV[3]));
        }
    }
}

Point<Inexact> ipeVectorToPoint(const ipe::Vector& v) {
    return {v.x, v.y};
}
