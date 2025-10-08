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
			splines.emplace_back(customConvertPathToSpline(*subpath, matrix));
		}
	}

	return splines;
}

CubicBezierSpline customConvertPathToSpline(const ipe::SubPath& path, const ipe::Matrix& matrix) {
    CubicBezierSpline spline;
    if (path.type() == ipe::SubPath::EClosedSpline) {
        std::vector<ipe::Bezier> beziers;
        path.asClosedSpline()->beziers(beziers);
        for (auto bezier : beziers) {
            spline.appendCurve(
                    Point<Inexact>(bezier.iV[0].x, bezier.iV[0].y), Point<Inexact>(bezier.iV[1].x, bezier.iV[1].y),
                    Point<Inexact>(bezier.iV[2].x, bezier.iV[2].y), Point<Inexact>(bezier.iV[3].x, bezier.iV[3].y));
        }
    } else if (path.type() == ipe::SubPath::ECurve) {
        std::vector<ipe::Bezier> bzs;

        const ipe::Curve* curve = path.asCurve();

        if (curve->segment(0).type() == ipe::CurveSegment::ESegment) {
            std::vector<ipe::Vector> vs(curve->countSegments() + 1);
            for (int i = 0; i < curve->countSegments(); ++i) {
                vs[i] = curve->segment(i).cp(0);
            }
            vs[curve->countSegments()] = curve->segment(curve->countSegments() - 1).cp(1);

            std::vector<ipe::Vector> toRemove;
            for (int i = 1; i < curve->countSegments() - 2; ++i) {
                auto p0 = vs[i-1];
                auto p1 = vs[i];
                auto p2 = vs[i+1];
                auto p3 = vs[i+2];
                if ((p2-p1).len() / ((p1-p0).len() + (p3-p2).len()) < 1.0 / 10.0) {
                    // remove p1 or p2
                    toRemove.push_back(p1);
                }
            }
            auto rit = std::remove_if(vs.begin(), vs.end(), [&toRemove](const auto& v) { return std::find(toRemove.begin(), toRemove.end(), v) != toRemove.end(); });
            vs.erase(rit, vs.end());

            ipe::Bezier::cardinalSpline(vs.size(), &vs[0], 0.5, bzs);
        } else {
            for (int i = 0; i < curve->countSegments(); ++i) {
                curve->segment(i).beziers(bzs);
            }
        }

        for (const auto& bz : bzs) {
            spline.appendCurve(
                    Point<Inexact>(bz.iV[0].x, bz.iV[0].y), Point<Inexact>(bz.iV[1].x, bz.iV[1].y),
                    Point<Inexact>(bz.iV[2].x, bz.iV[2].y), Point<Inexact>(bz.iV[3].x, bz.iV[3].y));
        }
    } else {
        throw std::runtime_error("Ipe SubPath is neither a ClosedSpline nor a Curve.");
    }
    return spline;
}