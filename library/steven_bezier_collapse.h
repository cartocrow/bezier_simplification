#pragma once

#include "bezier_collapse.h"
#include "fit_cubic.h"

#include <cartocrow/core/arrangement_helpers.h>

// for debugging only; todo: remove
#include <cartocrow/renderer/ipe_renderer.h>

namespace cartocrow::curved_simplification {
Number<Inexact> getAreaD0(Number<Inexact> a, Point<Inexact> p0, Direction<Inexact> t0, Direction<Inexact> t1, Number<Inexact> d1, Point<Inexact> p3);
Number<Inexact> getAreaD1(Number<Inexact> a, Point<Inexact> p0, Direction<Inexact> t0, Number<Inexact> d0, Direction<Inexact> t1, Point<Inexact> p3);
CubicBezierCurve
createCubicBezierFromPolar(Point<Inexact> p0, Direction<Inexact> t0Dir, Number<Inexact> d0,
						   Direction<Inexact> t1Dir, Number<Inexact> d1, Point<Inexact> p3);
Number<Inexact> length(const Polyline<Inexact>& pl);
bool isStraight(const CubicBezierCurve& curve);
// Whether curve c1 connects smoothly to curve c2: whether the end tangent of c1 aligns approximately with the
// start tangents of c2
bool connectsSmoothlyTo(const CubicBezierCurve& c1, const CubicBezierCurve& c2);

template <typename BG> struct StevenBCTraits {
	int tSteps;
	int nSegs;
	bool debug;

	StevenBCTraits(bool debug = false, int tSteps = 10,	int nSegs = 100) :
	      tSteps(tSteps), nSegs(nSegs), debug(debug) {};

	void determineCollapse(typename BG::Edge_handle e) {
		auto& edata = e->data();
		if (e->source()->degree() != 2 || e->target()->degree() != 2) {
			edata.collapse = std::nullopt;
			return;
		}
		const auto& c0 = e->prev()->curve();
		const auto& c1 = e->curve();
		const auto& c2 = e->next()->curve();

		// Determine cost, point, before, after for edge e.
		CubicBezierSpline beforeSpline;
		beforeSpline.appendCurve(c0);
		beforeSpline.appendCurve(c1);
		beforeSpline.appendCurve(c2);
		auto beforePl = beforeSpline.polyline(nSegs);

		// We keep the source and target.
		// We keep the tangents at source and target to ensure G^1 continuity
		auto p0 = c0.source();
		auto p1 = c0.sourceControl();
		auto p5 = c2.targetControl();
		auto p6 = c2.target();
		// Goal: find p2, p3, p4

		double minErr = std::numeric_limits<double>::infinity();
		std::optional<std::pair<CubicBezierCurve, CubicBezierCurve>> best = std::nullopt;

		renderer::IpeRenderer ipeRenderer;

		// choose a point on c1 to use as p3.
		// for now just equidistant t; should this be made adaptive (todo)?
		for (int it = 0; it < tSteps; ++it) {
			double t = tSteps == 1 ? 0.5 : static_cast<double>(it) / (tSteps - 1);
			auto p3 = c1.evaluate(t);
			auto v = c1.tangent(t);
            CubicBezierSpline firstPart;
            CubicBezierSpline secondPart;
            if (t < M_EPSILON) {
                firstPart.appendCurve(c0);
                secondPart.appendCurve(c1);
                secondPart.appendCurve(c2);
            } else if (t > 1 - M_EPSILON) {
                firstPart.appendCurve(c0);
                firstPart.appendCurve(c1);
                secondPart.appendCurve(c2);
            } else {
                auto [c11, c12] = c1.split(t);
                firstPart.appendCurve(c0);
                firstPart.appendCurve(c11);
                secondPart.appendCurve(c12);
                secondPart.appendCurve(c2);
            }
            auto firstPartPl = firstPart.polyline(nSegs);
            std::vector<Point<Inexact>> firstPartPts(firstPartPl.vertices_begin(), firstPartPl.vertices_end());
            auto secondPartPl = secondPart.polyline(nSegs);
            std::vector<Point<Inexact>> secondPartPts(secondPartPl.vertices_begin(), secondPartPl.vertices_end());
            auto firstPartFit = firstPart.numCurves() == 1 ? firstPart.curve(0) : fitCurve(firstPartPts, c0.sourceControl() - p0, -v);
            auto secondPartFit = secondPart.numCurves() == 1 ? secondPart.curve(0) : fitCurve(secondPartPts, v, c2.targetControl() - p6);
            p1 = firstPartFit.sourceControl();
            p5 = secondPartFit.targetControl();

			auto a0 = -firstPart.signedArea();
			auto a1 = -secondPart.signedArea();

            std::optional<CubicBezierCurve> c0_maybe;
            if (a0 == 0) {
                c0_maybe = CubicBezierCurve(p0, p3);
            } else {
                auto v0 = p1 - p0;
                auto t00 = v0.direction();
                auto d00 = sqrt(v0.squared_length());
                auto t01Wrong = (-v).direction();
                Direction<Inexact> t01(-t01Wrong.dx(), t01Wrong.dy());
                auto d01 = getAreaD1(a0, p0, t00, d00, t01, p3);
                if (d01 > 0) {
                    c0_maybe = createCubicBezierFromPolar(p0, t00, d00, t01, d01, p3);
                }
            }

            if (!c0_maybe.has_value()) continue;

            std::optional<CubicBezierCurve> c1_maybe;
            if (a1 == 0) {
                c1_maybe = CubicBezierCurve(p3, p6);
            } else {
                auto v1 = p6 - p5;
                auto t11Wrong = v1.direction();
                Direction<Inexact> t11(t11Wrong.dx(), -t11Wrong.dy());
                auto d11 = sqrt(v1.squared_length());
                auto t10 = v.direction();
                auto d10 = getAreaD0(a1, p3, t10, t11, d11, p6);
                if (d10 > 0) {
                    c1_maybe = createCubicBezierFromPolar(p3, t10, d10, t11, d11, p6);
                }
            }

            if (!c1_maybe.has_value()) continue;
            auto& c0_ = *c0_maybe;
            auto& c1_ = *c1_maybe;

			CubicBezierSpline afterSpline;
			afterSpline.appendCurve(c0_);
			afterSpline.appendCurve(c1_);
			auto afterPl = afterSpline.polyline(10);
            auto afterPlE = pretendExact(afterPl);

			double maxKappaBefore = 0;
			for (int i = 0; i <= nSegs; ++i) {
				double ti = static_cast<double>(i) / nSegs;
				double kappa = std::max(std::max(abs(c0.curvature(ti)), abs(c1.curvature(ti))), abs(c2.curvature(ti)));
				if (kappa > maxKappaBefore) {
					maxKappaBefore = kappa;
				}
			}

			double maxKappaAfter = 0;
            // Ignore start and end curvature
			for (int i = nSegs/10; i <= 9 * nSegs / 10; ++i) {
//			for (int i = 0; i <= nSegs; ++i) {
				double ti = static_cast<double>(i) / nSegs;
				double kappa = std::max(abs(c0_.curvature(ti)), abs(c1_.curvature(ti)));
				if (kappa > maxKappaAfter) {
					maxKappaAfter = kappa;
				}
			}

			double err, symDiffErr, curvatureErr = 0;
            bool createsSmoothConnection = false;
            createsSmoothConnection = !connectsSmoothlyTo(c0, c1) && !connectsSmoothlyTo(c1, c2) && connectsSmoothlyTo(c0_, c1_);
			if ((createsSmoothConnection || maxKappaAfter < maxKappaBefore) && !CGAL::do_curves_intersect(afterPlE.edges_begin(), afterPlE.edges_end())) {
				Arrangement<Exact> arr;
				std::vector<Arrangement<Exact>::X_monotone_curve_2> beforePlXMCurves;
				for (auto eit = beforePl.edges_begin(); eit != beforePl.edges_end(); ++eit) {
					beforePlXMCurves.emplace_back(pretendExact(*eit));
				}
				// todo use insert_non_intersecting_curves once it is guaranteed that the input curve is planar
				CGAL::insert(arr, beforePlXMCurves.begin(), beforePlXMCurves.end());
				CGAL::insert(arr, afterPlE.edges_begin(), afterPlE.edges_end());
                symDiffErr = 0;
				for (auto fit = arr.faces_begin(); fit != arr.faces_end(); ++fit) {
					if (fit->is_unbounded()) continue;
					auto pwh = approximate(face_to_polygon_with_holes<Exact>(fit));
					symDiffErr += abs(pwh.outer_boundary().area());
					for (const auto& h : pwh.holes()) {
						symDiffErr -= abs(h.area());
					}
				}

//                curvatureErr = - maxKappaBefore / maxKappaAfter;
                if (isStraight(c0) || isStraight(c1) || isStraight(c2)) {
                    curvatureErr = 0;
                }
                if (debug) {
                    std::cerr << "Sym diff err (sqrt): " << sqrt(symDiffErr) << "\n";
                    std::cerr << "Curvature err: " << curvatureErr << "; before: " << maxKappaBefore << "; after: " << maxKappaAfter << "\n";
                }
				err = sqrt(symDiffErr);// + 5 * curvatureErr;

                if (createsSmoothConnection) {
                    err /= 10;
                }
			} else {
				err = std::numeric_limits<double>::infinity();
			}

			if (err < minErr) {
				minErr = err;
				best = {c0_, c1_};
			}

			if (debug) { //&& err < std::numeric_limits<double>::infinity()) {
				ipeRenderer.addPainting([beforeSpline](renderer::GeometryRenderer& renderer) {
					renderer.setMode(renderer::GeometryRenderer::vertices | renderer::GeometryRenderer::stroke);
					renderer.setStroke(Color(200, 200, 200), 3.0);
					renderer.draw(beforeSpline);
				}, "Before");
				ipeRenderer.addPainting([afterSpline, err, symDiffErr, curvatureErr, createsSmoothConnection, maxKappaBefore, maxKappaAfter](renderer::GeometryRenderer& renderer) {
					renderer.setStroke(Color(0, 0, 0), 3.0);
					renderer.draw(afterSpline);
					renderer.drawText({10, 10}, "Err: " + std::to_string(err));
                    renderer.drawText({10, 20}, "Sym diff err: " + std::to_string(sqrt(symDiffErr)));
                    renderer.drawText({10, 30}, "Max. curvature before: " + std::to_string(maxKappaBefore));
                    renderer.drawText({10, 40}, "Max. curvature after: " + std::to_string(maxKappaAfter));
                    if (createsSmoothConnection) {
                        renderer.drawText({10, 50}, "Creates smooth connection!");
                    }
				}, "After");
				ipeRenderer.nextPage();
			}
		}

		if (debug) {
			ipeRenderer.save("debugging.ipe");
		}

        if (best.has_value()) {
		    edata.collapse = { .cost=minErr, .before=best->first, .after=best->second };
        } else {
            edata.collapse = std::nullopt;
        }
	}
};
}
