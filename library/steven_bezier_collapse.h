#pragma once

#include "bezier_collapse.h"
#include "fit_cubic.h"
#include "spiro_helpers/spiro_helpers.h"

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
        edata.collapse = std::nullopt;
		if (e->source()->degree() != 2 || e->target()->degree() != 2) {
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

        // Degenerate cases
        if (c0.source() == c0.target()) {
            edata.collapse = { .cost=0, .before=c1, .after=c2 };
            return;
        }
        if (c1.source() == c1.target()) {
            edata.collapse = { .cost=0, .before=c0, .after=c2 };
            return;
        }
        if (c2.source() == c2.target()) {
            edata.collapse = { .cost=0, .before=c0, .after=c1 };
            return;
        }

		double minErr = std::numeric_limits<double>::infinity();
		std::optional<std::pair<CubicBezierCurve, CubicBezierCurve>> best = std::nullopt;

		renderer::IpeRenderer ipeRenderer;

        auto evaluate = [&](const CubicBezierCurve& c0_, const CubicBezierCurve& c1_) {
            CubicBezierSpline afterSpline;
            afterSpline.appendCurve(c0_);
            afterSpline.appendCurve(c1_);
            auto afterPl = afterSpline.polyline(10);
            auto afterPlE = pretendExact(afterPl);

            double maxKappaBefore = 0;
            for (int curveIndex = 0; curveIndex < beforeSpline.numCurves(); ++curveIndex) {
                for (int i = 0; i <= nSegs; ++i) {
                    double ti = static_cast<double>(i) / nSegs;
                    double kappa = abs(beforeSpline.curvature({curveIndex, ti}));
                    if (kappa > maxKappaBefore) {
                        maxKappaBefore = kappa;
                    }
                }
            }

            double maxKappaAfter = 0;
            double avgKappaAfter = 0;
            // Ignore start and end curvature
//			for (int i = nSegs/10; i <= 9 * nSegs / 10; ++i) {

            std::vector<CubicBezierSpline::SplineParameter> afterSplineInflectionsT;
            afterSpline.inflectionsT(std::back_inserter(afterSplineInflectionsT));
            std::vector<double> curvatureTotals(afterSplineInflectionsT.size() + 1);
            for (int curveIndex = 0; curveIndex < afterSpline.numCurves(); ++curveIndex) {
                for (int i = 0; i <= nSegs; ++i) {
                    double ti = static_cast<double>(i) / nSegs;
                    double kappa = afterSpline.curvature({curveIndex, ti});
                    int bin = std::distance(afterSplineInflectionsT.begin(), std::lower_bound(afterSplineInflectionsT.begin(), afterSplineInflectionsT.end(),
                                                                                              CubicBezierSpline::SplineParameter{curveIndex, ti}));
                    curvatureTotals[bin] += kappa;
                    avgKappaAfter += abs(kappa);
                    if (abs(kappa) > maxKappaAfter) {
                        maxKappaAfter = abs(kappa);
                    }
                }
            }
            avgKappaAfter /= (nSegs + 1);

            double err, symDiffErr, curvatureErr = 0;
            bool createsSmoothConnection = false;
            createsSmoothConnection = !connectsSmoothlyTo(c0, c1) && !connectsSmoothlyTo(c1, c2) && connectsSmoothlyTo(c0_, c1_);
            if ((createsSmoothConnection || maxKappaAfter < 1.1 * maxKappaBefore) && !CGAL::do_curves_intersect(afterPlE.edges_begin(), afterPlE.edges_end())) {
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
                ipeRenderer.addPainting([afterSpline, err, symDiffErr, curvatureErr, createsSmoothConnection, maxKappaBefore, maxKappaAfter, avgKappaAfter, curvatureTotals, this](renderer::GeometryRenderer& renderer) {
                    renderer.setStroke(Color(0, 0, 0), 3.0);

                    for (int curveIndex = 0; curveIndex < afterSpline.numCurves(); ++curveIndex) {
                        for (int tStep = 0; tStep <= nSegs; ++tStep) {
                            double t = static_cast<double>(tStep) / nSegs;
                            auto c = afterSpline.curvature({curveIndex, t});
                            auto n = afterSpline.normal({curveIndex, t});
                            n /= sqrt(n.squared_length());
                            auto p = afterSpline.position({curveIndex, t});
                            renderer.draw(Segment<Inexact>(p, p + 500 * c * n));
                        }
                    }

                    renderer.draw(afterSpline);
                    std::vector<CubicBezierSpline::SplinePoint> infls;
                    afterSpline.inflections(std::back_inserter(infls));
                    for (const auto& infl : infls) {
                        renderer.draw(infl.point);
                    }
                    renderer.setHorizontalTextAlignment(renderer::GeometryRenderer::HorizontalTextAlignment::AlignLeft);
                    renderer.drawText({10, 10}, "Err: " + std::to_string(err));
                    renderer.drawText({10, 20}, "Sym diff err: " + std::to_string(sqrt(symDiffErr)));
                    renderer.drawText({10, 30}, "Max. curvature before: " + std::to_string(maxKappaBefore));
                    renderer.drawText({10, 40}, "Max. curvature after: " + std::to_string(maxKappaAfter));
                    renderer.drawText({10, 50}, "#Inflections: " + std::to_string(infls.size()));
                    renderer.drawText({10, 60}, "Avg. curvature: " + std::to_string(avgKappaAfter));
                    std::stringstream ssBins;
                    ssBins << "[";
                    for (int it = 0; it < curvatureTotals.size(); ++it) {
                        ssBins << curvatureTotals[it];
                        if (it != curvatureTotals.size() - 1) ssBins << ", ";
                    }
                    ssBins << "]";
                    renderer.drawText({10, 70}, "Curvature bins: " + ssBins.str());
                    if (createsSmoothConnection) {
                        renderer.drawText({10, 50}, "Creates smooth connection!");
                    }
                }, "After");
                ipeRenderer.nextPage();
            }
        };


		// choose a point on c1 to use as p3.
		// for now just equidistant t; should this be made adaptive (todo)?
		for (int it = tSteps / 2; it < tSteps * 2 - tSteps / 2; ++it) {
            const CubicBezierCurve& curveToSample = (it < tSteps) ? c0 : ((it < 2 * tSteps) ? c1 : c2);
            int curveToSampleIndex = it / tSteps;
			double t = tSteps == 1 ? 0.5 : static_cast<double>(it % tSteps) / (tSteps - 1);
			auto p3 = curveToSample.evaluate(t);
			auto v = curveToSample.tangent(t);
            CubicBezierSpline firstPart;
            CubicBezierSpline secondPart;
            if (t < M_EPSILON) {
                if (curveToSampleIndex == 0) {
                    secondPart.appendCurve(c0);
                    secondPart.appendCurve(c1);
                    secondPart.appendCurve(c2);
                }
                if (curveToSampleIndex == 1) {
                    firstPart.appendCurve(c0);
                    secondPart.appendCurve(c1);
                    secondPart.appendCurve(c2);
                }
                if (curveToSampleIndex == 2) {
                    firstPart.appendCurve(c0);
                    firstPart.appendCurve(c1);
                    secondPart.appendCurve(c2);
                }
            } else if (t > 1 - M_EPSILON) {
                if (curveToSampleIndex == 0) {
                    firstPart.appendCurve(c0);
                    secondPart.appendCurve(c1);
                    secondPart.appendCurve(c2);
                }
                if (curveToSampleIndex == 1) {
                    firstPart.appendCurve(c0);
                    firstPart.appendCurve(c1);
                    secondPart.appendCurve(c2);
                }
                if (curveToSampleIndex == 2) {
                    firstPart.appendCurve(c0);
                    firstPart.appendCurve(c1);
                    firstPart.appendCurve(c2);
                }
            } else {
                if (curveToSampleIndex >= 1) {
                    firstPart.appendCurve(c0);
                }
                if (curveToSampleIndex == 2) {
                    firstPart.appendCurve(c1);
                }
                auto [part1, part2] = curveToSample.split(t);

                firstPart.appendCurve(part1);
                secondPart.appendCurve(part2);
                if (curveToSampleIndex == 0) {
                    secondPart.appendCurve(c1);
                }
                if (curveToSampleIndex <= 1) {
                    secondPart.appendCurve(c2);
                }
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
            if (abs(a0) < M_EPSILON) {
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
            if (abs(a1) < M_EPSILON) {
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

            evaluate(c0_, c1_);
		}

        // Try other collapse approach: fit a smooth curve
        auto otherCollapse = [&]() {
            auto threshold = 0.001;

            Vector<Inexact> v0 = (p1 - p0);
            v0 /= v0.squared_length();
            Vector<Inexact> v6 = (p6 - p5);
            v6 /= v6.squared_length();
            std::vector<Point<Inexact>> points({p0, p0 + M_EPSILON * v0, p6 - M_EPSILON * v6, p6});
            auto spline = spiro(points.begin(), points.end(), false);
            auto l = CGAL::bisector(p0, p6);
            std::vector<CubicBezierSpline::SplinePoint> inters;
            spline.intersections(l, std::back_inserter(inters));
            auto currentControl = inters.back().point;
            auto targetArea = beforeSpline.signedArea();
            auto currentArea = spline.signedArea();

            if (debug) {
                std::cout << "Target area: " << targetArea << std::endl;
                std::cout << "Initial spiro spline has area: " << currentArea << std::endl;
            }

            if (abs(targetArea - currentArea) < threshold) {
                if (debug) {
                    std::cout << "Done" << std::endl;
                }
                // done
            } else {
                // exponential search to find interval that contains the control point that achieves the target area.
                auto midpoint = CGAL::midpoint(p0, p6);
                auto v = currentControl - midpoint;
                // try currentControl + v and currentControl - v = midpoint
                auto option1 = currentControl + v;
                auto option2 = midpoint;
                auto splineAndArea = [p0, p6, v0, v6](const Point<Inexact> &control) {
                    std::vector<Point<Inexact>> controlPoints(
                            {p0, p0 + M_EPSILON * v0, control, p6 - M_EPSILON * v6, p6});
                    auto candidateSpline = spiro(controlPoints.begin(), controlPoints.end(), false);
                    return std::pair(candidateSpline, candidateSpline.signedArea());
                };
                auto [s1, a1] = splineAndArea(option1);
                auto [s2, a2] = splineAndArea(option2);

                // To ease the logic, ensure that option 1 has smaller area.
                if (a1 > a2) {
                    std::swap(a1, a2);
                    std::swap(s1, s2);
                    std::swap(option1, option2);
                }

                if (debug) {
                    std::cout << "Tried control points: " << option1 << "  and  " << option2 << std::endl;
                    std::cout << "Splines have area: " << a1 << "  and  " << a2 << std::endl;
                }

                std::tuple<Point<Inexact>, CubicBezierSpline, double> lower, upper;
                if (a1 <= targetArea && targetArea <= a2) {
                    if (debug) {
                        std::cout << "Found interval" << std::endl;
                    }
                    // found interval
                    // must lie in either [a1, currentArea) or [currentArea, a2]
                    if (targetArea < currentArea) {
                        lower = {option1, s1, a1};
                        upper = {currentControl, spline, currentArea};
                    } else {
                        lower = {currentControl, spline, currentArea};
                        upper = {option2, s2, a2};
                    }
                } else {
                    // need to search further
                    if (targetArea < a1) {
                        // search in direction of s1
                        auto candidateControl = currentControl + 2 * (option1 - currentControl);
                        auto [candidateSpline, candidateArea] = splineAndArea(candidateControl);
                        while (targetArea < candidateArea) {
                            if (debug) {
                                std::cout << "Candidate area: " <<  candidateArea << std::endl;
                            }
                            if (candidateArea >= a1) {
                                // cannot find right control point
                                return;
                            }
                            assert(candidateArea < a1); // make sure area is going in the right direction
                            currentControl = option1;
                            s1 = candidateSpline;
                            a1 = candidateArea;
                            option1 = candidateControl;
                            candidateControl = currentControl + 2 * (option1 - currentControl);
                            std::tie(candidateSpline, candidateArea) = splineAndArea(candidateControl);
                        }
                        // now targetArea >= candidateArea
                        // also, targetArea < a1
                        // so we have an interval
                        lower = {candidateControl, candidateSpline, candidateArea};
                        upper = {option1, s1, a1};
                    } else { // targetArea > a2
                        // search in direction of s2
                        auto candidateControl = currentControl + 2 * (option2 - currentControl);
                        auto [candidateSpline, candidateArea] = splineAndArea(candidateControl);
                        while (targetArea > candidateArea) {
                            if (debug) {
                                std::cout << "Candidate area: " <<  candidateArea << std::endl;
                            }
//                            assert(candidateArea > a2); // make sure area is going in the right direction
                            if (candidateArea <= a2) {
                                // cannot find right control point
                                return;
                            }
                            currentControl = option2;
                            s2 = candidateSpline;
                            a2 = candidateArea;
                            option2 = candidateControl;
                            candidateControl = currentControl + 2 * (option2 - currentControl);
                            std::tie(candidateSpline, candidateArea) = splineAndArea(candidateControl);
                        }
                        // now targetArea <= candidateArea
                        // also, targetArea > a2
                        // so we have an interval
                        lower = {option2, s2, a2};
                        upper = {candidateControl, candidateSpline, candidateArea};
                    }
                }

                int remainingIterations = 100;
                while (std::get<2>(lower) < std::get<2>(upper) - threshold && remainingIterations >= 0) {
                    if (debug) {
                        std::cout << "[" << std::get<2>(lower) << ", " << std::get<2>(upper) << "]" << std::endl;
                    }
                    // invariant: targetArea lies within the [lower, upper] interval, which has size larger than the threshold.
                    // cut the interval in twice
                    auto &[lowerC, lowerS, lowerA] = lower;
                    auto &[upperC, upperS, upperA] = upper;
//                    auto cutControl = lowerC + (targetArea - lowerA) / (upperA - lowerA) * (upperC - lowerC);
                    auto cutControl = CGAL::midpoint(lowerC, upperC);
                    auto [cutSpline, cutArea] = splineAndArea(cutControl);
                    if (targetArea < cutArea) {
                        upper = {cutControl, cutSpline, cutArea};
                    } else {
                        lower = {cutControl, cutSpline, cutArea};
                    }
                    if (abs(targetArea - cutArea) < threshold) break;
                    --remainingIterations;
                }

                auto final = abs(std::get<2>(lower) - targetArea) < abs(std::get<2>(upper) - targetArea) ?
                             lower : upper;
                if (abs(std::get<2>(final) - targetArea) > threshold) return;
                spline = std::get<1>(final);
                currentControl = std::get<0>(final);
            }

            // spline contains the spiro spline.
            // currentControl contains the spiro spline control point on the bisector line.
            // we need to approximate the spiro spline with two cubic Bézier curves, with exact tangents and exact area.
            // we will cut the spline at the currentControl point, and approximate the two parts separately.

            // cut spline
            auto cutParam = spline.nearest(currentControl).param;
            auto [spline1, spline2] = spline.split(cutParam);
            auto spline1Pl = spline1.polyline(nSegs);
            auto spline2Pl = spline2.polyline(nSegs);
            std::vector<Point<Inexact>> spline1Pts(spline1Pl.vertices_begin(), spline1Pl.vertices_end());
            std::vector<Point<Inexact>> spline2Pts(spline2Pl.vertices_begin(), spline2Pl.vertices_end());

            auto v = spline.tangent(cutParam);
            auto curve1 = fitCurve(spline1Pts, v0, -v);
            auto curve2 = fitCurve(spline2Pts, v, -v6);
            p1 = curve1.sourceControl();
            p5 = curve2.targetControl();
            auto p3 = currentControl;
            auto a0 = -spline1.signedArea();
            auto a1 = -spline2.signedArea();

            std::optional<CubicBezierCurve> c0_maybe;
            if (a0 == 0) {
                c0_maybe = CubicBezierCurve(p0, p3);
            } else {
                auto t00 = v0.direction();
                auto d00 = sqrt((p1-p0).squared_length());
                auto t01Wrong = (-v).direction();
                Direction<Inexact> t01(-t01Wrong.dx(), t01Wrong.dy());
                auto d01 = getAreaD1(a0, p0, t00, d00, t01, p3);
                if (d01 > 0) {
                    c0_maybe = createCubicBezierFromPolar(p0, t00, d00, t01, d01, p3);
                }
            }

            if (!c0_maybe.has_value()) return;

            std::optional<CubicBezierCurve> c1_maybe;
            if (a1 == 0) {
                c1_maybe = CubicBezierCurve(p3, p6);
            } else {
                auto t11Wrong = v6.direction();
                Direction<Inexact> t11(t11Wrong.dx(), -t11Wrong.dy());
                auto d11 = sqrt((p6-p5).squared_length());
                auto t10 = v.direction();
                auto d10 = getAreaD0(a1, p3, t10, t11, d11, p6);
                if (d10 > 0) {
                    c1_maybe = createCubicBezierFromPolar(p3, t10, d10, t11, d11, p6);
                }
            }

            if (!c1_maybe.has_value()) return;

            if (debug) {
                std::cout << "a0: " << a0 << " a1: " << a1 << std::endl;

                ipeRenderer.addPainting([spline, spline1, spline2, curve1, curve2, c0_maybe, c1_maybe](renderer::GeometryRenderer &renderer) {
                    renderer.setMode(renderer::GeometryRenderer::stroke);
                    renderer.setStroke(Color(0, 0, 0), 3.0);
                    renderer.draw(spline);
                    renderer.setStroke(Color(255, 0, 0), 3.0);
                    renderer.draw(spline1);
                    renderer.setStroke(Color(255, 255, 0), 3.0);
                    renderer.draw(curve1);
                    renderer.setStroke(Color(255, 0, 255), 3.0);
                    renderer.draw(*c0_maybe);

                    renderer.setStroke(Color(0, 255, 0), 3.0);
                    renderer.draw(spline2);
                    renderer.setStroke(Color(0, 255, 255), 3.0);
                    renderer.draw(curve2);
                    renderer.setStroke(Color(0, 0, 255), 3.0);
                    renderer.draw(*c1_maybe);
                });
                ipeRenderer.nextPage();
            }

            evaluate(*c0_maybe, *c1_maybe);
//        }
        };

        otherCollapse();

		if (debug) {
			ipeRenderer.save("debugging.ipe");
		}

        if (best.has_value()) {
		    edata.collapse = { .cost=minErr, .before=best->first, .after=best->second };
        }
	}
};
}
