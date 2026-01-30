#pragma once

#include "bezier_collapse.h"
#include "schneider.h"
#include "spiro_helpers/spiro_helpers.h"

#include <cartocrow/core/arrangement_helpers.h>

#define DEBUG 0

#if DEBUG
#include <cartocrow/renderer/ipe_renderer.h>
#endif

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
    int symDiffSegs;
	bool debug;
#if DEBUG
    renderer::IpeRenderer ipeRenderer;
#endif

	StevenBCTraits(bool debug = false, int tSteps = 10,	int nSegs = 100, int symDiffSegs = 10) :
	      tSteps(tSteps), nSegs(nSegs), symDiffSegs(symDiffSegs), debug(debug) {};

private:
    /// Returns two cubic Bézier curves that approximate a spiro spline that in turn approximates the input spline.
    /// The returned curves match the tangents and signed area of the input.
    /// Precondition: the input spline is G^0 continuous and consists of three cubic Bézier curves.
    std::optional<std::pair<CubicBezierCurve, CubicBezierCurve>>
    findSpiro(const CubicBezierSpline& beforeSpline) {
        auto p0 = beforeSpline.controlPoint(0);
        auto p1 = beforeSpline.controlPoint(1);
        auto p5 = beforeSpline.controlPoint(8);
        auto p6 = beforeSpline.controlPoint(9);

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
        if (inters.empty()) {
            std::cerr << "Empty intersection with line in spiro; should be impossible!" << std::endl;
            std::cerr << "Line: " << l << std::endl;
            std::cerr << "Spiro control points; start: " << p0 << " with tangent " << v0 << "  to target: " << p6 << " with tangent: " << v6 << std::endl;
            return std::nullopt;
        }
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
                    int remainingIterations = 20;
                    while (targetArea < candidateArea && remainingIterations > 0) {
                        if (debug) {
                            std::cout << "Candidate area: " <<  candidateArea << std::endl;
                        }
                        if (candidateArea >= a1) {
                            // cannot find right control point
                            return std::nullopt;
                        }
                        assert(candidateArea < a1); // make sure area is going in the right direction
                        currentControl = option1;
                        s1 = candidateSpline;
                        a1 = candidateArea;
                        option1 = candidateControl;
                        candidateControl = currentControl + 2 * (option1 - currentControl);
                        std::tie(candidateSpline, candidateArea) = splineAndArea(candidateControl);
                        --remainingIterations;
                    }
                    if (targetArea < candidateArea) {
                        return std::nullopt;
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
                    int remainingIterations = 20;
                    while (targetArea > candidateArea && remainingIterations > 0) {
                        if (debug) {
                            std::cout << "Candidate area: " <<  candidateArea << std::endl;
                        }
                        if (candidateArea <= a2) {
                            // cannot find right control point
                            return std::nullopt;
                        }
                        currentControl = option2;
                        s2 = candidateSpline;
                        a2 = candidateArea;
                        option2 = candidateControl;
                        candidateControl = currentControl + 2 * (option2 - currentControl);
                        std::tie(candidateSpline, candidateArea) = splineAndArea(candidateControl);
                        --remainingIterations;
                    }
                    if (targetArea > candidateArea) {
                        return std::nullopt;
                    }
                    // now targetArea <= candidateArea
                    // also, targetArea > a2
                    // so we have an interval
                    lower = {option2, s2, a2};
                    upper = {candidateControl, candidateSpline, candidateArea};
                }
            }

            int remainingIterations = 10;
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
//                if (abs(std::get<2>(final) - targetArea) > threshold) return;
            spline = std::get<1>(final);
            currentControl = std::get<0>(final);
        }

        // spline contains the spiro spline.
        // currentControl contains the spiro spline control point on the bisector line.
        // we need to approximate the spiro spline with two cubic Bézier curves, with exact tangents and exact area.
        // we will cut the spline at the currentControl point, and approximate the two parts separately.

        // cut spline
        auto cutParam = spline.nearest(currentControl).param;
        Polygon<Inexact> triangle;
        triangle.push_back(spline.source());
        triangle.push_back(currentControl);
        triangle.push_back(spline.target());
        double willBeAdded = triangle.area();
        auto [spline1, spline2] = spline.split(cutParam);
        if (spline1.empty() || spline2.empty()) {
            std::cout << "Spline empty!" << std::endl;
            return std::nullopt;
        }
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
        auto s1a = spline1.signedArea();
        auto s2a = spline2.signedArea();
        auto theTargetSum = targetArea - willBeAdded;
        auto a0 = -s1a * theTargetSum / (s1a + s2a);
//        auto a1 = -spline2.signedArea();
        auto a1 = -s2a * theTargetSum / (s1a + s2a);

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

        if (!c0_maybe.has_value()) return std::nullopt;

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

        if (!c1_maybe.has_value()) return std::nullopt;

#if DEBUG
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
#endif

        return std::pair(*c0_maybe, *c1_maybe);
    }

    /// Returns two cubic Bézier curves that approximate curves c0, c1, and c2 by splitting the spline at parameter param.
    /// The returned curves match the tangents and signed area of the input.
    std::optional<std::pair<CubicBezierCurve, CubicBezierCurve>>
    closeFit(const CubicBezierSpline& beforeSpline, const CubicBezierSpline::SplineParameter& param) {
        auto p0 = beforeSpline.controlPoint(0);
        auto p1 = beforeSpline.controlPoint(1);
        auto p5 = beforeSpline.controlPoint(8);
        auto p6 = beforeSpline.controlPoint(9);
        const auto& c0 = beforeSpline.curve(0);
        const auto& c1 = beforeSpline.curve(1);
        const auto& c2 = beforeSpline.curve(2);

        auto p3 = beforeSpline.evaluate(param);
        auto v = beforeSpline.tangent(param);

        auto [firstPart, secondPart] = beforeSpline.split(param);
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
            if (d01 > 0 && isfinite(d01)) {
                c0_maybe = createCubicBezierFromPolar(p0, t00, d00, t01, d01, p3);
            }
        }

        if (!c0_maybe.has_value()) return std::nullopt;

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
            if (d10 > 0 && isfinite(d10)) {
                c1_maybe = createCubicBezierFromPolar(p3, t10, d10, t11, d11, p6);
            }
        }

        if (!c1_maybe.has_value()) return std::nullopt;

        return std::pair(*c0_maybe, *c1_maybe);
    }

    struct CurvatureInfo {
        std::vector<double> curvatureTotals;
        double maximumAbsoluteCurvature;
        double averageAbsoluteCurvature;

    };
    CurvatureInfo evaluateCurvature(CubicBezierSpline spline) {
        double maxKappa = 0;
        double avgKappa = 0;

        std::vector<CubicBezierSpline::SplineParameter> splineInflectionsT;
        spline.inflectionsT(std::back_inserter(splineInflectionsT));
        std::vector<double> curvatureTotals(splineInflectionsT.size() + 1);
        for (int curveIndex = 0; curveIndex < spline.numCurves(); ++curveIndex) {
            for (int i = 0; i <= nSegs; ++i) {
                double ti = static_cast<double>(i) / nSegs;
                double kappa = spline.curvature({curveIndex, ti});
                int bin = std::distance(splineInflectionsT.begin(), std::lower_bound(splineInflectionsT.begin(), splineInflectionsT.end(),
                                                                                     CubicBezierSpline::SplineParameter{curveIndex, ti}));
                curvatureTotals[bin] += kappa;
                avgKappa += abs(kappa);
                if (abs(kappa) > maxKappa) {
                    maxKappa = abs(kappa);
                }
            }
        }
        avgKappa /= (nSegs + 1);

        return {curvatureTotals, maxKappa, avgKappa};
    }

    double evaluateSymDiff(const CubicBezierSpline& beforeSpline, const CubicBezierSpline& afterSpline) {
        auto beforePlE = pretendExact(beforeSpline.polyline(symDiffSegs));
        auto afterPlE = pretendExact(afterSpline.polyline(symDiffSegs));
        Arrangement<Exact> arr;
        std::vector<Arrangement<Exact>::X_monotone_curve_2> beforePlXMCurves;
        for (auto eit = beforePlE.edges_begin(); eit != beforePlE.edges_end(); ++eit) {
            beforePlXMCurves.emplace_back(*eit);
        }
        CGAL::insert_non_intersecting_curves(arr, beforePlXMCurves.begin(), beforePlXMCurves.end());
        CGAL::insert(arr, afterPlE.edges_begin(), afterPlE.edges_end());
#if DEBUG
        if (debug) {
            ipeRenderer.addPainting([beforeSpline, afterSpline, arr](renderer::GeometryRenderer& renderer) {
                renderer.setMode(renderer::GeometryRenderer::stroke);
                renderer.draw(beforeSpline);
                renderer.draw(afterSpline);
                for (auto eit = arr.edges_begin(); eit != arr.edges_end(); ++eit) {
                    Segment<Exact> seg = eit->curve();
                    renderer.draw(seg);
                }
            });
            ipeRenderer.nextPage();
        }
#endif
        double symDiffErr = 0;
        for (auto fit = arr.faces_begin(); fit != arr.faces_end(); ++fit) {
            if (fit->is_unbounded()) continue;
            auto pwh = approximate(face_to_polygon_with_holes<Exact>(fit));
            symDiffErr += abs(pwh.outer_boundary().area());
            for (const auto& h : pwh.holes()) {
                symDiffErr -= abs(h.area());
            }
        }
        return symDiffErr;
    }

public:
    void determineCollapse(typename BG::Edge_handle e) {
#if DEBUG
        if (debug) {
            ipeRenderer = renderer::IpeRenderer();
        }
#endif

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

        // Closed spline that consists of three Béziers.
        // Current approaches do not work well.
        // For now, just do not simplify these further.
        if (c0.source() == c2.target()) {
            return;
        }

        double minErr = std::numeric_limits<double>::infinity();
        std::optional<std::pair<CubicBezierCurve, CubicBezierCurve>> best = std::nullopt;

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


        auto evaluate = [&](const CubicBezierCurve& c0_, const CubicBezierCurve& c1_) {
            CubicBezierSpline afterSpline;
            afterSpline.appendCurve(c0_);
            afterSpline.appendCurve(c1_);
            auto afterPl = afterSpline.polyline(10);

            auto curvatureInfo = evaluateCurvature(afterSpline);
            auto maxKappaAfter = curvatureInfo.maximumAbsoluteCurvature;
            auto avgKappaAfter = curvatureInfo.averageAbsoluteCurvature;
            auto curvatureTotals = curvatureInfo.curvatureTotals;

            double err = std::numeric_limits<double>::infinity();
            double symDiffErr = std::numeric_limits<double>::infinity();
            bool createsSmoothConnection = (!connectsSmoothlyTo(c0, c1) || !connectsSmoothlyTo(c1, c2)) && connectsSmoothlyTo(c0_, c1_);
            if (debug) {
                std::cout << "c0 connects smoothly to c1: " << connectsSmoothlyTo(c0, c1) << std::endl;
                std::cout << "c1 connects smoothly to c2: " << connectsSmoothlyTo(c1, c2) << std::endl;
                std::cout << "c0_ connects smoothly to c1_: " << connectsSmoothlyTo(c0_, c1_) << std::endl;
            }
            if ((createsSmoothConnection || maxKappaAfter < 1.1 * maxKappaBefore) && !afterSpline.selfIntersects(0.01)) {
                symDiffErr = evaluateSymDiff(beforeSpline, afterSpline);
                err = sqrt(symDiffErr);

                if (createsSmoothConnection) {
                    err /= 10;
                }
            }

            if (err < minErr) {
                minErr = err;
                best = {c0_, c1_};
            }

#if DEBUG
            if (debug) {
                ipeRenderer.addPainting([beforeSpline](renderer::GeometryRenderer& renderer) {
                    renderer.setMode(renderer::GeometryRenderer::vertices | renderer::GeometryRenderer::stroke);
                    renderer.setStroke(Color(200, 200, 200), 3.0);
                    renderer.draw(beforeSpline);
                }, "Before");
                ipeRenderer.addPainting([afterSpline, err, symDiffErr, createsSmoothConnection, maxKappaBefore, maxKappaAfter, avgKappaAfter, curvatureTotals, this](renderer::GeometryRenderer& renderer) {
                    renderer.setStroke(Color(200, 0, 200), 0.1);
                    renderer.setStrokeOpacity(100);

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

                    renderer.setStroke(Color(0, 0, 0), 3.0);
                    renderer.setStrokeOpacity(255);
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
                        renderer.drawText({10, 80}, "Creates smooth connection!");
                    }
                }, "After");
                ipeRenderer.nextPage();
            }
#endif
        };


        // choose a point on c1 to use as p3.
        // for now just equidistant t
        for (int it = tSteps / 2; it < tSteps * 2 - tSteps / 2; ++it) {
            int curveToSampleIndex = it / tSteps;
            double t = tSteps == 1 ? 0.5 : static_cast<double>(it % tSteps) / (tSteps - 1);
            auto closeFitResult = closeFit(beforeSpline, {curveToSampleIndex, t});
            if (closeFitResult.has_value()) {
                evaluate(closeFitResult->first, closeFitResult->second);
            }
        }

        // Try other collapse approach: fit a smooth curve
        auto spiroResult = findSpiro(beforeSpline);
        if (spiroResult.has_value()) {
            evaluate(spiroResult->first, spiroResult->second);
        }

#if DEBUG
        if (debug) {
            ipeRenderer.save("debugging.ipe");
        }
#endif

        if (best.has_value()) {
            edata.collapse = { .cost=minErr, .before=best->first, .after=best->second };
        }
    }
};
}
