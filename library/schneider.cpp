#include "schneider.h"
#include <cmath>
#include <cassert>

namespace cartocrow::curved_simplification {
inline double length(const Vector<Inexact>& v) { return std::sqrt(v.squared_length()); }
inline Vector<Inexact> normalize(const Vector<Inexact>& v) {
    double len = length(v);
    if (len == 0.0) return Vector<Inexact>(0.0, 0.0); // assume Vector ctor takes two doubles
    return v / len;
}

// Bernstein basis for cubic
inline double B0(double u) { double t = 1.0 - u; return t*t*t; }
inline double B1(double u) { double t = 1.0 - u; return 3.0 * u * t * t; }
inline double B2(double u) { double t = 1.0 - u; return 3.0 * u * u * t; }
inline double B3(double u) { return u*u*u; }

// Parameterize points by chord length
std::vector<double> chordLengthParameterize(const std::vector<Point<Inexact>>& pts, int first, int last) {
    int n = last - first + 1;
    std::vector<double> u(n);
    u[0] = 0.0;
    for (int i=first+1;i<=last;++i) {
        Vector<Inexact> diff = pts[i] - pts[i-1];
        u[i-first] = u[i-first-1] + length(diff);
    }
    double denom = u[n-1];
    if (denom == 0.0) {
        for (int i=0;i<n;++i) u[i] = 0.0;
        return u;
    }
    for (int i=1;i<n;++i) u[i] /= denom;
    return u;
}

// Compute unit tangents
Vector<Inexact> computeLeftTangent(const std::vector<Point<Inexact>>& pts, int index) {
    Vector<Inexact> t = pts[index+1] - pts[index];
    return normalize(t);
}
Vector<Inexact> computeRightTangent(const std::vector<Point<Inexact>>& pts, int index) {
    Vector<Inexact> t = pts[index-1] - pts[index];
    return normalize(t);
}
Vector<Inexact> computeCenterTangent(const std::vector<Point<Inexact>>& pts, int center) {
    Vector<Inexact> v1 = pts[center-1] - pts[center];
    Vector<Inexact> v2 = pts[center] - pts[center+1];
    Vector<Inexact> t((v1.x() + v2.x())/2.0, (v1.y() + v2.y())/2.0);
    return normalize(t);
}

// Compute max squared error and split point (returns maxSquaredDist)
double computeMaxError(const std::vector<Point<Inexact>>& pts, int first, int last,
                       const CubicBezierCurve& bez, const std::vector<double>& u, int& splitPoint) {
    splitPoint = first + (last - first)/2;
    double maxDist = 0.0;
    for (int i=first+1;i<last;++i) {
        Point<Inexact> P = bez.evaluate(u[i-first]);
        Vector<Inexact> diff = Vector<Inexact>(P.x() - pts[i].x(), P.y() - pts[i].y());
        double d2 = squared_length(diff);
        if (d2 >= maxDist) {
            maxDist = d2;
            splitPoint = i;
        }
    }
    return maxDist;
}

// Newton-Raphson root finding for parameter refinement
double newtonRaphsonRootFind(const CubicBezierCurve& Q, const Point<Inexact>& P, double u) {
    Point<Inexact> Q_u = Q.evaluate(u);
    Vector<Inexact> Q1_u = Q.derivative(u);
    Vector<Inexact> Q2_u = Q.derivative2(u);

    double numerator = (Q_u.x() - P.x()) * Q1_u.x() + (Q_u.y() - P.y()) * Q1_u.y();
    double denominator = Q1_u.x()*Q1_u.x() + Q1_u.y()*Q1_u.y()
                       + (Q_u.x() - P.x()) * Q2_u.x() + (Q_u.y() - P.y()) * Q2_u.y();
    if (denominator == 0.0) return u;
    double uPrime = u - numerator / denominator;
    if (uPrime < 0.0) uPrime = 0.0;
    else if (uPrime > 1.0) uPrime = 1.0;
    return uPrime;
}

// Reparameterize
std::vector<double> reparameterize(const std::vector<Point<Inexact>>& pts, int first, int last,
                                   const std::vector<double>& u, const CubicBezierCurve& bez) {
    int n = last - first + 1;
    std::vector<double> uPrime(n);
    for (int i=first;i<=last;++i) uPrime[i-first] = newtonRaphsonRootFind(bez, pts[i], u[i-first]);
    return uPrime;
}

// Generate Bezier control points by least-squares (Schneider)
CubicBezierCurve generateBezier(const std::vector<Point<Inexact>>& pts, int first, int last,
                                    const std::vector<double>& uPrime,
                                    const Vector<Inexact>& tHat1, const Vector<Inexact>& tHat2) {
    int nPts = last - first + 1;
    if (nPts == 3) {
        auto& a = pts[first];
        auto& b = pts[first+1];
        auto& c = pts[last];
        std::vector<Point<Inexact>> newPts({a, CGAL::midpoint(a, b), b, CGAL::midpoint(b, c), c});
        auto& ua = uPrime[first];
        auto& ub = uPrime[first+1];
        auto& uc = uPrime[last];
        std::vector<double> uAdapted({ua, (ua + ub) / 2, ub, (ub + uc) / 2, uc});
        return generateBezier(newPts, 0, 4, uAdapted, tHat1, tHat2);
    }
    std::vector<std::array<Vector<Inexact>,2>> A(nPts);
    for (int i=0;i<nPts;++i) {
        A[i][0] = tHat1 * B1(uPrime[i]);
        A[i][1] = tHat2 * B2(uPrime[i]);
    }
    double C00 = 0, C01 = 0, C11 = 0;
    double X0 = 0, X1 = 0;
    for (int i=0;i<nPts;++i) {
        C00 += A[i][0] * A[i][0];
        C01 += A[i][0] * A[i][1];
        C11 += A[i][1] * A[i][1];
        // tmp = P_i - (P0 * B0 + P0 * B1 + Pn * B2 + Pn * B3)
        double b0 = B0(uPrime[i]);
        double b1 = B1(uPrime[i]);
        double b2 = B2(uPrime[i]);
        double b3 = B3(uPrime[i]);
        Point<Inexact> P0 = pts[first];
        Point<Inexact> Pn = pts[last];
        Point<Inexact> tmp( pts[first + i].x() - (P0.x()*b0 + P0.x()*b1 + Pn.x()*b2 + Pn.x()*b3),
                    pts[first + i].y() - (P0.y()*b0 + P0.y()*b1 + Pn.y()*b2 + Pn.y()*b3));
        X0 += A[i][0] * Vector<Inexact>(tmp.x(), tmp.y());
        X1 += A[i][1] * Vector<Inexact>(tmp.x(), tmp.y());
    }
    double det = C00 * C11 - C01 * C01;
    double alpha_l = 0.0, alpha_r = 0.0;
    if (det != 0.0) {
        alpha_l = (X0 * C11 - X1 * C01) / det;
        alpha_r = (C00 * X1 - C01 * X0) / det;
    }
    double segLength = length(pts[last] - pts[first]);
    double epsilon = 1.0e-6 * segLength;
    auto& p0 = pts[first];
    auto& p3 = pts[last];
    if (alpha_l < epsilon || alpha_r < epsilon) {
        double dist = segLength / 3.0;
        Point<Inexact> p1( p0.x() + tHat1.x() * dist, p0.y() + tHat1.y() * dist );
        Point<Inexact> p2( p3.x() + tHat2.x() * dist, p3.y() + tHat2.y() * dist );
        return {p0, p1, p2, p3};
    }
    Point<Inexact> p1( p0.x() + tHat1.x() * alpha_l, p0.y() + tHat1.y() * alpha_l );
    Point<Inexact> p2( p3.x() + tHat2.x() * alpha_r, p3.y() + tHat2.y() * alpha_r );
    return {p0, p1, p2, p3};
}

// The recursive fitter
void fitCubicRecursive(const std::vector<Point<Inexact>>& pts, int first, int last,
                       const Vector<Inexact>& tHat1, const Vector<Inexact>& tHat2, double allowedError,
                       CubicBezierSpline& outSpline, int maxRecursion=10, int reparameterizationIterations=10) {
    int nPts = last - first + 1;
    if (nPts == 2) {
        double dist = length(pts[last] - pts[first]) / 3.0;
        Point<Inexact> p1( pts[first].x() + tHat1.x()*dist, pts[first].y() + tHat1.y()*dist );
        Point<Inexact> p2( pts[last].x() + tHat2.x()*dist, pts[last].y() + tHat2.y()*dist );
        outSpline.appendCurve(pts[first], p1, p2, pts[last]);
        return;
    }
    auto u = chordLengthParameterize(pts, first, last);
    auto bez = generateBezier(pts, first, last, u, tHat1, tHat2);
    int splitPoint = 0;
    double maxError = computeMaxError(pts, first, last, bez, u, splitPoint);
    if (maxError < allowedError) {
        outSpline.appendCurve(bez);
        return;
    }
    for (int i=0;i<reparameterizationIterations;++i) {
        auto uPrime = reparameterize(pts, first, last, u, bez);
        bez = generateBezier(pts, first, last, uPrime, tHat1, tHat2);
        maxError = computeMaxError(pts, first, last, bez, uPrime, splitPoint);
        if (maxError < allowedError) {
            outSpline.appendCurve(bez);
            return;
        }
        u = std::move(uPrime);
    }
    if (maxRecursion == 0) {
        outSpline.appendCurve(bez);
        return;
    }
    // Split and fit recursively
    Vector<Inexact> tHatCenter = computeCenterTangent(pts, splitPoint);
    fitCubicRecursive(pts, first, splitPoint, tHat1, tHatCenter, allowedError, outSpline, maxRecursion - 1);
    Vector<Inexact> negCenter(-tHatCenter.x(), -tHatCenter.y());
    fitCubicRecursive(pts, splitPoint, last, negCenter, tHat2, allowedError, outSpline, maxRecursion - 1);
}

CubicBezierSpline fitSpline(const std::vector<Point<Inexact>>& points, double maxSquaredError, int reparameterizationIterations) {
    assert(points.size() >= 2);
    CubicBezierSpline spline;
    Vector<Inexact> tHat1 = computeLeftTangent(points, 0);
    Vector<Inexact> tHat2 = computeRightTangent(points, static_cast<int>(points.size()) - 1);
    fitCubicRecursive(points, 0, static_cast<int>(points.size()) - 1, tHat1, tHat2, maxSquaredError, spline, 10, reparameterizationIterations);
    return spline;
}

CubicBezierSpline fitSpline(const std::vector<Point<Inexact>>& points, double maxSquaredError, const Vector<Inexact>& startTangent, const Vector<Inexact>& endTangent, int reparameterizationIterations) {
    assert(points.size() >= 2);
    CubicBezierSpline spline;
    fitCubicRecursive(points, 0, static_cast<int>(points.size()) - 1, normalize(startTangent), normalize(endTangent), maxSquaredError, spline, 10, reparameterizationIterations);
    return spline;
}

CubicBezierCurve fitCurve(const std::vector<Point<Inexact>>& points, int reparameterizationIterations) {
    assert(points.size() >= 2);
    if (points.size() == 2) return CubicBezierCurve(points[0], points[1]);
    CubicBezierSpline spline;
    Vector<Inexact> tHat1 = computeLeftTangent(points, 0);
    Vector<Inexact> tHat2 = computeRightTangent(points, static_cast<int>(points.size()) - 1);
    fitCubicRecursive(points, 0, static_cast<int>(points.size()) - 1, tHat1, tHat2, 0, spline, 0, reparameterizationIterations);
    return spline.curve(0);
}

CubicBezierCurve fitCurve(const std::vector<Point<Inexact>>& points, const Vector<Inexact>& startTangent, const Vector<Inexact>& endTangent, int reparameterizationIterations) {
    assert(points.size() >= 2);
    CubicBezierSpline spline;
    fitCubicRecursive(points, 0, static_cast<int>(points.size()) - 1, normalize(startTangent), normalize(endTangent), 0.001, spline, 0, reparameterizationIterations);
    return spline.curve(0);
}

CubicBezierSpline fitTwoCurves(const std::vector<Point<Inexact>>& points, int reparameterizationIterations) {
    assert(points.size() >= 2);
    CubicBezierSpline spline;
    Vector<Inexact> tHat1 = computeLeftTangent(points, 0);
    Vector<Inexact> tHat2 = computeRightTangent(points, static_cast<int>(points.size()) - 1);
    fitCubicRecursive(points, 0, static_cast<int>(points.size()) - 1, tHat1, tHat2, 0, spline, 1, reparameterizationIterations);
    return spline;
}
}
