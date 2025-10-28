#include "steven_bezier_collapse.h"
#include <cartocrow/core/vector_helpers.h>

namespace cartocrow::curved_simplification {
Number<Inexact> getAreaD0(Number<Inexact> a, Point<Inexact> p0, Direction<Inexact> t0, Direction<Inexact> t1, Number<Inexact> d1, Point<Inexact> p3) {
	CGAL::Aff_transformation_2<Inexact> move(CGAL::TRANSLATION, CGAL::ORIGIN - p0);
	double hyp = sqrt(CGAL::squared_distance(p0, p3));
	CGAL::Aff_transformation_2<Inexact> rotate(CGAL::ROTATION, -(p3.y() - p0.y()) / hyp, (p3.x() - p0.x()) / hyp);
	CGAL::Aff_transformation_2<Inexact> scale(CGAL::SCALING, 1.0 / hyp);

	auto trans = scale * rotate * move;

	Point<Inexact> t0p(p0 + t0.vector());
	auto t0p_ = t0p.transform(trans);
	auto t0_ = (t0p_ - CGAL::ORIGIN).direction();
	Point<Inexact> t1p(p3.x() - t1.dx(), p3.y() + t1.dy());
	auto t1p_ = t1p.transform(trans);
	auto t1_ = Direction<Inexact>(1 - t1p_.x(), t1p_.y());
	auto d1_ = d1 / hyp;
	double a_ = a / (hyp * hyp);

	auto t0hyp = sqrt(t0_.dx() * t0_.dx() + t0_.dy() * t0_.dy());
	auto cost0 = t0_.dx() / t0hyp;
	auto sint0 = t0_.dy() / t0hyp;
	auto t1hyp = sqrt(t1_.dx() * t1_.dx() + t1_.dy() * t1_.dy());
	auto cost1 = t1_.dx() / t1hyp;
	auto sint1 = t1_.dy() / t1hyp;
	// Use sine of sum rule: sin(theta0 + theta1) = sin(theta0)cos(theta1) + cos(theta0)sin(theta1)
	auto sint0t1 = sint0*cost1 + cost0*sint1;

	double d0 = ((20.0/3.0)*(a_ - (3.0/10.0) * d1_ * sint1)) / (2 * sint0 - d1_ * sint0t1);
	return d0 * hyp;
}

Number<Inexact> getAreaD1(Number<Inexact> a, Point<Inexact> p0, Direction<Inexact> t0, Number<Inexact> d0, Direction<Inexact> t1, Point<Inexact> p3) {
	CGAL::Aff_transformation_2<Inexact> move(CGAL::TRANSLATION, CGAL::ORIGIN - p0);
	double hyp = sqrt(CGAL::squared_distance(p0, p3));
	CGAL::Aff_transformation_2<Inexact> rotate(CGAL::ROTATION, -(p3.y() - p0.y()) / hyp, (p3.x() - p0.x()) / hyp);
	CGAL::Aff_transformation_2<Inexact> scale(CGAL::SCALING, 1.0 / hyp);

	auto trans = scale * rotate * move;

	Point<Inexact> t0p(p0 + t0.vector());
	auto t0p_ = t0p.transform(trans);
	auto t0_ = (t0p_ - CGAL::ORIGIN).direction();
	Point<Inexact> t1p(p3.x() - t1.dx(), p3.y() + t1.dy());
	auto t1p_ = t1p.transform(trans);
	auto t1_ = Direction<Inexact>(1 - t1p_.x(), t1p_.y());
	auto d0_ = d0 / hyp;
	double a_ = a / (hyp * hyp);

	auto t0hyp = sqrt(t0_.dx() * t0_.dx() + t0_.dy() * t0_.dy());
	auto cost0 = t0_.dx() / t0hyp;
	auto sint0 = t0_.dy() / t0hyp;
	auto t1hyp = sqrt(t1_.dx() * t1_.dx() + t1_.dy() * t1_.dy());
	auto cost1 = t1_.dx() / t1hyp;
	auto sint1 = t1_.dy() / t1hyp;
	// Use sine of sum rule: sin(theta0 + theta1) = sin(theta0)cos(theta1) + cos(theta0)sin(theta1)
	auto sint0t1 = sint0*cost1 + cost0*sint1;

	double d1 = ((20.0/3.0)*(a_ - (3.0/10.0) * d0_ * sint0)) / (2 * sint1 - d0_ * sint0t1);
	return d1 * hyp;
}

CubicBezierCurve
createCubicBezierFromPolar(Point<Inexact> p0, Direction<Inexact> t0Dir, Number<Inexact> d0,
						   Direction<Inexact> t1Dir, Number<Inexact> d1, Point<Inexact> p3) {
	auto t0v = t0Dir.to_vector();
	auto t0 = t0v / sqrt(t0v.squared_length());
	auto t1v = t1Dir.to_vector();
	auto t1 = t1v / sqrt(t1v.squared_length());
	auto p1 = p0 + t0 * d0;
	Point<Inexact> p2 = { p3.x() - d1 * t1.x(), p3.y() + d1 * t1.y() };
	return {p0, p1, p2, p3};
}

Number<Inexact> length(const Polyline<Inexact>& pl) {
	double total = 0;
	for (auto eit = pl.edges_begin(); eit != pl.edges_end(); ++eit) {
		total += sqrt(eit->squared_length());
	}
	return total;
}

bool isStraight(const CubicBezierCurve& curve) {
    auto a1 = abs(CGAL::area(curve.control(0), curve.control(1), curve.control(2)));
    auto a2 = abs(CGAL::area(curve.control(1), curve.control(2), curve.control(3)));
    return a1 + a2 < M_EPSILON;
}

bool connectsSmoothlyTo(const CubicBezierCurve& c1, const CubicBezierCurve& c2) {
    auto control1 = CGAL::squared_distance(c1.targetControl(), c1.target()) > M_EPSILON ? c1.targetControl() :
                    (CGAL::squared_distance(c1.sourceControl(), c1.target()) > M_EPSILON ? c1.sourceControl() : c1.source());
    auto control2 = CGAL::squared_distance(c2.source(), c2.sourceControl()) > M_EPSILON ? c2.sourceControl() :
                    (CGAL::squared_distance(c2.source(), c2.targetControl()) > M_EPSILON ? c2.targetControl() : c2.target());
    return smallestAngleBetween(c1.target() - control1, control2 - c1.target()) < M_EPSILON;
}
}