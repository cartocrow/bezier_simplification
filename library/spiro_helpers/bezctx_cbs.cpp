#include "bezctx_cbs.h"

namespace cartocrow {
static void
bezctx_cbs_moveto(bezctx *z, double x, double y, int is_open) {
    auto *bc = (bezctx_cbs *)z;
    bc->last = Point<Inexact>(x, y);
}

void
bezctx_cbs_lineto(bezctx *z, double x, double y) {
    auto *bc = (bezctx_cbs *)z;
    Point<Inexact> nextPoint(x, y);
    bc->ss.appendCurve(bc->last, nextPoint);
    bc->last = nextPoint;
}

void
bezctx_cbs_quadto(bezctx *z, double xm, double ym, double x3, double y3)
{
    auto *bc = (bezctx_cbs *)z;
    Point<Inexact> m(xm, ym);
    Point<Inexact> t(x3, y3);
    bc->ss.appendCurve(bc->last, m, t);
    bc->last = t;
}

void
bezctx_cbs_curveto(bezctx *z, double x1, double y1, double x2, double y2,
                   double x3, double y3)
{
    auto *bc = (bezctx_cbs *)z;
    Point<Inexact> p1(x1, y1);
    Point<Inexact> p2(x2, y2);
    Point<Inexact> p3(x3, y3);
    bc->ss.appendCurve(bc->last, p1, p2, p3);
    bc->last = p3;
}

bezctx *
new_bezctx_cbs(void) {
    auto* result = new bezctx_cbs;

    result->base.moveto = bezctx_cbs_moveto;
    result->base.lineto = bezctx_cbs_lineto;
    result->base.quadto = bezctx_cbs_quadto;
    result->base.curveto = bezctx_cbs_curveto;
    result->base.mark_knot = nullptr;
    return &result->base;
}

CubicBezierSpline
bezctx_cbs_close(bezctx *z)
{
    auto *bc = (bezctx_cbs *)z;
    auto spline = bc->ss;
    delete bc;
    return spline;
}
}
