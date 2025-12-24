#pragma once

#include <cartocrow/core/cubic_bezier.h>
#include "spiro/spiroentrypoints.h"

namespace cartocrow {
typedef struct {
    bezctx base; /* This is a superclass of bezctx, and this is the entry for the base */
    Point<Inexact> last; // last point
    CubicBezierSpline ss;
} bezctx_cbs;

static void bezctx_cbs_moveto(bezctx *z, double x, double y, int is_open);

void bezctx_cbs_lineto(bezctx *z, double x, double y);

void bezctx_cbs_quadto(bezctx *z, double xm, double ym, double x3, double y3);

void bezctx_cbs_curveto(bezctx *z, double x1, double y1, double x2, double y2, double x3, double y3);

bezctx* new_bezctx_cbs();

CubicBezierSpline bezctx_cbs_close(bezctx *z);
}