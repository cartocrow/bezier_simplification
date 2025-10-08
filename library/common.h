#pragma once

#include <cartocrow/core/core.h>

namespace cartocrow::curved_simplification {

template<class Elt, typename Kernel>
struct GraphQueueTraits {
	using EltHandle = Elt;

	static void setIndex(Elt elt, int index) {
		elt->data().qid = index;
	}

	static int getIndex(Elt elt) {
		return elt->data().qid;
	}

	static int compare(Elt a, Elt b) {
        const auto& aclps = a->data().collapse;
        const auto& bclps = b->data().collapse;
		Number<Kernel> ac = aclps.has_value() ? aclps->cost : std::numeric_limits<double>::infinity();
		Number<Kernel> bc = bclps.has_value() ? bclps->cost : std::numeric_limits<double>::infinity();
		if (ac < bc) {
			return -1;
		}
		else if (ac > bc) {
			return 1;
		}
		else {
			return 0;
		}
	}
};
}