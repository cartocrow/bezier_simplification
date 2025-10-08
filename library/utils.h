#pragma once

#include <cartocrow/core/core.h>

namespace cartocrow::curved_simplification::utils {
template<typename T>
bool vectorRemove(T elt, std::vector<T>& vec) {
	auto pos = std::find(vec.begin(), vec.end(), elt);
	if (pos != vec.end()) {
		vec.erase(pos);
		return true;
	}
	else {
		return false;
	}
}
}
