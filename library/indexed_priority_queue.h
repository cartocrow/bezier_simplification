#pragma once

#include <vector>

namespace cartocrow::curved_simplification {
	template <class QT> concept QueueTraits = requires(typename QT::EltHandle elt, int i) {
		typename QT::EltHandle;

		{ QT::setIndex(elt, i) };

		{
			QT::getIndex(elt)
		} -> std::same_as<int>;

		{
			QT::compare(elt, elt)
		} -> std::same_as<int>; // negative if elt < elt2, positive if elt > elt2, zero if elt = elt2. Smallest value == highest priority (top of queue)
	};

	template <QueueTraits QT> class IndexedPriorityQueue {
	public:
		using EltHandle = QT::EltHandle;

	private:
		std::vector<EltHandle> queue;

		void siftUp(int k, EltHandle elt);
		void siftDown(int k, EltHandle elt);

	public:
		bool empty();

		void push(EltHandle elt);
		EltHandle pop();

		bool remove(EltHandle elt);
		bool contains(EltHandle elt);
		void update(EltHandle elt);
	};

} // namespace cartocrow::simplification

#include "indexed_priority_queue.hpp"