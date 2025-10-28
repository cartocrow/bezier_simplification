// -----------------------------------------------------------------------------
// IMPLEMENTATION OF TEMPLATE FUNCTIONS
// Do not include this file, but the .h file instead
// -----------------------------------------------------------------------------
 
namespace cartocrow::curved_simplification {

	template <QueueTraits QT>
	void IndexedPriorityQueue<QT>::siftUp(int k, EltHandle elt) {
		while (k > 0) {
			int parent = (k - 1) >> 1;
			EltHandle e = queue[parent];
			if (QT::compare(elt, e) >= 0) {
				break;
			}
			queue[k] = e;
			QT::setIndex(e, k);
			k = parent;
		}
		queue[k] = elt;
		QT::setIndex(elt, k);
	}

	template <QueueTraits QT>
	void IndexedPriorityQueue<QT>::siftDown(int k, EltHandle elt) {
		int half = queue.size() >> 1;
		while (k < half) {
			int child = (k << 1) + 1;
			EltHandle c = queue[child];
			int right = child + 1;
			if (right < queue.size() && QT::compare(c, queue[right]) > 0) {
				c = queue[child = right];
			}
			if (QT::compare(elt, c) <= 0) {
				break;
			}
			queue[k] = c;
			QT::setIndex(c, k);
			k = child;
		}
		queue[k] = elt;
		QT::setIndex(elt, k);
	}

	template <QueueTraits QT>
	bool IndexedPriorityQueue<QT>::empty() {
		return queue.empty();
	}

	template <QueueTraits QT>
	void IndexedPriorityQueue<QT>::push(EltHandle elt) {
		queue.push_back(elt);
		siftUp(queue.size() - 1, elt);
	}

	template <QueueTraits QT>
	QT::EltHandle IndexedPriorityQueue<QT>::pop() {
		EltHandle result = queue[0];
        for (auto eltH : queue) {
            if (QT::compare(eltH, result) < 0) {
                std::cout << "! Problem !" << std::endl;
                auto clps1 = result->data().collapse;
                auto clps2 = eltH->data().collapse;
                if (clps1.has_value() && clps2.has_value()) {
                    std::cout << clps1->cost << " is not smallest, but " << clps2->cost << " is." << std::endl;
                }
            }
        }
		QT::setIndex(result, -1);

		EltHandle last = queue[queue.size() - 1];
		queue.pop_back();
		if (!queue.empty()) {
			siftDown(0, last);
		}

		return result;
	}

	template <QueueTraits QT>
	bool IndexedPriorityQueue<QT>::remove(EltHandle elt) {
		int id = QT::getIndex(elt);
		if (id < 0 || id >= queue.size() || queue[id] != elt) {
			return false;
		}
		else {
			QT::setIndex(elt, -1);
			if (id == queue.size() - 1) {
				queue.pop_back();
			}
			else {
				EltHandle moved = queue[queue.size() - 1];
				queue.pop_back();
				siftDown(id, moved);
				if (queue[id] == moved) {
					siftDown(id, moved);
				}
			}

			return true;
		}
	}

	template <QueueTraits QT>
	bool IndexedPriorityQueue<QT>::contains(EltHandle elt) {
		int id = QT::getIndex(elt);
		if (id < 0 || id >= queue.size()) {
			return false;
		}
		else {
			return queue[id] == elt;
		}
	}

	template <QueueTraits QT>
	void IndexedPriorityQueue<QT>::update(EltHandle elt) {
		siftUp(QT::getIndex(elt), elt);
		siftDown(QT::getIndex(elt), elt);
	}

} // namespace cartocrow::simplification