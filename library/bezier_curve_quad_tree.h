#pragma once

#include <cartocrow/core/cubic_bezier.h>
#include <cartocrow/core/core.h>
#include <cartocrow/core/polyline.h>

namespace cartocrow::curved_simplification {
	namespace detail {
		template <typename PHandle, typename K> class PQTNode;
	}

	template <typename PHandle, typename K> class BezierCurveQuadTree {
	private:
		using Node = detail::PQTNode<PHandle, K>;

    public:
        Node* root;

    private:
		int maxdepth;
		Number<K> fuzziness;
        std::function<CubicBezierCurve(PHandle)> m_toBezierCurve;

		// Does the rectangle enclose the (possibly infinite) node?
		bool encloses(Rectangle<K>& rect, Node* node);

		// Is the (possibly infinite) node disjoint from the rectangle
		bool disjoint(Node* node, Rectangle<K>& rect);

		void findOverlappedRecursive(Node* n, Rectangle<K>& query, std::function<void(PHandle)> act);

		template <bool extend> Node* find(PHandle elt);

	public:
		BezierCurveQuadTree(Rectangle<K>& box, int depth, Number<K> fuzz, std::function<CubicBezierCurve(PHandle)> toBezierCurve);
		~BezierCurveQuadTree();

		void clear();
		void insert(PHandle elt);
		bool remove(PHandle elt);

		void findOverlapped(Rectangle<K>& query, std::function<void(PHandle)> act);
	};

} // namespace cartocrow::simplification

//#include "bezier_curve_quad_tree.hpp"

#include "utils.h"

namespace cartocrow::curved_simplification {

    namespace detail {
        template <typename PHandle, typename K> struct PQTNode {

            using Node = PQTNode<PHandle, K>;

            Node* parent;
            Node* lt = nullptr;
            Node* lb = nullptr;
            Node* rt = nullptr;
            Node* rb = nullptr;
            Rectangle<K>* rect;
            bool inf_left, inf_right, inf_bottom, inf_top;
            std::vector<PHandle>* elts;

            PQTNode(PQTNode<PHandle, K>* parent, Rectangle<K>* rect, bool inf_left, bool inf_right,
                    bool inf_bottom, bool inf_top)
                    : parent(parent), rect(rect), inf_left(inf_left), inf_right(inf_right),
                      inf_bottom(inf_bottom), inf_top(inf_top) {
                elts = new std::vector<PHandle>();
            }

            ~PQTNode() {
                delete rect;
                delete lt;
                delete lb;
                delete rt;
                delete rb;
                delete elts;
            }
        };
    }

    template <typename PHandle, typename K>
    bool BezierCurveQuadTree<PHandle, K>::encloses(Rectangle<K>& rect, Node* node) {
        // node extends further to the left
        if (node->inf_left || rect.xmin() > node->rect->xmin()) {
            return false;
        }

        // node extends further to the right
        if (node->inf_right || rect.xmax() < node->rect->xmax()) {
            return false;
        }

        // node extends further to the bottom
        if (node->inf_bottom || rect.ymin() > node->rect->ymin()) {
            return false;
        }

        // node extends further to the top
        if (node->inf_top || rect.ymax() < node->rect->ymax()) {
            return false;
        }

        return true;
    }

    template <typename PHandle, typename K>
    bool BezierCurveQuadTree<PHandle, K>::disjoint(Node* node, Rectangle<K>& rect) {
        // rect left of node
        if (!node->inf_left && rect.xmax() < node->rect->xmin()) {
            return true;
        }

        // rect right of node
        if (!node->inf_right && rect.xmin() > node->rect->xmax()) {
            return true;
        }

        // rect below node
        if (!node->inf_bottom && rect.ymax() < node->rect->ymin()) {
            return true;
        }

        // rect above node
        if (!node->inf_top && rect.ymin() > node->rect->ymax()) {
            return true;
        }

        return false;
    }

    template <typename PHandle, typename K>
    void BezierCurveQuadTree<PHandle, K>::findOverlappedRecursive(Node* n, Rectangle<K>& query, std::function<void(PHandle)> act) {

        if (n == nullptr) {
            return;
        }
        else if (disjoint(n, query)) {
            return;
        }

        if (encloses(query, n)) {
            for (PHandle elt : *(n->elts)) {
                act(elt);
            }
        }
        else {
            for (PHandle elt : *(n->elts)) {
                CubicBezierCurve curve = m_toBezierCurve(elt);
                if (utils::overlaps(query, curve)) {
                    act(elt);
                }
            }
        }

        findOverlappedRecursive(n->lb, query, act);
        findOverlappedRecursive(n->lt, query, act);
        findOverlappedRecursive(n->rb, query, act);
        findOverlappedRecursive(n->rt, query, act);
    }

    template <typename PHandle, typename K>
    template <bool extend>
    detail::PQTNode<PHandle, K>* BezierCurveQuadTree<PHandle, K>::find(PHandle elt) {
        int d = 0;
        Node* n = root;
        CubicBezierCurve curve = m_toBezierCurve(elt);

        auto bbox = curve.bbox();

        Number<K> seg_left = bbox.xmin();
        Number<K> seg_right = bbox.xmax();
        Number<K> seg_bottom = bbox.ymin();
        Number<K> seg_top = bbox.ymax();

        // NB: extending, n will never become null
        while (d < maxdepth) {
            d++;

            Number<K> xmid = (n->rect->xmin() + n->rect->xmax()) / 2;
            Number<K> ymid = (n->rect->ymin() + n->rect->ymax()) / 2;
            Number<K> fx = (n->rect->xmax() - n->rect->xmin()) / 2.0 * fuzziness;
            Number<K> fy = (n->rect->ymax() - n->rect->ymin()) / 2.0 * fuzziness;

            if (seg_right <= xmid + fx) {
                // left
                if (seg_top <= ymid + fy) {
                    // bottom
                    if (n->lb == nullptr) {
                        if constexpr (extend) {
                            n->lb =
                                    new Node(n,
                                             new Rectangle<K>(n->rect->xmin(), n->rect->ymin(),
                                                              xmid + fx, ymid + fx),
                                             n->inf_left, false, n->inf_bottom, false);
                        }
                        else {
                            return nullptr;
                        }
                    }
                    n = n->lb;
                }
                else if (seg_bottom >= ymid - fy) {
                    // top
                    if (n->lt == nullptr) {
                        if constexpr (extend) {
                            n->lt = new Node(n,
                                             new Rectangle<K>(n->rect->xmin(), ymid - fy,
                                                              xmid + fx, n->rect->ymax()),
                                             n->inf_left, false, false, n->inf_top);
                        }
                        else {
                            return nullptr;
                        }
                    }
                    n = n->lt;
                }
                else {
                    // here
                    return n;
                }
            }
            else if (seg_left >= xmid - fx) {
                // right
                if (seg_top <= ymid + fy) {
                    // bottom
                    if (n->rb == nullptr) {
                        if constexpr (extend) {
                            n->rb = new Node(n,
                                             new Rectangle<K>(xmid - fx, n->rect->ymin(),
                                                              n->rect->xmax(), ymid + fy),
                                             false, n->inf_right, n->inf_bottom, false);
                        }
                        else {
                            return nullptr;
                        }
                    }
                    n = n->rb;
                }
                else if (seg_bottom >= ymid - fy) {
                    // top
                    if (n->rt == nullptr) {
                        if constexpr (extend) {
                            n->rt =
                                    new Node(n,
                                             new Rectangle<K>(xmid - fx, ymid - fy,
                                                              n->rect->xmax(), n->rect->ymax()),
                                             false, n->inf_right, false, n->inf_top);
                        }
                        else {
                            return nullptr;
                        }
                    }
                    n = n->rt;
                }
                else {
                    // here
                    return n;
                }
            }
            else {
                // here
                return n;
            }
        }

        return n;
    }

    template <typename PHandle, typename K>
    BezierCurveQuadTree<PHandle, K>::BezierCurveQuadTree(Rectangle<K>& box, int depth, Number<K> fuzz, std::function<CubicBezierCurve(PHandle)> toBezier) {
        root = new Node(nullptr, new Rectangle<K>(box[0], box[2]), true, true, true, true);
        maxdepth = depth;
        fuzziness = fuzz;
        m_toBezierCurve = toBezier;
    }

    template <typename PHandle, typename K>
    BezierCurveQuadTree<PHandle, K>::~BezierCurveQuadTree() {
        delete root;
    }

    template <typename PHandle, typename K>
    void BezierCurveQuadTree<PHandle, K>::clear() {
        // just delete all nodes except the root
        delete root->lb;
        delete root->rb;
        delete root->lt;
        delete root->rt;
        root->lb = nullptr;
        root->rb = nullptr;
        root->lt = nullptr;
        root->rt = nullptr;

        // and empty the elements at the root
        root->elts->clear();
    }

    template <typename PHandle, typename K>
    void BezierCurveQuadTree<PHandle, K>::insert(PHandle elt) {
        Node* n = find<true>(elt);
        n->elts->push_back(elt);
    }

    template <typename PHandle, typename K>
    bool BezierCurveQuadTree<PHandle, K>::remove(PHandle elt) {
        Node* n = find<false>(elt);

        if (n == nullptr) {
            // it's supposed leaf doesn't exist
            return false;
        }

        auto position = std::find(n->elts->begin(), n->elts->end(), elt);
        if (position == n->elts->end()) {
            // the leaf doesn't actually contain the element
            return false;
        }

        // it's contained: erase and clean up the tree if possible
        n->elts->erase(position);

        while (n != root && n->elts->empty() && n->lb == nullptr && n->lt == nullptr &&
               n->rb == nullptr && n->rt == nullptr) {
            Node* p = n->parent;
            if (p->lb == n) {
                p->lb = nullptr;
            }
            else if (p->lt == n) {
                p->lt = nullptr;
            }
            else if (p->rb == n) {
                p->rb = nullptr;
            }
            else if (p->rt == n) {
                p->rt = nullptr;
            }
            delete n;
            n = p;
        }

        return true;
    }

    template <typename PHandle, typename K>
    void BezierCurveQuadTree<PHandle, K>::findOverlapped(Rectangle<K>& query, std::function<void(PHandle)> act) {
        findOverlappedRecursive(root, query, act);
    }

} // namespace cartocrow::simplification

