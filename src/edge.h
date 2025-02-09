#ifndef CWDT_EDGE_H
#define CWDT_EDGE_H

#include <array>
#include <iostream>

namespace CWDT {
	class edge {
	public:
		edge(const int v1, const int v2) {
			if (v1 < v2) {
				ev[0] = v1; ev[1] = v2;
			}
			else {
				ev[0] = v2; ev[1] = v1;
			}
		}

		bool operator==(edge const& e) const {
			return (ev[0] == e.ev[0] && ev[1] == e.ev[1]);
		}

		int& operator[](int i) {
			return ev[i];
		}

		const int& operator[](int i) const {
			return ev[i];
		}

	private:
		std::array<int, 2> ev; /* sorted vertex indices, v.first < v.second */
	};
}

template<>
struct std::hash<CWDT::edge> {
	size_t operator()(CWDT::edge const& key) const {
		return (size_t)key[0] * (size_t)key[1];
	}
};

#endif
