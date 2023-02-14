#include"spmtrx/spmutil.h"
#include"lprp/lprp.h"

using namespace spm;
using namespace std;

laplacian::laplacian() {
}

laplacian::laplacian(const char *tsv, option o, double threshold) {
	_w = spmutil::make_crs(tsv);
	if (o == option::npmi) {
		_w = spmutil::npmi(_w, threshold);
	}
	_w_t = _w.transpose();
	_diag_a = _w.sqrt_inv_diag_adjacency();
}

laplacian::~laplacian() {
}

laplacian::rst laplacian::labelprop(seed& s, int hop, double a) {
	if (s.empty()) {
		return make_pair(vector<instance>(), vector<pattern>());
	}
	spmcrs f_0 = _make_seed(s);
	spmcrs f_t = f_0;
	spmcrs f_p; // pattern
	f_0 *= (1.-a);
	for (int i = 0; i < hop; ++i) {
		f_t *= a;
		f_p = _w_t*(_diag_a*f_t);
		f_t = (_diag_a*(_w*f_p))-f_t+f_0;
	}
	return make_pair(_make_instance(f_t),_make_instance(f_p));
}

spmcrs laplacian::_make_seed(seed& s) {
	spmcrs v_t;
	sort(s.begin(), s.end());
	for (auto it = s.begin(); it != s.end(); ++it) {
		v_t.set(*it, 1);
	}
	v_t.cr();
	return v_t.transpose();
}

vector<laplacian::instance> laplacian::_make_instance(spmcrs& f) {
	vector<instance> i;
	int j = 0;
	for (auto it = f.begin(); it != f.end(); ++it, ++j) {
		for (auto k = *it; k < *(it+1); ++k) {
			i.push_back(make_pair(j, f.get_val(k)));
		}
	}
	return i;
}
