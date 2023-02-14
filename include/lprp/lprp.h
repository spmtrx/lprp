#ifndef LPRP_H
#define LPRP_H

#include"spmtrx/spmcrs.h"
namespace spm {
	class laplacian {
		public:
			enum class option {
				freq = 0,
				npmi
			};
			using seed = std::vector<unsigned int>;
			using instance = std::pair<unsigned int, double>;
			using pattern = std::pair<unsigned int, double>;
			using rst = std::pair<std::vector<instance>, std::vector<pattern> >;
			laplacian();
			laplacian(const char *tsv, option = option::freq, double threshold = 0.);
			~laplacian();
			virtual rst labelprop(seed& s, int hop, double a = 0.01);
		protected:
		private:
			spmcrs _w;
			spmcrs _w_t;
			spmcrs _diag_a;
			spmcrs _make_seed(seed& s);
			std::vector<instance> _make_instance(spmcrs& f);
			friend class boost::serialization::access;
			template<class Archive>
				void serialize(Archive& ar, unsigned int version) {
					ar & _w;
					ar & _w_t;
					ar & _diag_a;
				}
	};
}
#endif /** LPRP_H **/
