/***************************************************************************
 *  ThunderEgg, a library for solving Poisson's equation on adaptively
 *  refined block-structured Cartesian grids
 *
 *  Copyright (C) 2019-2020 ThunderEgg Developers. See AUTHORS.md file at the
 *  top-level directory.
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <https://www.gnu.org/licenses/>.
 ***************************************************************************/
#include <ThunderEgg/Domain.h>
#include <ThunderEgg/tpl/json.hpp>
#include <fstream>
#include <string>
template <int D> class DomainReader
{
	private:
	bool                                      neumann;
	std::array<int, D>                        ns;
	int                                       num_ghost;
	std::shared_ptr<ThunderEgg::Domain<D>>    coarser_domain;
	std::shared_ptr<ThunderEgg::Domain<D>>    finer_domain;
	std::shared_ptr<ThunderEgg::PatchInfo<D>> parsePatch(nlohmann::json &patch_j)
	{
		auto pinfo             = std::make_shared<ThunderEgg::PatchInfo<D>>();
		*pinfo                 = patch_j.get<ThunderEgg::PatchInfo<D>>();
		pinfo->num_ghost_cells = num_ghost;
		pinfo->ns              = ns;
		for (size_t d = 0; d < D; d++) {
			pinfo->spacings[d] /= ns[d];
		}

		for (ThunderEgg::Side<D> s : ThunderEgg::Side<D>::getValues()) {
			if (!pinfo->hasNbr(s)) {
				pinfo->neumann[s.getIndex()] = neumann;
			}
		}
		return pinfo;
	}

	public:
	DomainReader(std::string file_name, std::array<int, D> ns_in, int num_ghost_in,
	             bool neumann_in = false)
	: neumann(neumann_in), ns(ns_in), num_ghost(num_ghost_in)
	{
		int rank;
		MPI_Comm_rank(MPI_COMM_WORLD, &rank);

		nlohmann::json j;
		std::ifstream  input_stream(file_name);
		if (!input_stream.good()) {
			throw "could not open file";
		}
		input_stream >> j;
		input_stream.close();
		std::map<int, std::shared_ptr<ThunderEgg::PatchInfo<D>>> finer_map;
		for (nlohmann::json &patch_j : j.at("levels")[0]) {
			auto patch = parsePatch(patch_j);
			if (patch->rank == rank)
				finer_map[patch->id] = patch;
		}
		finer_domain = std::make_shared<ThunderEgg::Domain<D>>(finer_map, ns, num_ghost);
		std::map<int, std::shared_ptr<ThunderEgg::PatchInfo<D>>> coarser_map;
		for (nlohmann::json &patch_j : j.at("levels")[1]) {
			auto patch = parsePatch(patch_j);
			if (patch->rank == rank)
				coarser_map[patch->id] = patch;
		}
		coarser_domain = std::make_shared<ThunderEgg::Domain<D>>(coarser_map, ns, num_ghost);
	}
	std::shared_ptr<ThunderEgg::Domain<D>> getCoarserDomain()
	{
		return coarser_domain;
	}
	std::shared_ptr<ThunderEgg::Domain<D>> getFinerDomain()
	{
		return finer_domain;
	}
};
extern template class DomainReader<2>;