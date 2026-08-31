#include "Dictionary.h"
#include "Parser.h"
#include "StorageBin.h"
#include "Surface.h"
#include "SurfaceComp.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace
{
void require(bool condition, const std::string &message)
{
	if (!condition)
	{
		throw std::runtime_error(message);
	}
}

void require_near
(
	const std::map<int, double> &expected,
	const std::map<int, double> &actual,
	const std::string &operation
)
{
	require(expected.size() == actual.size(), operation + " changed map size");
	for (std::map<int, double>::const_iterator it = expected.begin();
		 it != expected.end(); ++it)
	{
		std::map<int, double>::const_iterator found = actual.find(it->first);
		require(found != actual.end(), operation + " lost a species key");
		double scale = std::max(std::fabs(it->second), 1e-30);
		require
		(
			std::fabs(found->second - it->second) <= 1e-12 * scale,
			operation + " changed a species concentration"
		);
	}
}

cxxSurface make_surface()
{
	cxxSurface surface;
	surface.Set_n_user(17);
	surface.Get_species_map()[7] = 1.23456789012345e-8;
	surface.Get_species_map()[42] = 9.87654321098765e-6;

	std::vector<cxxSurfaceComp> components(1);
	components[0].Set_formula("SurfOH");
	components[0].Set_charge_name("Surf");
	components[0].Set_master_element("Surf");
	surface.Set_surface_comps(components);
	return surface;
}

void test_raw_and_modify()
{
	cxxSurface source = make_surface();
	std::map<int, double> expected = source.Get_species_map();
	std::ostringstream raw;
	source.dump_raw(raw, 0);
	raw << "END\n";

	std::istringstream raw_input(raw.str());
	CParser raw_parser(raw_input);
	cxxStorageBin bin;
	bin.read_raw(raw_parser);
	cxxSurface *loaded = bin.Get_Surface(17);
	require(loaded != NULL, "SURFACE_RAW did not restore the surface");
	require_near(expected, loaded->Get_species_map(), "SURFACE_RAW");
	std::map<int, double> before_modify = loaded->Get_species_map();

	std::istringstream modify_input
	(
		"SURFACE_MODIFY 17\n"
		"    -transport 1\n"
		"END\n"
	);
	CParser modify_parser(modify_input);
	bin.read_raw(modify_parser);
	loaded = bin.Get_Surface(17);
	require(loaded != NULL, "SURFACE_MODIFY removed the surface");
	require(loaded->Get_transport(), "SURFACE_MODIFY did not update transport");
	require
	(
		loaded->Get_species_map() == before_modify,
		"unrelated SURFACE_MODIFY erased surface species"
	);
}

void test_serialization_and_mixing()
{
	cxxSurface source = make_surface();
	std::map<int, double> expected = source.Get_species_map();
	Dictionary dictionary;
	std::vector<int> ints;
	std::vector<double> doubles;
	source.Serialize(dictionary, ints, doubles);

	cxxSurface restored;
	int int_index = 0;
	int double_index = 0;
	restored.Deserialize(dictionary, ints, doubles, int_index, double_index);
	require
	(
		restored.Get_species_map() == expected,
		"serialization did not preserve surface species"
	);

	cxxSurface copied;
	copied.add(source, 3.0);
	require
	(
		copied.Get_species_map() == expected,
		"one-source copy did not preserve surface species"
	);

	cxxSurface second = source;
	second.Get_species_map()[7] = 5.0;
	copied.add(second, 2.0);
	require
	(
		copied.Get_species_map().empty(),
		"genuine mixture retained derived surface species"
	);
}
}

int main()
{
	try
	{
		test_raw_and_modify();
		test_serialization_and_mixing();
		std::cout << "Surface species restart state preserved\n";
		return 0;
	}
	catch (const std::exception &error)
	{
		std::cerr << error.what() << '\n';
		return 2;
	}
}
