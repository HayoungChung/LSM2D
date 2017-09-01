#include "slsm_Module.h"

#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>

slsm_Module::slsm_Module(){};

double slsm_Module::wrapper_optimise(std::vector<slsm::BoundaryPoint>& boundaryPoints_, const std::vector<double>& constraintDistances_,
        				std::vector<double>& lambdas_, std::vector <double>& Multipliers_, double moveLimit)
{
	double timeStep;
    slsm::Optimise optimise(boundaryPoints_, constraintDistances_, 
					lambdas_, Multipliers_, timeStep,moveLimit,false);
	optimise.solve();
	// std::cout << lambdas_[0] << "\t" << lambdas_[1] << "\n";
	// std::cout << Multipliers_[0] << "\t" << Multipliers_[1] << "\n";
	// std::cout << timeStep << "\n";

	return timeStep;
}

using namespace boost::python; 

class dummyA{};

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(cV1,slsm::LevelSet::computeVelocities, 1, 1);
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(cV2,slsm::LevelSet::computeVelocities, 4, 4);

BOOST_PYTHON_MODULE(slsm_Module){
    class_<std::vector<double>>("vector__double__")
		.def(vector_indexing_suite<std::vector<double>>())
        ;
		
	class_<slsm_Module>("slsm_Module") 
		.def("wrapper_optimise",&slsm_Module::wrapper_optimise)
		;


    {scope slsm = class_<dummyA>("slsm");
    class_<slsm::Hole>("Hole",init<double,double,double>())
        ;
	class_<std::vector<slsm::Hole>>("vector__Holes__")
		.def(vector_indexing_suite<std::vector<slsm::Hole>>())
		;

	class_<std::vector<slsm::Element>>("vector__Elements__")
		.def(vector_indexing_suite<std::vector<slsm::Element>>())
		;

	class_<std::vector<slsm::Node>>("vector__Nodes__")
		.def(vector_indexing_suite<std::vector<slsm::Node>>())
		;

    class_<slsm::Element>("Element")
        .add_property("area",&slsm::Element::area)
        ;

    class_<slsm::Mesh>("Mesh",init<unsigned int, unsigned int, bool>())
        .add_property("elements",&slsm::Mesh::elements)
        ;
    
    class_<slsm::LevelSet>("LevelSet",init<slsm::Mesh&, std::vector<slsm::Hole>&, double, unsigned int, bool>())
		.def("ComputeGradients",&slsm::LevelSet::computeGradients)
		.def("update",&slsm::LevelSet::update)
		.def("computeVelocities",static_cast< void(slsm::LevelSet::*)(const std::vector<slsm::BoundaryPoint>&)>
			(&slsm::LevelSet::computeVelocities),cV1())
        .def("computeVelocities",
			static_cast< double(slsm::LevelSet::*)(std::vector<slsm::BoundaryPoint>&, double&, const double, slsm::MersenneTwister&)>
			(&slsm::LevelSet::computeVelocities),cV2())  
		.def("reinitialise",&slsm::LevelSet::reinitialise)		
		.add_property("moveLimit",&slsm::LevelSet::moveLimit)
		;
		
    class_<slsm::Coord>("Coord")
        .add_property("x",&slsm::Coord::x)
        .add_property("y",&slsm::Coord::y)
        ;

	class_<slsm::BoundaryPoint>("BoundaryPoint") // STRUCT
		.add_property("sensitivities",&slsm::BoundaryPoint::sensitivities)
		.add_property("coord",&slsm::BoundaryPoint::coord)
		;
	class_<std::vector<slsm::BoundaryPoint>>("vector_BoundaryPoint_")
		.def(vector_indexing_suite<std::vector<slsm::BoundaryPoint>>())
		;

	class_<slsm::Boundary>("Boundary",init<slsm::LevelSet&>()) 
		.def("discretise",&slsm::Boundary::discretise)
		.add_property("points",&slsm::Boundary::points)
		.def("computeAreaFractions",&slsm::Boundary::computeAreaFractions)
		.def("ComputeNormalVectors",&slsm::Boundary::computeNormalVectors)
        .add_property("area",&slsm::Boundary::area)
		;

	class_<slsm::InputOutput>("InputOutput")
	// 	.def("saveLevelSetVTK",&slsm::InputOutput::saveLevelSetVTK)
	// 	.def("saveBoundarySegmentsTXT",&slsm::InputOutput::saveBoundarySegmentsTXT)
		;

    class_<slsm::Optimise>("Optimise",init<std::vector<slsm::BoundaryPoint>&, std::vector<double>&, std::vector<double>&, std::vector<double>&, double&>() )
        .def("solve",&slsm::Optimise::solve)
        ;
    

    }
}
