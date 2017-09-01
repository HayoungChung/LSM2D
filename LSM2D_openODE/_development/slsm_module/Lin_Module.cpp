#include "Lin_Module.h"

#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>

Lin_Module::Lin_Module(){
};

fea::Mesh* Lin_Module::get_CMesh(){
	return CMesh;
}

fea::ElasticitySensitivities<fea::BilinearQuad>* Lin_Module::get_CSensitivities(){
	return CSensitivities;
}

slsm::Mesh* Lin_Module::get_mesh(){
	return mesh;
}

slsm::LevelSet* Lin_Module::get_levelset(){
	return levelSet;
}

slsm::Boundary* Lin_Module::get_boundary(){
	return boundary;
}

std::vector<double>* Lin_Module::get_Field(){
	return &Field;
}

std::vector<double> Lin_Module::get_Field2(){
	return Field;
}

void Lin_Module::run (fea::Mesh& CMesh, slsm::Mesh& mesh, slsm::LevelSet& levelSet, 
								slsm::Boundary& boundary, slsm::InputOutput& io, std::vector<double>& Field)
{
    //     //Problem Dimension: 2D
    const unsigned int SpaceDimension = 2;
    double moveLimit = 0.9; // CFL limit
    double maxTime = 200; 
    double minArea = 0.5; // Vmin
    double sampleInterval = 1;
    double nextSample = 1;
        
    // //int NelX = 160, NelY = 80;        
    // double Xmin = 0.0, Xmax = 160.0, Ymin = 0.0, Ymax = 80.0;

    // //fea::Mesh CMesh; // Hayoung
    // CMesh.CreateMesh (Xmin, Xmax, Ymin, Ymax, NelX, NelY);
	
    //CMesh.PrintMesh(0);

    // Initialise a 160 x 80 non-periodic mesh.
    //slsm::Mesh mesh(NelX, NelY, false); // Hayoung
    double meshArea = mesh.width * mesh.height;

    // Create initial topology with holes.
    // std::vector <slsm::Hole> holes;
    // holes.push_back (slsm::Hole(16, 14, 5));
    // holes.push_back (slsm::Hole(32, 27, 5));
    // holes.push_back (slsm::Hole(48, 14, 5));
    // holes.push_back (slsm::Hole(64, 27, 5));
    // holes.push_back (slsm::Hole(80, 14, 5));
    // holes.push_back (slsm::Hole(96, 27, 5));
    // holes.push_back (slsm::Hole(112, 14, 5));
    // holes.push_back (slsm::Hole(128, 27, 5));
    // holes.push_back (slsm::Hole(144, 14, 5));
    // holes.push_back (slsm::Hole(16, 40, 5));
    // holes.push_back (slsm::Hole(32, 53, 5));
    // holes.push_back (slsm::Hole(48, 40, 5));
    // holes.push_back (slsm::Hole(64, 53, 5));
    // holes.push_back (slsm::Hole(80, 40, 5));
    // holes.push_back (slsm::Hole(96, 53, 5));
    // holes.push_back (slsm::Hole(112, 40, 5));
    // holes.push_back (slsm::Hole(128, 53, 5));
    // holes.push_back (slsm::Hole(144, 40, 5));
    // holes.push_back (slsm::Hole(16, 66, 5));
    // holes.push_back (slsm::Hole(48, 66, 5));
    // holes.push_back (slsm::Hole(80, 66, 5));
    // holes.push_back (slsm::Hole(112, 66, 5));
    // holes.push_back (slsm::Hole(144, 66, 5));
	
    // Initialise the level set object.
    //slsm::LevelSet levelSet(mesh, holes, moveLimit, 6, false); // Hayoung
	// levelSet(mesh,holes,moveLimit,6,false);

    // Initialise io object.
    // slsm::InputOutput io;

    // Reinitialise the level set to a signed distance function.
    levelSet.reinitialise();

    // Initialise the boundary object.
    // slsm::Boundary boundary(levelSet);// Hayoung

    // Perform initial boundary discretisation.
    boundary.discretise();	

    // Compute the initial element area fractions.
    boundary.computeAreaFractions();

    // Compute the initial boundary point normal vectors.
    boundary.computeNormalVectors();

	//Material - Elasticity	
    double E, v;
	E = 1.0; //N/mm^2
	v = 0.3;
	bool Isotropy = true;
	double Density = 7850; //kg/m^3
	fea::CHakElasticity CElasticityMaterial(Isotropy, Density, SpaceDimension);
	bool isPlaneStress = true;
	
	std::vector <double> BoundarySensitivities; // Hayoung
	
	//Computing the Elasticity tensor
	CElasticityMaterial.ComputeProperties (E, v, isPlaneStress); 
	std::vector <fea::CHakElasticity> Materials (1, CElasticityMaterial);

	//Integration Points
	fea::IntegrationPoint CIntegrationPoint (SpaceDimension);
	std::vector <fea::IntegrationPoint> IntegrationPoints (4, CIntegrationPoint);
	IntegrationPoints [0].Coordinate = { -0.577350269, -0.577350269 }; IntegrationPoints [0].weight = 1;
	IntegrationPoints [1].Coordinate = { 0.577350269, -0.577350269 }; IntegrationPoints [1].weight = 1; 
	IntegrationPoints [2].Coordinate = { 0.577350269, 0.577350269 }; IntegrationPoints [2].weight = 1;
	IntegrationPoints [3].Coordinate = { -0.577350269, 0.577350269 }; IntegrationPoints [3].weight = 1;
	std::string Interpolation = "BilinearQuad";
	unsigned int DoFperNode = 2;
	double Thickness = 1.0;
	
	//Modeling Void Material in FE
	double VoidMaterial = 1e-6;
	
	for (unsigned int i = 0; i < CMesh.nElements; i = i + 1)
	{
		CMesh.Elements [i].IntegrationPoints = &IntegrationPoints; //Assigning the address of one element's integration points to all the other elements to save memory
		CMesh.Elements [i].Interpolation = &Interpolation;
		CMesh.Elements [i].DoFperNode = &DoFperNode;
		CMesh.Elements [i].Thickness = &Thickness;
		if (mesh.elements [i].area < VoidMaterial)
			CMesh.Elements [i].AreaFraction = VoidMaterial; // From SLSM
		else
			CMesh.Elements [i].AreaFraction = mesh.elements [i].area; //From SLSM
	}
	CMesh.TotalDoF = 0;
	for (unsigned int i = 0; i < CMesh.nNodes; i = i + 1)
	{
		CMesh.Nodes [i].nDoF = &DoFperNode;
		CMesh.TotalDoF = CMesh.TotalDoF + *CMesh.Nodes [i].nDoF;
	}
	//Printing the base material properties
	Materials [0].PrintMaterial ();
	
	//Defining Dirichlet Boundary conditions
	double Tolerance = 0.001;

	CMesh.GetNodesByXCoordinate (0.0, Tolerance);

	fea::EssentialBoundaryCondition CDirichletBC;

	//In Cantilever, both x and y displacements are fixed
	std::vector <fea::EssentialBoundaryCondition> DirichletBC (CMesh.SelectedNodes.size () * 2); 
	unsigned int BCCounter = 0;
	for (unsigned int i = 0; i < CMesh.SelectedNodes.size (); i = i + 1)
	{
		CDirichletBC.NodeNumber = CMesh.SelectedNodes [i];
		for (unsigned int j = 0; j < DoFperNode; j = j + 1)
		{
			CDirichletBC.DoF = j;
			CDirichletBC.Value = 0.0;
			DirichletBC [BCCounter] = CDirichletBC;
			BCCounter = BCCounter + 1;
		}
	}
	CMesh.SelectedNodes.clear ();
	CMesh.CreateID (DirichletBC); //Create ID vector
	CMesh.CreateLM (DirichletBC); //Create LM arrays
	CMesh.nLoadCase = 1;
	// std::vector <double> Field (CMesh.TotalDoF*CMesh.nLoadCase);
	Field.resize(CMesh.TotalDoF*CMesh.nLoadCase);
	for (unsigned int i = 0; i < DirichletBC.size (); i = i + 1)
	{
		for (unsigned int j = 0; j < CMesh.nLoadCase; j = j + 1)
		{
			Field [CMesh.TotalDoF*j + DoFperNode*(DirichletBC [i].NodeNumber - 1) + DirichletBC [i].DoF] = DirichletBC [i].Value;
		}
	}

	//Defining Loads
	CMesh.DoFLoads.resize (CMesh.nLoadCase);
	CMesh.DoFLoads [0].resize (3); //Point Loads
	CMesh.GetNodesByXYCoordinates (160.0, 40.0, Tolerance);
	CMesh.DoFLoads [0] [0].resize (2); //Equation number and load value
	CMesh.DoFLoads [0] [0] [0] = CMesh.ID [CMesh.SelectedNodes [0] - 1] [1]; //Y direction load
	CMesh.DoFLoads [0] [0] [1] = -5; //N
	CMesh.SelectedNodes.clear ();
	
	CMesh.GetNodesByXYCoordinates (160.0, 39.0, Tolerance);
	CMesh.DoFLoads [0] [1].resize (2); //Equation number and load value
	CMesh.DoFLoads [0] [1] [0] = CMesh.ID [CMesh.SelectedNodes [0] - 1] [1]; //Y direction load
	CMesh.DoFLoads [0] [1] [1] = -2.5; //N
	CMesh.SelectedNodes.clear ();
	
	CMesh.GetNodesByXYCoordinates (160.0, 41.0, Tolerance);
	CMesh.DoFLoads [0] [2].resize (2); //Equation number and load value
	CMesh.DoFLoads [0] [2] [0] = CMesh.ID [CMesh.SelectedNodes [0] - 1] [1]; //Y direction load
	CMesh.DoFLoads [0] [2] [1] = -2.5; //N
	CMesh.SelectedNodes.clear ();

	//Material and Element Mappings
	CMesh.isStructured = true;
	for (unsigned int i = 0; i < CMesh.nElements; i = i + 1)
	{
		CMesh.Elements [i].MaterialNumber = 0;
	}
	Materials [0].ElementNumber = 0;
	//CMesh.PrintMesh (3);

	//Linear Elasticity with Bilinear Quadrilateral shape functions
	fea::CHakLinearElasticity<fea::BilinearQuad> CLinearElasticity (CMesh, Materials, SpaceDimension);

    CLinearElasticity.Assembly ();
    fea::CHakLinearSolver CSolve (CMesh, CLinearElasticity.irn, CLinearElasticity.jcn, CLinearElasticity.sparseK, CLinearElasticity.F, Field, 0);
    CSolve.MA57 ();

	io.saveLevelSetVTK(0, levelSet);	

	//fea::ElasticitySensitivities <fea::BilinearQuad> CSensitivities (CMesh, Field, SpaceDimension, Materials, 0.01);
};

using namespace boost::python;

class dummyA{};
class dummyB{};

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(cV1,slsm::LevelSet::computeVelocities, 1, 1);
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(cV2,slsm::LevelSet::computeVelocities, 4, 4);

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(cM1,fea::Mesh::CreateMesh,6,6);
// BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(cM2,fea::Mesh::CreateMesh,6,6);

BOOST_PYTHON_MODULE(Lin_Module){
	class_<std::vector<double>>("vector__Field__")
		.def(vector_indexing_suite<std::vector<double>>())
		;

	// class_<std::vector<fea::CHakElasticity>>("vector_Material_")
	// 	.def(vector_indexing_suite<std::vector<fea::CHakElasticity>>())
	// 	;
	
	
	class_<Lin_Module>("Lin_Module")//,init<unsigned int, unsigned int>())
		// .add_property("NelX",&Lin_Module::NelX)
		.def("run",&Lin_Module::run)
		.def("get_CMesh",&Lin_Module::get_CMesh, return_value_policy<reference_existing_object>())
		.def("get_mesh",&Lin_Module::get_mesh, return_value_policy<reference_existing_object>())
		.def("get_levelset",&Lin_Module::get_levelset, return_value_policy<reference_existing_object>())
		.def("get_boundary",&Lin_Module::get_boundary, return_value_policy<reference_existing_object>())
		.def("get_CSensitivities",&Lin_Module::get_CSensitivities, return_value_policy<reference_existing_object>())
		.def("get_Field",&Lin_Module::get_Field, return_value_policy<reference_existing_object>())
		.def("get_Field2",&Lin_Module::get_Field2)//, return_value_policy<copy_const_reference>())

		.add_property("Field",&Lin_Module::Field)
		// .add_property("Materials",&Lin_Module::Materials)
		;		
		
	{scope fea = class_<dummyA>("fea"); 

	class_<fea::Node>("Node",init<unsigned int>())
		.add_property("Coordinate",&fea::Node::Coordinate)
		;

	class_<std::vector<fea::Node>>("vector__Nodes__")
		.def(vector_indexing_suite<std::vector<fea::Node>>())
		;
	
	class_<fea::Mesh>("Mesh")
		.def("CreateMesh",static_cast<void(fea::Mesh::*)(double,double,double,double,unsigned int,unsigned int)>(&fea::Mesh::CreateMesh),cM1())
		.add_property("Nodes",&fea::Mesh::Nodes)
		;

	}
	// class_<fea::CHakElasticity>("CHakElasticity",init<bool,double,const unsigned int>())
	// 	;
	
	{scope slsm = class_<dummyB>("slsm");

	class_<slsm::Hole>("Hole",init<double,double,double>())
		;
	class_<std::vector<slsm::Hole>>("vector__Holes__")
		.def(vector_indexing_suite<std::vector<slsm::Hole>>())
		;
		
	class_<slsm::Mesh>("Mesh",init<unsigned int, unsigned int, bool>())
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
		;
	class_<fea::ElasticitySensitivities<fea::BilinearQuad>>("ElasticitySensitivities",
		init<fea::Mesh&, std::vector <double>&, unsigned int, std::vector <fea::CHakElasticity>&, double>())
		;
	
	class_<slsm::BoundaryPoint>("BoundaryPoint")
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
		;

	class_<slsm::InputOutput>("InputOutput")
		;
	}

}	
