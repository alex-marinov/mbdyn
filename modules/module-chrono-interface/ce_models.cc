// =============================================================================
// PROJECT CHRONO - http://projectchrono.org
//
// Copyright (c) 2014 projectchrono.org
// All rights reserved.
//
// Use of this source code is governed by a BSD-style license that can be found
// in the LICENSE file at the top level of the distribution and at
// http://projectchrono.org/license-chrono.txt.
//
//
//
// A C::E Models: sliding point and a box linked by a spring-dash device;
// The global reference frame has Y up.
//
// =============================================================================

// Modified Date: 2020/06/09

#include "chrono/ChConfig.h"
#include "chrono/utils/ChUtilsCreators.h"
#include "chrono/utils/ChUtilsGenerators.h"
#include "chrono/physics/ChBodyEasy.h"
#include "chrono_parallel/physics/ChSystemParallel.h"
#include "chrono_parallel/solver/ChIterativeSolverParallel.h"
#include "chrono/assets/ChPointPointDrawing.h"
#ifdef CHRONO_OPENGL
#include "chrono_opengl/ChOpenGLWindow.h"
#endif

#include "mbdyn_ce.h"

using namespace chrono;
using namespace chrono::collision;

// --------------------------------------------------------------------------


extern "C" void
MBDyn_CE_CEModel_Create(ChSystemParallelNSC * pMBDyn_CE_CEModel)
{
	GetLog() << "\n\nCopyright (c) 2017 projectchrono.org\nChrono version: " << CHRONO_VERSION << "\n\n";

	// ----------------
	// 1. Parameters
	// ----------------
	const double PI = 3.141592653589793238462643;
	double radius = 0.0040;	
	auto mat = std::make_shared<ChMaterialSurfaceNSC>();
	mat->SetFriction(0.52f);


	double spring2_coef = 1602.7;
	double damping2_coef = 0.116;
	double fre = 12.0;//the A of excitation

	

	int num_threads = 4;
	uint max_iteration_normal = 100;
	uint max_iteration_sliding =100;
	uint max_iteration_spinning = 0;
	uint max_iteration_bilateral = 100;

	// ------------------------
	// 2. Create the parallel system
	// --------------------------
	//ChSystemParallelNSC* pCEmodel_system;

	//ChSystemParallelNSC system;
	pMBDyn_CE_CEModel->Set_G_acc(ChVector<>(0, -9.81, 0));


	// Set number of threads

	 int max_threads = CHOMPfunctions::GetNumProcs();
	 if (num_threads < max_threads)
		 num_threads = max_threads-1;
	CHOMPfunctions::SetNumThreads(num_threads);

	// Set solver settings
	//system.SetTimestepperType(ChTimestepper::Type::NEWMARK);
	pMBDyn_CE_CEModel->ChangeSolverType(SolverType::APGDREF);

	pMBDyn_CE_CEModel->GetSettings()->perform_thread_tuning = true;

	pMBDyn_CE_CEModel->GetSettings()->solver.solver_mode = SolverMode::SLIDING;
	pMBDyn_CE_CEModel->GetSettings()->solver.max_iteration_normal = max_iteration_normal;
	pMBDyn_CE_CEModel->GetSettings()->solver.max_iteration_sliding = max_iteration_sliding;
	pMBDyn_CE_CEModel->GetSettings()->solver.max_iteration_spinning = max_iteration_spinning;
	pMBDyn_CE_CEModel->GetSettings()->solver.max_iteration_bilateral = max_iteration_bilateral;
	pMBDyn_CE_CEModel->GetSettings()->solver.alpha = 0;
	pMBDyn_CE_CEModel->GetSettings()->solver.contact_recovery_speed = 1e5;
	pMBDyn_CE_CEModel->GetSettings()->solver.use_full_inertia_tensor = false;
	pMBDyn_CE_CEModel->GetSettings()->solver.tolerance = 0.0;

	pMBDyn_CE_CEModel->GetSettings()->collision.collision_envelope = 0.01 * radius;
	pMBDyn_CE_CEModel->GetSettings()->collision.narrowphase_algorithm = NarrowPhaseType::NARROWPHASE_HYBRID_MPR;
	pMBDyn_CE_CEModel->GetSettings()->collision.bins_per_axis = vec3(10, 10, 10);

	// ------------------------
	// 3. Create the ground
	// --------------------------
	auto ground = std::shared_ptr<ChBody>(pMBDyn_CE_CEModel->NewBody());
	pMBDyn_CE_CEModel->Add(ground);
	ground->SetPos(ChVector<>(0, 0, 0));
	ground->SetBodyFixed(true);
	ground->SetIdentifier(0);
	ground->SetMass(1.0);
	ground->SetCollide(false);
	//ground->GetCollisionModel()->ClearModel();
	utils::AddBoxGeometry(ground.get(), ChVector<>(0.2,0.005,0.2), ChVector<>(0, -0.005, 0));
	
	// ------------------------
	// 4. Create the sliding point
	auto body_1 = std::shared_ptr<ChBody>(pMBDyn_CE_CEModel->NewBody());
	pMBDyn_CE_CEModel->AddBody(body_1);
	body_1->SetIdentifier(1);
	body_1->SetPos(ChVector<>(-0.1, 0.0, 0));
	//body_1->SetBodyFixed(true);
	//body_1->SetCollide(true);
	//body_1->GetCollisionModel()->ClearModel();
	utils::AddBoxGeometry(body_1.get(), ChVector<>(0.01, 0.05, 0.01), ChVector<>(0, 0.05, 0));
	//body_1->GetCollisionModel()->BuildModel();
	//body_1->GetMaterialSurfaceNSC()->SetFriction(0.0f);
	//body_1->GetMaterialSurfaceNSC()->SetRestitution(0.0f);

	// 4-1. Create the constraint between ground and sliding point. 

	auto linku = std::make_shared<ChLinkLockLock>();
	linku->Initialize(body_1, ground, ChCoordsys<>(ChVector<>(0, 0, 0)));
	auto mmotion_x = std::make_shared<ChFunction_Sequence>();
	auto mmotion_x1 = std::make_shared<ChFunction_Const>();
	mmotion_x1->Set_yconst(0.000);
	auto mmotion_x2=std::make_shared<ChFunction_Sine>(0, fre, 0.001);  // phase freq ampl
	mmotion_x->InsertFunct(mmotion_x1, 0.1, 1.0, true);
	mmotion_x->InsertFunct(mmotion_x2, 20, 1.0, true);
	linku->SetMotion_X(mmotion_x);
	pMBDyn_CE_CEModel->Add(linku);
	// --------------------------

	// ------------------------
	// 5. Create the container
	// ------------------------

	// Parameters for the containing bin
	int binId = 3;
	double thickness = 0.005;
	double width = 0.058/2.0+2*thickness;//L
	double length = 0.038/2.0 +2* thickness;//W
	double height = 0.038/2.0 + thickness ;//H	
	auto container = std::shared_ptr<ChBody>(pMBDyn_CE_CEModel->NewBody());
	pMBDyn_CE_CEModel->Add(container);
	container->SetPos(ChVector<>(0, 0, 0));
	container->SetBodyFixed(false);
	container->SetIdentifier(2);
	container->SetMass(0.293);
	container->SetCollide(true);
	container->GetCollisionModel()->ClearModel();
	utils::AddBoxGeometry(container.get(), ChVector<>(width, thickness, length), ChVector<>(0, 0+thickness, 0));//floor
	utils::AddBoxGeometry(container.get(), ChVector<>(thickness, height, length), ChVector<>(-width + thickness, height + thickness, 0));//left
	utils::AddBoxGeometry(container.get(), ChVector<>(thickness, height, length), ChVector<>(width - thickness, height + thickness, 0));//right
	utils::AddBoxGeometry(container.get(), ChVector<>(width, height, thickness), ChVector<>(0, height + thickness, -length + thickness));//back
	utils::AddBoxGeometry(container.get(), ChVector<>(width, height, thickness), ChVector<>(0, height + thickness, length - thickness), ChQuaternion<>(1, 0, 0, 0), false);//front
	utils::AddBoxGeometry(container.get(), ChVector<>(width, thickness, length), ChVector<>(0, 2 * height + thickness, 0));//ceil

	container->GetCollisionModel()->BuildModel();

	container->SetMaterialSurface(mat);
	//container->GetMaterialSurfaceNSC()->SetFriction(0.52f);
	//container->GetMaterialSurfaceNSC()->SetRestitution(0.0f);
	
	// ------------------------
	//  Create the container using translational joint
	// ------------------------
	
	auto link_CG = std::make_shared<ChLinkLockPointLine>();
	link_CG->Initialize(container,ground, ChCoordsys<>(ChVector<>(0, 0, 0)));
	pMBDyn_CE_CEModel->AddLink(link_CG);
	
	// Create the spring between body_1 and container. The spring end points are
	// specified in the body relative frames.
	auto spring_2 = std::make_shared<ChLinkTSDA>();
	spring_2->Initialize(body_1, container, true, ChVector<>(0, 0, 0), ChVector<>(0, 0, 0), true);
	spring_2->SetSpringCoefficient(spring2_coef);
	spring_2->SetDampingCoefficient(damping2_coef);
	pMBDyn_CE_CEModel->AddLink(spring_2);

	std::cout << "CE_model is created:\n";

	/*#ifdef CHRONO_OPENGL
	// -------------------------------
	// Create the visualization window
	// -------------------------------

	opengl::ChOpenGLWindow& gl_window = opengl::ChOpenGLWindow::getInstance();
	gl_window.Initialize(1280, 720, "Settling test", pMBDyn_CE_CEModel);
	gl_window.SetCamera(ChVector<>(-0.1, 0.1, 0.6), ChVector<>(-0.10, 0.1, 0.4), ChVector<>(0, 1.0, 0), 0.05f);
	gl_window.SetRenderMode(opengl::WIREFRAME);
#endif*/

	// ---------------
	// Simulate system
	// ---------------
	/*double time_step = 1e-4;
	double time_end = 5.0;
	double time_out = 0.1;
	bool output = false;
	double time = 0.0;
	//step: creat a file for data
	const std::string out_dir = GetChronoOutputPath() + "DEMO_GNUPLOT";
	// create a .dat file with three columns of demo data:
	std::string datafile = out_dir + "/particle.dat";
	ChStreamOutAsciiFile mdatafile(datafile.c_str());
	ChVector<> opos;
	ChVector<> ovel;
	ChVector<> oForce;
	//GetLog() << "Time: " << time << "\n";
	int out_i = 10;*/
	/*while (time < time_end)
	{
#ifdef CHRONO_OPENGL
		opengl::ChOpenGLWindow& gl_window = opengl::ChOpenGLWindow::getInstance();

		if (gl_window.Active()) {
			gl_window.DoStepDynamics(time_step);
			gl_window.Render();
			time += time_step;
		}
		else {
			break;
		}
#else
		// Run simulation for specified time
			system.DoStepDynamics(time_step);
			time += time_step;

		
#endif
	

	if (out_i == 10)
		{
			GetLog() << "Time: " << time << "\n";
			opos = body_1->GetPos();
			ovel = body_1->GetPos_dt();
			//oForce = ball->GetContactForce();
			mdatafile << time << "\t"
				<< container->GetPos().x() << "\t"
				<< container->GetPos_dt().x() << "\t"
				<< 0.1+opos.x() << "\t"
				<< ovel.x() << "\t"
				<< container->GetPos().y() << "\t"
				<< container->GetPos_dt().y() << "\n";
			//<< oForce.y() << "\n";
			out_i = 0;
		}
		out_i = out_i + 1;
	

#ifdef CHRONO_OPENGL


#endif
	}*/
	//return pMBDyn_CE_CEModel;
}