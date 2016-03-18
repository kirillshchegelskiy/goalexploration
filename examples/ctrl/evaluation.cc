#include "stage.hh"
#include <typeinfo>
#include <iostream>
#include <typeinfo>
#include <math.h>
using namespace Stg;

static const double cruisespeed = 0.4; 
static const double avoidspeed = 0.05; 
static const double avoidturn = 0.5;
static const double minfrontdistance = 0.5; // 0.6  
static const bool verbose = false;
static const double stopdist = 1.0;
static const int avoidduration = 12;

static const int placecells = 10;
static const int indim = 12;
//static const int lweightsnum = placecells*(placecells-1)/2;

static const double beta = 0.5;

typedef struct
{
	ModelPosition* pos;
	ModelRanger* sonar;
	double  w[placecells][indim];
	double w_lat[placecells][placecells];
	double a[placecells];
	double a_old[placecells];
	} robot_t;

int SonarUpdate( Model* mod, robot_t* robot );
int PositionUpdate( Model* mod, robot_t* robot );


double l2_norm(double const* u, int n) 
{
    double accum = 0.;
    for (int i = 0; i < n; ++i) 
    {
        accum += u[i] * u[i];
		}
    return sqrt(accum);
}


// Stage calls this when the model starts up
extern "C" int Init( Model* mod, CtrlArgs* args )
{

	robot_t* robot = new robot_t;
	
	double* norms = new double[placecells]();
	
	 // Initializing weights randomly
	
	for (int i=0; i<placecells; i++)
	{
		robot->a_old[i] = 0.0;
		for (int j=0; j<indim; j++)
		{
			robot->w[i][j] = static_cast<double>(rand() % 1000 - 500)/500;
			//std::cout << "w[" << i << "][" << j << "] = " << robot->w[i][j] << "\n";
			norms[i] += robot->w[i][j] * robot->w[i][j];
			}
		norms[i] = sqrt(norms[i]);
		}
	
	for (int i=0; i<placecells; i++)
	{
		for (int j=0; j<indim; j++)
		{
			robot->w[i][j] = robot->w[i][j]/norms[i];
			}
		
		}	
	
	
	for (int i=0; i<placecells; i++)
	{
		for (int j=0; j<placecells; j++)
		{
			robot->w_lat[i][j] = 0.0;
			}
		}
	
	robot->pos = (ModelPosition*)mod;
	robot->pos->AddCallback( Model::CB_UPDATE, (model_callback_t)PositionUpdate, robot );
	robot->pos->Subscribe(); // starts the position updates

	robot->sonar = (ModelRanger*)mod->GetChild( "ranger:0" );
	robot->sonar->AddCallback( Model::CB_UPDATE, (model_callback_t)SonarUpdate, robot );
	robot->sonar->Subscribe();
 
	return 0; 
}

int SonarUpdate( Model* mod, robot_t* robot )
{
	const std::vector<ModelRanger::Sensor>& sensors = robot->sonar->GetSensors();
	size_t sensor_count = sensors.size();
	//std::cout << "sensor_count type: " << typeid(sensor_count).name();
	
	double* ranges = new double[1]();
	
	int scount = static_cast<int>(sensor_count);
	//std::cout << scount << "\n";
	
	Pose pose = robot->pos->GetPose();
	double deg = pose.a*180/3.1415926;
	//std::cout << deg << "\n";
	int lfront = 0;
	int rfront = 0;
	
	if(scount>0)
	{
		ranges = new double[scount]();
		if (deg>=0) deg = deg+15;
		if (deg<0) deg = deg-15;
		
		for(int i=0; i<scount; i++)
		{
			//store all distances from sonar sensors in one array with orientation compensation
			if (deg>=0) 
				{
					lfront = (int)(deg/30) % 12;
					rfront = (11 + (int)(deg/30)) % 12;
					ranges[(i+(int)(deg/30)) % 12] = sensors[i].ranges[0]; 
				}
			else if (deg<0) 
				{
					lfront = (12-(int)(-deg/30)) % 12;
					rfront = (11+12-(int)(-deg/30)) % 12;
					ranges[(i+12-(int)(-deg/30)) % 12] = sensors[i].ranges[0];
				}
			
			}
		}
	
	//std::cout<<"rfront = "<<rfront<<", lfront = "<<lfront<<"\n";
		
	double* ranges_normd = new double[scount]();
	
	for (int i=0; i<scount; i++)
	{ 
		ranges_normd[i] = ranges[i]/(10.0*sqrt(indim));//normalize input distances array
		//std::cout<<"i = "<<i<<", range_normalized[i] = "<<ranges_normd[i]<<"\n";
		}	
	
	double tmp = 0;
	double tmpc = 0;
	
	for (int i=0; i<placecells; i++)
	{
		for (int j=0; j<indim; j++)
		{
			tmp += (robot->w[i][j])*ranges_normd[j];
			}
		for (int k=0; k<placecells; k++)
		{
			tmpc += (robot->w_lat[i][k])*(robot->a_old[k]);
			}
		//std::cout << "tmp[i] = " << tmp << "\n";
		//robot->a[i] = tanh(tmp + beta*robot->act_old*robot->w_lat[i][robot->winner_old]); // calculate PCs activation
		robot->a[i] = tanh(tmp + beta*tmpc); // calculate PCs activation with matrix multiplication
		tmp = 0;
		tmpc = 0;
		
		//std::cout << "a[" << i << "] = " << robot->a[i] << "\n";
		}
	
	int winner = std::distance(robot->a, std::max_element(robot->a, robot->a + placecells)); //get winner cell
	std::cout << "Winner cell # = " << winner << "\n";
	//robot->act_old = robot->a[winner];
	
	for (int i=0; i<placecells; i++)
	{
		robot->a_old[i] = robot->a[i];
	}
			
  return 0; // run again
}

int PositionUpdate( Model* mod, robot_t* robot )
{
  Pose pose = robot->pos->GetPose();

  printf( "Pose: [%.2f %.2f %.2f %.2f]\n",  pose.x, pose.y, pose.z, pose.a );
  
  robot->pos->GoTo(-5.0,-3.0,0.0);
  
  //double angle = pose.a*180/3.1415926;
  //std::cout << angle << "\n";

  return 0; // run again
}
