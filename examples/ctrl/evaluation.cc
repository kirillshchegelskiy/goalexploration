#include "stage.hh"
#include <typeinfo>
#include <iostream>
#include <typeinfo>
#include <math.h>
#include <fstream>
#include <sstream>
#include <string>
#include <stdexcept>
#include <cmath>
using namespace Stg;

static const double cruisespeed = 0.4; 
static const double avoidspeed = 0.05; 
static const double avoidturn = 0.5;
static const double minfrontdistance = 0.5; // 0.6  
static const bool verbose = false;
static const double stopdist = 1.0;
static const int avoidduration = 12;

static const int placecells = 100;
static const int indim = 12;
//static const int lweightsnum = placecells*(placecells-1)/2;
static const int gridnum = 127;

static const double beta = 0.0;

bool flagfinal = false;

typedef struct
{
	ModelPosition* pos;
	ModelRanger* sonar;
	double  w[placecells][indim];
	double w_lat[placecells][placecells];
	double a[placecells];
	double a_old[placecells];
	double locs[gridnum][gridnum][2];
	int cnti;
	int cntj;
	double winners[gridnum][gridnum];
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

class BadConversion : public std::runtime_error {
public:
  BadConversion(const std::string& s)
    : std::runtime_error(s)
    { }
};
inline double convertToDouble(const std::string& s)
{
  std::istringstream i(s);
  double x;
  if (!(i >> x))
    throw BadConversion("convertToDouble(\"" + s + "\")");
  return x;
}

// Stage calls this when the model starts up
extern "C" int Init( Model* mod, CtrlArgs* args )
{

	robot_t* robot = new robot_t;
	
	double* norms = new double[placecells]();
	
	robot->cnti=0;
	robot->cntj=0;
	
	std::ifstream f("weights.txt");
	std::string line;
	
	if (f.is_open())
	{
	for (int i=0; i<placecells; i++)
	{
		robot->a_old[i] = 0.0;
		for (int j=0; j<indim; j++)
		{
			if (getline (f, line)) robot->w[i][j] = convertToDouble(line);
			//robot->w[i][j] = static_cast<double>(rand() % 1000 - 500)/500;
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
			if (getline (f, line)) robot->w_lat[i][j] = convertToDouble(line);
			}
		}
	f.close();
	std::cout << "\n\nWeights successfully read from file!\n";
	}
	else std::cout << "Could not open weights file!\n";
	
	
	double* grid1 = new double[gridnum+2]();
	for (int i=0; i<gridnum+2; i++)
	{
		grid1[i] = 16.0*i/double(gridnum+1);
		//std::cout <<"grid1["<<i<<"] = "<<grid1[i]<<"\n";
		}
		
	for(int i=0; i<gridnum; i++)
	{	
		for(int j=0; j<gridnum; j++)
		{
			robot->locs[i][j][0] = grid1[i+1]-8.0;
			robot->locs[i][j][1] = grid1[j+1]-8.0;
			//std::cout <<"locs["<<i<<"]["<<j<<"] = ("<<robot->locs[i][j][0]<<", "<<robot->locs[i][j][1]<<")\n";
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
	Pose pose = robot->pos->GetPose();
	//printf( "Pose: [%.2f %.2f %.2f %.2f]\n",  pose.x, pose.y, pose.z, pose.a );
	//std::cout << "absx = " << std::abs(pose.x-robot->locs[robot->cnti][robot->cntj][0]) <<"\n";
	//std::cout << "absy = " << std::abs(pose.y-robot->locs[robot->cnti][robot->cntj][1]) <<"\n";
	
		//Calculate winners at desired locations
	if(std::abs(pose.x-robot->locs[robot->cnti][robot->cntj][0])<0.035 && std::abs(pose.y-robot->locs[robot->cnti][robot->cntj][1])<0.035)
	{
		//std::cout << "In Position\n";
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
		ranges_normd[i] = ranges[i]/(10.0*sqrt(indim)); // normalize input distances array
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
	
	int winner = std::distance(robot->a, std::max_element(robot->a, robot->a + placecells)); // get winner cell
	//std::cout << "Winner cell # at locs["<<robot->cnti<<"]["<<robot->cntj<<"] : " << winner << "\n";
	//robot->act_old = robot->a[winner];
	robot->winners[robot->cnti][robot->cntj]=winner; // store winner cell
	
		//At final location load all data to output file
	if(robot->cnti == gridnum-1 && robot->cntj == gridnum-1)
	{
		if(!flagfinal)
		{
			std::cout<<"Writing winners to file!\n";
			std::ofstream of;
			of.open("winners.txt");
			for (int i=0; i<gridnum; i++)
			{
				for (int j=0; j<gridnum; j++)
				{
					of << robot->winners[i][j] << "\n";
					}
				}
			of.close();
		flagfinal=true;
		}
	}
	
	for (int i=0; i<placecells; i++)
	{
		robot->a_old[i] = robot->a[i];
	}
		}
  return 0; // run again
}

int PositionUpdate( Model* mod, robot_t* robot )
{
  Pose pose = robot->pos->GetPose();

  //printf( "Pose: [%.2f %.2f %.2f %.2f]\n",  pose.x, pose.y, pose.z, pose.a );
  
  
  //std::cout <<"Current goal["<<robot->cnti<<"]["<<robot->cntj<<"]: ("<<robot->locs[robot->cnti][robot->cntj][0]<<", "<<robot->locs[robot->cnti][robot->cntj][1]<<")\n";
  //std::cout << "absx = " << std::abs(pose.x-robot->locs[robot->cnti][robot->cntj][0]) <<"\n";
  
	//Go through list of desired locations and visit them all one by one
  robot->pos->GoTo(robot->locs[robot->cnti][robot->cntj][0],robot->locs[robot->cnti][robot->cntj][1],1.57);
  if(std::abs(pose.x-robot->locs[robot->cnti][robot->cntj][0])<0.03 && std::abs(pose.y-robot->locs[robot->cnti][robot->cntj][1])<0.03)
  {
	  //std::cout<<"Reached locs["<<robot->cnti<<"]["<<robot->cntj<<"]\n";
	  if (robot->cntj == gridnum-1  && robot->cnti<gridnum-1)
	  {
		  robot->cnti++;
		  robot->cntj=0;
	  }
	  else if(robot->cntj<gridnum-1)
	  {
		  robot->cntj++;
	  }
  }
  
  //double angle = pose.a*180/3.1415926;
  //std::cout << angle << "\n";

  return 0; // run again
}
