#include "stage.hh"
#include <typeinfo>
#include <iostream>
#include <typeinfo>
#include <math.h>
using namespace Stg;

static const double cruisespeed = 0.4; 
static const double avoidspeed = 0.05; 
static const double avoidturn = 0.5;
static const double minfrontdistance = 1.0; // 0.6  
static const bool verbose = false;
static const double stopdist = 0.3;
static const int avoidduration = 10;

static const int placecells = 10;
static const int indim = 12;

static const double eta1 = 0.045;
static const double alpha = 0.62;
static const double eta2 = 0.5;
static const double delta = 0.998;


typedef struct
{
	ModelPosition* pos;
	ModelRanger* sonar;
	int avoidcount, randcount;
	double  w[placecells][indim];
	//double w_old[placecells][indim];
	//double wlat[placecells][placecells-1];
	double a[placecells];
	int winner_old;
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
 
	robot->avoidcount = 0;
	robot->randcount = 0;
	robot->winner_old = 0;
	
	double* norms = new double[placecells]();
	
	 // Initializing weights randomly
	
	for (int i=0; i<placecells; i++)
	{
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
	
	
	robot->pos = (ModelPosition*)mod;
	if( verbose ) robot->pos->AddCallback( Model::CB_UPDATE, (model_callback_t)PositionUpdate, robot );
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
	
	if(scount>0)
	{
		ranges = new double[scount]();
		
		for(int i=0; i<scount; i++)
		{
			ranges[i] = sensors[i].ranges[0]; //store all distances from sonar sensors in one array
			//std::cout<<"i = "<<i<<", range[i] = "<<ranges[i]<<"\n";
			}
		}
	
	double* ranges_normd = new double[scount]();
	
	for (int i=0; i<scount; i++)
	{
		ranges_normd[i] = ranges[i]/l2_norm(ranges, scount); //normalize input distances array
		//std::cout<<"i = "<<i<<", range_normalized[i] = "<<ranges_normd[i]<<"\n";
		}
	
	
	
	
	
	double tmp = 0;
	
	for (int i=0; i<placecells; i++)
	{
		for (int j=0; j<indim; j++)
		{
			tmp += (robot->w[i][j])*ranges_normd[j];
			
			}
		//std::cout << "tmp[i] = " << tmp << "\n";
		robot->a[i] = tanh(tmp); // calculate PCs activation
		tmp = 0;
		std::cout << "a[" << i << "] = " << robot->a[i] << "\n";
		}
	
	int winner = std::distance(robot->a, std::max_element(robot->a, robot->a + placecells)); //get winner cell
	std::cout << "Winner cell # = " << winner << "\n";
	for (int j=0; j<indim; j++)
	{
		robot->w[winner][j] += eta1*(ranges_normd[j] - robot->w[winner][j]); // Kohonen learning for curr winner
		std::cout << "dw[" << j << "] = " << eta1*(ranges_normd[j] - robot->w[winner][j]) << "\n";
		robot->w[robot->winner_old][j] += alpha*eta1*(ranges_normd[j] - robot->w[robot->winner_old][j]); // Kohonen learning for prev winner
		}
	
	//Normalizing the weights
	
	double* norms = new double[placecells]();
	
	for (int i=0; i<placecells; i++)
	{
		for (int j=0; j<indim; j++)
		{
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
			//std::cout << "w[" << i << "][" << j << "] = " << robot->w[i][j] << "\n";
			}
		}	
	
	robot->winner_old = winner; //saving winner cell number for next iteration
	
	
	
	
	bool obstruction = false;
	bool stop = false;
	
  // find the closest distance to the left and right and check if
  // there's anything in front
	double minleft = 1e6;
	double minright = 1e6;
  
  //printf ("%i \n", sample_count);
  
  
	if( ranges[0] < minfrontdistance || ranges[6] < minfrontdistance)
	{
		if( verbose ) puts( "  obstruction!" );
		obstruction = true;
		}
		
	if( ranges[6] < stopdist || ranges[0] < stopdist )
	{
		if( verbose ) puts( "  stopping!" );
		stop = true;
		}
      
    if( ranges[0] < stopdist ) minleft = std::min( minleft, ranges[0] );
    else if ( ranges[6] < stopdist)	minright = std::min( minright, ranges[6] );
      
	if( verbose ) 
    {
		puts( "" );
		printf( "minleft %.3f \n", minleft );
		printf( "minright %.3f\n ", minright );
		}

	if( obstruction || stop || (robot->avoidcount>0) )
    {
		if( verbose ) printf( "Avoid %d\n", robot->avoidcount );
	  		
		robot->pos->SetXSpeed( stop ? 0.0 : avoidspeed );      
      
      /* once we start avoiding, select a turn direction and stick
	 with it for a few iterations */
		if( robot->avoidcount < 1 )
        {
			if( verbose ) puts( "Avoid START" );
			robot->avoidcount = random() % avoidduration + avoidduration;
			 
			if( minleft < minright  )
			{
				robot->pos->SetTurnSpeed( -avoidturn );
				if( verbose ) printf( "turning right %.2f\n", -avoidturn );
				}
			else
			{
				robot->pos->SetTurnSpeed( +avoidturn );
				if( verbose ) printf( "turning left %2f\n", +avoidturn );
				}
			}
		
      robot->avoidcount--;
    } 
	else
    {
		if( verbose ) puts( "Cruise" );

		robot->avoidcount = 0;
		robot->pos->SetXSpeed( cruisespeed );	  
		robot->pos->SetTurnSpeed(  0 );
		}

  //  if( robot->pos->Stalled() )
  // 	 {
  // 		robot->pos->SetSpeed( 0,0,0 );
  // 		robot->pos->SetTurnSpeed( 0 );
  // }
			
  return 0; // run again
}

int PositionUpdate( Model* mod, robot_t* robot )
{
  Pose pose = robot->pos->GetPose();

  printf( "Pose: [%.2f %.2f %.2f %.2f]\n",
	  pose.x, pose.y, pose.z, pose.a );

  return 0; // run again
}

