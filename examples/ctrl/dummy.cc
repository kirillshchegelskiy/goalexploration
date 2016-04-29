#include "stage.hh"
#include <typeinfo>
#include <iostream>
#include <typeinfo>
#include <math.h>
#include <fstream>
using namespace Stg;

static const double cruisespeed = 0.4; 
static const double avoidspeed = 0.05; 
static const double avoidturn = 0.5;
static const double minfrontdistance = 0.01; // 0.6  
static const bool verbose = false;
static const double stopdist = 0.5;
static const int avoidduration = 12;

static const int placecells = 50;
static const int indim = 12;
//static const int lweightsnum = placecells*(placecells-1)/2;

static const double eta1 = 0.001;
static const double alpha = 0.5;
static const double eta2 = 0.05;
static const double delta = 0.998;
static const double beta = 0.5;

static const int trwindow = 10;
static const double training_criteria = 0.0003;


typedef struct
{
	ModelPosition* pos;
	ModelRanger* sonar;
	int avoidcount, randcount;
	double  w[placecells][indim];
	double w_old[placecells][indim];
	double w_lat[placecells][placecells];
	double a[placecells];
	int winner_old;
	double a_old[placecells];
	bool trained;
	double dw[trwindow];
	bool writeflag;
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
	robot->trained = false;
	robot->writeflag = true;
	//robot->act_old = 0.0;
	
	double* norms = new double[placecells]();
	
	 // Initializing weights randomly
	
	for (int i=0; i<placecells; i++)
	{
		robot->a_old[i] = 0.0;
		for (int j=0; j<indim; j++)
		{
			robot->w[i][j] = static_cast<double>(rand() % 1000)/1000;
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
			robot->w_old[i][j] = 0.0;
			}
		
		}	
	
	
	for (int i=0; i<placecells; i++)
	{
		for (int j=0; j<placecells; j++)
		{
			robot->w_lat[i][j] = 0.0;
			}
		}
	
	for (int i=0; i<trwindow; i++)
	{
		robot->dw[i] = 0.0;
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
	
	double* ranges_bit = new double[scount]();
	for (int i=0; i<scount; i++)
	{
		if (ranges[i]>3.0) ranges_bit[i]=0.0;
		else ranges_bit[i]=1.0;
		}
	
	double* ranges_normd = new double[scount]();
	
	for (int i=0; i<scount; i++)
	{
		//std::cout<<"ranges[" << i << "] = "<<ranges[i]<<"\n";
		ranges_normd[i] = ranges[i]/l2_norm(ranges, scount); //classic l2 normalization
		//ranges_normd[i] = ranges[i]/(10.0*sqrt(indim));//normalize input distances array - weird way w.r.t. max sensor range
		//ranges_normd[i] = ranges[i]/10.0; //normalize input distances array w.r.t. max sensor range
		
		//ranges_normd[i] = ranges_bit[i]/l2_norm(ranges_bit, scount); //bit input normalization
		
		//std::cout<<"i = "<<i<<", range_normalized[i] = "<<ranges_normd[i]<<"\n";
		}
	
	
	if (!robot->trained)
	{
	
	
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
	
	for (int j=0; j<indim; j++)
	{
		robot->w[winner][j] += eta1*(ranges_normd[j] - robot->w[winner][j]); // Kohonen learning for curr winner
		//std::cout << "dw[" << j << "] = " << eta1*(ranges_normd[j] - robot->w[winner][j]) << "\n";
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
		//std::cout << "norms[" << i << "] = " << norms[i] << "\n";
		}
	
	for (int i=0; i<placecells; i++)
	{
		for (int j=0; j<indim; j++)
		{
			robot->w[i][j] = robot->w[i][j]/norms[i];
			//std::cout << "w[" << i << "][" << j << "] = " << robot->w[i][j] << "\n";
			}
		}	
	
	//Registring weights change over a time window
	
	double dw_curr = 0.0;
	
	for (int i=0; i<placecells; i++)
	{
		for (int j=0; j<indim; j++)
		{
			dw_curr += (robot->w[i][j] - robot->w_old[i][j])*(robot->w[i][j] - robot->w_old[i][j]);		
			//std::cout << "dw_curr += " << (robot->w[i][j] - robot->w_old[i][j])*(robot->w[i][j] - robot->w_old[i][j]) << "\n";		
			}
		}
	dw_curr = sqrt(dw_curr);
	//std::cout << "dw_curr = " << dw_curr << "\n";
	
	for (int i=1; i<trwindow; i++)
	{
		robot->dw[i-1] = robot->dw[i];	
		}
	robot->dw[trwindow-1] = dw_curr;
	
	std::cout << "l2_norm of dw: " << l2_norm(robot->dw, trwindow) << "\n";
	
	if (l2_norm(robot->dw, trwindow) < training_criteria) //was 0.00002
	{
		robot->trained = true;
		if (robot->writeflag == true) 
		{
			std::ofstream f;
			f.open("weights.txt");
			for (int i=0; i<placecells; i++)
			{
				for (int j=0; j<indim; j++)
				{
					f << robot->w[i][j] << "\n";
					}
				}
			for (int i=0; i<placecells; i++)
			{
				for (int j=0; j<placecells; j++)
				{
					f << robot->w_lat[i][j] << "\n";
					}
				}
			f.close();
			robot->writeflag = false;
			std::cout << "Wrote to file!\n";
			}
		}
	std::cout << "Training status: " << robot->trained << "\n";
	
	//lateral weight learning
	robot->w_lat[winner][robot->winner_old] += eta2*robot->a[winner]*robot->a[robot->winner_old]*(1-robot->w_lat[winner][robot->winner_old]);
	robot->w_lat[robot->winner_old][winner] += eta2*robot->a[winner]*robot->a[robot->winner_old]*(1-robot->w_lat[winner][robot->winner_old]);
	//std::cout << "w_lat[" << winner << "][" << robot->winner_old << "] = " << robot->w_lat[winner][robot->winner_old] << "\n";
	
	for (int i=0; i<placecells; i++)
	{
		robot->a_old[i] = robot->a[i];
		for (int j=0; j<placecells; j++)
		{
			robot->w_lat[i][j] = delta*robot->w_lat[i][j]; //lateral weight decay
			if(i==j) robot->w_lat[i][j] = 0;
			//std::cout << "w_lat[" << i << "][" << j << "] = " << robot->w_lat[i][j] << "\n";
			}
		}
	//std::cout << "after decay w_lat[" << winner << "][" << robot->winner_old << "] = " << robot->w_lat[winner][robot->winner_old] << "\n";
	
	robot->winner_old = winner; //saving winner cell number for next iteration
	for (int i=0; i<placecells; i++) //saving old weights for next iteration
	{
		for (int j=0; j<indim; j++)
		{
			robot->w_old[i][j] = robot->w[i][j];
			}
		}
	
	
	} //end of training
	
	
	bool obstruction = false;
	bool stop = false;
	
  // find the closest distance to the left and right and check if
  // there's anything in front
	double minleft = 1e6;
	double minright = 1e6;
  
  //printf ("%i \n", sample_count);
  
  
	if( ranges[lfront] < minfrontdistance || ranges[rfront] < minfrontdistance)
	{
		if( verbose ) puts( "  obstruction!" );
		obstruction = true;
		}
		
	if( ranges[lfront] < stopdist || ranges[rfront] < stopdist )
	{
		if( verbose ) puts( "  stopping!" );
		stop = true;
		}
      
    if( ranges[lfront] < stopdist ) minleft = std::min( minleft, ranges[lfront] );
    else if ( ranges[rfront] < stopdist)	minright = std::min( minright, ranges[rfront] );
      
    //printf( "minleft %.3f, lfront %.3f \n", minleft, ranges[lfront] );
	//printf( "minright %.3f, rfront %.3f \n ", minright, ranges[rfront] );  
      
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
  //Pose pose = robot->pos->GetPose();

  //printf( "Pose: [%.2f %.2f %.2f %.2f]\n",  pose.x, pose.y, pose.z, pose.a );
  
  //double angle = pose.a*180/3.1415926;
  //std::cout << angle << "\n";

  return 0; // run again
}

