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
#include <stdio.h>
#include <limits.h>
#define DBL_MAX 1.7976931348623158e+308

using namespace Stg;

static const double cruisespeed = 0.4; 
static const double avoidspeed = 0.05; 
static const double avoidturn = 0.5;
static const double minfrontdistance = 0.2; // 0.6  
static const bool verbose = false;
static const double stopdist = 0.5;
static const int avoidduration = 12;
static const int minforward = 50;
static const int maxforward = 100;
static const int momentum = 30;

static const int placecells = 50;
static const int indim = 36;

static const double beta = 0.0;

static const int maxruns = 1000;

typedef struct
{
	ModelPosition* pos;
	ModelRanger* sonar;	
	int avoidcount, randcount, forwardcount, turncount;
	int forwarddist, turndist;
	double  w[placecells][indim];
	double w_lat[placecells][placecells][2];
	double a[placecells];
	double a_old[placecells];
	int itrs;
	int goals[placecells];
	int currgoal;
	int ttlgoals;
	int inter[placecells];
	int localgoal;
	int localind;
	bool turnedtogoal;
	bool moving;
	int currsrc;
	int goalrunnmbr;
	int itrstotal[maxruns];
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

int minDistance(double dist[], bool sptSet[], int ttlgoals)
{
   double min = DBL_MAX;
   int min_index;
 
   for (int v = 0; v < ttlgoals; v++)
     if (sptSet[v] == false && dist[v] <= min)
         min = dist[v], min_index = v;
 
   return min_index;
}

void dijkstra(double w_lat[placecells][placecells][2], int goals[placecells], int ttlgoals, int src, int currgoal, int inter[placecells])
{
	
	if (src==currgoal)
		{
			for (int i=0; i<ttlgoals; i++) inter[i] = -1;
		}
	else {
	double dist[ttlgoals];
	bool sptSet[ttlgoals];
	double graph[ttlgoals][ttlgoals];
	int parent[ttlgoals];
	
	//std::cout<<"Started Dijkstra\n";
	
	for (int i=0; i<ttlgoals; i++)
	{
		for (int j=0; j<ttlgoals; j++)
		{
			graph[i][j] = 1/w_lat[goals[i]][goals[j]][0];
		}
	}
	
	//std::cout<<"Stored 1/w_lat\n";
	
	for (int i=0; i<ttlgoals; i++) dist[i]=DBL_MAX, sptSet[i]=false, inter[i]=-1, parent[i]=-1;
	
	int src_ind = std::distance(goals, std::find(goals, goals + ttlgoals, src));
	//std::cout<<"src = "<<src<<", src_ind = "<<src_ind<<"\n";
	
	dist[src_ind] = 0;
	//dijkstra(robot->w_lat, robot->goals, robot->ttlgoals, robot->goals[0], robot->currgoal, robot->inter);
	//std::cout<<"Dijkstra initialization complete\n";
	
	for (int count = 0; count < ttlgoals-1; count++)
	{
       // Pick the minimum distance vertex from the set of vertices not
       // yet processed. u is always equal to src in first iteration.
        //std::cout<<"Started Dijkstra step "<<count<<"\n";
       
		int u = minDistance(dist, sptSet, ttlgoals);
		//std::cout<<"Found minDistance for vertex "<<u<<"("<<goals[u]<<")\n";
       // Mark the picked vertex as processed
		sptSet[u] = true;
 
       // Update dist value of the adjacent vertices of the picked vertex.
		for (int v = 0; v < ttlgoals; v++)
		{
         // Update dist[v] only if is not in sptSet, there is an edge from 
         // u to v, and total weight of path from src to  v through u is 
         // smaller than current value of dist[v]
			if (!sptSet[v] && graph[v][u] && dist[u] != DBL_MAX && dist[u]+graph[v][u] < dist[v]) 
			{
				//std::cout<<"Updating vertex "<<v<<"("<<goals[v]<<") with "<<graph[v][u]<<"\n";
				parent[v] = u;
				dist[v] = dist[u] + graph[v][u];
			}
		}
	}
	//std::cout<<"Dijkstra complete!\n";
	int tmp = std::distance(goals, std::find(goals, goals + ttlgoals, currgoal));;
	for (int i=0; i<ttlgoals-1; i++)
	{
		
		tmp = parent[tmp];
		if (tmp==src_ind) break;
		inter[i] = goals[tmp];
	}
	}
}


// Stage calls this when the model starts up
extern "C" int Init( Model* mod, CtrlArgs* args )
{

	robot_t* robot = new robot_t;
	
	double* norms = new double[placecells]();
	
	robot->itrs=0;
	robot->ttlgoals=0;
	
	robot->avoidcount = 0;
	robot->randcount = 0;
	robot->forwardcount = 0;
	robot->turncount = 0;
	robot->forwarddist = rand()%(maxforward-minforward) + minforward;
	robot->turndist = 2*(rand()%momentum) - momentum + 1;
	
	robot->turnedtogoal = false;
	robot->currsrc = 0;
	robot->goalrunnmbr = 0;
	
	std::ifstream f("weights.txt");
	std::ifstream fg("goals.txt");
	std::string line;
	std::string lineg;
	
	if (f.is_open())
	{
	for (int i=0; i<placecells; i++)
	{
		robot->a_old[i] = 0.0;
		robot->goals[i] = -1;
		robot->inter[i] = -1;
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
		if (getline (fg, lineg)) robot->goals[i] = convertToDouble(lineg);
		//std::cout << "goals[" << i << "] = " << robot->goals[i] << "\n";
		for (int j=0; j<placecells; j++)
		{
			if (getline (f, line)) robot->w_lat[i][j][0] = convertToDouble(line);
			if (getline (f, line)) robot->w_lat[i][j][1] = convertToDouble(line);
			}
		}
	f.close();
	std::cout << "\n\nWeights and goals successfully read from files!\n";
	}
	else std::cout << "Could not open weights/goals file!\n";
	
	for (int i=0; i<placecells; i++)
	{
		if (robot->goals[i]!=-1) robot->ttlgoals+=1;
	}
	//std::cout << "Total number of goals = " << robot->ttlgoals << "\n";
	robot->currgoal = robot->goals[rand()%(robot->ttlgoals)];
	//std::cout << "Next goal = " << robot->currgoal << "\n";
	/*
	for (int i=0; i<robot->ttlgoals; i++)
	{
		for (int j=0; j<robot->ttlgoals; j++)
		{
			std::cout << "1/w_lat["<<robot->goals[i]<<"]["<<robot->goals[j]<<"][0] = "<<1/robot->w_lat[robot->goals[i]][robot->goals[j]][0]<<", "<<"w_lat["<<robot->goals[i]<<"]["<<robot->goals[j]<<"][1] = "<<robot->w_lat[robot->goals[i]][robot->goals[j]][1]<<"\n";
		}
	}
	*/
	/*for (int i=0; i<robot->ttlgoals; i++)
	{
		for (int j=0; j<robot->ttlgoals; j++)
		{
			std::cout << "w_lat["<<robot->goals[i]<<"]["<<robot->goals[j]<<"][1] = "<<robot->w_lat[robot->goals[i]][robot->goals[j]][1]<<"\n";
		}
	}*/
	
	dijkstra(robot->w_lat, robot->goals, robot->ttlgoals, robot->goals[0], robot->currgoal, robot->inter);
	//std::cout << "Intermediate goals: \n";
	/*for (int i=0; i<robot->ttlgoals; i++)
	{
		if (robot->inter[i]!=-1) std::cout << robot->inter[i] << " ";
	}*/
	//std::cout<<"\n";
	for (int i=robot->ttlgoals-1; i>-1; i--)
	{
		if (robot->inter[i]!=-1) 
		{
			robot->localgoal = robot->inter[i];
			robot->localind = i;
			break;
		}
		if (i==0) robot->localgoal = robot->currgoal;
	}
	//std::cout << "Local goal: " << robot->localgoal << "\n";
	
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
		if (deg>=0) deg = deg+360/(indim*2);
		if (deg<0) deg = deg-360/(indim*2);
		
		int avoid_param = 0;
		
		if (indim==36) avoid_param = 2;
		if (indim==72) avoid_param = 4;
		
		for(int i=0; i<scount; i++)
		{
			//store all distances from sonar sensors in one array with orientation compensation
			if (deg>=0) 
				{
					lfront = (int)(deg/(360/indim) + avoid_param) % indim;
					rfront = (indim-1 + (int)(deg/(360/indim)) - avoid_param) % indim;
					ranges[(i+(int)(deg/(360/indim))) % indim] = sensors[i].ranges[0]; 
				}
			else if (deg<0) 
				{
					lfront = (indim-(int)(-deg/(360/indim)) + avoid_param) % indim;
					rfront = (indim-1+indim-(int)(-deg/(360/indim)) - avoid_param) % indim;
					ranges[(i+indim-(int)(-deg/(360/indim))) % indim] = sensors[i].ranges[0];
				}
			
			}
		}
	
	if(robot->goalrunnmbr < maxruns){
	
	
	double* ranges_bit = new double[scount]();
	for (int i=0; i<scount; i++)
	{
		if (ranges[i]>3.0) ranges_bit[i]=0.0;
		else ranges_bit[i]=1.0;
		}
	double* ranges_normd = new double[scount]();
	
	for (int i=0; i<scount; i++)
	{ 
		ranges_normd[i] = ranges[i]/l2_norm(ranges, scount); //classic l2 normalization
		//ranges_normd[i] = ranges[i]/(10.0*sqrt(indim)); // normalize input distances array
		
		//ranges_normd[i] = ranges_bit[i]/l2_norm(ranges_bit, scount); //bit input normalization
		
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
			tmpc += (robot->w_lat[i][k][0])*(robot->a_old[k]);
			}
		//std::cout << "tmp[i] = " << tmp << "\n";
		//robot->a[i] = tanh(tmp + beta*robot->act_old*robot->w_lat[i][robot->winner_old][0]); // calculate PCs activation
		robot->a[i] = tanh(tmp + beta*tmpc); // calculate PCs activation with matrix multiplication
		tmp = 0;
		tmpc = 0;
		
		//std::cout << "a[" << i << "] = " << robot->a[i] << "\n";
		}
	
	int winner = std::distance(robot->a, std::max_element(robot->a, robot->a + placecells)); // get winner cell
	
	if (robot->itrs==0) robot->currsrc = winner;
	
	if (winner==robot->currgoal) 
	{
		//std::cout << "\nGoal " << winner << " reached! It took " << robot->itrs << " iterations\n";
		
		if (robot->itrs > 1) 
		{
			std::cout << "Run " << robot->goalrunnmbr+1 << " complete \n";
			robot->itrstotal[robot->goalrunnmbr] = robot->itrs;
			robot->goalrunnmbr += 1;
		}
		
		robot->itrs = 0;
		robot->currgoal = robot->goals[rand()%(robot->ttlgoals)];
		//std::cout << "Next goal = " << robot->currgoal << "\n";
		
		//if(winner==robot->currgoal) std::cout << "WINNER == NEXT GOAL\n";
		dijkstra(robot->w_lat, robot->goals, robot->ttlgoals, winner, robot->currgoal, robot->inter);
		//std::cout << "Intermediate goals: ";
		/*sfor (int i=0; i<robot->ttlgoals; i++)
		{
			if (robot->inter[i]!=-1) std::cout << robot->inter[i] << " ";
		}*/
		for (int i=robot->ttlgoals-1; i>-1; i--)
		{
			if (robot->inter[i]!=-1) 
			{
				robot->localgoal = robot->inter[i];
				robot->localind = i;
				break;
			}
			if (i==0) robot->localgoal = robot->currgoal;
		}
		//std::cout << "Local goal: " << robot->localgoal << ", localind = " << robot->localind << "\n";
		robot->turnedtogoal = false;
		robot->currsrc = winner;
	}
	robot->itrs+=1;
		
	if (winner!=robot->currgoal && winner==robot->localgoal && robot->localind!=0)
	{
		//std::cout << "Local goal " << winner << " reached!\n";
		for (int i=robot->localind-1; i>-1; i--)
		{
			if (robot->inter[i]!=-1) 
			{
				robot->localgoal = robot->inter[i];
				robot->localind = i;
				break;
			}
		}
		//std::cout << "Next local goal: " << robot->localgoal << ", localind = " << robot->localind << "\n";
		robot->turnedtogoal = false;
		robot->currsrc = winner;
	}
	
	if (winner!=robot->currgoal && winner==robot->localgoal && robot->localind==0) 
	{
		//std::cout <<"Local goal "<<robot->localgoal<<" reached! Next local goal is current goal\n";
		robot->localgoal = robot->currgoal;	
		robot->turnedtogoal = false;
		robot->currsrc = winner;
	}
	
	for (int i=0; i<placecells; i++)
	{
		robot->a_old[i] = robot->a[i];
	}
		
		
	if(robot->goalrunnmbr == maxruns)
	{
		std::cout<<"Writing iterations numbers to file!\n";
		std::ofstream of;
		of.open("naviters.txt");
		for (int i=0; i<maxruns; i++)
		{
			of << robot->itrstotal[i] << "\n";
		}
		of.close();
	}
	
	
		
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
      
	if(!robot->turnedtogoal)
	{
		if (robot->moving == true)
		{
			//std::cout << "Stopping ROBOT movement\n";
			robot->pos->SetXSpeed(0.0);
			robot->pos->SetTurnSpeed(0.0);
			robot->moving = false;
		}
		//std::cout << "POSE" << pose.x<<" "<<pose.y<<" "<<pose.a<<"\n";
		//std::cout << "GoTo"<<pose.x<<" "<<pose.y<<" "<<robot->w_lat[robot->localgoal][robot->currsrc][1]<<"\n";
		robot->pos->GoTo(pose.x,pose.y,robot->w_lat[robot->localgoal][robot->currsrc][1]);
		if (std::abs(pose.a - robot->w_lat[robot->localgoal][robot->currsrc][1])<0.1) 
			{
				robot->turnedtogoal = true;
				//std::cout << "ROBOT Turned to local goal "<<robot->localgoal<<"\n";
			}
	}

	if (robot->turnedtogoal) {
	robot->moving=true;
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
	else //issues true random walk
    {
		if( verbose ) puts( "Cruise" );

		robot->avoidcount = 0;
		
		if (robot->forwardcount==robot->forwarddist)
		{
			//robot->forwardcount = 0;
			//robot->forwarddist = rand()%(maxforward-minforward) + minforward;
			robot->pos->SetXSpeed(0);
			if (robot->turndist > 0) 
			{
				robot->pos->SetTurnSpeed(+avoidturn);
				robot->turncount++;
			}
			if (robot->turndist < 0) 
			{
				robot->pos->SetTurnSpeed(-avoidturn);
				robot->turncount--;
			}
			
			if (robot->turncount==robot->turndist)
			{
				robot->forwardcount = 0;
				robot->forwarddist = rand()%(maxforward-minforward) + minforward;
				robot->turncount = 0;
				robot->turndist = 2*(rand()%momentum) - momentum + 1;
			}
		}
		else
		{
			robot->pos->SetXSpeed(cruisespeed);	  
			robot->pos->SetTurnSpeed(0);
			robot->forwardcount++;
		}
	}
	}
  }
  return 0; // run again
}

int PositionUpdate( Model* mod, robot_t* robot )
{
  Pose pose = robot->pos->GetPose();

  //printf( "Pose: [%.2f %.2f %.2f %.2f]\n",  pose.x, pose.y, pose.z, pose.a );

  return 0; // run again
}
