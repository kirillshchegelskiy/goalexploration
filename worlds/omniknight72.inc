# The robot with omnidirectional distance sensors

define knightsonar sensor
(
  # define the size of each transducer [xsize ysize zsize] in meters
	size [0.001 0.001 0.001 ] 
  # define the range bounds [min max]
	range [0 10.0]
  # define the angular field of view in degrees
	fov 5
  # define the number of samples spread over the fov
	samples 1
)

define knight_sonar ranger
(
	knightsonar( pose [ 0.1 0.01 0 2.5 ] )
	knightsonar( pose [ 0.1 0.02 0 7.5 ] )
	knightsonar( pose [ 0.1 0.03 0 12.5 ] )
	knightsonar( pose [ 0.1 0.04 0 17.5 ] )
	knightsonar( pose [ 0.1 0.05 0 22.5 ] )
	knightsonar( pose [ 0.1 0.06 0 27.5 ] )
	knightsonar( pose [ 0.1 0.07 0 32.5 ] )
	knightsonar( pose [ 0.1 0.08 0 37.5 ] )
	knightsonar( pose [ 0.1 0.09 0 42.5 ] )
	knightsonar( pose [ 0.09 0.1 0 47.5 ] )
	knightsonar( pose [ 0.08 0.1 0 52.5 ] )
	knightsonar( pose [ 0.07 0.1 0 57.5 ] )
	knightsonar( pose [ 0.06 0.1 0 62.5 ] )
	knightsonar( pose [ 0.05 0.1 0 67.5 ] )
	knightsonar( pose [ 0.04 0.1 0 72.5 ] )
	knightsonar( pose [ 0.03 0.1 0 77.5 ] )
	knightsonar( pose [ 0.02 0.1 0 82.5 ] )
	knightsonar( pose [ 0.01 0.1 0 87.5 ] )
	
	
	knightsonar( pose [ -0.01 0.1 0 92.5] )
	knightsonar( pose [ -0.02 0.1 0 97.5] )
	knightsonar( pose [ -0.03 0.1 0 102.5] )
	knightsonar( pose [ -0.04 0.1 0 107.5] )
	knightsonar( pose [ -0.05 0.1 0 112.5] )
	knightsonar( pose [ -0.06 0.1 0 117.5] )
	knightsonar( pose [ -0.07 0.1 0 122.5] )
	knightsonar( pose [ -0.08 0.1 0 127.5] )
	knightsonar( pose [ -0.09 0.1 0 132.5] )
	knightsonar( pose [ -0.1 0.09 0 137.5] )
	knightsonar( pose [ -0.1 0.08 0 142.5] )
	knightsonar( pose [ -0.1 0.07 0 147.5] )
	knightsonar( pose [ -0.1 0.06 0 152.5] )
	knightsonar( pose [ -0.1 0.05 0 157.5] )
	knightsonar( pose [ -0.1 0.04 0 162.5] )
	knightsonar( pose [ -0.1 0.03 0 167.5] )
	knightsonar( pose [ -0.1 0.02 0 172.5] )
	knightsonar( pose [ -0.1 0.01 0 177.5] )
	
	
	knightsonar( pose [ -0.1 -0.01 0 -177.5 ] )
	knightsonar( pose [ -0.1 -0.02 0 -172.5 ] )
	knightsonar( pose [ -0.1 -0.03 0 -167.5 ] )
	knightsonar( pose [ -0.1 -0.04 0 -162.5 ] )
	knightsonar( pose [ -0.1 -0.05 0 -157.5 ] )
	knightsonar( pose [ -0.1 -0.06 0 -152.5 ] )
	knightsonar( pose [ -0.1 -0.07 0 -147.5 ] )
	knightsonar( pose [ -0.1 -0.08 0 -142.5 ] )
	knightsonar( pose [ -0.1 -0.09 0 -137.5 ] )
	knightsonar( pose [ -0.09 -0.1 0 -132.5 ] )
	knightsonar( pose [ -0.08 -0.1 0 -127.5 ] )
	knightsonar( pose [ -0.07 -0.1 0 -122.5 ] )
	knightsonar( pose [ -0.06 -0.1 0 -117.5 ] )
	knightsonar( pose [ -0.05 -0.1 0 -112.5 ] )
	knightsonar( pose [ -0.04 -0.1 0 -107.5 ] )
	knightsonar( pose [ -0.03 -0.1 0 -102.5 ] )
	knightsonar( pose [ -0.02 -0.1 0 -97.5 ] )
	knightsonar( pose [ -0.01 -0.1 0 -92.5 ] )
	
	knightsonar( pose [ 0.01 -0.1 0 -87.5] )
	knightsonar( pose [ 0.02 -0.1 0 -82.5] )
	knightsonar( pose [ 0.03 -0.1 0 -77.5] )
	knightsonar( pose [ 0.04 -0.1 0 -72.5] )
	knightsonar( pose [ 0.05 -0.1 0 -67.5] )
	knightsonar( pose [ 0.06 -0.1 0 -62.5] )
	knightsonar( pose [ 0.07 -0.1 0 -57.5] )
	knightsonar( pose [ 0.08 -0.1 0 -52.5] )
	knightsonar( pose [ 0.09 -0.1 0 -47.5] )
	knightsonar( pose [ 0.1 -0.09 0 -42.5] )
	knightsonar( pose [ 0.1 -0.08 0 -37.5] )
	knightsonar( pose [ 0.1 -0.07 0 -32.5] )
	knightsonar( pose [ 0.1 -0.06 0 -27.5] )
	knightsonar( pose [ 0.1 -0.05 0 -22.5] )
	knightsonar( pose [ 0.1 -0.04 0 -17.5] )
	knightsonar( pose [ 0.1 -0.03 0 -12.5] )
	knightsonar( pose [ 0.1 -0.02 0 -7.5] )
	knightsonar( pose [ 0.1 -0.01 0 -2.5] )
)

define knight_base position 
(
  color "blue"			# Default color.
  drive "diff"		 	# Differential steering model.
  gui_nose 1                  	# Draw a nose on the robot so we can see which way it points
  obstacle_return 1           	# Can hit things.
  ranger_return 0.5            	# reflects sonar beams
  blob_return 1               	# Seen by blobfinders  
  fiducial_return 1           	# Seen as "1" fiducial finders

  localization "gps"         	  
  localization_origin [0 0 0 0] 	# Start odometry at (0, 0, 0).  

  # four DOF kinematics limits
  # [ xmin xmax ymin ymax zmin zmax amin amax ]				
  velocity_bounds [-0.5 0.5 0 0 0 0 -90.0 90.0 ]					
  acceleration_bounds [-0.5 0.5 0 0 0 0 -90 90.0 ]
)

define knight_base_no_sonar knight_base
(
  size [0.44 0.38 0.22]

  origin [0 0 0 0]

  mass 10.0 
 )


define omniknight_base knight_base_no_sonar
(
  knight_sonar( pose [0 0 -0.03 0] ) 
)

define omniknight omniknight_base
(  
  block( 
    points 8
    point[0] [-0.2 0.12]
    point[1] [-0.2 -0.12]
    point[2] [-0.12 -0.2555 ]
    point[3] [0.12 -0.2555]
    point[4] [0.2 -0.12]
    point[5] [0.2 0.12]
    point[6] [0.12 0.2555]
    point[7] [-0.12 0.2555]
    z [0 0.22]
  )
)
