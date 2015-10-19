# The robot with omnidirectional distance sensors

define knightsonar sensor
(
  # define the size of each transducer [xsize ysize zsize] in meters
	size [0.01 0.01 0.01 ] 
  # define the range bounds [min max]
	range [0 5.0]
  # define the angular field of view in degrees
	fov 15
  # define the number of samples spread over the fov
	samples 1
)

define knight_sonar ranger
(
	knightsonar( pose [ 0.2 0.1 0 15 ] )
	knightsonar( pose [ 0.1 0.1 0 45] )
	knightsonar( pose [ 0.05 0.1 0 75 ] )
	knightsonar( pose [ -0.05 0.1 0 105 ] )
	knightsonar( pose [ -0.1 0.1 0 135 ] )
	knightsonar( pose [ -0.2 0.1 0 165 ] )
	
	
	knightsonar( pose [ 0.2 -0.1 0 -15 ] )
	knightsonar( pose [ 0.1 -0.1 0 -45] )
	knightsonar( pose [ 0.05 -0.1 0 -75 ] )
	knightsonar( pose [ -0.05 -0.1 0 -105 ] )
	knightsonar( pose [ -0.1 -0.1 0 -135 ] )
	knightsonar( pose [ -0.2 -0.1 0 -165 ] )
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
