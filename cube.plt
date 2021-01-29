unset key
set border 0
unset tics
unset colorbox
set view 58,27,2
set view equal xyz
load 'cube.fct'

set pm3d depthorder hidden3d
set pm3d implicit
unset hidden3d

set lmargin 2
set rmargin 0
set bmargin 0
set tmargin 0
# get cube positions from file
add_cube(x,y,z,c) = sprintf('cube(%f,%f,%f,%i) w l ls %i,',x,y,z,c,c)
CMD = ''
stats 'cube_positions.txt' u 1:(CMD = CMD.add_cube($1,$2,$3,$4)) nooutput
set xrange [0:13.6]
set yrange [0:11.8]
set zrange [0:13]
CMD = 'splot '.CMD.'1/0 w l ls 2'
# plot cubes
eval(CMD)