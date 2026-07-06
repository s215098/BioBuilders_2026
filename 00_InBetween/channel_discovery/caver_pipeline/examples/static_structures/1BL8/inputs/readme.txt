#Read me file for analysis of 1BL8

Following parameters were set to valeus other than the default ones:

probe_radius 0.6 - to enable identification of even partially closed tunnels
shell_radius 5 - to allow proper construction of molecular surface excluding inner parts of the channel
frame_clustering_threshold 2 - to remove redundant tunnels which occurs more frequently with lower probe_radius
seed 1 - to disable randomization, and thus provide comparable results 
average_surface_frame no - averaging of surfaces is suitable rather for globular proteins than for channels
average_surface_global no - averaging of surfaces is suitable rather for globular proteins than for channels
number_of_approximating_balls 4 - number of approximating balls was reduced from 12 to 4 to allow calculation within 1200mb available for 32-bit Java. If you want to perform calculation using more precise approximation, you will have to use 64-bit version of Java, which is able to allocate more memory.
