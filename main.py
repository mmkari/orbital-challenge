'''
Created on 12.5.2016

@author: Miko Kari
'''
import sys, math
import numpy as np

# Get minimum distance from A to B using Dijkstra's algorithm.
def get_min_distance(graph, start, goal):
	best_distance = {id: None for id in graph.keys()}	# dictionary of current distances to frontier nodes
	previous_node = {}	# dictionary of previous nodes in best paths
	closed = []		# indicates closed nodes
	# initialize with start node
	current = start
	current_distance = 0
	best_distance[current] = current_distance
	# Start search
	while best_distance:
		if current == goal: break	# reached the goal
		
		# add new nodes to the frontier (defined by nodes with specified distance)
		for neighbor, distance in graph[current].items():
			# deal with this node if not yet closed
			if neighbor in closed: continue
			# calculate new distance from start to this node
			new_distance = current_distance + distance
			# replace or update previous distance if appropriate
			if best_distance[neighbor] is None or best_distance[neighbor] > new_distance:
				best_distance[neighbor] = new_distance
				previous_node[neighbor] = current	# mark source node for this currently best distance
		
		# remove current node from frontier and mark it closed
		closed.append(current)
		del best_distance[current]
				
		# determine best node in frontier to be the next current node
		try:
			# find frontier node at minimum distance from start node
			min_distance = min([value for value in best_distance.values() if not value == None])
			# get first item with that minimum value (if multiple with draw)
			current, current_distance = [(k,v) for k,v in best_distance.items() if v==min_distance][0]
		except:
			print('No path found!')
			break
	# end while loop, return dictionary with best transitions
	return previous_node
	

def get_graph(cartesian_data, satellite_ids):
	# Determine point pairs without ray intersection with Earth
	mat = np.ones((len(cartesian_data),len(cartesian_data)))*-9
	for i in range(len(cartesian_data)-1):
		for j in range(i+1,len(cartesian_data)):
			if check_ray(cartesian_data[i],cartesian_data[j],r) == 1:
				# the two points are connected
				d = np.linalg.norm(cartesian_data[i]-cartesian_data[j])
				mat[i][j] = d
				mat[j][i] = d
	with open('mat.txt','w') as outf:
		outf.writelines(','.join(str(item) for item in inner_list)+"\n" for inner_list in mat)
	
	# Build the graph as a dictionary of dictionaries containing neighbors and their distances
	graph = {}
	for i in range(len(satellite_ids)):
		neighbors_dict = {satellite_ids[index]: dist for index, dist in enumerate(mat[i]) if not dist == -9}
		graph[satellite_ids[i]] = neighbors_dict
	return graph

# Checks whether there's a line-of-sight b/w points p1 and p2
def check_ray(p1,p2, radius):
	epsilon = 0.05	# error margin in comparisons of floats
	rounding_digits = 2	# round floats to nearest 10 meters before comparisons
	unit_vector = (p2-p1)/np.linalg.norm(p2-p1)
	discriminant = math.pow(np.dot(unit_vector,p1),2) - math.pow(np.linalg.norm(p1),2) + radius**2
	
	# test for intersection between line and sphere
	if round(discriminant) < 0 - epsilon:	# make sure it really is negative
		# line did not intersect sphere, neither can line segment
		return 1	# no intersection
	else:
		# Still not clear whether SEGMENT and sphere have intersected
		# Intersections given by |p1|**2 - r**2 + 2d*p1*i + d**2 = 0, where i is the unit vector
		# check whether intersection occurs on line segment between p1 and p2 by solving for d.
		
		# get distance between p1 and p2, this is the upper limit for intersection d-value in x = p1+d(p2-p1)/|p2-p1|
		d_limit = round(np.linalg.norm(p1-p2),rounding_digits)
		# solve above quadratic equation for d: d = (-b +- sqrt(b**2 - 4a*c)) / 2a
		b = np.dot(unit_vector,p1)
		s = math.sqrt(discriminant)
		d1 = round(-b + s, rounding_digits)
		d2 = round(-b - s, rounding_digits)
		
		# There is intersection only if both d1 and d2 are within [0,d_limit]
		if 0 <= d1 and d1 <= d_limit and 0 <= d2 and d2 <= d_limit:
			return 0	# Intersection
		else:
			return 1	# no intersection

# Converts spherical coordinates to Cartesian coordinates
def spherical_to_cartesian(list,radius):
	# convert degrees to radians first
	lat = math.radians(float(list[0]))
	long = math.radians(float(list[1]))
	alt = float(list[2])
	k = alt + radius	# distance from center of sphere
	return np.array((k*math.cos(lat)*math.cos(long), k*math.cos(lat)*math.sin(long), k*math.sin(lat)))

if __name__ == '__main__':
	# read name of data file given as command line parameter
	data_file = sys.argv[1]	
	r = 6371	# radius of Earth
	
	# Read data from file
	with open(data_file,'r') as fin:
		# read first line containing SEED
		fin.readline()
		# read data lines
		line = fin.readline()
		satellite_ids = []
		spherical_data = []
		
		while 'ROUTE' not in line:
			# append new satellite info to list
			parts = line.split(',')	# gives: ID, latitude, longitude, altitude
			satellite_ids.append(parts[0])
			spherical_data.append(parts[1:])
			line = fin.readline()
		# handle ROUTE line
		route_data = line.split(',')
		start_point = route_data[1:3]
		start_point.append(0)	# add altitude field
		end_point = route_data[3:5]
		end_point.append(0)	# add altitude field
	
	# convert start and end points to Cartesian coordinates
	start_point = spherical_to_cartesian(start_point, r)
	end_point = spherical_to_cartesian(end_point, r)
	# Convert satellites' spherical coordinates to Cartesian
	cartesian_data = []
	for item in spherical_data:
		cartesian_data.append(spherical_to_cartesian(item,r))
	# add start and end points to the list of Cartesian coordinates
	cartesian_data.extend([start_point,end_point])
	# add start and end point IDs 'A' and 'B' to the list of satellite IDs
	satellite_ids.extend(['A','B'])
	
	# Determine graph formed by satellites and start and end points, accounting for lines-of-sight
	graph = get_graph(cartesian_data, satellite_ids)
	
	# Find shortest path in resulting graph (minimum length, minimum number of edges)
	# Dijkstra's algorithm
	start = 'A'
	goal = 'B'
	source = get_min_distance(graph, start, goal)
		
	# determine path from start to goal, given one was found
	if 'B' in source:
		# Path was found to goal, trace it back to start
		current = goal
		path = [current]
		while current != start:
			current = source[current]
			path.append(current)
		path.reverse()
		# save path to file
		with open('route.txt','w') as f:
			path = ','.join(str(item) for item in path[1:-1])
			f.write(path)
			print('Path stored: %s' % path)