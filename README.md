# orbital-challenge

From https://reaktor.com/orbital-challenge/:

>Your task is to create an algorithm that can route phone calls through space across a network of satellites, much like the Iridium satellite constellation. Due to unfortunate circumstances that took place in the launch, our satellite constellation did not turn out as planned, but instead the satellites are scattered randomly across the globe at altitudes between 300-700km. Your algorithm should return the intermediate hops across satellites needed to transmit a signal from a starting point on the ground to an end point (also on the ground) in valid order (e.g. SAT10,SAT22,SAT7). No signal can travel through Earth, so all hops must be made between satellites that have an unobstructed line of sight between each other. It is possible, albeit very unlikely that a working route cannot be found for a given data file. The route you generate does not need to be the optimal one (i.e. least amount of hops or shortest), but our engineers will appreciate such a solution, as it’ll reduce the overall bandwidth needed.

>This CSV formatted randomized data file contains a snapshot with the current position of the satellites in orbit expressed as:

> *ID,latitude,longitude,altitude*

>The last line in the file will contain a route between two random points on the Earth surface expressed as:

> *ROUTE,latitude of starting point,longitude of starting point,latitude of end point,longitude of end point*

>The file will also contain a random seed used to generate it, which you will need to submit for verification purposes along with your answer. To make things simple, Earth is assumed to be perfectly round and its radius is 6371km. All altitudes are expressed as kilometers above the surface.
