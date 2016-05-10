import math
import sys

def main():
  f = open('../data/data.txt', 'r')

  lines = f.readlines()
  satellites = read_satellites(lines)
  graph = create_graph(satellites)
  print(find_shortest_path(graph, satellites))

def read_satellites(lines):
  cartesian_coords = {}

  for l in lines[1:-1]:
    s = l.split(',')
    label = s[0]
    polar = float(s[1])
    azim = float(s[2])
    r = float(s[3])
    cartesian_coords[label] = spherical_to_cartesian(polar, azim, r)

  s = lines[-1].split(',')

  cartesian_coords["START"] = spherical_to_cartesian(float(s[1]), float(s[2]), 1e-3)
  cartesian_coords["STOP"] = spherical_to_cartesian(float(s[3]), float(s[4]), 1e-3)

  return cartesian_coords;

def spherical_to_cartesian(polar, azim, r):
  polar_in_rad = math.radians(polar + 90.0)
  azim_in_rad = math.radians(azim + 180.0)
  r_in_earth = r/1000.0 + 6.371

  sin_polar = math.sin(polar_in_rad)
  cos_azim = math.cos(azim_in_rad)
  cos_polar = math.cos(polar_in_rad)
  sin_azim = math.sin(azim_in_rad)

  return (r_in_earth * sin_polar * cos_azim, r_in_earth * sin_polar * sin_azim, r_in_earth * cos_polar)

def create_graph(satellite_list):
  graph = {}
  for key1, value1 in satellite_list.items():
    l = []
    for key2, value2 in satellite_list.items():
      if key1 != key2:
        if visible(value1, value2):
          l.append(key2)
    graph[key1] = l

  return graph;

def sub(v1, v2):
  return (v1[0] - v2[0], v1[1] - v2[1], v1[2] - v2[2]);

def add(v1, v2):
  return (v1[0] + v2[0], v1[1] + v2[1], v1[2] + v2[2]);

def dot(v1, v2):
  return v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2];

def length(v):
  return math.sqrt(dot(v, v))

def distance(v1, v2):
  return length(sub(v2, v1))

def normalize(v):
  l = length(v)
  if l < 1e-15:
    print("Degenerate vector")
    return (1.0, 0.0, 0.0)
  return (v[0] / l, v[1] / l, v[2] / l)

# ray sphere intersection
def visible(v1, v2):
  diff = sub(v2, v1)
  ray = normalize(diff)
  vdotray = dot(ray, v1)
  discriminant = vdotray**2 - dot(v1, v1) + 6.371*6.371
  if discriminant < 0: return True;
  sqrt = math.sqrt(discriminant)
  tplus = -vdotray + sqrt;
  tminus = -vdotray - sqrt;
  if tplus < 0 and tminus < 0: return True;
  distance = length(diff)
  if tplus > distance and tminus > distance: return True;
  return False;  

# djikstra
def find_shortest_path(graph, satellites):
  dist = {}
  prev = {}
  workingSet = []
  for key, val in graph.items():
    dist[key] = 1e30
    prev[key] = "NULL"
    workingSet.append(key)

  dist["START"] = 0
  while len(workingSet) > 0:
    minimum = min([v for k,v in dist.items() if k in workingSet])
    current = [k for k in workingSet if dist[k] == minimum][0]
    workingSet.remove(current)
    if current == "STOP":
      break;
    for neighbor in graph[current]:
      alt = dist[current] + distance(satellites[neighbor], satellites[current])
      if alt < dist[neighbor]:
        dist[neighbor] = alt
        prev[neighbor] = current

  ret = []
  current = "STOP"
  ret.append(current)
  while prev[current] != "NULL":
    current = prev[current]
    ret.append(current)

  return list(reversed(ret))

def test():
  assertEqual.passed = 0
  assertEqual.total = 0

  assertEqual(dot((0.0, 0.0, 0.0), (0.0, 0.0, 0.0)), 0.0, "dot(0,0) == 0")
  assertEqual(add((1.0, 2.0, 3.0), (3.0, 5.0, 7.0)), (4.0, 7.0, 10.0), "vector sum")
  assertEqual(sub((1.0, 1.0, 1.0), (1.0, 1.0, 1.0)), (0.0, 0.0, 0.0), "sub")
  assertEqual(length((1.0, 0.0, 0.0)), 1.0, "length")
  assertEqual(length((0.0, 2.0, 0.0)), 2.0, "length")
  assertVisible((6.500, 0.0, 0.0), (6.500, 1.0, 0.0), True)
  assertVisible((6.400, 0.0, 0.0), (-6.400, 0.0, 0.0), False)
  assertVisible((0.0, 6.400, 0.0), (0.0, -6.400, 0.0), False)
  assertVisible((0.0, 0.0, 6.400), (0.0, 0.0,-6.400), False)
  assertVisible((-6.400, 0.0, 0.0), (6.400, 0.0, 0.0), False)
  assertVisible((0.0, -6.400, 0.0), (0.0, 6.400, 0.0), False)
  assertVisible((0.0, 0.0, -6.400), (0.0, 0.0, 6.400), False)

  f = open('../data/data.txt', 'r')
  lines = f.readlines()

  satellites = read_satellites(lines)
  for key, val in satellites.items():
    if key != "START" and key != "STOP":
      assertEqual(length(val) >= 6.371 + 0.3, True, "Satellite at correct height")
      assertEqual(length(val) <= 6.371 + 0.7, True, "Satellite at correct height")
    else:
      assertEqual(length(val) >= 6.371, True, "Start/stop pos at correct height")
      assertEqual(length(val) <= 6.372, True, "Start/stop pos at correct height")

  graph = create_graph(satellites)
  for key, val in graph.items():
    for v in val:
      assertEqual(key in graph[v], True, "Graph is symmetric (%s %s)" % (key, v))

  shortest_path = find_shortest_path(graph, satellites)

  current = shortest_path[0]
  for v in shortest_path[1:]:
    assertVisible(satellites[v], satellites[current], True)
    current = v;

  print(str(assertEqual.passed) + "/" + str(assertEqual.total) + " tests passed.")

def assertEqual(v1, v2, message):
  print("Testing '"+ message +"'. Status: " + ("pass" if v1 == v2 else "******fail*******"))
  if (v1 == v2): assertEqual.passed +=  1
  assertEqual.total += 1

def assertVisible(v1, v2, is_visible):
  assertEqual(visible(v1, v2), is_visible, "is visible")

if __name__=='__main__':
  if "test" in sys.argv:
    test()
  else:
    main()
