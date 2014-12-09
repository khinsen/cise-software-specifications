from math import sqrt

class Vector(object):

    def __init__(self, x, y, z):
        self.x = x
        self.y = y
        self.z = z

class Cell(object):

    def check_configuration(self, configuration):
        raise NotImplementedError

    def distance(self, p1, p2):
        raise NotImplementedError

class InfiniteCell(object):

    def check_configuration(self, configuration):
        pass

    def distance(self, p1, p2):
        return sqrt((p2.x-p1.x)**2 + (p2.y-p1.y)**2 + (p2.z-p1.z)**2)

class OrthorhombicCell(object):

    def __init__(self, edges):
        self.edges = edges

    def check_configuration(self, configuration):
        for p in configuration:
            assert p.x >= 0. and p.x < self.edges[0]
            assert p.y >= 0. and p.y < self.edges[1]
            assert p.z >= 0. and p.z < self.edges[2]

    def distance(self, p1, p2):
        lx, ly, lz = self.edges
        return sqrt(self.minimum_image(p2.x-p1.x, lx)**2
                    + self.minimum_image(p2.y-p1.y, ly)**2
                    + self.minimum_image(p2.z-p1.z, lz)**2)

    def minimum_image(self, d, l):
        if d > 0.5*l:
            d -= l
        elif d < -0.5*l:
            d += l
        return d

def pair_energy(r):
    LJEnergy = 1.   # kJ/mol
    LJRadius = 0.34 # nm
    LJCutoff = 1.5  # nm
    sr6 = LJRadius**6 / r**6 if (r < LJCutoff) else 0
    return 4 * LJEnergy * (sr6 * sr6 - sr6)

def potential_energy(cell, configuration):
    cell.check_configuration(configuration)
    n = len(configuration)
    e = 0.
    for i in range(n):
        for j in range(i+1, n):
            r = cell.distance(configuration[i],
                              configuration[j])
            pe = pair_energy(r)
            e += pe
    return e

# Two test cases: a triangle and a cubic lattice

def triangle(h):
    cell = InfiniteCell()
    configuration = [Vector(0, 0, 0),
                     Vector(h, 0, 0),
                     Vector(0.5*h, 0.5*sqrt(3)*h, 0)]
    return potential_energy(cell, configuration)

def cubic_lattice(n, h):
    cell = OrthorhombicCell(np.array([n*h, n*h, n*h]))
    configuration = [Vector(h*x, h*y, h*z)
                     for x in range(n)
                     for y in range(n)
                     for z in range(n)]
    return potential_energy(cell, configuration)


print triangle(0.3)
print cubic_lattice(2, 0.375)
print cubic_lattice(4, 0.375)
