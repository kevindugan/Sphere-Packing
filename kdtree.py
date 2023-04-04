

class KDTree:
    def __init__(self, spheres):
        self.k = len(spheres[0].c)
        self.root = self.build(spheres)

    def __str__(self):
        return str(self.root)

    def build(self, spheres, depth=0):
        if len(spheres) == 0:
            return None
        
        sorted_spheres = sorted(spheres, key=lambda x: x.c[depth % self.k])
        middle = len(sorted_spheres) // 2

        return KDNode(value=sorted_spheres[middle],
                      left=self.build(spheres=sorted_spheres[:middle], depth=depth+1),
                      right=self.build(spheres=sorted_spheres[middle+1:], depth=depth+1)
               )

class KDNode:
    def __init__(self, value=None, left=None, right=None):
        self.value = value
        self.left = left
        self.right = right

    def __str__(self):
        return f"[{self.value}, {self.left}, {self.right}]"
    
from numpy import zeros
from matplotlib.patches import Circle, CirclePolygon
class Sphere:
    def __init__(self, center=zeros(2), radius=1, color='b'):
        assert center.shape == (2,)
        self.center = center
        self.radius = radius
        self.color = color

    def __str__(self) -> str:
        return f"({self.x:.3f}, {self.y:.3f}), {self.r:.3f}"

    @property
    def r(self):
        return self.radius

    @property
    def x(self):
        return self.center[0]

    @property
    def y(self):
        return self.center[1]
    
    @property
    def c(self):
        return self.center
    
    def draw(self):
        # return Circle(self.center, self.radius, color=self.color)
        return CirclePolygon(self.center, self.radius, color=self.color, resolution=10)
    
    def setColor(self, color):
        self.color = color