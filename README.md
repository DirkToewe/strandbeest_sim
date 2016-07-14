Python Strandbeest Simulation
=============================
My second Python Project is an application that allows the 2D Kinematics simulation of simple rod and ball
joint frameworks. As an example the movement of a [Strandbeest](https://en.wikipedia.org/wiki/Strandbeest)
leg is simulated. The app is an ideal introduction to Python in mechanical engineering as it covers two key
aspects: Visualization and Numerics.

<img src="https://github.com/DirkToewe/strandbeest_sim/blob/master/strandbeest_sim.gif" width="60%">

Model
=====
A ```Model``` consists of:
  * a set of nodes. Each node has a unique name and a current position (x,y)
  * a set edges (or rods), specified by the name of the two edges it connects
  * a set of boundary constraints that constrain the movement of the nodes

Each node introduces two degrees of freedom to the model (x and y). A rod enforces a constant distance
between two edges. Each rod constrains one degree of freedom. The boundary constraints constrain the
movement of a node. The NodeMoveBC constrains the speed of a node alongside a direction vector. The
velocity orthogonal to the direction is not constrained. Each NodeMoveBC there constrains one degree
of freedom. The RotationBC constrains a node to move in a circular motion around a pivot point.
It constraings two degrees of freedom. If there is as many (rendundancy free) constraints as
there are degrees of freedom, the system is terminate and can be simulated.

Nodes are added via the ```addNode()``` method, edges via ```addEdge()```. Boundary constraints
are added directly to the ```bcs``` list.

```python
import Strandbeest

leg = Strandbeest.Model()
leg.addNode('p0',   0.,             0.)
leg.addNode('p1', -35.0353991829,  19.5071987762)
leg.addNode('p2',  17.7855809088,  37.4956412365)
leg.addNode('p3',  50.7538954493,  -0.0954512769988)
leg.addNode('p4',  2.77183874468, -39.2021288959)
leg.addNode('p5', -17.0054756067, -84.0335669394)
leg.addNode('p6', -28.0419102143, -19.2671627536)
leg.addEdge('p0','p1')
leg.addEdge('p0','p2')
leg.addEdge('p0','p4')
leg.addEdge('p1','p2')
leg.addEdge('p1','p6')
leg.addEdge('p2','p3')
leg.addEdge('p3','p4')
leg.addEdge('p4','p5')
leg.addEdge('p4','p6')
leg.addEdge('p5','p6')
leg.bcs += [
  Strandbeest.NodeMoveBC(leg,'p0',(1,0),0),
  Strandbeest.NodeMoveBC(leg,'p0',(0,1),0),
  Strandbeest.RotationBC(leg,'p3',(38,7.8),-1)
]
```

The ```Model::increment()``` method advances the simulation by the specified time step ```dt```. The node
positions change accordingly. The [Runge Kutta Method (RK4)](https://en.wikipedia.org/wiki/Runge%E2%80%93Kutta_methods)
is used to numerically integrate the node movements. Generally speaking, the simulation precision increases
with smaller time steps.

Widget
======
The ```Strandbeest.Widget``` allows the visualization and inspection of the ```Model```'s simulation results.
Once a ```Model``` submitted to the ```Widget``` it must not be modified. The ```Widget``` is handled like
a regular [PyQt4](https://wiki.python.org/moin/PyQt) [QWidget](http://pyqt.sourceforge.net/Docs/PyQt4/qwidget.html),
with the difference of taking a ```Model``` as constructor argument:

```python
from PyQt4.QtGui import QApplication
import sys
import Strandbeest

app = QApplication(sys.argv)
wdgt = Strandbeest.Widget(model)
wdgt.setGeometry(200,200,800,600)
wdgt.show()
sys.exit( app.exec_() )
```

The Start/Stop button on the upper left hand of the window allows to start and stop the animation. The slider
changes the animation speed. Moving the mouse over the nodes displays the name of the nodes. A secondary mouse
button click on a node opens a context menu that allows to display a movement curve for the node and allows
to display the nodes values in a separate diagram.
