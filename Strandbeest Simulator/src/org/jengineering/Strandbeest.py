#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright 2016 Dirk Toewe
#
# This file is part of Game of Pyth.
#
# Game of Pyth is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# Game of Pyth is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with Game of Pyth. If not, see <http://www.gnu.org/licenses/>.
'''
Created on Jul 8, 2016

@author: Dirk Toewe
'''
from PyQt4.Qt import QColor, QComboBox
from PyQt4.QtCore import QTimer, Qt, QRectF, QPointF, QLineF
from PyQt4.QtGui import QWidget, QGridLayout, QLabel, QPushButton, QSlider, \
  QGraphicsView, QGraphicsScene, QPen, QGraphicsItem, QGraphicsEllipseItem,\
  QMenu, QPainterPath, QPixmap, QGraphicsItemGroup, QGraphicsLineItem,\
  QGraphicsPathItem
from collections import defaultdict
from itertools import chain
from sortedcontainers import SortedDict
from numpy import zeros, fromiter, linalg
from numpy.linalg import norm
from timeit import itertools

from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg       as FigureCanvas
from matplotlib.backends.backend_qt4agg import NavigationToolbar2QTAgg as NavigationToolbar
from matplotlib import pyplot

class _Record(object):
  '''
  Used to keep track of the simulation results. Allows changeListeners to be notified when
  new simulation results arrive.
  '''
  def __init__(self):
    self.states = []
    self.velocities = []
    self.times = []
    self.changeListeners = []

class Widget(QWidget):
  '''
  The widget is used to simulate and view the Strandbeest's movement. The Strandbeest.Widget
  assumes that the model is not changed once it's submitted to the __init__() method. It
  is also assumed that the model's velocities are time-invariant.
  '''

  def __init__(self,model):
    super(QWidget,self).__init__()
    if not isinstance(model,Model):
      raise Exception('Constructor argument of type Model expected.')
    self.setWindowTitle('Strandbeest simulator')

    # lazily simulates on full movement cycle of the strandbeest
    def iterateStates():
      '''
      Iterates over the first full movement cycle of the Strandbeest.
      In order to determine when the movement repeats, a heuristic is
      used, that stops simulating once the state vector becomes similar
      enough to the initial state vector. This heuristic assumes that
      the movement of the Strandbeest is time independent.
      '''
      x0 = model.state()
      yield x0
      for _ in range(2):
        model.increment(0.01)
        xi = model.state()
        yield xi
      limit = norm(xi-x0)
      while True:
        model.increment(0.01)
        xi = model.state()
        yield xi
        if norm(xi-x0) < limit: # <- heuristic condition for a full simulation cycle
          break

    # record simulation results
    record = _Record()
    view = _View(model,record)

    def recordStates(states):
      '''
      The _Canvas uses the states of the simulation to display movement curves of
      the nodes. For that purposes the simulated states are stored to the model in
      the record.states list.
      In order to notify the _View of newly simulated states, all callables in
      the record.changeListeners list are notified, with the new state vector as
      argument.
      '''
      for x in states:
        record.states.append(x)
        record.velocities.append( model.v(x,t=None) )
        record.times.append( model.t )
        for listener in record.changeListeners:
          listener(x)
        yield x
      del record.changeListeners

    def recordSceenshots(statesLoop):
      '''
      A little utility i created to be able to generate GIFs from
      the simulation. Saves screenshots from the current view of
      every 10th simulation step.
      '''
      i,j = 0,0
      for x in statesLoop:
        if i < 2*943 and 0 == i % 10:
          img = QPixmap.grabWidget(self)
          img = img.scaled(640,500,transformMode = Qt.SmoothTransformation)
          img.save('data/img%(j)d.jpg' % locals(),'jpg')
          j+=1
        i += 1
        yield x

    statesLoop = iterateStates()
    statesLoop = recordStates(statesLoop)
    statesLoop = itertools.cycle(statesLoop)
#     statesLoop = recordSceenshots(statesLoop)

    def tick():
      '''
      Advances the view by one tick (not neccessarily the simulation)
      '''
      model.setState( next(statesLoop) ) # <- in the first cycle this statement is redundant
    tick()
    timer = QTimer()
    timer.timeout.connect(tick)

    speedLabel = QLabel('? * Realtime')
    speedLabel.resize(speedLabel.sizeHint())
    speedLabel.setToolTip('Simulation speed [ticks/sec.].')
    
    speedSlider = QSlider(Qt.Horizontal)
    speedSlider.resize(speedSlider.sizeHint())
    speedSlider.setMinimum(0)
    speedSlider.setMaximum(6)
    speedSlider.setTickPosition(QSlider.TicksBelow)
    def tpsChange(value): # <- called whenever ticks/sec. change
      '''
      Changes the number of simulation ticks per second.
      '''
      value = 2**(value-2)
      speedLabel.setText( '%(value)6.2f * Realtime' % locals() )
      timer.setInterval(10/value)
    speedSlider.valueChanged[int].connect(tpsChange)
    speedSlider.setValue(4)

    def startStopClick():
      startStopClick.running ^= True 
      if startStopClick.running:
        timer.start()
        startStop.setText('Stop')
      else:
        timer.stop()
        startStop.setText('Start')
    startStop = QPushButton('Start/Stop')
    startStop.resize(startStop.sizeHint())
    startStop.setToolTip('Start the simulation.')
    startStop.clicked.connect(startStopClick)
    startStopClick.running = True
    startStopClick()

    grid = QGridLayout(self)
    grid.addWidget(startStop,0,0)
    grid.addWidget(speedSlider,0,1)
    grid.addWidget(speedLabel,0,2)
    grid.addWidget(view,1,0,1,3)
    

class _View(QGraphicsView):
  '''
  Displays the Strandbeest's nodes and edges. Allows to interactively display
  further data such as movement curves of each node via right-click.
  '''
  def __init__(self,model,record):
    def test(rect):
      x,y,w,h = rect.x(), rect.y(), rect.width(), rect.height()
      self.setSceneRect( QRectF(x-w/2,y-h/2,2*w,2*h) )
    _scene = QGraphicsScene()
    _scene.sceneRectChanged.connect(test)
    # organize the graphics items in layers so that the curves are
    # always in background and the nodes are always in foreground
    groupCurves = QGraphicsItemGroup()
    groupEdges  = QGraphicsItemGroup()
    groupNodes  = QGraphicsItemGroup()
    for group in groupCurves,groupEdges,groupNodes:
      group.setHandlesChildEvents(False)
      _scene.addItem(group)
    super(_View,self).__init__(_scene)
    
    self.setDragMode(QGraphicsView.ScrollHandDrag)

    nodeToScene = {}
    nodeToLineP1 = {} # <- keeps track of lines an their p1, in order to change it during an onNodeMove event
    nodeToLineP2 = {} # <- keeps track of lines an their p2, in order to change it during an onNodeMove event
    class NodeItem(QGraphicsEllipseItem):
      '''
      This nested class handles the visual representation of a node in the _View.
      The NodeItem allows to interactively display information about the node.
      If you move over the node, the name of the node is displayed. On a rightclick
      you can open a context menu which allows you to show/hide a movement curve.
      You can also open a separate window with a graph of the nodes data over time.
      '''
      def __init__(self,node):
        super(NodeItem,self).__init__(0,0,0,0,parent = groupNodes)
        self.setPen( QPen( QColor(0,0,0), 0 ) )
        self.setBrush( QColor(0,0,0) )
        self.setCursor(Qt.ArrowCursor)
        self.__node = node

      def _addPoint(self,state):
        iNode = 2 * model._nodes.index(self.__node)
        x,y = state[iNode:iNode+2]
        assert iNode == 2*model._nodes.index(self.__node)
        item = self._path
        path = item.path()
        path.lineTo(x,-y)
        item.setPath(path)

      def _addPath(self):
        iNode = 2 * model._nodes.index(self.__node)
        x0,y0 = record.states[0][iNode:iNode+2]
        path = QPainterPath()
        path.moveTo(x0,-y0)
        item = QGraphicsPathItem()
        item.setPen( QPen( QColor(0,0,255), 0 ) )
        item.setPath(path)
        self._path = item
        for state in record.states[1:]:
          self._addPoint(state)
        if hasattr(record,'changeListeners'):
          record.changeListeners.append(self._addPoint)
        groupCurves.addToGroup(item)

      def _delPath(self):
        path = self._path
        del self._path
        groupCurves.removeItem(path)
        if hasattr(record,'changeListeners'):
          record.changeListeners.remove(self._addPoint)

      def _addChart(self):
        fig = pyplot.figure()
        canvas = FigureCanvas(fig)

        iNode = 2 * model._nodes.index(self.__node)
        data = {
          't': list(record.times),
          'pos_x': [ x[iNode+0] for x in record.states],
          'pos_y': [ x[iNode+1] for x in record.states],
          'v_x': [ v[iNode+0] for v in record.velocities],
          'v_y': [ v[iNode+1] for v in record.velocities]
        }

        xBox,yBox = QComboBox(), QComboBox()
        keys = list( data.keys() )
        xBox.addItems(keys)
        yBox.addItems(keys)
        xBox.setCurrentIndex( keys.index('pos_y') )
        yBox.setCurrentIndex( keys.index('v_x') )

        plot = fig.add_subplot(1,1,1)
        chart, = plot.plot([],[])

        def replot(*args):
          xSelection = str( xBox.currentText() )
          ySelection = str( yBox.currentText() )
          chart.set_xdata( data[xSelection] )
          chart.set_ydata( data[ySelection] )
          plot.relim()
          plot.autoscale_view(True,True,True)
          canvas.draw()
        replot()

        def addState(x):
          assert iNode == 2 * model._nodes.index(self.__node)
          data['t'].append(record.times[-1])
          pos_x,pos_y = record.states[-1][iNode:iNode+2]
          data['pos_x'].append(pos_x)
          data['pos_y'].append(pos_y)
          v_x, v_y = record.velocities[-1][iNode:iNode+2]
          data['v_x'].append(v_x)
          data['v_y'].append(v_y) 
          replot()

        xBox.currentIndexChanged.connect(replot)
        yBox.currentIndexChanged.connect(replot)

        if hasattr(record,'changeListeners'):
          record.changeListeners.append(addState)

        class WidgetWithCloseEvent(QWidget):
          def __init__(this):
            super(WidgetWithCloseEvent,this).__init__()
          def closeEvent(this,event):
            super(WidgetWithCloseEvent,this).closeEvent(event)
            if event.isAccepted():
              del self._chart
              if hasattr(record,'changeListeners'):
                record.changeListeners.remove(addState)

        wdgt = WidgetWithCloseEvent()
        wdgt.setWindowTitle( "Line chart for node '"+self.__node+"'")

        toolbar = NavigationToolbar(canvas,wdgt)

        layout = QGridLayout()
        layout.addWidget( toolbar, 0,0, 1,4 )
        layout.addWidget(xBox,1,1)
        layout.addWidget(yBox,1,3)
        layout.addWidget( QLabel('x-Axis:'), 1,0, 1,1, Qt.AlignRight )
        layout.addWidget( QLabel('y-Axis:'), 1,2, 1,1, Qt.AlignRight )
        layout.addWidget( canvas, 2,0, 1,4 )

        wdgt.setLayout(layout)
        wdgt.show()
        self._chart = wdgt

      def contextMenuEvent(self,event):
        menu = QMenu()
        if hasattr(self,'_path'):
          hidePath = menu.addAction('Hide Path')
          hidePath.perform = self._delPath
        else:
          showPath = menu.addAction('Show Path')
          showPath.perform = self._addPath
        if not hasattr(self,'_chart'):
          showChart = menu.addAction('Show Chart')
          showChart.perform = self._addChart
        action = menu.exec_( event.screenPos() )
        if None != action:
          action.perform()

    def onNodeMove(name,x,y): # <- FIXME moving will currently not work as intended as the curve would not move and the record would become invalid? :(
      '''
      Whenever a node in model is moved, this method is to be notified
      '''
      r = 2
      nodeToScene[name].setRect(x-r,-y-r,2*r,2*r)
      for line in nodeToLineP1[name]:
        l = line.line()
        l.setP1( QPointF(x,-y) )
        line.setLine(l)
      for line in nodeToLineP2[name]:
        l = line.line()
        l.setP2( QPointF(x,-y) )
        line.setLine(l)

    def onNodeAdd(name,x,y):
      '''
      Whenever a new node is added to model, this method is to be notified.
      '''
      assert not name in nodeToScene
      item = NodeItem(name)
      item.setToolTip( "Node '%(name)s'" % locals() )
      item.setFlags( item.flags() | QGraphicsItem.ItemSendsScenePositionChanges | QGraphicsItem.ItemSendsGeometryChanges )
      nodeToScene[name] = item
      nodeToLineP1[name] = set()
      nodeToLineP2[name] = set()
      onNodeMove(name,x,y)
      
    def onEdgeAdd(start,end):
      '''
      Whenever a new edge is added to model, this method is to be notified.
      '''
      p1 = nodeToScene[start].rect().center()
      p2 = nodeToScene[end]  .rect().center()
      line = QGraphicsLineItem( QLineF(p1,p2) )
      line.setPen( QPen( QColor(32,32,32), 1 ) )
      groupEdges.addToGroup(line)
      nodeToLineP1[start].add(line)
      nodeToLineP2  [end].add(line)

    for name,(x,y) in model.nodes():
      onNodeAdd(name,x,y)
    for start,end in model.edges():
      onEdgeAdd(start,end)
    model.onNodeMoveListeners.add(onNodeMove)  
    model.onNodeAddListeners.add(onNodeAdd)
    model.onEdgeAddListeners.add(onEdgeAdd)

  def wheelEvent(self,event):
    if event.delta() > 0: scale = 1.1
    else                : scale = 0.9
    mx,my = event.x(), event.y()
    mouseInScene = self.mapToScene(mx,my)
    self.scale(scale,scale)
    mouseInScene -= self.mapToScene(mx,my)
    viewport = self.viewport() # <- FIXME the viewport's x and y offset in the QGraphicsView parent may have to be considered as well 
    centerInScene = self.mapToScene(
      viewport.width ()/2,
      viewport.height()/2
    )
    self.centerOn(centerInScene+mouseInScene ) # <- the point in the scene, the mouse is pointing at should not change while zooming

class Model(object):
  '''
  The model of a Stranbeest. The Model consists of a set of nodes, edges and boundary
  conditions. Each node has a unique name and a x and y position which may change
  whenever the simuation is incremented. Each node introduces two degrees of freedom.
  The edges are specified by the nodes they are connecting. The edges are the push/pull
  rods which connect the edges whith one another. An edges keeps the distances between
  two nodes constant and therefore constrains exactly one degree of freedom in the system.
  '''

  def __init__(self):
    '''
    Constructor
    '''
    self._nodes = SortedDict()
    self._edges = defaultdict(set)

  def addNode(self,name,x,y):
    if not isinstance(name,str  ): raise Exception("The 1st argument must be the node's name as str.")
    if not isinstance(x   ,float): raise Exception("The 2nd argument must be the node's x position as float.")
    if not isinstance(y   ,float): raise Exception("The 2nd argument must be the node's y position as float.")
    if name in self._nodes: raise Exception( 'There already exists a node by the name of "%(name)s"' % locals() )
    self._nodes[name] = x,y
    self.__t = 0.0
    for listener in self.onNodeAddListeners:
      listener(name,x,y)

  def addEdge(self,node1,node2):
    if node1 == node2:
      raise Exception('"node1" cannot be equal to "node2".')
    self._edges[node1].add(node2)
    self._edges[node2].add(node1)
    for listener in self.onEdgeAddListeners:
      listener( min(node1,node2), max(node1,node2) )

  def pos(self,name):
    return self._nodes[name]

  def move(self,name,x,y):
    self._nodes[name] = x,y
    for listener in self.onNodeMoveListeners:
      listener(name,x,y)

  def state(self):
    return fromiter( chain.from_iterable( self._nodes.values() ), float )

  def setState(self,state):
    for i,(x,y) in enumerate( zip(state[::2],state[1::2]) ):
      self.move(self._nodes.keys()[i],x,y)

  @property
  def t(self):
    return self.__t

  def increment(self,dt):
    v = self.v
    t0 = self.__t
    x0 = self.state()
    # https://en.wikipedia.org/wiki/Runge%E2%80%93Kutta_methods#The_Runge.E2.80.93Kutta_method
    k0 = v(x0,           t0)
    k1 = v(x0+k0*(dt/2), t0+dt/2)
    k2 = v(x0+k1*(dt/2), t0+dt/2)
    k3 = v(x0+k2*(dt),   t0+dt)
    self.setState( x0 + dt/6 * (k0+k1+k2+k3) )
    self.__t += dt

  def v(self,x,t):
    lhs = zeros( 2*[len(x)] )
    rhs = zeros( len(x) )
    iRows = iter( range( len(x) ) )
    for start,end in self.edges():
      iStart = 2*self._nodes.index(start)
      iEnd   = 2*self._nodes.index(end)
      iRow = next(iRows)
      dx = x[iEnd+0] - x[iStart+0] 
      dy = x[iEnd+1] - x[iStart+1]
      lhs[iRow,iStart+0] = dx; lhs[iRow,iEnd+0] = -dx
      lhs[iRow,iStart+1] = dy; lhs[iRow,iEnd+1] = -dy
      rhs[iRow] = 0
    for bc in self.bcs:
      bc.addEquations(x,t,iRows,lhs,rhs)
    return linalg.solve(lhs,rhs)

  def nodes(self):
    return self._nodes.iteritems()

  def edges(self):
    for node1,neighbors in self._edges.items():
      for node2 in neighbors:
        if node1 < node2:
          yield node1,node2

  bcs = []

  onEdgeAddListeners = set() # <- FIXME should be a multiset
  onNodeAddListeners = set() # <- FIXME should be a multiset
  onNodeMoveListeners = set() # <- FIXME should be a multiset

class RotationBC(object):
  '''
  A boundary constraint that restricts both degrees of freedom of a node in a way that
  the node moves in a circle around the given `pivot` with an angular speed of `speed`.
  '''
  def __init__(self,model,node,pivot,speed):
    super(RotationBC,self).__init__()
    self.__model = model
    self.__node = node
    self.__pivot = pivot
    self.__speed = speed

  def addEquations(self,x,t,iRows,lhs,rhs):
    '''
    Called by a Model several times during increment() in order to retrieve the
    equations for this boundary condition. The equations are written into the
    equation system represented by the left hand matrix lhs and the right hand
    vector rhs. One or more unused row indices in the equation system can be
    retrieved from the iRows vector.

    The RotationBC adds exactly two equations as it constrains exactly two degrees
    of freedom of a node.
    '''
    iNode = 2*self.__model._nodes.index(self.__node)
    d = x[iNode:iNode+2] - self.__pivot
    iRow0 = next(iRows)
    iRow1 = next(iRows)
    lhs[iRow0][iNode+0] = 1 # +d[1]*self.__speed
    lhs[iRow1][iNode+1] = 1 # -d[0]*self.__speed
    rhs[iRow0] = -d[1]*self.__speed
    rhs[iRow1] = +d[0]*self.__speed

class NodeMoveBC(object):
  '''
  Enforces that a `node` moves with the constant velocity `vel` along in the specified
  `direction`.
  '''
  def __init__(self,model,node,direction,vel):
    super(NodeMoveBC,self).__init__()
    self.__model = model
    self.__node = node
    self.__direction = direction
    self.__vel = vel

  def addEquations(self,x,t,iRows,lhs,rhs):
    '''
    Called by a Model several times during increment() in order to retrieve the
    equations for this boundary condition. The equations are written into the
    equation system represented by the left hand matrix lhs and the right hand
    vector rhs. One or more unused row indices in the equation system can be
    retrieved from the iRows vector.

    The NodeMoveBC adds exactly one equation as it constrains exactly one degree
    of freedom of a node.
    '''
    iNode = 2 * self.__model._nodes.index(self.__node)
    iRow = next(iRows)
    direction = self.__direction
    lhs[iRow,iNode+0] = direction[0]
    lhs[iRow,iNode+1] = direction[1]
    rhs[iRow] = self.__vel