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

from PyQt4.QtGui import QApplication
import sys

import Strandbeest


def createModel():
  model = Strandbeest.Model()
  # LEFT LEG
  model.addNode('p0',   0.,             0.)
  model.addNode('p1', -35.0353991829,  19.5071987762)
  model.addNode('p2',  17.7855809088,  37.4956412365)
  model.addNode('p3',  50.7538954493,  -0.0954512769988)
  model.addNode('p4',  2.77183874468, -39.2021288959)
  model.addNode('p5', -17.0054756067, -84.0335669394)
  model.addNode('p6', -28.0419102143, -19.2671627536)
  model.addEdge('p0','p1')
  model.addEdge('p0','p2')
  model.addEdge('p0','p4')
  model.addEdge('p1','p2')
  model.addEdge('p1','p6')
  model.addEdge('p2','p3')
  model.addEdge('p3','p4')
  model.addEdge('p4','p5')
  model.addEdge('p4','p6')
  model.addEdge('p5','p6')
  model.bcs += [
    Strandbeest.NodeMoveBC(model,'p0',(1,0),0),
    Strandbeest.NodeMoveBC(model,'p0',(0,1),0),
    Strandbeest.RotationBC(model,'p3',(38,7.8),-1)
  ]
  # RIGHT LEG
  for name,(x,y) in list( model.nodes() ):
    model.addNode(name+"'",2*38-x,y)
  for start,end in list( model.edges() ):
    model.addEdge(start+"'",end+"'")
  model.bcs += [
    Strandbeest.NodeMoveBC(model,"p0'",(1,0),0),
    Strandbeest.NodeMoveBC(model,"p0'",(0,1),0),
    Strandbeest.RotationBC(model,"p3'",(38,7.8),-1)
  ]
  # DISPLAY PIVOT AS FIXED NODE
  model.addNode('pivot', 38., 7.8)
  model.bcs += [
    Strandbeest.NodeMoveBC(model,'pivot',(1,0),0),
    Strandbeest.NodeMoveBC(model,'pivot',(0,1),0)
  ]
  return model

if __name__ == '__main__':
  app = QApplication(sys.argv)

  wdgt = Strandbeest.Widget( createModel() )
  wdgt.setGeometry(200,200,800,600)
  wdgt.show()

  sys.exit( app.exec_() )