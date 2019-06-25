# -*- coding: utf-8 -*-
"""
/***************************************************************************
 PROPAGACAO_ERRO
                                 A QGIS plugin
 PROPAGACAO_ERRO
                              -------------------
        begin                : 2018-08-06
        git sha              : $Format:%H$
        copyright            : (C) 2018 by Livia Faria Sampaio
        email                : livias_93@hotmail.com
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
"""

import numpy as np
from math import *
from PyQt4.QtCore import *
from PyQt4.QtGui import *
import processing
import os
import sys
from buddy import *
from qgis.core import *
import psycopg2
import os.path
import osgeo.ogr
import random
from sympy import Symbol

# Initialize Qt resources from file resources.py
import resources

# Import the code for the dialog
from PROPAGACAO_ERRO_dialog import PROPAGACAO_ERRODialog
import os.path


class PROPAGACAO_ERRO:
    """QGIS Plugin Implementation."""

    def __init__(self, iface):
        """Constructor.

        :param iface: An interface instance that will be passed to this class
            which provides the hook by which you can manipulate the QGIS
            application at run time.
        :type iface: QgisInterface
        """
        # Save reference to the QGIS interface
        self.iface = iface
        # initialize plugin directory
        self.plugin_dir = os.path.dirname(__file__)
        # initialize locale
        locale = QSettings().value('locale/userLocale')[0:2]
        locale_path = os.path.join(
            self.plugin_dir,
            'i18n',
            'PROPAGACAO_ERRO_{}.qm'.format(locale))

        if os.path.exists(locale_path):
            self.translator = QTranslator()
            self.translator.load(locale_path)

            if qVersion() > '4.3.3':
                QCoreApplication.installTranslator(self.translator)


        # Declare instance attributes
        self.actions = []
        self.menu = self.tr(u'&PROPAGACAO_ERRO')
        # TODO: We are going to let the user set this up in a future iteration
        self.toolbar = self.iface.addToolBar(u'PROPAGACAO_ERRO')
        self.toolbar.setObjectName(u'PROPAGACAO_ERRO')

    # noinspection PyMethodMayBeStatic
    def tr(self, message):
        """Get the translation for a string using Qt translation API.

        We implement this ourselves since we do not inherit QObject.

        :param message: String for translation.
        :type message: str, QString

        :returns: Translated version of message.
        :rtype: QString
        """
        # noinspection PyTypeChecker,PyArgumentList,PyCallByClass
        return QCoreApplication.translate('PROPAGACAO_ERRO', message)


    def add_action(
        self,
        icon_path,
        text,
        callback,
        enabled_flag=True,
        add_to_menu=True,
        add_to_toolbar=True,
        status_tip=None,
        whats_this=None,
        parent=None):
        """Add a toolbar icon to the toolbar.

        :param icon_path: Path to the icon for this action. Can be a resource
            path (e.g. ':/plugins/foo/bar.png') or a normal file system path.
        :type icon_path: str

        :param text: Text that should be shown in menu items for this action.
        :type text: str

        :param callback: Function to be called when the action is triggered.
        :type callback: function

        :param enabled_flag: A flag indicating if the action should be enabled
            by default. Defaults to True.
        :type enabled_flag: bool

        :param add_to_menu: Flag indicating whether the action should also
            be added to the menu. Defaults to True.
        :type add_to_menu: bool

        :param add_to_toolbar: Flag indicating whether the action should also
            be added to the toolbar. Defaults to True.
        :type add_to_toolbar: bool

        :param status_tip: Optional text to show in a popup when mouse pointer
            hovers over the action.
        :type status_tip: str

        :param parent: Parent widget for the new action. Defaults None.
        :type parent: QWidget

        :param whats_this: Optional text to show in the status bar when the
            mouse pointer hovers over the action.

        :returns: The action that was created. Note that the action is also
            added to self.actions list.
        :rtype: QAction
        """

        # Create the dialog (after translation) and keep reference
        self.dlg = PROPAGACAO_ERRODialog()

        icon = QIcon(icon_path)
        action = QAction(icon, text, parent)
        action.triggered.connect(callback)
        action.setEnabled(enabled_flag)

        if status_tip is not None:
            action.setStatusTip(status_tip)

        if whats_this is not None:
            action.setWhatsThis(whats_this)

        if add_to_toolbar:
            self.toolbar.addAction(action)

        if add_to_menu:
            self.iface.addPluginToMenu(
                self.menu,
                action)

        self.actions.append(action)

        return action

    def initGui(self):
        """Create the menu entries and toolbar icons inside the QGIS GUI."""

        icon_path = ':/plugins/PROPAGACAO_ERRO/icon.png'
        self.add_action(
            icon_path,
            text=self.tr(u'PROPAGACAO_ERRO'),
            callback=self.run,
            parent=self.iface.mainWindow())

        
        self.dlg.toolButton.clicked.connect(self.OpenLayer)
        self.loadLayer()  
        self.dlg.comboBox.clear()
		
    def loadLayer(self):
        self.dlg.comboBox.clear()
        layers = self.iface.legendInterface().layers()
        layer_list = []
        for layer in layers:
            if layer.type() == QgsMapLayer.VectorLayer:
                layer_list.append(layer.name())
        self.dlg.comboBox.addItems(layer_list)
            
    def OpenLayer(self):
        infile = str(QFileDialog.getOpenFileName(filter = 'shapefiles file (*.shp)'))
        if infile is not None:
            nome = os.path.splitext(os.path.basename(infile))[0]
            
            self.iface.mapCanvas().freeze()
            self.mylayer=self.iface.addVectorLayer(infile,nome,'ogr')
            self.ProjInstance  = QgsProject.instance()
            self.root=self.ProjInstance.layerTreeRoot()
            QgsMapLayerRegistry.instance().addMapLayer(self.mylayer)
            self.iface.mapCanvas().freeze(False)
            
            self.dlg.toolButton.clicked.connect(self.OpenLayer)
            self.loadLayer()
            nome = nome
        
        

    def unload(self):
        """Removes the plugin menu item and icon from QGIS GUI."""
        for action in self.actions:
            self.iface.removePluginMenu(
                self.tr(u'&PROPAGACAO_ERRO'),
                action)
            self.iface.removeToolBarIcon(action)
        # remove the toolbar
        del self.toolbar


    def run(self):
        """Run method that performs all the real work"""
        # show the dialog
        self.dlg.show()
        # Run the dialog event loop
        result = self.dlg.exec_()

        if result == 1:
            Xa = str(self.dlg.Xa.text())
            Ya = str(self.dlg.Ya.text())
            Za = str(self.dlg.Za.text())
            Xb = str(self.dlg.Xb.text())
            Yb = str(self.dlg.Yb.text())
            Zb = str(self.dlg.Zb.text())
            Xc = str(self.dlg.Xc.text())
            Yc = str(self.dlg.Yc.text())
            Zc = str(self.dlg.Zc.text())
            n = str(self.dlg.nseries.text())
            PrecAng = int(self.dlg.PrecAngular.text())
            TL2 = int(self.dlg.PrecLinear2.text())
            TL1 = int(self.dlg.PrecLinear1.text())
            coord_x = str(self.dlg.coord_X.text())
            coord_y = str(self.dlg.coord_Y.text())
            coord_z = str(self.dlg.coord_Z.text())
            
            x = coord_x
            y = coord_y
            z = coord_z
            
            Xa1 = Xa
            Ya1 = Ya
            Za1 = Za
            
            Xb1 = Xb
            Yb1 = Yb
            Zb1 = Zb
            
            Xc1 = Xc
            Yc1 = Yc
            Zc1 = Zc
            
        ### ENTRAR COM AS COORDENADAS DA ESTACAO
        conn = psycopg2.connect(database="PROPAGACAO_ERRO", host="localhost", user="postgres", password="postgres", port="5432")
        cur = conn.cursor()
        
        cur.execute ("DROP TABLE IF EXISTS estacao")
        cur.execute ("""CREATE TABLE estacao (id SERIAL PRIMARY KEY,nome TEXT,coord_x NUMERIC,coord_y NUMERIC,coord_z NUMERIC)""")
        cur.execute ("""SELECT AddGeometryColumn('public','estacao','geom','4326','POINT', 3)""")

        cur.execute ("""INSERT INTO estacao (nome,coord_x,coord_y,coord_z,geom) VALUES ('1','%s','%s','%s',ST_GeometryFromText('POINT(%s %s %s)','4326'))"""%(x,y,z,coord_x,coord_y,coord_z))

        conn.commit()
        cur.close()
        conn.close()
            
        #ADD VIEW TO QGIS estacao
        uri = QgsDataSourceURI()
        uri.setConnection("localhost", "5432", "PROPAGACAO_ERRO", "postgres", "postgres")
        uri.setDataSource("public", "estacao", "geom", "", "4326")
        teste = QgsVectorLayer(uri.uri(), "estacao", "postgres")
        QgsMapLayerRegistry.instance().addMapLayer(teste)
		
        ### ENTRAR COM AS COORDENADAS DO PONTO A
        conn = psycopg2.connect(database="PROPAGACAO_ERRO", host="localhost", user="postgres", password="postgres", port="5432")
        cur = conn.cursor()
        
        cur.execute ("DROP TABLE IF EXISTS pontoa")
        cur.execute ("""CREATE TABLE pontoa (id SERIAL PRIMARY KEY,nome TEXT,coord_x NUMERIC,coord_y NUMERIC,coord_z NUMERIC)""")
        cur.execute ("""SELECT AddGeometryColumn('public','pontoa','geom','4326','POINT', 3)""")

        cur.execute ("""INSERT INTO pontoa (nome,coord_x,coord_y,coord_z,geom) VALUES ('1','%s','%s','%s',ST_GeometryFromText('POINT(%s %s %s)','4326'))"""%(Xa,Ya,Za,Xa1,Ya1,Za1))

        conn.commit()
        cur.close()
        conn.close()
            
        #ADD VIEW TO QGIS ponto A
        uri = QgsDataSourceURI()
        uri.setConnection("localhost", "5432", "PROPAGACAO_ERRO", "postgres", "postgres")
        uri.setDataSource("public", "pontoa", "geom", "", "4326")
        teste = QgsVectorLayer(uri.uri(), "pontoa", "postgres")
        QgsMapLayerRegistry.instance().addMapLayer(teste)
		
        ### ENTRAR COM AS COORDENADAS DO PONTO B
        conn = psycopg2.connect(database="PROPAGACAO_ERRO", host="localhost", user="postgres", password="postgres", port="5432")
        cur = conn.cursor()
        
        cur.execute ("DROP TABLE IF EXISTS pontob")
        cur.execute ("""CREATE TABLE pontob (id SERIAL PRIMARY KEY,nome TEXT,coord_x NUMERIC,coord_y NUMERIC,coord_z NUMERIC)""")
        cur.execute ("""SELECT AddGeometryColumn('public','pontob','geom','4326','POINT', 3)""")

        cur.execute ("""INSERT INTO pontob (nome,coord_x,coord_y,coord_z,geom) VALUES ('1','%s','%s','%s',ST_GeometryFromText('POINT(%s %s %s)','4326'))"""%(Xb,Yb,Zb,Xb1,Yb1,Zb1))

        conn.commit()
        cur.close()
        conn.close()
            
        #ADD VIEW TO QGIS pontob
        uri = QgsDataSourceURI()
        uri.setConnection("localhost", "5432", "PROPAGACAO_ERRO", "postgres", "postgres")
        uri.setDataSource("public", "pontob", "geom", "", "4326")
        teste = QgsVectorLayer(uri.uri(), "pontob", "postgres")
        QgsMapLayerRegistry.instance().addMapLayer(teste)
		
        
        ### ENTRAR COM AS COORDENADAS DO PONTO C
        conn = psycopg2.connect(database="PROPAGACAO_ERRO", host="localhost", user="postgres", password="postgres", port="5432")
        cur = conn.cursor()
        
        cur.execute ("DROP TABLE IF EXISTS pontoc")
        cur.execute ("""CREATE TABLE pontoc (id SERIAL PRIMARY KEY,nome TEXT,coord_x NUMERIC,coord_y NUMERIC,coord_z NUMERIC)""")
        cur.execute ("""SELECT AddGeometryColumn('public','pontoc','geom','4326','POINT', 3)""")

        cur.execute ("""INSERT INTO pontoc (nome,coord_x,coord_y,coord_z,geom) VALUES ('1','%s','%s','%s',ST_GeometryFromText('POINT(%s %s %s)','4326'))"""%(Xc,Yc,Zc,Xc1,Yc1,Zc1))

        conn.commit()
        cur.close()
        conn.close()
            
        #ADD VIEW TO QGIS pontob
        uri = QgsDataSourceURI()
        uri.setConnection("localhost", "5432", "PROPAGACAO_ERRO", "postgres", "postgres")
        uri.setDataSource("public", "pontoc", "geom", "", "4326")
        teste = QgsVectorLayer(uri.uri(), "pontoc", "postgres")
        QgsMapLayerRegistry.instance().addMapLayer(teste)
        
        
        
        #CRIANDO LISTA VAZIA PARA INSERIR OS DADOS DOS PONTOS
        AX = []
        AY = []
        AZ = []
        TX = []
        TY = []
        TZ = []
        aX = []
        aY = []
        aZ = []
        bX = []
        bY = []
        bZ = []
        cX = []
        cY = []
        cZ = []
        Ponto = []
        
        
        if result:
            layers = QgsMapLayerRegistry.instance().mapLayers().values()
            for layer in layers:
                if layer.name() == 'PONTOS_ENCOSTA':
                    for feature in layer.getFeatures():
                        AX.append(feature['coord_x'])
                        AY.append(feature['coord_y'])
                        AZ.append(feature['coord_z'])
                        Ponto.append(feature['ponto'])
                if layer.name() == 'estacao':
                    for feature in layer.getFeatures():
                        TX.append(feature['coord_x'])
                        TY.append(feature['coord_y'])
                        TZ.append(feature['coord_z'])
                if layer.name() == 'pontoa':
                    for feature in layer.getFeatures():
                        aX.append(feature['coord_x'])
                        aY.append(feature['coord_y'])
                        aZ.append(feature['coord_z'])
                if layer.name() == 'pontob':
                    for feature in layer.getFeatures():
                        bX.append(feature['coord_x'])
                        bY.append(feature['coord_y'])
                        bZ.append(feature['coord_z'])
                if layer.name() == 'pontoc':
                    for feature in layer.getFeatures():
                        cX.append(feature['coord_x'])
                        cY.append(feature['coord_y'])
                        cZ.append(feature['coord_z'])
                        
            print Ponto
            print AX
            print AY
            print AZ
            print TX
            print TY
            print TZ
            print aX
            print aY
            print aZ
            print bX
            print bY
            print bZ
            print cX
            print cY
            print cZ
            print PrecAng 
            print TL2
            print TL1 
            
        ### Funcao para calculo de distancia
        distancia = [0 for x in range(len(AX))]
        def dist(AX1,TX,AY1,TY,AZ1,TZ):
            Distancia = sqrt(((AX1-TX)**2)+((AY1-TY)**2+((AZ1-TZ)**2)))
            return Distancia
        
        # Calcula a distancia entre a estacao e os pontos A, B, C:
        dTA = dist(TX[0],aX[0],TY[0],aY[0],TZ[0],aZ[0])
        dTB = dist(TX[0],bX[0],TY[0],bY[0],TZ[0],bZ[0])
        dTC = dist(TX[0],cX[0],TY[0],cY[0],TZ[0],cZ[0])
        print "as distancias entre os pontos T e A,B,C sao", dTA,dTB,dTC 
        
        
        def azimute(xa, ya, xb, yb):
            if (xa == xb) and (ya == yb):
                print('ERRO: ponto A coincide com ponto B.')
                return None
            else:
                alfa = np.arctan2(xb-xa, yb-ya)
                if alfa < 0:
                    alfa += 2*np.pi
                return alfa

    
        ### Calculo dos Azimutes TA, TB, TC:
        AzTA = azimute(aX[0],aY[0],TX[0],TY[0])
        AzTB = azimute(bX[0],bY[0],TX[0],TY[0])
        AzTC = azimute(cX[0],cY[0],TX[0],TY[0])
        print "os azimutes sao", AzTA,AzTB,AzTC
        


        def calcAng(xa, ya, xb, yb, xc, yc):
            alfaBA = azimute(xb, yb, xa, ya)
            alfaBC = azimute(xb, yb, xc, yc)
            alfa = alfaBA - alfaBC
            if alfa < 0:
                alfa += 2*np.pi
            return alfa

        def calcDrvAng(xa, ya, xb, yb, xc, yc):
            dAB2 = (xa-xb)**2 + (ya-yb)**2
            if (ya == yb):
                drv = [0, 0, 0, 0, 0, 0]
            else:
                drv = [(ya-yb)/dAB2, (xb-xa)/dAB2, (yb-ya)/dAB2, (xa-xb)/dAB2, 0, 0]
            dBC2 = (xb-xc)**2 + (yb-yc)**2
            if (yb != yc):
                drv[2] -= (yb-yc)/dBC2
                drv[3] -= (xc-xb)/dBC2
                drv[4] -= (yc-yb)/dBC2
                drv[5] -= (xb-xc)/dBC2
            return drv
        
        ###Funcao que calcula o angulo Horizontal 
        AngHorizontal = [0 for x in range(len(AX))]
        def angHorizontal(a,b):
            HTA2 = b-a
            if HTA2<0:
                return HTA2+360
            else:
                return HTA2
            
        
        # Calcula os angulos A, B, C:
        BTA = calcAng(bX[0],bY[0],TX[0],TY[0],aX[0],aY[0])
        CTB = calcAng(cX[0],cY[0],TX[0],TY[0],bX[0],bY[0])
        print "os angulos BTA e CTB sao", BTA, CTB
        
        
        ### Funcao para calculo do angulo zenital
        AngZenital = [0 for x in range(len(AX))]
        def ang (AX1,TX,AY1,TY,AZ1,TZ):
            angul = (AZ1-TZ)/(sqrt(((AX1-TX)**2)+((AY1-TY)**2)+((AZ1-TZ)**2)))
            angulo = degrees(acos(angul))
            return angulo

        ### Calculo dos angulos Zenitais:
        AngZAT = ang(aX[0],TX[0],aY[0],TY[0],aZ[0],TZ[0])
        AngZBT = ang(bX[0],TX[0],bY[0],TY[0],bZ[0],TZ[0])
        AngZCT = ang(cX[0],TX[0],cY[0],TY[0],cZ[0],TZ[0])
        print "os angulos zenitais sao", AngZAT,AngZBT,AngZCT
        
        
        ### Calculo da precisão angular:
        Prec_ang = ((PrecAng*sqrt(2))/3600)*pi/180
        print "a precisao angular calculada e ", Prec_ang
        
        ### Funcao para calculo dos desvios-padrao das distancias:
        DesvDist = [0 for x in range(len(AX))]
        def desvDist(Dist):
            SigD = sqrt(TL2**2+(TL1*Dist/1000)**2)/1000    
            return SigD

        ### Calculo dos desvios-padrao das distancias:
        DesvdTA = desvDist(dTA)
        DesvdTB = desvDist(dTB)
        DesvdTC = desvDist(dTC)
        print "o desvio padrao das distancias sao ", DesvdTA,DesvdTB,DesvdTC
        
        ### Calculo dos valores possivelmente medidos:
        
        # Calcula a distancia com erro:
        dTA1 = (np.random.standard_normal(1)*DesvdTA)+dTA
        dTB1 = (np.random.standard_normal(1)*DesvdTB)+dTB
        dTC1 = (np.random.standard_normal(1)*DesvdTC)+dTC
        print "as dists com erro sao ", dTA1,dTB1,dTC1
        
        # Calcula os angulos zenitais com erro:
        AngZAT1 = (np.random.standard_normal(1)*Prec_ang)+AngZAT
        AngZBT1 = (np.random.standard_normal(1)*Prec_ang)+AngZBT
        AngZCT1 = (np.random.standard_normal(1)*Prec_ang)+AngZCT
        BTA1 = (np.random.standard_normal(1)*Prec_ang)+BTA
        CTB1 = (np.random.standard_normal(1)*Prec_ang)+CTB
        print "os angulos com erro sao ", AngZAT1,AngZBT1,AngZCT1,BTA1,CTB1
        
        
        # Calculo das derivadas
        dapX = (TX[0]-aX[0])/dTA1;
        dapY = (TY[0]-aY[0])/dTA1;
        dapZ = (TZ[0]-aZ[0])/dTA1;
        dbpX = (TX[0]-bX[0])/dTB1;
        dbpY = (TY[0]-bY[0])/dTB1;
        dbpZ = (TZ[0]-bZ[0])/dTB1;
        dcpX = (TX[0]-cX[0])/dTC1;
        dcpY = (TY[0]-cY[0])/dTC1;
        dcpZ = (TZ[0]-cZ[0])/dTC1;
        dalfa = calcDrvAng(bX[0],bY[0],TX[0],TY[0],aX[0],aY[0]);
        dbeta = calcDrvAng(cX[0],cY[0],TX[0],TY[0],bX[0],bY[0]);
        dalfaX = dalfa[2];
        dalfaY = dalfa[3];
        dalfaZ = 0;
        dbetaX = dbeta[2];
        dbetaY = dbeta[3];
        dbetaZ = 0;
        dVapX = 0;
        dVapY = 0;
        dVapZ = -1/(dTA1*((1-(TZ[0]/dTA1-aZ[0]/dTA1)**2)**(1/2)))
        dVbpX = 0;
        dVbpY = 0;
        dVbpZ = -1/(dTB1*((1-(TZ[0]/dTB1-bZ[0]/dTB1)**2)**(1/2)))
        dVcpX = 0;
        dVcpY = 0;
        dVcpZ = -1/(dTC1*((1-(TZ[0]/dTC1-cZ[0]/dTC1)**2)**(1/2)))
        print "os valores das derivadas", dalfaX,dalfaY,dalfaZ,dbetaX,dbetaY,dbetaZ,dVapX,dVapY,dVapZ,dVbpX,dVbpY,dVbpZ,dVcpX,dVcpY,dVcpZ 
        
        
        def ajustaParametrico(Lb,P,Xo,Lo,A):
            n = np.size(A,0)
            u = np.size(A,1)
            # Check dimensional
            if (np.size(Lb,0) != n) or (np.size(Lb,1) != 1):
                print('ERRO: Dimensionamento do vetor Lb.')
                return None
            if (np.size(P,0) != n) or (np.size(P,1) != n):
                print('ERRO: Dimensionamento da matriz dos Pesos.')
                return None
            if (np.size(Xo,0) != u) or (np.size(Xo,1) != 1):
                print('ERRO: Dimensionamento do vetor Xo.')
                return None
            if (np.size(Lo,0) != n) or (np.size(Lo,1) != 1):
                print('ERRO: Dimensionamento do vetor Lo.')
                return None
            if (n <= u):
                print('ERRO: Dimensionamento da matriz A.')
                return None
            # Ajustamento
            L = Lo - Lb
            At = A.T
            N = At*P*A
            U = At*P*L
            Ninv = N.I
            X = -Ninv*U
            Xa = Xo + X
            V = A*X + L
            vp = V.T*P*V/(n-u)
            mvcXa = np.multiply(vp, Ninv)
            sigmaXa = np.sqrt(np.diagonal(mvcXa))
            return (Xa,sigmaXa)
        
        # Ajustamento Metodo Parametrico:
        Lb =np.matrix([dTA1,dTB1,dTC1,BTA1,CTB1,AngZAT1,AngZBT1,AngZCT1])
        Lo =(np.matrix([dTA,dTB,dTC,BTA,CTB,AngZAT,AngZBT,AngZCT])).T

        P = np.diag([1/DesvdTA**2, 1/DesvdTB**2, 1/DesvdTC**2, 1/Prec_ang**2, 1/Prec_ang**2, 1/Prec_ang**2, 1/Prec_ang**2, 1/Prec_ang**2])    
        print "a matriz Peso e", P

        A1 = np.matrix([[dapX,dapY,dapZ],[dbpX,dbpY,dapZ],[dcpX,dcpY,dcpZ],[dalfaX,dalfaY,dalfaZ],[dbetaX,dbetaY,dbetaZ],[dVapX,dVapY,dVapZ],[dVbpX,dVbpY,dVcpZ],[dVcpX,dVbpY,dVcpZ]])
        A = A1.astype(np.float32)
        print "a matriz A e", A

        Xo = np.matrix([TX,TY,TZ])
        print "a matriz Xo e", Xo

        R = ajustaParametrico(Lb,P,Xo,Lo,A)
        print "o elemento R e", R
        Xa1 = R[0]
        print "o elemento 1 de R e", Xa1 
        sigmaXa1 = R[1]
        print "o elemento 2 de R e", sigmaXa1

        PrecX = sigmaXa1[0]
        print "o elemento PrecX e", PrecX
        PrecY = sigmaXa1[1]
        print "o elemento PrecY e", PrecY
        PrecZ = sigmaXa1[2]
        print "o elemento PrecZ e", PrecZ

        
        #PROPAGACAO
        # Calcula a distancia entre a estacao e os pontos:
        for i in range(len(AX)):
            distancia[i] = dist(AX[i],TX[0],AY[i],TY[0],AZ[i],TZ[0])
        print "as distancias entre os pontos T e Ai sao", distancia
        
        ### Calculo dos desvios-padrao das distancias:
        for i in range(len(distancia)):
            DesvDist[i] = desvDist(distancia[i])
        print "o desvio padrao da distancia entre AX1 e TX e  ", DesvDist
        
        ### Calculo dos angulos Zenitais:
        for i in range(len(AX)):
            AngZenital[i] = ang(AX[i],TX[0],AY[i],TY[0],AZ[i],TZ[0])
        print "os angulos zenitais entre os pontos T e Ai sao", AngZenital
        
        ### Calculo dos Azimutes:
        Azimutes = [0 for x in range(len(AX))]
        for i in range(len(AX)):
            Azimutes[i] = azimute(AX[i],AY[i],TX[0],TY[0])
        print "os azimutes sao", Azimutes
        
        
        ### Calculo dos angulos Horizontais:
        for i in range(len(AX)):
            AngHorizontal[i] = angHorizontal(Azimutes[0],Azimutes[i])
            print "os angulos horizontais entre os pontos T e Ai sao", AngHorizontal

        ###PROPAGACAO DE ERRO
        ###MONTAGEM DA MATRIZ G DAS DERIVADAS
        
        ### Calculo da derivada dxdd
        dxdd = [0 for x in range(len(AX))]
        for i in range(len(AX)):
            dxdd[i] = sin(pi/180*AngZenital[i])*sin(pi/180*Azimutes[i])
        print "derivadas dxdd", dxdd
        
        ### Calculo da derivada dxdv
        dxdv = [0 for x in range(len(AX))]
        for i in range(len(AX)):
            dxdv[i] = distancia[i]*cos(pi/180*AngZenital[i])*sin(pi/180*Azimutes[i])
        print "derivadas dxdv", dxdv
        
        ### Calculo da derivada dxda
        dxda = [0 for x in range(len(AX))]
        for i in range(len(AX)):
            dxda[i] = distancia[i]*sin(pi/180*AngZenital[i])*cos(pi/180*Azimutes[i])
        print "derivadas dxda", dxda 
        
        ### Calculo da derivada dxdx
        dxdx = [0 for x in range(len(AX))]
        for i in range(len(AX)):
            dxdx[i] = 1
        print "derivadas dxdx", dxdx 
        
        ### Calculo da derivada dxdy
        dxdy = [0 for x in range(len(AX))]
        for i in range(len(AX)):
            dxdy[i] = 0
        print "derivadas dxdy", dxdy
        
        ### Calculo da derivada dxdz
        dxdz = [0 for x in range(len(AX))]
        for i in range(len(AX)):
            dxdz[i] = 0
        print "derivadas dxdz", dxdz
        
        ### Calculo da derivada dydd
        dydd = [0 for x in range(len(AX))]
        for i in range(len(AX)):
            dydd[i] = sin(pi/180*AngZenital[i])*cos(pi/180*Azimutes[i])
        print "derivadas dydd", dydd
        
        ### Calculo da derivada dydv
        dydv = [0 for x in range(len(AX))]
        for i in range(len(AX)):
            dydv[i] = distancia[i]*cos(pi/180*AngZenital[i])*cos(pi/180*Azimutes[i])
        print "derivadas dydv", dydv
        
        ### Calculo da derivada dyda
        dyda = [0 for x in range(len(AX))]
        for i in range(len(AX)):
            dyda[i] = -distancia[i]*sin(pi/180*AngZenital[i])*sin(pi/180*Azimutes[i])
        print "derivadas dyda", dyda
        
        
        ### Calculo da derivada dydx
        dydx = [0 for x in range(len(AX))]
        for i in range(len(AX)):
            dydx[i] = 0
        print "derivadas dydx", dydx
        
        ### Calculo da derivada dydy
        dydy = [0 for x in range(len(AX))]
        for i in range(len(AX)):
            dydy[i] = 1
        print "derivadas dydy", dydy
        
        ### Calculo da derivada dydz
        dydz = [0 for x in range(len(AX))]
        for i in range(len(AX)):
            dydz[i] = 0
        print "derivadas dydz", dydz
        
        ### Calculo da derivada dzdd
        dzdd = [0 for x in range(len(AX))]
        for i in range(len(AX)):
            dzdd[i] = cos(pi/180*AngZenital[i])
        print "derivadas dzdd", dzdd
        
        ### Calculo da derivada dzdv
        dzdv = [0 for x in range(len(AX))]
        for i in range(len(AX)):
            dzdv[i] = -distancia[i]*sin(pi/180*AngZenital[i])
        print "derivadas dzdv", dzdv
        
        ### Calculo da derivada dzda
        dzda = [0 for x in range(len(AX))]
        for i in range(len(AX)):
            dzda[i] = 0
        print "derivadas dzda", dzda
        
        ### Calculo da derivada dzdx
        dzdx = [0 for x in range(len(AX))]
        for i in range(len(AX)):
            dzdx[i] = 0
        print "derivadas dzdx", dzdx
        
        ### Calculo da derivada dzdy
        dzdy = [0 for x in range(len(AX))]
        for i in range(len(AX)):
            dzdy[i] = 0
        print "derivadas dzdy", dzdy
        
        ### Calculo da derivada dzdz
        dzdz = [0 for x in range(len(AX))]
        for i in range(len(AX)):
            dzdz[i] = 1
        print "derivadas dzdz", dzdz
        
        n1 = float(n)
        ### Calculo das MVCs ----essa MVC que eu divido por n ??
        MVCs = [0 for x in range(len(DesvDist))]
        for i in range(len(DesvDist)):
                    MVCs[i] = np.array([[((DesvDist[i]**2)/n1),0,0,0,0,0],[0,((Prec_ang**2)/n1),0,0,0,0],[0,0,((Prec_ang**2)/n1),0,0,0],[0,0,0,((PrecX**2)/n1),0,0],[0,0,0,0,((PrecY**2)/n1),0],[0,0,0,0,0,((PrecZ**2)/n1)]])
        print "as MVCs sao", MVCs
        
        
        ### Calculo da MATRIZ G
        Gs = [0 for x in range(len(DesvDist))]
        for i in range(len(DesvDist)):
                    Gs[i] = np.array([[dxdd[i],dxdv[i],dxda[i],dxdx[i],dxdy[i],dxdz[i]],[dydd[i],dydv[i],dyda[i],dydx[i],dydy[i],dydz[i]],[dzdd[i],dzdv[i],dzda[i],dzdx[i],dzdy[i],dzdz[i]]])
        print "as Gs sao", Gs
        
        
        ### Calculo da MATRIZ MVCY
        GMVC = [0 for x in range(len(DesvDist))]
        MVCYs = [0 for x in range(len(DesvDist))]
        for i in range(len(DesvDist)):
            GMVC[i] = np.dot(Gs[i],MVCs[i])
            MVCYs[i] = np.dot(GMVC[i],Gs[i].T)
        print "as MVCYs sao", MVCYs
        
        
        ### Calculo do Dp
        Dps = [0 for x in range(len(DesvDist))]
        for i in range(len(DesvDist)):
            Dps[i] = [sqrt(MVCYs[i][0,0]),sqrt(MVCYs[i][1,1]),sqrt(MVCYs[i][2,2])]
        print "as Dps sao", Dps
        
        
        # Calcula a RESULTANTE 3D entre a estacao e os pontos:
        dpM = [0 for x in range(len(DesvDist))]
        for i in range(len(DesvDist)):
            dpM[i] = [sqrt((Dps[i][0]**2)+(Dps[i][1]**2)+(Dps[i][2]**2))]
        print "os valores médios entre os dps sao", dpM
        
        
        # Calcula a RESULTANTE 2D entre a estacao e os pontos:
        dpM2D = [0 for x in range(len(DesvDist))]
        for i in range(len(DesvDist)):
            dpM2D[i] = [sqrt((Dps[i][0]**2)+(Dps[i][1]**2))]
        print "os valores dpM2D sao", dpM2D
        
        
        # Cria uma lista com os valores dos desvio padrao em x:
        dpx = [0 for x in range(len(DesvDist))]
        for i in range(len(DesvDist)):
            dpx[i] = [Dps[i][0]]
        print "os valores dos desvios padrao em x sao", dpx

        # Cria uma lista com os valores dos desvio padrao em y:
        dpy = [0 for x in range(len(DesvDist))]
        for i in range(len(DesvDist)):
            dpy[i] = [Dps[i][1]]
        print "os valores dos desvios padrao em y sao", dpy
        
        # Cria uma lista com os valores dos desvio padrao em z:
        dpz = [0 for x in range(len(DesvDist))]
        for i in range(len(DesvDist)):
            dpz[i] = [Dps[i][2]]
        print "os valores dos desvios padrao em z sao", dpz
        
        
        conn = psycopg2.connect(database="PROPAGACAO_ERRO", host="localhost", user="postgres", password="postgres", port="5432")
        cur = conn.cursor()
            
        for i in range(len(dpx)):
            cur.execute ("""UPDATE PONTOS_ENCOSTA SET dpx = '%s' WHERE ponto = '%s'"""%(str(dpx[i]).strip('[]'),Ponto[i]))
            
            
        for i in range(len(dpx)):
            cur.execute ("""UPDATE PONTOS_ENCOSTA SET dpy = '%s' WHERE ponto = '%s'"""%(str(dpy[i]).strip('[]'),Ponto[i]))
            
            
        for i in range(len(dpx)):
            cur.execute ("""UPDATE PONTOS_ENCOSTA SET dpz = '%s' WHERE ponto = '%s'"""%(str(dpz[i]).strip('[]'),Ponto[i]))
            
        conn.commit()
        cur.close()
        conn.close()
        
        
        
