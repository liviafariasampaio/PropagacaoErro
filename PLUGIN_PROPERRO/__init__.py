# -*- coding: utf-8 -*-
"""
/***************************************************************************
 PROPAGACAO_ERRO
                                 A QGIS plugin
 PROPAGACAO_ERRO
                             -------------------
        begin                : 2018-08-06
        copyright            : (C) 2018 by Livia Faria Sampaio
        email                : livias_93@hotmail.com
        git sha              : $Format:%H$
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
 This script initializes the plugin, making it known to QGIS.
"""


# noinspection PyPep8Naming
def classFactory(iface):  # pylint: disable=invalid-name
    """Load PROPAGACAO_ERRO class from file PROPAGACAO_ERRO.

    :param iface: A QGIS interface instance.
    :type iface: QgisInterface
    """
    #
    from .PROPAGACAO_ERRO import PROPAGACAO_ERRO
    return PROPAGACAO_ERRO(iface)
