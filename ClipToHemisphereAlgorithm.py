# -*- coding: utf-8 -*-

"""
Date                 : July 2019
Copyright            : (C) 2019 by Woravich Kumthonkittikul
Email                : woravich@hotmail.com
Adapted from         : Juernjakob Dugge (2016)
***************************************************************************
*                                                                         *
*   This program is free software; you can redistribute it and/or modify  *
*   it under the terms of the GNU General Public License as published by  *
*   the Free Software Foundation; either version 2 of the License, or     *
*   (at your option) any later version.                                   *
*                                                                         *
***************************************************************************
"""

__author__ = 'Woravich Kumthonkittikul'
__date__ = 'July 2019'
__copyright__ = '(C) 2019, Woravich Kumthonkittikul'
__developed_from = 'Juernjakob Dugge (2016)'

from PyQt5.QtCore import QCoreApplication
from qgis.core import (QgsProcessing,
                       QgsProcessingException,
                       QgsProcessingAlgorithm,
                       QgsProcessingParameterFeatureSource,
                       QgsProcessingParameterNumber,
                       QgsProcessingParameterVectorDestination,
                       QgsCoordinateReferenceSystem,
                       QgsCoordinateTransform,
                       QgsProject,
                       QgsVectorLayer,
                       QgsPointXY,
                       QgsFeature,
                       QgsGeometry,
                       QgsProject)
import processing
import numpy as np
import cmath

class ClipToHemisphereAlgorithm(QgsProcessingAlgorithm):
    """
    This is an example algorithm that takes a vector layer and
    creates a new identical one.

    It is meant to be used as an example of how to create your own
    algorithms and explain methods and variables used to do it. An
    algorithm like this will be available in all elements, and there
    is not need for additional work.

    All Processing algorithms should extend the QgsProcessingAlgorithm
    class.
    """

    # Constants used to refer to parameters and outputs. They will be
    # used when calling the algorithm from another algorithm, or when
    # calling from the QGIS console.

    INPUT = 'INPUT'
    OUTPUT = 'OUTPUT'
    CENTER_LATITUDE = 'CENTER_LATITUDE'
    CENTER_LONGITUDE = 'CENTER_LONGITUDE'
    SEGMENTS = 'SEGMENTS'

    def tr(self, string):
        """
        Returns a translatable string with the self.tr() function.
        """
        return QCoreApplication.translate('Processing', string)

    def createInstance(self):
        return ClipToHemisphereAlgorithm()

    def name(self):
        """
        Returns the algorithm name, used for identifying the algorithm. This
        string should be fixed for the algorithm, and must not be localised.
        The name should be unique within each provider. Names should contain
        lowercase alphanumeric characters only and no spaces or other
        formatting characters.
        """
        return 'clip_vector_to_hemishere'

    def displayName(self):
        """
        Returns the translated algorithm name, which should be used for any
        user-visible display of the algorithm name.
        """
        return self.tr('Clip a vector layer to the hemisphere centred on a user specified point')

    def group(self):
        """
        Returns the name of the group this algorithm belongs to. This string
        should be localised.
        """
        return self.tr('Clip to Hemisphere')

    def groupId(self):
        """
        Returns the unique ID of the group this algorithm belongs to. This
        string should be fixed for the algorithm, and must not be localised.
        The group id should be unique within each provider. Group id should
        contain lowercase alphanumeric characters only and no spaces or other
        formatting characters.
        """
        return 'clip_to_hemisphere'

    def shortHelpString(self):
        """
        Returns a localised short helper string for the algorithm. This string
        should provide a basic description about what the algorithm does and the
        parameters and outputs associated with it..
        """
        return self.tr("Clip a vector layer for being used in a custom orthographic projection (world from space) defining the centre of the projection")

    def initAlgorithm(self, config=None):
        """
        Here we define the inputs and output of the algorithm, along
        with some other properties.
        """

        # We add the input vector features source. It can have any kind of
        # geometry.
        self.addParameter(
            QgsProcessingParameterFeatureSource(
                self.INPUT,
                self.tr('Input layer'),
                [QgsProcessing.TypeVectorAnyGeometry]
            )
        )
        self.addParameter(
            QgsProcessingParameterNumber(
                self.CENTER_LATITUDE,
                self.tr('Latitude of center of hemisphere'),
                QgsProcessingParameterNumber.Double,
				defaultValue = 0.0
            )
        )
        self.addParameter(
            QgsProcessingParameterNumber(
                self.CENTER_LONGITUDE,
                self.tr('Longitude of center of hemisphere'),
                QgsProcessingParameterNumber.Double,
				defaultValue = 0.0
            )
        )
        self.addParameter(
            QgsProcessingParameterNumber(
                self.SEGMENTS,
                self.tr('Number of segments for approximating the hemisphere'),
				defaultValue = 500
            )
        )
        

        # We add a feature sink in which to store our processed features (this
        # usually takes the form of a newly created vector layer when the
        # algorithm is run in QGIS).
        self.addParameter(
            QgsProcessingParameterVectorDestination(
                self.OUTPUT,
                self.tr('Output layer')
            )
        )

    def processAlgorithm(self, parameters, context, feedback):
        """
        Here is where the processing itself takes place.
        """
        earthRadius = 6370997
        
        # Retrieve the feature source and sink. The 'dest_id' variable is used
        # to uniquely identify the feature sink, and must be included in the
        # dictionary returned by the processAlgorithm function.
        source = self.parameterAsVectorLayer(
            parameters,
            self.INPUT,
            context
        )
        
        #lat
        centerLatitude = self.parameterAsDouble(
            parameters,
            self.CENTER_LATITUDE,
            context
        )
        
        #lon
        centerLongitude = self.parameterAsDouble(
            parameters,
            self.CENTER_LONGITUDE,
            context
        )
        
        #segment
        segments = self.parameterAsInt(
            parameters,
            self.SEGMENTS,
            context
        )
        

        # If source was not found, throw an exception to indicate that the algorithm
        # encountered a fatal error. The exception text can be any string, but in this
        # case we use the pre-built invalidSourceError method to return a standard
        # helper text for when a source cannot be evaluated
        if source is None:
            raise QgsProcessingException(self.invalidSourceError(parameters, self.INPUT))
        
        if centerLatitude is None:
            raise QgsProcessingException(self.invalidSourceError(parameters, self.CENTER_LATITUDE))
        
        if centerLongitude is None:
            raise QgsProcessingException(self.invalidSourceError(parameters, self.CENTER_LATITUDE))
        
        if segments is None:
            raise QgsProcessingException(self.invalidSourceError(parameters, self.SEGMENTS))

        output = self.parameterAsOutputLayer(
            parameters,
            self.OUTPUT,
            context
        )
        
        # Send some information to the user
        if source.sourceCrs().authid() != 'EPSG:4326':
            feedback.pushDebugInfo('Input layer for "Clip to Hemisphere" does not use the WGS84 (EPSG:4326) CRS. This can cause unexpected results.')
        else:
            feedback.pushInfo('CRS is {}'.format(source.sourceCrs().authid()))
            
        sourceCrs = source.sourceCrs()
        targetProjString = "+proj=ortho +lat_0=" + str(centerLatitude) + \
            " +lon_0=" + str(centerLongitude) + \
            " +x_0=0 +y_0=0 +a=" + str(earthRadius) + \
            " +b=" + str(earthRadius) + \
            " +units=m +no_defs"
        targetCrs = QgsCoordinateReferenceSystem()
        targetCrs.createFromProj4(targetProjString)
        
        transformTargetToSrc = QgsCoordinateTransform(targetCrs,sourceCrs,QgsProject.instance())
        transformSrcToTarget = QgsCoordinateTransform(sourceCrs,targetCrs,QgsProject.instance())
        
        clipLayer = QgsVectorLayer("MultiPolygon?crs=epsg:4326", "clipLayer", "memory")
        pr = clipLayer.dataProvider()
        
        """
        This part was adapted from (C) 2016 by Juernjakob Dugge
        Source: https://plugins.qgis.org/plugins/ClipToHemisphere/
        """
        # Handle edge cases:
        # Hemisphere centered on the equator
        if centerLatitude == 0:
            # Hemisphere centered on the equator and including the antimeridian
            if abs(centerLongitude) >= 90:
                edgeEast = -180 - np.sign(centerLongitude) * \
                        (180 - abs(centerLongitude)) + 90
                edgeWest = 180 - np.sign(centerLongitude) * \
                        (180 - abs(centerLongitude)) - 90
                circlePoints = [[
                    [QgsPointXY(-180.01, latitude) for
                        latitude in np.linspace(90, -90, segments / 8)] +
                    [QgsPointXY(longitude, -90) for longitude in
                        np.linspace(-180, edgeEast, segments / 8)] +
                    [QgsPointXY(edgeEast, latitude) for latitude in
                        np.linspace(-90, 90, segments / 8)] +
                    [QgsPointXY(longitude, 90) for longitude in
                        np.linspace(edgeEast, -180, segments / 8)]
                    ],
                    [
                    [QgsPointXY(edgeWest, latitude) for latitude in
                        np.linspace(90, -90, segments / 8)] +
                    [QgsPointXY(longitude, -90) for longitude in
                        np.linspace(edgeWest, 180, segments / 8)] +
                    [QgsPointXY(180.01, latitude) for
                        latitude in np.linspace(-90, 90, segments / 8)] +
                    [QgsPointXY(longitude, 90) for longitude in
                        np.linspace(180, edgeWest, segments / 8)]
                    ]]
            # Hemisphere centered on the equator not including the antimeridian
            else:
                edgeWest = centerLongitude - 90
                edgeEast = centerLongitude + 90
                circlePoints = [[
                    [QgsPointXY(edgeWest, latitude) for latitude in
                        np.linspace(90, -90, segments / 4)] +
                    [QgsPointXY(longitude, -90) for longitude in
                        np.linspace(edgeWest, edgeEast, segments / 4)] +
                    [QgsPointXY(edgeEast, latitude) for
                        latitude in np.linspace(-90, 90, segments / 4)] +
                    [QgsPointXY(longitude, 90) for longitude in
                        np.linspace(edgeEast, edgeWest, segments / 4)]
                    ]]
        # Hemisphere centered on one of the poles
        elif abs(centerLatitude) == 90:
            circlePoints = [[
                [QgsPointXY(-180.01, latitude) for latitude in
                        np.linspace(45 + 0.5 * centerLatitude,
                                   -45 + 0.5 * centerLatitude,
                                   segments / 4)] +
                [QgsPointXY(longitude, -45 + 0.5 * centerLatitude)
                        for longitude in
                        np.linspace(-180, 180, segments / 4)] +
                [QgsPointXY(180.01, latitude) for latitude in
                        np.linspace(-45 + 0.5 * centerLatitude,
                                     45 + 0.5 * centerLatitude,
                                   segments / 4)] +
                [QgsPointXY(longitude, 45 + 0.5 * centerLatitude) for longitude in
                        np.linspace(180, -180, segments / 4)]
                ]]
        # All other hemispheres
        else:
            # Create a circle in the orthographic projection, convert the
            # circle coordinates to the source CRS
            angles = np.linspace(0, 2 * np.pi, segments, endpoint=False)
            circlePoints = np.array([
                transformTargetToSrc.transform(
                    QgsPointXY(np.cos(angle) * earthRadius * 0.9999,
                             np.sin(angle) * earthRadius * 0.9999)
                            ) for angle in angles
            ])

            # Sort the projected circle coordinates from west to east
            sortIdx = np.argsort(circlePoints[:, 0])
            circlePoints = circlePoints[sortIdx, :]
            circlePoints = [[[QgsPointXY(point[0], point[1])
                for point in circlePoints]]]

            # Find the circle point in the orthographic projection that lies
            # on the antimeridian by linearly interpolating the angles of the
            # first and last coordinates
            startGap = 180 + circlePoints[0][0][0][0]
            endGap = 180 - circlePoints[0][0][-1][0]
            totalGap = startGap + endGap
            startCoordinates = transformSrcToTarget.transform(circlePoints[0][0][0])
            endCoordinates = transformSrcToTarget.transform(circlePoints[0][0][-1])
            startAngle = np.arctan2(startCoordinates[0], startCoordinates[1])
            endAngle = np.arctan2(endCoordinates[0], endCoordinates[1])
            antimeridianAngle = cmath.phase(
                endGap / totalGap * cmath.rect(1, startAngle) +
                startGap / totalGap * cmath.rect(1, endAngle))
            antimeridianPoint = transformTargetToSrc.transform(QgsPointXY(
                np.sin(antimeridianAngle) * earthRadius * 0.9999,
                np.cos(antimeridianAngle) * earthRadius * 0.9999
                ))

            # Close the polygon
            circlePoints[0][0].extend(
                [QgsPointXY(180.01, latitude) for latitude in
                        np.linspace(antimeridianPoint[1],
                            np.sign(centerLatitude) * 90, segments / 4)] +
                [QgsPointXY(-180.01, latitude) for latitude in
                        np.linspace(np.sign(centerLatitude) * 90,
                            antimeridianPoint[1], segments / 4)]
                )
        
        # Create the feature and add it to the layer
        circle = QgsFeature()
        circle.setGeometry(QgsGeometry.fromMultiPolygonXY(circlePoints))

        pr.addFeatures([circle])
        pr.updateExtents()
        
        # We need to add the clipping layer to the layer list in order to be
        # able to use them with processing.runalg()
        clipLayerReg = QgsProject.instance().addMapLayer(clipLayer)
        
        # If sink was not created, throw an exception to indicate that the algorithm
        # encountered a fatal error. The exception text can be any string, but in this
        # case we use the pre-built invalidSinkError method to return a standard
        # helper text for when a sink cannot be evaluated
        if output is None:
            raise QgsProcessingException(self.invalidSinkError(parameters, self.OUTPUT))
        
        #clip the layer
        processing.run("qgis:clip", 
            {'INPUT': source,
             'OVERLAY': clipLayerReg, 
             'OUTPUT':output},
            feedback = feedback,
            context = context)
        

        


        return {self.OUTPUT: output}
