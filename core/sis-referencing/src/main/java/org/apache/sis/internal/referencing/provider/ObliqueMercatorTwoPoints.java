/*
 * Licensed to the Apache Software Foundation (ASF) under one or more
 * contributor license agreements.  See the NOTICE file distributed with
 * this work for additional information regarding copyright ownership.
 * The ASF licenses this file to You under the Apache License, Version 2.0
 * (the "License"); you may not use this file except in compliance with
 * the License.  You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
package org.apache.sis.internal.referencing.provider;

import javax.xml.bind.annotation.XmlTransient;
import org.opengis.parameter.ParameterDescriptor;
import org.apache.sis.parameter.ParameterBuilder;
import org.apache.sis.metadata.iso.citation.Citations;


/**
 * The provider for "<cite>Oblique Mercator</cite>" projection specified by two points on the central line.
 * This is different than the classical {@linkplain ObliqueMercator Oblique Mercator}, which uses a central
 * point and azimuth. This projection is specific to ESRI software.
 *
 * @author  Rueben Schulz (UBC)
 * @author  Martin Desruisseaux (Geomatys)
 * @version 1.0
 * @since   1.0
 * @module
 */
@XmlTransient
public class ObliqueMercatorTwoPoints extends ObliqueMercator {
    /**
     * For compatibility with different versions during deserialization.
     */
    private static final long serialVersionUID = 7202109026784761711L;

    /**
     * The operation parameter descriptor for the latitude / longitude of one of the two points.
     */
    public static final ParameterDescriptor<Double> LAT_OF_1ST_POINT, LONG_OF_1ST_POINT,
                                                    LAT_OF_2ND_POINT, LONG_OF_2ND_POINT;

    static {
        final ParameterBuilder builder = builder().setCodeSpace(Citations.ESRI, "ESRI");
        LAT_OF_1ST_POINT  = createLatitude (builder.addName("Latitude_Of_1st_Point"), false);
        LAT_OF_2ND_POINT  = createLatitude (builder.addName("Latitude_Of_2nd_Point"), false);
        LONG_OF_1ST_POINT = createLongitude(builder.addName("Longitude_Of_1st_Point"));
        LONG_OF_2ND_POINT = createLongitude(builder.addName("Longitude_Of_2nd_Point"));
    }

    /**
     * Constructs a new provider for a projection with false easting/northing at coordinate system natural origin.
     */
    public ObliqueMercatorTwoPoints() {
        this("Hotine_Oblique_Mercator_Two_Point_Natural_Origin", FALSE_EASTING, FALSE_NORTHING);
    }

    /**
     * Constructs a new provider for a projection of the given name. Used for the two variants:
     * with easting/northing values at coordinate system natural origin or at projection center.
     */
    ObliqueMercatorTwoPoints(final String name,
                             final ParameterDescriptor<Double> easting,
                             final ParameterDescriptor<Double> northing)
    {
        super(builder().setCodeSpace(Citations.ESRI, "ESRI").addName(name)
                .createGroupForMapProjection(
                        LAT_OF_1ST_POINT,    LONG_OF_1ST_POINT,
                        LAT_OF_2ND_POINT,    LONG_OF_2ND_POINT,
                        LATITUDE_OF_CENTRE,  SCALE_FACTOR,
                        easting, northing));
    }
}
