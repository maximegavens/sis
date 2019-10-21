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
import org.opengis.parameter.ParameterDescriptorGroup;
import org.opengis.parameter.ParameterNotFoundException;
import org.apache.sis.metadata.iso.citation.Citations;
import org.apache.sis.parameter.Parameters;
import org.apache.sis.referencing.operation.projection.NormalizedProjection;


/**
 * The provider for <cite>"Oblique Stereographic"</cite> projection (EPSG:9809).
 *
 * @author  Rueben Schulz (UBC)
 * @author  Martin Desruisseaux (Geomatys)
 * @version 0.8
 *
 * @see <a href="http://geotiff.maptools.org/proj_list/oblique_stereographic.html">GeoTIFF parameters for Oblique Stereographic</a>
 *
 * @since 0.6
 * @module
 */
@XmlTransient
public final class ObliqueStereographic extends AbstractStereographic {
    /**
     * For cross-version compatibility.
     */
    private static final long serialVersionUID = 6505988910141381354L;

    /**
     * The operation parameter descriptor for the <cite>Latitude of natural origin</cite> (φ₀) parameter value.
     * Valid values range is [-90 … 90]° and default value is 0°.
     *
     * <!-- Generated by ParameterNameTableGenerator -->
     * <table class="sis">
     *   <caption>Parameter names</caption>
     *   <tr><td> EPSG:    </td><td> Latitude of natural origin </td></tr>
     *   <tr><td> OGC:     </td><td> latitude_of_origin </td></tr>
     *   <tr><td> ESRI:    </td><td> Latitude_Of_Origin </td></tr>
     *   <tr><td> NetCDF:  </td><td> latitude_of_projection_origin </td></tr>
     *   <tr><td> GeoTIFF: </td><td> NatOriginLat </td></tr>
     *   <tr><td> Proj4:   </td><td> lat_0 </td></tr>
     * </table>
     */
    public static final ParameterDescriptor<Double> LATITUDE_OF_ORIGIN = TransverseMercator.LATITUDE_OF_ORIGIN;

    /**
     * The operation parameter descriptor for the <cite>Longitude of natural origin</cite> (λ₀) parameter value.
     * Valid values range is [-180 … 180]° and default value is 0°.
     *
     * <!-- Generated by ParameterNameTableGenerator -->
     * <table class="sis">
     *   <caption>Parameter names</caption>
     *   <tr><td> EPSG:    </td><td> Longitude of natural origin </td></tr>
     *   <tr><td> OGC:     </td><td> central_meridian </td></tr>
     *   <tr><td> ESRI:    </td><td> Central_Meridian </td></tr>
     *   <tr><td> NetCDF:  </td><td> longitude_of_projection_origin </td></tr>
     *   <tr><td> GeoTIFF: </td><td> NatOriginLong </td></tr>
     *   <tr><td> Proj4:   </td><td> lon_0 </td></tr>
     * </table>
     */
    public static final ParameterDescriptor<Double> LONGITUDE_OF_ORIGIN = Mercator1SP.LONGITUDE_OF_ORIGIN;

    /**
     * The operation parameter descriptor for the <cite>Scale factor at natural origin</cite> (k₀) parameter value.
     * Valid values range is (0 … ∞) and default value is 1.
     *
     * <!-- Generated by ParameterNameTableGenerator -->
     * <table class="sis">
     *   <caption>Parameter names</caption>
     *   <tr><td> EPSG:    </td><td> Scale factor at natural origin </td></tr>
     *   <tr><td> OGC:     </td><td> scale_factor </td></tr>
     *   <tr><td> ESRI:    </td><td> Scale_Factor </td></tr>
     *   <tr><td> NetCDF:  </td><td> scale_factor_at_projection_origin </td></tr>
     *   <tr><td> GeoTIFF: </td><td> ScaleAtNatOrigin </td></tr>
     *   <tr><td> Proj4:   </td><td> k </td></tr>
     * </table>
     */
    public static final ParameterDescriptor<Double> SCALE_FACTOR = Mercator1SP.SCALE_FACTOR;    // Same as PolarStereographicA.

    /**
     * The group of all parameters expected by this coordinate operation.
     */
    private static final ParameterDescriptorGroup PARAMETERS;
    static {
        PARAMETERS = builder()
                .addIdentifier(             "9809")
                .addName(                   "Oblique Stereographic")
                .addName(Citations.OGC,     "Oblique_Stereographic")
                .addName(Citations.ESRI,    "Double_Stereographic")
                .addName(Citations.GEOTIFF, "CT_ObliqueStereographic")
                .addName(Citations.S57,     "Oblique stereographic")
                .addName(Citations.S57,     "OST")
                .addName(Citations.PROJ4,   "sterea")
                .addName(                   "Roussilhe")
                .addIdentifier(Citations.GEOTIFF, "16")
                .addIdentifier(Citations.S57,     "14")
                .createGroupForMapProjection(
                        LATITUDE_OF_ORIGIN,
                        LONGITUDE_OF_ORIGIN,
                        SCALE_FACTOR,
                        FALSE_EASTING,
                        FALSE_NORTHING);
    }

    /**
     * Constructs a new provider.
     */
    public ObliqueStereographic() {
        super(PARAMETERS);
    }

    /**
     * {@inheritDoc}
     *
     * @return the map projection created from the given parameter values.
     */
    @Override
    protected NormalizedProjection createProjection(Parameters parameters) throws ParameterNotFoundException {
        return new org.apache.sis.referencing.operation.projection.ObliqueStereographic(this, parameters);
    }
}
