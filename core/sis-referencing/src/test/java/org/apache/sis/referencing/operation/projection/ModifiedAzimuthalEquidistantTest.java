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
package org.apache.sis.referencing.operation.projection;

import org.apache.sis.internal.referencing.Formulas;
import org.opengis.util.FactoryException;
import org.opengis.referencing.operation.TransformException;
import org.apache.sis.internal.referencing.provider.MapProjection;
import org.apache.sis.test.DependsOn;
import org.junit.Test;

import static java.lang.StrictMath.toRadians;


/**
 * Tests the {@link ModifiedAzimuthalEquidistant} class.
 *
 * @author  Martin Desruisseaux (Geomatys)
 * @author  Maxime Gavens (Geomatys)
 * @version 1.1
 * @since   1.1
 * @module
 */
@DependsOn(NormalizedProjectionTest.class)
public final strictfp class ModifiedAzimuthalEquidistantTest extends AzimuthalEquidistantTest {
    /**
     * Returns the method to be tested.
     */
    @Override
    MapProjection method() {
        return new org.apache.sis.internal.referencing.provider.ModifiedAzimuthalEquidistant();
    }

    /**
     * Tests the projection on a sphere. We override the method provides in parent class
     * because the point provided in Snyder is too far from projection centre.
     *
     * @throws FactoryException if an error occurred while creating the projection.
     * @throws TransformException if an error occurred while projecting the test point.
     */
    @Test
    @Override
    public void testSpherical() throws FactoryException, TransformException {
        tolerance = 20;                     // Same tolerance than in parent class.
        final double r = 6357767.51;        // Conformal sphere radius at the latitude being tested.
        testWithEPSG(r, r);
    }

    /**
     * Tests with the point published in EPSG guidance note.
     *
     * @throws FactoryException if an error occurred while creating the projection.
     * @throws TransformException if an error occurred while projecting the test point.
     */
    @Test
    @Override
    public void testWithEPSG() throws FactoryException, TransformException {
        tolerance = 0.01;
        testWithEPSG(CLARKE_A, CLARKE_B);
    }

    /**
     * Tests the <cite>"Modified Azimuthal Equidistant"</cite> (EPSG:9832) projection method.
     * This test is defined in GeoAPI conformance test suite.
     *
     * @throws FactoryException if an error occurred while creating the map projection.
     * @throws TransformException if an error occurred while projecting a coordinate.
     *
     * @see org.opengis.test.referencing.ParameterizedTransformTest#testModifiedAzimuthalEquidistant()
     */
    @Test
    public void runGeoapiTest() throws FactoryException, TransformException {
        createGeoApiTestNoDerivatives(method()).testModifiedAzimuthalEquidistant();
    }

    /**
     * Tests the derivatives at a few points on an ellipsoid. This method compares the derivatives computed
     * by the projection with an estimation of derivatives computed by the finite differences method.
     *
     * @throws FactoryException if an error occurred while creating the map projection.
     * @throws TransformException if an error occurred while projecting a point.
     */
    @Test
    @Override
    public void testDerivative() throws FactoryException, TransformException {
        createCompleteProjection(method(),
                CLARKE_A,
                CLARKE_B,
                40,                 // Longitude of natural origin (central-meridian)
                25,                 // Latitude of natural origin
                Double.NaN,         // Standard parallel 1
                Double.NaN,         // Standard parallel 2
                Double.NaN,         // Scale factor
                40000,              // False easting
                60000);             // False Northing
        final double delta = (1.0 / 60) / 1852;                 // Approximately 1 metre.
        derivativeDeltas = new double[] {delta, delta};
        tolerance = Formulas.LINEAR_TOLERANCE / 10;
        verifyDerivative(toRadians(-4), toRadians(40));
        verifyDerivative(toRadians(10), toRadians(60));
        verifyInverse(toRadians(15), toRadians(25));
        verifyDerivative(30, 27);
        verifyDerivative(27, 20);
        verifyDerivative(105,  30);
        verifyDerivative(100, -60);
        verifyDerivative(-100,  3);
        verifyDerivative( -56, 50);
        verifyDerivative( -20, 47);
    }
}
