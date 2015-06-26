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
package org.apache.sis.io.wkt;

import org.opengis.referencing.crs.*;
import org.opengis.test.wkt.CRSParserTest;
import org.apache.sis.test.DependsOn;
import org.junit.runner.RunWith;
import org.junit.runners.JUnit4;


/**
 * Tests Well-Known Text parser using the tests defined in GeoAPI. Those tests use the
 * {@link org.apache.sis.referencing.factory.GeodeticObjectFactory#createFromWKT(String)} method.
 *
 * @author  Martin Desruisseaux (IRD, Geomatys)
 * @since   0.6
 * @version 0.6
 * @module
 */
@RunWith(JUnit4.class)
@DependsOn(GeodeticObjectParserTest.class)
public class WKTParserTest extends CRSParserTest {
    /**
     * Creates a new test case using the default {@code CRSFactory} implementation.
     */
    public WKTParserTest() {
        super(org.apache.sis.internal.system.DefaultFactories.forClass(CRSFactory.class));
    }
}