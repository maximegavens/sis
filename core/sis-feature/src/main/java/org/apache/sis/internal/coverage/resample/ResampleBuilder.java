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
package org.apache.sis.internal.coverage.resample;

import java.awt.image.RenderedImage;
import java.awt.image.WritableRenderedImage;
import java.util.function.BiFunction;
import org.apache.sis.coverage.grid.GridCoverage;
import org.apache.sis.coverage.grid.GridGeometry;
import org.opengis.referencing.crs.CoordinateReferenceSystem;

/**
 *
 * @author Martin Desruisseaux (Geomatys)
 * @author Johann Sorel (Geomatys)
 */
public class ResampleBuilder {

    public ResampleBuilder setSource(GridCoverage coverage) {
        throw new UnsupportedOperationException("Not implemented");
    }

    public ResampleBuilder setSource(GridGeometry grid, RenderedImage image) {
        throw new UnsupportedOperationException("Not implemented");
    }

    public ResampleBuilder setTarget(GridCoverage coverage) {
        throw new UnsupportedOperationException("Not implemented");
    }

    public ResampleBuilder setTarget(CoordinateReferenceSystem crs) {
        throw new UnsupportedOperationException("Not implemented");
    }

    public ResampleBuilder setTarget(GridGeometry grid) {
        throw new UnsupportedOperationException("Not implemented");
    }

    /**
     * Resample to target coverage image.
     *
     * @param grid
     * @param image
     * @return
     */
    public ResampleBuilder setTarget(GridGeometry grid, WritableRenderedImage image) {
        throw new UnsupportedOperationException("Not implemented");
    }

    public ResampleBuilder setInterpolation(Interpolator interpolator) {
        throw new UnsupportedOperationException("Not implemented");
    }

    /**
     * Set target value used to fill samples when source datas are out of coverage extent.
     *
     * @param fill
     * @return
     */
    public ResampleBuilder setFillValue(double ... fill) {
        throw new UnsupportedOperationException("Not implemented");
    }

    /**
     * When resampling multiple coverage onto a single image, how values should be merge
     * may depend on the algorithm.
     * For example if a value is already defined in the target image and the source
     * image provide NaN values we could want to preserve the value in the image instead of the NaN.
     *
     * By default source samples will override any value in the target image.
     *
     * @param merger
     * @return
     */
    public ResampleBuilder setMerger(BiFunction<double[], double[], double[]> merger) {
        throw new UnsupportedOperationException("Not implemented");
    }

    public Resample build() {
        throw new UnsupportedOperationException("Not implemented");
    }

}
