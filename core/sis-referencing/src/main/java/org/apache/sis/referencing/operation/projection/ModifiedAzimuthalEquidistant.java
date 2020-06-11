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

import java.util.EnumMap;

import org.apache.sis.referencing.operation.matrix.Matrix2;
import org.opengis.util.FactoryException;
import org.opengis.referencing.operation.Matrix;
import org.opengis.referencing.operation.MathTransform;
import org.opengis.referencing.operation.MathTransformFactory;
import org.opengis.referencing.operation.OperationMethod;
import org.opengis.parameter.ParameterDescriptor;
import org.apache.sis.parameter.Parameters;
import org.apache.sis.referencing.operation.transform.ContextualParameters;
import org.apache.sis.referencing.operation.matrix.MatrixSIS;
import org.apache.sis.util.Workaround;

import static java.lang.Math.*;
import static org.apache.sis.internal.referencing.provider.ModifiedAzimuthalEquidistant.*;


/**
 * <cite>Modified Azimuthal Equidistant</cite> projection (EPSG:9832).
 * This is an approximation of the oblique form of the <cite>Azimuthal Equidistant</cite> projection.
 * For distances under 800 kilometres this modification introduces no significant error.
 *
 * @author  Martin Desruisseaux (Geomatys)
 * @author  Maxime Gavens (Geomatys)
 * @version 1.1
 *
 * @see AzimuthalEquidistant
 *
 * @since 1.1
 * @module
 */
public class ModifiedAzimuthalEquidistant extends AzimuthalEquidistant {
    /**
     * For compatibility with different versions during deserialization.
     */
    private static final long serialVersionUID = 96569177715708509L;

    /**
     * A term involving radius of curvature ν₀, the latitude of origin φ₀ and the eccentricity.
     * The semi-major axis length <var>a</var> is omitted since it is handled outside this class.
     */
    private final double ℯ2_ν0_sinφ0;

    /**
     * The ℯ⋅sin(φ₀)/√(1 − ℯ²) term, used in direct projection.
     */
    private final double G;

    /**
     * The ℯ⋅cos(φ₀)/√(1 − ℯ²) term. This is the <var>H</var> term in EPSG guidance notes
     * but without the cos(α) term (omitted because α depends on the point to project).
     *
     * <p>Note that during inverse projection, EPSG guidance notes has a <var>A</var> as:
     * −ℯ²⋅cos²φ₀/(1 − ℯ²)⋅cos²α. We opportunistically use Hp² for that purpose.</p>
     */
    private final double Hp;

    /**
     * The 3⋅ℯ²⋅sin(φ₀)⋅cos(φ₀)/(1 − ℯ²) term. This is the <var>B</var> term in EPSG guidance notes
     * for inverse projection but without the terms that depend on coordinates of transformed point.
     */
    private final double Bp;

    /**
     * Work around for RFE #4093999 in Sun's bug database
     * ("Relax constraint on placement of this()/super() call in constructors").
     */
    @Workaround(library="JDK", version="1.8")
    private static Initializer initializer(final OperationMethod method, final Parameters parameters) {
        final EnumMap<ParameterRole, ParameterDescriptor<Double>> roles = new EnumMap<>(ParameterRole.class);
        roles.put(ParameterRole.CENTRAL_MERIDIAN, LONGITUDE_OF_ORIGIN);
        roles.put(ParameterRole.FALSE_EASTING,    FALSE_EASTING);
        roles.put(ParameterRole.FALSE_NORTHING,   FALSE_NORTHING);
        return new Initializer(method, parameters, roles, STANDARD_VARIANT);
    }

    /**
     * Creates a Modified Azimuthal Equidistant projection from the given parameters.
     * The {@code method} argument can be the description of one of the following:
     *
     * <ul>
     *   <li><cite>"Modified Azimuthal Equidistant"</cite>.</li>
     * </ul>
     *
     * @param method      description of the projection parameters.
     * @param parameters  the parameter values of the projection to create.
     */
    public ModifiedAzimuthalEquidistant(final OperationMethod method, final Parameters parameters) {
        this(initializer(method, parameters));
    }

    /**
     * Work around for RFE #4093999 in Sun's bug database
     * ("Relax constraint on placement of this()/super() call in constructors").
     */
    @Workaround(library="JDK", version="1.8")
    private ModifiedAzimuthalEquidistant(final Initializer initializer) {
        super(initializer);
        final double axisRatio, ν0, f;
        axisRatio   = initializer.axisLengthRatio().doubleValue();
        ν0          = initializer.radiusOfCurvature(sinφ0);
        ℯ2_ν0_sinφ0 = eccentricitySquared * ν0 * sinφ0;
        f           = eccentricity / axisRatio;                 // √(1 - ℯ²) = b/a
        G           = f * sinφ0;
        Hp          = f * cosφ0;
        Bp          = 3*eccentricitySquared * (sinφ0*cosφ0) / (1 - eccentricitySquared);

        final MatrixSIS denormalize = context.getMatrix(ContextualParameters.MatrixRole.DENORMALIZATION);
        denormalize.convertBefore(0, ν0, null);
        denormalize.convertBefore(1, ν0, null);
    }

    /**
     * Returns the sequence of <cite>normalization</cite> → {@code this} → <cite>denormalization</cite> transforms
     * as a whole. The transform returned by this method expects (<var>longitude</var>, <var>latitude</var>)
     * coordinates in <em>degrees</em> and returns (<var>x</var>,<var>y</var>) coordinates in <em>metres</em>.
     *
     * <p>The non-linear part of the returned transform will be {@code this} transform, except if the ellipsoid
     * is spherical. In the later case, {@code this} transform will be replaced by a simplified implementation.</p>
     *
     * @param  factory  the factory to use for creating the transform.
     * @return the map projection from (λ,φ) to (<var>x</var>,<var>y</var>) coordinates.
     * @throws FactoryException if an error occurred while creating a transform.
     */
    @Override
    public MathTransform createMapProjection(final MathTransformFactory factory) throws FactoryException {
        AzimuthalEquidistant kernel = this;
        if (eccentricity == 0) {
            kernel = new AzimuthalEquidistant(this);
        }
        return context.completeTransform(factory, kernel);
    }

    /**
     * Converts the specified (λ,φ) coordinate and stores the (<var>x</var>,<var>y</var>) result in {@code dstPts}.
     *
     * @return the matrix of the projection derivative at the given source position,
     *         or {@code null} if the {@code derivate} argument is {@code false}.
     * @throws ProjectionException if the coordinate can not be converted.
     */
    @Override
    public Matrix transform(final double[] srcPts, final int srcOff,
                            final double[] dstPts, final int dstOff,
                            final boolean derivate) throws ProjectionException
    {
        final double λ     = srcPts[srcOff  ];
        final double φ     = srcPts[srcOff+1];
        final double cosλ  = cos(λ);
        final double sinλ  = sin(λ);
        final double sinλ2 = sinλ * sinλ;
        final double cosφ  = cos(φ);
        final double sinφ  = sin(φ);
        final double rν    = sqrt(1 - eccentricitySquared*(sinφ*sinφ));
        final double tanΨ  = ((1 - eccentricitySquared)*sinφ + ℯ2_ν0_sinφ0*rν) / cosφ;
        final double rcosΨ = sqrt(1 + tanΨ*tanΨ);
        final double α     = atan2(sinλ, cosφ0*tanΨ - sinφ0*cosλ);
        final double sinα  = sin(α);
        final double cosα  = cos(α);
        final double H     = cosα * Hp;

        /*
         * Equations are:    s  =  asin(cos(φ₀)⋅sin(Ψ) − sin(φ₀)⋅cos(Ψ)) ⋅ signum(cos(α))     for small α
         *                   s  =  asin(sin(λ)⋅cos(Ψ) / sin(α))                               for other α
         *
         * Using identity:   sin(atan(x))  =  x / √(1 + x²)
         * Rewrite as:       sin(Ψ)  =   tan(Ψ) / √(1 + tan²Ψ)
         */
        double c;
        if (abs(sinα) < ANGULAR_TOLERANCE) {
            c = (cosφ0*tanΨ - sinφ0) / rcosΨ;
            if (cosα < 0) c = -c;
        } else {
            c = sinλ / (rcosΨ * sinα);
        }
        c = asin(c);                    // After this line this is the `s` value in EPSG guidance notes.
        final double s = c;
        final double s2 = c  * c;
        final double s3 = s2 * c;
        final double s4 = s2 * s2;
        final double s5 = s4 * c;
        final double H2 = H*H;
        final double GH = G*H;
        c *= 1 - (s2/6   *  H2*(1 -   H2))
               + (s3/8   *  GH*(1 - 2*H2))
               + (s4/120 * (H2*(4 - 7*H2) - 3*(G*G)*(1 - 7*H2)))
               - (s5/48  * GH);

        if (dstPts != null) {
            dstPts[dstOff  ] = c * sinα;
            dstPts[dstOff+1] = c * cosα;
        }
        if (!derivate) {
            return null;
        }

        /*
         * End of map projection. Now compute the derivative, if requested.
         */
        final double Ψ      = atan(tanΨ);
        final double cosΨ   = cos(Ψ);
        final double sinΨ   = sin(Ψ);
        final double secΨ   = 1 / cosΨ;
        final double secΨ2  = secΨ * secΨ;

        final double tanφ   = sinφ / cosφ;
        final double secφ   = 1 / cosφ;
        final double secφ2  = secφ * secφ;

        final double tanα   = tan(α);
        final double cotα   = 1 / tanα;
        final double cscα   = 1 / sinα;

        final double dν_dφ      = eccentricitySquared * cosφ * sinφ / (rν * rν * rν);
        final double dΨ_dφ      = ((1  - eccentricitySquared) * secφ2  +  ℯ2_ν0_sinφ0 * secφ * tanφ * rν  -  ℯ2_ν0_sinφ0 * secφ * dν_dφ * rν*rν)
                / (1 + tanΨ * tanΨ);
        final double utilvar1   = cosφ0 * tanΨ  -  sinφ0 * cosλ;
        final double dα_dλ      = (cosλ * utilvar1  -  sinφ0 * sinλ2) / (sinλ2  +  utilvar1 * utilvar1);
        final double dα_dφ      = - cosφ0 * secΨ2 * sinλ * dΨ_dφ / (sinλ2  +  utilvar1 * utilvar1);

        double ds_dφ;
        double ds_dλ;
        double utilvar2;
        if (abs(sinα) < ANGULAR_TOLERANCE) {
            utilvar2 = cosΨ * sinφ0  -  cosφ0 * sinΨ;
            ds_dλ = 0;
            ds_dφ = (cosφ0 * cosΨ * dΨ_dφ  +  sinφ0 * sinΨ * dΨ_dφ) / sqrt(1  -  utilvar2 * utilvar2);

            if (cosα < 0) ds_dφ = -ds_dφ;
        } else {
            utilvar2 = cosΨ * cscα * sinλ;
            ds_dλ = (cosλ * cosΨ * cscα  -  utilvar2 * cotα * dα_dλ) / sqrt(1  -  utilvar2 * utilvar2);
            ds_dφ = - (cscα * sinλ * sinΨ * dΨ_dφ  +  utilvar2 * cotα * dα_dφ) / sqrt(1  -  utilvar2 * utilvar2);
        }

        final double H3 = H2 * H;
        final double G2 = G * G;

        final double dH_dλ = - eccentricity * cosφ0 * sinα * dα_dλ / sqrt(1 - eccentricitySquared);
        final double dH_dφ = - eccentricity * cosφ0 * sinα * dα_dφ / sqrt(1 - eccentricitySquared);

        final double h1 = H2*(1 - H2);
        final double h2 = GH*(1 - 2*H2);
        final double h3 = H2*(4 - 7*H2) - 3*G2*(1 - 7*H2);
        final double h4 = GH;

        final double dh1_dφ = -2*H3*dH_dφ + 2*H*(1 - H2)*dH_dφ;
        final double dh2_dφ = -4*G*H2*dH_dφ + G*(1 - 2*H2)*dH_dφ;
        final double dh3_dφ = 42*G2*H*dH_dφ - 14*H3*dH_dφ + 2*H*(4 - 7*H2)*dH_dφ;
        final double dh4_dφ = G*dH_dφ;

        final double dh1_dλ = -2*H3*dH_dλ + 2*H*(1 - H2)*dH_dλ;
        final double dh2_dλ = -4*G*H2*dH_dλ + G*(1 - 2*H2)*dH_dλ;
        final double dh3_dλ = 42*G2*H*dH_dλ - 14*H3*dH_dλ + 2*H*(4 - 7*H2)*dH_dλ;
        final double dh4_dλ = G*dH_dλ;

        final double dc_λ = ds_dλ * (
                1
                        - s2/6   * h1
                        + s3/8   * h2
                        + s4/120 * h3
                        - s5/48  * h4
        )
                + s * (
                        - s2/6   * dh1_dλ
                        + s3/8   * dh2_dλ
                        + s4/120 * dh3_dλ
                        - s5/48  * dh4_dλ

                        - s/3     * h1 * ds_dλ
                        + 3*s2/8  * h2 * ds_dλ
                        + s3/30   * h3 * ds_dλ
                        - 5*s4/48 * h4 * ds_dλ
        );

        final double dc_φ = ds_dφ * (
                1
                        - s2/6   * h1
                        + s3/8   * h2
                        + s4/120 * h3
                        - s5/48  * h4
        )
                + s * (
                        - s2/6   * dh1_dφ
                        + s3/8   * dh2_dφ
                        + s4/120 * dh3_dφ
                        - s5/48  * dh4_dφ

                        - s/3     * h1 * ds_dφ
                        + 3*s2/8  * h2 * ds_dφ
                        + s3/30   * h3 * ds_dφ
                        - 5*s4/48 * h4 * ds_dφ
        );

        final double dx_dλ = sinα * dc_λ  +  c * cosα * dα_dλ;
        final double dx_dφ = sinα * dc_φ  +  c * cosα * dα_dφ;
        final double dy_dλ = cosα * dc_λ  -  c * sinα * dα_dλ;
        final double dy_dφ = cosα * dc_φ  -  c * sinα * dα_dφ;

        return new Matrix2(dx_dλ, dx_dφ,      // ∂x/∂λ , ∂x/∂φ
                           dy_dλ, dy_dφ);     // ∂y/∂λ , ∂y/∂φ
    }

    /**
     * Converts the specified (<var>x</var>,<var>y</var>) coordinates
     * and stores the result in {@code dstPts} (angles in radians).
     */
    @Override
    protected void inverseTransform(final double[] srcPts, final int srcOff,
                                    final double[] dstPts, final int dstOff)
            throws ProjectionException
    {
        final double x    = srcPts[srcOff  ];
        final double y    = srcPts[srcOff+1];
        final double α    = atan2(x, y);                // Actually α′ in EPSG guidance notes.
        final double sinα = sin(α);
        final double cosα = cos(α);
              double negA = Hp * cosα; negA *= negA;    // mA = −A  compared to EPSG guidance note.
        final double B    = Bp * (1 + negA) * cosα;
        final double D2   = x*x + y*y;
        final double D    = sqrt(D2);                   // D = c′/ν₀, but division by ν₀ is already done here.
        final double J    = D + (negA*(1 -   negA)*(D2*D )/6)
                              - (   B*(1 - 3*negA)*(D2*D2)/24);
        final double J2   = J*J;
        final double K    = 1 + (negA*J2/2) - (B*(J2*J)/6);
        final double sinJ = sin(J);
        final double sinΨ = sinφ0*cos(J) + cosφ0*sinJ*cosα;
        final double cosΨ = sqrt(1 - sinΨ*sinΨ);
        dstPts[dstOff  ]  = asin(sinα*sinJ / cosΨ);
        dstPts[dstOff+1]  = atan((1 - eccentricitySquared*sinφ0*K / sinΨ) * (sinΨ/cosΨ)
                               / (1 - eccentricitySquared));
    }
}
