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
package org.apache.sis.geometry;

import java.text.Format;
import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.text.NumberFormat;
import java.text.DecimalFormat;
import java.text.FieldPosition;
import java.text.ParsePosition;
import java.text.ParseException;
import java.util.Arrays;
import java.util.Date;
import java.util.Locale;
import java.util.TimeZone;
import java.io.IOException;
import java.io.UncheckedIOException;
import javax.measure.Unit;
import javax.measure.UnitConverter;
import javax.measure.quantity.Time;
import javax.measure.IncommensurableException;
import org.opengis.geometry.DirectPosition;
import org.opengis.referencing.cs.AxisDirection;
import org.opengis.referencing.cs.CoordinateSystem;
import org.opengis.referencing.cs.CoordinateSystemAxis;
import org.opengis.referencing.crs.CoordinateReferenceSystem;
import org.opengis.referencing.crs.TemporalCRS;
import org.opengis.referencing.datum.Ellipsoid;
import org.apache.sis.internal.referencing.Formulas;
import org.apache.sis.internal.referencing.ReferencingUtilities;
import org.apache.sis.internal.system.Loggers;
import org.apache.sis.internal.util.LocalizedParseException;
import org.apache.sis.internal.util.Numerics;
import org.apache.sis.math.DecimalFunctions;
import org.apache.sis.util.logging.Logging;
import org.apache.sis.util.resources.Errors;
import org.apache.sis.util.ArgumentChecks;
import org.apache.sis.util.CharSequences;
import org.apache.sis.util.Characters;
import org.apache.sis.math.MathFunctions;
import org.apache.sis.measure.Angle;
import org.apache.sis.measure.AngleFormat;
import org.apache.sis.measure.Latitude;
import org.apache.sis.measure.Longitude;
import org.apache.sis.measure.Units;
import org.apache.sis.referencing.CRS;
import org.apache.sis.io.CompoundFormat;


/**
 * Formats spatiotemporal coordinates using number, angle and date formats inferred from the coordinate system.
 * The format for each coordinate is inferred from the
 * {@linkplain org.apache.sis.referencing.cs.DefaultCoordinateSystemAxis#getUnit() coordinate system units}
 * using the following rules:
 *
 * <ul>
 *   <li>Coordinate values in angular units are formatted as angles using {@link AngleFormat}.</li>
 *   <li>Coordinate values in temporal units are formatted as dates using {@link DateFormat}.</li>
 *   <li>Other values are formatted as numbers using {@link NumberFormat} followed by the unit symbol
 *       formatted by {@link org.apache.sis.measure.UnitFormat}.</li>
 * </ul>
 *
 * The format can be controlled by invoking the {@link #applyPattern(Class, String)} public method,
 * or by overriding the {@link #createFormat(Class)} protected method.
 *
 * <h2>Coordinate reference system</h2>
 * If the Coordinate Reference System (CRS) of the positions to format is known, it should
 * {@linkplain #setDefaultCRS(CoordinateReferenceSystem) be specified} before to invoke other methods in this
 * class because that CRS may impact calculation of properties like the number of fraction digits achieving a
 * {@linkplain #setPrecision(double, Unit) specified precision}.
 * That default CRS is used only if the {@link DirectPosition} does not specify its own CRS.
 *
 * <p>{@code CoordinateFormat} does <strong>not</strong> transform the given coordinates in a unique CRS.
 * If the coordinates need to be formatted in a specific CRS, then the caller should
 * {@linkplain org.apache.sis.referencing.operation.transform.AbstractMathTransform#transform(DirectPosition, DirectPosition)
 * transform the position} before to format it.</p>
 *
 * @author  Martin Desruisseaux (MPO, IRD, Geomatys)
 * @version 1.1
 *
 * @see AngleFormat
 * @see org.apache.sis.measure.UnitFormat
 * @see GeneralDirectPosition
 *
 * @since 0.8
 * @module
 */
public class CoordinateFormat extends CompoundFormat<DirectPosition> {
    /**
     * Serial number for cross-version compatibility.
     */
    private static final long serialVersionUID = 8324486673169133932L;

    /**
     * Maximal number of characters to convert to {@link String} if the text to parse is not a string instance.
     * This is an arbitrary limit that may change (or be removed) in any future SIS version.
     */
    private static final int READ_AHEAD_LIMIT = 256;

    /**
     * Maximal number of dimensions to use when parsing a coordinate without {@link #defaultCRS}.
     * This is an arbitrary limit that may change (or be removed) in any future SIS version.
     * To avoid this limitation, users are encouraged to specify a default CRS.
     */
    private static final int DEFAULT_DIMENSION = 4;

    /**
     * The separator between each coordinate values to be formatted.
     * The default value is a space.
     */
    private String separator;

    /**
     * The separator without spaces, used at parsing time.
     */
    private String parseSeparator;

    /**
     * The desired precisions for each coordinate, or {@code null} if unspecified.
     * The length of this array does not need to be equal to the number of dimensions;
     * extraneous values are ignored and missing values are assumed equal to 0.
     * A value of 0 means to use the default precision for that dimension.
     *
     * <p>Note that this is the precision specified by the user, which may differ from
     * the precision returned by {@link #getPrecisions()}.</p>
     *
     * @see #setPrecisions(double...)
     * @see #getPrecisions()
     */
    private double[] desiredPrecisions;

    /**
     * The coordinate reference system to assume if no CRS is attached to the position to format.
     * May be {@code null}.
     */
    private CoordinateReferenceSystem defaultCRS;

    /**
     * The coordinate reference system of the last {@link DirectPosition} that we formatted.
     * This is used for determining if we need to recompute all other transient fields in this class.
     */
    private transient CoordinateReferenceSystem lastCRS;

    /**
     * Constants for the {@link #types} array.
     */
    private static final byte LONGITUDE=1, LATITUDE=2, ANGLE=3, DATE=4, TIME=5;

    /**
     * The type for each value in the {@code formats} array, or {@code null} if not yet computed.
     * Types are: 0=number, 1=longitude, 2=latitude, 3=other angle, 4=date, 5=elapsed time.
     */
    private transient byte[] types;

    /**
     * The format instances given by {@link #getFormat(Class)}, to use by default when we have
     * not been able to configure the precision. This is the same array than {@link #formats},
     * unless {@link #setPrecisions(double...)} has been invoked.
     * Values at different indices may reference the same {@link Format} instance.
     */
    private transient Format[] sharedFormats;

    /**
     * The formats to use for formatting each coordinate value, or {@code null} if not yet computed.
     * The length of this array should be equal to the number of dimensions in {@link #lastCRS}.
     * Values at different indices may reference the same {@link Format} instance.
     */
    private transient Format[] formats;

    /**
     * The units for each dimension to be formatted as number.
     * We do not store this information for dimensions to be formatted as angle or date.
     */
    private transient Unit<?>[] units;

    /**
     * Conversions from arbitrary units to the unit used by formatter, or {@code null} if none.
     * For example in the case of dates, this is the conversions from temporal axis units to milliseconds.
     */
    private transient UnitConverter[] toFormatUnit;

    /**
     * Units symbols. Used only for coordinate to be formatted as ordinary numbers.
     * Non-null only if at least one coordinate is to be formatted that way.
     */
    private transient String[] unitSymbols;

    /**
     * Flags the coordinate values that need to be inverted before to be formatted.
     * This is needed for example if the axis is oriented toward past instead than future,
     * or toward west instead than east.
     *
     * @see #negate(int)
     */
    private transient long negate;

    /**
     * The time epochs. Non-null only if the at least on coordinate is to be formatted as a date.
     */
    private transient long[] epochs;

    /**
     * Dummy field position.
     */
    private transient FieldPosition dummy;

    /**
     * Temporary buffer to use if the {@code toAppendTo} argument given to {@link #format(DirectPosition, Appendable)}
     * is not an instance of {@code StringBuffer}.
     */
    private transient StringBuffer buffer;

    /**
     * Constructs a new coordinate format with default locale and timezone.
     */
    public CoordinateFormat() {
        this(Locale.getDefault(Locale.Category.FORMAT), TimeZone.getDefault());
    }

    /**
     * Constructs a new coordinate format for the specified locale and timezone.
     *
     * @param  locale    the locale for the new {@code Format}, or {@code null} for {@code Locale.ROOT}.
     * @param  timezone  the timezone, or {@code null} for UTC.
     */
    public CoordinateFormat(final Locale locale, final TimeZone timezone) {
        super(locale, timezone);
        parseSeparator = separator = " ";
    }

    /**
     * Returns the separator between each coordinate (number, angle or date).
     * The default value is a single space.
     *
     * @return the current coordinate separator.
     */
    public String getSeparator() {
        return separator;
    }

    /**
     * Sets the separator between each coordinate.
     * The default value is a single space.
     *
     * @param  separator  the new coordinate separator.
     */
    public void setSeparator(final String separator) {
        ArgumentChecks.ensureNonEmpty("separator", separator);
        this.separator = separator;
        parseSeparator = CharSequences.trimWhitespaces(separator);
        if (parseSeparator.isEmpty()) {
            parseSeparator = separator;
        }
    }

    /**
     * Returns the coordinate reference system to use if no CRS is explicitly associated to a given {@code DirectPosition}.
     * This CRS determines the type of format to use for each coordinate (number, angle or date) and the number of fraction
     * digits to use for achieving a {@linkplain #setPrecision(double, Unit) specified precision}.
     *
     * @return the default coordinate reference system, or {@code null} if none.
     */
    public CoordinateReferenceSystem getDefaultCRS() {
        return defaultCRS;
    }

    /**
     * Sets the coordinate reference system to use if no CRS is explicitly associated to a given {@code DirectPosition}.
     * This CRS is only a default; positions given in another CRS are <strong>not</strong> automatically transformed to
     * that CRS before formatting.
     *
     * @param  crs  the default coordinate reference system, or {@code null} if none.
     */
    public void setDefaultCRS(final CoordinateReferenceSystem crs) {
        defaultCRS = crs;
    }

    /**
     * Computes the value of transient fields from the given CRS.
     */
    private void configure(final CoordinateReferenceSystem crs) {
        types         = null;
        formats       = null;
        sharedFormats = null;
        units         = null;
        toFormatUnit  = null;
        unitSymbols   = null;
        epochs        = null;
        negate        = 0;
        lastCRS       = crs;
        if (crs == null) {
            return;
        }
        /*
         * If no CRS were specified, we will format everything as numbers. Working with null CRS
         * is sometime useful because null CRS are allowed in DirectPosition according ISO 19107.
         * Otherwise (if a CRS is given), infer the format subclasses from the axes.
         */
        final CoordinateSystem cs = crs.getCoordinateSystem();
        if (cs == null) {
            return;                                    // Paranoiac check (should never be null).
        }
        final int dimension = cs.getDimension();
        final byte[]   types   = new byte  [dimension];
        final Format[] formats = new Format[dimension];
        for (int i=0; i<dimension; i++) {
            final CoordinateSystemAxis axis = cs.getAxis(i);
            if (axis == null) {                                               // Paranoiac check.
                formats[i] = getDefaultFormat();
                continue;
            }
            final Unit<?> unit = axis.getUnit();
            /*
             * Formatter for angular units. Target unit is DEGREE_ANGLE.
             * Type is LONGITUDE, LATITUDE or ANGLE depending on axis direction.
             */
            if (Units.isAngular(unit)) {
                byte type = ANGLE;
                final AxisDirection dir = axis.getDirection();
                if      (AxisDirection.NORTH.equals(dir)) {type = LATITUDE;}
                else if (AxisDirection.EAST .equals(dir)) {type = LONGITUDE;}
                else if (AxisDirection.SOUTH.equals(dir)) {type = LATITUDE;  negate(i);}
                else if (AxisDirection.WEST .equals(dir)) {type = LONGITUDE; negate(i);}
                types  [i] = type;
                formats[i] = getFormat(Angle.class);
                setConverter(dimension, i, unit.asType(javax.measure.quantity.Angle.class).getConverterTo(Units.DEGREE));
                continue;
            }
            /*
             * Formatter for temporal units. Target unit is MILLISECONDS.
             * Type is DATE.
             */
            if (Units.isTemporal(unit)) {
                final CoordinateReferenceSystem t = CRS.getComponentAt(crs, i, i+1);
                if (t instanceof TemporalCRS) {
                    if (epochs == null) {
                        epochs = new long[dimension];
                    }
                    types  [i] = DATE;
                    formats[i] = getFormat(Date.class);
                    epochs [i] = ((TemporalCRS) t).getDatum().getOrigin().getTime();
                    setConverter(dimension, i, unit.asType(Time.class).getConverterTo(Units.MILLISECOND));
                    if (AxisDirection.PAST.equals(axis.getDirection())) {
                        negate(i);
                    }
                    continue;
                }
                types[i] = TIME;
                // Fallthrough: formatted as number.
            }
            /*
             * Formatter for all other units. Do NOT set types[i] since it may have been set
             * to a non-zero value by previous case. If not, the default value (zero) is the
             * one we want.
             */
            formats[i] = getFormat(Number.class);
            if (unit != null) {
                if (units == null) {
                    units = new Unit<?>[dimension];
                }
                units[i] = unit;
                final String symbol = getFormat(Unit.class).format(unit);
                if (!symbol.isEmpty()) {
                    if (unitSymbols == null) {
                        unitSymbols = new String[dimension];
                    }
                    unitSymbols[i] = symbol;
                }
            }
        }
        this.types    = types;         // Assign only on success.
        this.formats  = formats;
        sharedFormats = formats;
    }

    /**
     * Returns a clone of the format at the specified dimension. Format instances are cloned only when first needed.
     * The clones are needed when we want to change the format pattern (number of fraction digits, <i>etc.</i>) for
     * only one dimension, without impacting other dimensions that may use the same format.
     */
    private Format getFormatClone(final int dim) {
        if (formats == sharedFormats) {
            formats = formats.clone();
        }
        Format format = formats[dim];
        if (format == sharedFormats[dim]) {
            formats[dim] = format = (Format) format.clone();
        }
        return format;
    }

    /**
     * The default format to use when no CRS or no axis information is available. The coordinate type
     * could be anything (a date, an angle, …), but since we have no information we assume a number.
     * This method is defined for making clearer when such fallback is used.
     */
    private Format getDefaultFormat() {
        return getFormat(Number.class);
    }

    /**
     * Sets the unit converter at the given index.
     */
    private void setConverter(final int dimension, final int i, final UnitConverter c) {
        if (!c.isIdentity()) {
            if (toFormatUnit == null) {
                toFormatUnit = new UnitConverter[dimension];
            }
            toFormatUnit[i] = c;
        }
    }

    /**
     * Remembers that coordinate values at the given dimension will need to have their sign reverted.
     */
    private void negate(final int dimension) {
        if (dimension >= Long.SIZE) {
            throw new ArithmeticException(Errors.format(Errors.Keys.ExcessiveNumberOfDimensions_1, dimension + 1));
        }
        negate |= (1L << dimension);
    }

    /**
     * Returns {@code true} if the value at the given dimension needs to have its sign reversed.
     */
    private boolean isNegative(final int dimension) {
        return (negate & Numerics.bitmask(dimension)) != 0;
    }

    /**
     * Returns the precisions at which coordinate values are formatted in each dimension.
     * For example if coordinates in dimension <var>i</var> are formatted with two fraction digits,
     * then the precision reported in {@code precisions[i]} will be 0.01. If the precision can not
     * be determined for some dimensions, the corresponding values in the returned array will be 0.
     *
     * <p>The values returned by this method are not necessarily equal to the values specified in the last
     * call to {@link #setPrecisions(double...)}. For example if a precision of 0.03 has been requested for
     * a dimension whose coordinates are formatted as decimal numbers, then the actual precision returned
     * by this method for that dimension will be 0.01.</p>
     *
     * @return precision of coordinate values in each dimension (may contain 0 values for unknown precisions).
     *
     * @see AngleFormat#getPrecision()
     * @see DecimalFormat#getMaximumFractionDigits()
     *
     * @since 1.1
     */
    public double[] getPrecisions() {
        if (lastCRS != defaultCRS) {
            configure(defaultCRS);
        }
        Format[] cf = formats;
        if (cf == null) {
            cf = new Format[DEFAULT_DIMENSION];
            Arrays.fill(cf, getDefaultFormat());
        }
        final double[] precisions = new double[cf.length];
        for (int i=0; i<precisions.length; i++) {
            final Format f = cf[i];
            if (f instanceof DecimalFormat) {
                /*
                 * Intentionally check the DecimalFormat subtype, not the more generic NumberFormat type,
                 * because the calculation below assumes base 10 and assumes that fraction digits are for
                 * fractions of 1 (by contrast, CompactNumberFormat may apply fraction to larger values).
                 */
                precisions[i] = MathFunctions.pow10(-((DecimalFormat) f).getMaximumFractionDigits());
            } else if (f instanceof AngleFormat) {
                precisions[i] = ((AngleFormat) f).getPrecision();
            }
        }
        return precisions;
    }

    /**
     * Sets the desired precisions at which to format coordinate values in each dimension.
     * For example if {@code precisions[i]} is 0.05, then coordinates in dimension <var>i</var>
     * will be shown with two fraction digits when formatted as decimal numbers, or with "D°MM"
     * pattern when formatted as angles.
     *
     * <p>This precision does not have a clear relationship to the precision on the ground.
     * For example a precision of 0.01 could be one centimeter or 10 meters, depending if
     * the units of measurement in that dimension is meter or kilometer.
     * For a precision related to the ground, use {@link #setPrecision(double, Unit)} instead.</p>
     *
     * <p>If any value in the given array is 0 or {@link Double#NaN}, then there is a choice:
     * if {@link #setPrecision(double, Unit)} has been invoked, the ground precision specified to that
     * method will apply (if possible). Otherwise an implementation-specific default precision is used.
     * So it is possible to use {@link #setPrecision(double, Unit)} for specifying a precision in "real
     * world" units and to use this {@code setPrecisions(double...)} method for adjusting the precision
     * of only the vertical axis for example.</p>
     *
     * @param  precisions  desired precision at which to format coordinate values in each dimension
     *                     (may have 0 values for unspecified precision in some of those dimensions),
     *                     or {@code null} for restoring the default values.
     *
     * @see AngleFormat#setPrecision(double, boolean)
     * @see DecimalFormat#setMaximumFractionDigits(int)
     *
     * @since 1.1
     */
    public void setPrecisions(final double... precisions) {
        if (precisions == null) {
            desiredPrecisions = null;
            formats = sharedFormats;
        } else {
            if (desiredPrecisions == null || desiredPrecisions.length != 0) {
                desiredPrecisions = new double[precisions.length];
                Arrays.fill(desiredPrecisions, Double.NaN);
            }
            for (int i=0; i<precisions.length; i++) {
                double p = Math.abs(precisions[i]);
                if (!(p < Double.POSITIVE_INFINITY)) p = 0;                 // Use ! for replacing NaN.
                if (desiredPrecisions[i] != (desiredPrecisions[i] = p)) {
                    // Precision changed. Keep format up to date.
                    if (formats != null && i < formats.length) {
                        applyPrecision(i);
                    }
                }
            }
        }
    }

    /**
     * Sets the pattern of the format for the specified dimension according the desired precision.
     * The format to configure is {@code formats[dim]} and the pattern will be constructed from the
     * {@code desiredPrecisions[dim]} value. Caller must ensure that the given dimension is valid
     * for both {@link #formats} and {@link #desiredPrecisions} arrays.
     */
    private void applyPrecision(final int dim) {
        final double precision = desiredPrecisions[dim];
        if (precision > 0) {
            final Format format = formats[dim];
            /*
             * Intentionally check the DecimalFormat subtype, not the more generic NumberFormat type,
             * because the calculation below assumes base 10 and assumes that fraction digits are for
             * fractions of 1 (by contrast, CompactNumberFormat may apply fraction to larger values).
             */
            if (format instanceof DecimalFormat) {
                final int c = Math.max(DecimalFunctions.fractionDigitsForDelta(precision, false), 0);
                final DecimalFormat nf = (DecimalFormat) getFormatClone(dim);
                nf.setMinimumFractionDigits(c);
                nf.setMaximumFractionDigits(c);
            } else if (format instanceof AngleFormat) {
                ((AngleFormat) getFormatClone(dim)).setPrecision(precision, true);
            }
        }
    }

    /**
     * Adjusts the number of fraction digits to show in coordinates for achieving the given precision.
     * The {@link NumberFormat} and {@link AngleFormat} are configured for coordinates expressed in the
     * {@linkplain #getDefaultCRS() default coordinate reference system} defined at the moment this method is invoked.
     * The number of fraction digits is <em>not</em> updated if a different CRS is specified after this method call
     * or if the coordinates to format are associated to a different CRS.
     *
     * <p>The given resolution will be converted to the units used by coordinate system axes. For example if a 10 metres
     * resolution is specified but the {@linkplain #getDefaultCRS() default CRS} axes use kilometres, then this method
     * converts the resolution to 0.01 kilometre and uses that value for inferring that coordinates should be formatted
     * with 2 fraction digits. If the resolution is specified in an angular units such as degrees, this method uses the
     * {@linkplain org.apache.sis.referencing.datum.DefaultEllipsoid#getAuthalicRadius() ellipsoid authalic radius} for
     * computing an equivalent resolution in linear units. For example if the ellipsoid of default CRS is WGS84,
     * then this method considers a resolution of 1 second of angle as equivalent to a resolution of about 31 meters.
     * Conversions work also in the opposite direction (from linear to angular units) and are also used for choosing
     * which angle fields (degrees, minutes or seconds) to show.</p>
     *
     * @param  resolution  the desired resolution.
     * @param  unit        unit of the desired resolution.
     *
     * @see DecimalFormat#setMaximumFractionDigits(int)
     * @see AngleFormat#setPrecision(double, boolean)
     *
     * @since 1.0
     */
    @SuppressWarnings("null")
    public void setPrecision(double resolution, Unit<?> unit) {
        ArgumentChecks.ensureFinite("resolution", resolution);
        ArgumentChecks.ensureNonNull("unit", unit);
        resolution = Math.abs(resolution);
        if (Units.isTemporal(unit)) {
            return;                                 // Setting temporal resolution is not yet implemented.
        }
        /*
         * If the given resolution is linear (for example in metres), compute an equivalent resolution in degrees
         * assuming a sphere of radius computed from the CRS.  Conversely if the resolution is angular (typically
         * in degrees), computes an equivalent linear resolution. For all other kind of units, do nothing.
         */
        Resolution specified = new Resolution(resolution, unit, Units.isAngular(unit));
        Resolution related   = null;
        IncommensurableException error = null;
        if (specified.isAngular || Units.isLinear(unit)) try {
            related = specified.related(ReferencingUtilities.getEllipsoid(defaultCRS));
        } catch (IncommensurableException e) {
            error = e;
        }
        /*
         * We now have the requested resolution in both linear and angular units, if equivalence has been established.
         * Convert those resolutions to the units actually used by the CRS. If some axes use different units, keep the
         * units which result in the finest resolution.
         */
        boolean relatedUsed = false;
        if (defaultCRS != null) {
            final CoordinateSystem cs = defaultCRS.getCoordinateSystem();
            if (cs != null) {                                                   // Paranoiac check (should never be null).
                final int dimension = cs.getDimension();
                for (int i=0; i<dimension; i++) {
                    final CoordinateSystemAxis axis = cs.getAxis(i);
                    if (axis != null) {                                         // Paranoiac check.
                        final Unit<?> axisUnit = axis.getUnit();
                        if (axisUnit != null) try {
                            final double maxValue = Math.max(Math.abs(axis.getMinimumValue()),
                                                             Math.abs(axis.getMaximumValue()));
                            if (!specified.forAxis(maxValue, axisUnit) && related != null) {
                                relatedUsed |= related.forAxis(maxValue, axisUnit);
                            }
                        } catch (IncommensurableException e) {
                            if (error == null) error = e;
                            else error.addSuppressed(e);
                        }
                    }
                }
            }
        }
        if (error != null) {
            Logging.unexpectedException(Logging.getLogger(Loggers.MEASURE), CoordinateFormat.class, "setPrecision", error);
        }
        specified.setPrecision(this);
        if (relatedUsed) {
            related.setPrecision(this);
        }
    }

    /**
     * Desired resolution in a given units, together with methods for converting to the units of a coordinate system axis.
     * This is a helper class for {@link CoordinateFormat#setPrecision(double, Unit)} implementation. An execution of that
     * method typically creates two instances of this {@code Resolution} class: one for the resolution in metres and another
     * one for the resolution in degrees.
     */
    private static final class Resolution {
        /** The desired resolution in the unit of measurement given by {@link #unit}. */
        private double resolution;

        /** Maximal absolute value that we may format in unit of measurement given by {@link #unit}. */
        private double magnitude;

        /** Unit of measurement of {@link #resolution} or {@link #magnitude}. */
        private Unit<?> unit;

        /** Whether {@link #unit} is an angular unit. */
        final boolean isAngular;

        /** Creates a new instance initialized to the given resolution. */
        Resolution(final double resolution, final Unit<?> unit, final boolean isAngular) {
            this.resolution = resolution;
            this.unit       = unit;
            this.isAngular  = isAngular;
        }

        /**
         * If this resolution is in metres, returns equivalent resolution in degrees. Or conversely if this resolution
         * is in degrees, returns an equivalent resolution in metres. Other linear and angular units are accepted too;
         * they will be converted as needed.
         *
         * @param  ellipsoid  the ellipsoid, or {@code null} if none.
         * @return the related resolution, or {@code null} if none.
         */
        Resolution related(final Ellipsoid ellipsoid) throws IncommensurableException {
            final double radius = Formulas.getAuthalicRadius(ellipsoid);
            if (radius > 0) {                                       // Indirectly filter null ellipsoid.
                Unit<?> relatedUnit = ellipsoid.getAxisUnit();      // Angular if `unit` is linear, or linear if `unit` is angular.
                if (relatedUnit != null) {                          // Paranoiac check (should never be null).
                    double related;
                    if (isAngular) {
                        // Linear resolution  =  angular resolution in radians  ×  radius.
                        related = unit.getConverterToAny(Units.RADIAN).convert(resolution) * radius;
                    } else {
                        // Angular resolution in radians  =  linear resolution  /  radius
                        related = Math.toDegrees(unit.getConverterToAny(relatedUnit).convert(resolution) / radius);
                        relatedUnit = Units.DEGREE;
                    }
                    return new Resolution(related, relatedUnit, !isAngular);
                }
            }
            return null;
        }

        /**
         * Adjusts the resolution units for the given coordinate system axis. This methods select the units which
         * result in the smallest absolute value of {@link #resolution}.
         *
         * @param  maxValue  the maximal absolute value that a coordinate on the axis may have.
         * @param  axisUnit  {@link CoordinateSystemAxis#getUnit()}.
         * @return whether the given axis unit is compatible with the expected unit.
         */
        boolean forAxis(double maxValue, final Unit<?> axisUnit) throws IncommensurableException {
            if (!axisUnit.isCompatible(unit)) {
                return false;
            }
            final UnitConverter c = unit.getConverterToAny(axisUnit);
            final double r = Math.abs(c.convert(resolution));
            if (r < resolution) {
                resolution = r;                                         // To units producing the smallest value.
                unit = axisUnit;
            } else {
                maxValue = Math.abs(c.inverse().convert(maxValue));     // From axis units to selected units.
            }
            if (maxValue > magnitude) {
                magnitude = maxValue;
            }
            return true;
        }

        /**
         * Configures the {@link NumberFormat} or {@link AngleFormat} for a number of fraction digits
         * sufficient for the given resolution.
         */
        void setPrecision(final CoordinateFormat owner) {
            final Format format = owner.getFormat(isAngular ? Angle.class : Number.class);
            if (format instanceof DecimalFormat) {
                /*
                 * Intentionally check the DecimalFormat subtype, not the more generic NumberFormat type,
                 * because the calculation below assumes base 10 and assumes that fraction digits are for
                 * fractions of 1 (by contrast, CompactNumberFormat may apply fraction to larger values).
                 */
                if (resolution == 0) resolution = 1E-6;                     // Arbitrary value.
                final int p = Math.max(0, DecimalFunctions.fractionDigitsForDelta(resolution, true));
                final int m = Math.max(0, DecimalFunctions.fractionDigitsForDelta(Math.ulp(magnitude), false));
                ((DecimalFormat) format).setMinimumFractionDigits(Math.min(p, m));
                ((DecimalFormat) format).setMaximumFractionDigits(p);
            } else if (format instanceof AngleFormat) {
                ((AngleFormat) format).setPrecision(resolution, true);
            }
        }
    }

    /**
     * Returns the pattern for number, angle or date fields. The given {@code valueType} should be
     * {@code Number.class}, {@code Angle.class}, {@code Date.class} or a sub-type of the above.
     * This method may return {@code null} if the underlying format can not provide a pattern.
     *
     * <table class="sis">
     *   <caption>Pattern availability for type of value</caption>
     *   <tr><th>Value type</th>     <th>Base format class</th>    <th>Format with pattern</th></tr>
     *   <tr><td>{@link Number}</td> <td>{@link NumberFormat}</td> <td>{@link DecimalFormat}</td></tr>
     *   <tr><td>{@link Angle}</td>  <td>{@link AngleFormat}</td>  <td>{@link AngleFormat}</td></tr>
     *   <tr><td>{@link Date}</td>   <td>{@link DateFormat}</td>   <td>{@link SimpleDateFormat}</td></tr>
     * </table>
     *
     * @param  valueType  the base type of coordinate values to parse and format:
     *                    {@code Number.class}, {@code Angle.class} or {@code Date.class}.
     * @return the pattern for fields of the given type, or {@code null} if not applicable.
     *
     * @see #getFormat(Class)
     */
    public String getPattern(final Class<?> valueType) {
        final Format format = getFormat(valueType);
        if (format instanceof AngleFormat) {
            return ((AngleFormat) format).toPattern();
        } else if (format instanceof DecimalFormat) {
            return ((DecimalFormat) format).toPattern();
        } else if (format instanceof SimpleDateFormat) {
            return ((SimpleDateFormat) format).toPattern();
        } else {
            return null;
        }
    }

    /**
     * Sets the pattern for number, angle or date fields.
     * The pattern syntax depends on the {@code valueType} argument:
     *
     * <ul>
     *   <li>If {@code valueType} is {@code Number.class}, then the pattern syntax shall be as described in the
     *     {@link DecimalFormat} class. This pattern may be used for any coordinate to be formatted as plain number,
     *     for example in {@linkplain org.apache.sis.referencing.cs.DefaultCartesianCS Cartesian coordinate system}.</li>
     *   <li>If {@code valueType} is {@code Angle.class}, then the pattern syntax shall be as described in the
     *     {@link AngleFormat} class. This pattern may be used for any coordinate to be formatted as latitude or longitude,
     *     for example in {@linkplain org.apache.sis.referencing.cs.DefaultEllipsoidalCS ellipsoidal coordinate system}.</li>
     *   <li>If {@code valueType} is {@code Date.class}, then the pattern syntax shall be as described in the
     *     {@link SimpleDateFormat} class. This pattern may be used for any coordinate to be formatted as date and time,
     *     for example in {@linkplain org.apache.sis.referencing.cs.DefaultTimeCS time coordinate system}.</li>
     * </ul>
     *
     * @param  valueType  the base type of coordinate values to parse and format:
     *                    {@code Number.class}, {@code Angle.class} or {@code Date.class}.
     * @param  pattern    the pattern as specified in {@link DecimalFormat}, {@link AngleFormat}
     *                    or {@link SimpleDateFormat} javadoc.
     * @return {@code true} if the pattern has been applied, or {@code false} if {@code valueType} does not
     *         specify a known type or if the format associated to that type does not support patterns.
     * @throws IllegalArgumentException if the given pattern is invalid.
     */
    public boolean applyPattern(final Class<?> valueType, final String pattern) {
        ArgumentChecks.ensureNonNull("pattern", pattern);
        final Format format = getFormat(valueType);
        if (format instanceof DecimalFormat) {
            ((DecimalFormat) format).applyPattern(pattern);
        } else if (format instanceof SimpleDateFormat) {
            ((SimpleDateFormat) format).applyPattern(pattern);
        } else if (format instanceof AngleFormat) {
            ((AngleFormat) format).applyPattern(pattern);
        } else {
            return false;
        }
        formats = sharedFormats;            // For forcing an update of `formats` when needed.
        return true;
    }

    /**
     * Returns the base type of values parsed and formatted by this {@code Format} instance.
     *
     * @return {@code DirectPosition.class}.
     */
    @Override
    public final Class<DirectPosition> getValueType() {
        return DirectPosition.class;
    }

    /**
     * Formats the given coordinate.
     * The type of each coordinate value (number, angle or date) is determined by the CRS of the given
     * position if such CRS is defined, or from the {@linkplain #getDefaultCRS() default CRS} otherwise.
     *
     * @param  position  the coordinate to format.
     * @return the formatted position.
     */
    public String format(final DirectPosition position) {
        if (buffer == null) {
            buffer = new StringBuffer();
        }
        buffer.setLength(0);
        try {
            format(position, buffer);
        } catch (IOException e) {
            /*
             * Should never happen when writing into a StringBuffer, unless the user override the
             * format(…) method. We do not rethrow an AssertionError because of this possibility.
             */
            throw new UncheckedIOException(e);
        }
        return buffer.toString();
    }

    /**
     * Formats the given coordinate and appends the resulting text to the given stream or buffer.
     * The type of each coordinate value (number, angle or date) is determined by the CRS of the given
     * position if such CRS is defined, or from the {@linkplain #getDefaultCRS() default CRS} otherwise.
     *
     * @param  position    the coordinate to format.
     * @param  toAppendTo  where the text is to be appended.
     * @throws IOException if an error occurred while writing to the given appendable.
     */
    @Override
    @SuppressWarnings({"UnnecessaryBoxing", "null"})
    public void format(final DirectPosition position, final Appendable toAppendTo) throws IOException {
        ArgumentChecks.ensureNonNull("position",   position);
        ArgumentChecks.ensureNonNull("toAppendTo", toAppendTo);
        CoordinateReferenceSystem crs = position.getCoordinateReferenceSystem();
        if (crs == null) {
            crs = defaultCRS;                           // May still be null.
        }
        if (crs != lastCRS) {
            configure(crs);
        }
        /*
         * Unconditionally configure the formatters for the desired precisions because those precisions
         * may change for every point. Note that the formatters may not have been created if the CRS is
         * null (because `configure(…)` does not know which format to use), in which case generic number
         * formats will be used.
         */
        if (desiredPrecisions != null) {
            if (sharedFormats == null) {
                formats = sharedFormats = new Format[desiredPrecisions.length];
                Arrays.fill(formats, getDefaultFormat());
                types = new byte[formats.length];
            }
            final int n = Math.min(desiredPrecisions.length, formats.length);
            for (int i=0; i<n; i++) {
                applyPrecision(i);
            }
        }
        /*
         * Standard java.text.Format API can only write into a StringBuffer. If the given Appendable is not a
         * StringBuffer, then we will need to format in a temporary buffer before to copy to the Appendable.
         */
        final StringBuffer destination;
        if (toAppendTo instanceof StringBuffer) {
            destination = (StringBuffer) toAppendTo;
        } else {
            if (buffer == null) {
                buffer = new StringBuffer();
            }
            destination = buffer;
            destination.setLength(0);
        }
        if (dummy == null) {
            dummy = new FieldPosition(0);
        }
        /*
         * The format to use for each coordinate has been computed by `configure`. The format array length
         * should match the number of dimensions in the given position if the DirectPosition is consistent
         * with its CRS, but we will nevertheless verify has a paranoiac check.  If there is no CRS, or if
         * the DirectPosition dimension is (illegally) greater than the CRS dimension, then we will format
         * the coordinate as a number.
         */
        final int dimension = position.getDimension();
        for (int i=0; i < dimension; i++) {
            double value = position.getOrdinate(i);
            final Object object;
            final Format f;
            if (formats != null && i < formats.length) {
                f = formats[i];
                if (isNegative(i)) {
                    value = -value;
                }
                if (toFormatUnit != null) {
                    final UnitConverter c = toFormatUnit[i];
                    if (c != null) {
                        value = c.convert(value);
                    }
                }
                switch (types[i]) {
                    default:        object = Double.valueOf(value); break;
                    case LONGITUDE: object = new Longitude (value); break;
                    case LATITUDE:  object = new Latitude  (value); break;
                    case ANGLE:     object = new Angle     (value); break;
                    case DATE:      object = new Date(Math.round(value) + epochs[i]); break;
                }
            } else {
                object = value;
                f = getDefaultFormat();
            }
            /*
             * At this point we got the value to format together with the Format instance to use.
             */
            if (i != 0) {
                toAppendTo.append(separator);
            }
            if (f.format(object, destination, dummy) != toAppendTo) {
                toAppendTo.append(destination);
                destination.setLength(0);
            }
            if (unitSymbols != null && i < unitSymbols.length) {
                final String symbol = unitSymbols[i];
                if (symbol != null) {
                    toAppendTo.append(Characters.NO_BREAK_SPACE).append(symbol);
                }
            }
        }
    }

    /**
     * Parses a coordinate from the given character sequence.
     * This method presumes that the coordinate reference system is the {@linkplain #getDefaultCRS() default CRS}.
     * The parsing begins at the {@linkplain ParsePosition#getIndex() index} given by the {@code pos} argument.
     * If parsing succeeds, then the {@code pos} index is updated to the index after the last coordinate value and
     * the parsed coordinate is returned. Otherwise (if parsing fails), the {@code pos} index is left unchanged,
     * the {@code pos} {@linkplain ParsePosition#getErrorIndex() error index} is set to the index of the first
     * unparsable character and an exception is thrown with a similar {@linkplain ParseException#getErrorOffset()
     * error index}.
     *
     * @param  text  the character sequence for the coordinate to parse.
     * @param  pos   the index where to start the parsing.
     * @return the parsed coordinate (never {@code null}).
     * @throws ParseException if an error occurred while parsing the coordinate.
     */
    @Override
    public DirectPosition parse(final CharSequence text, final ParsePosition pos) throws ParseException {
        ArgumentChecks.ensureNonNull("text", text);
        ArgumentChecks.ensureNonNull("pos",  pos);
        final int start  = pos.getIndex();
        final int length = text.length();
        /*
         * The NumberFormat, DateFormat and AngleFormat work only on String values, not on CharSequence.
         * If the given text is not a String, we will convert an arbitrarily small section of the given
         * text. Note that this will require to adjust the ParsePosition indices.
         */
        final int offset;
        final String asString;
        final ParsePosition subPos;
        if (text instanceof String) {
            offset   = 0;
            subPos   = pos;
            asString = (String) text;
        } else {
            offset   = start;
            subPos   = new ParsePosition(0);
            asString = text.subSequence(start, Math.min(start + READ_AHEAD_LIMIT, length)).toString();
        }
        /*
         * The Format instances to be used for each coordinate values is determined by the default CRS.
         * If no such CRS has been specified, then we will parse everything as plain numbers.
         */
        if (lastCRS != defaultCRS) {
            configure(defaultCRS);
        }
        final double[] coordinates;
        Format format;
        final Format[] formats = this.formats;
        if (formats != null) {
            format      = null;
            coordinates = new double[formats.length];
        } else {
            format      = getDefaultFormat();
            coordinates = new double[DEFAULT_DIMENSION];
        }
        /*
         * For each coordinate value except the first one, we need to skip the separator.
         * If we do not find the separator, we may consider that we reached the coordinate
         * end ahead of time. We currently allow that only for coordinate without CRS.
         */
        for (int i=0; i < coordinates.length; i++) {
            if (i != 0) {
                final int end = subPos.getIndex();
                int index = offset + end;
                while (!CharSequences.regionMatches(text, index, parseSeparator)) {
                    if (index < length) {
                        final int c = Character.codePointAt(text, index);
                        if (Character.isSpaceChar(c)) {
                            index += Character.charCount(c);
                            continue;
                        }
                    }
                    if (formats == null) {
                        pos.setIndex(index);
                        return new GeneralDirectPosition(Arrays.copyOf(coordinates, i));
                    }
                    pos.setIndex(start);
                    pos.setErrorIndex(index);
                    throw new LocalizedParseException(getLocale(), Errors.Keys.UnexpectedCharactersAfter_2,
                            new CharSequence[] {text.subSequence(start, end), CharSequences.token(text, index)}, index);
                }
                subPos.setIndex(index + parseSeparator.length() - offset);
            }
            /*
             * At this point 'subPos' is set to the beginning of the next coordinate to parse in 'asString'.
             * Parse the value as a number, angle or date, as determined from the coordinate system axis.
             */
            if (formats != null) {
                format = formats[i];
            }
            @SuppressWarnings("null")
            final Object object = format.parseObject(asString, subPos);
            if (object == null) {
                /*
                 * If we failed to parse, build an error message with the type that was expected for that coordinate.
                 * If the given CharSequence was not a String, we may need to update the error index since we tried
                 * to parse only a substring.
                 */
                Class<?> type = Number.class;
                if (types != null) {
                    switch (types[i]) {
                        case LONGITUDE: type = Longitude.class; break;
                        case LATITUDE:  type = Latitude.class;  break;
                        case ANGLE:     type = Angle.class;     break;
                        case DATE:      type = Date.class;      break;
                    }
                }
                pos.setIndex(start);
                if (subPos != pos) {
                    pos.setErrorIndex(offset + subPos.getErrorIndex());
                }
                throw new LocalizedParseException(getLocale(), type, text, pos);
            }
            double value;
            if (object instanceof Angle) {
                value = ((Angle) object).degrees();
            } else if (object instanceof Date) {
                value = ((Date) object).getTime() - epochs[i];
            } else {
                value = ((Number) object).doubleValue();
            }
            /*
             * The conversions and sign reversal applied below shall be in exact reverse order than
             * in the 'format(…)' method. However we have one additional step compared to format(…):
             * the unit written after the coordinate value may not be the same than the unit declared
             * in the CRS axis, so we have to parse the unit and convert the value before to apply
             * the reverse of 'format(…)' steps.
             */
            if (units != null) {
                final Unit<?> target = units[i];
                if (target != null) {
                    final int base = subPos.getIndex();
                    int index = base;
                    /*
                     * Skip whitespaces using Character.isSpaceChar(…), not Character.isWhitespace(…),
                     * because we need to skip also the non-breaking space (Characters.NO_BREAK_SPACE).
                     * If we can not parse the unit after those spaces, we will revert to the original
                     * position (absence of unit will not be considered an error).
                     */
                    while (index < asString.length()) {
                        final int c = asString.codePointAt(index);
                        if (Character.isSpaceChar(c)) {
                            index += Character.charCount(c);
                            continue;
                        }
                        subPos.setIndex(index);
                        final Object unit = getFormat(Unit.class).parseObject(asString, subPos);
                        if (unit == null) {
                            subPos.setIndex(base);
                            subPos.setErrorIndex(-1);
                        } else try {
                            value = ((Unit<?>) unit).getConverterToAny(target).convert(value);
                        } catch (IncommensurableException e) {
                            index += offset;
                            pos.setIndex(start);
                            pos.setErrorIndex(index);
                            throw (ParseException) new ParseException(e.getMessage(), index).initCause(e);
                        }
                        break;
                    }
                }
            }
            if (toFormatUnit != null) {
                final UnitConverter c = toFormatUnit[i];
                if (c != null) {
                    value = c.inverse().convert(value);
                }
            }
            if (isNegative(i)) {
                value = -value;
            }
            coordinates[i] = value;
        }
        final GeneralDirectPosition position = new GeneralDirectPosition(coordinates);
        position.setCoordinateReferenceSystem(defaultCRS);
        return position;
    }

    /**
     * Returns a clone of this format.
     *
     * @return a clone of this format.
     */
    @Override
    public CoordinateFormat clone() {
        final CoordinateFormat clone = (CoordinateFormat) super.clone();
        clone.dummy  = null;
        clone.buffer = null;
        clone.configure(null);
        if (desiredPrecisions != null) {
            clone.desiredPrecisions = desiredPrecisions.clone();
        }
        return clone;
    }
}
