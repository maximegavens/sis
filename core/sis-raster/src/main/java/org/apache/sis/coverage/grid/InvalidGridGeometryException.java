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
package org.apache.sis.coverage.grid;


/**
 * Thrown when the argument specified to a method or constructor would result in an invalid {@link GridGeometry}.
 *
 * @author  Martin Desruisseaux (Geomatys)
 * @version 1.0
 * @since   1.0
 * @module
 */
public class InvalidGridGeometryException extends IllegalArgumentException {
    /**
     * Serial number for inter-operability with different versions.
     */
    private static final long serialVersionUID = -7849140502096470380L;

    /**
     * Constructs an exception with no detail message.
     */
    public InvalidGridGeometryException() {
    }

    /**
     * Constructs an exception with the specified detail message.
     *
     * @param  message  the detail message.
     */
    public InvalidGridGeometryException(final String message) {
        super(message);
    }

    /**
     * Constructs an exception with the specified cause.
     *
     * @param  cause  the cause for this exception.
     */
    public InvalidGridGeometryException(final Throwable cause) {
        super(cause);
    }

    /**
     * Constructs an exception with the specified detail message and cause.
     *
     * @param  message  the detail message.
     * @param  cause    the cause for this exception.
     */
    public InvalidGridGeometryException(final String message, final Throwable cause) {
        super(message, cause);
    }
}
