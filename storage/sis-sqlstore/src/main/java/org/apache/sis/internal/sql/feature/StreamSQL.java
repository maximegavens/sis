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
package org.apache.sis.internal.sql.feature;

import java.sql.Connection;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.sql.Statement;
import java.util.concurrent.Callable;
import java.util.concurrent.atomic.AtomicReference;
import java.util.function.Consumer;
import java.util.function.DoubleConsumer;
import java.util.function.DoubleFunction;
import java.util.function.DoubleUnaryOperator;
import java.util.function.Function;
import java.util.function.Predicate;
import java.util.function.Supplier;
import java.util.function.ToDoubleFunction;
import java.util.function.ToIntFunction;
import java.util.function.ToLongFunction;
import java.util.stream.DoubleStream;
import java.util.stream.IntStream;
import java.util.stream.LongStream;
import java.util.stream.Stream;
import java.util.stream.StreamSupport;

import org.opengis.feature.Feature;

import org.apache.sis.internal.util.DoubleStreamDecoration;
import org.apache.sis.internal.util.StreamDecoration;
import org.apache.sis.storage.DataStoreException;
import org.apache.sis.util.collection.BackingStoreException;

import static org.apache.sis.util.ArgumentChecks.ensureNonNull;

/**
 * Manages query lifecycle and optimizations. Operations like {@link #count()}, {@link #distinct()}, {@link #skip(long)}
 * and {@link #limit(long)} are delegated to underlying SQL database. This class consistently propagate optimisation
 * strategies through streams obtained using {@link #map(Function)}, {@link #mapToDouble(ToDoubleFunction)} and
 * {@link #peek(Consumer)} operations. However, for result consistency, no optimization is stacked once either
 * {@link #filter(Predicate)} or {@link #flatMap(Function)} operations are called, as they modify browing flow (the
 * count of stream elements is not bound 1 to 1 to query result rows).
 *
 * @since 1.0
 *
 * @author Alexis Manin (Geomatys)
 *
 */
class StreamSQL extends StreamDecoration<Feature> {

    final Features.Builder queryBuilder;
    boolean parallel;

    private Consumer<? super Feature> peekAction;

    StreamSQL(final Table source) {
        this(new Features.Builder(source));
    }

    StreamSQL(Features.Builder builder) {
        this.queryBuilder = builder;
    }

    @Override
    public <R> Stream<R> map(Function<? super Feature, ? extends R> mapper) {
        return new MappedStream<>(mapper, this);
    }

    @Override
    public IntStream mapToInt(ToIntFunction<? super Feature> mapper) {
        return super.mapToInt(mapper);
    }

    @Override
    public LongStream mapToLong(ToLongFunction<? super Feature> mapper) {
        return super.mapToLong(mapper);
    }

    @Override
    public DoubleStream mapToDouble(ToDoubleFunction<? super Feature> mapper) {
        return super.mapToDouble(mapper);
    }

    @Override
    public Stream<Feature> parallel() {
        parallel = true;
        return this;
    }

    @Override
    public Stream<Feature> sequential() {
        parallel = false;
        return this;
    }

    @Override
    public Stream<Feature> distinct() {
        queryBuilder.distinct = true;
        return this;
    }

    @Override
    @SuppressWarnings("unchecked")
    public Stream<Feature> peek(Consumer<? super Feature> action) {
        if (peekAction == null) {
            peekAction = action;
        } else {
            // Safe cast, because Stream values are strongly typed to O.
            peekAction = peekAction.andThen((Consumer) action);
        }

        return this;
    }

    @Override
    public Stream<Feature> limit(long maxSize) {
        if (queryBuilder.limit < 1) queryBuilder.limit = maxSize;
        else queryBuilder.limit = Math.min(queryBuilder.limit, maxSize);
        return this;
    }

    @Override
    public Stream<Feature> skip(long n) {
        queryBuilder.offset += n;
        return this;
    }

    @Override
    public long count() {
        // Avoid opening a connection if sql text cannot be evaluated.
        final String sql;
        try {
            sql = queryBuilder.getSnapshot(true);
        } catch (SQLException e) {
            throw new BackingStoreException("Cannot create SQL COUNT query", e);
        }
        try (Connection conn = queryBuilder.parent.source.getConnection()) {
            try (Statement st = conn.createStatement();
                 ResultSet rs = st.executeQuery(sql)) {
                if (rs.next()) {
                    return rs.getLong(1);
                } else return 0;
            }
        } catch (SQLException e) {
            throw new BackingStoreException("Cannot estimate feature set size using SQL COUNT query", e);
        }
    }

    @Override
    protected synchronized Stream<Feature> createDecoratedStream() {
        final AtomicReference<Connection> connectionRef = new AtomicReference<>();
        Stream<Feature> featureStream = Stream.of(uncheck(() -> queryBuilder.parent.source.getConnection()))
                .map(Supplier::get)
                .peek(connectionRef::set)
                .flatMap(conn -> {
                    try {
                        final Features iter = queryBuilder.build(conn);
                        return StreamSupport.stream(iter, parallel).onClose(iter);
                    } catch (SQLException | DataStoreException e) {
                        throw new BackingStoreException(e);
                    }
                })
                .onClose(() -> queryBuilder.parent.closeRef(connectionRef));
        if (peekAction != null) featureStream = featureStream.peek(peekAction);
        return featureStream;
    }

    /**
     * Transform a callable into supplier by catching any potential verified exception and rethrowing it as a {@link BackingStoreException}.
     * @param generator The callable to use in a non-verified error context. Must not be null.
     * @param <T> The return type of input callable.
     * @return A supplier that delegates work to given callable, wrapping any verified exception in the process.
     */
    private static <T> Supplier<T> uncheck(final Callable<T> generator) {
        ensureNonNull("Generator", generator);
        return () -> {
            try {
                return generator.call();
            } catch (RuntimeException e) {
                throw e;
            } catch (Exception e) {
                throw new BackingStoreException(e);
            }
        };
    }

    /**
     * Describes a stream on which a {@link Stream#map(Function) mapping operation} has been set. It serves to delegate
     * optimizable operation to underlying sql stream (which could be an indirect parent).
     *
     * @param <I> Type of object received as input of mapping operation.
     * @param <O> Return type of mapping operation.
     */
    private static final class MappedStream<I, O> extends StreamDecoration<O> {
        private Function<? super I, ? extends O> mapper;
        private Stream<I> source;

        private MappedStream(Function<? super I, ? extends O> mapper, Stream<I> source) {
            this.mapper = mapper;
            this.source = source;
        }

        @Override
        public Stream<O> peek(Consumer<? super O> action) {
            mapper = concatenate(mapper, action);
            return this;
        }

        @Override
        public Stream<O> distinct() {
            source = source.distinct();
            return this;
        }

        @Override
        public Stream<O> limit(long maxSize) {
            source = source.limit(maxSize);
            return this;
        }

        @Override
        public Stream<O> skip(long n) {
            source = source.skip(n);
            return this;
        }

        @Override
        public long count() {
            return source.count();
        }

        @Override
        public boolean isParallel() {
            return source.isParallel();
        }

        @Override
        public Stream<O> sequential() {
            source = source.sequential();
            return this;
        }

        @Override
        public Stream<O> parallel() {
            source = source.parallel();
            return this;
        }

        @Override
        public <R> Stream<R> map(Function<? super O, ? extends R> mapper) {
            return new MappedStream<>(this.mapper.andThen(mapper), source);
        }

        @Override
        public DoubleStream mapToDouble(ToDoubleFunction<? super O> mapper) {
            return new ToDoubleStream<I>(source, i -> mapper.applyAsDouble(this.mapper.apply(i)));
        }

        @Override
        protected Stream<O> createDecoratedStream() {
            // Break possible infinite loop by sinking source content through its spliterator (terminal op).
            final Stream<I> sink = StreamSupport.stream(source.spliterator(), source.isParallel());
            sink.onClose(source::close);
            return sink.map(mapper);
        }
    }

    /**
     * Same purpose as {@link MappedStream}, but specialized for {@link Stream#mapToDouble(ToDoubleFunction) double mapping}
     * operations.
     *
     * @param <T> Type of objects contained in source stream (before double mapping).
     */
    private static final class ToDoubleStream<T> extends DoubleStreamDecoration {

        Stream<T> source;
        ToDoubleFunction<T> toDouble;

        private ToDoubleStream(Stream<T> source, ToDoubleFunction<T> toDouble) {
            this.source = source;
            this.toDouble = toDouble;
        }

        @Override
        public DoubleStream peek(DoubleConsumer action) {
            final ToDoubleFunction<T> toDoubleFixedRef = toDouble;
            toDouble = t -> {
                final double value = toDoubleFixedRef.applyAsDouble(t);
                action.accept(value);
                return value;
            };
            return this;
        }

        @Override
        public DoubleStream map(DoubleUnaryOperator mapper) {
            return new ToDoubleStream<T>(source, t -> mapper.applyAsDouble(toDouble.applyAsDouble(t)));
        }

        @Override
        public <U> Stream<U> mapToObj(DoubleFunction<? extends U> mapper) {
            return new MappedStream<T, U>(t -> mapper.apply(toDouble.applyAsDouble(t)), source);
        }

        @Override
        public DoubleStream distinct() {
            source = source.distinct();
            return this;
        }

        @Override
        public DoubleStream limit(long maxSize) {
            source = source.limit(maxSize);
            return this;
        }

        @Override
        public DoubleStream skip(long n) {
            source = source.skip(n);
            return this;
        }

        @Override
        public long count() {
            return source.count();
        }

        @Override
        public Stream<Double> boxed() {
            return new MappedStream<>(t -> Double.valueOf(toDouble.applyAsDouble(t)), source);
        }

        @Override
        public boolean isParallel() {
            return source.isParallel();
        }

        @Override
        public DoubleStream sequential() {
            source = source.sequential();
            return this;
        }

        @Override
        public DoubleStream parallel() {
            source = source.parallel();
            return this;
        }

        @Override
        protected DoubleStream createDecoratedStream() {
            // Break possible cycle by sinking source content through its spliterator (terminal op).
            final Stream<T> sink = StreamSupport.stream(source.spliterator(), source.isParallel());
            sink.onClose(source::close);
            return sink.mapToDouble(toDouble);
        }
    }

    private static <I, O> Function<? super I, ? extends O> concatenate(final Function<? super I, ? extends O> function, final Consumer<? super O> consumer) {
        return i -> {
            final O o = function.apply(i);
            consumer.accept(o);
            return o;
        };
    }
}
