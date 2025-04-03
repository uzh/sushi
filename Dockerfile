FROM debian:bookworm-slim

# Install required packages
RUN apt-get update && apt-get install -y \
    build-essential \
    libssl-dev \
    zlib1g-dev \
    libreadline-dev \
    libyaml-dev \
    libgdbm-dev \
    libncurses5-dev \
    libffi-dev \
    libgmp-dev \
    curl \
    git \
    default-libmysqlclient-dev \
    libsqlite3-dev \
    && rm -rf /var/lib/apt/lists/*

# Download and compile Ruby
RUN cd /tmp && \
    curl -O https://cache.ruby-lang.org/pub/ruby/3.3/ruby-3.3.7.tar.gz && \
    tar -xzf ruby-3.3.7.tar.gz && \
    cd ruby-3.3.7 && \
    ./configure --prefix=/usr/local --enable-shared --disable-install-doc && \
    make -j $(nproc) && \
    make install && \
    cd / && \
    rm -rf /tmp/ruby-3.3.7 /tmp/ruby-3.3.7.tar.gz

# Install bundler and rubygems-bundler
RUN gem install bundler rubygems-bundler && \
    gem regenerate_binstubs

# Create a non-root user
RUN useradd -m sushi

# Set working directory
WORKDIR /app

# Copy entrypoint script
COPY master/entrypoint.sh /usr/bin/
RUN chmod +x /usr/bin/entrypoint.sh

COPY master/ ./
RUN chown -R sushi:sushi /app && \
    su sushi -c "bundle config set --local path 'vendor/bundle' && bundle install"

# Expose port 3000
EXPOSE 3000

# Set the command to start the Rails server
CMD ["bundle", "exec", "rails", "server", "-b", "0.0.0.0"] 

