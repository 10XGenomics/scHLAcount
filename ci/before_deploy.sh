#!/bin/bash

# package the build artifacts

set -ex

. "$(dirname $0)/utils.sh"

# Generate artifacts for release
mk_artifacts() {
    cargo build --target "$TARGET" --release
}

mk_tarball() {
    # When cross-compiling, use the right `strip` tool on the binary.
    local gcc_prefix="$(gcc_prefix)"
    # Create a temporary dir that contains our staging area.
    # $tmpdir/$name is what eventually ends up as the deployed archive.
    local tmpdir="$(mktemp -d)"
    local name="${PROJECT_NAME}-${TRAVIS_TAG}-${TARGET}"
    local staging="$tmpdir/$name"
    mkdir -p "$staging"
    # The deployment directory is where the final archive will reside.
    # This path is known by the .travis.yml configuration.
    local out_dir="$(pwd)/deployment"
    mkdir -p "$out_dir"
    # Find the correct (most recent) Cargo "out" directory. The out directory
    # contains shell completion files and the man page.
    local cargo_out_dir="$(cargo_out_dir "target/$TARGET")"

    # Copy the scHLAcount binary and strip it.
    cp "target/$TARGET/release/sc_hla_count" "$staging/sc_hla_count"
    #"${gcc_prefix}strip" "$staging/sc_hla_count"
    # Copy the licenses and README.
    cp {README.md,LICENSE} "$staging/"

    (cd "$tmpdir" && tar czf "$out_dir/$name.tar.gz" "$name")
    rm -rf "$tmpdir"
}

main() {
    mk_artifacts
    mk_tarball
}

main
