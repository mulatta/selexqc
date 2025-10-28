{
  perSystem = {
    self',
    pkgs,
    lib,
    ...
  }: let
    inherit (pkgs.stdenv) isLinux;
    linkerArgs = lib.optionalString isLinux "-C link-arg=-fuse-ld=mold";
  in {
    packages = {
      slxqc = pkgs.rustPlatform.buildRustPackage {
        pname = "slxqc";
        version = "0.1.0";

        src = lib.cleanSourceWith {
          src = ../.;
          filter = path: _type: let
            baseName = baseNameOf path;
          in
            !(lib.hasSuffix ".nix" baseName)
            && baseName != "target"
            && baseName != "result"
            && baseName != ".git"
            && baseName != ".github";
        };

        cargoLock = {
          lockFile = ../Cargo.lock;
        };

        nativeBuildInputs = with pkgs;
          [
            pkg-config
          ]
          ++ lib.optionals isLinux [mold];

        RUSTFLAGS = lib.concatStringsSep " " (
          lib.filter (x: x != "") [
            linkerArgs
          ]
        );

        # Parallel jobs
        buildPhase = ''
          runHook preBuild
          export CARGO_BUILD_JOBS=$NIX_BUILD_CORES
          cargo build --release
          runHook postBuild
        '';

        installPhase = ''
          runHook preInstall
          install -Dm755 target/release/slxqc $out/bin/slxqc
          runHook postInstall
        '';

        # Tests
        checkPhase = ''
          runHook preCheck
          cargo test --release
          runHook postCheck
        '';

        doCheck = true;

        meta = with lib; {
          description = "High-performance parallel FASTA/FASTQ sequence counter";
          longDescription = ''
            slxqc is a fast, parallel quality control and filtering tool for RNA Capture-SELEX
            NGS libraries. It validates library structure, detects constant regions, and filters
            sequences based on configurable criteria with support for FASTA, FASTQ, and
            compressed formats. Features include:

            - Fast parallel processing with zero-copy parsing
            - Structure-aware validation (upstream/downstream pairs)
            - Flexible AND/OR validation logic
            - Sequence filtering with multiple output formats
            - MultiQC integration for pipeline workflows
            - Comprehensive reporting (TXT/CSV/JSON)
          '';
          homepage = "https://github.com/mulatta/slxqc";
          license = licenses.mit;
          maintainers = [];
          platforms = platforms.linux ++ platforms.darwin;
          mainProgram = "slxqc";
        };
      };

      default = self'.packages.slxqc;
    };
  };
}
