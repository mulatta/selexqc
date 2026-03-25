{inputs, ...}: {
  perSystem = {
    self',
    pkgs,
    lib,
    ...
  }: let
    rustToolchain = pkgs.rust-bin.stable.latest.default;
    craneLib = (inputs.crane.mkLib pkgs).overrideToolchain rustToolchain;

    src = lib.cleanSourceWith {
      src = ../.;
      filter = path: type:
        (craneLib.filterCargoSources path type)
        || (lib.hasInfix "/tests/fixtures/" path && lib.hasSuffix ".fastq" (baseNameOf path));
    };

    commonArgs = {
      inherit src;
      pname = "selexqc";
      version = "0.2.0";
      nativeBuildInputs = with pkgs; [pkg-config];
    };

    cargoArtifacts = craneLib.buildDepsOnly commonArgs;
  in {
    packages = {
      selexqc = craneLib.buildPackage (commonArgs
        // {
          inherit cargoArtifacts;
          doCheck = pkgs.stdenv.buildPlatform.canExecute pkgs.stdenv.hostPlatform;
          meta = with lib; {
            description = "SELEX analysis toolkit: construct QC, sequence counting, and enrichment analysis";
            license = licenses.mit;
            platforms = ["x86_64-linux" "aarch64-linux" "aarch64-darwin"];
            mainProgram = "selexqc";
          };
        });

      default = self'.packages.selexqc;
    };
  };
}
