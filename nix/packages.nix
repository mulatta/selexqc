{inputs, ...}: {
  perSystem = {
    lib,
    pkgs,
    self',
    ...
  }: {
    packages = {
      selexqc = pkgs.rustPlatform.buildRustPackage {
        pname = "selexqc";
        version = "0.1.0";

        src = inputs.gitignore.lib.gitignoreSource ../.;

        cargoLock.lockFile = ../Cargo.lock;

        nativeBuildInputs = with pkgs; [
          pkg-config
        ];

        meta = with lib; {
          description = "High-performance quality control tool for RNA Capture-SELEX NGS libraries";
          homepage = "https://github.com/mulatta/selexqc";
          license = licenses.mit;
          maintainers = [];
          platforms = platforms.linux ++ platforms.darwin;
          mainProgram = "selexqc";
        };
      };

      default = self'.packages.selexqc;
    };
  };
}
