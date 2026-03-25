{inputs, ...}: {
  imports = [inputs.treefmt-nix.flakeModule];
  perSystem = {
    treefmt = {
      projectRootFile = ".git/config";

      programs = {
        # Nix formatters & linters
        alejandra.enable = true;
        deadnix.enable = true;
        statix.enable = true;

        # Rust formatters & linters
        rustfmt.enable = true;

        # Shell formatters & linter
        shellcheck.enable = true;
        shfmt.enable = true;

        # Other formatters
        keep-sorted.enable = true;
        yamlfmt.enable = true;
        taplo.enable = true;
      };

      settings.global.excludes = [
        "*.lock"
        ".gitignore"
        "results/**"
      ];
    };
  };
}
