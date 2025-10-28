{self, ...}: {
  flake.overlays.default = final: _prev: {
    inherit (self.packages.${final.system}) selexqc;
  };
}
